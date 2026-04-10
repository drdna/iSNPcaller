#!/usr/bin/perl

##############################
#
# written by Mark L. Farman
#
# Purpose: Mask genome sequences, perform reciprocal pairwise blast searches against each other 
#
# Note: works in incremental fasion: only newly added genomes are processed when main program is run
#
# Usage: perl iSNPcaller_MT.pl <Project_name> [follow runtime instructions]
#
##############################

# use strict;

use warnings;

#use File::Find::Rule;

use File::Copy;

use File::Path;

use constant MAX_PROCESSES => 48;

use Parallel::ForkManager;

use RepeatmaskerfastMT;

use Snpcallerfast;

use SimpleFastaHeaders;

##############################

# Declare global variables;

my @To_be_BLASTed;

my %NewGenomesHash;

##############################


my $GenomesDialog =     "There are currently no genomes to be analyzed in the GENOMES directory.".
                        "\nOpen another terminal window, copy your genome .fasta files into the ".
                        "directory and then type: run\n\n";

my $GenomeDirDialog =   "Cannot find expected GENOMES directory ".
                        "- either this is not an iSNPcaller project ".
                        ", or the required directory has been moved/deleted. Quitting program.\n";

my $NoProjectDir =      "Cannot find a directory with the specified project name:\n";

my $NoProjectDialog =   "Cannot find a project with the specified name in the current directory. ".
                        "Type another name or type quit.\n";


# Read command line argument(s), if any

if(@ARGV) {

  my $Project = $ARGV[0];
  $Project =~ s/\W$//;

  if(-d "$Project/GENOMES") {
    RUN_PROJECT($Project)
  }
  else {
    my $ValidProject = PROJECT_CHECK($Project);
    MAKE_DIRECTORIES($ValidProject)
  }
}

else {
  & INITIAL_DIALOG
}


sub INITIAL_DIALOG {

  print "Is this the first time you have run the Incremental_Multipairwise_Blast script?\n\n";

  my $Answer = <STDIN>;

  chomp($Answer);

  if($Answer =~ /^[Y|y]/) {

    print "Please provide a name for the project (no spaces allowed)\n";

    my $Project = <STDIN>;

    $Project =~ s/\W$//;

    my $ValidProject = PROJECT_CHECK($Project);

    MAKE_DIRECTORIES($ValidProject)

  }

  elsif($Answer =~ /^[N|n]/) {

    print "Please type the project name\n";

    my $Project = <STDIN>;

    chomp($Project);

    $Project =~ s/\W$//;

    RUN_PROJECT($Project);
   
  }

}


# Create directories for inputs and results

sub MAKE_DIRECTORIES {

  my $Project = $_[0];

  print "New project name is: $Project\n";

  mkpath( $Project."/GENOMES/PROCESSED") || die "Cannot create project - check that there is not already a directory with the same name in the current working directory\n";

  mkpath( $Project."/NEW_GENOMES/FASTA");

  mkpath( $Project."/NEW_GENOMES/BLAST_RESULTS");

  mkpath( $Project."/PROCESSED/FASTA");

  mkpath( $Project."/PROCESSED/BLAST_RESULTS");

  print "Directories used by the program have been created under the $Project top directory. ".
        "Open another terminal window (or use the GUI interface) and copy your genome ".
        ".fasta files into the $Project/GENOMES/ directory and then ".
        "type: run\n\n";;

  chdir($Project);

  & GET_RUN_INPUT;

}


sub PROJECT_CHECK {

  my $Project = $_[0];

  my $response = $Project;

  while(-d $Project || -f $Project) {

    print "Re-type a different project name that does not match a folder/file in you current working directory\n";

    $response = <STDIN>;

    chomp($response);

      $response =~ s/\W+$//;

    }

    while($Project =~ /[^A-Za-z0-9-_]/) {

      print "Re-type a different project name that contains only alphanumerics, hyphens, or underscores\n";

      $response = <STDIN>;

      chomp($response);

      $response =~ s/\W+$//;

    }

  return($response)

}



# Start analyzing genomes in project directory

sub RUN_PROJECT {

  my $Project = $_[0];

  # Look for directory matching provided project ID

  if(-d $Project) {

    chdir($Project);

    if(-d "GENOMES") {

      opendir(GENOMES, "GENOMES") || die $GenomeDirDialog;

      @FastaFiles = readdir(GENOMES);

      if( grep( /fasta$|fa$|fna$|fsa$/, @FastaFiles)) {

        & RUN_PROGRAM

      }

      else {

	print $GenomesDialog;

        my $response = <STDIN>;

        if($response =~ /^R/i) {

          & RUN_PROGRAM

        }

      }

    }

  }

  else {

    print $NoProjectDialog;

    my $response = <STDIN>;

    chomp($response);

    die "Quitting program.\n" if $response =~ /^Q/i;

  }

}

##############################

sub GET_RUN_INPUT {

  my $Answer = <STDIN>;

  chomp($Answer);

  if($Answer =~ /^[R|r]/) {

    & RUN_PROGRAM

  }

  else {

    while($Answer !~ /^[R|r]/) {

    print "Unrecognized entry.  Please try again\n\n";

    & GET_RUN_INPUT;

  }

}

sub RUN_PROGRAM {

  $FastaDir = 'GENOMES';

  SimpleFastaHeaders::rehead($FastaDir);

  & MASK_GENOMES;

  & ADD_NEW_GENOMES;

  & BLAST_NEW_GENOMES;                    # Blast the new genomes against one another

  & BLAST_NEW_V_OLD;		               # Blast the new genomes against the ones that have already been processed

  & MOVE_NEW_FASTAs_TO_PROCESSED_DIR;     # Move new .fasta files to FASTA folder in PROCESSED directory

  print   "\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n\n".
            "Running the SNP caller\n";

  Snpcallerfast::SNPs;

  & MOVE_OLD_BLAST_RESULTS;               # move old BLAST results from NEW_GENOME
                                                # directory into PROCESSED directory
}

sub MASK_GENOMES {

  RepeatmaskerfastMT::MASK();

  opendir(GENOMES, "GENOMES");

  my @FileList = readdir(GENOMES);

  foreach my $Masked_file (@FileList) {

    if($Masked_file =~ /masked.+fasta$/) {

      move("GENOMES/$Masked_file", "NEW_GENOMES/FASTA/$Masked_file")

    }

    elsif($Masked_file =~ /fasta$/) {

      move("GENOMES/$Masked_file", "GENOMES/PROCESSED")

    }

  }

}


sub ADD_NEW_GENOMES {

  my @New;

  my $New_genome;

  # open directory containing the new genome sequences

  opendir(NEW, "NEW_GENOMES/FASTA");

  @New = readdir(NEW);

  foreach $New_genome (@New) {

    if($New_genome =~ /\.fasta$|\.fsa$/) {
			
      $NewGenomesHash{$New_genome} = 1;

      push @To_be_BLASTed, $New_genome;               # add genome to the "To_be_BLASTed" list
			
    }

  }

}


sub BLAST_NEW_GENOMES {
	
  my %New_genome;

# Perform reciprocal pairwise BLASTs between all new genomes
	
# Takes each genome in To_be_BLASTed list and BLASTs against all others

  my $Num_new_genomes = @To_be_BLASTed;

  if($Num_new_genomes > 1) {
	
    print "\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n\n".
          "BLASTing new genomes against one another\n\n";
	
    my $i = 0;

    my $j = 0;

    for($i = 0; $i <= $Num_new_genomes-2; $i++) {

      $New_genome{$i} = 1;

      my $pm = Parallel::ForkManager->new(MAX_PROCESSES);

      for ($j = $i + 1; $j <= $Num_new_genomes-1; $j++) {

        my $pid = $pm->start and next;

        my ($Q, $SUB) = ($To_be_BLASTed[$i], $To_be_BLASTed[$j]);

        $Q =~ s/\.fasta//;

        $SUB =~ s/\.fasta//;

        $New_genome{$j} = 1;
				
	print "BLASTing query $To_be_BLASTed[$i] against subject $To_be_BLASTed[$j]\n";
				
        exec("blastn -subject NEW_GENOMES/FASTA/$To_be_BLASTed[$j] -query NEW_GENOMES/FASTA/$To_be_BLASTed[$i]".
				 " -out NEW_GENOMES/BLAST_RESULTS/$Q.$SUB.BLAST".
				 " -max_target_seqs 20000 -evalue 1e-20 -outfmt '6 qseqid sseqid qstart qend sstart send btop'".
                 	        " 2>/dev/null >>exceptions.txt") || die $!;
	
        $pm -> finish		
      }

      $pm->wait_all_children()

    }

  }

}


sub MOVE_NEW_FASTAs_TO_PROCESSED_DIR {
	
  # move new genome .fasta files to the PROCESSED directory
	
  foreach my $Genome (@To_be_BLASTed) {

    move("NEW_GENOMES/FASTA/$Genome", "PROCESSED/FASTA/$Genome")

  }

}

sub BLAST_NEW_V_OLD {

  print 	"\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n\n".
		"BLASTing new genome(s) against old ones\n\n";

# perform pairwise BLASTs between new and old genomes

# Reads PROCESSED/FASTA directory, reads To_be_BLASTed list
# performs blasts between "new" and "already processed" genomes

  opendir(PROCESSED, "PROCESSED/FASTA");
	
  my @Processed = readdir(PROCESSED);

  my $Fasta_files = "no";
		
  foreach my $Processed_genome (@Processed) {
			
    if(exists($NewGenomesHash{$Processed_genome})) {

      next

    }

    elsif($Processed_genome =~ /\.fasta$/) {

      $Fasta_files = "yes";

      my $SUB = $Processed_genome;

      $SUB =~ s/\.fasta$//;
			
        my $pm = Parallel::ForkManager->new(MAX_PROCESSES);

	foreach my $New_genome (@To_be_BLASTed) {

          my $pid = $pm->start and next;

          my $Query = $New_genome;

	  $Query =~ s/\.fasta$//;

	  $pm->finish if $New_genome eq $Processed_genome;
					
	  print "BLASTing query $New_genome against subject $Processed_genome\n";

	  system("blastn -subject PROCESSED/FASTA/$Processed_genome -query NEW_GENOMES/FASTA/$New_genome".
	         " -out NEW_GENOMES/BLAST_RESULTS/$Query.$SUB.BLAST".
                 " -max_target_seqs 20000 -evalue 1e-20 -outfmt '6 qseqid sseqid qstart qend sstart send btop'".
                 " 2>/dev/null >>exceptions.txt");
        
          $pm -> finish

        }

        $pm->wait_all_children()

      }

    }

    if($Fasta_files eq "no") {

      print "\nNo old genomes to analyze yet\n"

     }

  } 

}

sub MOVE_OLD_BLAST_RESULTS {

  system('mv NEW_GENOMES/BLAST_RESULTS/* PROCESSED/BLAST_RESULTS/')

}
