#!/usr/bin/perl

##############################
#
# written by Mark L. Farman
#
# Purpose: Mask genome sequences, perform pairwise blast searches, call SNPs in uniquely aligned regions 
#
# Note: Uses "fastest" SNPcaller algorithm
#
# Note: works in incremental fasion: only newly added genomes are processed when main program is run
#
# Usage: perl iSNPcaller.pl <Project_name> [hit run when prompted]
#
##############################

use strict;

use warnings;
use File::Find::Rule;
use File::Copy;
use File::Path;
use Repeatmaskerfast_rr;
use Snpcallerfaster_rerun;
use Megamatrix;

##############################

# Declare global variables;

my $Directory;
my @To_be_BLASTed;
my %NewGenomesHash;
my %BlastTempHash;
my $Border = "@" x 40;
my $Separator = '-' x 40;
my $Dialog_bar = "\n"."-" x 40;

##############################

# Read command line arguments

if(@ARGV) {
  if(@ARGV != 2) {
    die "\n$Dialog_bar\n".
	"Usage:- 1st run for specified project: perl SNPcaller.pl\n\n".
	"Subsequent runs: perl SNPcaller.pl <PROJECT_NAME> run\n\n"
  }

  my($DirName, $Run) = @ARGV;

  if($Run =~ /[RUN|run]/) {
    chdir("$DirName") || die "\n$Dialog_bar\n".
			     "CANNOT FIND A PROJECT OF THAT NAME IN THE WORKING DIRECTORY\n".
			     "Check project name and re-start analysis\n\n";
    & RUN_PROGRAM
  }

  else {
    & INITIAL_DIALOG
  }

}

else {
  & INITIAL_DIALOG;
}


##############################

## Initiate user dialog


sub INITIAL_DIALOG {

    print "\n$Border\n\nIS THIS THE FIRST TIME YOU HAVE RUN THE SNPCALLER PROGRAM FOR THIS PROJECT?\nNOTE: Select NO if you need to complete a failed run or if you want to add more genomes to analyze\n\n\t";
    my $Response = <STDIN>;
    chomp($Response);
    if($Response =~ /^[Y|y]/) {
        & CREATE_PROJECT
    }

    elsif($Response =~ /^[N|n]/) {
       print "\n$Dialog_bar\n".
	     "ENTER THE NAME OF THE PROJECT TO WHICH THE NEW GENOMES BELONG\n\n\t";
       my $Answer = <STDIN>;
       chomp($Answer);
       PROJECT_DIRECTORY($Answer);
    }

    else {
        while($Response !~ /^[YyNn]/) {
            print "\n$Dialog_bar\n".
	 	  "INVALID ENTRY.  Please try again\n\n";
            & INITIAL_DIALOG;
        }
    }
}


## Provide a name for the PROJECT DIRECTORY that will contain all input files and results

sub CREATE_PROJECT {

    print "$Dialog_bar\n".
	  "ENTER NAME OF PROJECT (alphanumerics only, with underscores/hyphens for spaces):\n\n\t";
    my $Directory = <STDIN>;
    chomp($Directory);
    if($Directory !~ /\s|\.|\#|\:|\$/) {
        if(-d $Directory) {
            DIRECTORY_ERROR($Directory)
        }
        else {
            MAKE_SUB_DIRECTORIES($Directory)
        }
    }

    else {
        PROJECT_ERROR($Directory)
    }
}


sub DIRECTORY_ERROR {

    my $Directory = $_[0];
    while(-d $Directory) {
        print "\n$Dialog_bar\n".
              "THERE IS ALREADY A DIRECTORY WITH THE CHOSEN PROJECT NAME - $Directory. PLEASE TRY AGAIN\n\n\t";
        $Directory = <STDIN>;
        chomp($Directory);
        while($Directory =~ /(\W)|(\s)|(\.)|(\#)|(:)|(\$)/) {
            print "\n$Dialog_bar\n".
                  "PROJECT NAME CONTACTS AN ILLEGAL CHARACTER - use only alphanumerics (A-Z, a-z, 0-9) and underscores. PLEASE TRY AGAIN:\n\n\t";
            $Directory = <STDIN>;
            chomp($Directory);
        }
    }
    MAKE_SUB_DIRECTORIES($Directory)
}

## Check for illegal characters in PROJECT NAME 

sub PROJECT_ERROR {
    
    my $Directory = $_[0];

    while($Directory =~ /(\W)|(\s)|(\.)|(\#)|(:)|(\$)/) {
        print "\n$Dialog_bar\n".
	      "PROJECT NAME CONTACTS AN ILLEGAL CHARACTER - use only alphanumerics (A-Z, a-z, 0-9) and underscores\n\n";
        & PROJECT_NAME
    }
}


## Check for existence of the specified PROJECT DIRECTORY in current working directory and define this PROJECT DIRECTORY as the new working directory

sub PROJECT_DIRECTORY {

  my $Directory = $_[0];

  if(-d $Directory) { 	
    chdir($Directory);
    my $WorkingDir = qx{pwd};
    print "\n$Dialog_bar\n".
	  "SUCCESSFULLY CHANGED DIRECTORY TO: $WorkingDir\n".
	  "If necessary, copy new genome .fasta files into the $Directory/NEW_GENOMES/".
	  "directory and then type: run\n\n\t";
    & GET_RUN_INPUT;
  }

  else {
    my $WorkingDir = `pwd`;
    chomp($WorkingDir);
    print "\n$Dialog_bar\n".
	  "CANNOT FIND A PROJECT OF THAT NAME IN THE WORKING DIRECTORY: $WorkingDir.\nChoose from the following options:\n".
	  "[1] Re-enter project name\n".
	  "[2] Quit program and change working directory before re-starting\n\n";
    my $Answer = <STDIN>;
    chomp($Answer);

    while($Answer !~ /1|2/) {
      print "\n$Dialog_bar\n".
	    "INVALID ENTRY. Select option 1 or 2\n\n\t";
      $Answer = <STDIN>;
      chomp($Answer);
    }

    if($Answer == 1) {
      print "\n$Dialog_bar\n".
	    "RE-ENTER PROJECT NAME\n\n\t";
      $Directory = <STDIN>;
      chomp($Directory);
      PROJECT_DIRECTORY($Directory)
    }

    exit if($Answer == 2)
  }    
}


# Create directories for inputs and results

sub MAKE_SUB_DIRECTORIES {

        my $Directory = $_[0];
        mkpath( "$Directory/NEW_GENOMES");	
        mkpath( "$Directory/NEW_GENOMES/DBs");
        mkpath( "$Directory/NEW_GENOMES/BLAST_RESULTS");
        mkpath( "$Directory/PROCESSED/FASTA/MASKED");
        mkpath( "$Directory/PROCESSED/DBs");
        mkpath( "$Directory/PROCESSED/BLAST_RESULTS");
        mkpath( "$Directory/LOGS");
	mkpath( "$Directory/SNP_COUNTS/SNP_FILES");
	print 	"$Dialog_bar\n".'DIRECTORIES USED BY THE PROGRAM HAVE BEEN CREATED UNDER THE TOP (PROJECT) DIRECTORY NAMED: '.'"'."$Directory".'"'."\n".
		"Copy your genome .fasta files into the $Directory/NEW_GENOMES/ directory and then ".
		"type: run\n\n\t";
        chdir("$Directory"); 
        GET_RUN_INPUT($Directory);
}

sub GET_RUN_INPUT {
    my $Directory = $_[0];
    my $Answer = <STDIN>;
    chomp($Answer);

    if($Answer =~ /^[R|r]/) {
	print "$Dialog_bar\n".'NOTE: NEXT TIME YOU RUN THIS PROGRAM, YOU CAN SIMPLY ADD YOUR NEW GENOMES TO THE '.'"PROJECT/NEW_GENOMES"'.' DIRECTORY AND BYPASS THE INITIAL DIALOG BY TYPING: '.
              "perl IMPB.pl $Directory run\n\n\t";
        & RUN_PROGRAM
    }

    else {
        while($Answer !~ /^[R|r]/) {
	    print "$Dialog_bar\nINVALID ENTRY. Please try again\n\n\t";
	    & GET_RUN_INPUT;
	}
    }
}

sub RUN_PROGRAM {

    & MASK_GENOMES;

    & CREATE_DBs;					# create blast databases

    & READ_BLASTTEMPFILE;

    & BLAST_NEW_GENOMES;					# Blast the new genomes against one another

    & BLAST_NEW_V_OLD;						# Blast the new genomes against the ones that have already been processed
    
    my $DateTime = Snpcallerfaster_rerun::SNPs();

    & MOVE_FASTAs_DBs_TO_PROCESSED_DIR;  				# Move new .fasta files to FASTA folder in PROCESSED directory

    & MOVE_BLAST_RESULTS_TO_PROCESSED_DIR;			# move recently-generated BLAST results from NEW_GENOME
                                         	                # directory into PROCESSED directory
    Megamatrix::MATRIX($DateTime);

    print "$Separator\nANALYSIS SUCCESSFULLY COMPLETED\n"

}
 

sub MASK_GENOMES {

    Repeatmaskerfast_rr::MASK();	# creates masked genomes in same directory as original fasta files
    my @Masked = glob("NEW_GENOMES/*.fasta NEW_GENOMES/*.fsa");
}

sub READ_BLASTDBTEMPFILE {

  my $BlastDBFile;
  my %BlastDBTempHash;

  if(-f "BLASTDB_tempfile") {
    open(BLASTDBTEMPFILE, "BLASTDB_tempfile");
    my $LineCount = 0;

    while(my $L = <BLASTDBTEMPFILE>) {
      chomp($L);
      $BlastDBFile = $L;
      $BlastDBTempHash{$BlastDBFile} = 1;
    }

    close BLASTDBTEMPFILE;
  }

  open(BLASTDBTEMPFILE, '>>', "BLASTDB_tempfile");

}

## Create blastable DBs
  
sub CREATE_DBs {
    my $pwd = `pwd`;
    my @New;
    my $New_genome;

# Create a BLASTable database for each new masked genome

    opendir(NEW, "NEW_GENOMES") || die "Cannot find NEW_GENOMES dir: $pwd\n";
    @New = readdir(NEW);

    foreach $New_genome (@New) {
        $NewGenomesHash{$New_genome} = 1;
        push @To_be_BLASTed, $New_genome if $New_genome =~ /masked\.fasta$|masked\.fsa$/;               # add genome to list of "new" sequences to be BLASTed against "old" genomes
	system("makeblastdb -in NEW_GENOMES/$New_genome -dbtype nucl 2>/dev/null >exceptions.txt");	# create database
    }
}


## Check for tempfile from a previous interrupted run 

sub READ_BLASTTEMPFILE {

  my $BlastFile;

  if(-f "BLAST_tempfile") {
    open(BLASTTEMPFILE, "BLAST_tempfile");
    my $LineCount = 0;
    
    while(my $L = <BLASTTEMPFILE>) {
      chomp($L);
      $BlastFile = $L;
      my($Q, $S) = split(/\t/, $BlastFile);
      $BlastTempHash{$Q}{$S} = 1;
    }

    close BLASTTEMPFILE;
  }

  open(BLASTTEMPFILE, '>>', "BLAST_tempfile");

}

## Blast new genomes against one another in pairwise fashion

sub BLAST_NEW_GENOMES {
	
	my %New_genome;
 	my $Num_new_genomes = @To_be_BLASTed;

        if($Num_new_genomes > 1) {	
		print "\n$Dialog_bar\n\n".
			"BLASTing new genomes against one another\n\n".
			"$Dialog_bar\n\n";
		my $i = 0;
	        my $j = 0;
		for($i = 0; $i <= $Num_new_genomes-2; $i++) {
			$New_genome{$i} = 1;
	                for ($j = $i + 1; $j <= $Num_new_genomes-1; $j++) {
				my ($Q, $DB) = ($To_be_BLASTed[$i], $To_be_BLASTed[$j]);
		                $Q =~ s/\.fasta//;
		                $DB =~ s/\.fasta//;
                                next if(exists($BlastTempHash{$Q}{$DB}));
				$New_genome{$j} = 1;
				print "BLASTing query $To_be_BLASTed[$i] against DB $To_be_BLASTed[$j]\n";
				
				if(system("blastn -db NEW_GENOMES/$To_be_BLASTed[$j] -query NEW_GENOMES/$To_be_BLASTed[$i]".
				 " -out NEW_GENOMES/BLAST_RESULTS/$Q.$DB.BLAST".
				 " -max_target_seqs 20000 -evalue 1e-200 -outfmt '6 qseqid sseqid qstart qend sstart send btop'".
                 	        " 2>/dev/null >>exceptions.txt")) {exit};

				print BLASTTEMPFILE "$Q\t$DB\n"
			}
		}
	}
}

## Reads PROCESSED/FASTA directory, reads To_be_BLASTed list
## performs pairwise blasts between "new" and "already processed" genomes

sub BLAST_NEW_V_OLD {
        my $FsaFileCount = 0;
	opendir(PROCESSED, "PROCESSED/FASTA/MASKED");
	my @Processed = readdir(PROCESSED);
	my $Fasta_files = "no";
	foreach my $Processed_genome (@Processed) {
		if(exists($NewGenomesHash{$Processed_genome})) {
			next
		}
		elsif($Processed_genome =~ /\.fasta$|\.fsa$/) {
			$FsaFileCount ++;

			if($FsaFileCount == 1) {
				print 	"\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n\n".
					"BLASTing new genome(s) against old ones\n\n";
			}

                        my $DB = $Processed_genome;
	                $DB =~ s/\.fasta$|\.fsa$//;

			foreach my $New_genome (@To_be_BLASTed) {
	           		my $Query = $New_genome;
	                        $Query =~ s/\.fasta$|\.fsa$//;
				next if($New_genome eq $Processed_genome);				
				next if(exists($BlastTempHash{$New_genome}{$Processed_genome}));

				print "BLASTing query $New_genome against DB $Processed_genome\n";

	                	if(system("blastn -db PROCESSED/DBs/$Processed_genome -query NEW_GENOMES/$New_genome".
					" -out 'NEW_GENOMES/BLAST_RESULTS/$Query.$DB.BLAST'".
					" -max_target_seqs 20000 -evalue 1e-200 -outfmt '6 qseqid sseqid qstart qend sstart send btop'".
                               		" 2>/dev/null >>exceptions.txt")) {exit};

           	               print BLASTTEMPFILE "$Query\t$DB\n"

            		}

	        }

    	}

	$FsaFileCount = 0

}


# move masked .fasta files into PROCESSED/FASTA directory

sub MOVE_FASTAs_DBs_TO_PROCESSED_DIR {

    my @Fastas_and_DBs = glob("NEW_GENOMES/*fasta* NEW_GENOMES/*fsa*");
    foreach my $FilePath (@Fastas_and_DBs) {       
      (my $File = $FilePath) =~ s/NEW_GENOMES\///;
      copy("NEW_GENOMES/$File", "PROCESSED/FASTA/MASKED/$File") if $File =~ /masked\.fasta$|masked\.fsa$/;
      move("NEW_GENOMES/$File", "PROCESSED/DBs/$File") if $File =~ /masked/;
      move("NEW_GENOMES/$File", "PROCESSED/FASTA/$File") if $File =~ /fasta$|fsa$/;
      unlink("$FilePath")
    } 
}

# move BLAST results files into PROCESSED/BLAST_RESULTS directory

sub MOVE_BLAST_RESULTS_TO_PROCESSED_DIR {
  my @Blast_results = glob("NEW_GENOMES/BLAST_RESULTS/*BLAST");

  foreach my $BlastPath (@Blast_results) {
    (my $BlastFile = $BlastPath) =~ s/NEW_GENOMES\/BLAST_RESULTS\///;
    move("$BlastPath", "PROCESSED/BLAST_RESULTS/$BlastFile")
  }
}
