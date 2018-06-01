package Repeatmaskerfast_rr;

# Masks repeated sequences using a custom algorithm that identifies and masks all sequences that occur in multiple 

# high scoring pairs in genome self-blast searches

 
#use strict;

#use warnings;

use File::Copy;

use utf8;

use open qw( :encoding(cp866) :std );

#############################
# INTERNAL ARGUMENTS
#############################

my $dir = "NEW_GENOMES";

my $help;

my $copies = 1;

my $mask_char = 'N';

my $Border = "@" x 40;

my $Separator = "-" x 40;

sub MASK { 

#############################
# DECLARE GLOBAL VARIABLES
#############################

  my $HomeDir = qx{pwd};

  my %SeqHash;

  my %SeqLength;

  my $prevqseqid;

  my %NumAlignments;

  my %MaskedHash;

  my $MaskTempHash;  

  my @MFAFilesList;

#############################
# RUN PROGRAM
#############################

  ## Check for a tempfile from an interrupted previous run. Read file and delete in preparation for new one.

  if(-f "Masking_tempfile") {

    open(MASKTEMPFILE, "Masking_tempfile");

    while(my $L = <MASKTEMPFILE>) {

      chomp($L);

      my ($Genome, $BlastStatus, $MaskStatus) = split(/\t/, $L);

      %{$MaskTempHash{$Genome}} = ('Blast' => $BlastStatus, 'Mask' => $MaskStatus);

    }

    close MASKTEMPFILE;

  }

  open(MASKTEMPFILE, '>>', "Masking_tempfile");

  & READ_DIRECTORY;

  & CHECK_FASTA_HEADERS;

  & PROCESS_FILES;

  close(MASKTEMPFILE, OUT);

  chomp($HomeDir);

  chdir($HomeDir)

}

##############################
# SUBROUTINES
##############################

sub READ_DIRECTORY {

	my $mfa;
	my @Masked_mfa;

        opendir(DIR, $dir) || die "Can't open directory. Please provide complete path\n";

	@MFAFilesList = readdir(DIR);

	chdir(DIR);

        print "\n$Border\n".
              "\nPROCESSING NEW GENOMES:\n\n";
        
        foreach $mfa (@MFAFilesList) {

		next if $mfa =~ /nin|nsq|nhr/;				# skip over blastdb files
		next if $MaskTempHash{$mfa}{Mask} eq "complete";	# skip over completed files from interrupted run
                if($mfa =~ /mask/) {					
			push @Masked_mfa, $mfa;
			next
		}
		if(($mfa !~ /^\./) && ($mfa =~ /fasta$|fsa$/)) {
                        push @FilesToProcess, $mfa;
                }
        }
	print "\n$Separator\nUnmasked genomes:\n";
        if(@FilesToProcess == 0) {
		print "\tnone in this run\n"
	}
	else {
		foreach my $Unmasked_mfa (@FilesToProcess) {
			print "\t$Unmasked_mfa\n"
        	}
	}
        print "\n$Separator\nMasked genomes:\n";
	if(@Masked_mfa == 0) {
                print "\tnone\n"
        }
        else {
		foreach my $Masked_mfa (@Masked_mfa) {
			print "\t$Masked_mfa\n"
		}
	}

}

sub CHECK_FASTA_HEADERS {

	print "\n$Border\n\nCHECKING GENOMES FOR CORRECT FASTA HEADERS\n\n";
	my $Line;
        my $Problem;
        my $Fasta_warning =     "Sequence header lines should identify the genome name and contig number.".
				"They cannot have spaces or special characters except for - (hyphen) and _ (underscore)\n".
                                "Recommended formats are: >Genome-name_contig1 or >Genome-name_contig00001, etc.\n".
                                "Note: 'contig' can be replaced by the following words (scaffold, NODE, mitochondrial, mt).\n".
				"Sample errors are shown below.\n\n";
	foreach my $mfaFile (@FilesToProcess) {
                my ($Error, $Warning, $Warning1, $Warning2, $Warning3, $Warning4) = (0, 0, 0, 0, 0, 0);
		$Problem = "no";
		print "$mfaFile\n";
		open(MFAFILE, $mfaFile);
		while($Line = <MFAFILE>) {
			if($Line =~ /^>|\d+/) {
				chomp($Line);			
				#$Line =~ s/\s+//;

                                # No genome identifier

				if($Line !~ /_/ && $Warning1 < 1) {
                                        if($Warning < 1) {
                                                print "Multifasta has one or more invalid FASTA headers.\n$Fasta_warning";
                                                $Warning ++;
                                        }
                                        print "\t$Genome-name identifier appears to be missing.\n\n";
                                        $Warning1 ++;
                                        $Problem = yes;
                                        $Error++
                                }

				# Characters following scaffold numbers

				if($Line !~ /chromosome\d+$|contig\d+$|mt\d+$|mitochondrial\d+$|scaffold\d+$|NODE\d+$/ && $Warning2 < 1) {
                                        if($Warning < 1) {
						print "Multifasta has one or more invalid FASTA headers.\n$Fasta_warning";
						$Warning ++;
					}
					print "\t$Line, contig/scaffold number should come last.\n\n";
                        		$Warning2 ++;
					$Problem = yes;
					$Error++
				}

		                # > symbol missing

				if($Line !~ /^>/ && $Line =~ /\d/ && $Warning3 < 1) {
                			if($Warning < 1) {
                                		print "Genome $mfaFile has an invalid FASTA file. $Fasta_warning";
                                		$Warning ++;
					}
					print "\t$Line, no > symbol.\n\n";             
                        		$Warning3 ++;
					$Problem = "yes";
					$Error ++
				}

		                # Illegal characters

				if($Line =~ /\s|\||\:|\./ && $Warning4 < 1) {
                			if($Warning < 1) {
                                		print "Genome $mfaFile has an invalid FASTA file. $Fasta_warning";
                                		$Warning ++;
                        		}
					print "\t$Line, only alphanumerics, hyphens and underscores allowed.\n\n";             
                        		$Warning4 ++;
					$Problem = "yes";
					$Error++
				}
			}
		}
		if($Error == 0) {
			print "passed fasta format check\n"
		}	
	}
	if($Problem eq "yes") {
		die "Please re-format fasta headers and restart the program. Note: you can skip the initial dialog by typing:\n".
		    "perl SNPcaller <PROJECT_NAME> run\n\n"
	}
}

sub PROCESS_FILES {

	if(@FilesToProcess < 1) {
		print "No new genomes to process. Checking for unfinished files";
        }

	else {
		 foreach my $mfa (@FilesToProcess) {
	                next if $mfa =~ /masked/;
	                if(($mfa !~ /^\./) && ($mfa =~ /fasta$|fsa$/)) {
				print MASKTEMPFILE "$mfa\t";		# write genome name to tempfile 
				print "\n$Separator\n";
				print "Masking genome: $mfa\n\n";

				GENOME_SELF_BLAST($mfa) unless $MaskTempHash{$mfa}{Blast} eq "complete";

				print MASKTEMPFILE "complete\t";

				READ_BLAST_FILE($mfa) unless $MaskTempHash{$mfa}{Mask} eq "complete";

				MASK_FASTA($mfa, \%RepeatSegments) unless $MaskTempHash{$mfa}{Mask} eq "complete";

			}
		}
	}

}


## BLAST each genome against itself

sub GENOME_SELF_BLAST {

	my ($mfa) = @_;

#	Generate blastable database.  Use backticks to prevent blast STDOUT from being seen by user

	print "Making blast database\n";
 
	my $blast = `makeblastdb -in $mfa -dbtype nucl -out $mfa 2>/dev/null`; # add a test to check if command failed

	print "Executing Genome Self BLAST.  This may take a while.  Please be patient.\n";

	if(system("blastn -query $mfa -db $mfa -max_target_seqs 20000 -evalue 1e-100 -out '$mfa.selfBLASTN_tempfile' -outfmt '6 qseqid qstart qend' 2>/dev/null" ) !=0 ) { 

		die "Error executing BLASTN search\t$!\n"

	}
				
}


## Read self BLAST file and count number of times each base position occurs in an alignment
 
sub READ_BLAST_FILE {

        my $CountString = "";

	my ($mfa) = @_;

	my $outfile = $mfa;

	my $TheLine;

	my($qseqid, $qstart, $qend);

	my $k;

	my @RepeatList;

        my %RepeatSegments;

        my $EOL;


# Read in selfBLAST tempfile

  open(BLAST, "$mfa.selfBLASTN_tempfile") || die "Can't open BLASTN results file\n";

  print "Analyzing BLASTN report\n";

  while($TheLine = <BLAST>) {

    chomp($TheLine);

    $prevqseqid = $qseqid;

    ($qseqid, $qstart, $qend) = split(/\t/, $TheLine);

    if(($prevqseqid ne $qseqid) || eof ){	# if current BLAST line is the first for the next query sequence (or end of file has been reached), create list of repeat regions in previous query  
              
      $prevqseqid = $qseqid if eof;

      $RepeatSegments{$prevqseqid} = [()];	#??? Clear out Repeat segments hash

      while($CountString =~ /(2+)/g) {		# examine count string for runs of 2's. Record lengths and positions and add to end of a list

        my $RepeatLength = length $1;

        push @RepeatList, $-[0];

        push @RepeatList, $RepeatLength;

      }

      $CountString = "";

      $RepeatSegments{$prevqseqid} = [@RepeatList];	# Hash the repeat lists based on query sequence id

      @RepeatList = ()

    }

    # Store alignment information for the current HSP

    $CountString .= "0" x ($qend - length $CountString) if $qend > length($CountString);	# Start a new count string if it doesn't exist. Zeroes indicate no alignment

    substr($CountString, $qstart - 1, 1 + $qend - $qstart) =~ tr/01/12/;			# Increment counts in repeated regions

  }

  close BLAST;

  unlink("$mfa.selfBLASTN_tempfile")

}


## Read in .fasta files and mask bases that occurred in multiple blast alignments

sub MASK_FASTA {

  my($mfa, $RepeatSegs_ref) = @_;		

  my $outfile = $mfa;

  if($outfile =~ /(\.+)/) {

    $outfile =~ s/(\.+)/_masked$1/

  }

  else {

    $outfile = $outfile."_masked"

  }


  open(OUT, '>', $outfile) || die "Can't open outfile $outfile for writing\n";

  %RepeatSegments = %$RepeatSegs_ref;

  my $EOL = $/;

  $/ = undef;

  # open genome sequence file

  open(GENOME, "$mfa") || die "Can't find genome sequence file\n";

  my $Genome = <GENOME>;

  # open database file for writing

  # read genome sequence

  my @Fasta_records = split(/>/, $Genome);

  shift @Fasta_records;

  print "Writing masked version of genome\n";

  foreach my $SeqEntry (@Fasta_records) {

    my ($Header, $Seq) = split(/\n/, $SeqEntry, 2);

    $Header =~ s/ .+//;

    $Seq =~ s/\W+//g;

    my @RepeatSeqs = @{$RepeatSegments{$Header}};

    while (my($Index, $Length) = splice(@RepeatSeqs, 0, 2)) {

      substr($Seq, $Index, $Length) =~ tr/gatcGATC/XXXXXXXX/;

    }

    if(length($Seq) > 100) {					# Print masked sequences if length > 100 nt

      print OUT ">$Header\n$Seq\n";

    }

  }

  close OUT;

  print "Masking complete\n";

  print MASKTEMPFILE "complete\n";

  $/ = $EOL

}


1;
