package RepeatmaskerfastMT;

#use strict;

#use warnings;

use File::Copy;

use constant MAX_PROCESSES => 48;

use Parallel::ForkManager;

use utf8;

use open qw( :encoding(cp866) :std );

#############################
# INTERNAL ARGUMENTS
#############################

my $help;

my $dir = "GENOMES";

my $copies = 1;

my $mask_char = 'n';

my $save = 'NO';

my $blastdb = 'NO';

sub MASK { 


#############################
# DECLARE GLOBAL VARIABLES
#############################

  my %SeqHash;

  my %SeqLength;

  my $prevqseqid;

  my %NumAlignments;

#############################
# RUN PROGRAM
#############################

  my $Home_dir = qx{pwd};

  chomp($Home_dir);

  & READ_DIRECTORY;

  close OUT;

  chdir($Home_dir);

  my $get = `ls`;

  if($get =~ /temp_genome_blastn_0/) {

    system('rm temp_genome_blastn_0')

  }

}

##############################
# SUBROUTINES
##############################

sub READ_DIRECTORY {

	my $mfa;

        opendir(DIR, $dir) || die "Can't open directory.  Please provide complete path\t$!\n";

	my @FilesList = readdir(DIR);

	chdir(DIR);


        print   "\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n".
        	"\nProcessing new genomes:\n\n";


	my $pm2 = Parallel::ForkManager->new(MAX_PROCESSES);

        foreach $mfa (@FilesList) {

                next if $MaskedHash{$mfa} eq "yes";

                if(($mfa !~ /^\./) && ($mfa =~ /fasta|fsa$/)) {

                        print "$mfa\n"

                }

        }

        print "\n";

	foreach $mfa (@FilesList) {

	        my $pid = $pm2->start and next;

		if($mfa =~ /masked.+fasta/) { # && $MaskedHash{$mfa} eq "yes") {

                        move("$mfa", "NEW_GENOMES/$mfa");

                        next

                }

                if(($mfa !~ /^\./) && ($mfa =~ /fasta|fsa$/)) {

			next if $MaskedHash{$mfa} eq "yes";

                        print "\nMasking genome: $mfa\n\n";

                        GENOME_SELF_BLAST($mfa);

                        READ_BLAST_FILE($mfa);

			MASK_FASTA($mfa, \%RepeatSegments)

		}

		$pm2 -> finish

	}

	$pm2->wait_all_children();

}


sub GENOME_SELF_BLAST {

	my ($mfa) = @_;

	print "Executing Genome Self BLAST.  This may take a while.  Please be patient.\n";

#	add an option to be verbose and print the following line

#	print "BLAST: query = $mfa db = $mfa\n";

	if(system("blastn -query $mfa -subject $mfa -max_target_seqs 20000 -evalue 1e-100 -out '$mfa.selfBLASTN_tempfile' -outfmt '6 qseqid qstart qend' 2>/dev/null" ) !=0 ) { 

		die "Error executing BLASTN search\n"

	}
				
}


sub READ_BLAST_FILE {

        my $CountString = "";

	my ($mfa) = @_;

	my $outfile = $mfa;

	my $TheLine;

	my($qseqid, $qstart, $qend);

	my $k;

	my $EOL;


# Read in selfBLAST tempfile

  open(BLAST, "$mfa.selfBLASTN_tempfile") || die "Can't open BLASTN results file\n";

  print "Analyzing BLASTN report\n\n";

  while($TheLine = <BLAST>) {

    chomp($TheLine);

    $prevqseqid = $qseqid;

    ($qseqid, $qstart, $qend) = split(/\t/, $TheLine);

    if(($prevqseqid ne $qseqid) || eof ){
              
      $prevqseqid = $qseqid if eof;

      $RepeatSegments{$prevqseqid} = [()];

      while($CountString =~ /(2+)/g) {

        $RepeatLength = length $1;

        push @RepeatList, $-[0];

        push @RepeatList, $RepeatLength;

      }

      $CountString = "";

      $RepeatSegments{$prevqseqid} = [@RepeatList];

      @RepeatList = ()

    }

    # Store alignment information for the current HSP

    $CountString .= "0" x ($qend - length $CountString) if $qend > length($CountString);

    substr($CountString, $qstart - 1, 1 + $qend - $qstart) =~ tr/01/12/;

  }

  close BLAST;

}

sub MASK_FASTA {

  ($mfa, $RepeatSegs_ref) = @_;

  $outfile = $mfa;

  if($outfile =~ /(\.+)/) {

    $outfile =~ s/(\.+)/_masked$1/

  }

  else {

    $outfile = $outfile."_masked"

  }


  open(OUT, '>', $outfile) || die "Can't open outfile for writing\n";

  %RepeatSegments = %$RepeatSegs_ref;

  $EOL = $/;

  $/ = undef;

# Initialize database hash to ensure it is empty

  # open genome sequence file

  open(GENOME, "$mfa") || die "Can't find genome sequence file\n";

  $Genome = <GENOME>;

  # open database file for writing

  # read genome sequence and create a database entries

  @Entries = split(/>/, $Genome);

  shift @Entries;

  foreach $SeqEntry (@Entries) {

    ($Header, $Seq) = split(/\n/, $SeqEntry, 2);

    $Header =~ s/ .+//;

    $Seq =~ s/\W+//g;

    $Counter ++;

    if($Counter == 1) {

      print "Started writing masked version of genome\n";

    }

    @RepeatSeqs = @{$RepeatSegments{$Header}};

    while (($Index, $Length) = splice(@RepeatSeqs, 0, 2)) {

      substr($Seq, $Index, $Length) =~ tr/A-Za-z/nnnnnnnn/;

    }

    if(length($Seq) > 100) {

      print OUT ">$Header\n$Seq\n";

    }

  }

  unlink "$mfa.selfBLASTN_tempfile";

  $Counter = 0;

  %RepeatSegments = undef;

  $/ = $EOL

}


1;
