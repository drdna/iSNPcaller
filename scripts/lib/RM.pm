package RMSO;

#use strict;

#use warnings;

use File::Copy;

use utf8;

#use open qw( :encoding(cp866) :std );

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

sub RUN {

  $mfa = $_[0];

  print "\n########################";
  print "\nMasking genome: $mfa\n\n";

  GENOME_SELF_BLAST($mfa);

  $RepeatSegmentsRef = READ_BLAST_FILE($mfa);

  MASK_FASTA($mfa, $RepeatSegmentsRef)

}

sub GENOME_SELF_BLAST {

  $mfa = @_[0];

# Generate blastable database.  Use backticks to prevent blast STDOUT from being seen by user

 
  print "\n###########################";
  print "\nExecuting Genome Self BLAST.  This may take a while.  Please be patient.\n";

  if(system("blastn -query $mfa -subject $mfa -max_target_seqs 20000 -out '$mfa.selfBLASTN_tempfile' -outfmt '6 qseqid qstart qend'")) { # 2>/dev/null" ) !=0 ) { 

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

  print "\n#######################";
  print "\nAnalyzing BLASTN report\n\n";

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

  system("rm '$mfa.selfBLASTN_tempfile'");

  return(\%RepeatSegments)

}

sub MASK_FASTA {

  my ($mfa, $RepeatSegs_ref) = @_;

  my $outfile = $mfa;
  if ($outfile =~ /(\..+)/) {
    $outfile =~ s/(\..+)/_masked$1/;
  } else {
    $outfile = $outfile . "_masked";
  }

  open(OUT, '>', $outfile) || die "Can't open outfile for writing\n";

  my %RepeatSegments = %$RepeatSegs_ref;

  my $EOL = $/;
  $/ = undef;

  # open genome sequence file
  open(GENOME, "$mfa") || die "Can't find genome sequence file\n";
  my $Genome = <GENOME>;
  close GENOME;

  my @Entries = split(/>/, $Genome);
  shift @Entries;

  my $Counter = 0;

  foreach my $SeqEntry (@Entries) {

    my ($Header, $Seq) = split(/\n/, $SeqEntry, 2);
    $Header =~ s/ .+//;
    $Seq =~ s/[^A-Za-z]//g;

    $Counter++;
    print "Started writing masked version of genome\n" if $Counter == 1;

    my @RepeatSeqs = @{$RepeatSegments{$Header}};
    while (my ($Index, $Length) = splice(@RepeatSeqs, 0, 2)) {
      substr($Seq, $Index, $Length) =~ tr/A-Za-z/n/;
    }

    print "Writing masked version of $Header\n" if ($Counter % 100) == 0;

    # Check if there is any contiguous unmasked stretch >200 nt
    if ($Seq =~ /[ACGT]{201,}/i) {
      print OUT ">$Header\n$Seq\n";
    }
  }

  close OUT;

  $/ = $EOL;

  print "\nDONE\n";
}


1;
