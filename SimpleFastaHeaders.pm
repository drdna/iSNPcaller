package SimpleFastaHeaders;

sub rehead {

  $indir = $_[0];

  HEADER($indir) if -f $indir;

  & READ_DIR if -d $indir;

}

sub READ_DIR {

  $indir =~ s/\/$//;

  opendir(FASTADIR, $indir) || die "There is not a directory of that name in the path you specified. Please try again\n\n";

  @FASTA = readdir(FASTADIR);

  $Success = 'no';

  foreach $Fasta (@FASTA) {

    if($Fasta =~ /fasta|fsa|fa|fna|mfa$/) {

      next if $Fasta =~ /_nh/;

      HEADER($Fasta);

      $Success='yes'

    }

  }

  close FASTADIR;

}

  
sub HEADER {

  my $SeqNo = 0;

  my $Fasta = @_[0];

  ($oldGenome_ID = $Fasta) =~ s/_.+|\..+//;
  
  if(@_ == 2) {

    $Genome_ID = $ARGV[1];

    $outfile = $Genome_ID."_nh.fasta";

  }

  else {

    ($Genome_ID = $Fasta) =~ s/_.+|\..+//;

    ($outfile = $Fasta) =~ s/\.\w+$/_nh\.fasta/;

  }


  if(-d $indir) {

    open(FASTAOUT, '>', "$indir/$outfile") || die "Problem creating corrected genome file: $!\n";

    open(MAPOUT, '>', "$indir/$Genome_ID"."_contig_map.txt") || die "Problem creating contig mapping file\n" ;

  }

  elsif(-f $indir) {

    open(FASTAOUT, '>', "$outfile") || die "Problem creating corrected genome file: $!\n";

    open(MAPOUT, '>', "$Genome_ID"."_contig_map.txt") || die "Problem creating contig mapping file\n";

  }

  print MAPOUT "Old-genomeID\tNew-genomeID\n$oldGenome_ID\t$Genome_ID\nNew_ID\t\tOld_ID\n$oldGenome_ID\t$Genome_ID\n";


  open(FASTA, "$indir/$Fasta");

  while($Line = <FASTA>) {

    if($Line =~ /^>/) {

      $SeqNo ++;

      print FASTAOUT ">$Genome_ID"."_contig$SeqNo\n";

      $Line =~ s/^>//;

      print MAPOUT "$Genome_ID"."_contig$SeqNo\t$Line"

    }

    else {

      print FASTAOUT "$Line"

    }

  }

  close FASTA;

  unlink("$indir/$Fasta");

  close FASTAOUT;

  close MAPOUT;

  print   "\n########################################################\n\n".
        "Sequence headers have been converted and written to the file: $outfile\n".
	"Mapping between new and old sequence IDs was written to the file: $Genome_ID"."_contig_map.txt\n\n"

}


1;
