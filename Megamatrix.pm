package Megamatrix;

sub MATRIX {

  my $pwd = `pwd`;
  my $DateTime = $_[0];

  open(SNP_SUMMARY, "SNP_counts_cumulative.txt") || die "Cannot open SNPfile: $pwd\n";

  open(OUT, '>', "SNP_counts_$DateTime.meg");

  while($L = <SNP_SUMMARY>) {

    chomp($L);

    ($Q, $S, $SNPs, $Aligned, $Weighted) = split(/\t/, $L);

    $Q =~ s/_.+//;

    $S =~ s/_.+//;

    $Value = sprintf("%.0f", $Weighted);

    $Distances{$Q}{$S} = $Value;

    $Counting{$Q}{$S} ++;

    $Strains{$Q} = 1;

    $Strains{$S} = 1

  }

  close SNP_summary;

  @Strains = keys %Strains;

  $ListLen = @Strains;

  print OUT "#mega\n!Title:SNP counts $DateTime;\n!Format DataType=Distance DataFormat=UpperRight NTaxa=$ListLen;\n\n";

  foreach $Strain (@Strains) {

    $StrainNum ++;

    if($StrainNum < 10) {

      print OUT "[ $StrainNum] #$Strain\n";

    }

    else {

      print OUT "[$StrainNum] #$Strain\n";

    }

    print OUT "\n"

  }

  print OUT "[";

  for($i = 1; $i <= $ListLen; $i++) {

    print OUT "\t$i"

  }

  print OUT " ]\n";

  $StrainNum = 0;

  $Gap = "\t";

  while($ListLen > 1) {

    $Strain1 = shift(@Strains);

    $ListLen = @Strains;

    $StrainNum ++;

    if($StrainNum < 10) {

      print OUT "[ $StrainNum]";

    }

    else {

      print OUT "[$StrainNum]";

    }

    print OUT "$Gap";

    for($i = 0; $i <= $ListLen - 1; $i++) {

        $Strain2 = $Strains[$i];

        $Weighted_SNP1 = $Distances{$Strain1}{$Strain2};

        $Weighted_SNP2 = $Distances{$Strain2}{$Strain1};

        $Distance = sprintf("%.0f", $Weighted_SNP1+$Weighted_SNP2);

        if($i == 0) {

            print OUT "$Distance"

        }

        else {

            print OUT "\t$Distance"

        }

    }

    $Gap .= "\t";

    print OUT "\n";

  }         
    
  $StrainNum ++;

  print OUT "[$StrainNum]\n";

  close OUT

}

1;
