#!/usr/bin/perl

# version: 1.1

die "Usage: perl Pairwise_matrix_MEGA.pl <SNPsummaryFile>\n" if @ARGV != 1;

open(DISALLOWED, $ARGV[1]);

while($D = <DISALLOWED>) {

  chomp($D);

  $Disallowed{$D} = 1

}

open(SNP_SUMMARY, $ARGV[0]);

while($L = <SNP_SUMMARY>) {

    chomp($L);

    if($L =~ /masked/) {

        @LineList = split(/\t/, $L);

        ($Q, $S) = @LineList;

	$Q =~ s/Query: //;

        $Q =~ s/_.+//;

        $S =~ s/Subject: //;

        $S =~ s/_.+//;

    }

    if($L =~ /Weighted SNPs = (\d+)/) {

        $Value = sprintf("%.0f",$1);

        next if exists($Disallowed{$Q}) || exists($Disallowed{$S});

        next if exists($Distances{$Q}{$S}) || exists($flipDistances{$S}{$Q});

        $Distances{$Q}{$S} = $Value;

        $flipDistances{$S}{$Q} = 1;

        $Counting{$Q}{$S} ++;

        $Strains{$Q} = 1;

        $Strains{$S} = 1

    }

}

close SNP_summary;

@Strains = keys %Strains;

$ListLen = @Strains;

$ARGLen = @ARGV;

#print "$ARGLen\n";

for ($t = 1; $t <= @ARGV-1; $t++) {

    $Title .= " $ARGV[$t]"

}

print "#mega\n!Title:$Title;\n!Format DataType=Distance DataFormat=UpperRight NTaxa=$ListLen;\n\n";

foreach $Strain (@Strains) {

  $StrainNum ++;

  if($StrainNum < 10) {

    print "[ $StrainNum] #$Strain\n";

  }

  else {

    print "[$StrainNum] #$Strain\n";

  }

  print "\n"

}

print "[";

for($i = 1; $i <= $ListLen; $i++) {

  print "\t$i"

}

print " ]\n";

$StrainNum = 0;

$Gap = "\t";

while($ListLen > 1) {

    $Strain1 = shift(@Strains);

    $ListLen = @Strains;

    $StrainNum ++;

    if($StrainNum < 10) {

      print "[ $StrainNum]";

    }

    else {

      print "[$StrainNum]";

    }

    print "$Gap";

    for($i = 0; $i <= $ListLen - 1; $i++) {

        $Strain2 = $Strains[$i];

#        if($Counting{$Strain1}{$Strain2} > 1) {

#            "Result missing: $Strain1 2: $Strain2\n"

#        }

#        if($Counting{$Strain2}{$Strain1} > 1) {

#            print "1: $Strain2 2: $Strain1\n"

#        }
 

        $Weighted_SNP1 = $Distances{$Strain1}{$Strain2};

        $Weighted_SNP2 = $Distances{$Strain2}{$Strain1};

	if($Weighted_SNP1 eq "") {

            print OUT "Result_missing: $Strain1 $Strain2\n"

        }

        if($Weighted_SNP2 eq "") {

             print OUT "Result_missing: $Strain2 $Strain1\n" 

        }
        

        $Average = sprintf("%.0f", $Weighted_SNP1+$Weighted_SNP2);

        if($i == 0) {

            print "$Average"

        }

        else {

            print "\t$Average"

        }

    }

    $Gap .= "\t";

    print "\n";

}         
    
$StrainNum ++;

print "[$StrainNum]\n";

close OUT
