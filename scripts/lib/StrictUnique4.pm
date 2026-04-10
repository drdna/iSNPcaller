package StrictUnique4;

use strict;
use warnings;
use Cwd qw(cwd abs_path);

# Global variables
my ($indir, $outdir, $cwd, $DateString);
our (%countString, %SNPsHash, %SNPsTempHash);
our (@QArray, @SArray);

sub SNPs {
    my ($arg1, $arg2) = @_;
    my $raw_in  = ref($arg1) eq 'ARRAY' ? $arg1->[0] : $arg1;
    my $raw_out = ref($arg2) eq 'ARRAY' ? $arg2->[0] : $arg2;

    $cwd = cwd();
    
    # Resolve absolute paths before changing directories
    $indir  = abs_path($raw_in);
    $outdir = abs_path($raw_out);

    die "ERROR: Input directory $raw_in not found.\n" if (!$indir || !-d $indir);

    unless (-d "$outdir") { mkdir "$outdir" or die $!; }
    unless (-d "$outdir/LOGS") { mkdir "$outdir/LOGS" or die $!; }

    &RUN_PROGRAM;
}

sub RUN_PROGRAM {
    chdir("$outdir") || die "Can't cd to $outdir\n";

    if (-f "SNPcall_tempfile") {
        open(TEMP_IN, '<', "SNPcall_tempfile");
        while (<TEMP_IN>) { chomp; $SNPsTempHash{$_} = 1; }
        close TEMP_IN;
    }
    
    open(SNPTEMPFILE, '>>', "SNPcall_tempfile");
    unless (-d "SNP_COUNTS") { mkdir "SNP_COUNTS" }

    &OPEN_OUTFILES;

    opendir(DIR, $indir) || die "Can't opendir $indir\n";
    my @FilesList = grep { /BLAST$/ && !/^\./ } readdir(DIR);
    closedir(DIR);

    foreach my $BlastFile (sort @FilesList) {
        next if exists $SNPsTempHash{$BlastFile};
        
        my @parts = split(/\./, $BlastFile);
        my $Q = $parts[0];
        my $S = $parts[1];
        
        my $Q_id = $Q; $Q_id =~ s/_.*//;
        my $S_id = $S; $S_id =~ s/_.*//;

        print "Processing: $BlastFile\n";
        print MAINREPORT "$Q_id\t$S_id\t";
        print CUMULATIVE "$Q_id\t$S_id\t";

        PRE_SCREEN($BlastFile);
        my ($SNPs, $MAligned, $MMAligned) = COUNT_SNPS($BlastFile, $Q_id, $S_id);
        
        # PRINT_SNPs is now called inside COUNT_SNPS to ensure OUT is open
        RESULTS($SNPs, $MAligned, $MMAligned, $BlastFile);
    }

    close MAINREPORT; close CUMULATIVE; close LOGFILE; close SNPTEMPFILE;
    chdir($cwd);
}

sub PRE_SCREEN {
    my ($File) = @_;
    undef @QArray; undef @SArray;
    open(BLAST, '<', "$indir/$File") || die "Can't open $indir/$File\n";
    while (<BLAST>) {
        chomp;
        my ($qid, $sid, $qpos, $qend, $spos, $send) = split(/\t/);
        my $qnum = ($qid =~ /(\d+)$/) ? $1 : 0;
        my $snum = ($sid =~ /(\d+)$/) ? $1 : 0;
        for my $k ($qpos .. $qend) { $QArray[$qnum][$k]++ }
        if ($spos < $send) { for my $j ($spos .. $send) { $SArray[$snum][$j]++ } }
        else { for my $j ($send .. $spos) { $SArray[$snum][$j]++ } }
    }
    close BLAST;
}

sub COUNT_SNPS {
    my ($File, $Q, $S) = @_;
    my ($SNPs, $MAligned, $MMAligned) = (0, 0, 0);

    open(BLAST, '<', "$indir/$File") || die $!;
    open(OUT, '>', "${Q}_v_${S}_out") || die $!;

    while (<BLAST>) {
        chomp;
        my ($qid, $sid, $qpos, $qend, $spos, $send, $btop) = split(/\t/);
        ($SNPs, $MAligned, $MMAligned) = ANALYZE_ALIGNMENT($qid, $sid, $qpos, $qend, $spos, $send, $btop, $SNPs, $MAligned, $MMAligned);
    }
    
    # Restored: Print SNPs to the OUT filehandle BEFORE closing it
    &PRINT_SNPs;

    close BLAST; 
    close OUT;
    return ($SNPs, $MAligned, $MMAligned);
}

sub ANALYZE_ALIGNMENT {
    my ($qid, $sid, $qpos, $qend, $spos, $send, $AlignSumm, $SNPs, $MAligned, $MMAligned) = @_;
    my @Align = split(/(\d+)/, $AlignSumm);
    shift @Align;
    my $AlignedSegment = 0;
    foreach my $Alignment (@Align) {
        $AlignedSegment ++;
        if($Alignment =~ /\d+/) {
            ($qpos, $spos, $MAligned) = COUNT_ALIGNED_MATCHES($qid, $sid, $qpos, $spos, $send, $Alignment, $MAligned);
            if($AlignedSegment == 1) {
                $qpos --;
                if($spos < $send) { $spos -- } else { $spos ++ }
            }
        }
        if($Alignment =~ /\w+/) {
            ($qpos, $spos, $SNPs, $MMAligned) = ANALYZE_ALIGNED_MISMATCHES($qid, $sid, $qpos, $spos, $send, $Alignment, $SNPs, $MMAligned)
        }
    }
    return ($SNPs, $MAligned, $MMAligned)
}

sub COUNT_ALIGNED_MATCHES {
    my($qid, $sid, $qpos, $spos, $send, $Alignment, $MAligned) = @_;
    my $qnum = ($qid =~ /(\d+)$/) ? $1 : 0;
    my $snum = ($sid =~ /(\d+)$/) ? $1 : 0;
    for(my $m = 0; $m <= $Alignment - 1; $m ++) {
        if(($QArray[$qnum][$qpos] == 1) && ($SArray[$snum][$spos] == 1 )) { $MAligned ++; }
        $qpos ++;
        if($spos < $send) { $spos ++ } else { $spos -- }
    }
    return($qpos, $spos, $MAligned)
}

sub ANALYZE_ALIGNED_MISMATCHES {
    my ($qid, $sid, $qpos, $spos, $send, $Alignment, $SNPs, $MMAligned) = @_;
    my $qnum = ($qid =~ /(\d+)$/) ? $1 : 0;
    my $snum = ($sid =~ /(\d+)$/) ? $1 : 0;
    my @MutsList = split(//, $Alignment);

    for(my $j = 0; $j <= $#MutsList-1; $j += 2) {
        my ($qnucl, $snucl) = ($MutsList[$j], $MutsList[$j+1]);

        # 1. Handle Mismatches (including Ns)
        if(($qnucl =~ /[A-Za-z]/) && ($snucl =~ /[A-Za-z]/)) {
            $qpos ++;
            if($spos < $send) { $spos ++ } else { $spos -- }

            # ONLY hash and count if they are actual nucleotides
            if(($qnucl =~ /[AGTCagtc]/) && ($snucl =~ /[AGTCagtc]/)) {
                
                my $strand = ($spos < $send) ? 'f' : 'r'; # Helper for readability
                $SNPsHash{$qid}{$qpos} = [($sid, $spos, $qnucl, $snucl, $strand)];

                if(($QArray[$qnum][$qpos] < 2) && ($SArray[$snum][$spos] < 2)) {
                    $SNPs ++; 
                    $MMAligned ++;
                }
            }
        }
        # 2. Handle Insertions in Subject (Gap in Query)
        elsif(($qnucl eq "-") && ($snucl =~ /[A-Za-z]/)) { 
            if($spos < $send) { $spos ++ } else { $spos -- } 
        }
        # 3. Handle Deletions in Subject (Gap in Subject)
        elsif(($qnucl =~ /[A-Za-z]/) && ($snucl eq "-")) { 
            $qpos ++; 
        }
    }
    return ($qpos, $spos, $SNPs, $MMAligned)
}

sub PRINT_SNPs {
    foreach my $qid (sort keys %SNPsHash) {
        foreach my $qpos (sort {$a <=> $b} keys %{$SNPsHash{$qid}}) {
            my @list = @{$SNPsHash{$qid}{$qpos}};
            my $sid = shift @list;
            print OUT "$qid\t$sid\t$qpos\t" . join("\t", @list) . "\n";
        }
    }
    undef %SNPsHash;
}

sub OPEN_OUTFILES {
    my $DateTime = localtime();
    my @dt = split(/\s+/, $DateTime);
    $DateString = join("_", $dt[1], $dt[2], $dt[4]);
    $DateString =~ s/\:/-/g;
    open(MAINREPORT, '>>', "SNP_counts_out_$DateString") || die $!;
    open(LOGFILE, '>>', "LOGS/logfile_$DateString") || die $!;
    open(CUMULATIVE, '>>', "SNP_counts_cumulative.txt") || die $!;
    my $old = select(MAINREPORT); $|=1; select(LOGFILE); $|=1; select($old);
}

sub RESULTS {
    my ($SNPs, $MAligned, $MMAligned, $BlastFile) = @_;
    my $TotalAligned = $MAligned + $MMAligned;
    if ($TotalAligned > 0) {
        my $Weighted = sprintf("%.0f", ($SNPs / $TotalAligned) * 1000000);
        print MAINREPORT "$SNPs\t$TotalAligned\t$Weighted\n";
        print CUMULATIVE "$SNPs\t$TotalAligned\t$Weighted\n";
        print SNPTEMPFILE "$BlastFile\n";
    } else {
        print MAINREPORT "0\t0\t0\n";
        print CUMULATIVE "0\t0\t0\n";
    }
}

1;
