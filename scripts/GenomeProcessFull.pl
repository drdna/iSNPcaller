#!/usr/bin/perl

use strict;
use warnings;
use RM;
use Parallel::ForkManager;
use File::Spec;  # Essential for cross-platform path handling

# --- 1. Automatic Resource Detection ---
sub get_core_count {
    my $cores = 2; 
    if ($^O eq 'linux') {
        $cores = `grep -c ^processor /proc/cpuinfo`;
    } elsif ($^O eq 'darwin') {
        $cores = `sysctl -n hw.ncpu`;
    } elsif ($^O eq 'MSWin32') {
        $cores = $ENV{NUMBER_OF_PROCESSORS};
    }
    chomp($cores);
    return $cores > 0 ? $cores : 2;
}

# --- 2. Input Handling & Setup ---
die "Usage: perl MasterScript.pl <dirname/filename> <optional: genome-ID> <optional: min-length>\n" if @ARGV < 1;

my $input    = $ARGV[0];
my $id_arg   = $ARGV[1] || ""; 
my $min_len  = $ARGV[2] || 0;  
my $max_proc = int(get_core_count() / 2) || 1;

my $pm = Parallel::ForkManager->new($max_proc);

my @queue;
if (-d $input) {
    opendir(my $dh, $input) || die "Can't open directory: $!\n";
    @queue = map { File::Spec->rel2abs(File::Spec->catfile($input, $_)) } grep { /\.(fasta|fsa|fa|fna|mfa)$/i } readdir($dh);
    closedir($dh);
} elsif (-f $input) {
    push(@queue, File::Spec->rel2abs($input));
} else {
    die "Error: '$input' not found.\n";
}

print "Starting pipeline on " . scalar(@queue) . " file(s) using $max_proc cores...\n";

# --- 3. The Execution Loop ---
foreach my $infile (@queue) {
    my $pid = (scalar(@queue) > 1) ? $pm->start : 0;
    
    if ($pid == 0) { 
        process_fasta($infile, $id_arg, $min_len);
        $pm->finish if (scalar(@queue) > 1); 
    }
}

$pm->wait_all_children;
print "\nAll tasks finished successfully.\n";

# --- 4. The Unified Processing Subroutine ---
sub process_fasta {
    my ($file_path, $requested_id, $threshold) = @_;
    
    # Split the path into directory, and filename
    my ($volume, $directories, $filename) = File::Spec->splitpath($file_path);

    # --- Targeted Assembler Metadata Removal ---
    my $genome_id = ($requested_id ne "") ? $requested_id : $filename;

    # 1. Strip file extension first so it doesn't interfere
    $genome_id =~ s/\.[^.]+$//;

    # 2. Strip specific assembler keywords that appear at the END of the string
    # This matches: _scaffolds, _contigs, _masked, _nh, etc.
    # We use 'i' for case-insensitivity and '$' to ensure it's at the end.
    $genome_id =~ s/_(scaffolds|scaf|contigs|ctg|assembly)$//i;

    # 3. Repeat once more in case they are stacked (e.g., _scaffolds_masked)
    $genome_id =~ s/_(scaffolds|scaf|contigs|ctg|assembly)$//i;

    # 4. Now convert the USER'S meaningful underscores to hyphens
    $genome_id =~ s/_/-/g;

    # 5. Final Alphanumeric Polish
    $genome_id =~ s/[^a-zA-Z0-9-]//g; # Remove punctuation
    $genome_id =~ s/-+/-/g;          # Collapse hyphens
    $genome_id =~ s/^-|-$//g;        # Trim edges

    # Create full paths for output files in the SAME directory as the input
    my $outfile_name = $genome_id . "_nh.fasta";
    my $mapfile_name = $genome_id . "_contig_map.txt";
    
    my $out_path = File::Spec->catpath($volume, $directories, $outfile_name);
    my $map_path = File::Spec->catpath($volume, $directories, $mapfile_name);

    # Open handles
    open(my $IN,  '<', $file_path) or die "Cannot open $file_path: $!\n";
    open(my $OUT, '>', $out_path)  or die "Cannot create $out_path: $!\n";
    open(my $MAP, '>', $map_path)  or die "Cannot create $map_path: $!\n";

    print $MAP "New_ID\tOld_Header\tLength\n";

    my $count = 0;
    my ($h, $s) = ("", "");

    my $write_seq = sub {
        my ($hdr, $seq) = @_;
        if (length($seq) >= $threshold) {
            $count++;
            my $new_id = "${genome_id}_contig$count";
            print $OUT ">$new_id\n$seq\n";
            print $MAP "$new_id\t$hdr\t" . length($seq) . "\n";
        }
    };

    while (my $line = <$IN>) {
        chomp $line;
        if ($line =~ /^>(.*)/) {
            $write_seq->($h, $s) if $h;
            $h = $1; $s = "";
        } else {
            $s .= $line;
        }
    }
    $write_seq->($h, $s) if $h; 

    close $IN; close $OUT; close $MAP;

    # --- Step 2: Run RM on the output ---
    if ($count > 0) {
        print "[PID $$] Saved output to: $directories\n";
        RM::RUN($out_path);
    }
}
