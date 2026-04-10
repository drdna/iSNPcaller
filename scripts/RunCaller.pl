#!/usr/bin/perl
use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/lib";

use UniqueVariants;

die "Usage: $0 <blast-dir> <snps-dir>\n" unless @ARGV == 2;

UniqueVariants::SNPs(\@ARGV);
