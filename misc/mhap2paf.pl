#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Std;

my %opts = ();
getopts("2", \%opts);
my $is_dbl = defined($opts{2});

die("Usage: mhap2paf.pl [-1] <in.mhap>\n") if (@ARGV == 0 && -t STDIN);

while (<>) {
	chomp;
	my @t = split;
	my $bl = $t[6] - $t[5] > $t[10] - $t[9]? $t[6] - $t[5] : $t[10] - $t[9];
	my $r = $t[2];
	my $ml = int($bl * ($r <= 1.? $r : .01 * $r) + .499);
	my $cm = "cm:i:" . int($t[3] + .499);
	my $rev = $t[4] == $t[8]? '+' : '-';
	print(join("\t", @t[0,7,5,6], $rev, @t[1,11,9,10], $ml, $bl, 255, $cm), "\n");
	print(join("\t", @t[1,11,9,10], $rev, @t[0,7,5,6], $ml, $bl, 255, $cm), "\n") if ($is_dbl);
}
