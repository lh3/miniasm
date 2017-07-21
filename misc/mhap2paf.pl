#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Std;

my %opts = ();
getopts("2f:l:", \%opts);
die("Usage: mhap2paf.pl [-2] [-f name_list] [-l min_len] <in.mhap>\n") if (@ARGV == 0 && -t STDIN);

my $is_dbl = defined($opts{2});
my $min_blen = defined($opts{l})? $opts{l} : 0;
my @a = ();
if (defined $opts{f}) {
	open(FH, $opts{f} =~ /\.gz$/? "gzip -dc $opts{f} |" : $opts{f}) || die;
	while (<FH>) {
		chomp;
		my @t = split;
		push(@a, $t[0]);
	}
	close(FH);
}

while (<>) {
	chomp;
	my @t = split;
	my $bl = $t[6] - $t[5] > $t[10] - $t[9]? $t[6] - $t[5] : $t[10] - $t[9];
	my $r = $t[2];
	my $ml = int($bl * ($r <= 1.? $r : .01 * $r) + .499);
	my $cm = "cm:i:" . int($t[3] + .499);
	my $rev = $t[4] == $t[8]? '+' : '-';
	next if $bl < $min_blen;
	if (@a) {
		$t[0] = $a[$t[0]-1];
		$t[1] = $a[$t[1]-1];
	}
	print(join("\t", @t[0,7,5,6], $rev, @t[1,11,9,10], $ml, $bl, 255, $cm), "\n");
	print(join("\t", @t[1,11,9,10], $rev, @t[0,7,5,6], $ml, $bl, 255, $cm), "\n") if ($is_dbl);
}
