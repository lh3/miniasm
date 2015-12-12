#!/usr/bin/env perl

use strict;
use warnings;

while (<>) {
	chomp;
	my @t = split("\t");
	if ($t[4] eq '-') {
		@t[3,4] = ($t[2] - $t[4], $t[2] - $t[3]);
	}
	if ($t[6] eq '-') {
		@t[8,9] = ($t[7] - $t[9], $t[7] - $t[8]);
	}
	my $bl = $t[12] + $t[13] + $t[14] + $t[15];
	print join("\t", @t[0,2..4], $t[1] eq $t[6]? '+' : '-', @t[5,7..9,12], $bl, 255), "\n";
}
