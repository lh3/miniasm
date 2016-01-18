#!/usr/bin/perl

use strict;
use warnings;

my $fn = shift(@ARGV);
open(FH, $fn =~ /\.gz$/? "gzip -dc $fn|" : $fn) || die;
my %h;
while (<FH>) {
	my @t = split;
	$h{"$t[0]\t$t[1]"} = 1;
}
close(FH);

while (<>) {
	my @t = split;
	$h{"$t[0]\t$t[5]"} = 2 if ($h{"$t[0]\t$t[5]"});
	$h{"$t[5]\t$t[0]"} = 2 if ($h{"$t[5]\t$t[0]"});
}

my @cnt = (0, 0);
for my $x (keys %h) {
	++$cnt[$h{$x}-1];
}
print("$cnt[0]\t$cnt[1]\t", $cnt[1]/($cnt[0]+$cnt[1]), "\n");
