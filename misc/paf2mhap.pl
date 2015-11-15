#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Std;

my %opts = ();
getopts("2p", \%opts);
my $is_dbl = defined($opts{2});
my $is_100 = defined($opts{p});

die("Usage: paf2mhap.pl [-2p] <in.fa> <in.paf>\n") if (@ARGV == 0);

warn("Parsing FASTA to create the name<=>id table...\n");
my %hash;
my $fn = shift(@ARGV);
open(FH, $fn =~ /\.gz$/? "gzip -dc {} |" : $fn) || die;
my $cnt = 0;
while (<FH>) {
	if (/^>(\S+)/) {
		$hash{$1} = ++$cnt unless defined($hash{$1});
	}
}
close(FH);

warn("Converting PAF to MHAP format...\n");
while (<>) {
	chomp;
	my @t = split;
	next if ($t[0] eq $t[5]); # NB: ignore self matches
	next if (!$is_dbl && $t[0] gt $t[5]);
	my $cnt = /cm:i:(\d+)/? $1 : 0;
	my $r = $t[9] / $t[10];
	$r = sprintf("%.4f", $is_100? 100. * $r : $r);
	die if !defined($hash{$t[0]}) || !defined($hash{$t[5]});
	print(join(" ", $hash{$t[0]}, $hash{$t[5]}, $r, $cnt, 0, @t[2,3,1], $t[4] eq '+'? 0 : 1, @t[7,8,6]), "\n");
}
