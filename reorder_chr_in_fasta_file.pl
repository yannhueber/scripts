#!/usr/bin/perl 

use strict;
use warnings;

my @order = ("10","11","12","1","2","3","4","5","6","7","8","9");

my $fasta = $ARGV[0];
my %memory;
open(FAS,"<$fasta") or die ("Cannot open $fasta...\n");

my $first = 1;
my $current_chr;
while(my $line = <FAS>) {
	if ($line =~ m/^>(\d+)/) {
		$current_chr = $1;
		print"$current_chr\n";
		$memory{$current_chr} = "";
	} else {
		$memory{$current_chr} .= $line;
	}
}


open(OUT,">$fasta.reordered") or die ("Cannot open $fasta.reordered...\n");
foreach my $chr (@order) {
	print OUT ">" . $chr . "\n" . $memory{$chr};
}

close(FAS);
close(OUT);









