#!/usr/bin/perl -w

use strict;
use warnings;

my $fastq_to_split = $ARGV[0];
my $barcode_file = $ARGV[1];

open(BARCODE,"<$barcode_file") or die ("Cannot open $barcode_file : $!");

while(<BARCODE>) {
	next if ($_ =~ m/^#/);
	my @infos = split(m/\t/,$_);
	my $identifier = $infos[0];
	my $barcode = $infos[1];
	open(TMP,">barcode_tmp") or die ("Cannot open barcode_tmp : $!");
	print TMP $identifier . "\t" . $barcode;
	close(TMP);
	my $launch_split = "zcat $fastq_to_split | fastx_barcode_splitter.pl --bcfile barcode_tmp --bol --exact --prefix output_barcode_split/ --suffix \".fastq\"";
	system("$launch_split");
	my $rm_tmp = "rm barcode_tmp";
	system($rm_tmp);
} 

