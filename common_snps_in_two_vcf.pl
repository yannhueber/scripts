#!/usr/bin/perl -w

use strict;
use warnings;

my $vcf_light = $ARGV[0]; #VCF with the smallest size
my $vcf_big = $ARGV[1];

open(VCF1,"<$vcf_light") or die ("Cannot open $vcf_light : $!\n");
my %inmind;
my $number_snp1=0;
while(<VCF1>) {
	next if ($_ =~ m/^#/);
	$number_snp1++;
	my @infos_line = split(m/\t/,$_);
	my $chr = $infos_line[0];
	my $pos = $infos_line[1];
	my $id = $chr . "" . $pos;
	$inmind{$id} = 1;
}
close(VCF1);



my $number_snp2=0;
my $common=0;
open(VCF2,"<$vcf_big") or die ("Cannot open $vcf_big : $!\n");
while(<VCF2>) {
	next if ($_ =~ m/^#/);
	$number_snp2++;
	my @infos_line = split(m/\t/,$_);
	my $chr = $infos_line[0];
	my $pos = $infos_line[1];
	my $id2 = $chr . "" . $pos;
	if ($inmind{$id2}) {
		$common++;
	}
}
close(VCF2);

print "Number of snps in file $vcf_light : $number_snp1\n";
print "Number of snps in file $vcf_big : $number_snp2\n";
print "SNPs in common : $common\n";


