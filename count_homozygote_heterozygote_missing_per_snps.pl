#!/usr/bin/perl -w

#Take an input VCF file and count homozygous, heterozygous and missing positions per snp

use warnings;
use strict;
use Getopt::Long;


if (@ARGV < 4)
{
	print "\n\nUsage: perl $0 -i vcf_in -o results_output\n\n\n";
	exit;
}


#variables
my ($input_vcf, $output);

GetOptions ("i|vcf_in=s" => \$input_vcf,	#string
	    "o|output=s" => \$output)		#string
or die ("Error in command line arguments\n");



open(VCF,"<$input_vcf") or die("Cannot open file $input_vcf : $!");
open(OUT,">$output") or die("Cannot open file $output : $!");

print OUT "CHR\tPOSITION\tHETEROZYGOUS\tHOMOZYGOUS\tMISSING\tTOTAL\n";


#Running through vcf
while (my $vcf_line = <VCF>) {
	next if ($vcf_line =~ m/^#/);
	chomp($vcf_line);
	my @infos_line = split(m/\t/,$vcf_line);

	my $homozygous_number=0;
	my $heterozygous_number=0;
	my $missing_geno_number=0;

	#snps
	my $indiv_count=0;
	my $chr = $infos_line[0];
	my $position = $infos_line[1];
	my @individuals_infos = @infos_line[9..$#infos_line];
	foreach my $indiv_infos (@individuals_infos) {
		$indiv_count++;
		my @indiv_infos_split = split(m/:/,$indiv_infos);
		my $genotype = $indiv_infos_split[0];
		if ($genotype =~ m/\.\/\./) {
			$missing_geno_number++;
		} elsif ($genotype =~ m/(\d{1})\/\1/) {
			$homozygous_number++;
		} else {
			$heterozygous_number++;
		}
	}

	#print output
	print OUT "$chr\t$position\t$heterozygous_number\t$homozygous_number\t$missing_geno_number\t$indiv_count\n";

}	


close(VCF);
close(OUT);
