#!/usr/bin/perl -w

#Take an input VCF file and count homozygous, heterozygous and missing positions per individual

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

my %indv_hash;
my @individuals_name;

GetOptions ("i|vcf_in=s" => \$input_vcf,	#string
	    "o|output=s" => \$output)		#string
or die ("Error in command line arguments\n");



open(VCF,"<$input_vcf") or die("Cannot open file $input_vcf : $!");
open(OUT,">$output") or die("Cannot open file $output : $!");

print OUT "INDIV\tHETEROZYGOUS\tHOMOZYGOUS\tMISSING\tTOTAL\n";


#Running through vcf
while (my $vcf_line = <VCF>) {
	next if ($vcf_line =~ m/^##/);
	
	chomp($vcf_line);
	my @infos_line = split(m/\t/,$vcf_line);
	
	#individuals
	my @indiv_infos = @infos_line[9..$#infos_line];

	if ($vcf_line =~ m/^#CHROM/) {
		@individuals_name = @indiv_infos;
	} else {
		my @genotypes = @indiv_infos;
		my $i=0;
		foreach my $genotype (@genotypes) {
			my $indv = $individuals_name[$i];
			$indv_hash{$indv}{"total_count"}++;

			if ($genotype =~ m/\.\/\./) {
				$indv_hash{$indv}{"missing"}++;
			} elsif ($genotype =~ m/(\d{1})\/\1/) {
				$indv_hash{$indv}{"hom"}++;
			} else {
				$indv_hash{$indv}{"het"}++;
			}
			$i++;
		}
	}
}


#print output
foreach my $indiv (keys(%indv_hash)) {
	if (! $indv_hash{$indiv}{"missing"}) {
		$indv_hash{$indiv}{"missing"} =0;
	}
	print OUT $indiv . "\t" . $indv_hash{$indiv}{"het"} . "\t" . $indv_hash{$indiv}{"hom"} . "\t" . $indv_hash{$indiv}{"missing"} . "\t" . $indv_hash{$indiv}{"total_count"} . "\n";
}


close(VCF);
close(OUT);
