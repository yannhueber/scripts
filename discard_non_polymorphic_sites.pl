#!/usr/bin/perl -w

#Take an input VCF file and discard non polymorphic sites.

use warnings;
use strict;
use Getopt::Long;

#variables
my ($input_vcf, $out_vcf);

if (scalar(@ARGV) < 4)
{
	print "\n\nUsage : perl $0 -i vcf_in -o vcf_out\n\n\n";
	exit;
}

GetOptions ("i|vcf_in=s" => \$input_vcf,		#string
	    "o|out=s" => \$out_vcf)			#string
or die ("Error in command line arguments\n");



#Open files
open(VCF_IN,"<$input_vcf") or die("Cannot open file $input_vcf : $!");
open(VCF_OUT,">$out_vcf") or die("Cannot open file $out_vcf : $!");

#Running through vcf file
while (defined (my $vcf_line = <VCF_IN>)) {
	if ($vcf_line =~ m/^#/) {
		print VCF_OUT $vcf_line;
		next;
	}
	chomp($vcf_line);
	my @infos_line = split(m/\t/,$vcf_line);
	my @individuals_infos = @infos_line[9..$#infos_line];
	
	my $i =0;
	my %hash;
	foreach my $indiv_infos (@individuals_infos) {
		my @indiv_infos_split = split(m/:/,$indiv_infos);
		my $genotype = $indiv_infos_split[0];
		next if ($genotype =~ m/\.\/\./);
		my @genotype_split = split(m/\//,$genotype);
		foreach my $genotype (@genotype_split) {
			if (! defined ($hash{$genotype})) {
				$i++;
				$hash{$genotype}++;
			}
		}
		if ($i > 1) {
			print VCF_OUT $vcf_line . "\n";
			last;
		}
	}
}
close(VCF_IN);
close(VCF_OUT);
