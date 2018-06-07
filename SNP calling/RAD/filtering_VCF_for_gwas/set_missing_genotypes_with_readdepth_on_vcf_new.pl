#!/usr/bin/perl -w


use warnings;
use strict;
use Getopt::Long;
use Data::Dumper;

#variables
my ($input_vcf, $read_depth_limit_min, $read_depth_limit_max, $read_depth_limit_per_allele, $max_indiv_nb_with_missing_genotype, $help);
my $output_vcf = "output.vcf";
my $missing_data_genotype = "./.:0,0:0";
sub usage();
sub read_depth_allele_ok($$$);


#parse command line
usage() if (scalar @ARGV==0);

GetOptions ("i|vcfin=s" => \$input_vcf,                #string
            "l|readdepthlimit_min=i" => \$read_depth_limit_min, #int
	    "m|readdepthlimit_max=i" => \$read_depth_limit_max, #int
	    "n|readdepthmin_allele=i" => \$read_depth_limit_per_allele, #int
            "d|maxindivwithmissingdata=i" => \$max_indiv_nb_with_missing_genotype, #int
	    "o|vcfout=s" => \$output_vcf) 		#string
or die ("Error in command line arguments\n");

usage() if ($help);

die "Error vcf file not specified (use '-i or --vcfin [FILENAME]')\n" unless defined $input_vcf;

die "Error invalid read depth limit min to output missing genotype (from 1 to infinite)\n" if ($read_depth_limit_min < 1 or $read_depth_limit_min > $read_depth_limit_max);

die "Error invalid read depth limit max to output missing genotype \n" if ($read_depth_limit_max <= $read_depth_limit_min);

die "Error invalid read depth limit min per allele (from 1 to infinite)\n" if ($read_depth_limit_per_allele < 1);

die "Error invalid maximum individuals with a missing genotype (from 0 to infinite)\n" if ($max_indiv_nb_with_missing_genotype < 0);



#Open vcf file
open(VCF_IN,"<$input_vcf") or die("Cannot open file $input_vcf : $!");
open(VCF_OUT,">$output_vcf") or die ("Cannot open file $output_vcf : $!");


#Running through vcf file
while (defined (my $vcf_line = <VCF_IN>)) {
	if ($vcf_line =~ m/^#/) {
		print VCF_OUT $vcf_line;
		next;
	}
	chomp($vcf_line);
	my @infos_line = split(m/\t/,$vcf_line);
	my $ref = $infos_line[3];
	my $alt = $infos_line[4];
	my @individuals_infos = @infos_line[9..$#infos_line];

	#New indiv genotypes
	my $indiv_with_missing_data_count=0;
	my @new_individuals_genotypes =();

	foreach my $indiv_infos (@individuals_infos) { #a loop for all the individuals
		my @indiv_infos_split = split(m/:/,$indiv_infos);
		my $genotype = $indiv_infos_split[0];
		my $total_read_depth = $indiv_infos_split[2];
		my $alleles_read_depth = $indiv_infos_split[1];

		my $rd_alleles_ok = read_depth_allele_ok($genotype, $alleles_read_depth, $read_depth_limit_per_allele);
	
		if ($genotype =~ m/\.\/\./) {
			push(@new_individuals_genotypes,$missing_data_genotype);
			$indiv_with_missing_data_count++;
		}
		elsif ( (!$rd_alleles_ok) or ($total_read_depth < $read_depth_limit_min) or ($total_read_depth > $read_depth_limit_max) ) {
			push(@new_individuals_genotypes,$missing_data_genotype);
			$indiv_with_missing_data_count++;
		} else {
			push(@new_individuals_genotypes,$indiv_infos);
		}
	}
	
	print "coucou le nb dindiv missing est de $indiv_with_missing_data_count et la limite est fixe a $max_indiv_nb_with_missing_genotype\n";

	#Discard sites with too much missing genotypes
	if ($indiv_with_missing_data_count <= $max_indiv_nb_with_missing_genotype) {
		#print into VCF File
		my $line_to_print="";
		foreach my $info (@infos_line[0..8]) {
			$line_to_print .= $info . "\t";
		}
		foreach my $indiv_genotype (@new_individuals_genotypes) {
			$line_to_print .= $indiv_genotype . "\t";
		}
		chomp($line_to_print);
		print VCF_OUT $line_to_print . "\n";
	}
}

close(VCF_IN);
close(VCF_OUT);



sub read_depth_allele_ok($$$) {
	my $geno = shift;
	my $allele_rd = shift;
	my $rd_limit_per_al = shift;
	my @split_al_rd = split(m/,/,$allele_rd);
	my $ok=0;

	if ($geno =~ m/\.\/\./) {
		$ok=0;
	} elsif ($geno =~ m/(\d)\/\1/) { #homozygote
		my $match = $1;
		my $depth = $split_al_rd[$match];
		if ($depth >= $rd_limit_per_al) {
			$ok=1;
		}
	} else { #heterozygote
		my @split_geno = split(m/\//,$geno);
		my $one = $split_geno[0];
		my $two = $split_geno[1];
		my $depth1 = $split_al_rd[$one];
		my $depth2 = $split_al_rd[$two];
		if ($depth1 >= $rd_limit_per_al and $depth2 >= $rd_limit_per_al) {
			$ok=1;
		}
	}
	return $ok;
}




sub usage() {
print<<EOF;


This script reads a VCF file and set to missing genotype each genotype that has a read depth below a defined parameter. It removes also the sites having a number of missing genotypes (for all the individuals) above another defined parameter.

usage: perl $0 -i VCF_FILNAME_IN -l MIN_READ_DEPTH -d MAX_MISSING_GENO [-o VCF_FILENAME_OUT]

Arguments:

-i|--vcfin FILE			- VCF_IN filename

-l|--readdepthlimit_min		- Set to missing genotype ("./.") all genotype having a total read depth below this parameter

-m|--readdepthlimit_max         - Set to missing genotype ("./.") all genotype having a total read depth above this parameter

-n|--readdepthmin_allele         - Set to missing genotype ("./.") all genotype having an allele read depth below this parameter

-d|--maxindivwithmissingdata	- Set the maximum number of individuals that can have a missing genotype for a giver variant site

-o|--vcfout			- VCF_OUT filename

--help				- This helpful help screen

EOF

exit 1;
}



