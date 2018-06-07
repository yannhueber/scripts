#!/usr/bin/perl -w

#Script written by Yann Hueber
#Take an input VCF file (from CORNELL)and print fastPHASE input file

use warnings;
use strict;
use Getopt::Long;
use Data::Dumper;

#variables
my ($input_vcf, $output_prefix, $read_depth_limit, $limit_allele_minor_count, $max_indiv_nb_with_missing_data);


if (scalar(@ARGV) < 6)
{
	print "\n\nUsage : perl $0 -i vcf_in -l read_depth_limit_int -a limit_allele_minor_count -d max_indiv_nb_with_missing_data -o outputs_prefix\n\n\n";
	exit;
}

GetOptions ("i|vcf_in=s" => \$input_vcf,		#string
	    "l|readdepthlimit=s" => \$read_depth_limit,	#string
	    "a|alleleminorcountlimit=s" => \$limit_allele_minor_count, #string
	    "d|maxindivwithmissingdata=s" => \$max_indiv_nb_with_missing_data, #string
	    "o|out=s" => \$output_prefix)		#string
or die ("Error in command line arguments\n");



#Open vcf file
open(VCF_IN,"<$input_vcf") or die("Cannot open file $input_vcf : $!");



#Running through vcf file
my @individuals = ();
my %genotypes;
my %print_site;
my $global_count_sites;
my %chr_count_sites = ("1" => 0, "2" => 0, "3" => 0, "4" => 0, "5" => 0, "6" => 0, "7" => 0, "8" => 0, "9" => 0, "10" => 0, "11" => 0, "12" => 0);


while (defined (my $vcf_line = <VCF_IN>)) {
	next if ($vcf_line =~ m/^##/);
	chomp($vcf_line);
	my @infos_line = split(m/\t/,$vcf_line);
	my @individuals_infos = @infos_line[9..$#infos_line];

	#Keep in mind individuals name
	if ($vcf_line =~ m/^#CHROM/) {
		@individuals = @individuals_infos;
		next;
	}
	

	#For one site
	my $chr = $infos_line[0];
	my $pos = $infos_line[1];
	my $ref = $infos_line[3];
	my $alt = $infos_line[4];
	
	#replace - by D for deletion
	$ref =~ s/-/D/;
	$alt =~ s/-/D/;
	my @ref_alt_list = ($ref,$alt);
	my $increment = 0;
	my %alleles_count;
	my ($first_allele,$second_allele);
	my $indiv_with_missing_data_count=0;

	foreach my $indiv_infos (@individuals_infos) { #a loop for each indiv
		my @indiv_infos_split = split(m/:/,$indiv_infos);
		my $genotype = $indiv_infos_split[0];
		my $total_read_depth = $indiv_infos_split[2];
		my $indv = $individuals[$increment];

		if ($total_read_depth < $read_depth_limit) {
			$indiv_with_missing_data_count++;
			$genotypes{$chr}{$pos}{$indv} = "?\t?";
		} else {
			my @genotype_split = split(m/\//,$genotype);
			$first_allele = $ref_alt_list[$genotype_split[0]];
			$second_allele = $ref_alt_list[$genotype_split[1]];
			$genotypes{$chr}{$pos}{$indv} = "$first_allele\t$second_allele";
			$alleles_count{$first_allele}++;
			$alleles_count{$second_allele}++;
		}
		$increment++;
	}
	
	#We print site only if the minor allele is present at least : limit_allele_minor_count
	my $bool = 1;
	my $number_of_alleles=0;
	foreach my $allele (keys(%alleles_count)) {
		$number_of_alleles++;
		if ($alleles_count{$allele} < $limit_allele_minor_count) {
			$bool = 0;
		}
	}
	if ($bool and $number_of_alleles > 1 and $indiv_with_missing_data_count <= $max_indiv_nb_with_missing_data) {
		$print_site{$chr}{$pos} = 1;
		$chr_count_sites{$chr} ++; #count number of sites times 2 (because there is 2 alleles)
		if ($chr ne "12") {
			$global_count_sites ++;
		}
	} else {
		$print_site{$chr}{$pos} = 0;
	}
}

close(VCF_IN);



#print output

open(OUT2,">$output_prefix\_fastphase_all.inp") or die ("Cannot open file $output_prefix\_fastphase_all.inp : $!");
print OUT2 scalar(@individuals) . "\n" . $global_count_sites . "\n" . "P";
my %out2_lines;



foreach my $chr (sort {$a <=> $b} keys(%genotypes)) {
	
	#one file per chr
	open(OUT1,">$output_prefix\_fastphase_$chr.inp") or die ("Cannot open file $output_prefix\_fastphase_$chr.inp : $!");
	print OUT1 scalar(@individuals) . "\n" . $chr_count_sites{$chr} . "\n" . "P";


	my $pos_ref = $genotypes{$chr};
	my %positions = %$pos_ref;
	
	#print first line
	foreach my $pos (sort {$a <=> $b} keys(%positions)) {
		if ($print_site{$chr}{$pos} == 1) {
			print OUT1 "\t" . $pos;
			if ($chr ne "12") {
				print OUT2 "\t" . $chr. "." . $pos;
			}
		}
	}
	print OUT1 "\n";


	#print three lines per individual
	foreach my $indiv (@individuals) {
		
		print OUT1 $indiv . "\n";
		my $first_genotype_line ="";
		my $second_genotype_line ="";
		foreach my $position (sort {$a <=> $b} keys(%positions)) {
			if ($print_site{$chr}{$position} == 1) {
				my @geno_split = split(m/\t/,$genotypes{$chr}{$position}{$indiv});
				$first_genotype_line .= $geno_split[0];
				$second_genotype_line .= $geno_split[1];
				if ($chr ne "12") {
					$out2_lines{$indiv}{"one"} .=  $geno_split[0];
					$out2_lines{$indiv}{"two"} .= $geno_split[1];
				}
			}
		}
		print OUT1 $first_genotype_line . "\n";
		print OUT1 $second_genotype_line. "\n";
	}
	close(OUT1);
}
print OUT2 "\n";

#print file for all chr but not the 12
foreach my $indiv (@individuals) {
	print OUT2 $indiv . "\n";
	print OUT2 $out2_lines{$indiv}{"one"} . "\n";
	print OUT2 $out2_lines{$indiv}{"two"} . "\n";
}


close(OUT2);
