#!/usr/bin/perl -w

#Script written by Yann Hueber
#Take an input VCF file and print severals genotype file 

use warnings;
use strict;
use Getopt::Long;
use Data::Dumper;

#variables
my ($input_vcf, $output_prefix, $read_depth_limit, $limit_allele_minor_count, $max_indiv_nb_with_missing_data);

my %transfo_ATGC = ("A" => "1", "T" => "2", "C" => "3", "G" => "4", "D" => 5);

if (scalar(@ARGV) < 6)
{
	print "\n\nUsage : perl $0 -i vcf_in -l read_depth_limit_int -a limit_allele_minor_count -d max_indiv_nb_with_missing_data -o outputs_prefix\n\n\n";
	exit;
}

GetOptions ("i|vcf_in=s" => \$input_vcf,		#string
	    "l|readdepthlimit=s" => \$read_depth_limit,	#string
	    "a|alleleminorcountlimit=s" => \$limit_allele_minor_count, #string
	    "d|maxindivwithmissingdata=s" => \$max_indiv_nb_with_missing_data,
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
	
	my @ref_alt_list = ();
	if ($alt =~ m/,/) {
		my @alt_split = split (m/,/,$alt);
		@ref_alt_list = ($ref,$alt_split[0],$alt_split[1]);
	} else {
		@ref_alt_list = ($ref,$alt);
	}
	my $increment = 0;
	my %alleles_count;
	my ($first_allele,$second_allele);
	my $indiv_with_missing_data_count=0;

	foreach my $indiv_infos (@individuals_infos) { #a loop for each indiv
		
		my @indv_info_split = split(m/:/,$indiv_infos);
		my $genotype = $indv_info_split[0];


		my $indv = $individuals[$increment];

		
		if ($genotype =~ m/\.\/\./) {
			$indiv_with_missing_data_count++;
			$genotypes{$chr}{$pos}{$indv} = "9\t9";
		} else {
			
			
			my @genotype_split = split(m/\//,$genotype);
			$first_allele = $ref_alt_list[$genotype_split[0]];
			$second_allele = $ref_alt_list[$genotype_split[1]];
	
			#print "HERE $first_allele $second_allele\n";
	
			#transform alleles in numbers
			$first_allele = $transfo_ATGC{$first_allele};
			$second_allele = $transfo_ATGC{$second_allele};
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
		$chr_count_sites{$chr} += 2; #count number of sites times 2 (because there is 2 alleles)
		if ($chr ne "12") {
			$global_count_sites +=2;
		}
	} else {
		$print_site{$chr}{$pos} = 0;
	}
}

close(VCF_IN);



#print output
#OUT1 and OUT2 are for darwin file .var
#OUT3 and OUT4 are for structure file .txt


my $first_line_to_print = "\@DARwin 5.0 - ALLELIC - 2";
open(OUT2,">$output_prefix\_darwin_all.var") or die ("Cannot open file $output_prefix\_darwin_all.var : $!");
print OUT2 $first_line_to_print . "\n" . scalar(@individuals) . "\t" . $global_count_sites . "\n" . "Unit";
my %out2_lines;


open(OUT4,">$output_prefix\_structure_all.txt") or die ("Cannot open file $output_prefix\_structure_all.txt : $!");


open(LINK_ALL,">$output_prefix\_position_number_link_all.txt") or die ("Cannot open file $output_prefix\_position_number_link_all.txt : $!");
print LINK_ALL "NUMBER\tCHR\tPOSITION\n";


my $j=1;
foreach my $chr (sort {$a cmp $b} keys(%genotypes)) {
	
	#one file per chr
	open(OUT1,">$output_prefix\_darwin_$chr.var") or die ("Cannot open file $output_prefix\_darwin_$chr.var : $!");
	print OUT1 $first_line_to_print . "\n" . scalar(@individuals) . "\t" . $chr_count_sites{$chr} . "\n" . "Unit";

	open(OUT3,">$output_prefix\_structure_$chr.txt") or die ("Cannot open file $output_prefix\_structure_$chr.txt : $!");

	open(LINK_CHR,">$output_prefix\_position_number_link_$chr.txt") or die ("Cannot open file $output_prefix\_position_number_link_$chr.txt : $!");
	print LINK_CHR "NUMBER\tPOSITION\n";

	my $pos_ref = $genotypes{$chr};
	my %positions = %$pos_ref;
	
	#print first line
	my $i=1;
	foreach my $pos (sort {$a <=> $b} keys(%positions)) {
		if ($print_site{$chr}{$pos} == 1) {
			print LINK_CHR "$i\t$pos\n";
			print OUT1 "\t" . $i . "\t" . ++$i;
			if ($chr ne "12") {
				print LINK_ALL "$j\t$chr\t$pos\n";
				print OUT2 "\t" . $j . "\t" . ++$j;
				$j++;
			}
			$i++;
		}
	}
	print OUT1 "\n";


	#print one line per individual
	my $incr = 1;
	foreach my $indiv (@individuals) {

		print OUT1 $incr;
		print "$incr\t$indiv\n";
		print OUT3 $incr;
		foreach my $position (sort {$a <=> $b} keys(%positions)) {
			if ($print_site{$chr}{$position} == 1) {
				print OUT1 "\t" . $genotypes{$chr}{$position}{$indiv};
				print OUT3 "\t" . $genotypes{$chr}{$position}{$indiv};
				if ($chr ne "12") {
					$out2_lines{$indiv} .=  "\t" . $genotypes{$chr}{$position}{$indiv};
				}
			}
		}
		print OUT1 "\n";
		print OUT3 "\n";
		$incr++;
	}
	close(OUT1);
	close(OUT3);
	close(LINK_CHR);
}
print OUT2 "\n";

#print file for all chr but not the 12
my $i = 1;
foreach my $indiv (@individuals) {
	print OUT2 $i . $out2_lines{$indiv} . "\n";
	print OUT4 $i . $out2_lines{$indiv} . "\n";
	$i++;
}


close(OUT2);
close(OUT4);
close(LINK_ALL);
