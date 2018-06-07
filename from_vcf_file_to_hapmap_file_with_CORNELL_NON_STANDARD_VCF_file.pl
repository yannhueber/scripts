#!/usr/bin/perl -w

#Script written by Yann Hueber
#Take an input VCF file and print a hapmap file
#Check IUPAC hash to see how genotypes are coded

use warnings;
use strict;
use Getopt::Long;

#variables
my ($input_vcf, $out_hapmap, $read_depth_limit);


#Modified IUPAC Code
my %IUPAC = ("AA" => "A",
	     "TT" => "T",
	     "CC" => "C",
	     "GG" => "G",
	     "AG" => "R",
	     "GA" => "R",
	     "CT" => "Y",
	     "TC" => "Y",
	     "AC" => "M",
	     "CA" => "M",
	     "TG" => "K",
	     "GT" => "K",
	     "AT" => "W",
	     "TA" => "W",
	     "CG" => "S",
	     "GC" => "S",
	     "A-" => "W",
	     "-A" => "W",
	     "T-" => "W",
	     "-T" => "W",
	     "G-" => "S",
	     "-G" => "S",
	     "C-" => "S",
	     "-C" => "S");

my %compl = ("A" => "T",
	     "T" => "A",
	     "C" => "G",
	     "G" => "C");


if (scalar(@ARGV) < 6) {
	print "\n\nUsage : perl from_vcf_file_to_hapmap_file_with_CORNELL_NON_STANDARD_VCF_file.pl -i vcf_in -l read_depth_limit_int -o hapmap_output\n\n\n";
	exit;
}

GetOptions ("i|vcf_in=s" => \$input_vcf,		#string
	    "l|readdepthlimit=s" => \$read_depth_limit,	#string
	    "o|out=s" => \$out_hapmap)			#string
or die ("Error in command line arguments\n");



#Open files
open(VCF_IN,"<$input_vcf") or die("Cannot open file $input_vcf : $!");
open(HAP,">$out_hapmap") or die("Cannot open file $out_hapmap : $!");

open(DEL,">how_deletion_was_coded.txt") or die ("Cannot open file how_deletion_was_coded.txt : $!");


#head line
print HAP "rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode";

my @individuals;

#Running through vcf file
while (defined (my $vcf_line = <VCF_IN>)) {
	next if ($vcf_line =~ m/^##/);
	chomp($vcf_line);
	my @infos_line = split(m/\t/,$vcf_line);
	my @individuals_infos = @infos_line[9..$#infos_line];

	#Print individuals name in hapmap file
	my $individuals_string="";
	if ($vcf_line =~ m/^#CHROM/) {
		@individuals = @individuals_infos;
		foreach (@individuals) {
			$individuals_string .= "\t$_";
		}
		print HAP $individuals_string . "\n";
		next;
	}
	
	#For one site
	my $chr = $infos_line[0];
	my $pos = $infos_line[1];
	my $major = $infos_line[3];
	my $minor = $infos_line[4];
	my @major_minor_list = ($major,$minor);

	#How deletion was coded --> print in file
	my $complementary;
	my $alleles="";
	if ($major eq "-" or $minor eq "-") {
		my $i=0;
		my $where;
		foreach (@major_minor_list) {
			if ($_ eq "-") {
				$where = $i;
				next;
			}
			$complementary = $compl{$_};
			print DEL "$chr\t$pos\t" . $complementary . "\n";
			$i++;
		}
		if ($where == 0) {
			$alleles .= $complementary . "/" .  $compl{$complementary};
		} else {
			$alleles .= $compl{$complementary} . "/" . $complementary;
		}
	}

	#line_begins
	my $rs = "S" . $chr . "_" . $pos;
	if (! $alleles) {
		$alleles = $major . "/" . $minor;
	}
	my $strand = "+";
	my $hapmap_line = "$rs\t$alleles\t$chr\t$pos\t$strand\tNA\tNA\tNA\tNA\tNA\tNA";

	foreach my $indiv_infos (@individuals_infos) { #a loop for each indiv
		my @indiv_infos_split = split(m/:/,$indiv_infos);
		my $genotype = $indiv_infos_split[0];
		my $read_depth = $indiv_infos_split[2];

		my $hapmap_genotype;

		if ($read_depth < $read_depth_limit) {
			$hapmap_genotype = "N";
		} else {
			my @genotype_split = split(m/\//,$genotype);
			my $genotype_string = $major_minor_list[$genotype_split[0]] . $major_minor_list[$genotype_split[1]];
			if ($genotype_string eq "--") {
				$genotype_string = $complementary . $complementary;
			}
			$hapmap_genotype = $IUPAC{$genotype_string};
		}
		$hapmap_line .= "\t$hapmap_genotype";
	}

	#print line 
	print HAP $hapmap_line . "\n";
}

close(VCF_IN);
close(HAP);
close(DEL);
