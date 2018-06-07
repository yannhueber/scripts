#! /usr/bin/perl -w


use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use File::Path;


#vcf file variable
my $file_to_treat = $ARGV[0];

#Usage
if (scalar(@ARGV) == 0) {
	print "\n\n\nperl $0 vcf_file\n\n\n";
	exit 1;
}


#Filename without extension
my $filename = fileparse($file_to_treat, qr/\.vcf$/);

#open file
open(VCF,"<$file_to_treat") or die("Cannot open $file_to_treat : $!");
open(PED,">$filename.ped") or die ("Cannot open $filename.ped : $!");
open(INFO,">$filename.info") or die ("Cannot open $filename.info");


my %genotypes;
my @individuals =();
my @positions = ();

while (defined (my $line = <VCF>)) {
	next if ($line =~ m/^##/);
	chomp($line);
	my @infos_line = split(m/\t/,$line);
	if ($line =~ m/^#CHROM/) {
		@individuals = @infos_line[9..$#infos_line];
		next;
	}

	my $rs = $infos_line[2];
	my $ref = $infos_line[3];
	my $alt = $infos_line[4];
	push(@positions,$rs);

	my @indiv_genotypes = @infos_line[11..$#infos_line];



	foreach (my $i = 0; $i < @individuals; $i++) {
		my $indv = $individuals[$i];
		my $indv_genotype = $indiv_genotypes[$i];
		$indv_genotype =~ s/0/$ref/g;
		$indv_genotype =~ s/1/$alt/g;
		$indv_genotype =~ s/\./N/g;
		$indv_genotype =~ s/\///;
		$genotypes{$indv}{$rs} = $indv_genotype;
	}
}


#print ped file
#header
#no header
#body
my $incr =1;
foreach my $indiv (sort {$a cmp $b } keys(%genotypes)) {
	my $ref = $genotypes{$indiv};
	my %rs = %$ref;
	my $line_to_print = "$incr\t$indiv\t0\t0\t2\taffected_status\t";	

	foreach my $rs (@positions) {
		my $indv_genotype = $genotypes{$indiv}{$rs};
		$indv_genotype =~ s/NN/00/;
		my $first_allele = substr($indv_genotype,0,1);
		my $second_allele =substr($indv_genotype,1,1);
		$line_to_print .= "$first_allele $second_allele\t";
	}
	chop($line_to_print);
	$line_to_print .= "\n";
	print PED $line_to_print;
	$incr++;
}



foreach my $rs (@positions) {
	my @rs_split = split(m/_/,$rs);
	my $position = $rs_split[1];
	print INFO "$rs\t$position\n";
}

close(HAPMAP);
close(PED);
close(INFO);

exit 1;
