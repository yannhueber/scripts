#!/usr/bin/perl -w


use warnings;
use strict;
use Getopt::Long;
use Data::Dumper;
use List::Util qw( min max );


#variables
my ($input_vcf);
my $output = "results.txt";
my $window;
my $step;

if (scalar (@ARGV)==0) {
	print ("\n\nperl $0 -i vcf_in -o result.txt -w window -s step\n\n");
	exit 1;
}


GetOptions ("i|vcfin=s" => \$input_vcf,         #string
	    "o|output=s" => \$output, 		#string
	    "s|step=i" => \$step,           #string
	    "w|window=i" => \$window)       	#string
or die ("Error in command line arguments\n");


die "Error vcf file not specified (use '-i or --vcfin [FILENAME]')\n" unless defined $input_vcf;

#Open vcf file
open(VCF_IN,"<$input_vcf") or die("Cannot open file $input_vcf : $!");
open(OUT,">$output") or die ("Cannot open file $output : $!");


my @individuals_infos;
my %mind;


#Running through vcf file
while (defined (my $vcf_line = <VCF_IN>)) {
	
	if ($vcf_line =~ m/^##/) {
		next;
	} else {
		chomp($vcf_line);
		my @infos_line = split(m/\t/,$vcf_line);

		if ($vcf_line =~ m/^#CHROM/) {
			@individuals_infos = @infos_line[9..$#infos_line];
		} else {

			my $chr = $infos_line[0];
			my $position = $infos_line[1];
			my @genotype_infos = @infos_line[9..$#infos_line];
			
			my $i = 0;
			foreach my $geno_infos (@genotype_infos) {
				my @split_geno = split(m/:/,$geno_infos);
				my $genotype = $split_geno[0];
				my $allele_depth = $split_geno[1];

				my $indv = $individuals_infos[$i];
				$mind{$indv}{$chr}{$position}{"geno"} = $genotype;
				$mind{$indv}{$chr}{$position}{"ad"} = $allele_depth;
				$i++;
			}
		}
	}
}




print OUT "INDIVIDUALS\tCHROMOSOME\tWINDOWS\tNUMBER_OF_SNPS\tHP\n";

my $rightlimit_window;




#print results

my $step_init = 0;


while ($step_init < $window) {
	foreach my $indiv (keys (%mind)) {
	
		#print "Mon indiv est : $indiv\n";
	
		my $chr_ref = $mind{$indiv};
		my %chr = %$chr_ref;

		foreach my $chrom (sort keys (%chr)) {
			
			my $position_ref = $chr{$chrom};
			my %positions = %$position_ref;
	
			$rightlimit_window=$window + $step_init;

			my @chunk =();
	
			#print "mon chromosome est le $chrom\n";
	
			foreach my $pos (sort { $a <=> $b } keys (%positions)) {

				my $genotype = $mind{$indiv}{$chrom}{$pos}{"geno"};
				my $allele_depth = $mind{$indiv}{$chrom}{$pos}{"ad"};
			

				if ($pos > $step_init and $pos <= $rightlimit_window) {
					push(@chunk,$allele_depth);
				} else {
					while ($pos > $rightlimit_window) {
						my $hp_to_print = calcul_hp(@chunk);
						print OUT "$indiv\t$chrom\t" . $hp_to_print;
						$rightlimit_window+=$window;
						@chunk =();
					}
					push(@chunk,$allele_depth);
				}
			}
			
			#print "jai fini mon chrom et je vais print le dernier chunk\n";
	
			my $hp_to_print = calcul_hp(@chunk);
			print OUT "$indiv\t$chrom\t" . $hp_to_print;
		}
	}
	$step_init += $step;
}


sub calcul_hp {
	my @chunk = @_;
	my $string_to_print = "";

	my $windowleft = $rightlimit_window - $window;

	my $num_snps = scalar (@chunk);

	my $max_sum=0;
	my $min_sum=0;
	foreach my $allele_depth (@chunk) {
		my @split_ad = split(m/,/,$allele_depth);
		my $min = min(@split_ad);
		my $max = max(@split_ad);
		$min_sum+=$min;
		$max_sum+=$max;
	}
	
	my $hp;
	if ($num_snps > 0) {
		$hp = (2*$max_sum*$min_sum) / ( ($max_sum+$min_sum)**2 );
	} else {
		$hp = "NA";	
	}

	$string_to_print .= "$windowleft-$rightlimit_window\t$num_snps\t$hp\n";

	return $string_to_print;
}


close(VCF_IN);
close(OUT);










