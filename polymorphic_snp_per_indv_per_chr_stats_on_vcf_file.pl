#Script written by Yann Hueber
#Take an input VCF file and a list of individuals (as defined in the VCF headerline) on a seperate line in an other file and output a new VCF recoded with only indivduals mentionned

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

#variables
my @indiv_list = ();
my $total_site_count = 0;
my %count_polymorphic_per_indv_per_chr;
my $end_list_number;


#Running through vcf
while (my $vcf_line = <VCF>) {
	next if ($vcf_line =~ m/^##/);
	chomp($vcf_line);
	my @infos_line = split(m/\t/,$vcf_line);

	#indiv line
	if ($vcf_line =~ m/^#CHROM/) {
		@indiv_list = @infos_line;
		$end_list_number = $#infos_line;
	} else { #snps
		my $chr = $infos_line[0];
		$total_site_count++;

		for (my $i=9; $i<=$end_list_number; $i++) {
			
			if (! defined($count_polymorphic_per_indv_per_chr{$indiv_list[$i]}{$chr})) {
				$count_polymorphic_per_indv_per_chr{$indiv_list[$i]}{$chr} = 0;
			}
			my $genotype_infos = $infos_line[$i];
			my @genotype_infos_split = split(m/:/,$genotype_infos);
			my $genotype = $genotype_infos_split[0];
			if ($genotype !~ m/(\d{1})\/\1/) {
				$count_polymorphic_per_indv_per_chr{$indiv_list[$i]}{$chr}++;
			}
		}
	}
}

close(VCF);


open(OUT,">$output") or die("Cannot open $output : $!");
print OUT "INDIV\tCHR\tPOLYMORPHIC_SNP_PERCENT\n";


#print OUTPUT
foreach my $indiv (sort {$a cmp $b} keys(%count_polymorphic_per_indv_per_chr)) {
	my $chr_ref = $count_polymorphic_per_indv_per_chr{$indiv};
	my %chr_hash = %$chr_ref;
	foreach my $chr(sort {$a cmp $b} keys (%chr_hash)) {
		my $polymorphic_snp_number = $count_polymorphic_per_indv_per_chr{$indiv}{$chr};
		my $polymorphic_snp_percent = $polymorphic_snp_number/$total_site_count*100;
		$polymorphic_snp_percent = sprintf("%.2f", $polymorphic_snp_percent);
		print OUT $indiv . "\t" . $chr . "\t" . $polymorphic_snp_percent . "\n";
	}
}


close(OUT);
