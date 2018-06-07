#!/usr/bin/perl -w

#This script takes an input VCF file and discard individuals with missingness percentage superior to INT (in parameter)
#Yann Hueber, july 2014

use warnings;
use strict;
use Getopt::Long;
use Data::Dumper;
use Carp;


# Forward declarations
sub parse_command_line();
sub load_vcf_in_file($);
sub indv_with_high_missingness();
sub launch_selectvariants(@);
sub usage();

#Global arguments set by command line
my $vcf_in;
my $vcf_out = "output.vcf";
my $reference_fasta;
my $indv_out = "individuals_deleted.txt";
my $max_missingness_perc = 100;
my $ploidy;

#Gobal variables
my @individuals;
my @indv_to_del; #remove those indv from the vcf file
my $filehandle;
#my %snps;
my $count_snps = 0;
my %count_missing;


################### Start of program #####################

parse_command_line();

load_vcf_in_file ( $vcf_in );

@indv_to_del = indv_with_high_missingness();



open(TXT,">$indv_out") or die("Cannot open $indv_out : $!");
if (scalar (@indv_to_del) == 0) {
	#print ("No individuals with missingness superior to $max_missingness_perc\n");
	my $cp_vcf_in = "cp $vcf_in $vcf_out";
	system("$cp_vcf_in");
}
else {
	#print in output file which individuals will be discarded
	foreach ( @indv_to_del ) {
		print TXT "Indiv to delete $individuals[$_] \n";
	}
	#launch GATK SelectVariants
	launch_selectvariants( @indv_to_del );
}
close(TXT);

#################### End of program ######################




sub parse_command_line() {
	my $help;

	usage() if (scalar @ARGV==0);

	my $result = GetOptions ( "i|vcfin=s" => \$vcf_in,
				  "r|reference=s" => \$reference_fasta,
				  "o|vcfout=s" => \$vcf_out,
				  "v|indvout=s" => \$indv_out,
				  "m|max_missingness_percentage=i" => \$max_missingness_perc,
				  "p|ploidy=i" => \$ploidy,
				  "h|help" => \$help);

	usage() if ($help);

	die "Error: vcf file not specified (use '-i or --vcfin [FILENAME]')\n" unless defined $vcf_in;
	
	die "Error invalid maximum missingness percentage (valid values are from 0 to 100)\n" if ($max_missingness_perc < 0 or $max_missingness_perc > 100);

	die "Error invalid ploidy (valid values are 2 or 3)\n" if (! ($ploidy == 2 or $ploidy == 3));

	die "Error in command line aguments\n" unless $result;
}






sub load_vcf_in_file($) {
	my $filename_in = shift or croak("Missing vcf_in file name");
	
	open(VCF_IN,"<$filename_in") or die("Cannot open vcf file $filename_in : $!\n");

	while (my $vcf_line = <VCF_IN>) {
		

		chomp($vcf_line);
		my @infos_line = split(m/\t/,$vcf_line);
		
		#File infos
		if ($vcf_line =~ m/^#CHROM/) {
			@individuals = @infos_line[9..$#infos_line];
			next;
		}
		
		#Snps
		$count_snps++;
		
		my @individuals_infos = @infos_line[9..$#infos_line];
		
		#Count missing per indv
		my $i=0;
		foreach my $indv_info (@individuals_infos) {
			my @info_split = split(m/:/,$indv_info);
			my $genotype = $info_split[0];
			if ($genotype =~ m/\.\/\./) {
				$count_missing{$individuals[$i]}++;
			}
			$i++;
		}
	}
	close(VCF_IN);
}



sub indv_with_high_missingness() {
	
	my @indv_to_del;

	my $num_indiv=0;
	foreach my $indv (@individuals) {
		my $perc_miss = ($count_missing{$indv} / $count_snps) * 100;
		if ($perc_miss > $max_missingness_perc) {
			push(@indv_to_del,$num_indiv);
		}
		$num_indiv++;
	}
	return @indv_to_del;
}



sub launch_selectvariants(@) {

	my @indv_to_delete = @_;
	my $indvs_string = "";

	foreach (@indv_to_delete) {
		my $indiv = $individuals[$_];
		$indvs_string.= "--exclude_sample_name $indiv ";
	}

	#launch GATK
	my $gakt_cmd = "/usr/local/java/jre8/bin/java -Xmx4G -jar /usr/local/bioinfo/GenomeAnalysisTK/3.4-46/GenomeAnalysisTK.jar -T SelectVariants -R $reference_fasta $indvs_string --excludeNonVariants -V $vcf_in -o $vcf_out";
	print ("/usr/local/java/jre8/bin/java -Xmx4G -jar /usr/local/bioinfo/GenomeAnalysisTK/3.4-46/GenomeAnalysisTK.jar -T SelectVariants -R $reference_fasta $indvs_string --excludeNonVariants -V $vcf_in -o $vcf_out");
	system("$gakt_cmd");
}




sub usage() {
print<<EOF;

This program reads a VCF file and remove individuals based 
on missingness provided in the parameters

usage: perl $0 -i VCF_FILENAME_IN -r REF_FASTA_FILENAME -m MISSINGNESS -v INDV_OUT_LIST_FILENAME -p PLOIDY [-o VCF_FILENAME_OUT]

Arguments:

-i|--vcfin FILE	                      - VCF_IN filename

-r|--reference FILE		      - FASTA_REF filename

-m|--max_missingness_percentage INT   - Set the maximum missingness percentage allowed to keep individuals

-o|--vcfout FILENAME		      - VCF_OUT filename

-v|--indvout FILENAME                 - INDV_OUT filename

-p|--ploidy INT			      - Set the ploidy (2 or 3)

--help 				      - This helpful help screen

EOF

exit 1;
}



