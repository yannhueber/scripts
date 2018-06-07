#!/usr/bin/perl -w

#This script takes an input VCF file and make a flat file/database


use warnings;
use strict;
use Getopt::Long;
use Data::Dumper;
use Carp;


# Forward declarations
sub parse_command_line();
sub load_vcf_in_file($$);
sub whos_rna($);
sub usage();

#Global arguments set by command line
my $vcf_in;
my $rna_acc_name_in;
my $out = "databaselike.txt";

#Gobal variables
my @individuals;
my %rna_acc;

################### Start of program #####################

parse_command_line();

whos_rna( $rna_acc_name_in );

load_vcf_in_file ( $vcf_in, $out );

#################### End of program ######################




sub parse_command_line() {
	my $help;

	usage() if (scalar @ARGV==0);

	my $result = GetOptions ( "i|vcfin=s" => \$vcf_in,
				  "r|rna=s" => \$rna_acc_name_in,
				  "o|out=s" => \$out,
				  "h|help" => \$help);

	usage() if ($help);

	die "Error: vcf file not specified (use '-i or --vcfin [FILENAME]')\n" unless defined $vcf_in;

	 die "Error: rna filename file not specified (use '-r or --rna [FILENAME]')\n" unless defined $rna_acc_name_in;
	
	die "Error: out file not specified (use '-o or --out [FILENAME]')\n" unless defined $out;

	die "Error in command line aguments\n" unless $result;
}



sub whos_rna($) {
	
	my $rna_in = shift or croak("Missing rna");

	open(RNA_IN,"<$rna_in") or die("Cannot open vcf file $rna_in : $!\n");

	while (my $rna_filename = <RNA_IN>) {
			
		chop($rna_filename);
		$rna_acc{$rna_filename}++;
	}
	
	close(RNA_IN);
}





sub load_vcf_in_file($$) {
	my $filename_in = shift or croak("Missing vcf_in file name");
	my $filename_out = shift or croak("Cannot create out file name");

	open(VCF_IN,"<$filename_in") or die("Cannot open vcf file $filename_in : $!\n");
	open(OUT,">$filename_out") or die("Cannot open vcf file $filename_out : $!\n");

	while (my $vcf_line = <VCF_IN>) {
		

		chomp($vcf_line);
		my @infos_line = split(m/\t/,$vcf_line);
		
		next if ($vcf_line =~ m/^##/);

		#File infos
		if ($vcf_line =~ m/^#CHROM/) {
			@individuals = @infos_line[9..$#infos_line];
			next;
		}
		
		my @individuals_infos = @infos_line[9..$#infos_line];
		

		my $chr = $infos_line[0];
		my $pos = $infos_line[1];	
		my $ref = $infos_line[3];
		my $alt = $infos_line[4];
		my @alleles = split(m/,/,$alt);
		unshift (@alleles,$ref); 

		my %position_stats;
		my %geno_from_rna;
		my $i=0;
		foreach my $indv_info (@individuals_infos) {
			my @info_split = split(m/:/,$indv_info);
			my $genotype = $info_split[0];
			if ($genotype =~ m/\.\/\./) {
				$position_stats{"missing"}++;
			} else {
				$position_stats{$genotype}++;
				my $ad = $info_split[2];
				my $indiv = $individuals[$i];
				if (defined $rna_acc{$indiv} and $ad >= 10) {
					$geno_from_rna{$genotype}++;
				}
			}
			$i++;
		}

		#print position infos
		my @inmind;
		my $nbre_indv = scalar(@individuals);
		foreach my $geno (keys (%position_stats)) {
			my $seen_nb = $position_stats{$geno};
			my $seen_perc = sprintf("%.2f", $seen_nb/$nbre_indv);
			
			if ($geno ne "missing") {
				my @geno_split = split (m/\//,$geno);
				my $all1 = $geno_split[0];
				my $all2 = $geno_split[1];
				my $allele1 = $alleles[$all1];
				my $allele2 = $alleles[$all2];
				my $geno_atgc = "$allele1/$allele2";
				my $inpush = $geno_atgc . ":" . $seen_perc;
				if (defined $geno_from_rna{$geno}) {
					$inpush .= ":RNA" . $geno_from_rna{$geno};
				}
				push(@inmind,$inpush);	
			} else {
				unshift(@inmind,"MISSING:$seen_perc");
			}
		}


		print OUT "$chr\t$pos\t$ref\t$alt";
		foreach (@inmind) {
			print OUT "\t$_";
		}
		print OUT "\n";
	}
	close(VCF_IN);
	close(OUT);
}




sub usage() {
print<<EOF;

usage: perl $0 -i VCF_FILENAME_IN -o OUT_DATABASE -r RNA_acc_names

Arguments:

-i|--vcfin FILE	                      - VCF_IN filename

-o|--out FILENAME		      - OUT filename

-r|--rna FILENAME                     - RNA ids filename

--help 				      - This helpful help screen

EOF

exit 1;
}



