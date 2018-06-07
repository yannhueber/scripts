#!/usr/bin/perl -w

#This script takes an input VCF file and two database (parents) and look for homeosnp alleles

use warnings;
use strict;
use Getopt::Long;
use Data::Dumper;
use Carp;


# Forward declarations
sub parse_command_line();
sub load_database($);
sub runthrough_vcf_in_file($$$$);
sub seek_alleles($$);
sub usage();


#Global arguments set by command line
my $vcf_in;
my $db_a;
my $db_b;
my $out="homeosnpallele.txt";
my $maf=0.05;


#Gobal variables
my @individuals;


################### Start of program #####################

parse_command_line();

my $hashref_first = load_database ( $db_a );

my $hashref_second = load_database ( $db_b );

runthrough_vcf_in_file ( $vcf_in, $out, $hashref_first, $hashref_second );

#################### End of program ######################




sub parse_command_line() {
	my $help;

	usage() if (scalar @ARGV==0);

	my $result = GetOptions ( "i|vcfin=s" => \$vcf_in,
				  "a|dba=s" => \$db_a,
				  "b|dbb=s" => \$db_b,
				  "m|maf=f" => \$maf,
				  "o|out=s" => \$out,
				  "h|help" => \$help);

	usage() if ($help);

	die "Error: vcf file not specified (use '-i or --vcfin [FILENAME]')\n" unless defined $vcf_in;
	
	die "Error: first database file not specified (use '-a or --dba [FILENAME]')\n" unless defined $db_a;
		
	die "Error: second database file not specified (use '-b or --dbb [FILENAME]')\n" unless defined $db_b;

	die "Error: maf not specified or not comprised between 0 and 1 (use '-m or --maf FLOAT')\n" unless $maf >= 0 and $maf <= 1;

	die "Error: out file not specified (use '-o or --out [FILENAME]')\n" unless defined $out;

	die "Error in command line aguments\n" unless $result;
}



sub load_database($) {
	my $db_file = shift or croak("Missing db file");
		
	open(DB,"<$db_file") or die("Cannot open db_file $db_file : $!\n");
	
	my %db_hash;
	
	while (my $db_line = <DB>) {
		chomp($db_line);
		my @infos = split(m/\t/,$db_line);
		my $chr = shift(@infos);
		my $pos = shift(@infos);
		my $chrpos = "$chr" . "-" . "$pos";
		shift(@infos);
		shift(@infos);
		foreach my $geno_perc (@infos) {
			push(@{$db_hash{$chrpos}}, $geno_perc);
		}
	}
		
	return \%db_hash;

}


sub runthrough_vcf_in_file($$$$) {
	
	my $filename_in = shift or croak("Missing input filename");
	my $filename_out = shift or croak("Missing output filename");
	my $dba = shift;
	my $dbb = shift;

	open(VCF_IN,"<$filename_in") or die("Cannot open vcf file $filename_in : $!\n");
	open(OUT,">$filename_out") or die("Cannot open vcf file $filename_out : $!\n");


	#print results
	print OUT "CHR\tPOS\tREF\tALT";


	while (my $vcf_line = <VCF_IN>) {
		
		chomp($vcf_line);
		my @infos_line = split(m/\t/,$vcf_line);
		
		next if ($vcf_line =~ m/^##/);

		#File infos
		if ($vcf_line =~ m/^#CHROM/) {
			@individuals = @infos_line[9..$#infos_line];
			foreach(@individuals) {
                        	print OUT "\t$_";
                	}
                	print OUT "\tA_ALLELES\tB_ALLELES\tCASE\n";
			next;
		}
		
		my @individuals_infos = @infos_line[9..$#infos_line];


		my $chr = $infos_line[0];
		my $pos = $infos_line[1];
		my $chrpos = $chr . "-" . $pos;	
		my $ref = $infos_line[3];
		my $alt = $infos_line[4];
		next if (length($ref) > 1 or length($alt) > 1); #biallelic
		my @alleles = split(m/,/,$alt);
		unshift (@alleles,$ref); 
		
		
		

		#Check case first (variant in dbA? variant in dbB?)
		my $case=0;
		my @common_a_b;

		my %a_alleles = seek_alleles($dba, $chrpos);
		my %b_alleles = seek_alleles($dbb, $chrpos);

		if (! %a_alleles) {
			$a_alleles{$ref}++;
		}
		if (! %b_alleles) {
			$b_alleles{$ref}++;
		}
		

		foreach my $a_allele (keys(%a_alleles)) {
			if (defined ($b_alleles{$a_allele})) {
				push(@common_a_b, $a_allele);
			}
		}

		if (scalar(@common_a_b) == 0) {
			if (scalar(keys (%a_alleles)) == 1 and scalar(keys (%b_alleles)) == 1) {
				$case =1;
			} elsif ( scalar(keys (%a_alleles)) == 1 and scalar(keys (%b_alleles)) > 1) {  
				$case =6; 
			} elsif ( scalar(keys (%a_alleles)) > 1 and scalar(keys (%b_alleles)) == 1) {
				$case =5;
			} elsif ( scalar(keys (%a_alleles)) > 1 and scalar(keys (%b_alleles)) > 1) {
				$case =7;
			}
		}
		
		next if ($case == 0); 

			
		print OUT "$chr\t$pos\t$ref\t$alt\t";



		my $indv_info_to_print="";
		foreach my $acc_info (@individuals_infos) {
			my %indv_geno_ad;
			my @acc_info_split = split(m/:/,$acc_info);
			my $genotype = $acc_info_split[0];
			if ($genotype =~ m/\.\/\./) {
				$indv_info_to_print .= "NA:NA\t";
			} else {
				my $ad = $acc_info_split[1];
				my @ad_split = split(m/,/,$ad);
				my @geno_split = split(m/\//,$genotype);
				foreach my $all (@geno_split) {
					my $atgc_allele = $alleles[$all];
					if (! defined $indv_geno_ad{$atgc_allele}) {
						my $allele_ad = $ad_split[$all];
						if (! defined $allele_ad) {
							$allele_ad = "NA";
						}
						$indv_info_to_print .=  $atgc_allele . ":" . $allele_ad . ";"; 
					}	$indv_geno_ad{$atgc_allele}++;
				}
				chop($indv_info_to_print);
				$indv_info_to_print .= "\t";
			}
		}

		chop($indv_info_to_print);
		print OUT  $indv_info_to_print;
		my $alleles_to_print="";
		foreach my $a_all (keys(%a_alleles)) {
			$alleles_to_print.=  $a_all . ","; 
		}
		chop($alleles_to_print);
		$alleles_to_print.= "\t";
		foreach my $b_all (keys(%b_alleles)) {
                        $alleles_to_print.=  $b_all . ",";
                }
		chop($alleles_to_print);

		print OUT "\t" . $alleles_to_print . "\t" . $case . "\n";
		
	}
	
	close(VCF_IN);
	close(OUT);
}




sub seek_alleles($$) {
	
	my $db = shift;
	my $chrpos = shift;

	my %alleles=();
	
	if (defined $db->{$chrpos}) {

		my $list_ref = $db->{$chrpos};
		foreach my $geno_perc (@$list_ref) {
			my @geno_perc_infos = split(m/:/,$geno_perc);
			my $geno = $geno_perc_infos[0];
			my $perc = $geno_perc_infos[1];
			my $rna_associated=0;
			if ($geno_perc =~ m/RNA/) {
				$rna_associated=1;
			}
			next if ($geno eq "MISSING");
			if ($perc >= $maf or $rna_associated) {
				my $all1 = substr($geno,0,1);
				my $all2 = substr($geno,-1,1);
				$alleles{$all1}++;
				$alleles{$all2}++;
			}
		}
	}
	return(%alleles);
}











sub usage() {
print<<EOF;

usage: perl $0 -i VCF_FILENAME_IN -a DB_A -b DB_B -m MAF-o OUT_FILENAME 

Arguments:

-i|--vcfin FILE	                      - VCF_IN filename

-a|--dba FILE	                      - DB_IN filename

-b|--dbb FILE	                      - DB_IN filename

-m|--maf FLOAT	                      - Minor allele frequency

-o|--out FILENAME		      		  - OUT filename

--help 				    			  - This helpful help screen

EOF

exit 1;
}



