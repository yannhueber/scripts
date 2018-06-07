#!/usr/bin/perl -w

#Output 
#Yann Hueber, august 2014

use warnings;
use strict;
use Getopt::Long;
use Data::Dumper;
use Carp;


# Forward declarations
sub parse_command_line();
sub open_out_file($);
sub load_genes_list($);
sub runthrough_vcf_file($);
sub usage();

#Global arguments set by command line
my $vcf_in;
my $genes_list_file = "ALL";
my $file_out = "allele_count.output";
my $min_read_depth = 10;

#Gobal variables
my $filehandle;
my %genes_name;
my @individuals;
my %stats_on_indv;

################### Start of program #####################

parse_command_line();

open_out_file( $file_out );

if (! ($genes_list_file eq "ALL")) {
	load_genes_list ( $genes_list_file );
}

print "LA\n";
print Dumper(%genes_name);
print "LO\n";

runthrough_vcf_file ( $vcf_in );

close($filehandle);

#################### End of program ######################


sub uniq {
	my %seen;
	grep !$seen{$_}++, @_;
}



sub parse_command_line() {
	my $help;

	usage() if (scalar @ARGV==0);

	my $result = GetOptions ( "i|vcfin=s" => \$vcf_in,
				  "g|geneslist=s" => \$genes_list_file,
				  "o|out=s" => \$file_out,
				  "m|min_read_depth=i" => \$min_read_depth,
				  "h|help" => \$help);

	usage() if ($help);

	die "Error: vcf file not specified (use '-i or --vcfin [FILENAME]')\n" unless defined $vcf_in;

	die "Error genes list file not specified (use '-g or --geneslist [FILENAME]')\n" unless defined $genes_list_file;

	die "Error invalid minimum read depth (valid values are from above 0)\n" if ($min_read_depth < 0);

	die "Error in command line aguments\n" unless $result;
}





sub open_out_file($) {
	my $filename_out = shift or croak("Missing vcf_out file name");
	open($filehandle,">$file_out") or die("Cannot open file $file_out : $!");
}




sub load_genes_list($) {
	my $genes_list = shift or croak("Missing genes_list_file filename");
	open(GENES,"<$genes_list") or die("Cannor open file $genes_list : $!");
	
	while(my $gene = <GENES>) {
		$gene =~ s/\s+$//;                
		$genes_name{$gene} = 1;
	}
}




sub runthrough_vcf_file($) {
	my $filename_in = shift or croak("Missing vcf_in filename");
	
	open(VCF_IN,"<$filename_in") or die("Cannot open vcf file $filename_in : $!\n");
	
	croak("Bad file handle") unless(defined($filehandle));

	while (my $vcf_line = <VCF_IN>) {
		

		chomp($vcf_line);
		my @infos_line = split(m/\t/,$vcf_line);
		
		#FIndiv infos
		if ($vcf_line =~ m/^#/) {

			if ($vcf_line =~ m/^#CHROM/) {

				@individuals = @infos_line[9..$#infos_line];

				#print headers
				my $header = "GENE\tEFFECT\tCHR\tPOS\tREF\tALLELE1\tALLELE2";
				foreach(@individuals) {
					$header.= "\t$_.allele1\t$_.allele2";
				}
				$header.="\tSNPEFF\n";
				print $filehandle $header;
			}
			next;
		}
		
		#Snps
		my $infos = $infos_line[7];
		my @info_split = split(m/;/,$infos);
		my $eff_string;
		foreach(@info_split) {
			if ($_ =~ m/^ANN=/) {
				$eff_string = $_
			}
		}
		
		my $gene_match = 0;
		my @genes;
		my @effects;
		

		my @eff_string_split = split(m/,/,$eff_string);
		foreach my $effect_infos (@eff_string_split) {
			print "my effect egal $effect_infos\n";

			my @effect_infos_split = split(m/\|/,$effect_infos);

			my $effect = $effect_infos_split[1];
			my $transcript_name = $effect_infos_split[6];

			print "my effect : $effect\n";

			next if (! defined ($transcript_name));

			$transcript_name =~ s/\.\d//;

			if (! %genes_name or $genes_name{$transcript_name}) {
				$gene_match++;
				push(@effects,$effect);
				push(@genes,$transcript_name);
			}
		}
		print "gene match number : $gene_match \n" if ($gene_match > 0);

		next if ($gene_match == 0);
		

		my $chr = $infos_line[0];
		my $pos = $infos_line[1];
		my $ref = $infos_line[3];
		my $alt = $infos_line[4];
		
		print "CHR AND POSITION $chr and $pos \n";

		#next if (length($ref) > 1 or length($alt) > 1); #keep only biallelic


		my $depth_string="";
		my %polymorphic;
		my @genotypes = @infos_line[9..$#infos_line];

		my $i=0;

		foreach my $genotype_info (@genotypes) {
			my $indv_name = $individuals[$i];
	
			print "my individual is $indv_name \n";

			my @splitted = split(m/:/,$genotype_info);
			my $total_read_depth = $splitted[2];
			if ($genotype_info =~ m/\.\/\./ or $total_read_depth eq "." or $total_read_depth < $min_read_depth) {
				

				print "Case genotype unknow or inf to 10 \n";

				#stats
				if ($genotype_info =~ m/\.\/\./) {
					$stats_on_indv{$indv_name}{"unknown_genotype"}++;
				} elsif ($total_read_depth eq "." or $total_read_depth < $min_read_depth) {
					$stats_on_indv{$indv_name}{"read_depth_inf_to_min"}++;
				}
				
				$depth_string .= "NA\t";
			} else {
				
				print "Case genotype OK\n";

				$stats_on_indv{$indv_name}{"read_depth_sup_to_min"}++;

				my $genotype = $splitted[0];
				my @genotype_split = split(m/\//,$genotype);
				foreach(@genotype_split) {
					$polymorphic{$_}++;
				}

				my $alleles_read_depth = $splitted[1];
				$depth_string.= "$alleles_read_depth\t";
			}
			$i++;
		}
		chop($depth_string);




		#Only biallelic (ref is not taken into account)
		if (scalar(keys(%polymorphic)) == 2) { #biallelic only


			my $allele1;
			my $allele2;
			my $new_depth_string="";
			my @depth_string_split = split(m/\t/,$depth_string);
			

			if ($alt =~ m/,/) {
				my @alt_split = split(m/,/,$alt);
				unshift(@alt_split,$ref);
				my ($one, $two) = sort {$a <=> $b} keys(%polymorphic);
				$allele1 = $alt_split[$one];
				$allele2 = $alt_split[$two];

				foreach my $depth_string (@depth_string_split) {
					if ($depth_string eq "NA") {
						$new_depth_string .= "NA\tNA\t";
						next;
					}
					my @depth_split = split(m/,/,$depth_string);
					$new_depth_string .= $depth_split[$one] . "\t" . $depth_split[$two] . "\t";
				}
			} else {
				$allele1 = $ref;
				$allele2 = $alt;
				foreach my $depth_string (@depth_string_split) {
					if ($depth_string eq "NA") {
						$new_depth_string .= "NA\tNA\t";
						next;
					}
					my @depth_split = split(m/,/,$depth_string);
					$new_depth_string .= $depth_split[0] . "\t" . $depth_split[1] . "\t";
				}
			}
			chop($new_depth_string);

			my $gene_string="";
			my @uniq_genes = uniq(@genes);
			foreach my $gene (@uniq_genes) {
				$gene_string.=$gene . "|";
			}
			chop($gene_string);

			my $effect_string="";
			my @uniq_effects = uniq(@effects);
			foreach my $effect (@uniq_effects) {
				$effect_string.=$effect . "|";
			}
			chop($effect_string);

			my $line_to_print = "$gene_string\t$effect_string\t$chr\t$pos\t$ref\t$allele1\t$allele2\t$new_depth_string\t$eff_string\n";
			print $filehandle $line_to_print;
		}
	}
}


#print stats
foreach my $individual (keys(%stats_on_indv)) {
        print "My indiv : $individual ---> missing genotypes : " . $stats_on_indv{$individual}{"unknown_genotype"} . " ; genotypes with read depth inf to min : " . $stats_on_indv{$individual}{"read_depth_inf_to_min"} .  " ; genotypes with read depth sup to min : " . $stats_on_indv{$individual}{"read_depth_sup_to_min"} . "\n";
}





sub usage() {
print<<EOF;
Simply print allele read depth from VCF file for genes provided as input


usage: perl $0 -i VCF_FILENAME_IN [-g GENES_LIST] -m MIN_READ_DEPTH [-o ALLELE_COUNT_OUT]

Arguments:

-i|--vcfin FILENAME                   - VCF_IN filename

-g|--geneslist FILENAME		      - GENE_LIST filename (one gene's name per line) (default is all genes)

-m|--min_read_depth 		      - Set the minimum read depth allowed 

-o|--out FILENAME		      - FILE_OUT filename

--help 				      - This helpful help screen

EOF

exit 1;
}



