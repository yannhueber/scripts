#!/usr/bin/perl -w
#$ -q normal.q
#$ -cwd
#$ -V
#$ -t 1-1
#$ -S /usr/bin/perl
#$ -pe parallel_fill 1

my $java_path = "/usr/local/java/jre8/bin";

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use File::Path;

if (@ARGV == 0)
{
        print "\n\nUsage: perl $0 -n number_of_accession -r reference_file_path -x extension_file_to_treat\n\n\n";
        exit 1;
}


#Variables
my ($extension, $reference, $nb_ind, $min_read_depth);

#Getoptions
GetOptions("x|extension=s" => \$extension, "r|ref=s" => \$reference, "n|nb=i" => \$nb_ind, "m|minrd=i" => \$min_read_depth) or die ("Error in command line arguments\n");

#List of all file in directory with specific extension
opendir(CWD,".") or die("Cannot open current directory\n");
my @files = grep { m/\.$extension$/ } readdir(CWD);
closedir(CWD);


#relative file path
my $sge_id = `echo \$SGE_TASK_ID` or die ("Cannot catch SGE_TASK_ID\n");
chomp($sge_id);
$sge_id--; #perl is in basis 0
my $filename = $files[$sge_id];
my @filenm_split = split (m/_/,$filename);
my $itc_name = $filenm_split[0];
my $filename_without_ext = fileparse($filename, qr/.$extension$/);


#Programs


#SnpEff
my $snpeff = "$java_path/java -jar ~/bin/snpEff/snpEff.jar -v macuminata_v2_04_2015 $filename > $filename_without_ext.snpeff.vcf";
print ("\nRunning Snpeff\n");
system("$snpeff");

my $an_gatk = $nb_ind*2; #*2 for diploids;
#GATK no missing
my $gatk = "$java_path/java -jar /usr/local/bioinfo/GenomeAnalysisTK/3.5-0/GenomeAnalysisTK.jar -T SelectVariants -R $reference -V $filename_without_ext.snpeff.vcf -select 'AN==$an_gatk' -o $filename_without_ext.snpeff.nomissing.vcf";
print ("\nRunning SelectVariants from GATK\n");
system("$gatk");


#SnpSift filtering
my $snpsift = "cat $filename_without_ext.snpeff.nomissing.vcf | $java_path/java -jar /home/hueber/bin/snpEff/SnpSift.jar filter 'isRef(GEN[ITC0084]) & isRef(GEN[ITC1586]) & isVariant(GEN[$itc_name])' > $filename_without_ext.snpeff.nomissing.AAAequaltoRefand$itc_name\_isvariant.vcf";
print ("\nRunning SnpSift\n");
system("$snpsift");


#Output allele count
my $alldepth = "perl output_allele_count_from_vcf.pl -i $filename_without_ext.snpeff.nomissing.AAAequaltoRefand$itc_name\_isvariant.vcf -m $min_read_depth -o $filename_without_ext.allelic_depth_on_all_genes";
print ("\nRunning output_allele_count\n");
system("$alldepth");



exit 1;


