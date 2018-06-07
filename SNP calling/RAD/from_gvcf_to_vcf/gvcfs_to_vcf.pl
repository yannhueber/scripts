#! /usr/bin/perl -w
#$ -q bigmem.q
#$ -cwd
#$ -V
#$ -S /usr/bin/perl
#$ -pe parallel_fill 32

#charge modules
#use "bioinfo/bwa/0.7.12";
#use "system/java/jre8";
#use "bioinfo/FastQC/0.11.3";
#use "compiler/gcc/4.9.2";
#use "bioinfo/bwa/0.7.12";

#my $gatk_path = "/usr/local/bioinfo/GenomeAnalysisTK/3.4-46";
my $gatk_path = "/usr/local/bioinfo/GenomeAnalysisTK/3.7-0";
my $java_path = "/usr/local/java/jre8/bin";
my $picard_path = "/usr/local/bioinfo/picard-tools/1.130";
my $samtools_path = "/usr/local/bioinfo/samtools/1.2/bin";

###################
#
# Licencied under CeCill-C (http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html)
#
# Intellectual property belongs to Bioversity
#
# Written by Yann Hueber 
#
#################### 


use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use File::Path;


############### MAIN ##############


if (@ARGV < 4)
{
        print "\n\nUsage: perl $0 -r reference_fasta -p output_prefix -x extension_file_to_treat\n\n\n";
        exit 1;
}


#Variables
my ($extension, $reference, $output_prefix);

#Getoptions
GetOptions("x|extension=s" => \$extension, "r|reference=s" => \$reference, "p|outputprefix=s" => \$output_prefix) or die ("Error in command line arguments\n");


#List of all file in directory with specific extension
opendir(CWD,".") or die("Cannot open current directory\n");
my @files = grep { m/\.$extension$/ } readdir(CWD);
closedir(CWD);


#Concat bam filenames
my $igvcfs = "";
foreach (@files) {
	$igvcfs .= "-V $_ ";
}


#Programs 

#CombineGVCFs
my $combinegvcfs = "$java_path/java -Xmx128g -XX:-UseGCOverheadLimit -jar $gatk_path/GenomeAnalysisTK.jar -T GenotypeGVCFs -R $reference -nt 32 $igvcfs -drf DuplicateRead -o $output_prefix.vcf";
print ("Call SNPs with GenotypeGVCFs on files $igvcfs...\n");
system("$combinegvcfs");

#VariantFiltration (Filter SNPs before applying BaseRecalibrator)
#my $variant_filtration  = "$java_path/java -Xmx4g -jar $gatk_path/GenomeAnalysisTK.jar -T VariantFiltration -R $reference -V $output_prefix.vcf --clusterSize 3 --clusterWindowSize 10 --filterExpression 'QD < 1.5' --maskExtension 0 --filterName QD_FILTER -maskName Mask -o $output_prefix.VF.vcf";
#print ("Variant filtration --filterExpression 'MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)' --filterExpression 'QD < 1.5' on $output_prefix.vcf...\n");
#system("$variant_filtration");

#SelectVariant to filter out filtered variants
#my $select_variant = "$java_path/java -Xmx4g -jar $gatk_path/GenomeAnalysisTK.jar -T SelectVariants -R $reference -V $output_prefix.VF.vcf -select 'vc.isNotFiltered()' -o $output_prefix.VF.filtered.vcf";
#print ("Filter out filtered SNPs with SelectVariants -select 'vc.isNotFiltered()' -selectType SNP on file $output_prefix.VF.vcf...\n");
#system("$select_variant");

print "Finished\n";

exit 1;
