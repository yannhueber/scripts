#! /usr/bin/perl -w
#$ -q bigmem.q
#$ -cwd
#$ -V
#$ -S /usr/bin/perl
#$ -pe parallel_fill 4

#charge modules
#use "bioinfo/bwa/0.7.12";
#use "system/java/jre8";
#use "bioinfo/FastQC/0.11.3";
#use "compiler/gcc/4.9.2";
#use "bioinfo/bwa/0.7.12";

my $gatk_path = "/usr/local/bioinfo/GenomeAnalysisTK/3.4-46";
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
        print "\n\nUsage: perl $0 -r reference_fasta -p output_prefix -c HC|UG -x extension_file_to_treat\n\n\n";
        exit 1;
}


#Variables
my ($extension, $reference, $snp_caller, $output_prefix);

#Getoptions
GetOptions("x|extension=s" => \$extension, "r|reference=s" => \$reference, "c|snpcaller=s" => \$snp_caller, "p|outputprefix=s" => \$output_prefix) or die ("Error in command line arguments\n");


#List of all file in directory with specific extension
opendir(CWD,".") or die("Cannot open current directory\n");
my @files = grep { m/\.$extension$/ } readdir(CWD);
closedir(CWD);


#Concat bam filenames
my $ibams = "";
foreach (@files) {
	$ibams .= "-I $_ ";
}


#Programs 

#HaplotypeCaller or UnifiedGenotyper
if ($snp_caller eq "HC") {
	my $hc_caller = "$java_path/java -Xmx20g -XX:ParallelGCThreads=4 -jar $gatk_path/GenomeAnalysisTK.jar -T HaplotypeCaller -nct 4 -R $reference $ibams -ploidy 3 -dontUseSoftClippedBases -dt NONE -drf DuplicateRead -stand_call_conf 20.0 -stand_emit_conf 20.0 -o $output_prefix.$snp_caller.vcf";
	print ("Call SNPs with HaplotypeCaller on files $ibams...\n");
	system("$hc_caller");
} elsif ($snp_caller eq "UG") {
	my $ug_caller = "$java_path/java -Xmx20g -XX:ParallelGCThreads=4 -jar $gatk_path/GenomeAnalysisTK.jar -T UnifiedGenotyper -nct 4 -R $reference $ibams -glm BOTH -ploidy 3 -dt NONE -drf DuplicateRead -stand_call_conf 20.0 -stand_emit_conf 20.0 -o $output_prefix.$snp_caller.vcf";
        print ("Call SNPs with UnifiedGenotyper on files $ibams...\n");
        system("$ug_caller");
} else {
	die("Please specify a variant caller in the command line, '-c HC' for HaplotypeCaller or '-c UG' for UnifiedGenotyper\n");
}

#VariantAnnotator (add annotation: MappingQualityZero)
my $variant_annotation = "$java_path/java -Xmx5g -jar $gatk_path/GenomeAnalysisTK.jar -T VariantAnnotator -R $reference -V $output_prefix.$snp_caller.vcf $ibams -A MappingQualityZero -o $output_prefix.$snp_caller.ann.vcf";
print ("Adding MappingQualityZero annotation in file $output_prefix.$snp_caller.vcf...\n");
system("$variant_annotation");

#VariantFiltration (Filter SNPs before applying BaseRecalibrator)
my $variant_filtration  = "$java_path/java -Xmx5g -jar $gatk_path/GenomeAnalysisTK.jar -T VariantFiltration -R $reference -V $output_prefix.$snp_caller.ann.vcf --clusterSize 3 --clusterWindowSize 10 --filterExpression 'MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)' --filterExpression 'QD < 1.5' --maskExtension 0 --filterName HARD_TO_VALIDATE --filterName QD_FILTER -maskName Mask -o $output_prefix.$snp_caller.ann.VF.vcf";
print ("Variant filtration --filterExpression 'MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)' --filterExpression 'QD < 1.5' on $output_prefix.$snp_caller.ann.vcf...\n");
system("$variant_filtration");

#SelectVariant to filter out filtered variants
my $select_variant = "$java_path/java -Xmx5g -jar $gatk_path/GenomeAnalysisTK.jar -T SelectVariants -R $reference -V $output_prefix.$snp_caller.ann.VF.vcf -select 'vc.isNotFiltered()' -o $output_prefix.$snp_caller.ann.VF.filtered.vcf";
print ("Filter out filtered SNPs with SelectVariants -select 'vc.isNotFiltered()' -selectType SNP on file $output_prefix.$snp_caller.ann.VF.vcf...\n");
system("$select_variant");

print "Finished\n";

exit 1;
