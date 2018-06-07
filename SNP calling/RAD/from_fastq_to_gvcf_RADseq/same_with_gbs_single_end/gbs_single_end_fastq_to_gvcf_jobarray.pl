#! /usr/bin/perl -w
#$ -q normal.q
#$ -cwd
#$ -V
#$ -t 1-106
#$ -S /usr/bin/perl
#$ -pe parallel_fill 1

#charge modules
#use "bioinfo/bwa/0.7.12";
#use "system/java/jre8";
#use "bioinfo/FastQC/0.11.3";
#use "compiler/gcc/4.9.2";
#use "bioinfo/bwa/0.7.12";

my $gatk_path = "/usr/local/bioinfo/GenomeAnalysisTK/3.4-46";
my $java_path = "/usr/local/java/jre8/bin";
my $bwa_path = "/usr/local/bioinfo/bwa/0.7.12";
my $picard_path = "/usr/local/bioinfo/picard-tools/1.130";
my $fastqc_path = "/usr/local/bioinfo/FastQC/0.11.3";
my $cutadapt_path = "/usr/local/bioinfo/python/2.7.9_build2/bin";
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


if (@ARGV < 6)
{
        print "\n\nUsage: perl $0 -r reference_file_path -x extension_file_to_treat -cu cultivar\n\n\n";
        exit 1;
}


#Variables
my ($extension, $reference, $cultivar);

#Getoptions
GetOptions("x|extension=s" => \$extension, "r|reference=s" => \$reference, "c|cultivar=s" => \$cultivar) or die ("Error in command line arguments\n");



#List of all file in directory with specific extension
opendir(CWD,".") or die("Cannot open current directory\n");  
my @files = grep { m/\.$extension$/ } readdir(CWD);
closedir(CWD);



#relative file path 
my $sge_id = `echo \$SGE_TASK_ID` or die ("Cannot catch SGE_TASK_ID\n");
chomp($sge_id);
$sge_id--; #perl is in basis 0
my $file = $files[$sge_id];
my $filenm_root = fileparse($file, qr/.$extension$/);



########################################################################################################################################################
########################################################################################################################################################
#Programs 


#fastqc
my $fastqc = "$fastqc_path/fastqc -j $java_path/java $file";
print ("fastqc on $file (one of two)\n");
system("$fastqc");

#Cutadapt
my $rm_adapter = "$cutadapt_path/cutadapt -f fastq -g TGGCACAGA -b AGATCGGAAGAGC -O 7 -m 30 --quality-base=33 -q 20,20 $file | gzip > $filenm_root.cutadapt.fastq.gz";
print ("cutadapt on $file\n");
system("$rm_adapter");

#fastqc again
$fastqc = "$fastqc_path/fastqc -j $java_path/java $filenm_root.cutadapt.fastq.gz";
print ("fastqc on $filenm_root.cutadapt.fastq.gz (two of two)\n");
system("$fastqc");

#gunzip
my $gunzip = "gunzip $filenm_root.cutadapt.fastq.gz";
print ("gunzip $filenm_root.cutadapt.fastq.gz\n");
system("$gunzip");

#Check homogenieity of fastq files forward and reverse
#my $compare_fastq = "perl compare_fastq_paired_v5.pl -f $forward_without_ext.cutadapt.fastq -r $reverse_without_ext.cutadapt.fastq -of $forward_without_ext.cutadapt.cleaned.fastq -or $reverse_without_ext.cutadapt.cleaned.fastq -os $filenm_root.cutadapt.cleaned.single";
#print ("compare fastq v5 on $forward_without_ext.cutadapt.fastq (forward) and $reverse_without_ext.cutadapt.fastq (reverse)\n");
#system("$compare_fastq");

#bwa
my $bwa_mapping = "$bwa_path/bwa mem $reference $filenm_root.cutadapt.fastq > $filenm_root.sam";
print ("Mapping with BWA on file $filenm_root.cutadapt.fastq...\n");
system("$bwa_mapping");

#Add read group, and sort bam
my $condition = $filenm_root;
my $add_rg = "$java_path/java -Xmx6g -jar $picard_path/picard.jar AddOrReplaceReadGroups I=$filenm_root.sam O=$filenm_root.RG.sorted.bam SO=coordinate RGLB=$cultivar RGPL=illumina RGPU=run RGSM=$condition RGID=$condition";
print ("Adding Read Group and sorting on file $filenm_root.sam...\n");
system("$add_rg");

#Index bam
my $index_bam = "$samtools_path/samtools index $filenm_root.RG.sorted.bam";
print ("Indexing on bam file $filenm_root.RG.sorted.bam...\n");
system("$index_bam");

#RealignerTargetCreator
my $realigner_target = "$java_path/java -Xmx6g -jar $gatk_path/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $reference -I $filenm_root.RG.sorted.bam -o $filenm_root.forIndelRealigner.intervals";
print ("Create RealignTargetCreator intervals with $filenm_root.RG.sorted.bam...\n");
system("$realigner_target");

#IndelRealigner
my $indelrealigner = "$java_path/java -Xmx6g -jar $gatk_path/GenomeAnalysisTK.jar -T IndelRealigner -R $reference -I $filenm_root.RG.sorted.bam -targetIntervals $filenm_root.forIndelRealigner.intervals -o $filenm_root.indelrealigned.RG.sorted.bam";
print ("IndelRealigning on file $filenm_root.RG.sorted.bam...\n");
system("$indelrealigner");

#HaplotypeCaller (pour BaseRecalibrator)
my $haplotype_caller = "$java_path/java -Xmx6g -jar $gatk_path/GenomeAnalysisTK.jar -T HaplotypeCaller -R $reference -I $filenm_root.indelrealigned.RG.sorted.bam -o $filenm_root.indelrealigned.RG.sorted.HC.vcf";
print ("Call SNPs with HaplotypeCaller on file $filenm_root.indelrealigned.RG.sorted.bam...\n");
system("$haplotype_caller");

#VariantAnnotator (add annotation: MappingQualityZero)
my $variant_annotation = "$java_path/java -Xmx6g -jar $gatk_path/GenomeAnalysisTK.jar -T VariantAnnotator -R $reference -V $filenm_root.indelrealigned.RG.sorted.HC.vcf -I $filenm_root.indelrealigned.RG.sorted.bam -A MappingQualityZero -o $filenm_root.indelrealigned.RG.sorted.HC.ann.vcf";
print ("Adding MappingQualityZero annotation in file $filenm_root.indelrealigned.RG.sorted.HC.vcf...\n");
system("$variant_annotation");

#VariantFiltration (Filter SNPs before applying BaseRecalibrator)
my $variant_filtration  = "$java_path/java -Xmx6g -jar $gatk_path/GenomeAnalysisTK.jar -T VariantFiltration -R $reference -V $filenm_root.indelrealigned.RG.sorted.HC.ann.vcf --clusterSize 3 --clusterWindowSize 10 --filterExpression 'MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)' --filterExpression 'QD < 1.5' --filterExpression 'DP < 15' --maskExtension 0 --filterName HARD_TO_VALIDATE --filterName QD_FILTER --filterName DP_FILTER --maskName Mask -o $filenm_root.indelrealigned.RG.sorted.HC.ann.VF.vcf";
print ("Variant filtration --filterExpression 'MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)' --filterExpression 'QD < 1.5' --filterExpression 'DP < 15' on $filenm_root.indelrealigned.RG.sorted.HC.ann.vcf...\n");
system("$variant_filtration");

#SelectVariant to filter out filtered variants
my $select_variant = "$java_path/java -Xmx6g -jar $gatk_path/GenomeAnalysisTK.jar -T SelectVariants -R $reference -V $filenm_root.indelrealigned.RG.sorted.HC.ann.VF.vcf -select 'vc.isNotFiltered()' -selectType SNP -o $filenm_root.indelrealigned.RG.sorted.HC.ann.VF.filtered.vcf";
print ("Filter out filtered SNPs with SelectVariants -select 'vc.isNotFiltered()' -selectType SNP on file $filenm_root.indelrealigned.RG.sorted.HC.ann.VF.vcf...\n");
system("$select_variant");

#BaseRecalibrator
my $base_recalibrator = "$java_path/java -Xmx6g -jar $gatk_path/GenomeAnalysisTK.jar -T BaseRecalibrator -R $reference -I $filenm_root.indelrealigned.RG.sorted.bam --knownSites $filenm_root.indelrealigned.RG.sorted.HC.ann.VF.filtered.vcf -o $filenm_root.recal_data.table";
print ("Create a recalibration table file with BaseRecalibrator with knownSites in file $filenm_root.indelrealigned.RG.sorted.HC.ann.VF.filtered.vcf with BAM file $filenm_root.indelrealigned.RG.sorted.bam...\n");
system("$base_recalibrator");

#Apply recalibration with Printreads
my $apply_recalibration = "$java_path/java -Xmx6g -jar $gatk_path/GenomeAnalysisTK.jar -T PrintReads -R $reference -I $filenm_root.indelrealigned.RG.sorted.bam -BQSR $filenm_root.recal_data.table -o $filenm_root.indelrealigned.RG.sorted.recalibrated.bam";
print ("Apply recalibration with PrintReads on file $filenm_root.indelrealigned.RG.sorted.bam...\n");
system("$apply_recalibration");

#Obtain a gVCF with HaplotypeCaller
my $haplotype_caller_gvcf = "$java_path/java -Xmx6g -jar $gatk_path/GenomeAnalysisTK.jar -T HaplotypeCaller -R $reference -I $filenm_root.indelrealigned.RG.sorted.recalibrated.bam --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -o $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.gvcf";
print ("Call SNPs with HaplotypeCaller on file $filenm_root.indelrealigned.RG.sorted.bam and output a gvcf (--emitRefConfidence GVCF)...\n");
system("$haplotype_caller_gvcf");

print "Finished\n";

exit 1;
