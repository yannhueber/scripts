#! /usr/bin/perl -w
#$ -q normal.q
#$ -cwd
#$ -V
#$ -t 1-206
#$ -S /usr/bin/perl
#$ -pe parallel_fill 1

#charge modules
use "bioinfo/bwa/0.7.12";
use "system/java/jre8";
use "bioinfo/FastQC/0.11.3";
use "compiler/gcc/4.9.2";
use "bioinfo/bwa/0.7.12";


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


my %hash_file;
my $i=0;
my $j=1;
foreach my $file (sort( {$a cmp $b} @files)) {
	if ($j%2 != 0) {
		$hash_file{$i}=$file;
	} else {
		$hash_file{$i}.= "-" . $file;
		$i++;
	}
	$j++;
}


#relative file path 
my $sge_id = `echo \$SGE_TASK_ID` or die ("Cannot catch SGE_TASK_ID\n");
chomp($sge_id);
$sge_id--; #perl is in basis 0
my $files_to_treat = $hash_file{$sge_id};
my @files_to_treat_split = split(m/-/,$files_to_treat);
my $forward = $files_to_treat_split[0];
my $forward_without_ext = fileparse($forward, qr/.$extension$/);
my $reverse = $files_to_treat_split[1];
my $reverse_without_ext = fileparse($reverse, qr/.$extension$/);


#Filename without extension
my $filenm_root = fileparse($forward, qr/_\d\.$extension$/);


#Programs 

foreach my $file ( ($forward,$reverse) ) {
	
	my $filenm_without_ext=fileparse($file, qr/.$extension$/);

	#fastqc
	my $fastqc = "fastqc -j java $file";
	print ("fastqc on $file (one of two)\n");
	system("$fastqc");

	#Cutadapt
	my $rm_adapter = "cutadapt -b AGATCGGAAGAGC -O 7 -m 30 --quality-base=64 -q 20,20 $file | gzip > $filenm_without_ext.cutadapt.fastq.gz";
	print ("cutadapt on $file\n");
	system("$rm_adapter");

	#fastqc again
	$fastqc = "fastqc -j java $filenm_without_ext.cutadapt.fastq.gz";
	print ("fastqc on $filenm_without_ext.cutadapt.fastq.gz (two of two)\n");
	system("$fastqc");

	#gunzip
	my $gunzip = "gunzip $filenm_without_ext.cutadapt.fastq.gz";
	print ("gunzip $filenm_without_ext.cutadapt.fastq.gz\n");
	system("$gunzip");
}

#Check homogenieity of fastq files forward and reverse
my $compare_fastq = "perl compare_fastq_paired_v5.pl -f $forward_without_ext.cutadapt.fastq -r $reverse_without_ext.cutadapt.fastq -of $forward_without_ext.cutadapt.cleaned.fastq -or $reverse_without_ext.cutadapt.cleaned.fastq -os $filenm_root.cutadapt.cleaned.single";
print ("compare fastq v5 on $forward_without_ext.cutadapt.fastq (forward) and $reverse_without_ext.cutadapt.fastq (reverse)\n");
system("$compare_fastq");

#bwa
my $bwa_mapping = "bwa mem $reference $forward_without_ext.cutadapt.cleaned.fastq $reverse_without_ext.cutadapt.cleaned.fastq > $filenm_root.sam";
print ("Mapping with BWA on file $forward_without_ext.cutadapt.cleaned.fastq and $reverse_without_ext.cutadapt.cleaned.fastq...\n");
system("$bwa_mapping");

#Add read group, and sort bam
my $condition = $filenm_root;
my $add_rg = "java -Xmx4g -jar picard.jar AddOrReplaceReadGroups I=$filenm_root.sam O=$filenm_root.RG.sorted.bam SO=coordinate RGLB=$cultivar RGPL=illumina RGPU=run RGSM=$condition RGID=$condition";
print ("Adding Read Group and sorting on file $filenm_root.sam...\n");
system("$add_rg");

#Index bam
my $index_bam = "samtools index $filenm_root.RG.sorted.bam";
print ("Indexing on bam file $filenm_root.RG.sorted.bam...\n");
system("$index_bam");

#RealignerTargetCreator
my $realigner_target = "java -Xmx4g -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R $reference -I $filenm_root.RG.sorted.bam --fix_misencoded_quality_scores -o $filenm_root.forIndelRealigner.intervals";
print ("Create RealignTargetCreator intervals with $filenm_root.RG.sorted.bam...\n");
system("$realigner_target");

#IndelRealigner
my $indelrealigner = "java -Xmx4g -jar GenomeAnalysisTK.jar -T IndelRealigner -R $reference -I $filenm_root.RG.sorted.bam -targetIntervals $filenm_root.forIndelRealigner.intervals --fix_misencoded_quality_scores -o $filenm_root.indelrealigned.RG.sorted.bam";
print ("IndelRealigning on file $filenm_root.RG.sorted.bam...\n");
system("$indelrealigner");

#HaplotypeCaller (pour BaseRecalibrator)
my $haplotype_caller = "java -Xmx4g -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R $reference -I $filenm_root.indelrealigned.RG.sorted.bam -o $filenm_root.indelrealigned.RG.sorted.HC.vcf";
print ("Call SNPs with HaplotypeCaller on file $filenm_root.indelrealigned.RG.sorted.bam...\n");
system("$haplotype_caller");

#VariantAnnotator (add annotation: MappingQualityZero)
my $variant_annotation = "java -Xmx4g -jar GenomeAnalysisTK.jar -T VariantAnnotator -R $reference -V $filenm_root.indelrealigned.RG.sorted.HC.vcf -I $filenm_root.indelrealigned.RG.sorted.bam -A MappingQualityZero -o $filenm_root.indelrealigned.RG.sorted.HC.ann.vcf";
print ("Adding MappingQualityZero annotation in file $filenm_root.indelrealigned.RG.sorted.HC.vcf...\n");
system("$variant_annotation");

#VariantFiltration (Filter SNPs before applying BaseRecalibrator)
my $variant_filtration  = "java -Xmx4g -jar GenomeAnalysisTK.jar -T VariantFiltration -R $reference -V $filenm_root.indelrealigned.RG.sorted.HC.ann.vcf --clusterSize 3 --clusterWindowSize 10 --filterExpression 'MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)' --filterExpression 'QD < 1.5' --filterExpression 'DP < 15' --maskExtension 0 --filterName HARD_TO_VALIDATE --filterName QD_FILTER --filterName DP_FILTER --maskName Mask -o $filenm_root.indelrealigned.RG.sorted.HC.ann.VF.vcf";
print ("Variant filtration --filterExpression 'MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)' --filterExpression 'QD < 1.5' --filterExpression 'DP < 15' on $filenm_root.indelrealigned.RG.sorted.HC.ann.vcf...\n");
system("$variant_filtration");

#SelectVariant to filter out filtered variants
my $select_variant = "java -Xmx4g -jar GenomeAnalysisTK.jar -T SelectVariants -R $reference -V $filenm_root.indelrealigned.RG.sorted.HC.ann.VF.vcf -select 'vc.isNotFiltered()' -selectType SNP -o $filenm_root.indelrealigned.RG.sorted.HC.ann.VF.filtered.vcf";
print ("Filter out filtered SNPs with SelectVariants -select 'vc.isNotFiltered()' -selectType SNP on file $filenm_root.indelrealigned.RG.sorted.HC.ann.VF.vcf...\n");
system("$select_variant");

#BaseRecalibrator
my $base_recalibrator = "java -Xmx4g -jar GenomeAnalysisTK.jar -T BaseRecalibrator -R $reference -I $filenm_root.indelrealigned.RG.sorted.bam --knownSites $filenm_root.indelrealigned.RG.sorted.HC.ann.VF.filtered.vcf -o $filenm_root.recal_data.table";
print ("Create a recalibration table file with BaseRecalibrator with knownSites in file $filenm_root.indelrealigned.RG.sorted.HC.ann.VF.filtered.vcf with BAM file $filenm_root.indelrealigned.RG.sorted.bam...\n");
system("$base_recalibrator");

#Apply recalibration with Printreads
my $apply_recalibration = "java -Xmx4g -jar GenomeAnalysisTK.jar -T PrintReads -R $reference -I $filenm_root.indelrealigned.RG.sorted.bam -BQSR $filenm_root.recal_data.table -o $filenm_root.indelrealigned.RG.sorted.recalibrated.bam";
print ("Apply recalibration with PrintReads on file $filenm_root.indelrealigned.RG.sorted.bam...\n");
system("$apply_recalibration");

#Obtain a gVCF with HaplotypeCaller
my $haplotype_caller_gvcf = "java -Xmx4g -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R $reference -I $filenm_root.indelrealigned.RG.sorted.recalibrated.bam --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -o $filenm_root.indelrealigned.RG.sorted.recalibrated.HC.gvcf";
print ("Call SNPs with HaplotypeCaller on file $filenm_root.indelrealigned.RG.sorted.bam and output a gvcf (--emitRefConfidence GVCF)...\n");
system("$haplotype_caller_gvcf");

print "Finished\n";

exit 1;
