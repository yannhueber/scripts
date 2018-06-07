#! /usr/bin/perl -w
#$ -q bigmem.q
#$ -cwd
#$ -V
#$ -t 1-48
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
my $star_path = "/usr/local/bioinfo/STAR/2.5.0b/bin";

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
        print "\n\nUsage: perl $0 -g genome_dir_path -r reference_fasta -x extension_file_to_treat -cu cultivar\n\n\n";
        exit 1;
}


#Variables
my ($extension, $reference, $genome_indexes_dir, $cultivar);

#Getoptions
GetOptions("x|extension=s" => \$extension, "r|reference=s" => \$reference, "g|genomeindexesdir=s" => \$genome_indexes_dir, "c|cultivar=s" => \$cultivar) or die ("Error in command line arguments\n");



#List of all file in directory with specific extension
opendir(CWD,".") or die("Cannot open current directory\n");  
my @files = grep { m/\.$extension$/ } readdir(CWD);
closedir(CWD);


#relative file path 
my $sge_id = `echo \$SGE_TASK_ID` or die ("Cannot catch SGE_TASK_ID\n");
chomp($sge_id);
$sge_id--; #perl is in basis 0
my $input_se_fastq = $files[$sge_id]; #se means single end


#Filename without extension
my $filenm_root = fileparse($input_se_fastq, qr/\.$extension$/);


#Programs 

#Discard filtered reads 
#my $discard_filtered_reads = "zcat $input_se_fastq | grep -A 3 '^@.* [^:]*:N:[^:]*:' | grep -v \"^--\$\" | gzip > $filenm_root.RPF.fastq.gz"; #RPF = reads primary filtered
#print ("Filter out bad reads on file $input_se_fastq\n");
#system("$discard_filtered_reads");

#Cutadapt
my $rm_adapter = "$cutadapt_path/cutadapt -b AGATCGGAAGAGC -O 7 -m 30 -q 20,20 $input_se_fastq | gzip > $filenm_root.cutadapt.fastq.gz"; #Remove adaptator, Trim on quality, discard read with length inf to 30 
print ("cutadapt on $filenm_root.fastq.gz\n");
system("$rm_adapter");

#fastqc
my $fastqc = "$fastqc_path/fastqc -j $java_path/java $filenm_root.cutadapt.fastq.gz";
print ("fastqc on $filenm_root.cutadapt.fastq.gz\n");
system("$fastqc");

#STAR mapping
my $star_mapping = "$star_path/STAR --genomeDir $genome_indexes_dir --readFilesIn $filenm_root.cutadapt.fastq.gz --readFilesCommand zcat --outSAMunmapped Within --outFileNamePrefix $filenm_root. --outSAMmapqUnique 255 --twopassMode Basic";
print ("Mapping with STAR (2-pass mode) on file $filenm_root.cutadapt.fastq.gz...\n");
system("$star_mapping");

#Add read group, and sort bam
my $condition = $filenm_root;
my $add_rg = "$java_path/java -Xmx5g -jar $picard_path/picard.jar AddOrReplaceReadGroups I=$filenm_root.Aligned.out.sam O=$filenm_root.RG.sorted.bam SO=coordinate RGLB=$cultivar RGPL=illumina RGPU=run RGSM=$condition RGID=$condition";
print ("Adding Read Group and sorting on file $filenm_root.Aligned.out.sam...\n");
system("$add_rg");

#Markduplicates and create index
my $markduplicates = "$java_path/java -Xmx5g -jar $picard_path/picard.jar MarkDuplicates I=$filenm_root.RG.sorted.bam O=$filenm_root.dedupped.RG.sorted.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=$filenm_root.output.mectrics";
print ("Marking duplicates and creating index on file $filenm_root.RG.sorted.bam...\n");
system("$markduplicates");


#SplitNCigarReads ######## add -U ALLOW_N_CIGAR_READS in case...
my $split_and_trim = "$java_path/java -Xmx5g -jar $gatk_path/GenomeAnalysisTK.jar -T SplitNCigarReads -R $reference -I $filenm_root.dedupped.RG.sorted.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS -o $filenm_root.split.dedupped.RG.sorted.bam";
print ("Split N CigarReads on file $filenm_root.dedupped.RG.sorted.bam...\n");
system("$split_and_trim");

#RealignerTargetCreator ######Add --fix_misencoded_quality_scores option if needed
my $realigner_target = "$java_path/java -Xmx5g -jar $gatk_path/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $reference -I $filenm_root.split.dedupped.RG.sorted.bam -o $filenm_root.forIndelRealigner.intervals";
print ("Create RealignTargetCreator intervals with $filenm_root.split.dedupped.RG.sorted.bam...\n");
system("$realigner_target");

#IndelRealigner ######Add --fix_misencoded_quality_scores option if needed
my $indelrealigner = "$java_path/java -Xmx5g -jar $gatk_path/GenomeAnalysisTK.jar -T IndelRealigner -R $reference -I $filenm_root.split.dedupped.RG.sorted.bam -targetIntervals $filenm_root.forIndelRealigner.intervals -o $filenm_root.indelrealigned.split.dedupped.RG.sorted.bam";
print ("IndelRealigning on file $filenm_root.split.dedupped.RG.sorted.bam...\n");
system("$indelrealigner");

#HaplotypeCaller (for BaseRecalibrator)
my $haplotype_caller = "$java_path/java -Xmx5g -jar $gatk_path/GenomeAnalysisTK.jar -T HaplotypeCaller -R $reference -I $filenm_root.indelrealigned.split.dedupped.RG.sorted.bam -ploidy 3 -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 -o $filenm_root.indelrealigned.split.dedupped.RG.sorted.HC.vcf";
print ("Call SNPs with HaplotypeCaller on file $filenm_root.indelrealigned.split.dedupped.RG.sorted.bam...\n");
system("$haplotype_caller");

#VariantAnnotator (add annotation: MappingQualityZero)
my $variant_annotation = "$java_path/java -Xmx5g -jar $gatk_path/GenomeAnalysisTK.jar -T VariantAnnotator -R $reference -V $filenm_root.indelrealigned.split.dedupped.RG.sorted.HC.vcf -I $filenm_root.indelrealigned.split.dedupped.RG.sorted.bam -A MappingQualityZero -o $filenm_root.indelrealigned.split.dedupped.RG.sorted.HC.ann.vcf";
print ("Adding MappingQualityZero annotation in file $filenm_root.indelrealigned.split.dedupped.RG.sorted.HC.vcf...\n");
system("$variant_annotation");

#VariantFiltration (Filter SNPs before applying BaseRecalibrator)
my $variant_filtration  = "$java_path/java -Xmx5g -jar $gatk_path/GenomeAnalysisTK.jar -T VariantFiltration -R $reference -V $filenm_root.indelrealigned.split.dedupped.RG.sorted.HC.ann.vcf --clusterSize 3 --clusterWindowSize 10 --filterExpression 'MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)' --filterExpression 'QD < 1.5' --filterExpression 'DP < 15' --maskExtension 0 --filterName HARD_TO_VALIDATE --filterName QD_FILTER --filterName DP_FILTER --maskName Mask -o $filenm_root.indelrealigned.split.dedupped.RG.sorted.HC.ann.VF.vcf";
print ("Variant filtration --filterExpression 'MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)' --filterExpression 'QD < 1.5' --filterExpression 'DP < 15' on $filenm_root.indelrealigned.split.dedupped.RG.sorted.HC.ann.vcf...\n");
system("$variant_filtration");

#SelectVariant to filter out filtered variants
my $select_variant = "$java_path/java -Xmx5g -jar $gatk_path/GenomeAnalysisTK.jar -T SelectVariants -R $reference -V $filenm_root.indelrealigned.split.dedupped.RG.sorted.HC.ann.VF.vcf -select 'vc.isNotFiltered()' -selectType SNP -o $filenm_root.indelrealigned.split.dedupped.RG.sorted.HC.ann.VF.filtered.vcf";
print ("Filter out filtered SNPs with SelectVariants -select 'vc.isNotFiltered()' -selectType SNP on file $filenm_root.indelrealigned.split.dedupped.RG.sorted.HC.ann.VF.vcf...\n");
system("$select_variant");

#BaseRecalibrator
my $base_recalibrator = "$java_path/java -Xmx5g -jar $gatk_path/GenomeAnalysisTK.jar -T BaseRecalibrator -R $reference -I $filenm_root.indelrealigned.split.dedupped.RG.sorted.bam --knownSites $filenm_root.indelrealigned.split.dedupped.RG.sorted.HC.ann.VF.filtered.vcf -o $filenm_root.recal_data.table";
print ("Create a recalibration table file with BaseRecalibrator with knownSites in file $filenm_root.indelrealigned.split.dedupped.RG.sorted.HC.ann.VF.filtered.vcf with BAM file $filenm_root.indelrealigned.split.dedupped.RG.sorted.bam...\n");
system("$base_recalibrator");

#Apply recalibration with Printreads
my $apply_recalibration = "$java_path/java -Xmx5g -jar $gatk_path/GenomeAnalysisTK.jar -T PrintReads -R $reference -I $filenm_root.indelrealigned.split.dedupped.RG.sorted.bam -BQSR $filenm_root.recal_data.table -o $filenm_root.indelrealigned.split.dedupped.RG.sorted.recalibrated.bam";
print ("Apply recalibration with PrintReads on file $filenm_root.indelrealigned.split.dedupped.RG.sorted.bam...\n");
system("$apply_recalibration");

print "Finished\n";



exit 1;
