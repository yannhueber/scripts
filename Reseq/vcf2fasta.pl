#! /usr/bin/perl -w
#$ -q normal.q
#$ -cwd
#$ -V
#$ -t 1-113
#$ -S /usr/bin/perl
##$ -pe something 12



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


if (scalar(@ARGV) == 0) {

	print "\n\nperl $0 -x extension\n\n";
	exit 1;
}





#Variables
my ($extension);

#Getoptions
GetOptions("x|extension=s"=>\$extension);


#List of all file in directory with specific extension
opendir(CWD,".") or die("Cannot open current directory\n");  
my @files = grep { m/\.$extension$/ } readdir(CWD);
closedir(CWD);

#relative file path 
my $sge_id = `echo \$SGE_TASK_ID` or die ("Cannot catch SGE_TASK_ID\n");
chomp($sge_id);
$sge_id--; #perl is in basis 0
my $file_to_treat = $files[$sge_id];

#Filename without extension
my $filename = fileparse($file_to_treat, qr/\.$extension$/);


#Programs 

#barcode splitter
my $vcf2fasta = "/usr/local/java/jre8/bin/java -Xmx2G -jar /usr/local/bioinfo/GenomeAnalysisTK/3.7-0/GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker -R ~/musa_ref_v2/Musa_genomes_V2.fasta -o $filename.fasta -L interval_file2.interval_list -V $file_to_treat";
print "java -Xmx2G -jar GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker -R ~/musa_ref_v2/Musa_genomes_V2.fasta -o $filename.fasta -L interval_file2.interval_list -V $file_to_treat\n";
system("$vcf2fasta");


print "Finished\n";
exit 1;
