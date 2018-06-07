#! /usr/bin/perl -w
#$ -q bioinfo.q
#$ -cwd
#$ -V
#$ -t 1-4
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

#Variables
my ($extension, $barcode_filename);

#Getoptions
GetOptions("x|extension=s"=>\$extension,
	   "b|barcodefile=s"=>\$barcode_filename);

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
my $split = "zcat $file_to_treat | fastx_barcode_splitter.pl --bcfile $barcode_filename --bol --exact --prefix results/ --suffix \".$filename.fastq\"";
print "fastx_barcode_splitter.pl working with file : $file_to_treat...\n";
system("$split");


print "Finished\n";

exit;
