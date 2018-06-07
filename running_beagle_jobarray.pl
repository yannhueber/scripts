#! /usr/bin/perl -w
#$ -q bioinfo.q
#$ -cwd
#$ -V
#$ -t 1-12
#$ -S /usr/bin/perl
#$ -pe parallel_fill 1



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

#chr_number
#my $chr_number;
#if ($filename =~ m/(\d+)/) {
#	$chr_number = $1;
#} else {
#	$chr_number = "all";
#}


#Programs 

#Beagle
my $beagle = "/usr/local/jdk1.7.0_40/bin/java -Xmx2g -jar /home/hueber/bin/b4.r1230.jar gt=$file_to_treat phase-its=40 impute-its=10 out=$filename.out";
print "BEAGLE working with file : $file_to_treat...\n";
system("$beagle");


print "Finished\n";



exit;
