#! /usr/bin/perl -w
#$ -q bioinfo.q
#$ -cwd
#$ -V
#$ -t 1-13
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
my $chr_number;
if ($filename =~ m/(\d+)/) {
	$chr_number = $1;
} else {
	$chr_number = "all";
}


#Programs 

#fastphase
my $fastphase = "fastPHASE -q0.8 -T10 -C40 -H200 -i -KU50 -KL10 -KI10 -U20 -w -Pm -o\"results_chr$chr_number\_fastphase\" $file_to_treat";
print "fastPHASE working with file : $file_to_treat...\n";
system("$fastphase");


print "Finished\n";



exit;
