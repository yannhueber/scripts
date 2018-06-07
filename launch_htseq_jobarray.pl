#! /usr/bin/perl -w
#$ -q normal.q
#$ -cwd
#$ -V
#$ -t 1-60
#$ -S /usr/bin/perl
#$ -pe parallel_fill 1

my $python_path = "/usr/local/bioinfo/python/2.7.9_build2/bin";
my $htseq_path = "/homedir/hueber/.local/bin";

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


if (@ARGV ==0)
{
        print "\n\nUsage: perl $0 -g genes.gff3 -x extension_file_to_treat\n\n\n";
        exit 1;
}


#Variables
my ($extension, $gff);

#Getoptions
GetOptions("x|extension=s" => \$extension, "g|gfffile=s" => \$gff) or die ("Error in command line arguments\n");


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


#Programs 

#HTSeq
my $htseqcount = "$python_path/python $htseq_path/htseq-count --order=pos --idattr=Parent --format=bam $file $gff > $filenm_root.counts";
print ("$python_path/python $htseq_path/htseq-count --order=pos --idattr=Parent --format=bam $file $gff > $filenm_root.counts\n");
system("$htseqcount");

print "Finished\n";

exit 1;
