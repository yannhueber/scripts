#! /usr/bin/perl -w
#$ -q normal.q
#$ -cwd
#$ -V
#$ -t 1-2
#$ -S /usr/bin/perl
#$ -pe parallel_fill 1


use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use File::Path;


#Usage
if (scalar (@ARGV) == 0)
{
        print "\n\n\nUsage: perl $0 -x extension -g gff_file\n\n\n";
        exit 1;
}


#Variables
my $extension;
my $gff_file;

#Getoptions
GetOptions("x|ext=s" => \$extension, "g|gff=s" => \$gff_file) or die ("Error in command line arguments\n");


#List of all file in directory with specific extension
opendir(CWD,".") or die("Cannot open current directory\n");
my @files = grep { m/\.$extension$/ } readdir(CWD);
closedir(CWD);


#relative file path 
my $sge_id = `echo \$SGE_TASK_ID` or die ("Cannot catch SGE_TASK_ID\n");
chomp($sge_id);
my $files_to_treat = $files[$sge_id];



#Programs
my $htsseqcount = "python /homedir/hueber/.local/bin/htseq-count --order=pos --idattr=Parent --format=bam $files_to_treat $gff_file";
print "htseq-count --order=pos --idattr=Parent --format=bam $files_to_treat $gff_file\n";
system("$htsseqcount");

print "Finished\n";

exit 1;
