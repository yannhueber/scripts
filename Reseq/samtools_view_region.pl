#! /usr/bin/perl -w
#$ -q normal.q
#$ -cwd
#$ -V
#$ -t 1-113
#$ -S /usr/bin/perl
#$ -pe parallel_fill 1


use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use File::Path;



my @regions = ("chr06:12617481-12721481","chr03:28457932-28462132","chr07:28197511-28204544");




############### MAIN ##############

#Usage
if (scalar (@ARGV) == 0)
{
        print "\n\n\nUsage: perl $0 -x extension\n\n\n";
        exit 1;
}


#Variables
my $extension;


#Getoptions
GetOptions("x|ext=s" => \$extension) or die ("Error in command line arguments\n");


#List of all file in directory with specific extension
opendir(CWD,".") or die("Cannot open current directory\n");
my @files = grep { m/\.$extension$/ } readdir(CWD);
closedir(CWD);


#relative file path 
my $sge_id = `echo \$SGE_TASK_ID` or die ("Cannot catch SGE_TASK_ID\n");
chomp($sge_id);
$sge_id--;
my $files_to_treat = $files[$sge_id];



#Programs
foreach my $region (@regions) {
	
	my $output_suffixe = $region;
	$output_suffixe =~ s/:/-/;

	#samtools view to select reads in a region
	my $samtools_view = "samtools view -h $files_to_treat $region > $files_to_treat.$output_suffixe.sam";
	print "samtools view -h $files_to_treat $region > $files_to_treat.$output_suffixe.sam\n";
	system("$samtools_view");

	my $samtobam = "samtools view -b $files_to_treat.$output_suffixe.sam -o $files_to_treat.$output_suffixe.bam2";
        print "samtools view -b $files_to_treat.$output_suffixe.sam -o $files_to_treat.$output_suffixe.bam2\n";
        system("$samtobam");

	my $flagstat = "samtools flagstat $files_to_treat.$output_suffixe.bam2 > $files_to_treat.$output_suffixe.flagstat";
        print "samtools flagstat $files_to_treat.$output_suffixe.bam2 > $files_to_treat.$output_suffixe.flagstat\n";
    	system("$flagstat");
}



print "Finished\n";

exit 1;













