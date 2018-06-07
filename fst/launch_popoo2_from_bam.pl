#! /usr/bin/perl -wi
#$ -q normal.q
#$ -cwd
#$ -V
#$ -t 1-55
#$ -S /usr/bin/perl
#$ -pe parallel_fill 1

#charge modules
#use "bioinfo/bwa/0.7.12";
#use "system/java/jre8";
#use "bioinfo/FastQC/0.11.3";
#use "compiler/gcc/4.9.2";
#use "bioinfo/bwa/0.7.12";


my $java_path = "/usr/local/java/jre8/bin";
my $samtools_path = "/usr/local/bioinfo/samtools/1.3/bin";

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

my $extension;
my $pop_size_file;



if (@ARGV == 0)
{
        print "\n\nUsage: perl $0 -s pop_size_tab_file -x extension_file_to_treat\n\n\n";
        exit 1;
}



#Getoptions
GetOptions("x|extension=s" => \$extension, "s|popsize=s" => \$pop_size_file)  or die ("Error in command line arguments\n");


my %pop_size;
open(F,"<$pop_size_file") or die("Cannot open $pop_size_file : $!");

while (my $line = <F>) {

	chomp($line);
	my @infos = split(m/\t/,$line);
	my $pop = $infos[0];
	my $size = $infos[1];
	$pop_size{$pop} = $size;
}




#List of all file in directory with specific extension
opendir(CWD,".") or die("Cannot open current directory\n");  
my @files = grep { m/\.$extension$/ } readdir(CWD);
closedir(CWD);


my @combinaison;
my $i =0;
my $j = scalar(@files);

while ($i < $j) {

        my $k = $i+1;

        while ($k < $j) {

                my $str_f = $files[$i] . "-" . $files[$k];
                push(@combinaison, $str_f);
                $k++;
        }
        $i++;
}


#relative file path 
my $sge_id = `echo \$SGE_TASK_ID` or die ("Cannot catch SGE_TASK_ID\n");
chomp($sge_id);
$sge_id--; #perl is in basis 0
my $files_to_treat = $combinaison[$sge_id];


my $bam1;
my $bam2;
if ($files_to_treat =~ m/(.*)-(.*)/ ) {
	$bam1 = $1;
	$bam2 = $2;
} else {
	print "problem with $files_to_treat\n";
	exit 1;
}





#filename
my @split_filenm1 = split(m/\./,$bam1);
my $pop1 = $split_filenm1[0];
my @split_filenm2 = split(m/\./,$bam2);
my $pop2 = $split_filenm2[0];

my $pop_size1 = $pop_size{$pop1};
my $pop_size2 = $pop_size{$pop2};


#Programs 

#mpileup
my $samtools_mpileup = "$samtools_path/samtools mpileup -B $bam1 $bam2 -o $pop1-$pop2.mpileup";
print ("$samtools_path/samtools mpileup -B $bam1 $bam2 -o $pop1-$pop2.mpileup\n");
system("$samtools_mpileup");


#mpileup2sync
my $mpileup_sync = "$java_path/java -Xmx3g -jar popoolation2_1201/mpileup2sync.jar --input $pop1-$pop2.mpileup --output $pop1-$pop2.sync --fastq-type sanger --min-qual 20";
print ("$java_path/java -Xmx3g -jar popoolation2_1201/mpileup2sync.jar --input $pop1-$pop2.mpileup --output $pop1-$pop2.sync --fastq-type sanger --min-qual 20\n");
system("$mpileup_sync");


#fst sliding windows
my $fst = "perl popoolation2_1201/fst-sliding.pl --input $pop1-$pop2.sync --output $pop1-$pop2.fst --min-count 3 --min-coverage 3 --max-coverage 200 --min-covered-fraction 0 --window-size 100000 --step-size 50000 --pool-size $pop_size1:$pop_size2";
print ("perl popoolation2_1201/fst-sliding.pl --input $pop1-$pop2.sync --output $pop1-$pop2.fst --min-count 3 --min-coverage 3 --max-coverage 200 --min-covered-fraction 0 --window-size 100000 --step-size 50000 --pool-size $pop_size1:$pop_size2\n");
system("$fst");


print "Finished\n";

exit 1;






