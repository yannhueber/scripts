#! /usr/bin/perl -w
#$ -q bioinfo.q
#$ -cwd
#$ -V
#$ -t 1-12
#$ -S /usr/bin/perl
##$ -pe parallel_something 12



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
use Statistics::Descriptive;

############### MAIN ##############

#Variables
my ($extension,$lod_limit,$pdiseq_limit);

#Getoptions
GetOptions("x|extension=s"=>\$extension, "min_lod=i"=>\$lod_limit, "max_pdiseq=f"=>\$pdiseq_limit);

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

#open file
open(FILE,"<$file_to_treat") or die("Cannot open $file_to_treat : $!");#Can be a LD file from haploview or from tassel (the two files are slightly different)
open(OUT1,">$filename.reduced") or die ("Cannot open $filename.reduced : $!");
open(SPLIT1,">$filename.reduced.0kb-10kb") or die ("Cannot open $filename.reduced.0kb-10kb : $!");
open(SPLIT2,">$filename.reduced.10kb-50kb") or die ("Cannot open $filename.reduced.10kb-50kb : $!");
open(SPLIT3,">$filename.reduced.50kb-100kb") or die ("Cannot open $filename.reduced.50kb-100kb : $!");
open(SPLIT4,">$filename.reduced.100kb-500kb") or die ("Cannot open $filename.reduced.100kb-500kb : $!");
open(SPLIT5,">$filename.reduced.500kb-1mb") or die ("Cannot open $filename.reduced.500kb-1mb : $!");
open(SPLIT6,">$filename.reduced.1mb-5mb") or die ("Cannot open $filename.reduced.1mb-5mb : $!");
open(SPLIT7,">$filename.reduced.5mb-10mb") or die ("Cannot open $filename.reduced.5mb-10mb : $!");
open(SPLIT8,">$filename.reduced.10mb-100mb") or die ("Cannot open $filename.reduced.10mb-100mb : $!");
open(SPLIT9,">$filename.reduced.100mb-infinity") or die ("Cannot open $filename.reduced.10mb-100mb : $!");
open(SPLIT10,">$filename.reduced.unlinked") or die ("Cannot open $filename.reduced.unlinked : $!");
open(OUT2,">$filename.summary") or die ("Cannot open $filename.summary : $!");


print "HERE $pdiseq_limit\n";


my %frames = ("1" => "0kb-10kb", "2" => "10kb-50kb", "3" => "50kb-100kb", "4" => "100kb-500kb", "5" => "500kb-1mb", "6" => "1mb-5mb", "7" => "5mb-10mb", "8" => "10mb-100mb", "9" => "100mb-infinity", "10" => "unlinked");
my %counter;

while (defined (my $line = <FILE>)) {
	if ($line =~ m/^Locus1/ or $line =~ m/^L1/) {
		print OUT1 $line; print SPLIT1 $line; print SPLIT2 $line; print SPLIT3 $line; print SPLIT4 $line; print SPLIT5 $line; print SPLIT6 $line; print SPLIT7 $line; print SPLIT8 $line; print SPLIT9 $line; print SPLIT10 $line;
		next;
	}
	chomp($line);
	my @infos_line = split(m/\t/,$line);
	my $confidence_on_r2 = 0; #LOD value for haploview and pDiseq for TASSEL
	my $dist_between_markers;
	my $r2;
	my $dprime;
	if ($file_to_treat =~ m/Tassel/) {
		my $pdiseq = $infos_line[15];
		if ($pdiseq < $pdiseq_limit) {
			$confidence_on_r2 = 1;
		}
		$dist_between_markers = $infos_line[12];
		$r2 = $infos_line[13];
		$dprime = $infos_line[14];
	} elsif ($file_to_treat =~ m/Haploview/) {
		my $lod = $infos_line[3];
		if ($lod >= $lod_limit) {
			$confidence_on_r2 = 1;
		}
		$dist_between_markers = $infos_line[7];
		$r2 = $infos_line[4];
		$dprime = $infos_line[2];
	} else {
		die("Is it a Tassel or an Haploview file !?\n");
	}

	next if (! $confidence_on_r2);
	print OUT1 $line . "\n";

	my $case;
	if ($dist_between_markers eq "N/A") {
		print SPLIT10 $line . "\n";
		$case = 10;
	}
	elsif ($dist_between_markers <= 10000) {
		print SPLIT1 $line . "\n";
		$case = 1;
	} elsif ($dist_between_markers <= 50000) {
		print SPLIT2 $line . "\n";
		$case = 2;
	} elsif ($dist_between_markers <= 100000) {
		print SPLIT3 $line . "\n";
		$case = 3;
	} elsif ($dist_between_markers <= 500000) {
		print SPLIT4 $line . "\n";
		$case = 4;
	} elsif ($dist_between_markers <= 1000000) {
		print SPLIT5 $line . "\n";
		$case = 5;
	} elsif ($dist_between_markers <= 5000000) {
		print SPLIT6 $line . "\n";
		$case = 6;
	} elsif ($dist_between_markers <= 10000000) {
		print SPLIT7 $line . "\n";
		$case = 7;
	} elsif ($dist_between_markers <= 100000000) {
		print SPLIT8 $line . "\n";
		$case = 8;
	} else {
		print SPLIT9 $line . "\n";
		$case = 9;
	}

	#increment r2, dprime and number of comparison per file
	if (! $counter{$case}) {
		$counter{$case}{"count"} = 1;
	} else {
		$counter{$case}{"count"}++;
	}
	$counter{$case}{"sum_r2"} .= $r2 . ",";
	$counter{$case}{"sum_dprime"} .= $dprime . ",";
}



#print summary

#header
print OUT2 "FRAME\tCOMPARISON_NUMBER\tR2_AVERAGE\tR2_VAR\tDPRIME_AVERAGE\tDPRIME_VAR\n";
#body
foreach my $frame (sort {$a <=> $b} keys(%counter)) {
	print $frame . "\n";
	chop($counter{$frame}{"sum_r2"});
	chop($counter{$frame}{"sum_dprime"});
	my @split_r2 = split(m/,/,$counter{$frame}{"sum_r2"});
	my @split_dprime = split(m/,/,$counter{$frame}{"sum_dprime"});
	my $stat1 = Statistics::Descriptive::Full->new();
	my $stat2 = Statistics::Descriptive::Full->new();
	$stat1->add_data(\@split_r2);
	$stat2->add_data(\@split_dprime);

	print OUT2 $frames{$frame} . "\t" . $counter{$frame}{"count"} . "\t" . $stat1->mean() . "\t" . $stat1->variance() . "\t" . $stat2->mean() . "\t" . $stat2->variance() . "\n";
}


close(OUT1);
close(OUT2);
close(SPLIT1);
close(SPLIT2);
close(SPLIT3);
close(SPLIT4);
close(SPLIT5);
close(SPLIT6);
close(SPLIT7);
close(SPLIT8);
close(SPLIT9);
close(SPLIT10);

exit 1;
