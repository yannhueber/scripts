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

############### MAIN ##############



#Variables
my ($extension);

if (scalar (@ARGV)==0) {
	print ("\n\nperl $0 -x extension\n\n");
	exit 1;
}


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

#open file
open(HAPMAP,"<$file_to_treat") or die("Cannot open $file_to_treat : $!");
open(PED,">$filename.ped") or die ("Cannot open $filename.ped : $!");
open(INFO,">$filename.info") or die ("Cannot open $filename.info");


my %genotypes;
my @individuals =();
my @positions = ();

while (defined (my $line = <HAPMAP>)) {
	next if ($line =~ m/^#@/);
	chomp($line);
	chop($line);
	my @infos_line = split(m/\t/,$line);
	if ($line =~ m/^rs/) {
		@individuals = @infos_line[11..$#infos_line];
		next;
	}

	my $rs = $infos_line[0];
	push(@positions,$rs);

	my @indiv_genotypes = @infos_line[11..$#infos_line];

	foreach (my $i = 0; $i < @individuals; $i++) {
		my $indv = $individuals[$i];
		my $indv_genotype = $indiv_genotypes[$i];
		$genotypes{$indv}{$rs} = $indv_genotype;
	}
}


#print ped file
#header
#no header
#body
my $incr =1;
foreach my $indiv (sort {$a cmp $b } keys(%genotypes)) {
	my $ref = $genotypes{$indiv};
	my %rs = %$ref;
	my $line_to_print = "$incr\t$indiv\t0\t0\t2\taffected_status\t";	

	foreach my $rs (@positions) {
		my $indv_genotype = $genotypes{$indiv}{$rs};
		$indv_genotype =~ s/NN/00/;
		my $first_allele = substr($indv_genotype,0,1);
		my $second_allele =substr($indv_genotype,1,1);
		$line_to_print .= "$first_allele $second_allele\t";
	}
	chop($line_to_print);
	$line_to_print .= "\n";
	print PED $line_to_print;
	$incr++;
}



foreach my $rs (@positions) {
	my @rs_split = split(m/_/,$rs);
	my $position = $rs_split[1];
	print INFO "$rs\t$position\n";
}

close(HAPMAP);
close(PED);
close(INFO);

exit 1;
