#! /usr/bin/perl -w
#$ -q bioinfo.q
#$ -cwd
#$ -V
#$ -t 1-2
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
open(DAR,">$filename.var") or die ("Cannot open $filename.var : $!");


my %genotypes;
my @individuals =();
my %transfo_ATGC = ("A" => "1", "T" => "2", "C" => "3", "G" => "4", "D" => 5); 
my $nb_markers=0;
my @positions;


while (defined (my $line = <HAPMAP>)) {
	chop($line);
	my @infos_line = split(m/\t/,$line);
	my @indiv_infos = @infos_line[11..$#infos_line];

	if ($line =~ m/^rs/) {
		@individuals = @indiv_infos;
		next;
	}

	my $rs = $infos_line[0];
	push(@positions,$rs);
	my @indiv_genotypes = @indiv_infos;
	$nb_markers++;
	for (my $i = 0; $i < @individuals; $i++) {
		my $indv = $individuals[$i];
		my $indv_genotype = $indiv_genotypes[$i];
		$genotypes{$indv}{$rs} = $indv_genotype;
	}
}


#print var file (for darwin)
#header
print DAR "\@DARwin 5.0 - ALLELIC - 2\n";

my $nb_indivs = scalar(@individuals);
my $nb_markers_x2 = $nb_markers*2;
print DAR "$nb_indivs\t$nb_markers_x2\nUnit";

for (my $j = 1; $j <= $nb_markers_x2; $j++) {
	print DAR "\t$j";
}

print DAR "\n";


#body
my $incr_indv =1;

foreach my $indiv (@individuals) {
	print DAR $incr_indv;
	foreach my $rs (@positions) {
		my $indv_genotype = $genotypes{$indiv}{$rs};
		my $indv_genotype_number = $transfo_ATGC{$indv_genotype};
		print DAR "\t$indv_genotype_number\t$indv_genotype_number";
	}
	print DAR "\n";
	$incr_indv++;
}

close(HAPMAP);
close(DAR);
