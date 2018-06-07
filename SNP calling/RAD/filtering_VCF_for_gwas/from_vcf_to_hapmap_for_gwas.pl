#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use File::Path;
use POSIX qw(ceil);


#Global variables
my ($reference, $file_extension, $maf, $mininbreedingcoeff, $indv_count, $ploidy, $missingness, $readdepthlimitmin, $readdepthlimitmax, $readdepthmin_allele, $maxindivwithmissingdata, $allowed_missing_percent);
my $phaseits = "5";
my $imputeits = "5";

my $gatk_path = "/usr/local/bioinfo/GenomeAnalysisTK/3.4-46";
my $java_path = "/usr/local/java/jre8/bin";



#print usage if no arguments are specified 
if (scalar @ARGV==0) {
	print "\n\nUsage: perl $0 -r REF_PATH -n number_of_indv -f minor_allele_frequency -c min_inbreeding_coeff -a allowed_missing_percent -i min_missingness -s readdepthlimitmin -t readdepthlimitmax -u readdepthmin_allele -v maxindivwithmissingdata -p ploidy -x file_extension\n\n";
	exit;
}


#GetOptions
GetOptions("r|reference=s" => \$reference,
	   "x|extension=s" => \$file_extension,
	   "f|maf=s" => \$maf,
	   "a|allowedmissingpercent=i" => \$allowed_missing_percent,
	   "n|nb_indv=i" => \$indv_count,
	   "h|phaseits=i" => \$phaseits,
	   "m|imputeits=i" => \$imputeits,
	   "p|ploidy=i" => \$ploidy,
	   "i|missingness=f" => \$missingness,
	   "s|readdepthlimit_min=i" => \$readdepthlimitmin,
	   "t|readdepthlimit_max=i" => \$readdepthlimitmax,
	   "u|readdepthminallele=i" => \$readdepthmin_allele,
	   "v|maxindivwithmissingdata=i" => \$maxindivwithmissingdata,
	   "c|min_ic=f" => \$mininbreedingcoeff);


#List of all vcf file in the directory
opendir(CWD,".") or die("Cannot access current directory\n");
my @files = grep { m/\.$file_extension$/ } readdir(CWD);
closedir(CWD);
foreach my $filename (sort(@files)) {

	my $filenm_no_ext = fileparse($filename, qr/\.$file_extension/);

	#VariantFiltration (Filter SNPs before applying BaseRecalibrator)
	my $wdyd1 = "VF";
	my $variant_filtration  = "$java_path/java -Xmx4g -jar $gatk_path/GenomeAnalysisTK.jar -T VariantFiltration -R $reference -V $filename --clusterSize 3 --clusterWindowSize 10 --filterExpression 'QD < 1.5' --maskExtension 0 --filterName QD_FILTER -maskName Mask -o $filenm_no_ext.$wdyd1.vcf";
	print ("Variant filtration --clusterSize 3 --clusterWindowSize 10 --filterExpression 'QD < 1.5' --maskExtension 0 --filterName QD_FILTER -maskName Mask on file $filename...\n");
	system("$variant_filtration");

	#SelectVariant to filter out filtered variants
	my $wdyd2 = "$wdyd1.filtered";
	my $select_variant = "$java_path/java -Xmx4g -jar $gatk_path/GenomeAnalysisTK.jar -T SelectVariants -R $reference -V $filenm_no_ext.$wdyd1.vcf -select 'vc.isNotFiltered()' -selectType SNP -o $filenm_no_ext.$wdyd2.vcf";
	print ("Filter out filtered SNPs with SelectVariants -select 'vc.isNotFiltered()' and -selectType SNP on file $filenm_no_ext.$wdyd1.vcf...\n");
	system("$select_variant");

	#discard individuals based on missingness
	$wdyd1 = $wdyd2;
	$wdyd2 = "$wdyd1.missing$missingness";
	my $discard_indv = "perl discard_individuals_based_on_missingness_v3.pl -i $filenm_no_ext.$wdyd1.vcf -r $reference -m $missingness -p $ploidy -o $filenm_no_ext.$wdyd2.vcf -v $filenm_no_ext.list_of_deleted_indv";
	print ("Discard indv with missingness inf to $missingness % on file $filenm_no_ext.$wdyd1.vcf\n");
	system("$discard_indv");

	#how many individuals were deleted?
	my $line_count = `wc -l < $filenm_no_ext.list_of_deleted_indv`;
	my $an = (($indv_count - $line_count)*((100-$allowed_missing_percent)/100));
	$an = (ceil($an))*2;
	my $unmoinsmaf = 1 - $maf;

	#gatk SelectVariants
	$wdyd1 = $wdyd2;
	$wdyd2 = "$wdyd1.missing$allowed_missing_percent.maf$maf.biallelic.inbreedingcoeffsup$mininbreedingcoeff";
	my $gatk_filter = "$java_path/java -jar $gatk_path/GenomeAnalysisTK.jar -nt 4 -T SelectVariants -R $reference -V $filenm_no_ext.$wdyd1.vcf --restrictAllelesTo BIALLELIC --selectexpressions '(AF >= $maf && AF <=  $unmoinsmaf) && AN>=$an && InbreedingCoeff > $mininbreedingcoeff' --excludeNonVariants -o $filenm_no_ext.$wdyd2.vcf";
	print ("Select variants --selectexpressions '(AF >= $maf && AF <=  $unmoinsmaf) && AN >= $an && InbreedingCoeff > $mininbreedingcoeff' --excludeNonVariants $filenm_no_ext.$wdyd1.vcf...\n");
	system("$gatk_filter");


	#Set to missing genotype each genotype that has a read depth below a defined parameter. It removes also the sites having a number of missing genotypes (for all the individuals) above another defined parameter
	$wdyd1 = $wdyd2;
	$wdyd2 = "$wdyd1.rdmin$readdepthlimitmin.rdmax$readdepthlimitmax.rdallelemin$readdepthmin_allele";
	my $set_missing = "perl set_missing_genotypes_with_readdepth_on_vcf_new.pl -i $filenm_no_ext.$wdyd1.vcf -l $readdepthlimitmin -m $readdepthlimitmax -n $readdepthmin_allele -d $maxindivwithmissingdata -o $filenm_no_ext.$wdyd2.vcf";
	print ("perl set_missing_genotypes_with_readdepth_on_vcf.pl -i $filenm_no_ext.$wdyd1.vcf -l $readdepthlimitmin -m $readdepthlimitmax -n $readdepthmin_allele -d $maxindivwithmissingdata -o $filenm_no_ext.$wdyd2.vcf ...\n");
	system("$set_missing");

	#Beagle 4.0
	$wdyd1 = $wdyd2;
	$wdyd2 = "$wdyd1.beagle";
	my $beagle = "$java_path/java -Xmx4g -jar beagle.r1399.jar gt=$filenm_no_ext.$wdyd1.vcf phase-its=$phaseits impute-its=$imputeits out=$filenm_no_ext.$wdyd2";
	print ("java -Xmx4g -jar beagle.r1399.jar gt=$filenm_no_ext.$wdyd1.vcf phase-its=$phaseits impute-its=$imputeits out=$filenm_no_ext.$wdyd2...\n");
	system("$beagle");

	#unzip
	my $unzip = "gunzip $filenm_no_ext.$wdyd2.vcf.gz";
	print ("gunzip $filenm_no_ext.$wdyd2.vcf.gz...\n");
	system("$unzip");

	#Changed to hapmap
	my $hapmp = "perl from_beagle_phased_vcf_file_to_hapmap_file.pl -i $filenm_no_ext.$wdyd2.vcf -o $filenm_no_ext.$wdyd2.hapmap -d $filenm_no_ext.how_del_was_coded.txt";
	print ("Changed to hapmap file $filenm_no_ext.$wdyd2.vcf...\n");
	system("$hapmp");
}


print ("Done..\n");

exit 1;
