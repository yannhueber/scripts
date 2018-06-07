#! /usr/bin/perl -w
#$ -q bioinfo.q
#$ -cwd
#$ -V



opendir(REPERTOIRE_COURANT,".") or die("Cannot open directory\n");  
my @fichiers = grep { m/\.fastq$/ } readdir(REPERTOIRE_COURANT);
closedir(REPERTOIRE_COURANT);



open(ALL_IN,">concat_fastq") or die ("Cannot open file concat_fastq");

my $i=0;
foreach my $file (@fichiers) {
	open(IN,"<$file") or die ("Cannot open $file");
	print "Concat file number $i : $file\n";
	$i++;
	while(<IN>) {
		print ALL_IN $_;
	}
	close(IN);
}

exit;
