#!/usr/bin/perl
=pod

=head1 NAME

blast_cluster.pl 

=head1 SYNOPSIS

blast_cluster.pl --help

=head1 REQUIRES

Perl5

=head1 DESCRIPTION

=cut

use strict;
use Carp qw (cluck confess croak);
use warnings;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use Bio::SeqIO;

=pod

=head1 OPTIONS

Example :

perl blast_cluster.pl --input <file fasta> --program blastn --database <bank> --directory <path> -o <output> 

=head2 Parameters

=over 4

=item B<--input or -i> :

Input file name (<string>)

=item B<--directory> :

Directory output (<string>)

=item B<--program or -p> :

Input program name (<string>)

   blastn
   blastp
   blastx
   tblastn 
   tblastx

Default = `blastn' 

=item B<--database> : 

BLAST database name (<string>)

=item B<--output or -o> :

Output file name (<string>)

=item B<--evalue or -e> :

Expectation value (E) threshold for saving hits (<real>)

Default = `10'

=item B<--format or -f> :

Alignment view options:

   0 = pairwise,
   1 = query-anchored showing identities,
   2 = query-anchored no identities,
   3 = flat query-anchored, show identities,
   4 = flat query-anchored, no identities,
   5 = XML Blast output,
   6 = tabular,
   7 = tabular with comment lines,
   8 = Text ASN.1,
   9 = Binary ASN.1
   10 = Comma-separated values

Default = `0'

=item B<--max_target_seq or -m> : 
   
Maximum number of aligned sequences to keep (<integer>)

Default = `10'

=item B<--num_seq_by_batch or -n> : 

Default = `8'

=item B<--max_thread>

Maximum number of simultaneously running tasks (<integer>)

Default = 24 

=item B<-q> : 

Default = `normal.q'

	normal.q
	long.q
	gnpannot.q
	esttik.q
	greenphyl.q
	arcad.q
	bigmem.q
	
=item B<-dust>
 
 Filter query sequence with DUST (Format: 'yes','level window linker', or 'no' to disable)

=item B<--parse_deflines> :

Should the query and subject defline(s) be parsed?
  0 = no (Default)
  1 = yes

=item B<-num_descriptions <Integer, >=0>

   Number of database sequences to show one-line descriptions for
   Default = '500'

=item B<-num_alignments <Integer, >=0>

   Number of database sequences to show alignments for
   Default = '250'

=item B<-task>  <String, Permissible values: 'blastn' 'blastn-short' 'dc-megablast' 'megablast' 'vecscreen' >

   Task to execute
   Default = 'megablast'

=item B<-clean> <Boolean>
   
   Default = 0

=item B<-advanced_option> <String>

   Advanced option, if you want to add more specific option (<string>) 
   e.g : -advanced_option ' -seg yes '

=item B<--help or -h> : 

Prints a brief help message and exits.

=back

=cut

my $input            = "";
my $program          = "blastn";
my $database         = "";
my $output           = "blast.out";
my $evalue           = 10;
my $format           = 0;
my $directory        = "";
my $queue            = "normal.q";
my $hold_jid         = "";
my $help;
my $man;
my $max_target_seqs  = 8;
my $num_seq_by_batch = 8;
my $cpt_queries      = 0;
my $cpt_file_temp    = 0;
my $cptJob           = 0;
my $parse_deflines   = 0;
my $num_threads      = 4;
my $max_thread       = 16;
my $clean            = 0;
my $task             = 'megablast';
my $restart;
my $task_opt		 = "";
my $dust  = ""           ;
my $advanced_option  = "";
my @seq;
my $outseq;
my $option = "";
my $usage = q/
 blast_cluster.pl [--input <file_fasta>] [--database <db_file>] [--directory <path>]

Parameters
       -input                  Input file name (<string>) [required]
       -directory              Directory to store result [required]
       -database               Input program name (<string>) [required]
       -program                Input db name (<string>), default blastn
       -output		           Output file name (<string>), default blast.out
       -format	               Alignment view option, default 0
       -evalue                 Expectation value (E) threshold for saving hits, default 10
       -max_target_seqs        Maximum number of aligned sequences to keep (<integer>), default 10
       -max_thread             Maximum number of simultaneously running tasks (<integer>), default 24
       -hold_jid               Wait for this job_id before processing the job_array (<integer>)
       -dust 		       	   Filter query sequence with dust
       -task                   Task to execute, default megablast
       -q		          	   default normal.q
       -clean                  Remove repository. Default true (0) 
       -num_seq_by_batch
       -advanced_option        To add more specific option (e.g : -advanced_option ' -seg yes ')
       -parse_deflines
       -help
/;

GetOptions(
	'input=s'         	 => \$input,
	'program=s'       	 => \$program,
	'database=s'		 => \$database,
	'output=s'        	 => \$output,
	'evalue=f'           => \$evalue,
	'dust=s'	         => \$dust,
	'format=s{1,}'           => \$format,
	'directory=s'        => \$directory,
	'max_target_seqs=i'   => \$max_target_seqs,
	'max_thread=i'       => \$max_thread,
	'hold_jid=i'         => \$hold_jid,
	'num_seq_by_batch=i' => \$num_seq_by_batch,
	'task=s'             => \$task,
	'q=s'                => \$queue,
	'parse_deflines'     => \$parse_deflines,
	'advanced_option=s'  => \$advanced_option,
	'restart=s'          => \$restart,
	'clean=i'            => \$clean,
	'man'                => \$man,
	'help|h|?'           => \$help
)  or pod2usage(1);
if ($help) { pod2usage(1); }
if ($man) { pod2usage( -verbose => 2 ); }
if ($input eq "") {
    print "Warn :: --input is empty\nPlease specify a sequence file to be analyzed\n";
    print $usage;
    exit 0;
}

if ($database eq "") {
    print "Warn :: --database is empty\nPlease specify a database\n";
    print $usage;
    exit 0;
}
if ($directory eq "") {
    print "Warn :: --directory is empty\nPlease specify a directory to store result\n";
    print $usage;
    exit 0;
}
if ($parse_deflines == 1) {
    $option .= " -parse_deflines ";
}
if($hold_jid ne "")
{
  my $jid = $hold_jid;
  $hold_jid = "#\$ -hold_jid $jid";  
}
if ($task ne "" && $program eq "blastn") {
    $option .= " -task $task";
}
if ($dust) {
    $option .= " -dust $dust ";
}
my $type_format = (split(/ /,$format))[0];
if ($max_target_seqs){
	if ($type_format < 4) {
		$option .= "  -num_descriptions $max_target_seqs -num_alignments $max_target_seqs ";
	}
	else {
		$option .= " -max_target_seqs $max_target_seqs ";
	}
}
unless ($queue =~ /\.q$/){
	$queue = $queue.".q";
}
mkdir $directory unless -e $directory;
my $temp_dir;
if ($restart) {
    $temp_dir = $directory ."/".$restart;
}
else {
    $temp_dir = $directory ."/temp$$";
    mkdir $temp_dir;
}
my $inseq = Bio::SeqIO->new(
	-file   => $input,
	-format => 'fasta'
);
my $temp_file = $temp_dir ."/input.fasta.temp";
my $start = 1;
while (my $seq = $inseq->next_seq) {
    $cpt_queries++;
    if ($cpt_queries == 1) {
    	$cpt_file_temp++;
		my $file_out = $temp_file .".". $cpt_file_temp .".blast";
		if (-e $file_out) {$start++;}
		$outseq  =  Bio::SeqIO->new(
			-file =>">$temp_file".".$cpt_file_temp",
			-format =>'fasta'
		);
	}
	$outseq-> write_seq($seq);
	if ($cpt_queries == $num_seq_by_batch) {
		$cpt_queries = 0; 
	}	  
}

#\$ -pe parallel_smp 4n
#open Sav, ">$temp_dir/decoupeJobBlast.sh" or die "pas possible\n";
 

open(OUT,">$temp_dir/decoupeJobBlast.sh") or die "Unable to write file\n";
 
print OUT "
#!/bin/bash
 
#\$ -N decoupeJobBlast
#\$ -cwd
#\$ -e $temp_dir
#\$ -o $temp_dir
#\$ -q $queue
#\$ -t $start-$cpt_file_temp
#\$ -tc $max_thread
#\$ -S /bin/bash
#\$ -pe parallel_smp $num_threads
$hold_jid

/usr/local/bioinfo/ncbi-blast/2.2.30/bin/". $program ." -num_threads ". $num_threads ." -query ". $temp_dir."/input.fasta.temp.\$SGE_TASK_ID  -db ". $database ." -evalue ".  $evalue ." -outfmt '". $format ."' -out ". $temp_dir ."/input.fasta.temp.\$SGE_TASK_ID.blast ". $option ." " . $advanced_option .";";


close OUT;

system ("qsub -sync y $temp_dir/decoupeJobBlast.sh");

if ($type_format == 5){

	use XML::DOM;
	open Sav, ">$directory/$output" or die "pas possible\n";
	my $start =<<"start"; 
<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "NCBI_BlastOutput.dtd">
<BlastOutput>
  <BlastOutput_program>blastx</BlastOutput_program>
  <BlastOutput_version>BLASTX 2.2.24+</BlastOutput_version>
  <BlastOutput_reference>Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&amp;auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and D
avid J. Lipman (1997), &quot;Gapped BLAST and PSI-BLAST: a new generation of protein database search programs&quot;, Nucleic Acids Res. 25:338
9-3402.</BlastOutput_reference>
  <BlastOutput_db>/bank/uniprot/uniprot_sprot_XA.fasta</BlastOutput_db>
  <BlastOutput_query-ID>Query_1</BlastOutput_query-ID>
  <BlastOutput_query-def>Query</BlastOutput_query-def>
  <BlastOutput_query-len>667</BlastOutput_query-len>
  <BlastOutput_param>
    <Parameters>
      <Parameters_matrix>BLOSUM62</Parameters_matrix>
      <Parameters_expect>$evalue</Parameters_expect>
      <Parameters_gap-open>11</Parameters_gap-open>
      <Parameters_gap-extend>1</Parameters_gap-extend>
      <Parameters_filter>L;</Parameters_filter>
    </Parameters>
  </BlastOutput_param>
  <BlastOutput_iterations>
start

	print Sav $start;
	
    opendir(DIR, "$temp_dir") || die "can't open $temp_dir\n"; 
   	my @dir= grep { /^[^\.]/ && /\.blast$/ && -r "$temp_dir/$_" } readdir(DIR);
   	closedir(DIR);
   
   	for (@dir)   {
		my $file = $temp_dir."/".$_;
		open File, "$file";
		my @lignes = <File>;
		my $first_it = 0;
		foreach my $ligne(@lignes){
			if ($ligne =~/.Iteration./ ){
				$first_it = 1;
			}
			if ($first_it eq 1 && $ligne!~/\/BlastOutput_iterations/ && $ligne!~/\/BlastOutput/){
				print Sav $ligne;
			}
		}
	  
   	}
  	my $end ="  </BlastOutput_iterations>\n</BlastOutput>";
 	print Sav $end;
 	
	close Sav;  
  



}
else{
 	system ("cat $temp_dir/input.fasta.temp.*.blast > $directory/$output");
 	wait;
}

if ($clean == 0){
	system ("rm -Rf $temp_dir");
}
=pod

=head1 AUTHORS

Xavier ARGOUT (CIRAD), xavier.argout@cirad.fr
Gaetan DROC (CIRAD), gaetan.droc@cirad.fr

=head1 VERSION

Version [version.1.0]

Date [20/04/2015]

=head1 SEE ALSO

Not Available

=cut
