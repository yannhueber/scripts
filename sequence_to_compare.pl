#Script written by Yann Hueber
#compare two sequence and show differences


use warnings;
use strict;
use Getopt::Long;
use Bio::Seq;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::Tools::Run::StandAloneBlast;;
use Bio::SimpleAlign;
use Bio::SearchIO::blast;

my $first_seq = $ARGV[0];
my $second_seq = $ARGV[1];

$first_seq  =~ s/[\s\]\[]//g;
$second_seq =~ s/[\s\]\[]//g;

print "first : " . $first_seq . " and second : " . $second_seq . "\n";

#my $run_blast_2seq = "/usr/local/bioinfo/blast/bin/bl2seq -i $first_seq -j $second_seq -p blastn -o output_bl2seq";


my $seq1 = Bio::Seq->new(-seq => $first_seq, alphabet => 'dna' );
my $seq2 = Bio::Seq->new(-seq => $second_seq, alphabet => 'dna' );

print $seq1->subseq(2,2) . "\n"; 

my $factory = Bio::Tools::Run::StandAloneBlast->new('program' => 'blastn', 'outfile' => 'bl2seq.out');

my $bl2seq_report = $factory->bl2seq($seq1, $seq2);



#Use AlignIO.pm to create a SimpleAlign object from the bl2seq report
my $str = Bio::AlignIO->new(-file   => '<bl2seq.out',
                         -format => 'bl2seq');

my $aln = $str->next_aln();



