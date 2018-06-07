#Script written by Yann Hueber

use warnings;
use strict;
use Getopt::Long;
use Bio::Seq;
use Bio::SeqIO;

my $fasta_filename = $ARGV[0];
my $chr = $ARGV[1];
my $pos = $ARGV[2];

#retrieve seqIO object with reference file and charge it in memory
my %seq_hash;
my $seqio_obj = Bio::SeqIO->new(-file => "<$fasta_filename", -format => "fasta");

while (my $seq_obj = $seqio_obj->next_seq){
	my $seq_name = $seq_obj->display_id;
	$seq_hash{$seq_name} = $seq_obj;
}

my $ref_allele = $seq_hash{$chr}->subseq($pos,$pos);


print "Chr $chr position $pos ref allele : $ref_allele\n";


