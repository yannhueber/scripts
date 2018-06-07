#! /usr/bin/perl -w
#$ -q bigmem.q
#$ -cwd
#$ -V
#$ -b yes
#$ -N Compare_fastq_paired

#This script will compare two FASTQ files and check their homogeneity in term of forward-reverse.
#It is based on their name
#The non-paired sequences will be printed in a third file, and can be treated latter as single sequences.
#The files contain pairs but are not truely paired, ie the first sequence of the 1st file did not correspond to the first on the 2nd file...
#It will work for Fastq with one line for sequence, one line for qual ONLY


use strict;
use Getopt::Long;
use Time::HiRes qw( time );
use Data::Dumper;

#print the line when they are in buffer
$|=1;

my $start_time =time();
my $courriel="francois.sabot-at-ird.fr";
my ($nomprog) = $0 =~/([^\/]+)$/;
my $MessAbruti ="\nUsage:
\t$nomprog -f forward_file -r reverse_file -of output_forward_file -or ouput_reverse_file -os output_single_file
or
\t$nomprog -forward forward_file -reverse reverse_file

The script will output three files: forward_file.fastq, reverse_file.fastq and single.fastq, in the working dir.
By default, output file has the same name as the input withe 'paired' added before the extension

DESCRIPTION

This script will compare two FASTQ files and check their homogeneity in term of forward-reverse.
It is based on their name
The non-paired sequences will be printed in a third file, and can be treated latter as single sequences.
The files contain pairs but are not truely paired, ie the first sequence of the 1st file did not correspond to the first on the 2nd file...
It will work for Fastq with one line for sequence, one line for qual ONLY
The Illumina file can be in any type (Solexa, 1.3, 1.5, 1.8)

	contact: $courriel\n\n";
	

unless (@ARGV) 
	{
	print "\nType --help for more informations\n\n";
	exit;
	}

my ($forward,$reverse,$output_forward,$output_reverse,$output_single,$help);

GetOptions("help|?|h" => \$help,	
		"f|forward=s"=>\$forward,
		"r|reverse=s"=>\$reverse,
		"of|output_forward=s"=>\$output_forward,
		"or|output_reverse=s"=>\$output_reverse,
		"os|output_single=s"=>\$output_single
		);
			


if ($help){print $MessAbruti; exit;}
my ($forwardname,$reversename, $singlename);
#Creating output names:
if(!defined($output_forward))
{
  my @forwardlist=split/\//,$forward;
  $forwardname=$forwardlist[-1];
  $forwardname=~s/\.\w{1,}$//;
  $forwardname.="_paired.fastq";
}
else
{
  $forwardname=$output_forward;
}

if(!defined($output_reverse))
{
  my @reverselist=split/\//,$reverse;
  $reversename=$reverselist[-1];
  $reversename=~s/\.\w{1,}$//;
  $reversename.="_paired.fastq";
}
else
{
  $reversename=$output_reverse;
}

if(!defined($output_single))
{
  my @forwardlist=split/\//,$forwardname;
  $singlename=$forwardlist[-1];
  $singlename=~s/^(.+)_\d_.*/$1/;
  $singlename.="_single.fastq";
}
else
{
  $singlename=$output_single;
}

print("Creating hash \n");
#Gathering name and modifying them to be able to compare.
my $comheadline="head -n 1 $forward";
my $headline=`$comheadline`;
chomp $headline;
$headline=substr($headline,1,4);

my $comid="grep $headline $forward";
my $list_id_for= `$comid`;
chomp $list_id_for;
my @id_for=split/\n/,$list_id_for;

#free some perl memory
undef($list_id_for);

my %idhash;
foreach my $name (@id_for)
{
  $name=~s/\/\d$//;
  $name=~s/\s\d:.+$//;
  $idhash{$name}++;
  
}

my $duration=time() - $start_time;
print("Hash created. Duration:". getDuration($duration)." \n");

#free some perl memory
undef(@id_for);
print(system("ps aux | grep $$"));

#Reading and outputing
open SINGLE,">$singlename" or die("\nCannot create $singlename file: $!\n");
open REV, $reverse;
open ROK, ">$reversename" or die ("\nCannot create $reversename file: $!\n");
print "\nPrinting rev...\n";
my %id_both;
while (<REV>) #Writing forward correct
    {
    my $in=$_;
    chomp $in;
    my $idhere=$in;
    $idhere=~s/\/\d$//;
    $idhere=~s/\s\d:.+$//;
    $in.="\n".<REV>.<REV>.<REV>;#Pick up the 4 lines
    if (exists($idhash{$idhere})) #They are in the common list
        {
        print ROK $in;
        delete($idhash{$idhere});
        $id_both{$idhere}++;        
        }
    else #They are not in the common list
        {
        print SINGLE $in;
        }
    }
close(REV);
close(ROK);
$duration=time() - $start_time;
print("Reverse written. Duration: ". getDuration($duration)." \n");
print(system("ps aux | grep $$"));

open FOR, $forward;
open FOK, ">$forwardname" or die ("\nCannot create $forwardname file: $!\n");
print "\nPrinting for...\n";
while (<FOR>) #Writing reverse correct
    {
    my $in=$_;
    chomp $in;
    my $idhere=$in;
    $idhere=~s/\/\d$//;
    $idhere=~s/\s\d:.+$//;
    $in.="\n".<FOR>.<FOR>.<FOR>; #Pick up the 4 lines
    if ($id_both{$idhere}) #They are in the common list
        {
        print FOK $in;
        }
    else #They are not in the common list
        {
        if(!exists($idhash{$idhere}))
          {
           print "Check $idhere\n";
          }
        print SINGLE $in;
        }
    }
close(FOR);
close(FOK);
close(SINGLE);

$duration=time() - $start_time;
print("Forward written. Duration: ". getDuration($duration)." \n");
print(system("ps aux | grep $$"));


############################
#
#  SUB
#
############################


sub getDuration
{
    my $duration = shift;
    if (60 > $duration)
    {
        # less than a minute
        return sprintf('%4.2f second(s)', $duration);
    }
    elsif (3600 > $duration)
    {
        # less than an hour
        return sprintf('%d:%04.2f', int($duration/60), ($duration%60));
    }
    elsif (86400 > $duration)
    {
        # less than a day
        return sprintf('%d:%02d:%04.2f', int($duration/3600), int($duration/60)%60, ($duration%60));
    }
    else
    {
        # more than a day
        return sprintf('%d days and %d:%02d:%04.2f', int($duration/86400), int($duration/3600)%24, int($duration/60)%60, ($duration%60));
    }
}


sub modifname { # Remove the 1 or 2 (or any text separated by a space from the name) at the end of the name
    my @list=@_;
    my @out;
    foreach my $name (@list)
        {
        $name=~ s/\/\d$//; # Remove /1 or /2 for Illumina1.3 or 1.5
        $name =~ s/\s\d:.+$//;# Remove 1:..:...: or anything after a space, for Illumina1.8+
        push @out, $name;
        }
    return @out;
    }
    
exit(0);
