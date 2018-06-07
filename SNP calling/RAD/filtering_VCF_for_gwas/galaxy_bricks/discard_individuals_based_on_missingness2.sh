#!/bin/bash


vcf_input=$1
reference=$2
missingness=$3
vcf_output=$4
indv_del_list=$5



mkdir tmpdir$$

cp -rf $vcf_input tmpdir$$/input.vcf
cp -rf $reference tmpdir$$/ref.fa

/usr/bin/perl $HOME/galaxy_dist/tools/filters/discard_individuals_based_on_missingness.pl -i tmpdir$$/input.vcf -r tmpdir$$/ref.fa -m $missingness -o $vcf_output -v $indv_del_list

rm -rf tmpdir$$

