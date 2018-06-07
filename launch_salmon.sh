#!/bin/bash

for FILE in *.fastq.gz
do
/home/yhueber/bin/Salmon-0.8.2_linux_x86_64/bin/salmon quant -i ../cdna_index -l SF -g ../musa_acuminata_v2_2.gff3 -r $FILE -o "cdna_quant""$FILE"
done
