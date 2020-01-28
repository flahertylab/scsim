#!/bin/bash
source activate bam-readcount
# Input example-snv-locations.csv
pos=$1
# Input filenames of bam file
file=$2

cat $1 |head -100|tr ',' '\t'|awk ' { print "ref_source",$2,$2 }'|tr ' ' '\t' > ${3}pos.bed
for f in `cat ${file}`; do
bam-readcount -w 0 -f ${4}ref_source.fa -l ${3}pos.bed $f > temp.tab
cat temp.tab|awk ' { print $1":"$2,$3,$6,$7,$8,$9}'|tr ':' '\t'| tr ' ' '\t'|awk ' { print $1":"$2,$3,$5,$19,$33,$47}'|tr ' ' '\t' > temp.count
cat temp.count |awk ' { if($3!=0) print $1,$2,"A"; else if($4!=0) print $1,$2,"C"; else if($5!=0) print $1,$2,"G";else  if($6!=0) print $1,$2,"T"}'|awk ' {if ($2==$3) print $1,"1";else print $1,"0"}'|tr ' ' '\t' > $f".mat"; done

cat ${3}bulk*.bam.mat > ${3}BULK_SNPS.tab

source deactivate