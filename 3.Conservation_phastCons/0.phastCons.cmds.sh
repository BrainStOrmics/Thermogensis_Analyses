#!/bin/bash
# This script shows the phastCons analysis step by step
set -e
echo "This script shows the phastCons ansylsis step by step"
echo "It is not prepared as a run-able script, please run it line by line"
exit 1

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons470way/hg38.phastCons470way.bw
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigToBedGraph
chmod a+w bigWigToBedGraph
./bigWigToBedGraph hg38.phastCons470way.bw hg38.phastCons470way.bedGraph &

# sort ensembl canonical
# ensembl canonical gtf is prepared from extract_canonical_transcript.R
less ensembl_canonical.cds.clean.gtf 
cat ensembl_canonical.cds.clean.gtf |grep -v "^#"|sort -k1,1 -k2,2n > ensembl_canonical.cds.clean.sorted.gtf
less ensembl_canonical.cds.clean.sorted.gtf

# intersect
./bedtools intersect -a hg38.phastCons470way.bedGraph -b ensembl_canonical.cds.clean.sorted.gtf -sorted > phastCons.cds.bed
./bedtools intersect -b phastCons.cds.bed -a ensembl_canonical.cds.clean.sorted.gtf -sorted -wa -wb > ensembl_canonical.cds.phastCons.gtf

# combine phastCons scores to a list for each gene
python phastCons.cds.merge.py ensembl_canonical.cds.phastCons.gtf ensembl_canonical.cds.phastCons.list.gtf > merge.log
less -S ensembl_canonical.cds.phastCons.list.gtf|grep -- "-1"|wc -l

cat ensembl_canonical.cds.phastCons.list.gtf |sort -t $'\t' -k9,9 -k1,1 -k4,5n > ensembl_canonical.cds.phastCons.list.sort_gene.gtf
python phastCons.gene.merge.py ensembl_canonical.cds.phastCons.list.sort_gene.gtf ensembl_canonical.gene.cds.phastCons.list.tsv 2>&1 |tee gene.merge.log

# stat and plot
# use phastCons.plot.R
