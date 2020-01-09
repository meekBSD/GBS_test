


1.  python ./pick_Exons_CDS.py ./transcripts_IDs.txt ./refGene.txt > Test_genes_exons.xls

2.  python ./transform_bed.py ./transcripts_IDs.txt > test_exons.bed

## Use R

```grep 'ENST000_N' /data/annotation/database/gencode.v30lift37.annotation.gtf > testGene.bed

awk '{if($3=="exon")print}' testGene.bed > Gene_exons.tsv

samtools depth -d 100000 -b testGene.bed ../sample.dup.bam -a > sample_testGene.txt

python classify_regions.py Gene_exons.tsv sample_testGene.txt > result_dep.txt

Rscript plot_exons.R```
