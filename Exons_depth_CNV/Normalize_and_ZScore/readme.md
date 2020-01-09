
calculate depth of each chromosome

```
for inbam in `ls /data/test_bams/*.bam`
do
    sample=$(basename ${inbam} | cut -f1 -d'.')
    samtools idxstats ${inbam} 1>${sample}.idx.stats.xls
done
```

get Exon length

```
python get_E_len.py
```

calculate Uniformity
```
python test_Uniformity.py
```

Normalize or Not normalize, then calculate ZScore
```
python calc_Zscore.py

or

python normalize_depth.py
python calc_Zscore_normal.py


```

filter Exons by Z_score

```
python pick_z_exons.py
```
