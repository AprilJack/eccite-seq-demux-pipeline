#!/usr/bin/bash
for sample in sample1_FB_name sample2_FB_name;
do echo $sample;
cellranger count --id=${sample}_out --libraries=${sample}_library.csv --transcriptome=/gpfs/tools/cellranger/refdata-gex-mm10-2020-A  --feature-ref=hash_FB.csv --localcores=12 --localmem=96 --expect-cells 10000
done;
