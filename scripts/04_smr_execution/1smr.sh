#!/bin/bash
for ((i=1; i<=22; i++))
do
    ./smr  --beqtl-summary ./sx/eQTL_besd_lite/Whole_Blood.lite --query 5.0e-8 --snp-chr $i --make-besd --out Whole_Blood.lite_chr$i
done
