#文件夹里
#!/bin/bash


for ((i=1; i<=22; i++))
do
    ./smr --beqtl-summary ./sx/Brain-mMeta/Brain-mMeta_chr$i --bfile g1000_eur/g1000_eur --gwas-summary asd-gwas/gwas_chr$i.txt --out ./Brain-mMeta-asd/Brain-mMeta-asd_chr$i --thread-num 10
done
