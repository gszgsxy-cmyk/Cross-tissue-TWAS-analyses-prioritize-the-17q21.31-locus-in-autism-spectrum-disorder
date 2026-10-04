#!/bin/bash
set -euo pipefail

# Usage:
#   bash run_Whole_Blood_chr.sh <CHR>
# Example:
#   bash run_Whole_Blood_chr.sh 1

CHR="${1:-}"
TISSUE="Whole_Blood"

if [[ -z "$CHR" ]]; then
  echo "ERROR: chromosome number is required (1-22)."
  echo "Usage: bash run_Whole_Blood_chr.sh <CHR>"
  exit 1
fi

if ! [[ "$CHR" =~ ^([1-9]|1[0-9]|2[0-2])$ ]]; then
  echo "ERROR: CHR must be an integer from 1 to 22."
  exit 1
fi

cd "__SMR_WSL_ROOT__"

OUTDIR="./GTEX_nonbrain_completion_2026-09-18/${TISSUE}"
mkdir -p "$OUTDIR"

OUTPREFIX="${OUTDIR}/GTEX_${TISSUE}_chr${CHR}"
LOGFILE="${OUTPREFIX}.log"

echo "============================================================" | tee "$LOGFILE"
echo "GTEx non-brain SMR completion run" | tee -a "$LOGFILE"
echo "Tissue: $TISSUE" | tee -a "$LOGFILE"
echo "Chromosome: $CHR" | tee -a "$LOGFILE"
echo "eQTL: ./sx/eQTL_besd_lite/${TISSUE}.lite" | tee -a "$LOGFILE"
echo "GWAS: ./asd-gwas/gwas_chr${CHR}.txt" | tee -a "$LOGFILE"
echo "LD reference: ./g1000_eur/g1000_eur" | tee -a "$LOGFILE"
echo "Output prefix: $OUTPREFIX" | tee -a "$LOGFILE"
echo "============================================================" | tee -a "$LOGFILE"

./smr \
  --beqtl-summary "./sx/eQTL_besd_lite/${TISSUE}.lite" \
  --bfile "./g1000_eur/g1000_eur" \
  --gwas-summary "./asd-gwas/gwas_chr${CHR}.txt" \
  --out "$OUTPREFIX" \
  --thread-num 10 \
  2>&1 | tee -a "$LOGFILE"

echo "============================================================" | tee -a "$LOGFILE"
echo "Finished: $TISSUE chr$CHR" | tee -a "$LOGFILE"
echo "============================================================" | tee -a "$LOGFILE"
