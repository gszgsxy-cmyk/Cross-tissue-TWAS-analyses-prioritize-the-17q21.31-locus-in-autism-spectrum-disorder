#!/bin/bash

TISSUES=(
  Colon_Transverse
  Minor_Salivary_Gland
  Pancreas
  Small_Intestine_Terminal_Ileum
  Whole_Blood
)

for tissue in "${TISSUES[@]}"; do

  echo "============================================================"
  echo "STARTING TISSUE: ${tissue}"
  echo "============================================================"

  for chr in {1..22}; do

    echo "------------------------------------------------------------"
    echo "Running ${tissue} chr${chr}"
    echo "------------------------------------------------------------"

    bash "GTEx_nonbrain_SMR_completion_scripts/run_${tissue}_chr.sh" "${chr}"

    status=$?

    if [ ${status} -ne 0 ]; then
      echo "============================================================"
      echo "FAILED: ${tissue} chr${chr}"
      echo "Exit status: ${status}"
      echo "STOPPING MASTER RUN"
      echo "============================================================"
      exit ${status}
    fi

  done

  echo "============================================================"
  echo "FINISHED TISSUE: ${tissue}"
  echo "============================================================"

done

echo "============================================================"
echo "ALL REMAINING 5 GTEx NON-BRAIN TISSUES COMPLETED"
echo "============================================================"
