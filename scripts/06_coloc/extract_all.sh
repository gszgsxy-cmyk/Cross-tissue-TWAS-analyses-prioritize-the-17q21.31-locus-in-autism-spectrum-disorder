#!/usr/bin/env bash
set -uo pipefail
base='__WORKSPACE_WSL_ROOT__/outputs/BrainMeta_fullcis_coloc_final_20260924'
printf "Probe\tExit_code\n" > "$base/05_logs/query_status.tsv"
test ! -e "$base/02_BrainMeta_queries/BrainMeta_ENSG00000073969.18.txt" || exit 90
__SMR_WSL_ROOT__/smr --beqtl-summary __SMR_WSL_ROOT__/sx/BrainMeta_cis_2/BrainMeta_cis_eQTL_chr17 --query 1 --probe ENSG00000073969.18 --out "$base/02_BrainMeta_queries/BrainMeta_ENSG00000073969.18" > "$base/02_BrainMeta_queries/BrainMeta_ENSG00000073969.18.stdout.log" 2>&1
printf "ENSG00000073969.18\t%s\n" "$?" >> "$base/05_logs/query_status.tsv"
test ! -e "$base/02_BrainMeta_queries/BrainMeta_ENSG00000108379.10.txt" || exit 90
__SMR_WSL_ROOT__/smr --beqtl-summary __SMR_WSL_ROOT__/sx/BrainMeta_cis_2/BrainMeta_cis_eQTL_chr17 --query 1 --probe ENSG00000108379.10 --out "$base/02_BrainMeta_queries/BrainMeta_ENSG00000108379.10" > "$base/02_BrainMeta_queries/BrainMeta_ENSG00000108379.10.stdout.log" 2>&1
printf "ENSG00000108379.10\t%s\n" "$?" >> "$base/05_logs/query_status.tsv"
test ! -e "$base/02_BrainMeta_queries/BrainMeta_ENSG00000120088.15.txt" || exit 90
__SMR_WSL_ROOT__/smr --beqtl-summary __SMR_WSL_ROOT__/sx/BrainMeta_cis_2/BrainMeta_cis_eQTL_chr17 --query 1 --probe ENSG00000120088.15 --out "$base/02_BrainMeta_queries/BrainMeta_ENSG00000120088.15" > "$base/02_BrainMeta_queries/BrainMeta_ENSG00000120088.15.stdout.log" 2>&1
printf "ENSG00000120088.15\t%s\n" "$?" >> "$base/05_logs/query_status.tsv"
test ! -e "$base/02_BrainMeta_queries/BrainMeta_ENSG00000159314.11.txt" || exit 90
__SMR_WSL_ROOT__/smr --beqtl-summary __SMR_WSL_ROOT__/sx/BrainMeta_cis_2/BrainMeta_cis_eQTL_chr17 --query 1 --probe ENSG00000159314.11 --out "$base/02_BrainMeta_queries/BrainMeta_ENSG00000159314.11" > "$base/02_BrainMeta_queries/BrainMeta_ENSG00000159314.11.stdout.log" 2>&1
printf "ENSG00000159314.11\t%s\n" "$?" >> "$base/05_logs/query_status.tsv"
test ! -e "$base/02_BrainMeta_queries/BrainMeta_ENSG00000176681.14.txt" || exit 90
__SMR_WSL_ROOT__/smr --beqtl-summary __SMR_WSL_ROOT__/sx/BrainMeta_cis_2/BrainMeta_cis_eQTL_chr17 --query 1 --probe ENSG00000176681.14 --out "$base/02_BrainMeta_queries/BrainMeta_ENSG00000176681.14" > "$base/02_BrainMeta_queries/BrainMeta_ENSG00000176681.14.stdout.log" 2>&1
printf "ENSG00000176681.14\t%s\n" "$?" >> "$base/05_logs/query_status.tsv"
test ! -e "$base/02_BrainMeta_queries/BrainMeta_ENSG00000185294.7.txt" || exit 90
__SMR_WSL_ROOT__/smr --beqtl-summary __SMR_WSL_ROOT__/sx/BrainMeta_cis_2/BrainMeta_cis_eQTL_chr17 --query 1 --probe ENSG00000185294.7 --out "$base/02_BrainMeta_queries/BrainMeta_ENSG00000185294.7" > "$base/02_BrainMeta_queries/BrainMeta_ENSG00000185294.7.stdout.log" 2>&1
printf "ENSG00000185294.7\t%s\n" "$?" >> "$base/05_logs/query_status.tsv"
test ! -e "$base/02_BrainMeta_queries/BrainMeta_ENSG00000185829.18.txt" || exit 90
__SMR_WSL_ROOT__/smr --beqtl-summary __SMR_WSL_ROOT__/sx/BrainMeta_cis_2/BrainMeta_cis_eQTL_chr17 --query 1 --probe ENSG00000185829.18 --out "$base/02_BrainMeta_queries/BrainMeta_ENSG00000185829.18" > "$base/02_BrainMeta_queries/BrainMeta_ENSG00000185829.18.stdout.log" 2>&1
printf "ENSG00000185829.18\t%s\n" "$?" >> "$base/05_logs/query_status.tsv"
test ! -e "$base/02_BrainMeta_queries/BrainMeta_ENSG00000186868.17.txt" || exit 90
__SMR_WSL_ROOT__/smr --beqtl-summary __SMR_WSL_ROOT__/sx/BrainMeta_cis_2/BrainMeta_cis_eQTL_chr17 --query 1 --probe ENSG00000186868.17 --out "$base/02_BrainMeta_queries/BrainMeta_ENSG00000186868.17" > "$base/02_BrainMeta_queries/BrainMeta_ENSG00000186868.17.stdout.log" 2>&1
printf "ENSG00000186868.17\t%s\n" "$?" >> "$base/05_logs/query_status.tsv"
test ! -e "$base/02_BrainMeta_queries/BrainMeta_ENSG00000204650.14.txt" || exit 90
__SMR_WSL_ROOT__/smr --beqtl-summary __SMR_WSL_ROOT__/sx/BrainMeta_cis_2/BrainMeta_cis_eQTL_chr17 --query 1 --probe ENSG00000204650.14 --out "$base/02_BrainMeta_queries/BrainMeta_ENSG00000204650.14" > "$base/02_BrainMeta_queries/BrainMeta_ENSG00000204650.14.stdout.log" 2>&1
printf "ENSG00000204650.14\t%s\n" "$?" >> "$base/05_logs/query_status.tsv"
test ! -e "$base/02_BrainMeta_queries/BrainMeta_ENSG00000214401.4.txt" || exit 90
__SMR_WSL_ROOT__/smr --beqtl-summary __SMR_WSL_ROOT__/sx/BrainMeta_cis_2/BrainMeta_cis_eQTL_chr17 --query 1 --probe ENSG00000214401.4 --out "$base/02_BrainMeta_queries/BrainMeta_ENSG00000214401.4" > "$base/02_BrainMeta_queries/BrainMeta_ENSG00000214401.4.stdout.log" 2>&1
printf "ENSG00000214401.4\t%s\n" "$?" >> "$base/05_logs/query_status.tsv"
test ! -e "$base/02_BrainMeta_queries/BrainMeta_ENSG00000214425.8.txt" || exit 90
__SMR_WSL_ROOT__/smr --beqtl-summary __SMR_WSL_ROOT__/sx/BrainMeta_cis_2/BrainMeta_cis_eQTL_chr17 --query 1 --probe ENSG00000214425.8 --out "$base/02_BrainMeta_queries/BrainMeta_ENSG00000214425.8" > "$base/02_BrainMeta_queries/BrainMeta_ENSG00000214425.8.stdout.log" 2>&1
printf "ENSG00000214425.8\t%s\n" "$?" >> "$base/05_logs/query_status.tsv"
test ! -e "$base/02_BrainMeta_queries/BrainMeta_ENSG00000225190.11.txt" || exit 90
__SMR_WSL_ROOT__/smr --beqtl-summary __SMR_WSL_ROOT__/sx/BrainMeta_cis_2/BrainMeta_cis_eQTL_chr17 --query 1 --probe ENSG00000225190.11 --out "$base/02_BrainMeta_queries/BrainMeta_ENSG00000225190.11" > "$base/02_BrainMeta_queries/BrainMeta_ENSG00000225190.11.stdout.log" 2>&1
printf "ENSG00000225190.11\t%s\n" "$?" >> "$base/05_logs/query_status.tsv"
test ! -e "$base/02_BrainMeta_queries/BrainMeta_ENSG00000232300.1.txt" || exit 90
__SMR_WSL_ROOT__/smr --beqtl-summary __SMR_WSL_ROOT__/sx/BrainMeta_cis_2/BrainMeta_cis_eQTL_chr17 --query 1 --probe ENSG00000232300.1 --out "$base/02_BrainMeta_queries/BrainMeta_ENSG00000232300.1" > "$base/02_BrainMeta_queries/BrainMeta_ENSG00000232300.1.stdout.log" 2>&1
printf "ENSG00000232300.1\t%s\n" "$?" >> "$base/05_logs/query_status.tsv"
test ! -e "$base/02_BrainMeta_queries/BrainMeta_ENSG00000238083.8.txt" || exit 90
__SMR_WSL_ROOT__/smr --beqtl-summary __SMR_WSL_ROOT__/sx/BrainMeta_cis_2/BrainMeta_cis_eQTL_chr17 --query 1 --probe ENSG00000238083.8 --out "$base/02_BrainMeta_queries/BrainMeta_ENSG00000238083.8" > "$base/02_BrainMeta_queries/BrainMeta_ENSG00000238083.8.stdout.log" 2>&1
printf "ENSG00000238083.8\t%s\n" "$?" >> "$base/05_logs/query_status.tsv"
test ! -e "$base/02_BrainMeta_queries/BrainMeta_ENSG00000261575.2.txt" || exit 90
__SMR_WSL_ROOT__/smr --beqtl-summary __SMR_WSL_ROOT__/sx/BrainMeta_cis_2/BrainMeta_cis_eQTL_chr17 --query 1 --probe ENSG00000261575.2 --out "$base/02_BrainMeta_queries/BrainMeta_ENSG00000261575.2" > "$base/02_BrainMeta_queries/BrainMeta_ENSG00000261575.2.stdout.log" 2>&1
printf "ENSG00000261575.2\t%s\n" "$?" >> "$base/05_logs/query_status.tsv"
test ! -e "$base/02_BrainMeta_queries/BrainMeta_ENSG00000262500.1.txt" || exit 90
__SMR_WSL_ROOT__/smr --beqtl-summary __SMR_WSL_ROOT__/sx/BrainMeta_cis_2/BrainMeta_cis_eQTL_chr17 --query 1 --probe ENSG00000262500.1 --out "$base/02_BrainMeta_queries/BrainMeta_ENSG00000262500.1" > "$base/02_BrainMeta_queries/BrainMeta_ENSG00000262500.1.stdout.log" 2>&1
printf "ENSG00000262500.1\t%s\n" "$?" >> "$base/05_logs/query_status.tsv"
test ! -e "$base/02_BrainMeta_queries/BrainMeta_ENSG00000262539.1.txt" || exit 90
__SMR_WSL_ROOT__/smr --beqtl-summary __SMR_WSL_ROOT__/sx/BrainMeta_cis_2/BrainMeta_cis_eQTL_chr17 --query 1 --probe ENSG00000262539.1 --out "$base/02_BrainMeta_queries/BrainMeta_ENSG00000262539.1" > "$base/02_BrainMeta_queries/BrainMeta_ENSG00000262539.1.stdout.log" 2>&1
printf "ENSG00000262539.1\t%s\n" "$?" >> "$base/05_logs/query_status.tsv"
test ! -e "$base/02_BrainMeta_queries/BrainMeta_ENSG00000263503.1.txt" || exit 90
__SMR_WSL_ROOT__/smr --beqtl-summary __SMR_WSL_ROOT__/sx/BrainMeta_cis_2/BrainMeta_cis_eQTL_chr17 --query 1 --probe ENSG00000263503.1 --out "$base/02_BrainMeta_queries/BrainMeta_ENSG00000263503.1" > "$base/02_BrainMeta_queries/BrainMeta_ENSG00000263503.1.stdout.log" 2>&1
printf "ENSG00000263503.1\t%s\n" "$?" >> "$base/05_logs/query_status.tsv"
test ! -e "$base/02_BrainMeta_queries/BrainMeta_ENSG00000264070.1.txt" || exit 90
__SMR_WSL_ROOT__/smr --beqtl-summary __SMR_WSL_ROOT__/sx/BrainMeta_cis_2/BrainMeta_cis_eQTL_chr17 --query 1 --probe ENSG00000264070.1 --out "$base/02_BrainMeta_queries/BrainMeta_ENSG00000264070.1" > "$base/02_BrainMeta_queries/BrainMeta_ENSG00000264070.1.stdout.log" 2>&1
printf "ENSG00000264070.1\t%s\n" "$?" >> "$base/05_logs/query_status.tsv"
test ! -e "$base/02_BrainMeta_queries/BrainMeta_ENSG00000264589.4.txt" || exit 90
__SMR_WSL_ROOT__/smr --beqtl-summary __SMR_WSL_ROOT__/sx/BrainMeta_cis_2/BrainMeta_cis_eQTL_chr17 --query 1 --probe ENSG00000264589.4 --out "$base/02_BrainMeta_queries/BrainMeta_ENSG00000264589.4" > "$base/02_BrainMeta_queries/BrainMeta_ENSG00000264589.4.stdout.log" 2>&1
printf "ENSG00000264589.4\t%s\n" "$?" >> "$base/05_logs/query_status.tsv"
test ! -e "$base/02_BrainMeta_queries/BrainMeta_ENSG00000266918.1.txt" || exit 90
__SMR_WSL_ROOT__/smr --beqtl-summary __SMR_WSL_ROOT__/sx/BrainMeta_cis_2/BrainMeta_cis_eQTL_chr17 --query 1 --probe ENSG00000266918.1 --out "$base/02_BrainMeta_queries/BrainMeta_ENSG00000266918.1" > "$base/02_BrainMeta_queries/BrainMeta_ENSG00000266918.1.stdout.log" 2>&1
printf "ENSG00000266918.1\t%s\n" "$?" >> "$base/05_logs/query_status.tsv"
test ! -e "$base/02_BrainMeta_queries/BrainMeta_ENSG00000279685.2.txt" || exit 90
__SMR_WSL_ROOT__/smr --beqtl-summary __SMR_WSL_ROOT__/sx/BrainMeta_cis_2/BrainMeta_cis_eQTL_chr17 --query 1 --probe ENSG00000279685.2 --out "$base/02_BrainMeta_queries/BrainMeta_ENSG00000279685.2" > "$base/02_BrainMeta_queries/BrainMeta_ENSG00000279685.2.stdout.log" 2>&1
printf "ENSG00000279685.2\t%s\n" "$?" >> "$base/05_logs/query_status.tsv"
test ! -e "$base/02_BrainMeta_queries/BrainMeta_ENSG00000280022.1.txt" || exit 90
__SMR_WSL_ROOT__/smr --beqtl-summary __SMR_WSL_ROOT__/sx/BrainMeta_cis_2/BrainMeta_cis_eQTL_chr17 --query 1 --probe ENSG00000280022.1 --out "$base/02_BrainMeta_queries/BrainMeta_ENSG00000280022.1" > "$base/02_BrainMeta_queries/BrainMeta_ENSG00000280022.1.stdout.log" 2>&1
printf "ENSG00000280022.1\t%s\n" "$?" >> "$base/05_logs/query_status.tsv"
