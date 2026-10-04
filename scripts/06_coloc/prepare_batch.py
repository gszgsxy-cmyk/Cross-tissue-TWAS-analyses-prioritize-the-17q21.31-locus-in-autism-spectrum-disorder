from pathlib import Path
import csv,hashlib
ROOT=Path(__file__).resolve().parents[3]
OUT=Path(__file__).resolve().parents[1]
def read(p):
 with open(p,encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rows):
 with open(p,'w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def sha(p):
 with open(p,'rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
for d in ['01_inventory','02_BrainMeta_queries','03_probe_results','04_summary','05_logs']:(OUT/d).mkdir(exist_ok=True)
src=ROOT/'outputs/coloc_input_audit_20260923/01_BrainMeta_coloc_input_completeness.tsv'
raw=read(src)
rows=[dict(Gene=r['gene'],Probe=r['BrainMeta_probe_ID'],Probe_BP=r['probe_bp'],BrainMeta_source_status=r['classification'],Source_inventory_file=str(src),Expected_stored_variants=r['source_BESD_total_variants']) for r in raw]
write(OUT/'01_inventory/01_final_23probe_inventory.tsv',rows)
assert len(rows)==23 and len({r['Probe'] for r in rows})==23,'STOP: inventory not exactly 23 unique probes'
assert all(r['BrainMeta_source_status'] in ('A','PASS') for r in rows),'STOP: completeness conflict'
assert all(r['query_preserves_all_stored_variants'].lower()=='true' and r['query_beta_SE_P_alleles_match_source'].lower()=='true' for r in raw)
protected=[]
for folder in ['MAPT_AS1_ABF_validation_20260924','MAPT_AS1_gene_centric_coloc_pilot_20260924','MAPT_AS1_susie_exploratory_20260924','coloc_input_audit_20260923','coloc_core_fixedlead_20260923','coloc_expanded_20260919']:
 protected += [p for p in (ROOT/'outputs'/folder).rglob('*') if p.is_file()]
protected += [p for p in ROOT.rglob('*') if p.is_file() and p.suffix.lower() in ('.docx','.xlsx') and OUT not in p.parents]
brain=Path('__SMR_ROOT__/sx/BrainMeta_cis_2/BrainMeta_cis_eQTL_chr17')
protected += [Path(str(brain)+s) for s in ('.besd','.epi','.esi')]
protected += [Path('__SMR_ROOT__/asd-gwas/gwas_chr17.txt')]
write(OUT/'05_logs/protected_files_before.tsv',[dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(set(protected))])
base='__WORKSPACE_WSL_ROOT__/outputs/BrainMeta_fullcis_coloc_final_20260924'
lines=['#!/usr/bin/env bash','set -uo pipefail',f"base='{base}'",'printf "Probe\\tExit_code\\n" > "$base/05_logs/query_status.tsv"']
for r in rows:
 probe=r['Probe'];prefix=f'$base/02_BrainMeta_queries/BrainMeta_{probe}'
 lines += [f'test ! -e "{prefix}.txt" || exit 90',f'__SMR_WSL_ROOT__/smr --beqtl-summary __SMR_WSL_ROOT__/sx/BrainMeta_cis_2/BrainMeta_cis_eQTL_chr17 --query 1 --probe {probe} --out "{prefix}" > "{prefix}.stdout.log" 2>&1',f'printf "{probe}\\t%s\\n" "$?" >> "$base/05_logs/query_status.tsv"']
(OUT/'00_scripts/extract_all.sh').write_text('\n'.join(lines)+'\n',encoding='utf-8',newline='\n')
print('23 A-class probes validated; protected hashes saved; extraction script ready.')
