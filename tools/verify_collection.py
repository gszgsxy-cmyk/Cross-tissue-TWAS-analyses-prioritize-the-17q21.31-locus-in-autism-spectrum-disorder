"""Check collected-file integrity and small archived summaries, without rerunning analyses."""
from pathlib import Path
import csv,hashlib,json,math

def main():
    root=Path(__file__).resolve().parents[1]
    with (root/'provenance/source_manifest.tsv').open(encoding='utf8',newline='') as f: rows=list(csv.DictReader(f,delimiter='\t'))
    for row in rows:
        p=root/row['path']
        assert p.is_file(),row['path']
        assert hashlib.sha256(p.read_bytes()).hexdigest()==row['packaged_sha256'],row['path']
    with (root/'reference/expected_results/S10A_Coloc_primary.tsv').open(encoding='utf8',newline='') as f: coloc=list(csv.DictReader(f,delimiter='\t'))
    assert len(coloc)==24
    for row in coloc:
        assert math.isclose(sum(float(row['PP_H'+str(i)]) for i in range(5)),1.0,abs_tol=1e-10)
    with (root/'reference/expected_results/S11D_Tissue_summary.tsv').open(encoding='utf8',newline='') as f: tissue=list(csv.DictReader(f,delimiter='\t'))
    assert len(tissue)==16
    assert sum(int(r['N_total_SMR_tests']) for r in tissue)==59449
    assert sum(int(r['N_FDR_lt_0.05']) for r in tissue)==59
    for r in tissue:
        assert sum(int(r[k]) for k in ['N_HEIDI_no_heterogeneity','N_HEIDI_heterogeneity','N_HEIDI_unavailable'])==int(r['N_FDR_lt_0.05'])
    with (root/'reference/expected_results/S5_Test_counts.tsv').open(encoding='utf8',newline='') as f: counts=list(csv.DictReader(f,delimiter='\t'))
    for method,features,genes in [('JTI',264,64),('FUSION',561,128)]:
        row=next(r for r in counts if r['Analysis']==method and r['Scope']=='All')
        assert int(row['FDR_significant_features'])==features and int(float(row['Unique_significant_genes']))==genes
    print(json.dumps({'collected_files_verified':len(rows),'primary_coloc_rows':len(coloc),'nonbrain_tissues':len(tissue),'status':'PASS','scope':'integrity and archived-summary structure only; no numerical analysis rerun'},indent=2))
if __name__=='__main__':main()
