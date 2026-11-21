"""Read-only integrity and saved-result consistency checks; no scientific analysis.

Run from any directory: python tools/verify_collection.py --report /path/report.json
Only an explicitly supplied report path is written. Exit 0: no FAIL; exit 1: FAIL.
"""
from pathlib import Path
import argparse,collections,csv,hashlib,io,json,math,zipfile

ROOT=Path(__file__).resolve().parents[1]
ARCHIVE_SHA='f52792140f945eb3d124c367afabc166de5194f58fe3d864ad79de509106293b'
def rows(path):
    with path.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def truth(v):return str(v).lower() in ('true','1')
def number(v):
    if v in ('','NA','NaN',None):return None
    return float(v)
def require(condition,message):
    if not condition:raise ValueError(message)
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--report',type=Path);args=ap.parse_args()
    checks=[]
    def run(name,fn):
        try:detail=fn();checks.append({'check':name,'status':'PASS','detail':detail})
        except Exception as e:checks.append({'check':name,'status':'FAIL','detail':f'{type(e).__name__}: {e}'})
    def manifest():
        m=rows(ROOT/'File_Manifest.tsv')
        for r in m:
            p=(ROOT/r['Path']).resolve();require(p.is_relative_to(ROOT),'Unsafe manifest path')
            require(p.is_file(),r['Path']+' missing')
            require(p.stat().st_size==int(r['Bytes']),r['Path']+' size changed')
            require(hashlib.sha256(p.read_bytes()).hexdigest()==r['SHA256'],r['Path']+' checksum mismatch')
        actual={p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*') if p.is_file()
                and not set(p.relative_to(ROOT).parts)&{'.git','__pycache__','output'} and p.name!='File_Manifest.tsv'}
        expected={r['Path'] for r in m}
        require(actual==expected,'Manifest coverage mismatch: extra='+str(sorted(actual-expected))+'; missing='+str(sorted(expected-actual)))
        return {'verified_files':len(m),'scope':'listed repository bytes, excluding manifest itself'}
    run('Repository file manifest',manifest)
    def archive():
        p=ROOT/'reproducibility/FUSION_GRCh38_Reproducibility.zip'
        require(p.stat().st_size==810699,'Archive size changed')
        require(hashlib.sha256(p.read_bytes()).hexdigest()==ARCHIVE_SHA,'Archive checksum changed')
        with zipfile.ZipFile(p) as z:
            require(len(z.namelist())==58,'Expected 58 archived files')
            require(z.testzip() is None,'Archive CRC failed')
            m=list(csv.DictReader(io.StringIO(z.read('SHA256_manifest.tsv').decode('utf-8-sig')),delimiter='\t'))
            require(len(m)==57,'Expected 57 manifest entries')
            for r in m:
                b=z.read(r['Relative_file_path'])
                require(len(b)==int(r['File_size']),'Archived size mismatch')
                require(hashlib.sha256(b).hexdigest()==r['SHA256'],'Archived checksum mismatch')
        return {'archive_files':58,'manifest_entries':57,'sha256':ARCHIVE_SHA}
    run('Unchanged FUSION archive and internal manifest',archive)
    def summaries():
        counts=rows(ROOT/'data/final_summaries/S5_Test_counts.tsv')
        for method,tests,assoc,genes in [('JTI',296168,264,64),('FUSION',167057,561,128)]:
            r=next(r for r in counts if r['Analysis']==method and r['Scope']=='All 29 tissues')
            require([int(r[k]) for k in ['Valid tests / input units','FDR-significant / retained units','Genes represented']]==[tests,assoc,genes],method+' counts')
        g=rows(ROOT/'data/final_summaries/S6_Integrated.tsv')
        require(len(g)==len({r['Ensembl_ID'] for r in g})==152,'Union IDs/count')
        j={r['Ensembl_ID'] for r in g if truth(r['JTI_significant'])};f={r['Ensembl_ID'] for r in g if truth(r['FUSION_significant'])}
        require((len(j),len(f),len(j&f),len(j|f))==(64,128,40,152),'Gene sets')
        c=rows(ROOT/'data/final_summaries/S10A_Coloc_primary.tsv')
        require(len(c)==len({r['Probe'] for r in c})==23,'Expected 23 unique probes')
        for r in c:
            require(math.isclose(sum(float(r['PP_H'+str(i)]) for i in range(5)),1,abs_tol=1e-10),'Posterior sum')
            strong=float(r['PP_H4'])>=.8 and float(r['PP_H4'])>3*float(r['PP_H3'])
            require(strong==truth(r['Strong_colocalization']),'Stored strong-coloc flag')
        require(not any(truth(r['Strong_colocalization']) for r in c),'Unexpected strong coloc')
        n=rows(ROOT/'data/final_summaries/S11D_Tissue_summary.tsv')
        require(len(n)==16 and sum(int(r['N_FDR_lt_0.05']) for r in n)==59,'Nonbrain summaries')
        regional=rows(ROOT/'data/final_summaries/S9C_Sensitivity.tsv')
        # The source sheet also contains explanatory records below the four analyses.
        # Preserve those exports; select actual analyses by their original row labels.
        names={'Core','Extended','Core all-testable sensitivity','Extended all-testable sensitivity'}
        regional=[r for r in regional if r['Analysis'] in names]
        require([(int(r['Input_models']),int(r['Unique_genes'])) for r in regional]==[(353,28),(358,33),(559,32),(1011,82)],'Regional model/gene counts')
        return {'JTI':[264,64],'FUSION':[561,128],'shared':40,'union':152,'saved_coloc_probes':23,'scope':'saved workbook exports, not recalculated analyses'}
    run('Final saved TWAS/coloc/nonbrain summaries',summaries)
    def huang():
        d=ROOT/'data/external_evidence';a=rows(d/'S13A_Molecular_convergence.tsv');b=rows(d/'S13B_Evidence_details.tsv')
        trace=rows(d/'Huang_S14_Source_Trace.tsv');strict=rows(d/'Huang_S14_Strict5.tsv')
        require(len(a)==len(b)==152,'S13 row counts')
        require(len(trace)==22 and len({r['Gene'] for r in trace})==10,'S14 source rows')
        require(len({r['Original_S1_row'] for r in trace})==22,'Duplicate source row')
        require({r['Gene'] for r in trace}=={r['Gene'] for r in a if truth(r['Huang2021_support'])},'Huang broad overlap')
        actual=[sum(truth(r[c]) for r in a) for c in ['Huang2021_support','Matoba2020_support','Gandal2022_support']]
        require(actual==[10,67,22],'Three source overlaps')
        dist=collections.Counter(int(r['N_Huang_Matoba_Gandal_sources']) for r in a)
        for r in a:
            require(int(r['N_Huang_Matoba_Gandal_sources'])==sum(truth(r[c]) for c in ['Huang2021_support','Matoba2020_support','Gandal2022_support']),'Three-source flag count')
        require([dist[i] for i in range(4)]==[70,67,13,2],'Three-source distribution')
        require({r['Gene'] for r in a if int(r['N_Huang_Matoba_Gandal_sources'])==3}=={'MAPT','XRN2'},'Triple genes')
        require(len(strict)==5 and {r['Gene'] for r in strict}=={'POU3F2','MSRA','MAPT','KIZ','NKX2-2'},'Strict five')
        for r in trace:require(truth(r['Strict_combined_significant'])==(float(r['Combined_meta_P'])<4e-7),'Strict flag vs combined P')
        require(all(float(r['Combined_meta_P'])<4e-7 for r in strict),'Strict source P')
        strict_overlap={r['Gene'] for r in strict if truth(r['Overlap_152_gene_union'])}
        require(strict_overlap=={'MAPT','KIZ'},'Strict overlap')
        support=[int(r['Gene'] in strict_overlap)+int(truth(r['Matoba2020_support']))+int(truth(r['Gandal2022_support'])) for r in a]
        require([sum(x>=1 for x in support),sum(x>=2 for x in support),sum(x==3 for x in support)]==[80,10,1],'Strict sensitivity saved-flag counts')
        require({r['Gene'] for r,n in zip(a,support) if n==3}=={'MAPT'},'Strict triple gene')
        return {'broad_overlap':10,'source_rows':22,'at_least_one':82,'at_least_two':15,'triple':['MAPT','XRN2'],'strict':[80,10,1],'scope':'saved flags/source records only; no rematching or FDR adjustment'}
    run('Huang source trace and complementary evidence',huang)
    def grids():
        root=ROOT/'scripts/10_final_figures';out=[]
        for code,method,group,ng,nt,stars,missing in [('Figure2A','JTI','Brain',40,13,152,238),('Figure2B','JTI','Nonbrain',40,16,115,312),('Figure3A','FUSION','Brain',76,13,263,468),('Figure3B','FUSION','Nonbrain',85,16,295,755)]:
            d=root/'circular_heatmaps/data';v=rows(d/f'{code}_{method}_{group}_FINAL_plotdata.tsv');g=rows(d/f'{code}_gene_order.tsv');t=rows(d/f'{code}_tissue_order.tsv')
            require(len(v)==ng*nt and len({(r['ENSEMBL_ID'],r['TISSUE']) for r in v})==ng*nt,code+' grid')
            require(len(g)==ng and len(t)==nt and {r['gene'] for r in g}=={r['DISPLAY_GENE'] for r in v} and {r['tissue'] for r in t}=={r['TISSUE'] for r in v},code+' orders')
            fdr='FDR_BRAIN' if group=='Brain' else 'FDR_NONBRAIN'
            require(sum(r['STAR']=='*' for r in v)==stars and sum(number(r['TWAS_Z']) is None for r in v)==missing,code+' stars/NA')
            for r in v:require((r['STAR']=='*')==(number(r[fdr]) is not None and number(r[fdr])<.05),code+' group star')
            out.append({'panel':code,'cells':len(v),'stars':stars,'missing_Z':missing})
        for name,ng,zcol in [('Figure_2D_JTI_shared21_alltissueFDR_plotdata.tsv',21,'zscore'),('Figure_3D_FUSION_shared32_alltissueFDR_plotdata.tsv',32,'TWAS.Z')]:
            v=rows(root/'cross_tissue_D/data'/name)
            require(len(v)==ng*29 and len({(r['ENSEMBL_ID'],r['TISSUE']) for r in v})==ng*29,name+' grid')
            for r in v:require(truth(r['SIG_ALL'])==(number(r['FDR_ALL']) is not None and number(r['FDR_ALL'])<.05),name+' all-tissue star')
        return out
    run('Plot inputs and distinct FDR definitions',grids)
    checks.append({'check':'Reproduction boundary','status':'WARNING','detail':'No TWAS/SMR/FUSION/coloc analysis or new drawing was run. Archived-source scripts need external inputs. S2 revised-layout adoption remains unconfirmed.'})
    result={'status':'FAIL' if any(c['status']=='FAIL' for c in checks) else 'PASS_WITH_WARNINGS','checks':checks,'scientific_analyses_executed':False}
    if args.report:
        args.report.parent.mkdir(parents=True,exist_ok=True);args.report.write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2));return 1 if result['status']=='FAIL' else 0
if __name__=='__main__':raise SystemExit(main())
