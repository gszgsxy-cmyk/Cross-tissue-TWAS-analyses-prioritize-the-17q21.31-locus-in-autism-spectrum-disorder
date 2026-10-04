"""Materialize a separate path-configured copy. Does not run any analysis."""
from pathlib import Path
import argparse,json,re,shutil

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config',required=True)
    ap.add_argument('--output',required=True)
    args=ap.parse_args()
    root=Path(__file__).resolve().parents[1]
    out=Path(args.output).resolve()
    if out.exists(): raise SystemExit('Output already exists; choose a new directory.')
    values=json.loads(Path(args.config).read_text(encoding='utf-8-sig'))
    for key,value in values.items():
        if not isinstance(value,str) or not value.strip(): raise SystemExit('Missing path value: '+key)
        if any(c in value for c in ['"', "'", '$', '`', '\n', '\r']): raise SystemExit('Unsafe path characters: '+key)
        values[key]=value.replace('\\','/')
    sources=[p for p in root.rglob('*') if p.is_file() and not set(p.relative_to(root).parts)&{'.git','configured','run','__pycache__'} and not p.name.endswith('.local.json')]
    planned=[]
    for p in sources:
        raw=p.read_bytes()
        try: text=raw.decode('utf8')
        except UnicodeDecodeError: planned.append((p,raw));continue
        tokens=set(re.findall(r'__([A-Z][A-Z0-9_]*ROOT|USER_HOME)__',text))
        missing=tokens-values.keys()
        if missing: raise SystemExit('Unresolved roots: '+', '.join(sorted(missing)))
        for token in tokens: text=text.replace('__'+token+'__',values[token])
        planned.append((p,text.encode('utf8')))
    out.mkdir(parents=True)
    for p,data in planned:
        target=out/p.relative_to(root);target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(data)
    print('Configured',len(planned),'files in',out)
    print('No analysis executed. Inspect input/output paths before running any script.')
if __name__=='__main__':main()
