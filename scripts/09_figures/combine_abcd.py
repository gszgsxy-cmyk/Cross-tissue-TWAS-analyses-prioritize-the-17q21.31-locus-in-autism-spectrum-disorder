from pathlib import Path
import io,json,subprocess,shutil
import pdfplumber
from pypdf import PdfReader,PdfWriter,Transformation
from pypdf.generic import RectangleObject
from reportlab.pdfgen import canvas
from PIL import Image
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'outputs/Figure4_ABCD_final_20260928';OUT.mkdir(exist_ok=True)
src=ROOT/'outputs/figure_update_20260920/source_figures'
inputs=[src/'Figure3A_brain_eQTL.pdf',src/'Figure3B_brain_mQTL.pdf',ROOT/'outputs/Figure4C_Hatton_original_style_20260928/Figure_4C_Hatton_EUR_mQTL.pdf',ROOT/'outputs/Figure4D_Venn_20260928/Figure_4D_CpG_overlap.pdf']
titles=['A   BrainMeta brain eQTL-SMR','B   Brain-mMeta brain mQTL-SMR','C   Hatton EUR whole-blood mQTL-SMR']
items=[];width=900
for i,f in enumerate(inputs):
 p=PdfReader(f).pages[0]
 if i<3:
  with pdfplumber.open(f) as doc:
   chars=[z for z in doc.pages[0].chars if z['top']>60]
   bounds=(min(z['x0'] for z in chars)-9,min(z['y0'] for z in chars)-10,max(z['x1'] for z in chars)+9,max(z['y1'] for z in chars)+10)
 else:bounds=(0,0,float(p.mediabox.width),float(p.mediabox.height))
 p.cropbox=RectangleObject(bounds);p.trimbox=RectangleObject(bounds)
 scale=(width-36)/(bounds[2]-bounds[0]);height=(bounds[3]-bounds[1])*scale
 items.append((p,bounds,scale,height))
height=sum(z[3]+(36 if i<3 else 0) for i,z in enumerate(items))+36
w=PdfWriter();page=w.add_blank_page(width,height)
buf=io.BytesIO();c=canvas.Canvas(buf,pagesize=(width,height));y=height-18
for i,(p,b,s,h) in enumerate(items):
 if i<3:
  c.setFont('Helvetica-Bold',13);c.drawString(20,y-16,titles[i]);y-=36
 y-=h
 page.merge_transformed_page(p,Transformation().translate(-b[0],-b[1]).scale(s).translate(18,y))
c.save();page.merge_page(PdfReader(buf).pages[0])
target=OUT/'Figure_4_ABCD.pdf'
with target.open('wb') as f:w.write(f)
pop='__USER_HOME__/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/Library/bin/pdftoppm.exe'
subprocess.run([pop,'-singlefile','-r','300','-png',str(target),str(OUT/'render')],check=True)
Image.MAX_IMAGE_PIXELS=None
im=Image.open(OUT/'render.png').convert('RGB');im.save(OUT/'Figure_4_ABCD.tiff',compression='tiff_lzw',dpi=(300,300))
im.thumbnail((1500,3200));im.save(OUT/'Figure_4_ABCD_preview.png');(OUT/'render.png').unlink()
with pdfplumber.open(target) as doc:
 text=doc.pages[0].extract_text()
 assert all(x in text for x in ['GABBR1','AL133325.3','cg25635805','cg00968488','61,780','21 / 30'])
 for i,f in enumerate(inputs[:3]):
  with pdfplumber.open(f) as source:
   probes=[l.split()[0] for l in source.pages[0].extract_text().splitlines() if 'rs' in l]
   assert len(probes)==[22,21,30][i]
   assert all(x in text for x in probes)
shutil.copy2(__file__,OUT/'combine_abcd.py')
(OUT/'source_manifest.json').write_text(json.dumps({'sources':[str(f) for f in inputs],'layout':'A-B-C-D vertically stacked; vector content retained','PDF_size_pt':[width,height],'TIFF_dpi':300,'verified_rows':[22,21,30]},ensure_ascii=False,indent=2),encoding='utf-8')
print('Saved and verified all 22/21/30 forest rows and D summary.')
