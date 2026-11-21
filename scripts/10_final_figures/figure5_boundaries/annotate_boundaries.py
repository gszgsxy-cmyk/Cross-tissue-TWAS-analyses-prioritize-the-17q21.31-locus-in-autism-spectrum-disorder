"""Scheme C: annotate the verified vector PDF. No statistical analysis.
Run with Python (pypdf, pdfplumber, reportlab, Pillow, numpy) and Poppler.
Original R source: FUSION.post_process_GRCh38_plot_geneCollapsed.R:498-701.
The R plot depends on unsaved in-memory selection/covariance objects.
This script preserves the completed, glyph-compatible vector base unchanged.
The output corrects only its x-axis unit text from MB to Mb.
"""
from pathlib import Path
import io, json, hashlib, subprocess
import numpy as np
from PIL import Image, ImageChops
import pdfplumber
from pypdf import PdfReader, PdfWriter
from pypdf.generic import TextStringObject
from reportlab.pdfgen import canvas

import argparse
parser=argparse.ArgumentParser(description='Add predefined core boundaries to the archived vector Figure 5; no FUSION analysis.')
parser.add_argument('--output-dir',type=Path,default=Path(__file__).resolve().parent/'output')
parser.add_argument('--check-only',action='store_true')
parser.add_argument('--raster',action='store_true',help='Also render PNG/TIFF with PDFium; PDF remains vector')
args=parser.parse_args()
BASE=Path(__file__).resolve().parent/'data/Figure_5_Conditional_Joint_base.pdf'
OUT=args.output_dir;STEM='Figure_5_Core_Boundaries_Connected'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(BASE)=='8282284bd3f882761ffe35b713f079197b5e688beafaad3adfa9c0c97959e7a2'
with pdfplumber.open(BASE) as d:
 p=d.pages[0];width,height=p.width,p.height
 # Calibrate with the six actual x-axis tick marks, not the page margins.
 horizontal=[l for l in p.lines if abs(l['y1']-l['y0'])<1e-5 and l['x1']-l['x0']>350]
 axis=horizontal[0];axis_y=axis['y0']
 ticks=sorted(l['x0'] for l in p.lines if abs(l['x1']-l['x0'])<1e-5 and abs(l['y1']-axis_y)<.02 and 4<l['y1']-l['y0']<6)
 mb=np.array([45,45.5,46,46.5,47,47.5]);assert len(ticks)==len(mb)
 slope,origin=np.polyfit(mb-45,np.array(ticks),1)
 residual=float(np.max(np.abs(origin+slope*(mb-45)-ticks)))
 assert residual<.006
 xpos=[float(origin+slope*(v-45)) for v in [45.414716,46.836265]]
 # Drawing bounds follow the existing data tracks; the accepted v2 gap
 # connection is added below. Keep lines away from axis labels and margins.
 upper_chars=[c for c in p.chars if c['bottom']<140]
 upper=(min(r['y0'] for r in p.rects)-.5,max(c['y1'] for c in upper_chars)+.5)
 lower=(min(c['y0'] for c in p.curves),max(c['y1'] for c in p.curves))
 chars=p.chars;curve_signature=[(c['x0'],c['x1'],c['y0'],c['y1'],c.get('non_stroking_color')) for c in p.curves]
 char_signature=[(c['text'],c['fontname'],c['x0'],c['y0'],c['size']) for c in p.chars]
 assert len(curve_signature)==1056 and not any('ZapfDingbats' in c['fontname'] for c in chars)
 rect_signature=[(r['x0'],r['x1'],r['y0'],r['y1'],r.get('non_stroking_color')) for r in p.rects]
 unit_words=[g for g in p.extract_words() if g['text']=='(MB)']
 assert len(unit_words)==1,'Expected exactly one original MB axis unit'
 unit=unit_words[0]
 def outside_unit(g):
  return not (unit['x0']-.01<=g['x0'] and g['x1']<=unit['x1']+.01
              and unit['top']-.01<=g['top'] and g['bottom']<=unit['bottom']+.01)
 other_char_signature=[(g['text'],g['fontname'],g['x0'],g['y0'],g['size']) for g in chars if outside_unit(g)]
 original_character_text=''.join(g['text'] for g in chars)

# Edit the existing text-show literal using its existing font. No masking,
# replacement font, raster editing, coordinate change or statistics are used.
prepared_writer=PdfWriter(clone_from=str(BASE));page=prepared_writer.pages[0]
stream=page.get_contents();operations=stream.operations;unit_edits=[]
for index,(operands,operator) in enumerate(operations):
 if operator!=b'TJ':continue
 array=operands[0]
 label=''.join(str(v) for v in array if isinstance(v,TextStringObject))
 if label!='chr 17 physical position (MB)':continue
 for part_index,part in enumerate(array):
  if isinstance(part,TextStringObject) and '(MB)' in part:
   array[part_index]=TextStringObject(str(part).replace('(MB)','(Mb)'))
   unit_edits.append(index)
assert len(unit_edits)==1,'Expected one axis-unit text edit'
stream.operations=operations;page.replace_contents(stream)
original_operations=PdfReader(BASE).pages[0].get_contents().operations
assert len(original_operations)==len(operations)
assert [i for i,(a,b) in enumerate(zip(original_operations,operations)) if a!=b]==unit_edits

if args.check_only:
 print(json.dumps({'status':'PASS','scope':'base, coordinate calibration, vector primitives and in-memory axis-unit correction; no figure written',
 'axis_unit':'Mb','unit_text_edits':len(unit_edits),'base_sha256':sha(BASE),
 'core_mb':[45.414716,46.836265],'x_pdf_pt':xpos,'max_residual_pt':residual,'vector_circles':len(curve_signature)}));raise SystemExit(0)
OUT.mkdir(parents=True,exist_ok=True)
if (OUT/(STEM+'.pdf')).exists():raise SystemExit('Refusing to overwrite PDF')

def uncovered(lo,hi,excluded):
 cursor=lo
 for a,b in sorted(excluded):
  a=max(lo,a);b=min(hi,b)
  if b<=cursor or a>=hi:continue
  if a>cursor:yield cursor,a
  cursor=max(cursor,b)
 if cursor<hi:yield cursor,hi

buffer=io.BytesIO();c=canvas.Canvas(buffer,pagesize=(width,height))
c.setStrokeColorRGB(.5,.5,.5);c.setLineWidth(.6);c.setLineCap(0)
segments=[]
for x in xpos:
 for track,(lo,hi) in [('gene',upper),('SNP',lower)]:
  exclusions=[]
  if track=='gene':
   # Leave a small gap around every intersected transcript label; do not
   # paint white masks and do not shift any text or boundary coordinate.
   exclusions=[(g['y0']-.65,g['y1']+.65) for g in upper_chars if g['x0']-.7<x<g['x1']+.7]
  for y0,y1 in uncovered(lo,hi,exclusions):
   if y1-y0<.2:continue
   c.saveState();path=c.beginPath();path.rect(0,lo,width,hi-lo);c.clipPath(path,stroke=0,fill=0)
   c.setDash([2.5,2],(y0-lo)%4.5);c.line(x,y0,x,y1);c.restoreState()
   segments.append({'track':track,'x_pdf_pt':x,'y0':y0,'y1':y1})
# The accepted v2 continues boundaries through the inter-track blank interval.
for x in xpos:
 c.setDash([2.5,2],(lower[1]-lower[0])%4.5);c.line(x,lower[1],x,upper[0])
c.save()
# Underlay keeps all original rectangles, labels and 1056 circles on top.
page.merge_page(PdfReader(buffer).pages[0],over=False)
w=PdfWriter();w.add_page(page)
with (OUT/(STEM+'.pdf')).open('wb') as f:w.write(f)
with pdfplumber.open(OUT/(STEM+'.pdf')) as d:
 p=d.pages[0]
 assert (p.width,p.height)==(width,height)
 assert original_character_text.replace('(MB)','(Mb)')==''.join(g['text'] for g in p.chars)
 assert other_char_signature==[(g['text'],g['fontname'],g['x0'],g['y0'],g['size']) for g in p.chars if outside_unit(g)]
 assert any(g['text']=='(Mb)' for g in p.extract_words())
 assert curve_signature==[(g['x0'],g['x1'],g['y0'],g['y1'],g.get('non_stroking_color')) for g in p.curves]
 # The merge clipping path can be represented as an extra rectangle;
 # all original visible gene rectangles must still be present.
 assert all(g in [(v['x0'],v['x1'],v['y0'],v['y1'],v.get('non_stroking_color')) for v in p.rects] for g in rect_signature)

qc={'core_bp':[45414716,46836265],'core_mb':[45.414716,46.836265],'x_pdf_pt':xpos,
 'line_width_pt':.6,'dash_pt':[2.5,2],'connected_gap':True,'statistics_rerun':False,
 'base_sha256':sha(BASE),'axis_unit':'Mb','unit_text_edits':len(unit_edits),
 'only_text_change':'MB -> Mb','other_text_and_circles_preserved':True}
(OUT/'annotation_QC.json').write_text(json.dumps(qc,indent=2),encoding='utf-8')
if args.raster:
 import pypdfium2 as pdfium
 doc=pdfium.PdfDocument(str(OUT/(STEM+'.pdf')));p=doc[0]
 im=p.render(scale=4500/width).to_pil().convert('RGB')
 im.save(OUT/(STEM+'.tiff'),compression='tiff_lzw',dpi=(600,600))
 im.thumbnail((1600,1600));im.save(OUT/(STEM+'_preview.png'));p.close();doc.close()
print(json.dumps(qc,indent=2))
