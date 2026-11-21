"""Candidate S2 vector assembly. Author adoption is not yet confirmed."""
from pathlib import Path
from copy import deepcopy
import argparse,io,json,csv
from pypdf import PdfReader,PdfWriter,Transformation
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output-dir',type=Path,default=Path(__file__).resolve().parent/'output')
parser.add_argument('--bold-font',type=Path)
parser.add_argument('--check-only',action='store_true')
args=parser.parse_args();PAN=Path(__file__).resolve().parent/'data';OUT=args.output_dir
with (PAN/'S2A_plot_data.tsv').open(encoding='utf-8-sig',newline='') as f:rows=list(csv.DictReader(f,delimiter='\t'))
assert len(rows)==5 and len({r['probeID'] for r in rows})==5
assert all(r['p_HEIDI']=='NA' and r['nsnp_HEIDI']=='NA' for r in rows)
assert len(PdfReader(PAN/'S2A_Hannon_forest.pdf').pages)==1
assert len(PdfReader(PAN/'S2B_Hannon_Hatton_Venn.pdf').pages)==1
if args.check_only:
 print(json.dumps({'status':'PASS','scope':'candidate archived-panel inputs only','CpGs':5,'adoption':'unconfirmed'}));raise SystemExit(0)
OUT.mkdir(parents=True,exist_ok=True)
if (OUT/'Supplementary_Figure_S2_Layout_v1.pdf').exists():raise SystemExit('Refusing to overwrite candidate PDF')
WIDTH=900.; MARGIN=18.; BODY_X=54.; BODY_WIDTH=WIDTH-MARGIN-BODY_X
A_CROP=(32.,36.,776.,166.)  # table + forest + ticks + axis title; no source heading
B_CROP=(42.,12.,708.,371.)  # all B results and footnote; no source headings
SCALE_A=BODY_WIDTH/(A_CROP[2]-A_CROP[0])
SCALE_B=BODY_WIDTH/(B_CROP[2]-B_CROP[0])
A_HEIGHT=(A_CROP[3]-A_CROP[1])*SCALE_A
B_HEIGHT=(B_CROP[3]-B_CROP[1])*SCALE_B
GLOBAL_BAND=42.;HEADER_BAND=36.;PANEL_GAP=24.
HEIGHT=2*MARGIN+GLOBAL_BAND+2*HEADER_BAND+A_HEIGHT+PANEL_GAP+B_HEIGHT
CURSOR=HEIGHT-MARGIN-GLOBAL_BAND
A_LABEL_Y=CURSOR-22.;A_Y=CURSOR-HEADER_BAND-A_HEIGHT
CURSOR=A_Y-PANEL_GAP
B_LABEL_Y=CURSOR-22.;B_Y=MARGIN
FONT=args.bold_font
if FONT is None or not FONT.is_file():raise SystemExit('Supply --bold-font /path/to/Arial-Bold.ttf to retain the approved heading font')
pdfmetrics.registerFont(TTFont('S2ArialBold',str(FONT)))
writer=PdfWriter();page=writer.add_blank_page(WIDTH,HEIGHT)
def place(source,crop,scale,x,y):
 p=deepcopy(PdfReader(source).pages[0]);p.cropbox.lower_left=crop[:2];p.cropbox.upper_right=crop[2:]
 transform=Transformation().translate(-crop[0],-crop[1]).scale(scale).translate(x,y)
 page.merge_transformed_page(p,transform,over=True,expand=False)
place(PAN/'S2A_Hannon_forest.pdf',A_CROP,SCALE_A,BODY_X,A_Y)
place(PAN/'S2B_Hannon_Hatton_Venn.pdf',B_CROP,SCALE_B,BODY_X,B_Y)
overlay=io.BytesIO();c=canvas.Canvas(overlay,pagesize=(WIDTH,HEIGHT))
c.setFillColorRGB(0,0,0)
c.setFont('S2ArialBold',24);c.drawString(MARGIN,HEIGHT-MARGIN-22,'Supplementary Figure S2')
for letter,title,y in [('A','Hannon fetal-brain mQTL-SMR',A_LABEL_Y),('B','Cross-resource CpG overlap',B_LABEL_Y)]:
 c.setFont('S2ArialBold',22);c.drawString(MARGIN,y,letter)
 c.setFont('S2ArialBold',20);c.drawString(BODY_X,y,title)
c.setFont('S2ArialBold',20)
c.drawString(BODY_X+(500-B_CROP[0])*SCALE_B,B_LABEL_Y,'Fetal-brain sensitivity')
c.save();page.merge_page(PdfReader(overlay).pages[0])
writer.add_metadata({'/Title':'Supplementary Figure S2: aligned headings and compact layout','/Subject':'Original five-CpG forest and nested overlap vectors preserved; layout revision only'})
pdf=OUT/'Supplementary_Figure_S2_Layout_v1.pdf'
with pdf.open('wb') as stream:writer.write(stream)
