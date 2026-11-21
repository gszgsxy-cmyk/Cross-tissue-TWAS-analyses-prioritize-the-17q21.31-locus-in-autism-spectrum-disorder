from pathlib import Path
import math, io, subprocess, shutil
import pandas as pd
from reportlab.pdfgen import canvas
from reportlab.lib.colors import HexColor
from pypdf import PdfReader, PdfWriter, Transformation
from PIL import Image

ROOT=Path('__WORKSPACE_ROOT__')
OUT=ROOT/'outputs/Supplementary_Figures_FINAL_WORKING_20261003'
P=OUT/'S2_panels'
SRC=ROOT/'outputs/fetalbrain_20260922/Final_QC'
d=pd.read_csv(SRC/'03_shared_testable_CpG_universe.tsv',sep='\t')
f=set(d.loc[d.Fetal_significant,'probeID']);h=set(d.loc[d.EUR_significant,'probeID'])
assert (len(d),len(f),len(h),len(f&h),len(f-h),len(h-f))==(6842,5,7,5,0,2)
s=pd.read_csv(SRC/'06_Fetal_vs_EUR_overlap_enrichment.tsv',sep='\t').iloc[0]
assert abs(s.One_sided_Fisher_P/1.68313936993289e-16-1)<1e-10
assert set(pd.read_csv(P/'S2A_plot_data.tsv',sep='\t').probeID)==f

c=canvas.Canvas(str(P/'S2B_Hannon_Hatton_Venn.pdf'),pagesize=(720,405))
c.setTitle('Hannon fetal-brain and Hatton EUR CpG overlap')
ink='#302B33';gray='#6C6B70';blue='#608AA6';sand='#C7B59B'
def t(x,y,z,size=10,font='Helvetica',col=ink,align='left'):
    c.setFillColor(HexColor(col));c.setFont(font,size)
    {'left':c.drawString,'center':c.drawCentredString,'right':c.drawRightString}[align](x,y,z)
def line(x,y,x2,y2,col='#DCDCDC',width=.6):
    c.setStrokeColor(HexColor(col));c.setLineWidth(width);c.line(x,y,x2,y2)
t(14,379,'B',17,'Helvetica-Bold')
t(42,380,'Cross-resource CpG overlap',12,'Helvetica-Bold')
t(42,360,'Shared testable background: 6,842 CpGs',10,col=gray)
line(480,65,480,382)
t(174,334,'Hannon fetal-brain',11,'Helvetica-Bold',align='center')
t(174,319,'Brain mQTL, n = 5',9,col=gray,align='center')
t(367,334,'Hatton EUR',11,'Helvetica-Bold',align='center')
t(367,319,'Whole-blood mQTL, n = 7',9,col=gray,align='center')
# All five significant Hannon CpGs belong to the significant Hatton set.
# Complete containment is required; a partial-overlap Venn would be incorrect.
xh,y,rh=238,216,88
rf=rh*math.sqrt(5/7);xf=xh-11
assert abs(xh-xf)+rf<rh
for x,r,fill,edge in [(xh,rh,'#E8DFD2',sand),(xf,rf,'#D1DEE6',blue)]:
    c.setFillColor(HexColor(fill));c.setStrokeColor(HexColor(edge));c.setLineWidth(1.1);c.circle(x,y,r,fill=1,stroke=1)
t(xf,y+5,'5',27,'Helvetica-Bold',align='center')
t(xf,y-14,'Shared',10,col=gray,align='center')
line(xh+rh-6,y,350,y,col=sand,width=1)
t(381,y+5,'2',25,'Helvetica-Bold',align='center')
t(381,y-14,'Hatton only',9,col=gray,align='center')
t(63,139,'Hannon only: 0',9,col=gray)
label='Overlap enrichment: P = 1.68 x 10'
lw=c.stringWidth(label,'Helvetica',10);sx=245-(lw+12)/2
t(sx,112,label,10);t(sx+lw+1,117,'-16',7)
t(245,94,'One-sided Fisher exact test',9,col=gray,align='center')
t(42,74,'Shared CpGs',9,'Helvetica-Bold')
shared=sorted(f)
t(42,59,'   '.join(shared[:3]),9)
t(42,45,'   '.join(shared[3:]),9)
t(500,380,'Fetal-brain sensitivity',12,'Helvetica-Bold')
t(500,360,'Hannon significant CpGs',10,col=gray)
c.setFillColor(HexColor('#F3ECE0'));c.rect(495,170,209,169,fill=1,stroke=0)
t(511,306,'5 / 5',25,'Helvetica-Bold')
t(511,284,'Within the 17q21.31 core interval',9.5)
c.setFillColor(HexColor(sand));c.rect(511,253,176,14,fill=1,stroke=0)
t(511,230,'MAPT-related: 2 CpGs',9)
t(511,210,'HEIDI unavailable: 5 CpGs',9)
t(511,188,'7,610 valid SMR tests',9)
t(500,146,'GRCh37/hg19 core interval',9,col=gray)
t(500,130,'chr17:43,492,082-44,913,631',8.5,col=gray)
t(500,101,'Circle areas proportional to counts.',8.5,col=gray)
t(500,85,'Hannon is fully contained in Hatton.',8.5,col=gray)
line(42,32,704,32)
t(42,18,'Original resource-specific FDR retained; Hatton n = 7 in the shared background (30 in the full resource).',8,col=gray)
c.save()
d.loc[d.Fetal_significant|d.EUR_significant].to_csv(P/'S2B_significant_CpG_membership.tsv',sep='\t',index=False)

# Stack the original-style forest and Venn panels, preserving vector content.
readerA=PdfReader(P/'S2A_Hannon_forest.pdf');a=readerA.pages[0]
readerB=PdfReader(P/'S2B_Hannon_Hatton_Venn.pdf');b=readerB.pages[0]
width=900;margin=18;usable=width-2*margin
sa=usable/float(a.mediabox.width);sb=usable/float(b.mediabox.width)
ha=float(a.mediabox.height)*sa;hb=float(b.mediabox.height)*sb
height=ha+hb+margin*2+14
writer=PdfWriter();page=writer.add_blank_page(width,height)
page.merge_transformed_page(a,Transformation().scale(sa).translate(margin,margin+hb+14))
page.merge_transformed_page(b,Transformation().scale(sb).translate(margin,margin))
writer.add_metadata({'/Title':'Supplementary Figure S2. Hannon fetal-brain sensitivity analysis and cross-resource comparison'})
pdf=OUT/'Supplementary_Figure_S2.pdf'
with pdf.open('wb') as stream:writer.write(stream)

pop=shutil.which('pdftoppm')
assert pop, 'Poppler pdftoppm is required'
subprocess.run([pop,'-singlefile','-r','120','-png',str(pdf),str(OUT/'Supplementary_Figure_S2_preview')],check=True)
subprocess.run([pop,'-singlefile','-r','300','-png',str(pdf),str(OUT/'Supplementary_Figure_S2')],check=True)
im=Image.open(OUT/'Supplementary_Figure_S2.png').convert('RGB')
im.save(OUT/'Supplementary_Figure_S2.tiff',compression='tiff_lzw',dpi=(300,300));im.close()
text='\n'.join(pg.extract_text() or '' for pg in PdfReader(pdf).pages)
assert all(cpg in text for cpg in f)
assert all(z in text for z in ['6,842','7,610','Hannon','Hatton','5 / 5'])
shutil.copy2(ROOT/'work/plot_S2_forest_20261003.R',P/'plot_S2_forest.R')
shutil.copy2(ROOT/'work/plot_S2_venn_assemble_20261003.py',P/'plot_S2_venn_assemble.py')
shutil.copy2(ROOT/'outputs/2026-09-19_GTEx_nonbrain_SMR_rebuild/07_scripts/figure4_style.R',P/'figure4_style.R')
print('S2 saved: five forest rows; nested Venn 0/5/2; shared universe 6842. Vector PDF and 300-dpi TIFF.')
