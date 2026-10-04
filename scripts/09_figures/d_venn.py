from pathlib import Path
import math,json,subprocess,shutil
import pandas as pd
from reportlab.pdfgen import canvas
from reportlab.lib.colors import HexColor
from PIL import Image
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'outputs/Figure4D_Venn_20260928';OUT.mkdir(exist_ok=True)
E=ROOT/'outputs/2026-09-22_Hatton_EUR_mQTL_rebuild/Enrichment'
d=pd.read_csv(E/'03_shared_testable_universe.tsv',sep='\t')
b=set(d.loc[d.Brain_significant,'probeID']);h=set(d.loc[d.EUR_significant,'probeID'])
assert (len(d),len(b),len(h),len(b&h))==(61780,16,19,4)
s=pd.read_csv(E/'06_enrichment_statistics.tsv',sep='\t')
reg=pd.read_csv(E/'02_region_2x2.tsv',sep='\t')
assert list(reg.Inside_core)==[21,82] and list(reg.Outside_core)==[9,137198]
pdf=OUT/'Figure_4D_CpG_overlap.pdf'
c=canvas.Canvas(str(pdf),pagesize=(720,405))
c.setTitle('Figure 4D: CpG overlap and regional distribution')
ink='#302B33';gray='#6C6B70';blue='#608AA6';sand='#C7B59B'
def t(x,y,z,size=10,font='Helvetica',col=ink,align='left'):
 c.setFillColor(HexColor(col));c.setFont(font,size)
 (c.drawCentredString if align=='center' else c.drawString)(x,y,z)
def line(x,y,x2,y2,col='#DCDCDC'):
 c.setStrokeColor(HexColor(col));c.setLineWidth(.6);c.line(x,y,x2,y2)
t(14,379,'D',17,'Helvetica-Bold')
t(42,380,'Cross-resource CpG overlap',12,'Helvetica-Bold')
t(42,360,'Shared testable background: 61,780 CpGs',10,col=gray)
t(500,380,'17q21.31 distribution',12,'Helvetica-Bold')
t(500,360,'Hatton EUR significant CpGs',10,col=gray)
line(480,65,480,382)
# Circle areas and intersection area are proportional to 16, 19 and 4.
r1=76.;scale=math.pi*r1*r1/len(b);r2=math.sqrt(scale*len(h)/math.pi)
def intersection(dist):
 a=math.acos((dist*dist+r1*r1-r2*r2)/(2*dist*r1))
 z=math.acos((dist*dist+r2*r2-r1*r1)/(2*dist*r2))
 return r1*r1*a+r2*r2*z-.5*math.sqrt((-dist+r1+r2)*(dist+r1-r2)*(dist-r1+r2)*(dist+r1+r2))
low=abs(r2-r1)+1e-6;high=r1+r2-1e-6
for _ in range(90):
 mid=(low+high)/2
 if intersection(mid)>scale*len(b&h):low=mid
 else:high=mid
dist=(low+high)/2;x1=(450-(r1+dist+r2))/2+r1;x2=x1+dist;y=230
for x,r,color in [(x1,r1,'#D1DEE6'),(x2,r2,'#E8DFD2')]:
 c.setFillColor(HexColor(color));c.circle(x,y,r,fill=1,stroke=0)
c.saveState();clip=c.beginPath();clip.circle(x1,y,r1);c.clipPath(clip,stroke=0,fill=0)
c.setFillColor(HexColor('#B9C3C7'));c.circle(x2,y,r2,fill=1,stroke=0);c.restoreState()
for x,r,color in [(x1,r1,blue),(x2,r2,sand)]:
 c.setStrokeColor(HexColor(color));c.setLineWidth(1.1);c.circle(x,y,r,fill=0,stroke=1)
cross=(x1+r1+x2-r2)/2
t(x1-25,y-7,'12',26,'Helvetica-Bold',align='center')
t(cross,y-7,'4',26,'Helvetica-Bold',align='center')
t(x2+26,y-7,'15',26,'Helvetica-Bold',align='center')
t(x1-12,334,'Brain-mMeta',11,'Helvetica-Bold',align='center')
t(x1-12,319,'Brain mQTL | n = 16',9,col=gray,align='center')
t(x2+15,334,'Hatton EUR',11,'Helvetica-Bold',align='center')
t(x2+15,319,'Whole-blood mQTL | n = 19',9,col=gray,align='center')
label='Overlap enrichment: P = 1.16 x 10'
lw=c.stringWidth(label,'Helvetica',10);sx=225-(lw+12)/2
t(sx,128,label,10);t(sx+lw+1,133,'-11',7)
t(225,108,'One-sided Fisher exact test',9,col=gray,align='center')
t(42,78,'Shared CpGs',9,'Helvetica-Bold')
t(42,61,'cg00846647   cg01798261   cg07368061   cg15295732',9)
# Regional summary uses all 30 significant Hatton probes, not the shared universe.
c.setFillColor(HexColor('#F3ECE0'));c.rect(495,170,209,169,fill=1,stroke=0)
t(511,306,'21 / 30',25,'Helvetica-Bold')
t(511,284,'70% within the core interval',10)
c.setFillColor(HexColor('#E0DEDA'));c.rect(511,252,176,14,fill=1,stroke=0)
c.setFillColor(HexColor(sand));c.rect(511,252,176*.7,14,fill=1,stroke=0)
t(511,232,'Within: 21',9);t(616,232,'Outside: 9',9)
t(511,206,'Regional enrichment',9,'Helvetica-Bold')
t(511,187,'P = 3.81 x 10',10);t(573,192,'-60',7)
t(500,143,'Full Hatton testable background:',9,col=gray)
t(500,128,'103 / 137,310 within the interval',9)
t(500,111,'Two-sided Fisher exact test',9,col=gray)
t(500,82,'GRCh37/hg19 core interval',8.5,col=gray)
t(500,67,'chr17:43,492,082-44,913,631',8.5,col=gray)
line(42,45,704,45)
t(42,28,'Significance uses the original resource-specific results; FDR was not recalculated in the shared background.',8,col=gray)
c.save()
pop=Path('__USER_HOME__/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/Library/bin/pdftoppm.exe')
subprocess.run([str(pop),'-singlefile','-r','600','-png',str(pdf),str(OUT/'render600')],check=True)
im=Image.open(OUT/'render600.png').convert('RGB');im.save(OUT/'Figure_4D_CpG_overlap.tiff',compression='tiff_lzw',dpi=(600,600))
im.thumbnail((1600,1000));im.save(OUT/'Figure_4D_CpG_overlap_preview.png')
(OUT/'render600.png').unlink()
d.loc[d.Brain_significant|d.EUR_significant].to_csv(OUT/'Venn_significant_probe_membership.tsv',sep='\t',index=False)
shutil.copy2(__file__,OUT/'d_venn.py')
(OUT/'README.md').write_text('Figure 4D. Vector PDF and 600 dpi TIFF.\n\nVenn circle and overlap areas are proportional to counts (16,19,4). Neither significant: 61,749. Original significance sets retained. Regional background is the full Hatton dataset. Fisher P values are nominal.\n\nSources: '+str(E)+'\n',encoding='utf-8')
print('Verified: N=61780; Venn=12/4/15; region=21/30. PDF and TIFF saved.')
