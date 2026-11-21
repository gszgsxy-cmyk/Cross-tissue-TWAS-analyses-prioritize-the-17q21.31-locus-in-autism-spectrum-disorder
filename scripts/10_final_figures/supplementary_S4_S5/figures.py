import pathlib,json
from reportlab.pdfgen import canvas
from reportlab.lib.colors import HexColor
import pypdfium2 as pdfium
import argparse
R=pathlib.Path(__file__).resolve().parents[3]
parser=argparse.ArgumentParser(description='Draw saved S4/S5 counts; no enrichment/FDR analysis.')
parser.add_argument('--output-dir',type=pathlib.Path,default=pathlib.Path(__file__).resolve().parent/'output')
parser.add_argument('--check-only',action='store_true')
args=parser.parse_args();W=R/'data/external_evidence';O=args.output_dir
d=json.loads((W/'figure_summary.json').read_text(encoding='utf8'))
assert [r['Overlap number'] for r in d['Overlap_summary']]==[10,67,22]
assert [r['N_genes'] for r in d['Support_distribution']]==[70,67,13,2]
if args.check_only:
 print(json.dumps({'status':'PASS','scope':'saved S4/S5 inputs only','counts':d}));raise SystemExit(0)
if O.exists() and any(O.glob('*.pdf')):raise SystemExit('Refusing to overwrite existing PDF outputs')
O.mkdir(parents=True,exist_ok=True)

def figure(name,title,labels,counts,colors,foot):
 p=O/(name+'.pdf');c=canvas.Canvas(str(p),pagesize=(432,324));c.setTitle(title);c.setAuthor('ASD TWAS project')
 c.setFillColor(HexColor('#151515'));c.setFont('Helvetica-Bold',12);c.drawString(54,298,title)
 c.setFont('Helvetica',9);c.drawString(54,280,'JTI / FUSION significant union (n = 152)')
 left,bottom,width,height=54,73,355,185; ymax=80
 c.setLineWidth(.6)
 for y in range(0,81,20):
  py=bottom+height*y/ymax;c.setStrokeColor(HexColor('#DFE3E6'));c.line(left,py,left+width,py);c.setFillColor(HexColor('#333333'));c.setFont('Helvetica',9);c.drawRightString(left-8,py-3,str(y))
 c.setStrokeColor(HexColor('#222222'));c.line(left,bottom,left+width,bottom);c.line(left,bottom,left,bottom+height)
 c.saveState();c.translate(17,bottom+height/2);c.rotate(90);c.setFont('Helvetica',10);c.drawCentredString(0,0,'Number of genes');c.restoreState()
 gap=width/len(labels);bw=gap*.52
 for i,(label,n,color) in enumerate(zip(labels,counts,colors)):
  x=left+gap*(i+.5);bar=height*n/ymax;c.setFillColor(HexColor(color));c.rect(x-bw/2,bottom,bw,bar,fill=1,stroke=0)
  c.setFillColor(HexColor('#222222'));c.setFont('Helvetica-Bold',10);c.drawCentredString(x,bottom+bar+7,str(n));c.setFont('Helvetica',10);c.drawCentredString(x,bottom-16,label)
 c.setFont('Helvetica',8)
 for i,line in enumerate(foot):c.drawString(54,34-i*12,line)
 c.showPage();c.save()
 pdf=pdfium.PdfDocument(str(p));page=pdf[0];im=page.render(scale=300/72).to_pil();im.save(O/(name+'.png'),dpi=(300,300));page.close();pdf.close()
figure('Figure_external_evidence_overlap_barplot','Overlap with external ASD evidence',['Huang 2021','Matoba 2020','Gandal 2022'],[r['Overlap number'] for r in d['Overlap_summary']],['#486A89','#78917B','#9D7C62'],['Non-exclusive counts from three distinct evidence sources.','Descriptive overlap; bars do not represent independent replication.'])
figure('Figure_gene_support_distribution','Number of external literature sources',['0 sources','1 source','2 sources','3 sources'],[r['N_genes'] for r in d['Support_distribution']],['#B8C0C7','#8FA3B3','#607F98','#315877'],['Sources counted: Huang 2021, Matoba 2020 and Gandal 2022 only.','SFARI, Fu and Walker are excluded from this figure.'])
print('Two vector PDFs and two 300-dpi PNGs generated.')
