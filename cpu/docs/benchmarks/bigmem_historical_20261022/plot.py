from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator,FuncFormatter,NullFormatter
R=Path(__file__).resolve().parent;allrows=json.loads((R/'measurements.json').read_text());rows=[r for r in allrows if r['energy_gate_only'] and r['metadata_verified']]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.labelsize':10,'axes.titlesize':11,'axes.spines.top':False,'axes.spines.right':False,'axes.edgecolor':'#B5C2CF','text.color':'#1C2C3D','axes.labelcolor':'#425266','xtick.color':'#647184','ytick.color':'#647184','axes.linewidth':.7,'savefig.dpi':230,'pdf.fonttype':42,'svg.fonttype':'none'})
fig,axs=plt.subplots(1,2,figsize=(9.5,4.6),sharey=True);fig.patch.set_facecolor('#FAFBFD');fig.subplots_adjust(left=.08,right=.98,top=.73,bottom=.20,wspace=.18)
fig.text(.08,.94,'HOW SPEEDUP SCALES',fontsize=10,weight='bold',color='#637386');fig.text(.08,.87,'Larger systems, larger gains',fontsize=21,weight='bold')
fig.text(.08,.80,'Historical build · bigmem001–003 · 128 MPI × 1 · complete NSOC + one-step SOC',fontsize=9,color='#637386')
C=['#087F8C','#DD7D63','#536DBA'];selection={'D3':['GaTe','MoTe2','WSe2'],'C3':['Bi2Se3','PbI2','HfSe2']};pretty={'MoTe2':'MoTe₂','WSe2':'WSe₂','Bi2Se3':'Bi₂Se₃','PbI2':'PbI₂','HfSe2':'HfSe₂'}
for ax,group in zip(axs,['D3','C3']):
 rr=[r for r in rows if r['group']==group];materials=sorted(set(r['material']for r in rr))
 for mat in materials:
  if mat in selection[group]:continue
  ss=sorted([r for r in rr if r['material']==mat],key=lambda r:r['atoms']);ax.plot([r['atoms']for r in ss],[r['speedup']for r in ss],color='#C7D0D9',lw=.9,alpha=.7,marker='.',ms=3,zorder=1)
 for j,mat in enumerate(selection[group]):
  ss=sorted([r for r in rr if r['material']==mat],key=lambda r:r['atoms']);x=[r['atoms']for r in ss];y=[r['speedup']for r in ss];ax.plot(x,y,color=C[j],lw=1.8,marker=['o','s','^'][j],ms=5,mec='white',mew=.6,label=pretty.get(mat,mat),zorder=3)
  if mat in ['GaTe','MoTe2','Bi2Se3']:
   ax.annotate(f'{y[-1]:.2f}×',(x[-1],y[-1]),xytext=(7,5 if mat!='MoTe2' else -9),textcoords='offset points',fontsize=9,weight='bold',color=C[j])
 ax.set_facecolor('white');ax.set_title(group+f'  /  {len(materials)} materials',loc='left',pad=10,weight='bold');ax.axhline(1,color='#8B99A8',ls=(0,(3,3)),lw=.8);ax.set_xscale('log');ax.set_xlim(180,3000);ax.set_ylim(.9,4.95);ax.set_yticks([1,2,3,4]);ax.xaxis.set_major_locator(FixedLocator([200,400,800,1600]));ax.xaxis.set_minor_formatter(NullFormatter());ax.xaxis.set_major_formatter(FuncFormatter(lambda x,p:f'{x:,.0f}'));ax.grid(axis='y',color='#E4EAF0',lw=.6);ax.legend(loc='upper left',frameon=False,fontsize=8,handlelength=1.8);ax.set_xlabel('Number of atoms',labelpad=7)
axs[0].set_ylabel('Original / optimized total time',labelpad=8)
fig.text(.08,.078,'Colored curves: named examples. Gray curves: all other energy-qualified materials. No fit or extrapolation.',fontsize=8,color='#637386')
fig.text(.08,.035,'122 / 127 cases pass the historical 0.1 meV/atom energy gate; 5 WTe₂ cases excluded. Not a current-version or force-equivalence claim.',fontsize=7.3,color='#637386')
for ext in ['png','pdf','svg']:fig.savefig(R/f'historical_speedup_scaling.{ext}',facecolor=fig.get_facecolor())
p=R/'historical_speedup_scaling.svg';p.write_text('\n'.join(l.rstrip() for l in p.read_text().splitlines())+'\n')
print('Plotted',len(rows),'qualified cases in',len({r['material'] for r in rows}),'materials; maximum',max(r['speedup']for r in rows))
