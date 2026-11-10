"""Reproduce figures from the accompanying sanitized measurement summaries."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator,FuncFormatter,NullFormatter
R=Path(__file__).resolve().parent
S=json.loads((R/'summary.json').read_text());P=json.loads((R/'paired.json').read_text())
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.labelsize':9,'axes.titlesize':11,'axes.spines.top':False,'axes.spines.right':False,'axes.edgecolor':'#C6CED6','axes.labelcolor':'#374455','xtick.color':'#647184','ytick.color':'#647184','text.color':'#182637','axes.linewidth':.7,'xtick.labelsize':8,'ytick.labelsize':8,'savefig.dpi':230,'pdf.fonttype':42,'svg.fonttype':'none'})
C={'stock':'#7F8B99','automatic':'#087F8C','explicit_subgroup':'#DF7C62'}
fig,axs=plt.subplots(2,2,figsize=(9,6.65),sharex=True)
fig.subplots_adjust(left=.09,right=.97,top=.75,bottom=.11,hspace=.19,wspace=.28)
fig.patch.set_facecolor('#FAFBFD')
fig.text(.09,.953,'OPENMX  /  CPU PERFORMANCE',size=9,weight='bold',color='#637386')
fig.text(.09,.903,'Less compute. Less memory.',size=21,weight='bold')
fig.text(.09,.86,'Full NSOC + one-step SOC · single-node · 64 MPI × 1',size=9,color='#637386')
fig.text(.66,.927,f'{max(p["speedup"] for p in P):.2f}×',size=24,weight='bold',color=C['automatic'])
fig.text(.66,.885,'largest measured speedup',size=8,color='#637386')
fig.text(.85,.927,f'{max(p["memory_reduction_percent"] for p in P):.0f}%',size=24,weight='bold',color=C['explicit_subgroup'])
fig.text(.85,.885,'less peak memory*',size=8,color='#637386')
for col,(mat,name,group,sub,node) in enumerate([('MoTe2','MoTe₂','D3','C3','037'),('HgI2','HgI₂','D2','C2','039')]):
 for profile in ['stock','automatic','explicit_subgroup']:
  rr=sorted([r for r in S if r['material']==mat and r['profile']==profile],key=lambda r:r['atoms'])
  if not rr:continue
  x=np.array([r['atoms'] for r in rr]);label={'stock':'Original 4.0','automatic':group+' · automatic','explicit_subgroup':sub+' · subgroup'}[profile]
  for row,metric in enumerate(['core_hours','peak_GiB']):
   ax=axs[row,col];y=np.array([r[metric] for r in rr]);lo=np.array([r[metric+'_min'] for r in rr]);hi=np.array([r[metric+'_max'] for r in rr])
   ax.errorbar(x,y,yerr=[y-lo,hi-y],color=C[profile],label=label,marker={'stock':'s','automatic':'o','explicit_subgroup':'^'}[profile],lw=1.65,ms=5,mec='white',mew=.8,capsize=2,zorder=3)
 axs[0,col].set_title(name+'  /  node'+node,loc='left',pad=12,weight='bold');axs[0,col].legend(frameon=False,fontsize=8,loc='upper left',labelspacing=.4,handlelength=1.6)
 for row in [0,1]:
  ax=axs[row,col];ax.set_facecolor('white');ax.set_xscale('log');ax.set_xlim(180,2500);ax.xaxis.set_major_locator(FixedLocator([200,400,800,1600]));ax.xaxis.set_minor_formatter(NullFormatter());ax.xaxis.set_major_formatter(FuncFormatter(lambda x,p:f'{x:,.0f}'));ax.grid(axis='y',color='#E5EAF0',lw=.6);ax.set_axisbelow(True)
 axs[0,col].set_yscale('log');axs[0,col].set_ylim(1,48);axs[0,col].yaxis.set_major_locator(FixedLocator([1,2,4,8,16,32]));axs[0,col].yaxis.set_major_formatter(FuncFormatter(lambda x,p:f'{x:g}'));axs[0,col].yaxis.set_minor_formatter(NullFormatter())
 axs[1,col].set_ylim(0,210);axs[1,col].set_yticks([0,50,100,150,200]);axs[1,col].set_xlabel('Number of atoms')
axs[0,0].set_ylabel('Total cost (core·h)');axs[1,0].set_ylabel('Workflow peak memory (GiB)')
fig.text(.09,.035,'150 Ry · criterion 10⁻⁷ · History 30 / StartPulay 15 · completed, uninterrupted runs only',size=8,color='#637386')
fig.text(.09,.012,'*Maxima across four matched default-profile cases; not necessarily the same case. Lines join measurements; no fitted scaling law.',size=7,color='#637386')
for ext in ['png','pdf','svg']:fig.savefig(R/f'scaling.{ext}',facecolor=fig.get_facecolor())
plt.close(fig)
fig,axs=plt.subplots(1,2,figsize=(9,3.65),sharey=True)
fig.subplots_adjust(left=.25,right=.94,top=.76,bottom=.16,wspace=.3)
fig.patch.set_facecolor('#FAFBFD');fig.text(.075,.93,'Matched results, automatic symmetry',size=16,weight='bold');fig.text(.075,.86,'Original 4.0 → current CPU · same node and physical inputs · NSOC + one-step SOC',size=8,color='#637386')
y=np.arange(len(P))[::-1];labels=[f'{r["material"].replace("MoTe2","MoTe₂").replace("HgI2","HgI₂")}  {r["angle_deg"]:.2f}°\n{r["atoms"]:,} atoms · {r["group"]}' for r in P]
for ax,metric,title,color,limit,ticks in [(axs[0],'speedup','Workflow speedup ↑',C['automatic'],3.7,[1,2,3]),(axs[1],'memory_reduction_percent','Peak memory reduction ↑',C['explicit_subgroup'],85,[0,25,50,75])]:
 ax.set_facecolor('white');baseline=1 if metric=='speedup' else 0
 for pos,r in zip(y,P):
  val=r[metric];ax.plot([baseline,val],[pos,pos],lw=4,color=color,alpha=.22,solid_capstyle='round');ax.scatter([val],[pos],s=65,color=color,edgecolor='white',lw=1,zorder=3);ax.text(val+.09 if metric=='speedup' else val+2,pos,f'{val:.2f}×' if metric=='speedup' else f'{val:.1f}%',va='center',fontsize=9,weight='bold',color=color)
 ax.axvline(baseline,color='#AEBAC6',lw=.8,ls='--');ax.set_xlim(.85 if metric=='speedup' else -3,limit);ax.set_xticks(ticks);ax.grid(axis='x',color='#E5EAF0',lw=.6);ax.set_title(title,loc='left',size=10,pad=10);ax.set_ylim(-.5,3.5);ax.spines['left'].set_visible(False);ax.tick_params(axis='y',length=0)
axs[0].set_yticks(y,labels);axs[1].tick_params(labelleft=False)
fig.text(.075,.045,'Means of 1–2 runs per configuration. Symmetry projection changes the discretized operator; numerical differences are reported separately.',fontsize=7,color='#637386')
for ext in ['png','pdf','svg']:fig.savefig(R/f'paired_gains.{ext}',facecolor=fig.get_facecolor())
plt.close(fig)
print('Wrote scaling and matched-gain figures (PNG / SVG / PDF).')

# Keep generated SVG paths readable without trailing whitespace.
for path in R.glob('*.svg'):
    path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines()) + '\n')
