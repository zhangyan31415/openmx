#!/usr/bin/env python3
"""Plot saved same-node reference/latest measurements; never launch OpenMX."""
from pathlib import Path
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator

HERE = Path(__file__).resolve().parent


def main():
    with (HERE/'measurements.csv').open(newline='') as stream:
        rows = list(csv.DictReader(stream))
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 5.3))
    styles = [('reference', 'CPU 4.0 source-snapshot reference', '#737c89', 's'),
              ('latest', 'Latest source + scf.Symmetry on', '#2563a6', 'o')]
    for ax, key, title, unit in [
            (axes[0], 'total_seconds', 'Complete NSOC + one-step SOC', 'Native time (s)'),
            (axes[1], 'workflow_peak_GiB', 'Full workflow physical memory peak', 'Task cgroup peak (GiB)')]:
        for version, label, color, marker in styles:
            points = sorted((int(r['atoms']), float(r[key])) for r in rows if r['configuration']==version)
            assert len(points)==2
            ax.plot([p[0] for p in points], [p[1] for p in points], color=color,
                    marker=marker, ms=6, lw=1.7, label=label)
            for x, y in points:
                ax.annotate(f'{y:.2f}', (x,y), xytext=(0,10 if version=='reference' else -18),
                            textcoords='offset points', ha='center', fontsize=10, color=color)
        ax.set_xscale('log')
        ax.set_xlim(195,415)
        ax.xaxis.set_major_locator(FixedLocator([222,366]))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x,_:'222\n9.43 degrees' if x==222 else '366\n7.34 degrees'))
        ax.xaxis.set_minor_locator(NullLocator())
        ax.set_ylim(0,max(float(r[key]) for r in rows)*1.25)
        ax.set_xlabel('Number of atoms (log scale); twist angle below')
        ax.set_ylabel(unit)
        ax.set_title(title, fontsize=11.5, fontweight='bold')
        ax.grid(axis='y',alpha=.2)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.suptitle(r'Fresh MoTe$_2$ benchmark: reference source vs latest source',
                 fontsize=15.3,fontweight='bold',y=.98)
    fig.text(.5,.905,'Same node040 | 64 MPI x 1 | same oneAPI toolchain and CORE-AVX2 target',
             ha='center',fontsize=10.5,color='#444444')
    fig.legend(handles,labels,ncol=2,loc='upper center',bbox_to_anchor=(.5,.867),frameon=False,fontsize=10)
    fig.text(.075,.026,
             '150 Ry | 300 K | SCF 1e-7 | History30 | StartPulay15. Each SOC reads its own converged NSOC.\n'
             'One run per version and size; lines connect measurements, not a fitted scaling law.\n'
             'Symmetry constraints change the discretized operator/SCF trajectory; physical differences are reported separately.',
             fontsize=8.9,color='#555555',linespacing=1.5)
    fig.subplots_adjust(top=.75,bottom=.28,left=.075,right=.97,wspace=.29)
    for suffix in ('png','svg'):
        output=HERE/('latest_vs_reference.'+suffix)
        fig.savefig(output,dpi=180,facecolor='white')
        if suffix=='svg':
            output.write_text('\n'.join(line.rstrip() for line in output.read_text().splitlines())+'\n')
    plt.close(fig)


if __name__=='__main__':
    main()
