"""Apply a narrow, context-checked workspace patch to a pinned Set_XC_Grid.c.

Keep the stock mathematical loops intact. Reuse only caller-private temporary
gradients. External derivative outputs (including stress) retain stock storage.
"""
from pathlib import Path
import sys


def prepare(source):
    def once(old,new):
        nonlocal source
        if source.count(old)!=1:
            raise ValueError('Set_XC_Grid source context differs: '+old[:70])
        source=source.replace(old,new)

    once('#include <omp.h>', '#include <omp.h>\n#include "cpu_xc_workspace.h"')
    once('  dDen_Grid_NULL_flag = 0;', '''  int cpu_xc_workspace = CpuXcWorkspaceMode();
  int cpu_xc_private = cpu_xc_workspace && XC_switch==4 &&
      XC_P_switch>=0 && XC_P_switch<=2 && !dDen_Grid && !dEXC_dGD;
  int cpu_xc_spin_alias = cpu_xc_private && cpu_xc_workspace==2;
  if(cpu_xc_spin_alias)CpuXcWorkspaceCheckInputs(Den0,Den1);
  dDen_Grid_NULL_flag = 0;''')
    once('''    for (k=0; k<=1; k++){
      dDen_Grid[k] = (double**)malloc(sizeof(double*)*3);''', '''    for (k=0; k<=1; k++){
      if(cpu_xc_spin_alias && k==1){dDen_Grid[k]=dDen_Grid[0];continue;}
      dDen_Grid[k] = (double**)malloc(sizeof(double*)*3);''')
    once('  if ( XC_switch==4 && dEXC_dGD==NULL){', '''  if(cpu_xc_private)dEXC_dGD=dDen_Grid;
  if ( XC_switch==4 && dEXC_dGD==NULL){''')
    once('''    for (k=0; k<=1; k++){
      for (i=0; i<3; i++){
        free(dDen_Grid[k][i]);''', '''    for (k=0; k<=1; k++){
      if(cpu_xc_spin_alias && k==1)continue;
      for (i=0; i<3; i++){
        free(dDen_Grid[k][i]);''')

    # The divergence phase no longer consumes density gradients. Its stock
    # low-density zero writes must not mutate the aliased derivative field
    # while other points read neighboring derivatives.
    marker='        calculate the second part of XC potential'
    pos=source.index(marker)
    before,after=source[:pos],source[pos:]
    old='''\t  dDen_Grid[0][0][MN] = 0.0;
\t  dDen_Grid[0][1][MN] = 0.0;
\t  dDen_Grid[0][2][MN] = 0.0;
\t  dDen_Grid[1][0][MN] = 0.0;
\t  dDen_Grid[1][1][MN] = 0.0;
\t  dDen_Grid[1][2][MN] = 0.0;'''
    if after.count(old)!=1:
        raise ValueError('unexpected late density-gradient writes')
    after=after.replace(old,'if(!cpu_xc_private){\n'+old+'\n}')
    source=before+after
    once('  dDen_Grid_calc_method = 2;', '''  dDen_Grid_calc_method = 2;
  if(cpu_xc_private)CpuXcWorkspaceReport(cpu_xc_workspace,My_NumGridD);''')
    return source


if __name__=='__main__':
    Path(sys.argv[2]).write_text(prepare(Path(sys.argv[1]).read_text()))
