#!/usr/bin/env python3
"""One-time source recovery recipe; imported sources build without result trees."""
from pathlib import Path
import hashlib,json,re
HERE=Path(__file__).resolve().parent
REPO=HERE.parents[2]
RESULTS=REPO.parent/'results'
DST=HERE/'source'
FUSED=RESULTS/'openmx_fused_return_20261025_190500'
P=RESULTS/'openmx_time_followup_20261025/nsoc_elpa_omp/source'
COMPACT=RESULTS/'openmx_nsoc_compact_734_20261024_232200'
LIFE=RESULTS/'openmx_nsoc_soc_memory_port_20261024_190300/lifetime/source'
PAIRED=RESULTS/'openmx_equivalent_blocks_20261022_065042/variants/paired/source'
SHARED=RESULTS/'openmx_nsoc_symmetry_memory_20261024_143800/source'
HGRID=RESULTS/'openmx_multinode_efficiency_20261025/realspace_blas3/native_h'
NU=RESULTS/'openmx_multinode_734_20261025_010500'
SPARSE=RESULTS/'openmx_multinode_solver_20261025_074000/sparse_transform/source'
UNITS={
'DFT':P/'DFT.direct.c','Cluster_DFT_Col':P/'Cluster_DFT_Col.direct.c',
'cpu_nsoc_sym':P/'cpu_nsoc_sym.direct.c','cpu_realspace_d3':P/'cpu_realspace_d3.panel.c',
'cpu_nsoc_owner':RESULTS/'openmx_nsoc_history_owner_20261024_205400/source/cpu_nsoc_owner.c',
'cpu_nsoc_owner_openmx':P/'cpu_nsoc_owner_openmx.c','cpu_nsoc_rep_tail':P/'cpu_nsoc_rep_tail.c',
'rep_panel':P/'rep_panel.c','selected_left':FUSED/'selected_left.c',
'cpu_nsoc_sparse':PAIRED/'cpu_nsoc_sparse.c',
'cpu_nsoc_imag':COMPACT/'combined/source/cpu_nsoc_imag.c',
'cpu_symhistory_diis':COMPACT/'history/diis_kernel/cpu_symhistory_diis.c',
'truncation':COMPACT/'history/sym_source/truncation.symhistory.c',
'Free_Arrays':COMPACT/'history/sym_source/Free_Arrays.symhistory.c',
'DIIS_Mixing_DM':COMPACT/'history/sym_source/DIIS_Mixing_DM.symhistory.c',
'Mixing_H':COMPACT/'combined/source/Mixing_H.compact.c',
'Hamiltonian_Cluster_Hs':COMPACT/'combined/source/Hamiltonian_Cluster_Hs.compact.c',
'Mixing_DM':PAIRED/'Mixing_DM.nsym.c','SCF2File':PAIRED/'SCF2File.realspace.c',
'Force':LIFE/'Force.fields.c','Set_Density_Grid':LIFE/'Set_Density_Grid.fields.c',
'Set_Orbitals_Grid':LIFE/'Set_Orbitals_Grid.numa.c','Set_dOrbitals_Grid':LIFE/'Set_dOrbitals_Grid.numa.c',
'Set_Aden_Grid':LIFE/'Set_Aden_Grid.fields.c','Set_ProExpn_VNA':LIFE/'Set_ProExpn_VNA.vna_shared.c',
'Set_OLP_Kin':SHARED/'Set_OLP_Kin.tight.c',
'Set_Hamiltonian':HGRID/'Set_Hamiltonian.h_blas3.c','cpu_h_grid_blas3':HGRID/'cpu_h_grid_blas3.c',
'cpu_grid_shared':NU/'numa/source/cpu_grid_shared.c',
'cpu_grid_fields_shared':NU/'fields/nsoc/cpu_grid_fields_shared.c',
'cpu_vna_shared':RESULTS/'openmx_multinode_solver_20261025_074000/vna_raw0/source/cpu_vna_shared.c'}
build=json.loads((FUSED/'BUILD_SELECTED.json').read_text())
flags=next(x['command']for x in build['commands']if x['label']=='compile_selected_solver')
search=[Path(x[2:])for x in flags if x.startswith('-I')]
search=[COMPACT/'history/sym_source',COMPACT/'history/diis_kernel',NU/'fields/nsoc',NU/'numa/source',HGRID,SPARSE]+search
records={}
def sha(data):return hashlib.sha256(data).hexdigest()
def put(src,text=None):
    data=src.read_bytes() if text is None else text.encode()
    path=DST/src.name
    if path.name in records:
        if path.read_bytes()!=data:
            if src.name in ('cpu_nsoc_sym.h','cpu_vna_shared.h') and all(line.strip() in re.sub(r'/\*.*?\*/','',path.read_text(),flags=re.S) for line in re.sub(r'/\*.*?\*/','',src.read_text(),flags=re.S).splitlines() if line.strip()):
                records[path.name].setdefault('compatible_older_origins',[]).append(str(src))
                return
            raise RuntimeError(f'conflicting source {src}')
        return
    path.write_bytes(data)
    records[path.name]={'origin':str(src),'origin_sha256':sha(src.read_bytes()),'sha256':sha(data),'rewritten':data!=src.read_bytes()}
    # Recursively vendor only project-specific includes; common stock headers
    # and external library headers come from standalone vendor/library roots.
    for name in re.findall(r'^\s*#\s*include\s*"([^"]+)"',data.decode(),re.M):
        local=Path(name)
        if local.is_absolute():raise RuntimeError(f'absolute include {src}: {name}')
        if name in ['openmx_common.h','tran_prototypes.h','lapack_prototypes.h','Inputtools.h','flpq_dm.h','mpi.h','omp.h']:continue
        choices=[src.parent/name]+[d/name for d in search]
        found=next((p for p in choices if p.is_file()),None)
        if found is None:raise RuntimeError(f'missing include {src}: {name}')
        put(found)
DST.mkdir(exist_ok=True)
for obj,src in UNITS.items():
    text=None
    if obj=='cpu_nsoc_sym':
        text=src.read_text();needle='/* CHECK uses only sampled O(N) vectors.'
        assert text.count(needle)==1
        text=text.replace(needle,'#include "fused_return_native_v2.inc"\n'+needle)
        head,tail=text.split('void CpuNSym_Solve(int step',1)
        needle='int parallel=sector_parallel_request(step);if(parallel)sector_parallel_prepare(S,H,W);'
        assert tail.count(needle)==1
        tail=tail.replace(needle,needle+'\n int cfs_mode=CfsNativeRequest();(void)CfsNativeSelectedRequest(cfs_mode);')
        old='double tE=MPI_Wtime();backtransform(step,n,H,W,phase);\n double tBack=MPI_Wtime();reorder(eval,ko,maxn,W,H);representative_record(n,maxn,eval);double end=MPI_Wtime();'
        new='double tE=MPI_Wtime(),tBack;\n if(cfs_mode){CfsNativeTail(step,n,maxn,eval,ko,H,W,phase,cfs_mode);tBack=MPI_Wtime();}\n else{backtransform(step,n,H,W,phase);tBack=MPI_Wtime();reorder(eval,ko,maxn,W,H);}\n representative_record(n,maxn,eval);double end=MPI_Wtime();'
        assert tail.count(old)==1
        text=head+'void CpuNSym_Solve(int step'+tail.replace(old,new)
        assert sha(text.encode())==build['overlay_translation_sha256'],'exact selected translation mismatch'
    elif obj=='selected_left':
        text=src.read_text().replace(str(SPARSE/'cpu_nsoc_sparse.frozen.c'),'cpu_nsoc_sparse.frozen.c')
    put(src,text)
units=[{'object':obj+'.o','source':str((DST/src.name).relative_to(REPO)),'language':'c'}for obj,src in UNITS.items()]
manifest={'schema_version':1,'lineage_binary':build['binary'],'lineage_binary_sha256':build['binary_sha256'],'exact_selected_translation_sha256':build['overlay_translation_sha256'],'units':units,'include_dirs':[str(DST.relative_to(REPO))],'files':records,'note':'DFT.direct.c is merge input only; final build uses context-checked NSOC+SOC DFT. Files with duplicate object names must be reconciled with SOC overrides.'}
(HERE/'SOURCES.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(f'Recovered {len(units)} units, {len(records)} sources/includes')
