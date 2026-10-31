"""Read-only scientific/performance acceptance; never launches calculations.

H/S parsing reuses the locked OpenMX format12/15 reader from the already verified
benchmark. Every periodic block and every SOC field is compared. No S=I model.
"""
from pathlib import Path
import argparse,hashlib,importlib.util,json,math,re,sys
import numpy as np
READER=Path(__file__).resolve().with_name('hs_union.py')
_spec=importlib.util.spec_from_file_location('locked_hs_union',READER)
_hs=importlib.util.module_from_spec(_spec);_spec.loader.exec_module(_hs)
FLOAT=r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][-+]?\d+)?'
def _number(v):return float(v.replace('D','E').replace('d','e'))
def classify_log(text,stage,exit_code):
    normal='The calculation was normally finished.' in text
    norm=re.findall(r'NormRD\s*=\s*('+FLOAT+r')\s+Criterion\s*=\s*('+FLOAT+r')',text)
    du=re.findall(r'dUele\s*=\s*('+FLOAT+r')',text)
    final=re.findall(r'CPU_NSOC_SYM_FINAL[^\n]*\bconverged=(\d)',text)
    converged=bool(norm and du and abs(_number(norm[-1][0]))<=_number(norm[-1][1]) and abs(_number(du[-1]))<=_number(norm[-1][1]))
    representation='RAW';projected_norm=None
    finals=re.findall(r'CPU_NSOC_SYM_FINAL[^\n]*',text)
    if finals:
        values=dict(re.findall(r'(\w+)=([^\s]+)',finals[-1]))
        representation=values.get('representation','RAW')
        if representation=='PROJECTED'and all(k in values for k in['NormRD_projected','dU','criterion']):
            projected_norm=_number(values['NormRD_projected']);criterion=_number(values['criterion'])
            converged=values.get('converged')=='1'and abs(projected_norm)<=criterion and abs(_number(values['dU']))<=criterion
        else:converged=converged and values.get('converged')=='1'
    steps=re.findall(r'(?:SCF\s*=\s*|SCF=\s*)(\d+)',text)
    # One-step SOC is intentionally a c2n evaluation, not a self-consistency claim.
    passed=exit_code==0 and normal and (converged if stage=='nsoc' else stage=='soc')
    return dict(passed=passed,normal_finish=normal,exit_code=exit_code,scf_converged=converged if stage=='nsoc' else False,
                convergence_representation=representation,projected_NormRD=projected_norm,
                completed_steps=max(map(int,steps))if steps else None,
                final_NormRD=_number(norm[-1][0])if norm else None,criterion=_number(norm[-1][1])if norm else None,
                final_dUele=abs(_number(du[-1]))if du else None,
                category='complete converged NSOC'if stage=='nsoc'else'one-step SOC from own NSOC checkpoint')
def summarize_memory(samples,hosts):
    """Timestamp must be common sample ID; do not combine different-time peaks."""
    batches={};peaks={h:None for h in hosts}
    for row in samples:
        h=row['node']
        if h not in peaks:raise ValueError('unlisted memory node '+h)
        v=int(row['physical_bytes']);t=row.get('sample_id',row.get('timestamp'))
        if t is None or v<0:raise ValueError('invalid memory sample')
        if h in batches.setdefault(t,{}):raise ValueError('duplicate memory sample')
        batches[t][h]=v;peaks[h]=v if peaks[h] is None else max(v,peaks[h])
    complete=[sum(v.values())for v in batches.values()if set(v)==set(hosts)]
    return dict(per_node_peak_bytes=peaks,cluster_synchronized_peak_bytes=max(complete)if complete else None,
                sum_node_peaks_bytes=sum(peaks.values())if all(v is not None for v in peaks.values())else None,
                complete_samples=len(complete),scope='sampled task physical memory; no summed RSS; peak may be missed between samples')
def matched_comparison(reference,candidate,kind):
    required=['physics_sha256','group','stage']
    if kind=='implementation':required+=['layout']
    elif kind=='layout':required+=['binary_sha256']
    elif kind!='physical':raise ValueError('unknown comparison kind')
    differences=[k for k in required if k not in reference or k not in candidate or reference[k]!=candidate[k]]
    return dict(comparable=not differences,different_or_missing= differences,kind=kind,
                scope='historical timing is context only; no cross-hardware attribution')
def energy_force(path):
    text=Path(path).read_text(errors='replace')
    vals=re.findall(r'^\s*Utot\.\s+('+FLOAT+r')',text,re.M)
    if not vals:raise ValueError('missing final total energy '+str(path))
    chunks=re.findall(r'<coordinates.forces\s*\n(.*?)coordinates.forces>',text,re.S)
    if not chunks:raise ValueError('missing native force output '+str(path))
    lines=chunks[-1].strip().splitlines();n=int(lines[0]);force=np.array([[_number(v)for v in l.split()[-3:]]for l in lines[1:]])
    if force.shape!=(n,3)or not np.isfinite(force).all():raise ValueError('invalid force payload')
    return _number(vals[-1]),force

def export_density(path):
    """Complete native DM fields only. EDM is not exported in .scfout."""
    e=_hs.Export(path);count=6 if e.header[1]==15 else 3
    index=len(e.labels)+3*(e.order[0]+1)
    dm=np.memmap(e.path,'<f8',mode='r',offset=e.offset+index*e.element_count*8,shape=(count,e.element_count))
    return e,dm

def compare_density(left,right,electrons):
    a,A=export_density(left);b,B=export_density(right)
    try:
        if a.labels!=b.labels or not np.array_equal(a.geometry,b.geometry)or not np.array_equal(a.orbitals,b.orbitals):raise ValueError('DM representation differs')
        union=a.blocks.keys()|b.blocks.keys();maxdiff=0.;trace_a=trace_b=0.;herm=0.;missing_reverse=0
        for key in union:
            la=a.blocks.get(key);lb=b.blocks.get(key);shape=(la or lb)[1];n=math.prod(shape)
            x=A[:,la[0]:la[0]+n]if la else np.zeros((A.shape[0],n));y=B[:,lb[0]:lb[0]+n]if lb else np.zeros((B.shape[0],n))
            if not np.isfinite(x).all()or not np.isfinite(y).all():raise ValueError('nonfinite DM')
            maxdiff=max(maxdiff,float(abs(x-y).max()))
            if la:trace_a+=float((x[0]+x[1]if A.shape[0]==6 else 2*x[0])@a.arrays[-1,la[0]:la[0]+n])
            if lb:trace_b+=float((y[0]+y[1]if B.shape[0]==6 else 2*y[0])@b.arrays[-1,lb[0]:lb[0]+n])
        for key,(off,shape)in b.blocks.items():
            aa,bb,xx,yy,zz=key;rev=b.blocks.get((bb,aa,-xx,-yy,-zz))
            if rev is None:missing_reverse+=1;continue
            ro,rs=rev;n=math.prod(shape)
            for field in ([0,1,4,5]if B.shape[0]==6 else[0,1,2]):
                x=B[field,off:off+n].reshape(shape);y=B[field,ro:ro+n].reshape(rs).T
                imaginary=field in([4,5]if B.shape[0]==6 else[1,2])
                herm=max(herm,float(abs(x+y if imaginary else x-y).max()))
        return dict(DM_all_exported_fields_max_difference=maxdiff,TrDS_reference=trace_a,TrDS_candidate=trace_b,
                    electron_error_reference=abs(trace_a-electrons),electron_error_candidate=abs(trace_b-electrons),
                    DM_diagonal_spin_Hermiticity_max=herm,reverse_blocks_missing=missing_reverse,
                    EDM_scope='not in scfout; requires native gsoc_density_rank*.bin audit',finite=True)
    finally:A._mmap.close();B._mmap.close();a.close();b.close()

def compare_outputs(reference,candidate,electrons,physical=False):
    ref,can=map(Path,[reference,candidate]);a=next(ref.glob('*.scfout'));b=next(can.glob('*.scfout'))
    ea,fa=energy_force(a.with_suffix('.out'));eb,fb=energy_force(b.with_suffix('.out'))
    if fa.shape!=fb.shape:raise ValueError('force atom count differs')
    hs=_hs.compare(a,b);density=compare_density(a,b,electrons)
    errors=dict(energy_delta_Ha=eb-ea,energy_error_meV_atom=abs(eb-ea)*27211.386245988/len(fa),
                force_max_vector_meV_A=float(np.linalg.norm(fb-fa,axis=1).max())*51422.0674763)
    gates=dict(energy=errors['energy_error_meV_atom']<=0.1,force=errors['force_max_vector_meV_A']<=1.0)if physical else dict(energy=abs(eb-ea)<=1e-8,force=errors['force_max_vector_meV_A']<=1e-3,H=hs['max_H']<=1e-10,S=hs['max_S']<=1e-10,DM=density['DM_all_exported_fields_max_difference']<=1e-10)
    gates.update(electrons=max(density['electron_error_reference'],density['electron_error_candidate'])<=1e-5,
                 Hermiticity=density['DM_diagonal_spin_Hermiticity_max']<=1e-9 and density['reverse_blocks_missing']==0)
    return dict(passed=all(gates.values()),comparison='physical projection accuracy'if physical else'implementation equivalence',gates=gates,errors=errors,HS=hs,density=density,reader_sha256=hashlib.sha256(READER.read_bytes()).hexdigest())

def inspect_workflow(directory):
    directory=Path(directory);wf=json.loads((directory/'workflow.json').read_text());rows={}
    for stage in['nsoc','soc']:
        d=directory/stage;r=json.loads((d/'run.json').read_text());log=(d/'stdout.log').read_text(errors='replace')
        rows[stage]=classify_log(log,stage,r.get('exit_code'))
        rows[stage].update(wall_seconds=r.get('wall_seconds'),launch_wall_seconds=r.get('launch_wall_seconds'),memory=r.get('memory'),HS_output=r.get('HS_output'))
        if not r.get('HS_output')or not Path(r['HS_output']).is_file():rows[stage]['passed']=False;rows[stage]['missing_HS']=True
        if stage=='soc'and rows[stage]['completed_steps']not in(None,1):rows[stage]['passed']=False;rows[stage]['extra_soc_steps']=True
    return dict(passed=wf.get('state')=='complete'and all(r['passed']for r in rows.values()),stages=rows,workflow=wf,
                spectrum_scope='native residual/S-orthonormality only if explicit audit artifacts present; never inferred from energy equality')
def phase_statistics(text):
    """Rank maxima/minima per actual GSOC sector invocation; no summed wall."""
    phases={}
    for line in text.splitlines():
        if not line.startswith('CPU_GSOC_PHASE '):continue
        v=dict(re.findall(r'(\w+)=([^\s]+)',line));rank=int(v['rank']);sector=int(v['sector'])
        # Multiple SCF calls can log same sector/rank; preserve invocation index.
        key=(sector,rank);call=sum(1 for k in phases if k[0]==sector and rank in phases[k]['ranks'])
        group=phases.setdefault((sector,call),{'sector':sector,'invocation':call,'dtype':v['dtype'],'n':int(v['n']),'ranks':{}})
        group['ranks'][rank]={p:float(v[p])for p in['assemble','transform','eigen','back']}
    rows=[]
    for group in phases.values():
        row={k:v for k,v in group.items()if k!='ranks'};row['rank_count']=len(group['ranks'])
        for phase in['assemble','transform','eigen','back']:
            values=[v[phase]for v in group['ranks'].values()]
            row[phase]=dict(rank_max=max(values),rank_min=min(values),rank_spread=max(values)-min(values))
        rows.append(row)
    return dict(GSOC=rows,scope='critical path stage rank maxima and spread; overlapping phase timings are not additive complete wall',
                NSOC_sector_parallel_lines=[l for l in text.splitlines()if l.startswith('CPU_NSOC_SECTOR_PARALLEL_PHASE ')])

def _read_density_audit(path):
    with Path(path).open('rb')as f:
        head=np.fromfile(f,'<i4',8)
        if len(head)!=8 or head[0]!=0x47534441 or head[1]!=1:raise ValueError('unsupported native audit')
        thermal=np.fromfile(f,'<f8',2);eigen=np.fromfile(f,'<f8',int(head[7])*2).reshape(-1,2);blocks={}
        for _ in range(int(head[4])):
            ga,na,count=np.fromfile(f,'<i4',3).tolist()
            for h in range(count):
                gb,nb=np.fromfile(f,'<i4',2).tolist();data=np.fromfile(f,'<f8',na*nb*11).reshape(na*nb,11)
                if not np.isfinite(data).all():raise ValueError('nonfinite native density audit')
                blocks[ga,h,gb,na,nb]=data
        traces=np.fromfile(f,'<f8',2)
        if len(traces)!=2 or f.read(1):raise ValueError('malformed audit trailer')
    return head,thermal,eigen,blocks,traces

def compare_native_audit(reference,candidate):
    left=sorted(Path(reference).glob('gsoc_density_rank*.bin'));right=sorted(Path(candidate).glob('gsoc_density_rank*.bin'))
    if not left or [p.name for p in left]!=[p.name for p in right]:raise ValueError('native audit rank files differ or missing')
    maximum=np.zeros(11);trace=[np.zeros(2),np.zeros(2)];target=[]
    for a,b in zip(left,right):
        H,T,E,A,TR=_read_density_audit(a);h,t,e,B,tr=_read_density_audit(b)
        if not np.array_equal(H,h)or not np.array_equal(T,t)or A.keys()!=B.keys():raise ValueError('native audit layout differs')
        for key in A:maximum=np.maximum(maximum,np.max(abs(A[key]-B[key]),axis=0))
        trace[0]+=TR;trace[1]+=tr;target=[np.array([E[:,1].sum(),E[:,0]@E[:,1]]),np.array([e[:,1].sum(),e[:,0]@e[:,1]])]
    errors=[float(abs(x-y).max())for x,y in zip(trace,target)]
    return dict(passed=float(maximum.max())<=1e-10 and max(errors)<=1e-7,DM_max=float(maximum[:6].max()),EDM_max=float(maximum[6:10].max()),
                S_max=float(maximum[10]),traces=[v.tolist()for v in trace],targets=[v.tolist()for v in target],trace_identity_error=errors,
                scope='complete tiny native 6DM/4EDM components plus overlap, TrDS=sumf,TrEDMS=sumf_epsilon; not application timing')

def main(argv=None):
    p=argparse.ArgumentParser();sub=p.add_subparsers(dest='cmd',required=True)
    s=sub.add_parser('compare');s.add_argument('reference');s.add_argument('candidate');s.add_argument('--electrons',type=float,required=True);s.add_argument('--physical',action='store_true')
    s=sub.add_parser('inspect');s.add_argument('directory')
    args=p.parse_args(argv);result=compare_outputs(args.reference,args.candidate,args.electrons,args.physical)if args.cmd=='compare'else inspect_workflow(args.directory)
    print(json.dumps(result,indent=2));return 0 if result['passed']else 1

if __name__=='__main__':raise SystemExit(main())
