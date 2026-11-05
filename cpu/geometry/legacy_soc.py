"""Explicit adapter for the verified legacy Gamma C3/D3 one-step SOC path.

prepare(..., export_policy='raw') never enables native export projection.
The historical TR/real fast path requires export_policy='symmetrize'; this is
an explicit physical representation choice, not an implied side effect.
"""
from pathlib import Path
import hashlib
import importlib.util
import json
import struct

import numpy as np

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[1]
def _module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod

_math=_module('legacy_spinor_math',HERE/'legacy_helpers/soc_spinor.py')
Q3=_math.Q3;U3=_math.U3
build_spinor_basis=_math.build_spinor_basis
maxerr=_math.maxerr

def _sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def _require(value,message):
    if not value:raise ValueError(message)

def _input_guard(text):
    values={}
    for line in text.splitlines():
        fields=line.split('#',1)[0].split()
        if len(fields)>1 and not fields[0].startswith('<'):
            values[fields[0].lower()]=' '.join(fields[1:]).lower()
    for key,val in [('scf.spinpolarization','nc'),('scf.spinorbit.coupling','on'),
                    ('scf.restart','c2n'),('scf.maxiter','1'),('scf.eigenvaluesolver','cluster')]:
        _require(values.get(key)==val,'legacy SOC requires fixed nonmagnetic Gamma one-step c2n Cluster SOC: '+key)
    _require(values.get('scf.kgrid','1 1 1').split()==['1','1','1'],'legacy SOC requires Gamma')
    _require(values.get('md.type','nomd') in ('nomd','none'),'legacy SOC requires fixed geometry')
    for key in ('scf.hubbard.u','scf.constraint.nc.spin','scf.zeeman.nc.spin','scf.zeeman.nc.orbital'):
        _require(values.get(key,'off')=='off','legacy SOC does not support '+key)

def _ordered_operations(raw,selected,group):
    """Revalidate mappings, then retain the historical +120deg-z gauge."""
    detector=_module('legacy_metadata_detector',HERE/'metadata/symmetry_detect.py')
    entries=selected.get('operations',[])
    _require(selected.get('valid') and entries,'validated symmetry operations required')
    original_entries=entries
    tolerance=float(selected.get('symprec_Ang',1e-5))
    bound=float(selected.get('map_tol_Ang',tolerance*1.25))
    _require(np.isfinite(bound) and 0<bound<=.02,'invalid geometry mapping bound')
    checked=detector._build_operation_set(raw,[o['W']for o in entries],
        [o.get('t',o.get('translation_frac',[0,0,0]))for o in entries],tolerance,bound/tolerance)
    _require(checked.get('valid'),'SOC operations fail geometry or closure validation: '+str(checked.get('error')))
    entries=checked['operations']
    def rotation(op):return np.asarray(op['R_cart_orthogonal'])
    def nearest(q):
        op=min(entries,key=lambda o:np.max(abs(rotation(o)-q)))
        _require(np.max(abs(rotation(op)-q))<1e-3,'legacy path requires z-oriented hexagonal C3/D3; use general GSOC for another orientation')
        return op
    qs=[np.eye(3),Q3,Q3@Q3];us=[np.eye(2,dtype=complex),U3,U3@U3]
    if group=='D3':
        # Preserve the chosen source C2 gauge, including the historical verified
        # ordering. The detector sorts its revalidated operation set internally.
        position={tuple(np.asarray(o['W']).ravel()):i for i,o in enumerate(original_entries)}
        choices=[o for o in entries if np.linalg.det(rotation(o))>.99 and rotation(o)[2,2]<-.99]
        choices.sort(key=lambda o:position[tuple(np.asarray(o['W']).ravel())])
        _require(choices,'D3 requires a proper in-plane C2 generator')
        source=original_entries[position[tuple(np.asarray(choices[0]['W']).ravel())]]
        qraw=np.asarray(source.get('rotation_cart_nominal',rotation(choices[0])))
        _require(np.max(abs(qraw-rotation(choices[0])))<1e-3,'source nominal C2 differs from validated operation')
        _,v=_math.linalg.eigh((qraw+qraw.T)/2);axis=v[:,-1]
        axis[2]=0;axis/=np.linalg.norm(axis)
        if axis[np.argmax(abs(axis))]<0:axis=-axis
        q2=2*np.outer(axis,axis)-np.eye(3)
        u2=-1j*np.array([[0,axis[0]-1j*axis[1]],[axis[0]+1j*axis[1],0]])
        qs.extend([np.linalg.matrix_power(Q3,k)@q2 for k in range(3)])
        us.extend([np.linalg.matrix_power(U3,k)@u2 for k in range(3)])
    ordered=[nearest(q)for q in qs]
    _require(len({tuple(np.asarray(o['W']).ravel())for o in ordered})==len(qs),'operation correspondence not bijective')
    # Native solver also tests the slab and the special ordered hexagonal cell.
    cell=raw['cell'];flat=abs(cell[0,2])+abs(cell[1,2])+abs(cell[2,0])+abs(cell[2,1])
    _require(flat<=4*bound and np.linalg.norm(Q3@cell[0]-cell[1])<=bound,
             'legacy path requires the established ordered hexagonal slab cell; use general GSOC')
    largest=0.
    for op,q in zip(ordered,qs):
        xyz=raw['cart'];perm=np.asarray(op['atom_permutation']);t=np.asarray(op['t'])
        delta=(xyz@q.T+t@cell-xyz[perm])@np.linalg.inv(cell);delta-=np.rint(delta)
        residual=float(np.max(np.linalg.norm(delta@cell,axis=1)))
        largest=max(largest,residual)
    _require(largest<=bound,'exact legacy angular convention exceeds accepted geometry mapping bound')
    return ordered,qs,us,bound,largest

def _write_basis(path,raw,ops,qs,result):
    dims=result['dimensions'];u=result['U'];ng=len(ops);n=result['N']
    with path.open('xb')as f:
        f.write(struct.pack('<'+('8i'if ng==3 else'9i'),0x43334253,1 if ng==3 else 2,
                            raw['natoms'],n,len(dims),*dims))
        f.write(struct.pack('<q',u.nnz))
        f.write(np.asarray(ops[1]['atom_permutation'],'<i4').tobytes())
        f.write(np.asarray(ops[1]['t'],'<f8').tobytes())
        if ng==6:
            f.write(np.asarray(ops[3]['atom_permutation'],'<i4').tobytes())
            f.write(np.asarray(qs[3],'<f8').tobytes());f.write(np.asarray(ops[3]['t'],'<f8').tobytes())
        f.write(np.asarray(u.indptr,'<i8').tobytes());f.write(np.asarray(u.indices,'<i4').tobytes())
        f.write(np.asarray(u.data,'<c16').tobytes())

def _write_map(path,raw,ops,qs,us,n):
    species=list(dict.fromkeys(raw['species']))
    blocks=[_math.orbital_block(raw['species_basis'][s],np.eye(3))for s in species]
    chunks=[struct.pack('<8i',0x52534433,2,raw['natoms'],n,len(ops),len(species),0,0),
            np.asarray(raw['cell'],'<f8').tobytes(),np.asarray(raw['cart'],'<f8').tobytes(),
            np.asarray([species.index(s)for s in raw['species']],'<i4').tobytes(),
            np.asarray([b.shape[0]for b in blocks],'<i4').tobytes()]
    angular=[]
    for op,q in zip(ops,qs):
        for key in ('W','atom_permutation','image_shift'):chunks.append(np.asarray(op[key],'<i4').tobytes())
        ds=[_math.orbital_block(raw['species_basis'][s],q)for s in species];angular.append(ds)
        chunks.extend(np.asarray(d,'<f8').tobytes()for d in ds)
    chunks.append(np.stack([np.stack([u.real,u.imag],axis=-1)for u in us]).astype('<f8').tobytes())
    pauli=[np.array([[0,1],[1,0]],complex),np.array([[0,-1j],[1j,0]],complex),np.diag([1.,-1.])]
    error=0.
    for q,u in zip(qs,us):
        for j,s in enumerate(pauli):error=max(error,float(np.max(abs(u@s@u.conj().T-sum(q[i,j]*pauli[i]for i in range(3))))))
    for i,qi in enumerate(qs):
        for j,qj in enumerate(qs):
            k=int(np.argmin([np.max(abs(q-qi@qj))for q in qs]))
            _require(np.max(abs(qs[k]-qi@qj))<1e-12,'nominal angular group not closed')
            for a,b,c in zip(angular[i],angular[j],angular[k]):error=max(error,float(np.max(abs(a@b-c))))
    _require(error<1e-12,'orbital or SU2 adjoint representation check failed')
    with path.open('xb')as f:f.write(b''.join(chunks))
    return error

def prepare(input_soc,symmetry_json,output_dir,group,*,export_policy='raw'):
    """Return paths and environment additions; never launch or edit input files.

    Default raw exports disable the legacy TR/real acceleration because its
    existing native validation hook necessarily symmetrizes exported H/S.
    Explicit ``export_policy='symmetrize'`` reproduces that verified fast path.
    A full runtime profile and independent converged NSOC checkpoint remain
    required. Runtime guards remain authoritative for magnetic/auxiliary modes.
    """
    _require(group in ('D3','C3'),'legacy SOC supports only explicit D3 or C3')
    _require(export_policy in ('raw','symmetrize'),'export policy must be raw or symmetrize')
    inp=Path(input_soc).resolve();meta=Path(symmetry_json).resolve();out=Path(output_dir).resolve()
    _require(not out.exists(),'output must be a new metadata directory')
    _input_guard(inp.read_text())
    detector=_module('legacy_input_detector',HERE/'metadata/symmetry_detect.py')
    raw=detector.parse_openmx(inp);metadata=json.loads(meta.read_text());selected=metadata.get('selected',metadata)
    ops,qs,us,bound,residual=_ordered_operations(raw,selected,group)
    result=build_spinor_basis(raw['species'],raw['species_basis'],ops[1]['atom_permutation'],
        None if group=='C3'else ops[3]['atom_permutation'],None if group=='C3'else qs[3])
    out.mkdir(parents=True)
    basis=out/'spinor.basis';map_path=out/'operations.socrs'
    _write_basis(basis,raw,ops,qs,result)
    error=_write_map(map_path,raw,ops,qs,us,result['N'])
    project=export_policy=='symmetrize'
    env={'OPENMX_CPU_GSOC':'0','OPENMX_CPU_GSOC_COMPACT':'0','OPENMX_CPU_C3':'3',
         'OPENMX_CPU_D3':'1'if group=='D3'else'0','OPENMX_CPU_C3_BASIS':str(basis),
         'OPENMX_CPU_C3_GEOMETRY_TOL_ANG':format(bound,'.17g'),
         'OPENMX_CPU_SOC_TR':'1'if project and group=='D3'else'0',
         'OPENMX_CPU_SOC_REAL':'1'if project and group=='D3'else'0'}
    if project:
        env.update(OPENMX_CPU_SOC_NATIVE_GROUP=group,OPENMX_CPU_SOC_NATIVE_MAP=str(map_path))
    info={'schema':'legacy-soc-adapter-v1','status':'prepared','group':group,'basis':str(basis),'map':str(map_path),
          'export_HS':'symmetrized'if project else'raw','env_additions':env,
          'files':{'soc_basis':str(basis),'soc_native_map':str(map_path)},
          'input_soc':str(inp),'input_soc_sha256':_sha(inp),'symmetry_json':str(meta),'symmetry_json_sha256':_sha(meta),
          'basis_sha256':_sha(basis),'map_sha256':_sha(map_path),'geometry_modified':False,'basis_modified':False,
          'basis_sector_dimensions':result['dimensions'],'N_scalar':result['N'],'N_spinor':2*result['N'],
          'geometry_mapping_bound_Ang':bound,'geometry_residual_Ang':residual,
          'checks':{**result['checks'],'map_SU2_angular_error':error},
          'full_runtime_profile_required':True,'runtime_HS_covariance_checked':False,
          'scope':'fixed nonmagnetic Gamma NSOC -> one SOC; legacy z-oriented hexagonal cell only; metadata validation'}
    (out/'LEGACY_SOC.json').write_text(json.dumps(info,indent=2)+'\n')
    return info
