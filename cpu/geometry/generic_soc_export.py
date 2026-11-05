"""Finite-range GSOC export map: polar AO R and axial spin SU2(det(R) R).

The payload uses the same bounded owner-neighbor projector as legacy SOC,
with version 3 declaring the general (including improper) spin convention.
No full-system dense matrix is allocated.
"""
from pathlib import Path
import hashlib,importlib.util,json,struct
import numpy as np
HERE=Path(__file__).resolve().parent;REPO=HERE.parents[1]
def _load(name,path):
 s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
spinor=_load('gsoc_export_spinor',HERE/'prepare_soc_basis.py')
def spin_lifts(rotations):
 return [spinor.spin_rotation(np.linalg.det(q)*np.asarray(q))for q in rotations]
def build(input_path,soc_metadata,output,provenance):
 inp=Path(input_path).resolve();output=Path(output).resolve();raw=spinor.detector.parse_openmx(inp)
 ops=soc_metadata['operations'];tol=float(soc_metadata['validation'].get('symprec_Ang',1e-5))
 checked=spinor.detector._build_operation_set(raw,[o['W']for o in ops],[o.get('t',[0,0,0])for o in ops],tol,1.25)
 if not checked.get('valid'):raise ValueError('generic SOC export operations fail independent revalidation')
 ops=checked['operations'];ops.sort(key=lambda o:0 if np.array_equal(o['W'],np.eye(3))and np.max(abs(np.asarray(o['t'])-np.rint(o['t'])))<1e-10 else 1)
 if not 1<=len(ops)<=192:raise ValueError('generic SOC export supports 1..192 operations')
 species=list(raw['species_basis']);qs=[np.asarray(o['R_cart_orthogonal'])for o in ops];us=spin_lifts(qs)
 angular=[[np.asarray(spinor.detector._orbital_rotation_block(raw['species_basis'][s],q)).real for s in species]for q in qs]
 counts=[len(d)for d in angular[0]];n=sum(counts[species.index(s)]for s in raw['species'])
 error=0.
 for i,q in enumerate(qs):
  for j,r in enumerate(qs):
   k=int(np.argmin([np.max(abs(v-q@r))for v in qs]))
   if np.max(abs(qs[k]-q@r))>1e-10:raise ValueError('Cartesian export group is not closed')
   # Spin lifts may differ by a sign, but their adjoint action must close.
   error=max(error,float(np.max(abs(np.kron(us[i]@us[j],(us[i]@us[j]).conj())-np.kron(us[k],us[k].conj())))))
   for a,b,c in zip(angular[i],angular[j],angular[k]):
    if np.max(abs(a@b-c))>1e-11:raise ValueError('AO export group representation does not close')
 if error>1e-11:raise ValueError('spin adjoint export group representation does not close')
 chunks=[struct.pack('<8i',0x52534433,3,raw['natoms'],n,len(ops),len(species),0,0),np.asarray(raw['cell'],'<f8').tobytes(),np.asarray(raw['cart'],'<f8').tobytes(),np.asarray([species.index(s)for s in raw['species']],'<i4').tobytes(),np.asarray(counts,'<i4').tobytes()]
 for op,ds in zip(ops,angular):
  for key in ('W','atom_permutation','image_shift'):chunks.append(np.asarray(op[key],'<i4').tobytes())
  chunks.extend(np.asarray(d,'<f8').tobytes()for d in ds)
 chunks.append(np.stack([np.stack([u.real,u.imag],axis=-1)for u in us]).astype('<f8').tobytes())
 if output.exists():raise ValueError('generic SOC map output already exists')
 output.parent.mkdir(parents=True,exist_ok=True);output.write_bytes(b''.join(chunks))
 result={'input':str(inp),'input_sha256':hashlib.sha256(inp.read_bytes()).hexdigest(),'source_basis_sha256':soc_metadata.get('basis_sha256'),'source_symmetry_metadata_sha256':soc_metadata.get('symmetry_metadata_sha256'),'schema':'RSD3-v3-general-SU2','spin_lift':'axial_detR_R','operations':len(ops),'natoms':raw['natoms'],'N_scalar':n,'max_spin_adjoint_group_error':error,'output':str(output),'output_sha256':hashlib.sha256(output.read_bytes()).hexdigest(),'geometry_modified':False,'runtime_operator_equivalence':'required sampled comparison against actual solved sector operator; not certified by metadata'}
 Path(provenance).write_text(json.dumps(result,indent=2)+'\n');return result
