#!/usr/bin/env python3
"""Build sparse NSYM v4 canonical C3 pairs, requiring a complex native solver.

The first two real column groups X,Y are NOT separate invariant sectors.
For a commuting real symmetric H they form [[A,-B],[B,A]], B antisymmetric.
Only the third group is an independent real invariant sector. Older binaries
reject version 4, so this format cannot silently invoke their real-half path.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
import struct
from pathlib import Path
import numpy as np
_pair_spec = importlib.util.spec_from_file_location('openmx_scalar_c3_pair_math', Path(__file__).with_name('c3_pair_math.py'))
_pair_module = importlib.util.module_from_spec(_pair_spec)
_pair_spec.loader.exec_module(_pair_module)
paired_basis = _pair_module.paired_basis

DEFAULT_HELPER = Path(__file__).resolve().parents[4]/'local-records/cpu_d3/historical/geometry_tools/prepare_nsoc_auto.py'


def _load(path):
    spec = importlib.util.spec_from_file_location('paired_c3_generic_helper',path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_v4(output, natoms, n, sizes, columns):
    """Same sparse storage/header length as v2; only version and meaning differ."""
    if len(sizes)!=3 or sizes[0]!=sizes[1] or min(sizes)<1 or sum(sizes)!=n:
        raise ValueError('v4 requires positive [m,m,nA] dimensions')
    ptr,rows,values = [0],[],[]
    for row,value in columns:
        mask = abs(value)>1e-13
        rows.extend(np.asarray(row[mask],dtype='<i4').tolist())
        values.extend(np.asarray(value[mask],dtype='<f8').tolist())
        ptr.append(len(rows))
    if len(ptr)!=n+1:
        raise ValueError('number of sparse columns differs from N')
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('wb') as f:
        f.write(struct.pack('<9i',0x4E534433,4,natoms,n,3,*sizes,0))
        f.write(struct.pack('<q',len(rows)))
        ident = np.arange(natoms,dtype='<i4').tobytes()
        f.write(ident)
        f.write(np.zeros(3,dtype='<f8').tobytes())
        f.write(ident)
        f.write(np.eye(3,dtype='<f8').tobytes())
        f.write(np.zeros(3,dtype='<f8').tobytes())
        f.write(np.asarray(ptr,dtype='<i8').tobytes())
        f.write(np.asarray(rows,dtype='<i4').tobytes())
        f.write(np.asarray(values,dtype='<f8').tobytes())
    return len(rows)


def build(input_path, symmetry_json, output, provenance=None, helper_path=DEFAULT_HELPER):
    input_path,symmetry_json,output=map(Path,(input_path,symmetry_json,output))
    helper = _load(Path(helper_path))
    detector = helper.detector
    structure = detector.parse_openmx(input_path)
    metadata = json.loads(symmetry_json.read_text())
    if metadata.get('input_sha256') and metadata['input_sha256']!=sha(input_path):
        raise ValueError('symmetry metadata was validated for a different input')
    selected = metadata.get('selected',metadata)
    ops = selected.get('operations',[])
    if (not selected.get('valid') or len(ops)!=3 or
        not selected.get('closure',{}).get('closed') or
        not selected.get('orbital_representation',{}).get('checked')):
        raise ValueError('requires validated closed three-operation scalar C3 metadata')
    if structure['spin_mode']!='none':
        raise ValueError('C3 real conjugate-pair optimization requires scalar NSOC')
    if any(np.linalg.det(np.asarray(op['R_cart_orthogonal']))<.999999 for op in ops):
        raise ValueError('C3 operations must be proper rotations')
    _,identity,multiply,_,_ = helper._group_tables(ops)
    generators=[i for i in range(3) if i!=identity and multiply(multiply(i,i),i)==identity]
    if len(generators)!=2:
        raise ValueError('operation set is not cyclic C3')
    gi=generators[0]
    inverse=multiply(gi,gi)
    perms=[np.asarray(op['atom_permutation'],dtype=int) for op in ops]
    natoms=structure['natoms']
    if any(sorted(p.tolist())!=list(range(natoms)) for p in perms):
        raise ValueError('invalid atom-permutation metadata')
    if not np.array_equal(perms[gi][perms[gi]],perms[inverse]):
        raise ValueError('C3 permutation multiplication failed')
    shells,offsets,n = helper._atom_offsets(structure)
    provider=detector._load_orbital_rotation_provider()
    if provider is None:
        raise RuntimeError('OpenMX real orbital rotation provider unavailable')
    orbital={}
    for species,species_shells in shells.items():
        for l,nr in species_shells:
            d=2*'spdf'.index(l)+1
            for opid in (gi,inverse):
                rot=np.asarray(provider(l,np.asarray(ops[opid]['R_cart_orthogonal'])))
                if rot.shape!=(d,d) or np.max(abs(rot.imag),initial=0)>1e-12:
                    raise ValueError('orbital representation is not real')
                orbital[species,l,opid]=rot.real
    columns=[[],[],[]]
    seen=set()
    checks=[]
    for atom in range(natoms):
        if atom in seen:
            continue
        orbit=helper._orbit(atom,perms)
        seen.update(orbit)
        species=structure['species'][atom]
        where={a:i for i,a in enumerate(orbit)}
        base=0
        for l,nr in shells[species]:
            d=2*'spdf'.index(l)+1
            size=len(orbit)*d
            reps=[]
            for opid in (gi,inverse):
                r=np.zeros((size,size))
                for a in orbit:
                    dest=int(perms[opid][a]); i,j=where[dest]*d,where[a]*d
                    r[i:i+d,j:j+d]=orbital[species,l,opid]
                reps.append(r)
            if np.max(abs(reps[0]@reps[0]-reps[1]),initial=0)>2e-11:
                raise ValueError('AO C3 generator/inverse multiplication failed')
            x,y,trivial=paired_basis(reps[0])
            u=np.column_stack([x,y,trivial])
            m=x.shape[1];n0=trivial.shape[1]
            c,s=-.5,np.sqrt(3)/2
            target=np.block([[c*np.eye(m),-s*np.eye(m),np.zeros((m,n0))],
                             [s*np.eye(m),c*np.eye(m),np.zeros((m,n0))],
                             [np.zeros((n0,m)),np.zeros((n0,m)),np.eye(n0)]])
            residual=float(np.max(abs(u.T@reps[0]@u-target),initial=0))
            orth=float(np.max(abs(u.T@u-np.eye(size)),initial=0))
            if max(residual,orth)>2e-11:
                raise ValueError('local canonical C3 pair check failed')
            for radial in range(nr):
                rows=np.asarray([offsets[a]+base+radial*d+mu for a in orbit for mu in range(d)])
                for bucket,part in zip(columns,(x,y,trivial)):
                    bucket.extend((rows,part[:,j]) for j in range(part.shape[1]))
                checks.append(dict(orbit=orbit,shell=l,radial=radial,m=m,nA=n0,
                                   canonical_residual=residual,orthogonality=orth))
            base+=nr*d
    sizes=[len(part) for part in columns]
    if sum(sizes)!=n:
        raise ValueError('local orbital pairs do not span the complete basis')
    nnz=write_v4(output,natoms,n,sizes,[entry for part in columns for entry in part])
    info=dict(schema='openmx-paired-c3-nsym-v4',format_version=4,N=n,natoms=natoms,
        point_group='C3',operation_count=3,generator_operation=gi,
        invariant_sector_sizes=[2*sizes[0],sizes[2]],storage_column_groups=sizes,nnz=nnz,
        basis_kind='canonical_real_conjugate_pairs',
        canonical_rotation='[[cos(theta)I,-sin(theta)I],[sin(theta)I,cos(theta)I]], theta=2pi/3',
        commuting_symmetric_block='[[A,-B],[B,A]], A symmetric, B antisymmetric',
        native_solver_requirement='complex m-by-m generalized Hermitian solve, both real partners restored',
        input=str(input_path.resolve()),input_sha256=sha(input_path),
        symmetry_metadata=str(symmetry_json.resolve()),symmetry_metadata_sha256=sha(symmetry_json),
        basis=str(output.resolve()),basis_sha256=sha(output),
        generator_sha256=sha(__file__),math_sha256=sha(Path(__file__).with_name('c3_pair_math.py')),
        helper=str(Path(helper_path).resolve()),helper_sha256=sha(helper_path),
        local_checks=checks,
        validation=dict(group_closed=True,orbital_representation_checked=True,
            maximum_canonical_rotation_residual=max(x['canonical_residual'] for x in checks),
            maximum_orthogonality_error=max(x['orthogonality'] for x in checks),
            actual_H_S_covariance_checked=False,native_solver_validated=False),
        final_process_state='sparse C3 pair basis generated; native H/S validation remains required')
    if provenance:
        Path(provenance).parent.mkdir(parents=True,exist_ok=True)
        Path(provenance).write_text(json.dumps(info,indent=2)+'\n')
    print(json.dumps({k:v for k,v in info.items() if k!='local_checks'},indent=2),flush=True)
    return info


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path,required=True)
    p.add_argument('--symmetry-json',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--provenance',type=Path)
    p.add_argument('--helper',type=Path,default=DEFAULT_HELPER)
    a=p.parse_args()
    build(a.input,a.symmetry_json,a.output,a.provenance,a.helper)


if __name__=='__main__':
    main()
