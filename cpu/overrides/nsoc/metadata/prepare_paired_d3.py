#!/usr/bin/env python3
"""Validated, axis-independent D3 Ex/Ey/A1/A2 basis in explicit NSYM v3.

Consume a previously validated operation set without detecting new symmetry,
changing its tolerance, or moving atoms. Dense algebra is restricted to one
atom orbit and angular shell (at most 6*(2*l+1) rows). In particular, the two
E partner columns have a common multiplicity orientation, so a commuting real
symmetric H/S has identical Ex and Ey blocks and zero inter-partner block.
H/S covariance remains a mandatory native runtime check.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import struct
import sys

import numpy as np
from scipy import linalg, sparse

PARENT_TOOLS = Path(__file__).resolve().parents[3]/'geometry/metadata'
MAGIC = 0x4e534433
TOL = 1e-11


def _load_parent():
    spec = importlib.util.spec_from_file_location('paired_d3_parent_auto', PARENT_TOOLS/'prepare_nsoc_auto.py')
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _err(x):
    return float(np.max(np.abs(x), initial=0.))


def local_basis(rotation, c2, tol=TOL):
    """Return dense local U, Ex/Ey/A1/A2 sizes, and strict algebra errors."""
    r=np.asarray(rotation,dtype=float);c=np.asarray(c2,dtype=float)
    if r.ndim!=2 or r.shape[0]!=r.shape[1] or c.shape!=r.shape or not np.isfinite(r).all() or not np.isfinite(c).all():
        raise ValueError('D3 representation must be finite square real matrices')
    n=len(r);eye=np.eye(n)
    errors={'R_orthogonality':_err(r.T@r-eye),'C_orthogonality':_err(c.T@c-eye),
            'C3_cube':_err(r@r@r-eye),'C2_square':_err(c@c-eye),'D3_relation':_err(c@r@c-r.T)}
    if max(errors.values())>tol:
        raise ValueError(f'Invalid orthogonal real D3 representation: {errors}')
    p0=(eye+r+r.T)/3
    def image(p):
        values,vectors=linalg.eigh((p+p.T)*.5)
        if np.max(np.minimum(abs(values),abs(values-1)),initial=0)>tol:
            raise ValueError('D3 projector is not idempotent')
        return vectors[:,values>.5]
    ex=image((eye-p0)@(eye+c)*.5)
    ey=(r-r.T)@ex/np.sqrt(3)
    a1=image(p0@(eye+c)*.5);a2=image(p0@(eye-c)*.5)
    u=np.column_stack((ex,ey,a1,a2));dims=[ex.shape[1],ey.shape[1],a1.shape[1],a2.shape[1]]
    if u.shape!=(n,n):
        raise ValueError(f'Incomplete D3 decomposition: {dims} vs {n}')
    errors.update(orthogonality=_err(u.T@u-eye),E_pair=_err(ey-(r-r.T)@ex/np.sqrt(3)),
                  C2_Ex=_err(c@ex-ex),C2_Ey=_err(c@ey+ey),C3_A=_err(r@np.column_stack((a1,a2))-np.column_stack((a1,a2))),
                  C2_A1=_err(c@a1-a1),C2_A2=_err(c@a2+a2),
                  C3_Ex=_err(r@ex+.5*ex-np.sqrt(3)/2*ey),C3_Ey=_err(r@ey+.5*ey+np.sqrt(3)/2*ex))
    if max(errors.values())>tol:
        raise ValueError(f'Unaligned or nonorthogonal paired D3 basis: {errors}')
    return u,dims,errors


def _validated_group(operations, natoms, species, helper):
    if len(operations)!=6:
        raise ValueError('Paired D3 requires exactly six validated operations')
    # Reject non-integral W before the inherited lookup coerces to integers.
    for op in operations:
        w=np.asarray(op['W'],dtype=float);t=np.asarray(op['t'],dtype=float)
        if w.shape!=(3,3) or not np.isfinite(w).all() or _err(w-np.rint(w))>TOL or t.shape!=(3,) or not np.isfinite(t).all():
            raise ValueError('Invalid fractional group operation')
    _,identity,multiply,_,_=helper._group_tables(operations)
    mult=np.asarray([[multiply(i,j) for j in range(6)] for i in range(6)],dtype=int)
    permutations=[];cart=[]
    for op in operations:
        raw=np.asarray(op['atom_permutation'])
        p=np.asarray(raw,dtype=int)
        if p.shape!=(natoms,) or not np.array_equal(p,raw) or sorted(p.tolist())!=list(range(natoms)):
            raise ValueError('Operation atom map must be a permutation')
        if any(species[a]!=species[p[a]] for a in range(natoms)):
            raise ValueError('Operation changes atomic species')
        q=np.asarray(op['R_cart_orthogonal'],dtype=float)
        if q.shape!=(3,3) or not np.isfinite(q).all() or _err(q.T@q-np.eye(3))>TOL or abs(np.linalg.det(q)-1)>TOL:
            raise ValueError('D3 Cartesian operation is not a proper orthogonal rotation')
        permutations.append(p);cart.append(q)
    if not np.array_equal(permutations[identity],np.arange(natoms)) or _err(cart[identity]-np.eye(3))>TOL:
        raise ValueError('Identity operation representation is invalid')
    errors={'cartesian_composition':0.0,'permutation_composition':0.0}
    for i in range(6):
        for j in range(6):
            k=mult[i,j]
            if not np.array_equal(permutations[i][permutations[j]],permutations[k]):
                raise ValueError(f'Atom permutations violate group multiplication ({i},{j})')
            errors['cartesian_composition']=max(errors['cartesian_composition'],_err(cart[i]@cart[j]-cart[k]))
    if errors['cartesian_composition']>TOL:
        raise ValueError(f'Cartesian rotations violate group multiplication: {errors}')
    order3=[i for i in range(6) if i!=identity and mult[mult[i,i],i]==identity]
    order2=[i for i in range(6) if i!=identity and mult[i,i]==identity]
    if len(order3)!=2 or len(order2)!=3:
        raise ValueError('Operation set does not have D3 element orders')
    r=order3[0]
    candidates=[c for c in order2 if mult[mult[c,r],c]==mult[r,r]]
    if len(candidates)!=3:
        raise ValueError('Operation set violates C R C = inverse(R)')
    return permutations,cart,mult,r,candidates[0],errors


def _write_v3(output,natoms,n,dims,columns,operations,r,c):
    ptr=[0];rows=[];values=[]
    for ri,va in columns:
        mask=abs(va)>1e-13
        rows.extend(np.asarray(ri[mask],dtype='<i4').tolist())
        values.extend(np.asarray(va[mask],dtype='<f8').tolist());ptr.append(len(rows))
    u=sparse.csc_matrix((values,rows,ptr),shape=(n,n))
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('wb') as f:
        f.write(struct.pack('<9i',MAGIC,3,natoms,n,4,*dims));f.write(struct.pack('<q',u.nnz))
        for data,dtype in [(operations[r]['atom_permutation'],'<i4'),(operations[r]['t'],'<f8'),
                           (operations[c]['atom_permutation'],'<i4'),(operations[c]['R_cart_orthogonal'],'<f8'),
                           (operations[c]['t'],'<f8'),(u.indptr,'<i8'),(u.indices,'<i4'),(u.data,'<f8')]:
            f.write(np.asarray(data,dtype=dtype).tobytes())
    return u


def build(input_path, symmetry_json, output, provenance=None):
    """Build from fixed metadata; return {'U','sizes','metadata'} (no dense N²)."""
    inp=Path(input_path).resolve();meta=Path(symmetry_json).resolve();out=Path(output).resolve()
    helper=_load_parent();structure=helper.detector.parse_openmx(inp)
    data=json.loads(meta.read_text());selected=data.get('selected',data)
    if 'input_sha256' in data and data['input_sha256']!=hashlib.sha256(inp.read_bytes()).hexdigest():
        raise ValueError('Metadata input SHA256 does not match supplied input')
    if not selected.get('valid',False) or not selected.get('closure',{}).get('closed',False) or not selected.get('orbital_representation',{}).get('checked',False):
        raise ValueError('D3 basis requires validated closed geometric and orbital metadata')
    if structure['spin_mode']!='none':
        raise ValueError('Paired real D3 basis supports scalar NSOC only')
    if data.get('natoms',structure['natoms'])!=structure['natoms']:
        raise ValueError('Metadata atom count does not match input')
    if 'species' in data and data['species']!=structure['species']:
        raise ValueError('Metadata atomic species order does not match input')
    if 'species_basis' in data and data['species_basis']!=structure['species_basis']:
        raise ValueError('Metadata AO basis order does not match input')
    operations=selected['operations'];perms,cart,mult,ri,ci,group_errors=_validated_group(operations,structure['natoms'],structure['species'],helper)
    shells,offsets,n=helper._atom_offsets(structure)
    provider=helper.detector._load_orbital_rotation_provider()
    if provider is None:
        raise RuntimeError('OpenMX real harmonic rotation provider is unavailable')
    angular={}
    for shell in {l for sh in shells.values() for l,_ in sh}:
        d=2*'spdf'.index(shell)+1;reps=[]
        for q in cart:
            block=np.asarray(provider(shell,q),dtype=complex)
            if block.shape!=(d,d) or not np.isfinite(block).all() or _err(block.imag)>TOL:
                raise ValueError(f'Invalid real AO rotation for {shell}')
            reps.append(block.real)
        for i in range(6):
            if _err(reps[i].T@reps[i]-np.eye(d))>TOL:
                raise ValueError(f'AO rotation is not orthogonal for {shell}')
            for j in range(6):
                if _err(reps[i]@reps[j]-reps[mult[i,j]])>TOL:
                    raise ValueError(f'AO rotations violate group multiplication for {shell}')
        angular[shell]=reps
    seen=set();sectors=[[],[],[],[]];local_checks=[];orbit_sizes=[];max_local=0;basis_errors={}
    for atom in range(structure['natoms']):
        if atom in seen:continue
        orbit=helper._orbit(atom,perms);seen.update(orbit);orbit_sizes.append(len(orbit))
        where={a:i for i,a in enumerate(orbit)};local=0;species=structure['species'][atom]
        for shell,nrad in shells[species]:
            d=2*'spdf'.index(shell)+1;size=len(orbit)*d
            if size>6*d:raise ValueError('Atom orbit exceeds D3 group order')
            reps=[]
            for gi in range(6):
                m=np.zeros((size,size))
                for source in orbit:
                    i0=where[int(perms[gi][source])]*d;j0=where[source]*d
                    m[i0:i0+d,j0:j0+d]=angular[shell][gi]
                reps.append(m)
            u,dims,errors=local_basis(reps[ri],reps[ci])
            for key,value in errors.items():basis_errors[key]=max(basis_errors.get(key,0.),value)
            max_local=max(max_local,size);bounds=np.r_[0,np.cumsum(dims)]
            for radial in range(nrad):
                rows=np.asarray([offsets[a]+local+radial*d+mu for a in orbit for mu in range(d)],dtype=int)
                for sector in range(4):
                    for col in range(bounds[sector],bounds[sector+1]):sectors[sector].append((rows,u[:,col]))
                local_checks.append({'orbit_representative':atom,'orbit_size':len(orbit),'shell':shell,'radial':radial,'sizes':dims})
            local+=nrad*d
    dims=[len(s) for s in sectors]
    if sum(dims)!=n or dims[0]!=dims[1]:raise ValueError('Incomplete paired D3 global basis')
    columns=[column for sector in sectors for column in sector]
    # Validate the sparse payload after coefficient pruning, before writing it.
    ptr=[0];rows=[];values=[]
    for ids,vals in columns:
        mask=abs(vals)>1e-13;rows.extend(ids[mask]);values.extend(vals[mask]);ptr.append(len(rows))
    sparse_u=sparse.csc_matrix((values,rows,ptr),shape=(n,n))
    delta=sparse_u.T@sparse_u-sparse.eye(n,format='csc')
    basis_errors['sparse_payload_orthogonality']=_err(delta.data)
    if basis_errors['sparse_payload_orthogonality']>TOL:raise ValueError('Sparse paired payload lost orthogonality')
    u=_write_v3(out,structure['natoms'],n,dims,columns,operations,ri,ci)
    info={'schema':'openmx-paired-d3-nsoc-sector-basis-v3','format_version':3,'input':str(inp),
          'input_sha256':hashlib.sha256(inp.read_bytes()).hexdigest(),'symmetry_metadata':str(meta),
          'symmetry_metadata_sha256':hashlib.sha256(meta.read_bytes()).hexdigest(),
          'generator_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
          'operation_count':6,'point_group':selected.get('point_group'),'N':n,'natoms':structure['natoms'],
          'sector_sizes':dims,'dimensions':dims,'distinct_solve_dimensions':[dims[0],*dims[2:]],'nnz':u.nnz,
          'generator_operation_indices':{'C3':ri,'C2':ci},'C3_rotation':cart[ri].tolist(),'C2_rotation':cart[ci].tolist(),
          'C3_shift_fractional':operations[ri]['t'],'C2_shift_fractional':operations[ci]['t'],
          'basis_kind':'real_paired_d3_equivalent_rows','column_order':'Ex,Ey,A1,A2; Ey=(R-R.T)Ex/sqrt(3)',
          'basis_errors':basis_errors,'group_errors':group_errors,'max_local_block_dimension':max_local,
          'local_block_count':len(local_checks),'atom_orbit_sizes':orbit_sizes,'local_checks':local_checks,
          'geometry_modified':False,'symmetry_redetected':False,'symmetry_tolerance_modified':False,
          'validation':{'group_closed':True,'orbital_representation_checked':True,'canonical_d3_rows_checked':True,
                        'equivalent_E_orientation_checked':True,'h_s_covariance_checked':False},
          'payload':str(out),'payload_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'payload_bytes':out.stat().st_size,
          'environment_prefix':sys.prefix,'host':os.uname().nodename,'working_directory':os.getcwd(),
          'python_executable':sys.executable,'python_version':sys.version,'numpy_version':np.__version__,
          'command_argv':sys.argv,'process_threads':'one process, BLAS/OMP threads 1',
          'calculation_category':'basis construction and mathematical validation; not an SCF result',
          'final_process_state':'paired basis generated; native H/S covariance checks still required'}
    if provenance is not None:
        p=Path(provenance);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(info,indent=2)+'\n')
    return {'U':u,'sizes':dims,'metadata':info}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,required=True);parser.add_argument('--symmetry-json',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--provenance',type=Path)
    args=parser.parse_args();result=build(args.input,args.symmetry_json,args.output,args.provenance)
    print(json.dumps(result['metadata'],indent=2))

if __name__=='__main__':main()
