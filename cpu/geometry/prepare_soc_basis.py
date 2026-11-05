#!/usr/bin/env python3
"""Sparse Gamma nonmagnetic SOC basis from validated crystallographic operations.

Select proper cyclic Cn/dihedral Dn (n=2,3,4,6), or explicit Ci/Cs subgroups.
Dn fixes the real structure C2_perp*T within each rotation phase sector.
All dense algebra is confined to one atom orbit and one angular/radial shell;
there is never a dense full-system basis. Other point groups may use these
explicit subgroups; this does not claim their full irreducible decomposition.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
import struct
import sys
from pathlib import Path
import numpy as np
from scipy import sparse

HERE=Path(__file__).resolve().parent
GEOMETRY=HERE/'metadata'

def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    return mod

detector=load('generic_soc_detector',GEOMETRY/'symmetry_detect.py')
nsoc=load('generic_soc_scalar_tools',GEOMETRY/'prepare_nsoc_auto.py')
MAGIC=0x47534F43

def error(a):return float(np.max(np.abs(a),initial=0))

def spin_c2(axis):
    axis=np.asarray(axis,dtype=float)
    if axis.shape!=(3,) or not np.all(np.isfinite(axis)) or np.linalg.norm(axis)<1e-14:
        raise ValueError('C2 axis must be a nonzero finite Cartesian vector')
    x,y,z=axis/np.linalg.norm(axis)
    return -1j*np.array([[z,x-1j*y],[x+1j*y,-z]])

def time_reversal(n):
    return np.kron(np.array([[0.,1.],[-1.,0.]]),np.eye(n))

def c2_axis(rotation,tol=1e-9):
    r=np.asarray(rotation,dtype=float)
    if r.shape!=(3,3) or not np.all(np.isfinite(r)) or error(r.T@r-np.eye(3))>tol or abs(np.linalg.det(r)-1)>tol or abs(np.trace(r)+1)>tol:
        raise ValueError('expected a proper Cartesian pi rotation')
    _,v=np.linalg.eigh((r+r.T)*.5)
    axis=v[:,-1]
    k=int(np.argmax(abs(axis)))
    if axis[k]<0:axis=-axis
    return axis

def spin_rotation(rotation,tol=1e-9):
    """SU(2) lift of a proper Cartesian rotation with angle in [0,pi]."""
    r=np.asarray(rotation,dtype=float)
    if r.shape!=(3,3) or error(r.T@r-np.eye(3))>tol or abs(np.linalg.det(r)-1)>tol:
        raise ValueError('SU2 lift requires a proper Cartesian rotation')
    skew=np.array([r[2,1]-r[1,2],r[0,2]-r[2,0],r[1,0]-r[0,1]])*.5
    cosine=float(np.clip((np.trace(r)-1)*.5,-1,1));sine=float(np.linalg.norm(skew))
    angle=float(np.arctan2(sine,cosine))
    if angle<tol:return np.eye(2,dtype=complex)
    if abs(angle-np.pi)<tol:axis=c2_axis(r)
    else:axis=skew/sine
    return np.cos(angle*.5)*np.eye(2)+np.sin(angle*.5)*spin_c2(axis)

def _real_gauge(q,a,tol):
    """Fix an antiunitary a*K with square +I on an invariant subspace."""
    rank=q.shape[1]
    if rank==0:return q
    m=q.conj().T@a@q.conj()
    if error(a@a.conj()-np.eye(a.shape[0]))>tol or error(a@q.conj()-q@m)>tol:
        raise ValueError('antiunitary real structure does not preserve the sector')
    embedding=np.block([[m.real,m.imag],[m.imag,-m.real]])
    ev,w=np.linalg.eigh((embedding+embedding.T)*.5);fixed=w[:,ev>0]
    if fixed.shape!=(2*rank,rank):raise ValueError('real structure fixed-space rank mismatch')
    result=q@(fixed[:rank]+1j*fixed[rank:])
    if error(a@result.conj()-result)>tol:raise ValueError('failed antiunitary real gauge')
    return result

def local_cyclic_basis(r,order,perpendicular=None,tol=2e-10):
    """All cyclic double-group sectors, optional dihedral real structure.

    Real phases are full independent sectors: T squares to -I there and does
    not license solving half of that sector. Only distinct conjugate phases
    are marked as reconstructible partners.
    """
    if order not in (2,3,4,6):raise ValueError('supported cyclic orders are 2,3,4,6')
    r=np.asarray(r,dtype=complex)
    if r.ndim!=2 or r.shape[0]!=r.shape[1] or r.shape[0]%2:raise ValueError('invalid spinor rotation shape')
    dim=r.shape[0];t=time_reversal(dim//2)
    if error(r.conj().T@r-np.eye(dim))>tol or error(np.linalg.matrix_power(r,order)+np.eye(dim))>tol:
        raise ValueError('cyclic spinor generator must be unitary and its nth power -I')
    if error(r@t-t@r.conj())>tol:raise ValueError('generator must commute with time reversal')
    a=None
    if perpendicular is not None:
        b=np.asarray(perpendicular,dtype=complex)
        if b.shape!=r.shape or error(b@b+np.eye(dim))>tol or error(b@r@b.conj().T-r.conj().T)>tol or error(b.conj().T@b-np.eye(dim))>tol:
            raise ValueError('dihedral C2 must square to -I and invert the cyclic generator')
        a=b@t
    phases=[complex(np.exp(1j*(2*k+1)*np.pi/order)) for k in range(order)]
    independent=[k for k,z in enumerate(phases) if z.imag>tol or abs(z.imag)<=tol]
    independent.sort(key=lambda k:float(np.angle(phases[k])%(2*np.pi)))
    source_indices=independent+[k for k in range(order) if phases[k].imag < -tol]
    source_indices=list(dict.fromkeys(source_indices))
    powers=[np.linalg.matrix_power(r,k) for k in range(order)]
    basis=[];labels=[];partners=[];dtypes=[]
    by_index={}
    for k in source_indices:
        phase=phases[k]
        if phase.imag < -tol:
            source=min(independent,key=lambda j:abs(phases[j]-phase.conjugate()))
            if abs(phases[source]-phase.conjugate())>tol:raise ValueError('missing conjugate phase')
            p=by_index[source];q=t@basis[p].conj()
        else:
            projector=sum((phase**(-j))*power for j,power in enumerate(powers))/order
            if error(projector-projector.conj().T)>tol:raise ValueError('cyclic projector not Hermitian')
            vals,v=np.linalg.eigh((projector+projector.conj().T)*.5);q=v[:,vals>.5]
            if a is not None:q=_real_gauge(q,a,tol)
            p=-1
        if error(r@q-phase*q)>tol:raise ValueError('cyclic phase sector failed invariance')
        by_index[k]=len(basis);basis.append(q);labels.append(phase);partners.append(p);dtypes.append(1 if a is not None else 2)
    full=np.column_stack(basis)
    orth=error(full.conj().T@full-np.eye(dim))
    if full.shape!=(dim,dim) or orth>tol:raise ValueError('cyclic sectors do not span an orthonormal full spinor basis')
    return {'basis':basis,'phase':labels,'partner':partners,'dtype':dtypes,'orthogonality_error':orth}

def local_improper_basis(r,kind,tol=2e-10):
    """Mirror ±i partners or full inversion-parity self-conjugate sectors."""
    r=np.asarray(r,dtype=complex)
    if r.ndim!=2 or r.shape[0]!=r.shape[1] or r.shape[0]%2:raise ValueError('invalid improper spinor representation')
    dim=r.shape[0];t=time_reversal(dim//2)
    if error(r.conj().T@r-np.eye(dim))>tol or error(r@t-t@r.conj())>tol:
        raise ValueError('improper operation must be unitary and commute with T')
    if kind=='Cs':
        q,info=local_sector_basis(r,tol=tol)
        return {'basis':[q,t@q.conj()],'phase':[1j,-1j],'partner':[-1,0],'dtype':[2,2],'orthogonality_error':info['orthogonality_error']}
    if kind!='Ci' or error(r@r-np.eye(dim))>tol or error(r-r.conj().T)>tol:
        raise ValueError('inversion spinor representation must be a Hermitian involution')
    vals,q=np.linalg.eigh((r+r.conj().T)*.5)
    basis=[q[:,vals>0],q[:,vals<0]]
    if any(b.shape[1]%2 for b in basis):raise ValueError('parity sector must retain its full Kramers dimension')
    full=np.column_stack(basis);orth=error(full.conj().T@full-np.eye(dim))
    if orth>tol:raise ValueError('parity sectors are not a unitary complete basis')
    return {'basis':basis,'phase':[1.+0j,-1.+0j],'partner':[-1,-1],'dtype':[2,2],'orthogonality_error':orth}

def local_sector_basis(r,perpendicular=None,tol=2e-10):
    """Return one +i sector; optionally choose the A=Rperp*T fixed real gauge."""
    r=np.asarray(r,dtype=complex);dim=r.shape[0];n=dim//2
    if r.shape!=(dim,dim) or dim%2 or error(r.conj().T@r-np.eye(dim))>tol or error(r@r+np.eye(dim))>tol:
        raise ValueError('invalid spinor C2 representation: must be unitary and square to -I')
    t=time_reversal(n)
    if error(r@t-t@r.conj())>tol:raise ValueError('rotation does not commute with time reversal')
    vals,vecs=np.linalg.eigh((-1j*r+(-1j*r).conj().T)*.5)
    q=vecs[:,vals>0]
    if q.shape!=(dim,n) or error(r@q-1j*q)>tol:
        raise ValueError('C2 sectors must each contain exactly N scalar orbitals')
    real=perpendicular is not None
    if real:
        p=np.asarray(perpendicular,dtype=complex)
        if p.shape!=r.shape or error(p.conj().T@p-np.eye(dim))>tol or error(p@p+np.eye(dim))>tol or error(r@p+p@r)>tol:
            raise ValueError('second C2 must be unitary and anticommute in the spinor representation')
        a=p@t
        if error(a@a.conj()-np.eye(dim))>tol:raise ValueError('antiunitary real structure must square to +I')
        m=q.conj().T@a@q.conj()
        if error(a@q.conj()-q@m)>tol:raise ValueError('real structure does not preserve C2 sector')
        # Real embedding of m*K. Its +1 eigenspace is a real orthonormal
        # fixed basis; complex inner products are real within this space.
        embedding=np.block([[m.real,m.imag],[m.imag,-m.real]])
        ev,w=np.linalg.eigh((embedding+embedding.T)*.5)
        fixed=w[:,ev>0]
        if fixed.shape!=(2*n,n):raise ValueError('real structure fixed-space rank mismatch')
        q=q@(fixed[:n]+1j*fixed[n:])
        if error(a@q.conj()-q)>tol:raise ValueError('failed to fix antiunitary basis gauge')
    partner=t@q.conj();full=np.column_stack((q,partner))
    orth=error(full.conj().T@full-np.eye(dim))
    if orth>tol:raise ValueError(f'local full spinor basis is not unitary: {orth:g}')
    return q,{'reduced_real':real,'orthogonality_error':orth,'sector_error':error(r@q-1j*q)}

def _validated(structure,metadata):
    selected=metadata.get('selected',metadata);ops=selected.get('operations',[])
    if not selected.get('valid') or not ops:raise ValueError('no validated symmetry operations')
    # Recompute geometry and mappings against THIS input instead of trusting
    # an earlier valid flag or silently accepting stale atom permutations.
    check_structure=dict(structure,spin_mode='none')
    tol=float(selected.get('symprec_Ang',0))
    if not np.isfinite(tol) or tol<=0:raise ValueError('positive recorded geometry tolerance required')
    checked=detector._build_operation_set(check_structure,[o['W'] for o in ops],[o['t'] for o in ops],tol,1.0,source='SOC revalidation')
    if not checked.get('valid') or not checked.get('orbital_representation',{}).get('checked'):
        raise ValueError(f'operation/geometry/AO validation failed: {checked.get("error",checked)}')
    original={nsoc._key(o['W'],o['t']):o for o in ops}
    for o in checked['operations']:
        old=original[nsoc._key(o['W'],o['t'])]
        if old.get('atom_permutation')!=o['atom_permutation']:
            raise ValueError('stale or malformed atom permutation in supplied metadata')
        if old.get('image_shift')!=o['image_shift']:
            raise ValueError('stale or malformed periodic image shifts in supplied metadata')
    return checked

def _subgroup(checked,requested):
    ops=checked['operations'];_,identity,multiply,_,_=nsoc._group_tables(ops)
    candidates=[]
    for i,op in enumerate(ops):
        try:axis=c2_axis(op['R_cart_orthogonal'])
        except ValueError:continue
        if multiply(i,i)==identity:candidates.append((i,axis))
    if not candidates:raise ValueError('no supported proper C2 subgroup; no silent identity fallback')
    # Prefer a laboratory-z-aligned axis only to make output deterministic;
    # all construction and checks use the actual arbitrary Cartesian axis.
    candidates.sort(key=lambda v:(-abs(v[1][2]),-abs(v[1][0]),v[0]))
    pair=None
    for i,a in candidates:
        for j,b in candidates:
            if i!=j and abs(np.dot(a,b))<1e-9 and multiply(i,j)==multiply(j,i):
                pair=(i,j);break
        if pair:break
    if requested=='D2' and pair is None:raise ValueError('requested D2 has no perpendicular commuting proper C2 pair')
    if requested!='C2' and pair is not None:
        i,j=pair;ids=[identity,i,j,multiply(i,j)];kind='D2'
        if len(set(ids))!=4:raise ValueError('D2 subgroup does not have four distinct elements')
    else:i=candidates[0][0];j=None;ids=[identity,i];kind='C2'
    selected=[ops[k] for k in ids]
    nsoc._group_tables(selected)
    return kind,selected,1,(2 if j is not None else None)

def _build_c2_d2(input_path,symmetry_json,output,provenance=None,group='auto'):
    input_path=Path(input_path);symmetry_json=Path(symmetry_json);output=Path(output)
    structure=detector.parse_openmx(input_path)
    # Geometry basis itself is spin independent; input may be its NSOC parent
    # or nonmagnetic NC SOC child. Native runtime must verify occupations/TR.
    if structure['spin_mode']=='collinear':raise ValueError('collinear magnetism is outside this nonmagnetic contract')
    checked=_validated(structure,json.loads(symmetry_json.read_text()))
    kind,ops,primary,perp=_subgroup(checked,group)
    shells,offsets,n=nsoc._atom_offsets(structure);perms=[np.asarray(o['atom_permutation']) for o in ops]
    provider=detector._load_orbital_rotation_provider()
    if provider is None:raise RuntimeError('OpenMX AO rotation provider unavailable')
    spin=[np.eye(2,dtype=complex) if k==0 else spin_c2(c2_axis(o['R_cart_orthogonal'])) for k,o in enumerate(ops)]
    rotations={}
    for species,ls in shells.items():
        for l,nr in ls:
            for g,op in enumerate(ops):
                value=np.asarray(provider(l,np.asarray(op['R_cart_orthogonal'])),dtype=complex)
                if error(value.imag)>1e-10:raise ValueError('orbital rotation must be real in OpenMX real-harmonic convention')
                rotations[species,l,g]=value.real
    cols=[];seen=set();checks=[];cache={};max_local=0
    for atom in range(structure['natoms']):
        if atom in seen:continue
        orbit=nsoc._orbit(atom,perms);seen.update(orbit);where={a:i for i,a in enumerate(orbit)}
        species=structure['species'][atom];base=0
        for l,nr in shells[species]:
            d=2*'spdf'.index(l)+1;size=d*len(orbit);max_local=max(max_local,2*size)
            # Identical local representation layouts may share local algebra,
            # while each orbit retains its distinct global row indices.
            key=(species,l,tuple(tuple(where[int(p[a])] for a in orbit) for p in perms))
            if key not in cache:
                reps=[]
                for g in range(len(ops)):
                    orbital=np.zeros((size,size))
                    for a in orbit:
                        b=int(perms[g][a]);i=where[b]*d;j=where[a]*d
                        orbital[i:i+d,j:j+d]=rotations[species,l,g]
                    reps.append(np.kron(spin[g],orbital))
                q,info=local_sector_basis(reps[primary],None if perp is None else reps[perp])
                # Check every operation in the double cover up to its central
                # sign. Actual group matrices can differ by -I after a lift.
                for a in range(len(ops)):
                    for b in range(len(ops)):
                        product=np.asarray(ops[a]['R_cart_orthogonal'])@np.asarray(ops[b]['R_cart_orthogonal'])
                        matches=[k for k,o in enumerate(ops) if error(product-np.asarray(o['R_cart_orthogonal']))<1e-9]
                        if len(matches)!=1:raise ValueError('ambiguous spatial subgroup product')
                        v=reps[a]@reps[b];w=reps[matches[0]]
                        if min(error(v-w),error(v+w))>2e-10:raise ValueError('double-group representation product failed')
                cache[key]=(q,info)
            q,info=cache[key];checks.append(info)
            for radial in range(nr):
                scalar=np.array([offsets[a]+base+radial*d+mu for a in orbit for mu in range(d)])
                rows=np.r_[scalar,n+scalar]
                for j in range(size):
                    mask=abs(q[:,j])>1e-13
                    cols.append((rows[mask],q[mask,j]))
            base+=nr*d
    if len(cols)!=n:raise ValueError(f'constructed {len(cols)} columns, expected N={n}')
    full=cols+[(np.where(row<n,row+n,row-n),np.where(row<n,-val.conj(),val.conj())) for row,val in cols]
    ptr=np.r_[0,np.cumsum([len(row) for row,val in full])].astype('<i8')
    rows=np.concatenate([row for row,val in full]).astype('<i4')
    vals=np.concatenate([val for row,val in full]).astype('<c16')
    # Enforce sorted row indices in every column for deterministic consumers.
    for j in range(2*n):
        lo,hi=ptr[j:j+2];order=np.argsort(rows[lo:hi]);rows[lo:hi]=rows[lo:hi][order];vals[lo:hi]=vals[lo:hi][order]
    output.parent.mkdir(parents=True,exist_ok=True)
    species_order=list(structure['species_basis']);species_idx=[species_order.index(s) for s in structure['species']]
    dtype=1 if kind=='D2' else 2
    with output.open('wb') as f:
        f.write(struct.pack('<6iq',MAGIC,1,structure['natoms'],n,2,len(ops),len(vals)))
        for value in ([n,n],[dtype,dtype],[-1,0]):f.write(np.asarray(value,dtype='<i4').tobytes())
        f.write(np.asarray(structure['cell'],dtype='<f8').tobytes())
        f.write(np.asarray(np.diff(offsets),dtype='<i4').tobytes());f.write(np.asarray(species_idx,dtype='<i4').tobytes())
        for op in ops:
            f.write(np.asarray(op['R_cart_orthogonal'],dtype='<f8').tobytes());f.write(np.asarray(op['t'],dtype='<f8').tobytes());f.write(np.asarray(op['atom_permutation'],dtype='<i4').tobytes())
        f.write(ptr.tobytes());f.write(rows.tobytes());f.write(vals.tobytes())
    info={'schema':'openmx-generic-soc-basis-v1','binary_magic':MAGIC,'input':str(input_path.resolve()),'input_sha256':hashlib.sha256(input_path.read_bytes()).hexdigest(),
          'symmetry_metadata':str(symmetry_json.resolve()),'symmetry_metadata_sha256':hashlib.sha256(symmetry_json.read_bytes()).hexdigest(),
          'basis':str(output.resolve()),'basis_sha256':hashlib.sha256(output.read_bytes()).hexdigest(),'natoms':structure['natoms'],'N_scalar':n,'N_spinor':2*n,
          'selected_subgroup':kind,'source_point_group':checked['point_group'],'full_point_group_reduction':len(ops)==len(checked['operations']),
          'sector_sizes':[n,n],'dtype':[dtype,dtype],'partner':[-1,0],'nnz':len(vals),'operations':ops,'species_order':species_order,
          'spin_order':'all up scalar AOs, then all down','time_reversal':'J K, J=i sigma_y; U_minus=J U_plus*; partner reduced coefficients=conjugate source',
          'real_structure':None if perp is None else 'C2_perpendicular*T fixes U_plus','primary_operation':primary,'perpendicular_operation':perp,
          'validation':{'geometry_revalidated':True,'scalar_ao_representation_revalidated':True,'max_geometry_residual_Ang':checked['max_geometry_residual_Ang'],'symprec_Ang':checked['symprec_Ang'],
                        'spinor_double_group':True,'time_reversal_squared':-1,'real_structure_squared':None if perp is None else 1,
                        'max_local_orthogonality_error':max(c['orthogonality_error'] for c in checks),'max_local_sector_error':max(c['sector_error'] for c in checks),
                        'h_s_covariance_checked':False,'occupations_time_reversal_checked':False},
          'scope':'Gamma, fixed geometry, nonmagnetic SOC; actual H/S and TR covariance required at native runtime; no material-name cases',
          'complexity':'linear in atom count for fixed basis and subgroup; dense local blocks only','largest_local_dimension':max_local,'largest_column_support':max(len(row) for row,val in cols),
          'python':sys.executable,'python_version':sys.version.split()[0]}
    if provenance:Path(provenance).write_text(json.dumps(info,indent=2)+'\n')
    return info

def _cyclic_subgroup(checked,requested):
    ops=checked['operations'];_,identity,multiply,_,_=nsoc._group_tables(ops)
    candidates=[]
    for i,op in enumerate(ops):
        r=np.asarray(op['R_cart_orthogonal'])
        if abs(np.linalg.det(r)-1)>1e-9:continue
        angle=float(np.arccos(np.clip((np.trace(r)-1)*.5,-1,1)))
        for order in (6,4,3,2):
            if abs(angle-2*np.pi/order)>1e-9:continue
            if requested!='auto' and int(requested[1:])!=order:continue
            powers=[identity]
            for k in range(1,order):powers.append(multiply(i,powers[-1]))
            if len(set(powers))!=order or multiply(i,powers[-1])!=identity:continue
            if order==2:axis=c2_axis(r)
            else:axis=np.array([r[2,1]-r[1,2],r[0,2]-r[2,0],r[1,0]-r[0,1]])/(2*np.sin(angle))
            pair=None
            for j,b in enumerate(ops):
                try:perp=c2_axis(b['R_cart_orthogonal'])
                except ValueError:continue
                if abs(np.dot(axis,perp))>1e-9 or multiply(j,j)!=identity:continue
                if multiply(multiply(j,i),j)!=powers[-1]:continue
                ids=powers+[multiply(j,k) for k in powers]
                if len(set(ids))==2*order:pair=j;break
            if requested.startswith('D') and pair is None:continue
            if requested.startswith('C'):pair=None
            candidates.append((order,pair is not None,abs(axis[2]),abs(axis[0]),-i,powers,pair))
    if not candidates:raise ValueError(f'no supported proper cyclic/dihedral subgroup matching {requested}')
    order,_,_,_,_,powers,pair=max(candidates,key=lambda v:v[:5])
    ids=powers if pair is None else powers+[multiply(pair,k) for k in powers]
    selected=[ops[k] for k in ids];nsoc._group_tables(selected)
    return ('C' if pair is None else 'D')+str(order),selected,order,(None if pair is None else order)

def _improper_subgroup(checked,requested):
    ops=checked['operations'];_,identity,multiply,_,_=nsoc._group_tables(ops)
    kinds=['Cs','Ci'] if requested=='auto' else [requested]
    for kind in kinds:
        for i,op in enumerate(ops):
            r=np.asarray(op['R_cart_orthogonal'])
            if abs(np.linalg.det(r)+1)>1e-9 or multiply(i,i)!=identity:continue
            if kind=='Ci' and error(r+np.eye(3))>1e-9:continue
            if kind=='Cs' and (abs(np.trace(r)-1)>1e-9 or error(r@r-np.eye(3))>1e-9):continue
            return kind,[ops[identity],op],2,None
    raise ValueError(f'no actual inversion/mirror subgroup matching {requested}')

def build(input_path,symmetry_json,output,provenance=None,group='auto'):
    if group in ('C2','D2'):
        return _build_c2_d2(input_path,symmetry_json,output,provenance,group)
    if group not in ('auto','C3','D3','C4','D4','C6','D6','Ci','Cs'):raise ValueError('unsupported requested subgroup')
    input_path=Path(input_path);symmetry_json=Path(symmetry_json);output=Path(output)
    structure=detector.parse_openmx(input_path)
    if structure['spin_mode']=='collinear':raise ValueError('collinear magnetism is outside this nonmagnetic contract')
    checked=_validated(structure,json.loads(symmetry_json.read_text()))
    if group in ('Ci','Cs'):kind,ops,order,perp=_improper_subgroup(checked,group)
    else:
        try:kind,ops,order,perp=_cyclic_subgroup(checked,group)
        except ValueError:
            if group!='auto':raise
            kind,ops,order,perp=_improper_subgroup(checked,'auto')
    improper=kind in ('Ci','Cs')
    # Preserve the existing C2/D2 payload construction, including ordering,
    # for old commands and for auto selection on existing HgI2 metadata.
    if order==2 and not improper:return _build_c2_d2(input_path,symmetry_json,output,provenance,kind)
    shells,offsets,n=nsoc._atom_offsets(structure);perms=[np.asarray(o['atom_permutation']) for o in ops]
    provider=detector._load_orbital_rotation_provider()
    if provider is None:raise RuntimeError('OpenMX AO rotation provider unavailable')
    spins=[spin_rotation(np.linalg.det(o['R_cart_orthogonal'])*np.asarray(o['R_cart_orthogonal'])) if improper else spin_rotation(o['R_cart_orthogonal']) for o in ops]
    rotations={}
    for species,ls in shells.items():
        for l,nr in ls:
            for g,op in enumerate(ops):
                value=np.asarray(provider(l,np.asarray(op['R_cart_orthogonal'])),dtype=complex)
                if error(value.imag)>1e-10:raise ValueError('orbital action must be real in OpenMX convention')
                rotations[species,l,g]=value.real
    sector_columns=[[] for k in range(order)];seen=set();cache={};checks=[];max_local=0;layout=None
    for atom in range(structure['natoms']):
        if atom in seen:continue
        orbit=nsoc._orbit(atom,perms);seen.update(orbit);where={a:i for i,a in enumerate(orbit)}
        species=structure['species'][atom];base=0
        for l,nr in shells[species]:
            d=2*'spdf'.index(l)+1;size=d*len(orbit);max_local=max(max_local,2*size)
            key=(species,l,tuple(tuple(where[int(p[a])] for a in orbit) for p in perms))
            if key not in cache:
                reps=[]
                for g in range(len(ops)):
                    orbital=np.zeros((size,size))
                    for a in orbit:
                        b=int(perms[g][a]);i=where[b]*d;j=where[a]*d
                        orbital[i:i+d,j:j+d]=rotations[species,l,g]
                    reps.append(np.kron(spins[g],orbital))
                cache[key]=local_improper_basis(reps[1],kind) if improper else local_cyclic_basis(reps[1],order,None if perp is None else reps[perp])
            layout=cache[key];checks.append(layout['orthogonality_error'])
            for radial in range(nr):
                scalar=np.array([offsets[a]+base+radial*d+mu for a in orbit for mu in range(d)]);rows=np.r_[scalar,n+scalar]
                for k,q in enumerate(layout['basis']):
                    for j in range(q.shape[1]):
                        mask=abs(q[:,j])>1e-13;sector_columns[k].append((rows[mask],q[mask,j]))
            base+=nr*d
    keep=[k for k,cols in enumerate(sector_columns) if cols];remap={k:i for i,k in enumerate(keep)}
    sizes=[len(sector_columns[k]) for k in keep];dtype=[layout['dtype'][k] for k in keep]
    partners=[-1 if layout['partner'][k]<0 else remap[layout['partner'][k]] for k in keep]
    phases=[layout['phase'][k] for k in keep]
    if sum(sizes)!=2*n or len(sizes)>48:raise ValueError('invalid total spinor sector dimension')
    columns=[entry for k in keep for entry in sector_columns[k]]
    ptr=np.r_[0,np.cumsum([len(row) for row,val in columns])].astype('<i8')
    rows=np.concatenate([row for row,val in columns]).astype('<i4');vals=np.concatenate([val for row,val in columns]).astype('<c16')
    for j in range(2*n):
        lo,hi=ptr[j:j+2];idx=np.argsort(rows[lo:hi]);rows[lo:hi]=rows[lo:hi][idx];vals[lo:hi]=vals[lo:hi][idx]
    output.parent.mkdir(parents=True,exist_ok=True);species_order=list(structure['species_basis'])
    with output.open('wb') as f:
        f.write(struct.pack('<6iq',MAGIC,2 if improper else 1,structure['natoms'],n,len(sizes),len(ops),len(vals)))
        for value in (sizes,dtype,partners):f.write(np.asarray(value,dtype='<i4').tobytes())
        f.write(np.asarray(structure['cell'],dtype='<f8').tobytes());f.write(np.asarray(np.diff(offsets),dtype='<i4').tobytes())
        f.write(np.asarray([species_order.index(s) for s in structure['species']],dtype='<i4').tobytes())
        for op in ops:
            f.write(np.asarray(op['R_cart_orthogonal'],dtype='<f8').tobytes());f.write(np.asarray(op['t'],dtype='<f8').tobytes());f.write(np.asarray(op['atom_permutation'],dtype='<i4').tobytes())
        f.write(ptr.tobytes());f.write(rows.tobytes());f.write(vals.tobytes())
    info={'schema':'openmx-generic-soc-basis-v2' if improper else 'openmx-generic-soc-basis-v1','binary_version':2 if improper else 1,'binary_magic':MAGIC,'input':str(input_path.resolve()),'input_sha256':hashlib.sha256(input_path.read_bytes()).hexdigest(),
          'symmetry_metadata':str(symmetry_json.resolve()),'symmetry_metadata_sha256':hashlib.sha256(symmetry_json.read_bytes()).hexdigest(),
          'basis':str(output.resolve()),'basis_sha256':hashlib.sha256(output.read_bytes()).hexdigest(),'natoms':structure['natoms'],'N_scalar':n,'N_spinor':2*n,
          'selected_subgroup':kind,'cyclic_order':None if improper else order,'source_point_group':checked['point_group'],'full_point_group_reduction':len(ops)==len(checked['operations']),
          'sector_sizes':sizes,'dtype':dtype,'partner':partners,'sector_phases':[[float(z.real),float(z.imag)] for z in phases],
          'nnz':len(vals),'operations':ops,'species_order':species_order,'primary_operation':1,'perpendicular_operation':perp,
          'spin_order':'all up scalar AOs, then all down','spin_rotation_convention':'axial SU2(det(R)*R), AO uses physical polar R' if improper else 'proper SU2(R), AO uses physical R','time_reversal':'J K, J=i sigma_y; conjugate distinct phase partner=J U_source*; self-conjugate phase solved in full',
          'real_structure':None if perp is None else 'C2_perpendicular*T fixes every sector basis; self-conjugate sectors are NOT halved',
          'validation':{'geometry_revalidated':True,'scalar_ao_representation_revalidated':True,'max_geometry_residual_Ang':checked['max_geometry_residual_Ang'],'symprec_Ang':checked['symprec_Ang'],
                        'spinor_lift_power':1 if kind=='Ci' else -1,'time_reversal_squared':-1,'real_structure_squared':None if perp is None else 1,
                        'max_local_orthogonality_error':max(checks),'h_s_covariance_checked':False,'occupations_time_reversal_checked':False},
          'scope':'Gamma fixed geometry nonmagnetic SOC; explicit operation-driven subgroup, not complete arbitrary point-group irreps; native H/S and TR covariance checks required',
          'complexity':'linear in atom count at fixed subgroup and basis; no dense global U','largest_local_dimension':max_local,
          'python':sys.executable,'python_version':sys.version.split()[0]}
    if provenance:Path(provenance).write_text(json.dumps(info,indent=2)+'\n')
    return info

DEFAULT_SYMPREC_LADDER=(1e-6,3e-6,1e-5,2e-5,3e-5,5e-5,1e-4)

def prepare(input_path,output,symmetry_json=None,provenance=None,group='auto',
            symprec_ladder=DEFAULT_SYMPREC_LADDER,geometry='strict'):
    """One-command metadata preparation only; never enable native runtime.

    Strict means keeping the supplied coordinates: no geometry projection.
    The detector's strict_residual is its existing status classification,
    while actual geometric tolerances and residuals are recorded explicitly.
    """
    if geometry!='strict':raise ValueError('only strict, unmodified input geometry is supported')
    ladder=np.asarray(symprec_ladder,dtype=float)
    if ladder.ndim!=1 or len(ladder)<2 or not np.all(np.isfinite(ladder)) or np.any(ladder<=0) or np.any(np.diff(ladder)<=0):
        raise ValueError('symprec ladder must contain at least two finite positive strictly increasing tolerances')
    input_path=Path(input_path).resolve();output=Path(output).resolve()
    provenance=Path(provenance).resolve() if provenance else output.with_suffix('.json')
    automatic=symmetry_json is None
    symmetry_json=output.with_suffix('.symmetry.json') if automatic else Path(symmetry_json).resolve()
    paths=[output,provenance]+([symmetry_json] if automatic else [])
    if len(set(paths))!=len(paths) or input_path in paths or (not automatic and symmetry_json in paths):
        raise ValueError('input, basis, provenance and detection output paths must be distinct')
    for path in paths:
        if path.exists():raise FileExistsError(f'output already exists: {path}; choose a new result path')
    if automatic:
        metadata=detector.run(input_path,symmetry_json,ladder.tolist(),5e-6,1.0,mode='auto')
        if metadata.get('stable_pair') is None or not metadata['selected'].get('valid'):
            raise ValueError(f'no stable validated geometry symmetry platform; inspect {symmetry_json}')
    info=build(input_path,symmetry_json,output,None,group)
    info['preparation']={'auto_detected':automatic,'geometry_policy':'strict_unmodified_coordinates',
                         'symmetry_metadata':str(symmetry_json),'production_input_modified':False,
                         'native_runtime_enabled':False,'scf_Symmetry_Auto_parser_connected':False}
    provenance.parent.mkdir(parents=True,exist_ok=True);provenance.write_text(json.dumps(info,indent=2)+'\n')
    return info

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',type=Path,required=True);p.add_argument('--symmetry-json',type=Path,
                   help='validated metadata; if omitted, detect automatically beside output')
    p.add_argument('--output',type=Path,required=True);p.add_argument('--provenance',type=Path);p.add_argument('--group',choices=['auto','C2','D2','C3','D3','C4','D4','C6','D6','Ci','Cs'],default='auto')
    p.add_argument('--symprec-ladder',type=float,nargs='+',default=DEFAULT_SYMPREC_LADDER,
                   help='strictly increasing geometry tolerance scan in Angstrom')
    p.add_argument('--geometry',choices=['strict'],default='strict',
                   help='verify supplied coordinates; never project or modify input geometry')
    a=p.parse_args()
    try:info=prepare(a.input,a.output,a.symmetry_json,a.provenance,a.group,a.symprec_ladder,a.geometry)
    except (ValueError,OSError,RuntimeError) as exc:p.error(str(exc))
    print(json.dumps({k:info[k] for k in ['selected_subgroup','N_scalar','N_spinor','sector_sizes','dtype','partner','nnz','basis','validation']},indent=2))
if __name__=='__main__':main()
