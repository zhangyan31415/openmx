"""Real spatial OpenMX D3 basis; no spin rotations or geometry modification.

Rows are OpenMX atom/l/radial/real-harmonic order. Columns are Ex,Ey,A1,A2;
paired Ey=(R-R.T)Ex/sqrt(3). Each dense construction is restricted to one
atom orbit and one l/radial shell (at most 6*(2*l+1) rows).
"""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import os
import re
import struct
import sys

import numpy as np
from scipy import linalg, sparse
from scipy.spatial import cKDTree
import spglib

OLD_ANGULAR = Path(__file__).resolve().parents[1]/'openmx_orbitals.py'
ANGULAR_SOURCE = Path(__file__).resolve().parents[1]/'comp2real.json'


def block(text, name):
    return text.split('<'+name,1)[1].split(name+'>',1)[0].strip().splitlines()


def maxerr(matrix):
    return float(np.max(np.abs(matrix.data), initial=0))


def clean(matrix):
    matrix=matrix.tocsc()
    matrix.data[np.abs(matrix.data)<1e-13]=0
    matrix.eliminate_zeros()
    return matrix


def local_basis(rotation, c2):
    """Decompose an orthogonal real D3 representation, keeping E partners."""
    r=np.asarray(rotation,dtype=float); c=np.asarray(c2,dtype=float)
    n=len(r); eye=np.eye(n)
    if max(np.max(abs(r@r@r-eye)),np.max(abs(c@c-eye)),np.max(abs(c@r@c-r.T)))>1e-10:
        raise ValueError('Invalid real D3 representation')
    p0=(eye+r+r.T)/3
    def image(p):
        w,v=linalg.eigh((p+p.T)/2)
        return v[:,w>.5]
    ex=image((eye-p0)@(eye+c)/2)
    ey=(r-r.T)@ex/np.sqrt(3)
    a1=image(p0@(eye+c)/2)
    a2=image(p0@(eye-c)/2)
    u=np.column_stack((ex,ey,a1,a2))
    dims=[ex.shape[1],ey.shape[1],a1.shape[1],a2.shape[1]]
    if u.shape!=(n,n) or np.max(abs(u.T@u-eye),initial=0)>1e-11:
        raise ValueError('Incomplete or nonorthogonal D3 basis')
    return u,dims


def geometry(text):
    atoms=[line.split() for line in block(text,'Atoms.SpeciesAndCoordinates')]
    cell=np.array([[float(v) for v in line.split()[:3]] for line in block(text,'Atoms.UnitVectors')])
    xyz=np.array([[float(v) for v in atom[2:5]] for atom in atoms])
    # All production/smoke sources are Angstrom; reject other units explicitly.
    for key in ('Atoms.UnitVectors.Unit','Atoms.SpeciesAndCoordinates.Unit'):
        if not re.search(r'^\s*'+re.escape(key)+r'\s+Ang\b',text,re.M|re.I):
            raise ValueError(f'{key} must be Ang')
    frac=(xyz@np.linalg.inv(cell))%1
    labels={name:i+1 for i,name in enumerate(sorted({a[1] for a in atoms}))}
    types=np.array([labels[a[1]] for a in atoms])
    data=spglib.get_symmetry((cell,frac,types),symprec=1e-4)
    if data is None:
        raise ValueError('No geometry symmetry dataset')
    r=np.array([[-.5,-np.sqrt(3)/2,0],[np.sqrt(3)/2,-.5,0],[0,0,1.]])
    c3=[]; c2=[]
    for w,t in zip(data['rotations'],data['translations']):
        q=cell.T@w@np.linalg.inv(cell.T)
        if np.max(abs(q-r))<1e-6:
            c3.append((np.linalg.norm(t-np.rint(t)),w,t))
        if np.linalg.det(q)>.99 and abs(np.trace(q)+1)<1e-7 and q[2,2]<-.99:
            c2.append((np.linalg.norm(t-np.rint(t)),w,t,q))
    if not c3 or not c2:
        raise ValueError('C3 and in-plane C2 operations required')
    _,w3,t3=min(c3,key=lambda x:x[0])
    _,w2,t2,q2=min(c2,key=lambda x:x[0])
    _,ax=linalg.eigh((q2[:2,:2]+q2[:2,:2].T)/2)
    axis=np.r_[ax[:,-1],0.]
    q2=2*np.outer(axis,axis)-np.eye(3)
    def mapping(w,t):
        targets=(frac@w.T+t)%1
        perm=np.empty(len(atoms),dtype=np.int32)
        for typ in set(types):
            ids=np.flatnonzero(types==typ)
            _,j=cKDTree(frac[ids],boxsize=1).query(targets[ids])
            perm[ids]=ids[j]
        if len(set(perm))!=len(atoms):
            raise ValueError('Geometry operation not a species-preserving permutation')
        return perm
    perm3=mapping(w3,t3); perm2=mapping(w2,t2)
    ident=np.arange(len(atoms))
    if not (np.array_equal(perm3[perm3[perm3]],ident) and np.array_equal(perm2[perm2],ident)
            and np.array_equal(perm2[perm3[perm2]],np.argsort(perm3))):
        raise ValueError('Atom permutations violate D3')
    shift3=np.round(3*(t3-np.rint(t3)))/3
    shift2=np.round(6*(t2-np.rint(t2)))/6
    def residual(q,shift,perm):
        df=(xyz@q.T+shift@cell-xyz[perm])@np.linalg.inv(cell)
        df-=np.rint(df)
        return float(np.max(np.linalg.norm(df@cell,axis=1)))
    errors=[residual(r,shift3,perm3),residual(q2,shift2,perm2)]
    if max(errors)>5e-5:
        raise ValueError(f'Geometry residual too large: {errors}')
    return atoms,perm3,shift3,r,perm2,shift2,q2,axis,errors


def angular_rotations(axis):
    # Reuse the established OpenMX Comp2Real ordering, without spin factors.
    spec=importlib.util.spec_from_file_location('openmx_c3_angular',OLD_ANGULAR)
    old=importlib.util.module_from_spec(spec); spec.loader.exec_module(old)
    cc=old.angular_tables(); c3=[]; c2=[]
    for l,c in enumerate(cc):
        ms=np.arange(-l,l+1)
        lp=np.zeros((2*l+1,2*l+1),complex)
        for k,m in enumerate(ms[:-1]):
            lp[k+1,k]=np.sqrt(l*(l+1)-m*(m+1))
        jx=(lp+lp.conj().T)/2; jy=(lp-lp.conj().T)/(2j)
        mats=[c.conj()@np.diag(np.exp(-2j*np.pi*ms/3))@c.T,
              c.conj()@linalg.expm(-1j*np.pi*(axis[0]*jx+axis[1]*jy))@c.T]
        if max(np.max(abs(mat.imag)) for mat in mats)>1e-13:
            raise ValueError('Orbital representation not real')
        c3.append(mats[0].real); c2.append(mats[1].real)
    return c3,c2


def build(input_path, output_path=None, provenance_path=None):
    inp=Path(input_path).resolve(); text=inp.read_text()
    atoms,p3,t3,q3,p2,t2,q2,axis,geomerr=geometry(text)
    definitions={a[0]:a[1] for line in block(text,'Definition.of.Atomic.Species') if (a:=line.split())}
    shells={name:[('spdf'.index(ch),int(count)) for ch,count in re.findall(r'([spdf])(\d+)',value.split('-',1)[1])]
            for name,value in definitions.items()}
    counts=[sum((2*l+1)*nr for l,nr in shells[a[1]]) for a in atoms]
    offsets=np.r_[0,np.cumsum(counts)]; n=int(offsets[-1])
    d3,d2=angular_rotations(axis)
    if np.max(abs(d3[1]-q3))>1e-12 or np.max(abs(d2[1]-q2))>1e-12:
        raise ValueError('OpenMX p-shell and Cartesian rotations disagree')
    def operator(perm,angular):
        ri=[];ci=[];values=[]
        for at,atom in enumerate(atoms):
            local=0
            for l,nrad in shells[atom[1]]:
                d=2*l+1
                for radial in range(nrad):
                    i,j=np.nonzero(abs(angular[l])>1e-13)
                    ri.extend(offsets[perm[at]]+local+i)
                    ci.extend(offsets[at]+local+j)
                    values.extend(angular[l][i,j]); local+=d
        return clean(sparse.csc_matrix((values,(ri,ci)),shape=(n,n)))
    rot=operator(p3,d3); c2=operator(p2,d2)
    sectors=[[],[],[],[]]; seen=set(); orbit_sizes=[]; localdims=[]
    for at in range(len(atoms)):
        if at in seen: continue
        orbit={at}; todo=[at]
        while todo:
            a=todo.pop()
            for b in (int(p3[a]),int(p2[a])):
                if b not in orbit: orbit.add(b);todo.append(b)
        orbit=sorted(orbit); seen.update(orbit); orbit_sizes.append(len(orbit))
        local=0
        for l,nrad in shells[atoms[at][1]]:
            d=2*l+1
            for radial in range(nrad):
                ids=np.array([int(offsets[a])+local+mu for a in orbit for mu in range(d)])
                localdims.append(len(ids))
                assert len(ids)<=6*d
                rr=rot[ids,:][:,ids].toarray(); cc=c2[ids,:][:,ids].toarray()
                u,dims=local_basis(rr,cc)
                for k,(start,end) in enumerate(zip(np.r_[0,np.cumsum(dims)[:-1]],np.cumsum(dims))):
                    for col in range(start,end):
                        mask=abs(u[:,col])>1e-13
                        sectors[k].append((ids[mask],u[mask,col]))
                local+=d
    rows=[]; values=[]; ptr=[0]
    for sector in sectors:
        for ri,va in sector:
            rows.extend(ri);values.extend(va);ptr.append(len(rows))
    u=clean(sparse.csc_matrix((values,rows,ptr),shape=(n,n)))
    dims=[len(s) for s in sectors]; eye=sparse.eye(n,format='csc')
    e=dims[0];bounds=np.r_[0,np.cumsum(dims)]
    ex=u[:,:e];ey=u[:,e:2*e]
    errors={'orthogonality':maxerr(u.T@u-eye),'C3_cube':maxerr(rot@rot@rot-eye),
            'C2_square':maxerr(c2@c2-eye),'D3_relation':maxerr(c2@rot@c2-rot.T),
            'E_pair':maxerr(ey-(rot-rot.T)@ex/np.sqrt(3)),
            'C2_Ex':maxerr(c2@ex-ex),'C2_Ey':maxerr(c2@ey+ey),
            'C3_A':maxerr(rot@u[:,2*e:]-u[:,2*e:]),
            'C2_A1':maxerr(c2@u[:,bounds[2]:bounds[3]]-u[:,bounds[2]:bounds[3]]),
            'C2_A2':maxerr(c2@u[:,bounds[3]:]+u[:,bounds[3]:])}
    if max(errors.values())>1e-11: raise ValueError(errors)
    info={'source_input':str(inp),'source_input_sha256':hashlib.sha256(inp.read_bytes()).hexdigest(),
          'angular_source':str(ANGULAR_SOURCE),'angular_source_sha256':hashlib.sha256(ANGULAR_SOURCE.read_bytes()).hexdigest(),
          'generator_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
          'natoms':len(atoms),'N':n,'dimensions':dims,'distinct_solve_dimensions':[dims[0],*dims[2:]],
          'nnz':u.nnz,'basis_errors':errors,'max_local_block_dimension':max(localdims),
          'local_block_count':len(localdims),'atom_orbit_sizes':orbit_sizes,
          'C3_geometry_max_error_ang':geomerr[0],'C2_geometry_max_error_ang':geomerr[1],
          'C3_shift_fractional':t3.tolist(),'C2_shift_fractional':t2.tolist(),
          'C3_rotation':q3.tolist(),'C2_rotation':q2.tolist(),'C2_axis':axis.tolist(),
          'environment_prefix':sys.prefix,'host':os.uname().nodename,
          'requested_input':str(Path(input_path).absolute()),'working_directory':os.getcwd(),
          'python_executable':sys.executable,'python_version':sys.version,
          'numpy_version':np.__version__,'spglib_version':spglib.__version__,
          'command_argv':sys.argv,
          'angular_helper':str(OLD_ANGULAR),'angular_helper_sha256':hashlib.sha256(OLD_ANGULAR.read_bytes()).hexdigest(),
          'final_process_state':'basis generated and all algebraic checks passed',
          'process_threads':'one process, BLAS/OMP threads 1 recommended',
          'representation':'real orbital NSOC; no spin factors; atom/l/radial/OpenMX Comp2Real row order',
          'column_order':'Ex,Ey,A1,A2; Ey=(R-R.T)Ex/sqrt(3); C2Ex=Ex,C2Ey=-Ey',
          'group_action':'active Cartesian orbital rotation with source atom mapped to perm[source]; translations fractional row vectors times cell',
          'geometry_modified':False,'calculation_category':'basis construction and mathematical validation; not an SCF result',
          'smoke_source_spin_fields_ignored':bool(re.search(r'^scf.SpinPolarization\s+NC\b',text,re.M|re.I))}
    if output_path is not None:
        out=Path(output_path).resolve();out.parent.mkdir(parents=True,exist_ok=True)
        with out.open('wb') as f:
            f.write(struct.pack('<9i',0x4e534433,1,len(atoms),n,4,*dims));f.write(struct.pack('<q',u.nnz))
            for data,dtype in [(p3,'<i4'),(t3,'<f8'),(p2,'<i4'),(q2,'<f8'),(t2,'<f8'),
                               (u.indptr,'<i8'),(u.indices,'<i4'),(u.data,'<f8')]:
                f.write(np.asarray(data,dtype=dtype).tobytes())
        info.update(payload=str(out),payload_sha256=hashlib.sha256(out.read_bytes()).hexdigest(),payload_bytes=out.stat().st_size)
    if provenance_path is not None:
        path=Path(provenance_path);path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(json.dumps(info,indent=2)+'\n')
    return {'U':u,'R':rot,'C2':c2,'sizes':dims,'metadata':info}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--provenance',type=Path)
    args=parser.parse_args()
    result=build(args.input,args.output,args.provenance)
    print(json.dumps(result['metadata'],indent=2))
