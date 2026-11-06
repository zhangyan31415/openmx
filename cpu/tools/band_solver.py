#!/usr/bin/env python3
"""Internal, standalone SLEPc worker for native OpenMX band dispersion.

No moirekp dependency. Native code supplies owned real-space blocks once; the
worker constructs matching sparse Bloch H/S and returns only requested states.
"""
from pathlib import Path
import argparse
import json
import os
import struct
import sys
import time

EV_PER_HA = 27.211386245988


def reference_ha(value, fermi):
    import numpy as np
    result = float(fermi) if value.lower() == 'fermi' else float(value)/EV_PER_HA
    if not np.isfinite(result):raise ValueError('Band.ReferenceEnergy must be Fermi or a finite absolute energy in eV')
    return result


def nearest_indices(values, target, count):
    import numpy as np
    if count < 1 or count > len(values):raise ValueError('insufficient converged eigenvalues')
    indices=np.lexsort((values,abs(values-target)))[:count]
    return indices[np.argsort(values[indices],kind='stable')]


def load_parts(directory, parts):
    import numpy as np
    rr=[];cc=[];translations=[];values=[];schema=None
    for rank in range(parts):
        path=Path(directory)/f'part-{rank:06d}.bin'
        with path.open('rb') as stream:
            if stream.read(8)!=b'OMBAND01':raise ValueError('invalid native band snapshot magic')
            header=np.fromfile(stream,'<i4',5)
            if len(header)!=5:raise ValueError('short band snapshot header')
            n,spin,blocks,fields,version=map(int,header)
            if n<1 or spin not in (0,1,3) or blocks<0 or fields!={0:2,1:3,3:8}[spin] or version!=1:
                raise ValueError('invalid native band snapshot schema')
            if schema is not None and schema!=(n,spin,fields):raise ValueError('inconsistent native band snapshots')
            schema=n,spin,fields
            for block in range(blocks):
                h=np.fromfile(stream,'<i4',7)
                if len(h)!=7:raise ValueError('short block header')
                row,col,na,nb,*R=map(int,h)
                if min(row,col)<0 or min(na,nb)<1 or row+na>n or col+nb>n:raise ValueError('block outside native AO basis')
                data=np.fromfile(stream,'<f8',na*nb*fields)
                if len(data)!=na*nb*fields or not np.isfinite(data).all():raise ValueError('invalid block values')
                rr.append(np.repeat(np.arange(row,row+na,dtype=np.int32),nb))
                cc.append(np.tile(np.arange(col,col+nb,dtype=np.int32),na))
                translations.append(np.tile(np.array(R,dtype=np.int32),(na*nb,1)))
                values.append(data.reshape(fields,na*nb))
            if stream.read(1):raise ValueError('trailing native band snapshot payload')
    if not values:raise ValueError('no native H/S blocks')
    return dict(n=schema[0],spin=schema[1],rows=np.concatenate(rr),cols=np.concatenate(cc),
                R=np.concatenate(translations),values=np.concatenate(values,axis=1))


def bloch_matrices(data, k, channel):
    import numpy as np
    import scipy.sparse as sp
    n=data['n'];v=data['values'];phase=np.exp(2j*np.pi*(data['R']@k))
    def matrix(values):
        a=sp.coo_matrix((values*phase,(data['rows'],data['cols'])),shape=(n,n)).tocsr()
        a.sum_duplicates();a.eliminate_zeros();return a
    S=matrix(v[-1])
    if data['spin']!=3:return matrix(v[channel]),S
    up=matrix(v[0]+1j*v[4]);down=matrix(v[1]+1j*v[5]);ud=matrix(v[2]+1j*(v[3]+v[6]))
    return sp.bmat([[up,ud],[ud.getH(),down]],format='csr'),sp.block_diag((S,S),format='csr')


def native_hermitian(matrix):
    import numpy as np
    import scipy.sparse as sp
    from scipy.sparse.linalg import norm
    if not np.isfinite(matrix.data).all():raise ValueError('nonfinite Bloch matrix')
    defect=float(norm(matrix-matrix.getH())/max(float(norm(matrix)),1e-300))
    if not np.isfinite(defect) or defect>5e-7:raise ValueError('native Bloch matrix is not Hermitian: '+str(defect))
    # The original zheevx UPLO=L reads the lower triangle, including real diag.
    lower=sp.tril(matrix,k=-1,format='csr')
    return (lower+lower.getH()+sp.diags(matrix.diagonal().real)).tocsr(),defect


def libraries():
    import numpy as np
    from petsc4py import PETSc
    from slepc4py import SLEPc
    if np.dtype(PETSc.ScalarType).kind!='c':raise ValueError('targeted bands require complex PETSc/SLEPc')
    factor=next((name for name in ('mkl_pardiso','mumps','superlu') if PETSc.Sys.hasExternalPackage(name)),None)
    if factor is None:raise ValueError('PETSc needs a sequential sparse LU backend (PARDISO, MUMPS or SuperLU)')
    return PETSc,SLEPc,factor


def solve(H,S,target,count,tolerance=1e-9):
    import numpy as np
    from scipy.sparse.linalg import eigsh
    PETSc,SLEPc,factor=libraries()
    H,hdef=native_hermitian(H);S,sdef=native_hermitian(S);n=H.shape[0]
    if not 0<count<n-1:raise ValueError('Band.NumBands must be smaller than matrix dimension minus one')
    start=time.monotonic()
    # Reject bases for which stock OpenMX would invoke its overlap penalty.
    # This is one extremal overlap eigenvalue, not a full S eigendecomposition.
    smin=float(eigsh(S,k=1,which='SA',return_eigenvectors=False,tol=1e-7)[0])
    if not np.isfinite(smin) or smin<=1e-9:raise ValueError('overlap requires native regularization; targeted bands reject this basis')
    overlap_seconds=time.monotonic()-start
    def petsc_matrix(a):
        if a.nnz>np.iinfo(PETSc.IntType).max:raise ValueError('sparse index capacity exceeded')
        m=PETSc.Mat().createAIJ(size=a.shape,csr=(a.indptr.astype(PETSc.IntType),a.indices.astype(PETSc.IntType),a.data.astype(PETSc.ScalarType)),comm=PETSc.COMM_SELF)
        m.assemble();m.setOption(PETSc.Mat.Option.HERMITIAN,True);return m
    A=B=eps=vr=vi=None
    try:
        A=petsc_matrix(H);B=petsc_matrix(S);B.setOption(PETSc.Mat.Option.SPD,True)
        eps=SLEPc.EPS().create(PETSc.COMM_SELF);eps.setOperators(A,B);eps.setProblemType(SLEPc.EPS.ProblemType.GHEP)
        eps.setType('krylovschur');wanted=min(n-2,count+max(4,count//4))
        eps.setDimensions(wanted,min(n,max(2*wanted+12,40)))
        eps.setWhichEigenpairs(SLEPc.EPS.Which.TARGET_REAL);eps.setTarget(target)
        st=eps.getST();st.setType('sinvert');st.setShift(target)
        ksp=st.getKSP();ksp.setType('preonly');pc=ksp.getPC();pc.setType('lu');pc.setFactorSolverType(factor)
        eps.setTolerances(tolerance,10000)
        solve_start=time.monotonic();eps.solve();solve_seconds=time.monotonic()-solve_start
        converged=eps.getConverged()
        if eps.getConvergedReason()<=0 or converged<count:raise ValueError('SLEPc did not converge all requested bands')
        values=np.array([complex(eps.getEigenvalue(i)) for i in range(converged)])
        if not np.isfinite(values).all() or max(abs(values.imag))>1e-10:raise ValueError('nonreal/nonfinite eigenvalues')
        selected=nearest_indices(values.real,target,count);vectors=[];residual=0.
        vr=A.createVecRight();vi=A.createVecRight()
        for index in selected:
            eps.getEigenpair(int(index),vr,vi);x=vr.getArray(readonly=True).copy();e=values[index].real
            hx=H@x;sx=S@x;den=np.linalg.norm(hx)+abs(e)*np.linalg.norm(sx)
            residual=max(residual,float(np.linalg.norm(hx-e*sx)/max(den,1e-300)))
            vectors.append(x)
        C=np.column_stack(vectors);orth=float(abs(C.conj().T@(S@C)-np.eye(count)).max())
        if not np.isfinite(residual+orth) or residual>1e-8 or orth>1e-8:raise ValueError(f'generalized eigenpair audit failed: residual={residual} Sorth={orth}')
        return values.real[selected],dict(residual=residual,Sorth=orth,overlap_min=smin,H_hermiticity=hdef,S_hermiticity=sdef,
            solver='SLEPc KrylovSchur shift-invert',factor=factor,PETSc=PETSc.Sys.getVersion(),SLEPc=SLEPc.Sys.getVersion(),
            matrix_dimension=n,requested=count,converged=converged,overlap_check_seconds=overlap_seconds,solve_seconds=solve_seconds,
            S_diagonal_shift=0.,native_triangle='lower',full_spectrum=False)
    finally:
        for value in (vi,vr,eps,B,A):
            if value is not None:value.destroy()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check',action='store_true')
    parser.add_argument('--directory',type=Path);parser.add_argument('--rank',type=int,default=0)
    parser.add_argument('--workers',type=int,default=1);parser.add_argument('--parts',type=int,default=1)
    parser.add_argument('--count',type=int);parser.add_argument('--target-ha',type=float)
    parser.add_argument('--threads',type=int,default=1)
    args=parser.parse_args()
    for key in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
        os.environ[key]=str(args.threads)
    if args.check:
        libraries();return 0
    import numpy as np
    started=time.monotonic();data=load_parts(args.directory,args.parts)
    with (args.directory/'kpoints.bin').open('rb') as stream:
        header=np.fromfile(stream,'<i4',1)
        if len(header)!=1 or header[0]<1:raise ValueError('bad native kpoint count')
        points=np.fromfile(stream,'<f8',int(header[0])*3).reshape(-1,3)
        if stream.read(1) or not np.isfinite(points).all():raise ValueError('bad native kpoint payload')
    audits=[];channels=2 if data['spin']==1 else 1
    final=args.directory/f'values-{args.rank:06d}.bin';temporary=final.with_suffix('.writing')
    with temporary.open('xb') as out:
        for index in range(args.rank,len(points),args.workers):
            for channel in range(channels):
                before=time.monotonic();H,S=bloch_matrices(data,points[index],channel)
                values,audit=solve(H,S,args.target_ha,args.count)
                out.write(struct.pack('<2i',index,channel));out.write(values.astype('<f8').tobytes())
                audits.append(dict(kpoint=index,channel=channel,k=points[index].tolist(),total_seconds=time.monotonic()-before,**audit))
    temporary.replace(final)
    (args.directory/f'audit-{args.rank:06d}.json').write_text(json.dumps(dict(status='complete',python=sys.executable,worker=args.rank,
        threads=args.threads,seconds=time.monotonic()-started,points=audits),indent=2)+'\n')
    return 0


if __name__=='__main__':
    try:raise SystemExit(main())
    except Exception as error:
        print('OPENMX_TARGET_BANDS_ERROR: '+str(error),file=sys.stderr,flush=True)
        raise SystemExit(2)
