"""Preserved sparse spinor gauge from the verified legacy C3/D3 generator.

Only angular-table lookup was replaced by the repository-local table provider.
No full-system dense basis is constructed.
"""
from pathlib import Path
import importlib.util
import re
import numpy as np
from scipy import sparse,linalg
_spec=importlib.util.spec_from_file_location('legacy_openmx_orbitals',Path(__file__).resolve().parents[1]/'openmx_orbitals.py')
_orbitals=importlib.util.module_from_spec(_spec);_spec.loader.exec_module(_orbitals)
Q3=np.array([[-.5,-np.sqrt(3)/2,0],[np.sqrt(3)/2,-.5,0],[0,0,1.]])
U3=np.diag(np.exp([-1j*np.pi/3,1j*np.pi/3]))
def require(value,message):
    if not value:raise ValueError(message)
def clean(x):
    x=x.tocsc();x.data[abs(x.data)<1e-13]=0;x.eliminate_zeros();x.sort_indices();return x
def maxerr(x):return float(max(abs(x.data),default=0))
def angular():return _orbitals.angular_tables()
def angular_rotations(q):return [_orbitals.get_any_rot_orb_twostep(s,q)for s in 'spdf']
def shell_list(definition):return [('spdf'.index(ch),int(count)) for ch,count in re.findall(r'([spdf])(\d+)',definition.split('-',1)[1])]

def orbital_block(definition,q):
    d=angular_rotations(q)
    return linalg.block_diag(*[d[l] for l,nrad in shell_list(definition) for _ in range(nrad)])

def action(species,defs,perm,q,spin):
    shells={s:shell_list(defs[s]) for s in set(species)}
    counts=[sum((2*l+1)*n for l,n in shells[s]) for s in species];offset=np.r_[0,np.cumsum(counts)];n=int(offset[-1]);ds=angular_rotations(q)
    rows=[];cols=[];values=[]
    for a,s in enumerate(species):
        require(species[int(perm[a])]==s,'operation changes species');local=0
        for l,nrad in shells[s]:
            for _ in range(nrad):
                for so in range(2):
                    for si in range(2):
                        for i in range(2*l+1):
                            for j in range(2*l+1):
                                z=spin[so,si]*ds[l][i,j]
                                if abs(z)>1e-13:rows.append(so*n+int(offset[perm[a]])+local+i);cols.append(si*n+int(offset[a])+local+j);values.append(z)
                local+=2*l+1
    return clean(sparse.csc_matrix((values,(rows,cols)),shape=(2*n,2*n)))

def build_spinor_basis(species,defs,perm,c2perm=None,q2=None):
    perm=np.asarray(perm,dtype=int);nat=len(species);require(sorted(perm.tolist())==list(range(nat)),'C3 permutation is not bijective')
    require(np.array_equal(perm[perm[perm]],np.arange(nat)),'C3 permutation order is not three')
    shells={s:shell_list(defs[s]) for s in set(species)}
    counts=[sum((2*l+1)*n for l,n in shells[s]) for s in species];offset=np.r_[0,np.cumsum(counts)];n=int(offset[-1]);cc=angular();cols=[[],[],[]];seen=set()
    for a in range(nat):
        if a in seen:continue
        cycle=[a];b=int(perm[a])
        while b!=a:cycle.append(b);b=int(perm[b]);require(len(cycle)<=3,'invalid C3 cycle')
        seen.update(cycle)
        for spin in range(2):
            local=0
            for l,nrad in shells[species[a]]:
                for _ in range(nrad):
                    for mi,m in enumerate(range(-l,l+1)):
                        for k in (range(3) if len(cycle)==3 else [(-m)%3]):
                            rows=[];values=[]
                            for power,at in enumerate(cycle):
                                coef=np.exp(-2j*np.pi*(m+k)*power/3)/np.sqrt(len(cycle))
                                for mu in range(2*l+1):
                                    z=coef*cc[l][mu,mi].conjugate()
                                    if abs(z)>1e-13:rows.append(spin*n+int(offset[at])+local+mu);values.append(z)
                            cols[(k+spin)%3].append((rows,values))
                    local+=2*l+1
    ri=[];cv=[];cp=[0];dims=[len(v) for v in cols]
    for sector in cols:
        for rows,values in sector:ri.extend(rows);cv.extend(values);cp.append(len(ri))
    u=clean(sparse.csc_matrix((cv,ri,cp),shape=(2*n,2*n)));r=action(species,defs,perm,Q3,U3);eye=sparse.eye(2*n,format='csc')
    eig=np.repeat(np.exp(-1j*np.pi/3+2j*np.pi*np.arange(3)/3),dims)
    checks={'C3_spinor_cube_error':maxerr(r@r@r+eye),'C3_basis_generator_error':maxerr(r@u-u@sparse.diags(eig))}
    result={'U':u,'R':r,'dimensions':dims,'N':n,'checks':checks}
    if c2perm is not None:
        q2=np.asarray(q2);require(np.max(abs(q2@q2-np.eye(3)))<1e-12,'invalid C2 rotation')
        _,ax=linalg.eigh((q2+q2.T)/2);axis=ax[:,-1]
        if axis[np.argmax(abs(axis))]<0:axis=-axis
        spin=-1j*np.array([[axis[2],axis[0]-1j*axis[1]],[axis[0]+1j*axis[1],-axis[2]]])
        c2=action(species,defs,np.asarray(c2perm),q2,spin)
        require(dims[0]==dims[1] and dims[2]%2==0,'D3 sector dimensions incompatible')
        m=dims[0];v=u[:,2*m:];a=u[:,:m];b=clean(c2@a);j=clean(-1j*(v.conj().T@c2@v));half=v.shape[1]//2
        require(np.all(np.diff(j.indptr)==1),'C2 action is not sparse paired in C3=-1 sector')
        wr=[];wc=[];plus=[];minus=[];seen=set()
        for at in range(2*half):
            if at in seen:continue
            start=j.indptr[at];partner=int(j.indices[start]);alpha=j.data[start]
            require(partner!=at and partner not in seen and abs(abs(alpha)-1)<1e-12,'invalid C2 pairing')
            k=len(seen)//2;seen.update([at,partner]);wr.extend([at,partner]);wc.extend([k,k]);plus.extend([1/np.sqrt(2),alpha/np.sqrt(2)]);minus.extend([1/np.sqrt(2),-alpha/np.sqrt(2)])
        p=clean(v@sparse.csc_matrix((plus,(wr,wc)),shape=(2*half,half)));q=clean(v@sparse.csc_matrix((minus,(wr,wc)),shape=(2*half,half)))
        u=clean(sparse.hstack([a,b,p,q],format='csc'));dims=[m,m,half,half]
        checks.update(C2_spinor_square_error=maxerr(c2@c2+eye),D3_spinor_group_relation_error=maxerr(c2@r@c2.conj().T-r.conj().T),D3_partner_gauge_error=maxerr(c2@a-b),C2_plus_sector_error=maxerr(c2@p-1j*p),C2_minus_sector_error=maxerr(c2@q+1j*q))
        result.update(U=u,dimensions=dims,C2=c2,C2_axis=axis,spin2=spin)
    checks['basis_unitarity_error']=maxerr(u.conj().T@u-eye)
    require(max(checks.values())<1e-12,f'spinor representation checks failed: {checks}')
    return result
