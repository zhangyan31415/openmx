"""Read OpenMX header12/15 exports and compare every periodic H/S atom block.

SOC compares four real H fields, three imaginary HNL fields and overlap,
without discarding spinor components or replacing nonorthogonal S by I.
Only header/atom metadata and a single block-sized scratch are resident.
"""
import hashlib,math,pathlib
import numpy as np

class Export:
    def __init__(self,path):
        self.path=pathlib.Path(path)
        with self.path.open('rb')as f:
            def read(n,dtype):
                v=np.fromfile(f,dtype,n)
                if len(v)!=n:raise ValueError('truncated metadata '+str(self.path))
                return v
            self.header=read(6,'<i4').tolist();nat,spin,_,_,_,cells=self.header
            if nat<=0 or spin not in (12,15)or cells<0:raise ValueError('unsupported OpenMX H/S header')
            self.order=read(1,'<i4').tolist();read((cells+1)*4,'<f8')
            ijk=read((cells+1)*4,'<i4').reshape(-1,4)[:,1:]
            self.orbitals=read(nat,'<i4');nn=read(nat,'<i4')
            if np.any(self.orbitals<=0)or np.any(nn<0):raise ValueError('invalid orbital/neighbor counts')
            neighbors=[read(int(n)+1,'<i4')for n in nn];images=[read(int(n)+1,'<i4')for n in nn]
            self.geometry=read(24+4*nat,'<f8');self.offset=f.tell()
        if not np.isfinite(self.geometry).all():raise ValueError('nonfinite geometry')
        self.blocks={};count=0
        for a,(bs,cs)in enumerate(zip(neighbors,images),1):
            for b,c in zip(bs,cs):
                b,c=int(b),int(c)
                if not 1<=b<=nat or not 0<=c<=cells:raise ValueError('invalid neighbor/image index')
                key=(a,b,*map(int,ijk[c]));shape=(int(self.orbitals[a-1]),int(self.orbitals[b-1]))
                if key in self.blocks:raise ValueError('duplicate periodic atom block')
                self.blocks[key]=(count,shape);count+=math.prod(shape)
        self.labels=['H','S']if spin==12 else ['H.real.0','H.real.1','H.real.2','H.real.3','H.imag.0','H.imag.1','H.imag.2','S']
        if self.path.stat().st_size<self.offset+len(self.labels)*count*8:raise ValueError('truncated H/S payload')
        self.arrays=np.memmap(self.path,'<f8',mode='r',offset=self.offset,shape=(len(self.labels),count))
        self.atom_count=nat;self.element_count=count
        self.geometry_sha256=hashlib.sha256(self.geometry.tobytes()).hexdigest()
    def close(self):self.arrays._mmap.close()

def compare(left,right):
    a=Export(left);b=None
    try:
        b=Export(right)
        if a.labels!=b.labels:raise ValueError('spin/field convention differs')
        if not np.array_equal(a.orbitals,b.orbitals)or not np.array_equal(a.geometry,b.geometry):raise ValueError('basis/lattice/atom geometry differs')
        keys=a.blocks.keys()|b.blocks.keys();stats={}
        for fi,name in enumerate(a.labels):
            maximum=squares=norm=additional=missing=0.;elements=0
            for key in keys:
                ls=a.blocks.get(key);rs=b.blocks.get(key);shape=ls[1]if ls else rs[1]
                if ls and rs and ls[1]!=rs[1]:raise ValueError('periodic block AO shape differs')
                size=math.prod(shape);x=a.arrays[fi,ls[0]:ls[0]+size]if ls else np.zeros(size);y=b.arrays[fi,rs[0]:rs[0]+size]if rs else np.zeros(size)
                if not np.isfinite(x).all()or not np.isfinite(y).all():raise ValueError('nonfinite H/S field '+name)
                delta=x-y;maximum=max(maximum,float(np.max(abs(delta))));squares+=float(np.dot(delta,delta));norm+=float(np.dot(x,x));elements+=size
                if not ls:additional=max(additional,float(np.max(abs(y))))
                if not rs:missing=max(missing,float(np.max(abs(x))))
            stats[name]={'max_absolute_difference':maximum,'RMS_difference':math.sqrt(squares/max(elements,1)),'relative_Frobenius_difference':math.sqrt(squares/max(norm,1e-300)),'additional_candidate_blocks_max_absolute_value':additional,'missing_candidate_blocks_max_reference_value':missing}
        return {'reference':str(a.path),'candidate':str(b.path),'atom_count':a.atom_count,'scalar_AO_count':int(sum(a.orbitals)),'header_spin':a.header[1],'geometry_sha256':a.geometry_sha256,'reference_blocks':len(a.blocks),'candidate_blocks':len(b.blocks),'added_candidate_blocks':len(b.blocks.keys()-a.blocks.keys()),'removed_candidate_blocks':len(a.blocks.keys()-b.blocks.keys()),'fields':stats,'max_H':max(v['max_absolute_difference']for k,v in stats.items()if k!='S'),'max_S':stats['S']['max_absolute_difference'],'comparison':'explicit union of periodic atom blocks; absent blocks are zero for comparison only; full exported H spinor fields and nonorthogonal S retained'}
    finally:
        a.close()
        if b is not None:b.close()
