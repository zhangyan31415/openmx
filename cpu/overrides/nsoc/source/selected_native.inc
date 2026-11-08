/* Included by selected_left.c, reusing its independent sparse-U plan.
 * The existing summation order is retained. Whole sparse components are
 * batched for value-only native redistribution; no full AO matrix is made. */
int CfsSelectedNativeLeft(CfsSelectedPlan*p,MPI_Comm world,const double*A,int lda,
 const int*inv,int wanted,int first,int last,double*output,size_t capacity,
 const int*ends,int sectors,const size_t*packed_offsets,int tile,const double*audit,int profile,
 CfsSelectedTimes*t,double*audit_error,char*error,size_t ez)
{
 if(!p)return fail(error,ez,"native sparse plan missing");
 int rank,np;MPI_Comm_rank(world,&rank);MPI_Comm_size(world,&np);
 int nk=last>=first?last-first+1:0;
 int okay=A&&output&&A!=output&&inv&&ends&&sectors>0&&sectors<=48&&
  lda>=(p->nr?p->nr:1)&&wanted>0&&wanted<=p->n&&tile>=p->maxblock&&tile<=1024&&
  np==p->pr*p->pc&&(!nk||(first>=1&&last<=wanted))&&
  (size_t)nk<=SIZE_MAX/sizeof(double)/(size_t)p->n&&capacity>=(size_t)nk*p->n;
 int all;MPI_Allreduce(&okay,&all,1,MPI_INT,MPI_MIN,world);
 if(!all)return fail(error,ez,"invalid native sparse return layout/buffers");
 int *meta=calloc((size_t)np*4+(size_t)p->n*4+(size_t)p->nc+2*(size_t)np+2,sizeof(int));
 if(!meta){MPI_Abort(world,96);abort();}
 int *coords=meta,*ranges=coords+2*np,*owner_state=ranges+2*np;
 int *sector=owner_state+p->n,*seen=sector+p->n,*recv_states=seen+p->n;
 int *columns=recv_states+p->n,*sc=columns+p->nc,*rc=sc+np;
 int me[2]={p->row,p->col},range[2]={first,last};
 MPI_Allgather(me,2,MPI_INT,coords,2,MPI_INT,world);
 MPI_Allgather(range,2,MPI_INT,ranges,2,MPI_INT,world);
 for(int q=0;q<p->n;q++)owner_state[q]=-1;
 for(int r=0;r<np;r++){
  if(coords[2*r]<0||coords[2*r]>=p->pr||coords[2*r+1]<0||coords[2*r+1]>=p->pc)okay=0;
  for(int s=0;s<r;s++)if(coords[2*s]==coords[2*r]&&coords[2*s+1]==coords[2*r+1])okay=0;
  int a=ranges[2*r],b=ranges[2*r+1];
  if(b>=a){if(a<1||b>wanted)okay=0;else for(int k=a-1;k<b;k++){if(owner_state[k]!=-1)okay=0;owner_state[k]=r;}}
 }
 for(int k=0;k<wanted;k++)if(owner_state[k]<0)okay=0;
 int begin=0;unsigned long long hash=1469598103934665603ULL;
 for(int s=0;s<sectors;s++){
  if(ends[s]<=begin||ends[s]>p->n)okay=0;
  else for(int q=begin;q<ends[s];q++)sector[q]=s;
  begin=ends[s];hash^=(unsigned)ends[s];hash*=1099511628211ULL;
 }
 if(begin!=p->n)okay=0;
 for(int q=0;q<p->n;q++){
  if(inv[q]<0||inv[q]>=p->n)okay=0;else if(seen[inv[q]]++)okay=0;
  hash^=(unsigned)inv[q];hash*=1099511628211ULL;
 }
 hash^=(unsigned)wanted;hash*=1099511628211ULL;
 unsigned long long h[2]={hash,~hash},hh[2];MPI_Allreduce(h,hh,2,MPI_UNSIGNED_LONG_LONG,MPI_MIN,world);
 if(hh[0]!=hash||hh[1]!=~hash)okay=0;
 MPI_Allreduce(&okay,&all,1,MPI_INT,MPI_MIN,world);
 if(!all){free(meta);return fail(error,ez,"native sparse return ownership/permutation mismatch");}
 int local_row_start[48],local_col_start[48],local_row_count[48];
 int rs=0,cs=0,previous=0;
 for(int s=0;s<sectors;s++){
  local_row_start[s]=rs;local_col_start[s]=cs;
  for(int q=previous;q<ends[s];q++){rs+=owner(q,p->nb,p->pr)==p->row;cs+=owner(q,p->nb,p->pc)==p->col;}
  local_row_count[s]=rs-local_row_start[s];previous=ends[s];
 }
 int count=0;for(int j=0;j<p->nc;j++){
  int q=(j/p->nb)*p->nb*p->pc+p->col*p->nb+j%p->nb;
  if(q<p->n&&inv[q]<wanted){columns[count++]=j;sc[owner_state[inv[q]]]++;}
 }
 /* Received states are ordered by unsorted source column, just as sends. */
 int recv_offsets[1025];if(p->pc>1024){free(meta);return fail(error,ez,"native return grid too wide");}
 int at=0;for(int c=0;c<p->pc;c++){
  recv_offsets[c]=at;
  for(int q=0;q<p->n;q++)if(owner(q,p->nb,p->pc)==c&&inv[q]>=first-1&&inv[q]<last)recv_states[at++]=inv[q]-(first-1);
 }
 recv_offsets[p->pc]=at;if(at!=nk){free(meta);return fail(error,ez,"native return local state count mismatch");}
 int *counts=calloc((size_t)np*4+2,sizeof(int));
 size_t maxout=(size_t)tile*count,maxrecv=(size_t)tile*nk;
 double *batch=malloc((maxout?maxout:1)*sizeof(double));
 double *send=malloc((maxout?maxout:1)*sizeof(double)),*recv=malloc((maxrecv?maxrecv:1)*sizeof(double));
 if(!counts||!batch||!send||!recv){MPI_Abort(world,96);abort();}
 int *sendcounts=counts,*recvcounts=sendcounts+np,*sd=recvcounts+np,*rd=sd+np+1;
 rc[0]=0;for(int r=0;r<np;r++){rc[r+1]=rc[r]+sc[r];sd[r]=rc[r];}
 for(int j=0;j<count;j++){
  int lc=columns[j],q=(lc/p->nb)*p->nb*p->pc+p->col*p->nb+lc%p->nb;
  seen[sd[owner_state[inv[q]]]++]=j;
 }
 CfsSelectedTimes tt={0};tt.columns=count;
 tt.plan_bytes=p->panel_bytes+p->metadata_bytes+(maxout*2+maxrecv)*sizeof(double);
 double coeff[NSOC_SPARSE_MAX_BLOCK*NSOC_SPARSE_MAX_BLOCK],maxerror=0.;int bad=0;
 for(int component0=0,component1;component0<p->components;component0=component1){
  component1=component0+1;
  int row0=p->ptr[component0];while(component1<p->components&&p->ptr[component1+1]-row0<=tile)component1++;
  int row1=p->ptr[component1],localrows=0;tt.batches++;
  double clock=profile?MPI_Wtime():0.;
  for(int component=component0;component<component1;component++){
   int b=p->ptr[component],e=p->ptr[component+1],m=e-b;
   tt.logical_values+=(unsigned long long)m*count;
   memset(p->panel,0,(size_t)m*count*sizeof(double));
   for(int a=0;a<m;a++){
    int g=p->cols[b+a];if(owner(g,p->nb,p->pr)!=p->row)continue;
    int i=local_index(g,p->nb,p->pr);
    for(int j=0;j<count;j++){
     int lc=columns[j],q=(lc/p->nb)*p->nb*p->pc+p->col*p->nb+lc%p->nb;
     if(sector[g]==sector[q]){
      int s=sector[g];size_t offset=packed_offsets?
       packed_offsets[s]+(size_t)(lc-local_col_start[s])*local_row_count[s]+i-local_row_start[s]:
       (size_t)lc*lda+i;
      p->panel[(size_t)a*count+j]=A[offset];
     }
    }
   }
   if(profile){tt.pack+=MPI_Wtime()-clock;clock=MPI_Wtime();}
   if(count&&p->pr>1)MPI_Allreduce(MPI_IN_PLACE,p->panel,m*count,MPI_DOUBLE,MPI_SUM,p->colcomm);
   if(profile){tt.collective+=MPI_Wtime()-clock;clock=MPI_Wtime();}
   memset(coeff,0,(size_t)m*m*sizeof(double));
   for(int a=b;a<e;a++)for(int64_t z=p->cp[p->cols[a]];z<p->cp[p->cols[a]+1];z++){
    double v=p->cv[z];if(v!=0.)coeff[(size_t)p->row_slot[p->ri[z]]*m+a-b]+=v;
   }
   for(int a=b;a<e;a++){
    int g=p->rows[a];if(owner(g,p->nb,p->pr)!=p->row)continue;
    double*out=batch+(size_t)localrows++*count;memset(out,0,(size_t)count*sizeof(double));
    for(int k=0;k<m;k++){
     double v=coeff[(size_t)(a-b)*m+k];if(v==0.)continue;
     for(int j=0;j<count;j++)out[j]+=v*p->panel[(size_t)k*count+j];
    }
    for(int j=0;j<count;j++){
     if(!isfinite(out[j]))bad=1;
     if(audit){double reference=audit[(size_t)columns[j]*lda+local_index(g,p->nb,p->pr)];if(!isfinite(reference))bad=1;else if(fabs(reference-out[j])>maxerror)maxerror=fabs(reference-out[j]);}
    }
   }
   if(profile){tt.multiply+=MPI_Wtime()-clock;clock=MPI_Wtime();}
  }
  sd[0]=rd[0]=0;
  for(int r=0;r<np;r++){
   int rows=0;for(int a=row0;a<row1;a++)if(owner(p->rows[a],p->nb,p->pr)==coords[2*r])rows++;
   long long ns=(long long)localrows*sc[r],nr=(long long)rows*(recv_offsets[coords[2*r+1]+1]-recv_offsets[coords[2*r+1]]);
   if(ns>INT_MAX||nr>INT_MAX||(long long)sd[r]+ns>INT_MAX||(long long)rd[r]+nr>INT_MAX){MPI_Abort(world,96);abort();}
   sendcounts[r]=(int)ns;recvcounts[r]=(int)nr;sd[r+1]=sd[r]+(int)ns;rd[r+1]=rd[r]+(int)nr;
   int k=sd[r];for(int a=0;a<localrows;a++)for(int j=rc[r];j<rc[r+1];j++)send[k++]=batch[(size_t)a*count+seen[j]];
   if(k!=sd[r+1])bad=1;
  }
  if(profile){tt.pack+=MPI_Wtime()-clock;clock=MPI_Wtime();}
  MPI_Alltoallv(send,sendcounts,sd,MPI_DOUBLE,recv,recvcounts,rd,MPI_DOUBLE,world);
  if(profile){tt.collective+=MPI_Wtime()-clock;clock=MPI_Wtime();}
  for(int r=0;r<np;r++){
   int k=rd[r],c=coords[2*r+1];
   for(int a=row0;a<row1;a++)if(owner(p->rows[a],p->nb,p->pr)==coords[2*r]){
    double*out=output+(size_t)p->rows[a]*nk;
    for(int j=recv_offsets[c];j<recv_offsets[c+1];j++)out[recv_states[j]]=recv[k++];
   }
   if(k!=rd[r+1])bad=1;
  }
  if(profile)tt.pack+=MPI_Wtime()-clock;
 }
 int anybad;MPI_Allreduce(&bad,&anybad,1,MPI_INT,MPI_MAX,world);
 double globalerror;MPI_Allreduce(&maxerror,&globalerror,1,MPI_DOUBLE,MPI_MAX,world);
 if(audit_error)*audit_error=globalerror;if(t)*t=tt;
 free(counts);free(batch);free(send);free(recv);free(meta);
 if(anybad)return fail(error,ez,"native sparse return produced invalid values/counts");
 return 0;
}
