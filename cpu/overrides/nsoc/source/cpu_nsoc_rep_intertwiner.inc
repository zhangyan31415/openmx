/* Proof of the actual RS AO action on the solver's complete sparse basis.
 * Work is sparse-column action; storage is O(N), never an N*N group matrix. */
typedef struct {int nat,*ao,*atom,*orb,*rotation_offset,*seen,*touched,epoch,count;const int*types;const int*norb;const int*perm;const double*d;double*value;} RepAction;
static void rep_add(RepAction*a,int row,double v){
 if(a->seen[row]!=a->epoch){a->seen[row]=a->epoch;a->value[row]=0.;a->touched[a->count++]=row;}a->value[row]+=v;
}
static void rep_column_action(RepAction*a,int column){
 a->count=0;if(a->epoch==INT_MAX){memset(a->seen,0,(size_t)B.N*sizeof(int));a->epoch=0;}a->epoch++;
 for(int64_t z=B.cp[column];z<B.cp[column+1];z++){
  int row=B.ri[z],atom=a->atom[row],orb=a->orb[row],species=a->types[atom],width=a->norb[species],destination=a->perm[atom];
  const double*d=a->d+a->rotation_offset[species];
  for(int i=0;i<width;i++)rep_add(a,a->ao[destination]+i,d[i*width+orb]*B.cv[z]);
 }
}
static double rep_dot(RepAction*a,int column){double v=0.;for(int64_t z=B.cp[column];z<B.cp[column+1];z++){int row=B.ri[z];if(a->seen[row]==a->epoch)v+=B.cv[z]*a->value[row];}return v;}
static void rep_subtract(RepAction*a,int column,double coefficient){for(int64_t z=B.cp[column];z<B.cp[column+1];z++)rep_add(a,B.ri[z],-coefficient*B.cv[z]);}
static int representative_intertwiner(int ng,int nat,const int*const*perms,const int*types,const int*norb,const double*const*rotations){
 RepAction a={0};int ns=0,okay=1;double maximum=0.,form_maximum=0.,orth_maximum=0.;
 for(int atom=0;atom<nat;atom++){if(types[atom]<0)return 0;if(types[atom]>=ns)ns=types[atom]+1;}
 a.nat=nat;a.types=types;a.norb=norb;a.ao=mem(nat+1,sizeof(int));a.atom=mem(B.N,sizeof(int));a.orb=mem(B.N,sizeof(int));a.rotation_offset=mem(ns,sizeof(int));a.seen=mem(B.N,sizeof(int));a.touched=mem(B.N,sizeof(int));a.value=mem(B.N,sizeof(double));
 for(int s=1;s<ns;s++)a.rotation_offset[s]=a.rotation_offset[s-1]+norb[s-1]*norb[s-1];
 for(int atom=0;atom<nat;atom++){int width=norb[types[atom]];if(width<=0||width>128||a.ao[atom]>B.N-width){okay=0;break;}a.ao[atom+1]=a.ao[atom]+width;for(int i=0;i<width;i++){a.atom[a.ao[atom]+i]=atom;a.orb[a.ao[atom]+i]=i;}}
 if(a.ao[nat]!=B.N)okay=0;
 for(int g=0;g<ng&&okay;g++){
  a.perm=perms[g];a.d=rotations[g];
  for(int atom=0;atom<nat;atom++)if(a.perm[atom]<0||a.perm[atom]>=nat||types[a.perm[atom]]!=types[atom])okay=0;
  if(!okay)break;
  double e[4],character[NSYM_MAX_SECTORS]={0};
  rep_column_action(&a,B.sec[0].ofs);e[0]=rep_dot(&a,B.sec[0].ofs);e[2]=rep_dot(&a,B.sec[1].ofs);
  rep_column_action(&a,B.sec[1].ofs);e[1]=rep_dot(&a,B.sec[0].ofs);e[3]=rep_dot(&a,B.sec[1].ofs);
  for(int i=0;i<2;i++)for(int j=0;j<2;j++){double x=e[i]*e[j]+e[2+i]*e[2+j]-(i==j);if(!isfinite(x)||fabs(x)>1e-10)okay=0;if(isfinite(x))orth_maximum=fmax(orth_maximum,fabs(x));}
  if(c3_basis()){double form=fmax(fabs(e[0]-e[3]),fabs(e[1]+e[2]));form_maximum=fmax(form_maximum,form);if(!isfinite(form)||form>1e-10)okay=0;}
  for(int k=2;k<B.nsec;k++){rep_column_action(&a,B.sec[k].ofs);character[k]=rep_dot(&a,B.sec[k].ofs);if(!isfinite(character[k])||fabs(fabs(character[k])-1.)>1e-10)okay=0;}
  for(int k=0;k<B.nsec&&okay;k++)for(int offset=0;offset<B.sec[k].n;offset++){
   int column=B.sec[k].ofs+offset;rep_column_action(&a,column);
   if(k==0){rep_subtract(&a,column,e[0]);rep_subtract(&a,B.sec[1].ofs+offset,e[2]);}
   else if(k==1){rep_subtract(&a,B.sec[0].ofs+offset,e[1]);rep_subtract(&a,column,e[3]);}
   else rep_subtract(&a,column,character[k]);
   for(int z=0;z<a.count;z++){double v=a.value[a.touched[z]];if(!isfinite(v)){okay=0;break;}maximum=fmax(maximum,fabs(v));if(fabs(v)>1e-10){okay=0;break;}}
  }
 }
 free(a.ao);free(a.atom);free(a.orb);free(a.rotation_offset);free(a.seen);free(a.touched);free(a.value);
 if(rank0()){printf("CPU_NSOC_REP_INTERTWINER group_ops=%d all_columns=%d max_abs=%.12e orth_max_abs=%.12e form_max_abs=%.12e tolerance=1e-10 admitted=%d sparse_action=1 dense_N2_extra=0 gamma_periodic_phase=1\n",ng,B.N,maximum,orth_maximum,form_maximum,okay);fflush(stdout);}return okay;
}
