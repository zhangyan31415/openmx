/**********************************************************************
  Hamiltonian_Cluster.c:

     Hamiltonian_Cluster.c is a subroutine to make a Hamiltonian matrix
     for cluster or molecular systems.

  Log of Hamiltonian_Cluster.c:

     22/Nov/2001  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include "openmx_common.h"
#include "cpu_hassembly_workspace.h"
#include "mpi.h"
#include "cpu_nsoc_owner_openmx.h"

double *Hex;

/* Compact Work1 is valid only with the owner DM/EDM path. The DFT entry
 * performs this collective agreement before any global workspace allocation. */
int CpuHAssembly_Compact(void)
{
  const char *s=getenv("OPENMX_CPU_NSOC_COMPACT_WORK1");
  if(!s || !strcmp(s,"0"))return 0;
  if(!strcmp(s,"1"))return 1;
  fprintf(stderr,"CPU_NSOC_COMPACT_WORK1_REJECT flag must be0 or1\n");
  MPI_Abort(mpi_comm_level1,98);abort();
}
void CpuHAssembly_Validate(void)
{
  int local[2]={CpuHAssembly_Compact(),CpuNSOCOwnerMode()},lo[2],hi[2];
  MPI_Allreduce(local,lo,2,MPI_INT,MPI_MIN,mpi_comm_level1);
  MPI_Allreduce(local,hi,2,MPI_INT,MPI_MAX,mpi_comm_level1);
  if(lo[0]!=hi[0] || (hi[0]&&!lo[1])){
    fprintf(stderr,"CPU_NSOC_COMPACT_WORK1_REJECT requires agreed flag and validated owner DM mode\n");
    MPI_Abort(mpi_comm_level1,98);abort();
  }
}

/* Native broadcasts already visit owner segments in ascending rank order.
 * Consume each segment immediately, retaining the exact native atom / AO-row /
 * neighbor / AO-column order and += order for repeated periodic AO pairs.
 * No extra assembly collective and no cached matrix values are introduced. */
static void cpu_hassembly_stream(double ****RH,double *Hs,const int *MP,
  int spin,int myworld1,int rank,int np,const int *counts,
  const int *atoms,const int *first_atom,const int *order,int global_count)
{
  int max_count=0;
  for(int id=0;id<np;id++)if(counts[id]>max_count)max_count=counts[id];
  double *panel=malloc((size_t)(max_count?max_count:1)*sizeof(double));
  if(!panel){fprintf(stderr,"CPU_NSOC_COMPACT_WORK1_REJECT panel allocation\n");MPI_Abort(mpi_comm_level1,98);abort();}
  if(spin==myworld1)for(size_t q=0;q<(size_t)na_rows*(size_t)na_cols;q++)Hs[q]=0.;
  for(int id=0;id<np;id++){
    int at=0;
    if(rank==id){
      for(int m=1;m<=Matomnum;m++){
        int g=M2G[m],nr=Spe_Total_CNO[WhatSpecies[g]];
        for(int i=0;i<nr;i++)for(int h=0;h<=FNAN[g];h++){
          int nc=Spe_Total_CNO[WhatSpecies[natn[g][h]]];
          for(int j=0;j<nc;j++)panel[at++]=RH[m][h][i][j];
        }
      }
      if(at!=counts[id]){fprintf(stderr,"CPU_NSOC_COMPACT_WORK1_REJECT pack count\n");MPI_Abort(mpi_comm_level1,98);abort();}
    }
    MPI_Bcast(panel,counts[id],MPI_DOUBLE,id,mpi_comm_level1);
    if(spin==myworld1){
      at=0;
      for(int an=first_atom[id];an<first_atom[id]+atoms[id];an++){
        int g=order[an],nr=Spe_Total_CNO[WhatSpecies[g]],anum=MP[g];
        for(int i=0;i<nr;i++)for(int h=0;h<=FNAN[g];h++){
          int q=natn[g][h],nc=Spe_Total_CNO[WhatSpecies[q]],bnum=MP[q];
          for(int j=0;j<nc;j++){
            int ig=anum+i,jg=bnum+j,brow=(ig-1)/nblk,bcol=(jg-1)/nblk;
            int prow=brow%np_rows,pcol=bcol%np_cols;
            if(my_prow==prow&&my_pcol==pcol){
              int il=(brow/np_rows)*nblk+(ig-1)%nblk;
              int jl=(bcol/np_cols)*nblk+(jg-1)%nblk;
              Hs[(size_t)jl*(size_t)na_rows+il]+=panel[at];
            }
            at++;
          }
        }
      }
      if(at!=counts[id]){fprintf(stderr,"CPU_NSOC_COMPACT_WORK1_REJECT unpack count\n");MPI_Abort(mpi_comm_level1,98);abort();}
    }
  }
  static int reported=0;
  if(!reported&&!rank){
    printf("CPU_NSOC_COMPACT_WORK1 active=1 global_elements=%d max_owner_elements=%d per_rank_panel_bytes=%zu retained_Work1_elements=1 assembly_collectives_unchanged=1 native_order_unchanged=1\n",global_count,max_count,(size_t)(max_count?max_count:1)*sizeof(double));fflush(stdout);reported=1;
  }
  free(panel);
}



static void cpu_hassembly_impl(double ****RH, double *Hs, int *MP, int spin, int myworld1, double *workspace, size_t capacity)
{
  int i,j,k;
  int MA_AN,GA_AN,LB_AN,GB_AN,AN;
  int wanA,wanB,tnoA,tnoB,Anum,Bnum,NUM;
  int num,tnum,num_orbitals;
  int ID,myid,numprocs,tag=999;
  int *My_NZeros;
  int *is1,*ie1,*is2;
  int *My_Matomnum,*order_GA;
  double *H1,sum;
  MPI_Status stat;
  MPI_Request request;
  int ig,jg,il,jl,prow,pcol,brow,bcol;

  /* MPI */

  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);
  MPI_Barrier(mpi_comm_level1);

  /* allocation of arrays */

  My_NZeros = (int*)malloc(sizeof(int)*numprocs);
  My_Matomnum = (int*)malloc(sizeof(int)*numprocs);
  is1 = (int*)malloc(sizeof(int)*numprocs);
  ie1 = (int*)malloc(sizeof(int)*numprocs);
  is2 = (int*)malloc(sizeof(int)*numprocs);
  order_GA = (int*)malloc(sizeof(int)*(atomnum+2));

  /* find my total number of non-zero elements in myid */

  My_NZeros[myid] = 0;
  for (MA_AN=1; MA_AN<=Matomnum; MA_AN++){
    GA_AN = M2G[MA_AN];
    wanA = WhatSpecies[GA_AN];
    tnoA = Spe_Total_CNO[wanA];

    num = 0;      
    for (LB_AN=0; LB_AN<=FNAN[GA_AN]; LB_AN++){
      GB_AN = natn[GA_AN][LB_AN];
      wanB = WhatSpecies[GB_AN];
      tnoB = Spe_Total_CNO[wanB];
      num += tnoB;
    }

    My_NZeros[myid] += tnoA*num;
  }

  for (ID=0; ID<numprocs; ID++){
    MPI_Bcast(&My_NZeros[ID],1,MPI_INT,ID,mpi_comm_level1);
  }

  tnum = 0;
  for (ID=0; ID<numprocs; ID++){
    tnum += My_NZeros[ID];
  }  

  is1[0] = 0;
  ie1[0] = My_NZeros[0] - 1;

  for (ID=1; ID<numprocs; ID++){
    is1[ID] = ie1[ID-1] + 1;
    ie1[ID] = is1[ID] + My_NZeros[ID] - 1;
  }  

  /* set is2 and order_GA */

  My_Matomnum[myid] = Matomnum;
  for (ID=0; ID<numprocs; ID++){
    MPI_Bcast(&My_Matomnum[ID],1,MPI_INT,ID,mpi_comm_level1);
  }

  is2[0] = 1;
  for (ID=1; ID<numprocs; ID++){
    is2[ID] = is2[ID-1] + My_Matomnum[ID-1];
  }
  
  for (MA_AN=1; MA_AN<=Matomnum; MA_AN++){
    order_GA[is2[myid]+MA_AN-1] = M2G[MA_AN];
  }

  for (ID=0; ID<numprocs; ID++){
    MPI_Bcast(&order_GA[is2[ID]],My_Matomnum[ID],MPI_INT,ID,mpi_comm_level1);
  }

  /* set MP */

  Anum = 1;
  for (i=1; i<=atomnum; i++){
    MP[i] = Anum;
    wanA = WhatSpecies[i];
    Anum += Spe_Total_CNO[wanA];
  }
  NUM = Anum - 1;

  /* set H1 */

  if(CpuHAssembly_Compact()){
    if(!CpuNSOCOwnerMode()){
      fprintf(stderr,"CPU_NSOC_COMPACT_WORK1_REJECT requires owner DM mode\n");MPI_Abort(mpi_comm_level1,98);abort();
    }
    cpu_hassembly_stream(RH,Hs,MP,spin,myworld1,myid,numprocs,My_NZeros,My_Matomnum,is2,order_GA,tnum);
    free(order_GA);free(is2);free(ie1);free(is1);free(My_Matomnum);free(My_NZeros);
    return;
  }

  /* Work1 is unused until the later DM reduction. Every H1 element is
   * rewritten and all blocking broadcasts complete before returning. */
  if(workspace){
    if(tnum<0 || (size_t)tnum>capacity){
      fprintf(stderr,"CPU_NSOC_HASM_WORKSPACE insufficient capacity required=%d supplied=%zu\n",tnum,capacity);
      MPI_Abort(mpi_comm_level1,98);
    }
    H1=workspace;
    static int reported=0;
    if(!reported && myid==0){printf("CPU_NSOC_HASM_WORKSPACE borrowed=1 required_elements=%d capacity_elements=%zu added_payload_bytes=0\n",tnum,capacity);fflush(stdout);reported=1;}
  }
  else H1 = (double*)malloc(sizeof(double)*(tnum+1));

  k = is1[myid];
  for (MA_AN=1; MA_AN<=Matomnum; MA_AN++){
    GA_AN = M2G[MA_AN];
    wanA = WhatSpecies[GA_AN];
    tnoA = Spe_Total_CNO[wanA];
    for (i=0; i<tnoA; i++){
      for (LB_AN=0; LB_AN<=FNAN[GA_AN]; LB_AN++){
        GB_AN = natn[GA_AN][LB_AN];
        wanB = WhatSpecies[GB_AN];
        tnoB = Spe_Total_CNO[wanB];
        for (j=0; j<tnoB; j++){
          H1[k] = RH[MA_AN][LB_AN][i][j]; 
          k++;
	}
      }
    }
  }

  /* MPI H1 */
    
  for (ID=0; ID<numprocs; ID++){
    k = is1[ID];
    MPI_Bcast(&H1[k], My_NZeros[ID], MPI_DOUBLE, ID, mpi_comm_level1);
  }

  /* H1 -> H */

  if (spin==myworld1){

    for(i=0;i<na_rows*na_cols;i++){
      Hs[i] = 0.0;
    }

    k = 0;
    for (AN=1; AN<=atomnum; AN++){
      GA_AN = order_GA[AN];
      wanA = WhatSpecies[GA_AN];
      tnoA = Spe_Total_CNO[wanA];
      Anum = MP[GA_AN];

      for (i=0; i<tnoA; i++){

	for (LB_AN=0; LB_AN<=FNAN[GA_AN]; LB_AN++){
	  GB_AN = natn[GA_AN][LB_AN];
	  wanB = WhatSpecies[GB_AN];
	  tnoB = Spe_Total_CNO[wanB];
	  Bnum = MP[GB_AN];

	  for (j=0; j<tnoB; j++){
	    ig = Anum+i;
	    jg = Bnum+j;
	    
	    brow = (ig-1)/nblk;
	    bcol = (jg-1)/nblk;

	    prow = brow%np_rows;
	    pcol = bcol%np_cols;

	    if(my_prow==prow && my_pcol==pcol){
	      il = (brow/np_rows+1)*nblk+1;
	      jl = (bcol/np_cols+1)*nblk+1;
	      if(((my_prow+np_rows)%np_rows) >= (brow%np_rows)){
		if(my_prow==prow){
		  il = il+(ig-1)%nblk;
		}
		il = il-nblk;
	      }
	      if(((my_pcol+np_cols)%np_cols) >= (bcol%np_cols)){
		if(my_pcol==pcol){
		  jl = jl+(jg-1)%nblk;
		}
		jl = jl-nblk;
	      }
	      Hs[(jl-1)*na_rows+il-1] += H1[k]; 
	    }
	    
	    k++;
	  }
	}
      }
    }

    /*
    for (i=0; i<na_rows; i++){
      for (j=0; j<na_rows; j++){
	Hs[i*na_rows+j] += 0.25*Hex[i*na_rows+j];
      }
    }
    */


    /*
    printf("WWW1\n");
    for (i=0; i<na_rows; i++){
      for (j=0; j<na_rows; j++){
        printf("%7.4f ", Hs[i*na_rows+j]);
      }
      printf("\n");
    }
    */
     
  } // if (spin==myworld1)

  /* freeing of arrays */

  if(!workspace)free(H1);
  free(order_GA);
  free(is2);
  free(ie1);
  free(is1);
  free(My_Matomnum);
  free(My_NZeros);

}


void Hamiltonian_Cluster_Hs(double ****RH,double *Hs,int *MP,int spin,int myworld1){
 cpu_hassembly_impl(RH,Hs,MP,spin,myworld1,NULL,0);
}
void CpuHAssembly_WithScratch(double ****RH,double *Hs,int *MP,int spin,int myworld1,double *workspace,size_t capacity){
 cpu_hassembly_impl(RH,Hs,MP,spin,myworld1,workspace,capacity);
}
