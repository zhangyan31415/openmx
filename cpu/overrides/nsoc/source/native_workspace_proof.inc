/* These records come from the exact size variables used by DFT mallocs.
 * They never derive allocation size from an environment flag or log text. */
static double *direct_workspace[4];static size_t direct_workspace_count[4];
void CpuNSym_DirectWorkspaceRecord(double*w,size_t nw,double*c,size_t nc,double*ed,size_t ne,double*pd,size_t np){
 direct_workspace[0]=w;direct_workspace[1]=c;direct_workspace[2]=ed;direct_workspace[3]=pd;direct_workspace_count[0]=nw;direct_workspace_count[1]=nc;direct_workspace_count[2]=ne;direct_workspace_count[3]=np;
}
void CpuNSym_DirectWorkspaceProve(double*w,double*c,double*ed,double*pd){
 if(direct_mode!=1)return;int okay=w&&c&&ed&&w==direct_workspace[0]&&c==direct_workspace[1]&&ed==direct_workspace[2]&&pd==direct_workspace[3]&&direct_workspace_count[0]==1&&direct_workspace_count[1]==1&&direct_workspace_count[2]==1,all;
 MPI_Allreduce(&okay,&all,1,MPI_INT,MPI_MIN,mpi_comm_level1);if(!all)die("DIRECT_SECTORS actual DFT scalar Work1/CDM1/EDM1 allocation proof failed");int rank;MPI_Comm_rank(mpi_comm_level1,&rank);
 printf("CPU_NSOC_DIRECT_WORK1 rank=%d allocated_elements=%zu CDM_elements=%zu EDM_elements=%zu PDM_elements=%zu allocation_source=actual_DFT_malloc_counts pointer_identity=1\n",rank,direct_workspace_count[0],direct_workspace_count[1],direct_workspace_count[2],direct_workspace_count[3]);fflush(stdout);
}
