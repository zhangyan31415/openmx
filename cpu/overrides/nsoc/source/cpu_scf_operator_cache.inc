/* Native position/momentum operators are geometry-only. Compute them before
 * ending the orbital-grid lifetime and transfer ownership to SCF2File once. */
static double ******cpu_scf_cached_position=NULL;
static double *****cpu_scf_cached_momentum=NULL;
static void CpuSCFFileAllocateOperators(double *******position,double ******momentum)
{
  int Mc_AN,Gc_AN,tno0,tno1,Cwan,h_AN,Gh_AN,Hwan,i,direction,order;
  double ******OLPpo;
  double *****OLPmo;
  /* allocation of array */

  /* Modified by N. Yamaguchi ***/
  /* OLPpo: matrix for position operator */


  OLPpo=(double******)malloc(sizeof(double*****)*3);
  for (direction=0; direction<3; direction++){
    OLPpo[direction]=(double*****)malloc(sizeof(double****)*order_max);
    for (order=0; order<order_max; order++){
      OLPpo[direction][order]=(double****)malloc(sizeof(double***)*(Matomnum+1));
      FNAN[0] = 0;
      for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

	if (Mc_AN==0){
	  Gc_AN = 0;
	  tno0 = 1;
	}
	else{
	  Gc_AN = M2G[Mc_AN];
	  Cwan = WhatSpecies[Gc_AN];
	  tno0 = Spe_Total_CNO[Cwan];
	}

	OLPpo[direction][order][Mc_AN] = (double***)malloc(sizeof(double**)*(FNAN[Gc_AN]+1));
	for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	  if (Mc_AN==0){
	    tno1 = 1;
	  }
	  else{
	    Gh_AN = natn[Gc_AN][h_AN];
	    Hwan = WhatSpecies[Gh_AN];
	    tno1 = Spe_Total_CNO[Hwan];
	  }

	  OLPpo[direction][order][Mc_AN][h_AN] = (double**)malloc(sizeof(double*)*tno0);
	  for (i=0; i<tno0; i++){
	    OLPpo[direction][order][Mc_AN][h_AN][i] = (double*)malloc(sizeof(double)*tno1);
	  }
	}
      }
    }
  }

  /* OLPmo: matrix for momentum operator */

  OLPmo=(double*****)malloc(sizeof(double****)*3);
  for (direction=0; direction<3; direction++){

    OLPmo[direction]=(double****)malloc(sizeof(double***)*(Matomnum+1));
    FNAN[0] = 0;
    for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){

      if (Mc_AN==0){
	Gc_AN = 0;
	tno0 = 1;
      }
      else{
	Gc_AN = M2G[Mc_AN];
	Cwan = WhatSpecies[Gc_AN];
	tno0 = Spe_Total_CNO[Cwan];
      }

      OLPmo[direction][Mc_AN] = (double***)malloc(sizeof(double**)*(FNAN[Gc_AN]+1));
      for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	if (Mc_AN==0){
	  tno1 = 1;
	}
	else{
	  Gh_AN = natn[Gc_AN][h_AN];
	  Hwan = WhatSpecies[Gh_AN];
	  tno1 = Spe_Total_CNO[Hwan];
	}

	OLPmo[direction][Mc_AN][h_AN] = (double**)malloc(sizeof(double*)*tno0);
	for (i=0; i<tno0; i++){
	  OLPmo[direction][Mc_AN][h_AN][i] = (double*)malloc(sizeof(double)*tno1);
	}
      }
    }
  }

  /* end of allocation of OLPpo and OLPp */

  *position=OLPpo;*momentum=OLPmo;
}
void CpuSCFFileCacheOperators(void)
{
  int rank;
  if(HS_fileout!=1)return;
  if(cpu_scf_cached_position||cpu_scf_cached_momentum){
    fprintf(stderr,"CPU_SCF_OPERATOR_CACHE duplicate preparation\n");abort();
  }
  CpuSCFFileAllocateOperators(&cpu_scf_cached_position,&cpu_scf_cached_momentum);
  Calc_OLPpo(cpu_scf_cached_position);
  Calc_OLPmo(cpu_scf_cached_momentum);
  MPI_Comm_rank(mpi_comm_level1,&rank);
  if(rank==Host_ID){printf("CPU_SCF_OPERATOR_CACHE prepared=1 position_and_momentum=1 native_quadrature=unchanged additional_integrals=0\n");fflush(stdout);}
}
static int CpuSCFFileCacheTake(double *******position,double ******momentum)
{
  if(!cpu_scf_cached_position)return 0;
  *position=cpu_scf_cached_position;*momentum=cpu_scf_cached_momentum;
  cpu_scf_cached_position=NULL;cpu_scf_cached_momentum=NULL;return 1;
}
