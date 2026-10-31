static void CpuHNLRemoteAtom(int Mc_AN,int IDR,int Original_Mc_AN,double*****CDM0,double*****iDM0,CpuHNLWorkspace*w){
 int Gc_AN,Cwan,fan,h_AN,Gh_AN,Mh_AN,Hwan,ian,n,jg,j0,jg0,Mj_AN0,po2,q_AN,Gq_AN,Mq_AN,Qwan,jan,kl,kl1,i,j;
 double dEx,dEy,dEz,pref,Stime_atom,Etime_atom;
 dcomplex***Hx0=w->m[0],***Hy0=w->m[1],***Hz0=w->m[2],***Hx1=w->m[3],***Hy1=w->m[4],***Hz1=w->m[5];
 w->atoms++;


	  dtime(&Stime_atom);

	  dEx = 0.0; 
	  dEy = 0.0; 
	  dEz = 0.0; 

	  Gc_AN = M2G[Mc_AN];
	  Cwan = WhatSpecies[Gc_AN];
	  fan = FNAN[Gc_AN];

	  h_AN = 0;
	  Gh_AN = natn[Gc_AN][h_AN];
	  Mh_AN = F_G2M[Gh_AN];
	  Hwan = WhatSpecies[Gh_AN];
	  ian = Spe_Total_CNO[Hwan];

	  n = F_Rcv_Num_WK[IDR];
	  jg = Rcv_GAN[IDR][n];

	  for (j0=0; j0<=fan; j0++){

	    jg0 = natn[Gc_AN][j0];
	    Mj_AN0 = F_G2M[jg0];

	    po2 = 0;
	    if (Original_Mc_AN==Mj_AN0){
	      po2 = 1;
	      q_AN = j0;
	    }

	    if (po2==1){w->pairs++;

	      Gq_AN = natn[Gc_AN][q_AN];
	      Mq_AN = F_G2M[Gq_AN];
	      Qwan = WhatSpecies[Gq_AN];
	      jan = Spe_Total_CNO[Qwan];
	      kl = RMI1[Mc_AN][h_AN][q_AN];

 	      w->calls++;dHNL(0,Mc_AN,h_AN,q_AN,DS_NL,Hx0,Hy0,Hz0);

	      /* contribution of force = Trace(CDM0*dH) */
	      /* spin non-polarization */

	      if (SpinP_switch==0){

                if (q_AN==h_AN) pref = 2.0;
                else            pref = 4.0; 

		for (i=0; i<ian; i++){
		  for (j=0; j<jan; j++){

		    dEx += pref*CDM0[0][Mh_AN][kl][i][j]*Hx0[0][i][j].r;
		    dEy += pref*CDM0[0][Mh_AN][kl][i][j]*Hy0[0][i][j].r;
		    dEz += pref*CDM0[0][Mh_AN][kl][i][j]*Hz0[0][i][j].r;
		  }
		}
	      }

	      /* collinear spin polarized or non-colliear without SO and LDA+U */

	      else if (SpinP_switch==1 || (SpinP_switch==3 && SO_switch==0 && Hub_U_switch==0
		   && Constraint_NCS_switch==0 && Zeeman_NCS_switch==0 && Zeeman_NCO_switch==0)){ 

		if (q_AN==h_AN) pref = 1.0;
		else            pref = 2.0; 

		for (i=0; i<Spe_Total_CNO[Hwan]; i++){
		  for (j=0; j<Spe_Total_CNO[Qwan]; j++){

		    dEx += pref*(  CDM0[0][Mh_AN][kl][i][j]*Hx0[0][i][j].r
			         + CDM0[1][Mh_AN][kl][i][j]*Hx0[1][i][j].r);
		    dEy += pref*(  CDM0[0][Mh_AN][kl][i][j]*Hy0[0][i][j].r
				 + CDM0[1][Mh_AN][kl][i][j]*Hy0[1][i][j].r);
		    dEz += pref*(  CDM0[0][Mh_AN][kl][i][j]*Hz0[0][i][j].r
				 + CDM0[1][Mh_AN][kl][i][j]*Hz0[1][i][j].r);
		  }
		}
	      }

	      /* spin non-collinear with spin-orbit coupling or with LDA+U */

	      else if ( SpinP_switch==3 && (SO_switch==1 || (Hub_U_switch==1 && F_U_flag==1) 
                     || 1<=Constraint_NCS_switch || Zeeman_NCS_switch==1 || Zeeman_NCO_switch==1)){

                if (q_AN==h_AN){

		  for (i=0; i<Spe_Total_CNO[Hwan]; i++){
		    for (j=0; j<Spe_Total_CNO[Qwan]; j++){

		      dEx +=
			  CDM0[0][Mh_AN][kl][i][j]*Hx0[0][i][j].r
			- iDM0[0][Mh_AN][kl][i][j]*Hx0[0][i][j].i
			+ CDM0[1][Mh_AN][kl][i][j]*Hx0[1][i][j].r
			- iDM0[1][Mh_AN][kl][i][j]*Hx0[1][i][j].i
		    + 2.0*CDM0[2][Mh_AN][kl][i][j]*Hx0[2][i][j].r
		    - 2.0*CDM0[3][Mh_AN][kl][i][j]*Hx0[2][i][j].i;

		      dEy += 
			  CDM0[0][Mh_AN][kl][i][j]*Hy0[0][i][j].r
			- iDM0[0][Mh_AN][kl][i][j]*Hy0[0][i][j].i
			+ CDM0[1][Mh_AN][kl][i][j]*Hy0[1][i][j].r
			- iDM0[1][Mh_AN][kl][i][j]*Hy0[1][i][j].i
		    + 2.0*CDM0[2][Mh_AN][kl][i][j]*Hy0[2][i][j].r
		    - 2.0*CDM0[3][Mh_AN][kl][i][j]*Hy0[2][i][j].i; 

		      dEz +=  
			  CDM0[0][Mh_AN][kl][i][j]*Hz0[0][i][j].r
			- iDM0[0][Mh_AN][kl][i][j]*Hz0[0][i][j].i
			+ CDM0[1][Mh_AN][kl][i][j]*Hz0[1][i][j].r
			- iDM0[1][Mh_AN][kl][i][j]*Hz0[1][i][j].i
		    + 2.0*CDM0[2][Mh_AN][kl][i][j]*Hz0[2][i][j].r
		    - 2.0*CDM0[3][Mh_AN][kl][i][j]*Hz0[2][i][j].i; 

		    }
		  }
		}

                else {

		  for (i=0; i<Spe_Total_CNO[Hwan]; i++){  /* Hwan */
		    for (j=0; j<Spe_Total_CNO[Qwan]; j++){ /* Qwan  */

		      dEx += 
			  CDM0[0][Mh_AN][kl ][i][j]*Hx0[0][i][j].r
			- iDM0[0][Mh_AN][kl ][i][j]*Hx0[0][i][j].i
			+ CDM0[1][Mh_AN][kl ][i][j]*Hx0[1][i][j].r
			- iDM0[1][Mh_AN][kl ][i][j]*Hx0[1][i][j].i
		    + 2.0*CDM0[2][Mh_AN][kl ][i][j]*Hx0[2][i][j].r
		    - 2.0*CDM0[3][Mh_AN][kl ][i][j]*Hx0[2][i][j].i;

		      dEy += 
			  CDM0[0][Mh_AN][kl ][i][j]*Hy0[0][i][j].r
			- iDM0[0][Mh_AN][kl ][i][j]*Hy0[0][i][j].i
			+ CDM0[1][Mh_AN][kl ][i][j]*Hy0[1][i][j].r
			- iDM0[1][Mh_AN][kl ][i][j]*Hy0[1][i][j].i
	 	    + 2.0*CDM0[2][Mh_AN][kl ][i][j]*Hy0[2][i][j].r
		    - 2.0*CDM0[3][Mh_AN][kl ][i][j]*Hy0[2][i][j].i; 

		      dEz +=
			  CDM0[0][Mh_AN][kl ][i][j]*Hz0[0][i][j].r
			- iDM0[0][Mh_AN][kl ][i][j]*Hz0[0][i][j].i
			+ CDM0[1][Mh_AN][kl ][i][j]*Hz0[1][i][j].r
			- iDM0[1][Mh_AN][kl ][i][j]*Hz0[1][i][j].i
	 	    + 2.0*CDM0[2][Mh_AN][kl ][i][j]*Hz0[2][i][j].r
		    - 2.0*CDM0[3][Mh_AN][kl ][i][j]*Hz0[2][i][j].i;

		    } /* j */
		  } /* i */

		  w->calls++;dHNL(0,Mc_AN,q_AN,h_AN,DS_NL,Hx1,Hy1,Hz1);
		  kl1 = RMI1[Mc_AN][q_AN][h_AN];

		  for (i=0; i<Spe_Total_CNO[Qwan]; i++){ /* Qwan */
		    for (j=0; j<Spe_Total_CNO[Hwan]; j++){ /* Hwan */

		      dEx += 
			  CDM0[0][Mq_AN][kl1][i][j]*Hx1[0][i][j].r
			- iDM0[0][Mq_AN][kl1][i][j]*Hx1[0][i][j].i
			+ CDM0[1][Mq_AN][kl1][i][j]*Hx1[1][i][j].r
			- iDM0[1][Mq_AN][kl1][i][j]*Hx1[1][i][j].i
		    + 2.0*CDM0[2][Mq_AN][kl1][i][j]*Hx1[2][i][j].r
		    - 2.0*CDM0[3][Mq_AN][kl1][i][j]*Hx1[2][i][j].i;

		      dEy += 
			  CDM0[0][Mq_AN][kl1][i][j]*Hy1[0][i][j].r
			- iDM0[0][Mq_AN][kl1][i][j]*Hy1[0][i][j].i
			+ CDM0[1][Mq_AN][kl1][i][j]*Hy1[1][i][j].r
			- iDM0[1][Mq_AN][kl1][i][j]*Hy1[1][i][j].i
		    + 2.0*CDM0[2][Mq_AN][kl1][i][j]*Hy1[2][i][j].r
		    - 2.0*CDM0[3][Mq_AN][kl1][i][j]*Hy1[2][i][j].i; 

		      dEz +=
			  CDM0[0][Mq_AN][kl1][i][j]*Hz1[0][i][j].r
			- iDM0[0][Mq_AN][kl1][i][j]*Hz1[0][i][j].i
			+ CDM0[1][Mq_AN][kl1][i][j]*Hz1[1][i][j].r
			- iDM0[1][Mq_AN][kl1][i][j]*Hz1[1][i][j].i
		    + 2.0*CDM0[2][Mq_AN][kl1][i][j]*Hz1[2][i][j].r
		    - 2.0*CDM0[3][Mq_AN][kl1][i][j]*Hz1[2][i][j].i;

		    } /* j */
		  } /* i */

                }
	      }

	    } /* if (po2==1) */
	  } /* j0 */             

	  /* force from #4B */

	  Gxyz[Gc_AN][41] += dEx;
	  Gxyz[Gc_AN][42] += dEy;
	  Gxyz[Gc_AN][43] += dEz;

	  /* timing */
	  dtime(&Etime_atom);
	  time_per_atom[Gc_AN] += Etime_atom - Stime_atom;

	
}
