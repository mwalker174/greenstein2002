/*       ----------------------------------------------------

	 NOTICE OF COPYRIGHT AND OWNERSHIP OF SOFTWARE
	 Copyright 2003, The Johns Hopkins University
	 School of Medicine. All rights reserved.

	 Name of Program: Local Control Model 
	 Version: Documented Version, C
	 Date: November 2003

	 --------------------------------------------------

	 initialize_state.c - Initialize states globally and for FRUs

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "parameters.h"
//#include "parameters_fcn_fru.h" 


void initialize_state(double state[N],
		      double FRU_states[NFRU_sim_max][Nstates_FRU],
		      int LType_state[NFRU_sim_max][Nclefts_FRU][Nindepstates_LType],
		      int RyR_state[NFRU_sim_max][Nclefts_FRU][NRyRs_per_cleft],
		      int CaMKII_state[NFRU_sim_max][Nclefts_FRU][Nmon_per_holo], 
		      int LCCPhosph_state[NFRU_sim_max][Nclefts_FRU],
		      int RyRPhosph_state[NFRU_sim_max][Nclefts_FRU][NRyRs_per_cleft],
		      int Ito2_state[NFRU_sim_max][Nclefts_FRU],
		      unsigned long mt[NFRU_sim_max][mtN+1],int mti[NFRU_sim_max])
{
	double a1,a2;
	int k, iFRU, icleft;
	FILE *ic_file;
	char filepath[256];

	int rewind_FRU;
	const int rewind_FRU_max=250;

	time_of_last_change=0.0;
	rewind_FRU=min(NFRU_sim,rewind_FRU_max);

	if (ic_file_flag) {
	  int RyR_NFRU,FRU_NFRU,LCh_NFRU,CaMKII_NFRU,LCCPhosph_NFRU,RyRPhosph_NFRU,Ito2_NFRU;
	  
	  // if initial conditions are provided in a file, open the 
	  // file and load the data

	  sprintf(filepath,"%s/%s",ic_dir,ic_RyR_file);
	  if ((ic_file=fopen(filepath,"rb"))==NULL) {
	    fprintf(stderr,"Cannot open file '%s'!\n",filepath);
	    abort();
	  }
	  RyR_NFRU=read_next_int(ic_file);
	  printf("NFRU in RyR file =%d\n",RyR_NFRU);
	  for(iFRU=0;iFRU<NFRU_sim;iFRU++) {
	    for(icleft=0;icleft<Nclefts_FRU;icleft++) {
	      for(k=0;k<NRyRs_per_cleft;k++) {
		RyR_state[iFRU][icleft][k]=read_next_int(ic_file);
	      }
	    }
	    if (((iFRU+1)%rewind_FRU)==0) {
	      rewind(ic_file);
	      RyR_NFRU=read_next_int(ic_file);
	    }
	  }
	  fclose(ic_file);
	  
	  sprintf(filepath,"%s/%s",ic_dir,ic_CaMKII_file);
	  if((ic_file=fopen(filepath,"rb"))==NULL){
		fprintf(stderr,"Cannot open file '%s'!\n",filepath);
		abort();
	  }
	  CaMKII_NFRU=read_next_int(ic_file);
	  printf("NFRU in CaMKII file =%d\n",CaMKII_NFRU);
	  for(iFRU=0;iFRU<NFRU_sim;iFRU++){
		for(icleft=0;icleft<Nclefts_FRU;icleft++){
			for(k=0;k<Nmon_per_holo;k++){
				CaMKII_state[iFRU][icleft][k]=read_next_int(ic_file);
			}
		}
		if (((iFRU+1)%rewind_FRU)==0) {
			rewind(ic_file);
			CaMKII_NFRU=read_next_int(ic_file);
		}
	  }
	  fclose(ic_file);
	  
	  sprintf(filepath,"%s/%s",ic_dir,ic_LCh_file);
	  if ((ic_file=fopen(filepath,"rb"))==NULL) {
	    fprintf(stderr,"Cannot open file '%s'!\n",filepath);
	    abort();
	  }
	  LCh_NFRU=read_next_int(ic_file);
	  printf("NFRU in LCh file =%d\n",LCh_NFRU);
	  for(iFRU=0;iFRU<NFRU_sim;iFRU++) {
	    for(icleft=0;icleft<Nclefts_FRU;icleft++) {
	      for(k=0;k<Nindepstates_LType;k++) {
		LType_state[iFRU][icleft][k]=read_next_int(ic_file);
	      }
	    }
	    if (((iFRU+1)%rewind_FRU)==0) {
	      rewind(ic_file);
	      RyR_NFRU=read_next_int(ic_file);
	    }
	  }
	  fclose(ic_file);

	  sprintf(filepath,"%s/%s",ic_dir,ic_LCCPhosph_file);
	  if ((ic_file=fopen(filepath,"rb"))==NULL) {
	    fprintf(stderr,"Cannot open file '%s'!\n",filepath);
	    abort();
	  }
	  LCCPhosph_NFRU=read_next_int(ic_file);
	  printf("NFRU in LCCPhosph file =%d\n",LCCPhosph_NFRU);
	  for(iFRU=0;iFRU<NFRU_sim;iFRU++){
	    for(icleft=0;icleft<Nclefts_FRU;icleft++){
	      LCCPhosph_state[iFRU][icleft]=read_next_int(ic_file);
	    }
	    if (((iFRU+1)%rewind_FRU)==0) {
	      rewind(ic_file);
	      RyR_NFRU=read_next_int(ic_file);
	    }
	  }
	  fclose(ic_file);

	  sprintf(filepath,"%s/%s",ic_dir,ic_RyRPhosph_file);
	  if ((ic_file=fopen(filepath,"rb"))==NULL) {
	    fprintf(stderr,"Cannot open file '%s'!\n",filepath);
	    abort();
	  }
	  RyRPhosph_NFRU=read_next_int(ic_file);
	  printf("NFRU in RyRPhosph file =%d\n",RyRPhosph_NFRU);
	  for(iFRU=0;iFRU<NFRU_sim;iFRU++) {
	    for(icleft=0;icleft<Nclefts_FRU;icleft++) {
	      for(k=0;k<NRyRs_per_cleft;k++) {
		RyRPhosph_state[iFRU][icleft][k]=read_next_int(ic_file);
	      }
	    }
	    if (((iFRU+1)%rewind_FRU)==0) {
	      rewind(ic_file);
	      RyRPhosph_NFRU=read_next_int(ic_file);
	    }
	  }
	  fclose(ic_file);
	  
	  
	  sprintf(filepath,"%s/%s",ic_dir,ic_Ito2_file);
	  if ((ic_file=fopen(filepath,"rb"))==NULL) {
	    fprintf(stderr,"Cannot open file '%s'!\n",filepath);
	    abort();
	  }
	  Ito2_NFRU=read_next_int(ic_file);
	  printf("NFRU in Ito2 file =%d\n",Ito2_NFRU);
	  for(iFRU=0;iFRU<NFRU_sim;iFRU++) {
	    for(icleft=0;icleft<Nclefts_FRU;icleft++) {
	      Ito2_state[iFRU][icleft]=read_next_int(ic_file);
	    }
	    if (((iFRU+1)%rewind_FRU)==0) {
	      rewind(ic_file);
	      RyR_NFRU=read_next_int(ic_file);
	    }
	  }
	  fclose(ic_file);
	  
	  sprintf(filepath,"%s/%s",ic_dir,ic_FRU_file);
	  if ((ic_file=fopen(filepath,"rb"))==NULL) {
	    fprintf(stderr,"Cannot open file '%s'!\n",filepath);
	    abort();
	  }
	  FRU_NFRU=read_next_int(ic_file);
	  printf("NFRU in FRUstates file =%d\n",FRU_NFRU);
	  for(iFRU=0;iFRU<NFRU_sim;iFRU++) {
	    for(k=0;k<Nstates_FRU;k++) {
	      FRU_states[iFRU][k]=read_next_double(ic_file);
	    }
	    if (((iFRU+1)%rewind_FRU)==0) {
	      rewind(ic_file);
	      RyR_NFRU=read_next_int(ic_file);
	    }
	  }
	  fclose(ic_file);
	  
	  sprintf(filepath,"%s/%s",ic_dir,ic_states_file);
	  if ((ic_file=fopen(filepath,"rb"))==NULL) {
	    fprintf(stderr,"Cannot open file '%s'!\n",filepath);
	    abort();
	  }
	  for(k=0;k<N;k++) {
	    state[k]=read_next_double(ic_file);
	  }
	  fclose(ic_file);
	  
	  if (use_seeds_file_flag) {
	    unsigned long buf[mtN+2];
	    FILE *seeds_file;
	    int i,numread;
	    
	    sprintf(filepath,"%s/%s",ic_dir,ic_seeds_file);
	    if ((seeds_file=fopen(filepath,"rb"))==NULL) {
	      fprintf(stderr,"Problem opening file %s\n",filepath);
	      abort();
	    }
	    
	    for(iFRU=0;iFRU<NFRU_sim;iFRU++) {		
	      numread=fread(buf,sizeof(unsigned long),mtN+2,seeds_file);
	      if (numread!=mtN+2) {
			fprintf(stderr,"Unable to read seeds from file '%s'\n",filepath);
	      }
	      mti[iFRU]=buf[1];
	      for(i=0;i<mtN;i++) {
			mt[iFRU][i]=buf[i+2];
	      }
				// A simple consistency check
	      if ((int)buf[0]!=iFRU) {
		fprintf(stderr,"RNG seed file fails consistency check: %ld!=%d\n",buf[0],iFRU);
		abort();
	      }
	    }
	    fclose(seeds_file);
	  }

	} else {
		// Otherwise, set inital conditions to default values 

		state[index_V]= -0.9543106E+02;
                state[index_mNa]= 0.2654309E-03;
                state[index_hNa]= 0.9985519E+00;
                state[index_jNa]= 0.9987711E+00;
		state[index_Nai]= 0.1000000E+02;	
		state[index_Ki] = 0.1554636E+03;
		state[index_Cai]= 0.7954739E-04;	
		state[index_CaNSR]=	0.176e0;
		state[index_xKs]=0.1501720E-03;
		state[index_LTRPNCa]=0.7485195E-01;
		state[index_HTRPNCa]=0.9748885E+00;
		state[index_C0Kv43]=0.590958E+00;
		state[index_C1Kv43]=0.155217E+00;
		state[index_C2Kv43]=0.152881E-01;
		state[index_C3Kv43]=0.669242E-03;
		state[index_OKv43]= 0.109861E-04;
		state[index_CI0Kv43]=0.220999E-00;
		state[index_CI1Kv43]=0.143513E-01;
		state[index_CI2Kv43]=0.596808E-03;
		state[index_CI3Kv43]=0.642013E-04;
		state[index_OIKv43]=0.184528E-02;
		state[index_C0Kv14]=0.7646312E+00;
		state[index_C1Kv14]=0.7393727E-01;
		state[index_C2Kv14]=0.2681076E-02;
		state[index_C3Kv14]=0.4321218E-04;
		state[index_OKv14]= 0.2620394E-06;
		state[index_CI0Kv14]=0.1493312E+00;
		state[index_CI1Kv14]=0.6441895E-02;
		state[index_CI2Kv14]=0.2012634E-02;
		state[index_CI3Kv14]=0.6128624E-03;
		state[index_OIKv14]= 0.3084292E-03;
		state[index_CaTOT]= 0.7417606475E+01;
		state[index_C1Herg]= 0.99e+00;
		state[index_C2Herg]= 0.8e-02;
		state[index_C3Herg]= 0.2e-02;
		state[index_OHerg]= 0.0e+00;
		state[index_IHerg]= 0.0e+00;
                state[index_CytCaMKII_I]= 1.0e+00;
                state[index_CytCaMKII_B]= 0.0e+00;
                state[index_CytCaMKII_P]= 0.0e+00;
                state[index_CytCaMKII_T]= 0.0e+00;
                state[index_CytCaMKII_A]= 0.0e+00;
                state[index_CytCaMKII_C]= 0.0e+00;
                state[index_CytCaMKII_U]= 0.0e+00;
                state[index_FracPLB_P]= 0.0e+00;


		for(iFRU = 0;iFRU<NFRU_sim;iFRU++) {		
			for(icleft = 0;icleft<Nclefts_FRU;icleft++) {			
				LType_state[iFRU][icleft][index_LCC_states] = 1;			
				LType_state[iFRU][icleft][index_LCC_Vinact] = 2;
				for(k = 0;k<NRyRs_per_cleft;k++) {
					RyR_state[iFRU][icleft][k] = 1;
				}
				for(k=0; k<Nmon_per_holo;k++){
					CaMKII_state[iFRU][icleft][k]=1;
				}
				LCCPhosph_state[iFRU][icleft] = 1;
				for(k = 0;k<NRyRs_per_cleft;k++) {
					RyRPhosph_state[iFRU][icleft][k] = 1;
				}
				Ito2_state[iFRU][icleft] = 2;
			}
			FRU_states[iFRU][index_frustates_CaJSR] = 0.7;   //  local JSR  
			FRU_states[iFRU][1] = 0.12e-3; 
			FRU_states[iFRU][2] = 0.12e-3;
			FRU_states[iFRU][3] = 0.12e-3;
			FRU_states[iFRU][4] = 0.12e-3; 
		}
	}

	if (vclamp_flag) {
	  state[index_V]=vclamp_hold;
	}

	if (iv_flag) {
	  state[index_V]=iv_clamp_hold;
	}

	// This section initalizes the calculation of total call Ca, 
	// which is used as one of the model algorithm verification tests

	a1 = 0.0;
	a2 = 0.0;
	for(iFRU = 0;iFRU<NFRU_sim;iFRU++) {
		a1 = a1 + FRU_states[iFRU][index_frustates_CaJSR] * (1.0+CSQNtot/(KmCSQN+FRU_states[iFRU][index_frustates_CaJSR]) );
		for(icleft = 0;icleft<Nclefts_FRU;icleft++) {
			a2 = a2 + FRU_states[iFRU][icleft+1]
			* (1.0+BSRtot/(KBSR+FRU_states[iFRU][icleft+1])+BSLtot/(KBSL+FRU_states[iFRU][icleft+1]));
		}
	}
	state[index_CaTOT]= 1.e6*(  (a1*VJSR + a2*VSS)*NFRU_scale
		+ state[index_CaNSR]*VNSR  
		+ state[index_Cai]*Vmyo 
		* (1.0 + CMDNtot/(KmCMDN+state[index_Cai]) + EGTAtot/(KmEGTA+state[index_Cai])) 
		+ (state[index_LTRPNCa]*LTRPNtot + state[index_HTRPNCa]*HTRPNtot)*Vmyo);  //picomoles


	// This section can be used to override selected initial conditions, typically used to 
	// set all LCCs and/or all RyRs to the closed resting state, also can be used to change
	// initial SR Ca load
	if (reset_RyRs_flag) {
		for(iFRU = 0;iFRU<NFRU_sim;iFRU++) {
			for(icleft = 0;icleft<Nclefts_FRU;icleft++) {			
				LType_state[iFRU][icleft][index_LCC_states] = 1;
				LType_state[iFRU][icleft][index_LCC_Vinact] = Oy_LType;
				for(k = 0;k<NRyRs_per_cleft;k++) {
					RyR_state[iFRU][icleft][k] = 1;
				}
				for(k = 0;k<NRyRs_per_cleft;k++) {
					RyRPhosph_state[iFRU][icleft][k] = 1;
				}
			}
		}
	}
	if (0) {
	  //if (1) {
		for(iFRU = 0;iFRU<NFRU_sim;iFRU++) {
			FRU_states[iFRU][index_frustates_CaJSR] = 0.3e0;
		}
		state[index_CaNSR]=	0.3e0;
	}
	//
	// Compute weighting factor for error scaling calculation in RK4 algorithm.
	// These are approx. 1/(max_abs_value) for each state.

	errweight[index_V] = 0; //1.e-2; // not independent
        errweight[index_mNa] = 1.0;
        errweight[index_hNa] = 1.0;
        errweight[index_jNa] = 1.0;
	errweight[index_Nai] = 0.1;
	errweight[index_Ki] = 1/140.0;
	errweight[index_Cai] = 1.e3;    // Ok?
	errweight[index_CaNSR] = 2.5; // Ok?
	errweight[index_xKs] = 1.0;    // Ok?
	errweight[index_LTRPNCa] = 1.0;
	errweight[index_HTRPNCa] = 1.0;
	errweight[index_C0Kv43] = 1.0;
	errweight[index_C1Kv43] = 1.0;
	errweight[index_C2Kv43] = 1.0;
	errweight[index_C3Kv43] = 1.0;
	errweight[index_OKv43] = 1.0;
	errweight[index_CI0Kv43] = 1.0;
	errweight[index_CI1Kv43] = 1.0;
	errweight[index_CI2Kv43] = 1.0;
	errweight[index_CI3Kv43] = 1.0;
	errweight[index_OIKv43] = 0.0; // 1.0, not independent
	errweight[index_C0Kv14] = 1.0;
	errweight[index_C1Kv14] = 1.0;
	errweight[index_C2Kv14] = 1.0;
	errweight[index_C3Kv14] = 1.0;
	errweight[index_OKv14] = 1.0;
	errweight[index_CI0Kv14] = 1.0;
	errweight[index_CI1Kv14] = 1.0;
	errweight[index_CI2Kv14] = 1.0;
	errweight[index_CI3Kv14] = 1.0;
	errweight[index_OIKv14] = 0.0; // 1.0, not independent
	errweight[index_CaTOT] = 0.1;
	errweight[index_C1Herg]= 1.0;
	errweight[index_C2Herg]= 1.0;
	errweight[index_C3Herg]= 1.0;
	errweight[index_OHerg]= 1.0;
	errweight[index_IHerg]= 0.0; // 1.0, not independent
        errweight[index_CytCaMKII_I]= 1.0;
        errweight[index_CytCaMKII_B]= 1.0;
        errweight[index_CytCaMKII_P]= 1.0;
        errweight[index_CytCaMKII_T]= 1.0;
        errweight[index_CytCaMKII_A]= 1.0;
        errweight[index_CytCaMKII_C]= 1.0;
        errweight[index_CytCaMKII_U]= 0.0; //not independent
        errweight[index_FracPLB_P]= 1.0;

	if (vclamp_flag) {
	  state[index_V]=vclamp_hold;
	}
	
}









