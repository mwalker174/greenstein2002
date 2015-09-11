/*       ----------------------------------------------------

	 NOTICE OF COPYRIGHT AND OWNERSHIP OF SOFTWARE
	 Copyright 2003, The Johns Hopkins University
	 School of Medicine. All rights reserved.
	 
	 Name of Program: Local Control Model 
	 Version: Documented Version, C
	 Date: November 2003
	 
	 --------------------------------------------------
	 
	 simfru_local.c - This subroutine runs the Monte Carlo simulation 
	 from time0 to timef for all of the stochastic elements 
	 within a single release unit (indexed by iFRU).
	 It runs in incremental time steps determined based 
	 on the transition rates of currently occupied states.  
	 On each step, the generation of random numbers are 
	 used to determine into which state each stochastic 
	 element has transitioned, if a transition has occured.
	 Local Ca2+ concentrations are integrated within these
	 time steps.
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "parameters.h"
//#include "parameters_fcn_fru.h"

void simfru_local(const double time0,const double timef,
				  double FRUdep_states0[Nstates_FRUdep],
				  double FRUdep_statesf[Nstates_FRUdep],
				  int LType_state[Nclefts_FRU][Nindepstates_LType],
				  int RyR_state[Nclefts_FRU][NRyRs_per_cleft], 
				  int CaMKII_state[Nclefts_FRU][Nmon_per_holo],
				  int LCCPhosph_state[Nclefts_FRU],
				  int Ito2_state[Nclefts_FRU],
				  double FRU_states[Nstates_FRU],
				  int *mti_loc,
				  unsigned long mt_loc[mtN+1])
{
	double time_interval;
	double FRU_states1[Nstates_FRU];
	double dFRU_states[Nstates_FRU];
	double k1[Nstates_FRU],y_1[Nstates_FRU];
	int icleft;
	double FRUdep_states_interval[Nstates_FRUdep];
	double FRUdep_states[Nstates_FRUdep];

	double LType_trans_rates[Nclefts_FRU][4];
	double LType_trans_probs[Nclefts_FRU][4];
	int LType_index[Nclefts_FRU][4];
	int LType_length[Nclefts_FRU];
	double LType_Vinact_exitrate[Nclefts_FRU];
	double LType_Vinact_exitprob[Nclefts_FRU];
	double Ito2_exitrate[Nclefts_FRU];
	double Ito2_exitprob[Nclefts_FRU];
	double CaMKII_trans_rates[Nclefts_FRU][Nmon_per_holo][4];
	double CaMKII_trans_probs[Nclefts_FRU][Nmon_per_holo][4];
	int CaMKII_index[Nclefts_FRU][Nmon_per_holo][4];
	int CaMKII_length[Nclefts_FRU][Nmon_per_holo];
	double LCCPhosph_trans_rates[Nclefts_FRU][4];
	double LCCPhosph_trans_probs[Nclefts_FRU][4];
	int LCCPhosph_index[Nclefts_FRU][4];
	int LCCPhosph_length[Nclefts_FRU];
	double RyR_trans_rates[Nclefts_FRU][NRyRs_per_cleft][4];
	double RyR_trans_probs[Nclefts_FRU][NRyRs_per_cleft][4];
	int RyR_index[Nclefts_FRU][NRyRs_per_cleft][4];
	int RyR_length[Nclefts_FRU][NRyRs_per_cleft];
	double max_exit_rate,Accum_Prob;
	double time_FRU,time_stepFRU;

#if 0
	double unidev_arr[NRVseqs_per_FRU*16]; // for mt19937
	double *unidev;

	int cycles_left;
	int min_cycles;
#else
	double unidev[NRVseqs_per_FRU]; // for mt19937
#endif

	int i,j,count,base_ind;
	int signal1[Nclefts_FRU];
	int count_cycle=0;

	double LType_open[Nclefts_FRU];
	double LType_PCa[Nclefts_FRU];
	int NRyR_open[Nclefts_FRU];
	double CaSStemp;

	time_FRU = time0;
	for(i=0;i<Nstates_FRUdep;i++) {
		FRUdep_states_interval[i] = FRUdep_statesf[i] - FRUdep_states0[i];
	}
	time_interval = timef-time0;

	FRUdep_states[index_frudep_V] = FRUdep_states0[index_frudep_V];
	FRUdep_states[index_frudep_Cai] = FRUdep_states0[index_frudep_Cai];
	FRUdep_states[index_frudep_CaNSR] = FRUdep_states0[index_frudep_CaNSR];
	// V       FRUdep_states[0]
	// Cai     FRUdep_states[1]
	// CaNSR   FRUdep_states[2]
	// CaJSR   FRU_states[0]
	// CaSS1   FRU_states[1]
	// CaSS2   FRU_states[2]
	// CaSS3   FRU_states[3]
	// CaSS4   FRU_states[4)

	for(icleft = 0;icleft<Nclefts_FRU;icleft++) {
		// set the indicator variables for channels that are open
		  if ( ((LType_state[icleft][index_LCC_states]==O1_LType)||(LType_state[icleft][index_LCC_states]==O2_LType)) 
			&& (LType_state[icleft][index_LCC_Vinact]==Oy_LType))  
		  { 
			LType_open[icleft] = 1.0;
		  } else {
			LType_open[icleft] = 0.0;
		  }

		// set the indicator variables for mode 2 channels
		  if (LCCPhosph_state[icleft]==6)
		  {
			LType_PCa[icleft] = PCa2;
		  } else {
			LType_PCa[icleft] = PCa1;
		  }

		NRyR_open[icleft]=0;
		for(i=0;i<NRyRs_per_cleft;i++) {
			if ((RyR_state[icleft][i]==O1_RyR)||(RyR_state[icleft][i]==O2_RyR) 
				|| (RyR_state[icleft][i]==O3_RyR)) {
					NRyR_open[icleft] = NRyR_open[icleft] + 1;
				}
		}
		signal1[icleft] = 0;
	}

#if 0
	min_cycles=floor(1.0+(timef-time_FRU)/10.001e-3);
	cycles_left=min(min_cycles,15);
	MersenneTwister_local(mti_loc,mt_loc,NRVseqs_per_FRU*cycles_left,unidev_arr);
#endif

	while (time_FRU<timef) {

		// fru_rates returns the rate of exit from each currently occupied
		// state into all neighboring states for each channel
		max_exit_rate=fru_rates_local(LType_state,RyR_state,CaMKII_state,LCCPhosph_state,Ito2_state,FRUdep_states,FRU_states,
			LType_trans_rates,LType_index,LType_length, LType_Vinact_exitrate,
			RyR_trans_rates,RyR_index,RyR_length,CaMKII_trans_rates,CaMKII_index,CaMKII_length,LCCPhosph_trans_rates,LCCPhosph_index,LCCPhosph_length,Ito2_exitrate,mti_loc,mt_loc);

		time_stepFRU = 0.1/max_exit_rate;
		// 	time_stepFRU = 0.05/max_exit_rate
		time_stepFRU = min(time_stepFRU,10.001e-3);
		// 	time_stepFRU = min(time_stepFRU,1.001e-3)
		CaSStemp = 0.e-3;
		for(icleft = 0;icleft<Nclefts_FRU;icleft++) {
			CaSStemp = max(CaSStemp,FRU_states[icleft+1]);
		}
		if (CaSStemp>1.5e-3) {
			time_stepFRU = min(time_stepFRU, 5.001e-3);
		}
		time_stepFRU = min(time_stepFRU, timef-time_FRU);
		// set time step based on the maximum exit rate

		for(i=0;i<Nstates_FRU;i++) {
			FRU_states1[i] = FRU_states[i];
		}

		// local velocity field for Ca concentrations within the unit
		fcn_fru(time_FRU,FRU_states1,FRUdep_states,LType_open,LType_PCa,NRyR_open,dFRU_states);

		for(i=0;i<Nstates_FRU;i++) { // First intermediate integration step (Trapezoidal Method)
			k1[i] = time_stepFRU*dFRU_states[i];
			y_1[i] = FRU_states1[i] + k1[i];
		}

		time_FRU = time_FRU + time_stepFRU;
		for(i=0;i<Nstates_FRUdep;i++) { // Interpolate the values of global states which are needed locally
			FRUdep_states[i] = FRUdep_states0[i]
			+ (time_FRU-time0)/time_interval*FRUdep_states_interval[i];
		}

		// event occurs at delta t
		// local velocity field for Ca concentrations within the unit
		fcn_fru(time_FRU,y_1,FRUdep_states,LType_open,LType_PCa,NRyR_open,dFRU_states);

		for(i=0;i<Nstates_FRU;i++) {	// Second intermediate integration step (Trapezoidal Method)
			FRU_states[i] = FRU_states[i] + (k1[i]+time_stepFRU*dFRU_states[i])/2.0;
		}

		// Calculate transition probabilities based on rates and time step
		for(icleft = 0;icleft<Nclefts_FRU;icleft++) {
			// There are always at least two possible transitions
			LType_trans_probs[icleft][0] = 1.0 - time_stepFRU*LType_trans_rates[icleft][0];
			LType_trans_probs[icleft][1] = time_stepFRU*LType_trans_rates[icleft][1];
			for(i=2;i<LType_length[icleft];i++) {
				LType_trans_probs[icleft][i] = time_stepFRU*LType_trans_rates[icleft][i];
			}

			LType_Vinact_exitprob[icleft] = 1.0 - time_stepFRU*LType_Vinact_exitrate[icleft];

			for(j=0;j<NRyRs_per_cleft;j++) {
				// There are always at least two possible transitions
				RyR_trans_probs[icleft][j][0] = 1.0 - time_stepFRU*RyR_trans_rates[icleft][j][0];
				RyR_trans_probs[icleft][j][1] = time_stepFRU*RyR_trans_rates[icleft][j][1];
				for(i=2;i<RyR_length[icleft][j];i++) {
					RyR_trans_probs[icleft][j][i] = time_stepFRU*RyR_trans_rates[icleft][j][i];
				}
			}

			for(j=0;j<Nmon_per_holo;j++){
				//There are always at least two possible transitions
				CaMKII_trans_probs[icleft][j][0] = 1.0 - time_stepFRU*CaMKII_trans_rates[icleft][j][0];
				CaMKII_trans_probs[icleft][j][1] = time_stepFRU*CaMKII_trans_rates[icleft][j][1];
				for(i=2;i<CaMKII_length[icleft][j];i++){
					CaMKII_trans_probs[icleft][j][i] = time_stepFRU*CaMKII_trans_rates[icleft][j][i];
				}
			}

			LCCPhosph_trans_probs[icleft][0] = 1.0 - time_stepFRU*LCCPhosph_trans_rates[icleft][0];
			LCCPhosph_trans_probs[icleft][1] = time_stepFRU*LCCPhosph_trans_rates[icleft][1];
			for(i=2;i<LCCPhosph_length[icleft];i++){
				LCCPhosph_trans_probs[icleft][i] = time_stepFRU*LCCPhosph_trans_rates[icleft][i];
			}
			Ito2_exitprob[icleft] = 1.0 - time_stepFRU*Ito2_exitrate[icleft];
		} // for

		// Generate uniform random variables for each stochastic element

#if 0
		if (cycles_left<1) {
			//min_cycles=floor(1.0+(timef-time_FRU)/time_stepFRU); // incorrect, may affect results
			min_cycles=floor(1.0+(timef-time_FRU)/10.001e-3);
			cycles_left=min(min_cycles,15);
			MersenneTwister_local(mti_loc,mt_loc,NRVseqs_per_FRU*cycles_left,unidev_arr);
		}

		unidev=&unidev_arr[NRVseqs_per_FRU*(cycles_left-1)];
		cycles_left--;
#else
		//MersenneTwister_local(mti_loc,mt_loc,NRVseqs_per_FRU,unidev);
		MersenneTwister_fast(mti_loc,mt_loc,NRVseqs_per_FRU,unidev);
#endif

		// Determine the actual transition that took place, if any
		for(icleft = 0;icleft<Nclefts_FRU;icleft++) {
			base_ind=icleft*NRVseqs_per_cleft;
			Accum_Prob = LType_trans_probs[icleft][0];
			if (unidev[base_ind]>Accum_Prob) {
				count = 1; // We know that state has been changed
				Accum_Prob += LType_trans_probs[icleft][1];
				while (unidev[base_ind]>Accum_Prob) {
					count++;
					Accum_Prob += LType_trans_probs[icleft][count];
				}

				LType_state[icleft][index_LCC_states] = LType_index[icleft][count];
			} // if

			if (unidev[1+base_ind]>=LType_Vinact_exitprob[icleft]) {
				LType_state[icleft][index_LCC_Vinact] = 3 - LType_state[icleft][index_LCC_Vinact];
			}

			for(i=0;i<NRyRs_per_cleft;i++) {
				Accum_Prob = RyR_trans_probs[icleft][i][0];
				if (unidev[i+2+base_ind]>Accum_Prob) {
					count = 1; // we know that state has been changed
					Accum_Prob += RyR_trans_probs[icleft][i][1];
					while (unidev[i+2+base_ind]>Accum_Prob) {
						count++;
						Accum_Prob += RyR_trans_probs[icleft][i][count];
					}

					RyR_state[icleft][i] = RyR_index[icleft][i][count];
					//THAPSIGARGIN CASE
					signal1[icleft] = 1;
					//signal1[icleft] = 0;
				}
				//THAPSIGARGIN CASE
				//RyR_state[icleft][i] = 1;
			} // for

			for(i=0;i<Nmon_per_holo;i++){
				Accum_Prob=CaMKII_trans_probs[icleft][i][0];
				if(unidev[i+2+NRyRs_per_cleft+base_ind]>Accum_Prob){
					count = 1; //we know that state has been changed
					Accum_Prob +=CaMKII_trans_probs[icleft][i][1];
					while (unidev[i+2+NRyRs_per_cleft+base_ind]>Accum_Prob){
						count++;
						Accum_Prob += CaMKII_trans_probs[icleft][i][count];
					}
					CaMKII_state[icleft][i] = CaMKII_index[icleft][i][count];
					//THIS IS ONLY FOR DEBUGGING PURPOSES
					//CaMKII_state[icleft][i] = 2;
				}
			} // for

			Accum_Prob=LCCPhosph_trans_probs[icleft][0];
			if(unidev[2+NRyRs_per_cleft+Nmon_per_holo+base_ind]>Accum_Prob){
				count = 1; //we know that state has been changed
				Accum_Prob +=LCCPhosph_trans_probs[icleft][1];
				while (unidev[2+NRyRs_per_cleft+Nmon_per_holo+base_ind]>Accum_Prob){
					count++;
					Accum_Prob += LCCPhosph_trans_probs[icleft][count];
				}
				LCCPhosph_state[icleft] = LCCPhosph_index[icleft][count];
			}

			if (unidev[(icleft+1)*NRVseqs_per_cleft-1]>=Ito2_exitprob[icleft]) {
				Ito2_state[icleft] = 3 - Ito2_state[icleft];
			} // if
		} // for

		// Reset indicator variables for which channels are open
		for(icleft = 0;icleft<Nclefts_FRU;icleft++) {
			  if ( ((LType_state[icleft][index_LCC_states]==O1_LType)||(LType_state[icleft][index_LCC_states]==O2_LType)) 
			     && (LType_state[icleft][index_LCC_Vinact]==Oy_LType)) 
			  {
				LType_open[icleft] = 1.0;
			  } else {
				LType_open[icleft] = 0.0;
			  } //if
			
			if (LCCPhosph_state[icleft]==6)
			{
				LType_PCa[icleft] = PCa2;
			} else {
				LType_PCa[icleft] = PCa1;
			} //if

			if (signal1[icleft]==1) {
				signal1[icleft] = 0;
				NRyR_open[icleft]=0;
				for(i = 0;i<NRyRs_per_cleft;i++) {
					if ((RyR_state[icleft][i]==O1_RyR)
						||(RyR_state[icleft][i]==O2_RyR) 
						||(RyR_state[icleft][i]==O3_RyR)) 
					{
						NRyR_open[icleft] = NRyR_open[icleft] + 1;
					}
				}
			} // if
		} // for
	}
}






