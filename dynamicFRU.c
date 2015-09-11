/*       ----------------------------------------------------

NOTICE OF COPYRIGHT AND OWNERSHIP OF SOFTWARE
Copyright 2003, The Johns Hopkins University
School of Medicine. All rights reserved.

Name of Program: Local Control Model 
Version: Documented Version, C
Date: November 2003

--------------------------------------------------

dynamicFRU.c - dynamic switching of number of FRUs 
used in computations

*/

#include "parameters.h"
#include <math.h>

//    Subroutine dynamicFRU takes care of changing number FRU:s used in computations.
// It also checks if change is needed or not.
// 
// 
// 
void dynamicFRU(const double t,
				double state[N],
				double FRU_states[NFRU_sim_max][Nstates_FRU],
				int LType_state[NFRU_sim_max][Nclefts_FRU][Nindepstates_LType],
				int RyR_state[NFRU_sim_max][Nclefts_FRU][NRyRs_per_cleft],
				int CaMKII_state[NFRU_sim_max][Nclefts_FRU][Nmon_per_holo],
				int LCCPhosph_state[NFRU_sim_max][Nclefts_FRU],
				int Ito2_state[NFRU_sim_max][Nclefts_FRU],
				int mti[NFRU_sim_max],
				unsigned long mt[NFRU_sim_max][mtN+1],
				const double PRyROpen,
				const double LCaOpenP)
{	
	double Voltage,sum;
	int change_phase;
	int iFRU,i,j,icleft;

	//      Change number of FRUs to be calculated

	Voltage = state[index_V];

	//        we want to use ratio of open L-type Ca channels
	//        i.e. ICa_Numerator/NFRU_sim=P_open_LCa in here
	//        using V for test purposes

	change_phase=0;
	if (t-time_of_last_change>1.0) {
		double glob_time_from_stim;

		// time to next stimulus
		glob_time_to_stim = shift-fmod(t,period);
		if (glob_time_to_stim<0.0) { 
			glob_time_to_stim=glob_time_to_stim+period;
		}
		glob_time_from_stim = fmod(t,period)-shift;
		if (glob_time_from_stim<0.0) { 
			glob_time_from_stim=glob_time_from_stim+period;
		}

		if (current_phase==0) { // now in low FRU # phase

			// if stimulus is coming within 100ms, switch to high FRU number
			// Here we assume that tstep is less than 100ms
			if ((glob_time_to_stim>=0.0)&&(glob_time_to_stim<100.0)) {
				change_phase=1;
			}

			// if L-type Ca channel open probability > threshold
			//  { switch to high FRU number
			if ((LCaOpenP>2e-6)||(PRyROpen>6.e-3)) {
				change_phase=1;
			}

		} else {		// in high FRU # phase

			// if membrane potential lower than -80mV 
			//    and open probability of L-type Ca channel is less that 1e-6 
			//    and ratio of open RyRs if lower than 2e-3 
			//  { switch to low FRU number
			if ((Voltage < -80)&&(LCaOpenP<1.e-6)&&(PRyROpen<2.e-3)) {
				change_phase=1;
			}

			// if time from last change is less than 1ms { for(not switch
			if (t-time_of_last_change<1.0) {
				change_phase=0;
			}

			// if time to next stimulus is less than 100ms { do not switch
			if ((glob_time_to_stim>=0.0)&&(glob_time_to_stim<100.0)) {
				change_phase=0;
			}

			// if time from previous stimulus is less than 50ms { do not switch
			if ((glob_time_from_stim>=0.0)&&(glob_time_from_stim<50.0)) {
				change_phase=0;
			}
		}
	}

	if (change_phase) { // switch FRU number
		time_of_last_change=t;
		if (current_phase==0) { // switching from low FRU number to high
			if (NFRU_sim_factor>1) {
				parallel_get_FRUs(FRU_states,LType_state,RyR_state,CaMKII_state,LCCPhosph_state,Ito2_state,mt,mti);

				for(j=0;j<NFRU_sim_factor-1;j++) {
					for(iFRU=0;iFRU<NFRU_sim_low;iFRU++) {
						for(i=0;i<Nstates_FRU;i++) {
							FRU_states[j*NFRU_sim_low+iFRU][i] = FRU_states[iFRU][i];
						}
						for(icleft = 0;icleft<Nclefts_FRU;icleft++) {
							for(i=0;i<NRyRs_per_cleft;i++) {
								RyR_state[j*NFRU_sim_low+iFRU][icleft][i] = RyR_state[iFRU][icleft][i];
							}
							LType_state[j*NFRU_sim_low+iFRU][icleft][index_LCC_states] = LType_state[iFRU][icleft][index_LCC_states];
							LType_state[j*NFRU_sim_low+iFRU][icleft][index_LCC_Vinact] = LType_state[iFRU][icleft][index_LCC_Vinact];
							for(i=0;i<Nmon_per_holo;i++){
								CaMKII_state[j*NFRU_sim_low+iFRU][icleft][i] = CaMKII_state[iFRU][icleft][i];
							}
							LCCPhosph_state[j*NFRU_sim_low+iFRU][icleft] = LCCPhosph_state[iFRU][icleft];
							Ito2_state[j*NFRU_sim_low+iFRU][icleft] = Ito2_state[iFRU][icleft];
						}
					}
				}

				// Redistribute data to slave processes
				initialize_mpi_state(FRU_states,LType_state,RyR_state,CaMKII_state,LCCPhosph_state,Ito2_state,mti,mt);
			}

			NFRU_sim=NFRU_sim_high;
			NFRU_scale=NFRU_scale_high;
			current_phase=1; // phase high FRU number
		} else {		// switching from high FRU number to low

			// average FRU Ca states to keep total Ca constant
			for(iFRU=0;iFRU<NFRU_sim_low;iFRU++) {
				for(i=0;i<Nstates_FRU;i++) {
					sum=0.0; // calculate average
					for(j=0;j<NFRU_sim_factor-1;j++) { // ok?
						sum=sum+FRU_states[iFRU+NFRU_sim_low*j][i];
					}
					FRU_states[iFRU][i] = sum/(double)(NFRU_sim_factor);
				}
			}

			NFRU_sim=NFRU_sim_low;
			NFRU_scale=NFRU_scale_low;
			current_phase=0; // phase low FRU number
		}
	}

}

