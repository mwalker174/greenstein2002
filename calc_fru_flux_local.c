/*       ----------------------------------------------------

	 NOTICE OF COPYRIGHT AND OWNERSHIP OF SOFTWARE
	 Copyright 2003, The Johns Hopkins University
	 School of Medicine. All rights reserved.

	 Name of Program: Local Control Model 
	 Version: Documented Version, C
	 Date: November 2003

	 --------------------------------------------------

	 calc_fru_flux_local.c - This subroutine is called by rk4am preceding 
	                         each call to fcn, and returns the values of 
							 fluxes which cross release unit boundaries.  
							 This routine sums over all the units to calculate 
							 these fluxes (Jxfer,Jtr,ICa,Ito2)

	 A lot of work is moved to mpi_master
*/


#include <math.h>

#include "parameters.h"
//#include "parameters_fcn_fru.h"

void calc_fru_flux_local(int NFRU_loc,
			 int LType_state_loc[MAX_LOCAL_NFRU][Nclefts_FRU][Nindepstates_LType],
			 int LCCPhosph_state_loc[MAX_LOCAL_NFRU][Nclefts_FRU],
			 int Ito2_state_loc[MAX_LOCAL_NFRU][Nclefts_FRU],
			 double FRU_states_loc[MAX_LOCAL_NFRU][Nstates_FRU],
			 double num[7])
{
  double ICa_numerator1;
  double ICa_numerator2;
  double sum_CaSS_local;
  double sum_CaJSR_local;
  
  int iFRU, icleft;
  double OCa_numerator1;
  double OCa_numerator2;
  int NIto2_Open;

  ICa_numerator1 = 0.0;
  ICa_numerator2 = 0.0;

  OCa_numerator1 = 0;
  OCa_numerator2 = 0;
  NIto2_Open = 0;
  for(iFRU = 0;iFRU<NFRU_loc;iFRU++) { 
    for(icleft = 0;icleft<Nclefts_FRU;icleft++) {
        if ( ((LType_state_loc[iFRU][icleft][index_LCC_states]==O1_LType)||(LType_state_loc[iFRU][icleft][index_LCC_states]==O2_LType)) 
	   && (LType_state_loc[iFRU][icleft][index_LCC_Vinact]==Oy_LType)) 
	{
	    //ICa_numerator = ICa_numerator + FRU_states_loc[iFRU][icleft+1]*exp_2VFRT-Cao*0.341; 
	  if (LCCPhosph_state_loc[iFRU][icleft]==6)
	  {
	    ICa_numerator2 = ICa_numerator2 + FRU_states_loc[iFRU][icleft+1];
	    OCa_numerator2 = OCa_numerator2 + 1;
	  } else {
	    ICa_numerator1 = ICa_numerator1 + FRU_states_loc[iFRU][icleft+1];
	    OCa_numerator1 = OCa_numerator1 + 1;
	  }
	}
      
      if (Ito2_state_loc[iFRU][icleft]==O_Ito2) 
		{
		NIto2_Open = NIto2_Open + 1;
	}
    }
  }

  sum_CaSS_local = 0.0;
  sum_CaJSR_local = 0.0;
  
  for(iFRU = 0;iFRU<NFRU_loc;iFRU++) {
    sum_CaJSR_local += FRU_states_loc[iFRU][index_frustates_CaJSR];
    for(icleft = 0;icleft<Nclefts_FRU;icleft++) {
      sum_CaSS_local += FRU_states_loc[iFRU][icleft+1];
    }
  }

  num[0] = sum_CaSS_local;
  num[1] = sum_CaJSR_local;
  num[2] = ICa_numerator1;
  num[3] = OCa_numerator1;
  num[4] = NIto2_Open;
  num[5] = ICa_numerator2;
  num[6] = OCa_numerator2;
}








