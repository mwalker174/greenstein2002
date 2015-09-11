/*       ----------------------------------------------------

	 NOTICE OF COPYRIGHT AND OWNERSHIP OF SOFTWARE
	 Copyright 2003, The Johns Hopkins University
	 School of Medicine. All rights reserved.

	 Name of Program: Local Control Model 
	 Version: Documented Version, C
	 Date: November 2003

	 --------------------------------------------------

	 calc_fru_avg_local.c - compute average FRU properties

*/

#ifdef USE_MPI
#include <mpi.h>
#endif

#include "parameters.h"
//#include "parameters_fcn_fru.h"

// Calculate additional output quantities that are algebraically 
// related to the state variables and stored in 'otherstate'

void calc_fru_avg_local(double num_stat[18],
			int NFRU_loc,
			double FRU_states_loc[MAX_LOCAL_NFRU][Nstates_FRU],
			int LType_state_loc[MAX_LOCAL_NFRU][Nclefts_FRU][Nindepstates_LType],
			int RyR_state_loc[MAX_LOCAL_NFRU][Nclefts_FRU][NRyRs_per_cleft],
			int CaMKII_state_loc[MAX_LOCAL_NFRU][Nclefts_FRU][Nmon_per_holo],
			int LCCPhosph_state_loc[MAX_LOCAL_NFRU][Nclefts_FRU],
			int Ito2_state_loc[MAX_LOCAL_NFRU][Nclefts_FRU])
{
  double CaSS,CaJSR,CaTOT_SS,CaTOT_JSR;
  double JRyR,PRyR_Open,PNorm_Mode,PnotVinact,PLType_Open,PIto2_Open;
  double CaMKII_Act;
  int CaMKIIStateTemp;
  double Act_coeff[8];
  int i,iFRU,icleft,NRyR_open,NRyR_ready;
  double LCCP0, LCCP1, LCCP2, LCCP3;
  double Mode2Open;
  double CaMKII_Phosph;

  CaSS = 0.0;
  CaJSR = 0.0;
  CaTOT_SS = 0.0;
  CaTOT_JSR = 0.0;
  JRyR = 0.0;
  PRyR_Open = 0.0;
  PNorm_Mode = 0.0;
  PnotVinact = 0.0;
  PLType_Open = 0.0;
  PIto2_Open = 0.0;
  NRyR_ready=0;
  CaMKII_Act = 0.0;
  CaMKIIStateTemp = 1;
  LCCP0 = 0.0;
  LCCP1 = 0.0;
  LCCP2 = 0.0;
  LCCP3 = 0.0;
  Mode2Open = 0.0;
  CaMKII_Phosph = 0.0;

  //Note: this vector is also defined in fru_rates_local.c
  Act_coeff[0]=0.0;
  Act_coeff[1]=0.0;
  Act_coeff[2]=0.75;
  Act_coeff[3]=1.0;
  Act_coeff[4]=0.8;
  Act_coeff[5]=0.8;
  Act_coeff[6]=0.17;
  Act_coeff[7]=0.0;

  for(iFRU = 0;iFRU<NFRU_loc;iFRU++) {
    CaJSR = CaJSR + FRU_states_loc[iFRU][index_frustates_CaJSR];
    CaTOT_JSR = CaTOT_JSR + FRU_states_loc[iFRU][index_frustates_CaJSR] * (1.0+CSQNtot/(KmCSQN+FRU_states_loc[iFRU][index_frustates_CaJSR]) );
    for(icleft = 0;icleft<Nclefts_FRU;icleft++) {
      CaSS  = CaSS + FRU_states_loc[iFRU][icleft+1];
      CaTOT_SS = CaTOT_SS +
	FRU_states_loc[iFRU][icleft+1]* (1.0+BSRtot/(KBSR+FRU_states_loc[iFRU][icleft+1])+BSLtot/(KBSL+FRU_states_loc[iFRU][icleft+1]));
      NRyR_open=0;
      for(i=0;i<NRyRs_per_cleft;i++) {
	if ((RyR_state_loc[iFRU][icleft][i]==O1_RyR)||(RyR_state_loc[iFRU][icleft][i]==O2_RyR)
	    ||(RyR_state_loc[iFRU][icleft][i]==O3_RyR)) 
	  {
	    NRyR_open = NRyR_open + 1;
	  }
	if ((RyR_state_loc[iFRU][icleft][i]==1)||(RyR_state_loc[iFRU][icleft][i]==2)) 
	  {
	    NRyR_ready = NRyR_ready + 1;
	  }
      }
      JRyR = JRyR + JRyRmax*((double)NRyR_open)*(FRU_states_loc[iFRU][index_frustates_CaJSR]-FRU_states_loc[iFRU][icleft+1]);
      PRyR_Open = PRyR_Open + ((double)NRyR_open);
      if (LType_state_loc[iFRU][icleft][index_LCC_states]<=6) 
	{
	  PNorm_Mode = PNorm_Mode + 1.0;
	}
      if (LType_state_loc[iFRU][icleft][index_LCC_Vinact]==Oy_LType) 
	{
	  PnotVinact = PnotVinact + 1.0;
	}
        if (((LType_state_loc[iFRU][icleft][index_LCC_states]==O1_LType)||(LType_state_loc[iFRU][icleft][index_LCC_states]==O2_LType))
	  &&(LType_state_loc[iFRU][icleft][index_LCC_Vinact]==Oy_LType)) 
       	{
	    PLType_Open = PLType_Open + 1.0;
		if (LCCPhosph_state_loc[iFRU][icleft]==6)
		{
		  Mode2Open = Mode2Open + 1.0;
		}
	}

      if (Ito2_state_loc[iFRU][icleft]==O_Ito2) 
	{
	  PIto2_Open = PIto2_Open + 1.0;
	}
      for(i=0;i<Nmon_per_holo;i++) {
	CaMKIIStateTemp = CaMKII_state_loc[iFRU][icleft][i];
	CaMKII_Act = CaMKII_Act + Act_coeff[CaMKIIStateTemp];
	//Clear CaMKIIStateTemp
	CaMKIIStateTemp = 1;
      }
      for(i=0;i<Nmon_per_holo;i++){
	if((CaMKII_state_loc[iFRU][icleft][i]==3)||(CaMKII_state_loc[iFRU][icleft][i]==4)||(CaMKII_state_loc[iFRU][icleft][i]==5)||(CaMKII_state_loc[iFRU][icleft][i]==6)){
		CaMKII_Phosph = CaMKII_Phosph + 1.0;
	}
      }
      if (LCCPhosph_state_loc[iFRU][icleft]==1)
	{
	  LCCP0 = LCCP0 + 1.0;
 	}
      if ((LCCPhosph_state_loc[iFRU][icleft]==2)||(LCCPhosph_state_loc[iFRU][icleft]==3))
	{
	  LCCP1 = LCCP1 + 1.0;
	}
      if ((LCCPhosph_state_loc[iFRU][icleft]==4)||(LCCPhosph_state_loc[iFRU][icleft]==5))
	{
	  LCCP2 = LCCP2 + 1.0;
	}
      if (LCCPhosph_state_loc[iFRU][icleft]==6)
	{
	  LCCP3 = LCCP3 + 1.0;
	}
    }
  }

  // The following values are sums, not averages
  num_stat[0]=CaSS;
  num_stat[1]=CaJSR;
  num_stat[2]=JRyR;
  num_stat[3]=PRyR_Open;
  num_stat[4]=(double)NRyR_ready;
  num_stat[5]=CaTOT_SS;
  num_stat[6]=CaTOT_JSR;
  num_stat[7]=PNorm_Mode;
  num_stat[8]=PnotVinact;
  num_stat[9]=PLType_Open;
  num_stat[10]=PIto2_Open;
  num_stat[11]=CaMKII_Act;
  num_stat[12]=LCCP0;
  num_stat[13]=LCCP1;
  num_stat[14]=LCCP2;
  num_stat[15]=LCCP3;
  num_stat[16]=CaMKII_Phosph;
  num_stat[17]=Mode2Open;
}


