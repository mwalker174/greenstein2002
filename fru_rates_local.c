/* -------------------

	 NOTICE OF COPYRIGHT AND OWNERSHIP OF SOFTWARE
	 Copyright 2003, The Johns Hopkins University
	 School of Medicine. All rights reserved.

	 Name of Program: Local Control Model 
	 Version: Documented Version, C
	 Date: November 2003

	 --------------------------------------------------

	 fru_rates_local.c - This subroutine returns the individual rates
	                     for leaving the current state, the sum of those rates, 
			     and an index array that specifies the destination 
			     state that corresponds to each exit rate. All of 
			     these quantities are returned for each RyR and 
			     the LType channel and the V-dep inactivation gate 
			     of the LType channel and the Ito2 channel
			     Also returned is the maximum of the sum of rates
			     leaving any current state.
*/

#include <math.h>
#include <stdlib.h>

#include "parameters.h"


double fru_rates_local(int LType_state[Nclefts_FRU][Nindepstates_LType],
		       int RyR_state[Nclefts_FRU][NRyRs_per_cleft],
                       int CaMKII_state[Nclefts_FRU][Nmon_per_holo],
		       int LCCPhosph_state[Nclefts_FRU],
		       int Ito2_state[Nclefts_FRU],
		       const double FRUdep_states[Nstates_FRUdep],
		       const double FRU_states[Nstates_FRU],
		       double LType_rates[Nclefts_FRU][4],
		       int LType_index[Nclefts_FRU][4],
		       int LType_length[Nclefts_FRU], 
		       double LType_Vdep_exitrate[Nclefts_FRU],
		       double RyR_rates[Nclefts_FRU][NRyRs_per_cleft][4],
		       int RyR_index[Nclefts_FRU][NRyRs_per_cleft][4],
		       int RyR_length[Nclefts_FRU][NRyRs_per_cleft],
                       double CaMKII_rates[Nclefts_FRU][Nmon_per_holo][4],
                       int CaMKII_index[Nclefts_FRU][Nmon_per_holo][4],
                       int CaMKII_length[Nclefts_FRU][Nmon_per_holo],
		       double LCCPhosph_rates[Nclefts_FRU][4],
		       int LCCPhosph_index[Nclefts_FRU][4],
		       int LCCPhosph_length[Nclefts_FRU],
		       double Ito2_exitrate[Nclefts_FRU],
		       int *mti_loc,
		       unsigned long mt_loc[mtN])
{
	// real	rnum
	double dnum;
	double V,CaSS[Nclefts_FRU];
	double Act_coeff[8];
	Act_coeff[0]=0.0;
        // there is no state 0, so this is just a place holder
	Act_coeff[1]=0.0;
        // activity coefficient in state 1
	Act_coeff[2]=0.75;
        //activity coefficient in state 2
	Act_coeff[3]=1;
        //activity coefficient in state 3
	Act_coeff[4]=0.8;
        //activity coefficient in state 4
	Act_coeff[5]=0.8;
        //activity coefficient in state 5
	Act_coeff[6]=0.17;
        //activity coefficient in state 6
	Act_coeff[7]=0.0;
        //activity coefficient in state 7
	

	const double fL=0.85; // transition	rate into open state (1/ms)
	const double gL=2.0; //	transition rate	out	of open	state (1/ms)
	const double gPhosph = 0.049; //1; //0.049; //CHANGE MADE ON JAN 15, 2008 TO REFLECT CHANGE IN GATING WHEN CHANNEL IS PHOSPHORYLATED
	const double fLprime=0.005;	// transition rate into	Ca mode	open state (1/ms)
	const double gLprime=7.0; // transition	rate out of	Ca mode	open state (1/ms)
	const double bL=1.9356;	// mode	transition parameter
	const double bL2=bL*bL;
	const double bL3=bL*bL*bL;
	const double bL4=bL*bL*bL*bL;
	const double aL=2.0; //	mode transition	parameter
	const double aL2=aL*aL;
	const double aL3=aL*aL*aL;
	const double aL4=aL*aL*aL*aL;
	const double omega=0.83*2.0*1.3*0.01;  // mode transition parameter	(1/ms)

	const double alphacf=4.0*1.2*0.416;
	const double betacf=4.0*0.45*0.049;
	//change made on Jan 24 to eliminate CDI
	const double gammacf=0.83*1.9*1.3*0.31*7.5*0.09233;	// (ms-1 mM-1)

	const double LCCP_kalpha = 0.0008;//1.1*0.000729;
	const double LCCP_kbeta = 0.0008;//1.1*0.00026;

	double alpha, beta, alpha_prime, beta_prime, gamma_rate;
	double C0_to_C1, C1_to_C2, C2_to_C3, C3_to_C4;
	double C1_to_C0, C2_to_C1, C3_to_C2, C4_to_C3;
	double CCa0_to_CCa1, CCa1_to_CCa2;
	double CCa2_to_CCa3, CCa3_to_CCa4;
	double CCa1_to_CCa0, CCa2_to_CCa1;
	double CCa3_to_CCa2, CCa4_to_CCa3;
	double C0_to_CCa0, C1_to_CCa1, C2_to_CCa2;
	double C3_to_CCa3, C4_to_CCa4;

	const double CCa0_to_C0	= omega;		// = omega
	const double CCa1_to_C1	= omega/bL;	// = omega/bL
	const double CCa2_to_C2	= omega/bL2;	// = omega/bL^2
	const double CCa3_to_C3	= omega/bL3;	// = omega/bL^3
	const double CCa4_to_C4	= omega/bL4;	// = omega/bL^4


	double yCa_inf,	tau_yCa;
	const double yCa_frac=0.4;	// asymptotic value	for	fraction of	LCCs that
	// voltage-inactivate at depolarized potentials
	//	int	Ito2_state[Nclefts_FRU][NFRU_sim_max];
	//	double Ito2_exitrate[Nclefts_FRU];
	const double KdIto2=0.1502;	 //	(mM)
	const double kbIto2=2.0;  // (ms-1)
	const double kfIto2=kbIto2/KdIto2;	// (ms-1 mM-1)

	double max_rate;

	double k12,k23,k34,k54,k25,k56;	// (ms-1)
	const double k21=250.0,	k32=9.6, k43=0.07/0.06667*13.0,	k45=0.07;
	const double k52=0.001235, k65=30.0; //	at 10 CaNSR	rises
	const double k12cf=3000.0, k23cf=10.0*30000.0, k34cf=0.6*3000.0; //	(ms-1 mM-2)
	const double k54cf=0.6*0.198,k25cf=10.0*300.0,k56cf=2.0*4.0*3000.0;
	// parameters below	determine threshold	Ca levels for switching	between	different
	// RyR models based	on rapid equilibrium approximations
	double Sat_term;
	// const double	k_rate=12.0**4.0/16.3/(1000.0**2))	
	const double k_rate=0.00127215;	 
	double Sat_term2;
	// const double	k_rate2=43.0**4.0/(1000.0**2))
	const double k_rate2=3.4188;

	// const double	threshCa34to7=sqrt(k32*k_rate/(0.005*k34cf)) )	 //	Factor of 200 (0.005) difference in	rates 
	const double threshCa34to7=0.0368369379834969;	 //	Factor of 200 (0.005) difference in	rates 
	// const double	threshCa56to8=sqrt(k52*k_rate/(0.005*k56cf - k54cf)) )
	const double threshCa56to8=0.00011447933531005 ;

	const double threshMAX=2.0;	
	// const double	threshMAXCa	= sqrt(threshMAX*k_rate)) 
	const double threshMAXCa = 0.0504410547074504; 

	//CaMKII rate constants
	const double I_to_B = 0.00001; //nM^(-1)*ms^(-1)
	const double I_to_U = 0.0000006; //ms^(-1)
	//const double k_auto = 0.0005; //ms^(-1)
	//THIS WAS CHANGED ON JAN 14, 2008 TO RECTIFY PREVIOUS MISTAKE
	const double k_auto = 0.0008; //ms^(-1)
	const double B_to_I = 0.0008; //ms^(-1)
	const double P_to_T = 0.001; //ms^(-1)
	const double T_to_P = 0.001; //uM^(-4)*ms^(-1)
	const double T_to_A = 0.0000008; //B_to_I/1000
	const double A_to_T = 0.00001; //nM^(-1)*ms^(-1)
	const double A_to_C = 0.000146; //ms^(-1)

	int mon_on_R;
	int mon_on_L;
	int state_mon_on_R;
	int state_mon_on_L;
	double B_to_P_1;
	double B_to_P_2;
	double B_to_P;
	double P_to_B;
	double C_to_A;
	double C_to_U;
	double U_to_I; 
	double Ca4CaM[Nclefts_FRU];
	double Ca4[Nclefts_FRU];

	double LCCDephosph;
	double CaMKII_ActCleft;
	int CaMKIIStateTemp;

	double a1,a2,a3;
	int	i, icleft;

	V =	FRUdep_states[index_frudep_V];
	for(icleft = 0;icleft<Nclefts_FRU;icleft++) {
		CaSS[icleft] = FRU_states[icleft+1];
	}

	alpha =	alphacf	* exp(0.012*(V-35.0));
	beta = betacf *	exp(-0.05*(V-35.0));
	alpha_prime	= aL*alpha;
	beta_prime = beta/bL;

	// Computation of the rates	is moved into select case statement

	// C0_to_C1	= 4.0*alpha
	// C1_to_C2	= 3.0*alpha
	// C2_to_C3	= 2.0*alpha
	// C3_to_C4	=	   alpha

	// CCa0_to_CCa1	= 4.0*alpha_prime
	// CCa1_to_CCa2	= 3.0*alpha_prime
	// CCa2_to_CCa3	= 2.0*alpha_prime
	// CCa3_to_CCa4	=	   alpha_prime

	// C1_to_C0	=	   beta
	// C2_to_C1	= 2.0*beta
	// C3_to_C2	= 3.0*beta
	// C4_to_C3	= 4.0*beta

	// CCa1_to_CCa0	=	   beta_prime
	// CCa2_to_CCa1	= 2.0*beta_prime
	// CCa3_to_CCa2	= 3.0*beta_prime
	// CCa4_to_CCa3	= 4.0*beta_prime

	// CCa0_to_C0 =	omega		// = omega
	// CCa1_to_C1 =	omega/bL	// = omega/bL
	// CCa2_to_C2 =	omega/bL2	// = omega/bL^2
	// CCa3_to_C3 =	omega/bL3	// = omega/bL^3
	// CCa4_to_C4 =	omega/bL4	// = omega/bL^4

	for(icleft = 0;icleft<Nclefts_FRU; icleft++) {
	  gamma_rate =	  gammacf*CaSS[icleft];

	  // Computation of the rates	is moved into select case statement
	  //	   C0_to_CCa0 =	gamma_rate		// = gamma_rate
	  //	   C1_to_CCa1 =	aL*gamma_rate	// = gamma_rate*aL
	  //	   C2_to_CCa2 =	aL2*gamma_rate	// = gamma_rate*aL^2
	  //	   C3_to_CCa3 =	aL3*gamma_rate	// = gamma_rate*aL^3
	  //	   C4_to_CCa4 =	aL4*gamma_rate	// = gamma_rate*aL^4
	  
	  // The following goto statement	is used	as a CASE structure
	  
	  switch (LType_state[icleft][index_LCC_states]) {
	    
	  case 1:	// if LType	is in state	1
	    C0_to_C1 = 4.0*alpha;
	    C0_to_CCa0 = gamma_rate;
	    LType_rates[icleft][0] = (C0_to_C1+C0_to_CCa0);	// sum of rates	leaving	current	state
	    LType_rates[icleft][1] = C0_to_C1;		// rate	from 1 to 2
	    LType_rates[icleft][2] = C0_to_CCa0;		// rate	from 1 to 7
	    LType_length[icleft] = 3;			// length=1+#states	connected to current state
				//		LType_index[icleft][0] = 1			// index(1)	is not used
	    LType_index[icleft][1] = 2;			// index(2)	state into which rates(2) takes	you	= 2
	    LType_index[icleft][2] = 7;			// index(3)	state into which rates(3) takes	you	= 7
	    break;
	    
	  case 2:
	    C1_to_C2 = 3.0*alpha;
	    C1_to_C0 =		beta;
	    C1_to_CCa1 = aL*gamma_rate;	// = gamma_rate*aL
	    LType_rates[icleft][0] = (C1_to_C0+C1_to_C2+C1_to_CCa1);
	    LType_rates[icleft][1] = C1_to_C0;
	    LType_rates[icleft][2] = C1_to_C2;
	    LType_rates[icleft][3] = C1_to_CCa1;
	    LType_length[icleft] = 4;
				//		LType_index[icleft][0] = 2
	    LType_index[icleft][1] = 1;
	    LType_index[icleft][2] = 3;
	    LType_index[icleft][3] = 8;
	    break;
	    
	  case 3:
	    C2_to_C1 = 2.0*beta;
	    C2_to_C3 = 2.0*alpha;
	    C2_to_CCa2 = aL2*gamma_rate;
	    LType_rates[icleft][0] = (C2_to_C1+C2_to_C3+C2_to_CCa2);
	    LType_rates[icleft][1] = C2_to_C1;
	    LType_rates[icleft][2] = C2_to_C3;
	    LType_rates[icleft][3] = C2_to_CCa2;
	    LType_length[icleft] = 4;
	    //		LType_index[icleft][0] = 3
	    LType_index[icleft][1] = 2;
	    LType_index[icleft][2] = 4;
	    LType_index[icleft][3] = 9;
	    break;

	  case 4:
	    C3_to_C4 =		alpha;
	    C3_to_C2 = 3.0*beta;
	    C3_to_CCa3 = aL3*gamma_rate;
	    LType_rates[icleft][0] = (C3_to_C2+C3_to_C4+C3_to_CCa3);
	    LType_rates[icleft][1] = C3_to_C2;
	    LType_rates[icleft][2] = C3_to_C4;
	    LType_rates[icleft][3] = C3_to_CCa3;
	    LType_length[icleft] = 4;
				//		LType_index[icleft][0] = 4
	    LType_index[icleft][1] = 3;
	    LType_index[icleft][2] = 5;
	    LType_index[icleft][3] = 10;
	    break;
	    
	  case 5:
	    C4_to_C3 = 4.0*beta;
	    C4_to_CCa4 = aL4*gamma_rate;
	    LType_rates[icleft][0] = (C4_to_C3+fL+C4_to_CCa4);
	    LType_rates[icleft][1] = C4_to_C3;
	    LType_rates[icleft][2] = fL;
	    LType_rates[icleft][3] = C4_to_CCa4;
	    LType_length[icleft] = 4;
				//		LType_index[icleft][0] = 5
	    LType_index[icleft][1] = 4;
	    LType_index[icleft][2] = 6;
	    LType_index[icleft][3] = 11;
	    break;
	    
	  case 6:
	//CHANGES MADE ON JAN 15, 2008 TO REFLECT CHANGES IN GATING AFTER LCC PHOSPHORYLATION
	    if(LCCPhosph_state[icleft]==6)
		{
		LType_rates[icleft][0] = gL*gPhosph;
		LType_rates[icleft][1] = gL*gPhosph;
		}
	    else
		{
	    	LType_rates[icleft][0] = gL;
	    	LType_rates[icleft][1] = gL;
		}
	    LType_length[icleft] = 2;
				//		LType_index[icleft][0] = 6
	    LType_index[icleft][1] = 5;
	    break;
	    
	  case 7:
				//		CCa0_to_C0 = omega	// constant
	    CCa0_to_CCa1 = 4.0*alpha_prime;
	    LType_rates[icleft][0] = (CCa0_to_CCa1+CCa0_to_C0);
	    LType_rates[icleft][1] = CCa0_to_C0;
	    LType_rates[icleft][2] = CCa0_to_CCa1;
	    LType_length[icleft] = 3;
				//		LType_index[icleft][0] = 7
	    LType_index[icleft][1] = 1;
	    LType_index[icleft][2] = 8;
	    break;
	    
	  case 8:
				//		CCa1_to_C1 = omega/bL  // constant
	    CCa1_to_CCa2 = 3.0*alpha_prime;
	    CCa1_to_CCa0 =		beta_prime;
	    LType_rates[icleft][0] = (CCa1_to_CCa0+CCa1_to_CCa2+CCa1_to_C1);
	    LType_rates[icleft][1] = CCa1_to_CCa0;
	    LType_rates[icleft][2] = CCa1_to_C1;
	    LType_rates[icleft][3] = CCa1_to_CCa2;
	    LType_length[icleft] = 4;
				//		LType_index[icleft][0] = 8;
	    LType_index[icleft][1] = 7;
	    LType_index[icleft][2] = 2;
	    LType_index[icleft][3] = 9;
	    break;
	    
	  case 9:
				//		CCa2_to_C2 = omega/bL2	// constant
	    CCa2_to_CCa3 = 2.0*alpha_prime;
	    CCa2_to_CCa1 = 2.0*beta_prime;
	    LType_rates[icleft][0] = (CCa2_to_CCa1+CCa2_to_CCa3+CCa2_to_C2);
	    LType_rates[icleft][1] = CCa2_to_CCa1;
	    LType_rates[icleft][2] = CCa2_to_C2;
	    LType_rates[icleft][3] = CCa2_to_CCa3;
	    LType_length[icleft] = 4;
				//		LType_index[icleft][0] = 9;
	    LType_index[icleft][1] = 8;
	    LType_index[icleft][2] = 3;
	    LType_index[icleft][3] = 10;
	    break;
	    
	  case 10:		
				//		CCa3_to_C3 = omega/bL3	// constant
	    CCa3_to_CCa4 =		alpha_prime;
	    CCa3_to_CCa2 = 3.0*beta_prime;
	    LType_rates[icleft][0] = (CCa3_to_CCa2+CCa3_to_CCa4+CCa3_to_C3);
	    LType_rates[icleft][1] = CCa3_to_CCa2;
	    LType_rates[icleft][2] = CCa3_to_C3;
	    LType_rates[icleft][3] = CCa3_to_CCa4;
	    LType_length[icleft] = 4;
				//		LType_index[icleft][0] = 10
	    LType_index[icleft][1] = 9;
	    LType_index[icleft][2] = 4;
	    LType_index[icleft][3] = 11;
	    break;
	    
	  case 11:
				//		CCa4_to_C4 = omega/bL4	// constant
	    CCa4_to_CCa3 = 4.0*beta_prime;
	    LType_rates[icleft][0] = (CCa4_to_CCa3+CCa4_to_C4+fLprime);
	    LType_rates[icleft][1] = CCa4_to_CCa3;
	    LType_rates[icleft][2] = CCa4_to_C4;
	    LType_rates[icleft][3] = fLprime;
	    LType_length[icleft] = 4;
				//		LType_index[icleft][0] = 11
	    LType_index[icleft][1] = 10;
	    LType_index[icleft][2] = 5;
	    LType_index[icleft][3] = 12;
	    break;
	    
	  case 12:		
	    LType_rates[icleft][0] = gLprime;
	    LType_rates[icleft][1] = gLprime;
	    LType_length[icleft] = 2;
				//		LType_index[icleft][0] = 12
	    LType_index[icleft][1] = 11;
	    break;
	    
	  default: //	Unknown	state
	    fprintf(stderr,"Unknown state %d in LCC in cleft %d local\n",
		    LType_state[icleft][index_LCC_states],icleft);
	    break;
	  }
	} // for


	// Voltage dependent inactivation gate

	yCa_inf	= yCa_frac/(1.0+exp((V + 12.5)/5.0)) + (1.0-yCa_frac);
	tau_yCa	= 60.0 + 340.0/(1.0	+ exp((V+30.0)/12.0));

	for(icleft = 0;icleft<Nclefts_FRU;icleft++) {
	  switch (LType_state[icleft][index_LCC_Vinact]) {
	  case Oy_LType:
	    LType_Vdep_exitrate[icleft]	= (1.0-yCa_inf)/tau_yCa;
	    break;
	  case Cy_LType:
	    LType_Vdep_exitrate[icleft]	= yCa_inf/tau_yCa;
	    break;
	  default: //	Unknown	state
	    fprintf(stderr,"Unknown state %d in LCC Vinact in cleft %d local\n",
		    LType_state[icleft][index_LCC_Vinact],icleft);
	    break;
	  }
	  
	  /*		if (LType_state[icleft][index_LCC_Vinact]==Oy_LType)	{
			LType_Vdep_exitrate[icleft]	= (1.0-yCa_inf)/tau_yCa;
			} else {
			LType_Vdep_exitrate[icleft]	= yCa_inf/tau_yCa;
			} */
	}


	// RyR 

	for(icleft = 0;icleft<Nclefts_FRU;icleft++) {
		Sat_term = min(threshMAX,(CaSS[icleft]*CaSS[icleft])/k_rate);
		Sat_term2 =	min(threshMAX,(CaSS[icleft]*CaSS[icleft])/k_rate2);
		k12	= k12cf	* Sat_term2;
		k23	= k23cf	* Sat_term;
		k34	= k34cf	* Sat_term;
		k54	= k54cf	* Sat_term;
		k25	= k25cf	* Sat_term;
		k56	= k56cf	* Sat_term;

		if (min(CaSS[icleft],threshMAXCa)<threshCa56to8) { // use CaSS to determine	which model	to use for RyR

		  for(i=0;i<NRyRs_per_cleft;i++)	{
		    
		    switch (RyR_state[icleft][i])	{
		      
		    case 1:
		      RyR_rates[icleft][i][0]	= k12;
		      RyR_rates[icleft][i][1]	= k12;
		      RyR_length[icleft][i] = 2;
		      //				RyR_index[icleft][i][0]	= 1;
		      RyR_index[icleft][i][1]	= 2;
		      break;
		      
		    case 2:
		      RyR_rates[icleft][i][0]	= (k21+k23+k25);
		      RyR_rates[icleft][i][1]	= k21;
		      RyR_rates[icleft][i][2]	= k23;
		      RyR_rates[icleft][i][3]	= k25;
		      RyR_length[icleft][i] = 4;
		      //				RyR_index[icleft][i][0]	= 2;
		      RyR_index[icleft][i][1]	= 1;
		      RyR_index[icleft][i][2]	= 3;
		      RyR_index[icleft][i][3]	= 5;
		      break;
		      
		    case 3:
		      RyR_rates[icleft][i][0]	= (k32+k34);
		      RyR_rates[icleft][i][1]	= k32;
		      RyR_rates[icleft][i][2]	= k34;
		      RyR_length[icleft][i] = 3;
		      //				RyR_index[icleft][i][0]	= 3;
		      RyR_index[icleft][i][1]	= 2;
		      RyR_index[icleft][i][2]	= 4;
		      break;
						
		    case 4:
		      RyR_rates[icleft][i][0]	= (k43+k45);
		      RyR_rates[icleft][i][1]	= k43;
		      RyR_rates[icleft][i][2]	= k45;
		      RyR_length[icleft][i] = 3;
		      //				RyR_index[icleft][i][0]	= 4
		      RyR_index[icleft][i][1]	= 3;
		      RyR_index[icleft][i][2]	= 5;
		      break;
		      
		    case 5:
		      RyR_rates[icleft][i][0]	= (k52+k54+k56);
		      RyR_rates[icleft][i][1]	= k52;
		      RyR_rates[icleft][i][2]	= k54;
		      RyR_rates[icleft][i][3]	= k56;
		      RyR_length[icleft][i] = 4;
		      //				RyR_index[icleft][i][0]	= 5
		      RyR_index[icleft][i][1]	= 2;
		      RyR_index[icleft][i][2]	= 4;
		      RyR_index[icleft][i][3]	= 6;
		      break;
		      
		    case 6:
		      RyR_rates[icleft][i][0]	= k65;
		      RyR_rates[icleft][i][1]	= k65;
		      RyR_length[icleft][i] = 2;
		      //				RyR_index[icleft][i][0]	= 6
		      RyR_index[icleft][i][1]	= 5;
		      break;
		      
		    case 7:	// Not a real state, just transition to	state 4	or 3
		      //dnum=MersenneTwisterOne(iFRU);
		      dnum=MersenneTwisterOne_local(mti_loc,mt_loc);
		      if (dnum<(k34/(k34+k43))) {
			RyR_state[icleft][i] = 4;
			RyR_rates[icleft][i][0]	= (k43+k45);
			RyR_rates[icleft][i][1]	= k43;
			RyR_rates[icleft][i][2]	= k45;
			RyR_length[icleft][i] = 3;
			//						RyR_index[icleft][i][0]	= 4
			RyR_index[icleft][i][1]	= 3;
			RyR_index[icleft][i][2]	= 5;
		      } else {
			RyR_state[icleft][i] = 3;
			RyR_rates[icleft][i][0]	= (k32+k34);
			RyR_rates[icleft][i][1]	= k32;
			RyR_rates[icleft][i][2]	= k34;
			RyR_length[icleft][i] = 3;
			//						RyR_index[icleft][i][0]	= 3
			RyR_index[icleft][i][1]	= 2;
			RyR_index[icleft][i][2]	= 4;
		      }
		      break;
		      
		    case 8:	// Not a real state, just transition to	state 6	or 5
		      //dnum=MersenneTwisterOne(iFRU);
		      dnum=MersenneTwisterOne_local(mti_loc,mt_loc);
		      if (dnum< k56/(k56+k65)	) {
			RyR_state[icleft][i] = 6;
			RyR_rates[icleft][i][0]	= k65;
			RyR_rates[icleft][i][1]	= k65;
			RyR_length[icleft][i] = 2;
			//						RyR_index[icleft][i][0]	= 6
			RyR_index[icleft][i][1]	= 5;
		      } else {
			RyR_state[icleft][i] = 5;
			RyR_rates[icleft][i][0]	= (k52+k54+k56);
			RyR_rates[icleft][i][1]	= k52;
			RyR_rates[icleft][i][2]	= k54;
			RyR_rates[icleft][i][3]	= k56;
			RyR_length[icleft][i] = 4;
			//						RyR_index[icleft][i][0]	= 5
			RyR_index[icleft][i][1]	= 2;
			RyR_index[icleft][i][2]	= 4;
			RyR_index[icleft][i][3]	= 6;
		      }
		      break;
		      
		    default: //	Unknown	state
		      printf("Unknown state %d in RyR %d cleft % d\n",RyR_state[icleft][i],i,icleft);
		      break;
		    }
		  }

		} else {
		  if (min(CaSS[icleft],threshMAXCa)<threshCa34to7) { 
		    
		    for(i=0;i<NRyRs_per_cleft;i++)	{
		      
		      switch (RyR_state[icleft][i])	{
			
		      case 1:
			RyR_rates[icleft][i][0]	= k12;
			RyR_rates[icleft][i][1]	= k12;
			RyR_length[icleft][i] = 2;
			//				RyR_index[icleft][i][0]	= 1
			RyR_index[icleft][i][1]	= 2;
			break;
			
		      case 2:
			RyR_rates[icleft][i][0]	= (k21+k23+k25);
			RyR_rates[icleft][i][1]	= k21;
			RyR_rates[icleft][i][2]	= k23;
			RyR_rates[icleft][i][3]	= k25;
			RyR_length[icleft][i] = 4;
			//				RyR_index[icleft][i][0]	= 2
			RyR_index[icleft][i][1]	= 1;
			RyR_index[icleft][i][2]	= 3;
			RyR_index[icleft][i][3]	= 8;
			break;

		      case 3:
			RyR_rates[icleft][i][0]	= (k32+k34);
			RyR_rates[icleft][i][1]	= k32;
			RyR_rates[icleft][i][2]	= k34;
			RyR_length[icleft][i] = 3;
			//				RyR_index[icleft][i][0]	= 3
			RyR_index[icleft][i][1]	= 2;
			RyR_index[icleft][i][2]	= 4;
			break;
			
		      case 4:
			RyR_rates[icleft][i][0]	= (k43+k45);
			RyR_rates[icleft][i][1]	= k43;
			RyR_rates[icleft][i][2]	= k45;
			RyR_length[icleft][i] = 3;
			//				RyR_index[icleft][i][0]	= 4
			RyR_index[icleft][i][1]	= 3;
			RyR_index[icleft][i][2]	= 8;
			break;
			
		      case 5:	// Not a real state, just transition to	state 8
			RyR_state[icleft][i] = 8;		
			a1 = k65/(k56+k65);
			a2 = a1*k52;
			a3 = a1*k54;
			RyR_rates[icleft][i][0]	= (a2+a3);
			RyR_rates[icleft][i][1]	= a2;
			RyR_rates[icleft][i][2]	= a3;
			RyR_length[icleft][i] = 3;
			//				RyR_index[icleft][i][0]	= 8
			RyR_index[icleft][i][1]	= 2;
			RyR_index[icleft][i][2]	= 4;
			break;

		      case 6:	// Not a real state, just transition to	state 8
			RyR_state[icleft][i] = 8;
			a1 = k65/(k56+k65);
			a2 = a1*k52;
			a3 = a1*k54;
			RyR_rates[icleft][i][0]	= (a2+a3);
			RyR_rates[icleft][i][1]	= a2;
			RyR_rates[icleft][i][2]	= a3;
			RyR_length[icleft][i] = 3;
			//				RyR_index[icleft][i][0]	= 8;
			RyR_index[icleft][i][1]	= 2;
			RyR_index[icleft][i][2]	= 4;
			break;
			
		      case 7:	// Not a real state, just transition to	state 4	or state 3
			//dnum=MersenneTwisterOne(iFRU);
			dnum=MersenneTwisterOne_local(mti_loc,mt_loc);
			if (dnum< k34/(k34+k43)	) {
			  RyR_state[icleft][i] = 4;
			  RyR_rates[icleft][i][0]	= (k43+k45);
			  RyR_rates[icleft][i][1]	= k43;
			  RyR_rates[icleft][i][2]	= k45;
			  RyR_length[icleft][i] = 3;
			  //					RyR_index[icleft][i][1]	= 4;
			  RyR_index[icleft][i][1]	= 3;
			  RyR_index[icleft][i][2]	= 8;
			} else {
			  RyR_state[icleft][i] = 3;
			  RyR_rates[icleft][i][0]	= (k32+k34);
			  RyR_rates[icleft][i][1]	= k32;
			  RyR_rates[icleft][i][2]	= k34;
			  RyR_length[icleft][i] = 3;
			  //					RyR_index[icleft][i][1]	= 3;
			  RyR_index[icleft][i][1]	= 2;
			  RyR_index[icleft][i][2]	= 4;
			}
			break;
			
		      case 8:
			a1 = k65/(k56+k65);
			a2 = a1*k52;
			a3 = a1*k54;
			RyR_rates[icleft][i][0]	= (a2+a3);
			RyR_rates[icleft][i][1]	= a2;
			RyR_rates[icleft][i][2]	= a3;
			RyR_length[icleft][i] = 3;
			//				RyR_index[icleft][i][0]	= 8;
			RyR_index[icleft][i][1]	= 2;
			RyR_index[icleft][i][2]	= 4;
			break;
			
		      default: //	Unknown	state
			printf("Unknown state %d in RyR %d cleft % d local\n",RyR_state[icleft][i],i,icleft);
			break;
		      } // switch
		    } // for

		  } else {
		    
		    for(i=0;i<NRyRs_per_cleft;i++)	{
		      
		      switch (RyR_state[icleft][i])	{
		      case 1:
			RyR_rates[icleft][i][0]	= k12;
			RyR_rates[icleft][i][1]	= k12;
			RyR_length[icleft][i] = 2;
			//				RyR_index[icleft][i][0]	= 1;
			RyR_index[icleft][i][1]	= 2;
			break;
			
		      case 2:
			RyR_rates[icleft][i][0]	= (k21+k23+k25);
			RyR_rates[icleft][i][1]	= k21;
			RyR_rates[icleft][i][2]	= k23;
			RyR_rates[icleft][i][3]	= k25;
			RyR_length[icleft][i] = 4;
			//				RyR_index[icleft][i][0]	= 2;
			RyR_index[icleft][i][1]	= 1;
			RyR_index[icleft][i][2]	= 7;
			RyR_index[icleft][i][3]	= 8;
			break;
			
		      case 3:	// Not a real state, just transition to	state 8
			RyR_state[icleft][i] = 7;
			a1 = k34/(k43+k34);
			a2 = (1.0-a1)*k32;
			a3 = a1*k45;
			RyR_rates[icleft][i][0]	= (a2+a3);
			RyR_rates[icleft][i][1]	= a2;
			RyR_rates[icleft][i][2]	= a3;
			RyR_length[icleft][i] = 3;
			//				RyR_index[icleft][i][0]	= 7;
			RyR_index[icleft][i][1]	= 2;
			RyR_index[icleft][i][2]	= 8;
			break;
			
		      case 4:	// Not a real state, just transition to	state 7
			RyR_state[icleft][i] = 7;
			a1 = k34/(k43+k34);
			a2 = (1.0-a1)*k32;
			a3 = a1*k45;
			RyR_rates[icleft][i][0]	= (a2+a3);
			RyR_rates[icleft][i][1]	= a2;
			RyR_rates[icleft][i][2]	= a3;
			RyR_length[icleft][i] = 3;
			//				RyR_index[icleft][i][0]	= 7;
			RyR_index[icleft][i][1]	= 2;
			RyR_index[icleft][i][2]	= 8;
			break;
			
		      case 5:	// Not a real state, just transition to	state 8
			RyR_state[icleft][i] = 8;		
			a1 = k65/(k56+k65);
			a2 = a1*k52;
			a3 = a1*k54;
			RyR_rates[icleft][i][0]	= (a2+a3);
			RyR_rates[icleft][i][1]	= a2;
			RyR_rates[icleft][i][2]	= a3;
			RyR_length[icleft][i] = 3;
			//				RyR_index[icleft][i][0]	= 8;
			RyR_index[icleft][i][1]	= 2;
			RyR_index[icleft][i][2]	= 7;
			break;
			
		      case 6:	// Not a real state, just transition to	state 8
			RyR_state[icleft][i] = 8;
			a1 = k65/(k56+k65);
			a2 = a1*k52;
			a3 = a1*k54;
			RyR_rates[icleft][i][0]	= (a2+a3);
			RyR_rates[icleft][i][1]	= a2;
			RyR_rates[icleft][i][2]	= a3;
			RyR_length[icleft][i] = 3;
			//				RyR_index[icleft][i][0]	= 8;
			RyR_index[icleft][i][1]	= 2;
			RyR_index[icleft][i][2]	= 7;
			break;
			
		      case 7:
			a1 = k34/(k43+k34);
			a2 = (1.0-a1)*k32;
			a3 = a1*k45;
			RyR_rates[icleft][i][0]	= (a2+a3);
			RyR_rates[icleft][i][1]	= a2;
			RyR_rates[icleft][i][2]	= a3;
			RyR_length[icleft][i] = 3;
			//				RyR_index[icleft][i][0]	= 7;
			RyR_index[icleft][i][1]	= 2;
			RyR_index[icleft][i][2]	= 8;
			break;
			
		      case 8:
			a1 = k65/(k56+k65);
			a2 = a1*k52;
			a3 = a1*k54;
			RyR_rates[icleft][i][0]	= (a2+a3);
			RyR_rates[icleft][i][1]	= a2;
			RyR_rates[icleft][i][2]	= a3;
			RyR_length[icleft][i] = 3;
			//				RyR_index[icleft][i][0]	= 8;
			RyR_index[icleft][i][1]	= 2;
			RyR_index[icleft][i][2]	= 7;
			break;
			
		      default:
			fprintf(stderr,"Unknown state %d in RyR %d cleft % d local\n",RyR_state[icleft][i],i,icleft);
			break;
		      } // switch
		    } // for
		  } // if
		} // if
	} // for

	// CaMKII
	for(icleft = 0; icleft<Nclefts_FRU;icleft++){
	   Ca4[icleft] = CaSS[icleft]*CaSS[icleft]*CaSS[icleft]*CaSS[icleft]*(1.e+12);
	   Ca4CaM[icleft] = SSCaM*(1.e+6)*Ca4[icleft]/(Ca4[icleft]+ 1);

	   //Added on Mar 27, 2008
	   CaMKII_ActCleft = 0.0;
	   for(i=0; i<Nmon_per_holo; i++){
		CaMKIIStateTemp = CaMKII_state[icleft][i];
		CaMKII_ActCleft = CaMKII_ActCleft + Act_coeff[CaMKIIStateTemp];
		CaMKIIStateTemp = 1;
	   }
	   CaMKII_ActCleft = CaMKII_ActCleft/Nmon_per_holo;

	   for(i=0; i<Nmon_per_holo; i++){

		switch(CaMKII_state[icleft][i]){

		case 1: //if CaMKII is in state 1 (unbound)
			CaMKII_rates[icleft][i][0] = (I_to_B*Ca4CaM[icleft] + I_to_U); //sum of rates leaving current state
			CaMKII_rates[icleft][i][1] = I_to_B*Ca4CaM[icleft]; //rate from 1 to 2
			CaMKII_rates[icleft][i][2] = I_to_U; //rate from 1 to 7
			CaMKII_length[icleft][i] = 3;
				//CaMKII_index[icleft][i][0] = 1; //index(1) is not used
			CaMKII_index[icleft][i][1] = 2; //index (2) state into which rates (2) takes you = 2
			CaMKII_index[icleft][i][2] = 7; //index(3) state into which rates(3) takes you = 3
			break;

		case 2:
		 	if (i<(Nmon_per_holo/2)){
				mon_on_R = (i+1)%(Nmon_per_holo/2);
				mon_on_L = (i+Nmon_per_holo/2 -1)%(Nmon_per_holo/2);
			}else{
				mon_on_R = (i+1)%(Nmon_per_holo/2) + Nmon_per_holo/2;
				mon_on_L = (i + Nmon_per_holo/2-1)%(Nmon_per_holo/2) + Nmon_per_holo/2; 
			}  
			state_mon_on_R=CaMKII_state[icleft][mon_on_R];
			state_mon_on_L=CaMKII_state[icleft][mon_on_L];
			B_to_P_1=k_auto*Act_coeff[state_mon_on_R];
			B_to_P_2=k_auto*Act_coeff[state_mon_on_L];
			//CHANGE MADE ON MAR 27, 2008
			//B_to_P = B_to_P_1 + B_to_P_2;
			B_to_P =(CaMKII_ActCleft*CaMKII_ActCleft/(CaMKII_ActCleft*CaMKII_ActCleft + 1.3*1.3))*Act_coeff[2]*(B_to_P_1 + B_to_P_2);
			CaMKII_rates[icleft][i][0]=B_to_P + B_to_I;
			CaMKII_rates[icleft][i][1]=B_to_I;
			CaMKII_rates[icleft][i][2]=B_to_P;
			CaMKII_length[icleft][i]=3;
				//CaMKII_index[icleft][i][0]=2;
			CaMKII_index[icleft][i][1]=1;
			CaMKII_index[icleft][i][2]=3;
			break;
		case 3:
			P_to_B=VmaxPP1*PP1;
			CaMKII_rates[icleft][i][0]=P_to_T + P_to_B;
			CaMKII_rates[icleft][i][1]=P_to_B;
			CaMKII_rates[icleft][i][2]=P_to_T;
			CaMKII_length[icleft][i]=3;
				//CaMKII_index[icleft][i][0]=3;
			CaMKII_index[icleft][i][1]=2;
			CaMKII_index[icleft][i][2]=4;
			break;
		case 4:
			CaMKII_rates[icleft][i][0]=T_to_P*Ca4[icleft] + T_to_A;
			CaMKII_rates[icleft][i][1]=T_to_P*Ca4[icleft];
			CaMKII_rates[icleft][i][2]=T_to_A;
			CaMKII_length[icleft][i]=3;
				//CaMKII_index[icleft][i][0]=4;
			CaMKII_index[icleft][i][1]=3;
			CaMKII_index[icleft][i][2]=5;
			break;
		case 5:
			CaMKII_rates[icleft][i][0]=A_to_C + A_to_T*(SSCaM*(1.e+6)-Ca4CaM[icleft]);
			CaMKII_rates[icleft][i][1]=A_to_T*(SSCaM*(1.e+6)-Ca4CaM[icleft]);
			CaMKII_rates[icleft][i][2]=A_to_C;
			CaMKII_length[icleft][i]=3;
				//CaMKII_index[icleft][i][0]=5;
			CaMKII_index[icleft][i][1]=4;
			CaMKII_index[icleft][i][2]=6;
			break;
		case 6:
			C_to_A=VmaxPP2A*PP2A;
			C_to_U=VmaxPP1*PP1;
			CaMKII_rates[icleft][i][0]=C_to_A+C_to_U;
			CaMKII_rates[icleft][i][1]=C_to_A;
			CaMKII_rates[icleft][i][2]=C_to_U;
			CaMKII_length[icleft][i]=3;
				//CaMKII_index[icleft][i][0]=6;
			CaMKII_index[icleft][i][1]=5;
			CaMKII_index[icleft][i][2]=7;
			break;
		case 7:
			U_to_I=VmaxPP2A*PP2A;
			CaMKII_rates[icleft][i][0]=U_to_I;
			CaMKII_rates[icleft][i][1]=U_to_I;
			CaMKII_length[icleft][i]=2;
				//CaMKII_index[icleft][i][0]=7;
			CaMKII_index[icleft][i][1]=1;
			break;
		}
	   }
	}

	//LCCPhosph
	
	for(icleft = 0; icleft<Nclefts_FRU;icleft++){
		CaMKIIStateTemp = 1;
		CaMKII_ActCleft = 0.0;
		for(i=0; i<Nmon_per_holo; i++){
			CaMKIIStateTemp = CaMKII_state[icleft][i];
			CaMKII_ActCleft = CaMKII_ActCleft + Act_coeff[CaMKIIStateTemp];
			CaMKIIStateTemp = 1;
		}

		CaMKII_ActCleft = CaMKII_ActCleft/Nmon_per_holo;

		//LCCDephosph = VmaxPP2A*PP2A + VmaxPP1*PP1;
		LCCDephosph = VmaxPP2A*PP2A;

		switch (LCCPhosph_state[icleft]){
		
		case 1:
			LCCPhosph_rates[icleft][0] = CaMKII_ActCleft*(2*LCCP_kalpha + LCCP_kbeta);
			LCCPhosph_rates[icleft][1] = CaMKII_ActCleft*2*LCCP_kalpha;
			LCCPhosph_rates[icleft][2] = CaMKII_ActCleft*LCCP_kbeta;
			LCCPhosph_length[icleft]=3;
			LCCPhosph_index[icleft][1]=2;
			LCCPhosph_index[icleft][2]=3;
			break;
		case 2:
			LCCPhosph_rates[icleft][0] = LCCDephosph + CaMKII_ActCleft*(LCCP_kalpha + LCCP_kbeta);
			LCCPhosph_rates[icleft][1] = LCCDephosph;
			LCCPhosph_rates[icleft][2] = CaMKII_ActCleft*LCCP_kalpha;
			LCCPhosph_rates[icleft][3] = CaMKII_ActCleft*LCCP_kbeta;
			LCCPhosph_length[icleft] = 4;
			LCCPhosph_index[icleft][1] = 1;
			LCCPhosph_index[icleft][2] = 5;
			LCCPhosph_index[icleft][3] = 4;
			break;
		case 3:
			LCCPhosph_rates[icleft][0] = LCCDephosph + CaMKII_ActCleft*2*LCCP_kalpha;
			LCCPhosph_rates[icleft][1] = LCCDephosph;
			LCCPhosph_rates[icleft][2] = CaMKII_ActCleft*2*LCCP_kalpha;
			LCCPhosph_length[icleft] = 3;
			LCCPhosph_index[icleft][1] = 1;
			LCCPhosph_index[icleft][2] = 4;
			break;
		case 4:
			LCCPhosph_rates[icleft][0] = 2*LCCDephosph + CaMKII_ActCleft*LCCP_kalpha;
			LCCPhosph_rates[icleft][1] = LCCDephosph;
			LCCPhosph_rates[icleft][2] = LCCDephosph;
			LCCPhosph_rates[icleft][3] = CaMKII_ActCleft*LCCP_kalpha;
			LCCPhosph_length[icleft] = 4;
			LCCPhosph_index[icleft][1] = 3;
			LCCPhosph_index[icleft][2] = 2;
			LCCPhosph_index[icleft][3] = 6;
		case 5:
			LCCPhosph_rates[icleft][0] = 2*LCCDephosph + CaMKII_ActCleft*LCCP_kbeta;
			LCCPhosph_rates[icleft][1] = 2*LCCDephosph;
			LCCPhosph_rates[icleft][2] = CaMKII_ActCleft*LCCP_kbeta;
			LCCPhosph_length[icleft] = 3;
			LCCPhosph_index[icleft][1] = 2;
			LCCPhosph_index[icleft][2] = 6;
			break;
		case 6:
			LCCPhosph_rates[icleft][0] = 3*LCCDephosph;
			LCCPhosph_rates[icleft][1] = 2*LCCDephosph;
			LCCPhosph_rates[icleft][2] = LCCDephosph;
			LCCPhosph_length[icleft] = 3;
			LCCPhosph_index[icleft][1] = 4;
			LCCPhosph_index[icleft][2] = 5;
			break;
		}
	}
 



	// Ito2

	for(icleft = 0;icleft<Nclefts_FRU;icleft++) {
		switch (Ito2_state[icleft]) {
			case C_Ito2: // closed
				Ito2_exitrate[icleft] =	kfIto2*CaSS[icleft];
				break;
			case O_Ito2: // open
				Ito2_exitrate[icleft] =	kbIto2;
				break;
			default:
				fprintf(stderr,"Unknown state %d in Ito2 %d cleft local\n",Ito2_state[icleft],icleft);
				break;
		}

/*		if (Ito2_state[icleft]==1) {	// 1 is	closed,	2 is open
			Ito2_exitrate[icleft] =	kfIto2*CaSS[icleft];
		} else {
			Ito2_exitrate[icleft] =	kbIto2;
		} 
*/
	}

	max_rate = 0.0;
	for(icleft = 0;	icleft <Nclefts_FRU; icleft++)	{
		max_rate = max(max_rate,LType_rates[icleft][0]);	// Calc. max of	sum	of exit	rates
		max_rate = max(max_rate,LType_Vdep_exitrate[icleft]);
		max_rate = max(max_rate,Ito2_exitrate[icleft]);

		//	   for(i=1;i<=NRyRs_per_cleft
		//		max_rate = max(max_rate,RyR_rates[icleft][i][0])
		//	   }

		// The same	unrolled
		max_rate = max(max_rate,RyR_rates[icleft][0][0]);
		max_rate = max(max_rate,RyR_rates[icleft][1][0]);
		max_rate = max(max_rate,RyR_rates[icleft][2][0]);
		max_rate = max(max_rate,RyR_rates[icleft][3][0]);
		max_rate = max(max_rate,RyR_rates[icleft][4][0]);

		max_rate = max(max_rate,CaMKII_rates[icleft][0][0]);
		max_rate = max(max_rate,CaMKII_rates[icleft][1][0]);
		max_rate = max(max_rate,CaMKII_rates[icleft][2][0]);
		max_rate = max(max_rate,CaMKII_rates[icleft][3][0]);
		max_rate = max(max_rate,CaMKII_rates[icleft][4][0]);
		max_rate = max(max_rate,CaMKII_rates[icleft][5][0]);
		max_rate = max(max_rate,CaMKII_rates[icleft][6][0]);
		max_rate = max(max_rate,CaMKII_rates[icleft][7][0]);
		max_rate = max(max_rate,CaMKII_rates[icleft][8][0]);
		max_rate = max(max_rate,CaMKII_rates[icleft][9][0]);
		max_rate = max(max_rate,CaMKII_rates[icleft][10][0]);
		max_rate = max(max_rate,CaMKII_rates[icleft][11][0]);

		max_rate = max(max_rate,LCCPhosph_rates[icleft][0]);
	}

#if 0
	if (max_rate>10000) {
		printf("maxrate=%g\n",max_rate);
		for(icleft = 0;	icleft <Nclefts_FRU; icleft++)	{
			printf("%d: %g (%d) %g %g %g\n",icleft,LType_rates[icleft][0],
				LType_state[icleft][index_LCC_states],
				LType_Vdep_exitrate[icleft],Ito2_exitrate[icleft],CaSS[icleft]);
		}
		for(icleft = 0;	icleft <Nclefts_FRU; icleft++)	{
			printf("%d: %g %d %g %d %g %d %g %d %g %d\n",icleft,
				RyR_rates[icleft][0][0],RyR_state[icleft][0],
				RyR_rates[icleft][1][0],RyR_state[icleft][1],
				RyR_rates[icleft][2][0],RyR_state[icleft][2],
				RyR_rates[icleft][3][0],RyR_state[icleft][3],
				RyR_rates[icleft][4][0],RyR_state[icleft][4]);
		}
		for(icleft=0; icleft<Nclefts_FRU; icleft++){
			printf("%d: %g %d %g %d %g %d %g %d %g %d %g %d %g %d %g %d %g %d %g %d %g %d %g %d\n",icleft,
				CaMKII_rates[icleft][0][0],CaMKII_state[icleft][0],
				CaMKII_rates[icleft][1][0],CaMKII_state[icleft][1],
				CaMKII_rates[icleft][2][0],CaMKII_state[icleft][2],
				CaMKII_rates[icleft][3][0],CaMKII_state[icleft][3],
				CaMKII_rates[icleft][4][0],CaMKII_state[icleft][4],
				CaMKII_rates[icleft][5][0],CaMKII_state[icleft][5],
				CaMKII_rates[icleft][6][0],CaMKII_state[icleft][6],
				CaMKII_rates[icleft][7][0],CaMKII_state[icleft][7],
				CaMKII_rates[icleft][8][0],CaMKII_state[icleft][8],
				CaMKII_rates[icleft][9][0],CaMKII_state[icleft][9],
				CaMKII_rates[icleft][10][0],CaMKII_state[icleft][10],
				CaMKII_rates[icleft][11][0],CaMKII_state[icleft][11]);
		}
		for(icleft=0; icleft<Nclefts_FRU; icleft++){
			printf("%d: %g %d\n",icleft,
				LCCPhosph_rates[icleft][0],LCCPhosph_state[icleft]);
		}
	}
#endif
	return max_rate;
}
