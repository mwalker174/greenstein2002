/*       ----------------------------------------------------

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
		       double Ito2_exitrate[Nclefts_FRU],
		       int *mti_loc,
		       unsigned long mt_loc[mtN])
{
	// real	rnum
	double dnum;
	double V,CaSS[Nclefts_FRU];	

	const double fL=0.85; // transition	rate into open state (1/ms)
	const double gL=2.0; //transition rate	out	of open	state (1/ms)
	const double fLprime=0.005;	// transition rate into	Ca mode	open state (1/ms)
	const double gLprime=7.0; // transition	rate out of	Ca mode	open state (1/ms)
	const double bL=1.9356;		// mode	transition parameter
	const double bL2=bL*bL;
	const double bL3=bL*bL*bL;
	const double bL4=bL*bL*bL*bL;
	const double kappa=1;
	const double kappa2 = kappa*kappa;
	const double kappa3 = kappa*kappa*kappa;
	const double kappa4 = kappa*kappa*kappa*kappa;
	const double aL=2.0; //	mode transition	parameter
	const double aL2=aL*aL;
	const double aL3=aL*aL*aL;
	const double aL4=aL*aL*aL*aL;
	const double omega=0.83*2.0*1.3*0.01;  // mode transition parameter	(1/ms)

	const double alphacf=4.0*1.2*0.416;
	const double betacf=4.0*0.45*0.049;
	const double gammacf=0.83*1.9*1.3*0.31*7.5*0.09233;	// (ms-1 mM-1)

	double alpha, beta, alpha_prime, beta_prime, gamma_rate;
	double C0_to_C1, C1_to_C2, C2_to_C3, C3_to_C4;
	double C1_to_C0, C2_to_C1, C3_to_C2, C4_to_C3;
	double CCa0_to_CCa1, CCa1_to_CCa2;
	double CCa2_to_CCa3, CCa3_to_CCa4;
	double CCa1_to_CCa0, CCa2_to_CCa1;
	double CCa3_to_CCa2, CCa4_to_CCa3;
	double C0_to_CCa0, C1_to_CCa1, C2_to_CCa2;
	double C3_to_CCa3, C4_to_CCa4;

	const double CCa0_to_C0	= kappa4*omega;		// = omega
	const double CCa1_to_C1	= kappa3*omega/bL;	// = omega/bL
	const double CCa2_to_C2	= kappa2*omega/bL2;	// = omega/bL^2
	const double CCa3_to_C3	= kappa*omega/bL3;	// = omega/bL^3
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

	double a1,a2,a3;
	int	i, icleft;

	V =	FRUdep_states[index_frudep_V];
	for(icleft = 0;icleft<Nclefts_FRU;icleft++) {
		CaSS[icleft] = FRU_states[icleft+1];
	}

	alpha =	alphacf	* exp(0.012*(V-35.0));
	beta = betacf *	exp(-0.05*(V-35.0));
	alpha_prime	= aL*alpha;
	beta_prime = beta/(bL*kappa);

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
	    LType_rates[icleft][0] = gL;
	    LType_rates[icleft][1] = gL;
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
	double nu = 1; //1.5;
	yCa_inf	= nu*( yCa_frac/(1.0+exp((V + 12.5)/5.0)) + (1.0-yCa_frac));
	//ORIGINALLY, MULTIPLIED TAU BY NU 
	tau_yCa	=  60.0 + 340.0/(1.0	+ exp((V+30.0)/12.0));

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
	}
#endif
	return max_rate;
}
