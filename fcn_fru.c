/*       ----------------------------------------------------
 
	 NOTICE OF COPYRIGHT AND OWNERSHIP OF SOFTWARE
	 Copyright 2003, The Johns Hopkins University
	 School of Medicine. All rights reserved.

	 Name of Program: Local Control Model 
	 Version: Documented Version, C
	 Date: September 2003

	 --------------------------------------------------*/

// NB: Indices fixed

#include <math.h>

#include "parameters.h"
//#include "parameters_fcn_fru.h"

// Computes and returns the velocity field for CaSS and CaJSR 
// within the FRU

void fcn_fru(const double time_FRU,
	     const double FRU_states1[Nstates_FRU],
	     const double FRUdep_states[Nstates_FRUdep],
	     const double LType_open[Nclefts_FRU],
	     const int NRyR_open[Nclefts_FRU],
	     double dFRU_states1[Nstates_FRU])
{
	double V,Cai,CaNSR,CaJSR;
	double CaSS_1,CaSS_2,CaSS_3,CaSS_4;
	double VF_over_RT, VFsq_over_RT, exp_VFRT;
	double a1,a2;

	const double JDHconstant=1.0/(2.0*VSS*Faraday*1000.0);

	double Jtr,beta_JSR,JRyRtot;
	double Jxfer_1,Jxfer_2,Jxfer_3,Jxfer_4;
	double Jss2ss_1,Jss2ss_2,Jss2ss_3,Jss2ss_4;
	double JRyR_1,JRyR_2,JRyR_3,JRyR_4;
	double JDHPR_1,JDHPR_2,JDHPR_3,JDHPR_4;
	double beta_SS_1,beta_SS_2,beta_SS_3,beta_SS_4;

	V = FRUdep_states[index_frudep_V];
	Cai = FRUdep_states[index_frudep_Cai];
	CaNSR = FRUdep_states[index_frudep_CaNSR];
	CaJSR = FRU_states1[index_frustates_CaJSR];

	CaSS_1 = FRU_states1[1];
	CaSS_2 = FRU_states1[2];
	CaSS_3 = FRU_states1[3];
	CaSS_4 = FRU_states1[4];

	Jtr = (CaNSR - CaJSR)/tautr;
	// JRyRtot = 0.0
	// for(icleft = 1;icleft<=Nclefts_FRU
	//     JRyR[icleft] = JRyRmax*(double)(NRyR_open[icleft])*(CaJSR-CaSS[icleft])
	//     JRyRtot = JRyRtot + JRyR(icleft]
	//     Jxfer[icleft] = (CaSS[icleft]-Cai)/tauxfer
	// }

	JRyR_1 = JRyRmax*((double)NRyR_open[0])*(CaJSR-CaSS_1);
	JRyR_2 = JRyRmax*((double)NRyR_open[1])*(CaJSR-CaSS_2);
	JRyR_3 = JRyRmax*((double)NRyR_open[2])*(CaJSR-CaSS_3);
	JRyR_4 = JRyRmax*((double)NRyR_open[3])*(CaJSR-CaSS_4);
	JRyRtot = JRyR_1+JRyR_2+JRyR_3+JRyR_4;

	Jxfer_1 = (CaSS_1-Cai)/tauxfer;
	Jxfer_2 = (CaSS_2-Cai)/tauxfer;
	Jxfer_3 = (CaSS_3-Cai)/tauxfer;
	Jxfer_4 = (CaSS_4-Cai)/tauxfer;

	Jss2ss_1 = (CaSS_1 - CaSS_2)/tauss2ss;
	Jss2ss_2 = (CaSS_2 - CaSS_3)/tauss2ss;
	Jss2ss_3 = (CaSS_3 - CaSS_4)/tauss2ss;
	Jss2ss_4 = (CaSS_4 - CaSS_1)/tauss2ss;

	VF_over_RT=V/RT_over_F;
	VFsq_over_RT=(1000.0*Faraday)*VF_over_RT;
	exp_VFRT = exp(2.0*VF_over_RT);

	if (fabs(V)<1.e-6) {
	  //   for(icleft = 1;icleft<=Nclefts_FRU // First order Taylor expansion
	  //     a1 =  CaSS(icleft]-Cao*0.341 
	  //     ICa(icleft] = PCa*2.0*1000.0*Faraday*a1*LType_open[icleft] //uA
	  //     JDHPR(icleft] = -ICa(icleft]*JDHconstant
	  //   }
	  JDHPR_1 = -PCa*2.0*1000.0*Faraday*(CaSS_1-Cao*0.341)*LType_open[0]*JDHconstant;
	  JDHPR_2 = -PCa*2.0*1000.0*Faraday*(CaSS_2-Cao*0.341)*LType_open[1]*JDHconstant;
	  JDHPR_3 = -PCa*2.0*1000.0*Faraday*(CaSS_3-Cao*0.341)*LType_open[2]*JDHconstant;
	  JDHPR_4 = -PCa*2.0*1000.0*Faraday*(CaSS_4-Cao*0.341)*LType_open[3]*JDHconstant;
	} else {
	  //  	  a2 = exp_VFRT-1.0;
	  //   for(icleft = 1;icleft<=Nclefts_FRU
	  //     a1 =  CaSS(icleft]*exp_VFRT-Cao*0.341 
	  //     ICa(icleft] = PCa*4.0*VFsq_over_RT*(a1/a2)*(double)(LType_open[icleft]) //uA
	  //     JDHPR(icleft] = -ICa(icleft]*JDHconstant
	  //   }

	  // The same unrolled
	  a2 = PCa*4.0*VFsq_over_RT/(exp_VFRT-1.0)*JDHconstant;
	  JDHPR_1 = -(CaSS_1*exp_VFRT-Cao*0.341)*a2*LType_open[0];
	  JDHPR_2 = -(CaSS_2*exp_VFRT-Cao*0.341)*a2*LType_open[1];
	  JDHPR_3 = -(CaSS_3*exp_VFRT-Cao*0.341)*a2*LType_open[2];
	  JDHPR_4 = -(CaSS_4*exp_VFRT-Cao*0.341)*a2*LType_open[3];
	}

	// for(icleft = 0;icleft<Nclefts_FRU;icleft++) {
	//     a1 = BSRtot*KBSR/((CaSS[icleft]+KBSR)*(CaSS[icleft]+KBSR))
	//     a2 = BSLtot*KBSL/((CaSS[icleft]+KBSL)*(CaSS[icleft]+KBSL))
	//     beta_SS(icleft] = 1.0/(1.0+a1+a2) 
	// }

	// The same unrolled
	a1 = BSRtot*KBSR/((CaSS_1+KBSR)*(CaSS_1+KBSR));
	a2 = BSLtot*KBSL/((CaSS_1+KBSL)*(CaSS_1+KBSL));
	beta_SS_1 = 1.0/(1.0+a1+a2) ;
	a1 = BSRtot*KBSR/((CaSS_2+KBSR)*(CaSS_2+KBSR));
	a2 = BSLtot*KBSL/((CaSS_2+KBSL)*(CaSS_2+KBSL));
	beta_SS_2 = 1.0/(1.0+a1+a2) ;
	a1 = BSRtot*KBSR/((CaSS_3+KBSR)*(CaSS_3+KBSR));
	a2 = BSLtot*KBSL/((CaSS_3+KBSL)*(CaSS_3+KBSL));
	beta_SS_3 = 1.0/(1.0+a1+a2) ;
	a1 = BSRtot*KBSR/((CaSS_4+KBSR)*(CaSS_4+KBSR));
	a2 = BSLtot*KBSL/((CaSS_4+KBSL)*(CaSS_4+KBSL));
	beta_SS_4 = 1.0/(1.0+a1+a2) ;

	a1 = CSQNtot*KmCSQN/((CaJSR+KmCSQN)*(CaJSR+KmCSQN));
	beta_JSR = 1.0/(1.0+a1);

	dFRU_states1[index_frustates_CaJSR] = beta_JSR*(Jtr - VSS/VJSR*JRyRtot); // dCaJSR

	dFRU_states1[1] = beta_SS_1*(JDHPR_1 + JRyR_1 - Jxfer_1 - Jss2ss_1 + Jss2ss_4); // dCaSS1
	dFRU_states1[2] = beta_SS_2*(JDHPR_2 + JRyR_2 - Jxfer_2 - Jss2ss_2 + Jss2ss_1); // dCaSS2
	dFRU_states1[3] = beta_SS_3*(JDHPR_3 + JRyR_3 - Jxfer_3 - Jss2ss_3 + Jss2ss_2); // dCaSS3
	dFRU_states1[4] = beta_SS_4*(JDHPR_4 + JRyR_4 - Jxfer_4 - Jss2ss_4 + Jss2ss_3); // dCaSS4
}


