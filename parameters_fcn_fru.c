/*       ----------------------------------------------------

	 NOTICE OF COPYRIGHT AND OWNERSHIP OF SOFTWARE
	 Copyright 2003, The Johns Hopkins University
	 School of Medicine. All rights reserved.

	 Name of Program: Local Control Model 
	 Version: Documented Version, C
	 Date: November 2003

	 --------------------------------------------------

	 parameters_fcn_fru.h - constants for FRUs

*/

#include "parameters.h"

/*------STANDARD IONIC CONCENTRATIONS-----------------------------------------*/

const double  Ko=    4.0;   // extracellular K+   concentration (mM)
const double  Nao= 138.0;   // extracellular Na+  concentration (mM)
const double  Cao=   2.0;   // extracellular Ca++ concentration (mM)
const double  Clo= 150.0;   // extracellular Cl-  concentration (mM)
const double  Cli=  20.0;   // intracellular Cl-  concentration (mM)

/*------PHYSICAL constANTS----------------------------------------------------*/

const double  Faraday=  96.5;     // Faraday's constant (C/mmol)
const double  Temp=    310.0;     // absolute temperature (K)
const double  Rgas=      8.314;   // ideal gas constant (J/[mol*K])
const double  RT_over_F= 8.314*310.0/96.5; // (Rgas*Temp/Faraday); //  Rgas*Temp/Faraday (mV)


/*------CELL GEOMETRY PARAMETERS----------------------------------------------*/

/* NOTE: The assumption of 1uF/cm^2 holds for this model
	      therefore Acap in cm^2 is equal to whole cell
	      capacitance in uF.
*/

const double Acap= 1.534e-4; // capacitive membrane area (cm^2)
const double VNSR =(0.53*2.10e-6); // junctional SR volume (uL)
const double Vmyo= 25.84e-6; // myoplasmic volume (uL)
const double VJSR= (((double)Nclefts_FRU)*0.53*0.5*0.2*1.05e-10); // network SR volume (uL)
const double VSS= (0.5*0.2*2.03e-12); // subspace volume (uL)

const double  PCa= (1.5/2.8*0.2*0.9*0.9468e-11); //(cm/s) *uF  

const double  PCl= 2.65e-15; //(cm/s) *uF 
	 
const double  BSLtot=   1.124;  // (mM)
const double  CSQNtot= 13.5; 
const double  BSRtot=   0.047; 
const double  KBSL=     0.0087; 
const double  KmCSQN=   0.63; 
const double  KBSR=     0.00087; 

/*------BUFFERING PARAMETERS--------------------------------------------------*/

	// total troponin low affinity site conc. (mM)
const double  LTRPNtot= 70.0e-3;
	// total troponin high affinity site conc. (mM)
const double  HTRPNtot= 140.0e-3;
	// Ca++ on rate for troponin high affinity sites (1/[mM*ms])
const double  khtrpn_plus= 20.0;
	// Ca++ off rate for troponin high affinity sites (1/ms)
const double  khtrpn_minus= 0.066e-3;
	// Ca++ on rate for troponin low affinity sites (1/[mM*ms])
const double  kltrpn_plus= 40.0;
	// Ca++ off rate for troponin low affinity sites (1/ms)
const double  kltrpn_minus= 40.0e-3;
	// total myoplasmic calmodulin concentration (mM)
const double  CMDNtot= 50.0e-3;
	// total myoplasmic EGTA concentration (mM)
const double  EGTAtot= 0.0;
	// Ca++ half sat. constant for calmodulin (mM)
const double  KmCMDN= 2.38e-3;
	// Ca++ half sat. constant for EGTA (mM)
const double  KmEGTA= 1.5e-4;

	// JSR to subspace through a single RyR (1/ms)
const double  JRyRmax= (0.6*1.5*1.1*3.96);
	// NSR to JSR (1/ms)
const double  tautr= 3.0;
	// subspace to cytosol (1/ms)
const double  tauxfer= 0.005;
	// intersubspace transfer rate (1/ms)
//const double  tauss2ss= (10.0*tauxfer);
const double  tauss2ss= (10.0*0.005);

