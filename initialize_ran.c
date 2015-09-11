/*       ----------------------------------------------------

	 NOTICE OF COPYRIGHT AND OWNERSHIP OF SOFTWARE
	 Copyright 2003, The Johns Hopkins University
	 School of Medicine. All rights reserved.

	 Name of Program: Local Control Model 
	 Version: Documented Version, C
	 Date: November 2003

	 --------------------------------------------------

	 initialize_ran.c - Initialize pseudo random number generators

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "parameters.h"

void initialize_ran(unsigned long mt[NFRU_sim_max][mtN+1],int mti[NFRU_sim_max])
{
  int a,j;

/* Initializes the random number generators by passing
   a seeds for the first call. */
	
  a=17256550;
  for(j = 0; j<NFRU_sim;j++) {
    a=a+23*j;
    sgrnd_local(abs(a+4093*j*j)+1,&mti[j],&mt[j][0]);
  }

}
