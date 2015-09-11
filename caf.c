/*       ----------------------------------------------------

	 NOTICE OF COPYRIGHT AND OWNERSHIP OF SOFTWARE
	 Copyright 2003, The Johns Hopkins University
	 School of Medicine. All rights reserved.

	 Name of Program: Local Control Model 
	 Version: Documented Version, C
	 Date: November 2003

	 --------------------------------------------------

	 caf.c - main routine

*/

#include <stdio.h>
#include <time.h>
#include <math.h>

#include "parameters.h"

#if USE_MPI
#include <mpi.h>
#endif

double time_of_last_change,glob_time_to_stim;
int current_phase;
double surand_seed[NFRU_sim_max];
double errweight[N];
int iv_n;

// important variables 
double NFRU_scale;
double NFRU_sim_factor;
int NFRU_sim;
int my_rank;
int procs;

/* 
	Main routine
*/
int main(int argc,char **argv)
{
printf("started main\n");
#if USE_MPI
  printf("about to initialize MPI\n");
  initialize_mpi(&argc,argv,&my_rank,&procs);
  printf("initialized MPI\n");
#else
  procs=1;
  my_rank=0;
#endif

  if (my_rank==0) {
	// Master of the puppets

    double state[N];
    double current[Ncur];
    double otherstates[Nother];
    //double Act_coeff[8];

    double time_now, tstep, oldstepsize,iv_time;
    int  filenumber;
    int  n_restart;

	// Files for output
    FILE *otherstates_file;
    FILE *states_file;
    FILE *currents_file;
    FILE *fruprop_file[MAXFILES],*frupropavg_file[MAXFILES],*frupropavgtot_file;

    // Global variables used for initializing and saving values, 
    // local ones (ending with _local) are are used for actual computations
    double FRU_states[NFRU_sim_max][Nstates_FRU];
    int LType_state[NFRU_sim_max][Nclefts_FRU][Nindepstates_LType];
    int RyR_state[NFRU_sim_max][Nclefts_FRU][NRyRs_per_cleft];
    int CaMKII_state[NFRU_sim_max][Nclefts_FRU][Nmon_per_holo];
    int LCCPhosph_state[NFRU_sim_max][Nclefts_FRU];
    int RyRPhosph_state[NFRU_sim_max][Nclefts_FRU][NRyRs_per_cleft];
    int Ito2_state[NFRU_sim_max][Nclefts_FRU];
    unsigned long mt[NFRU_sim_max][mtN+1]; /* the array for the state vector  */
    int mti[NFRU_sim_max]; /* mti==mtN+1 means mt[mtN] is not initialized */

    // timing
    clock_t initial_clock,end_clock;
    double dt1;
    
    // Timing initialization
    initial_clock=clock();
    
    // Setting up FRU numbers
    NFRU_sim=NFRU_sim_high;
    current_phase=1;
    NFRU_scale=NFRU_scale_high;
	NFRU_sim_factor=NFRU_sim_high/NFRU_sim_low;
    
	// Checks to make sure certain parameters are ok
    if (fabs(NFRU_scale_low*NFRU_sim_low-NFRU_scale_high*NFRU_sim_high)>1.e-8) {
      puts("NFRU_scale_low*NFRU_sim_low differs from NFRU_scale_high*NFRU_sim_high!");
    }

    if (procs>MAXPROCS) {
      printf("procs (%d) > MAXPROCS (%d)\n",procs,MAXPROCS);
    }
    
    if (NFRU_sim_factor<1.0) {
      puts("NFRU_sim_factor less than one!");
    }

    if (vclamp_flag&&iv_flag) {
      puts("Warning: vclamp on at the same time as ivclamp on!");
    }

    if (((fru_end-fru_start)>MAXFILES)&&(write_fru_props_flag==1)) {
      puts("Warning, fru_end-fru_start>MAXFILES!");
    }
    
    // step size for integrator
    oldstepsize=step_max;
    
    // initial count of resets in filenumber
    filenumber = filenumber0;
    
    // Initialize random number generator seeds for each sequence
    initialize_ran(mt,mti);
    
    // Set initial conditions on state variables
    initialize_state(state,FRU_states,LType_state,RyR_state,CaMKII_state,LCCPhosph_state,RyRPhosph_state,Ito2_state,mt,mti);

	//initActCoeff(Act_coeff);
	//fprintf(stdout, "%g", Act_coeff[3]);
	// Set up MPI system and distribute data to slave processes
    initialize_mpi_state(FRU_states,LType_state,RyR_state,CaMKII_state,LCCPhosph_state,RyRPhosph_state,Ito2_state,mti,mt);
    
    // open output files
    states_file=open_states(filenumber);
    currents_file=open_currents(filenumber);
    otherstates_file=open_otherstates(filenumber);
    
    if (write_fru_props_flag) { // Open output files
      open_fru_props(filenumber,fruprop_file);
      open_fru_props_avg(filenumber,frupropavg_file);
      frupropavgtot_file=open_fru_props_avgtot(filenumber);
    }
    
    n_restart = 0;
    iv_time = 0.0;
    iv_n = 0;
    
    if ((time_step>100.0)&&(NFRU_sim_factor>1)) {
      puts("Warning: dynamic FRU number does not work with time_step>100ms.");
    }

    // Print simulation time to screen each time_step
    printf("%s %s %s %s\n","Time","Phase","NFRU","V");
    
	// Main loop
    for(time_now=time_start;time_now<=time_end-.0001*time_step;time_now+=time_step) {		// Start integration loop
      // Print simulation time to screen each time_step
      printf("%g %i %i %g\n",time_now,current_phase,NFRU_sim,state[index_V]);
      
      n_restart = n_restart + 1;
      tstep = time_now + time_step;
      iv_time = iv_time + time_step;
      
      // Integrate the whole system over 
      // time interval from 'time' to 'tstep'
#if USE_PD_INTEGRATOR
      oldstepsize=rk54pd(time_now,tstep,oldstepsize,state,current);
#else
      oldstepsize=rk4am(time_now,tstep,oldstepsize,state,current);
#endif 
      //	    if (apclamp_flag) state[index_V) = v_at_t;
      
      // Calculate additional output quantities that are algebraically 
      // related to the state variables and stored in 'otherstate'
      printf("about to send calc fru avg\n");
	send_calc_fru_avg(state,otherstates);
      printf("sent calc fru avg\n");
      if (time_now>time_start) { // Write data in 'current' to the output files
		
		// Close output files and open new output files for data in 'current'
		// if this is a restart (see below)
		if (n_restart==1) {
			fclose(currents_file);
			currents_file=open_currents(filenumber);
		} else {
			write_currents(currents_file,time_now,current);
		}
      }
      
      // Write data in 'otherstates' to the output files
      write_states(states_file,tstep,state);
      write_otherstates(otherstates_file,tstep,otherstates);
	printf("wrote other states\n");
      
      if (write_fru_props_flag) {
		// Write data which define the local events with the release unit 
		// to the output files only if the flag is set
		parallel_get_FRUs(FRU_states,LType_state,RyR_state,CaMKII_state,LCCPhosph_state,RyRPhosph_state,Ito2_state,mt,mti);
		write_fru_props(&fruprop_file[0],&frupropavg_file[0],frupropavgtot_file,tstep,
				state[index_V],state[index_Cai],state[index_CaNSR],
				FRU_states,LType_state,LCCPhosph_state,RyR_state,RyRPhosph_state,Ito2_state);
      }
	    
      if (n_restart >= Ndump) {
		n_restart=0;
		      
		// First collect all states from slave processes
		parallel_get_FRUs(FRU_states,LType_state,RyR_state,CaMKII_state,LCCPhosph_state,RyRPhosph_state,Ito2_state,mt,mti);
		printf("Got FRUs\n");
		// Then save data
		restart_data(filenumber,state,FRU_states,LType_state,RyR_state,CaMKII_state,LCCPhosph_state,RyRPhosph_state,Ito2_state,mt,mti);
		printf("Restarted data\n");
		filenumber=filenumber+1;
		
		// Close output files
		fclose(otherstates_file);
		fclose(states_file);
		
		if (write_fru_props_flag) {
			int j;
		  
			for(j = fru_start;j<=fru_end;j++) {
				fclose(fruprop_file[j-fru_start]);
				fclose(frupropavg_file[j-fru_start]);
			}
			fclose(frupropavgtot_file);
		}

		if (time_now < time_end-1.0001*time_step) {
			states_file=open_states(filenumber);
			otherstates_file=open_otherstates(filenumber);
			if (write_fru_props_flag) {
				open_fru_props(filenumber,fruprop_file);
				open_fru_props_avg(filenumber,frupropavg_file);
				frupropavgtot_file=open_fru_props_avgtot(filenumber);
			}
		}
		printf("Opened new files\n");
      }	// Close output files, write data to restart files, and open new 
      // output files for all data (except currents which are done seperately
      // above because they are calculated after the next time step)
	    
      if (iv_flag&&(iv_time>=iv_clamp_period)) {
		iv_time=0.0;
		iv_n=iv_n+1;
		      
		oldstepsize=step_max;
		// Initialize random number generator seeds for each sequence
		initialize_ran(mt,mti);	
		// Set initial conditions on state variables
		initialize_state(state,FRU_states,LType_state,RyR_state,CaMKII_state,LCCPhosph_state,RyRPhosph_state,Ito2_state,mt,mti);

		// Redistribute initialized states to slave processes
		initialize_mpi_state(FRU_states,LType_state,RyR_state,CaMKII_state,LCCPhosph_state,RyRPhosph_state,Ito2_state,mti,mt);	      
      }
      
      // Dynamic change of FRU number, this is probably not practical in MPI implementation
      if (NFRU_sim_high!=NFRU_sim_low) {
		dynamicFRU(tstep,state,FRU_states,LType_state,RyR_state,CaMKII_state,LCCPhosph_state,RyRPhosph_state,Ito2_state,mti,mt,
			otherstates[index_PRyR_Open],otherstates[index_PLType_Open]);
      }
	printf("Ran dynamicFRU\n");
    }
    
    // Calculate and write currents for the final time step
    lastcall(time_end,state,current);
    write_currents(currents_file,time_end,current);
    
    // Close output files
    fclose(currents_file);
    fclose(otherstates_file);
    fclose(states_file);
	
    if (write_fru_props_flag) {
	  int j;
  
	  for(j = fru_start;j<=fru_end;j++) {
	    fclose(fruprop_file[j-fru_start]);
	    fclose(frupropavg_file[j-fru_start]);
	  }
	  fclose(frupropavgtot_file);
    }
    
    end_mpi();
    
    // timings
    end_clock=clock();
    dt1=(double)(end_clock-initial_clock)/(double)CLOCKS_PER_SEC;
    printf("Total time: %fs\n",dt1);
  } else {
    // Puppet ie slave process

    parallel();
  }
  
  return 0;
}

