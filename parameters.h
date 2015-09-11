/*       ----------------------------------------------------

	 NOTICE OF COPYRIGHT AND OWNERSHIP OF SOFTWARE
	 Copyright 2003, The Johns Hopkins University
	 School of Medicine. All rights reserved.

	 Name of Program: Local Control Model 
	 Version: Documented Version, C
	 Date: November 2003

	 --------------------------------------------------

	 parameters.h - constants, definitions and prototypes

*/

#include <stdio.h>

#include "indices.h"

// definition describing the computer
#define USE_MPI 1
#define USE_PD_INTEGRATOR 1
#define SIMPLEDEBUG 0
#define NOMEMCPY 1
#define MAXPROCS 64
#define MAXFILES 200

// definitions of the model
#define Nclefts_FRU 4
#define Nstates_FRU (1+Nclefts_FRU)
#define Nstates_FRUdep 3 

#define Nstates_LType 12
#define Nstates_RyR 6
#define NRyRs_per_cleft 5
#define Nindepstates_LType 2

#define NRVseqs_per_cleft (2+NRyRs_per_cleft+1)
#define NRVseqs_per_FRU (NRVseqs_per_cleft*Nclefts_FRU)

// The number of FRUs simulated, in low FRU# phase and in high FRU# phase
#define NFRU_sim_low 250
#define NFRU_scale_low 50.0 // ratio of 12500/NFRU_sim_low
#define NFRU_sim_high 250
#define NFRU_scale_high 50.0 // ratio of 12500/NFRU_sim_high

// for array sizes, has to be NFRU_sim_max> NFRU_sim_high
// MAX_LOCAL_NFRU should probably be the same as NFRU_sim_max
#define NFRU_sim_max 250
#define MAX_LOCAL_NFRU 250

// Border between low and high FRU number phases

#define low_phase_V -60.0
#define high_phase_V -50.0

// Random number generators
// MUST BE A NEGATIVE INTEGER
//#define ran_seed -895976	

#define use_seeds_file_flag 0
#define save_seeds_file_flag 0

//	//mt19937 parameters
#define mtN 624

/*	
	Simulation runs from 'time_start' to 'time_end'
	time_step = 	Time_step is the time step at which data is 
			to be saved, the integrator will be called once
			each time_step (in which adaptive or many fixed 
			smaller steps are taken and will force a data 
			point to be aligned in time every time_step.
			(This is different from the previous version, in that
			data is saved every time_step, not every 
			Nmesh*time_step.   Nmesh is not used in this 
			version.  Time_step may be large relative to the 
			adaptive time steps taken within each call to the 
			integrator, see step_max
	step_max =  	maximum time step attempted by the integrator 
			when using an adaptive algorithm, this 
		    	should be .LE. time_step. (Differs from previous 
			version in that the integrator will not 
			automatically attempt a step at time_step, rather 
			it will use min(time_step,max_step
	step_min =  	minimum time step attempted by the integrator
	tolrk =		truncation error tolerance in adaptive integration 
			algorithm before scaling (normalization
*/

#define time_start 0.0
#define time_end 12000.0
//#define time_end 50.0
#define time_step 1.0
#define step_min (1.e-6)
#define step_max 0.1
#define tolrk (1.e-6)

//	ts_sec = 1 if output time in sec, 0 for ms
#define ts_sec 1

//	algebraic or differential method for the computation of V
#define use_algebraic_flag 0
	// Note: algebraic is not usable at the moment, since Cl is not dynamic

//	pulse_duration = total duration of current pulse in ms 
//	pulse_amplitude = amplitude of pulse in uA/uF
//	period = stimulus period in ms
//	shift = offset to set delay of stim. pulses from time = 0.
//		(all pulses shifted by same amount when repetative
//		 stimuli are used;

#define pulse_duration 0
#define pulse_amplitude 0
#define  period 0
#define shift 0 

// Define parameters for voltage clamp.
#define vclamp_flag 0 
#define vclamp_duration 200.0 
#define vclamp_set   0.0 
#define vclamp_shift 300.0
#define vclamp_hold -80.0 
#define vclamp_period 500.0 

// Define parameters for voltage clamp with a prepulse.
#define vppclamp_flag 1 
#define vppclamp_ppset  -40.0 
#define vppclamp_ppduration 100.0 
#define vppclamp_shift 10000.0 
#define vppclamp_set     0.0 
#define vppclamp_duration 300.0 
#define vppclamp_hold  -80.0 

//	I-V relationship studies
#define iv_flag 0
#define iv_clamp_duration 50.0 
#define iv_clamp_period 200.0 
#define iv_shift 10.0 
#define iv_clamp_set   -40.0 
#define iv_clamp_hold  -100.0 
#define iv_clamp_step  10.0 

// Define conditions for heart failure.
#define chf_flag 0 
#define chfsc_IK1 0.68 
#define chfsc_Jup 0.38 
#define chfsc_INaCa 1.75 
#define chfsc_IKv43 0.143   // When Kv43Frac = 0.77 (in fcn.f 
                            // then chfsc_IKv43=0.143 corresponds to
			    // scaling total Ito1 by 0.34

//	Ap clamp by Antonis A. Armoundas
//     Define action potential clamp parameters
#define apclamp_flag 0
#define rep_interval 0.32
#define APsamples 10
#define ap_clamp_file "ap_clamp.txt"

// Reset RyRs to closed state at the beginning of simulation
// 0 = no,  1 = yes
#define reset_RyRs_flag 0

/*

	Output filenames etc

*/

//       output_dir defines the path to the output data file

#define output_dir "./"
#define output_states_file "states"
#define output_currents_file "currents"
#define output_otherstates_file "otherstates"
#define output_fruprops_file "fruprops"
#define output_frupropsavg_file "fruprops_avg"
#define output_frupropsavgtot_file "fruprops_avg_tot"
#define output_seeds_file "r_randomseeds"

#define filenumber0 0

#define ic_file_flag 1
//#define ic_file_flag=0

#define ic_dir "ic/stable"
#define ic_states_file "ic_states_stable.txt"
#define ic_FRU_file "ic_FRU_stable.txt"
#define ic_LCh_file "ic_LCh_stable.txt"
#define ic_RyR_file "ic_RyR_stable.txt"
#define ic_Ito2_file "ic_Ito2_stable.txt"
#define ic_seeds_file "ic_randomseeds_stable.txt"

//	write_fru_props_flag controls, whether output is written on individual
//	RyR:s or not, if speed is needed set to FALSE
//	if0 no information on open RyR proportion is written to file
//
#define write_fru_props_flag 0
#define fru_start 0
#define fru_end 9

//	Ndump x time_step   time interval for dump restart files
//
#define Ndump 10000

// Number of columns
//	N   number of state variables 
//	Ncur   number of variables in "current" array
#define N 37
#define Ncur 24
#define Nother 12

/*

    Prototypes here

*/
double MersenneTwisterOne_local(int *mti,unsigned long[mtN+1]);
void MersenneTwister_local(int *,unsigned long[mtN+1],const int,double *);
void MersenneTwister_fast(int *,unsigned long [mtN+1],const int ,double *);
void initialize_state(double[N],double[NFRU_sim_max][Nstates_FRU],int[NFRU_sim_max][Nclefts_FRU][Nindepstates_LType],int[NFRU_sim_max][Nclefts_FRU][NRyRs_per_cleft],int[NFRU_sim_max][Nclefts_FRU],unsigned long[NFRU_sim_max][mtN+1],int[NFRU_sim_max]);
double rk54pd(double,double,double,double[N],double[Ncur]);
double pd5_predict_step_size(double,double,double,double[N],double[Ncur],double,double,double,double);
void calc_fru_avg(const double[N],const double[NFRU_sim_max][Nstates_FRU],const int[NFRU_sim_max][Nclefts_FRU][Nindepstates_LType],const int[NFRU_sim_max][Nclefts_FRU][NRyRs_per_cleft],const int[NFRU_sim_max][Nclefts_FRU],double[Nother]);
FILE *open_currents(const int);
FILE *open_otherstates(const int);
FILE *open_states(const int);
void close_currents(void);
void close_states_otherstates(void);
void close_fru_props(void);
void write_states(FILE *,const double,const double[N]);
void write_otherstates(FILE *,const double,const double[Nother]);
void write_currents(FILE *,const double,const double[Ncur]);
void open_fru_props_avg(const int filenumber,FILE *[]);
void open_fru_props(const int,FILE *[]);
FILE *open_fru_props_avgtot(const int);
void write_fru_props(FILE *[],FILE *[],FILE *,const double,const double,const double,const double,double [NFRU_sim_max][Nstates_FRU],int [NFRU_sim_max][Nclefts_FRU][Nindepstates_LType],int [NFRU_sim_max][Nclefts_FRU][NRyRs_per_cleft],int [NFRU_sim_max][Nclefts_FRU]);
double fru_rates_local(int[Nclefts_FRU][Nindepstates_LType],int[Nclefts_FRU][NRyRs_per_cleft],int[Nclefts_FRU],const double[Nstates_FRUdep],const double[Nstates_FRU],double[Nclefts_FRU][4],int[Nclefts_FRU][4],
		       int[Nclefts_FRU],double[Nclefts_FRU],double[Nclefts_FRU][NRyRs_per_cleft][4],int[Nclefts_FRU][NRyRs_per_cleft][4],int[Nclefts_FRU][NRyRs_per_cleft],double[Nclefts_FRU],int *,unsigned long[mtN+1]);
void fcn_fru(const double,const double[Nstates_FRU],const double[Nstates_FRUdep],const double[Nclefts_FRU],const int[Nclefts_FRU],double[Nstates_FRU]);
void fcn(const double,double[N],double[N],double[Ncur],const int,const double,const double,const double,const double);
void simfru_local(const double,const double,double[Nstates_FRUdep],double[Nstates_FRUdep],int[Nclefts_FRU][Nindepstates_LType],int[Nclefts_FRU][NRyRs_per_cleft],int[Nclefts_FRU],double[Nstates_FRU],int *,unsigned long [mtN+1]);
void dynamicFRU(const double,double[N],double[NFRU_sim_max][Nstates_FRU],int[NFRU_sim_max][Nclefts_FRU][Nindepstates_LType],int[NFRU_sim_max][Nclefts_FRU][NRyRs_per_cleft],int[NFRU_sim_max][Nclefts_FRU],int[NFRU_sim_max],unsigned long[NFRU_sim_max][mtN+1],const double,const double);
void lastcall(double,double[N],double[Ncur]);
void initialize_ran(unsigned long[NFRU_sim_max][mtN+1],int[NFRU_sim_max]);
void restart_data(const int,double[N],double[NFRU_sim_max][Nstates_FRU],int[NFRU_sim_max][Nclefts_FRU][Nindepstates_LType],int[NFRU_sim_max][Nclefts_FRU][NRyRs_per_cleft],int[NFRU_sim_max][Nclefts_FRU],unsigned long[NFRU_sim_max][mtN+1],int mti[NFRU_sim_max]);
void sgrnd_local(unsigned long,int *,unsigned long[mtN+1]);
double read_next_double(FILE *);
int read_next_int(FILE *);
void initialize_mpi(int *,char **,int *,int *);
void end_mpi(void);
void parallel(void);
void distrib_simFRU(double,double,double[Nstates_FRUdep],double[Nstates_FRUdep],double *,double *,double *,double *);
void parallel_get_FRUs(double[NFRU_sim_max][Nstates_FRU],int[NFRU_sim_max][Nclefts_FRU][Nindepstates_LType],int[NFRU_sim_max][Nclefts_FRU][NRyRs_per_cleft],int[NFRU_sim_max][Nclefts_FRU],unsigned long[NFRU_sim_max][mtN+1],int[NFRU_sim_max]);
void send_save_state(void);
void parallel_save_state(void);
void parallel_resume_state(void);
void parallel_compute_simfru(void);
void send_resume_state(void);
void send_calc_fru_flux(double[Nstates_FRUdep],double *, double *, double *, double *);
void send_calc_fru_avg(double[N],double[Nother]);
void initialize_mpi_state(double[NFRU_sim_max][Nstates_FRU],int[NFRU_sim_max][Nclefts_FRU][Nindepstates_LType],int[NFRU_sim_max][Nclefts_FRU][NRyRs_per_cleft],int[NFRU_sim_max][Nclefts_FRU],int[NFRU_sim_max],unsigned long[NFRU_sim_max][mtN+1]);
void calc_fru_flux_local(int,int[NFRU_sim_max][Nclefts_FRU][Nindepstates_LType],int[NFRU_sim_max][Nclefts_FRU],double[NFRU_sim_max][Nstates_FRU],double[4]);
void parallel_calc_fru_avg(void);
void parallel_calc_fru_flux(void);
void calc_fru_avg_local(double[11],int,double[MAX_LOCAL_NFRU][Nstates_FRU],int[MAX_LOCAL_NFRU][Nclefts_FRU][Nindepstates_LType],int[MAX_LOCAL_NFRU][Nclefts_FRU][NRyRs_per_cleft],int[MAX_LOCAL_NFRU][Nclefts_FRU]);

// external variables

extern int my_rank;
extern int procs;
extern int NFRU_local;
extern double NFRU_scale;
extern int NFRU_sim;
extern int current_phase;
extern double time_of_last_change;
extern double glob_time_to_stim;
extern double FRU_states_local[MAX_LOCAL_NFRU][Nstates_FRU];
extern int LType_state_local[MAX_LOCAL_NFRU][Nclefts_FRU][Nindepstates_LType];
extern int RyR_state_local[MAX_LOCAL_NFRU][Nclefts_FRU][NRyRs_per_cleft];
extern int Ito2_state_local[MAX_LOCAL_NFRU][Nclefts_FRU];
extern int iv_n;
extern unsigned long mt_local[MAX_LOCAL_NFRU][mtN+1];
extern int mti_local[MAX_LOCAL_NFRU];
extern double FRU_states_local_hold[MAX_LOCAL_NFRU][Nstates_FRU];
extern int LType_state_local_hold[MAX_LOCAL_NFRU][Nclefts_FRU][Nindepstates_LType];
extern int RyR_state_local_hold[MAX_LOCAL_NFRU][Nclefts_FRU][NRyRs_per_cleft];
extern int Ito2_state_local_hold[MAX_LOCAL_NFRU][Nclefts_FRU];
extern unsigned long mt_local_hold[MAX_LOCAL_NFRU][mtN+1];
extern int mti_local_hold[MAX_LOCAL_NFRU];
extern double errweight[N];
extern double NFRU_sim_factor;


//    extern double ap(APsamples;
//	extern double v_at_t, v_at_tstep

// Message definitions

#define MSG_INIT 0
#define MSG_EXIT 1
#define MSG_COMP_SIMFRU 2
#define MSG_RESUME_STATE 3
#define MSG_SAVE_STATE 4
#define MSG_INIT_FRU 5
#define MSG_RETURN_FRU 6
#define MSG_CALC_FRU_AVG 7
#define MSG_CALC_FRU_FLUX 8

// two macro definitions, in case they are needed
#ifndef min
#define min(a,b) ((a)<(b)?(a):(b))
#endif

#ifndef max
#define max(a,b) ((a)<(b)?(b):(a))
#endif

// Open and Closed states in various MC models

#define O1_LType 6
#define O2_LType 12
#define Oy_LType 2
#define Cy_LType 1
#define O1_RyR 3
#define O2_RyR 4
#define O3_RyR 7
#define O_Ito2 2
#define C_Ito2 1

// External physical constants

extern const double  Ko,Nao,Cao,Clo,Cli;
extern const double  Faraday,Temp,Rgas,RT_over_F;
extern const double  Acap,VNSR,Vmyo,VJSR,VSS;
extern const double  PCa,PCl,BSLtot,CSQNtot,BSRtot,KBSL,KmCSQN,KBSR;
extern const double  LTRPNtot,HTRPNtot,khtrpn_plus,khtrpn_minus,kltrpn_plus,kltrpn_minus;
extern const double  CMDNtot,EGTAtot,KmCMDN,KmEGTA;
extern const double  JRyRmax,tautr,tauxfer,tauss2ss;
