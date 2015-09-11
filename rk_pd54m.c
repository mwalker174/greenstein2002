/*
* ---------------------------------------------------- 
*
* Notice OF COPYRIGHT AND OWNERSHIP OF SOFTWARE Copyright 2003, The Johns
* Hopkins University School of Medicine. All rights reserved. 
*
* Name of Program: Local Control Model Version: Documented Version, C Date:
* September 2003 
*
* -------------------------------------------------- 
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

#include "parameters.h"


#if USE_MPI
#include <mpi.h>
#endif

//#include "parameters_fcn_fru.h"

#define min_tolvalue 1.e-10

double rk54pd(double t,double tstep,double oldstepsize,
	      double state[N],double current[Ncur])
{
  //Prince - Dormand version of Runge - Kutta, rk5(4) 6m, 5th / 4th Order ADAPTIVE Step Algorithm
  // Engeln - M \ "ullgess, Uhlig: Numerical Algorithms with C. Page 440. Springer, 1996.

  // Passed Variables:
  //t starting time of integration
  // tstep End time of integration
  // state Contains initial state at time t on entry.
  // Contains final state at time t + tstep on exit
  
  // This routine controls the solution of the system of differential
  // equations from time = t to time = tstep by monitoring the truncation
  // error on each incremental step, and adjusting the step_size based
  // on the error after each attempted step.
  
  //
  //
  double          F[N], F1[N];
  double y_1[N], ym[N]; // y1 is a Bessel function in math library, hence y1->y_1
  double k1[N], k2[N], k3[N], k4[N], k5[N], k6[N];
  double          start_time, st_time, step_size, tr_error, errtmp, ostepsize, extstep = 0;
  double          dummycurrent[2];
  
  int             i;
  //int             errmax;
  //int           j;
  int             notdone, success, keepc, DepFlag;
  const int       dummyfalse = 0;
  const double    huge_number = 1e99;

  //int             iFRU, icleft;
  double          FRUdep_states0[Nstates_FRUdep];
  double          FRUdep_statesf[Nstates_FRUdep];
  double          time_to_stim;
  
  int             notaccepted, stepsno, forcedaccept;
  
  double          Jxfer, Jtr, ICa, Ito2;
  double          end_time;
  
  step_size = min(oldstepsize, tstep - t);
  //start with previous stepsize
  step_size = min(step_max, step_size);
  ostepsize = step_size;
  start_time = t;
  notdone = 1;
  success = 1;
  keepc = 1;
  
  DepFlag = 1;
  stepsno = 0;
  notaccepted = 0;
  forcedaccept = 0;
  
  while (notdone) {

    stepsno = stepsno + 1;
    
    FRUdep_states0[index_frudep_V] = state[index_V];
    FRUdep_states0[index_frudep_Cai] = state[index_Cai];
    FRUdep_states0[index_frudep_CaNSR] = state[index_CaNSR];
    
    if (success) {
      
      //Calculate fluxes that cross functional unit boundaries
      // (Jxfer, Jtr, ICa, Ito2) so that they may be passed to fcn.f
      send_calc_fru_flux(FRUdep_states0,&Jxfer, &Jtr, &ICa, &Ito2);
      
      //Calculate global velocity field
      fcn(start_time, state, F1, current, keepc, Jxfer, Jtr, ICa, Ito2);
      
      keepc = 0;
      
      //Store all FRU states and the random number seeds before proceeding here.If
      //this time step is rejected all values will be reset to these stored values.
      send_save_state();
      
    } //(success)
    
    for (i = 0; i < N; i++) {
      //First intermediate RK step
      k1[i] = step_size * F1[i];
      y_1[i] = state[i] + k1[i] / 5.0;
    }
    
    //Set estimated values for V,Cai, CaNSR at the end
    // of the intermediate RK54PD time step so that these
    // quantities can be interpolated during the upcoming
    // Monte Carlo simulation (call to distrib_simfru)
    FRUdep_statesf[index_frudep_V] = y_1[index_V];
    FRUdep_statesf[index_frudep_Cai] = y_1[index_Cai];
    FRUdep_statesf[index_frudep_CaNSR] = y_1[index_CaNSR];

    end_time = start_time + step_size / 5.0;
    distrib_simFRU(start_time,end_time,FRUdep_states0,FRUdep_statesf,
		   &Jxfer, &Jtr, &ICa, &Ito2);
    
    FRUdep_states0[index_frudep_V] = state[index_V];
    FRUdep_states0[index_frudep_Cai] = state[index_Cai];
    FRUdep_states0[index_frudep_CaNSR] = state[index_CaNSR];
    
    fcn(start_time + step_size / 5.0, y_1, F, dummycurrent, dummyfalse, 
	Jxfer, Jtr, ICa, Ito2);
    
    for (i = 0; i < N; i++) {
      //Second intermediate RK4M step
      k2[i] = step_size * F[i];
      y_1[i] = state[i] + (3.0 / 40.0 * k1[i] + 9.0 / 40.0 * k2[i]);
    }
    
    FRUdep_statesf[index_frudep_V] = y_1[index_V];
    FRUdep_statesf[index_frudep_Cai] = y_1[index_Cai];
    FRUdep_statesf[index_frudep_CaNSR] = y_1[index_CaNSR];
    
    st_time = start_time + step_size / 5.0;
    end_time = start_time + step_size * 3.0 / 10.0;
    distrib_simFRU(st_time,end_time,FRUdep_states0,FRUdep_statesf,
		   &Jxfer, &Jtr, &ICa, &Ito2);

    FRUdep_states0[index_frudep_V] = state[index_V];
    FRUdep_states0[index_frudep_Cai] = state[index_Cai];
    FRUdep_states0[index_frudep_CaNSR] = state[index_CaNSR];

    fcn(start_time + step_size * 3.0 / 10.0, y_1, F, dummycurrent, dummyfalse, Jxfer, Jtr, ICa, Ito2);
    
    for (i = 0; i < N; i++) {
      //Third intermediate RK4M step
      k3[i] = step_size * F[i];
      y_1[i] = state[i] + (3.0 / 10.0 * k1[i] - 9.0 / 10.0 * k2[i] + 6.0 / 5.0 * k3[i]);
    }
    
    FRUdep_statesf[index_frudep_V] = y_1[index_V];
    FRUdep_statesf[index_frudep_Cai] = y_1[index_Cai];
    FRUdep_statesf[index_frudep_CaNSR] = y_1[index_CaNSR];
    
    st_time = start_time + step_size * 3.0 / 10.0;
    end_time = start_time + step_size * 3.0 / 5.0;
    distrib_simFRU(st_time,end_time,FRUdep_states0,FRUdep_statesf,
		   &Jxfer, &Jtr, &ICa, &Ito2);
    
    FRUdep_states0[index_frudep_V] = state[index_V];
    FRUdep_states0[index_frudep_Cai] = state[index_Cai];
    FRUdep_states0[index_frudep_CaNSR] = state[index_CaNSR];
    
    fcn(start_time + step_size * 3.0 / 5.0, y_1, F, dummycurrent, dummyfalse, Jxfer, Jtr, ICa, Ito2);
    
    for (i = 0; i < N; i++) {
      //Forth intermediate RK4M step
      k4[i] = F[i] * step_size;
      y_1[i] = state[i] + 226.0 / 729.0 * k1[i] - 25.0 / 27.0 * k2[i] + 880.0 / 729.0 * k3[i] + 55.0 / 729.0 * k4[i];
    }
    
    FRUdep_statesf[index_frudep_V] = y_1[index_V];
    FRUdep_statesf[index_frudep_Cai] = y_1[index_Cai];
    FRUdep_statesf[index_frudep_CaNSR] = y_1[index_CaNSR];
    
    st_time = start_time + step_size * 3.0 / 5.0;
    end_time = start_time + step_size * 2.0 / 3.0;
    distrib_simFRU(st_time,end_time,FRUdep_states0,FRUdep_statesf,
		   &Jxfer, &Jtr, &ICa, &Ito2);
    
    FRUdep_states0[index_frudep_V] = state[index_V];
    FRUdep_states0[index_frudep_Cai] = state[index_Cai];
    FRUdep_states0[index_frudep_CaNSR] = state[index_CaNSR];
    
    fcn(start_time + step_size * 2.0 / 3.0, y_1, F, dummycurrent, dummyfalse, Jxfer, Jtr, ICa, Ito2);
    
    for (i = 0; i < N; i++) {
      //Fifth intermediate RK4M step
      k5[i] = step_size * F[i];
      y_1[i] = state[i] - 181.0 / 270.0 * k1[i] + 5.0 / 2.0 * k2[i] - 266.0 / 297.0 * k3[i] - 91.0 / 27.0 * k4[i] + 189.0 / 55.0 * k5[i];
    }
    
    FRUdep_statesf[index_frudep_V] = y_1[index_V];
    FRUdep_statesf[index_frudep_Cai] = y_1[index_Cai];
    FRUdep_statesf[index_frudep_CaNSR] = y_1[index_CaNSR];
    
    st_time = start_time + step_size * 2.0 / 3.0;
    end_time = start_time + step_size;
    distrib_simFRU(st_time,end_time,FRUdep_states0,FRUdep_statesf,
		   &Jxfer, &Jtr, &ICa, &Ito2);
    
    fcn(start_time + step_size, y_1, F, dummycurrent, dummyfalse, Jxfer, Jtr, ICa, Ito2);

    for (i = 0; i < N; i++) {
      //Fifth intermediate RK4M step
      k6[i] = step_size * F[i];
    }
    DepFlag = 0;
    //for next loop
    
    // Calculation of truncation error.The value of tr_error is
    // the maximum normalized error(maximum over all states).
    
    tr_error = 0.0;
    //errmax -= 1;
    for (i = 0; i < N; i++) {
      y_1[i] = state[i] + 31.0 / 540.0 * k1[i] + 190.0 / 297.0 * k3[i] - 145.0 / 108.0 * k4[i] + 351.0 / 220.0 * k5[i] + 1.0 / 20.0 * k6[i];
      ym[i] = state[i] + 19.0 / 216.0 * k1[i] + 1000.0 / 2079.0 * k3[i] - 125.0 / 216.0 * k4[i] + 81.0 / 88.0 * k5[i] + 5.0 / 56.0 * k6[i];
      errtmp = fabs((ym[i] - y_1[i]) * 0.2 * errweight[i]);
      if (errtmp < huge_number) {
	//if (errtmp > tr_error)
	//errmax = i;
	tr_error = max(errtmp, tr_error);
      } else {
	tr_error = huge_number;
      }
    }
    //printf("max error index=%d, error=%f\n", errmax, tr_error);
    
    if (tr_error < tolrk) {
      for (i = 0; i < N; i++) {
	if (fabs(y_1[i]) < min_tolvalue) {
	  y_1[i] = 0.0;
	}
	state[i] = y_1[i];
      }
      start_time = start_time + step_size;
      //puts("accept");
      
      success = 1;
      
      if (start_time >= tstep) {
	notdone = 0;
	//ready, lets go
      } else {
	
	step_size = min(0.85 * step_size * pow(tolrk / tr_error, 0.2), 4.0 * step_size);

	//extstep = pd5_predict_step_size(0.9 * step_size * pow(tolrk / tr_error, 0.2),
	//				start_time, start_time + step_max, state, current,
	//				Jxfer, Jtr, ICa, Ito2);
	
	//printf("orig: pred %g < %g\n", extstep, step_size);
	
	//step_size = min(step_size, extstep);
	
	notdone = 1;
	ostepsize = step_size;
	
	//if (current_phase == 1) {
	//step_size = min(0.85 * step_size * (tolrk / tr_error) ** 0.2, 4.0 * step_size)
	//
	//} else {
	//step_size = min(0.85 * step_size * (tolrk / tr_error) ** 0.2, 5.0 * step_size)
	//
	//}

	time_to_stim = shift - fmod(start_time, period);
	if (time_to_stim < 0.0) {
	  time_to_stim = time_to_stim + period;
	}
	if (time_to_stim < step_size) {
	  step_size = max(step_min, time_to_stim - 0.5 * step_min);
	}
      }
    } else {
      //puts("reject");
      

      if (step_size <= step_min) {
				//accept step, this should be rare
	for (i = 0; i < N; i++) {
	  if (fabs(y_1[i]) < min_tolvalue) {
	    y_1[i] = 0.0;
	  }
	  state[i] = y_1[i];
	}
	forcedaccept++;
	start_time = start_time + step_size;
	success = 1;
	if (start_time >= tstep) {
	  notdone = 0;
	} else {
	  step_size = step_min;
	  notdone = 1;
	}
	ostepsize = step_size;
	
      } else {
	// do not accept step
	// step_size < stimulus - start_time
	// notaccepted++;
	
	step_size = 0.85 * step_size * pow(tolrk / tr_error, 0.2);
	
	time_to_stim = shift - fmod(start_time, period);
	if (time_to_stim < 0.0) {
	  time_to_stim = time_to_stim + period;
	}
	if (time_to_stim < step_size) {
	  step_size = max(step_min, time_to_stim - 0.5 * step_min);
	}
	success = 0;
	notdone = 1;
	
	// resume saved state for FRUs
	send_resume_state();
      }
      
      //limits to step size
      step_size = min(step_max, max(step_min, step_size));
      step_size = min(step_size, tstep - start_time);
      //if (step_size < 1.e-10) {
      //notdone = 0;
      //
    }
    
    //puts("out");

    oldstepsize = ostepsize;
  }
  
  return oldstepsize;
}

/*

  pd5_predict_step_size

  We try to predict the size of step, to avoid redoing the stochastic
  simulation several times in a row

  The next algorithms is rk4m without stochastic code
  The only thing it is used for is estimation of step size
  It is accurate only when Ca currents are small
*/
double pd5_predict_step_size(double previous, double t, double tstep, 
                             double state[N],double current[Ncur],
			     double Jxfer, double Jtr, double ICa, double Ito2) 
{
	double          pred, st, st1 = 100;
	double          F[N],F1[N],y_1[N],k1[N],k2[N],k3[N],k4[N],k5[N],k6[N],ym[N];
	double          start_time, step_size, tr_error, errtmp;
	double          dummycurrent[1];

	int             i;
	int             notdone, success, keepc;
	int             dummyfalse = 0;
	double          os;
	const double    huge_number = 1.e99;

	double          time_to_stim;
	double          been_here = 0;

	os = previous;
	step_size = previous;
	start_time = t;
	notdone = 1;
	success = 1;
	keepc = 1;

	while (notdone) {

		//write(*, *) 't=', start_time, ' step=', step_size

		if (success) {

			fcn(start_time, state, F1, current, keepc, Jxfer, Jtr, ICa, Ito2);
			//Calculate global velocity field
			keepc = 0;

		} //(success)

		for (i = 0; i < N; i++) {
			//First intermediate RK4M step
			k1[i] = step_size * F1[i];
			y_1[i] = state[i] + k1[i] / 5.0;
		}

		fcn(start_time + step_size / 5.0, y_1, F, dummycurrent, dummyfalse, Jxfer, Jtr, ICa, Ito2);

		for (i = 0; i < N; i++) {
			//Second intermediate RK4M step
			k2[i] = step_size * F[i];
			y_1[i] = state[i] + (3.0 / 40.0 * k1[i] + 9.0 / 40.0 * k2[i]);
		}

		fcn(start_time + step_size * 3.0 / 10.0, y_1, F, dummycurrent, dummyfalse, Jxfer, Jtr, ICa, Ito2);

		for (i = 0; i < N; i++) {
			//Third intermediate RK4M step
			k3[i] = step_size * F[i];
			y_1[i] = state[i] + (3.0 / 10.0 * k1[i] - 9.0 / 10.0 * k2[i] + 6.0 / 5.0 * k3[i]);
		}

		fcn(start_time + step_size * 3.0 / 5.0, y_1, F, dummycurrent, dummyfalse, Jxfer, Jtr, ICa, Ito2);

		for (i = 0; i < N; i++) {
			//Forth intermediate RK4M step
			k4[i] = F[i] * step_size;
			y_1[i] = state[i] + 226.0 / 729.0 * k1[i] - 25.0 / 27.0 * k2[i] + 880.0 / 729.0 * k3[i] + 55.0 / 729.0 * k4[i];
		}

		fcn(start_time + step_size * 2.0 / 3.0, y_1, F, dummycurrent, dummyfalse, Jxfer, Jtr, ICa, Ito2);

		for (i = 0; i < N; i++) {
			//Fifth intermediate RK4M step
			k5[i] = step_size * F[i];
			y_1[i] = state[i] - 181.0 / 270.0 * k1[i] + 5.0 / 2.0 * k2[i] - 266.0 / 297.0 * k3[i] - 91.0 / 27.0 * k4[i] + 189.0 / 55.0 * k5[i];
		}

		fcn(start_time + step_size, y_1, F, dummycurrent, dummyfalse, Jxfer, Jtr, ICa, Ito2);

		for (i = 0; i < N; i++) {
			//Fifth intermediate RK4M step
			k6[i] = step_size * F[i];
		}

		//Calculation of truncation error.The value of tr_error is
		// the maximum normalized error(maximum over all states).

		tr_error = 0.0;
		for (i = 0; i < N; i++) {
			y_1[i] = state[i] + 31.0 / 540.0 * k1[i] + 190.0 / 297.0 * k3[i] - 145.0 / 108.0 * k4[i] + 351.0 / 220.0 * k5[i] + 1.0 / 20.0 * k6[i];
			ym[i] = state[i] + 19.0 / 216.0 * k1[i] + 1000.0 / 2079.0 * k3[i] - 125.0 / 216.0 * k4[i] + 81.0 / 88.0 * k5[i] + 5.0 / 56.0 * k6[i];
			errtmp = fabs((ym[i] - y_1[i]) * 0.2 * errweight[i]);
			if (errtmp < huge_number) {
				tr_error = max(errtmp, tr_error);
			} else {
				tr_error = huge_number;
			}
		}

		if (tr_error < tolrk) {
			if (been_here > 1) {
				if (step_size < st1)
					step_size = st1;
				success = 1;
				notdone = 0;
			} else {
				st1 = step_size;
				been_here = 1;
				notdone = 1;
				success = 0;
				st = 0.9 * step_size * pow(tolrk / tr_error, 0.2);
				if (step_size < st) {
					step_size = st;
				} else {
					success = 1;
					notdone = 0;
				}
			}
		} else {
			if (step_size <= step_min) {
				//success = 1;
				notdone = 0;
				if (start_time + step_size < tstep) {
					step_size = step_min;
				}
				//step_size = min(step_max, max(step_min, step_size))
				// step_size = min(step_size, tstep - start_time)
				// if (step_size < min_tolvalue) {
				//notdone = 0;

			} else {
				if (been_here > 0)
					been_here = 2;

				//step_size < stimulus - start_time
				time_to_stim = shift - fmod(start_time, period);
				if (time_to_stim < 0.0) {
					time_to_stim = time_to_stim + period;
				}
				step_size = 0.9 * step_size * pow(tolrk / tr_error,.2);

				//step_size = 0.95 * step_size * (tolrk / tr_error) **.2
				// if (time_to_stim < step_size) {
				//step_size = max(step_min, time_to_stim - (1.0e-10))
				//
				
				notdone = 1;
			}
			//step_size = min(step_max, max(step_min, step_size))
			// step_size = min(step_size, tstep - start_time)

			success = 0;

			if (step_size < min_tolvalue)
				notdone = 0;
		}
		pred = step_size;
	}

	return pred;
}
