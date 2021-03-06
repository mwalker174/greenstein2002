/*       ----------------------------------------------------

	 NOTICE OF COPYRIGHT AND OWNERSHIP OF SOFTWARE
	 Copyright 2003, The Johns Hopkins University
	 School of Medicine. All rights reserved.

	 Name of Program: Local Control Model 
	 Version: Documented Version, C
	 Date: November 2003

	 --------------------------------------------------

	 mpi_slave - routines for MPI slave processes

*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>

#include "parameters.h"

#if USE_MPI
#include <mpi.h>
#endif

//#include "parameters_fcn_fru.h" 

int NFRU_local;

double FRU_states_local[MAX_LOCAL_NFRU][Nstates_FRU];
int LType_state_local[MAX_LOCAL_NFRU][Nclefts_FRU][Nindepstates_LType];
int RyR_state_local[MAX_LOCAL_NFRU][Nclefts_FRU][NRyRs_per_cleft];
int Ito2_state_local[MAX_LOCAL_NFRU][Nclefts_FRU];
unsigned long mt_local[MAX_LOCAL_NFRU][mtN+1];
int mti_local[MAX_LOCAL_NFRU];

double FRU_states_local_hold[MAX_LOCAL_NFRU][Nstates_FRU];
int LType_state_local_hold[MAX_LOCAL_NFRU][Nclefts_FRU][Nindepstates_LType];
int RyR_state_local_hold[MAX_LOCAL_NFRU][Nclefts_FRU][NRyRs_per_cleft];
int Ito2_state_local_hold[MAX_LOCAL_NFRU][Nclefts_FRU];
unsigned long mt_local_hold[MAX_LOCAL_NFRU][mtN+1];
int mti_local_hold[MAX_LOCAL_NFRU];

#if USE_MPI
/* 
   void parallel_init_FRU()

   initializes actual MPI processes

*/
void parallel_init_FRU(void)
{
  int source;
  MPI_Request request[100];
  MPI_Status status[100];
  int ok=1;
  int tag;
  int i,j,reqno;

  source=0; // receive from master
  tag=1;

  // first receive number of FRU:s
  MPI_Recv(&NFRU_local,1,MPI_INT,source,tag,MPI_COMM_WORLD,&status[0]); 
  // We need the result before we can proceed

  for(i=0;i<NFRU_local;i++) {
    
    tag=i*11+my_rank*1000;
    reqno=0;
  
    MPI_Irecv(&LType_state_local[i][0][0],Nclefts_FRU*Nindepstates_LType,MPI_INT,source,tag,MPI_COMM_WORLD,&request[reqno]);
    tag++; reqno++;
    MPI_Irecv(&RyR_state_local[i][0][0],Nclefts_FRU*NRyRs_per_cleft,MPI_INT,source,tag,MPI_COMM_WORLD,&request[reqno]);
    tag++; reqno++;
    MPI_Irecv(&Ito2_state_local[i][0],Nclefts_FRU,MPI_INT,source,tag,MPI_COMM_WORLD,&request[reqno]);
    tag++; reqno++;
    MPI_Irecv(&FRU_states_local[i][0],Nstates_FRU,MPI_DOUBLE,source,tag,MPI_COMM_WORLD,&request[reqno]);
    tag++; reqno++;
    MPI_Irecv(&mti_local[i],1,MPI_INT,source,tag,MPI_COMM_WORLD,&request[reqno]);
    tag++; reqno++;
    MPI_Irecv(&mt_local[i][0],mtN,MPI_UNSIGNED_LONG,source,tag,MPI_COMM_WORLD,&request[reqno]);
    MPI_Waitall(reqno+1,request,status);
  }

  //puts("receive ok");

  // send OK
  ok=1;
  tag=1;
  MPI_Send(&ok,1,MPI_INT,source,tag,MPI_COMM_WORLD);
}

/* 
   void parallel_send_FRU(void)

   sends FRUs to the master process

*/
void parallel_send_FRU(void)
{
  int dest;
  MPI_Request request[100];
  MPI_Status status[100];
  int tag,i;
  int reqno=0;

  dest=0; // receive from master
  tag=my_rank;

  // first send number of FRUs
  MPI_Isend(&NFRU_local,1,MPI_INT,dest,tag,MPI_COMM_WORLD,&request[reqno]); 

  // then send individual FRUs
  for(i=0;i<NFRU_local;i++) {
    tag=my_rank*1000+i;
    MPI_Isend(&LType_state_local[i][0][0],Nclefts_FRU*Nindepstates_LType,MPI_INT,dest,tag,MPI_COMM_WORLD,&request[reqno]);
    reqno++;
    tag++;
    MPI_Isend(&RyR_state_local[i][0][0],Nclefts_FRU*NRyRs_per_cleft,MPI_INT,dest,tag,MPI_COMM_WORLD,&request[reqno]);
    reqno++;
    tag++;
    MPI_Isend(&Ito2_state_local[i][0],Nclefts_FRU,MPI_INT,dest,tag,MPI_COMM_WORLD,&request[reqno]);
    reqno++;
    tag++;
    MPI_Isend(&FRU_states_local[i][0],Nstates_FRU,MPI_DOUBLE,dest,tag,MPI_COMM_WORLD,&request[reqno]);
    reqno++;
    tag++;
    MPI_Isend(&mti_local[i],1,MPI_INT,dest,tag,MPI_COMM_WORLD,&request[reqno]);
    reqno++;
    tag++;
    MPI_Isend(&mt_local[i][0],mtN,MPI_UNSIGNED_LONG,dest,tag,MPI_COMM_WORLD,&request[reqno]);
    reqno++;
    
    MPI_Waitall(reqno,request,status);
    reqno=0;
  }
}
#endif

/* 
   void parallel_resume_state(void)

   resumes states from safe

*/
void parallel_resume_state(void)
{
  int i;
  int iFRU;
  int icleft;

#if NOMEMCPY
  for(iFRU=0;iFRU<NFRU_local;iFRU++) {
    for (i = 0; i < Nstates_FRU; i++) {
      FRU_states_local[iFRU][i] = FRU_states_local_hold[iFRU][i];
    }
    for (icleft = 0; icleft < Nclefts_FRU; icleft++) {
      for (i = 0; i < NRyRs_per_cleft; i++) {
	RyR_state_local[iFRU][icleft][i] = RyR_state_local_hold[iFRU][icleft][i];
      }
      LType_state_local[iFRU][icleft][index_LCC_states] = LType_state_local_hold[iFRU][icleft][index_LCC_states];
      LType_state_local[iFRU][icleft][index_LCC_Vinact] = LType_state_local_hold[iFRU][icleft][index_LCC_Vinact];
      Ito2_state_local[iFRU][icleft] = Ito2_state_local_hold[iFRU][icleft];
    }

    mti_local[iFRU] = mti_local_hold[iFRU];
    for (i = 0; i < mtN; i++) {
      mt_local[iFRU][i] = mt_local_hold[iFRU][i];
    }
  }
#else
  // is this faster ? does not seem to be...
  memcpy(FRU_states_local, FRU_states_local_hold, sizeof(double) * NFRU_local * Nstates_FRU);
  memcpy(RyR_state_local, RyR_state_local_hold, sizeof(int) * NFRU_local * Nclefts_FRU * NRyRs_per_cleft);
  memcpy(LType_state_local, LType_state_local_hold, sizeof(int) * NFRU_local * Nclefts_FRU * Nindepstates_LType);
  memcpy(Ito2_state_local, Ito2_state_local_hold, sizeof(int) * NFRU_local * Nclefts_FRU);
  memcpy(mti_local, mti_local_hold, sizeof(int) * NFRU_local);
  memcpy(mt_local, mt_local_hold, sizeof(unsigned long) * NFRU_local * (mtN+1));
#endif  
}

/* 
   void parallel_save_state(void)

   saves current state of affairs, in case it will have to be resumed

*/
void parallel_save_state(void)
{
  int i;
  int icleft;
  int iFRU;

#if NOMEMCPY
  for(iFRU=0;iFRU<NFRU_local;iFRU++) {
    for (i = 0; i < Nstates_FRU; i++) {
      FRU_states_local_hold[iFRU][i] = FRU_states_local[iFRU][i];
    }
    for (icleft = 0; icleft < Nclefts_FRU; icleft++) {
      for (i = 0; i < NRyRs_per_cleft; i++) {
	RyR_state_local_hold[iFRU][icleft][i] = RyR_state_local[iFRU][icleft][i];
      }
      LType_state_local_hold[iFRU][icleft][index_LCC_states] = LType_state_local[iFRU][icleft][index_LCC_states];
      LType_state_local_hold[iFRU][icleft][index_LCC_Vinact] = LType_state_local[iFRU][icleft][index_LCC_Vinact];
      Ito2_state_local_hold[iFRU][icleft] = Ito2_state_local[iFRU][icleft];
    }

    mti_local_hold[iFRU] = mti_local[iFRU];
    for (i = 0; i < mtN; i++) {
      mt_local_hold[iFRU][i] = mt_local[iFRU][i];
    }
  }
#else
  // is this faster ? does not seem to be...
  memcpy(FRU_states_local_hold, FRU_states_local, sizeof(double) * NFRU_local * Nstates_FRU);
  memcpy(RyR_state_local_hold, RyR_state_local, sizeof(int) * NFRU_local * Nclefts_FRU * NRyRs_per_cleft);
  memcpy(LType_state_local_hold, LType_state_local, sizeof(int) * NFRU_local * Nclefts_FRU * Nindepstates_LType);
  memcpy(Ito2_state_local_hold, Ito2_state_local, sizeof(int) * NFRU_local * Nclefts_FRU);
  memcpy(mti_local_hold, mti_local, sizeof(int) * NFRU_local);
  memcpy(mt_local_hold, mt_local, sizeof(unsigned long) * NFRU_local * (mtN+1));
#endif
}

/* 
   void parallel_compute_simfru(void)
   
   does the actual computations

*/
void parallel_compute_simfru(void)
{
  int i,dest;
  double start_time_local,end_time_local;
  double FRUdep_states0_local[Nstates_FRUdep];
  double FRUdep_statesf_local[Nstates_FRUdep];
  double num_stat[6],total_num[6];
  double input[2+2*Nstates_FRUdep];

#if USE_MPI
  // Minimize number of messages
  MPI_Bcast(input,2+2*Nstates_FRUdep,MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

  start_time_local=input[0];
  end_time_local=input[1];
  for(i=0;i<Nstates_FRUdep;i++) {
    FRUdep_states0_local[i]=input[2+i];
    FRUdep_statesf_local[i]=input[2+Nstates_FRUdep+i];
  }

  //printf("MPI error %d with proc %d, tag %d\n",status.MPI_ERROR,my_rank,tag);

  for(i=0;i<NFRU_local;i++) {
    simfru_local(start_time_local,end_time_local,
		 FRUdep_states0_local,
		 FRUdep_statesf_local,
		 &LType_state_local[i][0][0],
		 &RyR_state_local[i][0][0],
		 &Ito2_state_local[i][0],
		 &FRU_states_local[i][0],
		 &mti_local[i],
		 &mt_local[i][0]);
  }

  calc_fru_flux_local(NFRU_local,LType_state_local,Ito2_state_local,
		      FRU_states_local,num_stat);

#if USE_MPI
  dest=0;
  for(i=0;i<5;i++) total_num[i]=0;

  MPI_Reduce(num_stat,total_num,5,MPI_DOUBLE,MPI_SUM,dest,MPI_COMM_WORLD);
#endif
}

#if USE_MPI
/* 
   void parallel_calc_fru_avg(void)

   computes calc_fru_avg

*/
void parallel_calc_fru_avg(void)
{
  double num_stat[11];
  double total_num[11];
  int i;

  if (NFRU_local>0) {
    calc_fru_avg_local(num_stat,NFRU_local,FRU_states_local,LType_state_local,RyR_state_local,Ito2_state_local);
  } else {
    int i;
    
    for(i=0;i<11;i++)
      num_stat[i]=0;
  }

  for(i=0;i<11;i++) total_num[i]=0;
  MPI_Reduce(num_stat,total_num,11,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
}
#endif

/* 
   void parallel_calc_fru_flux(void)

   computes calc_fru_flux

*/
void parallel_calc_fru_flux(void)
{
  int i,dest;
  //int DepFlag;
  //double FRUdep_states_local[Nstates_FRUdep];
  double num_stat[7],total_num[7];

  //printf("proc %d calc fru flux\n",my_rank);

#if USE_MPI
  //MPI_Bcast(FRUdep_states_local,Nstates_FRUdep,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif  

  //DepFlag = 1;
  if (NFRU_local>0) {
    calc_fru_flux_local(NFRU_local,LType_state_local,Ito2_state_local,FRU_states_local,num_stat);
  } else {
    for(i=0;i<6;i++)
      num_stat[i]=0;
  }

#if USE_MPI
  dest=0;
  for(i=0;i<5;i++) total_num[i]=0;
  MPI_Reduce(num_stat,total_num,5,MPI_DOUBLE,MPI_SUM,dest,MPI_COMM_WORLD);
#endif
}

/* 
   void parallel(void)

   main loop, where slave processes sit unless they have work to do
   ie this is message handler

*/
void parallel(void)
{
#if USE_MPI
  int tag,source,task,cont;
  
  cont=1;

  //printf("Proc %d enters parallel\n",my_rank);
  
  do {
    tag=my_rank*11;
    source=0;
    task=-1;
    //printf("Proc %d waiting for message\n",my_rank);
    // We broadcast since it is the fastest method
    MPI_Bcast(&task,1,MPI_INT,0,MPI_COMM_WORLD);
    //printf("MPI error %d with proc %d\n",status.MPI_ERROR,my_rank);
    //printf("Proc %d received %d\n",my_rank,task);

    switch (task) {
    case MSG_INIT: // init
      printf("Proc %d received init!\n",my_rank);
      break;
    case MSG_EXIT: // exit
      cont=0;
      break;
    case MSG_COMP_SIMFRU: // simfru computation
      parallel_compute_simfru();
      break;
    case MSG_RESUME_STATE: // go back to previous state
      parallel_resume_state();
      break;
    case MSG_SAVE_STATE: // save state
      parallel_save_state();
      break;
    case MSG_INIT_FRU: // init FRU
      parallel_init_FRU();
      break;
    case MSG_RETURN_FRU: // return FRU
      parallel_send_FRU();
      break;
    case MSG_CALC_FRU_AVG: // calc_fru_avg_local
      parallel_calc_fru_avg();
      break;
    case MSG_CALC_FRU_FLUX: // calc_fru_flux
      parallel_calc_fru_flux();
      break;
    default:
      printf("Unknown task %d in proc %d!\n",task,my_rank);
    }
  } while (cont);
  //printf("Proc %d exits.\n",my_rank);
  MPI_Finalize();
#endif
}

 
