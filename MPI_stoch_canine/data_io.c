/*       ----------------------------------------------------

	 NOTICE OF COPYRIGHT AND OWNERSHIP OF SOFTWARE
	 Copyright 2003, The Johns Hopkins University
	 School of Medicine. All rights reserved.

	 Name of Program: Local Control Model 
	 Version: Documented Version, C
	 Date: November 2003

	 --------------------------------------------------

	 data_io.c - Routines for input and output

*/

#include "parameters.h"
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

FILE *open_states(const int filenumber)
{
	FILE *filee;
	char filename[256];

	sprintf(filename,"%s%s.%d.txt",output_dir,output_states_file,filenumber);
	if ((filee=fopen(filename,"w+"))==NULL) {
		fprintf(stderr,"Cannot open file '%s'!\n",filename);
		exit(-2);
	}
	fprintf(filee,
		"%s %s %s %s %s %s %s %s %s %s " 
		"%s %s %s %s %s %s %s %s %s %s " 
		"%s %s %s %s %s %s %s %s %s %s "
		"%s %s %s %s %s %s %s %s\n",
		"Time","V","mNa","hNa","jNa","Nai","Ki","Cai","CaNSR","xKs",
		"LTRPNCa","HTRPNCa","C0Kv43","C1Kv43","C2Kv43","C3Kv43","OKv43","CI0Kv43","CI1Kv43","CI2Kv43",
		"CI3Kv43","OIKv43","C0Kv14","C1Kv14","C2Kv14","C3Kv14","OKv14","CI0Kv14","CI1Kv14","CI2Kv14",
		"CI3Kv14","OIKv14","CaToT","C1Herg","C2Herg","C3Herg","OHerg","IHerg");

	return filee;
}

FILE *open_otherstates(const int filenumber)
{
	// Open output files for data in 'states' and 'otherstate'
	FILE *filee;
	char filename[256];

	sprintf(filename,"%s%s.%d.txt",output_dir,output_otherstates_file,filenumber);
	if ((filee=fopen(filename,"w+"))==NULL) {
		fprintf(stderr,"Cannot open file '%s'!\n",filename);
		exit(-3);
	}
	fprintf(filee,
		"%s %s %s %s %s %s %s %s %s %s " 
		"%s %s %s \n",
		"Time","CaSSavg","CaJSRavg","JRyRtot","PRyR_open","PRyR_ready","PNorm_mode","PnotVinact","PLType_open","CaToT2",
		"PIto2_open","CaJSRtot","CaSStot");

	return filee;
}


FILE *open_currents(const int filenumber) // Open output files for data in 'current'
{
	FILE *filee;
	char filename[256];

	sprintf(filename,"%s%s.%d.txt",output_dir,output_currents_file,filenumber);
	if ((filee=fopen(filename,"w+"))==NULL) {
		fprintf(stderr,"Cannot open file '%s'!\n",filename);
		exit(-5);
	}

	fprintf(filee,
		"%s %s %s %s %s %s %s %s %s %s " 
		"%s %s %s %s %s %s %s %s %s %s " 
		"%s %s %s %s %s \n",
		"Time","INa","IKr","IKs","Ito1","IK1","IKp","INaCa","INaK","IpCa",
		"ICab","INab","ICa","JDHPR","Jup","Jtrpn","Jtr","Jxfer","IKv43","IKv14",
		"IKv14_K","IKv14_Na","Ito2","Istim","Itot");


	return filee;
}

// Write data to output files for data in 'states'
void write_states(FILE *file,const double time,const double state[N])
{
  double tprint;
  int i;
  
  tprint = time;
  if (ts_sec) tprint = tprint * 1.e-3;
  
  fprintf(file,"%g ",tprint);
  
  // Write states
  for(i=0;i<N;i++) {
    fprintf(file,"%.8g ",state[i]);
  }
  fprintf(file,"\n");
}

// Write data to output files for data in 'otherstate'
void write_otherstates(FILE *file,const double time,const double otherstates[Nother])
{
  double tprint;
  int i;
  
  tprint = time;
  if (ts_sec) tprint = tprint * 1.e-3;
  fprintf(file,"%g ",tprint);
  
  // Write otherst
  for(i=0;i<Nother;i++) {
    fprintf(file,"%g ",otherstates[i]);
  }
  fprintf(file,"\n");
}

// Write data to output files for data in 'current'
void write_currents(FILE *file,const double time,const double current[Ncur])
{
	double tprint;
	int i;

	tprint = time;
	if (ts_sec) tprint = tprint * 1.e-3;

	fprintf(file,"%g ",tprint);

	// Write states
	for(i=0;i<Ncur;i++) {
		fprintf(file,"%g ",current[i]);
	}
	fprintf(file,"\n");
}

double read_next_double(FILE *file)
{
	unsigned char buf[256];
	unsigned char ch;
	int i=0;
	double val;

	do { 
	  ch = fgetc( file );
	} while (isspace(ch)&&(!feof(file)));

	while (isdigit(ch)||(ch=='E')||(ch=='e')||(ch=='-')||(ch=='+')||(ch=='.')) {
	  buf[i]=ch;
	  i++;
	  ch = fgetc( file );
	}
	buf[i]=0;

	if (i<1) {
	  fprintf(stderr,"Error reading ic file: no decent input!\n");
	  val=0;
	} else {
	  val=strtod((const char *)buf,NULL);
	}
	return val;
}

int read_next_int(FILE *file)
{
	unsigned char buf[256];
	unsigned char ch;
	int i=0;
	int val;

	do { 
	  ch = fgetc( file );
	} while (isspace(ch)&&(!feof(file)));

	while ((isdigit(ch)||(ch=='-')||(ch=='+'))&&(i<255)) {
	  buf[i]=ch;
	  i++;
	  ch = fgetc( file );
	}
	buf[i]=0;

	if ((i<1)||(i>30)) {
	  fprintf(stderr,"Error reading ic file: no decent input!\n");
	  val=0;
	} else {
	  val=(int)strtol((const char *)buf,NULL,10);
	}
	return val;
}

// Write data to all output files for local release unit data
void restart_data(const int filenumber,
		  double state[N],
		  double FRU_states[NFRU_sim_max][Nstates_FRU],
		  int LType_state[NFRU_sim_max][Nclefts_FRU][Nindepstates_LType],
		  int RyR_state[NFRU_sim_max][Nclefts_FRU][NRyRs_per_cleft],
		  int Ito2_state[NFRU_sim_max][Nclefts_FRU],
		  unsigned long mt[NFRU_sim_max][mtN+1],int mti[NFRU_sim_max])
{
	int iFRU,icleft,k;
	char filename[256];
	FILE *states_file;
	FILE *RyR_file;
	FILE *LCh_file;
	FILE *FRU_file;
	FILE *Ito2_file;

	sprintf(filename,"%sr_states.%d.txt",output_dir,filenumber);
	if ((states_file=fopen(filename,"w+"))==NULL) {
		fprintf(stderr,"Problem opening file %s\n",filename);
	}
	for( k=0;k<N;k++) {
		fprintf(states_file,"%g\n",state[k]); 
	}
	fclose(states_file);

	sprintf(filename,"%sr_FRU.%d.txt",output_dir,filenumber);
	if ((FRU_file=fopen(filename,"w+"))==NULL) {
		fprintf(stderr,"Problem opening file %s\n",filename);
	}

	sprintf(filename,"%sr_LCh.%d.txt",output_dir,filenumber);
	if ((LCh_file=fopen(filename,"w+"))==NULL) {
		fprintf(stderr,"Problem opening file %s\n",filename);
	}

	sprintf(filename,"%sr_RyR.%d.txt",output_dir,filenumber);
	if ((RyR_file=fopen(filename,"w+"))==NULL) {
		fprintf(stderr,"Problem opening file %s\n",filename);
	}

	sprintf(filename,"%sr_Ito2.%d.txt",output_dir,filenumber);
	if ((Ito2_file=fopen(filename,"w+"))==NULL) {
		fprintf(stderr,"Problem opening file %s\n",filename);
	}

	// Add number of FRUs used at the beginning of each file
	// except states
	fprintf(FRU_file,"%d\n",NFRU_sim);
	fprintf(RyR_file,"%d\n",NFRU_sim);
	fprintf(LCh_file,"%d\n",NFRU_sim);
	fprintf(Ito2_file,"%d\n",NFRU_sim);

	for(iFRU=0;iFRU<NFRU_sim;iFRU++) {
	  for(k=0;k<Nstates_FRU;k++) {
	    fprintf(FRU_file,"%g ",FRU_states[iFRU][k]);
	  }
	  fprintf(FRU_file,"\n");
	  for( icleft = 0;icleft<Nclefts_FRU;icleft++) {
	    fprintf(LCh_file,"%d %d\n",LType_state[iFRU][icleft][index_LCC_states], LType_state[iFRU][icleft][index_LCC_Vinact]);
	    for(k=0;k<NRyRs_per_cleft;k++) {
	      fprintf(RyR_file,"%d ",RyR_state[iFRU][icleft][k]);
	    }
	    fprintf(RyR_file,"\n");
	    fprintf(Ito2_file,"%d\n",Ito2_state[iFRU][icleft]);
	  } 
	}

	fclose(FRU_file);
	fclose(states_file);
	fclose(RyR_file);
	fclose(LCh_file);
	fclose(Ito2_file);

	if (save_seeds_file_flag) {
	  unsigned long buf[mtN+2];
	  FILE *seeds_file;
	  int i,numwritten;
	  
	  sprintf(filename,"%s%s.%d.txt",output_dir,output_seeds_file,filenumber);
	  if ((seeds_file=fopen(filename,"wb"))==NULL) {
	    fprintf(stderr,"Problem opening file %s\n",filename);
	    abort();
	  }
	  
	  for(iFRU=0;iFRU<NFRU_sim;iFRU++) {		
	    buf[0]=(unsigned long)iFRU;
	    buf[1]=(unsigned long)mti[iFRU];
	    for(i=0;i<mtN;i++) {
	      buf[i+2]=mt[iFRU][i];
	    }
	    numwritten=fwrite((void *)buf,sizeof(unsigned long),mtN+2,seeds_file);
	    if (numwritten<mtN+2) {
	      fprintf(stderr,"Unable to save seeds\n");
	    }
	  }
	  fclose(seeds_file);
	}
}


