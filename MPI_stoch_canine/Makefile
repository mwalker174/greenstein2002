# .KEEP_STATE:
# .SUFFIXES: .o .c
# .SUFFIXES: .o .c .h

OBJS =  caf.o calc_fru_flux_local.o \
	data_io.o fcn.o calc_fru_avg_local.o \
        fru_rates_local.c initialize_ran.o initialize_state.o \
        simfru_local.o write_fru_props.o mt19937_local.o rk_pd54m.o \
	fcn_fru.o dynamicFRU.o mpi_slave.o mpi_master.o parameters_fcn_fru.o

# For Linux with GCC
#
#COMP = gcc
#LIBS = -lm
#OPT_FLAGS = -O3
#LARGE_FLAGS = 
#DEBUG = 

# For IBM AIX with xlc (mpcc is xlc with MPI)
COMP = mpcc
LIBS = -lm -lmass #-bnso -bI:/lib/syscalls.exp # static linking 
OPT_FLAGS=-O5 -qarch=auto -qtune=auto -qcache=auto -qstrict 
COMP_FLAGS = -qcpluscmt $(OPT_FLAGS) $(LARGE_FLAGS)
LARGE_FLAGS= 
PROFILING = 
DEBUG	=

# For IBM AIX with xlc (mpcc is xlc with MPI)
#COMP = mpcc
#LIBS = -lm -lmass #-bnso -bI:/lib/syscalls.exp # static linking 
#OPT_FLAGS=
#LARGE_FLAGS=-bmaxdata:0x80000000 -bmaxstack:0x80000000 
#COMP_FLAGS = -qcpluscmt $(OPT_FLAGS) $(LARGE_FLAGS)
#PROFILING = -pg
#DEBUG	= -g -C -qinfo=dcl:lan:par:ret:gen:ord:cnv:cnd:por:pro:tru:c99:eff:ext:gnr:trd -qextchk -qcheck=all -qfloat=fltint -qflttrap=inexact:underflow:overflow:invalid:imprecise:zerodivide $(PROFILING) -qfloat=rndsngl:fltint:nans #-bnoquiet

# The name of the executable
NAME=stoch

##########################################
#               main program             #
##########################################

all:	$(OBJS)
	$(COMP) $(COMP_FLAGS) $(DEBUG) $(OBJS) $(LIBS) -o $(NAME)


unix:
	dos2unix `find . -type f`

mt:	mttesti19937.o
	$(COMP) $(COMP_FLAGS) $(DEBUG) mttesti19937.o $(LIBS) -o MT.exe

makeit:	$(OBJS)
	$(COMP) $(COMP_FLAGS) $(DEBUG) $(OBJS) $(LIBS) -o $(NAME)

clean:
	rm -f *.o $(NAME)

cleanall:
	rm -f *.o $(NAME) current* state* other* r_* nohup.out fruprops_*

###########################################
#               subroutines               #
###########################################

%.o: %.c
	$(COMP) -c $(COMP_FLAGS) $(DEBUG) $<

caf.o: caf.c  parameters.h
	$(COMP) -c $(COMP_FLAGS) $(DEBUG) $<

fcn.o: fcn.c parameters.h
	$(COMP) -c $(COMP_FLAGS) $(DEBUG) $<

rk4am.o: rk4am.c  parameters.h
	$(COMP) -c $(COMP_FLAGS) $(DEBUG) $<

calc_fru_avg_local.o: calc_fru_avg_local.c parameters.h
	$(COMP) -c $(COMP_FLAGS) $(DEBUG) $<

calc_fru_flux_local.o: calc_fru_flux_local.c parameters.h
	$(COMP) -c $(COMP_FLAGS) $(DEBUG) $<

fru_rates.o: fru_rates.c parameters.h
	$(COMP) -c $(COMP_FLAGS) $(DEBUG) $<

fru_rates_local.o: fru_rates_local.c parameters.h
	$(COMP) -c $(COMP_FLAGS) $(DEBUG) $<

initialize_ran.o: initialize_ran.c parameters.h
	$(COMP) -c $(COMP_FLAGS) $(DEBUG) $<

initialize_state.o: initialize_state.c parameters.h
	$(COMP) -c $(COMP_FLAGS) $(DEBUG) $<

mpi_slave.o: mpi_slave.c parameters.h
	$(COMP) -c $(COMP_FLAGS) $(DEBUG) $<

write_fru_props.o: write_fru_props.c parameters.h
	$(COMP) -c $(COMP_FLAGS) $(DEBUG) $<

data_io.o: data_io.c parameters.h
	$(COMP) -c $(COMP_FLAGS) $(DEBUG) $<

simfru_local.o: simfru_local.c parameters.h
	$(COMP) -c $(COMP_FLAGS) $(DEBUG) $<

mt19937_local.o: mt19937_local.c parameters.h
	$(COMP) -c $(COMP_FLAGS) $(DEBUG) $<

rk_pd54m.o: rk_pd54m.c parameters.h
	$(COMP) -c $(COMP_FLAGS) $(DEBUG) $<

fcn_fru.o: fcn_fru.c parameters.h
	$(COMP) -c $(COMP_FLAGS) $(DEBUG) $<

dynamicFRU.o: dynamicFRU.c parameters.h
	$(COMP) -c $(COMP_FLAGS) $(DEBUG) $<

mpi_master.o: mpi_master.c parameters.h
	$(COMP) -c $(COMP_FLAGS) $(DEBUG) $<

parameters_fcn_fru.o: parameters_fcn_fru.c parameters.h
	$(COMP) -c $(COMP_FLAGS) $(DEBUG) $<
