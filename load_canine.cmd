#!/bin/ksh
#@ output = std.out
#@ error = std.err
#@ notification = never

#@ environment = MP_LABELIO=yes;
#@ environment = MP_INFOLEVEL=0;
#@ environment = MP_SHARED_MEMORY=yes;
#@ environment = MP_INTRDELAY=100;
#@ environment = MP_CPU_USAGE=unique;
#@ environment = MP_EUILIB = ip;
#@ environment = MP_CSS_INTERRUPT=yes;


#@ node = 1
#@ tasks_per_node=10
#@ class = par
#@ network.MPI=en2,shared,ip
#@ node_usage = shared
#@ wall_clock_limit = 100:00:00
#@ job_type = parallel
#@ checkpoint = no
#@ queue


/usr/bin/poe ./stoch


