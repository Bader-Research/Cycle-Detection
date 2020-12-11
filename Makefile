# $Id: Makefile,v 1.15 1999/09/25 16:37:00 dbader Exp $
#
# $Log: Makefile,v $
# Revision 1.15  1999/09/25 16:37:00  dbader
# Changed CYGWIN_NT-4.0 to CYGWIN_NT_40.
# Added sed to OS shell uname -s to convert - to _ and remove "."'s.
#
# Revision 1.14  1999/09/25 16:09:57  dbader
# Added support for CYGWIN_NT-4.0
#
# Revision 1.13  1999/09/23 15:38:13  dbader
# Added new CYGWIN32_NT include and lib information.
#
# Revision 1.12  1999/07/24 12:49:07  dbader
# Added portability to FreeBSD.
#
# Revision 1.11  1999/07/24 11:36:34  dbader
# Added Intel Paragon build instructions.
#
# Revision 1.10  1999/07/23 16:12:22  dbader
# Modified to use CYGWIN32_NT
#
# Revision 1.9  1999/07/15 12:15:36  dbader
# Added automatic OS detection.
#
# Revision 1.8  1999/07/15 01:30:19  dbader
# Added MALLOC_DEBUG support.
#
# Revision 1.7  1999/07/12 21:27:46  dbader
# Added cycleInterval.
#
# Revision 1.6  1999/07/04 20:24:00  dbader
# Fixed minor dependencies problem.
#
# Revision 1.5  1999/07/03 23:00:39  dbader
# Added cycleConvex.
#
# Revision 1.4  1999/06/05 21:38:37  dbader
# Added -Wall compilation flag.
#
# Revision 1.3  1999/06/04 12:58:16  dbader
# Added -DDEBUG_PRINT option and pretty-print
#
# Revision 1.2  1999/06/01 19:50:28  dbader
# Added timing.[ch]
#
# Revision 1.1  1999/05/27 19:09:08  dbader
# Initial revision
#
#

EXECS     = cycle cycleConvex cycleInterval

# Uncomment this next line to turn on automatic MALLOC debugging
# DEBUG_MALLOC = YES

# Uncomment this next line to turn on debugging print statements
# DEBUG_PRINT = YES

OBJLIST	    = queue misc mpi-printf timing
OBJS        = $(addsuffix .$(OBJSUFFIX), $(OBJLIST))
LIBS	    = -lm

OS          = $(shell uname -s | sed -e "s/\-/_/g" -e "s/\.//g")
CFLAGS      = $(OPT) -D_$(OS)

ifeq ($(shell hostname), siesta)
OS	    = Paragon
endif

ifeq ($(OS), SunOS)
CFLAGS     += -Wall
OPT         = -O2
CC          = /research/red/dbader/MPI/mpich/bin/mpicc
OBJSUFFIX   = o
endif

ifeq ($(OS), Linux)
OPT         = -fast
CC          = /usr/parallel/mpich-gm.pgi/bin/mpicc
OBJSUFFIX   = o
endif

ifeq ($(OS), FreeBSD)
CFLAGS     += -Wall
OPT         = -O2
CC          = /usr/local/mpich/bin/mpicc
OBJSUFFIX   = o
endif

ifeq ($(OS), CYGWIN32_NT)
CFLAGS     += -I "d:/apps/hpvm/myrinet/include" 
LIBS	   += mpi.lib mpicref.lib FM.lib myrilib.lib \
	      shell32.lib advapi32.lib wsock32.lib mkl_m.lib
OPT         = -O2
#CC          = gcc
CC 	    = cl
OBJSUFFIX   = obj
EXECS_MS    = $(addsuffix .exe, $(EXECS))
endif

ifeq ($(OS), CYGWIN_NT_40)
CFLAGS     += -I "d:/apps/hpvm/include" 
LIBS	   += mpi.lib fm.lib advapi32.lib kernel32.lib wsock32.lib
OPT         = -O2
#CC          = gcc
CC 	    = cl
OBJSUFFIX   = obj
EXECS_MS    = $(addsuffix .exe, $(EXECS))
endif

ifeq ($(OS), Paragon)
CFLAGS	   += -I/Net/local/mpi/include -D_NX
OPT         = -O2
CC	    = /Net/local/sunmos/current/bin/sicc
LIBS       += -L/Net/local/mpi/lib/paragon/sunmos -lmpi
endif

ifdef DEBUG_MALLOC
OBJS       += rmalloc.$(OBJSUFFIX)
CFLAGS     += -DMALLOC_DEBUG
endif

ifdef DEBUG_PRINT
CFLAGS     += -DDEBUG_PRINT
endif

all: $(EXECS)

# $(EXECS): $(OBJS)
#	$(CC) $(CFLAGS) -o $@ $^ $(LIBS) 

cycle: cycle.$(OBJSUFFIX) $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS) 

cycleConvex: cycleConvex.$(OBJSUFFIX) $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS) 

cycleInterval: cycleInterval.$(OBJSUFFIX) $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS) 

tar: cycle.tar.gz

cycle.tar :
	tar cvf $@ *.[ch] Makefile

cycle.tar.gz : cycle.tar
	rm -f $@
	gzip $<

.SUFFIXES: .c

.c.o : 
	$(CC) $(CFLAGS) -c $<

%.$(OBJSUFFIX) : %.c
	$(CC) $(CFLAGS) -c $<


clean:
	rm cycle.$(OBJSUFFIX) \
	cycleConvex.$(OBJSUFFIX) \
	cycleInterval.$(OBJSUFFIX) \
	$(OBJS) $(EXECS) $(EXECS_MS) *~ PI* *.clog 

prettyprint:
	c2ps -UDEBUG_PRINT *.[ch] Makefile | lpr

cycle.c: cycle.h misc.h queue.h mpi-printf.h timing.h
cycleConvex.c: cycleConvex.h misc.h queue.h mpi-printf.h timing.h
cycleInterval.c: cycleInterval.h misc.h queue.h mpi-printf.h timing.h
queue.c: queue.h
queue.h: misc.h
misc.c: misc.h
mpi-printf.c: mpi-printf.h
timing.c: timing.h
