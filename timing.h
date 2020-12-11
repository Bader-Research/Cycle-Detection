#ifndef _TIMING_H
#define _TIMING_H

/*
  $Id: timing.h,v 1.10 1999/09/25 16:37:00 dbader Exp $

  $Log: timing.h,v $
  Revision 1.10  1999/09/25 16:37:00  dbader
  Changed CYGWIN_NT-4.0 to CYGWIN_NT_40.

  Revision 1.9  1999/09/25 16:09:57  dbader
  Added support for CYGWIN_NT-4.0

  Revision 1.8  1999/07/24 12:50:14  dbader
  Changed logging header.

*/

#include <string.h>
#include <sys/time.h>           /* struct timeval */

#ifdef DEBUG_MALLOC
#include "rmalloc.h"
#endif

#if defined(_Linux)||defined(_SunOS)||defined(_FreeBSD)
#define TIME_GTOD
#endif

#if defined(_NX)
#define TIME_NX
#endif

#if defined(_CYGWIN32_NT)||defined(_CYGWIN_NT_40)
#define TIME_MPI
#endif

#ifdef TIME_GC
struct timespec tp;
#define get_seconds()   (getclock(TIMEOFDAY, &tp), \
                        (double)tp.tv_sec + (double)tp.tv_nsec / 1000000000.0)

#define get_nseconds()  (getclock(TIMEOFDAY, &tp), \
                        (double)(1000000000*tp.tv_sec) + (double)tp.tv_nsec)
#endif

#ifdef TIME_GTOD
double get_seconds();
#endif

#ifdef TIME_MPI
#define get_seconds() MPI_Wtime()
#endif

#ifdef TIME_NX
#include <nx.h>
#define get_seconds() dclock()
#endif

#define TIME_STEPS   25
#define STEP_MAXLEN  80
#define TSTRLEN      25

/******************************************************************/

static int _timer_level = 0;

#define timer_init() \
  int step; \
  double _secs[TIME_STEPS]; \
  double _tsec[TIME_STEPS]; \
  char   _stepname[TIME_STEPS][STEP_MAXLEN];
#define timer_start() \
  MPI_Barrier(MPI_COMM_WORLD); \
  _secs[step] = get_seconds(); {}
#define timer_reset() \
  step = 0; \
  _timer_level++; \
  MPI_Barrier(MPI_COMM_WORLD); \
  timer_start(); {}
#define timer_mark(msg) \
  _secs[step] = get_seconds() - _secs[step]; \
  strcpy(_stepname[step], msg); \
  step++; \
  _secs[step] = get_seconds(); {}
#define timer_report(strm, msg, n) { \
  int ts, sl; \
  for (ts = 0 ; ts<step ; ts++) \
    MPI_Reduce(_secs+ts, _tsec+ts, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD); \
  if (MYNODE==0) { \
    for (ts=0 ; ts<step ; ts++) { \
      fprintf(strm,"(%3d)PE%3d: n: %12d (L%2d) %2d %s: ", \
	      NODES,MYNODE,n,_timer_level,ts,_stepname[ts]); \
      sl = strlen(_stepname[ts]); \
      while (sl++<TSTRLEN) fprintf(strm," "); \
      fprintf(strm," %12.6f\n",_tsec[ts]); \
    } \
    for (ts=0 ; ts<step ; ts++) { \
      fprintf(strm,"TIMER %3d %12d %2d %2d %12.6f\n", \
	      NODES,n,_timer_level,ts,_tsec[ts]); \
    } \
    for (ts=1 ; ts<step ; ts++) \
      _tsec[0] += _tsec[ts]; \
    fprintf(strm,"(%3d)PE%3d:   +   %s: ",NODES,MYNODE,msg); \
    sl = strlen(msg); \
    while (sl++<TSTRLEN) fprintf(strm," "); \
    fprintf(strm," %12.6f\n",_tsec[0]); \
    fflush(strm); \
  } \
  _timer_level--; \
  MPI_Barrier(MPI_COMM_WORLD); \
  }

/******************************************************************/

#define INIT_STEP() \
  int step = 0; \
  double secs[TIME_STEPS],tsec[TIME_STEPS]; \
  if (MYNODE==0) { \
    fprintf(outfile,"PE%3d: \n",MYNODE); \
    fflush(outfile); \
  } {}

#define START_STEP() \
  MPI_Barrier(MPI_COMM_WORLD); \
  secs[step] = get_seconds(); {}

#define END_STEP(a) \
  secs[step] = get_seconds() - secs[step]; \
  MPI_Reduce(secs+step, tsec+step, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD); \
  if (MYNODE==0) { \
    fprintf(outfile,"PE%3d:     %s: %9.6f\n",MYNODE,a,tsec[step]); \
    fflush(outfile); \
  } \
  step++; {}

#define REPORT_STEP(a) { \
  int ts; \
  if (MYNODE==0) { \
    for (ts=1 ; ts<step ; ts++) \
      tsec[0] += tsec[ts]; \
    fprintf(outfile,"(%3d)PE%3d:     %s: %9.6f\n",NODES,MYNODE,a,tsec[0]); \
    fflush(outfile); \
  } }


#endif
