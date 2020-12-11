/*
   $Id: misc.h,v 1.9 1999/09/25 16:37:00 dbader Exp $
   
   $Log: misc.h,v $
   Revision 1.9  1999/09/25 16:37:00  dbader
   Changed CYGWIN_NT-4.0 to CYGWIN_NT_40.

   Revision 1.8  1999/09/25 16:09:57  dbader
   Added support for CYGWIN_NT-4.0

   Revision 1.7  1999/07/23 16:12:22  dbader
   Modified to use CYGWIN32_NT

   Revision 1.6  1999/07/15 01:30:19  dbader
   Added MALLOC_DEBUG support.

   Revision 1.5  1999/07/03 23:00:39  dbader
   Converted assert_malloc to SAFE_MALLOC

   Revision 1.4  1999/06/02 01:28:24  dbader
   Added new assert_malloc()

   Revision 1.3  1999/05/27 19:25:13  dbader
   Daily Update 990527

   Revision 1.2  1999/05/27 19:20:54  dbader
   Daily Update 990525

   Revision 1.1  1999/05/27 19:09:08  dbader
   Initial revision


*/

#ifndef MISC_H
#define MISC_H

#ifdef DEBUG_MALLOC
#include "rmalloc.h"
#endif

#if defined(_CYGWIN32_NT)||defined(_CYGWIN_NT_40)
#define srandom(x) srand(x)
#define random() rand()
#endif

typedef struct {
  int tail;      /* Tail vertex label */
  int head;      /* Head vertex label */
  int tailAssn;  /* Tail vertex processor assignment */
  int headAssn;  /* Head vertex processor assignment */
} transArc_t;

void transArcCopy(transArc_t *, transArc_t *);

void *SAFE_MALLOC(int, char *);

void assert_malloc(void *, char *);

typedef int BOOL;

#define TRUE  1
#define FALSE 0

#define min(a,b) ((a)<(b)?(a):(b))
#define max(a,b) ((a)>(b)?(a):(b))

int lastB(int, int);
int setB(int, int);
int clearB(int, int);
int testB(int, int);
int clearLastB(int, int);

int log2_i(int);

#endif
