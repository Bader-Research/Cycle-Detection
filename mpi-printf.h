/*
   $Id: mpi-printf.h,v 1.2 1999/07/15 01:30:19 dbader Exp $
   
   $Log: mpi-printf.h,v $
   Revision 1.2  1999/07/15 01:30:19  dbader
   Added MALLOC_DEBUG support.

   Revision 1.1  1999/05/27 19:09:08  dbader
   Initial revision


*/

#ifndef _MPI_PRINTF_H
#define _MPI_PRINTF_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>

#ifdef DEBUG_MALLOC
#include "rmalloc.h"
#endif

void MPI_printf(char *fmt, ...);
void MPI_fprintf(FILE *fp, char *fmt, ...);

#endif
