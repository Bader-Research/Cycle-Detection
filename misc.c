/*

   $Log: misc.c,v $
   Revision 1.8  1999/07/03 23:00:39  dbader
   Converted assert_malloc to SAFE_MALLOC

   Revision 1.7  1999/06/04 12:21:26  dbader
   Changed log2 function to return ceiling rather than floor, thus
   handling non-power-of-two machine sizes

   Revision 1.6  1999/06/02 01:28:24  dbader
   Added new assert_malloc()

   Revision 1.5  1999/05/27 19:25:13  dbader
   Daily Update 990527

   Revision 1.4  1999/05/27 19:22:55  dbader
   Daily Update 990526

   Revision 1.3  1999/05/27 19:20:54  dbader
   Daily Update 990525

   Revision 1.2  1999/05/27 19:15:54  dbader
   Daily Update 990523

   Revision 1.1  1999/05/27 19:09:08  dbader
   Initial revision


*/

static char rcsid[] = "$Id: misc.c,v 1.8 1999/07/03 23:00:39 dbader Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "misc.h"

void assert_malloc(void *ptr, char *str) {
  if (ptr == NULL) {
    fprintf(stderr,"ERROR: NULL ptr: %s\n",str);
    exit(-1);
  }
  return;
}

void *SAFE_MALLOC(int len, char *msg) {
  void *ptr;
  ptr = malloc(len);
  assert_malloc(ptr, msg);
  return(ptr);
}

int lastB(int z, int h) {
  return(z & ((1<<h)-1));
}

int setB(int z, int h) {
  return (z | (1<<h));
}

int clearB(int z, int h) {
  return (z & ~(1<<h));
}

int testB(int z, int h) {
  return ((z & (1<<h)) > 0);
}

int clearLastB(int z, int h) {
  return (z & ~((1<<(h+1))-1));
}

void transArcCopy(transArc_t *dst, transArc_t *src) {
  /* Deep copy */
  dst->tail = src->tail;
  dst->head = src->head;
  dst->tailAssn = src->tailAssn;
  dst->headAssn = src->headAssn;
  return;
}

#define log2(a) (log((a))/log(2.0))

int log2_i(int num) {
  return((int)ceil(log2((double)num)));
}
