/*

   $Log: queue.c,v $
   Revision 1.5  1999/07/03 23:00:39  dbader
   Converted assert_malloc to SAFE_MALLOC

   Revision 1.4  1999/06/02 01:28:24  dbader
   Added new assert_malloc()

   Revision 1.3  1999/06/02 01:10:19  dbader
   Changed a typecast and return void.

   Revision 1.2  1999/05/27 19:15:54  dbader
   Daily Update 990523

   Revision 1.1  1999/05/27 19:09:08  dbader
   Initial revision


*/

static char rcsid[] = "$Id: queue.c,v 1.5 1999/07/03 23:00:39 dbader Exp $";

#include <stdio.h>
#include <stdlib.h>
#include "queue.h"
#include "misc.h"

void QelemCopy(QELEM_T *dst, QELEM_T *src) {
  transArcCopy((transArc_t *)dst, (transArc_t *)src);
  return;
}
	       
/* Initializing the Queue */
void Qinit(Q_t *q_ptr, int size) {
   *q_ptr = (Q_t)SAFE_MALLOC(sizeof(struct queue_tag),
			     "(queue.c) *q_ptr");
   (*q_ptr)->count=0;
   (*q_ptr)->front=0;
   (*q_ptr)->rear=-1;
   (*q_ptr)->MAXQ = size;
   (*q_ptr)->entry = (QELEM_T *)SAFE_MALLOC((*q_ptr)->MAXQ * sizeof(QELEM_T),
					    "(queue.c) (*q_ptr)->entry");
   return;
}

void Qfree(Q_t q_ptr) {
   free(q_ptr->entry);
   free(q_ptr);
   return;
}

/* Adding nodes to the Queue */
void Qadd(Q_t q_ptr, QELEM_T *elem) {
   if(q_ptr->count >= q_ptr->MAXQ) {
      printf("count when Q is full:%d\n",q_ptr->count);
      printf("Queue is full\n");
      exit(1);
   }
   else {
      q_ptr->count++;
      q_ptr->rear=(q_ptr->rear+1)%(q_ptr->MAXQ);
      QelemCopy(q_ptr->entry + q_ptr->rear, elem);
   }
   return;
}

/* Deleting nodes already visited,from the Queue */
void Qdel(Q_t q_ptr, QELEM_T *elem) {
   if(q_ptr->count <= 0) {
      printf("Queue is empty\n");
      exit(1);
   }
   else{
      q_ptr->count--;
      QelemCopy(elem, q_ptr->entry + q_ptr->front);
      q_ptr->front=(q_ptr->front+1)%(q_ptr->MAXQ);
   }
   return;
}

/* Determine Queue Size */
int Qsize(Q_t q_ptr) {
   return q_ptr->count;
}

/* Checking if Queue is empty */
BOOL Qempty(Q_t q_ptr) {
   return(q_ptr->count <= 0);
}

/* Checking if Queue is full */
BOOL Qfull(Q_t q_ptr) {
   return(q_ptr->count >= q_ptr->MAXQ);
}
