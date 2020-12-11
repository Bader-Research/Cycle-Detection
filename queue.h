/*
   $Id: queue.h,v 1.1 1999/05/27 19:09:08 dbader Exp $
   
   $Log: queue.h,v $
   Revision 1.1  1999/05/27 19:09:08  dbader
   Initial revision


*/

#ifndef QUEUE_H
#define QUEUE_H

#include "misc.h"

typedef transArc_t QELEM_T;

/* Defining the Q data structure used in the BFS */
typedef struct queue_tag{
   int count;
   int front;
   int rear;
   int MAXQ;
   QELEM_T *entry;
} *Q_t;

void Qinit(Q_t *, int);
void Qfree(Q_t);
void Qadd(Q_t, QELEM_T *);
void Qdel(Q_t, QELEM_T *);
int Qsize(Q_t);
BOOL Qempty(Q_t);
BOOL Qfull(Q_t);

void QelemCopy(QELEM_T *, QELEM_T *);

#endif
