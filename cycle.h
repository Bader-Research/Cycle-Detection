/*
   $Id: cycle.h,v 1.7 1999/07/15 01:30:19 dbader Exp $
   
   $Log: cycle.h,v $
   Revision 1.7  1999/07/15 01:30:19  dbader
   Added MALLOC_DEBUG support.

   Revision 1.6  1999/07/05 01:44:06  dbader
   Split Express Graph LUT into entr and exit vertices, and
   removed entrLUT.

   Revision 1.5  1999/06/01 13:15:12  dbader
   Added new field to pExpGraph

   Revision 1.4  1999/05/29 15:52:28  dbader
   Added new auxiliary fields for the packed express graph.

   Revision 1.3  1999/05/27 19:22:55  dbader
   Daily Update 990526

   Revision 1.2  1999/05/27 19:20:54  dbader
   Daily Update 990525

   Revision 1.1  1999/05/27 19:09:08  dbader
   Initial revision


*/

#ifndef CYCLE_H
#define CYCLE_H

#ifdef DEBUG_MALLOC
#include "rmalloc.h"
#endif

#define BLACK 0
#define WHITE 1
#define GREEN 2
#define RED   3

typedef struct arc_s {
  int head;        /* vertex at head of arc */
  int assn;        /* node assignment of the head vertex */
} *arc_t;

typedef struct vertex_s {
  int label;       /* My vertex label */
  int arcs;        /* Number of adjacent vertices */
  arc_t alist;     /* Array of arcs */
} *vertex_t;
 
typedef struct reachlist_s {
  int idx;
  struct reachlist_s *next;
} *reachlist_t;

typedef struct vertexList_s {
  int num;              /* Number of vertices */
  vertex_t vlist;       /* Array of vertices */
  unsigned char *color; /* Color of each vertex */
  reachlist_t R;        /* list of outgoing transarc tail vertices */
} *vertexList_t;

typedef struct expVertex_s {
  transArc_t         *transArc;
  int                 expArcsNum;
  struct expArc_s    *expArcs;
  struct expVertex_s *next;
} *expVertex_t;

typedef struct expArc_s {
  expVertex_t headVertex;
  struct expArc_s *next;
} *expArc_t;

typedef struct labelLookup_s {
  int         label;
  expVertex_t expVertex;
} *labelLookup_t;

#define ENTR_LUT 0
typedef struct eGraph_s {
  int entrNum;           /* Number of entrance vertices */
  int exitNum;           /* Number of exit vertices */
  expVertex_t entrV;     /* Entrance vertices */
  expVertex_t exitV;     /* Exit vertices */
#if ENTR_LUT
  labelLookup_t entrLUT; /* Table from vertex label -> express vertex */
  int entrLUTNum;        /* Number of filled LUT entries */
#endif
  labelLookup_t exitLUT; /* Table from vertex label -> express vertex */
  int exitLUTNum;        /* Number of filled LUT entries */
} *eGraph_t;

/***************************************************************/
typedef struct pArc_s {
  int idx;
  int head;
  struct pArc_s *next;
} *pArc_t;

typedef struct pArcFlat_s {
  int idx;
  int head;
} *pArcFlat_t;

typedef struct pExpGraph_s {
  int entrNum;           /* Number of entrance vertices */
  int exitNum;           /* Number of exit vertices */
  int arcNum;            /* Number of express arcs */
  int *data;             /* Data of packed express graph */
  /********************************************************/
  struct pArc_s *newArcs; /* Additional express arcs to merge in */
  int deadEntr;
  int deadExit;
  int deadArc;
  int newArc;
  BOOL *deadEntrMask;
  BOOL *deadExitMask;
} *pExpGraph_t;

/***************************************************************/

#endif
