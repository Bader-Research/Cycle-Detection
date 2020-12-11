/*
   $Id: cycleInterval.h,v 1.5 1999/07/15 01:30:19 dbader Exp $
   
   $Log: cycleInterval.h,v $
   Revision 1.5  1999/07/15 01:30:19  dbader
   Added MALLOC_DEBUG support.

   Revision 1.4  1999/07/14 09:55:01  dbader
   Removed PEG Format.
   Changed expGraph format not to record every express arc,
   rather, just the convex notation of the interval of reachable
   exit vertices.

   Revision 1.3  1999/07/14 02:01:28  dbader
   Modified the Discovery routine to use convex notation rather than reachlists.

   Revision 1.2  1999/07/12 21:27:09  dbader
   Working copy of the new Interval Graph representation.

   Revision 1.1  1999/07/12 15:50:04  dbader
   Initial revision

   Revision 1.4  1999/07/11 13:57:34  dbader
   Removed old Packed Convex Graph (PCG) representation code.

   Revision 1.3  1999/07/11 03:56:50  dbader
   Added Packed Interval Graph (PIG) representation, analog to
   the Packed Express Graph (PEG). Instead of storing express arcs,
   intervals are saved instead. Intervals can be used since the input
   contains convex bipartite graphs on each processor. This
   drastically reduces the packed graph size.

   Revision 1.2  1999/07/05 01:44:06  dbader
   Split Express Graph LUT into entr and exit vertices, and
   removed entrLUT.

   Revision 1.1  1999/07/03 23:03:44  dbader
   Initial revision

   Revision 1.1  1999/07/03 23:02:53  dbader
   Initial revision

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
  int C0;          /* minlabel of outgoing transarc tail vertices */
  int C1;          /* maxlabel of outgoing transarc tail vertices */
} *vertex_t;
 
typedef struct vertexList_s {
  int num;              /* Number of vertices */
  vertex_t vlist;       /* Array of vertices */
  unsigned char *color; /* Color of each vertex */
} *vertexList_t;

typedef struct expVertex_s {
  transArc_t         *transArc;
  int C0;
  int C1;
  struct expVertex_s *next;
} *expVertex_t;

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

typedef struct pIntervalGraph_s {
  int entrNum;           /* Number of entrance vertices */
  int exitNum;           /* Number of exit vertices */
  int intervalNum;       /* Number of intervals */
  int *data;             /* Data of packed express graph */
  /********************************************************/
  int deadEntr;
  int deadExit;
  BOOL *deadEntrMask;
  BOOL *deadExitMask;
  int deadInt;
} *pIntervalGraph_t;

/***************************************************************/

#endif
