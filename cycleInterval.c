/*

   $Log: cycleInterval.c,v $
   Revision 1.6  1999/07/15 01:00:31  dbader
   Increased the memory allocation for initial transarcs.

   Revision 1.5  1999/07/14 09:55:01  dbader
   Removed PEG Format.
   Changed expGraph format not to record every express arc,
   rather, just the convex notation of the interval of reachable
   exit vertices.

   Revision 1.4  1999/07/14 09:24:09  dbader
   Modified Find_Local_Cycles to record the initial transarcs,
   thus, removing the second local DFS from the Discovery routine.

   Revision 1.3  1999/07/14 02:01:28  dbader
   Modified the Discovery routine to use convex notation rather than reachlists.

   Revision 1.2  1999/07/12 21:27:09  dbader
   Working copy of the new Interval Graph representation.

   Revision 1.1  1999/07/12 15:50:04  dbader
   Initial revision

   Revision 1.6  1999/07/11 15:42:13  dbader
   Fixed cases where a NULL ptr is freed.

   Revision 1.5  1999/07/11 13:57:34  dbader
   Removed old Packed Convex Graph (PCG) representation code.
   Fixed a minor bug with PIG --> an uninitialized pointer that was being
   freed.

   Revision 1.4  1999/07/11 03:56:50  dbader
   Added Packed Interval Graph (PIG) representation, analog to
   the Packed Express Graph (PEG). Instead of storing express arcs,
   intervals are saved instead. Intervals can be used since the input
   contains convex bipartite graphs on each processor. This
   drastically reduces the packed graph size.

   Revision 1.3  1999/07/05 01:44:06  dbader
   Split Express Graph LUT into entr and exit vertices, and
   removed entrLUT.

   Revision 1.2  1999/07/04 20:24:40  dbader
   Modified the code so that a Lattice Input Graph could be used
   in the current framework. This involved modifying several functions
   that map between index and label, so that it is not assumed that
   an index has a linear mapping to the vertices' labels.
   Also, since several index to label (and vice versa) mapping functions
   exist that are input-dependent, a new style for choosing inputs is
   given in the code preamble.

   Revision 1.1  1999/07/03 23:03:44  dbader
   Initial revision

   Revision 1.1  1999/07/03 23:02:53  dbader
   Initial revision

   Revision 1.27  1999/07/03 23:00:39  dbader
   Converted assert_malloc to SAFE_MALLOC

   Revision 1.26  1999/07/03 22:42:35  dbader
   PRINT_MERGESIZE defined

   Revision 1.25  1999/06/23 21:39:52  dbader
   Changed order of loop in Free_Input so that Linux wouldn't hang.

   Revision 1.24  1999/06/08 19:01:18  dbader
   Added -DUSE_RAND option to select rand()/srand() instead of
   random()/srandom().

   Revision 1.23  1999/06/08 18:57:30  dbader
   Added printing of input graph description

   Revision 1.22  1999/06/07 16:51:49  dbader
   Changed packed express graph express-arc lists on entrance
   vertices to be kept in sorted order per each entrance vertex.
   Linear searches for express arcs have been changed to binary
   searches.

   Revision 1.21  1999/06/05 21:38:37  dbader
   Fixed a color[a->index] -> color[aidx] problem in Phase 1.
   Added a new graph type with local cycles.
   Fixed several warning conditions (unused or unset variables)

   Revision 1.20  1999/06/04 12:58:16  dbader
   Changed DEBUG_PRINT to define.

   Revision 1.19  1999/06/04 12:21:26  dbader
   Fixed memory leak in CleanUp function.

   Revision 1.18  1999/06/03 12:45:09  dbader
   Changed reporting to include problem size.

   Revision 1.17  1999/06/02 13:26:48  dbader
   Discontinued detection if a cycle is detected during the local phase.

   Revision 1.16  1999/06/02 03:25:42  dbader
   Added argc/argv tests.

   Revision 1.15  1999/06/02 03:21:24  dbader
   Miscellaneous minor changes to N_DEFAULT

   Revision 1.14  1999/06/02 01:28:24  dbader
   Added new assert_malloc()

   Revision 1.13  1999/06/01 21:17:18  dbader
   Fixed misplaced RemoveVertex call.

   Revision 1.12  1999/06/01 19:50:28  dbader
   Added timing features.

   Revision 1.11  1999/06/01 13:15:12  dbader
   Detected multiple exit vertices with the same label.
   WORKING implementation at this point.
   Merge Express Packed graphs finished.

   Revision 1.10  1999/05/29 15:52:28  dbader
   Implementing the merge of packed express graphs with new arcs and deleted
   vertices. When a vertex is deleted, it's marked with a negative label.
   Arcs are added to a linked list, and then the CleanUp routine flattens
   these changes back to a packed express graph

   Revision 1.9  1999/05/28 23:06:21  dbader
   Added more merge testing, printing out new express arcs.

   Revision 1.8  1999/05/28 17:56:41  dbader
   Added check for cycles during merge.
   Correctly identifies predessors of exit0 and
   successors of entr1 during merge.

   Revision 1.7  1999/05/28 14:00:21  dbader
   Determined an error: after an express merge, an entrance
   vertex can have adjacent transarcs from two different
   exit vertices from different processors.

   Revision 1.6  1999/05/27 23:17:10  dbader
   Daily Update 990527a

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

static char rcsid[] = "$Id: cycleInterval.c,v 1.6 1999/07/15 01:00:31 dbader Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
#include "misc.h"
#include "cycleInterval.h"
#include "queue.h"
#include "mpi-printf.h"
#include "timing.h"

#define N_DEFAULT  (1<<10)
#define MAXARCS     4
#define INPUTCYCLE  1

#define PRINT_MERGESIZE 0

/*********************************************************/
#define INPUT_LATTICE     0
#define INPUT_LINEAR      1
#define INPUT_LOCALCYCLES 2

#define INPUT_NUM INPUT_LATTICE

/**********/

void Create_Input_Lattice(int, vertexList_t *);
void Create_Input_Linear(int, vertexList_t *);
void Create_Input_LocalCycles(int, vertexList_t *);

void (*Create_Input_func[])() = {
  Create_Input_Lattice,
  Create_Input_Linear,
  Create_Input_LocalCycles
};

#define Create_Input(n, verts) ((Create_Input_func[INPUT_NUM])((n), (verts)))

/**********/

int localIndex_Lattice(vertexList_t, int);
int localIndex_Default(vertexList_t, int);

int (*localIndex_func[])() = {
  localIndex_Lattice,
  localIndex_Default, /* Linear */
  localIndex_Default /* LocalCycles */
};

#define localIndex(verts, n) ((localIndex_func[INPUT_NUM])((verts), (n)))

/*********************************************************/

BOOL CYCLE_FOUND;

int MYNODE, NODES;
FILE *outfile;
FILE *errfile;

char *INPUT_GRAPH_NAME;

#define MPI_TAG_SIZES   0
#define MPI_TAG_DATA    1

#define ENTRV           0
#define EXITV           1

int randmax(int m) {
  /* Return a random integer (0 <= r < m) */
#ifdef USE_RAND
  return(rand() % m);
#else
  return(random() % m);
#endif
}

int randrange(int a, int b) {
  /* Return a random integer (a <= r < b) */
  int range;
  range = b-a;
  return (randmax(range) + a);
}

#if 0
static int intCompare(const void *a, const void *b) {
  return (*(int *)a - *(int *)b);
}
#endif

static int transArcCompareHead(const void *a, const void *b) {
  return (((transArc_t *)a)->headAssn - ((transArc_t *)b)->headAssn);
}

static int expGraphLabelCompare(const void *a, const void *b) {
  return (((labelLookup_t)a)->label - ((labelLookup_t)b)->label);
}


int calcLatticeLabel(int offset, int vLocRow, int vLocCol, int nRow, int nCol) {
  int lab;
  
  if ((vLocRow == 0) || (vLocRow == nRow-1) ||
      (vLocCol == 0) || (vLocCol == nCol-1)) {
    /* Boundary vertex */
    if (vLocRow==0) {
      /* Top Row, left to right */
      lab = vLocCol;
    }
    else {
      if (vLocCol == nCol-1) {
	/* Right Column, top to bottom */
	lab = nCol + vLocRow - 1;
      }
      else {
	if (vLocRow == nRow-1) {
	  /* Bottom Row, right to left  */
	  lab = 2*nCol + nRow - vLocCol - 3;
	}
	else /* (vLocCol == 0) */ {
	  /* Left Column, bottom to top */
	  lab = 2*(nCol+nRow) - vLocRow - 4;
	}
      }
    }
  }
  else {
    /* Interior vertex */
    lab  = 2*(nRow + nCol) - 4;
    lab += (vLocRow-1)*(nCol-2) + (vLocCol-1);
  }
  return (offset+lab);
}

void Create_Input_Lattice(int n, vertexList_t *verts) {
  int num;
  int offset;
  int i;
  vertex_t v;
  arc_t a;
  int pRow, pCol, myRow, myCol, nRow, nCol, vLoc, vLocRow, vLocCol;
  
  MPI_Barrier(MPI_COMM_WORLD);

#if INPUTCYCLE
  INPUT_GRAPH_NAME = "Lattice with cycle";
#else
  INPUT_GRAPH_NAME = "Lattice with no cycle";
#endif

  *verts = (vertexList_t)SAFE_MALLOC(sizeof(struct vertexList_s),
				     "(cycle.c) *verts");

  num = (*verts)->num = (int)floor((double)n/(double)NODES);
  (*verts)->vlist = (vertex_t)SAFE_MALLOC(num*sizeof(struct vertex_s),
					  "(cycle.c) (*verts)->vlist");

  pRow = (int)floor(sqrt((double)NODES));
  pCol = (int)floor(sqrt((double)NODES));
  if (pRow*pCol != NODES) {
    fprintf(errfile,"ERROR: Must use a square number of procs\n");
    /* exit(-1); */
  }

  myRow = MYNODE / pCol;
  myCol = MYNODE % pCol;

  nRow  = (int)floor(sqrt((double)num));
  nCol  = (int)floor(sqrt((double)num));
  if (nRow*nCol != num) {
    fprintf(errfile,"ERROR: Must use a square for vertices per proc\n");
    /* exit(-1); */
  }
    
  offset = num*MYNODE;
  for (i=0 ; i<num ; i++) {
    v = (*verts)->vlist + i;

    vLoc    = i;
    vLocRow = vLoc / nCol;
    vLocCol = vLoc % nCol;

    v->label = calcLatticeLabel(offset, vLocRow, vLocCol, nRow, nCol);
    v->C0 = -1;
    v->C1 = -1;

    v->arcs  = 0;
    if (vLocCol < nCol - 1)
      v->arcs++;
    if (vLocRow < nRow - 1)
      v->arcs++;
    if ((vLocCol == nCol-1) && (myCol < pCol-1))
      v->arcs++;
    if ((vLocRow == nRow-1) && (myRow < pRow-1))
      v->arcs++;
#if INPUTCYCLE
    /* Add arc from last vertex to first vertex */
#if 0
    if (v->label == (num*NODES)-1) {
      v->arcs++;
    }
#else
    if ((MYNODE == NODES-1) && (vLocRow == nRow-1) && (vLocCol == nCol-1)) {
      v->arcs++;
    }
#endif
#endif
    if (v->arcs > 0) 
      v->alist = (arc_t)SAFE_MALLOC((v->arcs)*sizeof(struct arc_s),
				    "(cycle.c) v->alist");
    else
      v->alist = (arc_t)NULL;
    
    a = v->alist;
    if (vLocCol < nCol - 1) {
      a->head = calcLatticeLabel(offset, vLocRow, vLocCol+1, nRow, nCol);
      a->assn = MYNODE;
      a++;
    }
    if (vLocRow < nRow - 1) {
      a->head = calcLatticeLabel(offset, vLocRow+1, vLocCol, nRow, nCol);
      a->assn = MYNODE;
      a++;
    }
    if ((vLocCol == nCol-1) && (myCol < pCol-1)) {
      a->head = calcLatticeLabel(offset+num, vLocRow, 0, nRow, nCol);
      a->assn = MYNODE + 1;
      a++;
    }
    if ((vLocRow == nRow-1) && (myRow < pRow-1)) {
      a->head = calcLatticeLabel(offset+(pCol*num), 0, vLocCol, nRow, nCol);
      a->assn = MYNODE + pCol;
      a++;
    }
#if INPUTCYCLE
    /* Add arc from last vertex to first vertex */
#if 0
    if (v->label == (num*NODES)-1)
#else
    if ((MYNODE == NODES-1) && (vLocRow == nRow-1) && (vLocCol == nCol-1)) 
#endif
      {
      a->head = 0;
      a->assn = 0;
      a++;
    }
#endif
  }
  
  (*verts)->color = (unsigned char *)SAFE_MALLOC(num*sizeof(unsigned char),
						 "(cycle.c) (*verts)->color");

  return;
}

void Create_Input_Linear(int n, vertexList_t *verts) {
  int num;
  int offset;
  int i;
  vertex_t v;
  arc_t a;
  
  MPI_Barrier(MPI_COMM_WORLD);

#if INPUTCYCLE
  INPUT_GRAPH_NAME = "Linear with cycle";
#else
  INPUT_GRAPH_NAME = "Linear with no cycle";
#endif

  *verts = (vertexList_t)SAFE_MALLOC(sizeof(struct vertexList_s),
				     "(cycle.c) *verts");

  num = (*verts)->num = (int)floor((double)n/(double)NODES);
  (*verts)->vlist = (vertex_t)SAFE_MALLOC(num*sizeof(struct vertex_s),
					  "(cycle.c) (*verts)->vlist");

  if (num*NODES != n) {
    fprintf(errfile,"ERROR: n does not divide p evenly\n");
    /* exit(-1); */
  }
    
  offset = num*MYNODE;
  for (i=0 ; i<num ; i++) {
    v = (*verts)->vlist + i;
    v->label = offset + i;
    v->C0    = -1;
    v->C1    = -1;
    v->arcs  = 1;
    if (v->label == n-1) {
#if INPUTCYCLE
      v->arcs = 1;
#else
      v->arcs = 0;
#endif
    }
    
    v->alist = (arc_t)SAFE_MALLOC((v->arcs)*sizeof(struct arc_s),
				  "(cycle.c) v->alist");
    a = v->alist;
    if (v->label < n-1) {
      a->head = v->label + 1;
      a->assn = MYNODE;
      if (i==num-1) 
	a->assn++;
    }
    else {
#if INPUTCYCLE
      a->head = 0;
      a->assn = 0;
#endif
    }
  }
  
  (*verts)->color = (unsigned char *)SAFE_MALLOC(num*sizeof(unsigned char),
						 "(cycle.c) (*verts)->color");

  return;
}

void Create_Input_LocalCycles(int n, vertexList_t *verts) {
  int num;
  int offset;
  int i;
  vertex_t v;
  arc_t a;
  
  MPI_Barrier(MPI_COMM_WORLD);

  INPUT_GRAPH_NAME = "LocalCycles";

  *verts = (vertexList_t)SAFE_MALLOC(sizeof(struct vertexList_s),
				     "(cycle.c) *verts");

  num = (*verts)->num = (int)floor((double)n/(double)NODES);
  (*verts)->vlist = (vertex_t)SAFE_MALLOC(num*sizeof(struct vertex_s),
					  "(cycle.c) (*verts)->vlist");

  if (num*NODES != n) {
    fprintf(errfile,"ERROR: n does not divide p evenly\n");
    /* exit(-1); */
  }
    
  offset = num*MYNODE;
  for (i=0 ; i<num ; i++) {
    v = (*verts)->vlist + i;
    v->label = offset + i;
    v->C0    = -1;
    v->C1    = -1;
    v->arcs  = 1;
    
    v->alist = (arc_t)SAFE_MALLOC((v->arcs)*sizeof(struct arc_s),
				  "(cycle.c) v->alist");
    a = v->alist;
    if (i < num-1) {
      a->head = v->label + 1;
      a->assn = MYNODE;
    }
    else {
      a->head = offset;
      a->assn = MYNODE;
    }
  }
  
  (*verts)->color = (unsigned char *)SAFE_MALLOC(num*sizeof(unsigned char),
						 "(cycle.c) (*verts)->color");

  return;
}

void Print_Input(vertexList_t verts) {
  int i, k, p;
  vertex_t v;
  arc_t a;
  
  MPI_Barrier(MPI_COMM_WORLD);

  for (p=0 ; p<NODES ; p++) {
    if (p==MYNODE) {

      fprintf(outfile,"PE%3d: vertices: %12d\n",MYNODE, verts->num);

      for (i=0 ; i<verts->num ; i++) {
	v = verts->vlist + i;
	fprintf(outfile,"PE%3d: vert[%12d]: label %12d arcs %2d\n",
		MYNODE, i,v->label,v->arcs);
	for (k=0 ; k<v->arcs ; k++) {
	  a = v->alist + k;
	  fprintf(outfile,"PE%3d: vert[%12d]:    arc[%2d]: (%12d, %3d)\n",
		  MYNODE, i, k,
		  a->head, a->assn);
	}
      }
      fflush(outfile);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  
  return;
}

int localIndex_Default(vertexList_t verts, int label) {
  int idx;
  /* Map the local vertex label to index */
  
  if ((label / verts->num) != MYNODE)
    fprintf(errfile,
	 "PE%3d: ERROR: localIndex of label off processor! l:%d n:%d d:%d\n",
	 MYNODE, label, verts->num, label/verts->num);

  idx = label % verts->num;
  if (verts->vlist[idx].label != label)
    fprintf(errfile, "PE%3d: ERROR: localIndex of label %d incorrect.\n",
	    MYNODE, label);

  return(idx);
}

int localIndex_Lattice(vertexList_t verts, int label) {
  int idx;
  int nRow, nCol;
  int vLocRow, vLocCol;
  int borderCnt, locLabel;
  /* Map the local vertex label to index */
  
  if ((label / verts->num) != MYNODE)
    fprintf(errfile,
	 "PE%3d: ERROR: localIndex of label off processor! l:%d n:%d d:%d\n",
	 MYNODE, label, verts->num, label/verts->num);

  nRow  = (int)floor(sqrt((double)verts->num));
  nCol  = (int)floor(sqrt((double)verts->num));
  locLabel = label - (verts->num * MYNODE);
  borderCnt = 2*(nRow + nCol) - 4;
  if (locLabel < borderCnt) {
    /* Boundary */
    if (locLabel < nCol) {
      /* Top Row */
      vLocRow = 0;
      vLocCol = locLabel;
    }
    else {
      if (locLabel < nCol + nRow - 1) {
	/* Right Column */
	vLocRow = locLabel - nCol + 1;
	vLocCol = nCol - 1;
      }
      else {
	if (locLabel < 2*nCol + nRow - 2) {
	  /* Bottom Row */
	  vLocRow = nRow - 1;
	  vLocCol = 2*nCol + nRow -locLabel - 3;
	}
	else {
	  /* Left Column */
	  vLocRow = 2*(nCol+nRow) - locLabel - 4;
	  vLocCol = 0;
	}
      }
    }
  }
  else {
    /* Interior */
    if (nCol >= 2) {
      locLabel -= borderCnt;
      vLocRow = locLabel / (nCol - 2);
      vLocRow++;
      vLocCol = locLabel % (nCol - 2);
      vLocCol++;
    }
    else {
      /* Single vertex */
      vLocRow =0;
      vLocCol =0;
    }
  }

  idx = vLocRow * nCol + vLocCol;

  if ((idx < 0) || (idx >= verts->num)) 
    fprintf(errfile, "PE%3d: ERROR: localIndex_Lattice  %d incorrect. idx OOB\n",
	    MYNODE, label);
    
  if (verts->vlist[idx].label != label)
    fprintf(errfile, "PE%3d: ERROR: localIndex of label %d incorrect.\n",
	    MYNODE, label);

  return(idx);
}

int globalLabel(vertexList_t verts, int idx) {
  /* Map the local index into the vertex label */
  int label;
  
  if ((idx < 0) || (idx >= verts->num))
    fprintf(errfile,"PE%3d: ERROR: globalLabel of idx bad!\n",MYNODE);

#if 0
  label = idx + (verts->num * MYNODE);
  if (label != verts->vlist[idx].label)
    fprintf(errfile,"PE%3d: ERROR: globalLabel (%d) of idx (%d) doesn't match!\n",
	    MYNODE, label, idx);
#else
  label = verts->vlist[idx].label;
#endif
  
  return(label);
}

void visitLocal(vertexList_t verts, int vidx,
		int *transArcNum, transArc_t *initTransArcs) {
  int i;
  vertex_t v;
  arc_t a;
  int aidx;
  transArc_t *tarc;
  
  verts->color[vidx] = RED;
  
  v = verts->vlist + vidx;
  for (i=0 ; i<v->arcs ; i++) {
    a = v->alist + i;
    if (a->assn == MYNODE) {  /* local-arc */
#if 1
      /* CHECK */
      aidx = localIndex(verts, a->head);
      if ((aidx < 0) || (aidx > verts->num)) {
	fprintf(errfile,"PE%3d: ERROR: aidx (%d) > verts->num (%d)\n",
		MYNODE, aidx, verts->num);
      }
#endif
      switch (verts->color[aidx]) {
      case WHITE:
	visitLocal(verts, aidx, transArcNum, initTransArcs);
	break;
      case RED:
	fprintf(outfile,"PE%3d: A local cycle (%d, %d) has been found!\n",
		MYNODE, v->label, a->head);
	CYCLE_FOUND = TRUE;
	break;
      case BLACK:
      default: ;
      }
#if 1
      if (verts->vlist[vidx].C0 < 0) {
	verts->vlist[vidx].C0 = verts->vlist[aidx].C0;
	verts->vlist[vidx].C1 = verts->vlist[aidx].C1;
      }
      else {
	verts->vlist[vidx].C0 = min(verts->vlist[vidx].C0, verts->vlist[aidx].C0);
	verts->vlist[vidx].C1 = max(verts->vlist[vidx].C1, verts->vlist[aidx].C1);
      }
#endif
    }
    else {
      tarc = initTransArcs + *transArcNum;
      tarc->tail = v->label;
      tarc->head = a->head;
      tarc->tailAssn = MYNODE;
      tarc->headAssn = a->assn;
      (*transArcNum)++;

      if (verts->vlist[vidx].C0 < 0) {
	verts->vlist[vidx].C0 = v->label;
	verts->vlist[vidx].C1 = v->label;
      }
      else {
	verts->vlist[vidx].C0 = min(verts->vlist[vidx].C0, v->label);
	verts->vlist[vidx].C1 = max(verts->vlist[vidx].C1, v->label);
      }
    }
  }

  verts->color[vidx] = BLACK;
  
  return;
}

void Find_Local_Cycles(vertexList_t verts, int *transArcNum,
		       transArc_t **transArcs) {
  int i, num;

  num = verts->num;
  *transArcNum = 0;
  *transArcs = (transArc_t *)SAFE_MALLOC(MAXARCS * num * sizeof(transArc_t),
					 "(cycle.c) transArcs");
  
  for (i=0 ; i<num ; i++)
    verts->color[i] = WHITE;

  for (i=0 ; i<num ; i++) {
    if ((verts->color[i] == WHITE) && (!CYCLE_FOUND)) {
      /* vertex i has not yet been visited */
      visitLocal(verts, i, transArcNum, *transArcs);
    }
  }

  return;
}


void Discovery(vertexList_t verts,
	       int initTransArcCount,
	       transArc_t *initTransArcs,
	       int *termTransArcCount,
	       transArc_t **termTransArcs)
{
  int i;
  int *send_cnt, *send_offset;
  int *recv_cnt, *recv_offset;
  MPI_Datatype MPI_TRANSARC;
  
#ifdef DEBUG_PRINT
  MPI_fprintf(outfile,"I have %d initial transarc%s\n",
	      initTransArcCount, (initTransArcCount==1)?"":"s");
#endif
  
  if (initTransArcs > 0) 
    qsort(initTransArcs, initTransArcCount,
	  sizeof(transArc_t), transArcCompareHead);

#ifdef DEBUG_PRINT
  {
    int k;
    MPI_Barrier(MPI_COMM_WORLD);
    for (k=0 ; k<NODES ; k++) {
      if (MYNODE==k) {
	for (i=0 ; i<initTransArcCount ; i++) {
	  fprintf(outfile,"PE%3d:  Send init arc[%2d]: (%3d, %3d, %3d, %3d)\n",
		  MYNODE, i,
		  initTransArcs[i].tail,
		  initTransArcs[i].head,
		  initTransArcs[i].tailAssn,
		  initTransArcs[i].headAssn);
	}
	fflush(outfile);
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }
#endif

  send_cnt = (int *)SAFE_MALLOC(NODES*sizeof(int),
				"(cycle.c) send_cnt");

  send_offset = (int *)SAFE_MALLOC(NODES*sizeof(int),
				   "(cycle.c) send_offset");

  recv_cnt = (int *)SAFE_MALLOC(NODES*sizeof(int),
				"(cycle.c) recv_cnt");

  recv_offset = (int *)SAFE_MALLOC(NODES*sizeof(int),
				   "(cycle.c) recv_offset");

  for (i=0 ; i<NODES ; i++)
    send_cnt[i] = 0;
  
  for (i=0 ; i<initTransArcCount ; i++)
    send_cnt[initTransArcs[i].headAssn]++;

  send_offset[0] = 0;
  for (i=1 ; i<NODES ; i++)
    send_offset[i] = send_offset[i-1] + send_cnt[i-1];


  MPI_Alltoall(send_cnt, 1, MPI_INT,
	       recv_cnt, 1, MPI_INT,
	       MPI_COMM_WORLD);

  recv_offset[0] = 0;
  for (i=1 ; i<NODES ; i++)
    recv_offset[i] = recv_offset[i-1] + recv_cnt[i-1];

  *termTransArcCount = recv_offset[NODES-1] + recv_cnt[NODES-1];

  *termTransArcs = (transArc_t *)SAFE_MALLOC(*termTransArcCount * sizeof(transArc_t),
					     "(cycle.c) *termTransArcs");

  MPI_Type_contiguous(sizeof(transArc_t)/sizeof(int),MPI_INT,&MPI_TRANSARC);
  MPI_Type_commit(&MPI_TRANSARC);

#ifdef DEBUG_PRINT
  {
    int k;
    MPI_Barrier(MPI_COMM_WORLD);
    for (k=0 ; k<NODES ; k++) {
      if (MYNODE==k) {
	for (i=0 ; i<NODES ; i++) {
	  fprintf(outfile,"PE%3d: [%d] send_cnt: %d recv_cnt: %d  send_offset: %d recv_offset: %d\n",
		  MYNODE, i,
		  send_cnt[i], recv_cnt[i],
		  send_offset[i], recv_offset[i]);
	}
	fflush(outfile);
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }
#endif

  MPI_Alltoallv(initTransArcs, send_cnt, send_offset, MPI_TRANSARC,
		*termTransArcs, recv_cnt, recv_offset, MPI_TRANSARC,
		MPI_COMM_WORLD);
  
  MPI_Type_free(&MPI_TRANSARC);

#ifdef DEBUG_PRINT
  {
    int k;
    MPI_Barrier(MPI_COMM_WORLD);
    for (k=0 ; k<NODES ; k++) {
      if (MYNODE==k) {
	for (i=0 ; i<*termTransArcCount ; i++) {
	  fprintf(outfile,"PE%3d: Rec'd term arc[%2d]: (%3d, %3d, %3d, %3d)\n",
		  MYNODE, i,
		  (*termTransArcs)[i].tail,
		  (*termTransArcs)[i].head,
		  (*termTransArcs)[i].tailAssn,
		  (*termTransArcs)[i].headAssn);
	}
	fflush(outfile);
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }
#endif

  MPI_Barrier(MPI_COMM_WORLD);

  free(recv_offset);
  free(recv_cnt);
  free(send_offset);
  free(send_cnt);
  return;
}

void AddEntrVertex(eGraph_t expGraph, transArc_t *tarc) {
  expVertex_t entrVertex;
  
  entrVertex = (expVertex_t)SAFE_MALLOC(sizeof(struct expVertex_s),
					"(cycle.c) entrVertex");

  entrVertex->transArc   = tarc;
  entrVertex->next       = expGraph->entrV;

  expGraph->entrV        = entrVertex;
  expGraph->entrNum++;

#if ENTR_LUT
  (expGraph->entrLUT + expGraph->entrLUTNum)->label = entrVertex->transArc->head;
  (expGraph->entrLUT + expGraph->entrLUTNum)->expVertex = entrVertex;
  expGraph->entrLUTNum++;
#endif
  
  return;
}

void AddExitVertex(eGraph_t expGraph, transArc_t *tarc) {
  expVertex_t exitVertex;
  
  exitVertex = (expVertex_t)SAFE_MALLOC(sizeof(struct expVertex_s),
					"(cycle.c) exitVertex");

  exitVertex->transArc   = tarc;
  exitVertex->next       = expGraph->exitV;

  expGraph->exitV        = exitVertex;
  expGraph->exitNum++;
  
  (expGraph->exitLUT + expGraph->exitLUTNum)->label = exitVertex->transArc->tail;
  (expGraph->exitLUT + expGraph->exitLUTNum)->expVertex = exitVertex;
  expGraph->exitLUTNum++;
  
  return;
}

int binSearchLUT(labelLookup_t lut, int minIdx, int maxIdx, int target) {
  /* Search searchList recursively using binary search to find target.
     Return -1 if it does not exist, or target's index if it does. */
  register int midIdx;

  midIdx = minIdx + (maxIdx - minIdx) / 2;
  if (maxIdx < minIdx)
    return(-1);
  else if (target < lut[midIdx].label)
    return binSearchLUT(lut, minIdx,     midIdx - 1, target);
  else if (target > lut[midIdx].label)
    return binSearchLUT(lut, midIdx + 1, maxIdx    , target);
  else
    return midIdx;
}

#if ENTR_LUT
expVertex_t getEntrVertexLUT(eGraph_t expGraph, int label) {
  int idx;
  expVertex_t v;

  idx = binSearchLUT(expGraph->entrLUT, 0, expGraph->entrLUTNum - 1, label);
  if (idx>=0)
    v = expGraph->entrLUT[idx].expVertex;
  else
    v = (expVertex_t)NULL;
  
  return v;
}
#endif

expVertex_t getExitVertexLUT(eGraph_t expGraph, int label) {
  int idx;
  expVertex_t v;

  idx = binSearchLUT(expGraph->exitLUT, 0, expGraph->exitLUTNum - 1, label);
  if (idx>=0)
    v = expGraph->exitLUT[idx].expVertex;
  else
    v = (expVertex_t)NULL;
  
  return v;
}

void AddExpressArcSet(vertexList_t verts, eGraph_t expGraph,
		      expVertex_t entrVertex) {
  int idx;
  int vlabel;
  
  vlabel = entrVertex->transArc->head;

  idx = localIndex(verts, vlabel);

  entrVertex->C0 = verts->vlist[idx].C0;
  entrVertex->C1 = verts->vlist[idx].C1;

  return;
}

void SortExpGraphLUT(eGraph_t expGraph) {
#if ENTR_LUT
  qsort(expGraph->entrLUT, expGraph->entrLUTNum, sizeof(struct labelLookup_s),
	expGraphLabelCompare);
#endif
  qsort(expGraph->exitLUT, expGraph->exitLUTNum, sizeof(struct labelLookup_s),
	expGraphLabelCompare);
  return;
}

void Create_Express(vertexList_t verts,
		    int initTransArcCount,
		    transArc_t *initTransArcs,
		    int termTransArcCount,
		    transArc_t *termTransArcs,
		    eGraph_t *expGraph)
{
  int i;
  expVertex_t entrVertex;
  int expVertNum;

  *expGraph = (eGraph_t)SAFE_MALLOC(sizeof(struct eGraph_s),
				    "(cycle.c) *expGraph");

  (*expGraph)->entrNum = 0;
  (*expGraph)->exitNum = 0;
  (*expGraph)->entrV   = (expVertex_t)NULL;
  (*expGraph)->exitV   = (expVertex_t)NULL;

  expVertNum           = initTransArcCount + termTransArcCount;
  
#if ENTR_LUT
  (*expGraph)->entrLUT = (labelLookup_t)
    SAFE_MALLOC(termTransArcCount * sizeof(struct labelLookup_s),
		"(cycle.c) (*expGraph)->entrLUT");
  (*expGraph)->entrLUTNum  = 0;
#endif
  
  (*expGraph)->exitLUT = (labelLookup_t)
    SAFE_MALLOC(initTransArcCount * sizeof(struct labelLookup_s),
		"(cycle.c) (*expGraph)->exitLUT");
  (*expGraph)->exitLUTNum  = 0;

  for (i=0 ; i<initTransArcCount ; i++) {
    AddExitVertex(*expGraph, initTransArcs+i);
  }

  if ((*expGraph)->exitLUTNum != initTransArcCount)
    fprintf(errfile,"PE%3d: ERROR: exitLUTNum (%d) != initTransArcCount (%d)\n",
	    MYNODE, (*expGraph)->exitLUTNum, initTransArcCount);

  for (i=0 ; i<termTransArcCount ; i++) {
    AddEntrVertex(*expGraph, termTransArcs+i);
  }

#if ENTR_LUT
  if ((*expGraph)->entrLUTNum != termTransArcCount)
    fprintf(errfile,"PE%3d: ERROR: entrLUTNum (%d) != termTransArcCount (%d)\n",
	    MYNODE, (*expGraph)->entrLUTNum, termTransArcCount);
#endif
  
  SortExpGraphLUT(*expGraph);
  
  entrVertex = (*expGraph)->entrV;
  while (entrVertex != (expVertex_t)NULL) {
    AddExpressArcSet(verts, *expGraph, entrVertex);
    entrVertex = entrVertex->next;
  }
    
  return;
}

void Free_Input(vertexList_t verts) {
  int i;

  free(verts->color);
  
  for (i=verts->num-1 ; i>=0 ; i--)
    if ((verts->vlist + i)->alist != (arc_t)NULL)
      free((verts->vlist + i)->alist);

  free(verts->vlist);
  free(verts);
  
  return;
}

void Free_expVertex(expVertex_t v) {
  free(v);
  return;
}

void Free_Express(eGraph_t expGraph) {
  expVertex_t v, nextv;

  free(expGraph->exitLUT);
#if ENTR_LUT
  free(expGraph->entrLUT);
#endif
  
  v = expGraph->entrV;
  while (v != (expVertex_t)NULL) {
    nextv = v->next;
    Free_expVertex(v);
    v = nextv;
  }
  
  v = expGraph->exitV;
  while (v != (expVertex_t)NULL) {
    nextv = v->next;
    Free_expVertex(v);
    v = nextv;
  }

  free(expGraph);
  return;
}

int getArcNum(eGraph_t expGraph) {
  int num;
  expVertex_t v;

  num = 0;
  v = expGraph->entrV;
  while (v != (expVertex_t)NULL) {
    num += (v->C1 - v->C0 + 1);
    v = v->next;
  }

  return num;
}

void Print_Express(eGraph_t expGraph) {
  int i, k;
  expVertex_t v;
  
  MPI_fprintf(outfile,"expGraph entr: %12d exit: %12d exarcs: %12d\n",
	      expGraph->entrNum, expGraph->exitNum, getArcNum(expGraph));
 
  for (i=0 ; i<NODES ; i++) {
    if (MYNODE==i) {
      fprintf(outfile,"PE%3d: Entrance vertices: %12d\n",
	      MYNODE, expGraph->entrNum);
      fflush(outfile);
      v = expGraph->entrV;
      k=0;
      while (v != (expVertex_t)NULL) {
	fprintf(outfile,"PE%3d: \t[%12d]: (%12d {%3d}, %12d {%3d})\n",
		MYNODE, k,
		v->transArc->tail, v->transArc->tailAssn,
		v->transArc->head, v->transArc->headAssn);
	fprintf(outfile,"PE%3d: \t\texpArcs: [%12d %12d]\n",
		MYNODE, v->C0, v->C1);
	if (v->transArc->headAssn != MYNODE)
	  fprintf(errfile,"PE%3d: ERROR: headAssn != MYNODE (%d, %d)\n", MYNODE,
		  v->transArc->headAssn, MYNODE);
	k++;
	v = v->next;
      }

      fprintf(outfile,"PE%3d:     Exit vertices: %12d\n",
	      MYNODE, expGraph->exitNum);
      fflush(outfile);
      v = expGraph->exitV;
      k=0;
      while (v != (expVertex_t)NULL) {
	fprintf(outfile,"PE%3d: \t[%12d]: (%12d {%3d}, %12d {%3d})\n",
		MYNODE, k,
		v->transArc->tail, v->transArc->tailAssn,
		v->transArc->head, v->transArc->headAssn);
	if (v->transArc->tailAssn != MYNODE)
	  fprintf(errfile,"PE%3d: ERROR: tailAssn != MYNODE (%d, %d)\n", MYNODE,
		  v->transArc->tailAssn, MYNODE);
	k++;
	v = v->next;
      }
      fflush(outfile);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  return;
}

/*************************************************************************/
   
int pig_GetVSize() {
  /* label, adj, adjassn, intidx */
  return 4;
}

int pig_GetISize() {
  /* C0, C1, nextIdx */
  return 3;
}

int pig_GetDataInInts_Var(int entrNum, int exitNum, int intervalNum) {
  return (pig_GetVSize()*(entrNum+exitNum) + pig_GetISize()*intervalNum);
}

int pig_GetDataInInts(pIntervalGraph_t pig) {
  return (pig_GetDataInInts_Var(pig->entrNum, pig->exitNum, pig->intervalNum));
}

int pig_GetDataInBytes(pIntervalGraph_t pig) {
  return (pig_GetDataInInts(pig)*sizeof(int));
}

int pig_GetLabel(pIntervalGraph_t pig, int idx, BOOL vType) {
  int x;
  switch (vType) {
  case ENTRV:
    x = *( pig->data + (idx*pig_GetVSize()) );
    break;
  case EXITV:
    x = *( pig->data + ((pig->entrNum + idx)*pig_GetVSize()) );
    break;
  default:
    x = -1;
    fprintf(errfile,"ERROR: pig_GetLabel()\n");
  }
  return x;
}

void pig_SetLabel(pIntervalGraph_t pig, int idx, BOOL vType, int val) {
  switch (vType) {
  case ENTRV:
    *( pig->data + (idx*pig_GetVSize()) ) = val;
    break;
  case EXITV:
    *( pig->data + ((pig->entrNum + idx)*pig_GetVSize()) ) = val;
    break;
  default:
    fprintf(errfile,"ERROR: pig_SetLabel()\n");
  }
  return;
}

int pig_GetAdj(pIntervalGraph_t pig, int idx, BOOL vType) {
  int x;
  switch (vType) {
  case ENTRV:
    x = *( pig->data + (idx*pig_GetVSize()) + 1);
    break;
  case EXITV:
    x = *( pig->data + ((pig->entrNum + idx)*pig_GetVSize()) + 1);
    break;
  default:
    x = -1;
    fprintf(errfile,"ERROR: pig_GetAdj()\n");
  }
  return x;
}

void pig_SetAdj(pIntervalGraph_t pig, int idx, BOOL vType, int val) {
  switch (vType) {
  case ENTRV:
    *( pig->data + (idx*pig_GetVSize()) + 1) = val;
    break;
  case EXITV:
    *( pig->data + ((pig->entrNum + idx)*pig_GetVSize()) + 1) = val;
    break;
  default:
    fprintf(errfile,"ERROR: pig_SetAdj()\n");
  }
  return;
}

int pig_GetAdjAssn(pIntervalGraph_t pig, int idx, BOOL vType) {
  int x;
  switch (vType) {
  case ENTRV:
    x = *( pig->data + (idx*pig_GetVSize()) + 2);
    break;
  case EXITV:
    x = *( pig->data + ((pig->entrNum + idx)*pig_GetVSize()) + 2);
    break;
  default:
    x = -1;
    fprintf(errfile,"ERROR: pig_GetAdjAssn()\n");
  }
  return x;
}

void pig_SetAdjAssn(pIntervalGraph_t pig, int idx, BOOL vType, int val) {
  switch (vType) {
  case ENTRV:
    *( pig->data + (idx*pig_GetVSize()) + 2) = val;
    break;
  case EXITV:
    *( pig->data + ((pig->entrNum + idx)*pig_GetVSize()) + 2) = val;
    break;
  default:
    fprintf(errfile,"ERROR: pig_SetAdjAssn()\n");
  }
  return;
}

int pig_GetIntervalIndex(pIntervalGraph_t pig, int idx, BOOL vType) {
  int x;
  switch (vType) {
  case ENTRV:
    x = *( pig->data + (idx*pig_GetVSize()) + 3);
    break;
  case EXITV:
    x = *( pig->data + ((pig->entrNum + idx)*pig_GetVSize()) + 3);
    if (x != -1)
      fprintf(errfile,"ERROR: pig_GetIntervalIndex() EXITV.\n");
    break;
  default:
    x = -1;
    fprintf(errfile,"ERROR: pig_GetIntervalIndex()\n");
  }
  return x;
}

void pig_SetIntervalIndex(pIntervalGraph_t pig, int idx, BOOL vType, int val) {
  switch (vType) {
  case ENTRV:
    *( pig->data + (idx*pig_GetVSize()) + 3) = val;
    break;
  case EXITV:
    *( pig->data + ((pig->entrNum + idx)*pig_GetVSize()) + 3) = val;
    if (val != -1)
      fprintf(errfile,"ERROR: pig_SetIntervalIndex() EXITV.\n");
    break;
  default:
    fprintf(errfile,"ERROR: pig_SetIntervalIndex()\n");
  }
  return;
}

void pig_GetInterval(pIntervalGraph_t pig, int intIdx, int *C0, int *C1, int *nextIdx) {
  int *offset;
  if (intIdx > pig->intervalNum)
    fprintf(errfile,"PE%3d: ERROR: pig_GetInterval intIdx: %d intervalNum: %d\n",
	    MYNODE, intIdx, pig->intervalNum);
  offset = pig->data + (pig->entrNum + pig->exitNum)*pig_GetVSize() +
    intIdx*pig_GetISize();
  *C0 = *offset;
  *C1 = *(offset + 1);
  *nextIdx = *(offset + 2);
  return;
}

int pig_GetIntervalC0(pIntervalGraph_t pig, int intIdx) {
  int x;
  if (intIdx > pig->intervalNum)
    fprintf(errfile,"PE%3d: ERROR: pig_GetIntervalC0 intIdx: %d intervalNum: %d\n",
	    MYNODE, intIdx, pig->intervalNum);
  x = *(pig->data + (pig->entrNum + pig->exitNum)*pig_GetVSize() +
	intIdx*pig_GetISize());
  return x;
}

int pig_GetIntervalC1(pIntervalGraph_t pig, int intIdx) {
  int x;
  if (intIdx > pig->intervalNum)
    fprintf(errfile,"PE%3d: ERROR: pig_GetIntervalC1 intIdx: %d intervalNum: %d\n",
	    MYNODE, intIdx, pig->intervalNum);
  x = *(pig->data + (pig->entrNum + pig->exitNum)*pig_GetVSize() +
	intIdx*pig_GetISize() + 1);
  return x;
}

int pig_GetIntervalNext(pIntervalGraph_t pig, int intIdx) {
  int x;
  if (intIdx > pig->intervalNum)
    fprintf(errfile,"PE%3d: ERROR: pig_GetIntervalNext intIdx: %d intervalNum: %d\n",
	    MYNODE, intIdx, pig->intervalNum);
  x = *(pig->data + (pig->entrNum + pig->exitNum)*pig_GetVSize() +
	intIdx*pig_GetISize() + 2);
  return x;
}


void pig_SetIntervalC0(pIntervalGraph_t pig, int intIdx, int val) {
  if (intIdx > pig->intervalNum)
    fprintf(errfile,"PE%3d: ERROR: pig_SetIntervalC0 intIdx: %d intervalNum: %d\n",
	    MYNODE, intIdx, pig->intervalNum);
  *(pig->data + (pig->entrNum + pig->exitNum)*pig_GetVSize() +
    intIdx*pig_GetISize()) = val;
  return;
}

void pig_SetIntervalC1(pIntervalGraph_t pig, int intIdx, int val) {
  if (intIdx > pig->intervalNum)
    fprintf(errfile,"PE%3d: ERROR: pig_SetIntervalC1 intIdx: %d intervalNum: %d\n",
	    MYNODE, intIdx, pig->intervalNum);
  *(pig->data + (pig->entrNum + pig->exitNum)*pig_GetVSize() +
    intIdx*pig_GetISize() + 1) = val;
  return;
}

void pig_SetIntervalNext(pIntervalGraph_t pig, int intIdx, int val) {
  if (intIdx > pig->intervalNum)
    fprintf(errfile,"PE%3d: ERROR: pig_SetIntervalNext intIdx: %d intervalNum: %d\n",
	    MYNODE, intIdx, pig->intervalNum);
  *(pig->data + (pig->entrNum + pig->exitNum)*pig_GetVSize() +
    intIdx*pig_GetISize() + 2) = val;
  return;
}

int pig_GetEntrIdx(pIntervalGraph_t pig, int label, int adj) {
  int i;
  
  for (i=0 ; i<pig->entrNum ; i++) {
    if ((pig_GetLabel(pig, i, ENTRV) == label) &&
	(pig_GetAdj(pig, i, ENTRV) == adj))
      return(i);
  }

  fprintf(errfile,"PE%3d: ERROR: pig_GetEntrIdx()\n",MYNODE);
  return (-1);
}

void pig_Init(pIntervalGraph_t pig) {
  int i;
  
  pig->data    = (int *)SAFE_MALLOC(pig_GetDataInBytes(pig),
				    "(cycle.c) pig->data");

  pig->deadEntr = 0;
  pig->deadExit = 0;

  pig->deadEntrMask = (BOOL *)SAFE_MALLOC(pig->entrNum * sizeof(BOOL),
					  "(cycle.c) pig->deadEntrMask");
  for (i=0 ; i<pig->entrNum ; i++)
    *(pig->deadEntrMask + i) = FALSE;

  pig->deadExitMask = (BOOL *)SAFE_MALLOC(pig->exitNum * sizeof(BOOL),
					  "(cycle.c) pig->deadExitMask");
  for (i=0 ; i<pig->exitNum ; i++)
    *(pig->deadExitMask + i) = FALSE;

  pig->deadInt = 0;
  
  return;
}

void Create_IntervalPack(eGraph_t expGraph, pIntervalGraph_t *pIntervalGraph) {
  expVertex_t v;
  pIntervalGraph_t pig;
  int k;
  BOOL setFirstInterval;

  *pIntervalGraph = (pIntervalGraph_t)SAFE_MALLOC(sizeof(struct pIntervalGraph_s),
					"(cycle.c) *pIntervalGraph");
  pig = *pIntervalGraph;

  pig->entrNum      = expGraph->entrNum;
  pig->exitNum      = expGraph->exitNum;
  pig->intervalNum  = 1;

  /* Check: Change this later to accept more than one interval */

  pig_Init(pig);

  setFirstInterval = FALSE;
  v = expGraph->entrV;
  k=0;
  while (v != (expVertex_t)NULL) {
    pig_SetLabel  (pig, k, ENTRV, v->transArc->head);
    pig_SetAdj    (pig, k, ENTRV, v->transArc->tail);
    pig_SetAdjAssn(pig, k, ENTRV, v->transArc->tailAssn);
    pig_SetIntervalIndex(pig, k, ENTRV, 0);
    
    if (!setFirstInterval) {
      setFirstInterval = TRUE;
      pig_SetIntervalC0(pig, 0, v->C0);
      pig_SetIntervalC1(pig, 0, v->C1);
      pig_SetIntervalNext(pig, 0, -1);
    }
    k++;
    v = v->next;
  }

  v = expGraph->exitV;
  k=0;
  while (v != (expVertex_t)NULL) {
    pig_SetLabel  (pig, k, EXITV, v->transArc->tail);
    pig_SetAdj    (pig, k, EXITV, v->transArc->head);
    pig_SetAdjAssn(pig, k, EXITV, v->transArc->headAssn);
    pig_SetIntervalIndex(pig, k, EXITV, -1);
    k++;
    v = v->next;
  }

  return;
}


void Print_myIntervalPack(pIntervalGraph_t pig) {
  int j;
  int intIdx, C0, C1, nextIdx;
  
  fprintf(outfile,"PE%3d: Packed Interval Graph with %12d Entr and %12d exit:\n",
	  MYNODE, pig->entrNum, pig->exitNum);
  for (j=0 ; j<pig->entrNum ; j++) {
    fprintf(outfile,"PE%3d: \tEntrVertex[%12d]: l:%12d adj:%12d(%3d)   IntIdx:%12d\n",
	    MYNODE, j,
	    pig_GetLabel  (pig, j, ENTRV),
	    pig_GetAdj    (pig, j, ENTRV),
	    pig_GetAdjAssn(pig, j, ENTRV),
	    pig_GetIntervalIndex(pig, j, ENTRV));
    intIdx = pig_GetIntervalIndex(pig, j, ENTRV);
    while (intIdx >= 0) {
      pig_GetInterval(pig, intIdx, &C0, &C1, &nextIdx);
      fprintf(outfile,"PE%3d: \t\tInterval[%12d]: [%12d, %12d]\n",
	      MYNODE, intIdx, C0, C1);
      intIdx = nextIdx;
    }
  }
  for (j=0 ; j<pig->exitNum ; j++) {
    fprintf(outfile,"PE%3d: \tExitVertex[%12d]: l:%12d adj:%12d(%3d)   IntIdx:%12d\n",
	    MYNODE, j,
	    pig_GetLabel  (pig, j, EXITV),
	    pig_GetAdj    (pig, j, EXITV),
	    pig_GetAdjAssn(pig, j, EXITV),
	    pig_GetIntervalIndex(pig, j, EXITV));
  }
  fflush(outfile);
  return;
}

void Print_IntervalPack(pIntervalGraph_t pig) {
  int i;

  for (i=0 ; i<NODES ; i++) {
    if (MYNODE==i) {
      Print_myIntervalPack(pig);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  
  return;
}

void Free_IntervalPack(pIntervalGraph_t pig) {
  free(pig->deadEntrMask);
  free(pig->deadExitMask);
  if (pig->data != (int *)NULL)
    free(pig->data);
  free(pig);
  return;
}

void Send_IntervalPack(pIntervalGraph_t pig, int toNode) {
  int sizeBuf[3];

  sizeBuf[0] = pig->entrNum;
  sizeBuf[1] = pig->exitNum;
  sizeBuf[2] = pig->intervalNum;

  MPI_Send(sizeBuf, 3, MPI_INT, toNode, MPI_TAG_SIZES, MPI_COMM_WORLD);
  MPI_Send(pig->data, pig_GetDataInInts(pig), MPI_INT, toNode, MPI_TAG_DATA,
	   MPI_COMM_WORLD);
  
#ifdef DEBUG_PRINT
  fprintf(outfile,"PE%3d: I SENT THIS TO PE%3d\n",MYNODE,toNode);
  Print_myIntervalPack(pig);
  fprintf(outfile,"PE%3d: I SENT THIS **********\n",MYNODE);
  fflush(outfile);
#endif
  return;
}

void Recv_IntervalPack(pIntervalGraph_t *pIntervalGraph, int fromNode) {
  pIntervalGraph_t pig;
  int         sizeBuf[3];
  MPI_Status  mpstat;
  
  *pIntervalGraph = (pIntervalGraph_t)SAFE_MALLOC(sizeof(struct pIntervalGraph_s),
						  "(cycle.c) *pIntervalGraph");
  pig = *pIntervalGraph;

  MPI_Recv(sizeBuf, 3, MPI_INT, fromNode, MPI_TAG_SIZES, MPI_COMM_WORLD,
	   &mpstat);

  pig->entrNum      = sizeBuf[0];
  pig->exitNum      = sizeBuf[1];
  pig->intervalNum  = sizeBuf[2];

  pig_Init(pig);
  
  MPI_Recv(pig->data, pig_GetDataInInts(pig), MPI_INT, fromNode, MPI_TAG_DATA,
	   MPI_COMM_WORLD, &mpstat);

#ifdef DEBUG_PRINT
  fprintf(outfile,"PE%3d: I RECEIVED THIS FROM PE%3d\n",MYNODE,fromNode);
  Print_myIntervalPack(pig);
  fprintf(outfile,"PE%3d: I RECEIVED THIS **********\n",MYNODE);
  fflush(outfile);
#endif
  return;
}


BOOL pig_CheckArcHead(pIntervalGraph_t pig, int idx, int target) {
  int minval, maxval;
  int intIdx, nextIdx;
  BOOL found;

  intIdx = pig_GetIntervalIndex(pig, idx, ENTRV);
  found = FALSE;
  while (!found && (intIdx >= 0)) {
    pig_GetInterval(pig, intIdx, &minval, &maxval, &nextIdx);
    found = ((minval <= target) && (target <= maxval));
    intIdx = nextIdx;
  }
  
  return(found);
}


void pig_GetEntrPred(pIntervalGraph_t pig, int label, int *num, int **listPred) {
  int i, n;

  *listPred = (int *)SAFE_MALLOC(pig->entrNum * sizeof(int),
				 "(cycle.c) *listPred");

  n=0;
  for (i=0 ; i<pig->entrNum ; i++) {
    if (pig_CheckArcHead(pig, i, label)) {
      *(*listPred + n) = i; 
      n++;
    }
  }
    
  if (n > pig->entrNum)
    fprintf(errfile,"PE%3d: ERROR: pig_GetEntrPred(%3d), n:%d > pig->entrNum: %d\n",
	    MYNODE, label, n, pig->entrNum);
  *num = n;
  
  return;
}


BOOL pig_CheckInterval(pIntervalGraph_t pig, int C0, int C1) {
  int s;
  BOOL found;
  
  found = FALSE;
  s = 0;
  while ((!found) && (s<pig->intervalNum)) {
    found = ((pig_GetIntervalC0(pig, s) == C0) &&
	     (pig_GetIntervalC1(pig, s) == C1));
    s++;
  }
  
  return (found);
}

void pig_RemoveVertex(pIntervalGraph_t pig, int idx, BOOL vType) {
  int i, lab_exit;
  int minval, maxval, nextIdx;
  BOOL found;
  
#ifdef DEBUG_PRINT
  fprintf(outfile,"PE%3d: Removing %s vertex idx: %2d label: %2d\n",
	  MYNODE, (vType==ENTRV)?"ENTR":"EXIT",
	  idx, pig_GetLabel(pig,idx,vType));
#endif

  switch (vType) {
  case ENTRV:
    pig->deadEntr++;
    pig->deadEntrMask[idx] = TRUE;
#ifdef DEBUG_PRINT
    fprintf(outfile,"PE%3d: Removing %s vertex\n",
	    MYNODE, (vType==ENTRV)?"ENTR":"EXIT");
#endif
    break;
  case EXITV:
    pig->deadExit++;
    pig->deadExitMask[idx] = TRUE;
    lab_exit = pig_GetLabel(pig, idx, vType);
    /* CHECK:  */

    /* If we are the last exit vertex with this label, remove Interval arcs
       aimed at us */

    found = FALSE;
    i=0;
    while ((i<pig->exitNum) && (!found)) {
      found = ((pig_GetLabel(pig, i, EXITV) == lab_exit) &&
	       (pig->deadExitMask[i] == FALSE));
      i++;
    }
    if (!found) {
      for (i=0 ; i<pig->intervalNum ; i++) {
	pig_GetInterval(pig, i, &minval, &maxval, &nextIdx);
	if ((minval <= lab_exit) && (lab_exit <= maxval)) {
	  if (minval == maxval) {
	    /* Remove this interval */
	    pig_SetIntervalC0(pig, i, -1);
	    pig_SetIntervalC1(pig, i, -1);
	    pig->deadInt++;
	  }
	  else {
	    if (minval == lab_exit)
	      pig_SetIntervalC0(pig, i, minval+1);
	    else {
	      if (maxval == lab_exit)
		pig_SetIntervalC1(pig, i, maxval-1);
	    }
	  }
	}
      }
    }
    break;
  default:
    fprintf(errfile,"ERROR: pig_RemoveVertex()\n");
  }
  return;
}



void pig_CleanUp(pIntervalGraph_t pig) {
  /* Merge in the new Intervals and delete all vertices with marked labels */
  pIntervalGraph_t pigNew;
  int vIdx, i;

  pigNew = (pIntervalGraph_t)SAFE_MALLOC(sizeof(struct pIntervalGraph_s),
					 "(cycle.c) pigNew");
  
  pigNew->entrNum = pig->entrNum - pig->deadEntr;
  pigNew->exitNum = pig->exitNum - pig->deadExit;
  pigNew->intervalNum = pig->intervalNum;

#ifdef DEBUG_PRINT
  fprintf(outfile,"PE%3d: dead entr: %3d dead exit: %3d dead int: %3d \n",
	  MYNODE, pig->deadEntr, pig->deadExit, pig->deadInt);
#endif

  pig_Init(pigNew);

  /* CHECK */

  /* Entrance Vertices */
  vIdx = 0;
#ifdef DEBUG_PRINT
  fprintf(outfile,"PE%3d: cu Entrance vertices (num %d)\n",MYNODE, pig->entrNum);
#endif
  for (i=0 ; i<pig->entrNum ; i++) {
    if (pig->deadEntrMask[i] == FALSE) {
#ifdef DEBUG_PRINT
  fprintf(outfile,"PE%3d: cu working on idx: %d (label: %d)\n",MYNODE,
	  i, pig_GetLabel(pig, i, ENTRV));
#endif
      pig_SetLabel  (pigNew, vIdx, ENTRV, pig_GetLabel  (pig, i, ENTRV));
      pig_SetAdj    (pigNew, vIdx, ENTRV, pig_GetAdj    (pig, i, ENTRV));
      pig_SetAdjAssn(pigNew, vIdx, ENTRV, pig_GetAdjAssn(pig, i, ENTRV));
      pig_SetIntervalIndex(pigNew, vIdx, ENTRV,
			   pig_GetIntervalIndex(pig, i, ENTRV));
      vIdx++;
    }
  }

  if (vIdx != pigNew->entrNum)
    fprintf(errfile,"PE%3d: ERROR: vIdx (%d) != pigNew->entrNum (%d)\n",
	    MYNODE, vIdx, pigNew->entrNum);

#ifdef DEBUG_PRINT
  fprintf(outfile,"PE%3d: cu entrance done\n",MYNODE);
#endif

  /* Exit Vertices */
  vIdx   = 0;

  for (i=0 ; i<pig->exitNum ; i++) {
    if (pig->deadExitMask[i] == FALSE) {
      pig_SetLabel  (pigNew, vIdx, EXITV, pig_GetLabel  (pig, i, EXITV));
      pig_SetAdj    (pigNew, vIdx, EXITV, pig_GetAdj    (pig, i, EXITV));
      pig_SetAdjAssn(pigNew, vIdx, EXITV, pig_GetAdjAssn(pig, i, EXITV));
      pig_SetIntervalIndex(pigNew, vIdx, EXITV,
			   pig_GetIntervalIndex(pig, i, EXITV));
      vIdx++;
    }
  }

  if (vIdx != pigNew->exitNum)
    fprintf(errfile,"PE%3d: ERROR: vIdx (%d) != pigNew->exitNum (%d)\n",
	    MYNODE, vIdx, pigNew->exitNum);

#ifdef DEBUG_PRINT
  fprintf(outfile,"PE%3d: cu exit done (num: %d) \n",MYNODE, pig->exitNum);
#endif

#ifdef DEBUG_PRINT
  fprintf(outfile,"PE%3d: Copying over %d intervals\n",MYNODE, pigNew->intervalNum);
#endif

  memcpy(pigNew->data + (pigNew->entrNum + pigNew->exitNum)*pig_GetVSize(),
	 pig->data + (pig->entrNum + pig->exitNum)*pig_GetVSize(),
	 pigNew->intervalNum * pig_GetISize() * sizeof(int));

  /* Swap in new data and parameters */
  free(pig->data);
  pig->entrNum  = pigNew->entrNum;
  pig->exitNum  = pigNew->exitNum;
  pig->intervalNum   = pigNew->intervalNum;
  pig->data     = pigNew->data;
#if 1
  pigNew->data  = (int *)NULL;
#endif
  pig->deadEntr = pigNew->deadEntr;
  pig->deadExit = pigNew->deadExit;
  pig->deadInt  = pigNew->deadInt;
  
  pig->deadEntrMask = (BOOL *)SAFE_MALLOC(pig->entrNum * sizeof(BOOL),
					  "(cycle.c) pig->deadEntrMask");
  for (i=0 ; i<pig->entrNum ; i++)
    *(pig->deadEntrMask + i) = FALSE;

  pig->deadExitMask = (BOOL *)SAFE_MALLOC(pig->exitNum * sizeof(BOOL),
					  "(cycle.c) pig->deadExitMask");
  for (i=0 ; i<pig->exitNum ; i++)
    *(pig->deadExitMask + i) = FALSE;

#ifdef DEBUG_PRINT
  fprintf(outfile,"PE%3d: cu swapped data structure\n",MYNODE);
#endif

  Free_IntervalPack(pigNew);
  
#ifdef DEBUG_PRINT
  fprintf(outfile,"PE%3d: cu Freed IntervalPack\n",MYNODE);
#endif

  return;
}


void IntervalMerge(pIntervalGraph_t *p0Ptr, pIntervalGraph_t p1, int h) {
  /* Merge p1 into p0 */
  pIntervalGraph_t pTemp;
  pIntervalGraph_t p0;
  int *ptr, *loc0, *loc1;
  int siz;
  int i;
  int curAdjAssn, origAdjAssn;
  int exit0, exit0Idx, entr1, entr1Idx;
  int intIdx, intOffset;

  BOOL *mark;
  BOOL seenLoop;
  int j, C0, C1, nextIdx, entr1intIdx;
    
  p0 = *p0Ptr;
  
#ifdef DEBUG_PRINT
  fprintf(outfile,"PE%3d: MY P0  MERGE AT THE START OF THINGS\n",MYNODE);
  Print_myIntervalPack(p0);
  fprintf(outfile,"PE%3d: MY P0  MERGE AT THE START OF THINGS DONE\n",MYNODE);
  fflush(outfile);
#endif

  pTemp = (pIntervalGraph_t)SAFE_MALLOC(sizeof(struct pIntervalGraph_s),
					"(cycle.c) pTemp");
  
  pTemp->entrNum      = p0->entrNum      + p1->entrNum;
  pTemp->exitNum      = p0->exitNum      + p1->exitNum;
  pTemp->intervalNum  = p0->intervalNum  + p1->intervalNum;

  pig_Init(pTemp);

  ptr  = pTemp->data;

  siz  = p0->entrNum*pig_GetVSize();
  loc0 = p0->data;
  memcpy(ptr, loc0, siz*sizeof(int));
  ptr += siz;

  siz  = p1->entrNum*pig_GetVSize();
  loc1 = p1->data;
  memcpy(ptr, loc1, siz*sizeof(int));
  ptr += siz;

  siz   = p0->exitNum*pig_GetVSize();
  loc0 += p0->entrNum*pig_GetVSize();
  memcpy(ptr, loc0 , siz*sizeof(int));
  ptr += siz;

  siz   = p1->exitNum*pig_GetVSize();
  loc1 += p1->entrNum*pig_GetVSize();
  memcpy(ptr, loc1, siz*sizeof(int));
  ptr += siz;
  
  siz   = p0->intervalNum*pig_GetISize();
  loc0 += p0->exitNum*pig_GetVSize();
  memcpy(ptr, loc0 , siz*sizeof(int));
  ptr += siz;

  siz   = p1->intervalNum*pig_GetISize();
  loc1 += p1->exitNum*pig_GetVSize();
  memcpy(ptr, loc1, siz*sizeof(int));
  ptr += siz;

  intOffset = p0->intervalNum;

  for (i=p0->entrNum ; i<p0->entrNum+p1->entrNum ; i++) {
    intIdx = pig_GetIntervalIndex(pTemp, i, ENTRV);
    if (intIdx >= 0) 
      pig_SetIntervalIndex(pTemp, i, ENTRV, intIdx + intOffset);
  }

  for (i=p0->intervalNum ; i<p0->intervalNum+p1->intervalNum ; i++) {
    intIdx = pig_GetIntervalNext(pTemp, i);
    if (intIdx >= 0) 
      pig_SetIntervalNext(pTemp, i, intIdx + intOffset);
  }


#ifdef DEBUG_PRINT
  fprintf(outfile,"PE%3d: MY P0  MERGE\n",MYNODE);
  Print_myIntervalPack(p0);
  fprintf(outfile,"PE%3d: MY P0  MERGE DONE\n",MYNODE);
  fflush(outfile);
  fprintf(outfile,"PE%3d: MY P1  MERGE\n",MYNODE);
  Print_myIntervalPack(p1);
  fprintf(outfile,"PE%3d: MY P1  MERGE DONE\n",MYNODE);
  fflush(outfile);

  fprintf(outfile,"PE%3d: MY P0+P1    MERGE\n",MYNODE);
  Print_myIntervalPack(pTemp);
  fprintf(outfile,"PE%3d: MY P0+P1    MERGE DONE\n",MYNODE);
  fflush(outfile);
#endif

  /****** MERGE HERE ********/
  mark = (BOOL *)SAFE_MALLOC(pTemp->intervalNum * sizeof(BOOL),
			     "(cycleInterval.c) mark");

  for (i=0 ; i<pTemp->exitNum ; i++) {
    exit0Idx = i;
    origAdjAssn = pig_GetAdjAssn(pTemp, exit0Idx, EXITV);
    curAdjAssn  = clearLastB(origAdjAssn, h);
#ifdef DEBUG_PRINT
    fprintf(outfile,"PE%3d: h: %d orig: %d cur: %d\n",
	    MYNODE, h, origAdjAssn, curAdjAssn);
#endif
    if (curAdjAssn == MYNODE) {
      /* we have a winner */

      /* must locate the set of entrance vertices pointing to us,
	 the entrance vertex we point to, and
	 the set of exit vertices that it points to */

      /* entr0 ------> exit0 ------> entr1 ------> exit1 */
      /*         exit0Idx,EXITV     entr1Idx,ENTRV       */
      
      /* exit0Idx,EXITV = index of the exit vertex */
      /* exit0 = label of the exit vertex */
      exit0 = pig_GetLabel(pTemp, exit0Idx, EXITV);

      /* entr1 = label of the entrance vertex we point to */
      entr1     = pig_GetAdj    (pTemp, exit0Idx, EXITV);

      /* entr1Idx,ENTRV = index of the entrance vertex */
      entr1Idx  = pig_GetEntrIdx(pTemp, entr1, exit0);
#ifdef DEBUG_PRINT
      fprintf(outfile,"PE%3d: h: %d orig: %d cur: %d exit0: %3d entr1: %3d entr1Idx: %3d\n",
	    MYNODE, h, origAdjAssn, curAdjAssn, exit0, entr1, entr1Idx);
#endif

      /* If the Interval arc exists from entr1 to exit0, we have a cycle! */

      if (pig_CheckArcHead(pTemp, entr1Idx, exit0)) {
	fprintf(outfile,"PE%3d: HALT: CYCLE DETECTED DURING MERGE (%3d, %3d)\n",
		MYNODE, exit0, entr1);
	CYCLE_FOUND = TRUE;
	return;
      }

      /* Must connect the intervals from entr1 to the intervals that contain exit0 */
      /* Trace each interval list that matches to its end. If we haven't seen
	 the new interval contained in one of these, then add it as new, and
	 point the last index to it */
	 
      entr1intIdx = pig_GetIntervalIndex(pTemp, entr1Idx, ENTRV);

      for (j=0 ; j<pTemp->intervalNum ; j++)
	mark[j] = FALSE;

      for (j=0 ; j<pTemp->intervalNum ; j++) {
	if (!mark[j]) {
	  mark[j] = TRUE;
	  intIdx = j;
	  pig_GetInterval(pTemp, intIdx, &C0, &C1, &nextIdx);
	  if ((C0 <= exit0) && (exit0 <= C1)) {

	    seenLoop = FALSE;
	    while ((intIdx != entr1intIdx)  && (nextIdx >= 0) && !seenLoop) {
	      intIdx = nextIdx;
	      if (mark[intIdx])
		seenLoop = TRUE;
	      else {
		mark[intIdx] = TRUE;
		pig_GetInterval(pTemp, intIdx, &C0, &C1, &nextIdx);
	      }
	      
	    }

	    if ((nextIdx < 0) && (intIdx != entr1intIdx) && !seenLoop) {
	      pig_SetIntervalNext(pTemp, intIdx, entr1intIdx);
#ifdef DEBUG_PRINT
	      fprintf(outfile,"PE%3d: (%d -> %d) Linking interval %3d -> %3d\n",
		      MYNODE, exit0, entr1, intIdx, entr1intIdx);
#endif
	    }
	  
	  }
	}
      }
      
      pig_RemoveVertex(pTemp, exit0Idx, EXITV);
      pig_RemoveVertex(pTemp, entr1Idx, ENTRV);
	 
    } /* (curAdjAssn  == MYNODE) */

  } /* foreach exit0 */

  pig_CleanUp(pTemp);
  
#ifdef DEBUG_PRINT
  Print_myIntervalPack(pTemp);
#endif

  free(mark);
  
  *p0Ptr = pTemp;

  Free_IntervalPack(p0);
  return;
}

void Merge_IntervalPack(pIntervalGraph_t *pig) {
  int h;
  int logp;
  pIntervalGraph_t pIntervalGraph_recv;
  int found, result;

  /* logp is the ceiling of log2(NODES) */
  logp = log2_i(NODES);
  for (h=0 ; h<logp ; h++) {
#if PRINT_MERGESIZE
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    if (lastB(MYNODE, h)==0) {
#if PRINT_MERGESIZE
      fprintf(outfile,"PE%3d: [%3d] entr: %12d exit: %12d exarcs: %12d\n",
	      MYNODE, h, (*pig)->entrNum, (*pig)->exitNum, (*pig)->arcNum);
#endif
      if (testB(MYNODE,h)==0) {
	/* Make sure that if we're not using a power-of-two nodes
	   that we just sit idle for a stage where our partner is missing */
	if (setB(MYNODE,h) < NODES) {
#ifdef DEBUG_PRINT
	  fprintf(outfile,"PE%3d: (h:%3d) I'm receiving from PE%3d\n",
		  MYNODE, h, setB(MYNODE,h));
	  fflush(outfile);
#endif
	  Recv_IntervalPack(&pIntervalGraph_recv, setB(MYNODE,h));
	  IntervalMerge(pig, pIntervalGraph_recv, h);

	  Free_IntervalPack(pIntervalGraph_recv);
	}
      }
      else {
#ifdef DEBUG_PRINT
	fprintf(outfile,"PE%3d: (h:%3d) I'm sending to PE%3d\n",
		MYNODE, h, clearB(MYNODE,h));
	fflush(outfile);
#endif
	Send_IntervalPack(*pig, clearB(MYNODE,h));
      }
    }
#if 1
    found = (CYCLE_FOUND==TRUE)?1:0;
    MPI_Allreduce(&found, &result, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if (result)
      return;
#else
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }
  return;
}


/*************************************************************************/

void Report_Cycles() {
  int found, result;
  
  MPI_Barrier(MPI_COMM_WORLD);

  found = ((CYCLE_FOUND==TRUE)?1:0);
  
  MPI_Reduce(&found, &result, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
  
  if (MYNODE==0) {
    if (result > 0)
      fprintf(outfile,"RESULT: A cycle was found!\n");
    else
      fprintf(outfile,"RESULT: No cycles were found.\n");
  }
  
  fflush(outfile);
  return;
}

void main(int argc, char **argv) {
  int n;
  vertexList_t myVerts;
  int initTransArcCount, termTransArcCount;
  transArc_t *initTransArcs, *termTransArcs;
  eGraph_t expGraph;
  pIntervalGraph_t pIntervalGraph;
  int found, result;

  timer_init();
  
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &MYNODE);
  MPI_Comm_size(MPI_COMM_WORLD, &NODES);

  outfile = stdout;
#if 0
  errfile = stderr;
#else
  errfile = outfile;
#endif

  switch (argc) {
  case 1:
    n = N_DEFAULT; break;
  case 2:
    n = atoi(argv[1]); break;
  default:
    if (MYNODE==0)
      fprintf(errfile,"Usage '%s n'\n",argv[0]);
    n = N_DEFAULT;
  }

  if (MYNODE==0) {
    fprintf(outfile,"NODES: %3d\n",NODES);
    fprintf(outfile,"n: %12d\n",n);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  CYCLE_FOUND = FALSE;
  
  timer_reset();
  
  Create_Input(n, &myVerts);

  if (MYNODE==0) {
    fprintf(outfile,"Input: %s and %12d vertices\n",INPUT_GRAPH_NAME,n);
  }

  timer_mark("Create_Input");

#ifdef DEBUG_PRINT
  Print_Input(myVerts);
  timer_mark("Print Input");
#endif

  Find_Local_Cycles(myVerts,
		    &initTransArcCount,
		    &initTransArcs);
  timer_mark("Find Local Cycles");

  found = (CYCLE_FOUND==TRUE)?1:0;
  MPI_Allreduce(&found, &result, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  if (!result) {
  
    Discovery(myVerts,
	      initTransArcCount, 
	      initTransArcs,
	      &termTransArcCount,
	      &termTransArcs);
    timer_mark("Discovery");

    Create_Express(myVerts,
		   initTransArcCount, 
		   initTransArcs,
		   termTransArcCount,
		   termTransArcs,
		   &expGraph);
    timer_mark("Create_Express");

#ifdef DEBUG_PRINT
    Print_Express(expGraph);
    timer_mark("Print_Express");
#endif
    
    Create_IntervalPack(expGraph, &pIntervalGraph);
    timer_mark("Create_IntervalPack");

#ifdef DEBUG_PRINT
    Print_IntervalPack(pIntervalGraph);
    timer_mark("Print_IntervalPack");
#endif
  
    Merge_IntervalPack(&pIntervalGraph);
    timer_mark("Merge_IntervalPack");
    
    Free_IntervalPack(pIntervalGraph);
    
    Free_Express(expGraph);
    free(termTransArcs);
  }

  Report_Cycles();
  timer_report(outfile,"CYCLE DETECTION", n);

  MPI_Barrier(MPI_COMM_WORLD);

  free(initTransArcs);
  Free_Input(myVerts);

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  
}
