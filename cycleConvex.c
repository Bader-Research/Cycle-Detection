/*

   $Log: cycleConvex.c,v $
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

static char rcsid[] = "$Id: cycleConvex.c,v 1.6 1999/07/11 15:42:13 dbader Exp $";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
#include "misc.h"
#include "cycleConvex.h"
#include "queue.h"
#include "mpi-printf.h"
#include "timing.h"

#define N_DEFAULT  (1<<10)
#define MAXARCS     4
#define INPUTCYCLE  1

#define PRINT_MERGESIZE 0

#define USE_PEG 0

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

static int intCompare(const void *a, const void *b) {
  return (*(int *)a - *(int *)b);
}

static int transArcCompareHead(const void *a, const void *b) {
  return (((transArc_t *)a)->headAssn - ((transArc_t *)b)->headAssn);
}

static int expGraphLabelCompare(const void *a, const void *b) {
  return (((labelLookup_t)a)->label - ((labelLookup_t)b)->label);
}

static int pArcFlatCompare(const void *a, const void *b) {
  return (((pArcFlat_t)a)->idx - ((pArcFlat_t)b)->idx);
}

static int pIntervalFlatCompare(const void *a, const void *b) {
  return (((pIntervalFlat_t)a)->idx - ((pIntervalFlat_t)b)->idx);
}

void Create_Input0(int n, vertexList_t *verts) {
  int num;
  int offset;
  int i, k;
  vertex_t v;
  arc_t a;
  
  MPI_Barrier(MPI_COMM_WORLD);

  INPUT_GRAPH_NAME = "Input0";

  *verts = (vertexList_t)SAFE_MALLOC(sizeof(struct vertexList_s),
				     "(cycle.c) *verts");

  num = (*verts)->num = (int)floor((double)n/(double)NODES);
  (*verts)->vlist = (vertex_t)SAFE_MALLOC(num*sizeof(struct vertex_s),
					  "(cycle.c) (*verts)->vlist");

#ifdef USE_RAND
  srand(time(0)+MYNODE);
#else
  srandom(time(0)+MYNODE);
#endif
  offset = num*MYNODE;
  for (i=0 ; i<num ; i++) {
    v = (*verts)->vlist + i;
    v->label = offset + i;
    v->arcs  = randmax(MAXARCS);
    v->alist = (arc_t)SAFE_MALLOC(MAXARCS*sizeof(struct arc_s),
				  "(cycle.c) v->alist");
    for (k=0 ; k<v->arcs ; k++) {
      a = v->alist + k;
#if 0
      a->head = randrange(offset,offset+num);
      a->assn = MYNODE;
#else
      if (MYNODE==0 && k==0) {
	a->head = offset+num+i;
	a->assn = MYNODE+1;
      }
      else {
	a->head = min(offset+num-1 , v->label + k + 1);
	if (a->head == v->label)
	  v->arcs = k;
	a->assn = MYNODE;
      }
#endif
    }
  }
  
  (*verts)->color = (unsigned char *)SAFE_MALLOC(num*sizeof(unsigned char),
						 "(cycle.c) (*verts)->color");

  (*verts)->R = (reachlist_t)SAFE_MALLOC(num*sizeof(struct reachlist_s),
					 "(cycle.c) (*verts)->R");
  
  return;
}

void Create_Input_Lattice_Orig(int n, vertexList_t *verts) {
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
    v->label = offset + i;

    vLoc    = i;
    vLocRow = vLoc / nCol;
    vLocCol = vLoc % nCol;
    
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
    if (v->label == (num*NODES)-1) {
      v->arcs++;
    }
#endif
    v->alist = (arc_t)SAFE_MALLOC((v->arcs)*sizeof(struct arc_s),
				  "(cycle.c) v->alist");
    a = v->alist;
    if (vLocCol < nCol - 1) {
      a->head = v->label + 1;
      a->assn = MYNODE;
      a++;
    }
    if (vLocRow < nRow - 1) {
      a->head = v->label + nCol;
      a->assn = MYNODE;
      a++;
    }
    if ((vLocCol == nCol-1) && (myCol < pCol-1)) {
      a->head = v->label + num - nCol + 1;
      a->assn = MYNODE + 1;
      a++;
    }
    if ((vLocRow == nRow-1) && (myRow < pRow-1)) {
      a->head = v->label + (pCol - 1)*num + nCol;
      a->assn = MYNODE + pCol;
      a++;
    }
#if INPUTCYCLE
    /* Add arc from last vertex to first vertex */
    if (v->label == (num*NODES)-1) {
      a->head = 0;
      a->assn = 0;
      a++;
    }
#endif
  }
  
  (*verts)->color = (unsigned char *)SAFE_MALLOC(num*sizeof(unsigned char),
						 "(cycle.c) (*verts)->color");

  (*verts)->R = (reachlist_t)SAFE_MALLOC(num*sizeof(struct reachlist_s),
					 "(cycle.c) (*verts)->R");
  
  return;
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

  (*verts)->R = (reachlist_t)SAFE_MALLOC(num*sizeof(struct reachlist_s),
					 "(cycle.c) (*verts)->R");
  
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

  (*verts)->R = (reachlist_t)SAFE_MALLOC(num*sizeof(struct reachlist_s),
					 "(cycle.c) (*verts)->R");
  
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

  (*verts)->R = (reachlist_t)SAFE_MALLOC(num*sizeof(struct reachlist_s),
					 "(cycle.c) (*verts)->R");
  
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

void visitLocal(vertexList_t verts, int vidx) {
  int i;
  vertex_t v;
  arc_t a;
  int aidx;
  
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
	visitLocal(verts, aidx);
	break;
      case RED:
	fprintf(outfile,"PE%3d: A local cycle (%d, %d) has been found!\n",
		MYNODE, v->label, a->head);
	CYCLE_FOUND = TRUE;
	break;
      case BLACK:
      default: ;
      }
    }
  }

  verts->color[vidx] = BLACK;
  
  return;
}

void Find_Local_Cycles(vertexList_t verts) {
  int i;
  int num;

  num = verts->num;
  
  for (i=0 ; i<num ; i++)
    verts->color[i] = WHITE;

  for (i=0 ; i<num ; i++) {
    if ((verts->color[i] == WHITE) && (!CYCLE_FOUND)) {
      /* vertex i has not yet been visited */
      visitLocal(verts, i);
    }
  }

  return;
}

void ReachList_Init(reachlist_t R) {
  R->next  = (reachlist_t)NULL;
  return;
}

void Free_ReachList(reachlist_t R) {
  reachlist_t link, lastlink;

  link = R->next;
  while (link != (reachlist_t)NULL) {
    lastlink = link;
    link = link->next;
    free(lastlink);
  }

  free(R);
  
  return;
}

BOOL ReachList_Empty(reachlist_t R) {
  return (R->next == (reachlist_t)NULL);
}

void ReachList_Add_NoCheck(reachlist_t R, int idx) {
  reachlist_t newLink;
  
  newLink = (reachlist_t)SAFE_MALLOC(sizeof(struct reachlist_s),
				     "(cycle.c) newLink");

  newLink->idx  = idx;
  newLink->next = R->next;
  R->next  = newLink;
  return;
}

void ReachList_Add(reachlist_t R, int idx) {
  reachlist_t RPtr;
  BOOL found;

  found = FALSE;
  
  if (R->next != (reachlist_t)NULL) {
    RPtr = R->next;
    do {
      if (RPtr->idx == idx) 
	found = TRUE;
      RPtr = RPtr->next;
    } while ((RPtr != (reachlist_t)NULL) && !found);
  }

  if (!found) 
    ReachList_Add_NoCheck(R, idx);
  
  return;
}

void ReachList_Merge(reachlist_t Ra, reachlist_t Rb) {
  /* Merge ReachList Rb into Ra */
  reachlist_t RPtr;
  if (Rb->next != NULL) {
    /* Rb is non-trivial */
    if (Ra->next == NULL) {
      /* Copy Rb to Ra */
      RPtr = Rb->next;
      while (RPtr != NULL) {
	ReachList_Add_NoCheck(Ra, RPtr->idx);
	RPtr = RPtr->next;
      }
    }
    else {
      /* Merge Rb into Ra */
      RPtr = Rb->next;
      while (RPtr != NULL) {
	ReachList_Add(Ra, RPtr->idx);
	RPtr = RPtr->next;
      }
    }
  }
  return;
}

reachlist_t visitDiscovery(vertexList_t verts, int vidx, Q_t Q) {
  int i;
  int widx;
  vertex_t v;
  arc_t a;
  reachlist_t Rv, Rw;
  QELEM_T tmpElem;

  Rv = verts->R + vidx;
  ReachList_Init(Rv);

  verts->color[vidx] = RED;

  v = verts->vlist + vidx;
  for (i=0 ; i<v->arcs ; i++) {
    a = v->alist + i;
    if (a->assn == MYNODE) {  /* local-arc incident from v */
      widx = localIndex(verts, a->head);
      switch (verts->color[widx]) {
      case WHITE:
	Rw = visitDiscovery(verts, widx, Q);
	ReachList_Merge(Rv, Rw);
	break;
      case RED: fprintf(outfile,"PE%3d: A local cycle has been found!\n",
			MYNODE);
	break;
      case BLACK:
	ReachList_Merge(Rv, verts->R + widx);
	break;
      case GREEN:
	break;
      default:
	fprintf(errfile,
		"PE%3d: ERROR! in default case of visitDiscovery() %d\n",
		MYNODE, verts->color[widx]);
	;
      }
    }
    else { /* trans-arc */
      tmpElem.tail = v->label;
      tmpElem.head = a->head;
      tmpElem.tailAssn = MYNODE;
      tmpElem.headAssn = a->assn;
      Qadd(Q, &tmpElem); 
      ReachList_Add(Rv, vidx);
    }
  }

  if (ReachList_Empty(Rv)) {
    verts->color[vidx] = GREEN;
  }
  else {
    verts->color[vidx] = BLACK;
  }
  
  return Rv;
}

void Discovery(vertexList_t verts,
	       int *initTransArcCount,
	       transArc_t **initTransArcs,
	       int *termTransArcCount,
	       transArc_t **termTransArcs)
{
  int i;
  int num;
  Q_t Q;
  transArc_t tarc;
  int *send_cnt, *send_offset;
  int *recv_cnt, *recv_offset;
  MPI_Datatype MPI_TRANSARC;
  
  num = verts->num;
  
  Qinit(&Q, MAXARCS*num);
  
  for (i=0 ; i<num ; i++) 
    verts->color[i] = WHITE;

  for (i=0 ; i<num ; i++) {
    if (verts->color[i] == WHITE) { /* vertex i has not yet been visited */
      visitDiscovery(verts, i, Q);
    }
  }

  *initTransArcCount = Qsize(Q);

#ifdef DEBUG_PRINT
  MPI_fprintf(outfile,"I have %d initial transarc%s\n", *initTransArcCount,
	      (*initTransArcCount==1)?"":"s");
#endif
  
  if (*initTransArcCount > 0) {

    /* Sort the Q */
    *initTransArcs =
      (transArc_t *)SAFE_MALLOC(*initTransArcCount * sizeof(transArc_t),
				"(cycle.c) *initTransArcs");
    i = 0;
    while (!Qempty(Q)) {
      Qdel(Q, &tarc);
      transArcCopy((*initTransArcs) + i , &tarc);
      i++;
    };
    if (i != *initTransArcCount)
      fprintf(errfile,"PE%3d: ERROR: (i != *initTransArcCount) \n",MYNODE);

    qsort(*initTransArcs, *initTransArcCount,
	  sizeof(transArc_t), transArcCompareHead);

  }
  else {
    *initTransArcs = (transArc_t *)NULL;
  }

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
  
  for (i=0 ; i<*initTransArcCount ; i++)
    send_cnt[((*initTransArcs)+i)->headAssn]++;

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

  MPI_Alltoallv(*initTransArcs, send_cnt, send_offset, MPI_TRANSARC,
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
  Qfree(Q);
  return;
}

void AddEntrVertex(eGraph_t expGraph, transArc_t *tarc) {
  expVertex_t entrVertex;
  
  entrVertex = (expVertex_t)SAFE_MALLOC(sizeof(struct expVertex_s),
					"(cycle.c) entrVertex");

  entrVertex->transArc   = tarc;
  entrVertex->expArcsNum = 0;
  entrVertex->expArcs    = (expArc_t)NULL;
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
  exitVertex->expArcsNum = 0;
  exitVertex->expArcs    = (expArc_t)NULL;
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

reachlist_t getReachList(vertexList_t verts, int label) {
  int idx;

  idx = localIndex(verts,label);
  
  return (verts->R + idx);
}

void AddExpressArc(expVertex_t tailVertex, expVertex_t headVertex) {
  expArc_t newArc;

  newArc = (expArc_t)SAFE_MALLOC(sizeof(struct expArc_s),
				 "(cycle.c) newArc");

  newArc->headVertex  = headVertex;
  newArc->next        = tailVertex->expArcs;

  tailVertex->expArcs = newArc;
  tailVertex->expArcsNum++;

  return;
}

void AddExpressArcSet(vertexList_t verts, eGraph_t expGraph,
		      expVertex_t entrVertex) {
  int vlabel;
  int reachLabelIdx;
  int reachLabel;
  expVertex_t exitVertex;
  reachlist_t rlist;
  
  vlabel = entrVertex->transArc->head;
  rlist = getReachList(verts, vlabel);

  if (rlist->next != (reachlist_t)NULL) {
    rlist = rlist->next;
  }
  else
    rlist = (reachlist_t)NULL;
  
  while (rlist != (reachlist_t)NULL) {
    reachLabelIdx = rlist->idx;
    reachLabel    = globalLabel(verts, reachLabelIdx);
    exitVertex    = getExitVertexLUT(expGraph, reachLabel);
    if (exitVertex == (expVertex_t)NULL)
      fprintf(errfile,"PE%3d: ERROR: getExitVertexLUT returned NULL\n", MYNODE);
    
    AddExpressArc(entrVertex, exitVertex);
    
    rlist = rlist->next;
  } 

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
  
  Free_ReachList(verts->R);
  /*  free(verts->R); */
  free(verts->color);
  
  for (i=verts->num-1 ; i>=0 ; i--)
    if ((verts->vlist + i)->alist != (arc_t)NULL)
      free((verts->vlist + i)->alist);

  free(verts->vlist);
  free(verts);
  
  return;
}

void Free_expVertex(expVertex_t v) {
  expArc_t a, nexta;

  a = v->expArcs;
  while (a != (expArc_t)NULL) {
    nexta = a->next;
    free(a);
    a = nexta;
  }
  
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
  expArc_t a;

  num = 0;
  v = expGraph->entrV;
  while (v != (expVertex_t)NULL) {
    a = v->expArcs;
    while (a != (expArc_t)NULL) {
      num++;
      a = a->next;
    }
    v = v->next;
  }

  return num;
}

void Print_Express(eGraph_t expGraph) {
  int i, j, k;
  expVertex_t v;
  expArc_t a;
  
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
	fprintf(outfile,"PE%3d: \t\texpArcs: %12d\n",
		MYNODE, v->expArcsNum);
	if (v->transArc->headAssn != MYNODE)
	  fprintf(errfile,"PE%3d: ERROR: headAssn != MYNODE (%d, %d)\n", MYNODE,
		  v->transArc->headAssn, MYNODE);
	a = v->expArcs;
	j=0;
	while (a != (expArc_t)NULL) {
	  fprintf(outfile,"PE%3d: \t\t  Arc [%12d]: (%12d, %12d)\n",
		  MYNODE, j, v->transArc->head, a->headVertex->transArc->tail);
	  j++;
	  a = a->next;
	}
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

int peg_GetVSize() {
  return 4;
}

int peg_GetASize() {
  return 1;
}

int peg_GetDataInInts_Var(int entrNum, int exitNum, int arcNum) {
  return (peg_GetVSize()*(entrNum+exitNum) + peg_GetASize()*arcNum);
}

int peg_GetDataInInts(pExpGraph_t peg) {
  return (peg_GetDataInInts_Var(peg->entrNum, peg->exitNum, peg->arcNum));
}

int peg_GetDataInBytes(pExpGraph_t peg) {
  return (peg_GetDataInInts(peg)*sizeof(int));
}

int peg_GetLabel(pExpGraph_t peg, int idx, BOOL vType) {
  int x;
  switch (vType) {
  case ENTRV:
    x = *( peg->data + (idx*peg_GetVSize()) );
    break;
  case EXITV:
    x = *( peg->data + ((peg->entrNum + idx)*peg_GetVSize()) );
    break;
  default:
    x = -1;
    fprintf(errfile,"ERROR: peg_GetLabel()\n");
  }
  return x;
}

void peg_SetLabel(pExpGraph_t peg, int idx, BOOL vType, int val) {
  switch (vType) {
  case ENTRV:
    *( peg->data + (idx*peg_GetVSize()) ) = val;
    break;
  case EXITV:
    *( peg->data + ((peg->entrNum + idx)*peg_GetVSize()) ) = val;
    break;
  default:
    fprintf(errfile,"ERROR: peg_SetLabel()\n");
  }
  return;
}

int peg_GetAdj(pExpGraph_t peg, int idx, BOOL vType) {
  int x;
  switch (vType) {
  case ENTRV:
    x = *( peg->data + (idx*peg_GetVSize()) + 1);
    break;
  case EXITV:
    x = *( peg->data + ((peg->entrNum + idx)*peg_GetVSize()) + 1);
    break;
  default:
    x = -1;
    fprintf(errfile,"ERROR: peg_GetAdj()\n");
  }
  return x;
}

void peg_SetAdj(pExpGraph_t peg, int idx, BOOL vType, int val) {
  switch (vType) {
  case ENTRV:
    *( peg->data + (idx*peg_GetVSize()) + 1) = val;
    break;
  case EXITV:
    *( peg->data + ((peg->entrNum + idx)*peg_GetVSize()) + 1) = val;
    break;
  default:
    fprintf(errfile,"ERROR: peg_SetAdj()\n");
  }
  return;
}

int peg_GetAdjAssn(pExpGraph_t peg, int idx, BOOL vType) {
  int x;
  switch (vType) {
  case ENTRV:
    x = *( peg->data + (idx*peg_GetVSize()) + 2);
    break;
  case EXITV:
    x = *( peg->data + ((peg->entrNum + idx)*peg_GetVSize()) + 2);
    break;
  default:
    x = -1;
    fprintf(errfile,"ERROR: peg_GetAdjAssn()\n");
  }
  return x;
}

void peg_SetAdjAssn(pExpGraph_t peg, int idx, BOOL vType, int val) {
  switch (vType) {
  case ENTRV:
    *( peg->data + (idx*peg_GetVSize()) + 2) = val;
    break;
  case EXITV:
    *( peg->data + ((peg->entrNum + idx)*peg_GetVSize()) + 2) = val;
    break;
  default:
    fprintf(errfile,"ERROR: peg_SetAdjAssn()\n");
  }
  return;
}

int peg_GetArcNum(pExpGraph_t peg, int idx, BOOL vType) {
  int x;
  switch (vType) {
  case ENTRV:
    x = *( peg->data + (idx*peg_GetVSize()) + 3);
    break;
  case EXITV:
    x = *( peg->data + ((peg->entrNum + idx)*peg_GetVSize()) + 3);
    if (x != 0)
      fprintf(errfile,"ERROR: peg_GetArcNum() EXITV shouldn't have arcs.\n");
    break;
  default:
    x = -1;
    fprintf(errfile,"ERROR: peg_GetArcNum()\n");
  }
  return x;
}

void peg_SetArcNum(pExpGraph_t peg, int idx, BOOL vType, int val) {
  switch (vType) {
  case ENTRV:
    *( peg->data + (idx*peg_GetVSize()) + 3) = val;
    break;
  case EXITV:
    *( peg->data + ((peg->entrNum + idx)*peg_GetVSize()) + 3) = val;
    if (val != 0)
      fprintf(errfile,"ERROR: peg_SetArcNum() EXITV shouldn't have arcs.\n");
    break;
  default:
    fprintf(errfile,"ERROR: peg_SetArcNum()\n");
  }
  return;
}

int peg_GetArcOffset(pExpGraph_t peg, int idxEntr) {
  int i, offset;
  offset=0;
  for (i=0 ; i<idxEntr ; i++)
    offset += peg_GetArcNum(peg, i, ENTRV);
  return (peg_GetVSize()*(peg->entrNum + peg->exitNum) + (peg_GetASize()*offset));
}

int peg_GetArcHead(pExpGraph_t peg, int idxEntr, int arcnum) {
  int x;
  x = *( peg->data + peg_GetArcOffset(peg, idxEntr) + arcnum);
  return x;
}

void peg_SetArcHead(pExpGraph_t peg, int idxEntr, int arcnum, int val) {
  *( peg->data + peg_GetArcOffset(peg, idxEntr) + arcnum) = val;
  return;
}

int peg_GetEntrIdx(pExpGraph_t peg, int label, int adj) {
  int i;
  
  for (i=0 ; i<peg->entrNum ; i++) {
    if ((peg_GetLabel(peg, i, ENTRV) == label) &&
	(peg_GetAdj(peg, i, ENTRV) == adj))
      return(i);
  }

  fprintf(errfile,"PE%3d: ERROR: peg_GetEntrIdx()\n",MYNODE);
  return (-1);
}

void peg_Init(pExpGraph_t peg) {
  int i;
  
  peg->data    = (int *)SAFE_MALLOC(peg_GetDataInBytes(peg),
				    "(cycle.c) peg->data");

  peg->newArcs  = (pArc_t)NULL;
  peg->deadEntr = 0;
  peg->deadExit = 0;
  peg->deadArc  = 0;
  peg->newArc   = 0;

  peg->deadEntrMask = (BOOL *)SAFE_MALLOC(peg->entrNum * sizeof(BOOL),
					  "(cycle.c) peg->deadEntrMask");
  for (i=0 ; i<peg->entrNum ; i++)
    *(peg->deadEntrMask + i) = FALSE;

  peg->deadExitMask = (BOOL *)SAFE_MALLOC(peg->exitNum * sizeof(BOOL),
					  "(cycle.c) peg->deadExitMask");
  for (i=0 ; i<peg->exitNum ; i++)
    *(peg->deadExitMask + i) = FALSE;

  return;
}

void peg_SortArcs(pExpGraph_t peg, int idx, BOOL vType) {
  
  switch (vType) {
  case ENTRV:
    qsort(peg->data + peg_GetArcOffset(peg, idx),
	  peg_GetArcNum(peg, idx, ENTRV),
	  sizeof(int), intCompare);
    break;
  case EXITV:
    fprintf(errfile,"ERROR: peg_SortArcs() EXITV shouldn't have arcs.\n");
    break;
  default:
    fprintf(errfile,"ERROR: peg_SortArcs()\n");
  }
  return;
}

void Create_ExpressPack(eGraph_t expGraph, pExpGraph_t *pExpGraph) {
  expVertex_t v;
  expArc_t a;
  pExpGraph_t peg;
  int j, k;

  *pExpGraph = (pExpGraph_t)SAFE_MALLOC(sizeof(struct pExpGraph_s),
					"(cycle.c) *pExpGraph");
  peg = *pExpGraph;

  peg->entrNum = expGraph->entrNum;
  peg->exitNum = expGraph->exitNum;
  peg->arcNum  = getArcNum(expGraph);

  peg_Init(peg);

  v = expGraph->entrV;
  k=0;
  while (v != (expVertex_t)NULL) {
    peg_SetLabel  (peg, k, ENTRV, v->transArc->head);
    peg_SetAdj    (peg, k, ENTRV, v->transArc->tail);
    peg_SetAdjAssn(peg, k, ENTRV, v->transArc->tailAssn);
    peg_SetArcNum (peg, k, ENTRV, v->expArcsNum);
    
    a = v->expArcs;
    j=0;
    while (a != (expArc_t)NULL) {
      peg_SetArcHead(peg, k, j, a->headVertex->transArc->tail);
      j++;
      a = a->next;
    }
    peg_SortArcs(peg, k, ENTRV);
    k++;
    v = v->next;
  }

  v = expGraph->exitV;
  k=0;
  while (v != (expVertex_t)NULL) {
    peg_SetLabel  (peg, k, EXITV, v->transArc->tail);
    peg_SetAdj    (peg, k, EXITV, v->transArc->head);
    peg_SetAdjAssn(peg, k, EXITV, v->transArc->headAssn);
    peg_SetArcNum (peg, k, EXITV, 0);
    k++;
    v = v->next;
  }

  return;
}


void Print_myExpressPack(pExpGraph_t peg) {
  int j, k;
  int arcnum;
  
  fprintf(outfile,"PE%3d: Packed Express Graph with %12d Entr and %12d exit:\n",
	  MYNODE, peg->entrNum, peg->exitNum);
  for (j=0 ; j<peg->entrNum ; j++) {
    arcnum = peg_GetArcNum (peg, j, ENTRV);
    fprintf(outfile,"PE%3d: \tEntrVertex[%12d]: l:%12d adj:%12d(%3d)   Arcs:%12d\n",
	    MYNODE, j,
	    peg_GetLabel  (peg, j, ENTRV),
	    peg_GetAdj    (peg, j, ENTRV),
	    peg_GetAdjAssn(peg, j, ENTRV),
	    arcnum);
    for (k=0 ; k<arcnum ; k++) {
      fprintf(outfile,"PE%3d: \t\tExpressArc[%12d]: to %12d\n",
	      MYNODE, k, peg_GetArcHead(peg, j, k));
    }
  }
  for (j=0 ; j<peg->exitNum ; j++) {
    fprintf(outfile,"PE%3d: \tExitVertex[%12d]: l:%12d adj:%12d(%3d)   Arcs:%12d\n",
	    MYNODE, j,
	    peg_GetLabel  (peg, j, EXITV),
	    peg_GetAdj    (peg, j, EXITV),
	    peg_GetAdjAssn(peg, j, EXITV),
	    peg_GetArcNum (peg, j, EXITV));
  }
  fflush(outfile);
  return;
}

void Print_ExpressPack(pExpGraph_t peg) {
  int i;

  for (i=0 ; i<NODES ; i++) {
    if (MYNODE==i) {
      Print_myExpressPack(peg);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  
  return;
}

void Free_ExpressPack(pExpGraph_t peg) {
  pArc_t a, nexta;

  a = peg->newArcs;
  while (a != (pArc_t)NULL) {
    nexta = a->next;
    free(a);
    a = nexta;
  }

  free(peg->deadEntrMask);
  free(peg->deadExitMask);
  free(peg->data);
  free(peg);
  return;
}

void Send_ExpressPack(pExpGraph_t peg, int toNode) {
  int sizeBuf[3];

  sizeBuf[0] = peg->entrNum;
  sizeBuf[1] = peg->exitNum;
  sizeBuf[2] = peg->arcNum;

  MPI_Send(sizeBuf, 3, MPI_INT, toNode, MPI_TAG_SIZES, MPI_COMM_WORLD);
  MPI_Send(peg->data, peg_GetDataInInts(peg), MPI_INT, toNode, MPI_TAG_DATA,
	   MPI_COMM_WORLD);
  
#ifdef DEBUG_PRINT
  fprintf(outfile,"PE%3d: I SENT THIS TO PE%3d\n",MYNODE,toNode);
  Print_myExpressPack(peg);
  fprintf(outfile,"PE%3d: I SENT THIS **********\n",MYNODE);
  fflush(outfile);
#endif
  return;
}

void Recv_ExpressPack(pExpGraph_t *pExpGraph, int fromNode) {
  pExpGraph_t peg;
  int         sizeBuf[3];
  MPI_Status  mpstat;
  
  *pExpGraph = (pExpGraph_t)SAFE_MALLOC(sizeof(struct pExpGraph_s),
					"(cycle.c) *pExpGraph");
  peg = *pExpGraph;

  MPI_Recv(sizeBuf, 3, MPI_INT, fromNode, MPI_TAG_SIZES, MPI_COMM_WORLD,
	   &mpstat);

  peg->entrNum = sizeBuf[0];
  peg->exitNum = sizeBuf[1];
  peg->arcNum  = sizeBuf[2];

  peg_Init(peg);
  
  MPI_Recv(peg->data, peg_GetDataInInts(peg), MPI_INT, fromNode, MPI_TAG_DATA,
	   MPI_COMM_WORLD, &mpstat);
  
#ifdef DEBUG_PRINT
  fprintf(outfile,"PE%3d: I RECEIVED THIS FROM PE%3d\n",MYNODE,fromNode);
  Print_myExpressPack(peg);
  fprintf(outfile,"PE%3d: I RECEIVED THIS **********\n",MYNODE);
  fflush(outfile);
#endif
  return;
}

int peg_CheckArcHead(pExpGraph_t peg, int idx, int min_idx, int max_idx, int target) {
  int mid_idx;
  int midhead;
  
  mid_idx = min_idx + (max_idx - min_idx) / 2;
  if (max_idx < min_idx)
    return (-1);
  else {
    midhead = peg_GetArcHead(peg, idx, mid_idx);
    if (target < midhead)
      return peg_CheckArcHead(peg, idx, min_idx, mid_idx - 1, target);
    else {
      if (target > midhead)
	return peg_CheckArcHead(peg, idx, mid_idx + 1, max_idx, target);
      else
	return (mid_idx);
    }
  }
}

void peg_GetEntrPred(pExpGraph_t peg, int label, int *num, int **listPred) {
  int i, n;
  int arcnum;
  int foundIdx;
  BOOL found;

  *listPred = (int *)SAFE_MALLOC(peg->entrNum * sizeof(int),
				 "(cycle.c) *listPred");

  n=0;
  for (i=0 ; i<peg->entrNum ; i++) {
    arcnum = peg_GetArcNum(peg, i, ENTRV);
    foundIdx = peg_CheckArcHead(peg, i, 0, arcnum-1, label);
    found = (foundIdx >= 0);
    if (found) {
      *(*listPred + n) = i; 
      n++;
    }
  }
    
  if (n > peg->entrNum)
    fprintf(errfile,"PE%3d: ERROR: peg_GetEntrPred(%3d), n:%d > peg->entrNum: %d\n",
	    MYNODE, label, n, peg->entrNum);
  *num = n;
  
  return;
}

BOOL peg_CheckArc(pExpGraph_t peg, int idx, int arcHeadLabel) {
  int arcnum;
  pArc_t parc;
  int foundIdx;
  BOOL found;

  arcnum = peg_GetArcNum(peg, idx, ENTRV);
  foundIdx = peg_CheckArcHead(peg, idx, 0, arcnum-1, arcHeadLabel);
  found = (foundIdx >= 0);
  
  parc = peg->newArcs;
  while ((parc != (pArc_t)NULL) && (!found)) {
    found = ((parc->idx == idx) && (parc->head == arcHeadLabel));
    parc = parc->next;
  }
  
  return (found);
}

void peg_RemoveVertex(pExpGraph_t peg, int idx, BOOL vType) {
  int i, arcnum, k, lab_exit;
  BOOL found;
#ifdef DEBUG_PRINT
  int lab_entr;
#endif
  
#ifdef DEBUG_PRINT
  fprintf(outfile,"PE%3d: Removing %s vertex idx: %2d label: %2d\n",
	  MYNODE, (vType==ENTRV)?"ENTR":"EXIT",
	  idx, peg_GetLabel(peg,idx,vType));
#endif

  switch (vType) {
  case ENTRV:
    peg->deadEntr++;
    peg->deadEntrMask[idx] = TRUE;
    peg->deadArc += peg_GetArcNum(peg, idx, vType);
#ifdef DEBUG_PRINT
    fprintf(outfile,"PE%3d: Removing %s vertex \t killing %2d arcs\n",
	  MYNODE, (vType==ENTRV)?"ENTR":"EXIT", peg_GetArcNum(peg,idx,vType));
#endif
    break;
  case EXITV:
    peg->deadExit++;
    peg->deadExitMask[idx] = TRUE;
    lab_exit = peg_GetLabel(peg, idx, vType);
    /* CHECK:  */

    /* If we are the last exit vertex with this label, remove express arcs
       aimed at us */

    found = FALSE;
    i=0;
    while ((i<peg->exitNum) && (!found)) {
      found = ((peg_GetLabel(peg, i, EXITV) == lab_exit) &&
	       (peg->deadExitMask[i] == FALSE));
      i++;
    }
    if (!found) {
      for (i=0 ; i<peg->entrNum ; i++) {
	if (peg->deadEntrMask[i] == FALSE) {
	  arcnum = peg_GetArcNum(peg, i, ENTRV);
#ifdef DEBUG_PRINT
	  lab_entr = peg_GetLabel(peg, i, ENTRV);
	  fprintf(outfile,"PE%3d: Removing %s vertex \t checking arcs from %2d\n",
		  MYNODE, (vType==ENTRV)?"ENTR":"EXIT",
		  lab_entr);
#endif
	  k = peg_CheckArcHead(peg, i, 0, arcnum-1, lab_exit);
	  if (k >= 0) {
	      peg_SetArcHead(peg, i, k, -lab_exit);
	      peg->deadArc++;
	  }
	}
      }
    }
    break;
  default:
    fprintf(errfile,"ERROR: peg_RemoveVertex()\n");
  }
  return;
}

void peg_AddArc(pExpGraph_t peg, int idx, int arcHeadLabel) {
  pArc_t parc;

  parc = (pArc_t)SAFE_MALLOC(sizeof(struct pArc_s),
			     "(cycle.c) parc");

  parc->idx  = idx;
  parc->head = arcHeadLabel;
  parc->next = peg->newArcs;

  peg->newArcs = parc;
  peg->newArc++;

  return;
}

void peg_CleanUp(pExpGraph_t peg) {
  /* Merge in the new Express Arcs and delete all vertices with marked labels */
  pArc_t a, nexta;
  pArcFlat_t newArcsArr;
  int k;
  pExpGraph_t pegNew;
  int arcIdx, vIdx, i, nextIdx;
  int totArc, curArc, newArc, idx0, j;
  int aIdx, alab;

  pegNew = (pExpGraph_t)SAFE_MALLOC(sizeof(struct pExpGraph_s),
				    "(cycle.c) pegNew");
  
  pegNew->entrNum = peg->entrNum - peg->deadEntr;
  pegNew->exitNum = peg->exitNum - peg->deadExit;
  pegNew->arcNum  = peg->arcNum  - peg->deadArc + peg->newArc;
#ifdef DEBUG_PRINT
  fprintf(outfile,"PE%3d: dead entr: %3d dead exit: %3d dead arc: %3d new arc: %3d\n",
	  MYNODE, peg->deadEntr, peg->deadExit, peg->deadArc, peg->newArc);
#endif

  peg_Init(pegNew);

  newArcsArr = (pArcFlat_t)SAFE_MALLOC(peg->newArc * sizeof(struct pArcFlat_s),
				       "(cycle.c) newArcsArr");

  k = 0;
  a = peg->newArcs;
  while (a != (pArc_t)NULL) {
    newArcsArr[k].idx  = a->idx;
    newArcsArr[k].head = a->head;
    nexta = a->next;
    free(a);
    a = nexta;
    k++;
  }

  if (k != peg->newArc)
    fprintf(errfile,"PE%3d: ERROR: k (%d) != peg->newArc (%d) \n",
	    MYNODE, k, peg->newArc);

  qsort(newArcsArr, k, sizeof(struct pArcFlat_s), pArcFlatCompare);

#ifdef DEBUG_PRINT
  for (i=0 ; i<k ; i++)
    fprintf(outfile,"PE%3d: newArcsArr[%3d]: (%3d, %3d)\n",
	    MYNODE, i, newArcsArr[i].idx, newArcsArr[i].head);
#endif
  /* CHECK */

  /* Entrance Vertices */
  arcIdx = 0;
  vIdx   = 0;
  
  for (i=0 ; i<peg->entrNum ; i++) {
    if (peg->deadEntrMask[i] == FALSE) {
      peg_SetLabel  (pegNew, vIdx, ENTRV, peg_GetLabel  (peg, i, ENTRV));
      peg_SetAdj    (pegNew, vIdx, ENTRV, peg_GetAdj    (peg, i, ENTRV));
      peg_SetAdjAssn(pegNew, vIdx, ENTRV, peg_GetAdjAssn(peg, i, ENTRV));
      totArc = peg_GetArcNum(peg, i, ENTRV);

      curArc = 0;
      for (j=0 ; j<totArc ; j++) {
	if (peg_GetArcHead(peg, i, j) >= 0)
	  curArc++;
      }
      
      idx0   = arcIdx;
      if (arcIdx < k) {
	nextIdx = newArcsArr[arcIdx].idx;
	if (nextIdx < i)
	  fprintf(errfile,"PE%3d: ERROR: newArcs[arcIdx(%d)].idx (%d) < i (%d)\n",
		  MYNODE, arcIdx, newArcsArr[arcIdx].idx, i);
	while ((arcIdx<k) && (nextIdx == i)) {
	  arcIdx++;
	  if (arcIdx<k)
	    nextIdx =  newArcsArr[arcIdx].idx;
	}
      }
      newArc = arcIdx - idx0;
      peg_SetArcNum(pegNew, vIdx, ENTRV, curArc + newArc);
      aIdx = 0;
      for (j=0 ; j<totArc ; j++) {
	alab = peg_GetArcHead(peg, i, j);
	if (alab >= 0) {
	  peg_SetArcHead(pegNew, vIdx, aIdx, alab);
	  aIdx++;
	}
      }
      if (aIdx != curArc)
	fprintf(errfile,"PE%3d: ERROR: aIdx (%d) != curArc (%d)\n",
		MYNODE, aIdx, curArc);
      for (j=0 ; j<newArc ; j++)
	peg_SetArcHead(pegNew, vIdx, curArc+j, newArcsArr[idx0+j].head);

      peg_SortArcs(pegNew, vIdx, ENTRV);
      
      vIdx++;
    }
  }

  if (vIdx != pegNew->entrNum)
    fprintf(errfile,"PE%3d: ERROR: vIdx (%d) != pegNew->entrNum (%d)\n",
	    MYNODE, vIdx, pegNew->entrNum);
  
  
  /* Exit Vertices */
  vIdx   = 0;

  for (i=0 ; i<peg->exitNum ; i++) {
    if (peg->deadExitMask[i] == FALSE) {
      peg_SetLabel  (pegNew, vIdx, EXITV, peg_GetLabel  (peg, i, EXITV));
      peg_SetAdj    (pegNew, vIdx, EXITV, peg_GetAdj    (peg, i, EXITV));
      peg_SetAdjAssn(pegNew, vIdx, EXITV, peg_GetAdjAssn(peg, i, EXITV));
      peg_SetArcNum (pegNew, vIdx, EXITV, peg_GetArcNum (peg, i, EXITV));
      vIdx++;
    }
  }

  if (vIdx != pegNew->exitNum)
    fprintf(errfile,"PE%3d: ERROR: vIdx (%d) != pegNew->exitNum (%d)\n",
	    MYNODE, vIdx, pegNew->exitNum);
  

  free(newArcsArr);

  /* Swap in new data and parameters */
  free(peg->data);
  peg->entrNum  = pegNew->entrNum;
  peg->exitNum  = pegNew->exitNum;
  peg->arcNum   = pegNew->arcNum;
  peg->data     = pegNew->data;
#if 1
  pegNew->data  = (int *)NULL;
#endif
  peg->newArcs  = pegNew->newArcs;
  peg->deadEntr = pegNew->deadEntr;
  peg->deadExit = pegNew->deadExit;
  peg->deadArc  = pegNew->deadArc;
  peg->newArc   = pegNew->newArc;
  
  peg->deadEntrMask = (BOOL *)SAFE_MALLOC(peg->entrNum * sizeof(BOOL),
					  "(cycle.c) peg->deadEntrMask");
  for (i=0 ; i<peg->entrNum ; i++)
    *(peg->deadEntrMask + i) = FALSE;

  peg->deadExitMask = (BOOL *)SAFE_MALLOC(peg->exitNum * sizeof(BOOL),
					  "(cycle.c) peg->deadExitMask");
  for (i=0 ; i<peg->exitNum ; i++)
    *(peg->deadExitMask + i) = FALSE;
    
  Free_ExpressPack(pegNew);
  return;
}

void ExpressMerge(pExpGraph_t *p0Ptr, pExpGraph_t p1, int h) {
  /* Merge p1 into p0 */
  pExpGraph_t pTemp;
  pExpGraph_t p0;
  int *ptr, *loc0, *loc1;
  int siz;
  int i, j, k;
  int curAdjAssn, origAdjAssn;
  int exit0, exit0Idx, entr1, entr1Idx, exit1;
  int entr1ArcNum;
  int p, s;

  int entr0Num;
  int *entr0List;

  int exit1Num;
  int *exit1List;
  int exit1Lab;

  p0 = *p0Ptr;
  
#ifdef DEBUG_PRINT
  fprintf(outfile,"PE%3d: MY P0  MERGE AT THE START OF THINGS\n",MYNODE);
  Print_myExpressPack(p0);
  fprintf(outfile,"PE%3d: MY P0  MERGE AT THE START OF THINGS DONE\n",MYNODE);
  fflush(outfile);
#endif

  pTemp = (pExpGraph_t)SAFE_MALLOC(sizeof(struct pExpGraph_s),
				   "(cycle.c) pTemp");
  
  pTemp->entrNum = p0->entrNum + p1->entrNum;
  pTemp->exitNum = p0->exitNum + p1->exitNum;
  pTemp->arcNum  = p0->arcNum  + p1->arcNum;

  peg_Init(pTemp);

  ptr  = pTemp->data;

  siz  = p0->entrNum*peg_GetVSize();
  loc0 = p0->data;
  memcpy(ptr, loc0, siz*sizeof(int));
  ptr += siz;

  siz  = p1->entrNum*peg_GetVSize();
  loc1 = p1->data;
  memcpy(ptr, loc1, siz*sizeof(int));
  ptr += siz;

  siz   = p0->exitNum*peg_GetVSize();
  loc0 += p0->entrNum*peg_GetVSize();
  memcpy(ptr, loc0 , siz*sizeof(int));
  ptr += siz;

  siz   = p1->exitNum*peg_GetVSize();
  loc1 += p1->entrNum*peg_GetVSize();
  memcpy(ptr, loc1, siz*sizeof(int));
  ptr += siz;
  
  siz   = p0->arcNum*peg_GetASize();
  loc0 += p0->exitNum*peg_GetVSize();
  memcpy(ptr, loc0 , siz*sizeof(int));
  ptr += siz;

  siz   = p1->arcNum*peg_GetASize();
  loc1 += p1->exitNum*peg_GetVSize();
  memcpy(ptr, loc1, siz*sizeof(int));
  ptr += siz;

#ifdef DEBUG_PRINT
  fprintf(outfile,"PE%3d: MY P0  MERGE\n",MYNODE);
  Print_myExpressPack(p0);
  fprintf(outfile,"PE%3d: MY P0  MERGE DONE\n",MYNODE);
  fflush(outfile);
  fprintf(outfile,"PE%3d: MY P1  MERGE\n",MYNODE);
  Print_myExpressPack(p1);
  fprintf(outfile,"PE%3d: MY P1  MERGE DONE\n",MYNODE);
  fflush(outfile);

  fprintf(outfile,"PE%3d: MY P0+P1    MERGE\n",MYNODE);
  Print_myExpressPack(pTemp);
  fprintf(outfile,"PE%3d: MY P0+P1    MERGE DONE\n",MYNODE);
  fflush(outfile);
#endif

  /****** MERGE HERE ********/

  for (i=0 ; i<pTemp->exitNum ; i++) {
    exit0Idx = i;
    origAdjAssn = peg_GetAdjAssn(pTemp, exit0Idx, EXITV);
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
      exit0 = peg_GetLabel(pTemp, exit0Idx, EXITV);

      /* entr1 = label of the entrance vertex we point to */
      entr1     = peg_GetAdj    (pTemp, exit0Idx, EXITV);

      /* entr1Idx,ENTRV = index of the entrance vertex */
      entr1Idx  = peg_GetEntrIdx(pTemp, entr1, exit0);
#ifdef DEBUG_PRINT
      fprintf(outfile,"PE%3d: h: %d orig: %d cur: %d exit0: %3d entr1: %3d entr1Idx: %3d\n",
	    MYNODE, h, origAdjAssn, curAdjAssn, exit0, entr1, entr1Idx);
#endif

      /* If the express arc exists from entr1 to exit0, we have a cycle! */

      if (peg_CheckArc(pTemp, entr1Idx, exit0)) {
	fprintf(outfile,"PE%3d: HALT: CYCLE DETECTED DURING MERGE (%3d, %3d)\n",
		MYNODE, exit0, entr1);
	CYCLE_FOUND = TRUE;
	return;
      }
      
      /* Get list of the entrance vertices that have express arcs to exit0 */
      peg_GetEntrPred(pTemp, exit0, &entr0Num, &entr0List);

      if (entr0Num > 0) {
	/* For each express arc of entr1, replace it with the new express from
	   entr0List to exit1 */
	entr1ArcNum = peg_GetArcNum(pTemp, entr1Idx, ENTRV);
	if (entr1ArcNum > 0) {
#ifdef DEBUG_PRINT
	  fprintf(outfile,"PE%3d: entr1Idx: %2d [%2d : %2d] entr1ArcNum: %2d\n",
		  MYNODE, entr1Idx, peg_GetLabel(pTemp, entr1Idx, ENTRV),
		  entr1, entr1ArcNum);
#endif
      
	  exit1List = (int *)SAFE_MALLOC(pTemp->exitNum * sizeof(int),
					 "(cycle.c) exit1List");
	  exit1Num = 0;
      
	  for (k=0 ; k<entr1ArcNum ; k++) {
	    exit1 = peg_GetArcHead(pTemp, entr1Idx, k);
	    if (exit1 >= 0) {
	      /* Add Express Arc from entr0 to exit1 */
	      /* There could be multiple exit1's with the same label */
	      for (j=0 ; j<pTemp->exitNum ; j++) {
		if (peg_GetLabel(pTemp, j, EXITV) == exit1) {
		  exit1List[exit1Num] = j;
		  exit1Num++;
		}
	      }
	    }
	  }

	  if (exit1Num > pTemp->exitNum) {
	    fprintf(errfile,"PE%3d: ERROR: exit1Num (%d) > pTemp->exitNum (%d)\n",
		    MYNODE, exit1Num, pTemp->exitNum);

	  }

	  if (exit1Num > 0) {
	    /* We have entr0 and exit1 vertices to link */
#ifdef DEBUG_PRINT
	    for (k=0 ; k < entr0Num ; k++) {
	      fprintf(outfile,
		      "PE%3d: \t (%2d) exit0: %3d n:%3d  entr0Idx: %3d  label: %3d\n",
		      MYNODE, k, exit0, entr0Num, entr0List[k],
		      peg_GetLabel(pTemp, entr0List[k], ENTRV));
	      fflush(outfile);
	    }
	    for (k=0 ; k < exit1Num ; k++) {
	      fprintf(outfile,
		      "PE%3d: \t (%2d) exit0: %3d n:%3d  exit1Idx: %3d  label: %3d\n",
		      MYNODE, k, exit0, exit1Num, exit1List[k],
		      peg_GetLabel(pTemp, exit1List[k], EXITV));
	      fflush(outfile);
	    }
#endif

	    for (p=0 ; p<entr0Num ; p++) {
	      for (s=0 ; s<exit1Num ; s++) {
#ifdef DEBUG_PRINT
		fprintf(outfile,"PE%3d: Testing express arc (%3d, %3d) idx %3d -> %3d\n",
			MYNODE,
			peg_GetLabel(pTemp, entr0List[p], ENTRV),
			peg_GetLabel(pTemp, exit1List[s], EXITV),
			entr0List[p], exit1List[s]);
		fflush(outfile);
#endif
		/* If the arc doesn't exit, add it */
		exit1Lab = peg_GetLabel(pTemp, exit1List[s], EXITV);
		if (!peg_CheckArc(pTemp, entr0List[p], exit1Lab)) {
#ifdef DEBUG_PRINT
		  fprintf(outfile,"PE%3d: Adding  express arc (%3d, %3d) idx %3d -> %3d\n",
			MYNODE,
			peg_GetLabel(pTemp, entr0List[p], ENTRV),
			peg_GetLabel(pTemp, exit1List[s], EXITV),
			entr0List[p], exit1List[s]);
		  fflush(outfile);
#endif
		  peg_AddArc(pTemp, entr0List[p], exit1Lab);
		}
#ifdef DEBUG_PRINT
		else {
		  fprintf(outfile,"PE%3d: Tossing express arc (%3d, %3d) idx %3d -> %3d\n",
			MYNODE,
			peg_GetLabel(pTemp, entr0List[p], ENTRV),
			peg_GetLabel(pTemp, exit1List[s], EXITV),
			entr0List[p], exit1List[s]);
		  fflush(outfile);
		}
#endif
		
	      } /* foreach exit1 */
	    } /* foreach entr0 */

	  } /* (exit1Num > 0) */
	  
	  free(exit1List);
	  
	} /* (entr1ArcNum > 0) */
	
      } /* (entr0Num > 0) */
      
      free(entr0List);
      
      peg_RemoveVertex(pTemp, exit0Idx, EXITV);
      peg_RemoveVertex(pTemp, entr1Idx, ENTRV);
	    
    } /* (curAdjAssn  == MYNODE) */
    
  } /* foreach exit0 */

  peg_CleanUp(pTemp);

#ifdef DEBUG_PRINT
  Print_myExpressPack(pTemp);
#endif
  
  *p0Ptr = pTemp;

  Free_ExpressPack(p0);
  return;
}

void Merge_ExpressPack(pExpGraph_t *peg) {
  int h;
  int logp;
  pExpGraph_t pExpGraph_recv;
  int found, result;

  /* logp is the ceiling of log2(NODES) */
  logp = log2_i(NODES);
#ifdef DEBUG_PRINT
  MPI_fprintf(outfile,"NODES: %d logp: %d\n",NODES,logp);
#endif
  for (h=0 ; h<logp ; h++) {
#if PRINT_MERGESIZE
    MPI_Barrier(MPI_COMM_WORLD);
#endif
#ifdef DEBUG_PRINT
    MPI_fprintf(outfile,"lastB(%3d): %3d testB(%3d): %3d \t setB(%3d): %3d clearB(%3d): %3d\n",h,lastB(MYNODE,h),
		h,testB(MYNODE,h),h,setB(MYNODE,h),h,clearB(MYNODE,h));
#endif
    if (lastB(MYNODE, h)==0) {
#if PRINT_MERGESIZE
      fprintf(outfile,"PE%3d: [%3d] entr: %12d exit: %12d exarcs: %12d\n",
	      MYNODE, h, (*peg)->entrNum, (*peg)->exitNum, (*peg)->arcNum);
#endif
      if (testB(MYNODE,h)==0) {
	/* Make sure that if we're not using a power-of-two nodes
	   that we just sit idle for a stage where our partner is missing */
	if (setB(MYNODE,h) < NODES) {
#ifdef DEBUG_PRINT
	  fprintf(outfile,"PE%3d: BEFORE MERGE\n",MYNODE);
	  Print_myExpressPack(*peg);
	  fprintf(outfile,"PE%3d: BEFORE MERGE DONE\n",MYNODE);
	  fflush(outfile);
#endif
#ifdef DEBUG_PRINT
	  fprintf(outfile,"PE%3d: (h:%3d) I'm receiving from PE%3d\n",
		  MYNODE, h, setB(MYNODE,h));
	  fflush(outfile);
#endif
	  Recv_ExpressPack(&pExpGraph_recv, setB(MYNODE,h));
	  ExpressMerge(peg, pExpGraph_recv, h);

#ifdef DEBUG_PRINT
	  fprintf(outfile,"PE%3d: AFTER MERGE\n",MYNODE);
	  Print_myExpressPack(*peg);
	  fprintf(outfile,"PE%3d: AFTER MERGE DONE\n",MYNODE);
	  fflush(outfile);
#endif
	  Free_ExpressPack(pExpGraph_recv);
#ifdef DEBUG_PRINT
	  fprintf(outfile,"PE%3d: AFTER FREEING EXPRESSPACK MERGE\n",MYNODE);
	  Print_myExpressPack(*peg);
	  fprintf(outfile,"PE%3d: AFTER FREEING EXPRESSPACK MERGE DONE\n",MYNODE);
	  fflush(outfile);
#endif
	}
      }
      else {
#ifdef DEBUG_PRINT
	fprintf(outfile,"PE%3d: (h:%3d) I'm sending to PE%3d\n",
		MYNODE, h, clearB(MYNODE,h));
	fflush(outfile);
#endif
	Send_ExpressPack(*peg, clearB(MYNODE,h));
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
   
int pig_GetVSize() {
  /* label, adj, adjassn, icount */
  return 4;
}

int pig_GetISize() {
  return 2;
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

int pig_GetIntervalCount(pIntervalGraph_t pig, int idx, BOOL vType) {
  int x;
  switch (vType) {
  case ENTRV:
    x = *( pig->data + (idx*pig_GetVSize()) + 3);
    break;
  case EXITV:
    x = *( pig->data + ((pig->entrNum + idx)*pig_GetVSize()) + 3);
    if (x != 0)
      fprintf(errfile,"ERROR: pig_GetIntervalCount() EXITV.\n");
    break;
  default:
    x = -1;
    fprintf(errfile,"ERROR: pig_GetIntervalCount()\n");
  }
  return x;
}

void pig_SetIntervalCount(pIntervalGraph_t pig, int idx, BOOL vType, int val) {
  switch (vType) {
  case ENTRV:
    *( pig->data + (idx*pig_GetVSize()) + 3) = val;
    break;
  case EXITV:
    *( pig->data + ((pig->entrNum + idx)*pig_GetVSize()) + 3) = val;
    if (val != 0)
      fprintf(errfile,"ERROR: pig_SetIntervalCount() EXITV.\n");
    break;
  default:
    fprintf(errfile,"ERROR: pig_SetIntervalCount()\n");
  }
  return;
}

void pig_SetIntervalOffsets(pIntervalGraph_t pig) {
  int i, *lastPtr;
  pig->intList = (int **)SAFE_MALLOC(pig->entrNum*sizeof(int *),
				     "(cycleConvex.c) pig->intList");

  lastPtr = pig->data + (pig->entrNum + pig->exitNum)*pig_GetVSize();
  
  for (i=0 ; i<pig->entrNum ; i++) {
    pig->intList[i] = lastPtr;
    lastPtr += pig_GetIntervalCount(pig, i, ENTRV)*pig_GetISize();
  }
  return;
}

int pig_GetIntervalC0(pIntervalGraph_t pig, int idxEntr, int intnum) {
  int x;
  if (intnum > pig_GetIntervalCount(pig, idxEntr, ENTRV))
    fprintf(errfile,"PE%3d: ERROR: pig_GetIntervalC0 idx: %d intnum: %d\n",
	    MYNODE, idxEntr, intnum);
  x = *(pig->intList[idxEntr] + 2*intnum);
  return x;
}

int pig_GetIntervalC1(pIntervalGraph_t pig, int idxEntr, int intnum) {
  int x;
  if (intnum > pig_GetIntervalCount(pig, idxEntr, ENTRV))
    fprintf(errfile,"PE%3d: ERROR: pig_GetIntervalC0 idx: %d intnum: %d\n",
	    MYNODE, idxEntr, intnum);
  x = *(pig->intList[idxEntr] + 2*intnum + 1);
  return x;
}

void pig_SetIntervalC0(pIntervalGraph_t pig, int idxEntr, int intnum, int val) {
  if (intnum > pig_GetIntervalCount(pig, idxEntr, ENTRV))
    fprintf(errfile,"PE%3d: ERROR: pig_SetIntervalC0 idx: %d intnum: %d\n",
	    MYNODE, idxEntr, intnum);
  *(pig->intList[idxEntr] + 2*intnum) = val;
  return;
}

void pig_SetIntervalC1(pIntervalGraph_t pig, int idxEntr, int intnum, int val) {
  if (intnum > pig_GetIntervalCount(pig, idxEntr, ENTRV))
    fprintf(errfile,"PE%3d: ERROR: pig_SetIntervalC0 idx: %d intnum: %d\n",
	    MYNODE, idxEntr, intnum);
  *(pig->intList[idxEntr] + 2*intnum + 1) = val;
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

  pig->newInt = 0;
  pig->deadInt = 0;
  pig->newIntervals = (pInterval_t)NULL;
  pig->intList = (int **)NULL;
  
  return;
}

void pig_SortIntervals(pIntervalGraph_t pig, int idx, BOOL vType) {
  
  switch (vType) {
  case ENTRV:
    qsort(pig->intList[idx], 
	  pig_GetIntervalCount(pig, idx, ENTRV),
	  pig_GetISize()*sizeof(int), intCompare);
    break;
  case EXITV:
    fprintf(errfile,"ERROR: pig_SortIntervals() EXITV.\n");
    break;
  default:
    fprintf(errfile,"ERROR: pig_SortIntervals()\n");
  }
  return;
}

void Create_IntervalPack(eGraph_t expGraph, pIntervalGraph_t *pIntervalGraph) {
  expVertex_t v;
  expArc_t a;
  pIntervalGraph_t pig;
  int k;
  int val, minval, maxval;

  *pIntervalGraph = (pIntervalGraph_t)SAFE_MALLOC(sizeof(struct pIntervalGraph_s),
					"(cycle.c) *pIntervalGraph");
  pig = *pIntervalGraph;

  pig->entrNum      = expGraph->entrNum;
  pig->exitNum      = expGraph->exitNum;
  pig->intervalNum  = pig->entrNum;

  pig_Init(pig);

  for (k=0 ; k<pig->entrNum ; k++)
    pig_SetIntervalCount(pig, k, ENTRV, 1);

  pig_SetIntervalOffsets(pig);

  v = expGraph->entrV;
  k=0;
  while (v != (expVertex_t)NULL) {
    pig_SetLabel  (pig, k, ENTRV, v->transArc->head);
    pig_SetAdj    (pig, k, ENTRV, v->transArc->tail);
    pig_SetAdjAssn(pig, k, ENTRV, v->transArc->tailAssn);
    
    a = v->expArcs;
    minval = maxval = -1;
    if (a != (expArc_t)NULL) {
      maxval = minval = a->headVertex->transArc->tail;
      a = a->next;
    }
    while (a != (expArc_t)NULL) {
      val = a->headVertex->transArc->tail;
      minval = min(minval, val);
      maxval = max(maxval, val);
      a = a->next;
    }
    if (maxval >= 0) {
      pig_SetIntervalC0(pig, k, 0, minval);
      pig_SetIntervalC1(pig, k, 0, maxval);
    }
#if 0
    pig_SortIntervals(pig, k, ENTRV);
#endif
    k++;
    v = v->next;
  }

  v = expGraph->exitV;
  k=0;
  while (v != (expVertex_t)NULL) {
    pig_SetLabel  (pig, k, EXITV, v->transArc->tail);
    pig_SetAdj    (pig, k, EXITV, v->transArc->head);
    pig_SetAdjAssn(pig, k, EXITV, v->transArc->headAssn);
    pig_SetIntervalCount(pig, k, EXITV, 0);
    k++;
    v = v->next;
  }

  return;
}


void Print_myIntervalPack(pIntervalGraph_t pig) {
  int j, k;
  int intervals;
  
  fprintf(outfile,"PE%3d: Packed Interval Graph with %12d Entr and %12d exit:\n",
	  MYNODE, pig->entrNum, pig->exitNum);
  for (j=0 ; j<pig->entrNum ; j++) {
    intervals = pig_GetIntervalCount (pig, j, ENTRV);
    fprintf(outfile,"PE%3d: \tEntrVertex[%12d]: l:%12d adj:%12d(%3d)   Ints:%12d\n",
	    MYNODE, j,
	    pig_GetLabel  (pig, j, ENTRV),
	    pig_GetAdj    (pig, j, ENTRV),
	    pig_GetAdjAssn(pig, j, ENTRV),
	    intervals);
    for (k=0 ; k<intervals ; k++) {
      fprintf(outfile,"PE%3d: \t\tInterval[%12d]: [%12d, %12d]\n",
	      MYNODE, k,
	      pig_GetIntervalC0(pig, j, k), pig_GetIntervalC1(pig, j, k));
    }
  }
  for (j=0 ; j<pig->exitNum ; j++) {
    fprintf(outfile,"PE%3d: \tExitVertex[%12d]: l:%12d adj:%12d(%3d)   Ints:%12d\n",
	    MYNODE, j,
	    pig_GetLabel  (pig, j, EXITV),
	    pig_GetAdj    (pig, j, EXITV),
	    pig_GetAdjAssn(pig, j, EXITV),
	    pig_GetIntervalCount(pig, j, EXITV));
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
  pInterval_t i, nexti;

  i = pig->newIntervals;
  while (i != (pInterval_t)NULL) {
    nexti = i->next;
    free(i);
    i = nexti;
  }
  
  free(pig->deadEntrMask);
  free(pig->deadExitMask);
  if (pig->intList != (int **)NULL)
    free(pig->intList);
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

  pig_SetIntervalOffsets(pig);
  
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
  int intervals;
  int s;
  BOOL found;

  s = 0;
  intervals = pig_GetIntervalCount(pig, idx, ENTRV);
  found = FALSE;
  while (!found && (s < intervals)) {
    minval = pig_GetIntervalC0(pig, idx, s);
    maxval = pig_GetIntervalC1(pig, idx, s);
    found = ((minval <= target) && (target <= maxval));
    s++;
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

BOOL pig_CheckArc(pIntervalGraph_t pig, int idx, int arcHeadLabel) {
  pInterval_t pint;
  BOOL found;

  found = pig_CheckArcHead(pig, idx, arcHeadLabel);
  
  pint = pig->newIntervals;
  while ((pint != (pInterval_t)NULL) && (!found)) {
    found = ((pint->idx == idx) &&
	     (pint->C0 <= arcHeadLabel) && (arcHeadLabel <= pint->C1));
    pint = pint->next;
  }

  return (found);
}

BOOL pig_CheckInterval(pIntervalGraph_t pig, int idx, int C0, int C1) {
  int s, intervals;
  pInterval_t pint;
  BOOL found;

  intervals = pig_GetIntervalCount(pig, idx, ENTRV);
  found = FALSE;
  s = 0;
  while ((!found) && (s<intervals)) {
    found = ((pig_GetIntervalC0(pig, idx, s) == C0) &&
	     (pig_GetIntervalC1(pig, idx, s) == C1));
    s++;
  }
  
  pint = pig->newIntervals;
  while ((pint != (pInterval_t)NULL) && (!found)) {
    found = ((pint->idx == idx) &&
	     (pint->C0 == C0) && (pint->C1 == C1));
    pint = pint->next;
  }

  return (found);
}

void pig_RemoveVertex(pIntervalGraph_t pig, int idx, BOOL vType) {
  int i, intervals, k, lab_exit;
  int minval, maxval;
  BOOL found;
#ifdef DEBUG_PRINT
  int lab_entr;
#endif
  
#ifdef DEBUG_PRINT
  fprintf(outfile,"PE%3d: Removing %s vertex idx: %2d label: %2d\n",
	  MYNODE, (vType==ENTRV)?"ENTR":"EXIT",
	  idx, pig_GetLabel(pig,idx,vType));
#endif

  switch (vType) {
  case ENTRV:
    pig->deadEntr++;
    pig->deadEntrMask[idx] = TRUE;
    pig->deadInt += pig_GetIntervalCount(pig, idx, vType);
#ifdef DEBUG_PRINT
    fprintf(outfile,"PE%3d: Removing %s vertex \t killing %2d intervals\n",
	    MYNODE, (vType==ENTRV)?"ENTR":"EXIT",
	    pig_GetIntervalCount(pig,idx,vType));
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
      for (i=0 ; i<pig->entrNum ; i++) {
	if (pig->deadEntrMask[i] == FALSE) {
#ifdef DEBUG_PRINT
	  lab_entr = pig_GetLabel(pig, i, ENTRV);
	  fprintf(outfile,
		  "PE%3d: Removing %s vertex \t checking intervals from %2d\n",
		  MYNODE, (vType==ENTRV)?"ENTR":"EXIT",
		  lab_entr);
#endif
	  intervals = pig_GetIntervalCount(pig, i, ENTRV);
	  for (k=0 ; k<intervals ; k++) {
	    minval = pig_GetIntervalC0(pig, i, k);
	    maxval = pig_GetIntervalC1(pig, i, k);
	    if ((minval <= lab_exit) && (lab_exit <= maxval)) {
	      if (minval == maxval) {
		/* Remove this interval */
		pig_SetIntervalC0(pig, i, k, -1);
		pig_SetIntervalC1(pig, i, k, -1);
		pig->deadInt++;
	      }
	      else {
		if (minval == lab_exit)
		  pig_SetIntervalC0(pig, i, k, minval+1);
		else {
		  if (maxval == lab_exit)
		    pig_SetIntervalC1(pig, i, k, maxval-1);
		}
	      }
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

void pig_AddInterval(pIntervalGraph_t pig, int idx, int C0, int C1) {
  pInterval_t pint;

  pint = (pInterval_t)SAFE_MALLOC(sizeof(struct pInterval_s),
				  "(cycle.c) pint");

  pint->idx  = idx;
  pint->C0   = C0;
  pint->C1   = C1;
  pint->next = pig->newIntervals;

  pig->newIntervals = pint;
  pig->newInt++;

  return;
}

void pig_CleanUp(pIntervalGraph_t pig) {
  /* Merge in the new Intervals and delete all vertices with marked labels */
  pInterval_t iv, nextiv;
  pIntervalFlat_t newIntsArr;
  int k;
  pIntervalGraph_t pigNew;
  int intIdx, vIdx, i, nextIdx;
  int totInt, curInt, newInt, idx0, j;
  int iIdx, minlab;
  int *intPtr, *intPtr0;

  pigNew = (pIntervalGraph_t)SAFE_MALLOC(sizeof(struct pIntervalGraph_s),
					 "(cycle.c) pigNew");
  
  pigNew->entrNum = pig->entrNum - pig->deadEntr;
  pigNew->exitNum = pig->exitNum - pig->deadExit;
  pigNew->intervalNum = pig->intervalNum - pig->deadInt + pig->newInt;
#ifdef DEBUG_PRINT
  fprintf(outfile,"PE%3d: dead entr: %3d dead exit: %3d dead int: %3d new int: %3d\n",
	  MYNODE, pig->deadEntr, pig->deadExit, pig->deadInt, pig->newInt);
#endif

  pig_Init(pigNew);

  newIntsArr = (pIntervalFlat_t)
    SAFE_MALLOC(pig->newInt * sizeof(struct pIntervalFlat_s),
		"(cycle.c) newIntsArr");

  k = 0;
  iv = pig->newIntervals;
  while (iv != (pInterval_t)NULL) {
    newIntsArr[k].idx = iv->idx;
    newIntsArr[k].C0  = iv->C0;
    newIntsArr[k].C1  = iv->C1;
    nextiv = iv->next;
    free(iv);
    iv = nextiv;
    k++;
  }

  if (k != pig->newInt)
    fprintf(errfile,"PE%3d: ERROR: k (%d) != pig->newInt (%d) \n",
	    MYNODE, k, pig->newInt);

  qsort(newIntsArr, k, sizeof(struct pIntervalFlat_s), pIntervalFlatCompare);

#ifdef DEBUG_PRINT
  for (i=0 ; i<k ; i++)
    fprintf(outfile,"PE%3d: newIntsArr[%3d]: (%3d, [%3d %3d])\n",
	    MYNODE, i, newIntsArr[i].idx, newIntsArr[i].C0, newIntsArr[i].C1);
#endif
  /* CHECK */

  /* Entrance Vertices */
  intIdx = 0;
  vIdx = 0;
  intPtr = pigNew->data + (pigNew->entrNum + pigNew->exitNum)*pig_GetVSize();
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

      totInt = pig_GetIntervalCount(pig, i, ENTRV);

      curInt = 0;
      for (j=0 ; j<totInt ; j++) {
	if (pig_GetIntervalC0(pig, i, j) >= 0)
	  curInt++;
      }
      
      idx0   = intIdx;
      if (intIdx < k) {
	nextIdx = newIntsArr[intIdx].idx;
	if (nextIdx < i)
	  fprintf(errfile,"PE%3d: ERROR: newInts[intIdx(%d)].idx (%d) < i (%d)\n",
		  MYNODE, intIdx, newIntsArr[intIdx].idx, i);
	while ((intIdx<k) && (nextIdx == i)) {
	  intIdx++;
	  if (intIdx<k)
	    nextIdx =  newIntsArr[intIdx].idx;
	}
      }
      newInt = intIdx - idx0;
      pig_SetIntervalCount(pigNew, vIdx, ENTRV, curInt + newInt);
      iIdx = 0;
      intPtr0 = intPtr;
      for (j=0 ; j<totInt ; j++) {
	minlab = pig_GetIntervalC0(pig, i, j);
	if (minlab >= 0) {
	  *intPtr = minlab;
	  intPtr++;
	  *intPtr = pig_GetIntervalC1(pig, i, j);
	  intPtr++;
	  iIdx++;
	}
      }
      if (iIdx != curInt)
	fprintf(errfile,"PE%3d: ERROR: iIdx (%d) != curArc (%d)\n",
		MYNODE, iIdx, curInt);
      for (j=0 ; j<newInt ; j++) {
	*intPtr = newIntsArr[idx0+j].C0;
	intPtr++;
	*intPtr = newIntsArr[idx0+j].C1;
	intPtr++;
      }

      if (curInt + newInt != (int)(intPtr - intPtr0) / 2)
	fprintf(errfile,
		"PE%3d: ERROR: intPtr (cur: %d new: %d [tot:%d]) (%p - %p = %d)\n",
		MYNODE, curInt, newInt, curInt+newInt, intPtr, intPtr0,
		(int)(intPtr - intPtr0)/2);

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
      pig_SetIntervalCount(pigNew, vIdx, EXITV,
			   pig_GetIntervalCount(pig, i, EXITV));
      vIdx++;
    }
  }

  if (vIdx != pigNew->exitNum)
    fprintf(errfile,"PE%3d: ERROR: vIdx (%d) != pigNew->exitNum (%d)\n",
	    MYNODE, vIdx, pigNew->exitNum);
  
#ifdef DEBUG_PRINT
  fprintf(outfile,"PE%3d: cu exit done (num: %d) \n",MYNODE, pig->exitNum);
#endif

  free(newIntsArr);

#ifdef DEBUG_PRINT
  fprintf(outfile,"PE%3d: cu Freed newIntsArr\n",MYNODE);
#endif
  /* Swap in new data and parameters */
  free(pig->data);
  pig->entrNum  = pigNew->entrNum;
  pig->exitNum  = pigNew->exitNum;
  pig->intervalNum   = pigNew->intervalNum;
  pig->data     = pigNew->data;
#if 1
  pigNew->data  = (int *)NULL;
#endif
  pig->newIntervals  = pigNew->newIntervals;
  pig->deadEntr = pigNew->deadEntr;
  pig->deadExit = pigNew->deadExit;
  pig->deadInt  = pigNew->deadInt;
  pig->newInt   = pigNew->newInt;
  
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

  free(pig->intList);
  pig_SetIntervalOffsets(pig);

  for (i=0 ; i<pig->entrNum ; i++)
    pig_SortIntervals(pig, i, ENTRV);
    
#ifdef DEBUG_PRINT
  fprintf(outfile,"PE%3d: cu set new interval offsets\n",MYNODE);
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
  int entr1IntervalNum;
  int p, s;
  int minval, maxval;
    
  int entr0Num;
  int *entr0List;

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

  pig_SetIntervalOffsets(pTemp);
  
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

      if (pig_CheckArc(pTemp, entr1Idx, exit0)) {
	fprintf(outfile,"PE%3d: HALT: CYCLE DETECTED DURING MERGE (%3d, %3d)\n",
		MYNODE, exit0, entr1);
	CYCLE_FOUND = TRUE;
	return;
      }
      
      /* Get list of the entrance vertices that have Interval arcs to exit0 */
      pig_GetEntrPred(pTemp, exit0, &entr0Num, &entr0List);

      if (entr0Num > 0) {
	/* For each Interval of entr1, replace it with the new Interval from
	   entr0List to exit1 */
	entr1IntervalNum = pig_GetIntervalCount(pTemp, entr1Idx, ENTRV);
	if (entr1IntervalNum > 0) {
#ifdef DEBUG_PRINT
	  fprintf(outfile,"PE%3d: entr1Idx: %2d [%2d : %2d] entr1IntervalNum: %2d\n",
		  MYNODE, entr1Idx, pig_GetLabel(pTemp, entr1Idx, ENTRV),
		  entr1, entr1IntervalNum);
#endif

	  for (p=0 ; p<entr0Num ; p++) {
	    for (s=0 ; s<entr1IntervalNum ; s++) {
	      minval = pig_GetIntervalC0(pTemp, entr1Idx, s);
	      maxval = pig_GetIntervalC1(pTemp, entr1Idx, s);
	      if (!pig_CheckInterval(pTemp, entr0List[p], minval, maxval)) {
#ifdef DEBUG_PRINT
		fprintf(outfile,"PE%3d: Adding  Interval (%3d, [%3d %3d]) idx %3d\n",
			MYNODE,
			pig_GetLabel(pTemp, entr0List[p], ENTRV),
			minval, maxval,
			entr0List[p]);
		fflush(outfile);
#endif
		pig_AddInterval(pTemp, entr0List[p], minval, maxval);
	      }
#ifdef DEBUG_PRINT
	      else {
		fprintf(outfile,"PE%3d: Tossing Interval (%3d, [%3d %3d]) idx %3d\n",
			MYNODE,
			pig_GetLabel(pTemp, entr0List[p], ENTRV),
			minval, maxval,
			entr0List[p]);
		fflush(outfile);
	      }
#endif
		
	    } /* foreach exit1 */
	  } /* foreach entr0 */
	  
	} /* (entr1ArcNum > 0) */
	
      } /* (entr0Num > 0) */
      
      free(entr0List);
      
      pig_RemoveVertex(pTemp, exit0Idx, EXITV);
      pig_RemoveVertex(pTemp, entr1Idx, ENTRV);
	    
    } /* (curAdjAssn  == MYNODE) */
    
  } /* foreach exit0 */

  pig_CleanUp(pTemp);

#ifdef DEBUG_PRINT
  Print_myIntervalPack(pTemp);
#endif
  
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
#if USE_PEG
  pExpGraph_t pExpGraph;
#else
  pIntervalGraph_t pIntervalGraph;
#endif
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

  Find_Local_Cycles(myVerts);
  timer_mark("Find Local Cycles");

  found = (CYCLE_FOUND==TRUE)?1:0;
  MPI_Allreduce(&found, &result, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  if (!result) {
  
    Discovery(myVerts,
	      &initTransArcCount, 
	      &initTransArcs,
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
    
#if USE_PEG
    Create_ExpressPack(expGraph, &pExpGraph);
    timer_mark("Create_ExpressPack");

#ifdef DEBUG_PRINT
    Print_ExpressPack(pExpGraph);
    timer_mark("Print_ExpressPack");
#endif
  
    Merge_ExpressPack(&pExpGraph);
    timer_mark("Merge_ExpressPack");

    Free_ExpressPack(pExpGraph);
#else /* USE_PIG */
    Create_IntervalPack(expGraph, &pIntervalGraph);
    timer_mark("Create_IntervalPack");

#ifdef DEBUG_PRINT
    Print_IntervalPack(pIntervalGraph);
    timer_mark("Print_IntervalPack");
#endif
  
    Merge_IntervalPack(&pIntervalGraph);
    timer_mark("Merge_IntervalPack");
    
    Free_IntervalPack(pIntervalGraph);
#endif
    
    Free_Express(expGraph);
    free(termTransArcs);
    free(initTransArcs);
  }

  Report_Cycles();
  timer_report(outfile,"CYCLE DETECTION", n);

  MPI_Barrier(MPI_COMM_WORLD);

  Free_Input(myVerts);

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  
}
