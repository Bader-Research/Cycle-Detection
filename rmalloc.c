/* =====================================================================
   File:	rmalloc.c
   Author:	Rammi <rammi@quincunx.escape.de>
   Date:	11/16/1995 (started)

   Disclaimer:	Use at your own risk.
		Released to the public domain.

   Content:	Debug wrapper functions for the malloc library.
   		For more information see rmalloc.h

   Changes:
   		04/11/1996 (Rammi)
		Changed to hashed table for faster access
   		04/15/1996 (Rammi)
		Included statistics
		08/16/1997 (Rammi)
		Automatic output of used memory on exit
		08/25/1997 (Rammi)
		Catching signals in critical situations (partly)
		02/18/1998 (Rammi)
		Testing memory before outputting statistic
		Overworked signal handling in IsPossibleFilePos()
		Made it unnecessary of changing all mallocs etc
		to upper case
		02/19/1998 (Rammi)
		Little changes to compile on Alpha (64bit) without
		warnings
		03/10/1998 (Rammi)
		Added comments.
		03/24/1998 (Rammi)
		Allowed compilation without WITH_FLAGS
		04/07/1998 (Rammi)
		All output goes to stderr.
		1.11beta is released as 1.11!
		05/28/1998 (Rammi)
		Changed comments to english for public release
		Added some more signal handling.
		06/01/1998 (Rammi)
		Added output of (flagged) strings in ELOQUENT mode.
		Changed all names from my_... to R...
		This is version 1.12!
		11/13/1998 (Rammi)
		Multiple defined label when using ELOQUENT mode now fixed.
		Added getcwd wrapper as a prototype how to handle library
		functions returning heap memory.
		This is version 1.13!
	        06/10/99 (Rammi)
		The getcwd wrapper had a bug as Greg Silverman pointed 
		out (I should better have tested it!). Also fixed a
		missing prototype and improved the signal handling in
		rtest to allow receiving more signals while handling 
		one. Version is now 1.14.

		07/06/99 (Rammi)
		Put everything in CVS. Started new efforts to expand the
		debug library to also include heap allocated in (standard)
		library functions.
   =====================================================================
   CVS STUFF
   ---------------------------------------------------------------------
   Name:	$Source: /research/red/dbader/MPI/cycle/RCS/rmalloc.c,v $

   Author:	$Author: dbader $

   Date:	$Date: 1999/07/15 01:30:19 $

   Revision:	$Revision: 1.1 $

   History:	$Log: rmalloc.c,v $
   History:	Revision 1.1  1999/07/15 01:30:19  dbader
   History:	Initial revision
   History:
   History:	Revision 1.8  1999/07/14 20:49:43  rammi
   History:	*** empty log message ***
   History:	
   History:	Revision 1.7  1999/07/14 20:45:26  rammi
   History:	*** empty log message ***
   History:	
   History:	Revision 1.6  1999/07/14 20:42:18  rammi
   History:	Changed tracked in traced in Statistics
   History:	
   History:	Revision 1.5  1999/07/12 18:35:16  rammi
   History:	Added public domain stuff to header.
   History:	Moved include of malloc.c and made it depending on TRACK_EXTERNAL_CALLS flag
   History:	Added TRACK_EXTERNAL_CALLS to allow easy switch-back to version 1.14
   History:	Renamed NO_LIB_ABORT to ABORT_ON_EXTERNAL_ERRORS with inversed meaning (to allow setting from cc command line).
   History:	Added Wrapper for Doug´s routines vs standard lib routines.
   History:	Added one constant tag for all external calls to allow easy compare of pointers (no strcmp necessary).
   History:	Added shortcut for external calls in IsPossibleFilePos()
   History:	Added special handling of errors in external calls depending on ABORT_ON_EXTERNAL_ERRORS
   History:	Moved Global.isInitialized setting in Initialize() because some implementations use malloc in fprintf/atexit
   History:	Added some lines to header shown on initialisation.
   History:	TestAll() now depends on previous initialisation.
   History:	Changed calls to standard malloc lib to use wrapper macros.
   History:	Added external call logic to Rmalloc_stat.
   History:	Changed Rmalloc_retag and Rmalloc_set to return the handled buffer (allowing nested calls).
   History:	
   History:	Revision 1.4  1999/07/06 15:30:21  rammi
   History:	Extended IsPossibleFilePos to accept the LIB_TAG.
   History:	
   History:	Revision 1.3  1999/07/06 15:10:09  rammi
   History:	Added support for NO_LIB_ABORT precompiler flag.
   History:	
   History:	Revision 1.2  1999/07/06 15:01:43  rammi
   History:	Added support for malloc in libraries.
   History:	Standard calls to malloc/realloc/free/strdup are now also
   History:	passed thru the malloc debug library (when rmalloc.o is
   History:	linked as last program before libc).
   History:	This is a beta version for testing.
   History:	
   ===================================================================== */


/* =========
   INCLUDEs:
   ========= */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <strings.h>
#include <assert.h>
#include <setjmp.h>
#include <signal.h>

#undef  MALLOC_DEBUG		/* here we need the correct malloc functions */
#define RM_NEED_PROTOTYPES	/* but we want to compare prototypes */
#include "rmalloc.h"

/* ========
   DEFINEs:
   ======== */

/* Actual version */
#define VERSION         "2.00beta0"

/* ================================================================== */
/* ============ Switch settings for different behaviours ============ */
/* ============        Please set as needed              ============ */
/* ================================================================== */

/* This switch sets, how and when the allocated blocks are tested 
 * on correctness. Each block is tested at least when 
 * reallocating/freeing it.
 * Possible values:
 * 	0:		Minimum testing. Uses less memory, but 
 *	                does not allow statistics.
 * 	1:		Extra testing possible by using RM_TEST
 *	                macro. Statistics possible.
 *   	2:		Testing ALL blocks on every malloc/free.
 *	                Statistics possible.
 */
#ifndef RM_TEST_DEPTH
#define RM_TEST_DEPTH 	1
#endif

/* Set this to track calls from non-prepared files (i.e. files
 * which did not include rmalloc.h = especially libraries).
 * This file has to be linked after all libraries but before
 * libc!
 * Setting thishas the advantage that every malloc call is passed
 * through the debug library but has the disadvantage that you must 
 * not link this file to release versions (besides not setting
 * MALLOC_DEBUG when compiling your sources)!  
 */
#define TRACK_EXTERNAL_CALLS

/* Switch on EXTENTED alloc information. (Makes sense if an actual 
 * error is observed such as a multiple free)
 */
/*#define ELOQUENT*/

/* Allows setting of special flags with the RM_SET_FLAGS nacro.
 * Needs more memory.
 */
#define WITH_FLAGS 

/* Switch on if you want to abort on errors in libraries.
 */
#if defined(WITH_FLAGS) && defined(TRACK_EXTERNAL_CALLS)
/*#define ABORT_ON_EXTERNAL_ERRORS*/
#endif /* WITH_FLAGS  &&  TRACK_EXTERNAL_CALLS */

/* Allows realloc(NULL, ...) 
 * Posix allows this but there are some old malloc libraries 
 * which crash on this. Switch on if you want to be compatible.
 */
/* #define ALLOW_REALLOC_NULL */

/* Allows free(NULL)
 * I still consider this an error in my progs because I use 
 * NULL always as a very special value.
 */
/* #define ALLOW_FREE_NULL */


/* ================================================================== */
/* ========================  Other Defines  ========================= */
/* ================================================================== */

/* Take Doug's or system's routines? */
#ifdef TRACK_EXTERNAL_CALLS
   /* Doug Lea's routines */
#  define MALLOC(x)	mALLOc(x)
#  define REALLOC(x,s)	rEALLOc(x,s)
#  define FREE(x)	fREe(x)
#else /* !TRACK_EXTERNAL_CALLS */
   /* System's routines */
#  define MALLOC(x)	malloc(x)
#  define REALLOC(x,s)	realloc(x,s)
#  define FREE(x)	free(x)
#endif /* TRACK_EXTERNAL_CALLS */

/* Alignment (8 bytes are ok for virtually all machines) */
#define ALIGNMENT	8

/* Output message header: */
#define HEAD		"<MALLOC_DEBUG>\t"

/* Statistics message header: */ 
#define STAT_HEAD	"<MALLOC_STATS>\t"

/* Alignment padding: */
#define ALIGN(s)	(((s+ALIGNMENT-1)/ALIGNMENT)*ALIGNMENT)

/* Magic marker for block start: */
#define PREV_STOP 	0x55555555

/* Additional space for block begin to keep alignment: */   
#define START_SPACE     ALIGN(sizeof(begin))

/* Additional space at end of block */
#define END_SPACE       (sizeof(End))

/* Overall additional space per block: */
#define EXTRA_SPACE     (START_SPACE+END_SPACE)

/* Hashtable size: */
#define HASHSIZE	257

/* Hash function: */
#define HASH(p)		((((unsigned long)(p))/ALIGNMENT)%HASHSIZE)


/* ==========================
   STRUCTs, TYPEDEFs & ENUMs:
   ========================== */

/* This Information is added to the beginning of every
 * allocated block.
 */
typedef struct _begin {
  unsigned       StpA;		/* Magic bytes */
#if RM_TEST_DEPTH > 0
  struct _begin	*Next,		/* for linking in forward */
	        *Prev;		/* and backward direction */
#endif
  const char	*File;		/* Filepos of allocation command */
  size_t	 Size;		/* Size demanded */
#ifdef WITH_FLAGS
  unsigned       Flags;		/* Special flags */
#endif
  unsigned 	 StpB;		/* Magic bytes */
} begin;


/*
 * Global data.
 */
typedef struct _global {
  unsigned	 isInitialized;	/* Flag: already initialized? */
  unsigned	 BlockCount;    /* Number of allocated blocks */
} global;



/* =======
   CONSTs: 
   ======= */

/* Magic block end: */
static unsigned char End[] = {
  0xA5, 0xA5, 0xA5, 0xA5,	/* odd */
  0x5B, 0x5B, 0x5B, 0x5B,	/* odd */
  0xAB, 0xAB, 0xAB, 0xAB,	/* odd */
  0xAA, 0x55, 0xAA, 0x55	/* odd */
};


/* Position tag for direct calls of malloc etc.: */
static char *EXTERNAL_CALL_TAG = "[external call]";


/* ========
   GLOBALs:
   ======== */



/* =======
   LOCALs:
   ======= */

#if RM_TEST_DEPTH > 0
/* Stop marker for linked list of allocated blocks: */
static begin ChainTempl = {
  PREV_STOP,
  &ChainTempl,
  &ChainTempl,
  "Special",
  0,
#ifdef WITH_FLAGS
  0,
#endif
  PREV_STOP
};

/* Hashed linked lists: */
static begin Chain[HASHSIZE];

/* Global data used: */
static global Global = {
  0,				/* is initialized?  */
  0				/* number of blocks */
};
#endif


/* ========
   FORWARD:
   ======== */

static int       FindBlk(const unsigned char *P);
static jmp_buf   errorbuf;


/* ===============================================================
   			IMPLEMENTATION
   =============================================================== */

#ifdef TRACK_EXTERNAL_CALLS
/* We need a special malloc implementation with changed names */
#include "malloc.c"
#endif


/* =============================================================================
   Function:		FatalSignal	// local //
   Author:		Rammi
   Date:		08/25/1997

   Return:		---

   Parameter:		signum    signal number

   Purpose:		Signal handler für fatal signals (SIGBUS, SIGSEGV)
   ============================================================================= */
static void FatalSignal(int signum) 
{
  /* --- jump to a save place --- */
  longjmp(errorbuf, signum);
}


/* =============================================================================
   Function:		IsPossibleFilePos	// local //
   Author:		Rammi
   Date:		11/30/1996

   Return:		!= 0    possibly ok
                        0       seems not so

   Parameter:		file	possible filename
			size    possible size

   Purpose:		Decide whether file could be a filename and 
   			size a block size.
   ============================================================================= */
static int IsPossibleFilePos(const char *file, int size)
{
  void   (*old_sigsegv_handler)(int) = SIG_DFL;
  void   (*old_sigbus_handler)(int)  = SIG_DFL;
  char    *dp;
  int      ret;

  if (file == EXTERNAL_CALL_TAG) {
    /* called by an unwrapped routine */
    return size >= 0;
  }
  
  if (setjmp(errorbuf)) {
    /* uh oh, we got a kick in the ass */
    signal(SIGSEGV, old_sigsegv_handler);
    signal(SIGBUS,  old_sigbus_handler);
    return 0;
  }
  
  /* --- the following is dangerous! So catch signals --- */
  old_sigsegv_handler = signal(SIGSEGV, FatalSignal);
  old_sigbus_handler  = signal(SIGBUS,  FatalSignal);
  
  dp  = strchr(file, ':');	/* file pos needs : */
  
  ret =  (dp   &&   dp-file > 3   &&   !strncmp(dp-2, ".c", 2)   &&   
	  atoi(dp+1) > 0   &&   size >= 0);
  
  /* --- danger is over! --- */
  signal(SIGSEGV, old_sigsegv_handler);
  signal(SIGBUS,  old_sigbus_handler);
  
  return ret;
}

/* =============================================================================
   Function:		ControlBlock	// lokal //
   Author:		Rammi
   Date:		11/16/1995

   Return:		---
   Parameter:		Bkl	Pos of allocated block (original)
   			file	file pos from where initial lib function
				was called

   Purpose:		Control integrity of block

   ============================================================================= */
static void ControlBlock(begin *B, const char *file)
{
  unsigned char *p = (((unsigned char *)B)+START_SPACE);
#if RM_TEST_DEPTH > 0
  int DoAbort = 0;
#endif
  /* === the very beginning === */
  if (B->StpA != PREV_STOP) {
#if RM_TEST_DEPTH > 0
    DoAbort = 1;
#endif
    fprintf(stderr, HEAD 
	    "Corrupted block begin (overwritten from elsewhere)\n"
	    "\tshould be: %08x\n"
	    "\tis:        %08x\n"
	    "\tblock was allocated in %s [%u Bytes]\n"
	    "\terror was detected in  %s\n",
	    PREV_STOP,
	    B->StpA,
	    B->File,
	    (unsigned) B->Size,
	    file);
  }
  
  /* === begin of user data === */
  if (B->StpB != PREV_STOP) {
#if RM_TEST_DEPTH > 0
    DoAbort = 1;
#endif
    fprintf(stderr, HEAD 
	    "Corrupted block begin (possibly written back)\n"
	    "\tshould be: %08x\n"
	    "\tis:        %08x\n"
	    "\tblock was allocated in %s [%u Bytes]\n"
	    "\terror was detected in  %s\n",
	    PREV_STOP,
	    B->StpB,
	    B->File,
	    (unsigned) B->Size,
	    file);
  }
  
  /* === end of user data === */
  if (memcmp(p+B->Size, End, END_SPACE) != 0) {
    unsigned char *E = (unsigned char *)(p+B->Size);
    int i;
    int found = 0;
#if RM_TEST_DEPTH > 0
    DoAbort = 1;
#endif
    fprintf(stderr, HEAD 
	    "Corrupted block end (possibly written past the end)\n"
	    "\tshould be:");
    for (i = 0;   i < END_SPACE;   i++) {
      fprintf(stderr, i%4 ? "%02x" : " %02x", End[i]);
    }
    
    fprintf(stderr,
	    "\n\tis:       ");
    for (i = 0;   i < END_SPACE;   i++) {
      fprintf(stderr, i%sizeof(int) ? "%02x" : " %02x", E[i]);
    }
    fprintf(stderr, "\n\tblock was allocated in %s [%u Bytes]\n"
	    "\terror was detected in  %s\n",
	    B->File,
	    (unsigned) B->Size,
	    file);
    
#if RM_TEST_DEPTH > 0
    if (!((unsigned long)E % sizeof(void *))   &&   
	!(*(unsigned long *)E % sizeof(void *))) {  /* because of alignment */
      /* Special service: look if memory was overwritten with pointer */
      if (FindBlk(*(unsigned char **)E)) {
	begin *b = (begin *)((*(unsigned char **)E)-START_SPACE);
	if (IsPossibleFilePos(b->File, b->Size)) {
	  fprintf(stderr, 
		  "\tFirst %d bytes of overwritten memory can be interpreted\n"
		  "\t\tas a pointer to a block allocated in:\n"
		  "\t\t%s [%u Bytes]\n",
		  sizeof(void *),
		  b->File,
		  (unsigned) b->Size);
	  found = 1;
	}
      }
    }
    if (!found)
#endif	
    {
      /* Look, what we can find... */
      int  j;
      
      for (j = END_SPACE-1;   j >= 0;   j--) {
	if (E[j] != End[j]) {
	  break;
	}
      }
      if (j >= 0   &&   !E[j]) {
	/* Ends with '\0', so it's possibly a string */
	if (j > 0) {
	  while (--j >= 0) {
	    if (!E[j]) {
	      break;
	    }
	  }
	  if (j < 0) {
	    fprintf(stderr, 
		    "\tLooks somewhat like a too long string,\n"
		    "\t\tending with \"%s\"\n",  
		    E);
	  }
	}
	else {
	  /* Off by one? */
	  fprintf(stderr,
		  "\tLooks like string allocated one byte too short\n"
		  "\t\t(missing the null byte)\n");
	}
      }
    }  
  }
    
#if RM_TEST_DEPTH > 0
  /* Die LOUD */
  if (DoAbort) {
#ifndef ABORT_ON_EXTERNAL_ERRORS
    if (!(B->Flags & RM_EXTERNAL)) {
      abort();
    }
#else
    abort();
#endif
  }
#endif
}



#if RM_TEST_DEPTH > 0
void Rmalloc_stat(const char *file);

/* =============================================================================
   Function:		Exit // local //
   Author:		Rammi
   Date:		08/19/1997

   Return:		---
   Parameter:		---

   Purpose:		Function called on exit
   ============================================================================= */
static void Exit(void) 
{
  Rmalloc_stat("[atexit]");	/* show statistics */
}
   

/* =============================================================================
   Function:		Initialize	// local //
   Author:		Rammi
   Date:		11.04.1996

   Return:		---
   Parameter:		---

   Purpose:		Necessary initializations

   ============================================================================= */
static void Initialize(void) 
{
  int i;
  
  Global.isInitialized = 1;	/* some implementations use malloc in */
				/* fprintf() or atexit() so this has  */
				/* to be set early! */

  fprintf(stderr,
	  HEAD "rmalloc -- malloc wrapper V " VERSION "\n"
	  "\tby Rammi <mailto:rammi@quincunx.escape.de>\n"
	  "\tInfo on http://www.escape.de/user/quincunx/rmdebug.html\n"
	  "\tCompiled with following options:\n"
#if RM_TEST_DEPTH==1
	  "\t\ttesting:\tonly actual block\n"
#elif RM_TEST_DEPTH==2
	  "\t\ttesting:\tall allocated blocks\n"
#else
	  "\t\ttesting:\tcomment missing in " __FILE__ ":" INT2STRING(__LINE__) "\n"
#endif
#ifdef TRACK_EXTERNAL_CALLS
	  "\t\texternal calls:\ttracked\n"
# if ABORT_ON_EXTERNAL_ERRORS
	  "\t\t   ext. errors:\tabort\n"
# else
	  "\t\t   ext. errors:\tdon't abort\n"
# endif
#else
	  "\t\texternal:\tnot tracked\n"
#endif
#ifdef ELOQUENT
	  "\t\teloquence:\tON\n"
#else
	  "\t\teloquence:\tOFF\n"
#endif
#ifdef ALLOW_REALLOC_NULL
	  "\t\trealloc(0):\tALLOWED\n"
#else
	  "\t\trealloc(0):\tNOT ALLOWED\n"
#endif
#ifdef ALLOW_FREE_NULL
	  "\t\tfree(0):\tALLOWED\n"
#else
	  "\t\tfree(0):\tNOT ALLOWED\n"
#endif
#ifdef WITH_FLAGS
	  "\t\tflags:  \tUSED\n"
#else
	  "\t\tflags:  \tUNUSED\n"
#endif
	  "\t\talignment:\t" INT2STRING(ALIGNMENT) "\n"
	  "\t\tpre space:\t%d\n"
	   "\t\tpost space:\t%d\n"
	  "\t\thash tab size:\t" INT2STRING(HASHSIZE) "\n\n",
	  START_SPACE, END_SPACE);

  /* --- init list heads --- */  
  for (i = 0;   i < HASHSIZE;   i++) {
    memcpy(Chain+i, &ChainTempl, sizeof(begin));
    Chain[i].Next = Chain[i].Prev = Chain+i;
  }
  
  /* --- show statistics at exit --- */
  (void)atexit(Exit);
  
}


/* =============================================================================
   Function:		TestAll		// local //
   Author:		Rammi
   Date:		16.11.1995

   Return:		---
   Parameter:		file		file pos where lib function was 
   					called

   Purpose:	        Test all allocated blocks for inconsistencies
   ============================================================================= */
static void TestAll(const char *file)
{

  if (Global.isInitialized) {
    begin *B;
    int    i;

    for (i = 0;   i < HASHSIZE;   i++) {
      B = Chain[i].Next;
      
      /* === Once around the circle === */
      while (B != &Chain[i]) {
	ControlBlock(B, file);
	B = B->Next;
      }
    }
  }

}


/* =============================================================================
   Function:		AddBlk		// local //
   Author:		Rammi
   Date:		16.11.1995

   Return:		---
   Parameter:		Blk		New block (original pos.)
   			file		called from 

   Purpose:		Add new block to the list
   ============================================================================= */
static void AddBlk(begin *Blk, const char *file)
{
  int hash = HASH(Blk);		/* hash val */
  
  if (!Global.isInitialized) {
    Initialize();
  }

#if RM_TEST_DEPTH > 1
  TestAll(file);
#endif
  /* --- insert it --- */
  Blk->Next = Chain[hash].Next;
  Blk->Prev = &Chain[hash];
  Chain[hash].Next->Prev = Blk;
  Chain[hash].Next = Blk;
  
  Global.BlockCount++;
  
}


/* =============================================================================
   Function:		DelBlk		// local //
   Author:		Rammi
   Date:		16.11.1995

   Return:		---

   Parameter:		Blk		block to remove
   			file		called from

   Purpose:		Remove block from list.
   			React angry if block is unknown
   ============================================================================= */
static void DelBlk(begin *Blk, const char *file)
{
  begin *B;			/* run var  */
  int    hash = HASH(Blk);	/* hash val */

  /* look if block is known */
  for (B = Chain[hash].Next;   B != &Chain[hash];   B = B->Next) {
    if (B == Blk) {
      goto found_actual_block;	/* friendly goto */
    }
  }
  /* not found */
  fprintf(stderr, HEAD
	  "Double or false delete\n"
	  "\tHeap adress of block: %p\n"
	  "\tDetected in %s\n",
	  ((char *)Blk)+START_SPACE, file);
  {
    void   (*old_sigsegv_handler)(int) = SIG_DFL;
    void   (*old_sigbus_handler)(int)  = SIG_DFL;
    
    if (setjmp(errorbuf)) {
      /* uh oh, we got a kick in the ass */
      signal(SIGSEGV, old_sigsegv_handler);
      signal(SIGBUS,  old_sigbus_handler);
    }
    else {
      /* --- the following is dangerous! So catch signals --- */
      old_sigsegv_handler = signal(SIGSEGV, FatalSignal);
      old_sigbus_handler  = signal(SIGBUS,  FatalSignal);

      if (IsPossibleFilePos(Blk->File, Blk->Size)) {
	fprintf(stderr, 
		"\tTrying identification (may be incorrect!):\n"
		"\t\tAllocated in %s [%u Bytes]\n",
		Blk->File, (unsigned) Blk->Size);
      }
      signal(SIGSEGV, old_sigsegv_handler);
      signal(SIGBUS,  old_sigbus_handler);
    }
  }
#ifndef ABORT_ON_EXTERNAL_ERRORS
  if (file != EXTERNAL_CALL_TAG) {
    abort();			/* die loud */
  }
#else
  abort();			/* die loud */
#endif
  
found_actual_block:
#if RM_TEST_DEPTH > 1
  /* check everything */
  TestAll(file);
#else 
  /* test integrity of actual block */
  ControlBlock(Blk, file);
#endif

  /* remove: */
  Blk->Next->Prev = Blk->Prev;
  Blk->Prev->Next = Blk->Next;
  
  Global.BlockCount--;
  
#ifdef ELOQUENT    
  fprintf(stderr,
	  HEAD "Delete: %d Bytes allocated in %s (from %s)\n", 
	  Blk->Size, Blk->File, file);
#ifdef WITH_FLAGS
  if (Blk->Flags & RM_STRING) {
    char *c;
    /* look for eos */
    for (c = (char *)Blk + START_SPACE;  
	 c - (char *)Blk + START_SPACE < Blk->Size;   
	 c++) {
      if (!*c) {
	fprintf(stderr,
		HEAD "\tContains string: \"%s\"\n",
		(char *)Blk + START_SPACE);
	goto found_old_block;
      }
    }
    /* not found */
    fprintf(stderr,
	    HEAD "\tContains string without null byte\n");
found_old_block:
    ;
  }
#endif /* WITH_FLAGS */
#endif /* ELOQUENT */

#ifdef WITH_FLAGS
  if (Blk->Flags & RM_STATIC) {
    fprintf(stderr,
	    HEAD "WARNING: freeing block marked as STATIC (in %s)\n", file);
  }
#endif /* WITH_FLAGS */
}



/* =============================================================================
   Function:		FindBlk		// local //
   Author:		Rammi
   Date:		11/30/1996

   Return:		0               not found
                        1               found

   Parameter:		P		block (user pos)

   Purpose:		look if block is known
   ============================================================================= */
static int FindBlk(const unsigned char *P)
{
  begin *B;
  const begin *Blk = (const begin *)(P - START_SPACE);
  int hash = HASH(Blk);
  
  /* look if block is known */
  for (B = Chain[hash].Next;   B != &Chain[hash];   B = B->Next) {
    if (B == Blk) {
      
      return 1;
    }
  }
  return 0;
}
#endif /* RM_TEST_DEPTH > 0 */



/* =============================================================================
   Function:		SetBlk		// local //
   Author:		Rammi
   Date:		11/16/1995

   Return:		pointer to block (user pos.)
   
   Parameter:		Blk		pointer to block (original pos.)
   			size		size (user)
   			file		called from
			flags           flags (when compiled WITH_FLAGS)

   Purpose:		Setzt unsere internen Informationen

   ============================================================================= */
#ifdef WITH_FLAGS
static void *SetBlk(void *Blk, size_t size, const char *file, unsigned flags)
#else
static void *SetBlk(void *Blk, size_t size, const char *file)
#endif
{

  ((begin *)Blk)->StpA  = PREV_STOP;
  ((begin *)Blk)->File  = file;
  ((begin *)Blk)->Size  = size;
#ifdef WITH_FLAGS
  ((begin *)Blk)->Flags = flags;
#endif
  ((begin *)Blk)->StpB  = PREV_STOP;
  memcpy(((char *)Blk)+START_SPACE+size, End, END_SPACE);
  
#if RM_TEST_DEPTH > 0
  AddBlk((begin *)Blk, file);
#endif
#ifdef ELOQUENT
  fprintf(stderr,
	  HEAD "Adding: %p, %d Bytes (from %s)\n", 
	  ((char *)Blk)+START_SPACE, size, file);
#endif /* ELOQUENT */
  return ((char *)Blk)+START_SPACE;
}



/* =============================================================================
   Function:		Rmalloc	// external //
   Author:		Rammi
   Date:		11/16/1995

   Return:		New prepared memory block with size size (user)
   
   Parameter:		size		demanded size
   			file		called from where?

   Purpose:		wrapper for malloc
   ============================================================================= */
void *Rmalloc(size_t size, const char *file)
{
  void *ret;			/* ret val */

  if (size == 0) {
    fprintf(stderr, HEAD "WARNING: malloc() demands 0 Bytes (in %s)\n", file);
  }


  ret = MALLOC(size+EXTRA_SPACE); /* get extended block  */

  if (ret) {
    /* initialize */
#ifdef WITH_FLAGS
    return SetBlk(ret, size, file, 0);
#else
    return SetBlk(ret, size, file);
#endif
  }
  else {
    fprintf(stderr,
	    HEAD "WARNING: Out of memory! Returning NULL (in %s)\n", file);
    return NULL;
  }
}



/* =============================================================================
   Function:		Rcalloc		// external //
   Author:		Rammi
   Date:		11/16/1995

   Return:		New (cleared) memory block of size nelem*size
   
   Parameter:		nelem		nr blocks (as stupid as calloc)
   			size		size of one block
   			file		called from

   Purpose:		Wrapper function for calloc
   ============================================================================= */
void *Rcalloc(size_t nelem, size_t size, const char *file) 
{
  void *ret;
  
  /* calculate correct size now */
  size *= nelem;
  
  if (size == 0) {
    fprintf(stderr,
	    HEAD "WARNING: calloc() demands 0 Bytes (in %s)\n", file);
  }
    
  /* Rmalloc makes nearly all the work */
  ret = Rmalloc(size, file);
    
  if (ret) {
    /* clear */
    memset(ret, 0, size);
    return ret;
  }
  else {
    fprintf(stderr,
	    HEAD "WARNING: Out of memory! Returning NULL (in %s)\n", file);
    return NULL;
  }
}



/* =============================================================================
   Function:		Rrealloc	// external //
   Author:		Rammi
   Date:		11/16/1995

   Return:		New block of size size (user pos.)
   
   Parameter:		p		previous pointer
   			size		new size
   			file		called from

   Purpose:		Wrapper function for realloc
   ============================================================================= */
void *Rrealloc(void *p, size_t size, const char *file)
{
  void     *ret;
#ifdef WITH_FLAGS
  unsigned  flags = 0;
#endif
  
  if (p == NULL) {
#ifndef ALLOW_REALLOC_NULL
    fprintf(stderr, HEAD "Realloc of NULL pointer (in %s)\n", file);
    abort();
#else /* ALLOW_REALLOC_NULL */
#ifdef ELOQUENT
    fprintf(stderr,
	    HEAD "WARNING: realloc of NULL pointer (in %s)\n", file);
#endif
    return Rmalloc(size, file);
#endif /* ALLOW_REALLOC_NULL */
  }
#ifdef WITH_FLAGS
  else {
    /* keep flags */
    flags = ((begin *)(((char *)p)-START_SPACE))->Flags;
  }
#endif /* WITH_FLAGS */
  
  if (size == 0) {
    fprintf(stderr,
	    HEAD "WARNING: realloc() demands 0 Bytes (in %s)\n", file);
  }
  
#if RM_TEST_DEPTH > 0
  /* remove old block from list */
  DelBlk((begin *)(((char *)p)-START_SPACE), file);
#endif
  /* get everything new */
  ret = REALLOC(((char *)p)-START_SPACE, size+EXTRA_SPACE);

  if (ret) {
    /* Initialize new block */
#ifdef WITH_FLAGS
    return SetBlk(ret, size, file, flags);
#else
    return SetBlk(ret, size, file);
#endif
  }
  else {
    fprintf(stderr,
	    HEAD "WARNING: Out of memory! Returning NULL (in %s)\n", file);
    return NULL;
  }
}



/* =============================================================================
   Function:		Rfree		// external //
   Author:		Rammi
   Date:		11/16/1995

   Return:		---
   
   Parameter:		p		block to free (user pos.)
   			file		called from

   Purpose:		Wrapper function for free()

   ============================================================================= */
void Rfree(void *p, const char *file)
{
#ifdef ELOQUENT
  fprintf(stderr, HEAD "Free: %p (called from: %s)\n", p, file);
#endif /* ELOQUENT */
  if (p == NULL) {
#ifdef ALLOW_FREE_NULL
#ifdef ELOQUENT
    fprintf(stderr, HEAD "WARNING: Freeing NULL pointer (in %s)\n", file);
#endif
    return;
#else /* !ALLOW_FREE_NULL */
    fprintf(stderr, HEAD "Trying to free NULL pointer (in %s)\n", file);
    abort();
#endif /* !ALLOW_FREE_NULL */
  }
#if RM_TEST_DEPTH > 0
  /* Remove block from list */
  DelBlk((begin *)(((char *)p)-START_SPACE), file);
#endif
  /* free block */
  FREE(((char *)p)-START_SPACE);
  
}    	



/* =============================================================================
   Function:		Rstrdup		// external //
   Author:		Rammi
   Date:		11/16/1995

   Return:		New memory with copied string.
   
   Parameter:		s		string to copy
   			file		called from

   Purpose:		Wrapper function for strdup()
   ============================================================================= */
char *Rstrdup(const char *s, const char *file) 
{
  size_t size;	/* needed memory */
  char *ret;
  
  if (s == NULL) {
    fprintf(stderr, HEAD "Calling strdup(NULL) (in %s)\n", file);
    abort();
  }
  size = strlen(s)+1;
    
  /* Rmalloc() does nearly all the work */
  ret = Rmalloc(size, file);
  
  if (ret) {
    /* copy string */
    strcpy(ret, s);
#ifdef WITH_FLAGS
    Rmalloc_set_flags(ret, RM_STRING, "<by strdup>");
#endif
    return ret;
  }
  else {
    fprintf(stderr,
	    HEAD "WARNING: Out of memory! Returning NULL (in %s)\n", file);
    return NULL;
  }    
    
}


/* =============================================================================
   Function:		Rgetcwd		// external //
   Author:		Rammi
   Date:		11/13/1998

   Return:		New memory with copied string depending on input (if
                        buffer == NULL)
   
   Parameter:	        buffer		buffer for write (or NULL)
                        size		buffer size
   			file		called from

   Purpose:		Wrapper function for getcwd() which sometimes returns
			memory from heap.
   ============================================================================= */
char *Rgetcwd(char *buffer, size_t size, const char *file) 
{
  char *ret = getcwd(buffer, size);

  if (ret && !buffer) {
    /* create new memory to get internals fixed */
    char *newret = Rstrdup(ret, file);
    FREE(ret);			/* free old stuff */
    ret = newret;		/* this was missing before 1.14 */
				/* thanks to Greg Silverman who discovered it! */
  }
  
  return ret;
}


/* =============================================================================
   Function:		Rmalloc_test		// external //
   Author:		Rammi
   Date:		04/11/1995

   Return:		---
   
   Parameter:		file		called from

   Purpose:		Explicitely test all blocks for integrity 

   ============================================================================= */
void Rmalloc_test(const char *file)
{
#if RM_TEST_DEPTH > 0
  TestAll(file);
#else
  fprintf(stderr, HEAD __FILE__ 
	  " not compiled with RM_TEST_DEPTH > 0, call in %s senseless.\n", file,);
#endif
}



/* =============================================================================
   Function:		BlockSort		// local //
   Author:		Rammi
   Date:		04/15/1995

   Return:		< 0		A < B
   			0		A == B
   			> 0		A > B
   
   Parameter:		A, B
   
   Purpose:		sort function for qsort

   ============================================================================= */
static int BlockSort(const begin **A, const begin **B)
{
  int   ret;
  
  /* Sort for adress of string (tricky!) */
  if ((ret = (*A)->File - (*B)->File)) {
    return ret;
  }
  
  /* sort for size */
  return ((int)(*A)->Size - (int)(*B)->Size);
}



/* =============================================================================
   Function:		Rmalloc_stat		// extern //
   Author:		Rammi
   Date:		04/15/1995

   Return:		---
   
   Parameter:		file		caled from

   Purpose:		Show statistic
   ============================================================================= */
void Rmalloc_stat(const char *file)
{
#if RM_TEST_DEPTH > 0
  TestAll(file);

  fprintf(stderr,
	  STAT_HEAD "============ STATISTICS (%s) =============\n", file);
  if (!Global.BlockCount) {
    fprintf(stderr, STAT_HEAD "Nothing allocated.\n");
  }
  else {
    const begin	**BlockVec;
    
    if ((BlockVec = (const begin **)MALLOC(Global.BlockCount*sizeof(begin *))) == NULL) {
      fprintf(stderr, STAT_HEAD "Couldn't allocate enough memory for statistics. Going on...\n");
    }
    else {
      int	 i = 0;
      int	 j;
      begin	*B;
      int	 count;
      size_t	 Mem = 0;
      int        nrBlocks;
#ifdef WITH_FLAGS
      size_t     StaticMem   = 0;
      size_t	 ExternalMem = 0;
#endif
      
      /* add all blocks to vector */
      for (j = 0;   j < HASHSIZE;   j++) {
	for (B = Chain[j].Next;   B != &Chain[j];   B = B->Next) {
#ifdef WITH_FLAGS
	  if (B->Flags & RM_STATIC) {
	    StaticMem += B->Size;
	  }
#ifdef TRACK_EXTERNAL_CALLS
	  else if (B->Flags & RM_EXTERNAL) {
	    ExternalMem += B->Size;
	  }
#endif /* TRACK_EXTERNAL_MEM */
	  else {
	    BlockVec[i++] = B;
	  }
#else  /* !WITH_FLAGS */
	  BlockVec[i++] = B;
#endif /* !WITH_FLAGS */
	}
      }
#ifdef WITH_FLAGS
      assert(i <= Global.BlockCount);
#else
      assert(i == Global.BlockCount);
#endif
      nrBlocks = i;
      
      qsort(BlockVec, nrBlocks, /* sort */
	    sizeof(begin *), 
	    (int (*)(const void *, const void *))BlockSort);
      
      for (i = 0;   i < nrBlocks;   i = j) {
	count = 1;
	for (j = i+1;   j < nrBlocks;   j++) {
	  if (BlockSort(BlockVec+i, BlockVec+j) != 0) {
	    break;
	  }
	  /* are equal */
	  count++;
	}
	fprintf(stderr,
		STAT_HEAD "%6d x %8u Bytes in %s\n",
		count,
		(unsigned) BlockVec[i]->Size,
		BlockVec[i]->File);
	Mem += count*BlockVec[i]->Size;
      }

      /* and give free */
      FREE(BlockVec);

      /* --- Summary --- */
#ifdef WITH_FLAGS
      fprintf(stderr, STAT_HEAD "*Variable*\t%12u Bytes\n",
	      (unsigned) Mem);
      fprintf(stderr, STAT_HEAD "*Static*  \t%12u Bytes\n",
	      (unsigned) StaticMem);
# ifdef TRACK_EXTERNAL_CALLS
      fprintf(stderr, STAT_HEAD "*External*\t%12u Bytes\n",
	      (unsigned) ExternalMem);
# endif
      fprintf(stderr, STAT_HEAD "---------------------------------------------\n");
      fprintf(stderr, STAT_HEAD "*TOTAL*   \t%12u Bytes\n",
	      (unsigned) (Mem+StaticMem+ExternalMem));
#else /* !WITH_FLAGS */
      fprintf(stderr, STAT_HEAD "*TOTAL*   \t%12u Bytes\n",
	      (unsigned) Mem);
#endif /* !WITH_FLAGS */
      fprintf(stderr, STAT_HEAD "=============================================\n");
      fprintf(stderr, STAT_HEAD "*Internal*\t%12u traced Blocks\n"
	              STAT_HEAD "         =\t%12u surplus Bytes\n",
	      (unsigned) (Global.BlockCount),
	      (unsigned) (Global.BlockCount*(START_SPACE+END_SPACE)));
    }
  }
  fprintf(stderr, STAT_HEAD "============= END OF STATISTICS =============\n");
#else
  fprintf(stderr, HEAD __FILE__ " not compiled with RM_TEST_DEPTH > 0, call in %s senseless.\n", file);
#endif
}



/* =============================================================================
   Function:		Rmalloc_retag		// external //
   Author:		Rammi
   Date:		12/12/1997

   Return:		The pointer parameter p.
   
   Parameter:		p		pointer to allocated block (user)
   			file		called from

   Purpose:		Change file position in header.
   ============================================================================= */
void *Rmalloc_retag(void *p, const char *file)
{
  if (p) {
    begin *info = (begin *)(((char *)p)-START_SPACE);
    
    /* --- test integrity --- */
    ControlBlock(info, file);
    
    /* --- change file pos --- */
    info->File = file;
  }
  return p;
}



/* =============================================================================
   Function:		Rmalloc_set_flags		// external //
   Author:		Rammi
   Date:		12/12/1997

   Return:		The pointer paramter p
   
   Parameter:		p               pointer to allocated block (user)
                        flags           Flags to set
                        file		called from

   Purpose:	        Set flags in header
   ============================================================================= */
void *Rmalloc_set_flags(void *p, unsigned flags, const char *file)
{
#ifdef WITH_FLAGS
  if (p) {
    begin *info = (begin *)(((char *)p)-START_SPACE);
    
    /* --- test integrity --- */
    ControlBlock(info, file);
    
    /* --- change flags --- */
    info->Flags |= flags;
  }
#endif
  return p;
}





/* =============================================================================
   Function:		Rmalloc_reinit		// external //
   Author:		Rammi
   Date:		05/28/1998

   Return:		---
   
   Parameter:		---

   Purpose:	        This reinits the lists. This is only for test purposes.
   			DON'T USE THIS FUNCTION!
   ============================================================================= */
void Rmalloc_reinit(void)
{
  int i;
  /* --- init list heads (discarding everything!) --- */  
  for (i = 0;   i < HASHSIZE;   i++) {
    memcpy(Chain+i, &ChainTempl, sizeof(begin));
    Chain[i].Next = Chain[i].Prev = Chain+i;
  }
}



#ifdef TRACK_EXTERNAL_CALLS
/* These are the routines called by external functions */


/* =============================================================================
   Function:		malloc			// external //
   Author:		Rammi
   Date:		07/06/1999

   Return:		newly allocated block
   
   Parameter:		size	block size

   Purpose:	        Overwriting standard malloc.
   ============================================================================= */
void *malloc(size_t size)
{
  void *block = Rmalloc(size, EXTERNAL_CALL_TAG);
  Rmalloc_set_flags(block, RM_EXTERNAL, EXTERNAL_CALL_TAG);
  return block;
}



/* =============================================================================
   Function:		realloc			// external //
   Author:		Rammi
   Date:		07/06/1999

   Return:		newly allocated block
   
   Parameter:		block	old block
			size	block size

   Purpose:	        Overwriting standard realloc.
   ============================================================================= */
void *realloc(void *block, size_t size)
{
  void *newblock = Rrealloc(block, size, EXTERNAL_CALL_TAG);
  Rmalloc_set_flags(newblock, RM_EXTERNAL, EXTERNAL_CALL_TAG);
  return newblock;
}



/* =============================================================================
   Function:		strdup			// external //
   Author:		Rammi
   Date:		07/06/1999

   Return:		newly allocated string
   
   Parameter:		str	string to copy

   Purpose:	        Overwriting standard strdup.
			(not really needed)
   ============================================================================= */
char *strdup(const char *str)
{
  char *ret = Rstrdup(str, EXTERNAL_CALL_TAG);
  Rmalloc_set_flags(ret, RM_EXTERNAL, EXTERNAL_CALL_TAG);
  return ret;
}



/* =============================================================================
   Function:		free			// external //
   Author:		Rammi
   Date:		07/06/1999

   Return:		---
   
   Parameter:		block	block to be freed

   Purpose:	        Overwriting standard free.
   ============================================================================= */
void free(void *block)
{
  Rfree(block, EXTERNAL_CALL_TAG);
}

#endif /* TRACK_EXTERNAL_CALLS */
