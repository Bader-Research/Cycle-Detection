/*

   $Log: mpi-printf.c,v $
   Revision 1.2  1999/05/27 19:15:54  dbader
   Daily Update 990523

   Revision 1.1  1999/05/27 19:09:08  dbader
   Initial revision


*/

static char rcsid[] = "$Id: mpi-printf.c,v 1.2 1999/05/27 19:15:54 dbader Exp $";

#include "mpi-printf.h"
#include "mpi.h"

void MPI_fprintf(FILE *fp, char *fmt, ...)
{
#define BUFSIZE 128
  int rank;
  int tag=10, p,np;
  va_list args;
  char str[BUFSIZE];
  MPI_Status status;
  va_start(args,fmt);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  if (rank==0) {
    fprintf(fp,"(PE%3d): ",rank);
    vfprintf(fp,fmt,args);
    for (p=1; p<np; p++) {
      MPI_Recv(str,BUFSIZE,MPI_CHAR,p,tag, MPI_COMM_WORLD, &status);
      fprintf(fp,"(PE%3d): ",p);
      fprintf(fp,"%s",str);
    }
    fflush(fp);
  } else {
    vsprintf(str,fmt,args);
    MPI_Send(str,BUFSIZE,MPI_CHAR,0,tag,MPI_COMM_WORLD);
  }
  va_end(args);
  return;
}

void MPI_printf(char *fmt, ...)
{
#if 0
  va_list args;
  va_start(args,fmt);
  MPI_fprintf(stdout,fmt,args);
  va_end(args);
  return;
#else
#define BUFSIZE 128
  int rank;
  int tag=10, p,np;
  va_list args;
  char str[BUFSIZE];
  MPI_Status status;
  va_start(args,fmt);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  if (rank==0) {
    printf("(PE%3d): ",rank);
    vprintf(fmt,args);
    for (p=1; p<np; p++) {
      MPI_Recv(str,BUFSIZE,MPI_CHAR,p,tag, MPI_COMM_WORLD, &status);
      printf("(PE%3d): ",p);
      printf("%s",str);
    }
    fflush(stdout);
  } else {
    vsprintf(str,fmt,args);
    MPI_Send(str,BUFSIZE,MPI_CHAR,0,tag,MPI_COMM_WORLD);
  }
  va_end(args);
  return;
#endif
}




