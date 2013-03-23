/***********************************************************
 ***********************************************************/
#ifndef MACROS_H
#define MACROS_H
/*
 * $Id: macros.h 1930 2010-02-16 14:46:27Z markus $
 *
 * defines a set of macros for time critical housekeeping tasks
 * Stolen from aln3nn
 *
 * Markus Riester, markus@bioinf.uni-leipzig.de */

/* allocate memory of given type and number of */
/* elements */
#define MALLOC(P,T,S) { if((P=(T*)malloc(sizeof(T)*(size_t)(S)))==NULL) { \
      FATAL("malloc failed"); \
    } memset(P,0,sizeof(T)*(size_t)(S)); \
  }

/* re-allocate memory of given type and number of */
/* elements */
#define REALLOC(P,T,S) { if((P=(T*)realloc(P,sizeof(T)*(size_t)(S)))==NULL) { \
      fprintf(stderr,"realloc failed.\n");				\
      exit(EXIT_FAILURE); \
    } }

/* allocates an array of chars */
#define MAKE1DCHAR(D,n1) { MALLOC(D,char,n1); }

/* allocates a two-dimensional array of chars */
#define MAKE2DCHAR(D,n1,n2,i) { MALLOC(D,char*,n1); \
    for(i=0;i<n1;i++) {	\
      MALLOC(D[i],char,n2); \
    } }

/* allocates a three-dimensional array of chars */
#define MAKE3DCHAR(D,n1,n2,n3,i,j) { MALLOC(D,char**,n1); \
    for(i=0;i<n1;i++) {	\
      MALLOC(D[i],char*,n2); \
      for (j=0;j<n2;j++) { \
	MALLOC(D[i][j],char,n3); \
      } } }

/* allocates an array of ints */
#define MAKE1DINT(D,n1)	{ MALLOC(D,int,n1); }

/* lower triangular matrices */
#define TLENGTH(N)   N*(N+1)/2
#define MAKETRIANGULAR(D,T,N) { MALLOC(D,T,TLENGTH(N)); }
/* when we know a is larger than b (in loops for example) */
#define TIDX(a, b)   a * ( a + 1 ) / 2  + b
/* else */
#define TIDXL(a, b)  ( a < b ) ? b*(b + 1)/2 + a : a*(a+1)/2 + b
#define FREETRIANGULAR(D) { FREE1D(D); } 

/* allocates an two-dimensional array of ints */
#define MAKE2DINT(D,n1,n2,i) { MALLOC(D,int*,n1); \
    for(i=0;i<n1;i++) {	\
      MALLOC(D[i],int,n2);			\
    } }
/* allocates a three-dimensional array of floats */
#define MAKE3DINT(D,n1,n2,n3,i,j) { MALLOC(D,int**,n1); \
    for(i=0;i<n1;i++) {	\
      MALLOC(D[i],int*,n2); \
      for (j=0;j<n2;j++) { \
	MALLOC(D[i][j],int,n3); \
      } } }


/* allocates an array of floats */
#define MAKE1DFLOAT(D,n1) { MALLOC(D,float,n1); }

/* allocates an two-dimensional array of floats */
#define MAKE2DFLOAT(D,n1,n2,i) { MALLOC(D,float*,n1); \
    for(i=0;i<n1;i++) {	\
      MALLOC(D[i],float,n2); \
    } }

/* allocates a three-dimensional array of floats */
#define MAKE3DFLOAT(D,n1,n2,n3,i,j) { MALLOC(D,float**,n1); \
    for(i=0;i<n1;i++) {	\
      MALLOC(D[i],float*,n2); \
      for (j=0;j<n2;j++) { \
	MALLOC(D[i][j],float,n3); \
      } } }

/* allocates an array of doubles */
#define MAKE1DDOUBLE(D,n1) { MALLOC(D,double,n1); }

/* allocates an two-dimensional array of doubles */
#define MAKE2DDOUBLE(D,n1,n2,i)	{ MALLOC(D,double*,n1); \
    for(i=0;i<n1;i++) {	\
      MALLOC(D[i],double,n2); \
    } }

#define MAKE1DBOOL(D,n1) { MALLOC(D,bool,n1); }
#define MAKE2DBOOL(D,n1,n2,i)	{ MALLOC(D,bool*,n1); \
    for(i=0;i<n1;i++) {	\
      MALLOC(D[i],bool,n2); \
    } }

/* frees an array of arbitrary type */
#define FREE(D)	{ free(D); D=NULL; }

/* frees an array of arbitrary type */
#define FREE1D(D) { free(D); D=NULL; }

/* frees an two-diemensional array of arbitrary type */
#define FREE2D(D,n1,i) { for(i=0;i<n1;i++) free(D[i]); \
    free (D); D=NULL; }

/* frees an three-diemensional array of arbitrary type */
#define FREE3D(D,n1,n2,i,j) { for(i=0;i<n1;i++) { \
      for (j=0;j<n2;++j) free(D[i][j]);	\
      free(D[i]); }			\
    free (D); D=NULL; }

/* returns the minimum of two values */
#define MIN(a,b) (((a)<=(b))?(a):(b))

/* returns the maximum of two values */
#define MAX(a,b) (((a)>=(b))?(a):(b))

/* returns the minimum of three values */
#define MIN3(a,b,c) (((a)<=(b))?(((a)<=(c))?(a):(c)):((b)<=(c))?(b):(c))

/* returns the maximum of three values */
#define MAX3(a,b,c) (((a)>=(b))?(((a)>=(c))?(a):(c)):((b)>=(c))?(b):(c))

/* returns the number of the minimum argument of two values */
#define MINARG(a,b) (((a)<=(b))?0:1)

/* returns the number of the maximum argument of two values */
#define MAXARG(a,b) (((a)>=(b))?0:1)

/* returns the number of the minimum of three values */
#define MINARG3(a,b,c) (((a)<=(b))?(((a)<=(c))?0:2):((b)<=(c))?1:2)

/* returns the number of the maximum of three values */
#define MAXARG3(a,b,c) (((a)>=(b))?(((a)>=(c))?0:2):((b)>=(c))?1:2)

/* prints an error message to stderr and terminate execution */
#include "vtprogressbar.h"
#define FATAL(msg) { fprintf(stderr,"Error: %s\n\n",msg); VTPROGRESSBARcursorVisible(); exit (EXIT_FAILURE); }

#if HAVE_CONFIG_H
#include "../config.h"
#else
#define DEBUG 0
#endif

/* error message for internal errors (bugs) */
#define FATALINT(msg) { fprintf(stderr,"Internal Error: %s\n\nPlease report bug to: %s\n\n",msg, PACKAGE_BUGREPORT); exit (EXIT_FAILURE); }

/* prints a warning message to stderr */
#define WARN(msg) { fprintf(stderr,"Warning: %s\n\n",msg); }

/* swaps two items */
#define SWAP(a,b,tmp) { tmp=a; a=b; b=tmp; }

/* swap two items in an array */
#define SWAP_TUPLE(t,e1,e2,tmp) { tmp=t[e1]; t[e1]=t[e2]; t[e2]=tmp; }

#include "../libdir/dcmt/include/dc.h"
#include "vtprogressbar.h"
/* init random number generator with seed A for N threads */
#define SRAND(A, I, N) \
    if (Options.Verbosity > 0) \
        VTPROGRESSBARupdate("Initializing Mersenne Twister", N, 0);\
    MALLOC(Options._MTS, mt_struct*, N); \
    for (I=0;I<N;I++) { \
        if (Options.Verbosity > 0) \
            VTPROGRESSBARupdate("Initializing Mersenne Twister", N, I);\
        Options._MTS[I] = get_mt_parameter_id_st(31,521,(I<9?I:I+1),A); \
        if (Options._MTS[I] == NULL) FATALINT("macros.h: SRAND"); \
        sgenrand_mt(I*A,Options._MTS[I]); \
    }\
        if (Options.Verbosity > 0) \
            VTPROGRESSBARcompleted("Initializing Mersenne Twister");

/* clean up mersenne twister */
#define RANDDESTROY(I, N) for (I=0; I<N; I++) free_mt_struct(Options._MTS[I]); \
    FREE(Options._MTS);

#define RAND_MAX_MT 2147483648
#ifdef HAVE_OPENMP
#include <omp.h>
/* return random int in interval [0, A-1] */
#define RANDINT(A) ((int)( A * ( genrand_mt( Options._MTS[ omp_get_thread_num() ] ) /(RAND_MAX_MT+1.0))))
/* return random double in interval [0, A) */
#define RANDDBL(A) ((double)(genrand_mt( Options._MTS[ omp_get_thread_num() ] ) / (((double)RAND_MAX_MT + 1) / (A))))
/* return random double in interval [0, 1) */
#define RANDDBLONE ((double)(genrand_mt( Options._MTS[ omp_get_thread_num() ] ) / ((double)RAND_MAX_MT + 1)))
#else
#define RANDINT(A) ((int)( A * ( genrand_mt( Options._MTS[ 0 ] ) /(RAND_MAX_MT+1.0))))
#define RANDDBL(A) ((double)(genrand_mt( Options._MTS[ 0 ] ) / (((double)RAND_MAX_MT + 1) / (A))))
#define RANDDBLONE ((double)(genrand_mt( Options._MTS[ 0 ] ) / ((double)RAND_MAX_MT + 1)))
#endif

#define RANDNORM(sig, p) do {   \
    double v1,v2,rsq,fac;         \
    do {                          \
        v1=RANDDBL(2.)-1.0;       \
        v2=RANDDBL(2.)-1.0;       \
        rsq=v1*v1+v2*v2;          \
    } while (rsq >= 1.0 || rsq == 0.0); \
    fac=sqrt(-2.0*log(rsq)/rsq);        \
    p = sig*v2*fac;                 \
} while(0);

/* opens a file for writing */
#define FOPENW(fp, fn) if ((fp = fopen(fn, "w+")) == NULL) \
            FATAL("Cannot open file for writing!");
/* opens a file for appending */
#define FOPENA(fp, fn) if ((fp = fopen(fn, "a")) == NULL) \
            FATAL("Cannot open file for appending!");
/* closes a filehandle, this is really 0, not NULL! */
#define FCLOSE(fp) if ((fclose(fp)) != 0) \
            FATAL("Cannot close file!");
#define FFLUSH(fp) if ((fflush(fp)) != 0) \
            FATAL("Cannot flush file!");

#include <float.h>
#include <math.h>

/* compare two floats */
#define CMP_FLOAT(fa, fb) ((fabsf((fa)-(fb)) < FLT_EPSILON)?1:0)
/* compare two doubles */
#define CMP_DBL(fa, fb) ((fabs((fa)-(fb)) < FLT_EPSILON)?1:0)

#if __INTEL_COMPILER
/* Disable ICC's remark #1419: external declaration in primary source file       *
 * This is legal ANSI C code so we disable the remark that is turned on with /W4 */
#pragma warning ( disable : 1419 )
#endif

/* for the banchmarking of hash functions in uthash.h when compiled with
 * HASH_EMIT_KEYS */
#ifdef HASH_EMIT_KEYS
#include <unistd.h>
#endif

/* 
  vim: ft=c sw=4 ts=4 expandtab
*/

#endif
