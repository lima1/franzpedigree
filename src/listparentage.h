#ifndef LISTPARENTAGE_H
#define LISTPARENTAGE_H

/*
 * $Id: listparentage.h 1648 2009-06-09 18:40:00Z markus $
 *
 * Markus Riester, markus@bioinf.uni-leipzig.de */

#include <stdio.h>
#include "ll.h"

typedef struct nst {
    int v;
    int w;
    double p;
    double lod;
    struct nst *next;
} LISTPARENTAGENODE;

void LISTPARENTAGEadd(LISTPARENTAGENODE **, int, int,double, double);
LISTPARENTAGENODE **LISTPARENTAGEsearch(LISTPARENTAGENODE **, int, int);

/* standard functions on lists */
LLAPI(LISTPARENTAGENODE, LISTPARENTAGE)

#endif

/* 
  vim: ft=c sw=4 ts=4 expandtab
*/
