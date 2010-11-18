#ifndef LISTDAG_H
#define LISTDAG_H

/*
 * $Id: listdag.h 1478 2009-04-21 11:31:28Z markus $
 *
 * Markus Riester, markus@bioinf.uni-leipzig.de */

#include <stdio.h>
#include "ll.h"

typedef struct nsdg {
    int v;
    struct nsdg *next;
} LISTDAGNODE;

void LISTDAGadd(LISTDAGNODE **, int);
LISTDAGNODE **LISTDAGsearch(LISTDAGNODE **, int);

/* standard functions on lists */
LLAPI(LISTDAGNODE, LISTDAG)

#endif

/* 
  vim: ft=c sw=4 ts=4 expandtab
*/
