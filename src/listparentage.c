/*
 * $Id: listparentage.c 1648 2009-06-09 18:40:00Z markus $
 *
 * Markus Riester, markus@bioinf.uni-leipzig.de */

#include <macros.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "listparentage.h"

void
LISTPARENTAGEadd(LISTPARENTAGENODE ** p, int v, int w, double prob, double lod)
{
    LISTPARENTAGENODE *n = (LISTPARENTAGENODE *) malloc(sizeof(LISTPARENTAGENODE));

    if (n == NULL) {
        FATAL("Malloc failed");
    }    
    n->next = *p;
    *p = n;
    n->v = v;
    n->w = w;
    n->p = prob;
    n->lod = lod;
}

/* for LISTprint */
static void
dumpNode(FILE *fp, LISTPARENTAGENODE * n)
{
    fprintf(fp,"%d %d (p: %f, LOD: %f);", n->v, n->w, n->p, n->lod);
}

/* for LISTsort */
static int
cmp(LISTPARENTAGENODE * a, LISTPARENTAGENODE * b)
{
    if (a->lod > b->lod) {
        return -1;
    } else if (a->lod < b->lod) {
        return 1;
    }
    return 0;
}

LISTPARENTAGENODE **
LISTPARENTAGEsearch(LISTPARENTAGENODE ** n, int v, int w)
{
    while (*n != NULL) {
        if (((*n)->v == v && (*n)->w == w) || ((*n)->v == w && (*n)->w == v)) {
            return n;
        }
        n = &(*n)->next;
    }
    return NULL;
}

/* generate standard functions on lists */
LLCODE(LISTPARENTAGE, LISTPARENTAGENODE)

/* 
  vim: ft=c sw=4 ts=4 expandtab
*/
