/*
 * $Id: listdag.c 1478 2009-04-21 11:31:28Z markus $
 *
 * Markus Riester, markus@bioinf.uni-leipzig.de */

#include "macros.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "listdag.h"

void
LISTDAGadd(LISTDAGNODE ** p, int v)
{
    LISTDAGNODE *n = (LISTDAGNODE *) malloc(sizeof(LISTDAGNODE));

    if (n == NULL) {
        FATAL("Malloc failed");
    }    

    n->next = *p;
    *p = n;
    n->v = v;
}

/* for LISTprint */
static void
dumpNode(FILE *fp, LISTDAGNODE * n)
{
    fprintf(fp,"%d; ", n->v);
}

/* for LISTsort */
static int
cmp(LISTDAGNODE * a, LISTDAGNODE * b)
{
    if (a->v > b->v) {
        return -1;
    } else if (a->v < b->v) {
        return 1;
    }
    return 0;
}

LISTDAGNODE **
LISTDAGsearch(LISTDAGNODE ** n, int v)
{
    while (*n != NULL) {
        if ((*n)->v == v) {
            return n;
        }
        n = &(*n)->next;
    }
    return NULL;
}

/* generate standard functions on lists */
LLCODE(LISTDAG, LISTDAGNODE)

/* 
  vim: ft=c sw=4 ts=4 expandtab
*/
