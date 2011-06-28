/*
 * $Id: ll_tests.c 1667 2009-06-16 11:03:01Z markus $
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tap.h"
#include "listdag.h"

int
main(int argc, char *argv[])
{
    LISTDAGNODE *t=NULL;

    plan_tests(7);
    ok(LISTDAGcnt(&t) == 0, "Empty list");
    LISTDAGadd(&t,1);
    ok(LISTDAGcnt(&t) == 1, "One element in list");
    LISTDAGadd(&t,4);
    LISTDAGadd(&t,2);
    ok(LISTDAGcnt(&t) == 3, "Three elements in list");
    t=LISTDAGsort(t);
    ok(t->v == 4, "4 largest value in list");
    LISTDAGremove(LISTDAGsearch(&t,2));
    ok(LISTDAGcnt(&t) == 2, "Two elements in list");
    ok(LISTDAGsearch(&t,2) == NULL, "2 not in list anymore");
    LISTDAGdestroy(&t);
    ok(t == NULL, "empty");
    return exit_status();       /* Return the correct exit code */
}

/* 
  vim: ft=c sw=4 ts=4 expandtab
*/
