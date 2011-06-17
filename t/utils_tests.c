/*
 * $Id: utils_tests.c 1678 2009-06-19 12:39:16Z markus $
 */

#include "macros.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "tap.h"
#include "utils.h"
#include "pvalue.h"

#define EPSILON 0.0001

int
main(int argc, char *argv[])
{
    double fret, var, avg;
    double farray[] = {0,1,2,3,4,5,5,5,4,3,2,1,0};
    int i, farray_2[] = {10,10,10,10,9,11};
    float farray_3[] = {1,1,1,1,1,1};
    double p[]          = { 0.5677791, 0.6668058, 0.1914607, 0.3691475, 0.3456122, 0.9614572, 0.9934841, 0.6397942, 0.02237115, 0.2144817 };
    double p2[]         = { 0.05677791, 0.06668058, 0.01914607, 0.03691475, 0.03456122, 0.09614572, 0.09934841, 0.06397942, 0.002237115, 0.02144817};
    double p_adjusted[] = { 0.5677791, 0.6668058, 0.1914607, 0.3691475, 0.3456122, 0.9614572, 0.9934841, 0.6397942, 0.02237115, 0.2144817 };

    char test_string_1[] = "the quick brown fox jumps over the lazy dog\n";
    char test_string_2[] = "the quick brown fox jumps over the lazy dog";
    char test_string_3[] = " 123    quick brown fox    ";
    char test_string_4[] = "               ";
    char test_string_5[] = "quick brown fox";
    char test_string_6[] = "_123____quick_brown_fox____";
    int *D;

    PVALUE_LOOKUP pv_lup;

    plan_tests(44);
    
    UTILSfactlnInit(2000);
    fret = UTILSfactln(0);
    ok(fabs(fret - 0.) <  EPSILON, "log(0!) == 0");
    fret = UTILSfactln(1);
    ok(fabs(fret - 0.) <  EPSILON, "log(1!) == 0");
    fret = UTILSfactln(2);
    ok(fabs(fret - 0.6931472) <  EPSILON, "log(2!) == 0.6931472");
    fret = UTILSfactln(100);
    ok(fabs(fret - 363.7394) <  EPSILON, "log(100!) == 363.7394");
   
    ok(MINARG(3,2) == 1, "minarg macro correct");
    ok(MINARG(-3,2) == 0, "minarg macro correct");

    /* test variance with just one element */
    UTILScalcVarianceArray(farray,2,2,&var, &avg);
    ok(CMP_DBL(avg, farray[2]), "avg = element");
    ok(CMP_DBL(var, 0), "var = 0");
    /* and with two elements*/
    UTILScalcVarianceArray(farray,2,3,&var, &avg);
    ok(CMP_DBL(avg, 2.5), "avg = 2.5");
    ok(CMP_DBL(var, 0.5), "mean = 0.5");

    MAKETRIANGULAR(D,int,10);
    D[TIDX(1,2)] = 3;
    D[TIDX(8,9)] = 1;
    ok(D[TIDXL(2,1)] = 3, "triangular ok");
    ok(D[TIDXL(8,9)] = 1, "triangular ok");
    FREE1D(D);

    /* and with more elements */
    UTILScalcVarianceArray(farray,0,12,&var, &avg);
    ok(fabs(var -  3.397436) < EPSILON, "var ==  3.397436");
    ok(fabs(avg -  2.692308) < EPSILON, "avg ==  2.692308");
    UTILScalcVarianceArray(farray,2,11,&var, &avg);
    ok(fabs(var -  2.044444) < EPSILON, "var ==  2.044444");
    ok(fabs(avg -  3.4) < EPSILON, "avg ==  3.4");

    fret = UTILSbeta(2,10);
    ok(CMP_DBL(fret, 0.009090909), "beta function seems to work");

    fret = UTILScalcDistance(51.28,12.35,51.3,12.38);
    ok(fabs(fret - 3049.5396) < EPSILON, "distance connewitz, cospudener correct");

    fret = UTILSchoose(100,2);
    ok(fabs(fret - 4950) < EPSILON, "100 over 2 is 45");

    fret = UTILSchoose(174,4);
    ok(CMP_DBL(fret, 36890001), "174 over 4 is 36890001");

    fret = exp(UTILSbinln(5,0,0.2));
    ok(fabs(fret - 0.32768) < EPSILON, "bin correct");
    fret = exp(UTILSbinln(5,1,0.2));
    ok(fabs(fret - 0.4096) < EPSILON, "bin correct");
    fret = exp(UTILSbinln(5,2,0.2));
    ok(fabs(fret - 0.20480) < EPSILON, "bin correct");
    fret = exp(UTILSbinln(5,3,0.2));
    ok(fabs(fret - 0.0512) < EPSILON, "bin correct");
    fret = exp(UTILSbinln(5,4,0.2));
    ok(fabs(fret - 0.0064) < EPSILON, "bin correct");
    fret = exp(UTILSbinln(5,5,0.2));
    ok(fabs(fret - 0.00032) < EPSILON, "bin correct");
    
    for (i=0; i<6; i++)
        farray_3[i] /= 6;
    
    fret = exp(UTILSmultiln(60,farray_2, farray_3,6));
    ok(fabs(fret - 0.000068) < EPSILON, "multi correct");
    farray_2[4] = 10;
    farray_2[5] = 10;
    fret = exp(UTILSmultiln(60,farray_2, farray_3,6));
    ok(fabs(fret - 0.000075) < EPSILON, "multi correct");

    /*fprintf(stderr, "DEBUG: %f\n", fret); */
    UTILSstripNewline(test_string_1, strlen(test_string_1));
    ok(strcmp(test_string_1, test_string_2) == 0,"newline removed");
    UTILScopyNonBlanks(test_string_3, test_string_4,4,strlen(test_string_3));
    ok(strcmp(test_string_4, test_string_5) == 0,"leading and trailing blanks removed");
/*    fprintf(stderr, "XXX%sXXX\n",test_string_1);
    fprintf(stderr, "XXX%sXXX\n",test_string_2);*/
    UTILSreplaceSpace(test_string_3, '_',strlen(test_string_3));
    ok(strcmp(test_string_3, test_string_6) == 0,"blanks replaced");
   
    UTILSfactlnDestroy();
    UTILSpvAdjustBH(p, p_adjusted, 10);
    ok(fabs(p_adjusted[0] - 0.9934841) < EPSILON, "BH correct");
    ok(fabs(p_adjusted[9] - 0.2237115) < EPSILON, "BH correct");
    ok(fabs(p_adjusted[4] - 0.8335072) < EPSILON, "BH correct");

    UTILSpvAdjustHolm(p2, p_adjusted, 10);
    ok(fabs(p_adjusted[0] - 0.02237115) < EPSILON, "holm correct");
    ok(fabs(p_adjusted[9] - 0.28388955) < EPSILON, "holm correct");
    ok(fabs(p_adjusted[4] - 0.24192854) < EPSILON, "holm correct");

    pv_lup = PVALUEinit(20, -7.1, 7.92);
    //PVALUEdump(stderr, pv_lup);
    PVALUEadd(pv_lup, 82, 0.001, 1.);
    PVALUEadd(pv_lup, -8, 0.9998, 1.);
    ok(CMP_DBL(PVALUEfind(pv_lup, 10), 0.001), "too large value ok");   
    ok(CMP_DBL(PVALUEfind(pv_lup, -855), 0.9998), "too small value ok");   
    PVALUEadd(pv_lup, 6, 0.01, 1.);
    ok(CMP_DBL(PVALUEfind(pv_lup, 6.005), 0.01), "in range value ok");   
    PVALUEadd(pv_lup, -6, 0.82, 1.);
    //PVALUEdump(stderr, pv_lup);
    ok(CMP_DBL(PVALUEfind(pv_lup, 7.005), 1.), "value wrong because not yet filled");   
    PVALUEfill(pv_lup);
    ok(CMP_DBL(PVALUEfind(pv_lup, 7.005), 0.01), "now correct after filling");   
    ok(CMP_DBL(PVALUEfindCritical(pv_lup, 0.005), 7.92), "find critical value correct");   
    ok(CMP_DBL(PVALUEfindCritical(pv_lup, 0.5), 5.548421), "find critical value correct");   
//    PVALUEdump(stderr, pv_lup);
    PVALUEdestroy(pv_lup);

    return exit_status();       /* Return the correct exit code */
}

/* 
  vim: ft=c sw=4 ts=4 expandtab
*/
