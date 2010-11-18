/*
 * $Id$
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "tap.h"
#include "global.h"
#include "macros.h"
#include "hwe.h"
#include "utils.h"

OPTIONS Options;
DATA    Data;
PROBS   Probs;

int
main(int argc, char *argv[])
{
    int i,j, ***genotypes, n; 
    int testdata [45][1][2] = {
        { { 2, 3} }, 
        { { 3, 2} }, 
        { { 2, 3} }, 
        { { 3, 3} }, 
        { { 2, 7} }, 
        { { 7, 2} }, 
        { { 2, 7} }, 
        { { 7, 2} }, 
        { { 2, 7} }, 
        { { 3, 7} }, 
        { { 3, 7} }, 
        { { 3, 7} }, 
        { { 3, 7} }, 
        { { 3, 7} }, 
        { { 3, 7} }, 
        { { 3, 7} }, 
        { { 3, 7} }, 
        { { 3, 7} }, 
        { { 7, 3} }, 
        { { 7, 3} }, 
        { { 7, 3} }, 
        { { 7, 3} }, 
        { { 7, 3} }, 
        { { 7, 3} }, 
        { { 7, 3} }, 
        { { 7, 3} }, 
        { { 7, 3} }, 
        { { 7, 7} }, 
        { { 8, 2} }, 
        { { 2, 8} }, 
        { { 8, 2} }, 
        { { 8, 3} }, 
        { { 8, 3} }, 
        { { 8, 3} }, 
        { { 8, 3} }, 
        { { 8, 3} }, 
        { { 3, 8} }, 
        { { 8, 3} }, 
        { { 8, 8} }, 
        { { 8, 8} }, 
        { { 8, 7} }, 
        { { 8, 7} }, 
        { { 8, 7} }, 
        { { 8, 7} }, 
        { { 8, 7} }, 
    },
    testdata2[8][1][2] = {
        { { 8, 7} }, 
        { { 8, 7} }, 
        { { 8, 7} }, 
        { { 8, 7} }, 
        { { 7, 8} }, 
        { { 7, 8} }, 
        { { 7, 8} }, 
        { { 7, 8} }, 
    }     
    ,
    testdata3[8][1][2] = {
        { { 7, 7} }, 
        { { 7, 7} }, 
        { { 8, 8} }, 
        { { 8, 8} }, 
        { { 7, 8} }, 
        { { 7, 8} }, 
        { { 7, 8} }, 
        { { 7, 8} }, 
    },
    testdata4[9][1][2] = {
        { { 1, 1} }, 
        { { 1, 2} }, 
        { { 1, 1} }, 
        { { 2, 2} }, 
        { { 1, 1} }, 
        { { 1, 1} }, 
        { { 1, 1} }, 
        { { 1, 1} }, 
        { { 2, 2} }, 
    }    
    ;

    FREQS freqs;
    double fret;

    HWETEST hwe;
    HWERESULT hwe_result;
    plan_tests(23);
#ifdef HAVE_OPENMP
    n = omp_get_max_threads();
#else
    n=1;
#endif
    
    SRAND(0,i, n);
    /* 2,3,7,8 */
    MAKE3DINT(genotypes, 45, 1, 2, i,j);
    /* normaly, this is done by freq.c */
    MAKE1DINT(freqs.allele_id_lookup, 8-2+1);
    for (i=0; i<(8-2+1); i++)
        freqs.allele_id_lookup[i] = -1;
    UTILSfactlnInit(2000);
    freqs.allele_id_lookup[2-2] = 0;
    freqs.allele_id_lookup[3-2] = 1;
    freqs.allele_id_lookup[7-2] = 2;
    freqs.allele_id_lookup[8-2] = 3;
    freqs.min = 2;
    freqs.max = 8;

    for (i=0; i<45; i++)
        for (j=0; j<2; j++)
            genotypes[i][0][j] = testdata[i][0][j];

    hwe = HWEinit(4);
    HWEcalcFreqs(hwe,genotypes,45,0,4,&freqs);
    ok(HWEgetFreq(hwe,2,2) == 0, "freq correct"); 
    ok(HWEgetFreq(hwe,8,8) == 2, "freq correct"); 
    ok(HWEgetFreq(hwe,3,7) == 18, "freq correct"); 
    ok(HWEgetFreq(hwe,7,3) == 18, "freq correct"); 
    ok(HWEgetAlleleFreq(hwe,2) == 11, "allele freq correct"); 
    ok(HWEgetAlleleFreq(hwe,3) == 30, "allele freq correct"); 
    ok(HWEgetAlleleFreq(hwe,7) == 30, "allele freq correct"); 
    ok(HWEgetAlleleFreq(hwe,8) == 19, "allele freq correct"); 
    //HWEdumpFreqs(stderr, hwe);
    HWEtest(hwe,10000,100,1700);
    hwe_result = HWEgetResult(hwe);
    ok(fabs(hwe_result.p_value - 0.017) < 0.001, "p-value correct");
    ok(fabs(hwe_result.se - 0.) < 0.001, "estimated error correct");
    //HWEdumpResults(stderr, hwe);
    HWEdestroy(hwe);
    hwe = HWEinit(4);
    HWEcalcFreqs(hwe,genotypes,45,0,4,&freqs);
    HWEtest(hwe,20000,250,2500);
    hwe_result = HWEgetResult(hwe);
    ok(fabs(hwe_result.p_value - 0.017) < 0.001, "p-value correct");
    ok(fabs(hwe_result.se - 0.) < 0.001, "estimated error correct");
    //HWEdumpResults(stderr, hwe);
    FREE1D(freqs.allele_id_lookup);
    createSumLog2(hwe);
    fret = exp( calcSampleProb(hwe, 0,8,0) );
    ok(fabs(fret-0.01989) < EPSILON,"sample probability correct");
    fret = exp( calcSampleProb(hwe, 1,6,1) );
    ok(fabs(fret-0.27848) < EPSILON,"sample probability correct");
    fret = exp( calcSampleProb(hwe, 2,4,2) );
    ok(fabs(fret-0.52215) < EPSILON,"sample probability correct");
    fret = exp( calcSampleProb(hwe, 3,2,3) );
    ok(fabs(fret-0.17404) < EPSILON,"sample probability correct");
    fret = exp( calcSampleProb(hwe, 4,0,4) );
    ok(fabs(fret-0.00544) < EPSILON,"sample probability correct");
    HWEdestroy(hwe);
    hwe = HWEinit(2);
    FREE3D(genotypes,45,1,i,j);
    MAKE3DINT(genotypes, 8, 1, 2, i,j);
    MAKE1DINT(freqs.allele_id_lookup, 2);
    freqs.allele_id_lookup[0] = 0;
    freqs.allele_id_lookup[1] = 1;
    freqs.min= 7;
    freqs.max = 8;

    for (i=0; i<8; i++)
        for (j=0; j<2; j++)
            genotypes[i][0][j] = testdata2[i][0][j];

    HWEcalcFreqs(hwe,genotypes,8,0,2,&freqs);
    HWEtest(hwe,20000,250,2500);

    hwe_result = HWEgetResult(hwe);

    ok(fabs(hwe_result.p_value - 0.02533) < EPSILON, "p-value correct");
    ok(fabs(hwe_result.se - 0.) < EPSILON, "estimated error correct");
    for (i=0; i<8; i++)
        for (j=0; j<2; j++)
            genotypes[i][0][j] = testdata3[i][0][j];

    HWEcalcFreqs(hwe,genotypes,8,0,2,&freqs);
    HWEtest(hwe,20000,250,2500);
    hwe_result = HWEgetResult(hwe);
    ok(fabs(hwe_result.p_value - 1.0) < EPSILON, "p-value correct");
    ok(fabs(hwe_result.se - 0.) < EPSILON, "estimated error correct");
    HWEdestroy(hwe);
    FREE1D(freqs.allele_id_lookup);
    FREE3D(genotypes,8,1,i,j);
    /* the HWE.exact R example */
    MAKE3DINT(genotypes, 9, 1, 2, i,j);
    /* normaly, this is done by freq.c */
    MAKE1DINT(freqs.allele_id_lookup, 2-1+1);
    for (i=0; i<(2-1+1); i++)
        freqs.allele_id_lookup[i] = -1;
    freqs.allele_id_lookup[0] = 0;
    freqs.allele_id_lookup[1] = 1;
    freqs.min = 1;
    freqs.max = 2;

    for (i=0; i<9; i++)
        for (j=0; j<2; j++)
            genotypes[i][0][j] = testdata4[i][0][j];

    hwe = HWEinit(2);
    HWEcalcFreqs(hwe,genotypes,9,0,2,&freqs);
    HWEtest(hwe,10000,100,1700);
    hwe_result = HWEgetResult(hwe);
    ok(fabs(hwe_result.p_value - 0.05882) < 0.0001, "p-value correct");
    ok(fabs(hwe_result.se - 0.) < EPSILON, "estimated error correct");
/*    HWEdumpFreqs(stderr,hwe);
    HWEdumpResults(stderr,hwe);*/

    HWEdestroy(hwe);
    FREE1D(freqs.allele_id_lookup);
    FREE3D(genotypes,9,1,i,j);

    UTILSfactlnDestroy();
    RANDDESTROY(i,n);
    return exit_status();       /* Return the correct exit code */
}

/* 
  vim: ft=c sw=4 ts=4 expandtab
*/
