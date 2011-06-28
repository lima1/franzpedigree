/*
 * $Id: utils.c 1977 2010-03-16 16:52:31Z markus $
 *
 * Some code taken from the Numerical Recipes bible.
 *
 * Copyright (C) 2008-2010 Universitaet Leipzig  
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Markus Riester, markus@bioinf.uni-leipzig.de */

/* just for FATAL */
#include "macros.h"
#include <math.h>
/* for UTILScopyNonBlanks */
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>

/* for UTILScalcDistance and UTILScalcLogFac*/

#include <string.h>

#include "utils.h"
#include "pvalue.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846 
#endif

/* calculates variance and mean of an array in the specified range (including
 * start, excluding end */
void
UTILScalcVarianceArray(double *arr, unsigned int start, unsigned int end, double *var,
                       double *avg)
{
    unsigned int i, n = 0;
    double mean = 0., M2 = 0., delta;

    if (start >= end) {
        *avg = arr[start];
        *var = 0;
        return;
    }

    for (i = start; i <= end; i++) {
        n++;
        delta = arr[i] - mean;
        mean += delta/(double)n;
        M2 += delta*(arr[i] -mean);
    }

    *avg = mean;
    *var = M2/(double)(n-1);
}

void
UTILScalcVarianceArrayInt(int *arr, unsigned int start, unsigned int end, double *var,
                       double *avg)
{
    unsigned int i, n = 0;
    double mean = 0., M2 = 0., delta;

    if (start >= end) {
        *avg = arr[start];
        *var = 0;
        return;
    }

    for (i = start; i <= end; i++) {
        n++;
        delta = arr[i] - mean;
        mean += delta/(double)n;
        M2 += delta*(arr[i] -mean);
    }

    *avg = mean;
    *var = M2/(double)(n-1);
}

int
UTILScompare_doubles_decr (const void *a, const void *b)
{
    const double *da = (const double *) a;
    const double *db = (const double *) b;

    return (*db > *da) - (*db < *da);
}

int
UTILScompare_doubles_incr (const void *a, const void *b)
{
    const double *da = (const double *) a;
    const double *db = (const double *) b;

    return (*da > *db) - (*da < *db);
}

int
UTILScompare_ints_incr (const void *a, const void *b)
{
    int *arg1 = (int *)a;
    int *arg2 = (int *)b;

    if (*arg1 < *arg2)
        return -1;
    else if (*arg1 == *arg2)
        return 0;
    else
        return 1;
}


void
UTILScmpDistributions(double *a1, double *a2, size_t n1, size_t n2, PVALUE_LOOKUP *pv_lup, int w1)
{
    size_t i, j, latest_j = 0, start_j = 0;
    bool latest_j_found = false;
    double ratio;
    qsort(a1, n1, sizeof(double), UTILScompare_doubles_decr);
    qsort(a2, n2, sizeof(double), UTILScompare_doubles_decr);

    /* no "bad" values */
    if (n1 == 0 && n2 > 0) { 
        PVALUEadd(*pv_lup, a2[n2 - 1], 0.,1.0);
        *pv_lup = PVALUEinit(20000, a2[n2 - 1], a2[0]);
    } else {
        if (n2 == 0) FATALINT("cmpDistributions");
        *pv_lup = PVALUEinit(20000, a1[n1 - 1], a2[0]);
    }

    for (i = 0; i < n1; i++) {
        
        ratio = 0.; 
        for (j=start_j; j < n2; j++) {
            /* what is the first value that is smaller than a1? */
            if (a2[j] < a1[i] || j == n2 - 1) {
                if (i==0 && j>0) { 
                    latest_j = j-1;
                    latest_j_found = true;
                }    
                ratio = (j + 1) / (double)(w1 * i + j + 2);
    //            fprintf(stderr, "DEBUG: %i %i %f %f ratio: %f\n", i,j, a1[i],a2[j],ratio);
                break;
            }
            /* where can we start to look for the first value smaller than
             * a1[i+1]? So we don't have O(n^2)? */ 
            if (i < (n1-1) && a2[j] >= a1[i+1]) 
               start_j = j;

        }
        PVALUEadd(*pv_lup, a1[i], 1. - ratio,((j+1)/(double)n2));
    }
    /* now, all values are correct (better than the best wrong (a1[n1-1])) */
    if (latest_j_found) {
        for (i=latest_j; i > 0; --i) {
           PVALUEadd(*pv_lup, a2[i], 0., (i+1)/(double)n2);
        }   
    }
    PVALUEfill(*pv_lup);
    //PVALUEdump(stderr,*pv_lup,10);
}

double
UTILSdeg2rad(double deg)
{
    return (deg * M_PI / 180.);
}

double
UTILScalcDistance(double lat1, double lon1, double lat2, double lon2)
{
    double dlon, dlat, a, slon, slat;

    lat1 = UTILSdeg2rad(lat1);
    lat2 = UTILSdeg2rad(lat2);
    lon1 = UTILSdeg2rad(lon1);
    lon2 = UTILSdeg2rad(lon2);

    dlon = lon1 - lon2;
    dlat = lat1 - lat2;
    slat = sin(dlat / 2.);
    slon = sin(dlon / 2.);

    a = slat * slat + ( cos(lat1) * cos(lat2) * slon * slon );
    /*return 2 * METER_RHO * atan2( sqrt(a), sqrt(1.-a) ); */
    return 2 * METER_RHO * asin(sqrt(a));
}

double
UTILSbinln(int N, int k, double p)
{
    if (N <= 0 || k < 0 || k > N || p < 0. || p > 1.)
        FATALINT("Bad argument in UTILSbinln()");
    return k*log(p)+(N-k)*log(1.-p)
            +UTILSfactln(N)-UTILSfactln(k)-UTILSfactln(N-k);
}

double
UTILSmultiln(int N, int *k, float *p, int q)
{
    int i, j;
    double sumpi = 0.;

    /* calc p_i^k[i]/k[i]! */
    for (i = 0; i < q; i++) {
        for (j = 0; j < k[i]; j++) {
            sumpi += log(p[i]);
        }
        sumpi -= UTILSfactln(k[i]);
    }
    return UTILSfactln(N) + sumpi;
}

double
UTILSgammln(double xx)
{
    int j;
    double x, tmp, y, ser;
    static const double cof[14] = { 57.1562356658629235, -59.5979603554754912,
        14.1360979747417471, -0.491913816097620199, .339946499848118887e-4,
        .465236289270485756e-4, -.983744753048795646e-4,
            .158088703224912494e-3,
        -.210264441724104883e-3, .217439618115212643e-3,
            -.164318106536763890e-3,
        .844182239838527433e-4, -.261908384015814087e-4,
            .368991826595316234e-5
    };
    if (xx <= 0)
        FATALINT("Bad arg in UTILSgammln()");
    y = x = xx;
    tmp = x + 5.24218750000000000;
    tmp = (x + 0.5) * log(tmp) - tmp;
    ser = 0.999999999999997092;
    for (j = 0; j < 14; j++)
        ser += cof[j] / ++y;
    return tmp + log(2.5066282746310005 * ser / x);
}

double
UTILSfactrl(int n)
{
    int i;
    static double a[171];
    static bool init = true;

    if (init) {
        init = false;
        a[0] = 1.;
        for (i = 1; i < 171; i++)
            a[i] = i * a[i - 1];
    }
    if (n < 0 || n > 170)
        FATALINT("Parameter in UTILSfactrl() out of range");
    return a[n];
}

static double * _fact_lookup;
static int _fact_lookup_size = -1;

void
UTILSfactlnInit(int n)
{
    int i;
    MAKE1DDOUBLE(_fact_lookup, n+1);
    for (i = 0; i <= n; i++)
        _fact_lookup[i] = UTILSgammln(i + 1.);
    _fact_lookup_size = n;
}

double
UTILSfactln(int n)
{
    if (n < 0 || _fact_lookup_size < 0)
        FATALINT("UTILSfactln()");
    if (n <= _fact_lookup_size)
        return _fact_lookup[n];
    WARN("value in factln for lookup table too large, please report bug");
    return UTILSgammln(n + 1.);
}

void
UTILSfactlnDestroy(void)
{
    FREE1D(_fact_lookup);
} 

double
UTILSbeta(double z, double w)
{
    return exp(UTILSgammln(z) + UTILSgammln(w) - UTILSgammln(z + w));
}

double
UTILSchoose(int n, int k)
{
    if (n < 0 || k < 0 || k > n)
        FATALINT("Bad args in UTILSchoose()");
    if (n < 171)
        return floor(0.5 +
                     UTILSfactrl(n) / (UTILSfactrl(k) * UTILSfactrl(n - k)));
    return floor(0.5 +
                 exp(UTILSfactln(n) - UTILSfactln(k) - UTILSfactln(n - k)));
}

void
UTILScopyNonBlanks(char *source, char *dest, unsigned int min, unsigned int max)
{
    unsigned int i, j;

    if (max == 0)
        FATALINT("max must be greater than 0");

    /* copies string from min to max in dest, ignores leading 
     * and trailing whitespaces */
    for (i = min; i < max; i++)
        if (!isspace(source[i]))
            for (j = max - 1; j >= i; j--) {
                /*fprintf(stderr, "%i %i %c %c\n", i, j, source[i], source[j]); */
                if (!isspace(source[j])) {
                    if ((strncpy(dest, (source + i), (size_t) (j - i + 1))) ==
                        NULL)
                        FATAL("Cannot copy string");
                    if ((j-i+1)>0)
                        dest[j - i + 1] = '\0';
                    return;
                }
            }
}

void
UTILSstripNewline(char *str, size_t size)
{
    size_t i;

    /* remove the null terminator */
    for (i = 0; i < size; ++i) {
        if (str[i] == '\n') {
            str[i] = '\0';
            return;
        }
    }
}

void 
UTILSreplaceSpace(char *str, char replace, size_t size) 
{
    size_t i;

    /* remove the null terminator */
    for (i = 0; i < size; ++i) 
        if (isspace(str[i])) 
            str[i] = replace;
}

void 
UTILSremoveSpecialChars(char *str, char *dest, size_t size) 
{
    size_t i,j=0;
    int ascii;

    /* remove the null terminator */
    for (i = 0; i < size; i++) {
        if (str[i] == '\0') {
            dest[j++] = str[i];
            return;
        }    
        ascii = (int)str[i];
        if (ascii < 33 || ascii > 125) continue;
        dest[j++] = str[i];
    }
}

void
UTILSpvAdjustBH(double *p, double *p_adjusted, int n)
{
    int i;
    if (n < 1)
        return;

    qsort(p, (size_t)n, sizeof(double), UTILScompare_doubles_decr);
    
    p_adjusted[0] = p[0];
    for (i=1; i<n;i++) {
        p_adjusted[i] = MIN(p_adjusted[i-1], (n * p[i] / (double)(n-i) ) );
        //fprintf(stderr, "DEBUG %i  %f: %f\n", i, p[i], p_adjusted[i]);
    }
}

void
UTILSpvAdjustHolm(double *p, double *p_adjusted, int n)
{
    int i;
    if (n < 1)
        return;

    qsort(p, (size_t)n, sizeof(double), UTILScompare_doubles_incr);
    

    p_adjusted[0] = MIN(1.0, ( n * p[0] ) );
    for (i=1; i<n;i++) {
        p_adjusted[i] = MIN(1.0, MAX(p_adjusted[i-1], ( (n-i) * p[i] )) );
        //fprintf(stderr, "DEBUG %i  %f: %f\n", i, p[i], p_adjusted[i]);
    }
}

void UTILSshuffleInt(int *array, size_t n)
{
    if (n > 1) {
        size_t i;
    for (i = 0; i < n - 1; i++) {
      size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
      int t = array[j];
      array[j] = array[i];
      array[i] = t;
    }
    }
}

/* 
  vim: ft=c sw=4 ts=4 expandtab
*/
