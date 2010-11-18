/*
    PVALUEdestroy(pv_lup);
 * $Id: pvalue.c 1977 2010-03-16 16:52:31Z markus $
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

#include "macros.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "pvalue.h"

#define D_MAX_VAL 9999.

typedef struct S_PVALUE
{
    double v;      /* the observed value */
    double pv;     /* the corresponding pvalue */
    double sensitivity;
} PVALUE;

struct S_PVALUE_LOOKUP
{
    int length;   /* how many values in lookup? */
    double min;   /* smalles and                */
    double max;   /*              largest value */
    double  _d;
    PVALUE *pvalues;
};

PVALUE_LOOKUP
PVALUEinit(int length, double min, double max)
{
    int i;
    PVALUE_LOOKUP pv_lup = malloc(sizeof *pv_lup);
    if (pv_lup == NULL)
        FATAL("malloc failed");

    if (length < 2)
        FATALINT("PVALUEinit: length");
    
//    fprintf(stderr, "min: %f, max %f MAX %f\n",min,max,D_MAX_VAL);
    /* to avoid the cases max is inf */
    max = MIN(max, D_MAX_VAL);

    if (min >= max)
        FATALINT("PVALUEinit: min/max");

    pv_lup->length = length;
    pv_lup->min    = min;
    pv_lup->max    = max;

    MALLOC(pv_lup->pvalues, PVALUE, length);
    pv_lup->_d = ( max - min ) / (double)(length - 1); 

    pv_lup->pvalues[0].v = min;
    for (i=1; i<length; i++) {
        pv_lup->pvalues[i].v  = pv_lup->pvalues[i-1].v + pv_lup->_d;
        pv_lup->pvalues[i].pv = 1.;
    }    

    assert(CMP_DBL(pv_lup->pvalues[length-1].v, max));

    return pv_lup;
}

static inline int
findElement(PVALUE_LOOKUP pv_lup, double v)
{    
    int idx;

    if (v < pv_lup->min)
        v = pv_lup->min;
    if (v > pv_lup->max)
        v = pv_lup->max;
    
    idx = (int)((v-pv_lup->min)/pv_lup->_d);

    assert(idx >= 0 && idx < pv_lup->length);

    return idx;
}

double 
PVALUEfind(PVALUE_LOOKUP pv_lup, double v)
{
    int idx = findElement(pv_lup, v);

    return pv_lup->pvalues[idx].pv;
}

void
PVALUEadd(PVALUE_LOOKUP pv_lup, double v, double pv, double sensitivity)
{    
    int idx = findElement(pv_lup, v);

    assert(pv >= 0. && pv <= 1.);

    pv_lup->pvalues[idx].pv = pv;
    pv_lup->pvalues[idx].sensitivity = sensitivity;
}

double
PVALUEfindCritical(PVALUE_LOOKUP pv_lup, double pv)
{
    int i;
    double ret = pv_lup->pvalues[pv_lup->length - 1].v;
    for (i=0; i<pv_lup->length; i++) {
        if ( pv_lup->pvalues[i].pv < pv) {
            ret = pv_lup->pvalues[i].v;
            break;
        }    
    }
    return ret;
}

void
PVALUEfill(PVALUE_LOOKUP pv_lup)
{
    int i;

    for (i=1; i<pv_lup->length; i++)
        if (pv_lup->pvalues[i].pv > pv_lup->pvalues[i-1].pv) {
            pv_lup->pvalues[i].pv = pv_lup->pvalues[i-1].pv;
            pv_lup->pvalues[i].sensitivity = pv_lup->pvalues[i-1].sensitivity;
        }            
}

void
PVALUEdump(FILE *fp, PVALUE_LOOKUP pv_lup, int grid)
{
    int i, inc = 1, length = 0;
//    fprintf(fp, "***\nPVALUE LOOPKUP, size: %i, grid: %f\n\n", pv_lup->length, pv_lup->_d);
    fprintf(fp, "%-10s  %-10s  %-11s\n","Delta", "p-Value", "Sensitivity");

    if (grid > 0) {
        for (i=0; i<pv_lup->length; i+=inc)
            if (!CMP_DBL(pv_lup->pvalues[i].pv, 0.)) length++;
        inc = MAX(1,length / grid);
    }
    else {
        inc = 1;
        length = pv_lup->length;
     }    

    for (i=0; i<length; i+=inc)
        fprintf(fp, "%10.5f  %10.5f  %11.5f\n", pv_lup->pvalues[i].v, pv_lup->pvalues[i].pv, pv_lup->pvalues[i].sensitivity);

    fprintf(fp, "\n");
}

void
PVALUEdestroy(PVALUE_LOOKUP pv_lup)
{
    FREE1D(pv_lup->pvalues);
    FREE(pv_lup);
}    
/* 
  vim: ft=c sw=4 ts=4 expandtab
*/
