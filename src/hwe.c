/*
 * $Id: hwe.c 2058 2010-05-12 15:19:58Z markus $
 *
 * Mostly taken from the hwe program, by Guo (now GPL in PyPop).
 *
 * Copyright (C) 1992. Sun-Wei Guo.
 * Copyright (C) 2008-2010. Universitaet Leipzig  
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

#define _POSIX_SOURCE
#include "macros.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "hwe.h"
#include "global.h"
#include "utils.h"

#define  RATIO(u, v)   ( (double) (u) ) / ( 1.0 + (double) (v) )
#define  TRANS(x)     (MIN(1.0, x))/2.0 /* transition probability */

extern OPTIONS Options;

struct S_HWETEST
{
    int *matrix;         /* array containing a lower nxn matrix              */
    int *matrix_orig;    /* array containing a lower nxn matrix, for output  */
    int matrix_length;   /* the length of the array                          */
    int *allele_freq;     
    int num_alleles;     /* the number of alleles at the locus               */
    int max_alleles;     /* the maximum number of alleles (when run not
                            in parallel, this allows use to malloc only once */
    int num_genotypes;
    int steps;
    int chunks;
    int chunksize;
    int *work;           /* temp array of size MAX_ALLELES                   */
    double *_sumlog2_lookup;
    HWERESULT result;
    FREQS *freqs;
};

typedef struct S_INDEX
{
    int i1;
    int i2;
    int j1;
    int j2;
    int type;
    double cst;
} Index;


static void
test_switch(int *a, Index index, int *switch_ind, int *switch_type,
            double *p1_rt, double *p2_rt)
{
    register int k11, k22, k12, k21;

    *switch_ind = 0;

    k11 = TIDXL(index.i1, index.j1);
    k22 = TIDXL(index.i2, index.j2);
    k12 = TIDXL(index.i1, index.j2);
    k21 = TIDXL(index.i2, index.j1);

    if (index.type <= 1) {      /* type = 0, 1 */
        if (a[k11] > 0 && a[k22] > 0) {
            *switch_ind = 1;
            *switch_type = 0;   /* D-switchable */
            *p1_rt =
                RATIO(a[k11], a[k12]) * RATIO(a[k22], a[k21]) * index.cst;
        }
        if (a[k12] > 0 && a[k21] > 0) {
            *switch_ind += 1;
            *switch_type = 1;   /* R-switchable */
            *p2_rt =
                RATIO(a[k12], a[k11]) * RATIO(a[k21], a[k22]) / index.cst;
        }

    } else {                    /* type = 2 */
        if (a[k11] > 0 && a[k22] > 0) {
            *switch_ind = 1;
            *switch_type = 0;   /* D-switchable */
            *p1_rt =
                RATIO(a[k11], a[k12] + 1.0) * RATIO(a[k22],
                                                    a[k12]) * index.cst;
        }
        if (a[k12] > 1) {
            *switch_ind += 1;
            *switch_type = 1;   /* R-switchable */
            *p2_rt =
                RATIO(a[k12], a[k11]) * RATIO(a[k12] - 1, a[k22]) / index.cst;
        }

    }
}

static inline void
do_switch(int *a, Index index, int type)
{
    register int k11, k22, k12, k21;

    k11 = TIDXL(index.i1, index.j1);
    k12 = TIDXL(index.i1, index.j2);
    k21 = TIDXL(index.i2, index.j1);
    k22 = TIDXL(index.i2, index.j2);


    if (type == 0) {            /* D-switch */
        --a[k11];
        --a[k22];
        ++a[k12];
        ++a[k21];
    } else {                    /* R-switch */
        ++a[k11];
        ++a[k22];
        --a[k12];
        --a[k21];
    }
}

static inline double
cal_prob(int *a, Index index, double ln_p_old, int *actual_switch)
{

    double p1_ratio = 0., p2_ratio = 0.;
    register double ln_p_new;
    double rand_num;
    int switch_ind, type;

    *actual_switch = 0;

/* determine the switchability and direction of switch for given face */

    test_switch(a, index, &switch_ind, &type, &p1_ratio, &p2_ratio);
/*  printf("%d %d %d swhind=%d\n",index.i1, index.i2,index.i3,switch_ind);*/

    switch (switch_ind) {
        case 0:                /* non-switchable */

            ln_p_new = ln_p_old;        /* retain the pattern, probability unchanged */
            break;

        case 1:                /* partially-switchable */

            if (type == 1)
                p1_ratio = p2_ratio;
            rand_num = RANDDBLONE;

            if (rand_num < TRANS(p1_ratio)) {   /* switch w/ transition P TRANS */
                do_switch(a, index, type);
                ln_p_new = ln_p_old + log(p1_ratio);    /* ln P_after-switch */
                *actual_switch = 1;
            } else              /* remain the same w/ P = 1 - TRANS */
                ln_p_new = ln_p_old;    /* probability unchanged */
            break;

        default:               /* fully switchable */
            rand_num = RANDDBLONE;

            if (rand_num <= TRANS(p1_ratio)) {
                do_switch(a, index, 0); /* D-switch */
                ln_p_new = ln_p_old + log(p1_ratio);    /* ln P_after-switch */
                *actual_switch = 2;
            } else if (rand_num <= TRANS(p1_ratio) + TRANS(p2_ratio)) {
                do_switch(a, index, 1); /* R-switch */
                ln_p_new = ln_p_old + log(p2_ratio);
                *actual_switch = 2;
            } else
                ln_p_new = ln_p_old;
            break;
    }

    return (ln_p_new);
}

static inline void
random_choose(HWETEST hwe, int *k1, int *k2)
{
    int temp, i, k = hwe->num_alleles, l;

    for (i = 0; i < k; ++i)
        hwe->work[i] = i;

    *k1 = RANDINT(k);

    --k;

    for (i = *k1; i < k; ++i)
        hwe->work[i] = i + 1;
   
    l = RANDINT(k);
#ifdef DEBUG    
    if(l < 0 || l >= hwe->max_alleles) 
        fprintf(stderr, "hwe.c: pre %i post %i\n", k, l);
#endif    
    assert(l >= 0 && l < hwe->max_alleles);
    *k2 = hwe->work[l];

    if (*k1 > *k2) {
        SWAP(*k1,*k2,temp);
    }
}


static void
select_index(HWETEST hwe, Index * index)
{

    int i1, i2, j1, j2;
    int k = 0;
    int l = 0;

/* generate row indices */

    random_choose(hwe, &i1, &i2);

    index->i1 = i1;
    index->i2 = i2;

/* generate column indices */

    random_choose(hwe, &j1, &j2);

    index->j1 = j1;
    index->j2 = j2;

/* calculate Delta = d(i1,j1) + d(i1,j2) + d(i2,j1) + d(i2,j2) */

    if (i1 == j1)
        ++k;

    if (i1 == j2)
        ++k;

    if (i2 == j1)
        ++k;

    if (i2 == j2)
        ++k;

    index->type = k;

    if ((i1 == j1) || (i2 == j2))
        ++l;

    index->cst = (l == 1) ? pow(2.0, (double)k) : pow(2.0, -(double)k);
}


static int
getMatrixId(HWETEST hwe, int a1, int a2)
{
    assert(a1 >= hwe->freqs->min && a2 >= hwe->freqs->min);
    return TIDXL(hwe->freqs->allele_id_lookup[a1 - hwe->freqs->min],
             hwe->freqs->allele_id_lookup[a2 - hwe->freqs->min]);
}

HWERESULT
HWEgetResult(HWETEST hwe)
{
    return hwe->result;
}

inline void
createSumLog2(HWETEST hwe) {
    int i, n1 = hwe->allele_freq[0];

    /* create the log(2^n) lookup table */
    MAKE1DDOUBLE(hwe->_sumlog2_lookup, n1+1);
    hwe->_sumlog2_lookup[0] = 0;
    hwe->_sumlog2_lookup[1] = log(2);
    for (i=2; i<n1+1; i++)
        hwe->_sumlog2_lookup[i] = hwe->_sumlog2_lookup[i-1] + hwe->_sumlog2_lookup[1];
}

inline double
calcSampleProb(HWETEST hwe, int n11, int n12, int n22)
{
    int n = n11 + n12 + n22, n1=2*n11+n12, n2=2*n22+n12;
    return hwe->_sumlog2_lookup[n12] + (UTILSfactln(n)-(UTILSfactln(n11)+UTILSfactln(n12)+UTILSfactln(n22)))-(UTILSfactln(n+n) - (UTILSfactln(n1)+UTILSfactln(n2)));
}

/* see Principles of population genetics, forth edition, p 59 
 *
 * TODO: I am not sure if this is numerically ok...
 * */
static void
calcTwoAllelesExactTest(HWETEST hwe)
{
    int n11, n11_observed = hwe->matrix[TIDX(0,0)], n22, n1 = hwe->allele_freq[0], n2 = hwe->allele_freq[1], combinations, n12, tmp; 
    double *fwork, p = 0., fcum = 0.;
    
    if (hwe->num_alleles != 2)
        FATALINT("exact test with wrong allele number");

    if (n1 > n2) {
        SWAP(n1,n2,tmp);
        n11_observed = hwe->matrix[TIDX(1,1)];
    }
    
    combinations = 1+n1/2;

    MAKE1DDOUBLE(fwork, combinations);
    createSumLog2(hwe);

    for(n11=0; n11<combinations; n11++) {
        n12 = n1 - 2 * n11;
        n22 = (n2 - n12) / 2;
        fwork[n11] = calcSampleProb(hwe, n11, n12, n22);
    //    fprintf(stderr, "DEBUG: n11: %i, n12 = %i, n22 = %i: %f\n",n11,n12,n22,exp(fwork[n11]));
        if (n11 == n11_observed)
            p = fwork[n11];
    } 
    FREE1D(hwe->_sumlog2_lookup);

    qsort(fwork, (size_t)combinations, sizeof(double), UTILScompare_doubles_incr);
    for(n11=0; n11<combinations; n11++) {
        fcum += exp( fwork[n11] );
      //  fprintf(stderr, "DEBUG: %i = %f   %f %i\n",n11,fwork[n11], p, n11_observed);
        if (fabs(fwork[n11] - p) < EPSILON) {
            hwe->result.p_value = fcum;
            break;
        }    
    }    
    
    FREE1D(fwork);
}

void
HWEtest(HWETEST hwe, int chunks, int chunksize, int steps)
{
    double ln_p_observed, ln_p_simulated, p_mean, p_square;
    double p_simulated;
    int counter, actual_switch;
    Index index;
    register int i, j;
    
    hwe->result.se = 0.;

    if (hwe->num_alleles == 2) {
        calcTwoAllelesExactTest(hwe);
        return;
    }
    else if (hwe->num_alleles <= 1) {
        hwe->result.p_value = 1.0;
        return;
    }

    hwe->steps = steps;
    hwe->chunks = chunks;
    hwe->chunksize = chunksize;
    ln_p_observed = 0.0;
    ln_p_simulated = ln_p_observed;

    p_mean = p_square = 0.;

    hwe->result.p_value = hwe->result.se = 0.;  /* initialization */

    hwe->result.swch_count[0] = hwe->result.swch_count[1] =
        hwe->result.swch_count[2] = 0;

    for (i = 0; i < steps; ++i) {       /* de-memorization for given steps */
        select_index(hwe, &index);

        ln_p_simulated =
            cal_prob(hwe->matrix, index, ln_p_simulated, &actual_switch);

        ++hwe->result.swch_count[actual_switch];
    }

    for (i = 0; i < chunks; ++i) {
        counter = 0;

        for (j = 0; j < chunksize; ++j) {
            select_index(hwe, &index);
            ln_p_simulated =
                cal_prob(hwe->matrix, index, ln_p_simulated, &actual_switch);

            if (ln_p_simulated <= ln_p_observed)
                ++counter;

            ++hwe->result.swch_count[actual_switch];
        }
        p_simulated = counter / (double)chunksize;
        p_mean += p_simulated;
        p_square += p_simulated * p_simulated;

    }

    p_mean /= chunks;
    hwe->result.p_value = p_mean;
    hwe->result.se = p_square / ((double)chunks) / (chunks - 1.0)
        - p_mean / (chunks - 1.0) * p_mean;
    hwe->result.se = sqrt(hwe->result.se);

}

int
HWEgetFreq(HWETEST hwe, int a1, int a2)
{
    return hwe->matrix[getMatrixId(hwe, a1, a2)];
}

int
HWEgetAlleleFreq(HWETEST hwe, int allele)
{
    return hwe->allele_freq[hwe->freqs->
                            allele_id_lookup[allele - hwe->freqs->min]];
}

void 
HWEdumpFreqs(FILE *fp, HWETEST hwe)
{
    register int i, j, k, l;
    char *line;

    MALLOC(line, char, hwe->num_alleles * 6+1);
    line[0] = '-';

    k = 1;

    fprintf(fp, "Observed genotype frequencies: \n\n");

    for (i = 0; i < hwe->num_alleles; ++i) {
        for (j = k; j < k + 5; ++j)
            line[j] = '-';

        line[j] = '\0';
        k = j;
        fprintf(fp, "     %s\n", line);
        l = -1;
        for (j = 0; j < hwe->freqs->max - hwe->freqs->min +1; j++)
            if (hwe->freqs->allele_id_lookup[j] == i)
                l = j+hwe->freqs->min;

        fprintf(fp, "%4i |",l);

        for (j = 0; j <= i; ++j) {
            l = TIDX(i, j);
            fprintf(fp, "%4d|", hwe->matrix_orig[l]);
        }

        fprintf(fp, "\n");
    }

    fprintf(fp, "     %s\n\n", line);
    FREE(line);
}

void
HWEdumpResults(FILE * fp, HWETEST hwe)
{

    int total_step;
    total_step = hwe->steps + hwe->chunks * hwe->chunksize;

    if (hwe->num_alleles == 2) {
        fprintf(fp, "Exact test P-value            : %7.4g\n", hwe->result.p_value);
    }
    else {
        fprintf(fp, "Randomization test P-value    : %7.4g  (%7.4g) \n",
            hwe->result.p_value, hwe->result.se);
        fprintf(fp, "Percentage of partial switches: %6.2f \n",
            hwe->result.swch_count[1] / (double) total_step * 100);
        fprintf(fp, "Percentage of full switches   : %6.2f \n",
            hwe->result.swch_count[2] / (double) total_step * 100);
        fprintf(fp, "Percentage of all switches    : %6.2f \n",
            (hwe->result.swch_count[1] + hwe->result.swch_count[2]) / (double) total_step * 100 );
    }
    fprintf(fp, "\n");
    FFLUSH(fp);
}

static void
calcAlleleFreqs(HWETEST hwe)
{
    int i, j, l;

    for (i = 0; i < hwe->num_alleles; ++i) {
        l = TIDX(i, i);
        hwe->allele_freq[i] = hwe->matrix[l];

        for (j = 0; j < hwe->num_alleles; ++j) {
            l = TIDXL(i, j);
            hwe->allele_freq[i] += hwe->matrix[l];
        }
    }
}

void
HWEcalcFreqs(HWETEST hwe, int ***genotypes, int num_genotypes, int locus,
             int num_alleles, FREQS *freqs)
{
    int i, a1, a2;
    
    assert(locus >= 0);

    /* clear matrix */
    for (i = 0; i < hwe->matrix_length; i++)
        hwe->matrix[i] = 0;

    hwe->freqs = freqs;
    assert(hwe->freqs->min >= 0);
    for (i = 0; i < num_genotypes; i++) {
        a1 = genotypes[i][locus][0];
        a2 = genotypes[i][locus][1];
        if (a1 < 0 || a2 < 0)
            continue;
        hwe->matrix[getMatrixId(hwe, a1, a2)]++;
    }
    for (i = 0; i < hwe->matrix_length; i++)
        hwe->matrix_orig[i] = hwe->matrix[i];

    hwe->num_alleles = num_alleles;
    hwe->num_genotypes = num_genotypes;
    calcAlleleFreqs(hwe);
}

HWETEST
HWEinit(int max_alleles)
{
    HWETEST hwe = malloc(sizeof *hwe);
    if (hwe == NULL)
        FATAL("malloc failed");

    hwe->matrix_length = TLENGTH(max_alleles);
    hwe->max_alleles = max_alleles;
    hwe->_sumlog2_lookup = NULL;
    MAKE1DINT(hwe->matrix, hwe->matrix_length);
    MAKE1DINT(hwe->matrix_orig, hwe->matrix_length);
    MAKE1DINT(hwe->allele_freq, max_alleles);
    MAKE1DINT(hwe->work, max_alleles);
    return hwe;
}

void
HWEdestroy(HWETEST hwe)
{
    FREETRIANGULAR(hwe->matrix);
    FREETRIANGULAR(hwe->matrix_orig);
    FREE1D(hwe->allele_freq);
    FREE1D(hwe->work);
    FREE(hwe);
}

/* 
  vim: ft=c sw=4 ts=4 expandtab
*/
