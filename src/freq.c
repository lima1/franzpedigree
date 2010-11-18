/*
 * $Id: freq.c 2058 2010-05-12 15:19:58Z markus $
 *
 * References:
 * 
 * PIC: Construction  of  a Genetic Linkage Map in Man Using  Restriction 
 * Fragment Length Polymorphisms. Botstein et al. 1980
 *
 * Exclusion Probabilities: Comparisons of three probability formulae for 
 * parentage exclusion. Jamieson & Taylor 1997.
 * Parentage and sibship exclusions: higher statistical power with
 * more family members. J. Wang. Heredity 2007.
 *
 * Identity Probabilities: Estimating the probability of  identity  
 * among  genotypes  in  natural  populations:  cautions and guidelines.
 * Waits et al. 2001.
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
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.*
 *
 * Markus Riester, markus@bioinf.uni-leipzig.de */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "freq.h"

#include "macros.h"
#include "hwe.h"
#include "dataio.h"
#include "null.h"
#include "global.h"
#include "dag.h"
#include "lod.h"
#include "vtprogressbar.h"
#include "utils.h"

extern DATA Data;
extern PROBS Probs;
extern OPTIONS Options;

/* defines the weighting of an observed allele in genotype id */
#define ALLELE_COUNTER(id) ( Options.RametAlleleFreqs ? Data.id_mapping[(id)].ramets : 1 )

void
FREQinit(void)
{
    int i;

    MAKE1DDOUBLE(Data.h_obs, Data.num_loci);
    MAKE1DDOUBLE(Data.h_exp, Data.num_loci);
    MAKE1DDOUBLE(Data.pic, Data.num_loci);
    MAKE1DDOUBLE(Data.pid, Data.num_loci);
    MAKE1DDOUBLE(Data.pid_sib, Data.num_loci);
    MAKE2DDOUBLE(Data.exclp_sp, D_EXCLP_OFFSPRING, Data.num_loci, i);
    MAKE2DDOUBLE(Data.exclp_mk, D_EXCLP_OFFSPRING, Data.num_loci, i);
    MAKE2DDOUBLE(Data.exclp_pp, D_EXCLP_OFFSPRING, Data.num_loci, i);
    /* only trios and quadruples */
    MAKE2DDOUBLE(Data.exclp_fullsibs, 2, Data.num_loci, i);
    MAKE1DDOUBLE(Data.p_null, Data.num_loci);
    MAKE1DDOUBLE(Options.ProportionTyped, Data.num_loci);

    MAKE1DINT(Data.num_alleles, Data.num_loci);
    MALLOC(Data.alleles, int *, Data.num_loci);

    MAKE1DINT(Data.num_typed, Data.num_loci);
    MALLOC(Data.hweresults, HWERESULT, Data.num_loci);
}

void
FREQgetAlleleRange(int locus, RANGE * r)
{
    int i, j, cnt;

    r->min = -1;
    r->max = -1;
    r->num_values = 0;
    for (i = 0; i < Data.num_samples; i++) {
        cnt = ALLELE_COUNTER(i);
        /* simply ignore missing data */
        if (Data.id_mapping[i].genotype_obs[locus][0] < 0 ||
            Data.id_mapping[i].genotype_obs[locus][1] < 0)
            continue;
        for (j = 0; j < 2; j++) {

            r->num_values += cnt;
            
            if (r->min < 0
                || Data.id_mapping[i].genotype_obs[locus][j] < r->min)
                r->min = Data.id_mapping[i].genotype_obs[locus][j];
            if (r->max < 0
                || Data.id_mapping[i].genotype_obs[locus][j] > r->max)
                r->max = Data.id_mapping[i].genotype_obs[locus][j];
        }
    }
}

static void
inRange(int id, int locus, int allele_id) 
{
    int allele = Data.id_mapping[id].genotype_obs[locus][allele_id]; 
    char msg[LINESIZE];
    if (allele >= 0 && ( allele < Data.allele_frequencies[locus].min || 
        allele > Data.allele_frequencies[locus].max ) ) {
        snprintf(msg, sizeof(msg), 
                "Partially missing genotype of %s at locus %i with unkown allele (%i).", Data.id_mapping[id].description,locus,allele);
        FATAL(msg);
    }
}    
static void
calcFrequencies(int locus)
{
    int i, j, cnt;
    RANGE r;

    FREQgetAlleleRange(locus, &r);
    MAKE1DDOUBLE(Data.allele_frequencies[locus].freqs, r.max - r.min + 1);
    MAKE1DINT(Data.allele_frequencies[locus].allele_id_lookup,
                  r.max - r.min + 1);
    Data.allele_frequencies[locus].min = r.min;
    Data.allele_frequencies[locus].max = r.max;

    for (i = 0; i < Data.num_samples; i++) {

        cnt = ALLELE_COUNTER(i);
        assert(cnt >= 1);

        if (Data.id_mapping[i].genotype_obs[locus][0] < 0 ||
            Data.id_mapping[i].genotype_obs[locus][1] < 0) {
            inRange(i,locus,0); 
            inRange(i,locus,1); 
            continue;
        }    

        for (j = 0; j < 2; j++)
            Data.allele_frequencies[locus].freqs[Data.id_mapping[i].
                                                    genotype_obs[locus]
                                                    [j] - r.min] +=
                cnt / (double)r.num_values;
    }      
}


/* update the frequencies given a pedigree */
void
FREQupdate(Dag D)
{
    int i, j, k, cnt, indegree, num_alleles;
    Edge incoming[2];

    /* reset */
    for (i = 0; i < Data.num_loci; i++) {
        Data.allele_frequencies[i].num_ind = 0;
        num_alleles = Data.allele_frequencies[i].max - Data.allele_frequencies[i].min +
             1;
        for (j = 0; j < num_alleles; j++) {
            /* make sure that every observed allele is considered */
            if (Data.allele_frequencies[i].freqs[j] > 0.) {
                Data.allele_frequencies[i].freqs[j] = 1.;
                Data.allele_frequencies[i].num_ind++;
            } else
                Data.allele_frequencies[i].freqs[j] = 0.;
        }
    }

    /* count alleles */
    for (i = 0; i < Data.num_samples; i++) {
        indegree = DAGincomingE(D, i, incoming);

        /* ignore offspring alleles */
        if (indegree >= 2)
            continue;

        cnt = ALLELE_COUNTER(i);
        for (j = 0; j < Data.num_loci; j++) {
            if (indegree == 0) {
                for (k = 0; k < 2; k++)
                    if (Data.id_mapping[i].genotype_obs[j][k] >= 0) {
                        Data.allele_frequencies[j].freqs[Data.id_mapping[i].
                                                         genotype_obs[j][k] -
                                                         Data.
                                                         allele_frequencies
                                                         [j].min] += cnt;
                        Data.allele_frequencies[j].num_ind += cnt;
                    }
            } else if (indegree == 1) {
                for (k = 0; k < 2; k++) {
                    if (Data.id_mapping[i].genotype_obs[j][k] >= 0
                        && Data.id_mapping[i].genotype_obs[j][k] !=
                        Data.id_mapping[incoming[0].v].genotype_obs[j][0]
                        && Data.id_mapping[i].genotype_obs[j][k] !=
                        Data.id_mapping[incoming[0].v].genotype_obs[j][1]) {
                        Data.allele_frequencies[j].freqs[Data.id_mapping[i].
                                                         genotype_obs[j][k] -
                                                         Data.
                                                         allele_frequencies
                                                         [j].min] += cnt;
                        Data.allele_frequencies[j].num_ind += cnt;
                        /*  fprintf(stderr, "DEBUG: %i %i %i %i %i %f\n",
                           Data.id_mapping[i].ramets, i,j,k, incoming[0].v,
                           Data.allele_frequencies[j].num_ind); */
                    }
                }
            }
        }
    }
    /* We counted at the beginning every allele once. Undo that, excluding the
     * alleles with still only one count. We don't wanna loose rare alleles. */
    for (i = 0; i < Data.num_loci; i++) {
        num_alleles = Data.allele_frequencies[i].max - Data.allele_frequencies[i].min +
             1;
        for (j = 0; j < num_alleles; j++) 
            if (Data.allele_frequencies[i].freqs[j] > 1.0) {
                Data.allele_frequencies[i].num_ind--;
                Data.allele_frequencies[i].freqs[j] -= 1.0;
            }
    }

    /* make probabilities */
    for (i = 0; i < Data.num_loci; i++)
        for (j = 0;
             j <
             Data.allele_frequencies[i].max - Data.allele_frequencies[i].min +
             1; j++)
            Data.allele_frequencies[i].freqs[j] /= Data.allele_frequencies[i].
                num_ind;
}

void
FREQcalcAlleleFrequencies(void)
{
    int i;
    
    MALLOC(Data.allele_frequencies, FREQS, Data.num_loci);

#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(guided)
#endif    
    for (i = 0; i < Data.num_loci; i++)
        calcFrequencies(i);
}

static void
calcHeterozygosity(int locus)
{
    int i, j, n_alleles, max_n_alleles = 0;
    double b, c, d, sumf1, sumf2, sumf3, sumf4, sumf5, sumf6, sumf7, sumf8,
        sumh, picf, pow2, n;
    FREQS f;

    sumh = 0;
    for (i = 0; i < Data.num_samples; i++) {
        if (Data.id_mapping[i].genotype_obs[locus][0] < 0 ||
            Data.id_mapping[i].genotype_obs[locus][1] < 0)
            continue;
        if (Data.id_mapping[i].genotype_obs[locus][0] !=
            Data.id_mapping[i].genotype_obs[locus][1])
            sumh++;
        Data.num_typed[locus]++;
    }
    /* if user defined a proportion, use this value. otherwise use the
     * data to estimate this values */
    if (Options.ProportionTypedDefined < 0)
        Options.ProportionTyped[locus] =
            Data.num_typed[locus] / (double)Data.num_samples;
    else
        Options.ProportionTyped[locus] = Options.ProportionTypedDefined;

    Data.h_obs[locus] = sumh / Data.num_typed[locus];
    f = Data.allele_frequencies[locus];

    sumf1 = sumf2 = sumf3 = sumf4 = sumf5 = sumf6 = sumf7 = sumf8 = 0.;

    n_alleles = 0;
    MAKE1DINT(Data.alleles[locus], (f.max - f.min + 1));
    for (i = 0; i < f.max - f.min + 1; i++) {
        if (f.freqs[i] > 0) {
            f.allele_id_lookup[i] = n_alleles;
            Data.alleles[locus][n_alleles] = i + f.min;
            n_alleles++;
        } else {
            f.allele_id_lookup[i] = -1;
            continue;
        }
        sumf1 += f.freqs[i];
        pow2 = f.freqs[i] * f.freqs[i];
        sumf2 += pow2;
        sumf3 += f.freqs[i] * pow2;
        sumf4 += pow2 * pow2;
        sumf5 += f.freqs[i] * pow2 * pow2;
        sumf6 += pow2 * pow2 * pow2;
        sumf7 += f.freqs[i] * pow2 * pow2 * pow2;
        sumf8 += pow2 * pow2 * pow2 * pow2;
    }
    Data.num_alleles[locus] = n_alleles;

    REALLOC(Data.alleles[locus], int, n_alleles);

    if (n_alleles > max_n_alleles)
        max_n_alleles = n_alleles;

    Data.h_exp[locus] = 1. - sumf2;
    Data.exclp_sp[0][locus] =
        1 - 4 * sumf2 + 2 * sumf2 * sumf2 + 4 * sumf3 - 3 * sumf4;
    Data.exclp_mk[0][locus] =
        1 - 2 * sumf2 + sumf3 + 2 * sumf4 - 3 * sumf5 - 2 * sumf2 * sumf2 +
        3 * sumf2 * sumf3;
    Data.exclp_pp[0][locus] =
        1 + 4 * sumf4 - 4 * sumf5 - 3 * sumf6 - 8 * sumf2 * sumf2 +
        8 * sumf2 * sumf3 + 2 * sumf3 * sumf3;

    b = 0.25;
    c = 0.5;
    d = 0.75;
    for (i = 1; i < D_EXCLP_OFFSPRING; i++) {
        b *= b;
        c *= c;
        d *= d;
        Data.exclp_mk[i][locus] =
            1 - 4 * c * sumf2 - 2 * (1 + 4 * b - 2 * c) * sumf2 * sumf2 +
            8 * (b + c - d) * sumf2 * sumf2 * sumf2 - 2 * (1 -
                                                           3 * c) *
            sumf3 - 2 * (6 * b + 3 * c - 4 * d) * sumf3 * sumf3 + (3 +
                                                                   8 * b -
                                                                   6 *
                                                                   c) *
            sumf4 - 4 * (9 * b + 6 * c - 7 * d) * sumf2 * sumf4 +
            4 * (7 * b + c - 2 * d) * (sumf2 * sumf3 - sumf5) +
            2 * (20 * b + 11 * c - 14 * d) * sumf6;

        Data.exclp_sp[i][locus] =
            1 - 8 * c * sumf2 - 4 * (1 - 6 * b) * sumf2 * sumf2 -
            8 * (3 * b - 3 * c + d) * sumf2 * sumf2 * sumf2 - 4 * (1 -
                                                                   6 * b -
                                                                   c) *
            sumf3 - 8 * (21 * b - 12 * c + d) * sumf2 * sumf3 + 2 * (3 +
                                                                     48 *
                                                                     b -
                                                                     36 *
                                                                     c +
                                                                     4 *
                                                                     d) *
            sumf3 * sumf3 + 6 * (1 - 14 * b + 4 * c) * sumf4 + 4 * (2 +
                                                                    45 *
                                                                    b -
                                                                    37 *
                                                                    c +
                                                                    7 *
                                                                    d) *
            sumf2 * sumf4 + 2 * (1 + 108 * b - 62 * c + 4 * d) * sumf5 -
            (15 + 264 * b - 204 * c + 28 * d) * sumf6;
        Data.exclp_pp[i][locus] =
            1 - 32 * b * sumf2 * sumf2 - 8 * (1 -
                                              2 * c) * sumf2 * sumf2 *
            sumf2 * sumf2 - 16 * (1 - 2 * b - c) * sumf2 * sumf2 * sumf3 -
            8 * (1 + 3 * b - 4 * c) * sumf3 * sumf3 - 64 * (3 * b -
                                                            3 * c +
                                                            d) * sumf2 *
            sumf3 * sumf3 + 16 * b * sumf4 + 8 * (1 - 24 * b +
                                                  10 * c) * sumf2 *
            sumf2 * sumf4 + 8 * (3 - 26 * b + 13 * c -
                                 4 * d) * sumf3 * sumf4 + 2 * (5 +
                                                               120 * b -
                                                               94 * c +
                                                               16 * d) *
            sumf4 * sumf4 - 16 * (3 * b - c) * (sumf5 - 2 * sumf2 * sumf3)
            - 16 * (2 * b - c) * sumf2 * (3 * sumf4 + 8 * sumf5 -
                                          2 * sumf2 * sumf2)
            + 16 * (1 + 31 * b - 31 * c + 9 * d) * sumf3 * sumf5 +
            4 * (1 + b - 4 * c) * sumf6 + 32 * (1 + 20 * b - 15 * c +
                                                2 * d) * sumf2 * sumf6 -
            4 * (1 - 114 * b + 67 * c - 8 * d) * sumf7 - (996 * b -
                                                          880 * c +
                                                          176 * d +
                                                          59) * sumf8;
    }
    if (Data.num_alleles[locus] < 3) {
        Data.exclp_fullsibs[0][locus] = 0.;
        Data.exclp_fullsibs[1][locus] = 0.;
    } else {
        Data.exclp_fullsibs[0][locus] =
            1 - 30 * sumf2 * sumf2 + 16 * sumf2 * sumf2 * sumf2 -
            22 * sumf3 * sumf3 + 72 * sumf2 * sumf3 - 60 * sumf2 * sumf4 +
            15 * sumf4 - 48 * sumf5 + 56 * sumf6;
        Data.exclp_fullsibs[1][locus] =
            1 - 60 * sumf2 * sumf2 * sumf2 * sumf2 -
            192 * sumf2 * sumf2 * sumf3 - 112 * sumf3 * sumf3 +
            288 * sumf2 * sumf3 * sumf3 + 396 * sumf2 * sumf2 * sumf4 +
            480 * sumf3 * sumf4 - 315 * sumf4 * sumf4 + 240 * sumf2 * sumf5 -
            552 * sumf3 * sumf5 + 84 * sumf6 - 696 * sumf2 * sumf6 -
            480 * sumf7 + 918 * sumf8;
    }
    n = (double)Data.num_typed[locus];
    Data.pid[locus] =
        (n * n * n * (2 * sumf2 * sumf2 - sumf4) -
         2 * n * n * (sumf3 + 2 * sumf3) + n * (9 * sumf2 + 2) -
         6) / (double)((n - 1) * (n - 2) * (n - 3));
    Data.pid_sib[locus] =
        0.25 + 0.5 * sumf2 + 0.5 * sumf2 * sumf2 - 0.25 * sumf4;

    picf = 0.;
    for (i = 0; i < f.max - f.min; i++)
        for (j = i + 1; j < f.max - f.min + 1; j++)
            picf += f.freqs[i] * f.freqs[i] * f.freqs[j] * f.freqs[j];
    Data.pic[locus] = Data.h_exp[locus] - 2 * picf;

    Data.h_exp[locus] *=
        2. * Data.num_typed[locus] / (double)(2. * Data.num_typed[locus] -
                                                1.);

}

static void
calcHWE(int locus, HWETEST * hwe, int ***genotypes)
{
    HWERESULT hweresult;

    /* Hardy-Weinberg */
    hwe[locus] = HWEinit(Data.num_alleles[locus]);
    HWEcalcFreqs(hwe[locus], genotypes, Data.num_samples, locus,
                 Data.num_alleles[locus], &Data.allele_frequencies[locus]);
    HWEtest(hwe[locus], Options.HWESteps, Options.HWEChunks,
            Options.HWEChunkSize);
    hweresult = HWEgetResult(hwe[locus]);

    Data.hweresults[locus].p_value = hweresult.p_value;
    Data.hweresults[locus].se = hweresult.se;
}

void
FREQcalcSummaryStatistics(void)
{
    int i, ***genotypes, loci_completed = 0;
    bool dumphweresults = false;
    double f,s = 0;
    HWETEST *hwe;
    NULLTEST *nt;
    FILE *hweout = NULL, *locisum = NULL;

    FREQinit();

    /* run the hwetest and null allele test on multiple cpus */
    MALLOC(hwe, HWETEST, Data.num_loci);
    MALLOC(nt, NULLTEST, Data.num_loci);

    /* for the wrapper around hwe.c */
    MALLOC(genotypes, int **, Data.num_samples);

    for (i = 0; i < Data.num_samples; i++)
        genotypes[i] = Data.id_mapping[i].genotype_obs;

    if (Options.Verbosity > 0) 
        VTPROGRESSBARupdate("Allele Frequency Analysis", Data.num_loci, 0);
    
#ifdef HAVE_OPENMP
    #pragma omp parallel for schedule(guided)
#endif
    for (i = 0; i < Data.num_loci; i++) {
        
#ifdef HAVE_OPENMP
        if (omp_get_thread_num() == 0 && Options.Verbosity >= 1)
#else            
        if (Options.Verbosity >= 1) 
#endif 
                VTPROGRESSBARupdate("Allele Frequency Analysis",
                                    Data.num_loci, loci_completed);
        calcHeterozygosity(i);
        calcHWE(i, hwe, genotypes);

        /* null alleles */
        nt[i] = NULLinit(i);
        Data.p_null[i] = NULLcalcPn(nt[i]);

#ifdef HAVE_OPENMP
#pragma omp atomic
#endif        
        loci_completed += 1;
    }
    
    /* which is the most polymorphic locus? */
    Data.max_num_alleles = -1;
    for (i = 0; i < Data.num_loci; i++) {
        f = 1 - Data.h_obs[i] / Data.h_exp[i];
        s += 1 - ( (1. - f) / (1. + f));
        if (Data.num_alleles[i] > Data.max_num_alleles)
            Data.max_num_alleles = Data.num_alleles[i];
    }
    Data.SelfingRateAvg = s / Data.num_loci;

    /* output detailed hwe test results */
    if (strlen(Options.HWETestOutfilename) > 0) {
        dumphweresults = true;
        FOPENW(hweout, Options.HWETestOutfilename);
        fprintf(hweout, "HWE Steps                   : %i\n",
                Options.HWESteps);
        fprintf(hweout, "HWE Chunks                  : %i\n",
                Options.HWEChunks);
        fprintf(hweout, "HWE Chunk Size              : %i\n",
                Options.HWEChunkSize);
        fprintf(hweout, "\n\n");
    }

    if (dumphweresults) {
        for (i = 0; i < Data.num_loci; i++) {
            fprintf(hweout, "*** Locus %s ***\n\n", Data.loci_ids[i]);
            HWEdumpFreqs(hweout, hwe[i]);

            HWEdumpResults(hweout, hwe[i]);
        }
        FCLOSE(hweout);
    }

    FREE(genotypes);

    /* output detailed summary statistics */
    FOPENW(locisum, Options.LociOutfilename);
    DATAIOdumpDatasetDetails(locisum);
    for (i = 0; i < Data.num_loci; i++) {
        HWEdestroy(hwe[i]);
        NULLdump(locisum, nt[i]);
        NULLdestroy(nt[i]);
    }
    FCLOSE(locisum);
    FREE1D(hwe);
    FREE(nt);

    if (Options.Verbosity >= 1)
        VTPROGRESSBARcompleted("Allele Frequency Analysis");
}

void
FREQdumpSummaryStatistics(FILE * fp)
{
    int i, j;
    double avgexclp_sp[D_EXCLP_OFFSPRING], avgexclp_mk[D_EXCLP_OFFSPRING],
        avgexclp_pp[D_EXCLP_OFFSPRING], avgexclp_fullsibs[2], avg, var,
        avgid = 1., avgid_sib = 1.;

    /* output the statistics in a compact format for summary.dat   */
    fprintf(fp, "*** Summary Statistics ***\n\n");
    fprintf(fp,
            "%-10s  %7s  %5s  %5s  %5s  %5s  %5s  %5s  %5s  %5s  %5s  %5s  %5s  %6s  %6s  %6s\n",
            "Locus", "Alleles", "Min", "Max", "N", "Hobs", "Hexp", "PIC",
            "EX 1P", "EX 2P", "EX PP", "ID", "IDsib", "P_NULL", "HWE PV",
            "HWE SE");

    for (j = 0; j < D_EXCLP_OFFSPRING; j++)
        avgexclp_sp[j] = avgexclp_mk[j] = avgexclp_pp[j] = 1.;

    avgexclp_fullsibs[0] = avgexclp_fullsibs[1] = 1.;

    for (i = 0; i < Data.num_loci; i++) {
        fprintf(fp,
                "%-10s  %7i  %5i  %5i   %5i  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.4f  %.4f  %.4f\n",
                Data.loci_ids[i], Data.num_alleles[i], Data.allele_frequencies[i].min,
                Data.allele_frequencies[i].max, Data.num_typed[i],
                Data.h_obs[i], Data.h_exp[i], Data.pic[i],
                Data.exclp_sp[0][i], Data.exclp_mk[0][i], Data.exclp_pp[0][i],
                Data.pid[i], Data.pid_sib[i], Data.p_null[i],
                Data.hweresults[i].p_value, Data.hweresults[i].se);
        for (j = 0; j < D_EXCLP_OFFSPRING; j++) {
            avgexclp_sp[j] *= (1 - Data.exclp_sp[j][i]);
            avgexclp_mk[j] *= (1 - Data.exclp_mk[j][i]);
            avgexclp_pp[j] *= (1 - Data.exclp_pp[j][i]);
        }
        for (j = 0; j < 2; j++)
            avgexclp_fullsibs[j] *= (1 - Data.exclp_fullsibs[j][i]);

        avgid *= Data.pid[i];
        avgid_sib *= Data.pid_sib[i];
        
    }

    for (i = 0; i < D_EXCLP_OFFSPRING; i++) {
        avgexclp_sp[i] = 1 - avgexclp_sp[i];
        avgexclp_mk[i] = 1 - avgexclp_mk[i];
        avgexclp_pp[i] = 1 - avgexclp_pp[i];
    }
    for (i = 0; i < 2; i++)
        avgexclp_fullsibs[i] = 1 - avgexclp_fullsibs[i];

    fprintf(fp, "\n");
    fprintf(fp, "Number of sampling locations: %i\n", Data.num_populations);
    if (Data.num_samples != Data.num_ramets) {
        fprintf(fp, "Number of genets            : %i\n", Data.num_samples);
        fprintf(fp, "Number of ramets            : %i\n", Data.num_ramets);
    }
    else
        fprintf(fp, "Number of individuals       : %i\n", Data.num_samples);
    fprintf(fp, "Number of loci              : %i\n", Data.num_loci);
    fprintf(fp, "\n");

    UTILScalcVarianceArrayInt(Data.num_alleles,0,Data.num_loci-1,&var, &avg);
    fprintf(fp, "Average number of alleles   : %7.3f (+- %.3f)\n", avg,sqrt(var));
    UTILScalcVarianceArray(Data.h_obs,0,Data.num_loci-1,&var, &avg);
    fprintf(fp, "Average observed heterozyg. : %7.3f (+- %.3f)\n", avg, sqrt(var));
    UTILScalcVarianceArray(Data.h_exp,0,Data.num_loci-1,&var, &avg);
    fprintf(fp, "Average expected heterozyg. : %7.3f (+- %.3f)\n", avg, sqrt(var));
    UTILScalcVarianceArray(Data.pic,0,Data.num_loci-1,&var, &avg);
    fprintf(fp, "Average PIC                 : %7.3f (+- %.3f)\n", avg, sqrt(var));
    fprintf(fp, "\n");
    fprintf(fp,
            "Cumulative exclusion probability when 1 to %i fullsibs are genotyped\n",
            D_EXCLP_OFFSPRING);
    fprintf(fp, "  First Parent              :");
    for (i = 0; i < D_EXCLP_OFFSPRING; i++)
        fprintf(fp, " %.7f ", avgexclp_sp[i]);
    fprintf(fp, "\n  Second Parent             :");
    for (i = 0; i < D_EXCLP_OFFSPRING; i++)
        fprintf(fp, " %.7f ", avgexclp_mk[i]);
    fprintf(fp, "\n  Parent Pair               :");
    for (i = 0; i < D_EXCLP_OFFSPRING; i++)
        fprintf(fp, " %.7f ", avgexclp_pp[i]);
    fprintf(fp, "\n");
    fprintf(fp, "\nCumulative probability of identity\n");
    fprintf(fp, "  2 unrelated individuals   : %.7f\n", avgid);
    fprintf(fp, "  Siblings                  : %.7f\n", avgid_sib);
    fprintf(fp, "\n");
    fprintf(fp, "Cumulative sibship exclusion probability\n");
    fprintf(fp, "  3 unrelated individuals   : %.7f\n", avgexclp_fullsibs[0]);
    fprintf(fp, "  4 unrelated individuals   : %.7f\n", avgexclp_fullsibs[1]);
    if (Options.Selfing) {
        fprintf(fp, "\n");
        fprintf(fp, "Estimated Selfing Rate (s)  : %.3f\n", Data.SelfingRateAvg);
    }
    fprintf(fp, "\n\n");
    FFLUSH(fp);
}

void
FREQdump(FILE * fp)
{
    int i, j, n, nnn;

    fprintf(fp, "%i\n", Data.num_loci);

    for (i = 0; i < Data.num_loci; i++) {
        n = (Data.allele_frequencies[i].max - Data.allele_frequencies[i].min +
             1);
        nnn = 0;
        for (j = 0; j < n; j++)
            if (Data.allele_frequencies[i].freqs[j] > 0.)
                nnn++;
        fprintf(fp, "%i %i %i\n", nnn, Data.allele_frequencies[i].min,
                Data.allele_frequencies[i].max);
        for (j = 0; j < n; j++)
            if (Data.allele_frequencies[i].freqs[j] > 0.)
                fprintf(fp, "%i %f\n", j + Data.allele_frequencies[i].min,
                        Data.allele_frequencies[i].freqs[j]);
    }
}

void 
FREQdestroyPreMCMC(void)
{
    int i;
    FREE(Data.hweresults);
    FREE(Data.h_obs);
    FREE(Data.h_exp);
    FREE2D(Data.exclp_mk, D_EXCLP_OFFSPRING, i);
    FREE2D(Data.exclp_pp, D_EXCLP_OFFSPRING, i);
    FREE2D(Data.exclp_fullsibs, 2, i);
    FREE(Data.p_null);
    FREE(Data.pic);
    FREE(Data.pid);
    FREE(Data.pid_sib);
    FREE(Options.ProportionTyped);
}

void
FREQdestroyPostMCMC(void)
{
    int i;

    for (i = 0; i < Data.num_loci; i++) {
        FREE1D(Data.allele_frequencies[i].freqs);
        FREE1D(Data.allele_frequencies[i].allele_id_lookup);
    }
    FREE2D(Data.exclp_sp, D_EXCLP_OFFSPRING, i);
    FREE(Data.allele_frequencies);
    FREE2D(Data.alleles, Data.num_loci, i);
    FREE(Data.num_alleles);
    FREE(Data.num_typed);
}

/* 
  vim: ft=c sw=4 ts=4 expandtab
*/
