/*
 * $Id: null.c 1885 2010-01-25 15:31:17Z markus $
 *
 * Estimates Null Allele Frequency with the simple Chakraborty formula
 * or the EM algorithm as described in:
 *
 * "Maximum likelihood estimation of the frequency of null alleles
 *  at microsatellite loci." Kalinowski & Taper 2006.
 *
 *  With known sub pedigrees, the homozygote/homozygote mismatches are
 *  incorporated (our method)
 *
 * Copyright (C) 2008 Universitaet Leipzig  
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "null.h"

#include "global.h"
#include "dataio.h"
#include "options.h"
#include "freq.h"

struct S_NULLTEST
{
    int locus;                  /* the locus to test                          */
    double p_n;                 /* the null allele frequency                  */
    double beta;                /* the genotyping failure rate, see K&T paper */
    double *p_i_it;             /* the altered allele frequencies             */
    int *ni;                    /* the number of i alleles                    */
    int *nii;                   /* number of i homozygotes                    */
    int *sum_nij;               /* number of i heterozygotes                  */
    int *niihm;                 /* number of i homozygotes with mismatches    */
    int **nij;                  /* number of i.j heterozygotes                */
    int nmm;                    /* number untyped genotypes                   */
};

extern OPTIONS Options;
extern DATA Data;

NULLTEST
NULLinit(int locus)
{
    int i;
    NULLTEST nt = malloc(sizeof *nt);
    if (nt == NULL)
        FATAL("malloc failed");

    nt->locus = locus;
    MAKE2DINT(nt->nij, Data.num_alleles[locus], Data.num_alleles[locus], i);
    MAKE1DDOUBLE(nt->p_i_it, Data.num_alleles[locus]);
    MAKE1DINT(nt->ni, Data.num_alleles[locus]);
    MAKE1DINT(nt->nii, Data.num_alleles[locus]);
    MAKE1DINT(nt->niihm, Data.num_alleles[locus]);
    MAKE1DINT(nt->sum_nij, Data.num_alleles[nt->locus]);
    return nt;
}

double
NULLcalcPnChakraborty(NULLTEST nt)
{
    nt->p_n =
        (Data.h_exp[nt->locus] -
         Data.h_obs[nt->locus]) / (Data.h_exp[nt->locus] +
                                   Data.h_obs[nt->locus]);
    return nt->p_n;
}

static void
getHomozygoteMismatches(int locus, int *has_hm)
{
    int i, j, k;

    /* has_hm is 1 if genotype is a null allele for sure */
    for (i = 0; i < Data.num_samples; i++)
        has_hm[i] = 0;
    /* we need a sub-pedigree for this */
    if (!Data.use_pedigree)
        return;

    /* now search for homozygote/homozygote mismatches */
    for (i = 0; i < Data.num_samples; i++)
        for (j = 0; j < 2; j++) {
            k = Data.in_pedigree[i][j];
            if (k < 0)
                continue;
            if (Data.in_pedigree_mismatches[i][j][locus] > 0 &&
                Data.id_mapping[i].genotype_obs[locus][0] ==
                Data.id_mapping[i].genotype_obs[locus][1]
                && Data.id_mapping[k].genotype_obs[locus][0] ==
                Data.id_mapping[k].genotype_obs[locus][1]) {

                has_hm[i] = 1;
                has_hm[k] = 1;
                /* otherwise in_pedigree_mismatches is wrong */
                assert(Data.id_mapping[i].genotype_obs[locus][0] !=
                       Data.id_mapping[k].genotype_obs[locus][0]);
            }
        }
}

/* TODO I should move this to freq.c */
void
NULLdump(FILE * fp, NULLTEST nt)
{
    int i, allele_id, sum_ni=0,sum_nij=0,sum_nii=0;
    double q;

    fprintf(fp, "*** Locus %s ***\n\n", Data.loci_ids[nt->locus]);
    fprintf(fp, " +--------+--------+------------+------------+-----------------------+--------------------------+\n");
    fprintf(fp, " | %6s | %6s | %10s | %10s | %10s %10s | %24s |\n", "Allele", "Count",
            "Heterozyg.", "Homozyg.", "Frequency", "SE",
            "Frequency (Null Alleles)");
    fprintf(fp, " +--------+--------+------------+------------+-----------------------+--------------------------+\n");

    /* QnD: calc. sum ni*/
    for (i = 0;i < FREQ_RANGE(nt->locus); i++) {
        if (Data.allele_frequencies[nt->locus].freqs[i] <= 0.)
            continue;
        allele_id =
            FREQ_ALLELE_ID(nt->locus,
                            i + Data.allele_frequencies[nt->locus].min);

        sum_ni  += nt->ni[allele_id];
    }
    for (i = 0;i < FREQ_RANGE(nt->locus); i++) {
        if (Data.allele_frequencies[nt->locus].freqs[i] <= 0.)
            continue;
        allele_id =
            FREQ_ALLELE_ID(nt->locus,
                            i + Data.allele_frequencies[nt->locus].min);
        q = Data.allele_frequencies[nt->locus].freqs[i];

        fprintf(fp, " | %6i | %6i | %10i | %10i | %10.5f %10.5f | %24.5f |\n",
                Data.allele_frequencies[nt->locus].min + i, nt->ni[allele_id],
                nt->sum_nij[allele_id], nt->nii[allele_id],q,
                sqrt((q*(1-q))/(double)sum_ni),
                nt->p_i_it[allele_id]);

        sum_nij += nt->sum_nij[allele_id];
        sum_nii += nt->nii[allele_id];
    }
    fprintf(fp, " +--------+--------+------------+------------+-----------------------+--------------------------+\n");
    fprintf(fp, "   %6s   %6i   %10i   %10i\n","TOTAL:", sum_ni, sum_nij, sum_nii);
    fprintf(fp, "\n");
    fprintf(fp, "Observed Heterozygosity     : %.4f\n",
            Data.h_obs[nt->locus]);
    fprintf(fp, "Expected Heterozygosity     : %.4f\n",
            Data.h_exp[nt->locus]);
    fprintf(fp, "Hardy-Weinberg Test\n");
    fprintf(fp, "   p-Value                  : %.4f\n",
            Data.hweresults[nt->locus].p_value);
    fprintf(fp, "   Standard Error           : %.4f\n",
            Data.hweresults[nt->locus].se);
    fprintf(fp, "Estimated Null Allele Freq. : %.4f\n", nt->p_n);
    fprintf(fp, "Estimated Genotyping Failure: %.4f (%i untyped)\n", nt->beta,
            nt->nmm);
    fprintf(fp, "Polymorphic Inform. Content : %.4f\n", Data.pic[nt->locus]);
    fprintf(fp, "\n");
    fprintf(fp, "Exclusion probability when 1 to %i fullsibs are genotyped\n", D_EXCLP_OFFSPRING);
    fprintf(fp, "  First Parent              :");
    for (i=0; i<D_EXCLP_OFFSPRING; i++)
        fprintf(fp, " %.7f ", Data.exclp_sp[i][nt->locus]);
    fprintf(fp, "\n  Second Parent             :");
    for (i=0; i<D_EXCLP_OFFSPRING; i++)
        fprintf(fp, " %.7f ", Data.exclp_mk[i][nt->locus]);

    fprintf(fp, "\n  Parent Pair               :");
    for (i=0; i<D_EXCLP_OFFSPRING; i++)
        fprintf(fp, " %.7f ", Data.exclp_pp[i][nt->locus]);
    fprintf(fp, "\n");
    fprintf(fp, "\nProbability of identity\n");
    fprintf(fp, "  2 unrelated individuals   : %.7f\n", Data.pid[nt->locus]);
    fprintf(fp, "  Siblings                  : %.7f\n",
            Data.pid_sib[nt->locus]);
    fprintf(fp, "\n");
    fprintf(fp, "Sibship exclusion probability\n");
    fprintf(fp, "  3 unrelated individuals   : %.7f\n", Data.exclp_fullsibs[0][nt->locus]);
    fprintf(fp, "  4 unrelated individuals   : %.7f\n", Data.exclp_fullsibs[1][nt->locus]);
    fprintf(fp, "\n\n");
}

double
NULLcalcPnKalinowskiTaper(NULLTEST nt)
{
    int i, j, *has_hm, a1, a1_id, a2, a2_id, atmp, N, conv_events = 0;

    /* initialize the iteraton variables */
    double p_n_it =
        NULLcalcPnChakraborty(nt), p_n_it_old, p_n_it_sqr, beta, sumt;

    MAKE1DINT(has_hm, Data.num_samples);
    for (i = 0; i < Data.num_alleles[nt->locus]; i++) {
        nt->p_i_it[i] = Data.allele_frequencies[nt->locus].freqs[i];
        nt->nii[i] = 0;
        nt->niihm[i] = 0;
    }

    N = Data.num_samples;

    nt->nmm = 0;
    /* now calculate the genotype numbers */
    for (i = 0; i < Data.num_samples; i++) {
        /* the two alleles */
        a1 = Data.id_mapping[i].genotype_obs[nt->locus][0];
        a2 = Data.id_mapping[i].genotype_obs[nt->locus][1];
        if (a2 < a1)
            SWAP(a1, a2, atmp);

        a1_id = FREQ_ALLELE_ID(nt->locus, a1);
        a2_id = FREQ_ALLELE_ID(nt->locus, a2);

        if (a1 < 0 && a2 < 0) {
            nt->nmm++;
            continue;
        }
        /* ignore all with one missing allele, a2 is >= 0 because of swap */
        if (a1 < 0) {
            assert(a2_id >= 0);
            N--;
            /* should we do this? */
            nt->ni[a2_id]++;
            continue;
        }
        nt->ni[a1_id]++;
        nt->ni[a2_id]++;
        if (a1 == a2)
            nt->nii[a1_id]++;
        else
            nt->nij[a1_id][a2_id]++;
    }

    /*printf("beta: %i, %f %f, \n", nt->nmm, p_n_it, N *p_n_it * p_n_it);*/

    for (i = 0; i < Data.num_alleles[nt->locus]; i++) {
        nt->sum_nij[i] = 0;
        for (j = 0; j < Data.num_alleles[nt->locus]; j++) {
            if (i == j)
                continue;
            nt->sum_nij[i] += nt->nij[MIN(i, j)][MAX(i, j)];
        }
    }
    getHomozygoteMismatches(nt->locus, has_hm);
    for (i = 0; i < Data.num_samples; i++) {
        if (!has_hm[i])
            continue;
        a1 = Data.id_mapping[i].genotype_obs[nt->locus][0];
        nt->niihm[FREQ_ALLELE_ID(nt->locus, a1)]++;
        nt->nii[FREQ_ALLELE_ID(nt->locus, a1)]--;
    }

    /* I don't think this is necessary, but negative p_n don't make much sense
     * either */
    if (p_n_it < 0.)
        p_n_it = EPSILON;

    /* use the initial guess of p_n to calculate the expected number of
     * nmm, then calculate the ratio to get the first guess of beta     */
    beta = MAX(0, (nt->nmm - N *p_n_it * p_n_it) / (double)nt->nmm);

    /* now start the iterations */
    for (i = 0; i < D_NULLSTEPS; i++) {
        /*fprintf(stderr, "DEBUG: locus %i iteration %i p_n %f, beta %f\n", nt->locus, i, p_n_it, beta); */
        /* update p_i */
        p_n_it_sqr = p_n_it * p_n_it;
        p_n_it_old = p_n_it;

        for (j = 0; j < Data.num_alleles[nt->locus]; j++) {

            nt->p_i_it[j] = 1 / (double)(2 * N) * ((2 * nt->nii[j] *
                                                    (nt->p_i_it[j] +
                                                     p_n_it)) /
                                                   (nt->p_i_it[j] +
                                                    2. * p_n_it) +
                                                   /* our method */
                                                   nt->niihm[j] +
                                                   nt->sum_nij[j] +
                                                   2. * nt->nmm * (beta /
                                                                   (double)
                                                                   (beta +
                                                                    p_n_it_sqr
                                                                    -
                                                                    beta *
                                                                    p_n_it_sqr)
                                                                   *
                                                                   nt->
                                                                   p_i_it[j])
                );

        }
        sumt = 0.;
        for (j = 0; j < Data.num_alleles[nt->locus]; j++)
            sumt +=
                nt->nii[j] * (2 * p_n_it / (nt->p_i_it[j] + 2 * p_n_it)) +
                nt->niihm[j];

        p_n_it = 1 / (double)(2 * N) * (sumt +
                                        2 * nt->nmm *
                                        ((p_n_it_sqr - beta * p_n_it_sqr +
                                          beta * p_n_it) / (double)(beta +
                                                                    p_n_it_sqr
                                                                    -
                                                                    beta *
                                                                    p_n_it_sqr))
            );

        beta =
            1 / (double)N *(nt->nmm *
                            (beta /
                             (double)(beta + p_n_it_sqr * (1. - beta))));
        if (fabs(p_n_it - p_n_it_old) < D_NULLEPSILON)
            conv_events++;
        else
            conv_events = 0;

        if (conv_events > D_NULLNEPSILON)
            break;
    }
    FREE1D(has_hm);
    nt->p_n = p_n_it;
    nt->beta = beta;
    return p_n_it;
}

void
NULLdestroy(NULLTEST nt)
{
    int i;

    FREE1D(nt->p_i_it);
    FREE1D(nt->ni);
    FREE1D(nt->nii);
    FREE1D(nt->sum_nij);
    FREE1D(nt->niihm);
    FREE2D(nt->nij, Data.num_alleles[nt->locus], i);
    FREE(nt);
}

/* 
  vim: ft=c sw=4 ts=4 expandtab
*/
