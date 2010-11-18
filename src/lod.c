/*
 * $Id: lod.c 1930 2010-02-16 14:46:27Z markus $
 * 
 * Formulas (with typos) see Appendix of
 *
 * S.T. Kalinowski, M.L. Taper, and T.C. Marshall. Revising how the 
 * computer program CERVUS accommodates genotyping error increases 
 * success in paternity assignment. Mol. Ecol., 16:1099--1106, Mar 2007.
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

#include "macros.h"

#include <stdlib.h>             /* for free */
#include <math.h>               /* for log */
#include <assert.h>

#include "lod.h"
#include "global.h"
#include "freq.h"
#include "prob.h"

extern DATA Data;
extern OPTIONS Options;

/* some constants for the likelihood formulas */
static short epsilonsmallocd = 0;       /* already calculated? */

/* for parent-pairs  (triples) */
static double *epsilon0;
static double *epsilon1;
static double *epsilon2;

/* for single-parents (dyads) */
static double *epsilond0;
static double *epsilond3;


#define IGNORE_LOCI { if ( ignore_ary != NULL) {                             \
                         if (ignore_ary[i]) continue;                        \
                      }                                                      \
                      else if (child[i][0] < 0 && child[i][1] < 0) continue;  }
#define IGNORE_LOCI_UPDATE if ( ignore_ary != NULL && ignore_ary[lc.locus]) return 0.;

/* for loci without missing alleles we can skip many if statements in the
 * transition prob. */
typedef double (*ptTransDyad)(int[],int[],int);
typedef double (*ptTransTriple)(int[],int[],int[], int);

static double calcTransDyad(int child[], int parent[], int locus);
static double calcTransDyadMissing(int child[], int parent[], int locus);

static double calcTransTriple(int child[], int p1[], int p2[], int locus);
static double calcTransTripleMissing(int child[], int p1[], int p2[], int locus);

static ptTransDyad   *transdyadfncs;
static ptTransTriple *transtriplefncs;

void
LODinit()
{
    int     i;

    if (epsilonsmallocd)
        return;
    epsilonsmallocd = 1;
    /* typingerror for each locus */
    MAKE1DDOUBLE(epsilon0, Data.num_loci);
    MAKE1DDOUBLE(epsilon1, Data.num_loci);
    MAKE1DDOUBLE(epsilon2, Data.num_loci);

    MAKE1DDOUBLE(epsilond0, Data.num_loci);
    MAKE1DDOUBLE(epsilond3, Data.num_loci);

    MALLOC(transdyadfncs,ptTransDyad,Data.num_loci);
    MALLOC(transtriplefncs,ptTransTriple,Data.num_loci);

    for (i = 0; i < Data.num_loci; i++) {
        if (1 || Data.locus_has_missing_data[i]) {
            transdyadfncs[i]   = &calcTransDyadMissing;
            transtriplefncs[i] = &calcTransTripleMissing;
        }    
        else {
            transdyadfncs[i]   = &calcTransDyad;
            transtriplefncs[i] = &calcTransTriple;
        }    
    }
    LODrecalcErrorConstants();
    /* ignore loci where both alleles are missing */
    PROBmakeIgnoreArrays(true);

}

void
LODrecalcErrorConstants() {
    int i;
    for (i = 0; i < Data.num_loci; i++) {
        epsilon0[i] = (1. - Data.TypingError[i]) * (1. - Data.TypingError[i]) *
            (1. - Data.TypingError[i]);
        epsilon1[i] =
            Data.TypingError[i] * (1. - Data.TypingError[i]) *
            (1. - Data.TypingError[i]);

        epsilon2[i] =
            3 * Data.TypingError[i] * Data.TypingError[i] * (1. -
                                                             Data.
                                                             TypingError[i]);

        epsilon2[i] =
             Data.TypingError[i] * Data.TypingError[i] * (3. - 2. * 
                                                             Data.
                                                             TypingError[i]);

        epsilond0[i] = (1. - Data.TypingError[i]) * (1. - Data.TypingError[i]);

        epsilond3[i] = Data.TypingError[i] * (2. - Data.TypingError[i]);
    }
}

static inline double
compareGenotypes(int g1[], int g2[], int locus)
{
    if ((g1[0] == g2[0] && g1[1] == g2[1]) ||
        (g1[0] == g2[1] && g1[1] == g2[0]))
        return 1.;

    if (g2[0] < 0) {
        if (g2[1] < 0)
            return LODcalcGenotypeProb(g1, locus);
        if (g1[0] == g2[1])
            return FREQ_ALLELE_PROB(locus, g1[1]);
        if (g1[1] == g2[1])
            return FREQ_ALLELE_PROB(locus, g1[0]);
    } else if (g2[1] < 0) {
        if (g1[0] == g2[0])
            return FREQ_ALLELE_PROB(locus, g1[1]);
        if (g1[1] == g2[0])
            return FREQ_ALLELE_PROB(locus, g1[0]);
    }
    return 0.;
}

static double
calcTransTripleMissing(int child[], int p1[], int p2[], int locus)
{    
    int     g[2], im, tmp;
    double  p;

    /* both alleles of an putative parent missing? */
    
    if (p1[0] < 0 && p1[1] < 0)
        return LODcalcTransProbDyad(child, p2, locus);

    if (p2[0] < 0 && p2[1] < 0)
        return LODcalcTransProbDyad(child, p1, locus);

    /* first one always missing */
    if (p1[1] < 0)
        SWAP(p1[0], p1[1], tmp);
    if (p2[1] < 0)
        SWAP(p2[0], p2[1], tmp);

    /* one allele of child genotype missing? */
    if (child[0] < 0 || child[1] < 0) {
        im = MAXARG(child[0], child[1]);
            p = FREQ_CMP_ALLELES_2(locus,child[im],
                                    p1[0]) + FREQ_CMP_ALLELES_2(locus,child[im],
                                                            p1[1]);

            p += FREQ_CMP_ALLELES_2(locus,child[im],
                                    p2[0]) + FREQ_CMP_ALLELES_2(locus,child[im],
                                                            p2[1]);

        return p / 4.;
    }

    /* loci with both missing values should be ignored */
    assert(child[0] >= 0 && child[1] >= 0);
    p = 0.;
    g[0] = p1[0];
    g[1] = p2[0];
    p += compareGenotypes(child, g, locus);
    g[0] = p1[1];
    g[1] = p2[0];
    p += compareGenotypes(child, g, locus);

    g[0] = p1[0];
    g[1] = p2[1];
    p += compareGenotypes(child, g, locus);

    /* no missing values anymore */
    g[0] = p1[1];
    g[1] = p2[1];
    if ((child[0] == g[0] && child[1] == g[1]) ||
        (child[0] == g[1] && child[1] == g[0]))
        p++;

    return p / 4;
}

static double
calcTransTriple(int child[], int p1[], int p2[], int locus)
{   
    int g[2];
    double p;
    p = 0.;
    g[0] = p1[0];
    g[1] = p2[0];
    if ((child[0] == g[0] && child[1] == g[1]) ||
        (child[0] == g[1] && child[1] == g[0]))
        p++;
    g[0] = p1[1];
    g[1] = p2[0];
    if ((child[0] == g[0] && child[1] == g[1]) ||
        (child[0] == g[1] && child[1] == g[0]))
        p++;

    g[0] = p1[0];
    g[1] = p2[1];
    if ((child[0] == g[0] && child[1] == g[1]) ||
        (child[0] == g[1] && child[1] == g[0]))
        p++;

    g[0] = p1[1];
    g[1] = p2[1];
    if ((child[0] == g[0] && child[1] == g[1]) ||
        (child[0] == g[1] && child[1] == g[0]))
        p++;

    return p / 4;
}

double
LODcalcTransProbTriple(int child[], int p1[], int p2[], int locus)
{
    return transtriplefncs[locus](child,p1,p2,locus);
}

inline int
LODcmpSetsDyad(int allele, int p[])
{
    int r = 0;
    if (allele == p[0]) r = 1;
    if (allele == p[1]) r++; 
    return r;
}

static double
calcTransDyad(int child[], int parent[], int locus)
{    
    int     ps[2] = { -1, -1};
    ps[0] = LODcmpSetsDyad(child[0], parent);

    switch (ps[0]) {
        case 2:
            /* ai.aj ai.ai ? p = p(aj) */
            return FREQ_ALLELE_PROB(locus, child[1]);
        case 1:
            if (child[0] == child[1]) {
                /* ai.ai ai.aj */

                return 0.5 * FREQ_ALLELE_PROB(locus, child[0]);
            } else {
                ps[1] = LODcmpSetsDyad(child[1], parent);
                /* ai.aj ai.ak */
                //fprintf(stderr, "DEBUG %d\n",ps[1]);
                return 0.5 * ( FREQ_ALLELE_PROB(locus, child[1]) +
                    ps[1] *  FREQ_ALLELE_PROB(locus, child[0]));
            }

        default:
            ps[1] = LODcmpSetsDyad(child[1], parent);
            return ps[1] * 0.5 * FREQ_ALLELE_PROB(locus, child[0]);
    }
}

static double
calcTransDyadMissing(int child[], int parent[], int locus)
{
    double  p;
    int     im, tp[2];


    if (parent[0] < 0 || parent[1] < 0) {
        /* ai.aj ?.? */
        if (parent[0] == parent[1])
            return LODcalcGenotypeProb(child, locus);

        /* ai.aj ak.? */
        if (( child[0] >= 0 && child[1] >= 0)) {
            im = MAXARG(parent[0], parent[1]);
            assert(parent[im] >= 0);
            tp[0] = child[0] == parent[im] ? 1 : 0;
            //tp[0] = FREQ_CMP_ALLELES_2(locus,child[0],parent[im]);

            if (child[0] == child[1]) {
                /* case 1: sampled parent allele is segregated */
                p = 0.5 * (tp[0] * FREQ_ALLELE_PROB(locus, child[1])) +
                    /* case 2: missing allele is the segregated one */
                    0.5 * LODcalcGenotypeProb(child, locus);
            } else {
                tp[1] = child[1] == parent[im] ? 1 : 0;
              //  tp[1] = FREQ_CMP_ALLELES_2(locus,child[1],parent[im]);
                
                /* case 1: sampled parent allele is segregated */
                p = 0.5 * ((tp[0] * FREQ_ALLELE_PROB(locus, child[1]) +
                            tp[1] * FREQ_ALLELE_PROB(locus, child[0])) +

                    /* case 2: missing allele is the segregated one */
                           LODcalcGenotypeProb(child, locus));
            }
            return p;
        }
   }

    if (child[0] < 0 || child[1] < 0) {
        im = MAXARG(child[0], child[1]);
        assert(child[im] >= 0);
        p = 0.5 * FREQ_ALLELE_PROB(locus, child[im]) +
            0.25 * 
                    (FREQ_CMP_ALLELES_2(locus,child[im], parent[0]) +
                    FREQ_CMP_ALLELES_2(locus,child[im], parent[1]));
        return p;
    }
    return calcTransDyad(child, parent, locus);
}

double
LODcalcTransProbDyad(int child[], int parent[], int locus)
{
    return transdyadfncs[locus](child,parent,locus);
}

double
LODcalcGenotypeProb(int genotype[], int locus)
{
    double  p[2];

    p[0] = FREQ_ALLELE_PROB(locus, genotype[0]);
    p[1] = FREQ_ALLELE_PROB(locus, genotype[1]);
    /* homozygotes (or one missing, p is then 1, see freq.h */
    if (FREQ_CMP_ALLELES(genotype[0], genotype[1]))
        return p[0] * p[1];

    /* heterozygotes */
    return 2 * p[0] * p[1];
}

/* The functions that calculate L(H1) and L(H2), new CERVUS paper */
inline static double
calcPchildTripleLocus(int locus, int **child, int **p1, int **p2, double tpt)
{
    double  tpd1, tpd2, pg;
    tpd1 = LODcalcTransProbDyad(child[locus], p1[locus], locus);
    tpd2 = LODcalcTransProbDyad(child[locus], p2[locus], locus);
    pg = LODcalcGenotypeProb(child[locus], locus);
    assert(epsilon0[locus] > 0);
    return log(epsilon0[locus] * tpt + epsilon1[locus] * (tpd1 + tpd2 + pg) +
               epsilon2[locus] * pg);
}

inline static double
calcPchildDyadLocus(int locus, int **child, double tpd)
{
    double  pg;
    pg = LODcalcGenotypeProb(child[locus], locus);
    return log(epsilond0[locus] * tpd + pg * epsilond3[locus]);
}

inline static double
calcDenominatorLocus(int locus, int **genotype)
{
    return log(LODcalcGenotypeProb(genotype[locus], locus));
}

inline static double
calcDenominatorTripleMotherKnownLocus(int locus, int **genotype, int **mother)
{
    double  tp, gp;
    tp = LODcalcTransProbDyad(genotype[locus], mother[locus], locus);
    gp = LODcalcGenotypeProb(genotype[locus], locus);
    return log(epsilon0[locus] * tp + epsilon1[locus] * (tp + 2 * gp) +
               epsilon2[locus] * gp);
}

/* The functions that calculate the multilocus genotype probs. */
double
LODcalcPchildTriple(int **child, int **p1, int **p2, int *mismatching,
                    bool known, bool * ignore_ary)
{
    int     i;
    double  tpt, pl = 0;

    *mismatching = 0;
    for (i = 0; i < Data.num_loci; i++) {
        IGNORE_LOCI;
        tpt = LODcalcTransProbTriple(child[i], p1[i], p2[i], i);
        if (CMP_DBL(tpt, 0.)) {
            (*mismatching)++;
            if (!known && *mismatching > Options.MaxMismatchingTriple)
                return MINUSINF;
        }
        pl += calcPchildTripleLocus(i, child, p1, p2, tpt);
    }
    return pl;
}

/* A fast variant that re-calcs this Prob. after change of one allele */
double
LODcalcPchildTripleDelta(int ***genotypes, int child_id, int p1_id, int p2_id,
                         LOCUS_COORD lc, int old_allele, bool * ignore_ary)
{
    int     new_allele = genotypes[lc.id][lc.locus][lc.allele];
    double  tpt, delta;
    IGNORE_LOCI_UPDATE;
    if (FREQ_CMP_ALLELES(old_allele, new_allele))
        return 0.;
    if (child_id != lc.id && p1_id != lc.id && p2_id != lc.id)
        return 0.;

    genotypes[lc.id][lc.locus][lc.allele] = old_allele;
    tpt =
        LODcalcTransProbTriple(genotypes[child_id][lc.locus],
                               genotypes[p1_id][lc.locus],
                               genotypes[p2_id][lc.locus], lc.locus);
    delta =
        calcPchildTripleLocus(lc.locus, genotypes[child_id], genotypes[p1_id],
                              genotypes[p2_id], tpt);

    genotypes[lc.id][lc.locus][lc.allele] = new_allele;
    tpt =
        LODcalcTransProbTriple(genotypes[child_id][lc.locus],
                               genotypes[p1_id][lc.locus],
                               genotypes[p2_id][lc.locus], lc.locus);
    return calcPchildTripleLocus(lc.locus, genotypes[child_id],
                                 genotypes[p1_id], genotypes[p2_id],
                                 tpt) - delta;
}

double
LODcalcPchildDyad(int **child, int **parent, int *mismatching, bool known,
                  bool * ignore_ary)
{
    int     i;
    double  pl = 0., tpd;

    *mismatching = 0;
    for (i = 0; i < Data.num_loci; i++) {
        IGNORE_LOCI;

        tpd = LODcalcTransProbDyad(child[i], parent[i], i);

        if (CMP_DBL(tpd, 0.)) {
            (*mismatching)++;
            /* skip the calculations when we have oberved too many
             * mismatches */
            if (!known && *mismatching > Options.MaxMismatching)
                return MINUSINF;
        }

        pl += calcPchildDyadLocus(i, child, tpd);
    }

    return pl;
}

double
LODcalcPchildDyadDelta(int ***genotypes, int child_id, int p_id,
                       LOCUS_COORD lc, int old_allele, bool * ignore_ary)
{
    int     new_allele = genotypes[lc.id][lc.locus][lc.allele];
    double  tpd, delta;
    IGNORE_LOCI_UPDATE;
    if (FREQ_CMP_ALLELES(old_allele, new_allele))
        return 0.;
    if (child_id != lc.id && p_id != lc.id)
        return 0.;

    genotypes[lc.id][lc.locus][lc.allele] = old_allele;
    tpd =
        LODcalcTransProbDyad(genotypes[child_id][lc.locus],
                             genotypes[p_id][lc.locus], lc.locus);
    delta = calcPchildDyadLocus(lc.locus, genotypes[child_id], tpd);

    genotypes[lc.id][lc.locus][lc.allele] = new_allele;
    tpd =
        LODcalcTransProbDyad(genotypes[child_id][lc.locus],
                             genotypes[p_id][lc.locus], lc.locus);
    return calcPchildDyadLocus(lc.locus, genotypes[child_id], tpd) - delta;
}


double
LODcalcDenominator(int **child, bool * ignore_ary)
{
    int     i;
    double  ld = 0.;

    for (i = 0; i < Data.num_loci; i++) {
        IGNORE_LOCI;
        ld += calcDenominatorLocus(i, child);
    }
    return ld;
}

double
LODcalcDenominatorTripleMotherKnownDelta(int ***genotypes, int child_id,
                                         int mother_id, LOCUS_COORD lc,
                                         int old_allele, bool * ignore_ary)
{
    int     new_allele = genotypes[lc.id][lc.locus][lc.allele];
    double  delta;
    IGNORE_LOCI_UPDATE;
    if (FREQ_CMP_ALLELES(old_allele, new_allele))
        return 0.;
    if (child_id != lc.id && mother_id != lc.id)
        return 0.;
    genotypes[lc.id][lc.locus][lc.allele] = old_allele;
    delta =
        calcDenominatorTripleMotherKnownLocus(lc.locus, genotypes[child_id],
                                              genotypes[mother_id]);

    genotypes[lc.id][lc.locus][lc.allele] = new_allele;
    return calcDenominatorTripleMotherKnownLocus(lc.locus, genotypes[child_id],
                                                 genotypes[mother_id]) - delta;
}

double
LODcalcDenominatorTripleMotherKnown(int **child, int **mother,
                                    bool * ignore_ary)
{
    int     i;
    double  ld = 0.;

    for (i = 0; i < Data.num_loci; i++) {
        IGNORE_LOCI;
        ld += calcDenominatorTripleMotherKnownLocus(i, child, mother);
    }
    return ld;
}

double
LODcalcDenominatorDelta(int ***genotypes, int child_id, LOCUS_COORD lc,
                        int old_allele, bool * ignore_ary)
{
    int     new_allele = genotypes[lc.id][lc.locus][lc.allele];
    double  delta;
    IGNORE_LOCI_UPDATE;
    if (FREQ_CMP_ALLELES(old_allele, new_allele) || child_id != lc.id)
        return 0.;
    genotypes[lc.id][lc.locus][lc.allele] = old_allele;
    delta = calcDenominatorLocus(lc.locus, genotypes[child_id]);

    genotypes[lc.id][lc.locus][lc.allele] = new_allele;
    return calcDenominatorLocus(lc.locus, genotypes[child_id]) - delta;
}

bool
LODhasMismatchLocus(int locus, int **child, int v, int w, bool *ignore_ary) {
    double tp;
    if ( ignore_ary != NULL) {            
        if (ignore_ary[locus]) return false;
    } 
    else if (child[locus][0] < 0 && child[locus][1] < 0) return false;  

    if (v < 0 && w < 0) return false;
    if (w < 0) {
        tp = LODcalcTransProbDyad(child[locus], Data.id_mapping[v].genotype_obs[locus], locus);
        if (CMP_DBL(tp,0)) return true;
    }
    else {
        tp = LODcalcTransProbTriple(child[locus], Data.id_mapping[v].genotype_obs[locus], Data.id_mapping[w].genotype_obs[locus], locus);
        if (CMP_DBL(tp,0)) return true;
    }        
    return false;
}

int LODcalcMismatches(int child_id, int v, int w, bool *ignore_ary)
{
    int i, mismatching = 0, **child = Data.id_mapping[child_id].genotype_obs;
    for (i = 0; i < Data.num_loci; i++) {
        if (LODhasMismatchLocus(i, child,v,w,ignore_ary)) mismatching++;
    }
    return mismatching;
}

void
LODdestroy()
{
    FREE1D(epsilon0);
    FREE1D(epsilon1);
    FREE1D(epsilon2);
    FREE1D(epsilond0);
    FREE1D(epsilond3);
    FREE1D(transdyadfncs);
    FREE1D(transtriplefncs);
    epsilonsmallocd = 0;
}

/* 
  vim: ft=c sw=4 ts=4 expandtab
*/
