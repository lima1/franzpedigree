/*
 * $Id: prob.c 2064 2010-05-26 11:20:32Z markus $
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
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "prob.h"

#include "global.h"
#include "pedigree.h"
#include "dataio.h"
#include "pvalue.h"
#include "utils.h"
#include "dag.h"
#include "lod.h"
#include "freq.h"
#include "listparentage.h"
#include "vtprogressbar.h"

#include "uthash.h"

extern DATA Data;
extern PROBS Probs;
extern OPTIONS Options;

static void leaveOnlyIntersection();

/* create a lookup where we can get the adjusted fullsib p-values */
typedef struct S_PV_ADJ_LUP
{
    double  pv;
    double  pv_adjusted;
    UT_hash_handle hh;
} PV_ADJ_LUP;

void
PROBinit(void)
{
    LODinit();
    Probs.simresults.pv_lup_po = NULL;
    Probs.simresults.pv_lup_hs = NULL;
    Probs.simresults.pv_lup_u = NULL;
}

/* sort posterior lists (containing the probs. of all PO pairs and 
 * triples) by different fields */
static int
comparePostsLOD(const void *a, const void *b)
{
    const POSTERIORS *da = (const POSTERIORS *)a;
    const POSTERIORS *db = (const POSTERIORS *)b;

    return (db->lod > da->lod) - (db->lod < da->lod);
}

static int
comparePostsPP(const void *a, const void *b)
{
    const POSTERIORS *da = (const POSTERIORS *)a;
    const POSTERIORS *db = (const POSTERIORS *)b;

    return (db->observed > da->observed) - (db->observed < da->observed);
}

static int
comparePostsPPNielsen(const void *a, const void *b)
{
    const POSTERIORS *da = (const POSTERIORS *)a;
    const POSTERIORS *db = (const POSTERIORS *)b;

    return (db->p_opt > da->p_opt) - (db->p_opt < da->p_opt);
}

int
PROBsearchPosteriors(POSTERIORS ** posteriors, int sample, int v, int w)
{
    int     i;

    if (v < w) SWAP(v,w,i);

    for (i = 0; i < Probs.num_posteriors[sample]; i++) {
        if (posteriors[sample][i].v == v &&
            posteriors[sample][i].w == w)
            return i;
    }
    return -1;
}

/* if user specified sampling location distances, then these two methods
 * filter pairs accordingly */
static  bool
isInMatingArea(int parent_id, int child_id)
{
    double  dist;

    if (!Data.use_distances)
        return true;

    dist =
        Data.population_distances[Data.id_mapping[parent_id].
                                  population->id][Data.id_mapping[child_id].
                                                  population->id];

    switch (Data.id_mapping[parent_id].sex) {
        case 1:
            if (Options.FemMaxDist < 0 || dist <= Options.FemMaxDist)
                return true;
            break;
        case 2:
            if (Options.MaleMaxDist < 0 || dist <= Options.MaleMaxDist)
                return true;
            break;
        default:
            if (Options.FemMaxDist < 0. || Options.MaleMaxDist < 0.
                || dist <= MAX(Options.FemMaxDist, Options.MaleMaxDist))
                return true;
            break;
    }

    return false;
}

static  bool
areParentsInMatingArea(int parent_id1, int parent_id2)
{
    if (!Data.use_distances || Options.ParentsMaxDist < 0)
        return true;

    if (Data.population_distances[Data.id_mapping[parent_id1].population->id]
        [Data.id_mapping[parent_id2].population->id]
        <= Options.ParentsMaxDist)
        return true;

    return false;
}

/* and here the age filter (the years of birth and death) */
static  bool
isReproductive(int sample, int child)
{
    int     min, max, age;

    /* without age data, we can't exclude sample as parent */
    if (Data.id_mapping[sample].birth < 0 || Data.id_mapping[child].birth < 0)
        return true;

    /* calculate the age at the birth of the child */
    age = Data.id_mapping[child].birth - Data.id_mapping[sample].birth;

    /* child would be older than sample */
    if (age < 0)
        return false;

    /* check the sex specific range */
    switch (Data.id_mapping[sample].sex) {
        case 1:
            min = Options.FemRepro.min;
            max = Options.FemRepro.max;
            break;
        case 2:
            min = Options.MaleRepro.min;
            max = Options.MaleRepro.max;
            break;
        default:
            min = MIN(Options.FemRepro.min, Options.MaleRepro.min);
            max = MAX(Options.FemRepro.max, Options.MaleRepro.max);
            break;
    }

    /* parent already dead at birth of child */
    if (Data.id_mapping[sample].death >= 0
        && age > Data.id_mapping[sample].death - Data.id_mapping[sample].birth)
        return false;

    if (min <= age && max >= age)
        return true;

    return false;
}

int
PROBcommonTypedLoci(int id1, int id2)
{
    int     i, common_typed_loci = 0;

    for (i = 0; i < Data.num_loci; i++)
        if (!(Data.id_mapping[id1].genotype_obs[i][0] < 0 ||
              Data.id_mapping[id1].genotype_obs[i][1] < 0 ||
              Data.id_mapping[id2].genotype_obs[i][0] < 0 ||
              Data.id_mapping[id2].genotype_obs[i][1] < 0))
            common_typed_loci++;

    return common_typed_loci;
}

/* here the public function that starts all tests */
bool
PROBisCandidateParent(int parent_id, int child_id)
{

    if (!isReproductive(parent_id, child_id))
        return false;
    if (!isInMatingArea(parent_id, child_id))
        return false;
    if (PROBcommonTypedLoci(parent_id, child_id) < Options.MinTyped)
        return false;

    return true;
}

/*for the sibling calculation, we do similar tests */
static  bool
siblingsAgeTest(int id1, int id2)
{
    int     age_diff =
        abs(Data.id_mapping[id1].birth - Data.id_mapping[id2].birth);

    if (Data.id_mapping[id1].birth < 0 || Data.id_mapping[id2].birth < 0)
        return true;

    /* _MinRepro contains the minimum repro range, typically FemRepro.max -
     * Femrepro.min. Set in cmdline.c. */
    if (age_diff <= Options._MinReproRange)
        return true;

    return false;
}

static  bool
siblingsDistTest(int id1, int id2)
{
    if (!Data.use_distances || Options.SibMaxDist < 0)
        return true;

    if (Data.population_distances[Data.id_mapping[id1].population->id]
        [Data.id_mapping[id2].population->id]
        > Options.SibMaxDist)
        return false;

    return true;
}

bool
PROBcouldBeSiblings(int id1, int id2)
{
    if (!siblingsAgeTest(id1, id2))
        return false;

    if (!siblingsDistTest(id1, id2))
        return false;

    if (PROBcommonTypedLoci(id1, id2) < Options.MinTyped)
        return false;

    return true;
}

/* ok, this is the monster function that determines all possible single parents 
 * and parent pairs of the child_id. It stores all triples and pairs in a
 * linked list, which we later (PROBcalcPosteriors) turn into a nice array... */
static void
calcParents(int child_id)
{
    int     i, j,               /* the two parent ids             */
            mismatching, mismatching_triple;    /* number of observed mismatches  */
    bool    known, both_known;  /* relationship known a priori?   */
    double  p, ldd, ldt;        /* likelihood of a relationship   */

    /* initialize structures for triples and dyads */
    Probs._dyads[child_id] = Probs._triples[child_id] = NULL;
    Probs._num_dyads[child_id] = Probs._num_triples[child_id] = 0;


    Probs.num_candidates[child_id] = 0;
    Probs.num_candidates_m[child_id] = Probs.num_candidates_f[child_id] = 0;

    /* the LOD denominators are always the same, so calc them only once */
    ldd = LODcalcDenominator(Data.id_mapping[child_id].genotype_obs, NULL);

    /* if user provided a pedigree with --pedigreein and child_id has an incoming
     * arc, then calculate H2. We skip the following test for mintyped
     * when we now the relationship  */
    if (HAS_KNOWN_MOTHER(child_id)) {
        i = Data.in_pedigree[child_id][0];
        ldt =
            LODcalcDenominatorTripleMotherKnown(Data.id_mapping
                                                [child_id].genotype_obs,
                                                Data.
                                                id_mapping[i].genotype_obs,
                                                NULL);
        goto mother;
    } else
        ldt = ldd;

    if (Data.id_mapping[child_id].typed_loci < Options.MinTyped)
        return;

    for (i = 0; i < Data.num_samples; i++) {

      mother:

        known = both_known = false;
        if (Data.use_pedigree) {
            if (Data.in_pedigree[child_id][0] >= 0) {
                if (Data.in_pedigree[child_id][0] != i)
                    break;
                else
                    known = true;
            }
            if (Data.in_pedigree[child_id][1] >= 0) {
                both_known = true;
                assert(known);
                j = Data.in_pedigree[child_id][1];
                goto father;
            }
        }

        if (child_id == i)      /* child cannot be its own parent */
            continue;

        if (!known) {
            if (Data.id_mapping[i].typed_loci < Options.MinTyped)
                continue;

            /* was putative parent at the time of birth reproductive? */
            if (!PROBisCandidateParent(i, child_id))
                continue;
        }

        p = LODcalcPchildDyad(Data.id_mapping[child_id].genotype_obs,
                              Data.id_mapping[i].genotype_obs,
                              &mismatching, known, NULL);

        if (!known && p > MINUSINF) {
            if (mismatching <= Options.MaxMismatchingDyad) {
                LISTPARENTAGEadd(&Probs._dyads[child_id], i, -1, p, p - ldd);
                Probs._num_dyads[child_id]++;
            } else {
                if (Data.id_mapping[i].sex == 1)
                    Probs._sum_filtered_likelihood_dyads_f[child_id] += exp(p);
                else
                    Probs._sum_filtered_likelihood_dyads_m[child_id] += exp(p);
            }
        } else if (known) {
            if (mismatching > Options.MaxMismatchingTriple)
                fprintf(stderr,
                        "Warning: many mismatches (%i) in known relationship (%s %s).\n",
                        mismatching, Data.id_mapping[i].description,
                        Data.id_mapping[child_id].description);
            LISTPARENTAGEadd(&Probs._dyads[child_id], i, -1, ldt, 0);
            Probs._num_dyads[child_id]++;
            j = 0;
            goto father;
        }

        Probs.candidates[child_id][Probs.num_candidates[child_id]++] = i;

        if (Data.has_sex_data && Data.id_mapping[i].sex > 0) {
            if (Data.id_mapping[i].sex == 1)
                Probs.num_candidates_f[child_id] += 2;
            else
                Probs.num_candidates_m[child_id] += 2;
        } else {
            Probs.num_candidates_f[child_id]++;
            Probs.num_candidates_m[child_id]++;
        }

        /* don't calculate triples with too many mismatches */
        if (mismatching > Options.MaxMismatchingTriple)
            continue;

        for (j = i; j < Data.num_samples; j++) {

          father:
            if (Data.use_pedigree && Data.in_pedigree[child_id][1] >= 0
                && Data.in_pedigree[child_id][1] != j)
                break;

            if (child_id == j)
                continue;
            if (!Options.Selfing && i == j)
                continue;
            if (!both_known) {
                if (Data.id_mapping[j].typed_loci < Options.MinTyped)
                    continue;
                if (!PROBisCandidateParent(j, child_id))
                    continue;

                /* distance father/mother too large? */
                if (!areParentsInMatingArea(i, j))
                    continue;

                if (Data.has_sex_data) {
                    /* if both parents have sex data and have the same sex,
                     * skip */
                    if (!(Options.Selfing && i == j)) {
                        if (Data.id_mapping[i].sex > 0
                            && Data.id_mapping[j].sex > 0
                            && Data.id_mapping[i].sex ==
                            Data.id_mapping[j].sex)
                            continue;
                    }
                }
            }

            p = LODcalcPchildTriple(Data.id_mapping[child_id].genotype_obs,
                                    Data.id_mapping[i].genotype_obs,
                                    Data.id_mapping[j].genotype_obs,
                                    &mismatching_triple, both_known, NULL);

            //printf("DEBUG: child_id = %s, i = %s, j = %s, %E, %i \n", Data.id_mapping[child_id].description,Data.id_mapping[i].description,Data.id_mapping[j].description,p, mismatching_triple);
            if (p > MINUSINF || both_known) {
                LISTPARENTAGEadd(&Probs._triples[child_id], i, j, p, p - ldt);
                Probs._num_triples[child_id]++;
            }

            if (known) {
                Probs.candidates[child_id][Probs.num_candidates[child_id]++] =
                    j;

                if (Data.has_sex_data && Data.id_mapping[j].sex > 0) {
                    if (Data.id_mapping[j].sex == 1)
                        Probs.num_candidates_f[child_id] += 2;
                    else
                        Probs.num_candidates_m[child_id] += 2;
                } else {
                    Probs.num_candidates_f[child_id]++;
                    Probs.num_candidates_m[child_id]++;
                }
            }
        }
    }

    if (!Options.Monoecious) {
        Probs.num_candidates_f[child_id] /= 2;
        Probs.num_candidates_m[child_id] /= 2;
    }

    if (Options.ndefined) {
        Probs.num_candidates_f[child_id] = Options.n;
        Probs.num_candidates_m[child_id] = Options.n;
    }
}

void
PROBcalcTriplesAndDyads(void)
{
    int     i, j, completed = 0;        /* the number of completed individuals, needed
                                           for parallel calculation                    */
    /* allocate space for the triples */
    MALLOC(Probs._triples, LISTPARENTAGENODE *, Data.num_samples);
    MAKE1DINT(Probs._num_triples, Data.num_samples);

    /* and for the dyads */
    MALLOC(Probs._dyads, LISTPARENTAGENODE *, Data.num_samples);
    MAKE1DINT(Probs._num_dyads, Data.num_samples);

    /* in contrast to CERVUS, we create the lists of candidate mothers     *
     * and fathers internally later. Here we allocated the space for them  */
    MAKE1DINT(Probs.num_candidates, Data.num_samples);
    MAKE1DINT(Probs.num_candidates_f, Data.num_samples);
    MAKE1DINT(Probs.num_candidates_m, Data.num_samples);
    MALLOC(Probs.candidates, int*, Data.num_samples);
    MAKE1DBOOL(Probs.is_candidate_parent, Data.num_samples);

    /* we will filter some parentages and we have to store their likelihood */
    MAKE1DDOUBLE(Probs._sum_filtered_likelihood_triples, Data.num_samples);
    MAKE1DDOUBLE(Probs._sum_filtered_likelihood_dyads_f, Data.num_samples);
    MAKE1DDOUBLE(Probs._sum_filtered_likelihood_dyads_m, Data.num_samples);

    Options.MaxMismatching =
        MAX(Options.MaxMismatchingDyad, Options.MaxMismatchingTriple) + 2;


#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(dynamic,4)
#endif
    for (i = 0; i < Data.num_samples; i++) {

#ifdef HAVE_OPENMP
        if (omp_get_thread_num() == 0 && Options.Verbosity >= 1)
#else
        if (Options.Verbosity >= 1)
#endif
            VTPROGRESSBARupdate("LOD Calculation",
                                Data.num_samples, completed);

        Probs._sum_filtered_likelihood_dyads_f[i] = 0.;
        Probs._sum_filtered_likelihood_dyads_m[i] = 0.;
        Probs._sum_filtered_likelihood_triples[i] = 0.;
        MAKE1DINT(Probs.candidates[i],Data.num_samples);
        calcParents(i);

        /* if there are no candidates or all others are candidates, we don't *
         * need this array and free the O(n^2) space                         */
        if (Probs.num_candidates[i] > 0 && Probs.num_candidates[i] < Data.num_samples - 1) {
            REALLOC(Probs.candidates[i], int, Probs.num_candidates[i]);
            assert(Probs.num_candidates[i] <= Data.num_samples);
        }
        else {
            FREE1D(Probs.candidates[i]);
        }

#ifdef HAVE_OPENMP
#pragma omp atomic
#endif
        completed++;
    }

    /* free the unused memory of the candidate lists */
    for (i = 0; i < Data.num_samples; i++) 
        if (Probs.num_candidates[i] == Data.num_samples - 1) {
            for (j = 0; j < Data.num_samples; j++) {
                if (i == j) continue;
                Probs.is_candidate_parent[j] = true;
            }
            /* everybody is a candidate now? then we can end this */
            if (Probs.is_candidate_parent[i]) 
                break;
        }
        else if (Probs.num_candidates[i] > 0) 
            for (j = 0; j < Probs.num_candidates[i]; j++)
                Probs.is_candidate_parent[Probs.candidates[i][j]] = true;

    if (Options.Verbosity >= 1)
        VTPROGRESSBARcompleted("LOD Calculation");
}

static void
intersectionPosts(LISTDAGNODE ** posts, int i, int j)
{
    int     found;
    LISTDAGNODE *tp, *tq, *intersection = NULL;

    for (tp = posts[i]; tp != NULL; tp = tp->next) {
        found = 0;
        for (tq = posts[j]; tq != NULL; tq = tq->next)
            if (Probs.posteriors[i][tp->v].v == Probs.posteriors[j][tq->v].v
                && Probs.posteriors[i][tp->v].w ==
                Probs.posteriors[j][tq->v].w) {
                found = 1;
                break;
            }
        if (found)
            LISTDAGadd(&intersection, tp->v);
    }
    LISTDAGdestroy(&posts[i]);
    posts[i] = intersection;
}

static void
recalcPosteriorProbs(int id, POSTERIORS ** posteriors, double Nf, double Nm)
{
    int     i;
    double  sum_l = 0., log_c_f, log_c_m, log_c_both, sum_r = 0.,max=MINUSINF;
    
    if (!Options.AssumeCompleteSample) {
    /* we need the number of unsampled individuals, not the total number */
    Nf = (double)MAX(1, (Nf - Probs.num_candidates_f[id]));
    Nm = (double)MAX(1, (Nm - Probs.num_candidates_m[id]));
    }
    else {
        Nf = Nm = 1;
    }

    log_c_f = log(Nf);
    log_c_m = log(Nm);
    
    if (Options.Monoecious)
        log_c_both = log((Nf * (Nf+1))/2);
    else
        log_c_both = log_c_f + log_c_m;
    
    /* first find best parentage */
    for (i = 0; i < Probs.num_posteriors[id]; i++) {
        posteriors[id][i].lw = posteriors[id][i].l;

        if (posteriors[id][i].v < 0 && posteriors[id][i].w < 0)
            posteriors[id][i].lw += log_c_both;
        else if (posteriors[id][i].w < 0)
            posteriors[id][i].lw +=
                Data.id_mapping[posteriors[id][i].v].sex == 1 ?
                log_c_m : log_c_f;

        if (!Options.IgnoreRamets && Data.num_ramets - Data.num_samples > 0) {
            if (posteriors[id][i].v >= 0) {
                posteriors[id][i].lw +=
                    log(Data.id_mapping[posteriors[id][i].v].ramets);
                sum_r += Data.id_mapping[posteriors[id][i].v].ramets;
            }    
            if (posteriors[id][i].w >= 0) {
                posteriors[id][i].lw +=
                    log(Data.id_mapping[posteriors[id][i].w].ramets);
                sum_r += Data.id_mapping[posteriors[id][i].w].ramets;
            }    
        }

        if (posteriors[id][i].lw > max) max = posteriors[id][i].lw;
    }
    
    /* now sum them up, relative to best parentage */
    for (i = 0; i < Probs.num_posteriors[id]; i++) {
        posteriors[id][i].lw = posteriors[id][i].l;

        if (posteriors[id][i].v < 0 && posteriors[id][i].w < 0)
            posteriors[id][i].lw += log_c_both;
        else if (posteriors[id][i].w < 0)
            posteriors[id][i].lw +=
                Data.id_mapping[posteriors[id][i].v].sex == 1 ?
                log_c_m : log_c_f;

        if (!Options.IgnoreRamets && Data.num_ramets - Data.num_samples > 0) {
            if (posteriors[id][i].v >= 0) {
                posteriors[id][i].lw +=
                    log(Data.id_mapping[posteriors[id][i].v].ramets);
                sum_r += Data.id_mapping[posteriors[id][i].v].ramets;
            }    
            if (posteriors[id][i].w >= 0) {
                posteriors[id][i].lw +=
                    log(Data.id_mapping[posteriors[id][i].w].ramets);
                sum_r += Data.id_mapping[posteriors[id][i].w].ramets;
            }    
        }

        sum_l += MAX(0,exp(posteriors[id][i].lw-max));
    }
    
    if (exp(max) > 0) {
        /* the likelihood of the filtered triples and dyads */
        sum_l += ( Probs._sum_filtered_likelihood_triples[id] +  Probs._sum_filtered_likelihood_dyads_f[id] * Nm +
                Probs._sum_filtered_likelihood_dyads_m[id] * Nf) / exp(max);
    }

    if (!Options.AssumeCompleteSample)
        for (i = 0; i < Probs.num_posteriors[id]; i++)
            posteriors[id][i].p_opt = log(exp(posteriors[id][i].lw - max) / sum_l);
    else 
        for (i = 0; i < Probs.num_posteriors[id]; i++)
            posteriors[id][i].p_opt = posteriors[id][i].l;

}

void
PROBrecalcPosteriorProbabilities(POSTERIORS ** posteriors, double Nf,
                                 double Nm)
{
    int     i;

    for (i = 0; i < Data.num_samples; i++)
        recalcPosteriorProbs(i, posteriors, Nf, Nm);
}

POSTERIORS **
PROBclonePosteriors(POSTERIORS ** posteriors)
{
    int     i, j;
    POSTERIORS **pcloned;

    MALLOC(pcloned, POSTERIORS *, Data.num_samples);
    for (i = 0; i < Data.num_samples; i++) {
        MALLOC(pcloned[i], POSTERIORS, Probs.num_posteriors[i]);
        for (j = 0; j < Probs.num_posteriors[i]; j++) {
            pcloned[i][j] = posteriors[i][j];
        }
    }
    return pcloned;
}

/* Fills the connections matrix. A 'true' indicates that there exists
 * a triple or pair with i and j being parent(s) and/or offspring.
 * this is cool to know when we change the genotype of one individual. Instead
 * of updating all LOD scores, we only update the "connected" LOD scores.
 */
static void
calcConnected()
{
    int     i, j, id_v, id_w;

    for (i = 0; i < Data.num_samples; i++) {

        Probs.are_connected[TIDX(i,i)] = true;

        for (j = 0; j < Probs.num_posteriors[i]; j++) {
            id_v = Probs.posteriors[i][j].v;
            id_w = Probs.posteriors[i][j].w;

            assert(id_v < Data.num_samples && id_w < Data.num_samples);

            if (id_v >= 0) 
                Probs.are_connected[TIDXL(i,id_v)] = true;
            if (id_w >= 0) {
                Probs.are_connected[TIDXL(i,id_w)] = true;
            }
        }
    }
}

/* this function filters out all parentages with negative LOD score AND a   *
 * very small posterior probability. They don't contribute much to the      *
 * posterior denominator so we can ignore them.                             */
static void
filterPosteriorProbabilities(void)
{
    int     i, j, k;
    for (i = 0; i < Data.num_samples; i++) {

        /* noting to filter ? */
        if (Probs.num_posteriors[i] == 0)
            continue;

        qsort(Probs.posteriors[i], (size_t) Probs.num_posteriors[i],
              sizeof(POSTERIORS), comparePostsLOD);

        /* filter all with negative LOD score */
        for (j = 0; j < Probs.num_posteriors[i]; j++) {
            if (Probs.posteriors[i][j].lod < -0.5) {
                /* sum the likelihood of the ones we will ignore */
                for (k = j; k < Probs.num_posteriors[i]; k++)
                    if (Probs.posteriors[i][k].v >= 0
                        && Probs.posteriors[i][k].w >= 0)
                        Probs._sum_filtered_likelihood_triples[i] +=
                            exp(Probs.posteriors[i][k].l);
                    else if (Probs.posteriors[i][k].v >= 0) {
                        if (Data.id_mapping[Probs.posteriors[i][k].v].sex == 1)
                            Probs._sum_filtered_likelihood_dyads_f[i] +=
                                exp(Probs.posteriors[i][k].l);
                        else
                            Probs._sum_filtered_likelihood_dyads_m[i] +=
                                exp(Probs.posteriors[i][k].l);
                    }
             /*   fprintf(stderr, "%s f: %E m: %E triples: %E\n",Data.id_mapping[i].description, Probs._sum_filtered_likelihood_dyads_f[i], Probs._sum_filtered_likelihood_dyads_m[i],
                        Probs._sum_filtered_likelihood_triples[i]);*/

                Probs.num_posteriors[i] = j;
                REALLOC(Probs.posteriors[i], POSTERIORS,
                        Probs.num_posteriors[i]);
                break;
            }
        }    
        qsort(Probs.posteriors[i], (size_t) Probs.num_posteriors[i],
              sizeof(POSTERIORS), comparePostsPPNielsen);
    }
}

void
PROBmakeIgnoreArrays(bool no_candidates)
{
    int     i, j, k;

    /* the array Data.ignore_ary[i][j] is true when the jth locus of the ith
     * individual should be ignored in all LOD calculations i (where i is the
     * offspring) */
#ifdef HAVE_OPENMP    
#pragma omp parallel for private(i,j,k)
#endif    
    for (i = 0; i < Data.num_samples; i++)
        for (j = 0; j < Data.num_loci; j++) {

            if (Data.num_alleles[j] < 2) {
                if (i == 0 && no_candidates)
                    fprintf(stderr, "Warning: Locus %s not informative.\n",
                            Data.loci_ids[j]);
                Data.ignore_ary[i][j] = true;
                continue;
            }
            Data.ignore_ary[i][j] = false;

            /* without missing data, use all */
            if (Data._num_missing_data == 0)
                continue;

            /* if offspring genotype lacks both alleles, ignore */
            if (Data.id_mapping[i].genotype_obs[j][0] < 0 &&
                Data.id_mapping[i].genotype_obs[j][1] < 0) {
                Data.ignore_ary[i][j] = true;
                continue;
            }

            /* i is a founder (and never offspring)? */
            if (no_candidates || Probs.num_candidates[i] < 1)
                continue;

            /* only use these loci, where at least one candidate father has no
             * missing alleles */
            Data.ignore_ary[i][j] = true;
            if (Probs.num_candidates[i] == Data.num_samples - 1) 
            for (k = 0; k < Data.num_samples; k++) {
                if (i != k && Data.id_mapping[k].genotype_obs[j][0] >= 0
                    && Data.id_mapping[k].genotype_obs[j][1] >= 0) {
                    Data.ignore_ary[i][j] = false;
                    break;
                }
            }
            else
            for (k = 0; k < Probs.num_candidates[i]; k++) {
                if (Data.id_mapping[Probs.candidates[i][k]].
                    genotype_obs[j][0] >= 0
                    && Data.id_mapping[Probs.candidates[i][k]].
                    genotype_obs[j][1] >= 0) {
                    Data.ignore_ary[i][j] = false;
                    break;
                }    
            }
        }
}

/* this will generate a memory efficient data structure out of the linked
 * lists with all the possible pairs and triples... */
void
PROBcalcPosteriors(void)
{
    int     i, j, h, it, id, ig, max_ci, num_triples, num_dyads, num_empty;
    double  cp[3], max_cp, cf, cm, ld;
    LISTPARENTAGENODE *nt;
    LISTPARENTAGENODE *nd;

    if (Options.Verbosity > 1 ) 
        fprintf(stderr, "\nSummarizing LOD results ...");

    MALLOC(Probs.posteriors, POSTERIORS *, Data.num_samples);
    MAKE1DINT(Probs.num_posteriors, Data.num_samples);
    MAKETRIANGULAR(Probs.are_connected, bool, Data.num_samples);

    for (i = 0; i < Data.num_samples; i++)
        for (j = 0; j < i; j++)
            Probs.are_connected[TIDX(i,j)] = false;

    for (i = 0; i < Data.num_samples; i++) {

        num_triples = Probs._num_triples[i];
        num_dyads = Probs._num_dyads[i];

        nt = Probs._triples[i];
        nd = Probs._dyads[i];

        /* we know that i has parents, so remove the no parents hypothesis */
        if (HAS_KNOWN_MOTHER(i))
            num_empty = 0;
        else
            num_empty = 1;

        Probs.num_posteriors[i] = num_triples + num_dyads + num_empty;
        MALLOC(Probs.posteriors[i], POSTERIORS, Probs.num_posteriors[i]);
        h = 0;
        it = 0;
        id = 0;
        ig = 0;
        max_ci = 0;
        ld = LODcalcDenominator(Data.id_mapping[i].genotype_obs, NULL);
	
	assert(ld > MINUSINF);

        /* merge all probabilities in one array, keep sorting */
        while (it < num_triples || id < num_dyads || ig < num_empty) {
            max_cp = cp[0] = cp[1] = cp[2] = MINUSINF;
            if (it < num_triples)
                cp[0] = nt->p;
            if (id < num_dyads)
                cp[1] = nd->p;
            if (ig < num_empty)
                cp[2] = ld;

            for (j = 0; j < 3; j++)
                if (cp[j] > max_cp) {
                    max_cp = cp[j];
                    max_ci = j;
                }
            /*printf("DEBUG: sample: %i, h: %i, max_ci: %i %f\n", i,h,max_ci, max_cp); */
            switch (max_ci) {
                case 0:
                    Probs.posteriors[i][h].l = nt->p;
                    if (nt->v < nt->w) {
                        Probs.posteriors[i][h].v = nt->w;
                        Probs.posteriors[i][h].w = nt->v;
                    } else {
                        Probs.posteriors[i][h].v = nt->v;
                        Probs.posteriors[i][h].w = nt->w;
                    }
                    Probs.posteriors[i][h].lod = nt->lod;
                    it++;
                    nt = nt->next;
                    break;
                case 1:
                    Probs.posteriors[i][h].l = nd->p;
                    Probs.posteriors[i][h].v = nd->v;
                    Probs.posteriors[i][h].w = -1;
                    Probs.posteriors[i][h].lod = nd->lod;
                    id++;
                    nd = nd->next;
                    break;
                case 2:
                    Probs.posteriors[i][h].l = cp[2];
                    Probs.posteriors[i][h].v = -1;
                    Probs.posteriors[i][h].w = -1;
                    Probs.posteriors[i][h].lod = 0;
                    ig++;
                    break;
            }
            h++;
        }
        LISTPARENTAGEdestroy(&Probs._triples[i]);
        LISTPARENTAGEdestroy(&Probs._dyads[i]);
    }

    /* now we don't need this linked lists anymore */
    FREE(Probs._triples);
    FREE(Probs._dyads);
    FREE(Probs._num_triples);
    FREE(Probs._num_dyads);

    if (Options.Verbosity > 1) 
        fprintf(stderr, " done.\n\nCalculating posterior probabilities ...");

    cf = cm = 0.;
    Probs.mean_candidates_f = Probs.mean_candidates_m = 0.;
    /* we will filter some parentages and we have to store their likelihood */

    for (i = 0; i < Data.num_samples; i++) {

        if (Probs.is_candidate_parent[i]) {
            if (Data.id_mapping[i].sex == 1)
                Probs.num_candidate_mothers += 2;
            else if (Data.id_mapping[i].sex == 2)
                Probs.num_candidate_fathers += 2;
            else {
                Probs.num_candidate_mothers++;
                Probs.num_candidate_fathers++;
            }
        }
        if (Probs.num_candidates_f[i] > 0) {
            cf++;
            Probs.mean_candidates_f += (double)Probs.num_candidates_f[i];
        }
        if (Probs.num_candidates_m[i] > 0) {
            cm++;
            Probs.mean_candidates_m += (double)Probs.num_candidates_m[i];
        }
    }

    if (cf > 0)
        Probs.mean_candidates_f /= cf;
    if (cm > 0)
        Probs.mean_candidates_m /= cm;

    Probs.num_candidate_mothers /= 2;
    Probs.num_candidate_fathers /= 2;

    if (Options.Verbosity > 1) 
        fprintf(stderr, " done.\n\nFinding informative loci ...");

    PROBmakeIgnoreArrays(false);

    if (Options.Verbosity > 1) 
        fprintf(stderr, " done.\n\n");

    PROBrecalcPosteriorProbabilities(Probs.posteriors, (double)Options.Nf,
                                     (double)Options.Nm);
    if (Options.Verbosity > 1) 
        fprintf(stderr, "Filtering posterior probabilities ...");
    filterPosteriorProbabilities();

    if (Options.Verbosity > 1) 
        fprintf(stderr, " done.\n\nCalculating connection matrix ...");

    calcConnected();
    if (Options.Fullsibs)
        PROBcalcFullsibs();

    if (Options.Fullsibs || strlen(Options.FullsibInfilename) > 0)
        leaveOnlyIntersection();

    PROBrecalcPosteriorProbabilities(Probs.posteriors, (double)Options.Nf,
                                     (double)Options.Nm);

    if (Options.Verbosity > 1) 
        fprintf(stderr, " done.\n\n");
}

static void
updateChild(int child_id, POSTERIORS ** posteriors, int ***genotypes)
{
    int     i, mismatching;
    double  lod;

    if (HAS_KNOWN_MOTHER(child_id)) {
        lod =
            LODcalcDenominatorTripleMotherKnown(genotypes[child_id],
                                                genotypes[Data.in_pedigree
                                                          [child_id][0]],
                                                Data.ignore_ary[child_id]);
    } else {
        lod =
            LODcalcDenominator(genotypes[child_id], Data.ignore_ary[child_id]);
    }
    for (i = 0; i < Probs.num_posteriors[child_id]; i++) {

        if (posteriors[child_id][i].v < 0)
            posteriors[child_id][i].l = lod;
        else {
            if (posteriors[child_id][i].w < 0) {
                posteriors[child_id][i].l =
                    LODcalcPchildDyad(genotypes[child_id],
                                      genotypes[posteriors[child_id][i].v],
                                      &mismatching, 1,
                                      Data.ignore_ary[child_id]);
                if (HAS_KNOWN_MOTHER(child_id))
                    posteriors[child_id][i].l = lod;

            } else {
                posteriors[child_id][i].l =
                    LODcalcPchildTriple(genotypes[child_id],
                                        genotypes[posteriors[child_id][i].v],
                                        genotypes[posteriors[child_id][i].w],
                                        &mismatching, 1,
                                        Data.ignore_ary[child_id]);
            }
        }
        posteriors[child_id][i].lod = posteriors[child_id][i].l - lod;
    }
}

static void
updateChildLocus(int child_id, POSTERIORS ** posteriors, int ***genotypes,
                 LOCUS_COORD lc, int old_allele)
{
    int     i;
#ifndef NDEBUG
    int     mismatching;
#endif
    double  lod;

    assert(Probs.num_posteriors[child_id] > 0);

    if (HAS_KNOWN_MOTHER(child_id)) {
        lod = posteriors[child_id][0].l - posteriors[child_id][0].lod +
            LODcalcDenominatorTripleMotherKnownDelta(genotypes, child_id,
                                                     Data.in_pedigree[child_id]
                                                     [0], lc, old_allele,
                                                     Data.ignore_ary
                                                     [child_id]);

        assert(CMP_DBL(lod,
                       LODcalcDenominatorTripleMotherKnown(genotypes[child_id],
                                                           genotypes
                                                           [Data.in_pedigree
                                                            [child_id][0]],
                                                           Data.ignore_ary
                                                           [child_id])));
    } else {
        /* TODO: one could also use delta functions for the denominators */
        lod = posteriors[child_id][0].l - posteriors[child_id][0].lod +
            LODcalcDenominatorDelta(genotypes, child_id, lc, old_allele,
                                    Data.ignore_ary[child_id]);
        assert(CMP_DBL
               (lod,
                LODcalcDenominator(genotypes[child_id],
                                   Data.ignore_ary[child_id])));
    }
    for (i = 0; i < Probs.num_posteriors[child_id]; i++) {
        if (posteriors[child_id][i].v < 0)
            posteriors[child_id][i].l = lod;
        else {
            if (posteriors[child_id][i].w < 0) {
                if (HAS_KNOWN_MOTHER(child_id))
                    posteriors[child_id][i].l = lod;
                else {
                    posteriors[child_id][i].l +=
                        LODcalcPchildDyadDelta(genotypes, child_id,
                                               posteriors[child_id][i].v, lc,
                                               old_allele,
                                               Data.ignore_ary[child_id]);
                    assert(CMP_DBL
                           (posteriors[child_id][i].l,
                            LODcalcPchildDyad(genotypes[child_id],
                                              genotypes[posteriors[child_id]
                                                        [i].v], &mismatching,
                                              1, Data.ignore_ary[child_id])));
                }

            } else {
                posteriors[child_id][i].l +=
                    LODcalcPchildTripleDelta(genotypes, child_id,
                                             posteriors[child_id][i].v,
                                             posteriors[child_id][i].w, lc,
                                             old_allele,
                                             Data.ignore_ary[child_id]);
                assert(CMP_DBL
                       (posteriors[child_id][i].l,
                        LODcalcPchildTriple(genotypes[child_id],
                                            genotypes[posteriors[child_id]
                                                      [i].v],
                                            genotypes[posteriors[child_id]
                                                      [i].w], &mismatching, 1,
                                            Data.ignore_ary[child_id])));
            }
        }
        posteriors[child_id][i].lod = posteriors[child_id][i].l - lod;
    }
}


void
PROBrecalcPosteriors(MCMC_STATE * state)
{
    int     i;

    for (i = 0; i < Data.num_samples; i++)
        updateChild(i, state->posteriors, state->genotypes);

    PROBrecalcPosteriorProbabilities(state->posteriors, state->Nf, state->Nm);
}

#ifndef NDEBUG
static  bool
updateGenotypeAssertion(MCMC_STATE * state)
{
    int     i, j;
    bool    ret = true;
    POSTERIORS **tmp;
    tmp = PROBclonePosteriors(state->posteriors);

    for (i = 0; i < Data.num_samples; i++)
        updateChild(i, tmp, state->genotypes);

    PROBrecalcPosteriorProbabilities(tmp, state->Nf, state->Nm);
    for (i = 0; i < Data.num_samples; i++) {
        for (j = 0; j < Probs.num_posteriors[i]; j++)
            if (tmp[i][j].v != state->posteriors[i][j].v ||
                tmp[i][j].w != state->posteriors[i][j].w ||
                !CMP_DBL(tmp[i][j].l, state->posteriors[i][j].l) ||
                !CMP_DBL(tmp[i][j].p_opt, state->posteriors[i][j].p_opt)
                ) {
                fprintf(stderr, "updateGenotypeAssertion failed: %i, %i\n", i,
                        j);
                ret = false;
            }
    }
    if (!ret) {
    PROBdumpPosteriors(stderr,tmp,NULL,NULL,2);
    PROBdumpPosteriors(stderr,state->posteriors,NULL,NULL,2);
    }    
    PROBdestroyPosteriors(tmp);
    return ret;
}
#endif

void
PROBupdateGenotype(MCMC_STATE * state, LOCUS_COORD lc, int old_allele)
{
    int     i;

    for (i = 0; i < Data.num_samples; i++)
        if (old_allele < 0) {
            updateChild(i, state->posteriors, state->genotypes);
            recalcPosteriorProbs(i, state->posteriors, state->Nf, state->Nm);
        }        
        else if (Probs.are_connected[TIDXL(i,lc.id)]) {
            updateChildLocus(i, state->posteriors, state->genotypes, lc,
                                 old_allele);
            recalcPosteriorProbs(i, state->posteriors, state->Nf, state->Nm);
        }
    /* we did the updating really fast, check if this was correct by
     * calculating it the easy way and then comparing the results      */
    assert(updateGenotypeAssertion(state));
}

static inline int
calcFullsibsOf(int i)
{
    int     j, k, l, count = 0, t;
    bool    common_posterior;
    char    msg[LINESIZE];
    IBD_LL *ibd;


    for (j = 0; j < i; j++) {
        t = TIDX(i,j);
        Data.fs_matrix[t].excluded = false;

        /* only detect siblings in offspring generation? */
        if (!Options.SiblingsInParentalGeneration
            && (Probs.num_candidates[i] <= 0
                && Probs.num_candidates[j] <= 0)) {
            Data.fs_matrix[t].excluded = true;
            continue;
        }

        /* check age data if i and j could be siblings */
        if (!PROBcouldBeSiblings(i, j)) {

            if (Data.fs_matrix[t].fullsibs) {
                snprintf(msg, sizeof(msg),
                         "Meta data conflict. %s and %s really fullsibs? Ignoring them now...",
                         Data.id_mapping[i].description,
                         Data.id_mapping[j].description);
                WARN(msg);
                Data.fs_matrix[t].fullsibs = false;
            }
            Data.fs_matrix[t].excluded = true;
            continue;
        }
        common_posterior = false;

        ibd = &Data.fs_matrix[t].ibd;

        IBDcalcRelationshipLikelihoods(Data.id_mapping[i].genotype_obs,
                                       Data.id_mapping[j].genotype_obs, ibd);
        /* get the p-Values */
        Data.fs_matrix[t].ibd_pv.po =
            PVALUEfind(Probs.simresults.pv_lup_po, ibd->fs - ibd->po);
        Data.fs_matrix[t].ibd_pv.hs =
            PVALUEfind(Probs.simresults.pv_lup_hs, ibd->fs - ibd->hs);
        Data.fs_matrix[t].ibd_pv.u =
            PVALUEfind(Probs.simresults.pv_lup_u, ibd->fs - ibd->u);

        if (Data.fs_matrix[t].fullsibs || ibd->fs < ibd->u
            || ibd->fs < ibd->po || ibd->fs < ibd->hs) {
            Data.fs_matrix[t].excluded = true;
            continue;
        }

        count++;


        for (k = 0; k < Probs.num_posteriors[i]; k++) {

            Data.fs_matrix[t].common_parents = 0;

            l = PROBsearchPosteriors(Probs.posteriors, j,
                                     Probs.posteriors[i][k].v,
                                     Probs.posteriors[i][k].w);
            if (l >= 0) {
                if (Probs.posteriors[i][k].v >= 0)
                    Data.fs_matrix[t].common_parents++;
                if (Probs.posteriors[i][k].w >= 0)
                    Data.fs_matrix[t].common_parents++;
                common_posterior = true;
                break;
            }
        }
        if (Data.fs_matrix[t].fullsibs && !common_posterior) {
            snprintf(msg, sizeof(msg),
                     "No common parents found. %s and %s really fullsibs? Ignoring them now...",
                     Data.id_mapping[i].description,
                     Data.id_mapping[j].description);
            WARN(msg);
            Data.fs_matrix[t].fullsibs = false;
        }
    }
    return count;
}

static void
dumpFullsibsCSV()
{
    int     i, j,t;
    IBD_LL *ibd;
    char    c;
    FILE   *fp = NULL;

    if (!Options.SiblingsOutfilename[2])
        return;
    FOPENW(fp, Options.SiblingsOutfilename[2]);
    fprintf(fp,
            "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n",
            "Genotype 1", "Genotype 2", "FS", "PO", "HS", "pV PO", "pV HS",
            "pV U", "Common Parents", "Common Loci", "Indirect");
    for (i = 0; i < Data.num_samples; i++)
        for (j = 0; j < i; j++) {
            t = TIDX(i,j);
            if (Data.fs_matrix[t].fullsibs) {
                ibd = &Data.fs_matrix[t].ibd;
                c = Data.fs_matrix[t].indirect ? '*' : ' ';
                /*        if (!Data.fs_matrix[i][j].excluded) {
                   if (ibd->fs < ibd->u || ibd->fs < ibd->po || ibd->fs < ibd->hs)
                   continue; */
                fprintf
                    (fp,
                     "%s,%s,%7.3f,%7.3f,%7.3f,%E,%E,%E,%i,%i,%c\n",
                     Data.id_mapping[i].description,
                     Data.id_mapping[j].description,
                     ibd->fs - ibd->u,
                     ibd->po - ibd->u,
                     ibd->hs - ibd->u,
                     Data.fs_matrix[t].ibd_pv_adjusted.po,
                     Data.fs_matrix[t].ibd_pv_adjusted.hs,
                     Data.fs_matrix[t].ibd_pv_adjusted.u,
                     Data.fs_matrix[t].common_parents,
                     PROBcommonTypedLoci(i, j), c);
            }
        }     
    FCLOSE(fp);
}

static void
dumpFullsibs()
{
    int     i, j, cnt = 0, t;
    IBD_LL *ibd;
    FILE   *fp = NULL;
    char    c;

    if (Options.SiblingsOutformat[0])
        DATAIOfullsibOut(Options.SiblingsOutfilename[0]);
    if (Options.SiblingsOutformat[2])
        dumpFullsibsCSV();
    if (!Options.SiblingsOutformat[1])
        return;
    FOPENW(fp, Options.SiblingsOutfilename[1]);
    DATAIOdumpDatasetDetails(fp);
    if (Options.SiblingsCorrectionMethod == 1)
        fprintf(fp, "p-Value Correction Method   : Benjamini-Hochberg");
    else
        fprintf(fp, "p-Value Correction Method   : Holm");
    fprintf(fp, "\np-Value Cutoff\n");
    fprintf(fp, "  No common parents         : %E\n", Options.THsiblings[0]);
    fprintf(fp, "  One common parent         : %E\n", Options.THsiblings[1]);
    fprintf(fp, "  Two common parents        : %E\n\n", Options.THsiblings[2]);
    fprintf(fp, "Null hypotheses:\n");
    fprintf(fp, "  Parent-Offspring (PO)     : %d\n", Options.FullsibAlternatives[0]);
    fprintf(fp, "  Halfsib (HS)              : %d\n", Options.FullsibAlternatives[1]);
    fprintf(fp, "  Unrelated (U)             : %d\n\n", Options.FullsibAlternatives[2]);

    fprintf(fp, "See %s for the typing error rates.\n\n", Options.Outfilename);
    fprintf(fp,
            "Pairs marked with a * are indirectly detected fullsibs (A/B,B/C FS -> A/C also FS)\n\n");

    fprintf(fp,
            " +------------+------------+---------+---------+---------+--------------+--------------+--------------+----------------+-------------+\n");
    fprintf(fp,
            " | %10s | %10s | %7s | %7s | %7s | %12s | %12s | %12s | %14s | %11s |\n",
            "Genotype 1", "Genotype 2", "FS", "PO", "HS", "pV PO", "pV HS",
            "pV U", "Common Parents", "Common Loci");
    fprintf(fp,
            " +------------+------------+---------+---------+---------+--------------+--------------+--------------+----------------+-------------+\n");
    for (i = 0; i < Data.num_samples; i++)
        for (j = i + 1; j < Data.num_samples; j++) {
            t = TIDX(j,i);
            if (Data.fs_matrix[t].fullsibs) {
                ibd = &Data.fs_matrix[t].ibd;
                cnt++;
                c = Data.fs_matrix[t].indirect ? '*' : ' ';
                /*        if (!Data.fs_matrix[i][j].excluded) {
                   if (ibd->fs < ibd->u || ibd->fs < ibd->po || ibd->fs < ibd->hs)
                   continue; */
                fprintf
                    (fp,
                     " |%c%10s | %10s | %7.3f | %7.3f | %7.3f | %E | %E | %E | %14i | %11i |\n",
                     c, Data.id_mapping[i].description,
                     Data.id_mapping[j].description, ibd->fs - ibd->u,
                     ibd->po - ibd->u, ibd->hs - ibd->u,
                     Data.fs_matrix[t].ibd_pv_adjusted.po,
                     Data.fs_matrix[t].ibd_pv_adjusted.hs,
                     Data.fs_matrix[t].ibd_pv_adjusted.u,
                     Data.fs_matrix[t].common_parents,
                     PROBcommonTypedLoci(i, j));
            }
        }    
    fprintf(fp,
            " +------------+------------+---------+---------+---------+--------------+--------------+--------------+----------------+-------------+\n");
    fprintf(fp, "\nTOTAL: %i\n\nRejected Dyads:\n", cnt);
    cnt = 0;
    fprintf(fp,
            " +------------+------------+---------+---------+---------+--------------+--------------+--------------+----------------+-------------+\n");
    for (i = 0; i < Data.num_samples; i++)
        for (j = i + 1; j < Data.num_samples; j++) {
            t = TIDX(j,i);
            if (!Data.fs_matrix[t].excluded
                && !Data.fs_matrix[t].fullsibs) {
                cnt++;
                ibd = &Data.fs_matrix[t].ibd;
                /*        if (!Data.fs_matrix[i][j].excluded) {
                   if (ibd->fs < ibd->u || ibd->fs < ibd->po || ibd->fs < ibd->hs)
                   continue; */
                c = Data.fs_matrix[t].indirect ? '*' : ' ';
                fprintf
                    (fp,
                     " |%c%10s | %10s | %7.3f | %7.3f | %7.3f | %E | %E | %E | %14i | %11i |\n",
                     c,
                     Data.id_mapping[i].description,
                     Data.id_mapping[j].description,
                     ibd->fs - ibd->u,
                     ibd->po - ibd->u,
                     ibd->hs - ibd->u,
                     Data.fs_matrix[t].ibd_pv_adjusted.po,
                     Data.fs_matrix[t].ibd_pv_adjusted.hs,
                     Data.fs_matrix[t].ibd_pv_adjusted.u,
                     Data.fs_matrix[t].common_parents,
                     PROBcommonTypedLoci(i, j));
            }
            }
    fprintf(fp,
            " +------------+------------+---------+---------+---------+--------------+--------------+--------------+----------------+-------------+\n");
    fprintf(fp, "\nTOTAL (REJECTED!): %i\n", cnt);
    FCLOSE(fp);
}

static PV_ADJ_LUP *
addPV(PV_ADJ_LUP * pv_adj_lup, double pv, double pv_adjusted, PV_ADJ_LUP * s)
{
    s->pv = pv;
    s->pv_adjusted = pv_adjusted;
    HASH_ADD(hh, pv_adj_lup, pv, sizeof(double), s);
    return pv_adj_lup;
}

static inline PV_ADJ_LUP *
findPV(PV_ADJ_LUP * pv_adj_lup, double pv)
{
    PV_ADJ_LUP *s;
    HASH_FIND(hh, pv_adj_lup, &pv, sizeof(double), s);
    if (s)
        return s;
    FATALINT("findPV");
}

static void
adjustPvalues(int n)
{
    int     i, j, k = 0, l[3],t;
    double **p, **p_adjusted, p_last;
    PV_ADJ_LUP *pv_adj_lup[3], *s, **ss;

    MAKE2DDOUBLE(p, 3, n, i);
    MAKE2DDOUBLE(p_adjusted, 3, n, i);

    for (i = 0; i < Data.num_samples; i++)
        for (j = 0; j < i; j++) {
            t = TIDX(i,j);
            if (!Data.fs_matrix[t].excluded) {
                p[0][k] = Data.fs_matrix[t].ibd_pv.po;
                p[1][k] = Data.fs_matrix[t].ibd_pv.hs;
                p[2][k] = Data.fs_matrix[t].ibd_pv.u;
                k++;
            }
        }    

    assert(k == n);

    MALLOC(ss, PV_ADJ_LUP *, 3);

#ifdef HAVE_OPENMP
#pragma omp parallel for private(i, j, p_last)
#endif
    for (i = 0; i < 3; i++) {
        pv_adj_lup[i] = NULL;
        p_last = -1.;
        l[i] = 0;

        if (Options.SiblingsCorrectionMethod == 1)
            UTILSpvAdjustBH(p[i], p_adjusted[i], n);
        else
            UTILSpvAdjustHolm(p[i], p_adjusted[i], n);

        MALLOC(ss[i], PV_ADJ_LUP, n);
        for (j = 0; j < n; j++) {
            if (p[i][j] == p_last)
                continue;
            pv_adj_lup[i] =
                addPV(pv_adj_lup[i], p[i][j], p_adjusted[i][j], &ss[i][l[i]]);
            p_last = p[i][j];
            l[i]++;
        }
        HASH_SORT(pv_adj_lup[i], UTILScompare_doubles_decr);
    }

    for (i = 0; i < Data.num_samples; i++)
        for (j = 0; j < i; j++) {
            t = TIDX(i,j);
            if (!Data.fs_matrix[t].excluded) {
                s = findPV(pv_adj_lup[0], Data.fs_matrix[t].ibd_pv.po);
                Data.fs_matrix[t].ibd_pv_adjusted.po = s->pv_adjusted;

                s = findPV(pv_adj_lup[1], Data.fs_matrix[t].ibd_pv.hs);
                Data.fs_matrix[t].ibd_pv_adjusted.hs = s->pv_adjusted;

                s = findPV(pv_adj_lup[2], Data.fs_matrix[t].ibd_pv.u);
                Data.fs_matrix[t].ibd_pv_adjusted.u = s->pv_adjusted;
            }
        }

    FREE2D(p, 3, i);
    FREE2D(p_adjusted, 3, i);
    for (i = 0; i < 3; i++)
        while (pv_adj_lup[i]) {
            s = pv_adj_lup[i];  /* grab pointer to first item */
            HASH_DEL(pv_adj_lup[i], s); /* delete it (users advances to next) */
        }
    FREE2D(ss, 3, i);
}

static void
checkIndirectFullsib(int i, int j, int k, bool remove, bool * changed)
{
    int     itm, jtm,ijtm;
    double  th;
    IBD_LL  ibd;
    itm = TIDXL(i,k);
    jtm = TIDXL(j,k);
    ijtm = TIDXL(i,j);
    if (Data.fs_matrix[jtm].fullsibs &&
        !Data.fs_matrix[itm].fullsibs) {
        if (remove && Data.fs_matrix[itm].excluded) {
            /*fprintf(stderr, "debug: %s %s %s %i %i\n",
                    Data.id_mapping[i].description,
                    Data.id_mapping[j].description,
                    Data.id_mapping[k].description, im, kim);*/
            ibd = Data.fs_matrix[ijtm].ibd_pv_adjusted;
            th = Options.THsiblings[Data.fs_matrix[ijtm].common_parents] /
                D_TH_SIBLING_IND;
           // fprintf(stderr, "debug: %f %f %f %f  %f\n", ibd.po, ibd.fs,ibd.hs,ibd.u,th);
            if (ibd.po > th || ibd.hs > th || ibd.u > th) {
                Data.fs_matrix[ijtm].fullsibs = false;
                Data.fs_matrix[ijtm].indirect = true;
                *changed = true;
                return;
            }
            ibd = Data.fs_matrix[jtm].ibd_pv_adjusted;
            th = Options.THsiblings[Data.fs_matrix[jtm].common_parents] /
                D_TH_SIBLING_IND;
            //fprintf(stderr, "debug: %f %f %f %f  %f\n", ibd.po, ibd.fs,ibd.hs,ibd.u,th);
            if (ibd.po > th || ibd.hs > th || ibd.u > th) {
                Data.fs_matrix[jtm].fullsibs = false;
                Data.fs_matrix[jtm].indirect = false;
                *changed = true;
                return;
            }
            /*
            ibd = Data.fs_matrix[im][kim].ibd;
            fprintf(stderr, "debug: %f %f %f %f\n", ibd.fs, ibd.po, ibd.hs, ibd.u);
            ibd = Data.fs_matrix[im][kim].ibd_pv;
            fprintf(stderr, "debug: %f %f %f %f\n", ibd.fs, ibd.po, ibd.hs, ibd.u);
            // Data.fs_matrix[i][j].excluded = true;*/
            return;
        }
        if (!remove) {
            Data.fs_matrix[itm].fullsibs = true;
            Data.fs_matrix[itm].indirect = true;
            *changed = true;
        }
    }
}


static void
calcIndirectFullsibs(void)
{
    int     i, j, k, l = 0,t, count = 0;
    bool    changed[2];
        for (i = 0; i < Data.num_samples; i++)
            for (j = 0; j < i; j++) {
                t = TIDX(i,j);
                if (Data.fs_matrix[t].fullsibs
                    && !Data.fs_matrix[t].indirect)
                    for (k = 0; k < Data.num_samples; k++) {
                        if (k == i || k == j) continue;
                        /* j and k also FS? then i and k FS */
                        checkIndirectFullsib(i, j, k, true, &changed[0]);
                        /* i and k also FS? then j and k FS */
                        checkIndirectFullsib(j, i, k, true, &changed[1]);
                    }

            }
    changed[0] = changed[1] = true;
    while ((changed[0] || changed[1]) && l < 100) {
        l++;
        changed[0] = changed[1] = false;
    //    fprintf(stderr, "%i iter\n", l);
        for (i = 0; i < Data.num_samples; i++)
            for (j = 0; j < i; j++) 
                if (Data.fs_matrix[TIDX(i,j)].fullsibs)
                    for (k = 0; k < Data.num_samples; k++) {
                        if (k == i || k == j) continue;
                        /* j and k also FS? then i and k FS */
                        checkIndirectFullsib(i, j, k, false, &changed[0]);
                        /* i and k also FS? then j and k FS */
                        checkIndirectFullsib(j, i, k, false, &changed[1]);
                    }
    }
    for (i = 0; i < Data.num_samples; i++)
        for (j = 0; j < i; j++) {
            t = TIDX(i,j);
            if (Data.fs_matrix[t].indirect) count++;
        }        
    if (count>Data.num_samples*0.1) WARN("Many indirect fullsib assignments. Try running without --fullsibtest");
}

void
PROBcalcFullsibs(void)
{
    int     i, j, n, completed = 0, count = 0, t;
    double  th;

    DATAIOinitFullsibMatrix();

    n = Data.num_samples * (Data.num_samples + 1) / 2;

    /* now calculate for each individual the sibling probabilities */
#ifdef HAVE_OPENMP
#pragma omp parallel for reduction(+:count)
#endif
    for (i = 0; i < Data.num_samples; i++) {
#ifdef HAVE_OPENMP
        if (omp_get_thread_num() == 0 && Options.Verbosity >= 1)
#else
        if (Options.Verbosity >= 1)
#endif
            VTPROGRESSBARupdate("IBD Calculation", n, completed);

        count = count + calcFullsibsOf(i);

#ifdef HAVE_OPENMP
#pragma omp atomic
#endif
        completed += Data.num_samples - i;
    }

    if (Options.Verbosity >= 1)
        VTPROGRESSBARcompleted("IBD Calculation");

    /* adjust the observed p-Values with the Benjamini-Hochberg method */
    adjustPvalues(count);

    /* for pairs we could not exclude as sibling pair (because of age for
     * example) we compaire the p-Values with the user defined 
     * (--fdrsiblings) thresholds                                         */
    for (i = 0; i < Data.num_samples; i++)
        for (j = 0; j < i; j++) {
            t = TIDX(i,j);
            if (!Data.fs_matrix[t].excluded) {
                /* this is one is an easy false positive we can fix */
                if ((HAS_KNOWN_MOTHER(i) || HAS_KNOWN_MOTHER(j)) &&
                    Data.in_pedigree[i][0] != Data.in_pedigree[j][0]
                    && Data.fs_matrix[t].common_parents == 0)
                    continue;

                th = Options.THsiblings[Data.fs_matrix[t].common_parents];

                if (
                       /* skip Parent-Offspring Test? */
                       ((!isReproductive(i, j) || !isReproductive(j, i)) ||
                        Data.fs_matrix[t].ibd_pv_adjusted.po <= th) &&
                       Data.fs_matrix[t].ibd_pv_adjusted.hs <= th &&
                       Data.fs_matrix[t].ibd_pv_adjusted.u <= th)
                    Data.fs_matrix[t].fullsibs = true;

                if (Data.fs_matrix[t].fullsibs)
                    Data.fs_matrix[t].halfsibs = false;
            }
        }
    calcIndirectFullsibs();
    /* output the results */
    dumpFullsibs();
}


/* This function removes candidate parent pairs of fullsibs that they don't
 * have in common. For example if A and B are fullsibs and C and D is a parent
 * pair only for A and not for B, we remove C/D from the posterior list of A 
 */
static void
leaveOnlyIntersection(void)
{
    int     i, j, k, l, max_posteriors = 0, *cposs, fs_parents;
    bool    found;
    POSTERIORS tmp;
    LISTDAGNODE **posts, *tp;

    MALLOC(posts, LISTDAGNODE *, Data.num_samples);
    for (i = 0; i < Data.num_samples; i++)
        posts[i] = NULL;

    for (i = 0; i < Data.num_samples; i++) {
        if (Probs.num_posteriors[i] > max_posteriors)
            max_posteriors = Probs.num_posteriors[i];
        for (j = 0; j < i; j++) {


            /* we are only interested in fullsibs here (allow only common    *
             * pairs)                                                        */
            if (Data.fs_matrix[TIDX(i,j)].fullsibs) {
                for (k = 0; k < Probs.num_posteriors[i]; k++) {
                    /* search posterior k of i in j's posterior list         */
                    l = PROBsearchPosteriors(Probs.posteriors, j,
                                             Probs.posteriors[i][k].v,
                                             Probs.posteriors[i][k].w);

                    /* found? ok, this could be the correct parentage        *
                     * (fullsibs have the same parents, you know...          */
                    if (l >= 0) {
                        if (LISTDAGsearch(&posts[i], k) == NULL)
                            LISTDAGadd(&posts[i], k);
                        if (LISTDAGsearch(&posts[j], l) == NULL)
                            LISTDAGadd(&posts[j], l);
                    }
                }
            } else if (Data.fs_matrix[TIDX(i,j)].halfsibs) {
                if (posts[i] == NULL)
                    for (k = 0; k < Probs.num_posteriors[i]; k++)
                        LISTDAGadd(&posts[i], k);

                for (k = 0; k < Probs.num_posteriors[i]; k++) {
                    found = false;
                    for (l = 0; l < Probs.num_posteriors[j]; l++) {

                        if ((Probs.posteriors[i][k].w < 0 ||
                             (Probs.posteriors[i][k].v ==
                              Probs.posteriors[j][l].v
                              || Probs.posteriors[i][k].v ==
                              Probs.posteriors[j][l].w)
                             || (Probs.posteriors[i][k].w ==
                                 Probs.posteriors[j][l].v
                                 || Probs.posteriors[i][k].w ==
                                 Probs.posteriors[j][l].w))) {
                            found = true;
                            continue;
                        }
                    }
                    if (!found && LISTDAGsearch(&posts[i], k) != NULL)
                        LISTDAGremove(LISTDAGsearch(&posts[i], k));
                }
            }
        }
    }

    /* if we know that A and B are fullsibs and B and C, then A and C are    *
     * also fullsibs, so we have to build these intersections. 6 iterations  *
     * should be more than enough...                                         */
    for (k = 0; k < 6; k++)
        for (i = 0; i < Data.num_samples; i++)
            for (j = 0; j < i; j++)
                if (Data.fs_matrix[TIDX(i,j)].fullsibs) {
                    intersectionPosts(posts, i, j);
                    intersectionPosts(posts, j, i);
                }

    MAKE1DINT(cposs, max_posteriors);
    for (i = 0; i < Data.num_samples; i++) {
        if (posts[i] == NULL)
            continue;

        for (j = 0; j < Probs.num_posteriors[i]; j++)
            cposs[j] = j;

        fs_parents = LISTDAGcnt(&posts[i]);
        Probs.num_posteriors[i] = fs_parents;

        for (tp = posts[i]; tp != NULL; tp = tp->next) {
            fs_parents--;
            /* we always swap posteriors so that the position of all in 
             * posts[i] is smaller than num_posteriors */
            if (cposs[tp->v] != fs_parents) {
                tmp = Probs.posteriors[i][fs_parents];
                Probs.posteriors[i][fs_parents] =
                    Probs.posteriors[i][cposs[tp->v]];
                Probs.posteriors[i][cposs[tp->v]] = tmp;
                cposs[fs_parents] = cposs[tp->v];
                cposs[tp->v] = fs_parents;
            }
        }
        LISTDAGdestroy(&posts[i]);
    }

    FREE(cposs);
    FREE(posts);
}

static void
printID(FILE * fp, int id)
{
    if (id < 0)
        fprintf(fp, ",");
    else
        fprintf(fp, "%s,%i", Data.id_mapping[id].description,
                Data.id_mapping[id].typed_loci);
}

double
parentPosterior(int i,int v, POSTERIORS ** posteriors,int k) {
    int j;
    double pp, sum = 0;

        for (j = 0; j < Probs.num_posteriors[i]; j++) {
            if (Probs.mh_sampled_pedigrees < 1)
                pp = exp(posteriors[i][j].p_opt);
            else
                pp = posteriors[i][j].observed /
                    (double)Probs.mh_sampled_pedigrees;
            if (v < 0) {
                if (k == 1 && Probs.posteriors[i][j].v == v) sum += pp;
                else if (k == 2 && Probs.posteriors[i][j].w == v) sum += pp;
            } else {
                if (Probs.posteriors[i][j].v == v) sum += pp;
                else if (Probs.posteriors[i][j].w == v) sum += pp;
            }
        }
        return sum;
}

void
PROBdumpPosteriors(FILE * fp, POSTERIORS ** posteriors, Dag D, Dag D_correct,
                   int outformat)
{
    int     i, j, k, post_id_dag, post_id_dag_correct, cnt_common, mismatching;
    bool    assigned, *ignore_ary;
    char    chosen;
    double  pp, d1, dg;

    fprintf(fp,
            "Offspring,Loci Typed,Parent 1,Loci Typed,Parent 2,Loci Typed,LOD,Posterior,Common Loci Typed,Mismatches,n_f,n_m,Pair LOD Parent 1,Pair LOD Parent 2,Posterior Parent 1,Posterior Parent 2\n");

    for (i = 0; i < Data.num_samples; i++) {
        if (Probs.num_candidates[i] <= 0)
            continue;

        if (Data.ignore_ary != NULL)
            ignore_ary = Data.ignore_ary[i];
        else
            ignore_ary = NULL;


        if (Probs.mh_sampled_pedigrees < 1)
            qsort(posteriors[i], (size_t) Probs.num_posteriors[i],
                  sizeof(POSTERIORS), comparePostsPPNielsen);
        else
            qsort(posteriors[i], (size_t) Probs.num_posteriors[i],
                  sizeof(POSTERIORS), comparePostsPP);

        post_id_dag = PEDIGREEgetPostId(posteriors, D, i);
        post_id_dag_correct = PEDIGREEgetPostId(posteriors, D_correct, i);
        if (D_correct != NULL && post_id_dag_correct < 0)
            fprintf(stderr,
                    "Warning: Correct parents of sample %i (%s) not considered!\n",
                    i, Data.id_mapping[i].description);

        assigned = false;
        for (j = 0; j < Probs.num_posteriors[i]; j++) {
            chosen = ' ';
            if (D != NULL && post_id_dag == j)
                chosen = '<';
            if (j >= Probs.num_posteriors[i])
                chosen = 'X';
            if (D_correct != NULL && post_id_dag_correct == j
                && post_id_dag != post_id_dag_correct)
                chosen = '!';

            /* no mcmc necessary ? */
            if (Probs.mh_sampled_pedigrees < 1)
                pp = exp(posteriors[i][j].p_opt);
            else
                pp = posteriors[i][j].observed /
                    (double)Probs.mh_sampled_pedigrees;


            if (j == 0 || ( outformat == OUTFORMAT_DETAIL && pp > 0.05)) {      /* sorted, so first one is the most-likeliest */
                assigned = true;
                cnt_common = 0;
                for (k = 0; k < Data.num_loci; k++) {
                    if (ignore_ary != NULL && ignore_ary[k])
                        continue;
                    cnt_common++;
                }
                printID(fp, i);
                fprintf(fp, ",");
                printID(fp, posteriors[i][j].v);
                fprintf(fp, ",");
                printID(fp, posteriors[i][j].w);
                mismatching =
                    LODcalcMismatches(i, posteriors[i][j].v,
                                      posteriors[i][j].w, ignore_ary);
                fprintf(fp, ",%E,%.4f,%i,%i,%i,%i,", posteriors[i][j].lod, pp,
                        cnt_common, mismatching, Probs.num_candidates_f[i],
                        Probs.num_candidates_m[i]);
                dg = LODcalcDenominator(Data.id_mapping[i].genotype_obs,
                                        ignore_ary);
                if (posteriors[i][j].v >= 0) {
                    d1 = LODcalcPchildDyad(Data.id_mapping[i].genotype_obs,
                                           Data.id_mapping[posteriors[i][j].v].
                                           genotype_obs, &mismatching, 1,
                                           ignore_ary);
                    fprintf(fp, "%E,", d1 - dg);
                } else
                    fprintf(fp, ",");
                if (posteriors[i][j].w >= 0) {
                    d1 = LODcalcPchildDyad(Data.id_mapping[i].genotype_obs,
                                           Data.id_mapping[posteriors[i][j].w].
                                           genotype_obs, &mismatching, 1,
                                           ignore_ary);
                    fprintf(fp, "%E", d1 - dg);
                }
                fprintf(fp, ",%.4f,%.4f", parentPosterior(i,posteriors[i][j].v,posteriors,1),
                                          parentPosterior(i,posteriors[i][j].w,posteriors,2));

                if (chosen != ' ')
                    fprintf(fp, ",%c", chosen);
                fprintf(fp, "\n");

            }
        }
        if (!assigned) {
            printID(fp, i);
            fprintf(fp, ",");

            if (Data.use_pedigree && Data.in_pedigree[i][0] > 0)
                printID(fp, Data.in_pedigree[i][0]);

            fprintf(fp, ",");

            if (Data.use_pedigree && Data.in_pedigree[i][1] > 0)
                printID(fp, Data.in_pedigree[i][1]);

            fprintf(fp, ",,,,%i,%i\n",
                    Probs.num_candidates_f[i], Probs.num_candidates_m[i]);
        }
    }
}


void
PROBdestroyPosteriors(POSTERIORS ** p)
{
    int     i;
    FREE2D(p, Data.num_samples, i);
}

void
PROBdestroyPreMCMC(void)
{
    int     i;
    /* some maybe already free'd */
    for (i = 0; i < Data.num_samples; i++)
        if (Probs.candidates[i] != NULL)
            FREE1D(Probs.candidates[i]);

    FREE1D(Probs.candidates);
    FREE1D(Probs.is_candidate_parent);
    if (Probs.simresults.pv_lup_po != NULL)
        PVALUEdestroy(Probs.simresults.pv_lup_po);
    if (Probs.simresults.pv_lup_hs != NULL)
        PVALUEdestroy(Probs.simresults.pv_lup_hs);
    if (Probs.simresults.pv_lup_u != NULL)
        PVALUEdestroy(Probs.simresults.pv_lup_u);
}

void
PROBdestroyPostMCMC(void)
{

    LODdestroy();
    PROBdestroyPosteriors(Probs.posteriors);
    FREE(Probs._sum_filtered_likelihood_dyads_f);
    FREE(Probs._sum_filtered_likelihood_dyads_m);
    FREE(Probs._sum_filtered_likelihood_triples);

    FREE(Probs.num_candidates);
    FREE(Probs.num_candidates_f);
    FREE(Probs.num_candidates_m);
    FREE(Probs.num_posteriors);
    FREETRIANGULAR(Probs.are_connected);
}

/* 
  vim: ft=c sw=4 ts=4 expandtab
*/
