/*
 * $Id: mcmc.c 2064 2010-05-26 11:20:32Z markus $
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

#define _POSIX_C_SOURCE 1
#include "macros.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>

#include "mcmc.h"
#include "genotype.h"
#include "dataio.h"
#include "freq.h"
#include "pedigree.h"
#include "prob.h"
#include "lod.h"
#include "utils.h"
#include "vtprogressbar.h"
#include "options.h"
#include "exact.h"

extern PROBS Probs;
extern DATA Data;
extern OPTIONS Options;

#define POSTERIORSCORE(i,j) state->posteriors[i][j].p_opt
static void gibbsChangeGenotype(MCMC_STATE * state, LOCUS_COORD lc);

static void fillMissingData(MCMC * mcmc, MCMC_STATE * state,
                            unsigned int *completed, unsigned int total);

static void createGibbsErrorRates(MCMC_STATE *state);

/* every thread stores its best pedigree found so far. this *
 * function initializes the data structure *best            */
static void
initBest(MCMC * mcmc, STATE_SNAPSHOT * best)
{
    best->p = MINUSINF;
    best->D = DAGinit(Data.num_samples, 2, mcmc->no_cycles);
    best->Nf = best->Nm = 1.0;

    if (mcmc->genotypes_fixed)
        best->genotypes = NULL;
    else
        best->genotypes = DATAIOcloneGenotypes();
}

/* and this function frees the memory                       */
static void
destroyBest(MCMC * mcmc, STATE_SNAPSHOT * best)
{
    DAGdestroy(best->D);
    if (!mcmc->genotypes_fixed)
        DATAIOdestroyGenotypes(best->genotypes);
}

MCMC
MCMCinit()
{
    int     i, j;
    unsigned int expected_observations;
    MCMC    mcmc;
    LOCUS_COORD lc;

    mcmc.D_true = NULL;
    mcmc.c0 = 1.;
    mcmc.mh_acceptr = 0.;
    mcmc.swaps = mcmc.swaps_accepted = mcmc.freqs_updated = 0;

    mcmc.genotypes_fixed = mcmc.posteriors_fixed = true;

    if (Options.GibbsMissingData && Data._num_missing_data > 0) {
        mcmc.genotypes_fixed = false;
        /* every thread gets its own genotypes */
        mcmc.genotypes = NULL;
    } else {
        /* only one central genotype array, best pedigree as the same
         * genotype */
        mcmc.genotypes = DATAIOcloneGenotypes();
    }

    /* if we change N or genotypes, every thread needs his own posteriors  */
    if (!(Options.Nmdefined && Options.Nfdefined) || !mcmc.genotypes_fixed)
        mcmc.posteriors_fixed = false;

    mcmc.no_cycles = true;

    /* without age data or if individuals with same year of birth could be *
     * parent and offspring, we have to check for cycles                   */
    if (Data.has_missing_birth_data
        || MIN(Options.FemRepro.min, Options.MaleRepro.min) < 1)
        mcmc.no_cycles = false;

    mcmc.p_true = MINUSINF;

    initBest(&mcmc, &mcmc.best);


    /* initialize the arrays that store the paramaters extracted from the 
     * sampled pedigrees */
    expected_observations = Options.MHIterations / Options.MHsamplefreq;
    MAKE1DDOUBLE(mcmc.mh_Nf, expected_observations);
    MAKE1DDOUBLE(mcmc.mh_Nm, expected_observations);
    MAKE1DDOUBLE(mcmc.mh_sr_f, expected_observations);
    MAKE1DDOUBLE(mcmc.mh_sr_m, expected_observations);
    MAKE1DDOUBLE(mcmc.mh_p, expected_observations);
    if (Options.GibbsTypingError) {
        MAKE2DDOUBLE(mcmc.mh_te, Data.num_loci,expected_observations,i);
    }
    else 
        mcmc.mh_te = NULL;
    
    if (Data.num_populations > 1) {
        MALLOC(mcmc.mh_pop, unsigned int*, expected_observations);
        for (i=0;i<expected_observations;i++) 
            MAKETRIANGULAR(mcmc.mh_pop[i], unsigned int, Data.num_populations);
    }
    else mcmc.mh_pop = NULL;

    if (Data.num_ramets != Data.num_samples) {
        MAKE1DDOUBLE(mcmc.mh_clonal, expected_observations);
    } else
        mcmc.mh_clonal = NULL;

    /* we also might log the sampled missing alleles */
    mcmc.log_missing_alleles = false;
    mcmc.mh_missing_alleles = NULL;
    if (strlen(Options.MissingAllelesOutfilename) > 0 && !mcmc.genotypes_fixed) {
        MALLOC(mcmc.mh_missing_alleles, int *, Data._num_missing_data);
        for (i = 0; i < Data._num_missing_data; i++) {
            lc = Data._missing_data[i];
            MAKE1DINT(mcmc.mh_missing_alleles[i], Data.num_alleles[lc.locus]);
            for (j = 0; j < Data.num_alleles[lc.locus]; j++)
                mcmc.mh_missing_alleles[i][j] = 0;
        }
        mcmc.log_missing_alleles = true;
    }

    if (Options.Selfing) {
        MAKE1DDOUBLE(mcmc.mh_selfing, expected_observations);
    } else
        mcmc.mh_selfing = NULL;

    /* we only need to sample new parents for offspring with unknown 
     * parents */
    MAKE1DINT(mcmc._to_change, Data.num_samples);

    mcmc._num_to_change = 0;
    mcmc.N = 0;
    for (i = 0; i < Data.num_samples; i++) {
        if (Probs.num_posteriors[i] > 1) {
            mcmc.N += (Probs.num_posteriors[i] - 1);
            mcmc._to_change[mcmc._num_to_change++] = i;
        }
    }
    mcmc.N_posteriors = mcmc.N;

    if (mcmc._num_to_change > 0)
        REALLOC(mcmc._to_change, int, mcmc._num_to_change);

    mcmc.N_fixed = true;
    if (!Options.Nfdefined || !Options.Nmdefined) {
        mcmc.N_fixed = false;
        if (Options.Nfdefined)
            mcmc.N += D_MCMCCPGRID;
        if (Options.Nmdefined)
            mcmc.N += D_MCMCCPGRID;
    }

    if (!mcmc.genotypes_fixed)
        mcmc.N += Data._num_missing_data;

    return mcmc;
}

inline static void
calcNfm(MCMC_STATE * state)
{

    if (Options.Nfdefined)
        state->Nf = Options.Nf;
    else
        state->Nf =
            MIN(Options.Nfmax,
                (Probs.mean_candidates_f / state->sr_f * state->x_f));

    if (Options.Nmdefined)
        state->Nm = Options.Nm;
    else
        state->Nm =
            MIN(Options.Nmmax,
                (Probs.mean_candidates_m / state->sr_m * state->x_m));
}

MCMC_STATE
MCMCinitState(MCMC * mcmc, Dag D_start)
{
    int     i, j;

    MCMC_STATE state;

    state.converged = state.unheated = false;
    state.iter = state.accepted = 0;
    /* D is our directed acylic graph we are working with */
    if (D_start == NULL)
        state.D = DAGinit(Data.num_samples, 2, mcmc->no_cycles);
    else
        state.D = DAGcopy(D_start);

    state._addd = state._rmvd = 0;

    state.x_f = state.x_m = state.sr_f = state.sr_m = 1.;

    state.c = 1.;
    state.all_p = NULL;
    state._p_chi_temp = -1;
    state._thread_id = 0;
    (void)time(&state.begin);
    if (mcmc->posteriors_fixed) {
        state.posteriors = Probs.posteriors;
    } else {
        state.posteriors = PROBclonePosteriors(Probs.posteriors);
        for (i = 0; i < Data.num_samples; i++)
            for (j = 0; j < Probs.num_posteriors[i]; j++)
                state.posteriors[i][j].observed = 0;
    }

    /* we also change genotypes, so give every thread its own genotypes */
    if (mcmc->genotypes_fixed)
        state.genotypes = mcmc->genotypes;
    else
        state.genotypes = DATAIOcloneGenotypes();

    /* some space we need to calculate stuff */
    MALLOC(state._edgesrmvd, Edge, 4);
    MALLOC(state._edgesaddd, Edge, 4);
    MALLOC(state._edgesall, Edge, Data.num_samples * 2);
    MAKE1DDOUBLE(state._p_gibbs_genotype,
                 (Data.max_num_alleles * (Data.max_num_alleles + 1) / 2));

    if (Options.GibbsTypingError) {
        MAKE2DDOUBLE(state._gibbs_error_rate, Data.num_loci, D_MCMCTEGRID, i);
        MAKE1DDOUBLE(state._p_gibbs_error, D_MCMCTEGRID);
        createGibbsErrorRates(&state);
    }

    if (D_start == NULL) {
        state.p = 0.;
        calcNfm(&state);
    } else {
        PEDIGREEestmateSamplingRate(state.D, &state.sr_f, &state.sr_m);
        calcNfm(&state);
        PROBrecalcPosteriors(&state);
        state.p = MCMCscoreDag(state.D, &state);
    }

    initBest(mcmc, &state.best);

    state.freqs_updated = 0;
    return state;
}

/* For each individual, Probs.posteriors[id] contains a list of valid incoming
 * arcs. This function looks at the incoming arcs and returns the id of this
 * in Probs.posteriors[id]. Returns -1 if it is not in the list (which
 * shouldn't happen .... */
static int inline
incomingHypo(Dag D, MCMC_STATE * state, int id)
{
    int     indegree, v = -1, w = -1;
    Edge    incoming[2];

    indegree = DAGincomingE(D, id, incoming);
    if (indegree > 0) {
        v = incoming[0].v;
        if (indegree == 2)
            w = incoming[1].v;
    }
    return PROBsearchPosteriors(state->posteriors, id, v, w);
}


/* This fills all missing data by Gibbs-Sampling */
static void
fillMissingData(MCMC * mcmc, MCMC_STATE * state, unsigned int *completed,
                unsigned int total)
{
    unsigned int i;
    bool    last_both = false;
    LOCUS_COORD lc;

    if (mcmc->genotypes_fixed)
        return;
    for (i = 0; i < Data._num_missing_data; i++) {
        lc = Data._missing_data[i];
        state->genotypes[lc.id][lc.locus][lc.allele] = GENOTYPErandAllele(lc.locus);
    }
    PROBrecalcPosteriors(state);

    /* now, use Gibbs-Sampling */
    for (i = 0; i < Data._num_missing_data; i++) {
#ifdef HAVE_OPENMP
#pragma omp master
#endif
        if (Options.Verbosity > 0)
            VTPROGRESSBARupdate("Initial Sampling Missing Data",
                                total, *completed);

        lc = Data._missing_data[i];

        if (!last_both) {
            gibbsChangeGenotype(state, lc);
            last_both = lc.both;
        } else
            last_both = false;

        /* missing values are smaller 0 */
        assert(state->genotypes[lc.id][lc.locus][lc.allele] >= 0);
#ifdef HAVE_OPENMP
#pragma omp atomic
#endif
        (*completed)++;

    }
}

/* This stores the current pedigree if it is better than the the pedigree    *
 * with the highest likelihood found so far...                               */
static void
storeCurrentAsBest(MCMC * mcmc, MCMC_STATE * state, STATE_SNAPSHOT * best)
{
    if (state->p > best->p) {
        best->p = state->p;
        best->Nf = state->Nf;
        best->Nm = state->Nm;

        if (!DAGareEqual(best->D, state->D))
            DAGcopyE(state->D, best->D);
        if (!mcmc->genotypes_fixed)
            DATAIOcopyGenotypes(state->genotypes, best->genotypes);
    }
}

/* This tries to insert the posterior h of sample into the pedigree. Returns *
 * true on success and takes care of the logging for fast undo later (so     *
 * we don't need to recalculate the pedigree likelihood). Assumes that the   *
 * old parents of sample are already removed (removeH).                      */
static  bool
insertH(MCMC_STATE * state, int sample, int h)
{
    int     v, w, ret, i;

    Edge    e[2];

    assert(h >= 0);

    v = state->posteriors[sample][h].v;
    w = state->posteriors[sample][h].w;

#ifdef DEBUG
    if (Options.Verbosity >= 2)
        printf("\n\t...Trying (%i, %i) for sample %s", v, w,
               Data.id_mapping[sample].description);
#endif

    if (v >= 0 && w >= 0) {

        ret = PEDIGREEinsertParents(state->D, sample, v, w);
        if (ret) {
            if (Options.MaxDepth > 0) {
                DAGcalcDepthV(state->D);
                for (i = 0; i < Data.num_samples; i++)
                    if (DAGgetDepthV(state->D, i) > Options.MaxDepth) {
                        PEDIGREEremoveParents(state->D, sample, e);
                        return false;
                    }
            }
#ifdef DEBUG
            if (Options.Verbosity >= 2)
                printf("\n\t...success %s (%s, %s) P: %f",
                       Data.id_mapping[sample].description,
                       Data.id_mapping[v].description,
                       Data.id_mapping[w].description, POSTERIORSCORE(sample,
                                                                      h));
#endif
            assert(state->_addd == 0);
            state->_edgesaddd[state->_addd].w = sample;
            state->_edgesaddd[state->_addd + 1].w = sample;
            state->_edgesaddd[state->_addd].v = v;
            state->_edgesaddd[state->_addd + 1].v = w;
            state->_addd += 2;
        } else
            return false;
    } else if (v >= 0) {
        ret = DAGinsertE(state->D, EDGE(v, sample));
        if (ret) {
            if (Options.MaxDepth > 0) {
                DAGcalcDepthV(state->D);
                for (i = 0; i < Data.num_samples; i++)
                    if (DAGgetDepthV(state->D, i) > Options.MaxDepth) {
                        DAGremoveE(state->D, EDGE(v, sample));
                        return false;
                    }
            }
#ifdef DEBUG
            if (Options.Verbosity >= 2)
                printf("\n\t...success %s (%s, --) P: %f",
                       Data.id_mapping[sample].description,
                       Data.id_mapping[v].description, POSTERIORSCORE(sample,
                                                                      h));
#endif
            assert(state->_addd == 0);
            state->_edgesaddd[state->_addd].w = sample;
            state->_edgesaddd[state->_addd].v = v;
            state->_addd++;
        } else
            return false;
    }
#ifdef DEBUG
    if (Options.Verbosity >= 2)
        printf("\n\t...success %s (--, --) P: %f",
               Data.id_mapping[sample].description, POSTERIORSCORE(sample, h));
#endif
    return true;
}

/* This function calculates the likelihood of the pedigree D */
double
MCMCscoreDag(Dag D, MCMC_STATE * state)
{
    int     i, ip;
    double  p = 0.;

    for (i = 0; i < Data.num_samples; i++) {

        ip = incomingHypo(D, state, i);

        /* we did not find the posterior? this should not happen (either input
         * pedigree is wrong (which we should check) or something else is fucked 
         * up */
        if (ip < 0) {
            fprintf(stderr,
                    "Invalid posterior in MCMCscoreDag(): %s\n",
                    Data.id_mapping[i].description);
        } else {
            p += POSTERIORSCORE(i, ip);
            if (p <= MINUSINF)
                break;
        }
    }
    return p;
}

static double
MCMCscoreDagGT(Dag D, MCMC_STATE * state)
{
    int     i, ip;
    double  p = 0.;

    for (i = 0; i < Data.num_samples; i++) {
        ip = incomingHypo(D, state, i);
        p += state->posteriors[i][ip].l;
    }
    return p;
}

/* For performance reasons, we don't recalculate the likelihood of the whole
 * pedigree. We calculate the new likelihood by looking at the likelihood of a
 * change. To make sure that this works, we can use this assertion. This makes
 * everything really slow! So don't use in production code!                  */
#ifndef NDEBUG
static int
scoreDagAssertion(MCMC_STATE * state, double p)
{
    double  score;

    /* recalculates the score of the pedigree and compares it with p. It
     * prints an error msg and returns 0 for assert if score is wrong    */

    score = MCMCscoreDag(state->D, state);

    if (CMP_DBL(score, p))
        return 1;

    fprintf(stderr,
            "ASSERT SCOREDAG FAILED! Calculated Score: %f  True Score: %f x_f: %f x_m: %f\n",
            p, score, state->x_f, state->x_m);

    return 0;
}
#endif

/* This recalculates the posterior probabilities when N or a genotype was    *
 * changed.                                                                  */
static inline void
recalcPP(MCMC_STATE * state)
{
    PROBrecalcPosteriorProbabilities(state->posteriors, state->Nf, state->Nm);
    state->p = MCMCscoreDag(state->D, state);
}

/* This just removes all recently added arcs, adds the removed ones and      *
 * resets the score                                                          */
static inline void
undoStepPedigreeChange(MCMC_STATE * state, double p_old)
{
    int     i;

    /* it seems we have added some edges, so remove them */
    for (i = 0; i < state->_addd; i++) {
#ifdef DEBUG
        if (Options.Verbosity > 1)
            fprintf(stderr, "UNDO: removing edge %i->%i\n",
                    state->_edgesaddd[i].v, state->_edgesaddd[i].w);
#endif
        assert(DAGhasE(state->D, state->_edgesaddd[i]));
        DAGremoveE(state->D, state->_edgesaddd[i]);
    }
    state->_addd = 0;

    /* and if we have removed some, we should re-add them */
    for (i = 0; i < state->_rmvd; i++) {
#ifdef DEBUG
        if (Options.Verbosity > 1)
            fprintf(stderr, "UNDO: inserting edge %i->%i\n",
                    state->_edgesrmvd[i].v, state->_edgesrmvd[i].w);
#endif
        (void)DAGinsertE(state->D, state->_edgesrmvd[i]);
    }
    state->_rmvd = 0;
    /* and reset the score */
    state->p = p_old;

    assert(scoreDagAssertion(state, state->p));
}


static inline void
undoStepCandidateParentsChange(MCMC_STATE * state,
                               double old_x_f, double old_x_m,
                               double old_sr_f, double old_sr_m)
{

    state->sr_f = old_sr_f;
    state->sr_m = old_sr_m;
    state->x_f = old_x_f;
    state->x_m = old_x_m;
    calcNfm(state);
    recalcPP(state);
    assert(scoreDagAssertion(state, state->p));
}

/* This function "samples" the current pedigree, meaning it stores various of
 * its attributes (parentages, N, likelihood, selfing-rate). */
static void
samplePedigree(MCMC * mcmc, MCMC_STATE * state)
{
    int     i, posterior_id;
    LOCUS_COORD lc;
    //FILE *out;
    
    assert(CMP_DBL(state->c,1.));

    for (i = 0; i < Data.num_samples; i++) {
        posterior_id = PEDIGREEgetPostId(state->posteriors, state->D, i);
        if (posterior_id < 0 || posterior_id >= Probs.num_posteriors[i])
            FATALINT("Invalid posterior in updateObservedCounts()");

        state->posteriors[i][posterior_id].observed++;
    }
    mcmc->mh_Nf[Probs.mh_sampled_pedigrees] = state->Nf;
    mcmc->mh_Nm[Probs.mh_sampled_pedigrees] = state->Nm;
    mcmc->mh_p[Probs.mh_sampled_pedigrees] = state->p;
    PEDIGREEestmateSamplingRate(state->D, &state->sr_f, &state->sr_m);
    mcmc->mh_sr_f[Probs.mh_sampled_pedigrees] = state->sr_f;
    mcmc->mh_sr_m[Probs.mh_sampled_pedigrees] = state->sr_m;
    if (Options.Selfing) {
        mcmc->mh_selfing[Probs.mh_sampled_pedigrees] =
            //PEDIGREEcalcSelfingRate(state->D,mcmc->selfing_lods,mcmc->selfing_pp);
            PEDIGREEcalcSelfingRate(state->D);
    }

    if (Options.GibbsTypingError) {
        for (i=0; i<Data.num_loci;i++) 
            mcmc->mh_te[i][Probs.mh_sampled_pedigrees] = Data.TypingError[i];
    }
    if (Data.num_ramets != Data.num_samples)
        mcmc->mh_clonal[Probs.mh_sampled_pedigrees] =
            PEDIGREEcalcClonalRate(state->D);

    if (mcmc->log_missing_alleles) {
        for (i = 0; i < Data._num_missing_data; i++) {
            lc = Data._missing_data[i];
            mcmc->mh_missing_alleles[i][FREQ_ALLELE_ID
                                        (lc.locus,
                                         state->genotypes[lc.id][lc.locus]
                                         [lc.allele])]++;
        }
    }
    if (Data.num_populations > 1) {
        PEDIGREEpoplist(state, mcmc->mh_pop[Probs.mh_sampled_pedigrees]);
        /*FOPENA(out, "populations.csv");
        for (i=0; i<TLENGTH(Data.num_populations); i++) {
            if (i!=0) fprintf(out, ",");
            fprintf(out, "%i",matrix[i]);
        }
        fprintf(out, "\n");
        FCLOSE(out);*/
    }
    Probs.mh_sampled_pedigrees++;
}

/* the different acceptance functions */
static  bool
acceptStepRandom(double p_old, double p_new, double temperature, double qr)
{
    if (p_new <= MINUSINF)
        return false;
    return true;
}

static  bool
acceptStepSA(double p_old, double p_new, double temperature, double qr)
{
    double  r = RANDDBLONE;

    if (p_old < p_new)
        return true;

    if (r < exp((p_new - p_old + qr) / temperature))
        return true;

    return false;
}

static  bool
acceptStepMH(double p_old, double p_new, double temperature, double qr)
{
    double  r = RANDDBLONE;

    if (r < exp(p_new - p_old +qr)) {
        return true;
    }
    return false;
}

static  bool
acceptStepMCMCMCSWAP(double pi, double pj, double ti, double tj)
{
    double  res;

    res = exp( (pj - pi)/ti +  (pi - pj)/tj);
    //fprintf(stderr,"%f %f %f %f %f\n",pi,pj,ti,tj,res);

    if (RANDDBLONE < res) {
        return true;
    }
    return false;
}

static inline void
performStepCandidateParents(MCMC_STATE * state)
{
#define MAXDEVEXP 4. 
    state->x_f =
        state->sr_f +
        ((MAXDEVEXP - state->sr_f) / D_MCMCCPGRID * RANDINT(D_MCMCCPGRID));
    if (!Data.has_sex_data) {
        state->x_m = state->x_f;
        assert(state->sr_f == state->sr_m);
    } else
        state->x_m =
            state->sr_m +
            ((MAXDEVEXP - state->sr_m) / D_MCMCCPGRID * RANDINT(D_MCMCCPGRID));

    assert(state->x_f >= state->sr_f && state->x_f <= MAXDEVEXP);
    assert(state->x_m >= state->sr_m && state->x_m <= MAXDEVEXP);
    PEDIGREEestmateSamplingRate(state->D, &state->sr_f, &state->sr_m);
    calcNfm(state);
    recalcPP(state);
}

static void createGibbsErrorRates(MCMC_STATE *state)
{
    int i,j;
    double min, max, inc,r;
    for (i=0;i<Data.num_loci;i++) {
        min = MIN(0.1,Data._TypingErrorEst[i] - Data._TypingErrorEst[i]/2);
        max = MIN(0.25,Data._TypingErrorEst[i] + Data._TypingErrorEst[i]/2);
        inc = (max-min) / (double)(D_MCMCTEGRID-1); 

        r = min;
        for (j=0;j<D_MCMCTEGRID;j++) {
            state->_gibbs_error_rate[i][j] = r;
            r += inc;
        }
        fprintf(stderr, "%f min %f max %f inc %f %i %i\n", state->_gibbs_error_rate[i][D_MCMCTEGRID-1],min,max,inc,i,j);
        assert(CMP_DBL(state->_gibbs_error_rate[i][D_MCMCTEGRID-1],max));
    }
}

inline static double
sampleGibbs(double *p, int cnt, double max) {
    int i;
    double p_sum = 0.;
    for (i = 0; i < cnt; i++) {
        /* the likelihood ration Pedigree/Pedigree_max */
        p[i] -= max;
        p[i] = exp(p[i]);
        /* the cumulative ratio */
        p_sum += p[i];
        //      if (lc.locus == 3) if (lc.both) fprintf(stderr, "GIBBS %i ID: %i Locus: %i p: %f %f\n", i, lc.id, lc.locus, state->_p_gibbs_genotype[i], p_sum);
        p[i] = p_sum;

    }
    return RANDDBL(p_sum);
}

void
gibbsChangeTypingError(MCMC_STATE *state) {
    int i,locus;
    double rd, max = MINUSINF;
    locus = RANDINT((double)Data.num_loci); 
    for (i=0;i<D_MCMCTEGRID;i++) {
        Data.TypingError[i] = state->_gibbs_error_rate[locus][i];
        LODrecalcErrorConstants();
        PROBrecalcPosteriors(state);
        state->_p_gibbs_error[i] = MCMCscoreDagGT(state->D, state);
        if (state->_p_gibbs_error[i] > max) max = state->_p_gibbs_error[i];
    }
    rd = sampleGibbs(state->_p_gibbs_error, D_MCMCTEGRID, max);
        for (i = 0; i < D_MCMCTEGRID; i++) {
            if (rd < state->_p_gibbs_error[i]) {
                Data.TypingError[i] = state->_gibbs_error_rate[locus][i];
                LODrecalcErrorConstants();
                PROBrecalcPosteriors(state);
                state->p = MCMCscoreDag(state->D, state);
                return;
            }
        }
}

/* This samples an allele at the specified locus with a simple Gibbs-Sampler. */
void
gibbsChangeGenotype(MCMC_STATE * state, LOCUS_COORD lc)
{

    int     i, j, k, last_allele, cnt = 0;
    double  rd, max = MINUSINF;

    /* first, store the probability of pedigree for each allele */
    if (lc.both)
        for (i = 0; i < Data.num_alleles[lc.locus]; i++)
            for (j = i; j < Data.num_alleles[lc.locus]; j++) {
                /* update the likelihoods by setting the two new alleles */
                for (k = 0; k < D_PLOIDY; k++) {
                    lc.allele = k;
                    last_allele = state->genotypes[lc.id][lc.locus][lc.allele];
                    if (lc.allele == 0)
                        state->genotypes[lc.id][lc.locus][lc.allele] =
                            Data.alleles[lc.locus][i];
                    else
                        state->genotypes[lc.id][lc.locus][lc.allele] =
                            Data.alleles[lc.locus][j];
                    PROBupdateGenotype(state, lc, last_allele);
                }

                state->p = MCMCscoreDag(state->D, state);
                state->_p_gibbs_genotype[cnt] =
                    MCMCscoreDagGT(state->D, state);

                //        if (lc.locus == 3)    fprintf(stderr, "GIBBS %i ID: %i Locus: %i %i.%i p: %f\n", cnt, lc.id, lc.locus, Data.alleles[lc.locus][i],Data.alleles[lc.locus][j], state->p);

                if (state->_p_gibbs_genotype[cnt] > max)
                    max = state->_p_gibbs_genotype[cnt];
                cnt++;
    } else
        for (i = 0; i < Data.num_alleles[lc.locus]; i++) {
            last_allele = state->genotypes[lc.id][lc.locus][lc.allele];
            state->genotypes[lc.id][lc.locus][lc.allele] =
                Data.alleles[lc.locus][i];

            PROBupdateGenotype(state, lc, last_allele);
            state->p = MCMCscoreDag(state->D, state);
            state->_p_gibbs_genotype[cnt] = MCMCscoreDagGT(state->D, state);

            //fprintf(stderr, "GIBBS %i ID: %i Locus: %i last: %i, new: %i p: %f\n", i, lc.id, lc.locus, last_allele, Data.alleles[lc.locus][i], state->p);

            if (state->_p_gibbs_genotype[i] > max)
                max = state->_p_gibbs_genotype[i];
            cnt++;
        }

    
    rd= sampleGibbs(state->_p_gibbs_genotype, cnt, max);

    if (lc.both)
        for (i = 0; i < cnt; i++) {
            if (rd < state->_p_gibbs_genotype[i]) {
                cnt = 0;
                for (j = 0; j < Data.num_alleles[lc.locus]; j++)
                    for (k = j; k < Data.num_alleles[lc.locus]; k++) {
                        if (cnt == i) {
//          if (lc.locus == 3)                 fprintf(stderr, "SELECTED %i ID: %i Locus: %i %i.%i %f\n", cnt, lc.id, lc.locus, Data.alleles[lc.locus][j],Data.alleles[lc.locus][k], rd);

                            lc.allele = 0;
                            last_allele =
                                state->genotypes[lc.id][lc.locus][lc.allele];
                            state->genotypes[lc.id][lc.locus][lc.allele] =
                                Data.alleles[lc.locus][j];
                            PROBupdateGenotype(state, lc, last_allele);
                            lc.allele = 1;
                            last_allele =
                                state->genotypes[lc.id][lc.locus][lc.allele];
                            state->genotypes[lc.id][lc.locus][lc.allele] =
                                Data.alleles[lc.locus][k];
                            PROBupdateGenotype(state, lc, last_allele);
                            state->p = MCMCscoreDag(state->D, state);
                            return;
                        }
                        cnt++;
                    }
            }
    } else
        for (i = 0; i < cnt; i++) {
            if (rd < state->_p_gibbs_genotype[i]) {
                last_allele = state->genotypes[lc.id][lc.locus][lc.allele];
                state->genotypes[lc.id][lc.locus][lc.allele] =
                    Data.alleles[lc.locus][i];
                PROBupdateGenotype(state, lc, last_allele);
                state->p = MCMCscoreDag(state->D, state);
                return;
            }
        }

    /* we shouldn't reach this */
    assert(0);
}

/* This removes the (indegree) incoming arcs of sample i and, again, it     *
 * takes care of the logging for fast undo.                                 */
static void
removeH(MCMC_STATE * state, int id, int indegree, int *old_hypothesis)
{

    if (indegree > 0) {
        DAGremoveIncomingE(state->D, id, &state->_edgesrmvd[state->_rmvd]);
        if (indegree == 1)
            state->_edgesrmvd[1].v = -1;

        /* now look how the probability has changed */
        *old_hypothesis =
            PROBsearchPosteriors(state->posteriors, id,
                                 state->_edgesrmvd[0].v,
                                 state->_edgesrmvd[1].v);
        assert(*old_hypothesis >= 0);

        state->_delta_score = POSTERIORSCORE(id, *old_hypothesis);
    } else {
        *old_hypothesis = PROBsearchPosteriors(state->posteriors, id, -1, -1);
        //delta_score = 0.;
        state->_delta_score = POSTERIORSCORE(id, *old_hypothesis);
    }
#ifdef DEBUG
    if (Options.Verbosity >= 2)
        printf
            ("\t...successfully removed parents of %s (%i,%i) P: %f delta: %f\n",
             Data.id_mapping[id].description,
             state->posteriors[id][*old_hypothesis].v,
             state->posteriors[id][*old_hypothesis].w,
             state->posteriors[id][*old_hypothesis].p_opt,
             state->_delta_score);
#endif
    state->_rmvd = indegree;
}


/* this function randomly picks one individual, removes its incoming arcs and
 * inserts new incoming arcs (a new, different parent combination)
 * */
static inline bool
pedigreeChangeHypo(MCMC * mcmc, MCMC_STATE * state, double p_old, double *qr)
{
    int     i, r1, r2, indegree, old_hypothesis = -1;
    double r, dsum = 0., dsum2 = 0.;
    bool    ret;
    
    assert(scoreDagAssertion(state, state->p));

    /* remove parents of random vertex */
    r1 = mcmc->_to_change[RANDINT((double)mcmc->_num_to_change)];
    assert(r1 >= 0 && r1 < Data.num_samples && Probs.num_posteriors[r1] > 1);
    indegree = DAGindegree(state->D, r1);
#ifdef DEBUG
    if (Options.Verbosity >= 2)
        printf
            ("\tTrying to remove parents of %s... in thread %i (old score: %f) %i \n",
             Data.id_mapping[r1].description, state->_thread_id, p_old,
             Probs.num_posteriors[r1]);
#endif
    removeH(state, r1, indegree, &old_hypothesis);
    /* the -1 is correct, we don't want to sample the old parents again */
    r2 = RANDINT((double)(Probs.num_posteriors[r1] - 1));
    if (r2 >= old_hypothesis)
      r2++;

    dsum=0.;
    for (i=0; i<Probs.num_posteriors[r1];i++) {
        if (i == old_hypothesis)  {
            dsum2 = exp(POSTERIORSCORE(r1,i));
            continue;
        }
        dsum += exp(POSTERIORSCORE(r1,i));
    }
    dsum2 += dsum;
    r = RANDDBL(dsum);
    dsum=0.;
    for (i=0; i<Probs.num_posteriors[r1];i++) {
        if (i == old_hypothesis) continue;
        dsum += exp(POSTERIORSCORE(r1,i));
        if (dsum >= r) {
            r2 = i;
            break;
       }    
    }        
   
    assert(r2 >= 0 && r2 < Probs.num_posteriors[r1] && r2 != old_hypothesis);

    ret = insertH(state, r1, r2);
    /* is the parentage possible? (no directed cycle)                */
    if (ret) {
        state->_delta_score = POSTERIORSCORE(r1, r2) - state->_delta_score;
        state->p += state->_delta_score;
        assert(scoreDagAssertion(state, state->p));
    }
    /* q(x_j,x_i)/q(x_i,x_j) = q(x_i) / q(x_j) */
    *qr = ( POSTERIORSCORE(r1,old_hypothesis) - log(dsum2 - exp(POSTERIORSCORE(r1,r2))) ) - ( POSTERIORSCORE(r1,r2) - log(dsum2 - exp(POSTERIORSCORE(r1,old_hypothesis))) );
//    fprintf(stderr, "%f   %f %f %f %i %i %i\n", *qr, exp(POSTERIORSCORE(r1,old_hypothesis)), exp(POSTERIORSCORE(r1,r2)), dsum2, old_hypothesis,r2, Probs.num_posteriors[r1]);
    return ret;
}

/* This function randomly picks one arc in the pedigree and changes its
 * direction if possible.                       */
static  bool
pedigreeChangeSwap(MCMC_STATE * state, double p_old)
{
    int     i, num_edges, r, changed[2], changed_hypo[2];
    double  delta = 0.;
    Edge    swapped;            /* the swapped arc */

    /* pick a random arc */
    num_edges = DAGedges(state->_edgesall, state->D);
    assert(num_edges > 0);
    r = RANDINT(num_edges);

    /* store the old scores to calcuate the new score later */
    delta -=
        POSTERIORSCORE(state->_edgesall[r].v,
                       incomingHypo(state->D, state, state->_edgesall[r].v));
    delta -=
        POSTERIORSCORE(state->_edgesall[r].w,
                       incomingHypo(state->D, state, state->_edgesall[r].w));

    /* now remove this arc, in the case of selfing, remove both arcs */
    assert(DAGhasE(state->D, state->_edgesall[r]));
    DAGremoveE(state->D, state->_edgesall[r]);

    /* for the step undo, log this change */
    state->_edgesrmvd[state->_rmvd] = state->_edgesall[r];
    state->_rmvd++;

    if (DAGhasE(state->D, state->_edgesall[r])) {
        /* if we don't allow selfing, mother and father must be different */
        assert(Options.Selfing);
        DAGremoveE(state->D, state->_edgesall[r]);
        state->_edgesrmvd[state->_rmvd] = state->_edgesall[r];
        state->_rmvd++;
    }

    /* now try to insert this edge with swapped direction */
    swapped.v = state->_edgesall[r].w;
    swapped.w = state->_edgesall[r].v;
#ifdef DEBUG
    if (Options.Verbosity >= 2)
        printf
            ("\tTrying to swap edge (%i,%i)... in thread %i (old score: %f)\n",
             swapped.w, swapped.v, state->_thread_id, p_old);
#endif
    if (!DAGinsertE(state->D, swapped)) {
#ifdef DEBUG
        if (Options.Verbosity >= 2)
            printf
                ("\tFAILED to swap edge (%i,%i)  in thread %i: cycle (old score: %f)\n",
                 swapped.w, swapped.v, state->_thread_id, p_old);
#endif
        return false;
    }
    /* again, log for step undo */
    state->_edgesaddd[state->_addd] = swapped;
    state->_addd++;

    /* is the new pedigree valid? */
    changed[0] = swapped.v;
    changed[1] = swapped.w;

    for (i = 0; i < 2; i++) {
        changed_hypo[i] = incomingHypo(state->D, state, changed[i]);
        if (changed_hypo[i] < 0) {
#ifdef DEBUG
            if (Options.Verbosity >= 2)
                printf
                    ("\tFAILED to swap edge (%i,%i)  in thread %i: not valid (old score: %f)\n",
                     swapped.w, swapped.v, state->_thread_id, p_old);
#endif
            return false;
        }
    }

    /* update the score */
    delta += POSTERIORSCORE(changed[0], changed_hypo[0]);
    delta += POSTERIORSCORE(changed[1], changed_hypo[1]);

    state->p += delta;

    assert(scoreDagAssertion(state, state->p));
    return true;
}

/* This is the holy function that changes the pedigree.                     */
static  bool
performSampling(MCMC * mcmc, MCMC_STATE * state, bool logmh,
                bool(*accfp) (double, double, double,double))
{
    int     r1, r2;
    bool    ret;

    double  p_old = state->p, old_x_f =
        state->x_f, old_x_m = state->x_m, old_sr_f = state->sr_f, old_sr_m =
        state->sr_m, qr;

    /* make sure that current pedigree likelihood is correct */
    assert(scoreDagAssertion(state, state->p));
    assert(DAGisAcyclic(state->D));

    state->_addd = state->_rmvd = 0;
    state->_delta_score = 0.;

    r1 = state->iter % (mcmc->_num_to_change + 3);

    if (!mcmc->N_fixed && r1 == 0) {

        performStepCandidateParents(state);

        if (!accfp(p_old, state->p, state->c,0.)) {
            undoStepCandidateParentsChange(state, old_x_f, old_x_m,
                                           old_sr_f, old_sr_m);
            return false;
        } else {
            if (logmh)
                samplePedigree(mcmc, state);
            return true;
        }

    } else if (!mcmc->genotypes_fixed && r1 == 1) {

        gibbsChangeGenotype(state, Data._missing_data[RANDINT((double)
                                                              (Data.
                                                               _num_missing_data))]);

        if (logmh)
            samplePedigree(mcmc, state);
        return true;
    } else if (state->_thread_id == 0 && Options.GibbsTypingError && r1 == 2) {

        gibbsChangeTypingError(state);

        if (logmh)
            samplePedigree(mcmc, state);
        return true;
    } 
    
    qr = 0.;

    if (Data.has_missing_birth_data && DAGgetEdgeNumber(state->D) > 0) {
        r2 = RANDINT((double)
                     (Data.num_samples + DAGgetEdgeNumber(state->D)));

        if (r2 > Data.num_samples)
            ret = pedigreeChangeSwap(state, p_old);
        else
            ret = pedigreeChangeHypo(mcmc, state, p_old, &qr);
    } else
        ret = pedigreeChangeHypo(mcmc, state, p_old, &qr);

    /* if not, or acceptace function returned false, discard change    */
    if (!ret || !accfp(p_old, state->p, state->c, qr)) {
#ifdef DEBUG
        if (Options.Verbosity >= 2)
            printf("\n===> REJECTED thread_id %i (%f)\n",
                   state->_thread_id, state->p);
#endif
        undoStepPedigreeChange(state, p_old);
        return false;
    }
#ifdef DEBUG
    if (Options.Verbosity >= 2)
        printf("\n===> ACCEPTED thread_id %i (%f)\n",
               state->_thread_id, state->p);
#endif
    /* save current pedigree in Metropolis sampling                  */
    if (logmh)
        samplePedigree(mcmc, state);
    return true;
}

static void
calcInitialSATemperature(MCMC * mcmc, MCMC_STATE * state)
{
    int     i, m1 = 0, m2 = 0;
    double  d[D_SAADJUSTRUNS], accept;

    /* user specified init. temp? then we skip this here... */
    if (Options.SAc0 > 0.) {
        mcmc->c0 = Options.SAc0;
        return;
    }

    state->c = 0;
    state->p = MCMCscoreDag(state->D, state);

    /* generate a few random transitions */
    for (i = 0; i < D_SAADJUSTRUNS; i++) {
        state->iter = i;
        (void)performSampling(mcmc, state, 0, acceptStepRandom);
        d[i] = state->p;
    }
    accept = 0.;
    while (state->c < D_SAMAXINITTEMP) {
        if (Options.Verbosity > 0)
            VTPROGRESSBARupdateNP("SA Optimization");
        for (i = 1; i < D_SAADJUSTRUNS; i++) {
            if (acceptStepSA(d[i], d[i - 1], state->c, 0.))
                m1++;
            else
                m2++;
        }

        accept = m1 / (double)(m1 + m2);
        if (accept < Options.SAchi)
            state->c++;
        else
            break;
    };
    mcmc->c0 = state->c;
    if (isnan(mcmc->c0))
        mcmc->c0 = 2.;
}

static void
startSAChain(MCMC * mcmc, MCMC_STATE * state, unsigned int n, FILE * logfile)
{
    unsigned int j, temp_stage = 0, temp_stage_begin = 0, convergence_events =
        0, accepted = 0;
    double  acceptance_rate, last_c, sigma_c0 = 0., sigma_ci, stage_avg =
        0., old_stage_avg = 0., convergence;
    int     temp_stage_size = (int)(Options.SAbeta * mcmc->N);

    if (mcmc->N == 0) {
        fprintf(logfile, "Nothing to optimize.");
        return;
    }
#ifdef HAVE_OPENMP
    state->_thread_id = omp_get_thread_num();
#else
    state->_thread_id = 0;
#endif

    MAKE1DDOUBLE(state->all_p, temp_stage_size);
    state->p = MCMCscoreDag(state->D, state);

    state->c = mcmc->c0;

    state->converged = false;
    (void)time(&state->begin);

    for (j = 0; j < n; j++) {
        state->iter = j;

        if (performSampling(mcmc, state, false, acceptStepSA))
            accepted++;

        storeCurrentAsBest(mcmc, state, &state->best);
        assert(scoreDagAssertion(state, state->p));

        /*store current score so that we can calculate the process variance */
        state->all_p[state->iter - temp_stage_begin] = state->p;
        if (j < 2 || (j % (temp_stage_size - 1)) != 0)
            continue;
        temp_stage++;

        if (state->_thread_id == 0 && Options.Verbosity > 0)
            VTPROGRESSBARupdateNP("SA Optimization");

        old_stage_avg = stage_avg;
        UTILScalcVarianceArray(state->all_p, 0,
                               state->iter - temp_stage_begin, &sigma_ci,
                               &stage_avg);

        sigma_ci = sqrt(sigma_ci);
        acceptance_rate = accepted / (double)temp_stage_size;
#ifdef HAVE_OPENMP        
        flockfile(logfile);
#endif        
        fprintf(logfile,
                " | %10u | %10u | %10.3f | %10.3f | %8.3f | %8.3f |  %.2f | %7.3f | %7.3f | %6i |\n",
                temp_stage, temp_stage_begin,
                state->best.p, stage_avg, sigma_ci, state->c,
                acceptance_rate, state->x_f, state->x_m, state->_thread_id);
        FFLUSH(logfile);
#ifdef HAVE_OPENMP        
        funlockfile(logfile);
#endif        
        accepted = 0;
        temp_stage_begin = state->iter + 1;
        if (acceptance_rate < D_SAPRIORCHI && state->_p_chi_temp < 0)
            state->_p_chi_temp = state->c;

        /* calculate new temperature */
        if (state->c > EPSILON) {
            last_c = state->c;
            state->c =
                state->c / (1.0 +
                            (state->c *
                             log(1.0 +
                                 Options.SAdelta * (1 +
                                                    state->_thread_id))) /
                            (3.0 * sigma_ci));
            if (last_c - state->c < 0.01)
                state->c = last_c - 0.01;
            if (state->c < 0.001)
                state->c = 0.;
        } else {
            state->converged = true;
            break;
        }

        /* check convergence after second stage */
        if (temp_stage == 1) {
            /* store process variance of first stage for convergence */
            sigma_c0 = sigma_ci;
            if (CMP_DBL(sigma_c0, 0))
                sigma_c0 = EPSILON;

        } else if (temp_stage > 1) {
            /* calculate converge criteria according Almudevar, Theoretical
             * Population Biology 63 (2003) 63-75                           */

            convergence = (stage_avg - old_stage_avg) / sigma_c0;

            if (fabs(convergence) <= Options.SAepsilon) {
                convergence_events++;
                /* convergence inferred after s->N_epsilon convergence 
                 * events                                                   */
                if (convergence_events >= Options.SANepsilon) {
                    state->converged = true;
                    /* one greedy stage */
                    if (state->c < EPSILON)
                        break;
                    state->c = 0.;
                    /*break; */
                }
            } else {
                convergence_events = 0;
            }
        }
        /*
           if (Options.UpdateAlleleFreqs && state->_thread_id == 0) {
           FREQupdate(state->D);
           PROBrecalcPosteriors(state);
           state->p = MCMCscoreDag(state->D, state);
           mcmc->p_best =
           MCMCscoreDag(mcmc->D_best,state);
           } */
    }
    assert(scoreDagAssertion(state, state->p));
    FREE(state->all_p);
}

static void
dumpMHparams(MCMC * mcmc)
{
    unsigned int i, j, k, sum;
    char    c;
    double  fsum, p;
    FILE   *out;
    LOCUS_COORD lc;

    FOPENW(out, Options.MHparamsfile);
    if (Options.Selfing) {
        fprintf(out, "pedigree,p,Nf,Nm,samplingratef,samplingratem,selfing\n");

        for (i = 0; i < Probs.mh_sampled_pedigrees; i++) {
            fprintf(out, "%i,%f,%f,%f,%f,%f,%f", i,mcmc->mh_p[i],
                    mcmc->mh_Nf[i], mcmc->mh_Nm[i], mcmc->mh_sr_f[i],
                    mcmc->mh_sr_m[i], mcmc->mh_selfing[i]);
            if (Options.GibbsMissingData)
                for (j=0; j<Data.num_loci;j++)
                    fprintf(out, ",%f",Data.TypingError[i]);
            fprintf(out, "\n");
        }    
    } else {
        fprintf(out, "pedigree,p,Nf,Nm,samplingratef,samplingratem\n");
        for (i = 0; i < Probs.mh_sampled_pedigrees; i++)
            fprintf(out, "%i,%f,%f,%f,%f,%f\n", i,mcmc->mh_p[i], mcmc->mh_Nf[i],
                    mcmc->mh_Nm[i], mcmc->mh_sr_f[i], mcmc->mh_sr_m[i]);
    }
    FCLOSE(out);

    if (Data.num_populations > 1)
        PEDIGREEpopgraphList(mcmc);

    if (!mcmc->log_missing_alleles)
        return;

    FOPENW(out, Options.MissingAllelesOutfilename);
    DATAIOdumpDatasetDetails(out);
    for (i = 0; i < Data.num_loci; i++) {
        if (!Data.locus_has_missing_data[i])
            continue;

        FFLUSH(out);
        sum = 0;
        for (j = 0; j < Data._num_missing_data; j++) {
            lc = Data._missing_data[j];
            if (lc.locus != i)
                continue;
            sum++;
            if (sum == 1) {
                fprintf(out, "*** Locus %s ***\n\n", Data.loci_ids[i]);
                fprintf(out, "%-10s ", "Genotype");
                for (k = 0; k < FREQ_RANGE(i); k++) {
                    if (Data.allele_frequencies[i].freqs[k] <= 0.)
                        continue;
                    fprintf(out, "  %6u", Data.allele_frequencies[i].min + k);
                }
#ifdef DEBUG
                fprintf(out, "%6s", "Sum");
#endif
                fprintf(out, "\n");
            }
            fprintf(out, "%-10s ", Data.id_mapping[lc.id].description);
            fsum = 0.;
            for (k = 0; k < Data.num_alleles[i]; k++) {
                if (FREQ_ALLELE_ID
                    (i, mcmc->best.genotypes[lc.id][lc.locus][lc.allele]) == k)
                    c = '*';
                else
                    c = ' ';
                p = mcmc->mh_missing_alleles[j][k] /
                    (double)Probs.mh_sampled_pedigrees;
                fsum += p;
                fprintf(out, " %c%.4f", c, p);
            }
#ifdef DEBUG
            fprintf(out, "  %.4f", fsum);
#endif
            fprintf(out, "\n");
        }
        if (sum > 0)
            fprintf(out, "\nTOTAL: %u\n\n", sum);
    }
    FCLOSE(out);
}


static void
startSimulatedAnnealing(MCMC * mcmc, MCMC_STATE * states)
{
    int     i, tid;
    double *chainlengths;
    FILE   *out = NULL;

    FOPENW(out, Options.MCMClogfile);
    if (Options.SAchains > 0) {
        if (Options.Verbosity > 0)
            VTPROGRESSBARupdateNP("SA Optimization");
        calcInitialSATemperature(mcmc, &states[0]);
#ifdef HAVE_OPENMP
        if (omp_get_max_threads() > Options.SAchains)
            Options.SAchains = omp_get_max_threads();
#endif
    }

    MCMCdumpSettings(out, mcmc);
    fprintf(out,
            "\n +------------+------------+------------+------------+----------+----------+-------+---------+---------+--------+\n");
    fprintf(out,
            " | %10s | %10s | %10s | %10s | %8s | %8s | %5s | %7s | %7s | %6s |\n",
            "Stage", "Iteration", "Best Score", "Average", "Sigma",
            "Temp", "Ratio", "x_f", "x_m", "Thread");
    fprintf(out,
            " +------------+------------+------------+------------+----------+----------+-------+---------+---------+--------+\n");

    MAKE1DDOUBLE(chainlengths, Options.SAchains);
#ifdef HAVE_OPENMP
#pragma omp parallel for private (tid)
#endif
    for (i = 0; i < Options.SAchains; i++) {
#ifdef HAVE_OPENMP
        tid = omp_get_thread_num();
#else
        tid = 0;
#endif
        if (tid == 0 && Options.Verbosity > 0)
            VTPROGRESSBARupdateNP("SA Optimization");

        startSAChain(mcmc, &states[tid], Options.SAmaxIterations, out);
        chainlengths[i] = difftime(time(NULL), states[tid].begin);

        if (!states[tid].converged)
            WARN("Chain not converged, try again with relaxed parameters");
    }
    fprintf(out,
            " +------------+------------+------------+------------+----------+----------+-------+---------+---------+--------+\n\n");
    for (i = 0; i < Options.SAchains; i++)
        fprintf(out, "Time Chain %i: %f\n", i, chainlengths[i]);
    FCLOSE(out);
    FREE1D(chainlengths);
    if (Options.SAchains > 0 && Options.Verbosity > 0)
        VTPROGRESSBARcompleted("SA Optimization");

}

void
MCMCstart(MCMC * mcmc)
{
    int     i, tmp_posterior_id, mh_accepted = 0, swap_id1, swap_id2, nthreads;

    unsigned int gibbs_loci_completed = 0, j, jj, k, l, last_sampled =
        0, batches, burnin_batches;
    bool    logmh = false, recalcposteriors = false, tmp_bool;
    double e_loglik;
    MCMC_STATE *states, tmp_state;
    EXACT e;


#ifdef HAVE_OPENMP
    nthreads = MAX(Options.MHchains, omp_get_max_threads());
#else
    nthreads = MAX(Options.MHchains, 1);
#endif


    Probs.mh_sampled_pedigrees = 0;
    if (mcmc->_num_to_change == 0) return;

    MALLOC(states, MCMC_STATE, nthreads);
    states[0] = MCMCinitState(mcmc, NULL);
    (void)MCMCgreedy(mcmc, &states[0]);
    
/*    if (Options.Selfing)
        Selfing(mcmc);*/
    /* no MCMC necessary ? */
    if (!Options.UpdateAlleleFreqs && mcmc->posteriors_fixed && mcmc->no_cycles) {
        MCMCdestroyState(mcmc, states[0]);
        FREE(states);
        return;
    }
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for (i = 0; i < nthreads; i++) {
        /* already initialized */
        if (i != 0)
            states[i] = MCMCinitState(mcmc, states[0].D);
        fillMissingData(mcmc, &states[i], &gibbs_loci_completed,
                        (nthreads * Data._num_missing_data));
        states[i].p = MCMCscoreDag(states[i].D, &states[i]);
        assert(scoreDagAssertion(&states[i], states[i].p));
    }
    if (gibbs_loci_completed > 0 && Options.Verbosity > 0)
        VTPROGRESSBARcompleted("Initial Sampling Missing Data");

    if (Data.num_samples <= Options.SAexactMax && mcmc->posteriors_fixed) {
        if (Options.Verbosity > 1)
            fprintf(stderr, "Less than %d individuals, doing exhausitive pedigree enumeration.\n", Options.SAexactMax + 1);
        e = EXACTinit();
        EXACTgetBestSinks(e);
        EXACTfindSinkOrdering(e);
        e_loglik = EXACTfindPedigree(e);
        DAGcopyE(e.D, states[0].D);
        states[0].p = MCMCscoreDag(states[0].D, &states[0]);
        assert(CMP_DBL(states[0].p,e_loglik));
        EXACTdestroy(e);
    }
    else startSimulatedAnnealing(mcmc, states);

    if (Options.MHchains < 0)
        Options.MHchains = nthreads;
    if (Options.MHNswaps < 0)
        Options.MHNswaps = nthreads-1;

    for (i = 0; i < Options.MHchains; i++) {
        states[i].c = (1. + i * Options.MHTemp);
        states[i].p = MCMCscoreDag(states[i].D, &states[i]);
        storeCurrentAsBest(mcmc, &states[i], &mcmc->best);
    }

    /* we later only sample from the unheated chain */
    states[0].unheated = true;

    batches = (Options.MHBurnin + Options.MHIterations) / Options.MHswapfreq;
    burnin_batches = Options.MHBurnin / Options.MHswapfreq;
    
    if (Options.MHchains > 0) {
    for (l = 0; l < batches; l++) {
        if (Options.Verbosity > 0) {
            if (l < burnin_batches)
                VTPROGRESSBARupdate("MCMC (Burnin)", batches, l);
            else
                VTPROGRESSBARupdate("MCMC (Sampling)", batches, l);
        }
#ifdef HAVE_OPENMP
#pragma omp parallel for private(i,j,jj)
#endif
        for (i = 0; i < Options.MHchains; i++) {

            /* after allele frequency updating, we need to recalculate      *
             * everything :(                                                */
            if (recalcposteriors) {
                /* with post. fixed, we already did this */
                if (!mcmc->posteriors_fixed)
                    PROBrecalcPosteriors(&states[i]);
                states[i].p = MCMCscoreDag(states[i].D, &states[i]);
            }

            assert(scoreDagAssertion(&states[i], states[i].p));

            for (jj = 0; jj < Options.MHswapfreq; jj++) {
                j = jj + (l * Options.MHswapfreq);
                
                states[i].iter = j;

                if (j < D_MHRANDOMITER) {
                    /* start at a random point (accept all valid changes) */
                    (void)performSampling(mcmc, &states[i], false,
                                          acceptStepRandom);
                } else if (states[i].unheated) {

                    assert(CMP_DBL(states[i].c, 1.));

                    /* Only log the unheated chain after burnin. It might be *
                     * that the next step is invalid (cycle). In this case,  *
                     * logmh is still true and we try to log the next        *
                     * pedigree. If it was successful, then                  *
                     * mh_sampled_pedigrees is increased and logmh is set to * 
                     * false.                                                */
                    if (j > Options.MHBurnin
                        && (j % Options.MHsamplefreq) == 0) {
                        logmh = true;
                        last_sampled = Probs.mh_sampled_pedigrees;
                    } else if (last_sampled != Probs.mh_sampled_pedigrees)
                        logmh = false;

                    if (performSampling(mcmc, &states[i], logmh, acceptStepMH))
                        states[i].accepted++;

                } else {
                    /* the heated chains */
                    (void)performSampling(mcmc, &states[i], false,
                                          acceptStepSA);
                }

                storeCurrentAsBest(mcmc, &states[i], &states[i].best);
            }
        }

        recalcposteriors = false;

        if (Options.MHchains > 1)
            for (k = 0; k < Options.MHNswaps; k++) {
                mcmc->swaps++;
                /* pick a random pair of threads... */
                swap_id1 = RANDINT((double)(Options.MHchains));
                swap_id2 = swap_id1;
                while (swap_id1 == swap_id2) {
                    swap_id2 = RANDINT((double)(Options.MHchains));
                }

                /* now accept with the standard probability */
                if (acceptStepMCMCMCSWAP
                    (states[swap_id1].p, states[swap_id2].p,
                     states[swap_id1].c, states[swap_id2].c)) {
                    mcmc->swaps_accepted++;

                    /* swap is simply a swap of the temps */
                    SWAP(states[swap_id1].c, states[swap_id2].c, tmp_state.c);
                    SWAP(states[swap_id1].unheated, states[swap_id2].unheated,
                         tmp_bool);
                }
            }
        /* use the unheated pedigree to update allele freqs */
        if (Options.UpdateAlleleFreqs && l < (batches - 1) ) {
            for (k = 0; k < Options.MHchains; k++) {
                if (!states[k].unheated) continue;
                FREQupdate(states[k].D);
                states[k].freqs_updated++;
                recalcposteriors = true;
                if (mcmc->posteriors_fixed)
                    PROBrecalcPosteriors(&states[0]);
            }
        }
    }
    if (Options.Verbosity > 0)
        VTPROGRESSBARcompleted("MCMC (Sampling)");
    }
    assert(scoreDagAssertion(&states[0], states[0].p));

    /* now collect the information from the threads together and store them
     * in our global data structures */
    if (!mcmc->posteriors_fixed)
        for (i = 0; i < Data.num_samples; i++)
            for (j = 0; j < Probs.num_posteriors[i]; j++) {
                Probs.posteriors[i][j].observed = 0;
                for (k = 0; k < nthreads; k++) {
                    tmp_posterior_id =
                        PROBsearchPosteriors(states[k].posteriors, i,
                                             Probs.posteriors[i][j].v,
                                             Probs.posteriors[i][j].w);
                    assert(tmp_posterior_id >= 0);
                    Probs.posteriors[i][j].observed +=
                        states[k].posteriors[i][tmp_posterior_id].observed;
                }
            }

    /* write the results in mhparams.dat */
    dumpMHparams(mcmc);

    for (i = 0; i < nthreads; i++) {
        /* choose the best of all threads as the ML pedigree */
        storeCurrentAsBest(mcmc, &states[i], &mcmc->best);
        mh_accepted += states[i].accepted;
        mcmc->freqs_updated += states[i].freqs_updated;
    }

    mcmc->mh_acceptr = mh_accepted / (double)Options.MHIterations;

    PROBrecalcPosteriorProbabilities(Probs.posteriors,
                                     mcmc->best.Nf, mcmc->best.Nm);
    for (i = 0; i < nthreads; i++)
        MCMCdestroyState(mcmc, states[i]);
    FREE(states);
}

/* This function tries all valid parents of an individual. It starts        *
 * with the one with the highest likelihood. If this fails, it tries the    *
 * second most likeliest one and so on...                                   */
static  bool
foundAParentage(MCMC_STATE * state, int child_id)
{
    int     i;
    bool    found = false;

/* the possible parents for sample i are sorted according their
    * probability in Probs.posteriors */
    for (i = 0; i < Probs.num_posteriors[child_id]; i++) {
        state->_addd = 0;
        if (insertH(state, child_id, i)) {
            /* Was this possible? then we can go to the next individual ... 
             * if not, we have to try the next putative parents             */
            state->p += POSTERIORSCORE(child_id, i);
            found = true;
            break;
        }
    }
    return found;
}

/* An error message if there is no valid pedigree. Something is wrong with  *
 * the data. A false positive fullsib could also be the problem             */
static void
noGreedyPedigree(int id)
{
    fprintf(stderr,
            "There is no possible pedigree. Is the inpedigree correct?\n");
    fprintf(stderr,
            "Try --nofullsibtest. If this is working, examine %s and specify a\n",
            Options.SiblingsOutfilename[1]);
    fprintf(stderr,
            "p-Value cutoff with --fdrsibling. If you think it is a bug, report it.");
    fprintf(stderr, "\n\nCritical ID: %s\n", Data.id_mapping[id].description);
    exit(1);
}

/* This function constructs a minimal pedigree, that means it only adds the *
 * fixed arcs (known relationships)                                         */
static void
minPedigree(MCMC_STATE * state, bool * inserted)
{
    int     ip, i;

    /* first insert the parentages where we have no alternative             */
    for (i = 0; i < Data.num_samples; i++) {

        if (Probs.num_posteriors[i] > 1)
            continue;
        if (!foundAParentage(state, i))
            noGreedyPedigree(i);

        inserted[i] = true;
    }

    /* then insert all known relationships (where we do not have the no      *
     * parent alternative when all other arcs would introduce a directed     *
     * cycle                                                                 */
    if (Data.use_pedigree) {
        for (i = 0; i < Data.num_samples; i++) {
            if (inserted[i] || Data.in_pedigree[i][0] < 0)
                continue;
            ip = PROBsearchPosteriors(state->posteriors, i,
                                      Data.in_pedigree[i][0],
                                      Data.in_pedigree[i][1]);
            state->_addd = 0;
            if (!insertH(state, i, ip))
                noGreedyPedigree(i);
            state->p += POSTERIORSCORE(i, ip);
            inserted[i] = true;
        }
    }
}

bool
MCMCgreedy(MCMC * mcmc, MCMC_STATE * state)
{
    int     i;
    bool   *inserted;

    state->p = 0.;
    MAKE1DBOOL(inserted, Data.num_samples);
    if (DAGgetEdgeNumber(state->D) > 0)
        DAGremoveAllE(state->D);

    for (i = 0; i < Data.num_samples; i++)
        inserted[i] = false;

    if (!mcmc->no_cycles)
        minPedigree(state, inserted);

    for (i = 0; i < Data.num_samples; i++) {
        if (inserted[i])
            continue;
        if (!foundAParentage(state, i))
            noGreedyPedigree(i);
    }

    FREE1D(inserted);
    assert(scoreDagAssertion(state, state->p));
    storeCurrentAsBest(mcmc, state, &mcmc->best);
    return true;
}

void
MCMCdumpSettings(FILE * fp, MCMC * mcmc)
{
    DATAIOdumpDatasetDetails(fp);
    fprintf(fp, "*** MCMC Parameters ***\n\n");
    fprintf(fp, "Seed                         : %u\n", Options.Seed);
    fprintf(fp, "Update Allele Frequencies    : ");
    OPTIONSprintBool(fp, Options.UpdateAlleleFreqs);
    fprintf(fp, "Gibbs sample missing data    : ");
    OPTIONSprintBool(fp, Options.GibbsMissingData);
    fprintf(fp, "Check for Pedigree cycles    : ");
    OPTIONSprintBool(fp, !mcmc->no_cycles);
    fprintf(fp, "\n\n");
    fprintf(fp, "*** Simulated Annealing specific Parameters ***\n\n");
    fprintf(fp, "Initial acceptance rate chi  : %.3f\n", Options.SAchi);
    fprintf(fp, "Neighborhood Size Factor beta: %.3f\n", Options.SAbeta);
    fprintf(fp, "Increment delta              : %f\n", Options.SAdelta);
    fprintf(fp, "Converge event number        : %u\n", Options.SANepsilon);
    fprintf(fp, "Converge tolerance epsilon   : %.3f\n", Options.SAepsilon);
    fprintf(fp, "\n");
    fprintf(fp, "Calculated Parameters     \n");
    fprintf(fp, "  Neighboorhood Size         : %u\n", mcmc->N);
    fprintf(fp, "  Initial temperature c0     : %f\n", mcmc->c0);
    fprintf(fp, "\n\n");
    fprintf(fp, "*** Metropolis Hastings specific Parameters ***\n\n");
    fprintf(fp, "Burnin Iterations            : %u\n", Options.MHBurnin);
    fprintf(fp, "Iterations                   : %u\n", Options.MHIterations);
    fprintf(fp, "Sample frequency             : %u\n", Options.MHsamplefreq);
    fprintf(fp, "MCMCMC Swap frequency        : %u\n", Options.MHswapfreq);
    fprintf(fp, "MCMCMC Number of swaps       : %i (-1 == number threads)\n",
            Options.MHNswaps);
    fprintf(fp, "MCMCMC Temperature           : %.3f\n", Options.MHTemp);
    fprintf(fp, "\n\n");
    fprintf(fp, "*** Data Parameters ***\n\n");
    fprintf(fp, "  Nfmax                      : %i\n", Options.Nfmax);
    fprintf(fp, "  Nmmax                      : %i\n", Options.Nmmax);
}

/* This appends some statistics to the summary.txt file */
void
MCMCdumpResults(FILE * fp, MCMC * mcmc)
{
    int     i, max_depth, *depths;
    double  cp_avg, cp_var;

    fprintf(fp, "*** Maximum Likelihood Pedigree ***\n\n");
    if (mcmc->D_true != NULL)
        fprintf(fp, "Likelihood true pedigree    : %.3f\n", mcmc->p_true);
    if (mcmc->best.p > MINUSINF)
        fprintf(fp, "Likelihood                  : %.3f\n", mcmc->best.p);
    fprintf(fp, "Number of Individuals\n");
    max_depth = DAGcalcDepthV(mcmc->best.D);
    MAKE1DINT(depths, max_depth + 1);
    for (i = 0; i < Data.num_samples; i++)
        depths[DAGgetDepthV(mcmc->best.D, i)]++;
    for (i = 0; i < max_depth + 1; i++)
        fprintf(fp, "  Generation %4i           : %5i\n", i - max_depth,
                depths[i]);
    FREE1D(depths);

    fprintf(fp, "\n\n*** MCMC ***\n\n");
    if (Probs.mh_sampled_pedigrees < 1) {
        fprintf(fp, "No MCMC sampling\n");
    } else {
        fprintf(fp, "Observed Pedigrees in MH    : %u\n",
                Probs.mh_sampled_pedigrees);
        fprintf(fp, "Acceptance Ratio in MH      : %.3f\n", mcmc->mh_acceptr);
        fprintf(fp, "Allele Freq. updates in MH  : %u\n", mcmc->freqs_updated);
        /*UTILScalcVarianceArray(mcmc->mh_cp, 0, Probs.mh_sampled_pedigrees-1,
           &cp_var, &cp_avg);
           fprintf(fp, "Candidate Parents in MH     : %.3f (+- %.3f)\n", cp_avg,
           sqrt(cp_var)); */
        UTILScalcVarianceArray(mcmc->mh_sr_f, 0,
                               Probs.mh_sampled_pedigrees - 1, &cp_var,
                               &cp_avg);
        fprintf(fp, "Estimated Sampling Rate Fem.: %.3f (+- %.3f)\n", cp_avg,
                sqrt(cp_var));
        UTILScalcVarianceArray(mcmc->mh_sr_m, 0,
                               Probs.mh_sampled_pedigrees - 1, &cp_var,
                               &cp_avg);
        fprintf(fp, "Estimated Sampling Rate Male: %.3f (+- %.3f)\n", cp_avg,
                sqrt(cp_var));
        fprintf(fp, "Swap Attempts in MCMCMC     : %u\n", mcmc->swaps);
        fprintf(fp, "Accepted Swaps in MCMCMC    : %u\n",
                mcmc->swaps_accepted);
    }
}


void
MCMCdestroyState(MCMC * mcmc, MCMC_STATE state)
{
    int i;
    FREE1D(state._edgesrmvd);
    FREE1D(state._edgesaddd);
    FREE1D(state._edgesall);
    FREE1D(state._p_gibbs_genotype);
    if (!mcmc->posteriors_fixed)
        PROBdestroyPosteriors(state.posteriors);
    if (!mcmc->genotypes_fixed)
        DATAIOdestroyGenotypes(state.genotypes);
    if (Options.GibbsTypingError) {
        FREE2D(state._gibbs_error_rate,Data.num_loci,i);
        FREE1D(state._p_gibbs_genotype);
    }    
    destroyBest(mcmc, &state.best);
    DAGdestroy(state.D);
}

void
MCMCdestroy(MCMC mcmc)
{
    int     i;
    FREE1D(mcmc.mh_Nf);
    FREE1D(mcmc.mh_Nm);
    FREE1D(mcmc.mh_sr_f);
    FREE1D(mcmc.mh_sr_m);
    FREE1D(mcmc.mh_p);
    FREE1D(mcmc._to_change);
    if (Options.Selfing) {
        FREE1D(mcmc.mh_selfing);
    }
    if (Data.num_ramets != Data.num_samples)
        FREE1D(mcmc.mh_clonal);

    if (mcmc.D_true != NULL)
        DAGdestroy(mcmc.D_true);

    if (mcmc.genotypes_fixed)
        DATAIOdestroyGenotypes(mcmc.genotypes);
    
    destroyBest(&mcmc, &mcmc.best);

    if (mcmc.log_missing_alleles) {
        FREE2D(mcmc.mh_missing_alleles, Data._num_missing_data, i);
    }

    if (Options.GibbsTypingError) {
        FREE2D(mcmc.mh_te,Data.num_loci,i);
    }    
}

void
MCMCcheckML(MCMC * mcmc)
{
    int     i, **pedigreeincheck = NULL;
    MCMC_STATE state;

    if (strlen(Options.PedigreeCheckInfilename) == 0)
        return;

    MAKE2DINT(pedigreeincheck, Data.num_samples, 2, i);
    DATAIOpedigreein(Options.PedigreeCheckInfilename, pedigreeincheck);
    mcmc->D_true = DAGinit(Data.num_samples, 2, mcmc->no_cycles);
    for (i = 0; i < Data.num_samples; i++) {
        if (pedigreeincheck[i][0] >= 0)
            DAGinsertE(mcmc->D_true, EDGE(pedigreeincheck[i][0], i));
        if (pedigreeincheck[i][1] >= 0)
            DAGinsertE(mcmc->D_true, EDGE(pedigreeincheck[i][1], i));
    }
    state = MCMCinitState(mcmc, NULL);
    state.Nf = mcmc->best.Nf;
    state.Nm = mcmc->best.Nm;
    mcmc->p_true = MCMCscoreDag(mcmc->D_true, &state);
    printf("%2i %f  ,%f,%f\n", Options.NumLoci,
           1. - DAGcalcDistance(mcmc->best.D, mcmc->D_true), mcmc->best.p,
           mcmc->p_true);
    MCMCdestroyState(mcmc, state);
    FREE2D(pedigreeincheck, Data.num_samples, i);
}

/* 
  vim: ft=c sw=4 ts=4 expandtab
*/
