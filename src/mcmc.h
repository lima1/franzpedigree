#ifndef MCMC_H
#define MCMC_H
/*
 * $Id: mcmc.h 2064 2010-05-26 11:20:32Z markus $
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

#include "dag.h"
#include "global.h"
#include "time.h"

typedef struct S_STATE_SNAPSHOT
{
    Dag  D;
    double p;
    double Nf;
    double Nm;
    int ***genotypes;
} STATE_SNAPSHOT;

typedef struct S_MCMC_STATE
{
    Dag D;                      /* the current pedigree                      */
    double p;                   /* current score of D                        */
    double *all_p;              /* all p to calculate the process variance   */
    unsigned int iter;          /* current MCMC/SA iteration                 */
    double c;                   /* current temperature                       */
    double x_f;                 /* number of candidateparents                */
    double x_m;                 /* number of candidateparents                */
    double sr_f;
    double sr_m;
    double Nf;
    double Nm;
    bool converged;             /* did the chain converge?                   */
    POSTERIORS **posteriors;    /* each thread has its own probabilities     */
    int ***genotypes;           /* the true genotypes                        */
    time_t begin;               /* start of the MCMC sampling for output     */
    unsigned int accepted;       /* number of accepted changes               */
    unsigned int freqs_updated;  /* number of allele frequency updates       */
    bool unheated;

    /* some temporary stuff we need for the calculations, only malloc'ed     *
     * once for performance reasons                                          */
    Edge *_edgesrmvd;           /* the added and removed edges of a iter...  */
    Edge *_edgesaddd;           /* ... so we can undo the changes easily)    */
    int _rmvd;                  /* the number of removed incoming edges      */
    int _addd;                  /* the number of added egdes                 */
    double _delta_score;        /* We have changed the pedigree, this is the *
                                   correspondig change of the score          */
    Edge *_edgesall;            /* array for all edges, allocate the max.
                                   number (vertices * 2) once                */
    double _p_chi_temp;
    int _thread_id;
    double *_p_gibbs_genotype;   /* temp array for genotype gibbs sampling   */
    double *_p_gibbs_error;
    double **_gibbs_error_rate;
    STATE_SNAPSHOT best;
} MCMC_STATE;


typedef struct S_MCMC
{

    Dag D_true;                 /* the true pedigree if known (testing)      */
    unsigned int N;             /* neighborhood of states                    */
    unsigned int N_posteriors;  /* neighborhood of posteriors (subset of N)  */
    double p_true;              /* score of the true pedigree                */
    double c0;                  /* starting temperature                      */
    double mh_acceptr;          /* the ratio accepted steps/step in MH       */
    double *mh_Nf;              /* candidateparents                          */
    double *mh_Nm;              /* candidateparents                          */
    double *mh_sr_f;            /* the estimated sampling rate               */
    double *mh_sr_m;            /* the estimated sampling rate               */
    double *mh_selfing;         /* the observed selfing rate                 */
    double *mh_clonal;          /* the observed rate of clonal reproduction  */
    double *mh_p;               /* the observed score                        */
    double **mh_te;             /* the observed typing error rate            */
    unsigned int **mh_pop;      /* the observed inter-pop parentages         */
    unsigned int swaps;         /* the number of swaps in MCMCMC             */ 
    unsigned int swaps_accepted;/* number of accepted MCMCMC swaps           */
    int _num_to_change;         /* number of non-founder genotypes           */
    int *_to_change;            /* ids of non-foundder genotypes             */
    unsigned int freqs_updated; /* total number of allele frequency updates  */
    bool genotypes_fixed;       /* do we change genotypes during sampling?   */
    int ***genotypes;           /* the true genotypes (used when we don't    *
                                 * sample genotypes, because then threads    *
                                 * don't need their own genotype array       */
    bool posteriors_fixed;      /* do we need posterior arrays for every     *
                                 * thread?                                   */
    bool N_fixed;               /* do we have to estimate N, Nf and/or Nm?   */ 
    bool no_cycles;             /* do we have to check for cycles?           */
    bool log_missing_alleles;
    int **mh_missing_alleles;
    STATE_SNAPSHOT best;
} MCMC;

MCMC MCMCinit();
MCMC_STATE MCMCinitState(MCMC *, Dag);
double MCMCscoreDag(Dag, MCMC_STATE *);
double MCMCscorePriors(Dag, int *);
void MCMCstart(MCMC *);
bool MCMCgreedy(MCMC *, MCMC_STATE *);
void MCMCdestroy(MCMC);
void MCMCdestroyState(MCMC *, MCMC_STATE);
void MCMCdumpSettings(FILE *, MCMC *);
void MCMCdumpResults(FILE *, MCMC *);
void MCMCcheckML(MCMC *);

/* 
  vim: ft=c sw=4 ts=4 expandtab
*/

#endif
