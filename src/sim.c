/*
 * $Id: sim.c 1885 2010-01-25 15:31:17Z markus $
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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "sim.h"

#include "macros.h"
#include "dataio.h"
#include "utils.h"
#include "global.h"
#include "ibd.h"
#include "lod.h"
#include "prob.h"
#include "genotype.h"
#include "vtprogressbar.h"

extern DATA Data;
extern PROBS Probs;
extern OPTIONS Options;

#define NUM_OFFSPRING 2

/* we generate N individuals in  initSimulation, this returns a random one */
#define RAND_CANDITATE sim->candidates[RANDINT((D_SIMNUMCANDIDATES))]

struct S_SIMULATION
{
    int ***candidates;        /* N randomly generated diploid genotypes      */

    /* mismatch result                                                       */
    int **all_mm;             /* the distributions of mismatches             */
    int max_mm_dyad;          /* the estimated max. number of mismatches for */
    int max_mm_triple;        /* dyads and triples                           */

    /* IDB result                                                            */
    double **delta_po;        /* delta po 0=for fullsibs, 1=parent-offspring */
    double **delta_hs;        /* delta hs 1=half-sib                         */
    double **delta_u;         /* delta hs 1=half-sib                         */
};

struct S_SIMULATION_STEP
{
    int **mother;             /* the genotype of the mother                  */
    int **father;             /* and father                                  */
    int ***offspring;         /* their offspring                             */
    int **halfsib;            /* a half-sib for delta-fs                     */
    bool selfing;
    IBD_LL ibd;               /* the calculated log-likelih. for deltapo/hs  */
};

typedef struct S_SIMULATION_STEP *SIMULATION_STEP;

static void
initSimulation(SIMULATION sim)
{
    int i, j;
    DATAIOestimateTypingError();
    MAKE3DINT(sim->candidates, D_SIMNUMCANDIDATES,
             Data.num_loci, 2, i, j);

    for (i = 0; i < D_SIMNUMCANDIDATES; i++)
        GENOTYPErand(sim->candidates[i]);

    MAKE2DINT(sim->all_mm, 4, Options.SimulationSteps * 2, i);
    MAKE2DDOUBLE(sim->delta_po, 2, Options.SimulationSteps, i);
    MAKE2DDOUBLE(sim->delta_hs, 2, Options.SimulationSteps, i);
    MAKE2DDOUBLE(sim->delta_u,  2, Options.SimulationSteps, i);
   
    if (Options.Selfing && Options.SimulationSelfingRate < 0.) 
        Options.SimulationSelfingRate = Data.SelfingRateAvg;
}


static SIMULATION_STEP
initSimulationStep(void)
{
    int i;
    int **affair;
    SIMULATION_STEP step = malloc(sizeof *step);
    if (step == NULL) FATAL("malloc failed");
    
    step->mother = GENOTYPEcreate(false);

    if (Options.Selfing && RANDDBLONE < Options.SimulationSelfingRate) {
        step->selfing = true;
        step->father = step->mother;
    }    
    else {
        step->selfing = false;
        step->father = GENOTYPEcreate(false);
    }    

    /* generate the offspring genotype(s) */
    MALLOC(step->offspring, int **, NUM_OFFSPRING);
    for (i = 0; i < NUM_OFFSPRING; i++) {
        step->offspring[i] = GENOTYPEcreate(true);
        if (step->selfing) 
            GENOTYPEselfing(step->offspring[i], step->mother);
        else
            GENOTYPEmate(step->offspring[i], step->mother, step->father);

        GENOTYPEmutate(step->offspring[i]);
    }

    /* create a second mating partner for half-sib */
    affair = GENOTYPEcreate(false);
    step->halfsib = GENOTYPEcreate(true);
    GENOTYPEmate(step->halfsib, step->mother, affair);

    GENOTYPEmutate(step->mother);
    GENOTYPEmutate(step->father);
    GENOTYPEmutate(step->halfsib);

    FREE2D(affair, Data.num_loci, i);

    return step;
}

static void
destroySimulationStep(SIMULATION_STEP step)
{
    int i, j;

    FREE2D(step->mother, Data.num_loci, i);
    if (!step->selfing) {
        FREE2D(step->father, Data.num_loci, i);
    }    
    FREE2D(step->halfsib, Data.num_loci, i);
    FREE3D(step->offspring, NUM_OFFSPRING, Data.num_loci, i, j);
    FREE(step);
}

static void
calcMismatchesChild(SIMULATION sim, SIMULATION_STEP step, int child_id,
                    int *mismatching, int *mismatching_triple,
                    int *mismatching_wrong, int *mismatching_wrong_triple)
{

    (void)LODcalcPchildDyad(step->offspring[child_id], step->mother,
                            mismatching, 1,NULL);
    (void)LODcalcPchildDyad(step->offspring[child_id], RAND_CANDITATE,
                            mismatching_wrong, 1, NULL);
    (void)LODcalcPchildTriple(step->offspring[child_id], step->mother,
                              step->father, mismatching_triple, 1, NULL);
    (void)LODcalcPchildTriple(step->offspring[child_id], step->mother,
                              RAND_CANDITATE, mismatching_wrong_triple, 1, NULL);
}

static void
performStep(SIMULATION sim, unsigned int i)
{
    SIMULATION_STEP step;

    step = initSimulationStep();

    /* how many mismatches can we observe? */
    calcMismatchesChild(sim, step, 0, &(sim->all_mm[0][i]),
                        &(sim->all_mm[1][i]), &(sim->all_mm[2][i]),
                        &(sim->all_mm[3][i]));
    calcMismatchesChild(sim, step, 1,
                        &(sim->all_mm[0][i + Options.SimulationSteps]),
                        &(sim->all_mm[1][i + Options.SimulationSteps]),
                        &(sim->all_mm[2][i + Options.SimulationSteps]),
                        &(sim->all_mm[3][i + Options.SimulationSteps])
        );
    
    if (Options.Fullsibs) {
        /* get the distribution of relationship probs of full-sibs */
        IBDcalcRelationshipLikelihoods(step->offspring[0], step->offspring[1],
                                    &step->ibd);

        sim->delta_po[0][i] = step->ibd.fs - step->ibd.po;
        sim->delta_hs[0][i] = step->ibd.fs - step->ibd.hs;
        sim->delta_u[0][i]  = step->ibd.fs - step->ibd.u;

        /* compare it to with the parent offspring distribution */
        IBDcalcRelationshipLikelihoods(step->offspring[0], step->mother,
                                    &step->ibd);
        sim->delta_po[1][i] = step->ibd.fs - step->ibd.po;

        IBDcalcRelationshipLikelihoods(step->offspring[0], step->halfsib,
                                    &step->ibd);
        sim->delta_hs[1][i] = step->ibd.fs - step->ibd.hs;

        IBDcalcRelationshipLikelihoods(step->offspring[0], RAND_CANDITATE,
                                    &step->ibd);

        sim->delta_u[1][i] = step->ibd.fs - step->ibd.u;
    }
    destroySimulationStep(step);
}

/* simulate sampling process */
SIMULATION
SIMstart(void)
{
    int i; /* old openmp compiler complain otherwise */
    unsigned int j, k, completed = 0;

    /* FILE *out; */
    SIMULATION sim = malloc(sizeof *sim);

    PROBinit();
    initSimulation(sim);

#ifdef HAVE_OPENMP
    #pragma omp parallel for private(i)
#endif
    for (i = 0; i < Options.SimulationSteps; i++) {
#ifdef HAVE_OPENMP
        if (omp_get_thread_num() == 0 && Options.Verbosity >= 1)
#else            
        if (Options.Verbosity >= 1) 
#endif 
            VTPROGRESSBARupdate("Simulation", Options.SimulationSteps,
                                completed);
        
        performStep(sim, i);
#ifdef HAVE_OPENMP
        #pragma omp atomic
#endif
        completed++;
    }
    if (Options.Verbosity > 0 ) 
        VTPROGRESSBARcompleted("Simulation");

    if (Options.Verbosity > 1 ) 
        fprintf(stderr, "\nSummarizing simulation results ...");

    if (!Options.FullsibAlternativesDefined && !Data.has_age_data) {
        Options.FullsibAlternatives[0] = 2;
        Options.FullsibAlternatives[1] = 4;
        Options.FullsibAlternatives[2] = 1;
    }

#ifdef HAVE_OPENMP
#pragma omp parallel sections
#endif
    {
#ifdef HAVE_OPENMP
#pragma omp section
#endif
        {
            if (Options.Fullsibs)
            UTILScmpDistributions(sim->delta_po[1], sim->delta_po[0],
                                  Options.SimulationSteps,
                                  Options.SimulationSteps, 
                                  &Probs.simresults.pv_lup_po, Options.FullsibAlternatives[0]);
        }
#ifdef HAVE_OPENMP
#pragma omp section
#endif
        {
            if (Options.Fullsibs)
            UTILScmpDistributions(sim->delta_hs[1], sim->delta_hs[0],
                                  Options.SimulationSteps,
                                  Options.SimulationSteps,
                                  &Probs.simresults.pv_lup_hs,
                                  Options.FullsibAlternatives[1]
                                  );
        }
#ifdef HAVE_OPENMP
#pragma omp section
#endif
        {
            if (Options.Fullsibs)
            UTILScmpDistributions(sim->delta_u[1], sim->delta_u[0],
                                  Options.SimulationSteps,
                                  Options.SimulationSteps, 
                                  &Probs.simresults.pv_lup_u,
                                  Options.FullsibAlternatives[2]
                                  );
        }
#ifdef HAVE_OPENMP
#pragma omp section
#endif
        {
        for (i = 0; i < 4; i++)
            qsort(sim->all_mm[i], Options.SimulationSteps * 2, sizeof(int),
                UTILScompare_ints_incr);
        }    
    }

    i = (int)(Options.SimulationSteps * 2. * 0.995);
    j = (unsigned int)(Options.SimulationSteps * 2. * 0.005);

    sim->max_mm_dyad = sim->all_mm[0][i];
    sim->max_mm_triple = sim->all_mm[1][i];

    j = (unsigned int)(Options.SimulationSteps * 2.);

    /* can we go higher in our significance level? */
    for (i = (int)(Options.SimulationSteps * 2. * 0.999); i< j; i++) {
        k = Options.SimulationSteps * 2 - i;
        if (sim->all_mm[0][i] < sim->all_mm[2][k])
            sim->max_mm_dyad = sim->all_mm[0][i];
        if (sim->all_mm[1][i] < sim->all_mm[3][k])
            sim->max_mm_triple = sim->all_mm[1][i];
    }
    
    if (sim->max_mm_dyad < sim->all_mm[0][Options.SimulationSteps * 2 - 1]) 
        sim->max_mm_dyad++;
    
    if (sim->max_mm_triple < sim->all_mm[1][Options.SimulationSteps * 2 - 1]) 
        sim->max_mm_triple++;

    /* this shouldn't happen and I've never observed it */
    if (sim->max_mm_dyad > sim->max_mm_triple)
        sim->max_mm_triple = sim->max_mm_dyad;
    if (!Options.MaxMismatchingDefined) {
        Options.MaxMismatchingDyad = sim->max_mm_dyad;
        Options.MaxMismatchingTriple = sim->max_mm_triple;
    }
    if (Options.Verbosity > 1) 
        fprintf(stderr, " done.\n\n");
    return sim;
}

void
SIMdump(FILE * fp, SIMULATION sim)
{
    fprintf(fp, "*** Simulation Results ***\n\n");
    fprintf(fp, "Max. mismatching loci Dyad  : %i\n",
            sim->max_mm_dyad);
    fprintf(fp, "Max. mismatching loci Triple: %i\n",
            sim->max_mm_triple);
    if (Options.Fullsibs) {
        fprintf(fp, "0.05 delta PO               : %6.3f\n",
                PVALUEfindCritical(Probs.simresults.pv_lup_po, 0.05));
        fprintf(fp, "0.05 delta HS               : %6.3f\n",
                PVALUEfindCritical(Probs.simresults.pv_lup_hs, 0.05));
        fprintf(fp, "0.05 delta U                : %6.3f\n",
                PVALUEfindCritical(Probs.simresults.pv_lup_u, 0.05));
    }
    fprintf(fp, "\n\n");
}

static void
dumpHistogram(FILE * fp, int *array, unsigned n)
{
    unsigned int i, last_i = 0;

    for (i = 0; i < n + 1; i++)
        if (i == n || array[last_i] != array[i]) {
            fprintf(fp, "%-28i: %u\n", array[last_i], (i - last_i));
            last_i = i;
        }
    fprintf(fp, "\n");
}

void
SIMdumpDetailed(FILE * fp, SIMULATION sim)
{
    int i;

    DATAIOdumpDatasetDetails(fp);

    fprintf(fp, "*** Simulation Settings ***\n\n");
    fprintf(fp, "Seed                        : %u\n", Options.Seed);
    fprintf(fp, "Simulation Steps            : %u\n",
            Options.SimulationSteps);
    if (Options.Selfing) 
        fprintf(fp, "Selfing Rate                : %.3f\n", Options.SimulationSelfingRate);
    fprintf(fp, "Rate of Typing Error\n");
    for (i = 0; i < Data.num_loci; i++)
        fprintf(fp, "  Locus %-10s          : %.3f\n", Data.loci_ids[i],
                Data.TypingError[i]);
    fprintf(fp, "Proportion Typed\n");
    for (i = 0; i < Data.num_loci; i++)
        fprintf(fp, "  Locus %-10s          : %.3f\n", Data.loci_ids[i],
                Options.ProportionTyped[i]);
    fprintf(fp, "\n\n");

    fprintf(fp, "*** Mismatching Loci ***\n\n");
    if (Options.MaxMismatchingDefined)
        fprintf(fp, "(Values not used because --maxmismatching defined!)\n");
    fprintf(fp, "Max. mismatching loci Dyad  : %i\n", sim->max_mm_dyad);
    fprintf(fp, "Max. mismatching loci Triple: %i\n", sim->max_mm_triple);
    fprintf(fp, "\n");
    fprintf(fp, "Histogram Mismatches Parent-Offspring pairs:\n");
    dumpHistogram(fp, sim->all_mm[0], (Options.SimulationSteps * 2));
    fprintf(fp, "Histogram Mismatches Unrelated pairs:\n");
    dumpHistogram(fp, sim->all_mm[2], (Options.SimulationSteps * 2));
    fprintf(fp, "Histogram Mismatches Offspring-Mother-Father triples:\n");
    dumpHistogram(fp, sim->all_mm[1], (Options.SimulationSteps * 2));
    fprintf(fp, "Histogram Mismatches Offspring-Mother-Unrelated triples:\n");
    dumpHistogram(fp, sim->all_mm[3], (Options.SimulationSteps * 2));
    fprintf(fp, "\n");
    if (Options.Fullsibs) {
        fprintf(fp, "*** Fullsib Delta Distributions ***\n\n");
        fprintf(fp, "Delta Parent-Offspring:\n\n");
        PVALUEdump(fp, Probs.simresults.pv_lup_po, 20);
        fprintf(fp, "Delta Half-sib:\n\n");
        PVALUEdump(fp, Probs.simresults.pv_lup_hs, 20);
        fprintf(fp, "Delta Unrelated:\n\n");
        PVALUEdump(fp, Probs.simresults.pv_lup_u, 20);
    }
}

void
SIMdestroy(SIMULATION sim)
{
    int i, j;

    FREE3D(sim->candidates, D_SIMNUMCANDIDATES,
            Data.num_loci, i, j);
    FREE2D(sim->all_mm, 4, i);
    FREE2D(sim->delta_po, 2, i);
    FREE2D(sim->delta_hs, 2, i);
    FREE2D(sim->delta_u, 2, i);
    FREE(sim);
}


/* 
  vim: ft=c sw=4 ts=4 expandtab
*/
