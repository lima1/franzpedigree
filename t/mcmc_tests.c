/*
 * $Id: mcmc_tests.c 1698 2009-06-22 13:22:32Z markus $
 */

#include "macros.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tap.h"

#include "global.h"
#include "dataio.h"
#include "prob.h"
#include "ibd.h"
#include "lod.h"
#include "freq.h"
#include "options.h"
#include "dag.h"
#include "mcmc.h"

DATA Data;
PROBS Probs;
OPTIONS Options;

int
main(int argc, char *argv[])
{
    int i,j,n;
    long observed;
    char* filename = FILENAMEMD;
    double sumd;
    MCMC mcmc;

    plan_tests(20);
#ifdef HAVE_OPENMP
    n = omp_get_max_threads();
#else
    n=1;
#endif
    
    SRAND(Options.Seed,i, n);
    OPTIONSinit();
    Options.Verbosity = 0;
    
    Options.FemRepro.min = Options.MaleRepro.min = 14;
    Options.FemRepro.max = Options.MaleRepro.max = 80;
    Options.MinTyped = 2;
    Options.MaxMismatchingDefined = 1;
    Options.MaxMismatchingTriple = 0;
    Options.MaxMismatchingDyad = 0;
    Options.GibbsMissingData = true;
    DATAIOinfile(filename);
    FREQcalcAlleleFrequencies();
    FREQcalcSummaryStatistics();
    PROBinit();
    PROBcalcTriplesAndDyads();
    PROBcalcPosteriors();
    /*PROBdumpPosteriors(stderr, Probs.posteriors, NULL, NULL,1);*/
    ok(Probs.are_connected[TIDXL(0,1)] == true,  "grampa homer");
    ok(Probs.are_connected[TIDXL(0,2)] == true,  "grampa bart");
    ok(Probs.are_connected[TIDXL(0,3)] == false, "grampa not involved");
    ok(Probs.are_connected[TIDXL(0,4)] == false, "grampa not involved");
    ok(Probs.are_connected[TIDXL(0,5)] == false, "grampa not involved");
    ok(Probs.are_connected[TIDXL(0,6)] == false, "grampa not involved");

    PROBdestroyPreMCMC();
    DATAIOdestroyPreMCMC();
    FREQdestroyPreMCMC();
    Options.MHIterations = 10000;
    mcmc = MCMCinit();
    MCMCstart(&mcmc);

    MCMCdestroy(mcmc);
    
    for(i=0; i<Data.num_samples; i++) {
        sumd =  0.;
        observed = 0;
        for (j=0; j<Probs.num_posteriors[i]; j++) {
            sumd += exp(Probs.posteriors[i][j].p_opt);
            observed += Probs.posteriors[i][j].observed;
        }    
        ok(sumd > 0.7 && sumd <= 1.0, "posteriors after mcmc sum to 1");
        ok(observed == Probs.mh_sampled_pedigrees, "sum sampled pedigrees correct");
    }
    PROBdestroyPostMCMC();
    DATAIOdestroyPostMCMC();
    FREQdestroyPostMCMC();
    remove(Options.SiblingsOutfilename[1]);
    remove(Options.LociOutfilename);
    remove(Options.MHparamsfile);
    remove(Options.MCMClogfile);
    return exit_status();       /* Return the correct exit code */
}

/* 
  vim: ft=c sw=4 ts=4 expandtab
*/
