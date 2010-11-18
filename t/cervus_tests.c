/*
 * A few tests that check if we get the same results as CERVUS 3.0.3
 *
 * $Id: cervus_tests.c 1672 2009-06-18 15:16:21Z markus $
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
#include "genotype.h"
#include "pedigree.h"
#include "dag.h"
#include "mcmc.h"
#include "sim.h"

DATA Data;
PROBS Probs;
OPTIONS Options;

int
main(int argc, char *argv[])
{
    int i,n,ids[3];
    char* filename = FILENAME; 
    char* filenamePedigree = FILENAMEPED; 
    DESCRIPTION_HASH *s;

    SIMULATION sim;       /* the sampling simulation data          */

    plan_tests(8);
    OPTIONSinit();
#ifdef HAVE_OPENMP
    n = omp_get_max_threads();
#else
    n=1;
#endif
    
    Options.Verbosity = 0;
    SRAND(123,i, n);
    Options.TypingErrorDefined = 0.01;
    Options.FemRepro.min = Options.FemRepro.max = 1;
    Options.MaleRepro.min = Options.MaleRepro.max = 1;
    DATAIOinfile(filename);
    FREQcalcAlleleFrequencies();
    FREQcalcSummaryStatistics();
    PROBinit();
    sim = SIMstart();
    
    PROBcalcTriplesAndDyads();
    PROBcalcPosteriors();

    /* parentage, no missing data */
    HASH_FIND_STR(Data.sample_ids, "1570", s);
    ids[0] = s->id;
    HASH_FIND_STR(Data.sample_ids, "221", s);
    ids[1] = s->id;
    HASH_FIND_STR(Data.sample_ids, "376", s);
    ids[2] = s->id;

    ok(CMP_DBL(Probs.posteriors[ids[0]][PROBsearchPosteriors(Probs.posteriors,ids[0],ids[1],ids[2])].lod, 1.98843824612410E+0001), "Triple LOD correct");
    ok(CMP_DBL(Probs.posteriors[ids[0]][PROBsearchPosteriors(Probs.posteriors,ids[0],ids[1],-1)].lod, 1.08164969253876E+0001), "Pair LOD correct");
    ok(CMP_DBL(Probs.posteriors[ids[0]][PROBsearchPosteriors(Probs.posteriors,ids[0],ids[2],-1)].lod, 6.94685569833075E+0000), "Pair LOD correct");

    /* parentage, missing data */
    HASH_FIND_STR(Data.sample_ids, "1810", s);
    ids[0] = s->id;
    HASH_FIND_STR(Data.sample_ids, "313", s);
    ids[1] = s->id;
    HASH_FIND_STR(Data.sample_ids, "222", s);
    ids[2] = s->id;

    ok(CMP_DBL(Probs.posteriors[ids[0]][PROBsearchPosteriors(Probs.posteriors,ids[0],ids[1],ids[2])].lod, 2.04404149122917E+0001), "Triple LOD correct (missing values)");
    ok(CMP_DBL(Probs.posteriors[ids[0]][PROBsearchPosteriors(Probs.posteriors,ids[0],ids[1],-1)].lod, 8.61834142652084E+0000 ), "Pair LOD correct (missing values)");
    ok(CMP_DBL(Probs.posteriors[ids[0]][PROBsearchPosteriors(Probs.posteriors,ids[0],ids[2],-1)].lod, 8.24729736134361E+0000), "Pair LOD correct (missing values)");

    SIMdestroy(sim);
    PROBdestroyPreMCMC();
    DATAIOdestroyPreMCMC();
    FREQdestroyPreMCMC();
    PROBdestroyPostMCMC();
    DATAIOdestroyPostMCMC();
    FREQdestroyPostMCMC();

    DATAIOinfile(filename);
    strncpy(Options.PedigreeInfilename,filenamePedigree,200);
    DATAIOinitInpedigree();
    FREQcalcAlleleFrequencies();
    FREQcalcSummaryStatistics();
    PROBinit();
    sim = SIMstart();
    
    PROBcalcTriplesAndDyads();
    PROBcalcPosteriors();
//    PROBdumpPosteriors(stderr,Probs.posteriors, NULL,NULL,2);

    HASH_FIND_STR(Data.sample_ids, "1570", s);
    ids[0] = s->id;
    HASH_FIND_STR(Data.sample_ids, "221", s);
    ids[1] = s->id;
    HASH_FIND_STR(Data.sample_ids, "376", s);
    ids[2] = s->id;
    /* paternity, no missing */
    ok(CMP_DBL(Probs.posteriors[ids[0]][PROBsearchPosteriors(Probs.posteriors,ids[0],ids[1],ids[2])].lod, 1.29375267629102E+0001), "Triple LOD correct (mother known)");

    HASH_FIND_STR(Data.sample_ids, "1810", s);
    ids[0] = s->id;
    HASH_FIND_STR(Data.sample_ids, "313", s);
    ids[1] = s->id;
    HASH_FIND_STR(Data.sample_ids, "222", s);
    ids[2] = s->id;

    ok(CMP_DBL(Probs.posteriors[ids[0]][PROBsearchPosteriors(Probs.posteriors,ids[0],ids[1],ids[2])].lod, 1.18220734857708E+0001), "Triple LOD correct (mother known, missing values)");

    RANDDESTROY(i, n);
    remove(Options.LociOutfilename);

    return exit_status();       /* Return the correct exit code */
}

/* 
  vim: ft=c sw=4 ts=4 expandtab
*/
