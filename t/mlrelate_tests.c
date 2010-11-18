/*
 * A few tests that check if we get the same results as ML-relate.
 *
 * ML-relate counts partially missing genotypes (e.g. 100.?) in the allele
 * freq. estimations, so we ignore such loci here. ML-relate also does not
 * incorporate typing errors, so we have to set this rate to 0. here.
 *
 * $Id: mlrelate_tests.c 1689 2009-06-20 07:23:37Z markus $
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
    int i,n,ids[3],t;
    double p;
    char* filename = FILENAME; 
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
    Options.Fullsibs = true;
    Options.IgnoreAge = true;
    Options.TypingErrorDefined = 0.;
    Options.NumLoci = 5;
    DATAIOinfile(filename);
    Options.SiblingsOutformat[0] = Options.SiblingsOutformat[2] = false;
    Options.SiblingsOutformat[1] = true;
    snprintf(Options.SiblingsOutfilename[1], PATHLEN, "%s.%s", "siblings",
             D_TEXTSUFFIX);

    FREQcalcAlleleFrequencies();
    FREQcalcSummaryStatistics();
    PROBinit();
    sim = SIMstart();
    
    PROBcalcTriplesAndDyads();
    PROBcalcPosteriors();

    /* parentage, no missing data */
    HASH_FIND_STR(Data.sample_ids, "1574", s);
    ids[0] = s->id;
    HASH_FIND_STR(Data.sample_ids, "1578", s);
    ids[1] = s->id;
    HASH_FIND_STR(Data.sample_ids, "376", s);
    ids[2] = s->id;
    t = TIDXL(ids[0],ids[1]);
    p = Data.fs_matrix[t].ibd.fs;
    ok(fabs(p + 24.67) < 0.01, "FS ibd correct");
    p = Data.fs_matrix[t].ibd.fs - Data.fs_matrix[t].ibd.u;
    ok(fabs(p - 9.66) < 0.01, "delta u ibd correct");
    p = Data.fs_matrix[t].ibd.fs - Data.fs_matrix[t].ibd.hs;
    ok(fabs(p - 5.16) < 0.01, "delta hs ibd correct");
    p = Data.fs_matrix[t].ibd.fs - Data.fs_matrix[t].ibd.po;
    ok(fabs(p - 3.19) < 0.01, "delta po ibd correct");


    t = TIDXL(ids[2],ids[1]);
    p = Data.fs_matrix[t].ibd.fs;
    ok(fabs(p + 28.15) < 0.01, "FS ibd correct");
    p = Data.fs_matrix[t].ibd.fs - Data.fs_matrix[t].ibd.u;
    ok(fabs(p - 6.66) < 0.01, "delta u ibd correct");
    p = Data.fs_matrix[t].ibd.fs - Data.fs_matrix[t].ibd.hs;
    ok(fabs(p - 2.61) < 0.01, "delta hs ibd correct");
    p = Data.fs_matrix[t].ibd.fs - Data.fs_matrix[t].ibd.po;
    ok(fabs(p - 0.48) < 0.01, "delta po ibd correct");

    SIMdestroy(sim);
    PROBdestroyPreMCMC();
    DATAIOdestroyPreMCMC();
    FREQdestroyPreMCMC();
    PROBdestroyPostMCMC();
    DATAIOdestroyPostMCMC();
    FREQdestroyPostMCMC();

    RANDDESTROY(i, n);
    remove(Options.LociOutfilename);
    remove(Options.SiblingsOutfilename[1]);
    remove(Options.MismatchOutfilename);

    return exit_status();       /* Return the correct exit code */
}

/* 
  vim: ft=c sw=4 ts=4 expandtab
*/
