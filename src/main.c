/*
 * $Id: main.c 1885 2010-01-25 15:31:17Z markus $
 *
 *        _/_/_/_/  _/_/_/      _/_/    _/      _/            
 *      _/        _/    _/  _/    _/  _/_/    _/  _/_/_/_/   
 *     _/_/_/    _/_/_/    _/_/_/_/  _/  _/  _/      _/      
 *    _/        _/    _/  _/    _/  _/    _/_/    _/         
 *   _/        _/    _/  _/    _/  _/      _/  _/_/_/_/      
 *   
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
 * Markus Riester, markus@bioinf.uni-leipzig.de                  */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <signal.h>
#include <time.h>

#include "macros.h"             /* some important housekeeping macros   */
#include "dag.h"                /* the graph library                    */
#include "prob.h"               /* calculate all kinds of probabilities */
#include "mcmc.h"               /* the MCMC sampling                    */
#include "freq.h"               /* calculate allele frequencies         */
#include "sim.h"                /* sampling simulations                 */
#include "global.h"             /* constants and structs                */
#include "cmdline.h"            /* parse command line arguments         */
#include "dataio.h"             /* read and write the data              */
#include "options.h"            /* init/load/store options              */

/* global variables (see global.h)                               */
DATA    Data;
PROBS   Probs;
OPTIONS Options;

static void
sighandler(int sig)
{
    VTPROGRESSBARcursorVisible();
    fprintf(stderr, "\nAaargh\n");
    exit(EXIT_FAILURE);
}

int
main(int argc, char *argv[])
{
    unsigned int i, n;
    SIMULATION sim;
    MCMC    mcmc;
    GraphvizCluster cluster;    /* for the GRAPHVIZ output            */
    FILE   *freqsout = NULL, *out = NULL, *pout = NULL;
    time_t now;

    /* make the cursor visible again when user killed FRANz        */
    (void)signal(SIGABRT, &sighandler);
    (void)signal(SIGTERM, &sighandler);
    (void)signal(SIGINT, &sighandler);

#ifndef NDEBUG
    WARN("Reconfigure without --enable-assertions for optimal performance!");
#endif
    
    (void)time(&Data.StartTime);

    /* Load default and specified options and data                 */
    OPTIONSinit();
    CMDLINEparse(argc, argv);
    DATAIOinfile(Options.Infilename);

    /* load coordinates of sampling locations if available         */
    if (strlen(Options.Coordfilename) > 0)
        DATAIOgeofile(Options.Coordfilename, false);

    /* and/or use user provided distances of sampling locations    */
    if (strlen(Options.Geofilename) > 0)
        DATAIOgeofile(Options.Geofilename, true);

    /* use frequency file, otherwise calculate allele frequencies  */
    if (strlen(Options.FreqInfilename) > 0)
        DATAIOfreqfile(Options.FreqInfilename, 0);
    else {
        FREQcalcAlleleFrequencies();
        /* and then dump the frequencies in our own format         */
        FOPENW(freqsout, Options.FreqOutfilename);
        FREQdump(freqsout);
        FCLOSE(freqsout);
    }

    /* load prior pedigree (known mothers, sub-pedigrees)          */
    if (strlen(Options.PedigreeInfilename) > 0)
        DATAIOinitInpedigree();

    /* load known fullsibs if available                            */
    if (strlen(Options.FullsibInfilename) > 0)
        DATAIOfullsibFile(Options.FullsibInfilename);

    /* create n random number streams, one for each thread         */
#ifdef HAVE_OPENMP
    n = omp_get_max_threads();
#else
    n = 1;
#endif

    VTPROGRESSBARinit(PROGRESSBARSIZE);
    SRAND(Options.Seed, i, n);

    /* calculate the loci summary statistics                       */
    FREQcalcSummaryStatistics();

    /* simulate the sampling (mismatch distributions, IBD deltas)  */
    sim = SIMstart();

    FOPENW(out, Options.SimulationOutfilename);
    SIMdumpDetailed(out, sim);
    FCLOSE(out);

    /* Calculate probabilities of all parents-offspring triples    *
     * and dyads with the allowed number of mismatches             */
    PROBcalcTriplesAndDyads();

    /* This now generates an efficient data structure with the     *
     * likelihoods of all possible arcs in the pedigree. It also   *
     * filters very unlikely parentages.                           */
    PROBcalcPosteriors();

    /* Now we know all candidate parents, so we can create the 
     * CERVUS output                                               */
    if (strlen(Options.GenotypeOutfilenameCERVUS) > 0)
        DATAIOdumpGenotypesCERVUS();
    if (strlen(Options.OffspringOutfilenameCERVUS) > 0)
        DATAIOdumpOffspringCERVUS();
    if (strlen(Options.GenotypeOutfilenameGenepop) > 0)
        DATAIOdumpGenotypesGenepop();
    if (strlen(Options.GenotypeOutfilenameRMES) > 0)
        DATAIOdumpGenotypesRMES();
    if (strlen(Options.GenotypeOutfilenamePARENTE) > 0)
        DATAIOdumpGenotypesParente();
    
    /* Create the summary.txt file                                 */
    FOPENW(out, Options.Outfilename);
    DATAIOdumpDatasetDetails(out);
    FREQdumpSummaryStatistics(out);
    OPTIONSdump(out);
    DATAIOcheckData(out);
    if (Data.num_populations > 1)
        DATAIOdumpSamplingLocations(out);
    SIMdump(out, sim);

    /* Cleanup memory  */
    SIMdestroy(sim);
    PROBdestroyPreMCMC();
    DATAIOdestroyPreMCMC();
    FREQdestroyPreMCMC();

    mcmc = MCMCinit();

    if (!Options.NoReconstruction) {
        /* now start MCMC chains                                           */
        MCMCstart(&mcmc);

        if (Options.PedigreeOutformat[1]) {
            /* output DAG, cluster sampling locations in graphviz output   */
            cluster = DATAIOgraphvizCluster();
            /* display meta data in graphviz output                        */
            DATAIOaddGVNodeStyle(mcmc.best.D);
            DAGgraphviz(mcmc.best.D, Options.PedigreeOutfilename[1], cluster, true);
            DATAIOgraphvizClusterDestroy(cluster);
        }
        /* --pedigreeincheck compares our results with the true pedigree   */
        MCMCcheckML(&mcmc);
        MCMCdumpResults(out, &mcmc);
    }
    /* write parentage files                                       */
    for (i = 0; i < 2; i++)
        if (Options.POutformat[i]) {
            FOPENW(pout, Options.POutfilename[i]);
            PROBdumpPosteriors(pout, Probs.posteriors, mcmc.best.D,
                               mcmc.D_true, i);
            FCLOSE(pout);
        }

    /* output pedigree in our input format                         */
    if (Options.PedigreeOutformat[0])
        DATAIOpedigreeout(Options.PedigreeOutfilename[0], mcmc.best.D);
    if (Options.PedigreeOutformat[2])
        DATAIOpedigreeoutText(Options.PedigreeOutfilename[2], mcmc.best.D);
    
    (void)time(&now);
    fprintf(out, "\nRunning Time                : %d sec.\n", (unsigned int)( now - Data.StartTime));
    FCLOSE(out);
    DATAIOdumpMismatches();
    MCMCdestroy(mcmc);
    PROBdestroyPostMCMC();
    DATAIOdestroyPostMCMC();
    FREQdestroyPostMCMC();
    RANDDESTROY(i, n);

    VTPROGRESSBARcursorVisible();
    exit(EXIT_SUCCESS);
}

/* 
  vim: ft=c sw=4 ts=4 expandtab
*/
