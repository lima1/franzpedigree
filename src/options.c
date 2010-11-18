/*
 * $Id: options.c 2065 2010-05-26 11:21:52Z markus $
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
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "options.h"

extern DATA Data;
extern OPTIONS Options;

void
OPTIONSinit(void)
{
    strcpy(Options.Infilename, "");
    strcpy(Options.Geofilename, "");
    strcpy(Options.Coordfilename, "");
    strcpy(Options.HWETestOutfilename, "");
    strcpy(Options.FullsibInfilename, "");
    Options.FullsibInfile = false; 
    strcpy(Options.FreqInfilename, "");
    strcpy(Options.PedigreeInfilename, "");
    strcpy(Options.PedigreeCheckInfilename, "");
    strcpy(Options.MissingAllelesOutfilename, "");
    strcpy(Options.MismatchOutfilename, D_MISMATCHOUT);
    strcpy(Options.FreqOutfilename, D_FREQOUT);
    strcpy(Options.LociOutfilename, D_LOCIOUT);
    strcpy(Options.Outfilename, D_OUT);
    strcpy(Options.MCMClogfile, D_MCMCLOGFILE);
    strcpy(Options.MHparamsfile, D_MHPARAMSOUT);
    strcpy(Options.HaplotypeOutfilename, D_HAPLOOUT);
    strcpy(Options.GenotypeOutfilename, D_GENOTYPEOUT);
    strcpy(Options.GenotypeOutfilenameCERVUS, "");
    strcpy(Options.OffspringOutfilenameCERVUS, "");
    strcpy(Options.GenotypeOutfilenameGenepop, "");
    strcpy(Options.GenotypeOutfilenameRMES, "");
    strcpy(Options.GenotypeOutfilenamePARENTE, "");
    strcpy(Options.SimulationOutfilename, D_SIMULATIONOUT);
    /*strcpy(Options.Paramfile, "paramfile"); */
    Options.Verbosity = D_VERBOSITY;
    Options.Selfing = D_SELFING;
    Options.MaxDepth = D_MAXDEPTH;
    Options.FemRepro.min = D_FEMREPROMIN;
    Options.FemRepro.max = D_FEMREPROMAX;
    Options.MaleRepro.min = D_MALEREPROMIN;
    Options.MaleRepro.max = D_MALEREPROMAX;
    Options.FemMaxDist = -1;
    Options.MaleMaxDist = -1;
    Options.ParentsMaxDist = -1;
    Options.SibMaxDist = -1;
    Options.MinTyped = D_MINTYPED;
    Options.UpdateAlleleFreqs = D_SAUPDATEFREQS;
    Options.MaxMismatchingDefined = false;
    Options.MaxMismatchingTriple = D_MAXMISMATCHING;
    Options.MaxMismatchingDyad = D_MAXMISMATCHING;
    Options.TypingErrorDefined = -1;
    Options.ProportionTypedDefined = -1.;
    Options.SimulationSteps = D_SIMNUM;
    Options.Nf = 1;
    Options.Nm = 1;
    Options.Nr = -1;
    Options.n = -1;
    Options.ndefined = false;
    Options.Nfdefined = false;
    Options.Nmdefined = false;
    Options.Nfmax = -1;
    Options.Nmmax = -1;
    Options.SimulationSelfingRate = -1;
    Options.SAchains = D_SACHAINS;
    Options.SAmaxIterations = D_SAMAXITERATION;
    Options.SAchi = D_SACHI;
    Options.SAbeta = D_SABETA;
    Options.SAdelta = D_SADELTA;
    Options.SAepsilon = D_SAEPSILON;
    Options.SANepsilon = D_SANEPSILON;
    Options.SAc0 = -1;
    Options.SAexactMax = D_EXACTMAXINDIVIDUALS;
    Options.MHchains = -1;
    Options.MHBurnin = D_MHBURNINITER;
    Options.MHIterations = D_MHITERATION;
    Options.MHsamplefreq = D_MHSAMPLEFREQ;
    Options.MHswapfreq = D_MHSWAPFREQ;
    Options.MHTemp = D_MHTEMP;
    Options.NumLoci = D_NUMLOCI;
    Options.MHNswaps = D_MHNSWAPS;
    Options.Seed = (unsigned int)time(NULL);
    Options.IgnoreAge = false;
    Options.IgnoreSex = false;
    Options.IgnoreRamets = false;

    Options.Fullsibs = false;
    Options.FullsibAlternativesDefined = false;
    Options.FullsibAlternatives[0] = 1;
    Options.FullsibAlternatives[1] = 2;
    Options.FullsibAlternatives[2] = 1;
    Options.SiblingsCorrectionMethod = D_PV_CORRECTION_METHOD;
    Options.THsiblings[0] = D_TH_SIBLING_0;
    Options.THsiblings[1] = D_TH_SIBLING_1;
    Options.THsiblings[2] = D_TH_SIBLING_2;

    Options.HWESteps = D_HWESTEPS;
    Options.HWEChunks = D_HWECHUNKS;
    Options.HWEChunkSize = D_HWECHUNKSIZE;
    Options.NoReconstruction = false;
    Options.RametAlleleFreqs = false;
    Options.GibbsMissingData = false;
    Options.GibbsTypingError = false;
    Options.CollapseRamets = false;
    Options.SiblingsInParentalGeneration = D_SIBLING_IN_PARENTAL_GENERATION;

    Options.MonoeciousDefined = false;
    Options.Monoecious = false;

    Options.AssumeCompleteSample = false;
}

void
OPTIONSprintBool(FILE *fp, bool b) {
    if (b)
        fprintf(fp, "Yes\n");
    else
        fprintf(fp, "No\n");
}

void
OPTIONSdump(FILE * fp)
{
    int i;

    fprintf(fp, "*** Files ***\n\n");
    fprintf(fp, "Input\n");
    fprintf(fp, "  Genotype                  : %s\n", Options.Infilename);
    if (strlen(Options.FreqInfilename) > 0)
        fprintf(fp, "  Allele frequency file     : %s\n",
                Options.FreqInfilename);
    else
        fprintf(fp,
                "  Allele frequency file     : none (calculates frequencies)\n");
    if (strlen(Options.Geofilename) > 0)
        fprintf(fp, "  Geodist file              : %s\n",
                Options.Geofilename);
    else
        fprintf(fp, "  Geodist file              : none\n");
    if (strlen(Options.Coordfilename) > 0)
        fprintf(fp, "  Coordinate file           : %s\n",
                Options.Coordfilename);
    else
        fprintf(fp, "  Coordinate file           : none\n");

    if (strlen(Options.PedigreeInfilename) > 0)
        fprintf(fp, "  Pedigree file             : %s\n",
                Options.PedigreeInfilename);
    else
        fprintf(fp, "  Pedigree file             : none\n");
    if (strlen(Options.PedigreeCheckInfilename) > 0)
        fprintf(fp, "  Pedigree check file       : %s\n",
                Options.PedigreeCheckInfilename);

    fprintf(fp, "\n");
    fprintf(fp, "Output\n");
    fprintf(fp, "  Summary                   : %s\n", Options.Outfilename);
    fprintf(fp, "  Loci summary file         : %s\n",
            Options.LociOutfilename);
    fprintf(fp, "  Allele frequency file     : %s\n",
            Options.FreqOutfilename);
    fprintf(fp, "  Simulation results        : %s\n",
            Options.SimulationOutfilename);
    fprintf(fp, "  Mismatches                : %s\n",
            Options.MismatchOutfilename);
    for (i=0; i<2;i++)
        if (Options.POutformat[i])
            fprintf(fp, "  Parentage file            : %s\n", Options.POutfilename[i]);

    for (i=0; i<2;i++)
        if (Options.SiblingsOutformat[i])
        fprintf(fp, "  Siblings file             : %s\n",
                Options.SiblingsOutfilename[i]);
    if (Options.PedigreeOutformat[0])
        fprintf(fp, "  Pedigree outfile FRANz    : %s\n",
            Options.PedigreeOutfilename[0]);
    if (Options.PedigreeOutformat[1])
        fprintf(fp, "  Pedigree outfile dot      : %s\n",
            Options.PedigreeOutfilename[1]);
    if (Options.PedigreeOutformat[2])
        fprintf(fp, "  Pedigree outfile text     : %s\n",
            Options.PedigreeOutfilename[2]);

    fprintf(fp, "  MCMC logfile              : %s\n", Options.MCMClogfile);
    /*
    fprintf(fp, "  Haplotype file            : %s\n",
            Options.HaplotypeOutfilename);
    */        
    if (strlen(Options.HWETestOutfilename) > 0)
        fprintf(fp, "  Detailed HWE test results : %s\n",
                Options.HWETestOutfilename);
    else
        fprintf(fp,
                "  Detailed HWE test results : none (specify --hwetestout)\n");

    if (strlen(Options.MissingAllelesOutfilename) > 0)
        fprintf(fp, "  Missing Data statistics   : %s\n",
                Options.MissingAllelesOutfilename);
    else
        fprintf(fp,
                "  Missing Data statistics   : none (specify --missingout)\n");

    fprintf(fp, "\n\n");
    fprintf(fp, "*** Parentage Settings ***\n\n");
    fprintf(fp, "Femrepro                    : %i:%i\n", Options.FemRepro.min,
            Options.FemRepro.max);
    fprintf(fp, "Malerepro                   : %i:%i\n",
            Options.MaleRepro.min, Options.MaleRepro.max);
    fprintf(fp, "Selfing                     : ");
    OPTIONSprintBool(fp, Options.Selfing);
    fprintf(fp, "Monoecious                  : ");
    OPTIONSprintBool(fp, Options.Monoecious);

    fprintf(fp, "Use distances               : ");
    OPTIONSprintBool(fp, Data.use_distances);

/*    fprintf(fp, "Lambda Females              : %f\n", Options.LambdaFem);
    fprintf(fp, "Lambda Males or unknown Sex : %f\n", Options.LambdaMale);
    */
    if (Options.FemMaxDist > 0.)
        fprintf(fp, "Max. Distance Mother-Child  : %f\n", Options.FemMaxDist);
    if (Options.MaleMaxDist > 0.)
        fprintf(fp, "Max. Distance Father-Child  : %f\n", Options.MaleMaxDist);
    if (Options.ParentsMaxDist > 0.)
        fprintf(fp, "Max. Distance Mother-Father : %f\n", Options.ParentsMaxDist);
    if (Options.SibMaxDist > 0.)
        fprintf(fp, "Max. Distance Siblings      : %f\n", Options.SibMaxDist);

    if (Options.NumLoci > 0)
        fprintf(fp, "Number Loci                 : %i\n", Options.NumLoci);
    else
        fprintf(fp, "Number Loci                 : use all.\n");
    fprintf(fp, "Minimum typed loci          : %i\n", Options.MinTyped);
    fprintf(fp, "Max. mismatching loci Dyad  : %i\n",
            Options.MaxMismatchingDyad);
    fprintf(fp, "Max. mismatching loci Triple: %i\n",
            Options.MaxMismatchingTriple);

    fprintf(fp, "Rate of typing error\n");
    fprintf(fp, "  Average                   : %.3f\n", Data.TypingErrorAvg);
    for (i = 0; i < Data.num_loci; i++)
        fprintf(fp, "  Locus %-10s          : %.3f\n", Data.loci_ids[i], Data.TypingError[i]); /*, sqrt(1.96 * 1.96 / (Data.num_samples)));*/
    fprintf(fp, "Update Allele frequencies   : ");
    OPTIONSprintBool(fp, Options.UpdateAlleleFreqs);
    fprintf(fp, "Detect fullsibs             : ");
    OPTIONSprintBool(fp, Options.Fullsibs);
    if (Options.Fullsibs) {
        fprintf(fp, "  in parental generation    : ");
        OPTIONSprintBool(fp, Options.SiblingsInParentalGeneration);
        fprintf(fp, "p-Value Threshold           : %.3E, %.3E, %.3E\n", Options.THsiblings[0], Options.THsiblings[1], Options.THsiblings[2]);
        if (Options.SiblingsCorrectionMethod == 1)
            fprintf(fp, "p-Value Correction Method   : Benjamini-Hochberg\n");
        else 
            fprintf(fp, "p-Value Correction Method   : Holm\n");
    }

    if (Options.IgnoreSex)
        fprintf(fp, "\nSex meta data will be ignored.\n");
    if (Options.IgnoreRamets)
        fprintf(fp, "\nRamet counts will be ignored.\n");
    fprintf(fp, "\n\n");
    if (Options.MaxDepth > 0) {
        fprintf(fp, "*** Pedigree Settings ***\n\n");
        fprintf(fp, "Max. Pedigree depth         : %i\n", Options.MaxDepth);
        fprintf(fp, "\n\n");
    }
    FFLUSH(fp);
}

void
OPTIONSdumpSim(FILE * fp)
{
    fprintf(fp, "*** Files ***\n\n");
    fprintf(fp, "Input\n");
    fprintf(fp, "  Allele frequency file    : %s\n", Options.FreqInfilename);
    fprintf(fp, "Output\n");
    fprintf(fp, "  Summary                  : %s\n", Options.SimOutfilename);
    fprintf(fp, "  Genotypes                : %s\n",
            Options.GenotypeOutfilename);
    fprintf(fp, "  Pedigree outfile         : %s\n",
            Options.PedigreeOutfilename[0]);
    fprintf(fp, "  Pedigree mothers only    : %s\n",
            Options.PedigreeMotherOutfilename);
    fprintf(fp, "  Pedigree Graphviz        : %s\n",
            Options.PedigreeOutfilename[1]);
}

/* 
  vim: ft=c sw=4 ts=4 expandtab
*/
