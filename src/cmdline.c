/*
 * $Id: cmdline.c 2062 2010-05-25 10:24:32Z markus $
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
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.*
 *
 * Markus Riester, markus@bioinf.uni-leipzig.de */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <math.h>

#include "cmdline.h"

#include "global.h"
#include "dataio.h"
#include "macros.h"

extern OPTIONS Options;

static struct option long_options[] = {
    {"mhburniniter", required_argument, 0, 0},         /*  0 */
    {"sachains", required_argument, 0, 0},             /*  1 */
    {"samaxiter", required_argument, 0, 0},            /*  2 */
    {"pedigreein", required_argument, 0, 0},           /*  3 */
    {"pedigreeincheck", required_argument, 0, 0},      /*  4 */
    {"geofile", required_argument, 0, 0},              /*  5 */
    {"selfing", no_argument, 0, 0},                    /*  6 */
    {"fullsibpvth", required_argument, 0, 0},          /*  7 */
    {"maxdepth", required_argument, 0, 0},             /*  8 */
    {"pedigreeoutformat", required_argument, 0, 0},    /*  9 */
    {"freqin", required_argument, 0, 0},               /* 10 */
    {"freqout", required_argument, 0, 0},              /* 11 */
    {"femrepro", required_argument, 0, 0},             /* 12 */
    {"malerepro", required_argument, 0, 0},            /* 13 */
    {"pedigreeout", required_argument, 0, 0},          /* 14 */
    {"mintyped", required_argument, 0, 0},             /* 15 */
    {"out", required_argument, 0, 0},                  /* 16 */
    {"mcmclog", required_argument, 0, 0},              /* 17 */
    {"maxmismatching", required_argument, 0, 0},       /* 18 */
    {"typingerror", required_argument, 0, 0},          /* 19 */
    {"simiter", required_argument, 0, 0},              /* 20 */
    {"n", required_argument, 0, 0},                    /* 21 */
    {"noupdatefreqs", no_argument, 0, 0},              /* 22 */
    {"fullsibparental", no_argument, 0, 0},            /* 23 */
    {"sachi", required_argument, 0, 0},                /* 24 */
    {"sadelta", required_argument, 0, 0},              /* 25 */
    {"sabeta", required_argument, 0, 0},               /* 26 */
    {"saepsilon", required_argument, 0, 0},            /* 27 */
    {"sanepsilon", required_argument, 0, 0},           /* 28 */
    {"numloci", required_argument, 0, 0},              /* 29 */
    {"seed", required_argument, 0, 0},                 /* 30 */
    {"ignoresex", no_argument, 0, 0},                  /* 31 */
    {"nofullsibtest", no_argument, 0, 0},              /* 32 */
    {"siblingsout", required_argument, 0, 0},          /* 33 */
    {"pout", required_argument, 0, 0},                 /* 34 */
    {"hwechunks", required_argument, 0, 0},            /* 35 */
    {"hwechunksize", required_argument, 0, 0},         /* 36 */
    {"hwesteps", required_argument, 0, 0},             /* 37 */
    {"hwetestout", required_argument, 0, 0},           /* 38 */
    {"simulationout", required_argument, 0, 0},        /* 39 */
    {"N", required_argument, 0, 0},                    /* 40 */
    {"proportiontyped", required_argument, 0, 0},      /* 41 */
    {"Nr", required_argument, 0, 0},                   /* 42 */
    {"parentsmaxdist", required_argument, 0, 0},       /* 43 */
    {"coordfile", required_argument, 0, 0},            /* 44 */
    {"lociout", required_argument, 0, 0},              /* 45 */
    {"poutformat", required_argument, 0, 0},           /* 46 */
    {"femmaxdist", required_argument, 0, 0},           /* 47 */
    {"malemaxdist", required_argument, 0, 0},          /* 48 */
    {"mhiter", required_argument, 0, 0},               /* 49 */
    {"mhsamplefreq", required_argument, 0, 0},         /* 50 */
    {"cervusgenotypeout", required_argument, 0, 0},    /* 51 */
    {"fullsibin", required_argument, 0, 0},            /* 52 */
    {"simselfingrate", required_argument, 0, 0},       /* 53 */
    {"sibmaxdist", required_argument, 0, 0},           /* 54 */
    {"nogibbsmissing", no_argument, 0, 0},             /* 55 */
    {"updatefreqs", no_argument, 0, 0},                /* 56 */
    {"sacstart", required_argument, 0, 0},             /* 57 */
    {"mhswapfreq", required_argument, 0, 0},           /* 58 */
    {"mhparamsout", required_argument, 0, 0},          /* 59 */
    {"mhtemp", required_argument, 0, 0},               /* 60 */
    {"Nmax", required_argument, 0, 0},                 /* 61 */
    {"ignoreage", no_argument, 0, 0},                  /* 62 */
    {"cervusoffspringout", required_argument, 0, 0},   /* 63 */
    {"noreconstruction", no_argument, 0, 0},           /* 64 */
    {"gibbsmissing", no_argument, 0, 0},               /* 65 */
    {"fullsibtest", no_argument, 0, 0},                /* 66 */
    {"missingout", required_argument, 0, 0},           /* 67 */
    {"fullsibpvmethod", required_argument, 0, 0},      /* 68 */
    {"Nfmax", required_argument, 0, 0},                /* 69 */
    {"Nmmax", required_argument, 0, 0},                /* 70 */
    {"Nf", required_argument, 0, 0},                   /* 71 */
    {"Nm", required_argument, 0, 0},                   /* 72 */
    {"mhchains", required_argument, 0, 0},             /* 73 */
    {"genepopout", required_argument, 0, 0},           /* 74 */
    {"halfsibin", required_argument, 0, 0},            /* 75 */
    {"siblingsoutformat", required_argument, 0, 0},    /* 76 */
    {"fullsibH0", required_argument, 0, 0},            /* 77 */
    {"mismatchout", required_argument, 0, 0},          /* 78 */
    {"rmesout", required_argument, 0, 0},              /* 79 */
    {"ignoreramets", no_argument, 0, 0},               /* 80 */
    {"parenteout", required_argument, 0, 0},           /* 81 */
    {"saexactmax", required_argument, 0, 0},           /* 82 */

    {"verbose", no_argument, 0, 'v'},
    {"quiet", no_argument, 0, 'q'},
    {"help", no_argument, 0, 'h'},
    {"helpall", no_argument, 0, 'i'},
    {0, 0, 0, 0}
};

/*********************************************************************
shows usage / help screen
*********************************************************************/
static void
usage(int all)
{
    printf("NAME\n");
    printf("     %s\n\n", PACKAGE);
    printf("VERSION\n");
    printf("     %s\n\n", PACKAGE_VERSION);
    printf("SYNOPSIS\n");
    printf("     %s [Options] infile\n\n", PACKAGE);
    printf("DESCRIPTION\n");
    printf
        ("     %s reconstructs pedigrees (family trees) using polymorphic,\n",
         PACKAGE);
    printf("     codominant markers.\n\n");
    printf("OPTIONS\n");
    printf("  Prior Information:\n");
    printf("     --femrepro i:i       Age range in which females can \n");
    printf("                          reproduce (default %i:%i)\n",
           D_FEMREPROMIN, D_FEMREPROMAX);
    printf
        ("     --malerepro i:i      Age range in which males can reproduce\n");
    printf("                          (default %i:%i)\n", D_MALEREPROMIN,
           D_MALEREPROMAX);

    printf
        ("     --N i [--Nf --Nm]    Number of candidate parents (default auto). For parentage\n");
    printf
        ("                          inference, you can set the numbers for males and females.\n");
    printf
        ("     --Nmax i [--Nfmax    Maximum number of candidate parents when N is not\n");
    printf("     --Nmmax]             defined.\n");
    if (all)
        printf
            ("     --n i                Number of sampled candidate parents (default auto)\n");
    printf("     --pedigreein FILE    Filename of the pedigree input file\n");
    printf("     --fullsibin FILE     Filename of the fullsib input file\n");
    printf
        ("     --halfsibin FILE     Filename of the half or fullsib input file\n");
    printf
        ("     --geofile FILE       Pairwise distances of sampling locations\n");
    printf("     --coordfile FILE     Coordinates of sampling locations\n");
    if (all) {
        printf("     --pedigreeincheck    Filename of the true pedigree\n");
        printf("         FILE             (optional, for testing purposes)\n");
        printf("     --ignoreage          Ignores available age data\n");
        printf("     --ignoresex          Ignores available sex data\n");
    }
    printf("\n");
    printf("  Parentage Options:\n");
    printf("     --selfing            Allows selfing\n");
    printf
        ("     --mintyped i         Minimum number of typed loci (default 1+numloci/2)\n");
    printf
        ("     --maxmismatching     Maximum number of mismatching loci for\n");
    printf
        ("         i,j              single parents & parent pairs (default auto)\n");

    if (all)
        printf("     --numloci i          Use only the first i loci\n");
    printf("     --typingerror f      Rate of typing error (default %.3f)\n",
           D_TYPINGERROR);
    printf
        ("     --noreconstruction   Don't reconstruct the pedigree (just analyze the data)\n");
    printf
        ("     --femmaxdist f       Maximal distance of sampling locations\n");
    printf("     --malemaxdist f\n");
    printf("     --parentsmaxdist f\n");
    printf("     --sibmaxdist f\n");
    printf("\n");
    printf("  Fullsib Options:\n");
    printf
        ("     --[no]fullsibtest    Detect siblings (default --nofullsibtest)\n");
    printf
        ("     --fullsibparental    Detect fullsibs also in parental generation?\n");
    if (all) {
        printf
            ("     --fullsibpvmethod    The p-Value correction method for multiple testing.\n");
        printf
            ("                          1 = Benjamini-Hochberg, 2 = Holm. Default 1.\n");
        printf
            ("     --fullsibpvth f,f,f  The p-Value threshold of the sibling filter.\n");
        printf
            ("                          For 0,1 or 2 common candidate parents.\n");
        printf("                          Default is %.4E,%.4E,%.4E.\n",
               D_TH_SIBLING_0, D_TH_SIBLING_1, D_TH_SIBLING_2);
        printf("     --fullsibH0 i,i,i    Defines the null hypotheses (PO,HS,U). Examples:\n");
        printf("                          0,0,1 : FS vs. Unrelated.\n");
        printf("                          1,2,1*: FS vs. PO, 2xHS(HS,Aunt/Uncle), U.\n");
        printf("                          2,4,1#: FS vs. 2xPO (PO and OP),\n");
        printf("                                  4xHS(HS,Aunt/Uncle,Grandparent/-child), U.\n\n");
        printf("                          *default with age data\n");
        printf("                          #default without age data\n");
    }
    printf("\n");
    if (all) {
        printf("  Allele frequency Options:\n");
        printf
            ("     --freqin FILE        Filename of allele frequency input file\n");
        printf("                          (optional)\n");
        printf
            ("     --[no]updatefreqs    Update Allele Frequencies (default --updatefreqs).\n");
        printf("\n");
        printf("  Simulation Options:\n");
        printf
            ("     --simiter i          Number of simulation iterations (default %i)\n",
             D_SIMNUM);
        printf
            ("     --proportiontyped f  Proportion of typed loci (default auto)\n");
        printf
            ("     --simselfingrate f   Proportion of self-feritilization\n");
        printf("\n");
        printf("  HWE exact test options:\n");
        printf
            ("     --hwesteps i         Number of steps (default %i)\n",
             D_HWESTEPS);
        printf("     --hwechunks i        Number of chunks (default %i)\n",
               D_HWECHUNKS);
        printf("     --hwechunksize i     The chunk size (default %i)\n",
               D_HWECHUNKSIZE);
        printf("\n");
        printf("  Pedigree Constraints:\n");
        printf
            ("     --maxdepth i         Max. pedigree depth (generations)\n");
        printf("                          (default %i)\n", D_MAXDEPTH);
        printf("\n");
        printf("  MCMC Parameters:\n");
        printf("       --[no]gibbsmissing Gibbs sampling of missing data\n");
        printf("                          (default --nogibbsmissing)\n");
        printf("\n");
        printf("  SA specific Parameters:\n");
        printf("       --sachains i       Number of chains (default %i)\n",
               D_SACHAINS);
        printf
            ("       --samaxiter i      Max. number of iterations (default %i)\n",
             D_SAMAXITERATION);
        printf
            ("       --sachi f          Initial acceptance probability (default %.3f)\n",
             D_SACHI);
        printf
            ("       --sacstart f       Initial temperature (disables --sachi, default auto)\n");
        printf
            ("       --sabeta f         Neighbourhood size factor (default %.3f)\n",
             D_SABETA);
        printf("       --sadelta f        Increment (default %.3f)\n",
               D_SADELTA);
        printf
            ("       --saepsilon f      The convergence tolerance (default %.3f)\n",
             D_SAEPSILON);
        printf("       --sanepsilon i     Convergence events (default %i)\n",
               D_SANEPSILON);
        printf("       --saexactmax i     Do exhaustive enumeration instead of SA if dataset\n");
        printf("                          contains less than i+1 individuals (default %i)\n",
               D_EXACTMAXINDIVIDUALS);
        printf("\n");
        printf("  Metropolis Hastings specific  Parameters:\n");
        printf
            ("       --mhchains i       Number of chains (default number of CPU cores)\n");
        printf
            ("       --mhburniniter i   Number of burnin iterations (default %i)\n",
             D_MHBURNINITER);
        printf
            ("       --mhiter i         Number of iterations (default %i)\n",
             D_MHITERATION);
        printf
            ("       --mhsamplefreq i   Sample every ith pedigree (default %i)\n",
             D_MHSAMPLEFREQ);
        printf
            ("       --mhswapfreq i     For MCMCMC: try to swap every ith iteration (default %i)\n",
             D_MHSWAPFREQ);
        printf
            ("       --mhtemp f         For MCMCMC: the temperature of the ith chain (default %.3f)\n",
             D_MHTEMP);
        printf("\n");
    }
    printf("   Output options:\n");
    printf("     --out FILE           Filename of the summary output file\n");
    printf("                          (default %s)\n", D_OUT);
    printf
        ("     --lociout FILE       Filename of the loci summary output file\n");
    printf("                          (default %s)\n", D_LOCIOUT);
    printf
        ("     --mismatchout FILE   Filename of the mismatches output file\n");
    printf("                          (default %s)\n", D_MISMATCHOUT);
    printf
        ("     --freqout FILE       Filename of allele frequency output file\n");
    printf("                          requires --updatefreqs 0 (default %s)\n",
           D_FREQOUT);
    printf
        ("     --pout FILE          Prefix of the parentage output file(s)\n");
    printf("                          (default %s)\n", D_POUT);
    printf
        ("     --poutformat         Format(s) of the parentage outfile(s)\n");
    printf("          i,i               1: Most likely parentages (.%s)\n",
           D_CSVSUFFIX);
    printf("                            2: All with positive LOD (.%s)\n",
           D_CSVSUFFIX);
    printf("                          (default %s)\n", D_POUTFORMAT);
    printf("     --simulationout      Filename of the simulation result\n");
    printf("          FILE            file (default %s)\n", D_SIMULATIONOUT);
    printf("     --siblingsout FILE   Prefix of the of siblings output files\n");
    printf("                          (default %s)\n", D_SIBLINGSOUT);
    printf("     --siblingsoutformat  Format(s) of the siblings outfile(s)\n");
    printf("          i,i               1: FRANz format (.%s)\n",
           D_FRANZSUFFIX);
    printf("                            2: Text format (.%s)\n",
           D_TEXTSUFFIX);
    printf("                            3: CSV format (.%s)\n",
           D_CSVSUFFIX);
    printf("                          (default %s)\n", D_SIBLINGSOUTFORMAT);
    printf("     --pedigreeout FILE   Prefix of the pedigree output files\n");
    printf("                          (default %s)\n", D_PEDIGREEOUT);
    printf("     --pedigreeoutformat  Format(s) of the pedigree outfile(s)\n");
    printf("          i,i               1: FRANz format (.%s)\n",
           D_FRANZSUFFIX);
    printf("                            2: Graphviz format (.%s)\n",
           D_GVSUFFIX);
    printf("                            3: Text format (Id Sire Dam) (.%s)\n",
           D_TEXTSUFFIX);
    printf("                          (default %s)\n", D_PEDIGREEOUTFORMAT);
    printf("     --mcmclog FILE       Filename of MCMC log file\n");
    printf("                          (default %s)\n", D_MCMCLOGFILE);
    printf
        ("     --mhparamsout FILE   The MH parameter file, containing the sampled params\n");
    printf("                          (default %s)\n", D_MHPARAMSOUT);
    printf
        ("     --hwetestout FILE    Filename of the detailed HWE test results\n");
    printf
        ("     --missingout FILE    Filename of the missing data Gibbs sampler results\n");
    printf("\n");
    if (all) {
        printf("  Data conversion options:\n");
        printf
            ("     --cervusgenotypeout  Output the genotypes in CERVUS (CSV) format.\n");
        printf("          FILE\n");
        printf("     --cervusoffspringout Output a CERVUS offspring file.\n");
        printf("          FILE\n");
        printf
            ("     --parenteout FILE    Output the genotypes in PARENTE format.\n");
        printf
            ("     --genepopout FILE    Output the genotypes in Genepop format.\n");
        printf
            ("     --rmesout FILE       Output the genotypes in RMES format.\n");
        printf("\n");
    }
    printf("  Program options:\n");
    if (all)
        printf
            ("       --seed i           seed for random numbers (default: time)\n");
    printf
        ("   -v  --verbose          increase verbosity level (standard level: %d)\n",
         D_VERBOSITY);
    printf
        ("   -q  --quiet            quiet mode, no output except errors and warnings is\n");
    printf("                          generated (=verb. level 0)\n");
    if (all) {
        printf("   -h  --help             the basic options\n");
        printf("       --helpall          this help screen\n");
    } else {
        printf("   -h  --help             this help screen\n");
        printf("       --helpall          shows all options\n");
    }
    printf("\n");
    if (all) {
        printf("AUTHOR\n");
        printf("     Markus Riester       (University of Leipzig)\n");
        printf
            ("     Peter F. Stadler     (University of Leipzig, University of Vienna,\n");
        printf("                           Santa Fe Institute)\n");
        printf("     Konstantin Klemm     (University of Leipzig)\n");
        printf("\n");
        printf("REPORTING BUGS\n");
        printf("     %s\n\n", PACKAGE_BUGREPORT);
        printf("COPYRIGHT\n");
        printf
            ("     This is free software; see the source for copying conditions. There is\n");
        printf
            ("     NO warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR\n");
        printf("     PURPOSE\n\n");
        printf("Build Details:\n");
        printf("     Debugging : %i\n     Assertions: ",
#ifdef DEBUG
               1
#else
               0
#endif
            );

#ifdef NDEBUG
        printf("0\n");
#else
        printf("1\n");
#endif
        printf("     OpenMP    : ");
#ifdef HAVE_OPENMP
        printf("1\n");
#else
        printf("0\n");
#endif
    }
    exit(EXIT_SUCCESS);
}

static int
checkRange(char *cp, RANGE * r)
{
    if (sscanf(cp, "%i:%i", &r->min, &r->max) != 2)
        return 0;
    if (r->min > r->max)
        return 0;
    return 1;
}

static void
checkProbability(double f, char *optionname)
{
    if (f < 0. || f > 1.) {
        fprintf(stderr, "Option --%s not in range 0.0-1.0\n\n", optionname);
        exit(EXIT_FAILURE);
    }
}

static  bool
checkFormat(char *s, int id, char *optionname)
{
    int     tmp, res;
    char   *p, s2[LINESIZE];
    strncpy(s2, s, PATHLEN);
    p = strtok(s2, ",");
    while (p != NULL) {
        res = sscanf(p, "%i", &tmp);
        //    fprintf(stderr, "DEBUG %s, %i, %s\n", s,id,p);
        if (res != 1) {
            fprintf(stderr,
                    "Option --%s not a list of numbers (e.g. \"1,2,3\")\n\n",
                    optionname);
            exit(0);
        }
        if (tmp == id)
            return true;
        p = strtok(NULL, ",");
    }
    return false;
}

/******************************************************************
parses command line
******************************************************************/
void
CMDLINEparse(int argc, char **argv)
{
    int     c, i, j;
    int     option_index;
    float   tmp;
    RANGE   r;
    char    pedigreeoutformat[LINESIZE], pedigreeout[LINESIZE],
        poutformat[LINESIZE], pout[LINESIZE],
        siblingsoutformat[LINESIZE], siblingsout[LINESIZE];

    strncpy(pedigreeoutformat, D_PEDIGREEOUTFORMAT, LINESIZE);
    strncpy(pedigreeout, D_PEDIGREEOUT, LINESIZE);
    strncpy(poutformat, D_POUTFORMAT, LINESIZE);
    strncpy(pout, D_POUT, LINESIZE);
    strncpy(siblingsout, D_SIBLINGSOUT, LINESIZE);
    strncpy(siblingsoutformat, D_SIBLINGSOUTFORMAT, LINESIZE);

    while (1) {
        option_index = 0;

        c = getopt_long(argc, argv, "fvqh", long_options, &option_index);
        if (c == -1)
            break;

        switch (c) {
            case 0:            /* long options */
                switch (option_index) {
                    case 0:
                        if (sscanf(optarg, "%u", &Options.MHBurnin) == 0)
                            FATAL("Option --mhburniniter requires argument.");
                        break;
                    case 1:
                        if (sscanf(optarg, "%d", &Options.SAchains) == 0)
                            FATAL("Option --sachains requires argument.");
                        break;
                    case 2:
                        if (sscanf(optarg, "%u", &Options.SAmaxIterations) ==
                            0)
                            FATAL("Option --samaxiter requires argument.");
                        break;
                    case 3:
                        strncpy(Options.PedigreeInfilename, optarg, PATHLEN);
                        break;
                    case 4:
                        strncpy(Options.PedigreeCheckInfilename, optarg,
                                PATHLEN);
                        break;
                    case 5:
                        strncpy(Options.Geofilename, optarg, PATHLEN);
                        break;
                    case 6:
                        Options.Selfing = true;
                        break;
                    case 7:
                        if (sscanf
                            (optarg, "%lf,%lf,%lf", &Options.THsiblings[0],
                             &Options.THsiblings[1],
                             &Options.THsiblings[2]) != 3)
                            FATAL
                                ("Option --fullsibpvth not valid (fdr0,fdr1,fdr2, e.g 0.05,0.01,0.005)");
                        break;
                    case 8:
                        if (sscanf(optarg, "%hd", &Options.MaxDepth) == 0)
                            FATAL("Option --maxdepth requires argument.");
                        break;
                    case 9:
                        strncpy(pedigreeoutformat, optarg, LINESIZE);
                        break;
                    case 10:
                        strncpy(Options.FreqInfilename, optarg, PATHLEN);
                        break;
                    case 11:
                        strncpy(Options.FreqOutfilename, optarg, PATHLEN);
                        break;
                    case 12:
                        if (checkRange(optarg, &r) == 0)
                            FATAL
                                ("Option --femrepro not a valid range (int:int)");
                        Options.FemRepro = r;
                        break;
                    case 13:
                        if (checkRange(optarg, &r) == 0)
                            FATAL
                                ("Option --malerepro not a valid range (int:int)");
                        Options.MaleRepro = r;
                        break;
                    case 14:
                        strncpy(pedigreeout, optarg, LINESIZE);
                        break;
                    case 15:
                        if (sscanf(optarg, "%d", &Options.MinTyped) == 0)
                            FATAL("Option --mintyped requires argument.");
                        break;
                    case 16:
                        strncpy(Options.Outfilename, optarg, PATHLEN);
                        break;
                    case 17:
                        strncpy(Options.MCMClogfile, optarg, PATHLEN);
                        break;
                    case 18:
                        if (sscanf(optarg, "%i,%i", &i, &j) == 0)
                            FATAL
                                ("Option --maxmismatching requires argument.");
                        if (i < 0)
                            FATAL
                                ("Option --maxmismatching single parent < 0.");
                        if (j < 0)
                            FATAL("Option --maxmismatching parent pair < 0.");
                        Options.MaxMismatchingDyad = i;
                        Options.MaxMismatchingTriple = j;
                        Options.MaxMismatchingDefined = true;
                        break;
                    case 19:
                        if (sscanf(optarg, "%f", &tmp) == 0)
                            FATAL("Option --typingerror requires argument.");
                        checkProbability(tmp, "typingerror");
                        Options.TypingErrorDefined = tmp;
                        break;
                    case 20:
                        if (sscanf(optarg, "%u", &Options.SimulationSteps) ==
                            0)
                            FATAL("Option --simiter requires argument.");
                        break;
                    case 21:
                        if (sscanf(optarg, "%i", &Options.n) == 0)
                            FATAL("Option --n requires argument.");
                        Options.ndefined = true;
                        break;
                    case 22:
                        Options.UpdateAlleleFreqs = false;
                        break;
                    case 23:
                        Options.SiblingsInParentalGeneration = true;
                        break;
                    case 24:
                        if (sscanf(optarg, "%f", &tmp) == 0)
                            FATAL("Option --sachi requires argument.");
                        checkProbability(tmp, "sachi");
                        Options.SAchi = tmp;
                        break;
                    case 25:
                        if (sscanf(optarg, "%f", &tmp) == 0)
                            FATAL("Option --sadelta requires argument.");
                        Options.SAdelta = tmp;
                        break;
                    case 26:
                        if (sscanf(optarg, "%lf", &Options.SAbeta) == 0)
                            FATAL("Option --sabeta requires argument.");
                        break;
                    case 27:
                        if (sscanf(optarg, "%lf", &Options.SAepsilon) == 0)
                            FATAL("Option --saepsilon requires argument.");
                        break;
                    case 28:
                        if (sscanf(optarg, "%u", &Options.SANepsilon) == 0)
                            FATAL("Option --sanepsilon requires argument.");
                        break;
                    case 29:
                        if (sscanf(optarg, "%i", &Options.NumLoci) == 0)
                            FATAL("Option --numloci requires argument.");
                        break;
                    case 30:
                        if (sscanf(optarg, "%u", &Options.Seed) == 0)
                            FATAL("Option --seed requires argument.");
                        break;
                    case 31:
                        Options.IgnoreSex = true;
                        break;
                    case 32:
                        Options.Fullsibs = false;
                        break;
                    case 33:
                        strncpy(siblingsout, optarg, PATHLEN);
                        break;
                    case 34:
                        strncpy(pout, optarg, LINESIZE);
                        break;
                    case 35:
                        if (sscanf(optarg, "%i", &Options.HWEChunks) == 0)
                            FATAL("Option --hwechunks requires argument.");
                        break;
                    case 36:
                        if (sscanf(optarg, "%i", &Options.HWEChunkSize) == 0)
                            FATAL("Option --hwechunksize requires argument.");
                        break;
                    case 37:
                        if (sscanf(optarg, "%i", &Options.HWESteps) == 0)
                            FATAL("Option --hwesteps requires argument.");
                        break;
                    case 38:
                        strncpy(Options.HWETestOutfilename, optarg, PATHLEN);
                        break;
                    case 39:
                        strncpy(Options.SimulationOutfilename, optarg,
                                PATHLEN);
                        break;
                    case 40:
                        if (sscanf(optarg, "%i", &Options.Nf) == 0)
                            FATAL("Option --N requires argument.");
                        Options.Nm = Options.Nf;
                        Options.Nfdefined = true;
                        Options.Nmdefined = true;
                        break;
                    case 41:
                        if (sscanf(optarg, "%f", &tmp) == 0)
                            FATAL
                                ("Option --proportiontyped requires argument.");
                        checkProbability(tmp, "proportiontyped");
                        Options.ProportionTypedDefined = tmp;
                        break;
                    case 42:
                        if (sscanf(optarg, "%i", &Options.Nr) == 0)
                            FATAL("Option --Nr requires argument.");
                        break;
                    case 43:
                        if (sscanf(optarg, "%lf", &Options.ParentsMaxDist) ==
                            0)
                            FATAL
                                ("Option --parentsmaxdist requires argument.");
                        break;
                    case 44:
                        strncpy(Options.Coordfilename, optarg, PATHLEN);
                        break;
                    case 45:
                        strncpy(Options.LociOutfilename, optarg, PATHLEN);
                        break;
                    case 46:
                        strncpy(poutformat, optarg, LINESIZE);
                        break;
                    case 47:
                        if (sscanf(optarg, "%lf", &Options.FemMaxDist) == 0)
                            FATAL("Option --femmaxdist requires argument.");
                        break;
                    case 48:
                        if (sscanf(optarg, "%lf", &Options.MaleMaxDist) == 0)
                            FATAL("Option --malemaxdist requires argument.");
                        break;
                    case 49:
                        if (sscanf(optarg, "%u", &Options.MHIterations) == 0)
                            FATAL("Option --mhiter requires argument.");
                        break;
                    case 50:
                        if (sscanf(optarg, "%u", &Options.MHsamplefreq) == 0)
                            FATAL("Option --mhsamplefreq requires argument.");
                        break;
                    case 51:
                        strncpy(Options.GenotypeOutfilenameCERVUS, optarg,
                                PATHLEN);
                        break;
                    case 52:
                        strncpy(Options.FullsibInfilename, optarg, PATHLEN);
                        Options.FullsibInfile = true;
                        break;
                    case 53:
                        if (sscanf(optarg, "%f", &tmp) == 0)
                            FATAL
                                ("Option --simselfingrate requires argument.");
                        checkProbability(tmp, "simselfingrate");
                        Options.SimulationSelfingRate = tmp;
                        break;
                    case 54:
                        if (sscanf(optarg, "%lf", &Options.SibMaxDist) == 0)
                            FATAL("Option --sibmaxdist requires argument.");
                        break;
                    case 55:
                        Options.GibbsMissingData = false;
                        break;
                    case 56:
                        Options.UpdateAlleleFreqs = true;
                        break;
                    case 57:
                        if (sscanf(optarg, "%f", &tmp) == 0 || tmp <= EPSILON)
                            FATAL("Option --sacstart requires argument.");
                        Options.SAc0 = tmp;
                        break;
                    case 58:
                        if (sscanf(optarg, "%u", &Options.MHswapfreq) == 0)
                            FATAL("Option --mhswapfreq requires argument.");
                        break;
                    case 59:
                        strncpy(Options.MHparamsfile, optarg, PATHLEN);
                        break;
                    case 60:
                        if (sscanf(optarg, "%f", &tmp) == 0)
                            FATAL("Option --mhtemp requires argument.");
                        checkProbability(tmp, "mhtemp");
                        Options.MHTemp = tmp;
                        break;
                    case 61:
                        if (sscanf(optarg, "%i", &Options.Nfmax) == 0)
                            FATAL("Option --Nmax requires argument.");
                        Options.Nmmax = Options.Nfmax;
                        break;
                    case 62:
                        Options.IgnoreAge = true;
                        break;
                    case 63:
                        strncpy(Options.OffspringOutfilenameCERVUS, optarg,
                                PATHLEN);
                        break;
                    case 64:
                        Options.NoReconstruction = true;
                        break;
                    case 65:
                        Options.GibbsMissingData = true;
                        break;
                    case 66:
                        Options.Fullsibs = true;
                        WARN("EXPERIMENTAL FEATURE --fullsibtest. Please read the manual and set parameters manually")
                        break;
                    case 67:
                        strncpy(Options.MissingAllelesOutfilename, optarg,
                                PATHLEN);
                        break;
                    case 68:
                        if (sscanf
                            (optarg, "%hd",
                             &Options.SiblingsCorrectionMethod) == 0)
                            FATAL
                                ("Option --fullsibpvmethod requires argument.");
                        break;
                    case 69:
                        if (sscanf(optarg, "%d", &Options.Nfmax) == 0)
                            FATAL("Option --Nfmax requires argument.");
                        break;
                    case 70:
                        if (sscanf(optarg, "%d", &Options.Nmmax) == 0)
                            FATAL("Option --Nmmax requires argument.");
                        break;
                    case 71:
                        if (sscanf(optarg, "%d", &Options.Nf) == 0)
                            FATAL("Option --Nf requires argument.");
                        Options.Nfdefined = true;
                        break;
                    case 72:
                        if (sscanf(optarg, "%d", &Options.Nm) == 0)
                            FATAL("Option --Nm requires argument.");
                        Options.Nmdefined = true;
                        break;
                    case 73:
                        if (sscanf(optarg, "%d", &Options.MHchains) == 0)
                            FATAL("Option --mhchains requires argument.");
                        break;
                    case 74:
                        strncpy(Options.GenotypeOutfilenameGenepop, optarg,
                                PATHLEN);
                        break;
                    case 75:
                        strncpy(Options.FullsibInfilename, optarg, PATHLEN);
                        break;
                    case 76:
                        strncpy(siblingsoutformat, optarg, LINESIZE);
                        break;
                    case 77:
                        if (sscanf(optarg, "%d,%d,%d", &Options.FullsibAlternatives[0],&Options.FullsibAlternatives[1], &Options.FullsibAlternatives[2]) < 3)
                            FATAL("Option --fullsibH0 requires format <i,i,i>.");
                        Options.FullsibAlternativesDefined = true;
                        break;
                    case 78:
                        strncpy(Options.MismatchOutfilename, optarg,
                                PATHLEN);
                        break;
                    case 79:
                        strncpy(Options.GenotypeOutfilenameRMES, optarg,
                                PATHLEN);
                        break;
                    case 80:
                        Options.IgnoreRamets = true;
                        break;
                    case 81:
                        strncpy(Options.GenotypeOutfilenamePARENTE, optarg,
                                PATHLEN);
                        break;
                    case 82:
                        if (sscanf(optarg, "%d", &Options.SAexactMax) == 0)
                            FATAL("Option --saexactmax requires argument.");
                        break;
                }
                break;
            case 'v':
                Options.Verbosity = 2;
                break;
            case 'q':
                Options.Verbosity = 0;
                break;
            case 'h':
                usage(0);
                break;
            case 'i':
                usage(1);
                break;
            default:
                FATAL("Type -h for options.");
        }
    }
    if (optind < argc)
        strcpy(Options.Infilename, argv[optind++]);
    if ((strlen(Options.Infilename) == 0))
        FATAL("Please specify input file.");
    if (Options.SimulationSelfingRate >= 0. && !Options.Selfing)
        FATAL("--simselfingrate require selfing.");

    if (strlen(Options.MissingAllelesOutfilename) > 0
        && !Options.GibbsMissingData)
        FATAL("--missingout requires --gibbsmissing.");

    if (strlen(Options.FreqInfilename) > 0 && Options.UpdateAlleleFreqs) {
        WARN("Will not update allele frequencies.");
        Options.UpdateAlleleFreqs = 0;
    }
    if (strlen(Options.Geofilename) <= 0 && strlen(Options.Coordfilename) <= 0
        && (Options.FemMaxDist >= 0 || Options.MaleMaxDist >= 0
            || Options.SibMaxDist >= 0))
        FATAL("(Fem|Male|Sib)MaxDist require distances.");

    if ((!Options.Nfdefined && Options.Nfmax < 0)
        || (!Options.Nmdefined && Options.Nmmax < 0)) {
        WARN
            ("Neither --N, --Nmax, --Nf or --Nfmax and  --Nm or --Nmmax defined.\nAssuming complete sampling!");
        Options.AssumeCompleteSample = true;
        Options.Nfdefined = true;
        Options.Nmdefined = true;
    }
    for (i = 0; i < 3; i++)
        Options.PedigreeOutformat[i] =
            checkFormat(pedigreeoutformat, i + 1, "pedigreeoutformat");
    for (i = 0; i < 2; i++)
        Options.POutformat[i] = checkFormat(poutformat, i + 1, "poutformat");
    for (i = 0; i < 3; i++)
        Options.SiblingsOutformat[i] = checkFormat(siblingsoutformat, i + 1, "siblingsoutformat");

    snprintf(Options.PedigreeOutfilename[0], PATHLEN, "%s.%s", pedigreeout,
             D_FRANZSUFFIX);
    snprintf(Options.PedigreeOutfilename[1], PATHLEN, "%s.%s", pedigreeout,
             D_GVSUFFIX);
    snprintf(Options.PedigreeOutfilename[2], PATHLEN, "%s.%s", pedigreeout,
             D_TEXTSUFFIX);
    snprintf(Options.POutfilename[OUTFORMAT_CSV], PATHLEN, "%s.%s", pout,
             D_CSVSUFFIX);
    snprintf(Options.POutfilename[OUTFORMAT_DETAIL], PATHLEN, "%s.%s", pout,
             D_CSVSUFFIX);
    snprintf(Options.SiblingsOutfilename[0], PATHLEN, "%s.%s", siblingsout,
             D_FRANZSUFFIX);
    snprintf(Options.SiblingsOutfilename[1], PATHLEN, "%s.%s", siblingsout,
             D_TEXTSUFFIX);
    snprintf(Options.SiblingsOutfilename[2], PATHLEN, "%s.%s", siblingsout,
             D_CSVSUFFIX);

    Options._MinReproRange = MIN( (Options.FemRepro.max-Options.FemRepro.min),
                                  (Options.MaleRepro.max-Options.MaleRepro.min) );

}

/* 
  vim: ft=c sw=4 ts=4 expandtab
*/
