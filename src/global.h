#ifndef GLOBAL_H
#define GLOBAL_H
/*
 * $Id: global.h 2062 2010-05-25 10:24:32Z markus $
 *
 * Global stuff...
 *
 * Copyright (C) 2008 Universitaet Leipzig  
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
 * Markus Riester, markus@bioinf.uni-leipzig.de */

#include <stdbool.h>
#include <stdio.h>
#include <float.h>
#include <limits.h>
#include <time.h> 
#include "listparentage.h"
#include "uthash.h"
#include "ibd.h"
#include "pvalue.h"
#include "../libdir/dcmt/include/dc.h"

#define LINESIZE 65536
#define PATHLEN  256
#define PROGRESSBARSIZE 20
#define OUTFORMAT_CSV    0
#define OUTFORMAT_DETAIL 1 
#define MINUSINF -FLT_MAX
#define EPSILON 0.00001
#define NUL '\0'

/* default values: */
#define D_VERBOSITY 1
#define D_PLOIDY 2   
#define D_SELFING 0
#define D_MINSAMPLINGRATE 0.001 /* minimum estimated sampling rate
                                   (to get always non-zero sampling rates) */

/* the relative FDR for the sibling calculations */
#define D_PV_CORRECTION_METHOD 1 /* 1 BH, 2 holm */
#define D_TH_SIBLING_2 0.05
#define D_TH_SIBLING_1 0.001
#define D_TH_SIBLING_0 0.001
#define D_TH_SIBLING_IND 1.66
#define D_SIBLING_IN_PARENTAL_GENERATION false  /* try to detect siblings also
                                                   in parental generation?*/

#define D_TYPINGERROR 0.01
#define D_TYPINGERRORMIN 0.005 /* for auto (calc by known mothers), min error */
#define D_MAXDEPTH -1
#define D_MINTYPED -1  /* < 0 means num loci / 2 +1 */
#define D_MINPOSTERIORPROB 0.0005
#define D_SAMPLINGRATE 0.25
#define D_CONFSTRICT 0.95
#define D_CONFRELAXED 0.80
#define D_MAXMISMATCHING 1
#define D_NUMLOCI -1
#define D_MISMATCHOUTTH 0.05 /* filter all mismatches < 0.05 % */

/* filenames*/
#define D_OUT "summary.txt"        
#define D_LOCIOUT "locisummary.txt"
#define D_SIBLINGSOUT "siblings"
#define D_SIBLINGSOUTFORMAT "2"
#define D_FREQOUT "allelefreqs.dat" 
#define D_POUT "parentage"
#define D_POUTFORMAT "1"
#define D_PEDIGREEOUT "pedigree"
#define D_PEDIGREEOUTFORMAT "1,2"


/* file suffices */
#define D_FRANZSUFFIX "dat"
#define D_GVSUFFIX "dot"
#define D_TEXTSUFFIX "txt"
#define D_CSVSUFFIX "csv"

#define D_MCMCLOGFILE "mcmc.log"
#define D_MHDISTSOUT "mhdists.dat"
#define D_MHPARAMSOUT "mhparams.dat"
#define D_HAPLOOUT "haplotypes.dat"
#define D_GENOTYPEOUT "simgenotypes.dat"
#define D_SIMULATIONOUT "simulation.txt"
#define D_MISMATCHOUT "mismatches.txt"

#define D_MALEREPROMAX 1000
#define D_FEMREPROMAX 1000
#define D_MALEREPROMIN 0
#define D_FEMREPROMIN 0

/* Parser settings */
#define D_PARSE_DESC_LENGTH 10            /* length of the descriptions */
#define D_PARSE_MAX_NUM_ALLELES 20000     /* max. number of alleles     */
#define D_PARSE_MAX_ALLELE_RANGE 20000    /* max - min + 1 repeat count */
#define D_PARSE_MAX_ALLELE_VALUE 99999 
#define D_PARSE_MAX_ALLELE_VALUE_STRLEN 5 /* the length of the max. allele value */
#define D_PARSE_MAX_SAMPLING_LOCATIONS 100000 
#define D_PARSE_MAX_LOCI 1000000
#define D_PARSE_MAX_INDIVIDUALS 1000000  /* per sampling location */
#define D_PARSE_MAX_RAMETS 2000  

/* simulation defaults */
#define D_SIMNUM 50000
#define D_SIMNUMCANDIDATES 5000

/* we integrate N over a grid from n to Nmax and this defines the grid size */
#define D_MCMCCPGRID 200.
#define D_MCMCTEGRID 25

#define D_EXACTMAXINDIVIDUALS 25   /* exhaustive enumeration */

/* SA defaults, see Almudevar paper */
#define D_SACHAINS 2
#define D_SAMAXITERATION 100000000
#define D_SAMAXINITTEMP 2000
#define D_SAADJUSTRUNS 10000
#define D_SACHI 0.95
#define D_SABETA 3.0 
#define D_SADELTA 0.1
#define D_SAEPSILON 0.001
#define D_SANEPSILON 3
#define D_SAPRIORRUNS 2
#define D_SAPRIORCHI 0.25
#define D_SAUPDATEFREQS false

/* MH defaults */
#define D_MHRANDOMITER 2000    /* accept the first n iterations every change */
#define D_MHBURNINITER 500000  /* do not sample pedigrees here               */
#define D_MHITERATION 3000000  /* stop after n iterations                    */
#define D_MHSAMPLEFREQ 10      /* thinning interval                          */
#define D_MHSWAPFREQ 25        /* swapping interval                          */
#define D_MHNSWAPS -1
#define D_MHTEMP 0.2

/* HWE test. See hwe.c and the paper */
#define D_HWESTEPS 1000
#define D_HWECHUNKS 200
#define D_HWECHUNKSIZE 1000

/* Null alleles test */
#define D_NULLSTEPS 50000         /* maximal number of iterations   */
#define D_NULLEPSILON 0.0000001   /* convergence criteria           */
#define D_NULLNEPSILON 10         

/* error rate estimation. this defines the minimum number of known  */
/* parent offspring relationships                                   */
#define D_ERRORZ 1.96       /* z-value of the error rate estimation */
#define D_ERRORE 0.2        /* with 10% error                       */

/* exclusion probabilities */
#define D_EXCLP_OFFSPRING 7 /* for how many offspring?              */

typedef struct S_RANGE
{
    int min;
    int max;
    int num_values;
} RANGE;

typedef struct S_POPULATION
{
    int id;
    char description[D_PARSE_DESC_LENGTH+1];
    double latitude;
    double longitude;
} POPULATION;

/* a lookup table for the sample ids */
typedef struct S_DESCRIPTION_HASH {
    char name[D_PARSE_DESC_LENGTH+1];
    int id;
    UT_hash_handle hh;
} DESCRIPTION_HASH;

typedef struct S_LOCUS_COORD {
    int id;
    int locus;
    int allele;              /* 0 = first allele, 1 = second allele */
    bool both;               /* both alleles missing? */
} LOCUS_COORD;

typedef struct S_SAMPLE
{
    int id;                  /* the internal id of the sample               */
    int **genotype_obs;      /* the bserved genotype                        */
    int ramets;              /* how often this genotype was sampled         */
    char description[D_PARSE_DESC_LENGTH+1];    /* a short description of the sample           */
    int sex;                 /* <= 0: unkown; 1: female; 2: male            */ 
    int birth;               /* < 0 unknown                                 */ 
    int death;              
    int typed_loci;          /* the number of completely typed loci         */
    POPULATION *population;  /*pointer to the population of this sample     */
} SAMPLE;

typedef struct S_HWERESULT
{
  double p_value;  /* mean p-value */
  double se;       /* standard error of the p-value */
  int swch_count[3];  /* switch counts for partial and full switch */
} HWERESULT;

typedef struct S_FREQS
{
    double *freqs;
    int *allele_id_lookup; /* lowest allele has id 0, highest id n_alleles-1 */
    int min;               /* value of "smallest" allele                     */
    int max;               /* of largest allele                              */ 
    int num_ind;           /* number of individuals contributing in est.     */
} FREQS;

typedef struct S_POSTERIORS
{
    int v;           /* first putative parent,  negative value for no parent */
    int w;           /* second putative parent, negative value for no parent */
    double l;        /* likelihood L_i = P(data|v,w)                         */
    double lw;       /* some above, but weighted with sampling rate          */
    double p_opt;    /* parentage probability to optimize                    */
    double lod;      /* the LOD score                                        */
    int observed;    /* observed count in metropolis sampling                */
} POSTERIORS;

typedef struct S_FULLSIB 
{
    IBD_LL ibd;
    IBD_LL ibd_pv;
    IBD_LL ibd_pv_adjusted;
    short common_parents;
    bool fullsibs;
    bool halfsibs;
    bool excluded;
    bool indirect;
} FULLSIB;

typedef struct S_DATA
{
    char description[D_PARSE_DESC_LENGTH + 1];    
    SAMPLE **samples;        /* one for each population */
    int num_populations;
    int num_samples;         /* samples is the number of different genotypes */
    int num_ramets;          /* .. and ramets */
    int num_loci;
    int *num_samples_in_population;
    int *num_ramets_in_population;
    int **alleles;          /* array of alleles at a locus */
    int *num_alleles;
    int max_num_alleles;    /* the number of alleles at the highest polymorphic loci */
    int *num_typed;         /* number of completely typed genotypes at locus */
    HWERESULT *hweresults;
    SAMPLE *id_mapping;
    POPULATION *populations;
    double **population_distances;
    double max_distance;
    bool use_distances;          /* Flag indicating that user loaded a 
                                   geodist file                              */
    bool use_pedigree;           /* Flag indicating that user loaded a 
                                   pedigree file                             */
    FREQS *allele_frequencies;  /* a pointer to the global allele 
                                   frequencies                               */
    FULLSIB *fs_matrix;         /* NxN matrix. 1=(x,y) are fullsibs          */
    int **in_pedigree;
    int ***in_pedigree_mismatches;
    int *in_pedigree_relationships;
    bool has_age_data;
    bool has_missing_birth_data;
    bool has_sex_data;
    /* loci summary statistics:                                              */
    double *h_obs;              /* The observed heteroyzgosity at locus i    */
    double *h_exp;              /* The expected heterozygosity at locus i    */
    double *pic;                /* The PIC                                   */
    double *pid;                /* The unbiased prob. of identity            */
    double *pid_sib;            /* The prob. of identity for siblings        */
    double **exclp_sp;          /* The probability of exclusion at locus i   */
    double **exclp_mk;          /* mother known                              */
    double **exclp_pp;          /* mother unknown                            */
    double **exclp_fullsibs;    /* n unrelated ind. excluded as fullsibs     */
    double *p_null;             /* expected null allele frequency            */
    double *TypingError;
    double *_TypingErrorEst;    /* with pedigreein, we can est. the error    */
    double TypingErrorAvg;
    double SelfingRateAvg;      /* the estimated sr based on Hobs/Hexp       */
    
    /* missing data */
    int *locus_has_missing_data; 
    LOCUS_COORD *_missing_data;  /* lookup table for  quickly finding missing alleles */
    unsigned int _num_missing_data;
    char **loci_ids;
    bool **ignore_ary;/* should we ignore loci?                               */
    time_t StartTime;
    DESCRIPTION_HASH *sample_ids;
} DATA;

typedef struct S_ARYSTATS {
    unsigned int num;
    double avg;
    double var;
} ARYSTATS;

typedef struct S_SIMRESULTS {
    PVALUE_LOOKUP pv_lup_po;
    PVALUE_LOOKUP pv_lup_hs;
    PVALUE_LOOKUP pv_lup_u;
} SIMRESULTS;

typedef struct S_PROBS
{
    int *num_candidates;        /* the total number of candidate parents for 
                                   offspring i */
    int *num_candidates_f;      /* total number of candidate mothers */
    int *num_candidates_m;      /* total number of candidate fathers */
    double mean_candidates_f;
    double mean_candidates_m;
    int **candidates;           /* the IDs of all candidate parents for 
                                   offspring i*/
    bool *is_candidate_parent;  /* true when ind. i is a candidate parents, 
                                   so false for offspring generation            */
    unsigned int num_candidate_mothers; 
                                /* the total number of females that could be a
                                 * mother   */ 
    unsigned int num_candidate_fathers;
                                /* the same for males */

    POSTERIORS **posteriors;  /* this is a data-structure that contains 
                                 the probabilities of all valid arcs in 
                                 the pedigree                             */

    int *num_posteriors;      /* all for the pedigree reconstruction      */
    bool *are_connected;  /* are the two in one putative parent-offspring 
                              relationship? for updating calculations      */
    unsigned int mh_sampled_pedigrees;
    SIMRESULTS simresults;

    /* temporary data-structures during LOD score calculation. We later 
     * work with arrays, not lists anymore for performance reasons.       */
    LISTPARENTAGENODE **_triples;  /* temporarily list of all valid       */
    LISTPARENTAGENODE **_dyads;    /* parent(s)-offspring pairs/triples   */
    int *_num_triples;
    int *_num_dyads;
    double *_sum_filtered_likelihood_dyads_f;
    double *_sum_filtered_likelihood_dyads_m;
    double *_sum_filtered_likelihood_triples;
} PROBS;

typedef struct S_OPTIONS
{
    char Infilename[PATHLEN];
    char Outfilename[PATHLEN];
    char SimOutfilename[PATHLEN];
    char LociOutfilename[PATHLEN];
    char SiblingsOutfilename[3][PATHLEN];
    char GenotypeOutfilename[PATHLEN];
    char GenotypeOutfilenameCERVUS[PATHLEN];
    char OffspringOutfilenameCERVUS[PATHLEN];
    char GenotypeOutfilenameGenepop[PATHLEN];
    char GenotypeOutfilenameRMES[PATHLEN];
    char GenotypeOutfilenamePARENTE[PATHLEN];
    char Geofilename[PATHLEN];
    char Coordfilename[PATHLEN];
    char FullsibInfilename[PATHLEN];
    char MCMClogfile[PATHLEN];
    char MHparamsfile[PATHLEN];
    char FreqInfilename[PATHLEN];
    char FreqOutfilename[PATHLEN];
    char POutfilename[2][PATHLEN];
    char PedigreeInfilename[PATHLEN];
    char PedigreeCheckInfilename[PATHLEN];
    char PedigreeMotherOutfilename[PATHLEN];
    char HaplotypeOutfilename[PATHLEN];
    char HWETestOutfilename[PATHLEN];
    char SimulationOutfilename[PATHLEN];
    char MissingAllelesOutfilename[PATHLEN];
    char MismatchOutfilename[PATHLEN];
    char PedigreeOutfilename[3][PATHLEN];
    bool FullsibInfile;
    bool PedigreeOutformat[3];
    bool POutformat[2];
    bool SiblingsOutformat[3];
    int Verbosity;
    bool Selfing;
    double FemMaxDist;
    double MaleMaxDist;
    double ParentsMaxDist;
    double SibMaxDist;
    double TypingErrorDefined;
    short MaxDepth;
    int MinTyped;
    double *ProportionTyped;
    double ProportionTypedDefined;
    int MaxMismatching;
    int MaxMismatchingDyad;
    int MaxMismatchingTriple;
    bool MaxMismatchingDefined;
    unsigned int SimulationSteps;
    int Nf;
    int Nm;
    int Nr;
    int n;
    bool Nfdefined;
    bool Nmdefined;
    bool ndefined;
    int Nfmax;
    int Nmmax;
    double SimulationSelfingRate;
    double SAchi;
    double SAbeta;
    double SAdelta;
    double SAepsilon;
    double SAc0;
    unsigned int SANepsilon;
    int SAexactMax;
    int SAchains;
    int MHchains;
    unsigned int MHBurnin;
    int MHNswaps;
    unsigned int MHIterations;
    unsigned int MHsamplefreq;
    unsigned int MHswapfreq;
    double MHTemp;
    unsigned int SAmaxIterations;
    int NumLoci;
    unsigned int Seed;
    bool UpdateAlleleFreqs;
    bool IgnoreAge;
    bool IgnoreSex;
    bool IgnoreRamets;
    bool Fullsibs;
    int FullsibAlternatives[3];
    bool FullsibAlternativesDefined;
    int HWESteps;
    int HWEChunks;
    int HWEChunkSize;
    bool NoReconstruction;
    bool RametAlleleFreqs;
    bool GibbsMissingData;
    bool GibbsTypingError;
    bool CollapseRamets;
    bool SiblingsInParentalGeneration;
    RANGE FemRepro;
    RANGE MaleRepro;
    int _MinReproRange;
    short SiblingsCorrectionMethod;
    double THsiblings[3];
    bool Monoecious;
    bool MonoeciousDefined;
    bool AssumeCompleteSample;
    mt_struct **_MTS;
} OPTIONS;

#define HAS_KNOWN_MOTHER(id) (Data.use_pedigree && Data.in_pedigree[(id)][0] >= 0)

/* 
  vim: ft=c sw=4 ts=4 expandtab
*/

#endif
