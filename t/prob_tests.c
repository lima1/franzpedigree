/*
 * $Id: prob_tests.c 1873 2009-12-21 13:59:20Z markus $
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
    int i,j,child[2], p1[2], p2[2],n;
    int randres[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    double p, pl[13], m0, m1, m2, pg, ptmp;
    double randdblres[2], dmin[2], dmax[2], sr_f, sr_m;
    char* filename = FILENAME; 
    char* filenameMD = FILENAMEMD;
    char* filenamePM = FILENAMEPM;
    char* filenamePMDIST = FILENAMEPMDIST;
#ifdef GENRANDFILES    
    FILE *fprandint, *fprandone, *fprand100, *fpnorm;
#endif    
    LOCUS_COORD lc;

    SIMULATION sim;       /* the sampling simulation data          */
    MCMC mcmc;
    MCMC_STATE state;

    plan_tests(192);
    OPTIONSinit();
#ifdef HAVE_OPENMP
    n = omp_get_max_threads();
#else
    n=1;
#endif
    
    Options.Verbosity = 0;
    SRAND(123,i, n);
    Options.RametAlleleFreqs = true;
    DATAIOinfile(filename);
    Data.TypingError[1] = 0.02;
    Data.TypingError[2] = 0.03;
    Data.TypingError[3] = 0.04;
    Data.TypingError[4] = 0.05;
    Data.TypingError[5] = 0.06;
    Data.TypingError[6] = 0.07;
    Data.locus_has_missing_data[0] = true;

    /*Options.MaxMismatching = 0;*/
    /* calculate dyads, they require allele frequencies */
    FREQcalcAlleleFrequencies();
    FREQcalcSummaryStatistics();
    PROBinit();
    sim = SIMstart();
    
    PROBcalcTriplesAndDyads();
    ok(Probs._triples[1]->v == 3, "parent a of sample 1 is 3");
    ok(Probs._triples[1]->w == 4, "parent a of sample 1 is 4");
    /*printf("%f %f\n",Probs.triples[1]->p,log(pow(0.5,7)));*/


    /* number of triples correct? */
    ok(Probs._num_triples[0] == 0, "no putative parents for sample 0");
    ok(Probs._num_triples[1] == 2, "putative parents for sample 1");
    ok(Probs._num_triples[2] == 6, "putative parents for sample 2");
    ok(Probs._num_triples[3] == 2, "putative parents for sample 3");
    ok(Probs._num_triples[4] == 5, "putative parents for sample 4");

    PROBcalcPosteriors();
    //PROBdumpPosteriors(stderr, Probs.posteriors, NULL, NULL, 1);
    /* some macro tests */

    ok(CMP_FLOAT(0.000001, 0.000001) == 1, "equal");
    ok(CMP_FLOAT(0.000002, 0.000001) == 0, "not equal");
    ok(CMP_DBL(0.000001, 0.000001) == 1, "equal");
    ok(CMP_DBL(0.000002, 0.000001) == 0, "not equal");

    child[0] = 10;
    child[1] = 20;
    p1[0] = p1[1] = 10;  
    p2[0] = p2[1] = 20;

    ok(CMP_DBL(LODcalcTransProbTriple(child,p1,p2,0), 1.), "T(child|p1,p2) correct");
    ok(CMP_DBL(LODcalcTransProbTriple(child,p2,p1,0), 1.), "T(child|p1,p2) correct");

    p1[0] = p1[1] = p2[0] = p2[1] = 10;

    ok(CMP_DBL(LODcalcTransProbTriple(child,p1,p2,0), 0.), "T(child|p1,p2) correct");
    ok(CMP_DBL(LODcalcTransProbTriple(child,p2,p1,0), 0.), "T(child|p1,p2) correct");

    p1[0] = p1[1] = 10;  
    p2[0] = 20;
    p2[1] = 30;
    ok(CMP_DBL(LODcalcTransProbTriple(child,p1,p2,0), 0.5), "T(child|p1,p2) correct");
    ok(CMP_DBL(LODcalcTransProbTriple(child,p2,p1,0), 0.5), "T(child|p1,p2) correct");

    p1[0] = 5;
    p1[1] = 10;
    ok(CMP_DBL(LODcalcTransProbTriple(child,p1,p2,0), 0.25), "T(child|p1,p2) correct");
    ok(CMP_DBL(LODcalcTransProbTriple(child,p2,p1,0), 0.25), "T(child|p1,p2) correct");
    p1[0] = 10;
    p1[1] = 20;
    p2[0] = 10;
    p2[1] = 20;
    ok(CMP_DBL(LODcalcTransProbTriple(child,p1,p2,0), 0.5), "T(child|p1,p2) correct");
    ok(CMP_DBL(LODcalcTransProbTriple(child,p2,p1,0), 0.5), "T(child|p1,p2) correct");
    p1[0] = 10;
    p1[1] = 20;
    p2[0] = 20;
    p2[1] = 20;
    ok(CMP_DBL(LODcalcTransProbTriple(child,p1,p2,0), 0.5), "T(child|p1,p2) correct");
    ok(CMP_DBL(LODcalcTransProbTriple(child,p2,p1,0), 0.5), "T(child|p1,p2) correct");
    p1[0] = 10;
    p1[1] = 30;
    p2[0] = 10;
    p2[1] = 20;
    ok(CMP_DBL(LODcalcTransProbTriple(child,p1,p2,0), 0.25), "T(child|p1,p2) correct");
    ok(CMP_DBL(LODcalcTransProbTriple(child,p2,p1,0), 0.25), "T(child|p1,p2) correct");

  /*  PROBdumpPosteriors(stderr, Probs.posteriors, NULL, NULL, 1);
    fprintf(stderr, "Nf: %f Nm %f\n", Options.Nf, Options.Nm);*/
    for (i=0; i<Data.num_samples; i++) {
        p = pg = 0.;
        for (j=0; j<Probs.num_posteriors[i]; j++) {
            p+=exp(Probs.posteriors[i][j].p_opt);
            pg += exp(Probs.posteriors[i][j].lw);
        }

        ok(p > 0.4 && p <= 1.0, "Posteriors sum to one");
       ptmp =  pg/(pg+Probs._sum_filtered_likelihood_dyads_f[i] +  Probs._sum_filtered_likelihood_triples[i]+ Probs._sum_filtered_likelihood_dyads_m[i]);
        fprintf(stderr, "%f %i p=%E  p2=%E %E %E %E\n",ptmp,i,p,pg, Probs._sum_filtered_likelihood_triples[i],
                Probs._sum_filtered_likelihood_dyads_m[i], Probs._sum_filtered_likelihood_dyads_f[i]
                );
        ok(CMP_DBL(ptmp,p),"unfiltered correct");
    }    

    /*
        locus: 0 , p 0.006250, pl 0.000000
        locus: 1 , p 0.468750, pl -5.075174
        locus: 2 , p 0.006250, pl -5.832860
        locus: 3 , p 0.075000, pl -10.908033
        locus: 4 , p 0.012500, pl -13.498301
        locus: 5 , p 0.012500, pl -17.880327
        locus: 6 , p 0.225000, pl -22.262354
    */
    pl[0] = FREQ_ALLELE_PROB(0,147);
    pl[8] = FREQ_ALLELE_PROB(0,151);
    pl[9] = FREQ_ALLELE_PROB(0,157);
    pl[1] = FREQ_ALLELE_PROB(1,182);
    pl[7] = FREQ_ALLELE_PROB(1,188);
    pl[2] = FREQ_ALLELE_PROB(2,218);
    pl[3] = FREQ_ALLELE_PROB(3,226);
    pl[4] = FREQ_ALLELE_PROB(4,165);
    pl[5] = FREQ_ALLELE_PROB(5,168);
    pl[6] = FREQ_ALLELE_PROB(6,178);
    pl[10] = FREQ_ALLELE_PROB(2,212);
    pl[11] = FREQ_ALLELE_PROB(2,214);
    pl[12] = FREQ_ALLELE_PROB(2,220);
    ok(fabs(pl[0] - 0.012500) < EPSILON, "Allele Freq. 147 (locus 0)");
    ok(fabs(pl[1] - 0.512500) < EPSILON, "Allele Freq. 182 (locus 1)");
    ok(fabs(pl[7] - 0.425000) < EPSILON, "Allele Freq. 188 (locus 1)");
    ok(fabs(pl[2] - 0.012500) < EPSILON, "Allele Freq. 218 (locus 2)");
    ok(fabs(pl[3] - 0.075000) < EPSILON, "Allele Freq. 226 (locus 3)");
    ok(fabs(pl[4] - 0.012500) < EPSILON, "Allele Freq. 165 (locus 4)");
    ok(fabs(pl[5] - 0.012500) < EPSILON, "Allele Freq. 168 (locus 5)");
    ok(fabs(pl[6] - 0.450000) < EPSILON, "Allele Freq. 178 (locus 6)");
    child[0] = 147;
    child[1] = 147;
    /* check IBD */
    IBDcalcGenotypeProbabilities(child,child, &m0, &m1, &m2, 0);
    ok(fabs(m0 - pow(pl[0],4))< EPSILON, "P(0) correct");
    ok(fabs(m1 - pow(pl[0],3))< EPSILON, "P(1) correct");
    ok(fabs(m2 - (pl[0] * pl[0])) < EPSILON, "P(2) correct");
    p1[0] = 147;
    p1[1] = 151;
    IBDcalcGenotypeProbabilities(child,p1, &m0, &m1, &m2, 0);
    ok(fabs(m0 - 2*pow(pl[0],3)*pl[8])< EPSILON, "P(0) correct");
    ok(fabs(m1 - pl[0]*pl[0]*pl[8])< EPSILON, "P(1) correct");
    ok(fabs(m2 - 0.) < EPSILON, "P(2) correct");
    IBDcalcGenotypeProbabilities(p1,child, &m0, &m1, &m2, 0);
    ok(fabs(m0 - 2*pow(pl[0],3)*pl[8])< EPSILON, "P(0) correct");
    ok(fabs(m1 - pl[0]*pl[0]*pl[8])< EPSILON, "P(1) correct");
    ok(fabs(m2 - 0.) < EPSILON, "P(2) correct");
    /* aiai ajaj */
    p1[0] = 151;
    p1[1] = 151;
    IBDcalcGenotypeProbabilities(child,p1, &m0, &m1, &m2, 0);
    ok(fabs(m0 - pl[0]*pl[0]*pl[8]*pl[8]) < EPSILON, "P(0) correct");
    ok(fabs(m1 - 0.)< EPSILON, "P(1) correct");
    ok(fabs(m2 - 0.) < EPSILON, "P(2) correct");
    IBDcalcGenotypeProbabilities(p1,child, &m0, &m1, &m2, 0);
    ok(fabs(m0 - pl[0]*pl[0]*pl[8]*pl[8]) < EPSILON, "P(0) correct");
    ok(fabs(m1 - 0.)< EPSILON, "P(1) correct");
    ok(fabs(m2 - 0.) < EPSILON, "P(2) correct");
    p1[0] = 147;
    p1[1] = 151;
    /*aiaj aiaj */
    IBDcalcGenotypeProbabilities(p1,p1, &m0, &m1, &m2, 0);
    ok(fabs(m0 - 4*pl[0]*pl[0]*pl[8]*pl[8]) < EPSILON, "P(0) correct");
    ok(fabs(m1 - pl[0]*pl[8]*(pl[0]+pl[8]))< EPSILON, "P(1) correct");
    ok(fabs(m2 - 2*pl[0]*pl[8]) < EPSILON, "P(2) correct");
    /*aiaj ajal */
    child[0] = 147;
    child[1] = 151;
    p1[0] = 151;
    p1[1] = 157;
    IBDcalcGenotypeProbabilities(child,p1, &m0, &m1, &m2, 0);
    ok(fabs(m0 - 4*pl[0]*pl[8]*pl[8]*pl[9]) < EPSILON, "P(0) correct");
    ok(fabs(m1 - pl[0]*pl[8]*pl[9])< EPSILON, "P(1) correct");
    ok(fabs(m2 - 0) < EPSILON, "P(2) correct");
    IBDcalcGenotypeProbabilities(p1,child, &m0, &m1, &m2, 0);
    ok(fabs(m0 - 4*pl[0]*pl[8]*pl[8]*pl[9]) < EPSILON, "P(0) correct");
    ok(fabs(m1 - pl[0]*pl[8]*pl[9])< EPSILON, "P(1) correct");
    ok(fabs(m2 - 0) < EPSILON, "P(2) correct");
    /* aiaj alal */
    p1[0] = 157;
    p1[1] = 157;
    IBDcalcGenotypeProbabilities(p1,child, &m0, &m1, &m2, 0);
    ok(fabs(m0 - 2*pl[0]*pl[8]*pl[9]*pl[9]) < EPSILON, "P(0) correct");
    ok(fabs(m1 - 0.)< EPSILON, "P(1) correct");
    ok(fabs(m2 - 0.) < EPSILON, "P(2) correct");
    IBDcalcGenotypeProbabilities(child, p1, &m0, &m1, &m2, 0);
    ok(fabs(m0 - 2*pl[0]*pl[8]*pl[9]*pl[9]) < EPSILON, "P(0) correct");
    ok(fabs(m1 - 0.)< EPSILON, "P(1) correct");
    ok(fabs(m2 - 0.) < EPSILON, "P(2) correct");
    child[0] = 151;
    child[1] = 147;
    IBDcalcGenotypeProbabilities(child, p1, &m0, &m1, &m2, 0);
    ok(fabs(m0 - 2*pl[0]*pl[8]*pl[9]*pl[9]) < EPSILON, "P(0) correct");
    ok(fabs(m1 - 0.)< EPSILON, "P(1) correct");
    ok(fabs(m2 - 0.) < EPSILON, "P(2) correct");
    IBDcalcGenotypeProbabilities(p1, child, &m0, &m1, &m2, 0);
    ok(fabs(m0 - 2*pl[0]*pl[8]*pl[9]*pl[9]) < EPSILON, "P(0) correct");
    ok(fabs(m1 - 0.)< EPSILON, "P(1) correct");
    ok(fabs(m2 - 0.) < EPSILON, "P(2) correct");

    /*aiaj akal*/
    child[0] = 212;
    child[1] = 218;
    p1[0] = 214;
    p1[1] = 220;
    IBDcalcGenotypeProbabilities(child, p1, &m0, &m1, &m2, 2);
    ok(fabs(m0 - 4*pl[2]*pl[10]*pl[11]*pl[12]) < EPSILON, "P(0) correct");
    ok(fabs(m1 - 0.)< EPSILON, "P(1) correct");
    ok(fabs(m2 - 0.) < EPSILON, "P(2) correct");

    /* now check Trans. Probabilities for Dyad (see Meagher&Thompson,1986,p.
     * 98)*/
    child[0] = 147;
    child[1] = 147;
    p =LODcalcTransProbDyad(child, child, 0);
    ok(fabs(p - FREQ_ALLELE_PROB(0,147)) < EPSILON,
      "T(a1a1|a1a1,?) correct");
    p1[0] = 147;
    p1[1] = 151;
    p =LODcalcTransProbDyad(child, p1, 0);
    ok(fabs(p - 0.5*FREQ_ALLELE_PROB(0,147)) < EPSILON,
      "T(a1a1|a1a2,?) correct");
    p =LODcalcTransProbDyad(p1,child, 0);
    ok(fabs(p - FREQ_ALLELE_PROB(0,151)) < EPSILON,
      "T(a1a2|a1a1,?) correct");
    p =LODcalcTransProbDyad(p1,p1, 0);
    ok(fabs(p - (0.5 * (FREQ_ALLELE_PROB(0,147)+FREQ_ALLELE_PROB(0,151)))) < EPSILON,
      "T(a1a2|a1a2,?) correct");
    /* genotype probabilities */  
    p=LODcalcGenotypeProb(p1,0);
    pg = 2.*FREQ_ALLELE_PROB(0,147) * FREQ_ALLELE_PROB(0,151);
    ok(fabs(p - pg) < EPSILON,"P(a1,a2) correct");
    child[0] = 147;
    child[1] = 147;
    p=LODcalcGenotypeProb(child,0);
    pg = FREQ_ALLELE_PROB(0,147) * FREQ_ALLELE_PROB(0,147);
    ok(fabs(p - pg) < EPSILON,"P(a1,a1) correct");

    /* missing alleles */
    child[0] = 147;
    child[1] = 147;
    p1[0] = 151;
    p1[1] = -1;
    p=LODcalcGenotypeProb(p1,0);
    ok(CMP_DBL(p,FREQ_ALLELE_PROB(0,151)), "P(a,?) correct");
    p1[0] = -1;
    p1[1] = 151;
    ok(CMP_DBL(p,FREQ_ALLELE_PROB(0,151)), "P(?,a) correct");

    p =LODcalcTransProbDyad(child,p1, 0);
    pg = 0.5 * 0. + 0.5 * FREQ_ALLELE_PROB(0,147) * FREQ_ALLELE_PROB(0,147);
    ok(CMP_DBL(p,pg), "T(a1.a1|a2.?,?) correct");
    p1[0] = 147;
    p1[1] = -1;
    p =LODcalcTransProbDyad(child,p1, 0);
    pg = 0.5 * FREQ_ALLELE_PROB(0,147) + 0.5 * FREQ_ALLELE_PROB(0,147) * FREQ_ALLELE_PROB(0,147);
    ok(CMP_DBL(p,pg), "T(a1.a1|a1.?,?) correct");

    child[0] = 151;
    child[1] = 147;
    p1[0] = 151;
    p1[1] = -1;
    p =LODcalcTransProbDyad(child,p1, 0);
    pg = 0.5 * FREQ_ALLELE_PROB(0,147) + 0.5 * 2 * FREQ_ALLELE_PROB(0,147) * FREQ_ALLELE_PROB(0,151);
    ok(CMP_DBL(p,pg), "T(a1.a2|a1.?,?) correct");

    p1[0] = -1;
    p1[1] = 151;
    p =LODcalcTransProbDyad(child,p1, 0);
    ok(CMP_DBL(p,pg), "T(a1.a2|?.a1,?) correct");

    child[0] = 147;
    child[1] = -1;
    p1[0] = p1[1] = 147;
    p =LODcalcTransProbDyad(child,p1, 0);
    pg = 0.5*FREQ_ALLELE_PROB(0,147) + 0.25 * (1.+1.);
    ok(CMP_DBL(p,pg), "T(a1.?|a1.a1,?) correct");

    p1[0] = p1[1] = 151;
    p =LODcalcTransProbDyad(child,p1, 0);
    pg = 0.5*FREQ_ALLELE_PROB(0,147) + 0.25 * (0.+0.);
    ok(CMP_DBL(p,pg), "T(a1.?|a2.a2,?) correct");

    p1[0] = 147;
    p =LODcalcTransProbDyad(child,p1, 0);
    pg = 0.5*FREQ_ALLELE_PROB(0,147) + 0.25 * (1.+0.);
    ok(CMP_DBL(p,pg), "T(a1.?|a1.a2,?) correct");

    p1[1] = -1;
    p =LODcalcTransProbDyad(child,p1, 0);
    pg = 0.5*FREQ_ALLELE_PROB(0,147) + 0.25 * (1.+ FREQ_ALLELE_PROB(0,147));
    ok(CMP_DBL(p,pg), "T(a1.?|a1.?,?) correct");

    p1[0] = 151;
    p =LODcalcTransProbDyad(child,p1, 0);
    pg = 0.5*FREQ_ALLELE_PROB(0,147) + 0.25 * (0.+ FREQ_ALLELE_PROB(0,147));
    ok(CMP_DBL(p,pg), "T(a1.?|a2.?,?) correct");

    /**************************/
    /* missing alleles triple */
    /**************************/

    /* missing child allele */
    child[0] = 151;
    child[1] = -1;
    p1[0] = p1[1] = 151;
    p2[0] = p2[1] = 151;
    p =LODcalcTransProbTriple(child,p1, p2, 0);
    ok(CMP_DBL(p,1.), "T(a1.?|a1.a1,a1.a1) correct");

    p1[0] = p1[1] = 151;
    p2[0] = p2[1] = 147;
    p =LODcalcTransProbTriple(child,p1, p2, 0);
    ok(CMP_DBL(p,0.5), "T(a1.?|a1.a1,a2.a2) correct");
    child[0] = -1;
    child[1] = 151;
    p =LODcalcTransProbTriple(child,p1, p2, 0);
    ok(CMP_DBL(p,0.5), "T(a1.?|a1.a1,a2.a2) correct");

    p1[0] = p1[1] = 147;
    p2[0] = p2[1] = 147;
    p =LODcalcTransProbTriple(child,p1, p2, 0);
    ok(CMP_DBL(p,0.), "T(a1.?|a2.a2,a2.a2) correct");

    child[0] = 151;
    child[1] = -1;
    p1[0] = p1[1] = 147;
    p2[0] = 151;
    p2[1] = 147;
    p =LODcalcTransProbTriple(child,p1, p2, 0);
    ok(CMP_DBL(p,0.25), "T(a1.?|a2.a2,a1.a2) correct");
    child[0] = -1;
    child[1] = 147;
    p =LODcalcTransProbTriple(child,p1, p2, 0);
    ok(CMP_DBL(p,0.75), "T(?.147|a2.a2,a1.a2) correct");
    p2[1] = -1;
    p =LODcalcTransProbTriple(child,p1, p2, 0);
    ok(CMP_DBL(p,0.5+0.25*FREQ_ALLELE_PROB(0,147)+0.), "T(?.a2|a2.a2,a1.?) correct");
    p1[0] = -1;
    p1[1] = 147;
    p2[0] = 151;
    p2[1] = -1;
    p =LODcalcTransProbTriple(child,p1, p2, 0);
    ok(CMP_DBL(p,0.25+2*0.25*FREQ_ALLELE_PROB(0,147)+0.), "T(?.a2|a?.a2,a1.?) correct");

    child[0] = child[1] = 147;
    p1[0] = p1[1] = 147;
    p2[0] = -1;
    p2[1] = 147;

    p =LODcalcTransProbTriple(child,p1, p2, 0);
    /* 100% from p1, 50% p2 + 50 % missing */
    pg = 1. * 0.5 * ( 1 + FREQ_ALLELE_PROB(0, 147) );
    ok(CMP_DBL(p,pg), "T(a1.a1|a1.a1,a1.?) correct");
    p2[1] = -1;
    p2[0] = 147;
    p =LODcalcTransProbTriple(child,p1, p2, 0);
    ok(CMP_DBL(p,pg), "T(a1.a1|a1.a1,?.a1) correct");

    p2[0] = p2[1] = 147;
    p1[1] = -1;
    p1[0] = 147;
    p =LODcalcTransProbTriple(child,p1, p2, 0);
    ok(CMP_DBL(p,pg), "T(a1.a1|?.a1,a1.a1) correct");
    p1[0] = -1;
    p1[1] = 147;
    p =LODcalcTransProbTriple(child,p1, p2, 0);
    ok(CMP_DBL(p,pg), "T(a1.a1|?.a1,a1.a1) correct");

    child[0] = 151;
    child[1] = 147;
    p1[0] = p1[1] = 147;
    p2[0] = -1;
    p2[1] = 147;
    p =LODcalcTransProbTriple(child,p1, p2, 0);
    /* 100% from p1, 50% p2 + 50 % missing */
    pg = 0.5 * (1. + 1.) * 0.5 * ( 0 + FREQ_ALLELE_PROB(0, 151) );
    ok(CMP_DBL(p,pg), "T(a1.a2|a2.a2,a2.?) correct");
    p2[1] = -1;
    p2[0] = 147;
    ok(CMP_DBL(p,pg), "T(a1.a2|a2.a2,?.a2) correct");

    child[0] = child[1] = 147;
    p1[0] = -1;
    p1[1] = 147;
    p2[0] = -1;
    p2[1] = 147;
    p =LODcalcTransProbTriple(child,p1, p2, 0);
    /* 100% from p1, 50% p2 + 50 % missing */
    pg = 0.5 * (1. + FREQ_ALLELE_PROB(0,147)) * 0.5 * ( 1. + FREQ_ALLELE_PROB(0, 147) );
    ok(CMP_DBL(p,pg), "T(a1.a1|?.a1,?.a1) correct");
    p2[0] = 147;
    p2[1] = -1;
    p =LODcalcTransProbTriple(child,p1, p2, 0);
    ok(CMP_DBL(p,pg), "T(a1.a1|?.a1,a1.?) correct");
    p1[0] = 147;
    p1[1] = -1;
    p =LODcalcTransProbTriple(child,p1, p2, 0);
    ok(CMP_DBL(p,pg), "T(a1.a1|a1.?,a1.?) correct");

    child[0] = 151;
    child[1] = 147;
    p1[0] = -1;
    p1[1] = 147;
    p2[0] = -1;
    p2[1] = 147;
    p =LODcalcTransProbTriple(child,p1, p2, 0);
    /* 100% from p1, 50% p2 + 50 % missing */
    pg =  2*(0.5 * FREQ_ALLELE_PROB(0,151) * 0.5) +(0.5 * FREQ_ALLELE_PROB(0,147)*FREQ_ALLELE_PROB(0,151));
    //fprintf(stderr, "%f %f %f %f\n", p,pg,FREQ_ALLELE_PROB(0,147), FREQ_ALLELE_PROB(0,151));
    ok(CMP_DBL(p,pg), "T(a1.a2|?.a1,?.a1) correct");
    p1[1] = -1;
    p1[0] = 147;
    p2[0] = -1;
    p2[1] = 147;
    p =LODcalcTransProbTriple(child,p1, p2, 0);
    ok(CMP_DBL(p,pg), "T(a1.a2|a1.?,?.a1) correct");
    p2[1] = -1;
    p2[0] = 147;
    p =LODcalcTransProbTriple(child,p1, p2, 0);
    ok(CMP_DBL(p,pg), "T(a1.a2|a1.?,a1.?) correct");

    /* missing parent allele */

    /*p = log(0.5*pl[0]) + log(0.5*(pl[1]+pl[7])) + log(0.5*pl[2]) + 
        log(pl[3]) + log(pl[4]) + log(pl[5]) + log(0.5*pl[6]);
    p += log(Data.id_mapping[4].ramets);
    ok(abs(Probs.dyads[0]->p - p) < 0.00001, "T(g_0|g_4) correct");*/
    
    p = LODcalcDenominator(Data.id_mapping[0].genotype_obs, NULL);
    //fprintf(stderr, "%f %f", p, EPSILON);

    ok(fabs(p + 20.651171)< EPSILON, "LOD denom correct");

    ok(Data.allele_frequencies[0].allele_id_lookup[147-147] == 0, "id_lookup correct");
    ok(Data.allele_frequencies[0].allele_id_lookup[148-147] == -1, "id_lookup correct");
    ok(Data.allele_frequencies[0].allele_id_lookup[151-147] == 1, "id_lookup correct");
    ok(Data.allele_frequencies[0].allele_id_lookup[157-147] == 2, "id_lookup correct");

    ok(FREQ_ALLELE_ID(0,148) == -1, "id_lookup correct");
    ok(FREQ_ALLELE_ID(0,157) == 2, "id_lookup correct");

    PROBdestroyPreMCMC();
    DATAIOdestroyPreMCMC();
    FREQdestroyPreMCMC();
    PROBdestroyPostMCMC();
    DATAIOdestroyPostMCMC();
    FREQdestroyPostMCMC();
    Options.MinTyped = 2;
    Options.MaxMismatchingDefined = 1;
    Options.MaxMismatchingTriple = 0;
    Options.MaxMismatchingDyad = 0;
    Options.Monoecious = false;
    DATAIOinfile(filenameMD);
    FREQcalcAlleleFrequencies();
    FREQcalcSummaryStatistics();
    PROBinit();
    PROBcalcTriplesAndDyads();
    PROBcalcPosteriors();
    PROBrecalcPosteriorProbabilities(Probs.posteriors, 10., 5.);
    i = 2;
    ok(LODcalcMismatches(i,-1,-1,NULL) == 0, "no mismatches Bart");
    ok(LODcalcMismatches(i,1,-1,NULL) == 0, "no mismatches Bart-Homer");
    ok(LODcalcMismatches(i,1,5,NULL) == 0, "no mismatches Bart-Homer-Marge");
    ok(LODcalcMismatches(i,3,-1,NULL) == 1, "one mismatches Bart-Lisa");
    ok(LODcalcMismatches(i,3,6,NULL) == 2, "two mismatches Bart-Lisa-Flanders");
    /* unsampled unsampled bart */
    j= PROBsearchPosteriors(Probs.posteriors,i,-1,-1);
    ok(CMP_DBL(exp(Probs.posteriors[i][j].lw-Probs.posteriors[i][j].l),( (10. - 2.) * ( 5. - 2.) )), "8x3 unsampled parents");
    /* homer marge bart */
    j= PROBsearchPosteriors(Probs.posteriors,i,1,5);
    ok(CMP_DBL(exp(Probs.posteriors[i][j].lw-Probs.posteriors[i][j].l),1), "1 x homer marge parents");
    /* homer marge bart */
    j= PROBsearchPosteriors(Probs.posteriors,i,0,-1);
    ok(CMP_DBL(exp(Probs.posteriors[i][j].lw-Probs.posteriors[i][j].l),(10. - 2.)), "8 x unsampled mothers");

    /* homer marge bart */
    j= PROBsearchPosteriors(Probs.posteriors,i,5,-1);
    ok(CMP_DBL(exp(Probs.posteriors[i][j].lw-Probs.posteriors[i][j].l),(5. - 2.)), "3 x unsampled fathers");
    PROBrecalcPosteriorProbabilities(Probs.posteriors, 1., 1.);

    /********* alleles lookup things working? ***********/
    ok(Data.num_alleles[2] == 4, "num alleles ok");

    ok(Data.num_alleles[0] == 6, "num alleles at locus 1 correct");
    ok(Data.max_num_alleles == 6, "max num alleles correct");

    ok(Data.alleles[0][0] == 100, "allele 1 at locus 1 correct");
    ok(Data.alleles[0][1] == 110, "allele 2 at locus 1 correct");
    ok(Data.alleles[0][2] == 120, "allele 3 at locus 1 correct");
    ok(Data.alleles[0][3] == 140, "allele 4 at locus 1 correct");
    ok(Data.alleles[0][4] == 160, "allele 5 at locus 1 correct");
    ok(Data.alleles[0][5] == 170, "allele 6 at locus 1 correct");

    for (i=0; i<50000; i++) {
        j = GENOTYPErandAlleleUnif(2, 1);
        if (j<0)
            randres[0]++;
        else 
            randres[j-Data.allele_frequencies[2].min+1]++;
    }
    ok(randres[0] < 13000 && randres[1] < 13000 && randres[3] < 13000 && randres[5] < 13000, "randalleleunif");
    ok(randres[0] > 8000 && randres[1] > 8000 && randres[3] > 8000 && randres[5] > 8000, "randalleleunif");
    for (i=0; i<8; i++) {
/*        fprintf(stderr, "randres unif %i: %i\n", i, randres[i]);*/
        randres[i] = 0;
    }
    
    /*FREQdumpSummaryStatistics(stderr);*/
    for (i=0; i<50000; i++) {
        j = GENOTYPErandAllele(2);
        randres[j-Data.allele_frequencies[2].min]++;
    }    
    ok(randres[0] > 10000 && randres[2] > 15000 && randres[4] > 10000 && randres[6] < 15000, "randallele");
    for (i=0; i<8; i++) {
    /*    fprintf(stderr, "randres %i: %i\n", i, randres[i]);*/
        randres[i] = 0;
    }        

    /********* are_connected matrix correct? ***********/

    ok(Probs.are_connected[TIDXL(0,1)] == true,  "grampa homer");
    ok(Probs.are_connected[TIDXL(0,2)] == true,  "grampa bart");
    ok(Probs.are_connected[TIDXL(0,3)] == false, "grampa not involved");
    ok(Probs.are_connected[TIDXL(0,4)] == false, "grampa not involved");
    ok(Probs.are_connected[TIDXL(0,5)] == false, "grampa not involved");
    ok(Probs.are_connected[TIDXL(0,6)] == false, "grampa not involved");
    for (i=0; i<5; i++) 
        ok(Probs.are_connected[TIDXL(1,i)] == true, "homer involved");
    ok(Probs.are_connected[TIDXL(1,5)] == false, "homer involved");
    ok(Probs.are_connected[TIDXL(1,6)] == false, "homer involved");
    for (i=0; i<6; i++) 
        ok(Probs.are_connected[TIDXL(6,i)] == false, "flanders not involved");
    ok(Probs.are_connected[TIDXL(6,6)] == true, "flanders involved");

    for (i=0; i<Data.num_loci; i++) {
        p = 0.;
        for (j=0; j<=(Data.allele_frequencies[i].max-Data.allele_frequencies[i].min); j++) 
            p += FREQ_ALLELE_PROB(i,j+Data.allele_frequencies[i].min);   
        ok(fabs(p - 1.0) < EPSILON, "Allele Freqs sum to one");
    }    
    mcmc = MCMCinit();
    state = MCMCinitState(&mcmc, NULL);
    Options.UpdateAlleleFreqs = 0;
    Options.Nf = 1;
    Options.Nm = 1;
    Options.Nfdefined  = true;
    Options.Nmdefined  = true;
    PROBrecalcPosteriorProbabilities(Probs.posteriors, Options.Nf, Options.Nm);
    //PROBdumpPosteriors(stderr, Probs.posteriors, NULL, NULL, 0);
    PEDIGREEestmateSamplingRateSex(state.D, &sr_f, &sr_m);
    ok(CMP_DBL(sr_f, D_MINSAMPLINGRATE), "sampling rate f min with empty pedigree");
    ok(CMP_DBL(sr_m, D_MINSAMPLINGRATE), "sampling rate m min with empty pedigree");
    (void)MCMCgreedy(&mcmc, &state);
    DATAIOaddGVNodeStyle(state.D);
    DAGinsertE(state.D, EDGE(0,1));
    /*DAGdump(stderr, state.D);
    FREQdump(stderr);*/
    PEDIGREEestmateSamplingRateSex(state.D, &sr_f, &sr_m);
    //fprintf(stderr,"DEBUG %f %f\n",sr_f,sr_m);
    ok(CMP_DBL(sr_f, 3./4.), "sampling rate f min correct");
    ok(CMP_DBL(sr_m, 4./4.), "sampling rate m min correct");
    //fprintf(stderr,"%f %f\n", sr_f, sr_m);
    FREQupdate(state.D);
    /*FREQdump(stderr);
    PROBdumpPosteriors(stderr, Probs.posteriors, state.D, NULL, 1);*/
    PROBrecalcPosteriorProbabilities(Probs.posteriors, Options.Nf, Options.Nm);
    lc.id = 0;
    lc.locus = 2;
    state.genotypes[lc.id][lc.locus][0] = 306;
    state.genotypes[lc.id][lc.locus][1] = 306;
    
    ok(fabs(Data.allele_frequencies[0].freqs[100-100] - 1./7.) < EPSILON, "freq after updating correct");
    ok(fabs(Data.allele_frequencies[0].freqs[110-100] - 1./7.) < EPSILON, "freq after updating correct");
    ok(fabs(Data.allele_frequencies[0].freqs[120-100] - 1./7.) < EPSILON, "freq after updating correct");
    ok(fabs(Data.allele_frequencies[0].freqs[140-100] - 1./7.) < EPSILON, "freq after updating correct");
    ok(fabs(Data.allele_frequencies[0].freqs[160-100] - 2./7.) < EPSILON, "freq after updating correct");
    ok(fabs(Data.allele_frequencies[0].freqs[170-100] - 1./7.) < EPSILON, "freq after updating correct");

    MCMCdestroyState(&mcmc,state);
    MCMCdestroy(mcmc);

    /*ok(FREQ_ALLELE_PROB(i,-1) == 1.0,"Probability of missing allele is 1");   */
    PROBdestroyPreMCMC();
    DATAIOdestroyPreMCMC();
    FREQdestroyPreMCMC();
    PROBdestroyPostMCMC();
    DATAIOdestroyPostMCMC();
    FREQdestroyPostMCMC();
    DATAIOinfile(filenamePM);
    Options.FemRepro.min = 14;
    Options.FemRepro.max = 45;
    Options.MaleRepro.min = 14;
    Options.MaleRepro.max = 45;
    FREQcalcAlleleFrequencies();
    FREQcalcSummaryStatistics();
    PROBinit();
    sim = SIMstart();
    
    PROBcalcTriplesAndDyads();
    PROBcalcPosteriors();
    //PROBdumpPosteriors(stderr, Probs.posteriors, NULL, NULL,0);
    strcpy(Options.OffspringOutfilenameCERVUS,"tmp.cervus.1");
    DATAIOdumpOffspringCERVUS();
    ok(Probs.candidates[25][Probs.num_candidates[25]-1] == 27, "last cand.");
    ok(Probs.num_candidates[25] == 6, "6 cand.");
    ok(Probs.num_candidates[2] == 9, "6 cand.");

    Options.FemMaxDist = 1;
    Options.MaleMaxDist = 1;
    PROBdestroyPreMCMC();
    DATAIOdestroyPreMCMC();
    FREQdestroyPreMCMC();
    PROBdestroyPostMCMC();
    DATAIOdestroyPostMCMC();
    FREQdestroyPostMCMC();
    SIMdestroy(sim);
    DATAIOinfile(filenamePM);
    DATAIOgeofile(filenamePMDIST, 1);
    FREQcalcAlleleFrequencies();
    FREQcalcSummaryStatistics();
    PROBinit();
    sim = SIMstart();
    
    PROBcalcTriplesAndDyads();
    PROBcalcPosteriors();
    ok(Probs.num_candidates[2] == 5, "5 cand.");
    strcpy(Options.OffspringOutfilenameCERVUS,"tmp.cervus.2");
    DATAIOdumpOffspringCERVUS();
    PROBcalcTriplesAndDyads();
    PROBcalcPosteriors();


//    remove(Options.OffspringOutfilenameCERVUS);

    remove(Options.SiblingsOutfilename[1]);
    remove(Options.LociOutfilename);
    remove("tmp.cervus.1");
    remove("tmp.cervus.2");
    dmin[0] = 1.; dmax[0] = 0.;
    dmin[1] = 3.; dmax[1] = 0.;

    /* this is not a serious random number generator test, it is just here  *
     * to catch bugs in the macros (which normally show up very loudly). It *
     * actually tests if the generated numbers are within the range.        */
    for (i=0; i<10000; i++) {
        j = RANDINT(5.);
        randres[j]++;
        randdblres[0] = RANDDBLONE;
        randdblres[1] = RANDDBL(3.);
        for (j=0; j<2;j++) {
            if (randdblres[j] > dmax[j])
                dmax[j] = randdblres[j];

            if (randdblres[j] < dmin[j])
                dmin[j] = randdblres[j];
        }
        /*printf("%i:%f\n",i,randdblres[i]);*/
    } 
    
    ok(randres[5] == 0 , "RANDINT seems to work (very simple check, see prob_tests.c!)");
    for (i=0;i<5;i++)
        ok(randres[i] > 0, "RANDINT seems to work (very simple check, see prob_tests.c!)");

    ok(dmin[0] >= 0. && dmin[0] < 0.01 &&  dmax[0] > 0.99 && dmax[0] < 1., "RANDDBLONE seems to work (very simple check, see prob_tests.c!)");
    ok(dmin[1] >= 0. && dmin[1] < 0.01 &&  dmax[1] > 2.99 && dmax[1] < 3., "RANDDBLONE seems to work (very simple check, see prob_tests.c!)");

    /* use this files to examine the numbers in R for example */
#ifdef GENRANDFILES
    /* one should plot these numbers in R: */
    FOPENW(fprandint, "tmp.randint");
    FOPENW(fprandone, "tmp.randone");
    FOPENW(fprand100, "tmp.rand100");
    FOPENW(fpnorm, "tmp.norm");
    for (i=0; i< 100000; i++) {
        fprintf(fprandint,"%i\n", RANDINT(100));
        fprintf(fprand100,"%lf\n", RANDDBL(100));
        fprintf(fprandone,"%lf\n", RANDDBLONE);
        RANDNORM(1.,p);
        fprintf(fpnorm,"%lf\n", p );
    }    
    fclose(fprand100);
    fclose(fprandone);
    fclose(fprandint);
#endif

    SIMdestroy(sim);

    return exit_status();       /* Return the correct exit code */
}

/* 
  vim: ft=c sw=4 ts=4 expandtab
*/
