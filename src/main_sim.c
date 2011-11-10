/*
 * $Id: main_sim.c 1431 2009-04-10 22:05:46Z markus $
 *
 * Markus Riester, markus@bioinf.uni-leipzig.de                  */

#define _POSIX_SOURCE
#include "macros.h"             /* some important housekeeping macros   */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "global.h"             /* constants and structs                */

/* load the important algorithms                                 */

#include "dag.h"                /* the graph library                    */
#include "prob.h"               /* calculate all kinds of probabilities */
#include "freq.h"               /* calculate allele frequencies         */
#include "sim.h"                /* human pedigree data                  */
#include "genotype.h"
#include "utils.h"

/* load boring utils                                             */

#include "dataio.h"             /* read and write the data              */
#include "options.h"            /* init/load/store options              */


DATA    Data;
PROBS   Probs;
OPTIONS Options;

int     SimOffspring;

#define MAX_BIRTH_CNT 1
#define NUMSTART 1000
#define NUMIND 20000
#define STARTYEAR 1979
#define NULLLOCUS -1
#define CLONALRATE 0.8
#define SELFINGRATE 0.1
#define MIGRATIONRATE 0.01
#define ALLELE_IN "allelefreqs.dat"
#define PDG_OUT "simpedigree.dat"
#define MOTHERS_OUT "simmothers.dat"
#define PROB_DEATH 0.0
#define GTP_OUT "simgenotype.dat"
#define SIM_OUT "simsum.dat"
#define SIBLINGS_OUT "simsiblings.dat"
#define MAX_RAMETS 10000

typedef struct S_MARRIAGE
{
    int     partner_id;
    int     year_marriage;
    int     year_divorce;
} MARRIAGE;

typedef struct S_PARENTS
{
    int     mother_id;
    int     father_id;
} PARENTS;

static int year = STARTYEAR;
static int living = 0;
static int ramets = 0;
static int ind = 0;
static int born[NUMIND], died[NUMIND];
static MARRIAGE married[NUMIND];
static PARENTS parents[NUMIND];
static int selfed = 0;
static int notselfed = 0;
static int clonal_rep = 0;
static int sex_rep = 0;
static int kill_r = 0;
static int kill_g = 0;

static void
createRandomInd(int sex)
{
    GENOTYPErand(Data.samples[0][ind].genotype_obs);
    Data.samples[0][ind].sex = sex;
    Data.samples[0][ind].sex = 1;
    Data.samples[0][ind].ramets = 1;
    born[ind] = year;
    living++;
    ramets++;
    ind++;
}

static void
createClone(int id, int sex)
{
    //GENOTYPEclone(Data.samples[0][ind].genotype_obs,Data.samples[0][id].genotype_obs);
    Data.samples[0][id].ramets++;
    ramets++;
    fprintf(stdout, " producing clone of %i (genet has now %i ramets)\n", id, Data.samples[0][id].ramets);
}

static inline bool 
notLiving(int id) {
    return (born[id] < 0 || died[id] >= 0) ? true : false;
}

static void
_ramets_correct(void) {
    int i, sum =0;
    for (i = 0; i < ind; i++) {
        if (notLiving(i)) 
            continue;
        if (Data.samples[0][i].ramets == 0) FATALINT("_ramets_correct 1");
        sum += Data.samples[0][i].ramets;
    }
    if (ramets != sum) FATALINT("_ramets_correct");
}


static int
randomGenet(int id, bool newborn)
{
    int     i, j,*ids, r=0, s;

    MAKE1DINT(ids,ramets);
    _ramets_correct();
    if (ind < 1)
        FATALINT("randomGenet");

    for (i = 0; i < ind; i++) {
        if (notLiving(i))
            continue;
        if (!newborn && born[i] ==  year)
            continue;
        if (i == id)
            continue;
        for (j=0; j< Data.samples[0][i].ramets; j++)
            ids[r++] = i;
    }
    s = ids[RANDINT((double)r)];
    FREE1D(ids);
    //fprintf(stderr, "randomGenet: %d, %d\n",r,s);
    return s;
}

/*make twins?
int
getRandomBirthCount()
{
    int i;
    double p, p_cum = 0.;
    p = RANDDBLONE;
    for (i=0; i<MAX_BIRTH_CNT; i++) {
        p_cum += PROB_CHILDREN_BIRTH[i];
        if (p < p_cum)
            return i+1;
    }    
    return 1;
}
*/

static void
createOffspring(int mother_id, int father_id, int mb)
{
    int     i;
    if (mb <= 0)
        mb = 1;
    //mb = getRandomBirthCount();

    if (mb > 1)
        fprintf(stdout, "MULTIPLE BIRTH: %i\n", mb);
    for (i = 0; i < mb; i++) {
        if (mother_id == father_id) {
            GENOTYPEselfing(Data.samples[0][ind].genotype_obs,
                            Data.samples[0][mother_id].genotype_obs);
            selfed++;
        } else {
            GENOTYPEmate(Data.samples[0][ind].genotype_obs,
                         Data.samples[0][mother_id].genotype_obs,
                         Data.samples[0][father_id].genotype_obs);
            notselfed++;
        }
        fprintf(stdout, " producing new genet with id %i (mother: %i father: %i)\n", ind, mother_id, father_id);
        /* random sex */
        Data.samples[0][ind].sex = RANDINT(2.) + 1;
        Data.samples[0][ind].sex = 1;
        Data.samples[0][ind].ramets = 1;
        born[ind] = year;
        parents[ind].mother_id = mother_id;
        parents[ind].father_id = father_id;
        living++;
        ind++;
        ramets++;
        if (living > ramets) FATALINT("there cannot be more genets than ramets");
        /* in case of multiple births, cancel it when we already reached the
         * limit */
        if (ind == NUMIND)
            return;
    }
}

static void
killRamet(int i)
{
    if (Data.samples[0][i].ramets < 2) {
        died[i] = year;
        living--;
        kill_g++;
        fprintf(stdout, " killing genet with id %d (born: %i died: %i)\n",i, born[i], died[i]);
    } else 
        fprintf(stdout, " killing ramet with id %d (born: %i ramets: %d)\n",i, born[i], Data.samples[0][i].ramets);

    Data.samples[0][i].ramets--;
    ramets--;
    kill_r++;
}

int
main(int argc, char *argv[])
{
    int     i, ii, j, k, kr, n, age, indegree1, indegree2, *ids, old_genets, old_ramets, surv;
    double  p;
    FILE   *out = NULL, *siblings = NULL;
    Dag     D, D_mothers;
    GraphvizCluster dummy;
    Edge    incoming1[2], incoming2[2];
    
    MAKE1DINT(ids,NUMIND);
    assert(0);
    /* Load default and specified options and data */
    OPTIONSinit();
    DATAIOfreqfile(ALLELE_IN, 1);
#ifdef HAVE_OPENMP
    n = omp_get_max_threads();
#else
    n = 1;
#endif

    if (Options.Verbosity > 0)
        fprintf(stderr, "\nInitializing Mersenne Twister for %i thread(s) .",
                n);

    SRAND(Options.Seed, i, n);

    if (Options.Verbosity > 0)
        fprintf(stderr, " done.\n\n");

    Data.num_populations = 1;
    Data.num_samples = NUMIND;
    DATAIOinit();
    Data.num_samples_in_population[0] = NUMIND;
    MALLOC(Data.samples, SAMPLE *, Data.num_populations);
    MALLOC(Data.samples[0], SAMPLE, Data.num_samples_in_population[0]);
    FOPENW(out, SIM_OUT);
    OPTIONSdumpSim(out);
    for (i = 0; i < NUMIND; i++) {
        born[i] = -1;
        died[i] = -1;
        parents[i].mother_id = -1;
        parents[i].father_id = -1;
        married[i].partner_id = -1;
        married[i].year_marriage = -1;
        married[i].year_divorce = -1;
        MAKE2DINT(Data.samples[0][i].genotype_obs, Data.num_loci, 2, j);
        Data.samples[0][i].id = i;
        Data.samples[0][i].ramets = 1;
        sprintf(Data.samples[0][i].description, "SIM%i", i);
    }
    DATAIOcreateMapping();
    FREQinit();
    for (i = 0; i < Data.num_loci; i++)
        Data.TypingError[i] = 0.01;
    for (i = 0; i < Data.num_loci; i++)
        Options.ProportionTyped[i] = 1.0;

    for (i=0;i<NUMIND;i++) ids[i] = i;
    /*create random babies */
    for (i = 0; i < NUMSTART; i++)
        createRandomInd(1 + (i < (NUMSTART / 2.)));

    while (ind < NUMIND && living > 0) {
        k = 0; kr =0; ii=0;
        for (i = 0; i < NUMIND; i++) {
            if (born[i] < 0 || died[i] > 0)
                continue;
            if (Data.samples[0][i].ramets == 1) k++;
            if (Data.samples[0][i].ramets == 2) kr++;
            if (Data.samples[0][i].ramets == 3) ii++;
        }    
        fprintf(out, "YEAR: %i; genets = %i; ind. = %i; ramets = %i, single=%d, do=%d, t=%d  r=%f %d %d %d %d sum: %d\n",
                year, living, ind, ramets,k,kr,ii,kill_g/(double)kill_r,kill_g, kill_r, kill_g - old_genets, kill_r - old_ramets,
                (kill_g - old_genets)+(kill_r-old_ramets)
                );
        fprintf(stdout, "YEAR: %i; genets = %i; ind. = %i; ramets = %i, single=%d, do=%d, t=%d  r=%f %d %d %d %d sum: %d\n",
                year, living, ind, ramets,k,kr,ii,kill_g/(double)kill_r,kill_g, kill_r, kill_g - old_genets, kill_r - old_ramets,
                (kill_g - old_genets)+(kill_r-old_ramets)
                );
        old_genets = kill_g;
        old_ramets = kill_r;
        
        UTILSshuffleInt(ids, NUMIND);
        surv = 0;
        for (ii = 0; ii < NUMIND; ii++) {
            i = ids[ii];
            if (notLiving(i))
                continue;
            age = year - born[i];
            if (age < 1)
                continue;
            surv++;
            kr = Data.samples[0][i].ramets;
            fprintf(stdout, "%i AGE: %i  (year %i %i) kr %i ---  generated %i TOTAL %i --- genets: %i ramets: %i surv: %d\n", i, age,
                    year, born[i], kr, ind, NUMIND, living, ramets, surv);
            /*
            p = RANDDBLONE;
            if (p < PROB_DEATH) {
                continue;
            }*/
            for (k = 0; k < kr; k++) {
                /* make babies */
                p = RANDDBLONE;
                if (ind < NUMIND && p < CLONALRATE) {
                    createClone(i, 1);
                    clonal_rep++;
                    if (ramets > MAX_RAMETS) killRamet(randomGenet(-1, true));
                    continue;
                }
                if (ind < NUMIND) {
                    p = RANDDBLONE;
                    if (p < SELFINGRATE)
                        createOffspring(i, i, -1);
                    else {
                        p = RANDDBLONE;
                        if (p < MIGRATIONRATE) {
                            createRandomInd(1);
                            createOffspring(i, ind-1, -1);
                            if (ramets > MAX_RAMETS) killRamet(randomGenet(-1,true));
                        } else {
                            createOffspring(i, randomGenet(i,false), -1);
                        }    
                    }    
                    sex_rep++;
                    if (ramets > MAX_RAMETS) killRamet(randomGenet(-1,true));
                    continue;
                } else
                    break;
            }
            if (ind >= NUMIND)
                break;
        }
        /* all already dead? */
        if (living == 0 && ind + 1 < NUMIND)
            createRandomInd(RANDINT(2.) + 1);
        while (ramets > MAX_RAMETS) {
        fprintf(stdout, "KILLING IN THE YEAR: %i; living = %i; ind. = %i; ramets = %i\n",
                year, living, ind, ramets);
            killRamet(randomGenet(-1,true));
        }
        year++;
    }

    if (NULLLOCUS >= 0)
        for (i = 0; i < NUMIND; i++) {
            for (j = 0; j < 2; j++) {
                /*if (Data.allele_frequencies[NULLLOCUS].
                   allele_id_lookup[Data.samples[0][i].
                   genotype_obs[NULLLOCUS][j]] == 10) */
                if (Data.samples[0][i].genotype_obs[NULLLOCUS][j] == 100)
                    Data.samples[0][i].genotype_obs[NULLLOCUS][j] = -1;
            }

            if (Data.samples[0][i].genotype_obs[NULLLOCUS][0] < 0)
                Data.samples[0][i].genotype_obs[NULLLOCUS][0] =
                    Data.samples[0][i].genotype_obs[NULLLOCUS][1];
            else if (Data.samples[0][i].genotype_obs[NULLLOCUS][1] < 0)
                Data.samples[0][i].genotype_obs[NULLLOCUS][1] =
                    Data.samples[0][i].genotype_obs[NULLLOCUS][0];
        }

    for (i = 0; i < NUMIND; i++)
        GENOTYPEmutate(Data.samples[0][i].genotype_obs);
    /* output the complete pedigree */
    D = DAGinit(NUMIND, 2, 0);
    /* and the one with mothers only */
    FOPENW(siblings, "ramets.txt");
    D_mothers = DAGinit(NUMIND, 2, 0);
    for (i = 0; i < NUMIND; i++) {
        fprintf
            (out, "%s (%i): born in %i, died in %i (age %i). ",
             Data.samples[0][i].description, Data.samples[0][i].sex, born[i],
             died[i], died[i] - born[i]);
        if (married[i].partner_id >= 0)
            fprintf(out, "Married with %s (%i-?, age %i).",
                    Data.samples[0][married[i].partner_id].description,
                    married[i].year_marriage,
                    married[i].year_marriage - born[i]);
        if (parents[i].mother_id >= 0) {
            fprintf(out, " Child of %s and %s.",
                    Data.samples[0][parents[i].mother_id].description,
                    Data.samples[0][parents[i].father_id].description);
            DAGinsertE(D, EDGE(parents[i].mother_id, i));
            DAGinsertE(D, EDGE(parents[i].father_id, i));
            DAGinsertE(D_mothers, EDGE(parents[i].mother_id, i));
        }
        fprintf(out, "\n");
        Data.samples[0][i].birth = born[i];
        Data.samples[0][i].death = died[i];
        fprintf(siblings, "%i\n", Data.samples[0][i].ramets);
    }
    FCLOSE(siblings);
    FOPENW(siblings, SIBLINGS_OUT);
    for (i = 0; i < NUMIND; i++) {
        for (j = i + 1; j < NUMIND; j++) {
            if (i == j)
                continue;
            indegree1 = DAGincomingE(D, i, incoming1);
            indegree2 = DAGincomingE(D, j, incoming2);
            /* siblings must have the 2 same sampled parents */
            if (indegree1 != indegree2 || indegree1 != 2)
                continue;

            assert(indegree1 == 2);
            assert(indegree2 == 2);
            if ((incoming1[0].v == incoming2[0].v
                 && incoming1[1].v == incoming2[1].v)
                || (incoming1[0].v == incoming2[1].v
                    && incoming1[1].v == incoming2[0].v)) {
                fprintf(siblings, "%10s%10s%10s%10s\n",
                        Data.samples[0][i].description,
                        Data.samples[0][j].description,
                        Data.samples[0][incoming1[0].v].description,
                        Data.samples[0][incoming1[1].v].description);
            }
        }

    }
    FCLOSE(siblings);
    DATAIOpedigreeout(PDG_OUT, D);
    DATAIOpedigreeout(MOTHERS_OUT, D_mothers);
    dummy.num_cluster = 0;

    DATAIOaddGVNodeStyle(D);
    DAGgraphviz(D, "simpedigree.dot", dummy, false);
    //DAGdump(stderr, D);
    for (i = STARTYEAR; i < year; i++) {
        k = 0;
        kr  = 0;
        for (j = 0; j < NUMIND; j++)
            if (born[j] < i && (died[j] < 0 || died[j] > i)) {
                k++;
                kr += Data.samples[0][j].ramets;
            }
        fprintf(out, "year: %i; living: %i; ramets: %i\n", i, k, kr);
    }
    fprintf(out, "living: %i SELFED %i, NOT SELFED %i (rate %f), RAMETS %i ",
            living, selfed, notselfed, selfed / (double)(selfed + notselfed),
            ramets);
    fprintf(out, "CLONAL %i, SEX %i (rate %f)\n", clonal_rep, sex_rep,
            clonal_rep / (double)(sex_rep + clonal_rep));
    DAGdestroy(D);

    DATAIOgenotypeout(GTP_OUT);
    FCLOSE(out);
    /* free malloc'd memory */

    DATAIOdestroyPreMCMC();
    DATAIOdestroyPostMCMC();
    FREQdestroyPreMCMC();
    FREQdestroyPostMCMC();

    exit(EXIT_SUCCESS);
}

/* 
  vim: ft=c sw=4 ts=4 expandtab
*/
