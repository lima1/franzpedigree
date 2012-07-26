/*
 * $Id: pedigree.c 1977 2010-03-16 16:52:31Z markus $
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
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Markus Riester, markus@bioinf.uni-leipzig.de */

#include "macros.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "global.h"
#include "dag.h"
#include "prob.h"
#include "dataio.h"
#include "pedigree.h"
#include "pvalue.h"
#include "utils.h"

extern DATA Data;
extern PROBS Probs;
extern OPTIONS Options;

bool
PEDIGREEinsertParents(Dag D, int child, int p1, int p2)
{
    int ret;
    Edge tmp[2];

    tmp[0] = EDGE(p1, child);
    tmp[1] = EDGE(p2, child);
    ret = DAGinsertE(D, tmp[0]);
    if (!ret)
        return false;

    ret = DAGinsertE(D, tmp[1]);
    if (!ret) {
        DAGremoveE(D, tmp[0]);
        return false;
    }

    assert(DAGisAcyclic(D));

    return true;
}

void
PEDIGREEremoveParents(Dag D, int child, Edge e[])
{
    DAGremoveIncomingE(D, child, e);
}

int
PEDIGREEgetPostId(POSTERIORS ** posteriors, Dag D, int sample_id)
{
    int ret = -1;
    Edge incoming[2];

    if (D != NULL) {
        incoming[0].v = -1;
        incoming[1].v = -1;
        (void)DAGincomingE(D, sample_id, incoming);
        ret = PROBsearchPosteriors(posteriors, sample_id, incoming[0].v, incoming[1].v);
    }
    return ret;
}

static inline int
cntHomozygousLoci(int i) {
    int j, homo=0;
    for (j=0; j< Data.num_loci;j++)
        if (Data.id_mapping[i].genotype_obs[j][0] !=
            Data.id_mapping[i].genotype_obs[j][1] ) homo++;
    return homo;
}

static void
calcPV(double *sf, double *nsf, int csf, int cnsf, PVALUE_LOOKUP *pv_lup) {
    int i, j, lastsf=sf[0]-1, lastnsf = nsf[0]-1, lastj = 0;
    for (i=0; i<=Data.num_loci; i++)
        PVALUEadd(*pv_lup,i, 0., 1.);

    for (i=0; i< csf; i++) {
        if (lastsf == sf[i]) continue;
        lastsf = sf[i];
        if (nsf[0] >= sf[i]) {
            for (j=0;j<=sf[i];j++) {
  //              fprintf(stderr, "init adding %i\n",j);
                PVALUEadd(*pv_lup,j, 1., 1.);   
            }
            continue;     
        }
    //    fprintf(stderr, "Xlast %i sf[i]: %f j: %i  %i\n",lastsf,sf[i],j,cnsf);
        for (j=0; j < cnsf; j++) {
            if (lastnsf == nsf[j]) continue;
            
            lastj = j;
//             fprintf(stderr, "last %i sf[i]: %f nsf[j] %f. j:%i\n",lastnsf,sf[i],nsf[j],j);
            lastnsf = nsf[j];
            if (sf[i] == nsf[j]) {
//                fprintf(stderr, "adding %i %i\n",i,j);
                PVALUEadd(*pv_lup,sf[i], (i)/(double)(i+j), 1.);
                continue;
            }
        }        
    }                

    //PVALUEfill(*pv_lup);
}

 double
 PEDIGREEcalcSelfingRate(Dag D) 
 {
    int i, indegree, v = DAGgetVertexNumber(D);
    double selfings = 0.,parentages = 0.;
     Edge incoming[2];
     for (i=0; i<v; i++) {
         indegree = DAGincomingE(D, i, incoming);
         parentages += indegree / 2.;
        if (indegree == 2 && incoming[0].v == incoming[1].v)
             selfings++;
    }

     if (CMP_DBL(parentages, 0))
         return 0.;
  
     return selfings / parentages;
 }

double
PEDIGREEcalcSelfingRateN(Dag D) 
{
    int i, indegree, v = DAGgetVertexNumber(D), homo, csf = 0, cnsf=0, cnt=0;
    double selfings = 0.,parentages = 0.,*sf, *nsf, pvrsf = EPSILON, pvrnsf = EPSILON, r;
    Edge incoming[2];
    PVALUE_LOOKUP pv_lup;// = PVALUEinit(Data.num_loci*2,0.,(double)Data.num_loci*2.);
    MAKE1DDOUBLE(sf,Data.num_samples);
    MAKE1DDOUBLE(nsf,Data.num_samples);
    for (i=0; i<v; i++) {
        indegree = DAGincomingE(D, i, incoming);
        parentages += indegree / 2.;
        if (indegree == 2 && incoming[0].v == incoming[1].v) {
            selfings++;
            homo = cntHomozygousLoci(i);
            sf[csf++] = homo;
//            fprintf(stderr, "DEBUG: Selfed: %d\n",homo); 
        } else if (indegree > 0) {
            homo = cntHomozygousLoci(i);
            nsf[cnsf++] = homo;
  //          fprintf(stderr, "DEBUG: Not-Selfed: %d\n",homo); 
        }    
    }    
    pv_lup = PVALUEinit(Data.num_loci+1,0.,(double)Data.num_loci);
    qsort(nsf, cnsf, sizeof(double), UTILScompare_doubles_incr);
    qsort(sf,   csf, sizeof(double), UTILScompare_doubles_incr);
    calcPV(sf,nsf, csf,cnsf,&pv_lup);
//    PVALUEdump(stderr,pv_lup,Data.num_loci+1);
//    fprintf(stderr,"cnting done\n");
    FREE1D(sf);
    FREE1D(nsf);
    for (i=0; i<v; i++) {
        indegree = DAGincomingE(D, i, incoming);
        if (indegree != 0) continue;
        homo = cntHomozygousLoci(i);
        r = PVALUEfind(pv_lup, homo);
//        if (r < 0.3 || r > (1. - 0.3)) {
            pvrnsf += (1. - r);
            pvrsf  += r;
            cnt++;
//           fprintf(stderr,"PV: %i = %f  %f %f %f\n",homo,r, pvrsf, pvrnsf, pvrsf/pvrnsf );
  //      }
    }    
//    fprintf(stderr," %f %f %f %d\n",(pvrsf+selfings)/(pvrnsf+parentages), selfings/parentages, pvrsf/pvrnsf, cnt );
    if (CMP_DBL(parentages, 0))
        return 0.;
 
    return selfings / parentages;
//    return (selfings+pvrsf) / (pvrnsf+parentages);
}
/*
double
PEDIGREEcalcSelfingRate(Dag D, double *selfing_lods, double *selfing_pp) 
{
    int i, indegree, v = DAGgetVertexNumber(D);
    double selfings = 0.,parentages = 0.;
    Edge incoming[2];
    for (i=0; i<v; i++) {
        indegree = DAGincomingE(D, i, incoming);
        parentages += indegree / 2.;
        if (indegree == 2 && incoming[0].v == incoming[1].v)
            selfings++;
        if (0 && indegree == 0 &&selfing_pp[i]>0.99) {
            parentages+= selfing_pp[i];
            //parentages++;
            if (selfing_lods[i] > 0)
                //selfings++;
                selfings += selfing_pp[i];
        }
    }

    if (CMP_DBL(parentages, 0))
        return 0.;
 
    return selfings / parentages;
}*/
/*
double
PEDIGREEcalcClonalRate(Dag D) {
    int clones = 0, i, j,indegree,  v = DAGgetVertexNumber(D), cnt = 0;
    double *a, t,var,avg, outdegree, srr = Data.num_ramets/(double)Options.Nr, c ;
    Edge e[2];

    MAKE1DDOUBLE(a, Data.num_ramets);
    for (i=0; i<v; i++) {
        //if (Data.id_mapping[i].ramets < 2) continue;
        clones += Data.id_mapping[i].ramets - 1;
        c= Data.id_mapping[i].ramets - 1;
        outdegree = 0.;
        for (j=0;j<v;j++) {
            indegree = DAGincomingE(D, j,e);
            if (indegree == 1) {
                if (e[0].v == i) outdegree++;
            }    
            else if (indegree == 2 &&  (e[0].v == i || e[1].v ==i)) outdegree++;
        }   
        if (c> 0)  {
            c *= 1./srr;
            outdegree *= 1/(2.*srr);
            t=c/(double)(c+outdegree);
            for (j=0;j<Data.id_mapping[i].ramets-1;j++) {
                a[cnt++]=t;
            }
        fprintf(stderr, "%f %i %f %f\n",c,indegree,outdegree,t );
        }

        //if (indegree > 0)
         //   cnt++;

    }
    UTILScalcVarianceArray(a,0,cnt,&var,&avg);
        fprintf(stderr, "%i %f %f  sr: %f\n",cnt,avg,sqrt(var), srr);
    exit(1);
    FREE1D(a);
    return clones/(double)(clones+cnt);
}
*/

double
PEDIGREEcalcClonalRate(Dag D) {
    int i, indegree,  v = DAGgetVertexNumber(D), cnt = 0;
    double avg, o, O, sr = Data.num_ramets/(double)Options.Nr, clones =0.,sro, sg,srf, srm;
    Edge e[2];
    FILE *out;
    return 0.;
    FOPENW(out, "test.csv");
    for (i=0; i<v; i++) {
        //if (Data.id_mapping[i].ramets < 2) continue;
        clones += Data.id_mapping[i].ramets-1;
        indegree = DAGincomingE(D, i,e);
        if (indegree > 0) 
            o++;
    }   
    fprintf(stderr, "clones: %f, o: %f sr: %f nr: %i  ng: %i\n",clones,o, sr, Data.num_ramets,Data.num_samples);
    sg = sr;

    sro = (1 - (1. - sg) * (1. - sg));
    O = o / sro;
    clones = Options.Nr - O;
    avg=clones/(double)(clones+O);

    //sg = Data.num_samples/(Options.Nr-clones);
    fprintf(stderr, "clones: %f  o: %f O: %f cr: %f sr: %f sro: %f sg: %f\n",clones,o, O, avg, sr, sro, sg);
    sg = (sr * Data.num_samples) / o;
    sro = (1 - (1. - sg) * (1. - sg));
    O = o / sro;
    clones = Options.Nr - O;
    avg=clones/(double)(clones+O);
    //sg = Data.num_samples/(Options.Nr-clones);
    fprintf(stderr, "clones: %f  o: %f O: %f cr: %f sr: %f sro: %f sg: %f\n",clones,o, O, avg, sr, sro, sg);
    sg = (sr * Data.num_samples) / o;

    Data.has_sex_data = false;
    PEDIGREEestmateSamplingRate(D,&srf, &srm);
    sro = (1 - (1. - sg) * (1. - sg));
    sro = (2*sr)*srf;
    O = o / sro;
    clones = Options.Nr - O;
    avg=clones/(double)(clones+O);
    //sg = Data.num_samples/(Options.Nr-clones);
    fprintf(stderr, "clones: %f  o: %f O: %f cr: %f sr: %f sro: %f sg: %f   %f %i %f %f %f %f\n",clones,o, O, avg, sr, sro, sg,   o,Data.num_samples, o/(double)Data.num_samples, sr/sro, srf,srm);

/*    
    for (i= 0; i< 30000; i++) {
    sg = i/300000.;
    sro = (1 - (1. - sg) * (1. - sg));
    O = o / sro;
    clones = Options.Nr - O;
    avg=clones/(double)(clones+O);
    //sg = Data.num_samples/(Options.Nr-clones);
    //fprintf(stderr, "%f %f (%f) %f sr: %f %f sg: %f\n",clones,offspring, o, avg, sr, sro, sg);
    if (avg > 0.)
    fprintf(out, "%f,%f,%f,%f\n",avg, sro, sg,sr);
    }
*/    
    FCLOSE(out);
    FATAL("");
    return clones/(double)(clones+cnt);
}
void
PEDIGREEestmateSamplingRate(Dag D, double *sr_f, double *sr_m)
{
    int i, indegree_cnt[3] = { 0, 0, 0};
    double resf;
    if (!Options.IgnoreSex && Data.has_sex_data) {
        PEDIGREEestmateSamplingRateSex(D,sr_f,sr_m);
        return;
    }
        
    for (i = 0; i < Data.num_samples; i++) 
        indegree_cnt[DAGindegree(D, i)]++;
    
    resf =  1. / (double) ( indegree_cnt[1] / (double) (2.*indegree_cnt[2]) + 1);
    if (isnan(resf))
        resf = D_MINSAMPLINGRATE;
    else {
        resf = MIN(resf, 1.);
        resf = MAX(D_MINSAMPLINGRATE,resf);
    }
    *sr_f = resf;
    *sr_m = resf;
}

void 
PEDIGREEestmateSamplingRateSex(Dag D, double *sr_f, double *sr_m) {
    int i, indegree, p1;
    double cf=0, cm=0;
    Edge incoming[2];
    for (i = 0; i < Data.num_samples; i++) {
        /* only use offspring generations */
        if (Probs.num_candidates[i] < 1) continue;
        indegree = DAGincomingE(D,i,incoming);
        p1 = incoming[0].v;
        if (indegree > 0) {
            /* first parent is mother */
            if (Data.id_mapping[p1].sex == 1) {
                cf++;
                if (indegree == 2)
                    cm++;
            }
            /* father */
            else if (Data.id_mapping[p1].sex == 2) {
                cm++;
                if (indegree == 2)
                    cf++;
            }
            else {
                if (indegree == 1) {
                    cf += 0.5;
                    cm += 0.5;
                }
                else {
                    cf++;
                    cm++;
                }
            }
        }
    }
    //ind_f /= 2; ind_m /= 2;
    *sr_f = MIN( 1., MAX(D_MINSAMPLINGRATE, ( cf / (double) Probs.num_candidate_mothers) ));
    *sr_m = MIN( 1., MAX(D_MINSAMPLINGRATE, ( cm / (double) Probs.num_candidate_fathers) ));
   
    assert(*sr_f <= 1.0 && *sr_m <= 1.0);
}

void
PEDIGREEpopgraph(MCMC * s)
{
    Dag D;
    GraphvizCluster cluster;
    int i, child_pop, p1_pop, p2_pop;
    double dcp1, dcp2, dp1p2, max = 1;

/*    char graphviz[LINESIZE];*/

    Edge *e;

    D = DAGinit(Data.num_populations, -1, 1);
    for (i = 0; i < Data.num_populations; i++) {
        DAGlabelV(D, i, Data.populations[i].description, 50);
        /*DAGgraphvizV(D, i, "pos=\"\"");*/
    }

    MALLOC(e, Edge, 2);

    for (i = 0; i < Data.num_samples; i++) {
        if (DAGindegree(s->best.D, i) == 0)
            continue;
        DAGincomingE(s->best.D, i, e);
        /* get the population ids */
        child_pop = Data.id_mapping[i].population->id;
        p1_pop = Data.id_mapping[e[0].v].population->id;
        p2_pop = Data.id_mapping[e[1].v].population->id;

        /* printf("%i (%i):\n", i,child_pop);
           printf("\t(%i %i) (%i %i)", e[0].v, e[0].w,e[1].v,e[1].w);
           printf("\t%i: %i %i: %i\n", e[0].v, p1_pop, e[1].v,p2_pop ); */
        if (Data.use_distances) {
            dcp1 = Data.population_distances[child_pop][p1_pop];
            dcp2 = Data.population_distances[child_pop][p2_pop];
            dp1p2 = Data.population_distances[p1_pop][p2_pop];
            if (dcp1 > max && dcp2 > max && dp1p2 > max)
                continue;
            /*printf("Dist C-P1: %f, C-P2: %f, P1-P2:%f\n",dcp1,dcp2,dp1p2); */
        }
        /* don't draw triples where all members are from the same
         * population */
        if (p1_pop != child_pop)
            DAGinsertE(D, EDGE(p1_pop, child_pop));
        if (p2_pop != child_pop)
            DAGinsertE(D, EDGE(p2_pop, child_pop));
    }
    FREE(e);
    cluster.num_cluster = 0;
    DAGgraphviz(D, "pop.dot", cluster, false);
    /*DAGdump(stderr,D); */
    DAGdestroy(D);
}

void
PEDIGREEpoplist(MCMC_STATE * s, unsigned int * matrix)
{
    int i, child_pop, p1_pop, p2_pop, indegree;
    Edge e[2];

    for (i = 0; i < Data.num_samples; i++) {

        indegree = DAGindegree(s->D, i);

        if (indegree == 0) continue;

        DAGincomingE(s->D, i, e);
        /* get the population ids */
        child_pop = Data.id_mapping[i].population->id;
        p1_pop = Data.id_mapping[e[0].v].population->id;
        matrix[TIDXL(child_pop, p1_pop)]++;
        if (indegree == 2) {
            p2_pop = Data.id_mapping[e[1].v].population->id;
            matrix[TIDXL(child_pop, p2_pop)]++;
        }
    }
}

void
PEDIGREEpopgraphList(MCMC * mcmc)
{
    Dag D;
    GraphvizCluster cluster;
    int i, j, k;
    double var, avg, *tmp;

/*    char graphviz[LINESIZE];*/

    MAKE1DDOUBLE(tmp, Probs.mh_sampled_pedigrees);

    D = DAGinit(Data.num_populations, -1, 1);
    for (i = 0; i < Data.num_populations; i++) {
        DAGlabelV(D, i, Data.populations[i].description, 10);
        /*DAGgraphvizV(D, i, "pos=\"\"");*/
    }

    for (i = 0; i < Data.num_populations; i++)
        for (j = 0; j < i; j++) {
            for (k = 0; k<Probs.mh_sampled_pedigrees; k++) 
                tmp[k] = (double)mcmc->mh_pop[k][TIDX(i,j)];
            UTILScalcVarianceArray(tmp,0,Probs.mh_sampled_pedigrees-1,&var,&avg);
            if (i!=j)
                for (k=0; k<(int)avg;k++)
                    DAGinsertE(D, EDGE(i, j));

        }
    cluster.num_cluster = 0;
    DAGgraphviz(D, "pop.dot", cluster, false);
    /*DAGdump(stderr,D); */
    DAGdestroy(D);
    FREE1D(tmp);
}

/* 
  vim: ft=c sw=4 ts=4 expandtab
*/
