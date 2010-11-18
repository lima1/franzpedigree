/*
 * $Id: exact.c 2062 2010-05-25 10:24:32Z markus $
 *
 * This is an implementation of the exact algorithm described in:
 *
 * Cowell R (2009) Efficient maximum likelihood pedigree reconstruction,
 * Theoretical Population Biology, 76, 285-291.
 *
 * Copyright (C) 2009 Cowell, Robert.
 * Copyright (C) 2010 Universitaet Leipzig.
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "exact.h"

#include "macros.h"
#include "pedigree.h"
#include "prob.h"

extern PROBS Probs;
extern DATA Data;
extern OPTIONS Options;

static inline bool 
IncSubset(EXACT e)
{
    int     k = 1;
    e.subset[k] += 1;
    while (k < Data.num_samples + 2 && e.subset[k] > 1) {
        e.subset[k++] = 0;
        e.subset[k] += 1;
    }
    return e.subset[Data.num_samples + 1] == 0;
}

static inline void
InitSubset(EXACT e)
{
    int i;
    e.subset[0] = 1;
    for (i = 1; i < Data.num_samples + 2; ++i)
        e.subset[i] = 0;
}

static inline int
GetSubsetIndex(EXACT e)
{
    int     pow = 1, k = 1, tot = 0;
    while (k < Data.num_samples + 2) {
        tot += pow * e.subset[k++];
        pow *= 2;
    }
    return tot;
}

EXACT
EXACTinit()
{
    EXACT   e;
    int     k, n = 1;

    MAKE1DINT(e.subset, Data.num_samples + 2);
    k = Data.num_samples;
    while (k > 0) {
        n *= 2;
        k--;
    }
    n++;
    e.sn = n;
    MALLOC(e.sinks, short, n);
    for (k = 0; k < n; ++k)
        e.sinks[k] = 0;
    MAKE1DDOUBLE(e.scores, n);
    for (k = 0; k < n; ++k)
        e.scores[k] = 0.0;
    MAKE1DINT(e.order, Data.num_samples + 1);
    e.D = DAGinit(Data.num_samples, 2, 0);
    return e;
}

void
EXACTdestroy(EXACT e)
{
    FREE1D(e.subset);
    FREE1D(e.sinks);
    FREE1D(e.order);
    DAGdestroy(e.D);
}

static inline double
FindBestScore(EXACT e, int child)
{
    int     i;
    child--;
    assert(child >= 0 && child < Data.num_samples);
    for (i = 0; i < Probs.num_posteriors[child]; i++) {
        if (e.subset[Probs.posteriors[child][i].v + 1] > 0 &&
            e.subset[Probs.posteriors[child][i].w + 1] > 0)
            return Probs.posteriors[child][i].p_opt;
    }
    FATALINT("In exact.c, FindBestScore(). Could not find valid parentage! Try --saexactmax 0 and report bug!");
    return 0;
}

void
EXACTgetBestSinks(EXACT e)
{
    int     j, sink, index = 1, upvars;
    double  skore;

    if (Options.Verbosity > 0)
        VTPROGRESSBARupdate("Pedigree Enumeration",e.sn,0);

    InitSubset(e);
    while (IncSubset(e)) {
        assert(index < e.sn);
        e.scores[index] = 0;
        e.sinks[index] = -1;
        if (Options.Verbosity > 0) 
            VTPROGRESSBARupdate("Pedigree Enumeration",e.sn/16,index/16);
        for (j = 1; j <= Data.num_samples; ++j)
            if (e.subset[j] > 0) {
                e.subset[j] = 0;
                sink = j;
                upvars = GetSubsetIndex(e);
                skore = e.scores[upvars];
                skore += FindBestScore(e, sink);
                if (e.sinks[index] == -1 || skore > e.scores[index]) {
                    e.scores[index] = skore;
                    e.sinks[index] = (short)sink;
                }
                e.subset[j] = 1;
            }
        index++;
    }
    FREE1D(e.scores);
    if (Options.Verbosity > 0) 
        VTPROGRESSBARcompleted("Pedigree Enumeration");
}

void
EXACTfindSinkOrdering(EXACT e)
{
    int     i, bestsink, left = 0;
    for (i = 1; i <= Data.num_samples; ++i)
        e.subset[i] = 1;
    e.subset[0] = 0;
    e.subset[Data.num_samples + 1] = 0;
    for (i = Data.num_samples; i > 0; --i) {
        left = GetSubsetIndex(e);
        bestsink = e.sinks[left];
        e.order[i] = bestsink;
        e.subset[bestsink] = 0;
    }

}

static double
GetBestParents(EXACT e, int child, int *par1, int *par2)
{
    int     i;
    *par1 = *par2 = 0;
    /* Invdiduals start here at 0, not 1 */
    child--;
    for (i = 0; i < Probs.num_posteriors[child]; i++) {
        if (e.subset[Probs.posteriors[child][i].v + 1] > 0 &&
            e.subset[Probs.posteriors[child][i].w + 1] > 0) {
            *par1 = Probs.posteriors[child][i].v + 1;
            *par2 = Probs.posteriors[child][i].w + 1;
            return Probs.posteriors[child][i].p_opt;
        }
    }
    FATALINT("In exact.c, GetBestParents(). Could not find valid parentage! Try --saexactmax 0 and report bug!");
    return 0;
}

double
EXACTfindPedigree(EXACT e)
{
    double  loglik = 0;
    int     i, child, par1, par2;
    bool ret;
    
    InitSubset(e);
    //PROBdumpPosteriors(stderr,Probs.posteriors,NULL,NULL,1);
    for (i = 1; i <= Data.num_samples; ++i) {
        child = e.order[i];
        loglik += GetBestParents(e, child, &par1, &par2);
        e.subset[child] = 1;
        if (par1 > 0) {
            ret = DAGinsertE(e.D,EDGE(par1 - 1, child - 1));
            assert(ret);
        }
        if (par2 > 0) {
            ret = DAGinsertE(e.D,EDGE(par2 - 1, child - 1));
            assert(ret);
        }
    }
    //DAGdump(out, e.D);
    return loglik;
}
