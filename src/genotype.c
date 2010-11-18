/*
 * $Id: genotype.c 1885 2010-01-25 15:31:17Z markus $
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

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "genotype.h"

#include "global.h"
#include "macros.h"
#include "freq.h"

extern DATA Data;       /* for num loci         */
extern OPTIONS Options; /* for the typing error */

int**
GENOTYPEcreate(bool onlymalloc)
{
    int i;
    int **g;

    MAKE2DINT(g, Data.num_loci, D_PLOIDY, i);
    if (onlymalloc)
        return g;
    GENOTYPErand(g);
    return g;
}

void
GENOTYPEmate(int **offspring, int **mother, int **father)
{
    int i;
    
    for (i = 0; i < Data.num_loci; i++) {
        if (RANDDBLONE < 0.5) {
            offspring[i][0] = mother[i][0];
            offspring[i][1] = father[i][1];
        }
        else {
            offspring[i][0] = mother[i][1];
            offspring[i][1] = father[i][0];
        }    
    }
}

void
GENOTYPEselfing(int **offspring, int **mother) 
{    
    int i;
    
    for (i = 0; i < Data.num_loci; i++) {
        if (RANDDBLONE < 0.5) 
            offspring[i][0] = mother[i][0];
        else 
            offspring[i][0] = mother[i][1];

        if (RANDDBLONE < 0.5) 
            offspring[i][1] = mother[i][0];
        else 
            offspring[i][1] = mother[i][1];
    }
}

/* sample allele according their frequency */
int
GENOTYPErandAllele(int locus)
{
    int i;
    double fsum, r;

    fsum = 0.;
    r = RANDDBLONE;
    for (i = 0;
         i <
         Data.allele_frequencies[locus].max -
         Data.allele_frequencies[locus].min + 1; i++) {
        fsum += Data.allele_frequencies[locus].freqs[i];
        if (fsum > r)
            return i + Data.allele_frequencies[locus].min;
    }
    assert(0); 
    return 0;
}

/* sample alleles uniformely. */
int
GENOTYPErandAlleleUnif(int locus, bool null_allele)
{
    int r;

    if (null_allele) {
        r = RANDINT((double)(Data.num_alleles[locus] + 1.));
        /* return null allele */            
        if (r == Data.num_alleles[locus])
            return -1;
    }
    else {
        r = RANDINT((double)(Data.num_alleles[locus]));
    }
    
    assert(r >= 0 && r < Data.num_alleles[locus]);

    return Data.alleles[locus][r];
}

void
GENOTYPErand(int **g)
{
    int i;

    for (i = 0; i < Data.num_loci; i++) {
        g[i][0] = GENOTYPErandAllele(i);
        g[i][1] = GENOTYPErandAllele(i);
    }
}

void
GENOTYPEmutate(int **g)
{
    int i, typed = Data.num_loci;
    double p;

    for (i = 0; i < Data.num_loci; i++) {
        p = RANDDBLONE;
        if (p < Data.TypingError[i]) {
            g[i][0] = GENOTYPErandAllele(i);
            p = RANDDBLONE;
            if (p < Data.TypingError[i])
                g[i][1] = GENOTYPErandAllele(i);
        }
        p = RANDDBLONE;
        if (p < 1. - Options.ProportionTyped[i]
            && typed > Options.MinTyped) {
            g[i][0] = -1;
            g[i][1] = -1;
            typed--;
        }
    }
}

void
GENOTYPEdump(FILE * fp, int **g)
{
    int i;

    for (i = 0; i < Data.num_loci; i++) {
        fprintf(fp, "%i.%i", g[i][0], g[i][1]);
        if (i < Data.num_loci - 1)
            fprintf(fp, " ");
    }
}


/* 
  vim: ft=c sw=4 ts=4 expandtab
*/
