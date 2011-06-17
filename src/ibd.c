/*
 * $Id: ibd.c 2064 2010-05-26 11:20:32Z markus $
 *
 * See DNA-based methods for pedigree reconstruction and kinship analysis in
 * natural populations, Michael S. Blouin.
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
#include <math.h> /* for log */

#include "ibd.h"

#include "freq.h"

extern DATA Data; /* for the typing error */

/* probabilities that a dyad has a pair of genotypes given that they share m
 * alleles identical by descent 
 * See DNA-based methods for pedigree reconstruction and kinship analysis in
 * natural populations, Michael S. Blouin, Box 3
 */

void
IBDcalcGenotypeProbabilities(int *g1, int *g2, double *m0, double *m1, double *m2, int locus)
{
    int i, j;
    double p1[2], p2[2];

    p1[0] = FREQ_ALLELE_PROB(locus, g1[0]);
    p1[1] = FREQ_ALLELE_PROB(locus, g1[1]);
    p2[0] = FREQ_ALLELE_PROB(locus, g2[0]);
    p2[1] = FREQ_ALLELE_PROB(locus, g2[1]);

    /* see blouin 2003 */
    *m0 = p1[0] * p1[1] * p2[0] * p2[1];
    *m1 = 0;
    *m2 = 0;
    if (g1[0] != g1[1])
        *m0 *= 2;
    if (g2[0] != g2[1])
        *m0 *= 2;

    if (g1[0] != g1[1]) {
        if ((g1[0] == g2[0] && g1[1] == g2[1]) ||
            (g1[0] == g2[1] && g1[1] == g2[0])) {
            /*aiaj aiaj */
            *m1 = p1[0] * p1[1] * (p1[0] + p1[1]);
            *m2 = 2 * p1[0] * p1[1];
            return;
        } else if (g2[0] != g2[1]) {
            for (i = 0; i < 2; i++)
                for (j = 0; j < 2; j++)
                    if (g1[i] == g2[j]) {
                        if (i == 0)
                            *m1 = p1[1] * p2[0] * p2[1];
                        else
                            *m1 = p1[0] * p2[0] * p2[1];
                        return;
                    }
        } else {
            /* aiaj aiai */
            if (g1[0] == g2[0]) {
                *m1 = p1[0] * p1[0] * p1[1];
                return;
            }
            /* aiaj ajaj */
            else if (g1[1] == g2[0]) {
                *m1 = p1[0] * p1[1] * p1[1];
                return;
            }
        }
    } else {
        /* aiai aiai */
        if (g1[0] == g2[0] && g2[0] == g2[1]) {
            *m1 = p1[0] * p1[0] * p1[0];
            *m2 = p1[0] * p1[0];
            return;
        }
        /* aiai ai ? */
        else if (g1[0] == g2[0]) {
            *m1 = p1[0] * p1[0] * p2[1];
            return;
        }
        /* aiai ? ai */
        else if (g1[0] == g2[1]) {
            *m1 = p1[0] * p1[0] * p2[0];
            return;
        }
    }
}

void
IBDcalcRelationshipLikelihoods(int **g1, int **g2, IBD_LL *ibd)
{
    int i;
    double m0, m1, m2;

    ibd->po = 0.;
    ibd->fs = 0.;
    ibd->hs = 0.;
    ibd->u  = 0.;
    for (i = 0; i < Data.num_loci; i++) {

        if (g1[i][0] < 0 || g1[i][1] < 0 || g2[i][0] < 0 || g2[i][1] < 0) continue;

        IBDcalcGenotypeProbabilities(g1[i], g2[i], &m0, &m1, &m2, i);
        /* typing errors */
        m1 = (1.-Data.TypingError[i]) * m1 + Data.TypingError[i]*m0; 
        m2 = (1.-Data.TypingError[i]) * m2 + Data.TypingError[i]*m0; 

        ibd->po += log(m1);
        ibd->fs += log(0.25 * m0 + 0.5 * m1 + 0.25 * m2);
        ibd->hs += log(0.5 * m0 + 0.5 * m1);
        ibd->u  += log(m0);
    }
}
/* 
  vim: ft=c sw=4 ts=4 expandtab
*/
