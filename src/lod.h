#ifndef LOD_H
#define LOD_H

/*
 * $Id: lod.h 2064 2010-05-26 11:20:32Z markus $
 *
 * Calculates the LOD scores and transition probabilities.
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

#include <stdbool.h>
#include "global.h"

/* calculate the constants in the LOD score calculations */
void LODinit();
void LODrecalcErrorConstants();

/* calculate the LOD score numerator */
double LODcalcPchildTriple(int **, int **, int **, int *, bool, bool*);
double LODcalcPchildDyad(int **, int **, int *, bool, bool *);
double LODcalcGenotypeProb(int[], int);               /* genotype locus */

double LODcalcTransProbTriple(int[], int[], int[], int); /* child, p1, p2  */
double LODcalcTransProbDyad(int[], int[], int);      /* child p2 locus */

bool LODhasMismatchLocus(int locus, int **child, int v, int w, bool *ignore_ary);
int LODcalcMismatches(int child_id, int v, int w, bool *ignore_ary);

int LODcmpSetsDyad(int allele, int p[]);
/* the log-likelihoods of H2 */
double LODcalcDenominator(int **, bool *);
double LODcalcDenominatorTripleMotherKnown(int **, int **, bool *);

double LODcalcPchildTripleDelta(int ***genotypes, int child_id, int p1_id, int p2_id, LOCUS_COORD lc, int old_allele, bool *ignore_ary);
double LODcalcPchildDyadDelta(int ***genotypes, int child_id, int p_id, LOCUS_COORD lc, int old_allele, bool *ignore_ary);

double LODcalcDenominatorTripleMotherKnownDelta(int ***genotypes, int child_id, int mother_id, LOCUS_COORD lc, int old_allele, bool *ignore_ary);

double LODcalcDenominatorDelta(int ***genotypes, int child_id, LOCUS_COORD lc, int old_allele, bool *ignore_ary);
/* free malloc'd constants */
void LODdestroy();

/* 
  vim: ft=c sw=4 ts=4 expandtab
*/

#endif
