#ifndef HWE_H
#define HWE_H

/*
 * $Id: hwe.h 2058 2010-05-12 15:19:58Z markus $
 *
 * Performing the Exact Test of Hardy-Weinberg Proportion for
 * Multiple Alleles. Guo, Thompson.
 *
 * Copyright (C) 1992. Sun-Wei Guo.
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

#include "global.h"

typedef struct S_HWETEST *HWETEST;

/* malloc matrix etc*/
  HWETEST HWEinit(int max_alleles);

     void HWEcalcFreqs(HWETEST hwe, int ***genotypes, int samples, int locus, int,
                  FREQS*);
      int HWEgetFreq(HWETEST hwe, int, int);
      int HWEgetAlleleFreq(HWETEST hwe, int);
     void HWEtest(HWETEST hwe, int chunks, int chunksize, int steps);
HWERESULT HWEgetResult(HWETEST hwe);

/* free matrix */
     void HWEdestroy(HWETEST hwe);
     void HWEdumpFreqs(FILE *, HWETEST);
     void HWEdumpResults(FILE *, HWETEST);

/* exported for test suite */     
double calcSampleProb(HWETEST, int, int, int);
void createSumLog2(HWETEST);

/* 
  vim: ft=c sw=4 ts=4 expandtab
*/

#endif
