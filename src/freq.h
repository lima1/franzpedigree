#ifndef FREQ_H
#define FREQ_H
/*
 * $Id: freq.h 1885 2010-01-25 15:31:17Z markus $
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
 *
 * Markus Riester, markus@bioinf.uni-leipzig.de */

#include "dag.h"
#include "global.h"

void FREQinit(void);
void FREQdestroyPreMCMC(void);
void FREQdestroyPostMCMC(void);

void FREQcalcAlleleFrequencies(void);
void FREQupdate(Dag);       /* updates the frequencies using the DAG */
void FREQcalcSummaryStatistics(void);
void FREQdump(FILE * fp);
void FREQdumpSummaryStatistics(FILE * fp);

void FREQgetAlleleRange(int locus, RANGE*);  

/* some macros for accessing the frequencies */
#define FREQ_ALLELE_ID(locus,allele) ( (allele < 0) ? -1 : \
    Data.allele_frequencies[locus].allele_id_lookup[allele \
        - Data.allele_frequencies[locus].min] )
#define FREQ_ALLELE_PROB(locus,allele) ( ( allele < 0 ) ? 1. : \
    Data.allele_frequencies[locus].freqs[allele - Data.allele_frequencies[locus].min] )

/* missing alleles never cause a mismatch */
#define FREQ_CMP_ALLELES(a1,a2) ( (a1 < 0 || a2 < 0 || a1 == a2) ? 1: 0 )
#define FREQ_CMP_ALLELES_2(locus,a1,a2) ( ( a1 <  0 ) ? 1 : \
                                          ( a2 <  0 ) ? FREQ_ALLELE_PROB(locus,a1) : \
                                          ( a1 == a2) ? 1 : 0 )

/* all alleles are stored in a array from 0 to FREQ_RANGE */
#define FREQ_RANGE(locus) ( Data.allele_frequencies[locus].max - Data.allele_frequencies[locus].min + 1 )

/* 
  vim: ft=c sw=4 ts=4 expandtab
*/

#endif
