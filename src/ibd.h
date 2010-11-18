#ifndef IBD_H
#define IBD_H

/*
 * $Id: ibd.h 2064 2010-05-26 11:20:32Z markus $
 *
 * calculates the IBD likelihood ratios.
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

/* the IBD log-likelihoods for parent-offspring, full-/half-sib and unrelated
 * dyads */
typedef struct S_IBD_LL {
    double po;
    double fs;
    double hs;
    double u;
} IBD_LL;

void IBDcalcGenotypeProbabilities(int *, int *, double *, double *, double *, int locus);
void IBDcalcRelationshipLikelihoods(int **g1, int **g2, IBD_LL *ibd);

/* 
  vim: ft=c sw=4 ts=4 expandtab
*/

#endif
