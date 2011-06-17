#ifndef PROB_H
#define PROB_H
/*
 * $Id: prob.h 2064 2010-05-26 11:20:32Z markus $
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

#include "dag.h"
#include "ibd.h"
#include "global.h"
#include "mcmc.h"

/* calculates constants in the likelihood formulas */
void PROBinit(void);

/* the two main functions that calculate all probabilities */
void PROBcalcTriplesAndDyads(void);

void PROBcalcFullsibs(void);

/* this function merges all hypotheses in one array (instead of a linked list)
 * and sorts them by their posterior probability */

void PROBcalcPosteriors(void);

/* recalc the transition and genotype probabilites */
void PROBrecalcPosteriors(MCMC_STATE*);

void PROBrecalcPosteriorProbabilities(POSTERIORS**, double, double);

/* update probabilities after changing the genotype of one individual */
void PROBupdateGenotype(MCMC_STATE *, LOCUS_COORD, int);

POSTERIORS** PROBclonePosteriors(POSTERIORS**);

bool PROBisCandidateParent(int, int);

/* use age data to determine if they could be siblings */
bool PROBcouldBeSiblings(int, int);

int PROBsearchPosteriors(POSTERIORS **, int, int, int);

/* helper functions */
int  PROBcommonTypedLoci(int id1, int id2);
void PROBdumpPosteriors(FILE *, POSTERIORS **, Dag, Dag, int);
void PROBmakeIgnoreArrays(bool);
void PROBdestroyPosteriors();
void PROBdestroyPreMCMC(void);
void PROBdestroyPostMCMC(void);

/* 
  vim: ft=c sw=4 ts=4 expandtab
*/

#endif
