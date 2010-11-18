#ifndef EXACT_H
#define EXACT_H
/*
 * $Id: exact.h 1919 2010-02-09 20:46:04Z markus $
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
 
#include "dag.h"
#include "global.h"

typedef struct S_EXACT
{
	int * subset;    // array 0..n+1, denoting presence or not in the subset {1,..n}. //
	Dag D;           // The Maximum Likelihood Pedigree                               //
	short * sinks;
	double * scores;
	int * order;
    int sn;
} EXACT;

EXACT EXACTinit();
void EXACTdestroy();
void EXACTclear();
void EXACTgetBestSinks();
void EXACTfindSinkOrdering();
double EXACTfindPedigree(EXACT e); // returns log likelihood of pedigree found.

#endif
