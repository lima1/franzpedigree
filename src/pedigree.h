#ifndef PEDIGREE_H
#define PEDIGREE_H

/*
 * $Id: pedigree.h 2065 2010-05-26 11:21:52Z markus $
 *
 * Extends dag.c with functions that we need for our pedigree stuff
 * but which are not general enough to include them in the DAG library
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
#include "mcmc.h"

bool   PEDIGREEinsertParents(Dag,int,int,int);
void   PEDIGREEremoveParents(Dag,int, Edge[]);
int    PEDIGREEgetPostId(POSTERIORS **, Dag, int);
double PEDIGREEcalcSelfingRate(Dag D);
double PEDIGREEcalcClonalRate(Dag D);

/*estimates the sampling rate of females and males*/
void   PEDIGREEestmateSamplingRate(Dag, double*, double*);
/* here the version when sex data is available. it is called by the function
 * above, so you shouldnt call this directly (public for tests) */
void   PEDIGREEestmateSamplingRateSex(Dag, double*, double*);

/* a graph that visualizes the flow between sampling locations */
void   PEDIGREEpopgraph(MCMC*);
void   PEDIGREEpopgraphList(MCMC*);
void   PEDIGREEpoplist(MCMC_STATE*,unsigned int*);

/* 
  vim: ft=c sw=4 ts=4 expandtab
*/

#endif
