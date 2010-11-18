#ifndef DATAIO_H
#define DATAIO_H

/*
 * $Id: dataio.h 1885 2010-01-25 15:31:17Z markus $
 *
 * The parser for the input files.
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
#include "global.h"

void DATAIOinit(void);
void DATAIOinitFullsibMatrix(void);

/* input files */
void DATAIOinfile(char *);
void DATAIOfullsibFile(char *);
void DATAIOfullsibOut(char *);
void DATAIOgeofile(char *, bool);
void DATAIOfreqfile(char *, int);
void DATAIOpedigreein(char *, int**);

void DATAIOinitInpedigree();

/* output files */
void DATAIOpedigreeout(char *, Dag);
void DATAIOpedigreeoutText(char *, Dag);
void DATAIOgenotypeout(char *);

void DATAIOestimateTypingError(void);
void DATAIOcheckData(FILE *);
void DATAIOdump(FILE *);
void DATAIOdumpSamplingLocations(FILE *);
void DATAIOdumpDatasetDetails(FILE *);
void DATAIOdumpGenotypesCERVUS(void);
void DATAIOdumpOffspringCERVUS(void);
void DATAIOdumpGenotypesGenepop(void);
void DATAIOdumpGenotypesRMES(void);
void DATAIOdumpGenotypesParente(void);
void DATAIOdumpMismatches(void);

void DATAIOaddGVNodeStyle(Dag);
GraphvizCluster DATAIOgraphvizCluster(void);
void DATAIOgraphvizClusterDestroy(GraphvizCluster);
void DATAIOcreateMapping(void);
void DATAIOdestroyMapping(void);

void DATAIOdestroyPreMCMC(void);
void DATAIOdestroyPostMCMC(void);

int *** DATAIOcloneGenotypes(void);
void DATAIOcopyGenotypes(int ***, int ***);
void    DATAIOdestroyGenotypes(int *** gcloned);

/* 
  vim: ft=c sw=4 ts=4 expandtab
*/
#endif
