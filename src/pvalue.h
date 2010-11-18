#ifndef PVALUE_H
#define PVALUE_H

/*
 * $Id: pvalue.h 2064 2010-05-26 11:20:32Z markus $
 *
 * This file provides a simple data structure for p-values. In a first step,
 * the data structure is initialized by defining the range of the distribution
 * and a grid size. For example:
 *
 *   *pv_lup = PVALUEinit(100, -1, 3);
 *
 * This will initialize a lookup table in the range [-1, 3]. This range is
 * devided in a grid of size 100. Values smaller than -1 will return the value
 * stored for -1, values > 3 the one stored for 3. You should make sure that
 * these ends have p-Values of 0.0 and 1.0, respectively. 
 * 
 * In the second step, you add the p-values to the distribution. For example
 *
 *   PVALUEadd(*pv_lup, 2.5, 0.001);
 *
 * This says a value of 2.5 has a p-value of 0.001.
 *
 * In a last step, you can "smooth" the distribution:
 *
 *   PVALUEfill(*pv_lup);
 *
 * This will fill unset grid-points with the values of their nearest 
 * neighbors. 
 *
 * Now you can lookup p-values in O(1), for example:
 *
 *  PVALUEfind(*pv_lup, 2.64);
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

typedef struct S_PVALUE_LOOKUP *PVALUE_LOOKUP;

PVALUE_LOOKUP 
       PVALUEinit(int length, double min, double max);
void   PVALUEdestroy(PVALUE_LOOKUP);
void   PVALUEdump(FILE *fp, PVALUE_LOOKUP pv_lup, int grid);
double PVALUEfind(PVALUE_LOOKUP pv_lup, double v);
double PVALUEfindCritical(PVALUE_LOOKUP pv_lup, double pv);

void   PVALUEadd(PVALUE_LOOKUP pv_lup, double v, double pv, double sensitivity);
void   PVALUEfill(PVALUE_LOOKUP);

/* 
  vim: ft=c sw=4 ts=4 expandtab
*/

#endif
