#ifndef NULL_H
#define NULL_H

/*
 * $Id: null.h 568 2008-05-18 20:07:42Z markus $
 *
 * estimates null allele frequency
 *
 * Copyright (C) 2008 Universitaet Leipzig  
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

/* the default method */
#define NULLcalcPn NULLcalcPnKalinowskiTaper

typedef struct S_NULLTEST *NULLTEST;

NULLTEST NULLinit(int locus);                   /* alloc required memory */
double   NULLcalcPnChakraborty(NULLTEST);       
double   NULLcalcPnKalinowskiTaper(NULLTEST);
void     NULLdestroy(NULLTEST);
void     NULLdump(FILE*, NULLTEST);             /* dump some statistics */

/* 
  vim: ft=c sw=4 ts=4 expandtab
*/

#endif
