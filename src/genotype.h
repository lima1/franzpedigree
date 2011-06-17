#ifndef GENOTYPE_H
#define GENOTYPE_H

/*
 * $Id: genotype.h 1885 2010-01-25 15:31:17Z markus $
 *
 * Simple functions for genotypes (int**), mainly for simulation
 *
 * Copyright (C) 2008-2010. Universitaet Leipzig
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

int** GENOTYPEcreate(bool onlymalloc); /* creates a genotype, uses allelefreq 
                                        unless onlymalloc                    */
void GENOTYPEmutate(int **);         
void GENOTYPErand(int **);          /* fill the genotype with random alleles */

/* generate offspring genotype, errors according specified rate              */
void GENOTYPEmate(int **offspring, int **mother, int **father);
void GENOTYPEselfing(int **offspring, int **mother);

void GENOTYPEdump(FILE *, int **);
int GENOTYPErandAllele(int locus);
int GENOTYPErandAlleleUnif(int locus, bool null_allele);

/* 
  vim: ft=c sw=4 ts=4 expandtab
*/

#endif
