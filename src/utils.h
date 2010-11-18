#ifndef UTILS_H
#define UTILS_H

/*
 * $Id: utils.h 1977 2010-03-16 16:52:31Z markus $
 *
 * misc useful methods, no dependencies
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

#include <stdbool.h>

/* mean earth radius for hsin */
#define METER_RHO 6371640

/* String Utilities */
void UTILScopyNonBlanks(char *source, char *dest, unsigned int min, unsigned int max);
void UTILSstripNewline(char *str, size_t size);
void UTILSreplaceSpace(char *str, char replace, size_t size);
void UTILSremoveSpecialChars(char *, char*,size_t);

/* Math             */
double UTILSfactrl(int n);
void UTILSfactlnInit(int max);
double UTILSfactln(int n);
void UTILSfactlnDestroy(void);

void UTILScalcVarianceArray(double *arr, unsigned int start, unsigned int end, double *var,
                            double *avg);

void UTILScalcVarianceArrayInt(int *arr, unsigned int start, unsigned int end, double *var,
                            double *avg);
double UTILSgammln(double xx);
double UTILSchoose(int N, int k);
double UTILSbeta(double z, double w);
double UTILSbinln(int N, int k, double p);
double UTILSmultiln(int N, int *k, float *p, int q);


/* Misc             */

double UTILSdeg2rad(double deg);
/* return distance in km for lat1,long1, lat2,long2*/
double UTILScalcDistance(double lat1, double lon1, double lat2, double lon2);

int UTILScompare_doubles_decr (const void *a, const void *b);
int UTILScompare_doubles_incr (const void *a, const void *b);
int UTILScompare_ints_incr (const void *a, const void *b);

/* benjamini-hochberg adjusted p-values. this function SORTS p !! */
void UTILSpvAdjustBH(double *p, double *p_adjusted, int n);
void UTILSpvAdjustHolm(double *p, double *p_adjusted, int n);

#include "pvalue.h"

void UTILScmpDistributions(double *a1, double *a2, size_t n1, size_t n2, PVALUE_LOOKUP *pv_lup, int w);
void UTILSshuffleInt(int *, size_t);

/* 
  vim: ft=c sw=4 ts=4 expandtab
*/

#endif
