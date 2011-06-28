#ifndef OPTIONS_H
#define OPTIONS_H
/*
 * $Id: options.h 2065 2010-05-26 11:21:52Z markus $
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

#include "global.h"

void OPTIONSinit(void);
void OPTIONSdump(FILE *);
void OPTIONSdumpSim(FILE *);
void OPTIONSprintBool(FILE *, bool);

/* 
  vim: ft=c sw=4 ts=4 expandtab
*/

#endif
