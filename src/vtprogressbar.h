#ifndef VTPROGRESSBAR_H
#define VTPROGRESSBAR_H

/*
 * $Id: vtprogressbar.h 1885 2010-01-25 15:31:17Z markus $
 *
 * Most code shamelessly stolen from Steve Hoffman
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

void VTPROGRESSBARupdate(char *message, unsigned int complete, unsigned int processed);
void VTPROGRESSBARupdateNP(char *message);
void VTPROGRESSBARinit(unsigned int);
void VTPROGRESSBARcompleted(char *);
void VTPROGRESSBARcursorVisible(void);
void VTPROGRESSBARcursorInvisible(void);

/* 
  vim: ft=c sw=4 ts=4 expandtab
*/

#endif
