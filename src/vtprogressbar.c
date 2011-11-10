/*
 * $Id: vtprogressbar.c 1885 2010-01-25 15:31:17Z markus $
 *
 * Most code shamelessly stolen from Steve Hoffman
 *
 * http://en.wikipedia.org/wiki/ANSI_escape_code
 *   ... for an explanation what these funny fprintf statements
 *   do.
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

#ifdef WIN32
#include "windows.h"
#endif
#include "macros.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include "vtprogressbar.h"

static unsigned int size;
static int lastval;
static time_t started;

void
VTPROGRESSBARcursorInvisible(void)
{
#ifdef WIN32
     CONSOLE_CURSOR_INFO cciCursor;
     HANDLE hStdOut = GetStdHandle(STD_OUTPUT_HANDLE);
     
     if (GetConsoleCursorInfo(hStdOut, &cciCursor)) {
          cciCursor.bVisible=FALSE;
          (void)SetConsoleCursorInfo(hStdOut, &cciCursor);
     }     
#else
    fprintf(stderr, "%c%c%c%d%c", 27, '[', '?', 25, 'l');
#endif
}

void
VTPROGRESSBARcursorVisible(void)
{
#ifdef WIN32
     CONSOLE_CURSOR_INFO cciCursor;
     HANDLE hStdOut = GetStdHandle(STD_OUTPUT_HANDLE);
     
     if (GetConsoleCursorInfo(hStdOut, &cciCursor)) {
          cciCursor.bVisible=TRUE;
          (void)SetConsoleCursorInfo(hStdOut, &cciCursor);
     }     
#else
    fprintf(stderr, "%c%c%c%d%c", 27, '[', '?', 25, 'h');
#endif
}

void
VTPROGRESSBARinit(unsigned int s)
{
#ifndef WIN32
    fprintf(stderr, "%c%c%c", 27, '[', 's');
    fprintf(stderr, "%c%c%c", 27, '[', 'K');
#endif
    size = s;
    lastval = -1;
    started = 0;
    VTPROGRESSBARcursorInvisible();
}

void
VTPROGRESSBARcompleted(char *message)
{
    VTPROGRESSBARupdate(message, 100, 100);
    fprintf(stderr, "\n");
    lastval = -1;
    started = 0;
}

inline static void
cursorUp(void)
{
#ifdef WIN32
    CONSOLE_SCREEN_BUFFER_INFO csbi;
    COORD newCP;
    HANDLE hStdOut = GetStdHandle(STD_OUTPUT_HANDLE);
    GetConsoleScreenBufferInfo(hStdOut, &csbi);
    newCP.X  = 0; 
    newCP.Y  = csbi.dwCursorPosition.Y - 1;
    SetConsoleCursorPosition(hStdOut, newCP); 
#else
    fprintf(stderr, "%c%c%c", 27, '[', 'A');
#endif
}

void
VTPROGRESSBARupdate(char *message, unsigned int complete, unsigned int processed)
{
    unsigned int i, percent, bar;
    unsigned long remaining;
    time_t seconds;

    percent = (processed * 100) / (complete);
    if (lastval >= (int)percent)
        return;
    lastval = (int)percent;
    
    bar = (size * processed) / (complete);

    fprintf(stderr, "[");

    for (i = 0; i < bar; i++) 
        fprintf(stderr, "=");
    for (i = bar; i < size; i++) 
        fprintf(stderr, " ");

    fprintf(stderr, "]  %3u%c  %-34s", percent, '%', message);

    if (started == 0)
        started = time(NULL);
    seconds = time(NULL) - started; 
    if (seconds > 3 && percent != 100) {
        remaining = (unsigned long)(seconds * (complete-processed)/processed);
        if (remaining < 60)
            fprintf(stderr, "ETA: %4lu sec\n", remaining);
        else if (remaining < 3600)
            fprintf(stderr, "ETA: %4lu min\n", (unsigned long)ceil( remaining / 60.));
        else 
            fprintf(stderr, "ETA: %4lu h\n", (unsigned long)MIN(9999., ceil( remaining / 3600.)));
    }    
    else 
        fprintf(stderr, "%13s\n"," ");

    cursorUp();
}

void
VTPROGRESSBARupdateNP(char *message)
{
    unsigned int i;
    time_t now = time(NULL);
    static int pos = 0;
    static int failed = 0;
    static bool forward = true;
    /* we don't need started here, but we update only every second */
    if (started == now && failed < 10) {
        failed++;
        return;
    }    
    started = now;
    failed = 0;

    pos = forward ? pos + 1 : pos - 1;
    if (pos < 1) {
        forward = true;
        pos = 2;
    }
    else if (pos >= size -1) {
        forward = false;
        pos = size - 3;
    }

    fprintf(stderr, "[%*s",pos-1,"");

    if (pos == 0)
        fprintf(stderr, "=>");
    else if (pos == size - 1)
        fprintf(stderr, "<=");
    else
        fprintf(stderr, "<=>");
    for (i = pos + 2; i < size; i++) 
            fprintf(stderr, " ");
    fprintf(stderr, "]  %4s  %-34s\n", " ",  message);
    cursorUp();
}

/* 
  vim: ft=c sw=4 ts=4 expandtab
*/
