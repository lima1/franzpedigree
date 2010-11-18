#ifndef LL_H
#define LL_H

/*
 * $Id: ll.h 700 2008-06-06 13:06:59Z markus $
 *
 * C from http://en.wikipedia.org/wiki/Linked_list
 * and http://www.chiark.greenend.org.uk/~sgtatham/algorithms/listsort.html
 *
 * The Preprocessor ugliness is mine... 
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

#define LLAPI(T,PREFIX)\
int    PREFIX##cnt(T **); \
void   PREFIX##remove(T **); \
void   PREFIX##destroy(T ** n); \
void   PREFIX##dump(FILE*,T *); \
T     *PREFIX##sort(T *); \

#define LLREMOVE(PREFIX,T,P) \
void PREFIX##remove(T ** P){\
    if (*P != NULL) {     \
        T *n = *P;        \
        *P = (*P)->next;  \
        free(n);          \
    }}

#define LLCNT(PREFIX,T,P) \
int PREFIX##cnt(T ** P){\
    T *n; \
    int i = 0; \
    for (n = *P; n != NULL; n=n->next) { \
        i++; \
    }\
   return i;\
}\

#define LLDESTROY(PREFIX,T,P) \
void PREFIX##destroy(T ** P){\
    while (*P != NULL) { \
        PREFIX##remove(P);   \
    }}

#define LLDUMP(PREFIX,T,P) \
void PREFIX##dump(FILE *fp, T * P){\
    if (P == NULL) {  \
        fprintf(fp,"list is empty"); \
    } \
    while (P != NULL) { \
        dumpNode(fp,P); \
        P = P->next; \
    }}\

#define LLSORT(PREFIX,T,L) \
T * PREFIX##sort(T * L) \
{ \
    T *p, *q, *e, *tail; int insize, nmerges, psize, qsize, i; \
    if (!L) return NULL; \
    insize = 1; \
    while (1) { \
        p = L;\
        L = NULL; tail = NULL; \
        nmerges = 0;\
        while (p) {\
            nmerges++; \
            q = p;   \
            psize = 0; \
            for (i = 0; i < insize; i++) { \
                psize++; \
                q = q->next; \
                if (!q)break; }\
            qsize = insize;\
            while (psize > 0 || (qsize > 0 && q)) { \
                if (psize == 0) { \
                    e = q;\
                    q = q->next;\
                    qsize--;\
                } else if (qsize == 0 || !q) {\
                    /* q is empty; e must come from p. */\
                    e = p;\
                    p = p->next;\
                    psize--;\
                } else if (cmp(p, q) <= 0) {\
                    /* First element of p is lower (or same);\
                     * e must come from p. */\
                    e = p;\
                    p = p->next;\
                    psize--;\
                } else {\
                    /* First element of q is lower; e must come from q. */\
                    e = q;\
                    q = q->next;\
                    qsize--;\
                }\
                /* add the next element to the merged list */\
                if (tail) {\
                    tail->next = e;\
                } else {\
                    L = e;\
                }\
                tail = e;\
            }\
            /* now p has stepped `insize' places along, and q has too */\
            p = q;\
        }\
        tail->next = NULL;\
        /* If we have done only one merge, we're finished. */\
        if (nmerges <= 1)       /* allow for nmerges==0, the empty list case */\
            return L;\
        /* Otherwise repeat, merging lists twice the size */\
        insize *= 2;\
    }\
}\

#define LLCODE(PREFIX,T) \
LLCNT(PREFIX,T, p) \
LLREMOVE(PREFIX,T, p) \
LLDESTROY(PREFIX,T, p) \
LLDUMP(PREFIX,T, p) \
LLSORT(PREFIX,T, list) \

/* 
  vim: ft=c sw=4 ts=4 expandtab
*/

#endif
