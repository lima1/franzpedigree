#ifndef DAG_H
#define DAG_H

/*
 * $Id: dag.h 1885 2010-01-25 15:31:17Z markus $
 *
 * - A minimalistic but very fast and well tested graph library for 
 *   directed, acyclic graphs (DAGS)
 * - Adjacency list based
 * - More features:  
 *    * supports fast cycle detection (depth first search and transitive 
 *      closure)
 *    * max indegree (2 for pedigrees)
 *    * for max indegree 2, we use a 2xN matrix to avoid malloc (recompile
 *      with ADJLIST to turn this off) 
 *    * dynamic transitive closure (beta)
 *    * connected components
 *    * graphviz output (quick and dirty)
 *    * NO support for adding or removing vertices!
 *
 * IMPORTANT: Compile production code with NDEBUG. Otherwise it will be 
 *            *really* slow...
 *
 * Mostly stolen from Algorithms in C, Part 5. Robert Sedgewick.
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

#ifndef LINESIZE
#define LINESIZE 65536
#endif
#define LABELSIZE 10

/* store either in a Nx2 matrix or in a linked list
 * the former is faster for pedigrees, the latter more flexible */
#ifndef ADJLIST
#define STORAGEARRAY 1
#endif
#include "listdag.h"

typedef struct {
    int v;
    int w;
} Edge;
Edge EDGE(int, int);

/* for the very ad hoc graphviz output */
typedef struct {
    int  num_cluster;
    int *num_vertices_in_cluster;
    char **labels;
    int  **ids;
} GraphvizCluster;

typedef struct {
    char shape[10];
    char color[10];
    char fillcolor[10];
    double height;
    double width;
} GraphvizVertex;

typedef struct dag *Dag;


/* initializes a Dag with n vertices. if we allow cycles, we might not have a
 * DAG anymore and the algorithms that assume a DAG won't give the correct
 * results anymore...  */
Dag  DAGinit(int n, int max_indegree, bool allow_cycles); 

void DAGdestroy(Dag);                   /* frees all malloc'd space        */
Dag  DAGcopy(Dag);                      /* clones a Dag (vertices and 
                                           edges, no other attributes to 
                                           make it as fast as possible     */
void DAGcopyE(Dag source, Dag dest);    /* an even faster variant that removes 
                                          dest's edges and then copies 
                                          the edges from source to dest */

/* edges */
bool DAGinsertE(Dag, Edge);             /* inserts the specified edge, 
                                           checks for cycles (with dfs or
                                           TC if the later is up2date) and 
                                           max indegree                    */
void DAGremoveE(Dag, Edge);             /* remove the edge                 */ 
bool DAGhasE(Dag, Edge);                /* has D such an edge?             */
void DAGremoveAllE(Dag);                /* removes all edges               */
int  DAGremoveIncomingE(Dag, int, Edge[]); /* removes all incoming edges of
                                              the vertex in O(1)           */
int  DAGincomingE(Dag, int, Edge[]);    /* returns all incoming edges of 
                                           the vertex in O(1)              */
int  DAGoutgoingE(Dag, int, Edge[]);    /* returns all outgoing edges of 
                                           the vertex in O(n)              */

int  DAGedges(Edge[], Dag);             
int  DAGgetEdgeNumber(Dag);             /* returns the number of edges in 
                                           O(1)                            */

/* vertices */
int  DAGgetVertexNumber(Dag);             /* returns the number of vertices in 
                                           O(1)                            */
int  DAGindegree(Dag, int);             /* returns deg+ in O(1)            */
int  DAGoutdegree(Dag, int);            /* returns deg- in O(n)            */
void DAGlabelV(Dag, int, char *,size_t);
void DAGgraphvizV(Dag, int, GraphvizVertex);
int  DAGgetMaxIndegree(Dag);            /* what is the max. indegree of the
                                           DAG (-1 = no limit)             */
int  DAGcalcDepthV(Dag);                /* calcs the depth of all vertices in 
                                           O(V), returns the max.depth     */
int  DAGgetDepthV(Dag, int);            /* returns the depth of the vertex 
                                           in O(1) after DAGcalcDepth
                                        */

/* cycle detection */
bool DAGisAcyclic(Dag);                 /* Is Dag a Dag?                   */ 
bool DAGcycleE(Dag, Edge);              /* Would Edge introduce a cycle?   */
bool DAGgetAllowCycles(Dag);            /* Do we allow cycles 
                                           (we shouldn't do that)          */

/* transitive closure. ASSUMES A DAG for the dynamic programming 
 * algorithm
 * Compile with NDEBUG! (the assertions calculate a new TC matrix and 
 * compare it with the updated one                                         */
void DAGtc(Dag);                        /* calculate TC V*V matrix 
                                           O(V*(V+X))                      */
void DAGtcUpdateAdd(Dag,Edge);          /* update TC matrix after adding   */ 
void DAGtcUpdateRemove(Dag,Edge*,int);  /* update TC matrix after removing */ 

bool DAGgetTcUpToDate(Dag);             /* is TC matrix up2date?           */
bool DAGreach(Dag, int, int);           /* is first vertex ancestor of
                                           the second vertex?              */
bool DAGreachNoTC(Dag, int, int);       /* the same, but does not require a TC */

/* writes for every vertex the component id in the specified array
 * (min length is number vertices) and returns the number of connected
 * components */
int DAGconnectedComponents(Dag, int *component_ids); 

/* a simple distance measure: simply calculate the number of different
 * vertices and normalize by the number of vertices. Requires same number 
 * of vertices                                                             */
double DAGcalcDistance(Dag,Dag);

/* return 1 when the two Dags are equal */
bool DAGareEqual(Dag,Dag);

/* output */
void DAGdump(FILE *,Dag);               /* output the DAG                  */
void DAGtcDump(FILE *,Dag);             /* output the transitive closure   */
void DAGgraphviz(Dag, char *, GraphvizCluster, bool married_vertex); 
void DAGsearch(Dag);


void insertE(Dag,Edge);                 /* insert edge without any checks  */

/* 
  vim: ft=c sw=4 ts=4 expandtab
*/

#endif
