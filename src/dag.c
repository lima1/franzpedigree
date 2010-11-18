/*
 * $Id: dag.c 1885 2010-01-25 15:31:17Z markus $
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

#include "macros.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dag.h"
#include <assert.h>

#ifdef STORAGEARRAY
  /* returns only valid ids (> 0) */
   #define INITERATION(D,t,i) for (t##t = 0; t##t < 2; t##t++) { \
       t = D->adj[i][t##t];\
       if (t < 0) continue;
   #define CURRENTIN(t) t
   #define DEFITER(t) int t,t##t
#else       
   #define INITERATION(D,t,i) for (t = D->adj[i]; t != NULL; t = t->next) { 
   #define CURRENTIN(t) t->v
   #define DEFITER(t) LISTDAGNODE *t
#endif

struct dag
{
    int V;             /* the number of vertices                             */
    int E;             /* the number of edges                                */
    int max_indegree;  /* max. indegree (otherwise DAGinsertE will fail)     
                          a value smaller 0 turns this check off             */
    bool allow_cycles;  /* when 1, DAGinsertE won't check for directed cycles */
#ifdef STORAGEARRAY    
    int **adj;         /* the adjacency list                                 */
#else    
    LISTDAGNODE **adj; /* the adjacency list                                 */
#endif    
    bool **tc;          /* transitive closure matrix                          */
    bool tc_up_to_date; /* flag that is set to 0 after edge ins. or del.      */
    char **vlabels;    /* array of vertex labels, default label is vertex id */
    float *vsizes;     /* array of vertex sizes (default 1) for graphviz     */
    GraphvizVertex *vgraphviz;  /* array of graphviz attributes                       */
    int *vdepths;      /* array of vertex depths                             */ 
    int _maxdepth;
    /* some internal variables */
    int _cnt;
    int _cntP;
    int *_pre;
    int *_post;
};

static int depth;           /* debugging only: output of dfs trace in show() */

Edge
EDGE(int v, int w)
{
    Edge e;

    e.v = v;
    e.w = w;
    return e;
}

static bool
isAncestor(Dag D, Edge e)
{
    DEFITER(t);

    D->_pre[e.w] = 1;
    INITERATION(D,t,e.w)  
        if (CURRENTIN(t) == e.v) {
            return true;
        }
        if (D->_pre[CURRENTIN(t)] == -1) {
            if (isAncestor(D, EDGE(e.v, CURRENTIN(t))))
                return true;
        }
    }
    return false;
}

static bool
hasbacklinkR(Dag D, Edge e)
{
    bool ret;
    DEFITER(t);

    D->_pre[e.w] = D->_cnt++;
    INITERATION(D,t,e.w)  
        if (D->_pre[CURRENTIN(t)] == -1) {
            ret = hasbacklinkR(D, EDGE(e.w, CURRENTIN(t)));
            if (ret)
                return true;
        } else {
            if (D->_post[CURRENTIN(t)] == -1)
                return true;
        }
    }
    D->_post[e.w] = D->_cntP++;
    return false;
}

static void
addDepthR(Dag D, int v)
{
    int max = -1;
    DEFITER(t);

    /* mark vertex as done */
    D->_pre[v] = 1;

    /* depth of vertex is the max. depth of its children + 1 */
    max = -1;
    INITERATION(D,t,v)    
        if (D->_pre[CURRENTIN(t)] == -1) 
            addDepthR(D, CURRENTIN(t));
        
        if (D->vdepths[CURRENTIN(t)] > max)
            max = D->vdepths[CURRENTIN(t)];
    }
    D->vdepths[v] = max+1;
    if (D->vdepths[v] > D->_maxdepth)
        D->_maxdepth = D->vdepths[v];

}

bool
DAGisAcyclic(Dag D)
{
    int v;

    D->_cnt = 0;
    D->_cntP = 0;
    for (v = 0; v < D->V; v++) {
        D->_pre[v] = -1;
        D->_post[v] = -1;
    }
    for (v = 0; v < D->V; v++)
        if (D->_pre[v] == -1) 
            if (hasbacklinkR(D, EDGE(v, v))) 
                return false;

    return true;
}

int
DAGcalcDepthV(Dag D)
{
    int v;

    for (v = 0; v < D->V; v++) 
        D->_pre[v] = -1;

    for (v = 0; v < D->V; v++)
        if (D->_pre[v] == -1) 
            addDepthR(D, v);

    return D->_maxdepth;
}

int
DAGgetDepthV(Dag D, int v) {
    return D->vdepths[v];
}

Dag
DAGinit(int V, int max_indegree, bool allow_cycles)
{
    int v;
    Dag D = malloc(sizeof *D);
    if (D == NULL) 
        FATAL("malloc failed");

    D->V = V;
    D->E = 0;
    D->_maxdepth = 0;
#ifdef STORAGEARRAY
    if (max_indegree < 0 || max_indegree > 2 )
        FATAL("max_indegree > 2 or undefined. recompile with -DADJLIST");
    MAKE2DINT(D->adj, D->V, 2,v)
#else        
    MALLOC(D->adj, LISTDAGNODE *, V);
#endif    

    MALLOC(D->vlabels, char *, V);
    MALLOC(D->vgraphviz, GraphvizVertex, V);

    MALLOC(D->vsizes, float, V);
    MAKE1DINT(D->vdepths, V);
    MAKE1DINT(D->_pre, V);
    MAKE1DINT(D->_post, V);

    MAKE2DBOOL(D->tc, D->V, D->V, v);

    D->max_indegree = max_indegree;
/*    MALLOC(incoming,Edge,max_indegree);*/
    D->allow_cycles = allow_cycles;
    D->tc_up_to_date = false;

    for (v = 0; v < V; v++) {
#ifdef STORAGEARRAY
        D->adj[v][0] = -1;
        D->adj[v][1] = -1;
#else
        D->adj[v] = NULL;
#endif        

        MALLOC(D->vlabels[v], char, LABELSIZE+1);
        snprintf(D->vlabels[v], LABELSIZE * sizeof(char), "%i", v);

        D->vgraphviz[v].height = 0.;
        D->vgraphviz[v].width = 0.;
        strcpy(D->vgraphviz[v].shape, "diamond");
    }
    return D;
}

bool
DAGgetTcUpToDate(Dag D)
{
    return D->tc_up_to_date;
}

bool 
DAGgetAllowCycles(Dag D)
{
    return D->allow_cycles;
}

int
DAGgetMaxIndegree(Dag D)
{
    return D->max_indegree;
}

void
DAGlabelV(Dag D, int i, char *s, size_t length)
{
    if (length > LABELSIZE)
        FATALINT("DAGlabelV()");
    strncpy(D->vlabels[i], s, length);
}

void
DAGgraphvizV(Dag D, int i, GraphvizVertex gv)
{
   strncpy(D->vgraphviz[i].shape,gv.shape,10);
}

void
insertE(Dag D, Edge e)
{
#ifdef STORAGEARRAY
    int i;
    /* use first free slot to store the edge */
    i = D->adj[e.w][0] < 0 ? 0 : 1;
    D->adj[e.w][i] = e.v;
#else    
    LISTDAGadd(&D->adj[e.w], e.v);
#endif    
    D->E++;
}

double
DAGcalcDistance(Dag D1, Dag D2)
{
    int i, sum_found, wrong = 0;
    bool found;
    DEFITER(t);
    DEFITER(u);

    if (D1->V != D2->V)
        FATAL("Number vertices not eqal.");
    for (i=0; i<D1->V;i++) {
        sum_found = 0;
        INITERATION(D1,t,i)
            found = false;
            INITERATION(D2,u,i)
                if (CURRENTIN(t) == CURRENTIN(u)) {
                    found = true;
                    break;    
                }
            }
            if (found) sum_found++;
        } 
        wrong += (DAGindegree(D1,i)+DAGindegree(D2,i)) - 2*sum_found;
    }
    return wrong / (double)( D1->E + D2->E);
} 

bool
DAGareEqual(Dag D1, Dag D2)
{
    int i;
    bool found;
    DEFITER(t);
    DEFITER(u);

    if (D1->V != D2->V || D1->E != D2->E)
        return false;

    for (i=0; i<D1->V;i++) {
        if (DAGindegree(D1,i)!=DAGindegree(D2,i))
            return false;
        INITERATION(D1,t,i)
            found = false;
            INITERATION(D2,u,i)
                if (CURRENTIN(t) == CURRENTIN(u)) {
                    found = true;
                    break;    
                }
            }
            if (!found) return false;
        } 
    }
    return true;
} 

bool
DAGinsertE(Dag D, Edge e)
{
    int w = e.w;

    /* we have a TC? then check for cycles in O(1) 
    if (D->tc_up_to_date && !D->allow_cycles && DAGreach(D, e.w, e.v))
        return false;
    */

    /* max indegree check */
    if (D->max_indegree > 0 && D->max_indegree <= DAGindegree(D, w))
        return false;

    insertE(D, e);

    D->tc_up_to_date = false;

    /* check for cycles with a dfs */
    if (!D->allow_cycles) {
        if (DAGcycleE(D, e)) {
            DAGremoveE(D, e);
            return false;
        }
    }
    return true;
}

void
DAGremoveE(Dag D, Edge e)
{
#ifdef STORAGEARRAY
    if (D->adj[e.w][0] == e.v)
        D->adj[e.w][0] = -1;
    else     
        D->adj[e.w][1] = -1;
#else    
    LISTDAGremove(LISTDAGsearch(&D->adj[e.w],e.v));
#endif    
    D->E--;
    D->tc_up_to_date = false;
}

bool
DAGhasE(Dag D, Edge e)
{
#ifdef STORAGEARRAY
    if (D->adj[e.w][0] == e.v)
        return true;
    if (D->adj[e.w][1] == e.v)
        return true;
#else    
    if (LISTDAGsearch(&D->adj[e.w],e.v))
        return true;
#endif    
    return false;
}

int  
DAGgetEdgeNumber(Dag D) 
{
    return D->E;
} 

int  
DAGgetVertexNumber(Dag D) 
{
    return D->V;
} 

int
DAGremoveIncomingE(Dag D, int v, Edge edges[])
{
    DEFITER(t);
    int i = 0;
    
    INITERATION(D,t,v)     
        edges[i].v = CURRENTIN(t);
        edges[i].w = v;
        i++;
    }
#ifdef STORAGEARRAY
    D->adj[v][0] = -1;
    D->adj[v][1] = -1;
#else    
    LISTDAGdestroy(&D->adj[v]);
    D->adj[v] = NULL;
#endif    
    D->E -= i;
    return i;
}

void
DAGremoveAllE(Dag D)
{
    int v;

    for (v = 0; v < D->V; v++) {
#ifdef STORAGEARRAY
        D->adj[v][0] = -1;
        D->adj[v][1] = -1;
#else    
        LISTDAGdestroy(&D->adj[v]);
        D->adj[v] = NULL;
#endif 
    }     
    D->E = 0;
}

int
DAGincomingE(Dag D, int v, Edge edges[])
{
    DEFITER(t);
    int i = 0;

    INITERATION(D,t,v)     
        edges[i].v = CURRENTIN(t);
        edges[i].w = v;
        i++;
    }
    return i;
}

int
DAGoutgoingE(Dag D, int v, Edge edges[])
{
    DEFITER(t);
    int i, cnt = 0;

    for (i=0; i<D->V;i++) {
        INITERATION(D,t,i)
            if (CURRENTIN(t) == v) {
                edges[cnt].v = v;
                edges[cnt].w = i;
                cnt++;
                break;
            }
        }
    }
    return cnt;
}

int
DAGoutdegree(Dag D, int v)
{
    DEFITER(t);
    int i, cnt = 0;

    for (i=0; i<D->V;i++) {
        INITERATION(D,t,i)
            if (CURRENTIN(t) == v) {
                cnt++;
                break;
            }
        }
    }
    return cnt;
}

int
DAGindegree(Dag D, int v)
{
    int indegree = 0;
    DEFITER(t);

    INITERATION(D,t,v)    
        indegree++;
    }
    return indegree;
}

int
DAGedges(Edge a[], Dag D)
{
    int v, E = 0;
    DEFITER(t);

    for (v = 0; v < D->V; v++)
         INITERATION(D,t,v)   
            a[E++] = EDGE(CURRENTIN(t), v);
         }
    return E;
}

static void
TCdfsR(Dag D, Edge e)
{
    DEFITER(t);
    int v = e.w, i;

    D->_pre[v] = D->_cnt++;

    INITERATION(D,t,v)    
        D->tc[CURRENTIN(t)][v] = true;
        if (D->_pre[CURRENTIN(t)] > D->_pre[v])
            continue;
        if (D->_pre[CURRENTIN(t)] == -1)
            TCdfsR(D, EDGE(CURRENTIN(t), CURRENTIN(t)));
        for (i = 0; i < D->V; i++)
            if (D->tc[i][CURRENTIN(t)]) {
                D->tc[i][v] = true;
            }
    }
}

#ifndef NDEBUG
/* Debug function that clones the transitive closure, recomputes the original
 * one and then compares both. Useful for assertions in DAGtcUpdate*.
 * Make sure that this is never executed in production code.                 */
static bool
updatingWorks(Dag D)
{
    int i, j;
    bool **tc_test;

    MAKE2DBOOL(tc_test, D->V, D->V, i);
    for (i = 0; i < D->V; i++)
        for (j = 0; j < D->V; j++)
            tc_test[i][j] = D->tc[i][j];
    /*DAGtcDump(stderr, D);*/
    DAGtc(D);
   /* DAGtcDump(stderr, D);*/
    for (i = 0; i < D->V; i++)
        for (j = 0; j < D->V; j++)
            if (tc_test[i][j] != D->tc[i][j]) {
                DAGdump(stderr, D);
                FREE2D(tc_test, D->V, i);
                return false;
            }

    FREE2D(tc_test, D->V, i);
    return true;
}
#endif

void
DAGtcUpdateAdd(Dag D, Edge e)
{
    int i, j;
    D->tc_up_to_date = true;
    /* we can skip this if this edge is already added */
    if (D->tc[e.v][e.w]) return;

/*    printf("Inserting: %i %i\n", e.v, e.w);*/

    for (i = 0; i < D->V; i++)
        if ((D->tc[i][e.v] && !D->tc[i][e.w]) || i == e.v) {
            for (j = 0; j < D->V; j++)
                if (D->tc[e.w][j] || j == e.w)
                    D->tc[i][j] = true;
        }    
    assert(updatingWorks(D));
}

void
DAGtcUpdateRemove(Dag D, Edge * e, int edge_cnt)
{
    int i, j, id;
    int toupdatecnt = 0;
    //GraphvizCluster cluster;

/*    cluster.num_cluster = 0;
    for (i = 0; i < edge_cnt; i++)
        printf("Deleting: %i %i\n", e[i].v, e[i].w);
    printf("TCDDDUUUMP\n");
    DAGtcDump(stderr, D);
    DAGdump(stderr, D);*/

    D->_cnt = 0;
    for (i = 0; i < D->V; i++) {
        D->_pre[i] = -1;
        for (j = 0; j < edge_cnt; j++) {
            if (i == e[j].v || i == e[j].w || D->tc[i][e[j].v]
                || D->tc[e[j].w][i]) {
                D->_post[toupdatecnt++] = i;
                break;
            }
        }
    }

    /*DAGgraphviz(D, "debug", cluster);
    DAGtcDump(stderr, D);*/

    for (i = 0; i < toupdatecnt; i++)
        for (j = 0; j < toupdatecnt; j++)
            D->tc[D->_post[i]][D->_post[j]] = false;

    for (i = 0; i < toupdatecnt; i++) {
        id = D->_post[i];
        /*printf("UPDATE: %i, %i, %i, PRE %i\n", i, id, toupdatecnt,
         * pre[id]);*/

        if (D->_pre[id] == -1) {
            /*for (t=D->adj[id]; t != NULL; t=t->next) 
               D->tc[t->v][id]=0; */
            TCdfsR(D, EDGE(id, id));
        }
    }
    D->tc_up_to_date = true;
    assert(updatingWorks(D));
}

/* calculates the transitive closure */
void
DAGtc(Dag D)
{
    int v, w;
    /* dynamic programming algorithm is DAG only */
    assert(DAGisAcyclic(D));
    D->_cnt = 0;
    for (v = 0; v < D->V; v++) {
        D->_pre[v] = -1;
        for (w = 0; w < D->V; w++)
            D->tc[v][w] = false;
    }
    for (v = 0; v < D->V; v++)
        if (D->_pre[v] == -1)
            TCdfsR(D, EDGE(v, v));
    D->tc_up_to_date = true;
}

/* dumps the transitive closure */
void
DAGtcDump(FILE * fp, Dag D)
{
    int i, j;

    for (i = 0; i < D->V; i++) {
        fprintf(fp, "%i:\n", i);
        for (j = 0; j < D->V; j++)
            if (D->tc[i][j])
                fprintf(fp, "\t[%i][%i]: 1\n", i, j);
    }
}

bool
DAGreach(Dag D, int s, int t)
{
    assert(DAGgetTcUpToDate(D));
    return D->tc[s][t];
}

bool 
DAGreachNoTC(Dag D, int s, int t)
{
    int i;
    for (i = 0; i < D->V; i++)
        D->_pre[i] = -1;
    return isAncestor(D, EDGE(s, t));
}

void
DAGdump(FILE * fp, Dag D)
{
    int i,j, *component_ids, components;
    DEFITER(t);
    DAGcalcDepthV(D);

    fprintf(fp, "%d vertices, %d edges\n", D->V, D->E);
    MAKE1DINT(component_ids, D->V);
    components = DAGconnectedComponents(D, component_ids);
    for (j = 0; j<components;j++) {
        fprintf(fp, "\nConnected Component %i:\n", j+1);
        for (i = 0; i < D->V; i++) {
            if (component_ids[i] != j)
                continue;
            fprintf(fp, "  %s (depth: %i):", D->vlabels[i], D->vdepths[i]);
            INITERATION(D,t,i)    
                fprintf(fp, " %s ", D->vlabels[CURRENTIN(t)]);
            }
            fprintf(fp, "\n");
        }
    }
    FREE1D(component_ids);
}

bool
DAGcycleE(Dag D, Edge e)
{
    int w;

    D->_cnt = 0;
    D->_cntP = 0;
    for (w = 0; w < D->V; w++) {
        D->_pre[w] = -1;
        D->_post[w] = -1;
    }
    return hasbacklinkR(D, e);
}

void
DAGgraphviz(Dag D, char *filename, GraphvizCluster cluster, bool married_vertex)
{
    /* quick n dirty graphviz export */
    /* maybe I should use gvc.h. dependencies, dependencies... */

    char line[LINESIZE];
    int i, j;
    float width;
    FILE *out = NULL;
    bool *already_connected;
    Edge *edges;
    
    FOPENW(out, filename);

    /* print graphviz header */
    fputs
        ("digraph GRAPH_0 {\n\tedge [ dir=none];\n\tgraph [ rankdir=TB ];\n",
         out);
    fputs
        ("\tratio=auto;\n\tmincross=2.0;\n",
         out);
    fputs
        ("\tnode [\n\t\tfontsize=11,\n\t\tfillcolor=white,\n\t\tstyle=filled,\n\n];\n",
         out);
    MAKE1DBOOL(already_connected, D->V);
    for (i = 0; i < D->V; i++) already_connected[i] = false;
    if (cluster.num_cluster > 0) {
        for (i = 0; i < cluster.num_cluster; i++) {
            snprintf(line, sizeof(line),"\tsubgraph cluster_%i {\n\t\tlabel = \"%s\"\n", i,
                    cluster.labels[i]);
            fputs(line, out);
            for (j = 0; j < cluster.num_vertices_in_cluster[i]; j++) {
                width = D->vsizes[cluster.ids[i][j]];
                snprintf(line,sizeof(line),
                        "\t\t%i [ label=\"%s\", shape=%s, width=%f, height=%f ]\n",
                        cluster.ids[i][j], D->vlabels[cluster.ids[i][j]], D->vgraphviz[cluster.ids[i][j]].shape,
                        width, width);
                fputs(line, out);
                if (DAGindegree(D,cluster.ids[i][j]) > 0) 
                    snprintf(line, sizeof(line),"\tmarr%i [shape=diamond,style=filled,label=\"\",height=.1,width=.1]\n", cluster.ids[i][j]);
                fputs(line, out);
            }
            fputs("\t}\n", out);
        }
    } else {
        for (i = 0; i < D->V; i++) {
            width = D->vsizes[i];
            snprintf(line, sizeof(line),
                    "\t\t%i [ label=\"%s\", shape=%s, width=%f, height=%f ]\n",
                    i, D->vlabels[i], D->vgraphviz[i].shape, width, width);
            fputs(line, out);
        }
    }

    MALLOC(edges, Edge, D->E);
    DAGedges(edges, D);
    
    for (i = 0; i < D->E; i++) {
        if (married_vertex) {
            snprintf(line, sizeof(line),"\t%i -> marr%i  [ color=\"#000000\" weight=1 ]\n", edges[i].v,
                    edges[i].w);
            fputs(line, out);
            if (!already_connected[edges[i].w]) {
                snprintf(line, sizeof(line),"\tmarr%i -> %i  [ color=\"#000000\" weight=2 ]\n", edges[i].w, edges[i].w);
                already_connected[edges[i].w] = true;
                fputs(line, out);
            }
        }
        else {
            snprintf(line, sizeof(line),"\t%i -> %i  [ color=\"#000000\" weight=1 ]\n", edges[i].v,
                    edges[i].w);
                fputs(line, out);
        }    
    }
    FREE(edges);
    FREE(already_connected);
    fputs("}\n", out);
    FCLOSE(out);
}

void
DAGdestroy(Dag D)
{
    int i;
    if (D != NULL) {
        /* remove edges */
#ifndef STORAGEARRAY
        for (i = 0; i < D->V; i++)
            LISTDAGdestroy(&D->adj[i]);
#endif
        FREE2D(D->adj, D->V, i);
        FREE2D(D->tc, D->V, i);

        /* remove vertices incl. attributes */
        FREE2D(D->vlabels, D->V, i);
        FREE(D->vgraphviz);
        FREE(D->vsizes);
        FREE(D->vdepths);
        FREE(D->_pre);
        FREE(D->_post);

        FREE(D);
    }
}


/* DEBUG: */

static void
show(char *msg, int v, int w)
{
    int i;

    for (i = 0; i < depth; i++)
        printf(" ");
    printf("%i-%i %s\n", v, w, msg);
}

static void
dfsR(Dag D, Edge e)
{
    DEFITER(t);
    int v = e.v;

    depth++;
    show("tree", e.v, e.w);

    D->_pre[e.v] = D->_cnt++;
    INITERATION(D,t,e.v)    
        if (D->_pre[CURRENTIN(t)] == -1)
            dfsR(D, EDGE(CURRENTIN(t), e.v));
        else {
            v = CURRENTIN(t);
            if (D->_post[v] == -1)
                show("back", e.v, v);
            else if (D->_pre[v] > D->_pre[e.w])
                show("down", e.w, v);
            else
                show("cross", e.w, v);
        }
    }
    D->_post[e.v] = D->_cntP++;
    depth--;
}

static void
setComponentID(Dag D, int v, int *component_ids, int id)
{    
    DEFITER(t);
    component_ids[v] = id;
    /* apply id to all ancestors with unset component id */
    INITERATION(D,t,v)
        if (component_ids[CURRENTIN(t)] != id)
            setComponentID(D,CURRENTIN(t),component_ids, id);
    }
}

static int
getComponentID(Dag D, int v, int *component_ids)
{
    /* TODO: this assumes a max indegree of 2 */
    int i, id[2] = { -1, -1}, cnt = 0, max;
    DEFITER(t);

    /* has a component id? */
    if (component_ids[v] >= 0)
        return component_ids[v];

    /* otherwise check the ancestors */
    INITERATION(D,t,v)
        id[cnt++] = getComponentID(D, CURRENTIN(t), component_ids);
    }
    if (id[0] < 0)
        return id[1];
    if (id[1] < 0)
        return id[0];

    max = MAX(id[0], id[1]);

    for (i=0; i<D->V; i++)
        if (component_ids[i] == max)
            component_ids[i] = MIN(id[0], id[1]);

    return MIN(id[0], id[1]);
}    

int 
DAGconnectedComponents(Dag D, int *component_ids)
{
    int i,j, id, next_id = 0;
    bool found;

    for (i=0; i< D->V; i++)
        component_ids[i] = -1;

    for (i=0; i< D->V; i++)
        if (component_ids[i] < 0) {
            /* has any ancestor already a component id? */
            id = getComponentID(D, i, component_ids);
            /* no? then take the next available ... */
            if (id < 0) {
                id = next_id;
                next_id++;
            }    
            /* and apply it to all recursively */
            setComponentID(D, i, component_ids, id);
        }
    
    /* the following code makes sure that ids are from
     * 0 to number components */
    for (i=0; i<next_id; i++) {
        found = false;
        for (j=0; j< D->V; j++) {
            if (component_ids[j] == i) {
                found = true;
                break;
            }
        }
        if (!found) {
            for (j=0; j< D->V; j++) 
                if (component_ids[j] > i) 
                    component_ids[j]--;
            next_id--;
            i--;
        }
    }
    return next_id;
}

void
DAGsearch(Dag D)
{
    int v;

    D->_cnt = 0;
    D->_cntP = 0;
    for (v = 0; v < D->V; v++) {
        D->_pre[v] = -1;
        D->_post[v] = -1;
    }
    for (v = 0; v < D->V; v++)
        if (D->_pre[v] == -1) {
            dfsR(D, EDGE(v, v));
        }
}

void 
DAGcopyE(Dag Ds, Dag Dd)
{
    int i;
    DEFITER(t);
    DAGremoveAllE(Dd);
    for (i=0; i<Ds->V; i++) {
        INITERATION(Ds,t,i)     
            insertE(Dd, EDGE(CURRENTIN(t),i));
        }
    }
    assert(DAGareEqual(Ds,Dd));
}

Dag
DAGcopy(Dag D) 
{
    int i;
    DEFITER(t);
    Dag cloned;
    cloned = DAGinit(D->V, D->max_indegree, D->allow_cycles);
    for (i=0; i<D->V; i++) {
        INITERATION(D,t,i)     
            insertE(cloned, EDGE(CURRENTIN(t),i));
        }
    }
    assert(DAGareEqual(D,cloned));
    return cloned;
}    
/* 
  vim: ft=c sw=4 ts=4 expandtab
*/
