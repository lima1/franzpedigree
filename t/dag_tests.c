/*
 * $Id: dag_tests.c 1698 2009-06-22 13:22:32Z markus $
 */

#include "macros.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tap.h"
#include "dag.h"
#include "pedigree.h"
#include "global.h"
#include "prob.h"
#include "options.h"
#include "dataio.h"

static char labels[7][10] = { "Grampa","Homer","Bart", "Lisa", "Maggie", "Marge", "Flanders" };
static Edge e[2];

DATA Data;
OPTIONS Options;
PROBS Probs;

static void 
homer_marge_flanders_tc_tests(Dag D) {
    /* 14 tests */
    int i;
    ok(DAGreach(D,1,5) == 0, "Homer no ancestor of Marge");
    ok(DAGreachNoTC(D,1,5) == 0, "Homer no ancestor of Marge");
    for (i=0;i<7;i++) {
        ok(DAGreach(D,6,i) == 0, "Flanders no ancestor");
        ok(DAGreachNoTC(D,6,i) == 0, "Flanders no ancestor");
    }    
    for (i=2;i<5;i++) {
        ok(DAGreach(D,1,i) == 1, "Homer ancestor of Bart-Maggie");
        ok(DAGreachNoTC(D,1,i) == 1, "Homer ancestor of Bart-Maggie");
    }    
    for (i=2;i<5;i++) {
        ok(DAGreach(D,5,i) == 1, "Marge ancestor of Bart-Maggie");
        ok(DAGreach(D,5,i) == 1, "Marge ancestor of Bart-Maggie");
    }    
}

static void 
all_tc_tests(Dag D) {
    ok(DAGreach(D,0,2) == 1, "Grampa ancestor of Bart");
    ok(DAGreachNoTC(D,0,2) == 1, "Grampa ancestor of Bart");
    ok(DAGreachNoTC(D,0,5) == 0, "Grampa no ancestor of Marge");
    ok(DAGreachNoTC(D,0,6) == 0, "Grampa no ancestor of Flanders");
    homer_marge_flanders_tc_tests(D);
}


static void 
triple_tests(Dag D) {
    int i,ret;
    Edge outgoing[14];
    for (i=0;i<7;i++)    
        DAGlabelV(D,i,labels[i], 10);
    
    ret = DAGoutgoingE(D,1,outgoing);
    ok(ret == 0, "Homer no children yet");
    ok(PEDIGREEinsertParents(D,2,1,5), "Bart child of Homer and Marge");
    ret = DAGoutgoingE(D,1,outgoing);
    ok(ret == 1, "Homer one child yet");
    ok(outgoing[0].v == 1 && outgoing[0].w == 2, "Bart child of Homer");
    ok(PEDIGREEinsertParents(D,3,1,5), "Lisa child of Homer and Marge");
    ret = DAGoutgoingE(D,1,outgoing);
    ok(ret == 2, "Homer two children yet");
    ok(outgoing[0].v == 1 && outgoing[0].w == 2, "Bart child of Homer");
    ok(outgoing[1].v == 1 && outgoing[1].w == 3, "Lisa child of Homer");
    ok(PEDIGREEinsertParents(D,4,1,5), "Maggie child of Homer and Marge");
    ret = DAGoutgoingE(D,1,outgoing);
    ok(ret == 3, "Homer three children");
    ok(outgoing[0].v == 1 && outgoing[0].w == 2, "Bart child of Homer");
    ok(outgoing[1].v == 1 && outgoing[1].w == 3, "Lisa child of Homer");
    ok(outgoing[2].v == 1 && outgoing[2].w == 4, "Maggie child of Homer");
    ret = DAGoutdegree(D,1);
    ok(ret == 3, "Homer three children");
    ret = DAGoutgoingE(D,5,outgoing);
    ok(ret == 3, "Marge three children");
    ok(!PEDIGREEinsertParents(D,2,5,6), "Bart not child of Flanders and Marge");
    ret = DAGoutgoingE(D,0,outgoing);
    ok(ret == 0, "Grampa no children yet");
    ok(PEDIGREEinsertParents(D,1,0,6), "Homer child of Flanders and Grampa");
    ret = DAGoutgoingE(D,0,outgoing);
    ok(ret == 1, "Grampa one child");
    ok(outgoing[0].v == 0 && outgoing[0].w == 1, "Homer child of Grampa");
    PEDIGREEremoveParents(D,1, e);
    ret = DAGoutgoingE(D,0,outgoing);
    ok(ret == 0, "Grampa no children yet");
}

int
main(void)
{
    Dag D, Dcloned, G, P1;
    int i, ret, component_ids[7];
    double retf;
#ifdef ADJLIST
    plan_tests(278);
#else
    plan_tests(272);
#endif    
    D = DAGinit(7,2,0);
    for (i=0;i<7;i++)    
        DAGlabelV(D,i,labels[i], 10);

    ret = DAGconnectedComponents(D, component_ids); 
    ok(ret == 7, "num. components in empty graph ok");
    for (i=0;i<7;i++)    
        ok(component_ids[i] == i, "component id in empty graph ok");

    ret = DAGinsertE(D,EDGE(0,1));
    ok(ret == 1, "Grampa father of Homer");
    /*ok(DAGreach(D,0,1),"DAGreach works without TC");*/
    
    ret = DAGconnectedComponents(D, component_ids); 
    ok(ret == 6, "num. components in (nearly) empty graph ok");
    ok(component_ids[0] == component_ids[1], "component ids of homer and grampa identical");
    for (i=2;i<7;i++)    
        ok(component_ids[i] == i-1, "component id in (nearly) empty graph ok");

    ret = DAGisAcyclic(D);
    ok(ret == 1, "is acyclic");

    ret = DAGinsertE(D,EDGE(1,0));
    ok(ret == 0, "Homer can't be father of Grampa");
/*    ok(!DAGreach(D,1,0),"DAGreach works without TC");*/

    ret = DAGisAcyclic(D);
    ok(ret == 1, "still acyclic");

    ret = DAGinsertE(D,EDGE(1,2));
    ok(ret == 1, "Homer father of Bart");
/*    ok(DAGreach(D,0,2),"DAGreach works without TC");*/

    ret = DAGconnectedComponents(D, component_ids); 
    ok(ret == 5, "num. components in graph ok");
    ok(component_ids[0] == component_ids[1], "component ids of homer and grampa identical");
    ok(component_ids[1] == component_ids[2], "component ids of homer and bart identical");

    for (i=3;i<7;i++)    
        ok(component_ids[i] == i-2, "component id in (nearly) empty graph ok");

    ret = DAGinsertE(D,EDGE(2,0));
    ok(ret == 0, "Bart can't be father of Grampa");

    ok(DAGhasE(D,EDGE(2,0)) == 0, "Bart can't be father of Grampa");

    ret = DAGinsertE(D,EDGE(1,3));
    ok(ret == 1, "Homer father of Lisa");
    ok(DAGhasE(D,EDGE(1,3)), "Homer father of Lisa");

    ret = DAGinsertE(D,EDGE(1,4));
    ok(ret == 1, "Homer father of Maggie");
    ok(DAGhasE(D,EDGE(1,4)), "Homer father of Maggie");
    ok(!DAGhasE(D,EDGE(1,1)), "Homer not father of Homer");
    ret = DAGconnectedComponents(D, component_ids); 
    ok(ret == 3, "num. components in graph ok");

    ret = DAGinsertE(D,EDGE(5,2));
    ok(ret == 1, "Marge mother of Bart");

    ret = DAGconnectedComponents(D, component_ids); 
    ok(ret == 2, "num. components in graph ok");
    ok(component_ids[0] == component_ids[5], "component id grampa and marge identical");
    ok(component_ids[5] != component_ids[6], "component id flanders and marge not identical");
    ok(component_ids[6] == 1, "component id flanders correct");

    ret = DAGinsertE(D,EDGE(5,3));
    ok(ret == 1, "Marge mother of Lisa");

    ret = DAGinsertE(D,EDGE(5,4));
    ok(ret == 1, "Marge mother of Lisa");
  
    ret = DAGgetEdgeNumber(D);
    ok(ret == 7, "7 edges so far");
    
    ret = DAGincomingE(D,2,e);
    ok(ret == 2, "Bart has two parents");
    ok(e[0].w == 2 && e[1].w == 2, "Bart is child of its parents");
    ok((e[0].v == 1 && e[1].v == 5) || (e[0].v == 5 && e[1].v == 1), "Bart is child of its parents");

    Dcloned = DAGcopy(D);
    ok(DAGgetEdgeNumber(Dcloned) == 7, "DAGcopy: edge number correct");
    ok(DAGareEqual(Dcloned, D) == 1, "DAGcopy: are equal after cloning");
    ret = DAGcalcDepthV(D); 
    ok(ret == 2,"max depth is 2");
    ok(DAGgetDepthV(D,0) == 0, "depth of Grampa 0");
    ok(DAGgetDepthV(D,5) == 0, "depth of Marge 0");
    ok(DAGgetDepthV(D,6) == 0, "depth of Flanders 0");
    ok(DAGgetDepthV(D,1) == 1, "depth of Homer 1");
    ok(DAGgetDepthV(D,2) == 2, "depth of Bart 2");
    ok(DAGgetDepthV(D,3) == 2, "depth of Lisa 2");
    ok(DAGgetDepthV(D,4) == 2, "depth of Maggie 2");
   
    ok(CMP_DBL(DAGcalcDistance(Dcloned, D), 0.), "distance zero");
    ok(CMP_DBL(DAGcalcDistance(D, Dcloned), 0.), "distance zero");
    /* test remove of incoming edges */
    DAGremoveIncomingE(D,2,e);
    ok(DAGareEqual(Dcloned, D) == 0, "DAGareEqual: not equal anymore");
    ok(fabs(DAGcalcDistance(Dcloned, D) - 2./(7.+5.)) < EPSILON, "distance correct");
    ok(fabs(DAGcalcDistance(D, Dcloned) - 2./(7.+5.)) < EPSILON, "distance correct");
    ret = DAGgetEdgeNumber(D);
    ok(ret == 5, "5 edges now (w/o Bart)");
    ret = DAGinsertE(D,EDGE(1,2));
    ok(ret == 1, "Homer father of Bart");
    ok(fabs(DAGcalcDistance(Dcloned, D) - 1./(7.+6.)) < EPSILON, "distance correct");
    ok(fabs(DAGcalcDistance(D, Dcloned) - 1./(7.+6.)) < EPSILON, "distance correct");
    ret = DAGinsertE(D,EDGE(5,2));
    ok(ret == 1, "Marge mother of Bart");
    ok(fabs(DAGcalcDistance(Dcloned, D) - 0) < EPSILON, "distance correct");
    ok(fabs(DAGcalcDistance(D, Dcloned) - 0) < EPSILON, "distance correct");
    
    DAGremoveE(Dcloned, EDGE(5,2));
    ok(!DAGareEqual(D, Dcloned), "dags not equal anymore");
    DAGcopyE(D, Dcloned);
    ok(DAGareEqual(D, Dcloned), "dags again equal");

    /* transitive-closure tests */
    DAGtc(D);
    all_tc_tests(D);
    DAGtc(Dcloned);
    all_tc_tests(Dcloned);

    /* update tests */
    e[0] = EDGE(0,1);
    DAGremoveE(D,e[0]);
    ret = DAGgetEdgeNumber(D);
    ok(ret == 6, "6 edges now");
    DAGtcUpdateRemove(D,e,1);
    for (i=0;i<5;i++)
        ok(DAGreach(D,0,i) == 0, "Grampa not ancestor of Homer-Maggie");
    homer_marge_flanders_tc_tests(D);
    DAGinsertE(D,e[0]);
    DAGtcUpdateAdd(D,e[0]);
    all_tc_tests(D);

    ret = DAGinsertE(D,EDGE(1,0));
    ok(ret == 0, "Homer can't be father of Grampa");
    ret = DAGinsertE(D,EDGE(2,0));
    ok(ret == 0, "Bart can't be father of Grampa");
    ret = DAGinsertE(D,EDGE(2,1));
    ok(ret == 0, "Bart can't be father of Homer");

    ret = DAGisAcyclic(D);
    ok(ret == 1, "still acyclic");
    
    /* small directed cycle test */
    DAGremoveAllE(D);
    ret = DAGgetEdgeNumber(D);
    ok(ret == 0, "no edges anymore");
    DAGtc(D);
    ret = DAGinsertE(D,EDGE(0,1));
    DAGtcUpdateAdd(D,EDGE(0,1));
    ok(ret == 1, "Grampa father of Homer");
    DAGinsertE(D,EDGE(1,3));
    DAGtcUpdateAdd(D,EDGE(1,3));
    ok(ret == 1, "Homer father of Lisa");

    DAGinsertE(D,EDGE(1,2));
    DAGtcUpdateAdd(D,EDGE(1,2));

    ok(ret == 1, "Homer father of Bart");
    ok(DAGreach(D,1,2), "Homer ancestor of Bart");

    DAGinsertE(D,EDGE(3,2));
    DAGtcUpdateAdd(D,EDGE(3,2));
    ok(ret == 1, "Lisa mother of Bart");

    DAGremoveE(D,EDGE(3,2));
    e[0] = EDGE(3,2);
    DAGtcUpdateRemove(D,e,1);
    ok(DAGreach(D,1,2), "updating works after removing small directed cycle");
    ok(!DAGreach(D,3,2), "updating works after removing small directed cycle");


    /* TODO: should fail if sexdata avaiable */
    ret = DAGinsertE(D,EDGE(0,3));
    ok(ret == 1, "Grampa father of Lisa");

    //DAGdump(stderr,D);
   

    /* triples */
    P1 = DAGinit(7,2,0);
    DAGtc(P1);
    //DAGtcDump(stderr,P1);
    triple_tests(P1);
    DAGtc(P1);
    homer_marge_flanders_tc_tests(P1);
    DAGdestroy(P1);

    ok(DAGgetAllowCycles(D) == 0,  "we check for cycles");
    ok(DAGgetMaxIndegree(D) == 2,  "two parents");

#ifdef ADJLIST    
    /* run tests in Graph-mode */
    G = DAGinit(7,-1,1);
    for (i=0;i<7;i++)    
        DAGlabelV(G,i,labels[i], 10);

    ret = DAGinsertE(G,EDGE(0,1));
    ok(ret == 1, "Grampa father of Homer");
    ret = DAGisAcyclic(G);
    ok(ret == 1, "G is acyclic");

    ret = DAGinsertE(G,EDGE(1,0));
    ok(ret == 1, "Homer can now be father of Grampa");

    ret = DAGisAcyclic(G);
    ok(ret == 0, "not a DAG anymore");
    ok(DAGgetAllowCycles(G) == 1,  "no check for cycles");
    ok(DAGgetMaxIndegree(G) == -1, "no check for indegree");
    DAGdestroy(G);
#endif

    DAGremoveAllE(D);
    retf = PEDIGREEcalcSelfingRate(D);
    ok(CMP_DBL(retf, 0.), "selfingrate correct");
    ret = DAGinsertE(D,EDGE(0,1));
    ret = DAGinsertE(D,EDGE(0,1));
    retf = PEDIGREEcalcSelfingRate(D);
    ok(CMP_DBL(retf, 1.0), "selfingrate correct");
    
    ret = DAGinsertE(D,EDGE(1,2));
    ret = DAGinsertE(D,EDGE(5,2));
    retf = PEDIGREEcalcSelfingRate(D);
    ok(CMP_DBL(retf, 0.5), "selfingrate correct");
    DAGdestroy(Dcloned);
    //DAGtcDump(stderr,D);
    DAGdestroy(D);
    G = DAGinit(8,2,0);
    D = DAGinit(7,2,0);
    ok(DAGareEqual(D,G) == 0 && DAGareEqual(G,D) == 0, "empty graphs with diff. v not equal"); 
    DAGdestroy(G);
    G = DAGinit(7,2,0);
    ok(DAGareEqual(D,G) == 1 && DAGareEqual(G,D) == 1, "empty graphs with equal. v equal"); 
    DAGinsertE(D, EDGE(1,2));
    DAGinsertE(D, EDGE(5,2));
    DAGinsertE(G, EDGE(1,2));
    DAGinsertE(G, EDGE(6,2));
    ok(DAGareEqual(D,G) == 0 && DAGareEqual(G,D) == 0, "graphs not equal"); 
    DAGremoveE(G, EDGE(6,2));
    DAGinsertE(G, EDGE(5,2));
    ok(DAGareEqual(D,G) == 1 && DAGareEqual(G,D) == 1, "graphs equal"); 
    DAGdestroy(D);
    D = DAGinit(8,2,0);
    ok(DAGareEqual(D,G) == 0 && DAGareEqual(G,D) == 0, "graphs not equal (different no. vertices)");
    DAGdestroy(D);
    DAGdestroy(G);
    return exit_status();       /* Return the correct exit code */
}

/* 
  vim: ft=c sw=4 ts=4 expandtab
*/
