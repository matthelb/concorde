/****************************************************************************/
/*                                                                          */
/*  This file is part of CONCORDE                                           */
/*                                                                          */
/*  (c) Copyright 1995--1999 by David Applegate, Robert Bixby,              */
/*  Vasek Chvatal, and William Cook                                         */
/*                                                                          */
/*  Permission is granted for academic research use.  For other uses,       */
/*  contact the authors for licensing options.                              */
/*                                                                          */
/*  Use at your own risk.  We make no guarantees about the                  */
/*  correctness or usefulness of this code.                                 */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/*                  Exact Seperation of Blossoms (Padberg-Rao)              */
/*                                                                          */
/*                              TSP CODE                                    */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: Spring 1989 (Bico)                                                */
/*        Modified January 11, 1999 (Bico)                                  */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCtsp_exactblossom (CCtsp_lpcut_in **cuts, int *cutcount,           */
/*      int ncount, int ecount, int *elist, double *x,                      */
/*      CCrandstate *rstate)                                                */
/*       RUNS Padberg-Rao to seperate over blossom inequalities.            */
/*        -cuts (new cutting plans will be added to front of this list)     */
/*        -cutcount will return the number of new cuts found (can be NULL)  */
/*        -ncount is the number of nodes                                    */
/*        -ecount is the number of edges                                    */
/*        -elist is the edge list in node node format                       */
/*        -x is an lp solution vector                                       */
/*    NOTES:                                                                */
/*      The exactblossom  code was written very early in our TSP project.   */
/*      In January 1999 it was updated to fit into the current concorde,    */
/*      but the guts of the code are still written in our old sloppy        */
/*      style.  This is a good candidate for a rewrite (big speedups are    */
/*      probably possible without too much effort).                         */
/*                                                                          */
/*  int CCtsp_fastblossom (CCtsp_lpcut_in **cuts, int *cutcount,            */
/*      int ncount, int ecount, int *elist, double *x)                      */
/*    FINDS blossoms by looking at 0 < x < 1 graph for connected comps      */
/*     meeting an odd number of 1 edges.                                    */
/*                                                                          */
/*  int CCtsp_ghfastblossom (CCtsp_lpcut_in **cuts, int *cutcount,          */
/*      int ncount, int ecount, int *elist, double *x)                      */
/*    FINDS blossoms using a heuristic described by Groetschel and          */
/*     Holland. It works with the 0 < x < 1-EPS (with EPS = .3) graph,      */
/*     builds components, and picks a greedy set of teeth.                  */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/*    The code is set up to allow blossoms with single teeth.  To forbid    */
/*    these, undef ALLOW_SINGLE_TEETH.                                      */
/*                                                                          */
/****************************************************************************/


#include "machdefs.h"
#include "util.h"
#include "cut.h"
#include "tsp.h"

#define ALLOW_SINGLE_TEETH
#define BLOTOLERANCE .01
#define OTHEREND(e,n) ((e)->ends[0] == (n) ? (e)->ends[1] \
                                           : (e)->ends[0])

typedef struct edge {
    struct node    *ends[2];
    struct node    *splitter;
    double          x;
    int             magiclabel;
    struct edge    *next;
} edge;

typedef struct edgeptr {
    struct edge    *this;
    struct edgeptr *next;
} edgeptr;

typedef struct node {
    edgeptr        *adj;
    int             magiclabel;
    struct node    *next, *prev;
    struct node    *oddnode;
    edge           *pe;
    int             mark;
    int             num;
    int             name;
} node;

typedef struct nodeptr {
    struct node    *this;
    struct nodeptr *next;
} nodeptr;

typedef struct graph {
    int             ncount;
    node           *nodelist;
    int             ecount;
    edge           *edgelist;
    node           *pseudonodelist;
    edge           *pseudoedgelist;
    int             magicnum;
    node            pseudonodedummy;
    edge            pseudoedgedummy;
    CCptrworld      edge_world;
    CCptrworld      edgeptr_world;
    CCptrworld      node_world;
    CCptrworld      nodeptr_world;
} graph;

typedef struct toothobj {
    int  in;
    int  out;
} toothobj;

CC_PTRWORLD_ROUTINES (edge, edgealloc, edge_bulkalloc, edgefree)
CC_PTRWORLD_LEAKS_ROUTINE (edge, edge_check_leaks, x, double)

CC_PTRWORLD_LIST_ROUTINES (edgeptr, edge *, edgeptralloc, edgeptr_bulkalloc,
        edgeptrfree, edgeptr_listadd, edgeptr_listfree)
CC_PTRWORLD_LEAKS_ROUTINE (edgeptr, edgeptr_check_leaks, this, edge *)

CC_PTRWORLD_ROUTINES (node, nodealloc, node_bulkalloc, nodefree)
CC_PTRWORLD_LEAKS_ROUTINE (node, node_check_leaks, name, int)

CC_PTRWORLD_LIST_ROUTINES (nodeptr, node *, nodeptralloc, nodeptr_bulkalloc,
        nodeptrfree, nodeptr_listadd, nodeptr_listfree)
CC_PTRWORLD_LEAKS_ROUTINE (nodeptr, nodeptr_check_leaks, this, node *)


static void
    blolink_init (graph *G),
    blolink_free (graph *G),
    removedegreezero (graph *G),
    splitem (graph *G, node *n),
    splitedge (graph *G, edge *e, node *n),
    destroysplitgraph (graph *G),
    free_adj (graph *G),
    freesplitedges (graph *G),
    freesplitters (graph *G),
    markcuttree_cut (CC_GHnode *n, int v, node **names),
    initgraph (graph *G),
    freegraph (graph *G);

static int
    buildadj_from_pseudoedgelist (graph *G),
    searchtree (graph *G, CC_GHnode *n, node **names, CCtsp_lpcut_in **cuts,
        int *cutcount),
    loadcuttree_blossom (graph *G, int v, CCtsp_lpcut_in **cuts,
        int *cutcount),
    work_blossom (graph *G, nodeptr *handle, int tcount, edgeptr *teeth,
        CCtsp_lpcut_in **cuts, int *cutcount),
    add_blossom (graph *G, int hcount, int *handle, int tcount,
        toothobj *teeth, CCtsp_lpcut_in **cuts, int *cutcount),
    cuttree_tooth (edge *e, int v),
    oneend (edge *e, int v),
    buildgraph (graph *G, int ncount, int ecount, int *elist, double *x),
    grab_component (graph *G, node *n, int label, nodeptr **comp, double lbd,
        double ubd),
    grow_teeth (graph *G, nodeptr *handle, CCtsp_lpcut_in **cuts, int *cutcount),
    grow_ghteeth (graph *G, nodeptr *handle, CCtsp_lpcut_in **cuts,
        int *cutcount);


#define ONEMINUS 0.999999
#define ZEROPLUS 0.000001


int CCtsp_exactblossom (CCtsp_lpcut_in **cuts, int *cutcount, int ncount,
        int ecount, int *elist, double *x, CCrandstate *rstate)
{
    int i, k;
    node *n;
    edge *e;
    edgeptr *ep;
    CC_GHtree T;
    int gncount, gecount, markcount;
    int    *marks  = (int *) NULL;
    int    *gelist = (int *) NULL;
    node **names  = (node **) NULL;
    double *gecap  = (double *) NULL;
    graph G;
    int rval = 0;

/*
    printf ("CCtsp_exactblossom (%d, %d)...\n", ncount, ecount);
    fflush (stdout);
*/

    *cutcount = 0;
    CCcut_GHtreeinit (&T);
    initgraph (&G);
    blolink_init (&G);

    rval = buildgraph (&G, ncount, ecount, elist, x);
    if (rval) {
        fprintf (stderr, "buildgraph failed\n"); goto CLEANUP;
    }

    for (i = G.ecount, e = G.edgelist; i; i--, e++) {
        if (e->x > ONEMINUS) {
            e->splitter = G.nodelist;     /* just to kill the edge */
            e->ends[0]->mark = 1 - e->ends[0]->mark;
            e->ends[1]->mark = 1 - e->ends[1]->mark;
        } else if (e->x < ZEROPLUS) {
            e->splitter = G.nodelist;
        } else {
            e->splitter = (node *) NULL;
        }
    }

    G.pseudoedgelist = &G.pseudoedgedummy;
    G.pseudoedgelist->next = (edge *) NULL;

    G.magicnum++;
    for (n = G.pseudonodelist->next; n; n = n->next) {
        if (n->magiclabel != G.magicnum) {
            splitem (&G, n);
        }
    }

    free_adj (&G);
    rval = buildadj_from_pseudoedgelist (&G);
    if (rval) goto CLEANUP;
    
    removedegreezero (&G);

    gncount = 0;
    gecount = 0;
    markcount = 0;
    for (n = G.pseudonodelist->next; n; n = n->next) {
        for (ep = n->adj; ep; ep = ep->next) {
            ep->this->magiclabel = 0;
            gecount++;
        }
        n->num = gncount++;
        if (n->mark) markcount++;
    }
    gecount /= 2;

    if (gecount == 0) {
        /* printf ("No edges in blossom graph\n");  fflush (stdout); */
        rval = 0; goto CLEANUP; 
    }

    names  = CC_SAFE_MALLOC (gncount, node *);
    gelist = CC_SAFE_MALLOC (2*gecount, int);
    gecap  = CC_SAFE_MALLOC (gecount, double);
    if (!names || !gelist || !gecap) {
        fprintf (stderr, "out of memory in Xblossom\n");
        rval = 1; goto CLEANUP;
    }
    if (markcount) {
        marks = CC_SAFE_MALLOC (markcount, int);
        if (!marks) {
            fprintf (stderr, "out of memory in Xblossom\n");
            rval = 1; goto CLEANUP;
        }
    }

    k = 0;
    markcount = 0;
    for (i = 0, n = G.pseudonodelist->next; n; n = n->next, i++) {
        names[i] = n;
        for (ep = n->adj; ep; ep = ep->next) {
            if (ep->this->magiclabel == 0) {
                gelist[2*k]     = ep->this->ends[0]->num;
                gelist[2*k + 1] = ep->this->ends[1]->num;
                gecap[k]        = ep->this->x;
                k++;
                ep->this->magiclabel = 1;
            }
        }
        if (n->mark) marks[markcount++] = n->num;
    }

    if (markcount > 0) {
        rval = CCcut_gomory_hu (&T, gncount, gecount, gelist, gecap, 
                                markcount, marks, rstate);
        if (rval) {
            fprintf (stderr, "CCcut_gomory_hu failed\n"); goto CLEANUP;
        }
    }

    CC_IFFREE (marks, int);
    CC_IFFREE (gelist, int);
    CC_IFFREE (gecap, double);

    if (T.root) {
        rval = searchtree (&G, T.root, names, cuts, cutcount);
        if (rval) {
            fprintf (stderr, "searchtree failed\n"); goto CLEANUP;
        }
    }

CLEANUP:

    CCcut_GHtreefree (&T);
    destroysplitgraph (&G);
    freegraph (&G);
    blolink_free (&G);

    CC_IFFREE (marks, int);
    CC_IFFREE (gelist, int);
    CC_IFFREE (gecap, double);
    CC_IFFREE (names, node *);

    return rval;
}

static void blolink_init (graph *G)
{
    CCptrworld_init (&G->edge_world);
    CCptrworld_init (&G->edgeptr_world);
    CCptrworld_init (&G->node_world);
    CCptrworld_init (&G->nodeptr_world);
}

static void blolink_free (graph *G)
{
    int total, onlist;

    if (edge_check_leaks (&G->edge_world, &total, &onlist)) {
        fprintf (stderr, "WARNING: %d outstanding BLOSSOM-edges\n",
                 total - onlist);
    }
    CCptrworld_delete (&G->edge_world);

    if (edgeptr_check_leaks (&G->edgeptr_world, &total, &onlist)) {
        fprintf (stderr, "WARNING: %d outstanding BLOSSOM-edgeptrs\n",
                 total - onlist);
    }
    CCptrworld_delete (&G->edgeptr_world);

    if (node_check_leaks (&G->node_world, &total, &onlist)) {
        fprintf (stderr, "WARNING: %d outstanding BLOSSOM-nodes\n",
                 total - onlist);
    }
    CCptrworld_delete (&G->node_world);

    if (nodeptr_check_leaks (&G->nodeptr_world, &total, &onlist)) {
        fprintf (stderr, "WARNING: %d outstanding BLOSSOM-nodeptrs\n",
                 total - onlist);
    }
    CCptrworld_delete (&G->nodeptr_world);
}

static int buildadj_from_pseudoedgelist (graph *G)
{
    edge *e;
    int rval;

    for (e = G->pseudoedgelist->next; e; e = e->next) {
        rval = edgeptr_listadd (&(e->ends[0]->adj), e, &G->edgeptr_world);
        if (rval) return rval;
        rval = edgeptr_listadd (&(e->ends[1]->adj), e, &G->edgeptr_world);
        if (rval) return rval;
    }
    return 0;
}

static void removedegreezero (graph *G)
{
    node *n;

    for (n = G->pseudonodelist->next; n; n = n->next) {
        if (n->adj == (edgeptr *) NULL) {
            if (n->next != (node *) NULL) {
                n->next->prev = n->prev;
                n->prev->next = n->next;
            } else {
                n->prev->next = (node *) NULL;
            }
        }
    }
}

static void splitem (graph *G, node *n)
{
    edgeptr *ep;
    edge *e, *last;
    node *child;

    n->magiclabel = G->magicnum;
    for (ep = n->adj; ep; ep = ep->next) {
        e = ep->this;
        if (!e->splitter) {
            child = OTHEREND (e, n);
            if (child->magiclabel != G->magicnum)
                splitem (G, child);
        }
    }
    for (ep = n->adj, last = (edge *) NULL; ep; ep = ep->next) {
        e = ep->this;
        if (!e->splitter) {
            if (last)
                splitedge (G, last, n);
            last = e;
        }
    }

    if (last) {
        if (n->mark)
            splitedge (G, last, n);
        else
            splitedge (G, last, OTHEREND (last, n));
    }
}

static void splitedge (graph *G, edge *e, node *n)
{
    node *n1;
    edge *e1;

    n->mark = 1 - n->mark;

    n1 = nodealloc (&G->node_world);
    n1->magiclabel = 0;
    e->splitter = n1;
    n1->pe = e;
    n1->oddnode = n;  
    n1->mark = 1;
    n1->adj = (edgeptr *) NULL;

    n1->next = G->pseudonodelist->next;
    n1->prev = G->pseudonodelist;
    G->pseudonodelist->next->prev = n1;
    G->pseudonodelist->next = n1;


    e1 = edgealloc (&G->edge_world);
    e1->ends[0] = n;
    e1->ends[1] = n1;
    e1->x = 1.0 - e->x;
    e1->next = G->pseudoedgelist->next;
    G->pseudoedgelist->next = e1;


    e1 = edgealloc (&G->edge_world);
    e1->ends[0] = OTHEREND (e, n);
    e1->ends[1] = n1;
    e1->x = e->x;
    e1->next = G->pseudoedgelist->next;
    G->pseudoedgelist->next = e1;
}

static void destroysplitgraph (graph *G)
{
    free_adj (G);
    freesplitedges (G);
    freesplitters (G);
}

static void free_adj (graph *G)
{
    node *n;

    for (n = G->pseudonodelist->next; n; n = n->next) {
        edgeptr_listfree (&G->edgeptr_world, n->adj);
        n->adj = (edgeptr *) NULL;
    }
}

static void freesplitedges (graph *G)
{
    edge *e, *enext;

    for (e = G->pseudoedgelist->next; e; e = enext) {
        enext = e->next;
        edgefree (&G->edge_world, e);
    }
}

static void freesplitters (graph *G)
{
    node *n, *nnext, *prev;

    for (n = G->pseudonodelist->next, prev = G->pseudonodelist; n; n = nnext) {
        nnext = n->next;
        if (n->pe) {
            prev->next = nnext;
            if (nnext) nnext->prev = prev;
            nodefree (&G->node_world, n);
        } else {
            prev = n;
        }
    }
}

static int searchtree (graph *G, CC_GHnode *n, node **names,
        CCtsp_lpcut_in **cuts, int *cutcount)
{
    CC_GHnode *c;
    int rval = 0;

    if (n->ndescendants % 2 == 1  &&  n->ndescendants > 1  ) {
        if (n->cutval < 1.0 - BLOTOLERANCE) {
            G->magicnum++;
            markcuttree_cut (n, G->magicnum, names);
            rval = loadcuttree_blossom (G,  G->magicnum, cuts, cutcount);
            if (rval) {
                fprintf (stderr, "loadcuttree_blossom failed\n");
                goto CLEANUP;
            }
        }
    }
    for (c = n->child; c; c = c->sibling) {
        rval = searchtree (G, c, names, cuts, cutcount);
        if (rval) {
            goto CLEANUP;
        }
    }

CLEANUP:

    return rval;
}

static void markcuttree_cut (CC_GHnode *n, int v, node **names)
{
    CC_GHnode *c;
    int i;

    for (i = 0; i < n->listcount; i++) {
        names[n->nlist[i]]->magiclabel = v;
    }
    for (c = n->child; c; c = c->sibling) markcuttree_cut (c, v, names);
}

static int loadcuttree_blossom (graph *G, int v, CCtsp_lpcut_in **cuts,
        int *cutcount)
{
    int i, rval = 0;
    edge *e;
    node *n;
    nodeptr *handle = (nodeptr *) NULL;
    edgeptr *teeth  = (edgeptr *) NULL;
    int tcount = 0;
    int hcount = 0;

    for (i = G->ecount, e = G->edgelist; i; i--, e++) {
        if (oneend (e, v) && cuttree_tooth (e, v)) {
            rval = edgeptr_listadd (&teeth, e, &G->edgeptr_world);
            if (rval) goto CLEANUP;
            tcount++;
        }
    }
    for (i = G->ncount, n = G->nodelist; i; i--, n++) {
        if (n->magiclabel == v) {
            rval = nodeptr_listadd (&handle, n, &G->nodeptr_world);
            if (rval) goto CLEANUP;
            hcount++;
        }
    }

    if (hcount >= 3 && (tcount % 2 == 1)) {
        rval = work_blossom (G, handle, tcount, teeth, cuts,
                             cutcount);
        if (rval) {
            fprintf (stderr, "work_blossom failed\n"); goto CLEANUP;
        }
    }

CLEANUP:

    edgeptr_listfree (&G->edgeptr_world, teeth);
    nodeptr_listfree (&G->nodeptr_world, handle);

    return rval;
}

static int work_blossom (graph *G, nodeptr *handle, int tcount,
        edgeptr *teeth, CCtsp_lpcut_in **cuts, int *cutcount)
{
    edge *e;
    edgeptr *ep;
    nodeptr *np;
    toothobj *t = (toothobj *) NULL;
    toothobj *newteeth = (toothobj *) NULL;
    int i, k, newhcount, newtcount, rval = 0;
    int *del = (int *) NULL;
    int *add = (int *) NULL;
    int *hit = (int *) NULL;
    int *newhandle = (int *) NULL;

    /* Clean up intersecting teeth */

    t   = CC_SAFE_MALLOC (tcount, toothobj);
    hit = CC_SAFE_MALLOC (G->ncount, int);
    del = CC_SAFE_MALLOC (G->ncount, int);
    add = CC_SAFE_MALLOC (G->ncount, int);
    if (!t || !hit || !del || !add) {
        fprintf (stderr, "out of memory in work_blossom\n");
    }

    G->magicnum++;
    for (np = handle; np; np = np->next) {
        np->this->magiclabel = G->magicnum;
    }
    for (ep = teeth, i = 0; ep; ep = ep->next, i++) {
        e = ep->this;
        if (e->ends[0]->magiclabel == e->ends[1]->magiclabel) goto CLEANUP;
        if (e->ends[0]->magiclabel == G->magicnum) {
            t[i].in  = e->ends[0]->name;
            t[i].out = e->ends[1]->name;
        } else {
            t[i].in  = e->ends[1]->name;
            t[i].out = e->ends[0]->name;
        }
    }

    for (i = 0;  i < tcount; i++) {
        hit[t[i].in] = 0;
        del[t[i].in] = 0;
    } 
    for (i = 0;  i < tcount; i++) {
        if (hit[t[i].in]) del[t[i].in] = 1;
        else               hit[t[i].in] = 1;
    }

    for (i = 0; i < tcount; i++) {
        hit[t[i].out] = 0;
        add[t[i].out] = 0;
    }
    for (i = 0; i < tcount; i++) {
        if (hit[t[i].out]) add[t[i].out] = 1;
        else              hit[t[i].out] = 1;
    }

    for (i = 0; i < G->ncount; i++) hit[i] = 0;

    for (np = handle; np; np = np->next) {
        hit[np->this->name] = 1;
    }
    for (i = 0; i < tcount; i++) {
        if (del[t[i].in]) hit[t[i].in] = 0;
    }
    for (i = 0; i < tcount; i++) {
        if (add[t[i].out]) hit[t[i].out] = 1;
    }
  
    for (i = 0, newhcount = 0; i < G->ncount; i++) {
        if (hit[i]) newhcount++;
    }
    for (i = 0, newtcount = 0; i < tcount; i++) {
        if (hit[t[i].in] != hit[t[i].out]) newtcount++;
    }

#ifdef  ALLOW_SINGLE_TEETH
    if (newhcount >= 3 && (newtcount % 2 == 1)) {
#else
    if (newhcount >= 3 && (newtcount % 2 == 1) && newtcount >= 3) {
#endif
        newhandle = CC_SAFE_MALLOC (newhcount, int);
        newteeth  = CC_SAFE_MALLOC (newtcount, toothobj);
        if (!newhandle || !newteeth) {
            fprintf (stderr, "out of memory in work_blossom\n");
            rval = 1; goto CLEANUP;
        }
        k = 0;
        for (i = 0; i < G->ncount; i++) {
            if (hit[i]) newhandle[k++] = i;
        }
        k = 0;
        for (i = 0; i < tcount; i++) {
            if (hit[t[i].in] != hit[t[i].out]) {
                newteeth[k].in  = t[i].in;
                newteeth[k].out = t[i].out;
                k++;
            }
        }
        rval = add_blossom (G, newhcount, newhandle, newtcount, newteeth, cuts,
                            cutcount);
        CCcheck_rval (rval, "add_blossom failed");
   }

CLEANUP:

    CC_IFFREE (t, toothobj);
    CC_IFFREE (hit, int);
    CC_IFFREE (del, int);
    CC_IFFREE (add, int);
    CC_IFFREE (newhandle, int);
    CC_IFFREE (newteeth, toothobj);

    return rval;
}

static int add_blossom (graph *G, int hcount, int *handle, int tcount,
        toothobj *teeth, CCtsp_lpcut_in **cuts, int *cutcount)
{
    int itooth[2];
    CCtsp_lpcut_in *lc;
    int i, rval = 0;

    lc = CC_SAFE_MALLOC (1, CCtsp_lpcut_in);
    CCcheck_NULL (lc, "out of memory in add_blossom");
    CCtsp_init_lpcut_in (lc);

    rval = CCtsp_create_lpcliques (lc, tcount+1);
    CCcheck_rval (rval, "CCtsp_creat_lpcliques_failed");

    rval = CCtsp_array_to_lpclique (handle, hcount, &(lc->cliques[0]));
    CCcheck_rval (rval, "CCtsp_array_to_lpclique_failed");

    for (i = 0; i < tcount; i++) {
        itooth[0] = teeth[i].in;
        itooth[1] = teeth[i].out;
        rval = CCtsp_array_to_lpclique (itooth, 2, &(lc->cliques[i+1]));
        CCcheck_rval (rval, "CCtsp_array_to_lpclique_failed");
    }

    lc->rhs         = CCtsp_COMBRHS (lc);
    lc->sense       = 'G';
    lc->branch      = 0;

    rval = CCtsp_construct_skeleton (lc, G->ncount);
    if (rval) {
        fprintf (stderr, "CCtsp_construct_skeleton failed\n"); goto CLEANUP;
    }

    lc->next = *cuts;
    *cuts = lc;
    /* CCtsp_print_lpcut_in (lc); */

    (*cutcount)++;

CLEANUP:

    if (rval) {
        CCtsp_free_lpcut_in (lc);
        CC_IFFREE (lc, CCtsp_lpcut_in);
    }

    return rval;
}

static int cuttree_tooth (edge *e, int v)
{
    if (e->x > ONEMINUS) return 1;
    if (e->x < ZEROPLUS) return 0;

    if (e->splitter->magiclabel == v) {
        if (e->splitter->oddnode->magiclabel == v) return 0;
        else                                       return 1;
    } else {
        if (e->splitter->oddnode->magiclabel == v) return 1;
        else                                       return 0;
    }
}

static int oneend (edge *e, int v)
{
    if ((e->ends[0]->magiclabel == v && e->ends[1]->magiclabel != v) ||
        (e->ends[0]->magiclabel != v && e->ends[1]->magiclabel == v)) {
        return 1;
    } else {
        return 0;
    }
}

static int buildgraph (graph *G, int ncount, int ecount, int *elist, double *x)
{
    node *n;
    edge *e;
    int i, k, n1, n2;
    int rval = 0;

    initgraph (G);

    G->nodelist = CC_SAFE_MALLOC (ncount, node);
    G->edgelist = CC_SAFE_MALLOC (ecount, edge);

    if (!G->nodelist || !G->edgelist) {
        fprintf (stderr, "out of memory in buildgraph\n");
        rval = 1; goto CLEANUP;
    }
    G->ncount = ncount;
    G->ecount = ecount;

    for (i = 0; i < ncount; i++) {
        G->nodelist[i].adj = (edgeptr *) NULL;
        G->nodelist[i].magiclabel = 0;
        G->nodelist[i].name = i;
        G->nodelist[i].pe = (edge *) NULL;
        G->nodelist[i].mark = 0;
    }
    for (i = 0, e = G->edgelist, k = 0; i < ecount; i++, e++) {
        n1 = elist[k++];
        n2 = elist[k++];
        e->ends[0] = G->nodelist + n1;
        e->ends[1] = G->nodelist + n2;
        e->x = x[i];
        e->magiclabel = 0;
    }

    for (i = ecount, e = G->edgelist; i; i--, e++) {
        rval = edgeptr_listadd (&(e->ends[0]->adj), e, &G->edgeptr_world);
        if (rval) goto CLEANUP;
        rval = edgeptr_listadd (&(e->ends[1]->adj), e, &G->edgeptr_world);
        if (rval) goto CLEANUP;
    }

    G->pseudonodelist = &G->pseudonodedummy;
    G->pseudonodedummy.prev = (node *) NULL;
    G->pseudonodedummy.next = G->nodelist;
    for (i = 0, n = G->nodelist; i < G->ncount; i++, n++) {
        n->prev = n - 1;
        n->next = n + 1;
    }
    G->nodelist->prev = G->pseudonodelist;
    G->nodelist[G->ncount - 1].next = (node *) NULL;

CLEANUP:

    if (rval) freegraph (G);
    return rval;
}

static void initgraph (graph *G)
{
    if (G) {
        G->ncount = 0;
        G->nodelist = (node *) NULL;
        G->ecount = 0;
        G->edgelist = (edge *) NULL;
        G->pseudonodelist = (node *) NULL;
        G->pseudoedgelist = (edge *) NULL;
        G->magicnum = 0;
    }
}

static void freegraph (graph *G)
{
    int i;
    node *n;

    if (G->nodelist) {
        for (i = G->ncount, n = G->nodelist; i; i--, n++) {
            edgeptr_listfree (&G->edgeptr_world, n->adj);
            n->adj = (edgeptr *) NULL;
        }
        CC_FREE (G->nodelist, node);
    }
    CC_IFFREE (G->edgelist, edge);
}

int CCtsp_fastblossom (CCtsp_lpcut_in **cuts, int *cutcount, int ncount,
        int ecount, int *elist, double *x)
{
    graph G;
    int rval = 0;
    nodeptr *handle;
    int i, k;

    *cutcount = 0;
    initgraph (&G);
    blolink_init (&G);

    rval = buildgraph (&G, ncount, ecount, elist, x);
    if (rval) {
        fprintf (stderr, "buildgraph failed\n"); goto CLEANUP;
    }

    k = 0;
    for (i = 0; i < G.ncount; i++) {
        if (G.nodelist[i].mark == 0) {
            handle = (nodeptr *) NULL;
            rval = grab_component (&G, &(G.nodelist[i]), ++k, &handle,
                                   ZEROPLUS, ONEMINUS);
            if (rval) {
                fprintf (stderr, "grab_component failed\n"); goto CLEANUP;
            }
            grow_teeth (&G, handle, cuts, cutcount);
            nodeptr_listfree (&G.nodeptr_world, handle);
        }
    }

CLEANUP:

    freegraph (&G);
    blolink_free (&G);

    return rval;
}

static int grab_component (graph *G, node *n, int label, nodeptr **comp,
        double lbd, double ubd)
{
    node *v;
    edgeptr *ep;
    int rval;

    n->mark = label;
    rval = nodeptr_listadd (comp, n, &G->nodeptr_world);
    if (rval) return 1;

    for (ep = n->adj; ep; ep = ep->next) {
        if (ep->this->x > lbd && ep->this->x < ubd) {
             v = OTHEREND (ep->this, n);
             if (v->mark == 0) {
                 if (grab_component (G, v, label, comp, lbd, ubd)) {
                     return 1;
                 }
             }
        }
    }
    return 0;
}

static int grow_teeth (graph *G, nodeptr *handle, CCtsp_lpcut_in **cuts,
        int *cutcount)
{
    edgeptr *ep, *teeth = (edgeptr *) NULL;
    nodeptr *np;
    node *n;
    int tcount = 0, hcount = 0, rval = 0;

    for (np = handle; np; np = np->next) {
        hcount++;
    }
    if (hcount < 3) goto CLEANUP;
 
    for (np = handle; np; np = np->next) {
        n = np->this;
        for (ep = n->adj; ep; ep = ep->next) {
            if (ep->this->x >= ONEMINUS) {
                if (ep->this->ends[0]->mark != ep->this->ends[1]->mark) {
                    rval = edgeptr_listadd (&teeth, ep->this, &G->edgeptr_world);
                    if (rval) goto CLEANUP;
                    tcount++;
                }
            }
        }
    }

    if (tcount % 2) {
        rval = work_blossom (G, handle, tcount, teeth, cuts, cutcount);
        if (rval) {
            fprintf (stderr, "work_blossom failed\n"); goto CLEANUP;
        }
    }

CLEANUP:

    edgeptr_listfree (&G->edgeptr_world, teeth);
    return rval;
}

#define GH_EPS 0.3

int CCtsp_ghfastblossom (CCtsp_lpcut_in **cuts, int *cutcount, int ncount,
        int ecount, int *elist, double *x)
{
    graph G;
    int rval = 0;
    nodeptr *handle;
    int i, k;

    /* NOTE: Groetchel and Holland use a lowerbound of GH_EPS for the  */
    /*       edges allowed in the graph, while we use ZEROPLUS (a much */
    /*       smaller number).                                          */

    *cutcount = 0;
    initgraph (&G);
    blolink_init (&G);

    rval = buildgraph (&G, ncount, ecount, elist, x);
    if (rval) {
        fprintf (stderr, "buildgraph failed\n"); goto CLEANUP;
    }

    k = 0;
    for (i = 0; i < G.ncount; i++) {
        if (G.nodelist[i].mark == 0) {
            handle = (nodeptr *) NULL;
            rval = grab_component (&G, &(G.nodelist[i]), ++k, &handle,
                                   ZEROPLUS, 1.0 - GH_EPS);
            if (rval) {
                fprintf (stderr, "grab_component failed\n"); goto CLEANUP;
            }
            grow_ghteeth (&G, handle, cuts, cutcount);
            nodeptr_listfree (&G.nodeptr_world, handle);
        }
    }

CLEANUP:

    freegraph (&G);
    blolink_free (&G);

    return rval;
}

static int grow_ghteeth (graph *G, nodeptr *handle, CCtsp_lpcut_in **cuts,
        int *cutcount)
{
    edgeptr *ep, *teeth = (edgeptr *) NULL, *pteeth = (edgeptr *) NULL;
    edge *emax = (edge *) NULL, **tlist = (edge **) NULL;
    nodeptr *np;
    node *n;
    int i, ptcount = 0, hcount = 0, rval = 0;
    int *tperm = (int *) NULL;
    double z = 0.0, xemax = 0.0;
    double *xtlist = (double *) NULL;

    for (np = handle; np; np = np->next) {
        hcount++;
    }
    if (hcount < 3) goto CLEANUP;
 
    for (np = handle; np; np = np->next) {
        n = np->this;
        for (ep = n->adj; ep; ep = ep->next) {
            if (ep->this->x > ZEROPLUS) {
                if (ep->this->ends[0]->mark != ep->this->ends[1]->mark) {
                    if (ep->this->x >= 1.0 - GH_EPS) {
                        rval = edgeptr_listadd (&pteeth, ep->this, &G->edgeptr_world);
                        if (rval) goto CLEANUP;
                        ptcount++;
                    } else if (ep->this->x <= ZEROPLUS) {
                        if (ep->this->x > xemax) {
                            xemax = ep->this->x;
                            emax  = ep->this;
                        }
                    }
                } else {
                    z += ep->this->x;
                }
            }
        }
    }
    z *= 0.5;

    if (ptcount == 0) {
         goto CLEANUP;
    } else if (ptcount % 2 == 0) {
        if (emax) {
            rval = edgeptr_listadd (&pteeth, emax, &G->edgeptr_world);
            if (rval) goto CLEANUP;
            ptcount++;
        } 
    }

    tlist  = CC_SAFE_MALLOC (ptcount, edge *);
    xtlist = CC_SAFE_MALLOC (ptcount, double);
    tperm  = CC_SAFE_MALLOC (ptcount, int);
    if (!tlist || !xtlist || !tperm) {
        fprintf (stderr, "out of memory in grow_ghteeth\n");
        rval = 1; goto CLEANUP;
    }

    for (i = 0, ep = pteeth; ep; ep = ep->next, i++) {
        tlist[i]  = ep->this;
        xtlist[i] = -ep->this->x;
        tperm[i]  = i;
    }

    CCutil_double_perm_quicksort (tperm, xtlist, ptcount);
    if (ptcount % 2 == 0) ptcount--;

    z += tlist[tperm[0]]->x;
    rval = edgeptr_listadd (&teeth, tlist[tperm[0]], &G->edgeptr_world);
    if (rval) goto CLEANUP;
    i = 1;

    if (z > (double) hcount + (double) ((i - 1)/2) + BLOTOLERANCE) {
        rval = work_blossom (G, handle, i, teeth, cuts, cutcount);
        if (rval) {
            fprintf (stderr, "work_blossom failed\n"); goto CLEANUP;
        }
    } else {
        while (i < ptcount) {
            z += tlist[tperm[i]]->x;
            rval = edgeptr_listadd (&teeth, tlist[tperm[i]], &G->edgeptr_world);
            if (rval) goto CLEANUP;
            i++;
            z += tlist[tperm[i]]->x;
            rval = edgeptr_listadd (&teeth, tlist[tperm[i]], &G->edgeptr_world);
            if (rval) goto CLEANUP;
            i++;

            if (z > (double) hcount + (double) ((i - 1)/2) + BLOTOLERANCE) {
                rval = work_blossom (G, handle, i, teeth, cuts, cutcount);
                if (rval) {
                    fprintf (stderr, "work_blossom failed\n"); goto CLEANUP;
                }
                break;
            }
        }
    }

CLEANUP:

    edgeptr_listfree (&G->edgeptr_world, pteeth);
    edgeptr_listfree (&G->edgeptr_world, teeth);
    CC_IFFREE (tlist, edge *);
    CC_IFFREE (xtlist, double);
    CC_IFFREE (tperm, int);
    return rval;
}

