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
/*                  Find Combs Based on Greedy Growing                      */
/*                                                                          */
/*                              TSP CODE                                    */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: October 31, 1997                                                  */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCtsp_edge_comb_grower (CCtsp_lpcut_in **cuts, int *cutcount,       */
/*      int ncount, int ecount, int *elist, double *x,                      */
/*      CCtsp_tighten_info *stats)                                          */
/*    BUILDS combs with 1 inverted tooth using greedy cut.                  */
/*     -cuts (new cutting plans will be added to front of this list)        */
/*     -cutcount will return the number of new cuts found (can be NULL)     */
/*     -ncount is the number of nodes                                       */
/*     -ecount is the number of edges                                       */
/*     -elist is the edge list in node node format                          */
/*     -x is an lp solution vector                                          */
/*     -stats is used to monitor tighten                                    */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "cut.h"
#include "tsp.h"
#include "combs.h"

#define DOUBLE_EDGE_CUTOFF (10)
#define CUT_BATCH          (1000)  
#define GROW_EPS           (0.001) 
#define X_FLUFF            (1e-10)

typedef struct intptr {
    int this;
    struct intptr *next;
} intptr;

typedef struct node {
    struct edge *adj;
    int degree;
} node;

typedef struct edge {
    double x;
    int to;
} edge;

typedef struct graph {
    node *nodelist;
    edge *adjspace;
    int   ncount;
    int   ecount;
} graph;



static int
    grow_the_comb (CCtsp_lpgraph *L, CCtsp_tighten_info *stats,
        double *x, CC_GCgraph *G, int marker, CC_SRKexpinfo *expand,
        int *set, int *wset, int n1, int n2, int m1, int m2,
        CCtsp_lpcut_in **newcut, double *viol, int ncount, int orig_ncount),
    shrink_ones (CC_SRKgraph *G, int ncount, int ecount, int *elist,
        double *dlen, int *oncount, int *oecount, int **olist, double **olen,
        CC_SRKexpinfo *expand),
    buildgraph (graph *G, int ncount, int ecount, int *elist, double *x),
    grab_nonzero_x (int ecount, int *elist, double *x, int *new_ecount,
        int **new_elist, double **new_x, double tol);

static void
    initgraph (graph *G),
    freegraph (graph *G),
    growcomb_free_world (CCptrworld *intptr_world);


CC_PTRWORLD_LIST_ROUTINES(intptr, int, intptralloc, intptr_bulkalloc,
        intptrfree, intptr_listadd, intptr_listfree)
CC_PTRWORLD_LEAKS_ROUTINE (intptr, intptr_check_leaks, this, int)

int CCtsp_edge_comb_grower (CCtsp_lpcut_in **cuts, int *cutcount, int ncount,
        int ecount, int *elist, double *x, CCtsp_tighten_info *stats)
{
    int rval    = 0;
    int nblocks = 0;
    int kcut    = 0;
    int marker  = 0;
    int newecount, oncount, oecount, i, j, k, n, listcount;
    int ncutnodes = 0;
    int *cutnodes = (int *) NULL;
    int *set      = (int *) NULL;
    int *wset     = (int *) NULL;
    int *marks    = (int *) NULL;
    int *blockcnt = (int *) NULL;
    int *oelist   = (int *) NULL;
    int *newelist = (int *) NULL;
    int **blocks  = (int **) NULL;
    intptr *ones  = (intptr *) NULL;
    CCtsp_lpgraph L; 
    CC_GCgraph G;
    CC_SRKgraph S;
    CC_SRKexpinfo expand;
    graph H;
    double viol, szeit;
    double maxviol = 0.0;
    double *ox = (double *) NULL;
    double *newx = (double *) NULL;
    CCptrworld intptr_world;

    /* try to grow a 3-tooth combs, starting with a pair of 1's  */
    /* this is an idea Naddef & Clochard, described in an unpublished report */

    CCtsp_init_lpgraph_struct (&L);
    CCcut_SRK_init_graph (&S);
    CCcut_SRK_init_expinfo (&expand);
    CCcombs_GC_init_graph (&G);
    initgraph (&H);
    CCptrworld_init (&intptr_world);

    if (cutcount) *cutcount = 0;

    szeit = CCutil_zeit ();

    rval = grab_nonzero_x (ecount, elist, x, &newecount, &newelist, &newx,
                           X_FLUFF);
    if (rval) {
        fprintf (stderr, "grab_nonzero_x failed\n"); goto CLEANUP;
    }

    CCtsp_init_lpgraph_struct (&L);

    rval = CCtsp_build_lpgraph (&L, ncount, newecount, newelist,
                                (int *) NULL);
    if (rval) {
        fprintf (stderr, "CCtsp_build_lpgraph failed\n"); goto CLEANUP;
    }
    rval = CCtsp_build_lpadj (&L, 0, newecount);
    if (rval) {
        fprintf (stderr, "CCtsp_build_lpadj failed\n"); goto CLEANUP;
    }

    rval = shrink_ones (&S, ncount, newecount, newelist, newx,
                        &oncount, &oecount, &oelist, &ox, &expand);
    if (rval) {
        fprintf (stderr, "shrink_ones failed\n"); goto CLEANUP;
    }
    CCcut_SRK_free_graph (&S);
    CC_FREE (newelist, int);

    if (oncount == 0 || oecount == 0) goto CLEANUP;

    rval = CCcombs_find_blocks (oncount, oecount, oelist, ox,
             &nblocks, &blockcnt, &blocks, &ncutnodes, &cutnodes);
    if (rval) {
        fprintf (stderr, "CCcombs_find_blocks failed\n"); goto CLEANUP;
    }

    rval = CCcombs_GC_build_graph (&G, oncount, oecount, oelist, ox);
    if (rval) {
        fprintf (stderr, "CCcombs_GC_build_graph failed\n"); goto CLEANUP;
    }

    set   = CC_SAFE_MALLOC (oncount, int);
    wset  = CC_SAFE_MALLOC (ncount, int);
    marks = CC_SAFE_MALLOC (oncount, int); 
    if (set == (int *) NULL || wset == (int *) NULL || marks == (int *) NULL) {
        fprintf (stderr, "out of memory in CCcombs_block_combs\n");
        rval = 1; goto CLEANUP;
    }

    rval = buildgraph (&H, oncount, oecount, oelist, ox);
    if (rval) {
        fprintf (stderr, "buildgraph failed\n"); goto CLEANUP;
    }

    for (i = 0; i < oncount; i++) {
        marks[i] = 0;
    }

    for (i = 0; i < nblocks; i++) {
        for (j = 0; j < blockcnt[i]; j++) {
            marks[blocks[i][j]] = 1; 
        }
        ones = (intptr *) NULL;
        listcount = 0;
        for (j = 0; j < blockcnt[i]; j++) {
            intptr *ip, *jp;
            int n1, n2, m1, m2;
            n = blocks[i][j];
            for (k = 0; k < H.nodelist[n].degree; k++) {
                if (marks[H.nodelist[n].adj[k].to] == 0 &&
                    H.nodelist[n].adj[k].x >= 1.0 - GROW_EPS) {
                    rval = intptr_listadd (&ones, H.nodelist[n].adj[k].to, &intptr_world);
                    if (rval) goto CLEANUP;
                    rval = intptr_listadd (&ones, n, &intptr_world);
                    if (rval) goto CLEANUP;
                    listcount++;
                }
            }
            if (listcount < DOUBLE_EDGE_CUTOFF) {
                for (ip = ones; ip; ip = ip->next->next) {
                    n1 = ip->this;
                    n2 = ip->next->this;
                    for (jp = ip->next->next; jp; jp = jp->next->next) {
                        CCtsp_lpcut_in *newcut = (CCtsp_lpcut_in *) NULL;
                        m1 = jp->this;
                        m2 = jp->next->this;
                        if (m1 == n1 || m1 == n2 || m2 == n1 || m2 == n2) {
                            continue;
                        }
                        marker++;
                        rval = grow_the_comb (&L, stats, newx, &G, marker,
                                 &expand, set, wset, n1, n2, m1, m2, &newcut,
                                 &viol, oncount, ncount);
                        if (rval) {
                            fprintf (stderr, "grow_the_comb failed\n");
                        }
                        if (newcut) {
                            if (viol > CCtsp_MIN_VIOL) {
                                newcut->next = *cuts;
                                *cuts = newcut;
                                kcut++;
                            } else {
                                CCtsp_free_lpcut_in (newcut);
                                CC_IFFREE (newcut, CCtsp_lpcut_in);
                            }
                            if (viol > maxviol) maxviol = viol;
                        }
                    }
                }
            } else {
                for (ip = ones; ip && ip->next->next; ip = ip->next->next) {
                    CCtsp_lpcut_in *newcut = (CCtsp_lpcut_in *) NULL;
                    n1 = ip->this;
                    n2 = ip->next->this;
                    m1 = ip->next->next->this;
                    m2 = ip->next->next->next->this;
                    if (m1 == n1 || m1 == n2 || m2 == n1 || m2 == n2) {
                        continue;
                    }
                    marker++;
                    rval = grow_the_comb (&L, stats, newx, &G, marker,
                             &expand, set, wset, n1, n2, m1, m2, &newcut,
                             &viol, oncount, ncount);
                    if (rval) {
                        fprintf (stderr, "grow_the_comb failed\n");
                    }
                    if (newcut) {
                        if (viol > CCtsp_MIN_VIOL) {
                            newcut->next = *cuts;
                            *cuts = newcut;
                            kcut++;
                        } else {
                            CCtsp_free_lpcut_in (newcut);
                            CC_IFFREE (newcut, CCtsp_lpcut_in);
                        }
                        if (viol > maxviol) maxviol = viol;
                    }
                }
            }
        }
        for (j = 0; j < blockcnt[i]; j++) {
            marks[blocks[i][j]] = 0; 
        }
        intptr_listfree (&intptr_world, ones);
    }

    if (cutcount) *cutcount = kcut;

CLEANUP:

    CCtsp_free_lpgraph (&L);
    CCcut_SRK_free_graph (&S);
    CCcut_SRK_free_expinfo (&expand);
    CCcombs_GC_free_graph (&G);
    freegraph (&H);

    CC_IFFREE (cutnodes, int);
    for (i = 0; i < nblocks; i++) {
        CC_IFFREE (blocks[i], int);
    }
    CC_IFFREE (blocks, int *);
    CC_IFFREE (blockcnt, int);

    CC_IFFREE (newelist, int);
    CC_IFFREE (newx, double);

    CC_IFFREE (oelist, int);
    CC_IFFREE (ox, double);

    CC_IFFREE (set, int);
    CC_IFFREE (wset, int);
    CC_IFFREE (marks, int);

    growcomb_free_world (&intptr_world);

    return rval;
}

static int grow_the_comb (CCtsp_lpgraph *L, CCtsp_tighten_info *stats,
        double *x, CC_GCgraph *G, int marker, CC_SRKexpinfo *expand,
        int *set, int *wset, int n1, int n2, int m1, int m2,
        CCtsp_lpcut_in **newcut, double *viol, int ncount, int orig_ncount)
{
    int rval = 0;
    int j, setsize;
    double val;
    CCtsp_lpcut_in cut;
    CCtsp_lpcut_in *tcut = (CCtsp_lpcut_in *) NULL;

    *newcut = (CCtsp_lpcut_in *) NULL;
    *viol = 0.0;
    CCtsp_init_lpcut_in (&cut);

    if (ncount < 6) goto CLEANUP;

    G->nodelist[n1].mark = marker;
    G->nodelist[n2].mark = marker;
    G->nodelist[m1].mark = marker;
    G->nodelist[m2].mark = marker;

    set[0] = n1;
    set[1] = m1;
    setsize = 2;

    rval = CCcombs_greedy_cut (G, &setsize, set, marker, 1, 0, 1,
                               (int *) NULL, &val);
    if (rval) {
        fprintf (stderr, "CCcombs_greedy_cut failed\n"); goto CLEANUP;
    }

    if (setsize > ncount/4 || setsize < 3) {
        goto CLEANUP;
    }

    cut.cliques = CC_SAFE_MALLOC (4, CCtsp_lpclique);
    if (!cut.cliques) {
        fprintf (stderr, "out of memory in build_star\n");
        rval = 1; goto CLEANUP;
    }

    rval = CCtsp_shrunk_set_to_lpclique (setsize, set, wset, expand,
                                   &cut.cliques[0]);
    if (rval) {
        fprintf (stderr, "CCtsp_shrunk_set_to_lpclique failed\n");
        CC_FREE (cut.cliques, CCtsp_lpclique);
        goto CLEANUP;
    }

    setsize = 4;
    set[0] = n1;
    set[1] = n2;
    set[2] = m1;
    set[3] = m2;

    rval = CCtsp_shrunk_set_to_lpclique (setsize, set, wset, expand,
                                   &cut.cliques[1]);
    if (rval) {
        fprintf (stderr, "CCtsp_shrunk_set_to_lpclique failed\n");
        for (j = 0; j < 1; j++) {
            CCtsp_free_lpclique (&cut.cliques[j]);
        }
        CC_FREE (cut.cliques, CCtsp_lpclique);
        goto CLEANUP;
    }

    setsize = 2;
    set[0] = n1;
    set[1] = n2;

    rval = CCtsp_shrunk_set_to_lpclique (setsize, set, wset, expand,
                                   &cut.cliques[2]);
    if (rval) {
        fprintf (stderr, "CCtsp_shrunk_set_to_lpclique failed\n");
        for (j = 0; j < 2; j++) {
            CCtsp_free_lpclique (&cut.cliques[j]);
        }
        CC_FREE (cut.cliques, CCtsp_lpclique);
        goto CLEANUP;
    }
    

    setsize = 2;
    set[0] = m1;
    set[1] = m2;

    rval = CCtsp_shrunk_set_to_lpclique (setsize, set, wset, expand,
                                   &cut.cliques[3]);
    if (rval) {
        fprintf (stderr, "CCtsp_shrunk_set_to_lpclique failed\n");
        for (j = 0; j < 3; j++) {
            CCtsp_free_lpclique (&cut.cliques[j]);
        }
        CC_FREE (cut.cliques, CCtsp_lpclique);
        goto CLEANUP;
    }


    cut.cliquecount = 4;
    cut.rhs = 10;
    cut.sense = 'G';

    rval = CCtsp_construct_skeleton (&cut, orig_ncount);
    if (rval) {
        fprintf (stderr, "CCtsp_construct_skeleton failed\n");
        CCtsp_free_lpcut_in (&cut);
        goto CLEANUP;
    }

    tcut = CC_SAFE_MALLOC (1, CCtsp_lpcut_in);
    if (!tcut) {
        fprintf (stderr, "out of memory in grow_the_comb\n"); 
        CCtsp_free_lpcut_in (&cut);
        rval = 1; goto CLEANUP;
    }
    CCtsp_init_lpcut_in (tcut);

    rval = CCtsp_tighten_lpcut_in (L, &cut, x, tcut, stats, &val);  
    if (rval) {
        fprintf (stderr, "CCtsp_tighten_lpcut_in failed\n");
        CCtsp_free_lpcut_in (&cut);
        CC_IFFREE (tcut, CCtsp_lpcut_in);
        goto CLEANUP;
    }
    CCtsp_free_lpcut_in (&cut);

    *viol = -CCtsp_cutprice (L, tcut, x);
    *newcut = tcut;

CLEANUP:

    return rval;
}

static int buildgraph (graph *G, int ncount, int ecount, int *elist, double *x)
{
    int rval = 0;
    node *n;
    edge *e;
    int i, j;

    if (ncount) {
        G->nodelist = CC_SAFE_MALLOC (ncount, node);
        if (G->nodelist == (node *) NULL) {
            fprintf (stderr, "out of memory in buildgraph\n");
            rval = 1; goto CLEANUP;
        }
    }
    if (ecount) {
        G->adjspace = CC_SAFE_MALLOC (2 * ecount, edge);
        if (G->adjspace == (edge *) NULL) {
            fprintf (stderr, "out of memory in buildgraph\n");
            CC_IFFREE (G->nodelist, node);
            rval = 1; goto CLEANUP;
        }
    }

    n = G->nodelist;
    for (i = 0; i < ncount; i++) {
        n[i].degree = 0;
        n[i].adj    = (edge *) NULL;
    }

    for (i = 0; i < ecount; i++) {
        n[elist[2*i]].degree++;
        n[elist[2*i+1]].degree++;
    }

    e = G->adjspace;
    for (i = 0; i < ncount; i++) {
        n[i].adj = e;
        e += n[i].degree;
        n[i].degree = 0;
    }
    for (i = 0; i < ecount; i++) {
        j = elist[2*i];
        n[j].adj[n[j].degree].to =  elist[2*i+1];
        n[j].adj[n[j].degree].x  = x[i];
        n[j].degree++;
        j = elist[2*i+1];
        n[j].adj[n[j].degree].to =  elist[2*i];
        n[j].adj[n[j].degree].x  = x[i];
        n[j].degree++;
    }

    G->ncount = ncount;
    G->ecount = ecount;

CLEANUP:

    return rval;
}

static void initgraph (graph *G)
{
    if (G) {
        G->nodelist = (node *) NULL;
        G->adjspace = (edge *) NULL;
        G->ncount   = 0;
        G->ecount   = 0;
    }
}

static void freegraph (graph *G)
{
    if (G) {
        CC_IFFREE (G->nodelist, node);
        CC_IFFREE (G->adjspace, edge);
    }
}

static void growcomb_free_world (CCptrworld *intptr_world)
{
    int total, onlist;

    if (intptr_check_leaks (intptr_world, &total, &onlist)) {
        fprintf (stderr, "WARNING: %d outstanding intptrs\n",
                 total - onlist);
    }
    CCptrworld_delete (intptr_world);
}

static int shrink_ones (CC_SRKgraph *G, int ncount, int ecount, int *elist,
        double *dlen, int *oncount, int *oecount, int **olist, double **olen,
        CC_SRKexpinfo *expand)
{
    int k;

    if (CCcut_SRK_buildgraph (G, ncount, ecount, elist, dlen)) {
        fprintf (stderr, "buildgraph failed in shrink_ones\n");
        return 1;
    }
    CCcut_SRK_increment_marker (G);

    if (CCcut_SRK_defluff (G)) {
        fprintf (stderr, "CCcut_SRK_defluff failed in shrink_ones\n");
        return 1;
    }

    CCcut_SRK_identify_paths_to_edges (G, &k, 0);
    if (CCcut_SRK_grab_edges (G, oncount, oecount, olist, olen, expand)) {
        fprintf (stderr, "grab edges failed in shrink_ones\n");
        return 1;
    }

    return 0;
}

static int grab_nonzero_x (int ecount, int *elist, double *x, int *new_ecount,
        int **new_elist, double **new_x, double tol)
{
    int i;
    int count;

    *new_ecount = 0;
    *new_elist = (int *) NULL;
    *new_x = (double *) NULL;

    for (i = 0, count = 0; i < ecount; i++) {
        if (x[i] > tol) {
            count++;
        }
    }

    *new_elist = CC_SAFE_MALLOC (2*count, int);
    *new_x = CC_SAFE_MALLOC (count, double);
    if (!(*new_elist) || !(*new_x)) {
        fprintf (stderr, "out of memory in grab_nonzero_x\n");
        CC_IFFREE (*new_elist, int);
        CC_IFFREE (*new_x, double);
        return 1;
    }

    for (i = 0, count = 0; i < ecount; i++) {
        if (x[i] > tol) {
            (*new_elist)[2*count] = elist[2*i];
            (*new_elist)[2*count+1] = elist[2*i+1];
            (*new_x)[count] = x[i];
            count++;
        }
    }
    *new_ecount = count;

    return 0;
}

