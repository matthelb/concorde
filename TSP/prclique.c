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
/*                A Padberg-Rinaldi Cliquetree Heuristic                    */
/*                                                                          */
/*                              TSP CODE                                    */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: November 4, 1997                                                  */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCtsp_pr_cliquetree (CCtsp_lpcut_in **cuts, int *cutcount,          */
/*      int ncount, int ecount, int *elist, double *x,                      */
/*      CCtsp_tighten_info *stats)                                          */
/*    BUILDS clique trees based on the conected compenents of x-1.          */
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
    int odd;
    int qhandle;
} node;

typedef struct edge {
    double x;
    int to;
    int num;
} edge;

typedef struct graph {
    node *nodelist;
    edge *adjspace;
    int   ncount;
    int   ecount;
} graph;



static int
    search_component (CCtsp_lpgraph *L, CCtsp_tighten_info *stats, double *x,
        CC_SRKexpinfo *expand, int fullncount, CCtsp_lpcut_in **newcut,
        double *viol, graph *G, int ccount, int *comp,
        CCptrworld *intptr_world),
    check_comp (CCtsp_lpgraph *L, CCtsp_tighten_info *stats,
        double *x, int *set, int setcount, intptr *onelist, graph *G,
        CCtsp_lpcut_in **newcut, double *viol, CC_SRKexpinfo *expand,
        int fullncount, CCptrworld *intptr_world),
    grab_cliquetree (intptr *T, intptr *u1, intptr *u2, int setcount, int *set,
        intptr *onelist, graph *G, CCtsp_lpcut_in **cut,
        CC_SRKexpinfo *expand, int fullncount, CCptrworld *intptr_world),
    build_blossom (int *set, int setcount, intptr *toothlist,
       CC_SRKexpinfo *expand, int fullncount, CCtsp_lpcut_in **newcut,
       double *viol),
    growT (intptr **T, graph *H, CCptrworld *intptr_world),
    growTprime (intptr **T, graph *H, CCptrworld *intptr_world),
    divide_comp (intptr *T, graph *H, intptr **u1, intptr **u2,
        CCptrworld *intptr_world),
    find_split (graph *H, intptr *T, int *marks, int *success),
    shrink_ones (int ncount, int ecount, int *elist, double *dlen,
        int *oncount, int *oecount, int **olist, double **olen,
        CC_SRKexpinfo *expand),
    buildgraph (graph *G, int ncount, int ecount, int *elist, double *x),
    grab_nonzero_x (int ecount, int *elist, double *x, int *new_ecount,
        int **new_elist, double **new_x, double tol),
    grab_nonone_edges (int ecount, int *elist, double *x, int *new_ecount,
        int **new_elist, double tol),
    find_components (int ncount, int ecount, int *elist, double *x,
        int *ncomp, int **compcnt, int **comps);

static void
    connect_search (graph *G, int n, int marker, int *dstack, int *marks),
    get_unscanned_edge (int *head, int *elist, char *scanned, graph *G,
        int *e),
    mark_scanned_edges (intptr *T, char *scanned, graph *H, int *marks),
    mark_set_and_neighbors (graph *G, int count, int *set, int *marks,
        int marker),
    mark_set (int count, int *set, int *marks, int marker),
    initgraph (graph *G),
    freegraph (graph *G),
    prclique_free_world (CCptrworld *intptr_world);


CC_PTRWORLD_LIST_ROUTINES(intptr, int, intptralloc, intptr_bulkalloc,
        intptrfree, intptr_listadd, intptr_listfree)
CC_PTRWORLD_LEAKS_ROUTINE(intptr, intptr_check_leaks, this, int)


int CCtsp_pr_cliquetree (CCtsp_lpcut_in **cuts, int *cutcount, int ncount,
        int ecount, int *elist, double *x, CCtsp_tighten_info *stats)
{
    int rval    = 0;
    int ncomp   = 0;
    int kcut    = 0;
    int newecount, oncount, oecount, i, k;
    int *compcnt  = (int *) NULL;
    int *comps    = (int *) NULL;
    int *oelist   = (int *) NULL;
    int *newelist = (int *) NULL;
    graph G;
    CCtsp_lpgraph L; 
    CC_SRKexpinfo expand;
    double viol, szeit;
    double maxviol = 0.0;
    double *ox = (double *) NULL;
    double *newx = (double *) NULL;
    CCptrworld intptr_world;

    printf ("CCtsp_pr_cliquetree ...\n"); fflush (stdout);

    /* a version of the Padberg-Rinaldi clique tree heuristic described */
    /* in Mathematical Programming 47 (1990) 219-257 (section 7.1)      */

    CCtsp_init_lpgraph_struct (&L);
    CCcut_SRK_init_expinfo (&expand);
    initgraph (&G);
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

    rval = shrink_ones (ncount, newecount, newelist, newx, &oncount,
                        &oecount, &oelist, &ox, &expand);
    if (rval) {
        fprintf (stderr, "shrink_ones failed\n"); goto CLEANUP;
    }
    CC_FREE (newelist, int);

    if (oncount == 0 || oecount == 0) goto CLEANUP;

    rval = buildgraph (&G, oncount, oecount, oelist, ox);
    if (rval) {
        fprintf (stderr, "buildgraph failed\n"); goto CLEANUP;
    }
    rval = find_components (oncount, oecount, oelist, ox, &ncomp,
                            &compcnt, &comps);
    if (rval) {
        fprintf (stderr, "find_components failed\n"); goto CLEANUP;
    }

    for (i = 0, k = 0; i < ncomp; i++) {
        if (compcnt[i] >= 7) {
            CCtsp_lpcut_in *newcut = (CCtsp_lpcut_in *) NULL;
            rval = search_component (&L, stats, newx, &expand, ncount,
                         &newcut, &viol, &G, compcnt[i], comps + k,
                         &intptr_world);
            if (rval) {
                fprintf (stderr, "search_component failed\n"); goto CLEANUP;
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
        k += compcnt[i]; 
    }

    if (cutcount) *cutcount = kcut;

    printf ("PR Cliquetrees : %d cuts  %.4f max violation  %.2f seconds\n",
                  kcut, maxviol, CCutil_zeit () - szeit);
    fflush (stdout);

CLEANUP:

    freegraph (&G);
    CCtsp_free_lpgraph (&L);
    CCcut_SRK_free_expinfo (&expand);

    CC_IFFREE (comps, int);
    CC_IFFREE (compcnt, int);

    CC_IFFREE (newelist, int);
    CC_IFFREE (newx, double);

    CC_IFFREE (oelist, int);
    CC_IFFREE (ox, double);

    prclique_free_world (&intptr_world);

    return rval;
}

static int search_component (CCtsp_lpgraph *L, CCtsp_tighten_info *stats,
         double *x, CC_SRKexpinfo *expand, int fullncount,
         CCtsp_lpcut_in **newcut, double *viol, graph *G, int ccount,
         int *comp, CCptrworld *intptr_world)
{
    int rval = 0;
    int oddcount = 0;
    int extracount = 0;
    int *marks = (int *) NULL;
    int *extraset = (int *) NULL;
    int *set;
    int i, j, k, n, setcount;
    intptr *onelist = (intptr *) NULL;
    intptr *extralist = (intptr *) NULL;
    intptr *ip;

    marks = CC_SAFE_MALLOC (G->ncount, int);
    if (marks == (int *) NULL) {
        fprintf (stderr, "out of memory in search_component\n");
        rval = 1; goto CLEANUP;
    }
    mark_set_and_neighbors (G, ccount, comp, marks, 0);
    mark_set (ccount, comp, marks, 3);

    for (i = 0; i < ccount; i++) {
        n = comp[i];
        for (j = 0; j < G->nodelist[n].degree; j++) {
            if (G->nodelist[n].adj[j].x >= 1.0 - X_FLUFF) {
                if (marks[G->nodelist[n].adj[j].to] < 3) {
                    marks[G->nodelist[n].adj[j].to]++;
                    if (marks[G->nodelist[n].adj[j].to] == 2) {
                        extracount++;
                        rval = intptr_listadd (&extralist,
                                G->nodelist[n].adj[j].to, intptr_world);
                        if (rval) {
                            fprintf (stderr, "intptr_listadd failed\n");
                            goto CLEANUP;
                        }
                    }
                }
            }
        }
    }
    for (i = 0; i < ccount; i++) {
        n = comp[i];
        for (j = 0; j < G->nodelist[n].degree; j++) {
            if (G->nodelist[n].adj[j].x >= 1.0 - X_FLUFF) {
                if (marks[G->nodelist[n].adj[j].to] < 2) {
                    marks[n]++;
                    fflush (stdout);
                    rval = intptr_listadd (&onelist, n, intptr_world);
                    if (rval) {
                        fprintf (stderr, "intptr_listadd failed\n"); goto CLEANUP;
                    }
                    rval = intptr_listadd (&onelist, G->nodelist[n].adj[j].to,
                                        intptr_world);
                    if (rval) {
                        fprintf (stderr, "intptr_listadd failed\n"); goto CLEANUP;
                    }
                }
            }
        }
    }

    for (i = 0; i < ccount; i++) {
        n = comp[i];
        if (marks[n] == 4) {
            oddcount++;
        }
    }
    CC_IFFREE (marks, int);

    if (oddcount == 0 || oddcount == 2) goto CLEANUP; 

    if (extracount) {
        extraset = CC_SAFE_MALLOC (ccount + extracount, int);
        if (extraset == (int *) NULL) {
            fprintf (stderr, "out of memory in search_component\n");
            rval = 1; goto CLEANUP;
        }
        for (k = 0; k < ccount; k++) {
            extraset[k] = comp[k];
        }
        for (ip = extralist; ip; ip = ip->next) {
            extraset[k++] = ip->this;
        }
        if (k != ccount + extracount) {
            fprintf (stderr, "error in search_component\n");
            rval = 1; goto CLEANUP;
        }
        set = extraset;
        setcount = ccount + extracount;
    } else {
        set = comp;
        setcount = ccount;
    }

    
    if (oddcount % 2) {
        rval = build_blossom (set, setcount, onelist, expand, fullncount,
                              newcut, viol);
        if (rval) {
            fprintf (stderr, "build_blossom failed\n"); goto CLEANUP;
        }
    } else if (oddcount >= 4) {
        rval = check_comp (L, stats, x, set, setcount, onelist, G, newcut,
                           viol, expand, fullncount, intptr_world);
        if (rval) {
            fprintf (stderr, "check_comp failed\n"); goto CLEANUP;
        }
    } 


CLEANUP:

    CC_IFFREE (extraset, int);
    CC_IFFREE (marks, int);
    intptr_listfree (intptr_world, onelist);
    intptr_listfree (intptr_world, extralist);
    return rval;
}

static int check_comp (CCtsp_lpgraph *L, CCtsp_tighten_info *stats,
        double *x, int *set, int setcount, intptr *onelist, graph *G,
        CCtsp_lpcut_in **newcut, double *viol, CC_SRKexpinfo *expand,
        int fullncount, CCptrworld *intptr_world)
{
    int rval = 0;
    int hcount = 0;
    int i, j, e, test, round, scanhead;
    int *marks = (int *) NULL;
    int *hlist = (int *) NULL;
    int *hmarks = (int *) NULL;
    int *invset = (int *) NULL;
    intptr *T  = (intptr *) NULL;
    intptr *u1 = (intptr *) NULL;
    intptr *u2 = (intptr *) NULL;
    intptr *ip;
    char *scanned = (char *) NULL;
    double val;
    double *hx = (double *) NULL;
    node *n;
    graph H;
    CCtsp_lpcut_in *cut  = (CCtsp_lpcut_in *) NULL;
    CCtsp_lpcut_in *tcut = (CCtsp_lpcut_in *) NULL;

    initgraph (&H);

    if (setcount == 0 || G->ncount == 0) goto CLEANUP;

    invset = CC_SAFE_MALLOC (G->ncount, int);
    marks = CC_SAFE_MALLOC (G->ncount, int);
    if (marks == (int *) NULL || invset == (int *) NULL) {
        fprintf (stderr, "out of memory in check_comp\n");
        rval = 1;  goto CLEANUP;
    }
    mark_set_and_neighbors (G, setcount, set, marks, 0);
    mark_set (setcount, set, marks, 1);

    for (i = 0; i < setcount; i++) {
        invset[set[i]] = i;
    }

    for (i = 0; i < setcount; i++) {
        n = &G->nodelist[set[i]];
        for (j = 0; j < n->degree; j++) {
            if (marks[n->adj[j].to] && n->adj[j].to < set[i]) {
                hcount++;
            }
        }
    }
    if (hcount == 0) goto CLEANUP;

    hlist = CC_SAFE_MALLOC (2 * hcount, int);
    hx    = CC_SAFE_MALLOC (hcount, double);
    if (hlist == (int *) NULL || hx == (double *) NULL) {
        fprintf (stderr, "out of memory in check_comp\n");
        rval = 1;  goto CLEANUP;
    }

    hcount = 0;
    for (i = 0; i < setcount; i++) {
        n = &G->nodelist[set[i]];
        for (j = 0; j < n->degree; j++) {
            if (marks[n->adj[j].to] && n->adj[j].to < set[i]) {
                hlist[2*hcount] = i;
                hlist[2*hcount+1] = invset[n->adj[j].to];
                hx[hcount] = n->adj[j].x;
                hcount++;
            }
        }
    }

    rval = buildgraph (&H, setcount, hcount, hlist, hx);
    if (rval) {
        fprintf (stderr, "buildgraph failed\n"); goto CLEANUP;
    }

    for (ip = onelist; ip; ip = ip->next) {
        if (marks[ip->this]) {
            H.nodelist[invset[ip->this]].odd =
                  1 - H.nodelist[invset[ip->this]].odd;
        }
    }
    CC_IFFREE (marks, int);

    scanned = CC_SAFE_MALLOC (hcount, char);
    hmarks = CC_SAFE_MALLOC (setcount, int);
    if (scanned == (char *) NULL || hmarks == (int *) NULL) {
        fprintf (stderr, "out of memory in check_comp\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < setcount; i++) {
        hmarks[i] = 0;
    }

    for (round = 0; round < 2 && cut == (CCtsp_lpcut_in *) NULL; round++) {
        for (i = 0; i < hcount; i++) {
            scanned[i] = 0;
        }
        scanhead = 0;
        while (scanhead < hcount && cut == (CCtsp_lpcut_in *) NULL) {
            get_unscanned_edge (&scanhead, hlist, scanned, &H, &e);
            if (e != -1) {
                rval = intptr_listadd (&T, hlist[2*e], intptr_world);
                if (rval) {
                    fprintf (stderr, "intptr_listadd failed\n"); goto CLEANUP;
                }
                rval = intptr_listadd (&T, hlist[2*e+1], intptr_world);
                if (rval) {
                    fprintf (stderr, "intptr_listadd failed\n"); goto CLEANUP;
                }
                if (round == 0) {
                    rval = growT (&T, &H, intptr_world);
                    if (rval) {
                        fprintf (stderr, "growT failed\n"); goto CLEANUP;
                    }
                } else {
                    rval = growTprime (&T, &H, intptr_world);
                    if (rval) {
                        fprintf (stderr, "growTprime failed\n"); goto CLEANUP;
                    }
                }
                rval = divide_comp (T, &H, &u1, &u2, intptr_world);
                if (rval) {
                    fprintf (stderr, "divide_comp failed\n"); goto CLEANUP;
                }
                if (u1 != (intptr *) NULL && u2 != (intptr *) NULL) {
                    rval = grab_cliquetree (T, u1, u2, setcount, set, onelist,
                                  G, &cut, expand, fullncount, intptr_world);
                    if (rval) {
                        fprintf (stderr, "grab_cliquetree failed\n");
                        goto CLEANUP;
                    }
                }

                mark_scanned_edges (T, scanned, &H, hmarks);
                intptr_listfree (intptr_world, T);
                T = (intptr *) NULL;
                intptr_listfree (intptr_world, u1);
                u1 = (intptr *) NULL;
                intptr_listfree (intptr_world, u2);
                u2 = (intptr *) NULL;
            }
        } 
    }

    if (cut) {
        printf ("have a cut\n"); fflush (stdout);
        rval = CCtsp_test_pure_simple_cliquetree (L->ncount, cut, &test);
        if (rval) {
            fprintf (stderr, "CCtsp_test_pure_simple_cliquetree failed\n");
            CCtsp_free_lpcut_in (cut);
            CC_IFFREE (cut, CCtsp_lpcut_in);
            goto CLEANUP;
        }
        if (test == 0) {
            fprintf (stderr, "cliquetree did not pass test\n");
            CCtsp_print_lpcut_in (cut);
            CCtsp_free_lpcut_in (cut);
            CC_IFFREE (cut, CCtsp_lpcut_in);
            rval = 1; goto CLEANUP;
        }

        tcut = CC_SAFE_MALLOC (1, CCtsp_lpcut_in);
        if (tcut == (CCtsp_lpcut_in *) NULL) {
            fprintf (stderr, "out of memory in check_comp\n");
            CCtsp_free_lpcut_in (cut);
            CC_IFFREE (cut, CCtsp_lpcut_in);
            rval = 1; goto CLEANUP;
        }
        CCtsp_init_lpcut_in (tcut);
        rval = CCtsp_tighten_lpcut_in (L, cut, x, tcut, stats, &val);
        if (rval) {
            fprintf (stderr, "CCtsp_tighten_lpcut_in failed\n");
            CCtsp_free_lpcut_in (cut);
            CC_IFFREE (cut, CCtsp_lpcut_in);
            CC_IFFREE (tcut, CCtsp_lpcut_in);
            goto CLEANUP;
        }
        CCtsp_free_lpcut_in (cut);
        CC_IFFREE (cut, CCtsp_lpcut_in);

        *viol = -CCtsp_cutprice (L, tcut, x);
        *newcut = tcut;
    }

CLEANUP:

    freegraph (&H);
    CC_IFFREE (marks, int);
    CC_IFFREE (invset, int);
    CC_IFFREE (hlist, int);
    CC_IFFREE (hx, double);
    CC_IFFREE (scanned, char);
    CC_IFFREE (hmarks, int);
    intptr_listfree (intptr_world, T);
    intptr_listfree (intptr_world, u1);
    intptr_listfree (intptr_world, u2);

    return rval;
}

static int grab_cliquetree (intptr *T, intptr *u1, intptr *u2, int setcount,
       int *set, intptr *onelist, graph *G, CCtsp_lpcut_in **cut,
       CC_SRKexpinfo *expand, int fullncount, CCptrworld *intptr_world)
{
    int rval = 0;
    int i, e1, e2, tcount;
    int toothcount = 0;
    int *marks = (int *) NULL;
    int *wset  = (int *) NULL;
    int *tset  = (int *) NULL;
    int tooth[2];
    intptr *ip;
    intptr *teeth = (intptr *) NULL;
    intptr *lineup[3];
    CCtsp_lpcut_in *tcut = (CCtsp_lpcut_in *) NULL;

    *cut = (CCtsp_lpcut_in *) NULL;

    marks = CC_SAFE_MALLOC (G->ncount, int);
    if (marks == (int *) NULL) {
        fprintf (stderr, "out of memory in grab_cliquetree\n");
        rval = 1; goto CLEANUP;
    }

    mark_set_and_neighbors (G, setcount, set, marks, 0);
    for (ip = u1; ip; ip = ip->next) {
        marks[set[ip->this]] = 2;
    }
    for (ip = u2; ip; ip = ip->next) {
        marks[set[ip->this]] = 2;
    }
    for (ip = T; ip; ip = ip->next) {
        marks[set[ip->this]] = 1;
    }

    for (ip = onelist; ip; ip = ip->next->next) {
        e1 = ip->this;
        e2 = ip->next->this;

        if ((marks[e1] == 0 && marks[e2] == 2) ||
            (marks[e1] == 2 && marks[e2] == 0)) {
            toothcount++;
            rval = intptr_listadd (&teeth, e1, intptr_world);
            if (rval) {
                fprintf (stderr, "intptr_listadd failed\n"); goto CLEANUP;
            }
            rval = intptr_listadd (&teeth, e2, intptr_world);
            if (rval) {
                fprintf (stderr, "intptr_listadd failed\n"); goto CLEANUP;
            }
        }
    } 
    if (toothcount % 2 == 1) {
        fprintf (stderr, "odd number of 1-teeth in cliquetree\n");
        rval = 1; goto CLEANUP;
    }
    if (toothcount < 4) {
        fprintf (stderr, "only %d 1-teeth in cliquetree\n", toothcount);
        rval = 1; goto CLEANUP;
    }

    /* Now build the clique tree from T, u1, u2 and teeth */

    tcut = CC_SAFE_MALLOC (1, CCtsp_lpcut_in);
    CCcheck_NULL (tcut, "out of memory in grab_cliquetree");
    CCtsp_init_lpcut_in (tcut);
 
    tcut->cliques = CC_SAFE_MALLOC (toothcount + 3, CCtsp_lpclique);
    if (tcut->cliques == (CCtsp_lpclique *) NULL) {
        fprintf (stderr, "out of memory in build_blossom\n");
        rval = 1; goto CLEANUP;
    }

    wset = CC_SAFE_MALLOC (fullncount, int);
    tset = CC_SAFE_MALLOC (setcount, int);
    if (wset == (int *) NULL || tset == (int *) NULL) {
        fprintf (stderr, "out of memory in grab_cliquetree\n");
        rval = 1; goto CLEANUP;
    }

    lineup[0] = u1;
    lineup[1] = u2;
    lineup[2] = T;
   
    for (i = 0; i < 3; i++) {
        tcount = 0;
        for (ip = lineup[i]; ip; ip = ip->next) {
            tset[tcount++] = set[ip->this];
        }
        rval = CCtsp_shrunk_set_to_lpclique (tcount, tset, wset, expand,
                                       &tcut->cliques[i]);
        if (rval) {
            fprintf (stderr, "CCtsp_shrunk_set_to_lpclique failed\n");
            goto CLEANUP;
        }
        tcut->cliquecount++;
    }

    for (i = 0, ip = teeth; i < toothcount; i++, ip = ip->next->next) {
        tooth[0] = ip->this;
        tooth[1] = ip->next->this;
        rval = CCtsp_shrunk_set_to_lpclique (2, tooth, wset, expand,
                                             &tcut->cliques[i+3]);
        if (rval) {
            fprintf (stderr, "CCtsp_shrunk_set_to_lpclique failed\n");
            goto CLEANUP;
        }
        tcut->cliquecount++;
    }
    
    tcut->rhs = 2 * (tcut->cliquecount) + (tcut->cliquecount - 3);
    tcut->sense = 'G';

    rval = CCtsp_construct_skeleton (tcut, fullncount);
    if (rval) {
        fprintf (stderr, "CCtsp_construct_skeleton failed\n");
        goto CLEANUP;
    }
    *cut = tcut;

    
CLEANUP:

    if (rval) {
        CCtsp_free_lpcut_in (tcut);
        CC_IFFREE (tcut, CCtsp_lpcut_in);
    }

    CC_IFFREE (marks, int);
    CC_IFFREE (wset, int);
    CC_IFFREE (tset, int);
    intptr_listfree (intptr_world, teeth);

    return rval;
}

static int build_blossom (int *set, int setcount, intptr *toothlist,
       CC_SRKexpinfo *expand, int fullncount, CCtsp_lpcut_in **newcut,
       double *viol)
{
    int rval = 0;
    int tcount = 0;
    int i;
    int tooth[2];
    int *wset = (int *) NULL;
    intptr *ip;
    CCtsp_lpcut_in *cut = (CCtsp_lpcut_in *) NULL;

    *newcut = (CCtsp_lpcut_in *) NULL;
    *viol = 0.0;

    cut = CC_SAFE_MALLOC (1, CCtsp_lpcut_in);
    CCcheck_NULL (cut, "out of memory in build_blossom");
    CCtsp_init_lpcut_in (cut);

    for (ip = toothlist; ip && ip->next; ip = ip->next->next) {
        tcount++;
    }
    if (tcount % 2 == 0) {
        printf ("Warning: Even number of teeth in build_blossom\n");
        fflush (stdout);
        goto CLEANUP;
    }

    cut->cliques = CC_SAFE_MALLOC (tcount + 1, CCtsp_lpclique);
    if (cut->cliques == (CCtsp_lpclique *) NULL) {
        fprintf (stderr, "out of memory in build_blossom\n");
        rval = 1; goto CLEANUP;
    }
    cut->cliquecount = 0;

    wset = CC_SAFE_MALLOC (fullncount, int);
    if (wset == (int *) NULL) {
        fprintf (stderr, "out of memory in build_blossom\n");
        rval = 1; goto CLEANUP;
    }

    rval = CCtsp_shrunk_set_to_lpclique (setcount, set, wset, expand,
                                   &cut->cliques[0]);
    if (rval) {
        fprintf (stderr, "CCtsp_shrunk_set_to_lpclique failed\n");
        goto CLEANUP;
    }
    cut->cliquecount++;

    for (i = 0, ip = toothlist; i < tcount; i++, ip = ip->next->next) {
        tooth[0] = ip->this;
        tooth[1] = ip->next->this;
        rval = CCtsp_shrunk_set_to_lpclique (2, tooth, wset, expand,
                                             &cut->cliques[i+1]);
        if (rval) {
            fprintf (stderr, "CCtsp_shrunk_set_to_lpclique failed\n");
            goto CLEANUP;
        }
        cut->cliquecount++;
    }
    
    cut->rhs = CCtsp_COMBRHS (cut);
    cut->sense = 'G';

    rval = CCtsp_construct_skeleton (cut, fullncount);
    if (rval) {
        fprintf (stderr, "CCtsp_construct_skeleton failed\n");
        goto CLEANUP;
    }
    *newcut = cut;
    *viol = 1.0;

CLEANUP:

    if (rval) {
        CCtsp_free_lpcut_in (cut);
        CC_IFFREE (cut, CCtsp_lpcut_in);
    }
    CC_IFFREE (wset, int);
    return rval;
}

static int growT (intptr **T, graph *H, CCptrworld *intptr_world)
{
    int rval = 0;
    int ncount = H->ncount;
    int i, k;
    int *marks = (int *) NULL;
    intptr *ip;
    double delta;
    double *vals = (double *) NULL;
    node *n;
    node *nodelist = H->nodelist;
    CCpriority q;

    rval =  CCutil_priority_init (&q, ncount);
    if (rval) {
        fprintf (stderr, "CCutil_priority_init failed\n"); goto CLEANUP;
    } 
    for (i = 0; i < ncount; i++) {
        nodelist[i].qhandle = -1;
    }

    vals = CC_SAFE_MALLOC (ncount, double);
    if (vals == (double *) NULL) {
        fprintf (stderr, "out of memory in growT\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < ncount; i++) {
        vals[i] = 0.0;
    }

    marks = CC_SAFE_MALLOC (ncount, int);
    if (marks == (int *) NULL) {
        fprintf (stderr, "out of memory in growT\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < ncount; i++) {
        marks[i] = 0;
    }
    for (ip = *T; ip; ip = ip->next) {
        marks[ip->this] = 1;
    }

    for (ip = *T; ip; ip = ip->next) {
        n = &nodelist[ip->this];
        for (i = 0; i < n->degree; i++) {
            k = n->adj[i].to;
            vals[k] += n->adj[i].x;
            if (marks[k] == 0 /* && vals[k] != n->adj[i].x */ ) {
                if (nodelist[k].qhandle < 0) {
                    nodelist[k].qhandle = CCutil_priority_insert (&q,
                                  (void *) (&nodelist[k]), -vals[k]); 
                } else {
                    CCutil_priority_changekey (&q, nodelist[k].qhandle,
                                  -vals[k]);
                }
            }
        }
    }
    while ((n = (node *) CCutil_priority_deletemin (&q, &delta)) !=
           (node *) NULL && delta <= -1.0) {
        rval = intptr_listadd (T, (int) (n - nodelist), intptr_world);
        if (rval) {
            fprintf (stderr, "intptr_listadd failed\n"); goto CLEANUP;
        }
        marks[n - nodelist] = 1;
        for (i = 0; i < n->degree; i++) {
            k = n->adj[i].to;
            vals[k] += n->adj[i].x;
            if (marks[k] == 0  /* && vals[k] != n->adj[i].x */) {
                if (nodelist[k].qhandle < 0) {
                    nodelist[k].qhandle = CCutil_priority_insert (&q,
                                  (void *) (&nodelist[k]), -vals[k]); 
                } else {
                    CCutil_priority_changekey (&q, nodelist[k].qhandle,
                                  -vals[k]);
                }
            }
        }
    }


CLEANUP:

    CC_IFFREE (vals, double);
    CC_IFFREE (marks, int);
    CCutil_priority_free (&q);

    return rval;
}

static int growTprime (intptr **T, graph *H, CCptrworld *intptr_world)
{
    int rval = 0;
    int i, k;
    int tcount = 0;
    int ncount = H->ncount;
    int *marks = (int *) NULL;
    intptr *ip;
    double delta;
    double sum = 0.0;
    double *vals = (double *) NULL;
    node *n;
    node *nodelist = H->nodelist;
    CCpriority q;

    rval =  CCutil_priority_init (&q, ncount);
    if (rval) {
        fprintf (stderr, "CCutil_priority_init failed\n"); goto CLEANUP;
    } 
    for (i = 0; i < ncount; i++) {
        nodelist[i].qhandle = -1;
    }

    vals = CC_SAFE_MALLOC (ncount, double);
    if (vals == (double *) NULL) {
        fprintf (stderr, "out of memory in growTprime\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < ncount; i++) {
        vals[i] = 0.0;
    }

    marks = CC_SAFE_MALLOC (ncount, int);
    if (marks == (int *) NULL) {
        fprintf (stderr, "out of memory in growTprime\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < ncount; i++) {
        marks[i] = 0;
    }
    for (ip = *T; ip; ip = ip->next) {
        marks[ip->this] = 1;
        tcount++;
    }

    for (ip = *T; ip; ip = ip->next) {
        n = &nodelist[ip->this];
        for (i = 0; i < n->degree; i++) {
            if (marks[n->adj[i].to]) {
                sum += n->adj[i].x;
            }
        }
    }
    sum  = sum / 2.0;

    for (ip = *T; ip; ip = ip->next) {
        n = &nodelist[ip->this];
        for (i = 0; i < n->degree; i++) {
            k = n->adj[i].to;
            vals[k] += n->adj[i].x;
            if (marks[k] == 0) {
                if (nodelist[k].qhandle < 0) {
                    nodelist[k].qhandle = CCutil_priority_insert (&q,
                                  (void *) (&nodelist[k]), -vals[k]); 
                } else {
                    CCutil_priority_changekey (&q, nodelist[k].qhandle,
                                  -vals[k]);
                }
            }
        }
    }
    while ((n = (node *) CCutil_priority_deletemin (&q, &delta)) !=
           (node *) NULL && sum - delta >= (double) tcount - 1.5) {
        rval = intptr_listadd (T, (int) (n - nodelist), intptr_world);
        if (rval) {
            fprintf (stderr, "intptr_listadd failed\n"); goto CLEANUP;
        }
        marks[n - nodelist] = 1;
        sum += vals[n - nodelist];
        tcount++;
        for (i = 0; i < n->degree; i++) {
            k = n->adj[i].to;
            vals[k] += n->adj[i].x;
            if (marks[k] == 0) {
                if (nodelist[k].qhandle < 0) {
                    nodelist[k].qhandle = CCutil_priority_insert (&q,
                                  (void *) (&nodelist[k]), -vals[k]); 
                } else {
                    CCutil_priority_changekey (&q, nodelist[k].qhandle,
                                  -vals[k]);
                }
            }
        }
    }


CLEANUP:

    CC_IFFREE (vals, double);
    CC_IFFREE (marks, int);
    CCutil_priority_free (&q);

    return rval;
}

static int divide_comp (intptr *T, graph *H, intptr **u1, intptr **u2,
        CCptrworld *intptr_world)
{
    int rval = 0;
    int tcount = 0;
    int i, test;
    int *marks = (int *) NULL;
    intptr *ip;
    node *n;

    *u1 = (intptr *) NULL;
    *u2 = (intptr *) NULL;

    for (ip = T; ip; ip = ip->next) {
        tcount++;
    }

    if (tcount < 3) goto CLEANUP;

    marks = CC_SAFE_MALLOC (H->ncount, int);
    if (marks == (int *) NULL) {
        fprintf (stderr, "out of memory in divide_comp\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < H->ncount; i++) {
        marks[i] = 0;
    }

    /* first try {v \in T: delta(v) \in T} */

    for (ip = T; ip; ip = ip->next) {
        marks[ip->this] = 1;
    }
    for (ip = T; ip; ip = ip->next) {
        for (i = 0; i < H->nodelist[ip->this].degree; i++) {
            if (!marks[H->nodelist[ip->this].adj[i].to]) {
                break;
            }
        }
        if (i == H->nodelist[ip->this].degree) {
            marks[ip->this] = 2;
        }
    }
    for (ip = T; ip; ip = ip->next) {
        marks[ip->this]--;
    }
    rval = find_split (H, T, marks, &test);
    if (rval) {
        fprintf (stderr, "find_split failed\n"); goto CLEANUP;
    }

    if (!test) {
        /* now try all of T */

        for (i = 0; i < H->ncount; i++) {
            marks[i] = 0;
        }
        for (ip = T; ip; ip = ip->next) {
            marks[ip->this] = 1;
        }
        rval = find_split (H, T, marks, &test);
        if (rval) {
            fprintf (stderr, "find_split failed\n"); goto CLEANUP;
        }
    }

    if (test) {
        double bestx, xadj;
        int u1best, u2best;

        for (i = 0; i <  H->ncount; i++) {
            if (marks[i] == 2) {
                rval = intptr_listadd (u1, i, intptr_world);
                if (rval) {
                    fprintf (stderr, "intptr_listadd failed\n"); goto CLEANUP;
                }
            } else if (marks[i] == 3) {
                rval = intptr_listadd (u2, i, intptr_world);
                if (rval) {
                    fprintf (stderr, "intptr_listadd failed\n"); goto CLEANUP;
                }
            }
        }

        u1best = -1;
        bestx = -1.0;
        for (ip = T; ip; ip = ip->next) {
            n = &(H->nodelist[ip->this]);
            xadj = 0.0;
            for (i = 0; i < n->degree; i++) {
                if (marks[n->adj[i].to] == 2) {
                    xadj += n->adj[i].x;
                }
            }
            if (xadj > bestx) {
                bestx = xadj;
                u1best = ip->this;
            }
        }
        rval = intptr_listadd (u1, u1best, intptr_world);
        if (rval) {
            fprintf (stderr, "intptr_listadd failed\n"); goto CLEANUP;
        }

        u2best = -1;
        bestx = -1.0;
        for (ip = T; ip; ip = ip->next) {
            if (ip->this != u1best) {
                n = &H->nodelist[ip->this];
                xadj = 0.0;
                for (i = 0; i < n->degree; i++) {
                    if (marks[n->adj[i].to] == 3) {
                        xadj += n->adj[i].x;
                    }
                }
                if (xadj > bestx) {
                    bestx = xadj;
                    u2best = ip->this;
                }
            }
        }
        rval = intptr_listadd (u2, u2best, intptr_world);
        if (rval) {
            fprintf (stderr, "intptr_listadd failed\n"); goto CLEANUP;
        }
    }

CLEANUP:

    if (rval) {
        intptr_listfree (intptr_world, *u1);
        intptr_listfree (intptr_world, *u2);
        *u1 = (intptr *) NULL;
        *u2 = (intptr *) NULL;
    }

    CC_IFFREE (marks, int);
    return rval;
}

static int find_split (graph *H, intptr *T, int *marks, int *success)
{
    int rval = 0;
    int counter = 0;
    int i;
    int *dstack = (int *) NULL;
    intptr *ip;

    /* label the components in the unmarked graph as 2 and 3 */

    *success = 0;

    dstack = CC_SAFE_MALLOC (H->ncount, int);
    if (dstack == (int *) NULL) {
        fprintf (stderr, "out of memory in find_split\n");
        rval = 1; goto CLEANUP;
    } 

    for (i = 0; i < H->ncount && counter < 2; i++) {
        if (marks[i] == 0) {
            counter++;
            connect_search (H, i, counter + 1, dstack, marks);
        }
    }

    if (counter == 2 && i < H->ncount) {
        int odd_u1 = 0;
        int odd_u2 = 0;

        for (ip = T; ip; ip = ip->next) {
            marks[ip->this] = 1;
        }
        
        for (i = 0; i < H->ncount; i++) {
            if (H->nodelist[i].odd) {
                if (marks[i] == 2) {
                    odd_u1++;
                } else if (marks[i] == 3) {
                    odd_u2++;
                }
            }
        }
        if ((odd_u1 >= 2) && (odd_u2 >= 2) &&
            (odd_u1 % 2 == 0) && (odd_u2 % 2 == 0)) {
            *success = 1;
        }
    } 

CLEANUP:

    CC_IFFREE (dstack, int);
    return rval;
}

static void connect_search (graph *G, int n, int marker, int *dstack,
        int *marks)
{
    int i, k, head = 0;
    node *nodelist = G->nodelist;

    marks[n] = marker;
    dstack[head++] = n;

    while (head > 0) {
        n = dstack[--head];
        for (i = 0; i < nodelist[n].degree; i++) {
            k = nodelist[n].adj[i].to;
            if (!marks[k]) {
                marks[k] = marker;
                dstack[head++] =  k;
            }
        }
    }
}

static void  get_unscanned_edge (int *head, int *elist, char *scanned, graph *G,
        int *e)
{
    int i;

    for (i = *head; i < G->ecount; i++) {
        if (scanned[i] == 0 && G->nodelist[elist[2*i]].odd == 0 &&
                               G->nodelist[elist[2*i+1]].odd == 0) {
            break;
        }
    }
    if (i < G->ecount) {
        *head = i+1;
        *e = i;
    } else {
        *head = G->ecount;
        *e = -1;
    }
}

static void mark_scanned_edges (intptr *T, char *scanned, graph *H, int *marks)
{
    int i;
    intptr *ip;
    node *n;

    for (ip = T; ip; ip = ip->next) {
        n = &H->nodelist[ip->this];
        for (i = 0; i < n->degree; i++) {
            if (marks[n->adj[i].to] == 1) {
                scanned[n->adj[i].num] = 1;
            }
        }
        marks[ip->this] = 1;
    }
    for (ip = T; ip; ip = ip->next) {
        marks[ip->this] = 0;
    }
}

static void mark_set_and_neighbors (graph *G, int count, int *set, 
       int *marks, int marker)
{
    int i, j, n;

    for (i = 0; i < count; i++) {
        n = set[i];
        marks[n] = marker;
        for (j = 0; j < G->nodelist[n].degree; j++) {
            marks[G->nodelist[n].adj[j].to] = marker; 
        }
    }
}

static void mark_set (int count, int *set, int *marks, int marker)
{
    int i;

    for (i = 0; i < count; i++) {
        marks[set[i]] = marker;
    }
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
        n[i].odd     = 0;
        n[i].qhandle = -1;
        n[i].degree  = 0;
        n[i].adj     = (edge *) NULL;
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
        n[j].adj[n[j].degree].to  = elist[2*i+1];
        n[j].adj[n[j].degree].x   = x[i];
        n[j].adj[n[j].degree].num = i;
        n[j].degree++;
        j = elist[2*i+1];
        n[j].adj[n[j].degree].to  = elist[2*i];
        n[j].adj[n[j].degree].x   = x[i];
        n[j].adj[n[j].degree].num = i;
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

static void prclique_free_world (CCptrworld *intptr_world)
{
    int total, onlist;

    if (intptr_check_leaks (intptr_world, &total, &onlist)) {
        fprintf (stderr, "WARNING: %d outstanding intptrs\n",
                 total - onlist);
    }
    CCptrworld_delete (intptr_world);
}

static int shrink_ones (int ncount, int ecount, int *elist, double *dlen,
        int *oncount, int *oecount, int **olist, double **olen,
        CC_SRKexpinfo *expand)
{
    int rval = 0;
    int k;
    CC_SRKgraph S;

    CCcut_SRK_init_graph (&S);

    rval = CCcut_SRK_buildgraph (&S, ncount, ecount, elist, dlen);
    if (rval) {
        fprintf (stderr, "buildgraph failed in shrink_ones\n");
        goto CLEANUP;
    }
    CCcut_SRK_increment_marker (&S);

    rval = CCcut_SRK_defluff (&S);
    if (rval) {
        fprintf (stderr, "CCcut_SRK_defluff failed in shrink_ones\n");
        goto CLEANUP;
    }

    CCcut_SRK_identify_paths_to_edges (&S, &k, 0);
    rval = CCcut_SRK_grab_edges (&S, oncount, oecount, olist, olen, expand);
    if (rval) {
        fprintf (stderr, "grab edges failed in shrink_ones\n");
        goto CLEANUP;
    }

CLEANUP:

    CCcut_SRK_free_graph (&S);
    return rval;
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
    if (*new_elist == (int *) NULL || *new_x == (double *) NULL) {
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

static int grab_nonone_edges (int ecount, int *elist, double *x,
        int *new_ecount, int **new_elist, double tol)
{
    int i;
    int count;
    double cutoff = 1.0 - tol;

    *new_ecount = 0;
    *new_elist = (int *) NULL;

    for (i = 0, count = 0; i < ecount; i++) {
        if (x[i] < cutoff) {
            count++;
        }
    }
    if (count == 0) return 0;

    *new_elist = CC_SAFE_MALLOC (2*count, int);
    if (*new_elist == (int *) NULL) {
        fprintf (stderr, "out of memory in grab_nonzero_x\n");
        return 1;
    }

    for (i = 0, count = 0; i < ecount; i++) {
        if (x[i] < cutoff) {
            (*new_elist)[2*count] = elist[2*i];
            (*new_elist)[2*count+1] = elist[2*i+1];
            count++;
        }
    }
    *new_ecount = count;

    return 0;
}

static int find_components (int ncount, int ecount, int *elist, double *x,
        int *ncomp, int **compcnt, int **comps)
{
    int rval = 0;
    int newecount;
    int *newelist = (int *) NULL;

    *compcnt = (int *) NULL;
    *comps   = (int *) NULL;

    rval = grab_nonone_edges (ecount, elist, x, &newecount, &newelist, X_FLUFF);
    if (rval) {
        fprintf (stderr, "grab_nonone_edges failed\n");
        goto CLEANUP;
    }

    rval = CCcut_connect_components (ncount, newecount, newelist,
               (double *) NULL, ncomp, compcnt, comps);
    if (rval) {
        fprintf (stderr, "CCcut_connect_components failed\n");
        goto CLEANUP;
    }

CLEANUP:

    CC_IFFREE (newelist, int);
    return rval;
}
