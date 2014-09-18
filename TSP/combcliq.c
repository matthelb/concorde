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
/*               BUILD SIMPLE CLIQUETREE FROM A COMB                        */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: July 3, 1997                                                      */
/*                                                                          */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCtsp_test_pure_simple_cliquetree (int ncount,                      */
/*      CCtsp_lpcut_in *c, int *yes_no)                                     */
/*    TEST if cut is a two handled cliquetree (assumes first two cliques    */
/*     cliques in the cut are the handles).                                 */
/*     -ncount is the number of nodes in the graph.                         */
/*     -c points  to the cut.                                               */
/*     -yes_no will be set to either 0 or 1, with 1 meaning yes.            */
/*                                                                          */
/*  int CCtsp_comb_to_cliquetree (CCtsp_lpgraph *g, CC_GCgraph *h,          */
/*      double *x, CCtsp_lpcut_in *c, CCtsp_lpcut_in **d)                   */
/*    ATTEMPTS to build a cliquetree from the comb by adding a bunny.       */
/*     -x is an lp vector (it should match the edge set of the graph g)     */
/*     -c is the comb (it will be tested)                                   */
/*     -d returns the double decker or NULL if none is found                */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "macrorus.h"
#include "util.h"
#include "tsp.h"
#include "combs.h"


static void
    mark_GCgraph_clique (CC_GCgraph *g, CCtsp_lpclique *c, int marker);

static int
    comb_to_cliquetree (CCtsp_lpgraph *g, CC_GCgraph *h, CCtsp_lpclique *handle,
        int nteeth, CCtsp_lpclique **teeth, double *x, CCtsp_lpcut_in **cuts),
    build_cliquetree (CCtsp_lpcut_in **cut, int nteeth, CCtsp_lpclique **teeth,
        CCtsp_lpclique *hand1, CCtsp_lpclique *hand2, int ncount);



int CCtsp_comb_to_cliquetree (CCtsp_lpgraph *g, CC_GCgraph *h, double *x,
        CCtsp_lpcut_in *c, CCtsp_lpcut_in **d)
{
    int rval = 0;
    int i, test, ihandle;
    int nteeth = 0;
    CCtsp_lpclique **teeth = (CCtsp_lpclique **) NULL;
    CCtsp_lpclique *handle;

    rval = CCtsp_test_pure_comb (g->ncount, c, &test, &ihandle);
    if (rval) {
        fprintf (stderr, "CCtsp_test_pure_comb failed\n"); goto CLEANUP;
    }
    if (!test) goto CLEANUP;

    handle = &c->cliques[ihandle];
    teeth = CC_SAFE_MALLOC (c->cliquecount - 1, CCtsp_lpclique *);
    if (teeth == (CCtsp_lpclique **) NULL) {
        fprintf (stderr, "out of memory in CCtsp_comb_to_cliquetree\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < c->cliquecount; i++) {
        if (i != ihandle) {
            teeth[nteeth++] = &c->cliques[i];
        }
    }

    rval = comb_to_cliquetree (g, h, handle, nteeth, teeth, x, d);
    if (rval) {
        fprintf (stderr, "comb_to_cliquetree failed\n");
        goto CLEANUP;
    }

CLEANUP:

    CC_IFFREE (teeth, CCtsp_lpclique *);

    return rval;
}

static int comb_to_cliquetree (CCtsp_lpgraph *g, CC_GCgraph *h,
        CCtsp_lpclique *handle, int nteeth, CCtsp_lpclique **teeth, double *x,
        CCtsp_lpcut_in **cuts)
{
    int rval = 0;
    int nbig = 0;
    int t1[2];
    int t2[2];
    double t1val, t2val;
    int i, j, k, e, n, test;
    CCtsp_lpclique **ctlist = (CCtsp_lpclique **) NULL;
    CCtsp_lpclique **bigteeth = (CCtsp_lpclique **) NULL;
    CCtsp_lpclique bunny;
    CCtsp_lpclique ear1;
    CCtsp_lpclique ear2;
    CCtsp_lpcut_in *dp;
    int *gset = (int *) NULL;
    int *marks = (int *) NULL;
    int *ar   = (int *) NULL;
    int gcount = 0;
    int acount = 0;
    double gval, gbestval;
    int gbest;

    gset  = CC_SAFE_MALLOC (g->ncount, int);
    marks = CC_SAFE_MALLOC (g->ncount, int);
    if (gset == (int *) NULL || marks == (int *) NULL) {
        fprintf (stderr, "out of memory in comb_to_cliquetree\n");
        rval = 1; goto CLEANUP;
    }
        
    CCtsp_mark_clique_and_neighbors (g, handle, marks, 0);
    for (i = 0; i < nteeth; i++) {
        CCtsp_mark_clique_and_neighbors (g, teeth[i], marks, 0);
    }

    bigteeth = CC_SAFE_MALLOC (nteeth, CCtsp_lpclique *);
    if (!bigteeth) {
        fprintf (stderr, "out of memory in comb_to_cliquetree\n");
        rval = 1; goto CLEANUP;
    }

    CCtsp_mark_clique (handle, marks, 1);
    for (i = 0; i < nteeth; i++) {
        CCtsp_clique_marked_count (teeth[i], marks, 0, &test);
        if (test > 1) {
            bigteeth[nbig++] = teeth[i];
        }
    }

    mark_GCgraph_clique (h, handle, 1);
    for (i = 0; i < nteeth; i++) {
        mark_GCgraph_clique (h, teeth[i], 1);
    }
        
    for (k = 0; k < nbig; k++) {
        gbestval = 1000.0;
        gbest = -1;

        CCtsp_mark_clique_and_neighbors (g, handle, marks, 0);
        for (i = 0; i < nteeth; i++) {
            CCtsp_mark_clique_and_neighbors (g, teeth[i], marks, 0);
        }
        for (i = 0; i < nteeth; i++) {
            CCtsp_mark_clique (teeth[i], marks, 1);
        }
        CCtsp_mark_clique (handle, marks, 2);

        rval = CCtsp_clique_to_array (bigteeth[k], &ar, &acount);
        if (rval) {
            fprintf (stderr, "CCtsp_clique_to_array failed\n");
            goto CLEANUP;
        }

        for (i = 0; i < acount; i++) {
            n = ar[i];
            if (marks[n] == 1) {
                gcount = 1;
                gset[0] = n;

                rval = CCcombs_greedy_cut (h, &gcount, gset, 1, 2, 0, 2,
                                           (int *) NULL, &gval);
                if (rval) {
                    fprintf (stderr, "CCcombs_greedy_cut failed\n");
                    goto CLEANUP;
                }
                if (gcount >= 3) {
                    if (gval < gbestval) {
                        gbestval = gval;
                        gbest = n;
                    }
                }
            }
        }
        CC_IFFREE (ar, int);

        if (gbest != -1 && gbestval < 2.5) {
            gcount = 1;
            gset[0] = gbest;

            rval = CCcombs_greedy_cut (h, &gcount, gset, 1, 2, 0, 2,
                                       (int *) NULL, &gval);
            if (rval) {
                fprintf (stderr, "CCcombs_greedy_cut failed\n");
                goto CLEANUP;
            }

            rval = CCtsp_array_to_lpclique (gset, gcount, &bunny);
            if (rval) {
                fprintf (stderr, "CCtsp_array_to_lpclique failed\n");
                goto CLEANUP;
            }

            t1[0] = t1[1] = -1;
            t2[0] = t2[1] = -1;
            t1val = t2val = -1.0;

            CCtsp_mark_clique_and_neighbors (g, &bunny, marks, 0);
            CCtsp_mark_clique (&bunny, marks, 1);
            CCtsp_mark_clique (handle, marks, 2);
            for (i = 0; i < nteeth; i++) {
                CCtsp_mark_clique (teeth[i], marks, 2);
            }

            for (i = 0; i < gcount; i++) {
                n = gset[i];
                if (marks[n] == 1) {
                    for (j = 0; j < g->nodes[n].deg; j++) {
                        e = g->nodes[n].adj[j].edge;
                        if (marks[g->nodes[n].adj[j].to] == 0) {
                            if (x[e] > t1val) {
                                t1val = x[e];
                                t1[0] = n;
                                t1[1] = g->nodes[n].adj[j].to;
                            }
                        }
                    }
                }
            }
            if (t1val >= 0.1) {
                marks[t1[0]] = 2;
                marks[t1[1]] = 2;
                for (i = 0; i < gcount; i++) {
                    n = gset[i];
                    if (marks[n] == 1) {
                        for (j = 0; j < g->nodes[n].deg; j++) {
                            e = g->nodes[n].adj[j].edge;
                            if (marks[g->nodes[n].adj[j].to] == 0) {
                                if (x[e] > t2val) {
                                    t2val = x[e];
                                    t2[0] = n;
                                    t2[1] = g->nodes[n].adj[j].to;
                                }
                            }
                        }
                    }
                }
            }

            if (t1val >= 0.1 && t2val >= 0.1) {
                rval = CCtsp_array_to_lpclique (t1, 2, &ear1);
                if (rval) {
                    fprintf (stderr, "CCtsp_array_to_lpclique failed\n");
                    goto CLEANUP;
                }
                rval = CCtsp_array_to_lpclique (t2, 2, &ear2);
                if (rval) {
                    fprintf (stderr, "CCtsp_array_to_lpclique failed\n");
                    goto CLEANUP;
                }

                ctlist = CC_SAFE_MALLOC (nteeth + 2, CCtsp_lpclique *);
                if (!ctlist) {
                    fprintf (stderr, "out of memory in comb2ddecker\n");
                    CCtsp_free_lpclique (&ear1);
                    CCtsp_free_lpclique (&ear2);
                    CCtsp_free_lpclique (&bunny);
                    goto CLEANUP;
                }
                for (i = 0; i < nteeth; i++) {
                    ctlist[i] = teeth[i];
                }
                ctlist[i++] = &ear1;
                ctlist[i]   = &ear2;

                rval = build_cliquetree (&dp, nteeth + 2, ctlist, handle,
                                         &bunny, g->ncount);
                if (rval) {
                    fprintf (stderr, "build_cliquegtee failed\n");
                    CCtsp_free_lpclique (&ear1);
                    CCtsp_free_lpclique (&ear2);
                    CCtsp_free_lpclique (&bunny);
                    goto CLEANUP;
                }

                rval = CCtsp_test_pure_simple_cliquetree (g->ncount, dp, &test);
                if (rval) {
                    fprintf (stderr, "test_pure_simple_cliquetree failed\n");
                    CCtsp_free_lpclique (&ear1);
                    CCtsp_free_lpclique (&ear2);
                    CCtsp_free_lpclique (&bunny);
                    CCtsp_free_lpcut_in (dp);
                    goto CLEANUP;
                }
                if (test == 0) {
                    fprintf (stderr, "clique tree did not pass test\n");
                    CCtsp_free_lpclique (&ear1);
                    CCtsp_free_lpclique (&ear2);
                    CCtsp_free_lpclique (&bunny);
                    CCtsp_free_lpcut_in (dp);
                    rval = 1; goto CLEANUP;
                }


                CCtsp_free_lpclique (&ear1);
                CCtsp_free_lpclique (&ear2);

                dp->next = *cuts;
                *cuts = dp;
                CC_IFFREE (ctlist, CCtsp_lpclique *);
            }
            CCtsp_free_lpclique (&bunny);
        }
    }

CLEANUP:

    mark_GCgraph_clique (h, handle, 0);
    for (i = 0; i < nteeth; i++) {
        mark_GCgraph_clique (h, teeth[i], 0);
    }

    CC_IFFREE (gset, int);
    CC_IFFREE (marks, int);
    CC_IFFREE (ar, int);
    CC_IFFREE (bigteeth, CCtsp_lpclique *);
    CC_IFFREE (ctlist, CCtsp_lpclique *);
    return rval;
}

static int build_cliquetree (CCtsp_lpcut_in **cut, int nteeth,
        CCtsp_lpclique **teeth, CCtsp_lpclique *hand1,
        CCtsp_lpclique *hand2, int ncount)
{
    int rval = 0;
    CCtsp_lpcut_in *dp;
    int i, j;

    *cut = (CCtsp_lpcut_in *) NULL;

    dp = CC_SAFE_MALLOC (1, CCtsp_lpcut_in);
    CCcheck_NULL (dp, "out of memory in build_cliquetree");
    CCtsp_init_lpcut_in (dp);


    dp->cliques = CC_SAFE_MALLOC (nteeth + 2, CCtsp_lpclique);
    if (!dp->cliques) {
        fprintf (stderr, "out of memory in build_cliquetree\n");
        CC_FREE (dp, CCtsp_lpcut_in);
        rval = 1; goto CLEANUP;
    }

    rval = CCtsp_copy_lpclique (hand1, &dp->cliques[0]);
    if (rval) {
        fprintf (stderr, "CCtsp_copy_lpclique failed\n");
        CC_FREE (dp, CCtsp_lpcut_in);
        goto CLEANUP;
    }
    rval = CCtsp_copy_lpclique (hand2, &dp->cliques[1]);
    if (rval) {
        fprintf (stderr, "CCtsp_copy_lpclique failed\n");
        CCtsp_free_lpclique (&dp->cliques[0]);
        CC_FREE (dp, CCtsp_lpcut_in);
        goto CLEANUP;
    }

    for (i = 0; i < nteeth; i++) {
        rval = CCtsp_copy_lpclique (teeth[i], &dp->cliques[i+2]);
        if (rval) {
            fprintf (stderr, "CCtsp_copy_lpclique failed\n");
            for (j = 0; j < i+2; j++) {
                CCtsp_free_lpclique (&dp->cliques[j]);
            }
            CC_FREE (dp->cliques, CCtsp_lpclique);
            CC_FREE (dp, CCtsp_lpcut_in);
            goto CLEANUP;
        }
    }


    dp->cliquecount = nteeth + 2;
    dp->rhs = (2 * (nteeth + 2)) + (nteeth - 1);
    dp->sense = 'G';

    rval = CCtsp_construct_skeleton (dp, ncount);
    if (rval) {
        fprintf (stderr, "CCtsp_construct_skeleton failed\n");
        CCtsp_free_lpcut_in (dp);
        goto CLEANUP;
    }

    *cut = dp;

CLEANUP:

    return rval;
}

int CCtsp_test_pure_simple_cliquetree (int ncount, CCtsp_lpcut_in *c,
       int *rtest)
{
    int rval = 0;
    int *marks = (int *) NULL;
    int i, j, k, test, test2, rhs;

    /* Assumes first two cliques are the handles */

    *rtest = 1;

    marks = CC_SAFE_MALLOC (ncount, int);
    if (marks == (int *) NULL) {
        fprintf (stderr, "out of memory in test_pure_simple_cliquetree\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < c->cliquecount; i++) {
        CCtsp_mark_clique (&c->cliques[i], marks, 0);
    }

    /* Check that teeth are disjoint */

    for (i = 2; i < c->cliquecount; i++) {
        CCtsp_clique_marked_count (&c->cliques[i], marks, 1, &test);
        if (test > 0) {
            fprintf (stderr, "teeth are not disjoint\n");
            *rtest = 0; goto CLEANUP;
        }
        CCtsp_mark_clique (&c->cliques[i], marks, 1);
    }
    for (i = 2; i < c->cliquecount; i++) {
        CCtsp_mark_clique (&c->cliques[i], marks, 0);
    }

    /* Check that handles are disjoint */

    CCtsp_mark_clique (&c->cliques[0], marks, 1);
    CCtsp_clique_marked_count (&c->cliques[1], marks, 1, &test);
    if (test > 0) {
        fprintf (stderr, "handles are not disjoint\n");
        *rtest = 0; goto CLEANUP;
    }
    CCtsp_mark_clique (&c->cliques[0], marks, 0);

    /* Check that each tooth has a cavity */

    CCtsp_mark_clique (&c->cliques[0], marks, 1);
    CCtsp_mark_clique (&c->cliques[1], marks, 1);
    for (i = 2; i < c->cliquecount; i++) {
        CCtsp_clique_marked_count (&c->cliques[i], marks, 0, &test);
        if (test == 0) {
            fprintf (stderr, "tooth has no cavity\n");
            *rtest = 0; goto CLEANUP;
        }
    }
    CCtsp_mark_clique (&c->cliques[0], marks, 0);
    CCtsp_mark_clique (&c->cliques[1], marks, 0);

    /* Check that each handle intersects an odd number of teeth */

    for (j = 0; j <= 1; j++) {
         CCtsp_mark_clique (&c->cliques[j], marks, 1);
         k = 0;
         for (i = 2; i < c->cliquecount; i++) {
             CCtsp_clique_marked_count (&c->cliques[i], marks, 0, &test);
             if (test > 0) k++;
         }
         if (k % 2 == 0) {
             fprintf (stderr, "handle meets even number of teeth\n");
             *rtest = 0; goto CLEANUP;
         }
         if (k < 3) {
             fprintf (stderr, "handle meets only %d teeth\n", k);
             *rtest = 0; goto CLEANUP;
         }
         CCtsp_mark_clique (&c->cliques[j], marks, 0);
    }

    /* Check that structure is connected (i.e. there is a nonpendent tooth) */

    CCtsp_mark_clique (&c->cliques[0], marks, 1);
    CCtsp_mark_clique (&c->cliques[1], marks, 2);
    k = 0;
    for (i = 2; i < c->cliquecount; i++) {
        CCtsp_clique_marked_count (&c->cliques[i], marks, 1, &test);
        CCtsp_clique_marked_count (&c->cliques[i], marks, 2, &test2);
        if (test > 0 && test2 > 0) k++;
    }
    if (k != 1) {
        fprintf (stderr, "%d nonpendent teeth\n", k);
        *rtest = 0; goto CLEANUP;
    }

    /* RHS should be 2*ncliques + (nteeth-1) */

    rhs = (2 * c->cliquecount) + (c->cliquecount - 3);
    if (rhs != c->rhs) {
        fprintf (stderr, "rhs value is wrong\n");
        *rtest = 0; goto CLEANUP;
    }

CLEANUP:

    if (rval) *rtest = 0;
    CC_IFFREE (marks, int);
    return rval;
}

static void mark_GCgraph_clique (CC_GCgraph *g, CCtsp_lpclique *c, int marker)
{
    int j, tmp;

    CC_FOREACH_NODE_IN_CLIQUE (j, *c, tmp) {
        g->nodelist[j].mark = marker;
    }
}
