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
/*                    SEARCH FOR DOUBLE DECKERS                             */
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
/*  int CCtsp_test_pure_double_decker (CCtsp_lpcut_in *c, int *yes_no,      */
/*      int *handle1, int *handle2)                                         */
/*    TESTS if cut is a pure double decker (no flips permitted)             */
/*     -yes_no will be set to either 0 or 1, with 1 meaning yes.            */
/*     -handle1 and handle2 will be set to the handles (they can be NULL)   */
/*                                                                          */
/*  int CCtsp_comb_to_double_decker (CCtsp_lpgraph *g, CC_GCgraph *h,       */
/*      double *x, CCtsp_lpcut_in *c, CCtsp_lpcut_in **d)                   */
/*    ATTEMPTS to build a double decker from the comb by adding a second    */
/*     handle.                                                              */
/*     -x is an lp vector (it should match the edge set of the graph g)     */
/*     -c is the comb (it will be tested)                                   */
/*     -d returns a NULL terminated list of any new double deckers that     */
/*      were found                                                          */
/*                                                                          */
/*  int CCtsp_comb_to_star (CCtsp_lpgraph *g, CC_GCgraph *h, double *x,     */
/*      CCtsp_lpcut_in *c, CCtsp_lpcut_in **d)                              */
/*    ATTEMPTS to build a star inequalite from the comb using the           */
/*     lcm method of Naddef and Thienel.                                    */
/*     -x is an lp vector (it should match the edge set of the graph g)     */
/*     -c is the comb (it will be tested)                                   */
/*     -d returns a NULL terminated list of any new stars that were found   */
/*                                                                          */
/*  int CCtsp_comb_handling (CCtsp_lpgraph *g, CC_GCgraph *h, double *x,    */
/*      CCtsp_lpcut_in *c, CCtsp_lpcut_in **d)                              */
/*    ATTEMPTS to build star inequalities by adding nested handles and      */
/*     possibly stretching the teeth.                                       */
/*     -arguments are the same as above.                                    */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "macrorus.h"
#include "util.h"
#include "tsp.h"
#include "combs.h"
#include "verify.h"

#undef  TRY_QUAD_DECKERS


static void
    add_x_to_clique_neighbors (CCtsp_lpgraph *g, double *x, CCtsp_lpclique *c,
        double *values),
    add_x_to_node_neighbors (CCtsp_lpgraph *g, double *x,
        int n, double *values),
    subtract_x_from_node_neighbors (CCtsp_lpgraph *g, double *x,
        int n, double *values),
    divide_values_in_clique (CCtsp_lpclique *c, double *values,
        double divider),
    find_smallest_marked (CCtsp_lpclique *c, double *values, int *marks,
        int marker, int *bestnode, double *bestval),
    find_largest_marked (CCtsp_lpclique *c, double *values, int *marks,
        int marker, int *bestnode, double *bestval),
    mark_GCgraph_clique (CC_GCgraph *g, CCtsp_lpclique *c, int marker);

static int
    comb_to_ddecker (CCtsp_lpgraph *g, CCtsp_lpclique *handle, int nteeth,
        CCtsp_lpclique **teeth, double *x, CCtsp_lpcut_in **cuts),
    comb_to_star (CCtsp_lpgraph *g, CCtsp_lpclique *handle, int nteeth,
        CCtsp_lpclique **teeth, double *x, CCtsp_lpcut_in **cuts),
    comb_handling (CCtsp_lpgraph *g, CC_GCgraph *h, CCtsp_lpclique *handle,
        int nteeth, CCtsp_lpclique **teeth, double *x, CCtsp_lpcut_in **cuts),
    build_star (CCtsp_lpcut_in **cut, int nhandles, CCtsp_lpclique **handles,
        int *alpha, int nteeth, CCtsp_lpclique **teeth, int *beta, int ncount),
    stretch_teeth (CCtsp_lpgraph *g, CC_GCgraph *h,
        CCtsp_lpclique *handle, int nteeth, CCtsp_lpclique **teeth,
        CCtsp_lpclique **bigcliques, int forced, int *targets, int *test),
    grow_handle (CCtsp_lpgraph *g, double *x, CCtsp_lpclique *handle,
        int nteeth, CCtsp_lpclique **teeth, CCtsp_lpclique **newhandle,
        int *marks),
    add_or_subtract_nodes (CCtsp_lpclique *old, CCtsp_lpclique *new,
        int icount, int *ind, int add_or_sub);



int CCtsp_test_pure_double_decker (CCtsp_lpcut_in *c, int *yes_no,
        int *handle1, int *handle2)
{
    int rval = 0;
    int ccount = c->cliquecount;
    CCtsp_lpclique *cliques = c->cliques;
    int *hashvalues = (int *) NULL;
    int *hperm = (int *) NULL;
    int *marks = (int *) NULL;
    int *rcliques = (int *) NULL;
    int *rmult = (int *) NULL;
    int i, j, k, hk, alpha, beta, test, marked, marked2, rhs, maxn, rcount;
    int hand0, hand1;

    *yes_no = 0;
    if (handle1) *handle1 = -1;
    if (handle2) *handle2 = -1;

    if (ccount < 5 || c->sense != 'G') {
        printf ("wrong ccount or sense in ddecker\n"); fflush (stdout);
        goto CLEANUP;
    }

    /* find the duplicated cliques (using a hash function) */

    hashvalues = CC_SAFE_MALLOC (ccount, int);
    hperm      = CC_SAFE_MALLOC (ccount, int);
    if (!hashvalues || !hperm) {
        fprintf (stderr, "out of memory in CCtsp_test_pure_double_decker\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < ccount; i++) {
        hashvalues[i] =  (int) CCtsp_hashclique (&cliques[i]);
        hperm[i] = i;
    }
    CCutil_int_perm_quicksort (hperm, hashvalues, ccount);

    rcliques = CC_SAFE_MALLOC (ccount, int);
    rmult    = CC_SAFE_MALLOC (ccount, int);
    if (!rcliques || !rmult) {
        fprintf (stderr, "out of memory in CCtsp_test_pure_double_decker\n");
        rval = 1; goto CLEANUP;
    }

    rcliques[0] = hperm[0];
    rmult[0] = 1;
    rcount = 0;
    for (i = 1; i < ccount; i++) {
        if (hashvalues[hperm[i]] != hashvalues[rcliques[rcount]]) {
            rcliques[++rcount] = hperm[i];
            rmult[rcount] = 1;
        } else {
            int tcount = rcount;
            do {
                CCtsp_clique_eq (&cliques[rcliques[tcount]],
                                 &cliques[hperm[i]], &test);
                if (test) break;
                tcount--;
            } while (tcount >= 0 && hashvalues[hperm[i]] ==
                                    hashvalues[rcliques[tcount]]);

            if (test) {
                rmult[tcount]++;
            } else {
                rcliques[++rcount] = hperm[i];
                rmult[rcount] = 1;
            }
        }
    }
    rcount++;

    /* must have an odd number of distinct cliques */

    if (rcount % 2 != 1 || rcount < 5) {
        printf ("either an even number of distinct cliques or too few: %d\n",
                 rcount);
        fflush (stdout);
        for (i = 0; i < rcount; i++) {
            printf ("Clique %d [%d]: ", i, hashvalues[rcliques[i]]);
            CCtsp_print_lpclique (&cliques[rcliques[i]]);
        }
        goto CLEANUP;
    }

    /* no clique can appear more than twice */

    for (i = 0; i < rcount; i++) {
        if (rmult[i] > 2) {
            printf ("some clique appears more than twice\n"); fflush (stdout);
            goto CLEANUP;
        }
    }

    /* get a marks array for the cliques */

    maxn = 0;
    for (i = 0; i < rcount; i++) {
        for (j = 0; j < cliques[rcliques[i]].segcount; j++) {
            if (cliques[rcliques[i]].nodes[j].hi > maxn) {
                maxn = cliques[rcliques[i]].nodes[j].hi;
            }
        }
    }
    marks = CC_SAFE_MALLOC (maxn + 1, int);
    if (!marks) {
        fprintf (stderr, "out of memory in CCtsp_test_pure_double_decker\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < rcount; i++) {
        CCtsp_mark_clique (&cliques[rcliques[i]], marks, 0);
    }

    /* find a potential handle */

    hand0 = -1;
    CCtsp_mark_clique (&cliques[rcliques[0]], marks, 1);
    CCtsp_is_clique_marked (&cliques[rcliques[1]], marks, 1, &test);
    if (test) {
        CCtsp_is_clique_marked (&cliques[rcliques[2]], marks, 1, &marked);
        CCtsp_is_clique_marked (&cliques[rcliques[3]], marks, 1, &marked2);
        if (marked && marked2) {
            hand0 = 0;
        } else {
            hand0 = 1;
        }
    } else {
        for (i = 2; i < rcount; i++) {
            CCtsp_is_clique_marked (&cliques[rcliques[i]], marks, 1, &marked);
            if (marked) {
                hand0 = i;
                break;
            }
        }
    }
    CCtsp_mark_clique (&cliques[rcliques[0]], marks, 0);
    if (hand0 == -1) {
        printf ("no potential handle in ddecker\n"); fflush (stdout);
        goto CLEANUP;
    }


    /* find a clique that either contains hand0 or is contained in hand0 */

    hand1 = -1;
    CCtsp_mark_clique (&cliques[rcliques[hand0]], marks, 1);
    CCtsp_clique_count (&cliques[rcliques[hand0]], &k);
    for (i = 0; i < rcount; i++) {
        if (i != hand0) {
            CCtsp_is_clique_marked (&cliques[rcliques[i]], marks, 0, &marked);
            if (!marked) {
                hand1 = hand0;
                hand0 = i;
                break;
            } else {
                CCtsp_clique_marked_count (&cliques[rcliques[i]], marks, 1,
                                           &hk);
                if (hk == k) {
                    hand1 = i;
                    break;
                }
            }
        }
    }
    if (hand1 == -1) {
        printf ("no second handle in ddecker\n"); fflush (stdout);
        goto CLEANUP;
    }
    CCtsp_mark_clique (&cliques[rcliques[hand1]], marks, 0);

    /* hand0 is inner handle, and hand1 is outer handle */

    /* check that each tooth meets the inner handle */

    CCtsp_mark_clique (&cliques[rcliques[hand0]], marks, 1);
    for (i = 0; i < rcount; i++) {
        if (i != hand0 && i != hand1) {
            CCtsp_is_clique_marked (&cliques[rcliques[i]], marks, 1, &marked);
            if (!marked) {
                printf ("tooth does not meet inner handle\n");
                fflush (stdout);
                goto CLEANUP;
            }
        }
    }
    CCtsp_mark_clique (&cliques[rcliques[hand0]], marks, 0);


    /* check that each tooth has node outside of outer handle */

    CCtsp_mark_clique (&cliques[rcliques[hand1]], marks, 1);
    for (i = 0; i < rcount; i++) {
        if (i != hand0 && i != hand1) {
            CCtsp_is_clique_marked (&cliques[rcliques[i]], marks, 0, &marked);
            if (!marked) {
                printf ("tooth does not have node outside handles\n");
                fflush (stdout);
                goto CLEANUP;
            }
        }
    }
    CCtsp_mark_clique (&cliques[rcliques[hand1]], marks, 0);


    /* check that teeth are disjoint */

    for (i = 0; i < rcount; i++) {
        if (i != hand0 && i != hand1) {
            CCtsp_is_clique_marked (&cliques[rcliques[i]], marks, 1, &marked);
            if (marked) {
                printf ("teeth are not disjoint\n");
                fflush (stdout);
                goto CLEANUP;
            }
            CCtsp_mark_clique (&cliques[rcliques[i]], marks, 1);
        }
    }
    for (i = 0; i < rcount; i++) {
        if (i != hand0 && i != hand1) {
            CCtsp_mark_clique (&cliques[rcliques[i]], marks, 0);
        }
    }


    /* check that each duplicated tooth has no extra cavity */

    CCtsp_mark_clique (&cliques[rcliques[hand1]], marks, 1);
    CCtsp_mark_clique (&cliques[rcliques[hand0]], marks, 0);
    for (i = 0; i < rcount; i++) {
        if (i != hand0 && i != hand1) {
            CCtsp_is_clique_marked (&cliques[rcliques[i]], marks, 1, &marked);
            if ((marked && rmult[i] != 1) || (!marked && rmult[i] != 2)) {
                printf ("duplicated tooth has an extra cavity\n");
                fflush (stdout);
                goto CLEANUP;
            }
        }
    }


    /* check rhs value */

    k = (rcount - 3) / 2;
    beta = 0;
    for (i = 0; i < rcount; i++) {
        if (i != hand0 && i != hand1) {
            beta += rmult[i];
        }
    }
    alpha = 2;
    rhs = (2 * (k + 1) * alpha) + (2 * beta);

    if (rhs != c->rhs) {
        printf ("bad rhs in ddecker: %d instead of %d\n", c->rhs, rhs);
        fflush (stdout);
        goto CLEANUP;
    }

    *yes_no = 1;
    if (handle1) *handle1 = hand0;
    if (handle2) *handle2 = hand1;


CLEANUP:

    CC_IFFREE (hashvalues, int);
    CC_IFFREE (hperm, int);
    CC_IFFREE (rcliques, int);
    CC_IFFREE (rmult, int);
    CC_IFFREE (marks, int);
    return rval;
}

#define DD_MOVE_TOL   (0.50) /* (1.00) */
#define DD_MOVE_TOL_X (0.10)

int CCtsp_comb_to_double_decker (CCtsp_lpgraph *g, CC_GCgraph *h, double *x,
        CCtsp_lpcut_in *c, CCtsp_lpcut_in **d)
{
    int rval = 0;
    int i, test, ihandle;
    int nteeth = 0;
    CCtsp_lpclique **teeth     = (CCtsp_lpclique **) NULL;
    CCtsp_lpclique **bigteeth  = (CCtsp_lpclique **) NULL;
    CCtsp_lpclique *bigcliques = (CCtsp_lpclique *) NULL;
    CCtsp_lpclique *handle;

    *d = (CCtsp_lpcut_in *) NULL;

    rval = CCtsp_test_pure_comb (g->ncount, c, &test, &ihandle);
    if (rval) {
        fprintf (stderr, "CCtsp_test_pure_comb failed\n"); goto CLEANUP;
    }
    if (!test) goto CLEANUP;

    handle = &c->cliques[ihandle];
    teeth = CC_SAFE_MALLOC (c->cliquecount - 1, CCtsp_lpclique *);
    if (teeth == (CCtsp_lpclique **) NULL) {
        fprintf (stderr, "out of memory in CCtsp_comb_to_double_decker\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < c->cliquecount; i++) {
        if (i != ihandle) {
            teeth[nteeth++] = &c->cliques[i];
        }
    }

    rval = comb_to_ddecker (g, handle, nteeth, teeth, x, d);
    if (rval) {
        fprintf (stderr, "comb_to_ddecker failed\n"); goto CLEANUP;
    }
        
    rval = stretch_teeth (g, h, handle, nteeth, teeth, &bigcliques, 2,
                         (int *) NULL, (int *) NULL);
    if (rval) {
        fprintf (stderr, "strech_teeth failed\n");
        goto CLEANUP;
    }

    bigteeth = CC_SAFE_MALLOC (nteeth, CCtsp_lpclique *);
    if (bigteeth == (CCtsp_lpclique **) NULL) {
        fprintf (stderr, "out of memory in CCtsp_comb_to_double_decker\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < nteeth; i++) {
        bigteeth[i] = &bigcliques[i];
    }

    rval = comb_to_ddecker (g, handle, nteeth, bigteeth, x, d);
    if (rval) {
        fprintf (stderr, "comb_to_ddecker failed\n"); goto CLEANUP;
    }

CLEANUP:

    if (bigcliques) {
        for (i = 0; i < nteeth; i++) {
            CCtsp_free_lpclique (&bigcliques[i]);
        }
        CC_FREE (bigcliques, CCtsp_lpclique);
    }
    CC_IFFREE (teeth, CCtsp_lpclique *);
    CC_IFFREE (bigteeth, CCtsp_lpclique *);

    return rval;
}

static int comb_to_ddecker (CCtsp_lpgraph *g, CCtsp_lpclique *handle,
        int nteeth, CCtsp_lpclique **teeth, double *x, CCtsp_lpcut_in **cuts)
{
    int rval = 0;
    int i, trial, test, tnode, snode;
    int addcount  = 0;
    int addcount2 = 0;
    int addcount_x = 0;
    int subcount  = 0;
    int subcount2 = 0;
    int subcount_x = 0;
    int *wnodes, *wnodes2;
    double tval, sval;
    int *addnodes  = (int *) NULL;
    int *addnodes2 = (int *) NULL;
    int *addnodes_x = (int *) NULL;
    int *subnodes  = (int *) NULL;
    int *subnodes2 = (int *) NULL;
    int *subnodes_x = (int *) NULL;
    int *marks     = (int *) NULL;
    int *beta      = (int *) NULL;
    int alpha[4];     
    double *values = (double *) NULL;
    CCtsp_lpclique newclique, newclique2;
    CCtsp_lpcut_in *dp;
    CCtsp_lpclique *harray[4];
#ifdef TRY_QUAD_DECKERS
    CCtsp_lpclique newclique3;
    double qval;
    int *addnodes3 = (int *) NULL;
    int addcount3 = 0;
    int rnode;
#endif

    for (i = 0; i < 4; i++) {
        alpha[i] = 1;
    }

    marks = CC_SAFE_MALLOC (g->ncount, int);
    if (!marks) {
        fprintf (stderr, "out of memory in comb_to_ddecker\n");
        rval = 1; goto CLEANUP;
    }
    CCtsp_mark_clique_and_neighbors (g, handle, marks, 0);
    for (i = 0; i < nteeth; i++) {
        CCtsp_mark_clique_and_neighbors (g, teeth[i], marks, 0);
    }

    if (nteeth % 2 == 0) {
        fprintf (stderr, "even number of teeth\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < nteeth; i++) {
        CCtsp_clique_marked_count (teeth[i], marks, 1, &test);
        if (test > 0) {
            fprintf (stderr, "teeth are not disjoint\n");
            rval = 1; goto CLEANUP;
        }
        CCtsp_mark_clique (teeth[i], marks, 1);
    }
    for (i = 0; i < nteeth; i++) {
        CCtsp_mark_clique (teeth[i], marks, 0);
    }

    CCtsp_mark_clique (handle, marks, 1);
    for (i = 0; i < nteeth; i++) {
        CCtsp_clique_marked_count (teeth[i], marks, 1, &test);
        if (test == 0) {
            fprintf (stderr, "tooth does not meet handle\n");
            printf ("Tooth %d: ", i); fflush (stdout);
            CCtsp_print_lpclique (teeth[i]);
            printf ("Handle: "); fflush (stdout);
            CCtsp_print_lpclique (handle);
            rval = 1; goto CLEANUP;
        }
        CCtsp_clique_marked_count (teeth[i], marks, 0, &test);
        if (test == 0) {
            fprintf (stderr, "tooth does not have a cavity\n");
            rval = 1; goto CLEANUP;
        }
    }
    CCtsp_mark_clique (handle, marks, 0);

    values = CC_SAFE_MALLOC (g->ncount, double);
    if (!values) {
        fprintf (stderr, "out of memory in comb_to_ddecker\n");
        rval = 1; goto CLEANUP;
    }
    CCtsp_mark_clique_and_neighbors_double (g, handle, values, 0.0);
    for (i = 0; i < nteeth; i++) {
        CCtsp_mark_clique_and_neighbors_double (g, teeth[i], values, 0.0);
    }
    add_x_to_clique_neighbors (g, x, handle, values);
    divide_values_in_clique (handle, values, 2.0);

    addnodes  = CC_SAFE_MALLOC (nteeth, int);
    addnodes2 = CC_SAFE_MALLOC (nteeth, int);
    addnodes_x = CC_SAFE_MALLOC (nteeth, int);
    subnodes  = CC_SAFE_MALLOC (nteeth, int);
    subnodes2 = CC_SAFE_MALLOC (nteeth, int);
    subnodes_x = CC_SAFE_MALLOC (nteeth, int);
    beta      = CC_SAFE_MALLOC (nteeth, int);
    if (!addnodes || !addnodes2 || !addnodes_x || !subnodes || !subnodes2 ||
        !subnodes_x || !beta) {
        fprintf (stderr, "out of memory in comb_to_ddecker\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < nteeth; i++) {
        addnodes[i] = addnodes2[i] = subnodes[i] = subnodes2[i] = -1;
        addnodes_x[i] = subnodes_x[i] = -1;
    }

#ifdef TRY_QUAD_DECKERS
    addnodes3 = CC_SAFE_MALLOC (nteeth, int);
    if (!addnodes3) {
        fprintf (stderr, "out of memory in comb_to_ddecker\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < nteeth; i++) {
        addnodes3[i] = -1;
    }
#endif

    CCtsp_mark_clique (handle, marks, 1);
    for (i = 0; i < nteeth; i++) {
        CCtsp_clique_marked_count (teeth[i], marks, 0, &test);
        if (test > 1) {
            find_largest_marked (teeth[i], values, marks, 0, &tnode, &tval);
            if (tval >= DD_MOVE_TOL) {
                addnodes[i] = tnode;
                addcount++;
                if (test > 2)  {
                    marks[tnode] = 1;
                    find_largest_marked (teeth[i], values, marks, 0, &snode,
                                         &sval);
                    if (sval >= DD_MOVE_TOL) {
                        addnodes2[i] = snode;
                        addcount2++;
#ifdef TRY_QUAD_DECKERS
                        if (test > 3) {
                            marks[snode] = 1; 
                            find_largest_marked (teeth[i], values, marks, 0,
                                                 &rnode, &sval);
                            if (qval >= DD_MOVE_TOL) {
                                addnodes3[i] = rnode;
                                addcount3++;
                            }
                            marks[snode] = 0; 
                        }
#endif
                    }
                    marks[tnode] = 0;
                }
            }
            if (tval >= DD_MOVE_TOL_X) {
                addnodes_x[i] = tnode;
                addcount_x++;
            }
            if (addnodes[i] != -1) {
                add_x_to_node_neighbors (g, x, addnodes[i], values);
            }
        }
    }
    for (i = 0; i < nteeth; i++) {
        if (addnodes[i] != -1) {
            subtract_x_from_node_neighbors (g, x, addnodes[i], values);
        }
    }

    for (i = 0; i < nteeth; i++) {
        CCtsp_clique_marked_count (teeth[i], marks, 1, &test);
        if (test > 1) {
            find_smallest_marked (teeth[i], values, marks, 1, &tnode, &tval);
            tval = 2.0 - tval;
            if (tval >= DD_MOVE_TOL) {
                subnodes[i] = tnode;
                subcount++;
                if (test > 2)  {
                    marks[tnode] = 0;
                    find_smallest_marked (teeth[i], values, marks, 1, &snode,
                                          &sval);
                    sval = 2.0 - sval;
                    if (sval >= DD_MOVE_TOL) {
                        subnodes2[i] = snode;
                        subcount2++;
                    }
                    marks[tnode] = 1;
                }
            }
            if (tval >= DD_MOVE_TOL_X) {
                subnodes_x[i] = tnode;
                subcount_x++;
            }
            if (subnodes[i] != -1) {
                subtract_x_from_node_neighbors (g, x, subnodes[i], values);
            }
        }
    }

    for (trial = 0; trial <= 3; trial++) {
        if (trial == 0) {
            if (addcount == 0) continue;
            rval = add_or_subtract_nodes (handle, &newclique, nteeth,
                                          addnodes, 0);
            if (rval) {
                fprintf (stderr, "add_or_subtract_nodes failed\n");
                goto CLEANUP;
            }
            wnodes = addnodes;
        } else if (trial == 1) {
            if (subcount == 0) continue;
            rval = add_or_subtract_nodes (handle, &newclique, nteeth,
                                          subnodes, 1);
            if (rval) {
                fprintf (stderr, "add_or_subtract_nodes failed\n");
                goto CLEANUP;
            }
            wnodes = subnodes;
        } else if (trial == 2) {
            if (addcount_x <= addcount) continue;
            rval = add_or_subtract_nodes (handle, &newclique, nteeth,
                                          addnodes_x, 0);
            if (rval) {
                fprintf (stderr, "add_or_subtract_nodes failed\n");
                goto CLEANUP;
            }
            wnodes = addnodes_x;
        } else {
            if (subcount_x <= subcount) continue;
            rval = add_or_subtract_nodes (handle, &newclique, nteeth,
                                          subnodes_x, 1);
            if (rval) {
                fprintf (stderr, "add_or_subtract_nodes failed\n");
                goto CLEANUP;
            }
            wnodes = subnodes_x;
        }

        for (i = 0; i < nteeth; i++) {
            beta[i] = (wnodes[i] == -1 ? 2  : 1);
        }

        harray[0] = handle;
        harray[1] = &newclique;
        rval = build_star (&dp, 2, harray, alpha, nteeth, teeth, beta,
                           g->ncount);
        if (rval) {
            fprintf (stderr, "build_star failed\n");
            CCtsp_free_lpclique (&newclique);
            goto CLEANUP;
        }
        CCtsp_free_lpclique (&newclique);

        if (dp) {
            rval = CCtsp_test_pure_double_decker (dp, &test,
                                                 (int *) NULL, (int *) NULL);
            if (rval) {
                fprintf (stderr, "CCtsp_test_pure_double_decker failed\n");
                CCtsp_free_lpcut_in (dp);
                goto CLEANUP;
            }
            if (!test) {
                fprintf (stderr, "ddecker is not valid\n");
                CCtsp_print_lpcut_in (dp);
                rval = 1; goto CLEANUP;
            }

            dp->next = *cuts;
            *cuts = dp;
        }
    }

    /* Try triple deckers */

    for (trial = 0; trial <= 1; trial++) {
        if (trial == 1) {
            if (addcount2 == 0) continue;
            rval = add_or_subtract_nodes (handle, &newclique, nteeth,
                                          addnodes, 0);
            if (rval) {
                fprintf (stderr, "add_or_subtract_nodes failed\n");
                goto CLEANUP;
            }
            rval = add_or_subtract_nodes (&newclique, &newclique2, nteeth,
                                          addnodes2, 0);
            if (rval) {
                fprintf (stderr, "add_or_subtract_nodes failed\n");
                CCtsp_free_lpclique (&newclique);
                goto CLEANUP;
            }
            wnodes = addnodes;
            wnodes2 = addnodes2;
        } else {
            if (subcount2 == 0) continue;
            rval = add_or_subtract_nodes (handle, &newclique, nteeth,
                                          subnodes, 1);
            if (rval) {
                fprintf (stderr, "add_or_subtract_nodes failed\n");
                goto CLEANUP;
            }
            rval = add_or_subtract_nodes (&newclique, &newclique2, nteeth,
                                          subnodes2, 1);
            if (rval) {
                fprintf (stderr, "add_or_subtract_nodes failed\n");
                CCtsp_free_lpclique (&newclique);
                goto CLEANUP;
            }
            wnodes = subnodes;
            wnodes2 = subnodes2;
        }

        for (i = 0; i < nteeth; i++) {
            if (wnodes[i] == -1) {
                beta[i] = 3;
            } else {
                if (wnodes2[i] == -1) {
                    beta[i] = 2;
                } else {
                    beta[i] = 1;
               }
            }
        }

        harray[0] = handle;
        harray[1] = &newclique;
        harray[2] = &newclique2;
        rval = build_star (&dp, 3, harray, alpha, nteeth, teeth, beta,
                           g->ncount);
        if (rval) {
            fprintf (stderr, "build_star failed\n");
            CCtsp_free_lpclique (&newclique);
            CCtsp_free_lpclique (&newclique2);
            goto CLEANUP;
        }
        CCtsp_free_lpclique (&newclique);
        CCtsp_free_lpclique (&newclique2);

        if (dp) {
            dp->next = *cuts;
            *cuts = dp;
        }
    }

#ifdef TRY_QUAD_DECKERS
    if (addcount3 > 0) {
        printf ("Try"); fflush (stdout);
        rval = add_or_subtract_nodes (handle, &newclique, nteeth, addnodes, 0);
        if (rval) {
            fprintf (stderr, "add_or_subtract_nodes failed\n");
            goto CLEANUP;
        }
        rval = add_or_subtract_nodes (&newclique, &newclique2, nteeth,
                                      addnodes2, 0);
        if (rval) {
            fprintf (stderr, "add_or_subtract_nodes failed\n");
            CCtsp_free_lpclique (&newclique);
            goto CLEANUP;
        }
        rval = add_or_subtract_nodes (&newclique2, &newclique3, nteeth,
                                      addnodes3, 0);
        if (rval) {
            fprintf (stderr, "add_or_subtract_nodes failed\n");
            CCtsp_free_lpclique (&newclique);
            CCtsp_free_lpclique (&newclique2);
            goto CLEANUP;
        }

        for (i = 0; i < nteeth; i++) {
            if (addnodes[i] == -1) {
                beta[i] = 4;
            } else {
                if (addnodes2[i] == -1) {
                    beta[i] = 3;
                } else {
                    if (addnodes3[i] == -1) {
                        beta[i] = 2;
                    } else {
                        beta[i] = 1;
                    }
                }
            }
        }

        harray[0] = handle;
        harray[1] = &newclique;
        harray[2] = &newclique2;
        harray[3] = &newclique3;
        rval = build_star (&dp, 4, harray, alpha, nteeth, teeth, beta,
                           g->ncount);
        if (rval) {
            fprintf (stderr, "build_star failed\n");
            CCtsp_free_lpclique (&newclique);
            CCtsp_free_lpclique (&newclique2);
            goto CLEANUP;
        }
        CCtsp_free_lpclique (&newclique);
        CCtsp_free_lpclique (&newclique2);
        CCtsp_free_lpclique (&newclique3);

        if (dp) {
            dp->next = *cuts;
            *cuts = dp;
        }
    }
#endif /* TRY_QUAD_DECKERS */

CLEANUP:

    CC_IFFREE (marks, int);
    CC_IFFREE (values, double);
    CC_IFFREE (addnodes, int);
    CC_IFFREE (addnodes2, int);
    CC_IFFREE (addnodes_x, int);
    CC_IFFREE (subnodes, int);
    CC_IFFREE (subnodes2, int);
    CC_IFFREE (subnodes_x, int);
    CC_IFFREE (beta, int);
#ifdef TRY_QUAD_DECKERS
    CC_IFFREE (addnodes3, int);
#endif

    return rval;
}

#define MAX_STAR_COUNT 100

static int build_star (CCtsp_lpcut_in **cut, int nhandles,
        CCtsp_lpclique **handles, int *alpha, int nteeth,
        CCtsp_lpclique **teeth, int *beta, int ncount)
{
    int rval = 0;
    int betasum = 0, alphasum = 0;
    int i, j, k;
    CCtsp_lpcut_in *dp;

    *cut = (CCtsp_lpcut_in *) NULL;

    for (i = 0; i < nhandles; i++) {
        alphasum += alpha[i];
    }
    for (i = 0; i < nteeth; i++) {
        betasum += beta[i];
    }

    if (alphasum + betasum > MAX_STAR_COUNT) goto CLEANUP;

    dp = CC_SAFE_MALLOC (1, CCtsp_lpcut_in);
    CCcheck_NULL (dp, "out of memory in build_star");
    CCtsp_init_lpcut_in (dp);

    dp->cliques = CC_SAFE_MALLOC (alphasum + betasum, CCtsp_lpclique);
    if (!dp->cliques) {
        fprintf (stderr, "out of memory in build_star\n");
        CC_FREE (dp, CCtsp_lpcut_in);
        rval = 1; goto CLEANUP;
    }

    for (k = 0, i = 0; i < nhandles; i++) {
        for (j = 0; j < alpha[i]; j++) {
            rval = CCtsp_copy_lpclique (handles[i], &dp->cliques[k]);
            if (rval) {
                fprintf (stderr, "CCtsp_copy_lpclique failed\n");
                for (j = 0; j < k; j++) {
                    CCtsp_free_lpclique (&dp->cliques[j]);
                }
                CC_FREE (dp->cliques, CCtsp_lpclique);
                CC_FREE (dp, CCtsp_lpcut_in);
                goto CLEANUP;
            }
            k++;
        }
    }

    for (i = 0; i < nteeth; i++) {
        for (j = 0; j < beta[i]; j++) {
            rval = CCtsp_copy_lpclique (teeth[i], &dp->cliques[k]);
            if (rval) {
                fprintf (stderr, "CCtsp_copy_lpclique failed\n");
                for (j = 0; j < k; j++) {
                    CCtsp_free_lpclique (&dp->cliques[j]);
                }
                CC_FREE (dp->cliques, CCtsp_lpclique);
                CC_FREE (dp, CCtsp_lpcut_in);
                goto CLEANUP;
            }
            k++;
        }
    }

    dp->cliquecount = alphasum + betasum;
    k = nteeth / 2;
    dp->rhs = (2 * (k + 1) * alphasum) + (2 * betasum);
    dp->sense = 'G';

    if (dp->dominocount != 0) {
        printf ("DDECKER Yipes %d\n", dp->dominocount); fflush (stdout);
        exit (1);
    }
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

static int stretch_teeth (CCtsp_lpgraph *g, CC_GCgraph *h,
        CCtsp_lpclique *handle, int nteeth, CCtsp_lpclique **teeth,
        CCtsp_lpclique **bigcliques, int forced, int *targets, int *test)
{
    int rval = 0;
    int i, j, acount, gcount;
    int *gset  = (int *) NULL;
    int *ar    = (int *) NULL;
    double gval;
    CCtsp_lpclique *newteeth = (CCtsp_lpclique *) NULL;

    *bigcliques =  (CCtsp_lpclique *) NULL;
    if (test) *test = 1;

    gset  = CC_SAFE_MALLOC (g->ncount, int);
    if (gset == (int *) NULL) {
        fprintf (stderr, "out of memory in stretch_teeth\n");
        rval = 1; goto CLEANUP;
    }

    mark_GCgraph_clique (h, handle, 1);
    for (i = 0; i < nteeth; i++) {
        mark_GCgraph_clique (h, teeth[i], 1);
    }
        
    newteeth = CC_SAFE_MALLOC (nteeth, CCtsp_lpclique);
    if (newteeth == (CCtsp_lpclique *) NULL) {
        fprintf (stderr, "out of memory in stretch_teeth\n");
        rval = 1; goto CLEANUP;
    }

    for (i = 0; i < nteeth; i++) {
        if (targets == (int *) NULL || targets[i] != 0) {
            rval = CCtsp_clique_to_array (teeth[i], &ar, &acount);
            if (rval) {
                fprintf (stderr, "CCtsp_clique_to_array failed\n");
                for (j = 0; j < i; j++) {
                    CCtsp_free_lpclique (&newteeth[j]);
                }
                CC_FREE (newteeth, CCtsp_lpclique);
                goto CLEANUP;
            }
            for (j = 0; j < acount; j++) {
                gset[j] = ar[j];
            }
            gcount = acount;
            CC_IFFREE (ar, int);
            
            rval = CCcombs_greedy_cut (h, &gcount, gset, 1, forced, 0, forced,
                                       test, &gval);
            if (rval) {
                fprintf (stderr, "CCcombs_greedy_cut failed\n");
                for (j = 0; j < i; j++) {
                    CCtsp_free_lpclique (&newteeth[j]);
                }
                CC_FREE (newteeth, CCtsp_lpclique);
                goto CLEANUP;
            }
            if (test && (*test == 0)) {
                /* not enough room to grow the tooth */
                for (j = 0; j < i; j++) {
                    CCtsp_free_lpclique (&newteeth[j]);
                }
                CC_FREE (newteeth, CCtsp_lpclique);
                goto CLEANUP;
            }

            rval = CCtsp_array_to_lpclique (gset, gcount, &newteeth[i]);
            if (rval) {
                fprintf (stderr, "CCtsp_array_to_lpclique failed\n");
                for (j = 0; j < i; j++) {
                    CCtsp_free_lpclique (&newteeth[j]);
                }
                CC_FREE (newteeth, CCtsp_lpclique);
                goto CLEANUP;
            }
            mark_GCgraph_clique (h, &newteeth[i], 1);
        } else {
            rval = CCtsp_copy_lpclique (teeth[i], &newteeth[i]);
            if (rval) {
                fprintf (stderr, "CCtsp_copy_lpclique failed\n");
                for (j = 0; j < i; j++) {
                    CCtsp_free_lpclique (&newteeth[j]);
                }
                CC_FREE (newteeth, CCtsp_lpclique);
                goto CLEANUP;
            }
        } 
    }
    
    mark_GCgraph_clique (h, handle, 0);
    for (i = 0; i < nteeth; i++) {
        mark_GCgraph_clique (h, &newteeth[i], 0);
    }

    *bigcliques = newteeth;

CLEANUP:

    CC_IFFREE (gset, int);
    CC_IFFREE (ar, int);

    return rval;
}

static int add_or_subtract_nodes (CCtsp_lpclique *old, CCtsp_lpclique *new,
        int icount, int *ind, int add_or_sub)
{
    int *ar = (int *) NULL;
    int count = 0;
    int i;
    int rval = 0;

    ar = CC_SAFE_MALLOC (icount, int);
    if (!ar) {
        fprintf (stderr, "out of memory in add_or_subtract_nodes\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < icount; i++) {
        if (ind[i] != -1) {
            ar[count++] = ind[i];
        }
    }
    if (add_or_sub == 0) {
        rval = CCtsp_add_nodes_to_lpclique (old, new, count, ar);
        if (rval) {
            fprintf (stderr, "CCtsp_add_nodes_to_lpclique failed\n");
            goto CLEANUP;
        }
    } else {
        rval = CCtsp_delete_nodes_from_lpclique (old, new, count, ar);
        if (rval) {
            fprintf (stderr, "CCtsp_delete_nodes_from_lpclique failed\n");
            goto CLEANUP;
        }
    }

CLEANUP:

    CC_IFFREE (ar, int);
    return rval;
}

static void add_x_to_clique_neighbors (CCtsp_lpgraph *g, double *x,
        CCtsp_lpclique *c, double *values)
{
    int j, k, tmp;

    CC_FOREACH_NODE_IN_CLIQUE (j, *c, tmp) {
        for (k = 0; k < g->nodes[j].deg; k++) {
            values[g->nodes[j].adj[k].to] += x[g->nodes[j].adj[k].edge];
        }
    }
}

static void add_x_to_node_neighbors (CCtsp_lpgraph *g, double *x,
        int n, double *values)
{
    int k;

    for (k = 0; k < g->nodes[n].deg; k++) {
        values[g->nodes[n].adj[k].to] += x[g->nodes[n].adj[k].edge];
    }
}

static void subtract_x_from_node_neighbors (CCtsp_lpgraph *g, double *x,
        int n, double *values)
{
    int k;

    for (k = 0; k < g->nodes[n].deg; k++) {
        values[g->nodes[n].adj[k].to] -= x[g->nodes[n].adj[k].edge];
    }
}

static void divide_values_in_clique (CCtsp_lpclique *c, double *values,
               double divider)
{
    int j, tmp;

    CC_FOREACH_NODE_IN_CLIQUE (j, *c, tmp) {
        values[j] = values[j] / divider;
    }
}

static void find_smallest_marked (CCtsp_lpclique *c, double *values,
        int *marks, int marker, int *bestnode, double *bestval)
{
    int j, tmp;

    *bestval = CCtsp_LP_MAXDOUBLE;
    *bestnode = -1;

    CC_FOREACH_NODE_IN_CLIQUE (j, *c, tmp) {
        if (marks[j] == marker) {
            if (values[j] < *bestval) {
                *bestval = values[j];
                *bestnode = j;
            }
        }
    }
}

static void find_largest_marked (CCtsp_lpclique *c, double *values,
        int *marks, int marker, int *bestnode, double *bestval)
{
    int j, tmp;

    *bestval = -CCtsp_LP_MAXDOUBLE;
    *bestnode = -1;

    CC_FOREACH_NODE_IN_CLIQUE (j, *c, tmp) {
        if (marks[j] == marker) {
            if (values[j] > *bestval) {
                *bestval = values[j];
                *bestnode = j;
            }
        }
    }
}

static void mark_GCgraph_clique (CC_GCgraph *g, CCtsp_lpclique *c, int marker)
{
    int j, tmp;

    CC_FOREACH_NODE_IN_CLIQUE (j, *c, tmp) {
        g->nodelist[j].mark = marker;
    }
}


/*****************************  STARS  **********************************/


#define MAX_LCM   (25)     /* don't accept stars with larger alpha sum */
#define DELTA_EPS (0.001)


int CCtsp_comb_to_star (CCtsp_lpgraph *g, CC_GCgraph *h, double *x,
        CCtsp_lpcut_in *c, CCtsp_lpcut_in **d)
{
    int rval = 0, nteeth = 0;
    int i, test, ihandle;
    CCtsp_lpclique **teeth     = (CCtsp_lpclique **) NULL;
    CCtsp_lpclique **bigteeth  = (CCtsp_lpclique **) NULL;
    CCtsp_lpclique *bigcliques = (CCtsp_lpclique *) NULL;
    CCtsp_lpclique *handle;

    *d = (CCtsp_lpcut_in *) NULL;

    rval = CCtsp_test_pure_comb (g->ncount, c, &test, &ihandle);
    if (rval) {
        fprintf (stderr, "CCtsp_test_pure_comb failed\n"); goto CLEANUP;
    }
    if (!test) goto CLEANUP;

    handle = &c->cliques[ihandle];
    teeth = CC_SAFE_MALLOC (c->cliquecount - 1, CCtsp_lpclique *);
    if (teeth == (CCtsp_lpclique **) NULL) {
        fprintf (stderr, "out of memory in CCtsp_star_to_double_decker\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < c->cliquecount; i++) {
        if (i != ihandle) {
            teeth[nteeth++] = &c->cliques[i];
        }
    }

    rval = comb_to_star (g, handle, nteeth, teeth, x, d);
    if (rval) {
        fprintf (stderr, "comb_to_star failed\n"); goto CLEANUP;
    }
        
    rval = stretch_teeth (g, h, handle, nteeth, teeth, &bigcliques, 2,
                         (int *) NULL, (int *) NULL);
    if (rval) {
        fprintf (stderr, "strech_teeth failed\n");
        goto CLEANUP;
    }

    bigteeth = CC_SAFE_MALLOC (nteeth, CCtsp_lpclique *);
    if (bigteeth == (CCtsp_lpclique **) NULL) {
        fprintf (stderr, "out of memory in CCtsp_star_to_double_decker\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < nteeth; i++) {
        bigteeth[i] = &bigcliques[i];
    }

    rval = comb_to_star (g, handle, nteeth, bigteeth, x, d);
    if (rval) {
        fprintf (stderr, "comb_to_star failed\n"); goto CLEANUP;
    }

CLEANUP:

    if (bigcliques) {
        for (i = 0; i < nteeth; i++) {
            CCtsp_free_lpclique (&bigcliques[i]);
        }
        CC_FREE (bigcliques, CCtsp_lpclique);
    }
    CC_IFFREE (teeth, CCtsp_lpclique *);
    CC_IFFREE (bigteeth, CCtsp_lpclique *);

    return rval;
}

static int comb_to_star (CCtsp_lpgraph *g, CCtsp_lpclique *handle,
        int nteeth, CCtsp_lpclique **teeth, double *x, CCtsp_lpcut_in **cuts)
{
    int i, j, k, bestnode, hits, al, acount;
    int nbig = 0, handlenum = 0, rval = 0, lcm = 0;
    int *marks     = (int *) NULL;
    int *csize     = (int *) NULL;
    int *bsize     = (int *) NULL;
    int *wsize     = (int *) NULL;
    int *beta      = (int *) NULL;
    int *abeta     = (int *) NULL;
    int *alpha     = (int *) NULL;
    int *tload     = (int *) NULL;
    int *tadd      = (int *) NULL;
    int *buse      = (int *) NULL;
    int *ar        = (int *) NULL;
    int **norder   = (int **) NULL;
    int **nwinners = (int **) NULL;
    double *values   = (double *) NULL;
    double *bigdelta = (double *) NULL;
    double **nvalues = (double **) NULL;
    double bestval;
    CCtsp_lpclique **bigteeth = (CCtsp_lpclique **) NULL;
    CCtsp_lpclique **ghandles = (CCtsp_lpclique **) NULL;
    CCtsp_lpclique *handles   = (CCtsp_lpclique *) NULL;
    CCtsp_lpcut_in *dp;


    marks  = CC_SAFE_MALLOC (g->ncount, int);
    values = CC_SAFE_MALLOC (g->ncount, double);
    if (!marks || !values) {
        fprintf (stderr, "out of memory in comb_to_star\n");
        rval = 1; goto CLEANUP;
    }
    CCtsp_mark_clique_and_neighbors (g, handle, marks, 0);
    CCtsp_mark_clique_and_neighbors_double (g, handle, values, 0.0);
    for (i = 0; i < nteeth; i++) {
        CCtsp_mark_clique_and_neighbors (g, teeth[i], marks, 0);
        CCtsp_mark_clique_and_neighbors_double (g, teeth[i], values, 0.0);
    }

    csize  = CC_SAFE_MALLOC (nteeth, int);
    if (!csize) {
        fprintf (stderr, "out of memory in comb_to_star\n");
        rval = 1; goto CLEANUP;
    }

    CCtsp_mark_clique (handle, marks, 1);
    for (i = 0; i < nteeth; i++) {
        CCtsp_clique_marked_count (teeth[i], marks, 0, &(csize[i]));
    }
    for (i = 0, nbig = 0; i < nteeth; i++) {
        if (csize[i] >= 2) nbig++;
    }
    if (nbig == 0) goto CLEANUP;

    bigteeth = CC_SAFE_MALLOC (nbig, CCtsp_lpclique *);
    bsize    = CC_SAFE_MALLOC (nbig, int);
    wsize    = CC_SAFE_MALLOC (nbig, int);
    beta     = CC_SAFE_MALLOC (nbig, int);
    buse     = CC_SAFE_MALLOC (nbig, int);
    tadd     = CC_SAFE_MALLOC (nbig, int);
    tload    = CC_SAFE_MALLOC (nbig, int);
    bigdelta = CC_SAFE_MALLOC (nbig, double);
    if (!bigteeth || !bsize || !beta || !buse || !tadd || !tload ||
        !bigdelta || !wsize) {
        fprintf (stderr, "out of memory in comb_to_star\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0, nbig = 0; i < nteeth; i++) {
        if (csize[i] >= 2) {
            bigteeth[nbig] = teeth[i];
            rval = CCtsp_clique_delta (g, x, bigteeth[nbig], &bigdelta[nbig]);
            if (rval) {
                 fprintf (stderr, "CCtsp_clique_delta failed\n"); goto CLEANUP;
            }
            bsize[nbig++]  = csize[i];
        }
    }

    norder   = CC_SAFE_MALLOC (nbig, int *);
    nwinners = CC_SAFE_MALLOC (nbig, int *);
    nvalues  = CC_SAFE_MALLOC (nbig, double *);
    if (!norder || !nwinners || !nvalues) {
        fprintf (stderr, "out of memory in comb_to_star\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < nbig; i++) {
        norder[i]   = (int *) NULL;
        nwinners[i] = (int *) NULL;
        nvalues[i]  = (double *) NULL;
    }
    for (i = 0; i < nbig; i++) {
        norder[i]   = CC_SAFE_MALLOC (bsize[i], int);
        nvalues[i]  = CC_SAFE_MALLOC (bsize[i], double);
        nwinners[i] = CC_SAFE_MALLOC (bsize[i], int);
        if (!norder[i] || !nvalues[i] || !nwinners[i]) {
            fprintf (stderr, "out of memory in comb_to_star\n");
            rval = 1; goto CLEANUP;
        }
    }

    add_x_to_clique_neighbors (g, x, handle, values);
    for (i = 0; i < nbig; i++) {
        double win, delta;
        for (j = 0; j < bsize[i]; j++) {
            find_largest_marked (bigteeth[i], values, marks, 0, &bestnode,
                                 &bestval);
            add_x_to_node_neighbors (g, x, bestnode, values);
            marks[bestnode] = 1;
            norder[i][j] = bestnode;
            nvalues[i][j] = bestval;
        }
        for (j = 0; j < bsize[i]; j++) {
            subtract_x_from_node_neighbors (g, x, norder[i][j], values);
            marks[norder[i][j]] = 0;
        }

        wsize[i] = 0;
        win = 0.0;
        if (bigdelta[i] > 2.0 + DELTA_EPS) delta = -DELTA_EPS - 0.5;
        else                               delta =  DELTA_EPS;

        for (j = 0; j < bsize[i] - 1; j++) {
            win += (nvalues[i][j] - 1.0);
            if (win > delta) {
                (wsize[i])++;
                nwinners[i][j] = 1;
                win = 0.0;
            } else {
                nwinners[i][j] = 0;
            }
        }
    }

    lcm = wsize[0] + 1;
    for (i = 1; i < nbig; i++) {
        lcm = CCutil_our_lcm (lcm, wsize[i] + 1);
        if (lcm > MAX_LCM) goto CLEANUP;
    }
    if (lcm == 1) goto CLEANUP;

    for (i = 0; i < nbig; i++) {
        beta[i] = lcm / (wsize[i] + 1);
    }

    handles = CC_SAFE_MALLOC (lcm, CCtsp_lpclique);
    alpha   = CC_SAFE_MALLOC (lcm, int);
    if (!handles || !alpha) {
        fprintf (stderr, "out of memory in comb_to_star\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < lcm; i++) {
        CCtsp_init_lpclique (&(handles[i]));
    }

    for (i = 0, j = 0; i < nbig; i++) {
        j += (bsize[i] - 1);
    }
    ar = CC_SAFE_MALLOC (j, int);
    if (!ar) {
        fprintf (stderr, "out of memory in comb_to_star\n");
        rval = 1; goto CLEANUP;
    }


    for (i = 0; i < nbig; i++) {
        tload[i] = tadd[i] = buse[i] = 0;
    }
    hits = 0;

    while (hits < lcm) {
        acount = 0;
        for (i = 0; i < nbig; i++) {
            if (tadd[i] == 1) {
                do {
                    ar[acount++] = norder[i][buse[i]];
                    (buse[i])++;
                } while (nwinners[i][buse[i]-1] != 1);
                tadd[i] = 0;
            }
        }
        al = lcm + 1;
        for (i = 0; i < nbig; i++) {
            if (beta[i] - tload[i] < al) {
                al = beta[i] - tload[i];
            }
        }

        if (handlenum == 0) {
            rval = CCtsp_copy_lpclique (handle, &(handles[0]));
            if (rval) {
                fprintf (stderr, "CCtsp_copy_lpclique failed\n"); goto CLEANUP;
            }
        } else {
            rval = CCtsp_add_nodes_to_lpclique (&(handles[handlenum-1]),
                                       &(handles[handlenum]), acount, ar);
            if (rval) {
                fprintf (stderr, "CCtsp_add_nodes_to_lpclique failed\n");
                goto CLEANUP;
            }
        }
        alpha[handlenum++] = al;

        for (i = 0; i < nbig; i++) {
            tload[i] += al;
            if (tload[i] == beta[i]) {
                tload[i] = 0;
                tadd[i]  = 1;
            }
        }
        hits += al;
    }

    abeta = CC_SAFE_MALLOC (nteeth, int);
    if (!abeta) {
        fprintf (stderr, "out of memory in comb_to_star\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0, k = 0; i < nteeth; i++) {
        if (csize[i] >= 2) {
            abeta[i] = beta[k++];
        } else {
            abeta[i] = lcm;
        }
    }

    ghandles = CC_SAFE_MALLOC (handlenum, CCtsp_lpclique *);
    if (!ghandles) {
        fprintf (stderr, "out of memory in comb_to_star\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < handlenum; i++) {
        ghandles[i] = &(handles[i]);
    }

    rval = build_star (&dp, handlenum, ghandles, alpha, nteeth, teeth, abeta,
                       g->ncount); 
    if (rval) {
        fprintf (stderr, "build_star failed\n"); goto CLEANUP;
    }

    if (dp) {
        rval = CCverify_cut (dp, CC_TYPE_STAR, (int *) NULL);
        if (rval) {
            printf ("Bad star\n"); fflush (stdout);
            rval = 0; goto CLEANUP;
        }
        dp->next = *cuts;
        *cuts = dp;
    }

CLEANUP:
 
    CC_IFFREE (bigteeth, CCtsp_lpclique *);
    CC_IFFREE (ghandles, CCtsp_lpclique *);
    CC_IFFREE (marks, int);
    CC_IFFREE (csize, int);
    CC_IFFREE (bsize, int);
    CC_IFFREE (wsize, int);
    CC_IFFREE (beta, int);
    CC_IFFREE (abeta, int);
    CC_IFFREE (alpha, int);
    CC_IFFREE (ar, int);
    CC_IFFREE (buse, int);
    CC_IFFREE (tadd, int);
    CC_IFFREE (tload, int);
    CC_IFFREE (ar, int);
    CC_IFFREE (bigdelta, double);
    CC_IFFREE (values, double);
    if (nwinners) {
        for (i = 0; i < nbig; i++) {
            CC_IFFREE (nwinners[i], int);
        }
        CC_IFFREE (nwinners, int *);
    }
    if (nvalues) {
        for (i = 0; i < nbig; i++) {
            CC_IFFREE (nvalues[i], double);
        }
        CC_IFFREE (nvalues, double *);
    }
    if (norder) {
        for (i = 0; i < nbig; i++) {
            CC_IFFREE (norder[i], int);
        }
        CC_IFFREE (norder, int *);
    }
    if (handles) {
        for (i = 0; i < lcm; i++) {
            CCtsp_free_lpclique (&(handles[i]));
        }
        CC_IFFREE (handles, CCtsp_lpclique);
    }

    return rval;
}

int CCtsp_comb_handling (CCtsp_lpgraph *g, CC_GCgraph *h, double *x,
        CCtsp_lpcut_in *c, CCtsp_lpcut_in **d)
{
    int rval = 0;
    int i, test, ihandle, nteeth;
    CCtsp_lpclique **teeth = (CCtsp_lpclique **) NULL;
    CCtsp_lpclique *handle;

    *d = (CCtsp_lpcut_in *) NULL;

    rval = CCtsp_test_pure_comb (g->ncount, c, &test, &ihandle);
    if (rval) {
        fprintf (stderr, "CCtsp_test_pure_comb failed\n"); goto CLEANUP;
    }
    if (!test) goto CLEANUP;

    handle = &c->cliques[ihandle];
    teeth = CC_SAFE_MALLOC (c->cliquecount - 1, CCtsp_lpclique *);
    if (teeth == (CCtsp_lpclique **) NULL) {
        fprintf (stderr, "out of memory in CCtsp_star_to_double_decker\n");
        rval = 1; goto CLEANUP;
    }

    nteeth = 0;
    for (i = 0; i < c->cliquecount; i++) {
        if (i != ihandle) {
            teeth[nteeth++] = &c->cliques[i];
        }
    }

    rval = comb_handling (g, h, handle, nteeth, teeth, x, d);
    if (rval) {
        fprintf (stderr, "comb_handling failed\n"); goto CLEANUP;
    }

CLEANUP:

    CC_IFFREE (teeth, CCtsp_lpclique *);
    return rval;
}

#define MAX_ROUNDS (20)   /* Number of rounds of handle growing */

static int comb_handling (CCtsp_lpgraph *g, CC_GCgraph *h,
        CCtsp_lpclique *handle, int nteeth, CCtsp_lpclique **teeth, double *x,
        CCtsp_lpcut_in **cuts)
{
    int rval = 0;
    int i, round, test, extranodes;
    int *marks   = (int *) NULL;
    int *targets = (int *) NULL;
    int *alpha   = (int *) NULL;
    int *beta    = (int *) NULL;
    CCtsp_lpclique  *newcliques   = (CCtsp_lpclique *) NULL;
    CCtsp_lpclique  *oldcliques   = (CCtsp_lpclique *) NULL;
    CCtsp_lpclique **workingteeth = (CCtsp_lpclique **) NULL;
    CCtsp_lpclique **handlearray  = (CCtsp_lpclique **) NULL;
    CCtsp_lpclique *newhandle;
    CCtsp_lpcut_in *dp;

    marks  = CC_SAFE_MALLOC (g->ncount, int);
    if (!marks) {
        fprintf (stderr, "out of memory in comb_handling\n");
        rval = 1; goto CLEANUP;
    }
    CCtsp_mark_clique_and_neighbors (g, handle, marks, 0);
    for (i = 0; i < nteeth; i++) {
        CCtsp_mark_clique_and_neighbors (g, teeth[i], marks, 0);
    }

    targets      = CC_SAFE_MALLOC (nteeth, int);
    alpha        = CC_SAFE_MALLOC (MAX_ROUNDS, int);
    beta         = CC_SAFE_MALLOC (nteeth, int);
    workingteeth = CC_SAFE_MALLOC (nteeth, CCtsp_lpclique *);
    if (!targets || !alpha || !beta || !workingteeth) {
        fprintf (stderr, "out of memory in comb_handling\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < nteeth; i++) {
        workingteeth[i] = teeth[i];
        beta[i] = 1;
    }
    for (i = 0; i < MAX_ROUNDS; i++) {
        alpha[i] = 1;
    }


    handlearray  = CC_SAFE_MALLOC (MAX_ROUNDS + 1, CCtsp_lpclique *);
    if (!handlearray) {
        fprintf (stderr, "out of memory in comb_handling\n");
        rval = 1; goto CLEANUP;
    }
    handlearray[0] = handle;
    for (i = 0; i < MAX_ROUNDS; i++) {
        handlearray[i+1] = (CCtsp_lpclique *) NULL;
    }

    for (round = 0; round < MAX_ROUNDS; round++) {
        CCtsp_mark_clique (handlearray[round], marks, 1);
        for (i = 0; i < nteeth; i++) {
            CCtsp_clique_marked_count (workingteeth[i], marks, 0, &extranodes);
            if (extranodes < 2) targets[i] = 1;
            else                targets[i] = 0;
        }
        CCtsp_mark_clique (handlearray[round], marks, 0);
        rval = stretch_teeth (g, h, handlearray[round], nteeth, workingteeth,
                              &newcliques, 1, targets, &test);

        if (round > 0) {
            for (i = 0; i < nteeth; i++) {
                CCtsp_free_lpclique (workingteeth[i]);
            }
            CC_IFFREE (oldcliques, CCtsp_lpclique);
        }

        if (rval) {
            fprintf (stderr, "stretch_teeth failed\n");
            goto CLEANUP;
        }
        if (test == 0) {
            /* no room to grow the teeth */
            goto CLEANUP;
        }

        oldcliques = newcliques;
        for (i = 0; i < nteeth; i++) {
            workingteeth[i] = &newcliques[i];
            CCtsp_mark_clique_and_neighbors (g, workingteeth[i], marks, 0);
        }

        rval = grow_handle (g, x, handlearray[round], nteeth, workingteeth,
                            &newhandle, marks);
        if (rval) {
            fprintf (stderr, "grow_handle failed\n");
            for (i = 0; i < nteeth; i++) {
                CCtsp_free_lpclique (workingteeth[i]);
            }
            CC_IFFREE (oldcliques, CCtsp_lpclique);
            goto CLEANUP;
        }
        handlearray[round+1] = newhandle;

        rval = build_star (&dp, round+1, handlearray, alpha, nteeth, 
                           workingteeth, beta, g->ncount);
        if (rval) {
            fprintf (stderr, "build_star failed\n"); goto CLEANUP;
        }

        if (dp) {
            rval = CCverify_cut (dp, CC_TYPE_STAR, (int *) NULL);
            if (rval) {
                printf ("Bad star\n"); fflush (stdout);
                printf ("Number Handles: %d   Number Teeth: %d\n",
                           round+1, nteeth); fflush (stdout);
                CCtsp_print_lpcut_in (dp);
                rval = 0; goto CLEANUP;
            }
            dp->next = *cuts;
            *cuts = dp;
        }
    }

    if (round > 0) {
        for (i = 0; i < nteeth; i++) {
            CCtsp_free_lpclique (workingteeth[i]);
        }
        CC_IFFREE (oldcliques, CCtsp_lpclique);
    }

CLEANUP:

    if (handlearray) {
        for (i = 0; i < MAX_ROUNDS; i++) {
            if (handlearray[i+1] != (CCtsp_lpclique *) NULL) {
                CCtsp_free_lpclique (handlearray[i+1]);
                CC_FREE (handlearray[i+1], CCtsp_lpclique);
            }
        }
        CC_FREE (handlearray, CCtsp_lpclique *);
    }

    CC_IFFREE (marks, int);
    CC_IFFREE (targets, int);
    CC_IFFREE (alpha, int);
    CC_IFFREE (beta, int);
    CC_IFFREE (workingteeth, CCtsp_lpclique *);
    CC_IFFREE (handlearray, CCtsp_lpclique *);

    return rval;
}

static int grow_handle (CCtsp_lpgraph *g, double *x, CCtsp_lpclique *handle,
        int nteeth, CCtsp_lpclique **teeth, CCtsp_lpclique **newhandle,
        int *marks)
{
    int rval = 0;
    int i, tnode;
    double tval;
    int *addnodes = (int *) NULL;
    double *values = (double *) NULL;
    CCtsp_lpclique *clique = (CCtsp_lpclique *) NULL;
   
    *newhandle = (CCtsp_lpclique *) NULL;

    values   = CC_SAFE_MALLOC (g->ncount, double);
    addnodes = CC_SAFE_MALLOC (nteeth, int);
    if (!values || !addnodes) {
        fprintf (stderr, "out of memory in grow_handle\n");
        rval = 1; goto CLEANUP;
    }

    CCtsp_mark_clique_and_neighbors_double (g, handle, values, 0.0);
    for (i = 0; i < nteeth; i++) {
        CCtsp_mark_clique_and_neighbors_double (g, teeth[i], values, 0.0);
    }
    add_x_to_clique_neighbors (g, x, handle, values);


    CCtsp_mark_clique (handle, marks, 1);
    for (i = 0; i < nteeth; i++) {
        find_largest_marked (teeth[i], values, marks, 0, &tnode, &tval);
        addnodes[i] = tnode;
    }
    CCtsp_mark_clique (handle, marks, 0);

    clique = CC_SAFE_MALLOC (1, CCtsp_lpclique);
    if (!clique) {
        fprintf (stderr, "out of memory in grow_handle\n");
        rval = 1; goto CLEANUP;
    }


    rval = CCtsp_add_nodes_to_lpclique (handle, clique, nteeth, addnodes);
    if (rval) {
        fprintf (stderr, "CCtsp_add_nodes_to_lpclique failed\n");
        CC_IFFREE (clique, CCtsp_lpclique);
        goto CLEANUP;
    }
    *newhandle = clique;

CLEANUP:

    CC_IFFREE (addnodes, int);
    CC_IFFREE (values, double);
    return rval;
}
