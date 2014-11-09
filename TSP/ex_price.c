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
/*            ROUTINES TO PRICE EDGES USING BIGGUY ARITHMETIC               */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: January 21, 1997                                                  */
/*                                                                          */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCtsp_exact_price (CCtsp_lp *lp, CCbigguy *bound,                   */
/*      int complete_price, int phase1, int silent)                         */
/*    RETURNS a bound that is valid for the entire edge set; the values     */
/*         of the dual variables will be stored in lp->exact_dual unless    */
/*         the existing lp->exact_dual's cutcount agrees with the           */
/*         cutcount for lp                                                  */
/*      -lp is a pointer to the tsp lp                                      */
/*      -bound returns the LP bound                                         */
/*      -if complete_price is nonzero, then price over the complete         */
/*       graph, even if a valid full edge set is present                    */
/*      -phase1 is either 0 or 1, with 1 indicating that the pricing        */
/*       should be to determine a Farkas' lemma bound to prove that the     */
/*       LP is infeasbile                                                   */
/*                                                                          */
/*  int CCtsp_edge_elimination (CCtsp_lp *lp, int eliminate_sparse,         */
/*      int  silent)                                                        */
/*    USES the bound information to elimination edges and set edges to 1;   */
/*         the remaining edges are placed into lp->fulladj (the old adj     */
/*         is freed) and the fixed edges are placed on the list             */
/*         lp->fixededges; the dual values are taken from lp->exact_dual    */
/*      -lp is a pointer to the tsp lp; lp->exact_lowerbound will be used   */
/*       together with lp->upperbound to determine the cutoff to elim       */
/*       and fix edges                                                      */
/*      -if eliminate_sparse is nonzero and lp->full_edges_valid, then      */
/*       the elimination is based on the lp->fulladj list.  Otherwise       */
/*       the elimination is based on the lp->dat.                           */
/*      NOTES: this function does not alter the LP or lp->graph             */
/*                                                                          */
/*  int CCtsp_exact_dual (CCtsp_lp *lp)                                     */
/*    RETURNS the dual values as bigguys (used to store the info used       */
/*        to establish the exact lower bound); the values will be           */
/*        stored in lp->exact_dual (and the existing values freed)          */
/*     -lp is the CCtsp_lp                                                  */
/*                                                                          */
/*  int CCtsp_verify_infeasible_lp (CCtsp_lp *lp, int *yesno, int silent)   */
/*    VERIFIES that the lp is infeasible using exact pricing.               */
/*     -yesno is set to 1 if the lp is infeasible and 0 otherwise.          */
/*                                                                          */
/*  int CCtsp_verify_lp_prune (CCtsp_lp *lp, int *yesno, int silent)        */
/*    VERIFIES that the lp bound is less than the lp upperbound - 1.        */
/*     -yesno is set to 1 if bound < upperbound - 1 and 0 otherwise.        */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "macrorus.h"
#include "bigguy.h"
#include "tsp.h"

#define BIG_PRICE_GEN 20000

typedef struct bigpredge {
    int ends[2];
    int len;
    CCbigguy rc;
} bigpredge;


static int
    big_pricing_duals (CCtsp_lp *lp, CCbigguy *node_pi, CCbigguy *node_piest,
        CCbigguy *node_domino, CCbigguy *cut_pi, CCbigguy *clique_pi,
        CCbigguy *rhs_sum),
    big_price_list (CCtsp_lp *lp, int ecount, bigpredge *elist,
        CCbigguy *node_pi, CCbigguy *clique_pi, CCbigguy *cut_pi),
    big_generate_edges (CCtsp_lp *lp, int use_full_edges, CCbigguy *node_piest,
        CCbigguy *node_domino, int nwant, int *gencount, bigpredge *genlist,
        int *n1, int *n2, int *finished, CCbigguy cutoff, int phase1),
    add_to_inlist (CCtsp_lp *lp, int use_full_edges, bigpredge *inlist,
        int *count, int end0, int end1, int phase1),
    add_to_adj (CCtsp_lp *lp, CCtsp_genadj *adj, int end0, int end1,
        int *count),
    test_edge (int end1, int end2, int len, CCbigguy *node_pi,
        CCbigguy *node_domino, CCbigguy cutoff);


int CCtsp_exact_price (CCtsp_lp *lp, CCbigguy *bound, int complete_price,
        int phase1, int silent)
{
    int incount;
    bigpredge *inlist = (bigpredge *) NULL;
    CCbigguy penalty, rhs_sum;
    CCbigguy *node_pi = (CCbigguy *) NULL;
    CCbigguy *node_piest = (CCbigguy *) NULL;
    CCbigguy *node_domino = (CCbigguy *) NULL;
    CCbigguy *clique_pi = (CCbigguy *) NULL;
    CCbigguy *cut_pi = (CCbigguy *) NULL;
    int nbranch = 0;
    int n1, n2, i, j, finished;
    int rval = 0;
    double szeit;
    int use_full_edges;

    szeit = CCutil_zeit ();
    *bound = CCbigguy_ZERO;

    if (!lp->dat && !lp->full_edges_valid) {
        fprintf (stderr, "must have dat file or full edge set\n");
        return 1;
    }

    if (!lp->dat && complete_price) {
        fprintf (stderr, "must have dat file for complete price\n");
        return 1;
    }

    if (phase1) {
        printf ("phase 1 pricing\n");
        fflush (stdout);
    }

    use_full_edges = lp->full_edges_valid;
    
    if (complete_price) {
        printf ("Pricing COMPLETE GRAPH\n");
        fflush (stdout);
        use_full_edges = 0;
    }

    if (!lp->exact_dual || lp->exact_dual->cutcount != lp->cuts.cutcount) {
        fflush (stdout);
        rval = CCtsp_exact_dual (lp);
        if (rval) {
            fprintf (stderr, "CCtsp_exact_dual failed\n");
            goto CLEANUP;
        }
    }

    for (i = 0; i < lp->branchdepth; i++) {
        if (lp->branchhistory[i].ends[0] != -1) {
            nbranch++;
        }
    }

    incount = 0;
    if (BIG_PRICE_GEN >= lp->nfixededges + nbranch) {
        inlist = CC_SAFE_MALLOC (BIG_PRICE_GEN, bigpredge);
    } else {
        inlist = CC_SAFE_MALLOC (lp->nfixededges + nbranch, bigpredge);
    }
    CCcheck_NULL (inlist, "out of memory in CCtsp_exact_price");

    node_pi    = CC_SAFE_MALLOC (lp->graph.ncount, CCbigguy);
    CCcheck_NULL (node_pi, "out of memory in CCtsp_exact_price");
    node_piest = CC_SAFE_MALLOC (lp->graph.ncount, CCbigguy);
    CCcheck_NULL (node_piest, "out of memory in CCtsp_exact_price");
  
    if (lp->cuts.dominoend) {
        node_domino = CC_SAFE_MALLOC (lp->graph.ncount, CCbigguy);
        CCcheck_NULL (node_domino, "out of memory in CCtsp_exact_price");
    }

    if (lp->cuts.cliqueend) {
        clique_pi = CC_SAFE_MALLOC (lp->cuts.cliqueend, CCbigguy);
        CCcheck_NULL (clique_pi, "out of memory in CCtsp_exact_price");
    }
    if (lp->cuts.cutcount) {
        cut_pi = CC_SAFE_MALLOC (lp->cuts.cutcount, CCbigguy);
        CCcheck_NULL (cut_pi, "out of memory in CCtsp_exact_price");
    }

    rval = big_pricing_duals (lp, node_pi, node_piest, node_domino, cut_pi,
                              clique_pi, &rhs_sum);
    CCcheck_rval (rval, "big_pricing_duals failed");

    finished = 0;
    n1 = 0;
    n2 = (use_full_edges ? 0 : 1);
    penalty = CCbigguy_ZERO;

    while (!finished) {
        rval = big_generate_edges (lp, use_full_edges, node_piest,
                  node_domino, BIG_PRICE_GEN, &incount, inlist, &n1, &n2,
                  &finished, CCbigguy_ZERO, phase1);
        CCcheck_rval (rval, "big_generate_edges failed");

        rval = big_price_list (lp, incount, inlist, node_pi, clique_pi, cut_pi);
        CCcheck_rval (rval, "big_price_list failed");

        for (i = 0; i < incount; i++) {
            if (CCbigguy_cmp (inlist[i].rc, CCbigguy_ZERO) < 0) {
                CCbigguy_add (&penalty, inlist[i].rc);
#if 0
                {
                    int q;
                    q = CCtsp_find_edge (&lp->graph, inlist[i].ends[0],
                                   inlist[i].ends[1]);
                    if (q == -1) {
                        printf ("YIPES: %f [%d, %d %d]: %f %f\n",
                          CCbigguy_bigguytod (inlist[i].rc), inlist[i].ends[0],
                          inlist[i].ends[1], inlist[i].len,
                          CCbigguy_bigguytod (node_pi[inlist[i].ends[0]]),
                          CCbigguy_bigguytod (node_pi[inlist[i].ends[1]]));
                        fflush (stdout);
                    }
                }
#endif
            }
        }
    }

    if (lp->nfixededges + nbranch) {
        int end0, end1;

        /* Adjust bound for fixed/branch edges                            */
        /* If a fixed or a branch=1 edge has positive rc, it can be added */
        /* If a branch=0 edge has negative rc, it should be subtracted    */
        /* from the penalty since it was earlier added.                   */

        incount = 0;
        for (i = 0; i < lp->nfixededges; i++) {
            end0 = lp->fixededges[2*i];
            end1 = lp->fixededges[2*i+1];
            rval = add_to_inlist (lp, use_full_edges, inlist, &incount,
                                  end0, end1, phase1);
            if (rval) {
                fprintf (stderr, "add_to_inlist failed\n");
                goto CLEANUP;
            }
        }
        for (i = 0; i < lp->branchdepth; i++) {
            if (lp->branchhistory[i].ends[0] != -1) {
                end0 = lp->branchhistory[i].ends[0];
                end1 = lp->branchhistory[i].ends[1];
                rval = add_to_inlist (lp, use_full_edges, inlist, &incount,
                                      end0, end1, phase1);
                if (rval) {
                    fprintf (stderr, "add_to_inlist failed\n");
                    goto CLEANUP;
                }
            }
        }

        rval = big_price_list (lp, incount, inlist, node_pi, clique_pi, cut_pi);
        CCcheck_rval (rval, "big_price_list failed");

        for (i = 0; i < lp->nfixededges; i++) {
            if (CCbigguy_cmp (inlist[i].rc, CCbigguy_ZERO) > 0) {
                CCbigguy_add (&penalty, inlist[i].rc);
            }
        }
        for (i = lp->nfixededges, j = 0; i < lp->nfixededges+nbranch; i++) {
            while (lp->branchhistory[j].ends[0] == -1)
                j++;
            if (lp->branchhistory[j].rhs == 0) {
                if (CCbigguy_cmp (inlist[i].rc, CCbigguy_ZERO) < 0) {
                    CCbigguy_sub (&penalty, inlist[i].rc);
                }
            } else {
                if (CCbigguy_cmp (inlist[i].rc, CCbigguy_ZERO) > 0) {
                    
                    CCbigguy_add (&penalty, inlist[i].rc);
                }
            }
            j++;
        }
    }

    *bound = rhs_sum;
    CCbigguy_add (bound, penalty);

    if (!silent) {
        printf ("Exact Price Time: %.2f seconds\n", CCutil_zeit () - szeit);
        fflush (stdout);
    }
    rval = 0;

CLEANUP:

    CC_IFFREE (cut_pi, CCbigguy);
    CC_IFFREE (clique_pi, CCbigguy);
    CC_IFFREE (node_pi, CCbigguy);
    CC_IFFREE (node_piest, CCbigguy);
    CC_IFFREE (node_domino, CCbigguy);
    CC_IFFREE (inlist, bigpredge);
    return rval;
}

static int add_to_inlist (CCtsp_lp *lp, int use_full_edges, bigpredge *inlist,
        int *count, int end0, int end1, int phase1)
{
    int rval = 0;
    int j;
    int len = 0;
    CCtsp_genadj *adj = lp->fulladj;

    if (end0 > end1) {
        CC_SWAP (end0, end1, j);
    }
    if (phase1) {
        len = 0;
    } else {
        if (use_full_edges) {
            for (j = 0; j < adj[end0].deg; j++) {
                if (adj[end0].list[j].end == end1) {
                    len = adj[end0].list[j].len;
                    break;
                }
            }
            if (j == adj[end0].deg) {
                fprintf (stderr, "ERROR: fixed edge not in fulladj\n");
                rval = 1; goto CLEANUP;
            }
        } else {
            len = CCutil_dat_edgelen (end0, end1, lp->dat);
        }
    }
    inlist[*count].ends[0] = end0;
    inlist[*count].ends[1] = end1;
    inlist[*count].len     = len;
    (*count)++;

CLEANUP:

    return rval;
}

int CCtsp_edge_elimination (CCtsp_lp *lp, int eliminate_sparse, int silent)
{
    int incount;
    bigpredge *inlist = (bigpredge *) NULL;
    CCbigguy rhs_sum;
    CCbigguy *node_pi = (CCbigguy *) NULL;
    CCbigguy *node_piest = (CCbigguy *) NULL;
    CCbigguy *node_domino = (CCbigguy *) NULL;
    CCbigguy *clique_pi = (CCbigguy *) NULL;
    CCbigguy *cut_pi = (CCbigguy *) NULL;
    CCtsp_genadj *adj = (CCtsp_genadj *) NULL;
    CCtsp_genadjobj *adjspace = (CCtsp_genadjobj *) NULL;
    CCtsp_genadjobj *pa;
    int i, n1, n2, finished, nremain, nfixed, ek;
    int nbranch = 0;
    int end0, end1;
    int oldnfixed = lp->nfixededges;
    int rval = 0;
    CCbigguy cutoff, negcutoff;
    double szeit;

    /* This code has not been tested after branching (at present, we only */
    /* call edge_elimination at the root LP).                             */

    if (CCbigguy_cmp (lp->exact_lowerbound, CCbigguy_MINBIGGUY) == 0) {
        fprintf (stderr, "need an exact lowerbound to run elimination\n");
        return 1;
    }
    if (lp->upperbound == CCtsp_LP_MAXDOUBLE) {
        fprintf (stderr, "need an exact upperbound to run elimination\n");
        return 1;
    }
    if (!lp->exact_dual || lp->exact_dual->cutcount != lp->cuts.cutcount) {
        fprintf (stderr, "no exact_dual in CCtsp_edge_elimination\n");
        return 1;
    }

    szeit = CCutil_zeit ();

    cutoff = CCbigguy_dtobigguy (lp->upperbound);
    CCbigguy_sub (&cutoff, lp->exact_lowerbound);
    CCbigguy_sub (&cutoff, CCbigguy_ONE);
    negcutoff = CCbigguy_ZERO;
    CCbigguy_sub (&negcutoff, cutoff);
    if (!silent) {
        printf ("Edge Elimination Cutoff: %f\n", CCbigguy_bigguytod (cutoff));
        fflush (stdout);
    }
    if (CCbigguy_cmp (cutoff, CCbigguy_ZERO) < 0) {
        printf ("Cutoff is less than ZERO, do not eliminate\n");
        fflush (stdout);
        return 1;
    }

    incount = 0;
    inlist     = CC_SAFE_MALLOC (BIG_PRICE_GEN, bigpredge);
    CCcheck_NULL (inlist, "out of memory in CCtsp_edge_elimination");

    node_pi    = CC_SAFE_MALLOC (lp->graph.ncount, CCbigguy);
    CCcheck_NULL (node_pi, "out of memory in CCtsp_edge_elimination");
    node_piest = CC_SAFE_MALLOC (lp->graph.ncount, CCbigguy);
    CCcheck_NULL (node_piest, "out of memory in CCtsp_edge_elimination");

    if (lp->cuts.dominoend) {
        node_domino = CC_SAFE_MALLOC (lp->graph.ncount, CCbigguy);
        CCcheck_NULL (node_domino,  "out of memory in CCtsp_edge_elimination");
    }

    if (lp->cuts.cliqueend) {
        clique_pi = CC_SAFE_MALLOC (lp->cuts.cliqueend, CCbigguy);
        CCcheck_NULL (clique_pi,  "out of memory in CCtsp_edge_elimination");
    }
    if (lp->cuts.cutcount) {
        cut_pi = CC_SAFE_MALLOC (lp->cuts.cutcount, CCbigguy);
        CCcheck_NULL (cut_pi,  "out of memory in CCtsp_edge_elimination");
    }

    rval = big_pricing_duals (lp, node_pi, node_piest, node_domino, cut_pi,
                              clique_pi, &rhs_sum);
    CCcheck_rval (rval, "big_pricing_duals failed");

    adj = CC_SAFE_MALLOC (lp->graph.ncount, CCtsp_genadj);
    CCcheck_NULL (adj,  "out of memory in CCtsp_edge_elimination");

    for (i = 0; i < lp->graph.ncount; i++) {
        adj[i].deg = 0;
    }

    if (eliminate_sparse == 0) lp->full_edges_valid = 0;
    finished = 0;
    nremain = 0;
    nfixed = 0;
    n1 = 0;
    n2 = (lp->full_edges_valid ? 0 : 1);

    while (!finished) {
        rval = big_generate_edges (lp, lp->full_edges_valid, node_piest,
                   node_domino, BIG_PRICE_GEN, &incount, inlist, &n1, &n2,
                   &finished, cutoff, 0);
        if (rval) {
            fprintf (stderr, "big_generate_edges failed\n");
            CC_FREE (adj, CCtsp_genadj);
            goto CLEANUP;
        }
        rval = big_price_list (lp, incount, inlist, node_pi, clique_pi, cut_pi);
        if (rval) {
            fprintf (stderr, "big_price_list failed\n");
            CC_FREE (adj, CCtsp_genadj);
            goto CLEANUP;
        }

        for (i = 0; i < incount; i++) {
            if (CCbigguy_cmp (inlist[i].rc, cutoff) <= 0) {
                adj[inlist[i].ends[0]].deg++;
                nremain++;
                if (CCbigguy_cmp (inlist[i].rc, negcutoff) < 0) {
                    ek = CCtsp_find_edge (&(lp->graph),inlist[i].ends[0],
                                                 inlist[i].ends[1]);
                    if (ek != -1 && (lp->graph.edges[ek].fixed  == 0 &&
                                     lp->graph.edges[ek].branch == 0)) {
/* Bico 9.12.99
                    if (ek == -1 || (lp->graph.edges[ek].fixed  == 0 &&
                                     lp->graph.edges[ek].branch == 0)) {
*/
                        /*
                        printf ("[%d, %d] ", inlist[i].ends[0],
                                             inlist[i].ends[1]);
                        fflush (stdout);
                        */
                        nfixed++;
                    }
                }
            }
        }
    }

    /* leave enough room for the old fixed edges and the branch edges */

    for (i = 0; i < lp->branchdepth; i++) {
        if (lp->branchhistory[i].ends[0] != -1) {
            end0 = lp->branchhistory[i].ends[0];
            end1 = lp->branchhistory[i].ends[1];
        if (end0 < end1)
            adj[end0].deg++;
        else
            adj[end1].deg++;
            nbranch++;
        }
    }
    for (i = 0; i < lp->nfixededges; i++) {
        end0 = lp->fixededges[2*i];
        end1 = lp->fixededges[2*i+1];
        if (end0 < end1)
            adj[end0].deg++;
        else
            adj[end1].deg++;
    }

    if (nremain + lp->nfixededges + nbranch) {
        adjspace = CC_SAFE_MALLOC (nremain + lp->nfixededges + nbranch,
                                CCtsp_genadjobj);
        if (!adjspace) {
            fprintf (stderr, "out of memory in CCtsp_edge_elimination\n");
            CC_FREE (adj, CCtsp_genadj);
            rval = 1;
            goto CLEANUP;
        }
    }
    if (nfixed) {
        rval = CCutil_reallocrus_count ((void **) &(lp->fixededges),
                  2 * (lp->nfixededges + nfixed), sizeof (int));
        if (rval) {
            fprintf (stderr, "out of memory in CCtsp_edge_elimination\n");
            CC_FREE (adj, CCtsp_genadj);
            CC_IFFREE (adjspace, CCtsp_genadjobj);
            goto CLEANUP;
        }
    }

    pa = adjspace;
    for (i = 0; i < lp->graph.ncount; i++) {
        adj[i].list = pa;
        pa += adj[i].deg;
        adj[i].deg = 0;
    }

    finished = 0;
    n1 = 0;
    n2 = (lp->full_edges_valid ? 0 : 1);

    while (!finished) {
        rval = big_generate_edges (lp, lp->full_edges_valid, node_piest,
                   node_domino, BIG_PRICE_GEN, &incount, inlist, &n1, &n2,
                   &finished, cutoff, 0);
        if (rval) {
            fprintf (stderr, "big_generate_edges failed\n");
            CC_FREE (adj, CCtsp_genadj);
            CC_IFFREE (adjspace, CCtsp_genadjobj);
            goto CLEANUP;
        }
        rval = big_price_list (lp, incount, inlist, node_pi, clique_pi, cut_pi);
        if (rval) {
            fprintf (stderr, "big_price_list failed\n");
            CC_FREE (adj, CCtsp_genadj);
            CC_IFFREE (adjspace, CCtsp_genadjobj);
            goto CLEANUP;
        }
        for (i = 0; i < incount; i++) {
            if (CCbigguy_cmp (inlist[i].rc, cutoff) <= 0) {
                adj[inlist[i].ends[0]].list[adj[inlist[i].ends[0]].deg].end =
                                              inlist[i].ends[1];
                adj[inlist[i].ends[0]].list[adj[inlist[i].ends[0]].deg].len =
                                              inlist[i].len;
                adj[inlist[i].ends[0]].deg++;
                if (CCbigguy_cmp (inlist[i].rc, negcutoff) < 0) {
                    ek = CCtsp_find_edge (&(lp->graph),inlist[i].ends[0],
                                                 inlist[i].ends[1]);
                    if (ek != -1 && (lp->graph.edges[ek].fixed  == 0 &&
                                     lp->graph.edges[ek].branch == 0)) {
/* Bico 9.12.99
                    if (ek == -1 || (lp->graph.edges[ek].fixed  == 0 &&
                                     lp->graph.edges[ek].branch == 0)) {
*/
                        lp->fixededges[2*lp->nfixededges] = inlist[i].ends[0];
                        lp->fixededges[2*lp->nfixededges+1] =
                                                            inlist[i].ends[1];
                        lp->nfixededges++;
                    }
                }
            }
        }
    }

    /* some of the old fixed and branched edges may not have gotten into adj */

    for (i = 0; i < lp->branchdepth; i++) {
        if (lp->branchhistory[i].ends[0] != -1) {
            end0 = lp->branchhistory[i].ends[0];
            end1 = lp->branchhistory[i].ends[1];
            rval = add_to_adj (lp, adj, end0, end1, &nremain);
            if (rval) {
                fprintf (stderr, "add_to_adj failed\n"); goto CLEANUP;
            }
        }
    }
    for (i = 0; i < oldnfixed; i++) {
        end0 = lp->fixededges[2*i];
        end1 = lp->fixededges[2*i+1];
        rval = add_to_adj (lp, adj, end0, end1, &nremain);
        if (rval) {
            fprintf (stderr, "add_to_adj failed\n"); goto CLEANUP;
        }
    }
    rval = 0;

    CC_IFFREE (lp->fulladjspace, CCtsp_genadjobj);
    CC_IFFREE (lp->fulladj, CCtsp_genadj);
    lp->fullcount = nremain;
    lp->fulladjspace = adjspace;
    lp->fulladj = adj;
    lp->full_edges_valid = 1;

    if (!silent) {
        printf ("Remaining Edges: %d (with %d new fixed)\n", nremain, nfixed);
        printf ("Edge Elimination Time: %.2f seconds\n",
                CCutil_zeit () - szeit);
        fflush (stdout);
    }

CLEANUP:

    CC_IFFREE (cut_pi, CCbigguy);
    CC_IFFREE (clique_pi, CCbigguy);
    CC_IFFREE (node_pi, CCbigguy);
    CC_IFFREE (node_piest, CCbigguy);
    CC_IFFREE (node_domino, CCbigguy);
    CC_IFFREE (inlist, bigpredge);
    return rval;
}

static int add_to_adj (CCtsp_lp *lp, CCtsp_genadj *adj, int end0, int end1,
        int *count)
{
    int rval = 0;
    int len  = 0;
    int j, k;

    if (end0 > end1) {
        CC_SWAP (end0, end1, k);
    }
    for (k = 0; k < adj[end0].deg; k++) {
        if (adj[end0].list[k].end == end1)
            break;
    }
    if (k == adj[end0].deg) {
        if (lp->full_edges_valid) {
            for (j = 0; j < lp->fulladj[end0].deg; j++) {
                if (lp->fulladj[end0].list[j].end == end1) {
                    len = lp->fulladj[end0].list[j].len;
                    break;
                }
            }
            if (j == lp->fulladj[end0].deg) {
                fprintf (stderr, "ERROR: fixed/branch edge not in fulladj\n");
                rval = 1; goto CLEANUP;
            }
        } else {
            len = CCutil_dat_edgelen (end0, end1, lp->dat);
        }
        adj[end0].list[adj[end0].deg].end = end1;
        adj[end0].list[adj[end0].deg].len = len;
        adj[end0].deg++;
        (*count)++;
    }

CLEANUP:

    return rval;
}

static int big_pricing_duals (CCtsp_lp *lp, CCbigguy *node_pi,
        CCbigguy *node_piest, CCbigguy *node_domino, CCbigguy *cut_pi,
        CCbigguy *clique_pi, CCbigguy *rhs_sum)
{
    CCbigguy x;
    int i, j, k, s, tmp;
    int rval = 0;
    CCtsp_lpcut *c;

    *rhs_sum = CCbigguy_ZERO;

    if (!lp->exact_dual || lp->exact_dual->cutcount != lp->cuts.cutcount) {
        fprintf (stderr, "no exact_dual in big_pricing_duals\n");
        rval = 1; goto CLEANUP;
    }

    for (i = 0; i < lp->graph.ncount; i++) {
        node_pi[i] = lp->exact_dual->node_pi[i];
    }
    for (i = 0; i < lp->cuts.cutcount; i++) {
        cut_pi[i] = lp->exact_dual->cut_pi[i];
    }

    for (i = 0; i < lp->graph.ncount; i++)
        CCbigguy_addmult (rhs_sum, node_pi[i], 2);

    for (i = 0; i < lp->cuts.cutcount; i++) {
        x = cut_pi[i];
        CCbigguy_addmult (rhs_sum, x, lp->cuts.cuts[i].rhs);
        for (j = 0; j < lp->cuts.cuts[i].modcount; j++) {
            CCbigguy_addmult (rhs_sum, x,
                          2 * (lp->cuts.cuts[i].mods[j].mult - 128));
        }
    }

    for (i = 0; i < lp->cuts.cliqueend; i++) {
        clique_pi[i] = CCbigguy_ZERO;
    }

    for (i = 0; i < lp->cuts.cutcount; i++) {
        x = cut_pi[i];
        for (j = 0; j < lp->cuts.cuts[i].modcount; j++) {
            CCbigguy_addmult (&(node_pi[lp->cuts.cuts[i].mods[j].node]), x,
                    lp->cuts.cuts[i].mods[j].mult - 128);
        }
        if (lp->cuts.cuts[i].dominocount <= 0) {
            for (j = 0; j < lp->cuts.cuts[i].cliquecount; j++) {
                CCbigguy_add (&(clique_pi[lp->cuts.cuts[i].cliques[j]]), x);
            }
        }
    }

    for (i = 0; i < lp->graph.ncount; i++) {
        node_piest[i] = node_pi[i];
    }

    for (i = 0; i < lp->cuts.cliqueend; i++) {
        x = clique_pi[i];
        if (CCbigguy_cmp (x, CCbigguy_ZERO) > 0) {
            CC_FOREACH_NODE_IN_CLIQUE (j, lp->cuts.cliques[i], tmp) {
                CCbigguy_add (&(node_pi[j]), x);
                CCbigguy_add (&(node_piest[j]), x);
            }
        } else if (CCbigguy_cmp (x, CCbigguy_ZERO) < 0) {
            CC_FOREACH_NODE_IN_CLIQUE (j, lp->cuts.cliques[i], tmp) {
                CCbigguy_add (&(node_pi[j]), x);
            }
        }
    }

    if (node_domino) {
        for (i = 0; i < lp->graph.ncount; i++) {
            node_domino[i] = CCbigguy_ZERO;
        }
        for (i = 0; i < lp->cuts.cutcount; i++) {
            c = &lp->cuts.cuts[i];
            if (c->dominocount > 0) {
                x = cut_pi[i];
                if (CCbigguy_cmp (x, CCbigguy_ZERO) > 0) {
                    if (c->cliquecount != 1) {
                        fprintf (stderr, "YIPES: No handle for domino\n");
                        rval = 1;  goto CLEANUP;
                    }
                    CC_FOREACH_NODE_IN_CLIQUE (j,
                         lp->cuts.cliques[c->cliques[0]], tmp) {
                        CCbigguy_add (&(node_domino[j]), x);
                    }
                    for (k = 0; k < c->dominocount; k++) {
                        for (s = 0; s < 2; s++) {
                            CC_FOREACH_NODE_IN_CLIQUE (j,
                                 lp->cuts.dominos[c->dominos[k]].sets[s], tmp) {
                                CCbigguy_addmult (&(node_domino[j]), x, 2);
                            }
                        }
                    }
                } else if (CCbigguy_cmp (x, CCbigguy_ZERO) < 0) {
                    fprintf (stderr, "YIPES: negative domino\n");
                    rval = 1;  goto CLEANUP;
                }
            }
        }
    }

CLEANUP:

    return rval;
}

static int big_price_list (CCtsp_lp *lp, int ecount, bigpredge *elist,
    CCbigguy *node_pi, CCbigguy *clique_pi, CCbigguy *cut_pi)
{
    CCtsp_lpadj *adjspace = (CCtsp_lpadj *) NULL;
    CCtsp_lpnode *n = (CCtsp_lpnode *) NULL;
    int *temp_elist = (int *) NULL;
    int i, j, tmp, l, nzlist, nznext;
    CCtsp_lpadj *a;
    int marker = 0;
    CCbigguy x;
    int ncount = lp->graph.ncount;
    int ccount = lp->cuts.cliqueend;
    CCtsp_lpclique *c = lp->cuts.cliques;
    CCtsp_lpcut *cut;
    CCtsp_lpgraph g;
    int rval = 0;

    CCtsp_init_lpgraph_struct (&g);

    if (ecount == 0) goto CLEANUP;

    n = CC_SAFE_MALLOC (ncount, CCtsp_lpnode);
    CCcheck_NULL (n, "out of memory in big_price_list");
    adjspace = CC_SAFE_MALLOC (2*ecount, CCtsp_lpadj);
    CCcheck_NULL (adjspace, "out of memory in big_price_list");

    for (i = 0; i < ncount; i++) {
        n[i].deg = 0;
        n[i].mark = 0;
    }
    for (i = 0; i < ecount; i++) {
        elist[i].rc = CCbigguy_itobigguy (elist[i].len);
        CCbigguy_sub (&(elist[i].rc), node_pi[elist[i].ends[0]]);
        CCbigguy_sub (&(elist[i].rc), node_pi[elist[i].ends[1]]);
        n[elist[i].ends[0]].deg++;
        n[elist[i].ends[1]].deg++;
    }
    a = adjspace;
    for (i = 0; i < ncount; i++) {
        n[i].adj = a;
        a += n[i].deg;
        n[i].deg = 0;
    }
    for (i = 0; i < ecount; i++) {
        j = elist[i].ends[0];
        n[j].adj[n[j].deg].to = elist[i].ends[1];
        n[j].adj[n[j].deg].edge = i;
        n[j].deg++;
        j = elist[i].ends[1];
        n[j].adj[n[j].deg].to = elist[i].ends[0];
        n[j].adj[n[j].deg].edge = i;
        n[j].deg++;
    }

    for (i = 0; i < ccount; i++) {
        if (CCbigguy_cmp (clique_pi[i], CCbigguy_ZERO)) {
            x = clique_pi[i];
            CCbigguy_add (&x, clique_pi[i]);
            marker++;
            CC_FOREACH_NODE_IN_CLIQUE (j, c[i], tmp) {
                a = n[j].adj;
                for (l = 0; l < n[j].deg; l++) {
                    if (n[a[l].to].mark == marker) {
                        CCbigguy_add (&(elist[a[l].edge].rc), x);
                    }
                }
                n[j].mark = marker;
            }
        }
    }

    /* Price dominos, using nzlist */

    if (lp->cuts.dominoend > 0) {
        temp_elist = CC_SAFE_MALLOC (2*ecount, int);
        CCcheck_NULL (temp_elist, "out of memory in price_list");
        for (i = 0; i < ecount; i++) {
            temp_elist[2*i]   = elist[i].ends[0];
            temp_elist[2*i+1] = elist[i].ends[1];
        }
        rval = CCtsp_build_lpgraph (&g, ncount, ecount, temp_elist,
                                   (int *) NULL);
        CCcheck_rval (rval, "CCtsp_build_lpgraph failed");
        rval = CCtsp_build_lpadj (&g, 0, ecount);
        CCcheck_rval (rval, "CCtsp_build_lpadj failed");
        CC_FREE (temp_elist, int);

        for (i = 0; i < lp->cuts.cutcount; i++) {
            cut = &lp->cuts.cuts[i];
            if (cut->dominocount > 0) {
                x = cut_pi[i];
                if (CCbigguy_cmp (x, CCbigguy_ZERO) > 0) {
                    if (cut->cliquecount != 1) {
                        fprintf (stderr, "YIPES: No handle for domino\n");
                        rval = 1;  goto CLEANUP;
                    }
                    nzlist = CCtsp_lpcut_nzlist (&g, cut, lp->cuts.cliques,
                                                 lp->cuts.dominos, 0);
                    while (nzlist != -1) {
                        nznext = g.edges[nzlist].coefnext;
                        g.edges[nzlist].coefnext = -2;
                        if (g.edges[nzlist].coef) {
                            CCbigguy_addmult (&(elist[nzlist].rc), x,
                                              -g.edges[nzlist].coef);
                            g.edges[nzlist].coef = 0;
                        }
                        nzlist = nznext;
                    }
                } else if (CCbigguy_cmp (x, CCbigguy_ZERO) < 0) {
                    fprintf (stderr, "YIPES: negative domino\n");
                    rval = 1;  goto CLEANUP;
                }
            }
        }
    }

CLEANUP:

    CC_IFFREE (n, CCtsp_lpnode);
    CC_IFFREE (adjspace, CCtsp_lpadj);
    CC_IFFREE (temp_elist, int);
    CCtsp_free_lpgraph (&g);

    return rval; 
}

static int big_generate_edges (CCtsp_lp *lp, int use_full_edges,
        CCbigguy *node_piest, CCbigguy *node_domino, int nwant, int *gencount,
        bigpredge *genlist, int *n1, int *n2, int *finished, CCbigguy cutoff,
        int phase1)
{
    int i = *n1;
    int j = *n2;
    int cnt = 0;
    int len;
    CCtsp_genadj *adj;
    CCdatagroup *dat;
    int ncount = lp->graph.ncount;

    *gencount = 0;
    *finished = 0;

    if (!lp->dat && !use_full_edges) {
        fprintf (stderr, "no source of edges in big_generate_edges\n");
        return 1;
    }

    if (i >= ncount) {
        *finished = 1;
        return 0;
    }

    if (use_full_edges) {
        if (!phase1) {
            adj = lp->fulladj;
            for (; j < adj[i].deg; j++) {
                if (test_edge (i, adj[i].list[j].end, adj[i].list[j].len,
                               node_piest, node_domino, cutoff)) {
                    genlist[cnt].ends[0] = i;
                    genlist[cnt].ends[1] = adj[i].list[j].end;
                    genlist[cnt].len = adj[i].list[j].len;
                    cnt++;
                    if (cnt == nwant) {
                        *finished = 0;
                        *gencount = cnt;
                        *n1 = i; *n2 = j + 1;
                        return 0;
                    }
                }
            }
            for (i++; i < ncount; i++) {
                for (j = 0; j < adj[i].deg; j++) {
                    if (test_edge (i, adj[i].list[j].end, adj[i].list[j].len,
                                   node_piest, node_domino, cutoff)) {
                        genlist[cnt].ends[0] = i;
                        genlist[cnt].ends[1] = adj[i].list[j].end;
                        genlist[cnt].len = adj[i].list[j].len;
                        cnt++;
                        if (cnt == nwant) {
                            *finished = 0;
                            *gencount = cnt;
                            *n1 = i; *n2 = j + 1;
                            return 0;
                        }
                    }
                }
            }
        } else {
            adj = lp->fulladj;
            for (; j < adj[i].deg; j++) {
                if (test_edge (i, adj[i].list[j].end, 0, node_piest,
                               node_domino, cutoff)) {
                    genlist[cnt].ends[0] = i;
                    genlist[cnt].ends[1] = adj[i].list[j].end;
                    genlist[cnt].len = 0;
                    cnt++;
                    if (cnt == nwant) {
                        *finished = 0;
                        *gencount = cnt;
                        *n1 = i; *n2 = j + 1;
                        return 0;
                    }
                }
            }
            for (i++; i < ncount; i++) {
                for (j = 0; j < adj[i].deg; j++) {
                    if (test_edge (i, adj[i].list[j].end, 0, node_piest,
                                   node_domino, cutoff)) {
                        genlist[cnt].ends[0] = i;
                        genlist[cnt].ends[1] = adj[i].list[j].end;
                        genlist[cnt].len = 0;
                        cnt++;
                        if (cnt == nwant) {
                            *finished = 0;
                            *gencount = cnt;
                            *n1 = i; *n2 = j + 1;
                            return 0;
                        }
                    }
                }
            }
        }
    } else {
        if (!phase1) {
            dat = lp->dat;
            for (; j < ncount; j++) {
                len = CCutil_dat_edgelen (i, j, dat);
                if (test_edge (i, j, len, node_piest, node_domino, cutoff)) {
                    genlist[cnt].ends[0] = i;
                    genlist[cnt].ends[1] = j;
                    genlist[cnt].len = len;
                    cnt++;
                    if (cnt == nwant) {
                        *finished = 0;
                        *gencount = cnt;
                        *n1 = i; *n2 = j + 1;
                        return 0;
                    }
                }
            }
            for (i++; i < ncount; i++) {
                for (j = i + 1; j < ncount; j++) {
                    len = CCutil_dat_edgelen (i, j, dat);
                    if (test_edge (i, j, len, node_piest, node_domino,
                                  cutoff)) {
                        genlist[cnt].ends[0] = i;
                        genlist[cnt].ends[1] = j;
                        genlist[cnt].len = len;
                        cnt++;
                        if (cnt == nwant) {
                            *finished = 0;
                            *gencount = cnt;
                            *n1 = i; *n2 = j + 1;
                            return 0;
                        }
                    }
                }
            }
        } else {
            for (; j < ncount; j++) {
                if (test_edge (i, j, 0, node_piest, node_domino, cutoff)) {
                    genlist[cnt].ends[0] = i;
                    genlist[cnt].ends[1] = j;
                    genlist[cnt].len = 0;
                    cnt++;
                    if (cnt == nwant) {
                        *finished = 0;
                        *gencount = cnt;
                        *n1 = i; *n2 = j + 1;
                        return 0;
                    }
                }
            }
            for (i++; i < ncount; i++) {
                for (j = i + 1; j < ncount; j++) {
                    if (test_edge (i, j, 0, node_piest, node_domino, cutoff)) {
                        genlist[cnt].ends[0] = i;
                        genlist[cnt].ends[1] = j;
                        genlist[cnt].len = 0;
                        cnt++;
                        if (cnt == nwant) {
                            *finished = 0;
                            *gencount = cnt;
                            *n1 = i; *n2 = j + 1;
                            return 0;
                        }
                    }
                }
            }
        }
    }

    *n1 = ncount;
    *n2 = ncount;
    *gencount = cnt;
    *finished = 1;

    return 0;
}

static int test_edge (int end1, int end2, int len, CCbigguy *node_pi,
        CCbigguy *node_domino, CCbigguy cutoff)
{
    CCbigguy rc;

    rc = CCbigguy_itobigguy (len);
    CCbigguy_sub (&rc, node_pi[end1]);
    CCbigguy_sub (&rc, node_pi[end2]);
    if (node_domino) {
        CCbigguy_sub (&rc, node_domino[end1]);
        CCbigguy_sub (&rc, node_domino[end2]);
    }
    if (CCbigguy_cmp (rc, cutoff) <= 0) {
        return 1;
    } else {
        return 0;
    }
}

int CCtsp_exact_dual (CCtsp_lp *lp)
{
    int i;
    int ncount = lp->graph.ncount;
    int cutcount = lp->cuts.cutcount;
    double *d_node_pi = (double *) NULL;
    double *d_cut_pi = (double *) NULL;
    CCtsp_bigdual *d = (CCtsp_bigdual *) NULL;
    int rval = 0;

    rval = CCtsp_get_lp_result (lp, (double *) NULL, (double *) NULL,
                    (int *) NULL, (int **) NULL, (double **) NULL,
                    (double **) NULL, &d_node_pi, &d_cut_pi);
    if (rval) {
        fprintf (stderr, "get_lp_result failed\n");
        fflush (stdout);
        goto CLEANUP;
    }

    d = CC_SAFE_MALLOC (1, CCtsp_bigdual);
    if (!(d)) {
        fprintf (stderr, "out of memory in CCtsp_exact_dual C\n");
        rval = 1; goto CLEANUP;
    }
    d->cutcount = cutcount;
    d->node_pi = (CCbigguy *) NULL;
    d->cut_pi = (CCbigguy *) NULL;

    d->node_pi = CC_SAFE_MALLOC (ncount, CCbigguy);
    if (!d->node_pi) {
        fprintf (stderr, "out of memory in CCtsp_exact_dual B\n");
        CC_FREE (d, CCtsp_bigdual);
        rval = 1; goto CLEANUP;
    }

    for (i = 0; i < ncount; i++) {
        d->node_pi[i] = CCbigguy_dtobigguy (d_node_pi[i]);
    }

    if (cutcount) {
        d->cut_pi = CC_SAFE_MALLOC (cutcount, CCbigguy);
        if (!d->cut_pi) {
            fprintf (stderr, "out of memory in CCtsp_exact_dual A\n");
            CC_FREE (d->node_pi, CCbigguy);
            CC_FREE (d, CCtsp_bigdual);
            rval = 1; goto CLEANUP;
        }
        for (i = 0; i < lp->cuts.cutcount; i++) {
            d->cut_pi[i] = CCbigguy_dtobigguy (d_cut_pi[i]);
        }
    }

    if (lp->exact_dual) {
        CC_IFFREE (lp->exact_dual->node_pi, CCbigguy);
        CC_IFFREE (lp->exact_dual->cut_pi, CCbigguy);
        CC_FREE (lp->exact_dual, CCtsp_bigdual);
    }
    lp->exact_dual = d;


CLEANUP:

    CC_IFFREE (d_node_pi, double);
    CC_IFFREE (d_cut_pi, double);
    return rval;
}

void CCtsp_free_bigdual (CCtsp_bigdual **d)
{
    if (d) {
        if (*d) {
            CC_IFFREE ((*d)->node_pi, CCbigguy);
            CC_IFFREE ((*d)->cut_pi,  CCbigguy);
            CC_IFFREE ((*d), CCtsp_bigdual);
        }
    }
}

int CCtsp_verify_infeasible_lp (CCtsp_lp *lp, int *yesno, int silent)
{
    int rval = 0;
    CCbigguy exactbound;

    *yesno = 0;
    rval = CCtsp_exact_price (lp, &exactbound, 0, 1, silent);
    if (rval) {
        fprintf (stderr, "CCtsp_exact_price_failed\n"); goto CLEANUP;
    }

    if (!silent) {
        printf ("Exactbound: %f\n", CCbigguy_bigguytod (exactbound));
        fflush (stdout);
    }

    if (CCbigguy_cmp (exactbound, CCbigguy_ZERO) > 0) {
        printf ("Problem is shown to be infeasible\n"); fflush (stdout);
        *yesno = 1;
        lp->infeasible = 1;
        lp->lowerbound = CCtsp_LP_MAXDOUBLE;
        lp->exact_lowerbound = CCbigguy_MAXBIGGUY;
    } else {
        printf ("Did not verify an infeasible LP\n"); fflush (stdout);
        *yesno = 0;
    }

CLEANUP:

    return rval;
}

int CCtsp_verify_lp_prune (CCtsp_lp *lp, int *yesno, int silent)
{
    int rval = 0;
    CCbigguy exactbound, bnd;

    *yesno = 0;
    rval = CCtsp_exact_price (lp, &exactbound, 0, 0, silent);
    if (rval) {
        fprintf (stderr, "CCtsp_exact_price_failed\n"); goto CLEANUP;
    }
  
    if (!silent) {
        printf ("Exact LP bound: %f\n", CCbigguy_bigguytod (exactbound));
        fflush (stdout);
    }

    bnd = CCbigguy_dtobigguy (lp->upperbound);
    CCbigguy_sub (&bnd, CCbigguy_ONE);

    if (CCbigguy_cmp (exactbound, bnd) > 0) {
        if (!silent) {
            printf ("Can prune lp.\n"); fflush (stdout);
        }
        *yesno = 1;
        lp->exact_lowerbound = exactbound;
    } else {
        if (!silent) {
            printf ("Cannot prune lp.\n"); fflush (stdout);
        }
        *yesno = 0;
    }

CLEANUP:

    return rval;
}
