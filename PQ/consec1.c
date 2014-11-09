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
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCpq_consecutiveones (CCtsp_lpcut_in **cuts, int *cutcount,         */
/*      CCtsp_cuttree *ctree, CCtsp_lpcuts *pool, int ecount,               */
/*      int *elist, double *x)                                              */
/*    NONE                                                                  */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "tsp.h"
#include "pq.h"
#include "cuttree.h"
#include "consec1.h"

#define MINVIOL 0.001
#undef DEBUG

struct conflict_data_struct {
    CCpq_tree *pqt;
    CCtsp_lpclique *conflict;
};


static int
    found_clique_conflict (CCtsp_lpcut_in **cuts, int *cutcount,
        CCtsp_cuttree *ctree, int viol, CCtsp_lpcuts *pool, int *cliquenum,
        int *perm, int clique_cnt),
    found_ctree_conflict (CCtsp_lpcut_in **cuts, int *cutcount,
        int *conflict_list, int conflict_cnt, CCtsp_cuttree *ctree,
        CCtsp_lpcuts *pool),
    conflict_callback (int *arr, int cnt, int *stop, void *u_data),
    found_comb_conflict (CCtsp_lpcut_in **cuts, int *cutcount, int nodecount,
        CCtsp_lpclique *conflict_cliques[4]),
    histok (unsigned int *hist),
    add_conflict_to_cuts (CCtsp_lpcut_in **cuts, int *cutcount,
        CCtsp_lpclique *conflict_cliques[4], int ncount);


int CCpq_consecutiveones (CCtsp_lpcut_in **cuts, int *cutcount,
        CCtsp_cuttree *ctree, CCtsp_lpcuts *pool, int ecount, int *elist,
        double *x)
{
    int *cliquenums = (int *) NULL;
    double *cliquevals = (double *) NULL;
    int *perm = (int *) NULL;
    CCpq_tree pqt;
    int nodecount = ctree->nodecount;
    int cliquecount;
    CCtsp_lpclique *c;
    int status;
    int bot;
    int top;
    int rval;
    int i;
    int cnt_bot_contrad = 0;
    int cnt_top_contrad = 0;
    int cnt_trivial = 0;
    int cnt_nontrivial = 0;
    double bot3;
/*    double szeit = CCutil_zeit ();*/

    CCpq_tree_init (&pqt);
    *cutcount = 0;

    rval = CCpq_cuttree_to_pq (ctree, &pqt);
    if (rval) {
        fprintf (stderr, "CCpq_cuttree_to_pq failed\n");
        goto CLEANUP;
    }

    rval = CCtsp_get_clique_prices (pool, &cliquenums, &cliquevals,
                                    0.0, 4.0 - MINVIOL, &cliquecount,
                                    nodecount, ecount, elist, x);
    if (rval) {
        fprintf (stderr, "CCtsp_get_clique_prices failed\n");
        goto CLEANUP;
    }

    if (cliquecount == 0) {
        rval = 0; goto CLEANUP;
    }

    perm = CC_SAFE_MALLOC (cliquecount, int);
    if (perm == (int *) NULL) {
        fprintf (stderr, "Out of memory in CCpq_consecutiveones\n");
        rval = 1; goto CLEANUP;
    }
    for (i=0; i<cliquecount; i++) perm[i] = i;

    CCutil_double_perm_quicksort (perm, cliquevals, cliquecount);

    bot = 0;
    top = cliquecount-1;
    bot3 = 0.0;
    /* [0, bot) is base, (top, cliquecount-1] is tested */
    /* bot3 is cliquevals[bot] + [bot-1] + [bot-2] */

    while (bot < top) {
        if (bot3 + cliquevals[perm[top]] < 10.0 - MINVIOL) {
            rval = CCtsp_get_clique (pool, cliquenums[perm[bot]], &c);
            if (rval) {
                fprintf (stderr, "CCtsp_get_clique failed\n");
                goto CLEANUP;
            }
            rval = CCpq_apply_clique (&pqt, c, &status);
            if (rval) {
                fprintf (stderr, "CCpq_apply_clique failed\n");
                goto CLEANUP;
            }
            if (status == CCpq_STATUS_NOSOL) {
                cnt_bot_contrad++;
                rval = found_clique_conflict (cuts, cutcount, ctree,
                                              cliquenums[perm[bot]], pool,
                                              cliquenums, perm, bot);
                if (rval) {
                    fprintf (stderr, "found_clique_conflict failed\n");
                    goto CLEANUP;
                }
                perm[bot] = -1;
            } else if (status == CCpq_STATUS_TRIVIAL) {
                cnt_trivial++;
                perm[bot] = -1;
            } else if (status == CCpq_STATUS_NONTRIVIAL) {
                cnt_nontrivial++;
            } else {
                fprintf (stderr, "Unknown PQ status %d\n", status);
                rval = 1; goto CLEANUP;
            }
            bot++;
            bot3 += cliquevals[bot];
            if (bot >= 3) bot3 -= cliquevals[bot-3];
        } else {
            rval = CCtsp_get_clique (pool, cliquenums[perm[top]], &c);
            if (rval) {
                fprintf (stderr, "CCtsp_get_clique failed\n");
                goto CLEANUP;
            }
            CCpq_check_clique (&pqt, c, &status);
            if (status == CCpq_STATUS_NOSOL) {
                cnt_top_contrad++;
                rval = found_clique_conflict (cuts, cutcount, ctree,
                                              cliquenums[perm[top]], pool,
                                              cliquenums, perm, bot);
                if (rval) {
                    fprintf (stderr, "found_clique_conflict failed\n");
                    goto CLEANUP;
                }
            }
            top--;
        }
    }

/*
    printf ("PQ: %d bottom: %d nontrivial, %d trivial, %d contrad\n",
            bot, cnt_nontrivial, cnt_trivial, cnt_bot_contrad);
    printf ("PQ: %d top: %d contrad\n", cliquecount-1 - top, cnt_top_contrad);
    printf ("CCpq_consecutiveones finished in %.2f seconds\n",
            CCutil_zeit() - szeit);
    fflush (stdout);
*/

    rval = 0;

  CLEANUP:
    CC_IFFREE (cliquenums, int);
    CC_IFFREE (cliquevals, double);
    CC_IFFREE (perm, int);
    CCpq_tree_free (&pqt);
    return rval;
}

/* clique viol added to cliquenums[perm[0..cnt)] results in violation.
   Convert to comb. */

static int found_clique_conflict (CCtsp_lpcut_in **cuts, int *cutcount,
        CCtsp_cuttree *ctree, int viol, CCtsp_lpcuts *pool, int *cliquenums,
        int *perm, int clique_cnt)
{
    int conflict_list[4];
    int conflict_cnt;
    CCpq_tree pqt;
    CCtsp_lpclique *c;
    int i;
    int conflict_found;
    int status;
    int rval;

    CCpq_tree_init (&pqt);

    conflict_list[0] = viol;
    conflict_cnt = 1;

    for (;;) {
        rval = CCpq_cuttree_to_pq (ctree, &pqt);
        if (rval) {
            fprintf (stderr, "CCpq_cuttree_to_pq failed\n");
            goto CLEANUP;
        }
        for (i = 0; i < conflict_cnt; i++) {
            rval = CCtsp_get_clique (pool, conflict_list[i], &c);
            if (rval) {
                fprintf (stderr, "CCtsp_get_clique failed\n");
                goto CLEANUP;
            }

            rval = CCpq_apply_clique (&pqt, c, &status);
            if (rval) {
                fprintf (stderr, "CCpq_apply_clique failed\n");
                goto CLEANUP;
            }

            if (status == CCpq_STATUS_NOSOL) {
                rval = found_ctree_conflict (cuts, cutcount, conflict_list,
                                             conflict_cnt, ctree, pool);
                if (rval) {
                    fprintf (stderr, "found_ctree_conflict failed\n");
                }
                goto CLEANUP;
            }
        }
        if (conflict_cnt == 4) {
#ifdef DEBUG
            printf ("violation uses > 4 + ctree\n");
#endif
            rval = 0;
            goto CLEANUP;
        }

        conflict_found = 0;
        for (i=0; i<clique_cnt; i++) {
            if (perm[i] != -1) {
                rval = CCtsp_get_clique (pool, cliquenums[perm[i]], &c);
                if (rval) {
                    fprintf (stderr, "CCtsp_get_clique failed\n");
                    goto CLEANUP;
                }

                rval = CCpq_apply_clique (&pqt, c, &status);
                if (rval) {
                    fprintf (stderr, "CCpq_apply_clique failed\n");
                    goto CLEANUP;
                }

                if (status == CCpq_STATUS_NOSOL) {
                    conflict_list[conflict_cnt] = cliquenums[perm[i]];
                    conflict_cnt++;
                    clique_cnt = i;
                    conflict_found = 1;
                    break;
                }
            }
        }
        if (conflict_found == 0) {
#ifdef DEBUG
            printf ("violation vanished\n");
#endif
            rval = 0;
            goto CLEANUP;
        }
    }

  CLEANUP:
    CCpq_tree_free (&pqt);
    return rval;
}

static int found_ctree_conflict (CCtsp_lpcut_in **cuts, int *cutcount,
        int *conflict_list, int conflict_cnt, CCtsp_cuttree *ctree,
        CCtsp_lpcuts *pool)
{
    int clique_cnt = 0;
    CCtsp_lpclique *conflict_cliques[4];
    int nodecount = ctree->nodecount;
    CCpq_tree pqt;
    struct conflict_data_struct conflict_data;
    int status;
    int i;
    int rval;

    CCpq_tree_init (&pqt);

    if (conflict_cnt > 4) {
        fprintf (stderr, "Whoa, bogus conflict count %d\n", conflict_cnt);
        rval = 1; goto CLEANUP;
    }

    for (i=0; i<conflict_cnt; i++) {
        rval = CCtsp_get_clique (pool, conflict_list[i], &conflict_cliques[i]);
        if (rval) {
            fprintf (stderr, "CCtsp_get_clique failed\n");
            goto CLEANUP;
        }
    }

    clique_cnt = conflict_cnt;

    for (;;) {
        rval = CCpq_tree_trivial (&pqt, nodecount, ctree->extern_node);
        if (rval) {
            fprintf (stderr, "CCpq_tree_trivial failed\n");
            goto CLEANUP;
        }
        for (i=0; i<clique_cnt; i++) {
            rval = CCpq_apply_clique (&pqt, conflict_cliques[i], &status);
            if (rval) {
                fprintf (stderr, "CCpq_apply_clique failed\n");
                goto CLEANUP;
            }
            if (status == CCpq_STATUS_NOSOL) {
                if (i == 3) {
                    rval = found_comb_conflict (cuts, cutcount, nodecount,
                                                conflict_cliques);
                    if (rval) {
                        fprintf (stderr, "found_comb_conflict failed\n");
                    }
                    goto CLEANUP;
                }
#ifdef DEBUG
                printf ("conflict with < 4 cliques\n");
#endif
                rval = 0; goto CLEANUP;
            }
        }
        if (clique_cnt >= 4) {
#ifdef DEBUG
            printf ("no conflict with 4 cliques\n");
#endif
            rval = 0; goto CLEANUP;
        }

        conflict_data.pqt = &pqt;
        conflict_data.conflict = (CCtsp_lpclique *) NULL;

        rval = CCpq_cuttree_gen_cliques (ctree, (void *) &conflict_data,
                                    conflict_callback);
        if (rval) {
            fprintf (stderr, "CCpq_cuttree_gen_cliques failed\n");
            goto CLEANUP;
        }

        if (conflict_data.conflict == (CCtsp_lpclique *) NULL) {
#ifdef DEBUG
            printf ("conflict vanished\n");
#endif
            rval = 0; goto CLEANUP;
        }

        conflict_cliques[clique_cnt++] = conflict_data.conflict;
        CCpq_tree_free (&pqt);
    }

  CLEANUP:
    CCpq_tree_free (&pqt);
    for (i=conflict_cnt; i<clique_cnt; i++) {
        CCtsp_free_lpclique (conflict_cliques[i]);
        CC_IFFREE (conflict_cliques[i], CCtsp_lpclique);
    }
    return rval;
}

static int conflict_callback (int *arr, int cnt, int *stop, void *u_data)
{
    struct conflict_data_struct *conflict_data =
        (struct conflict_data_struct *) u_data;
    CCpq_tree *pqt = conflict_data->pqt;
    int rval;
    int status;
    int i;

    CCpq_clear_leaflist (pqt);
    for (i=0; i<cnt; i++) {
        CCpq_add_leaflist (pqt, arr[i]);
    }
    rval = CCpq_apply (pqt, &status);
    if (rval) {
        fprintf (stderr, "CCpq_apply failed\n");
        goto CLEANUP;
    }

    if (status == CCpq_STATUS_NOSOL) {
        conflict_data->conflict = CC_SAFE_MALLOC (1, CCtsp_lpclique);
        if (conflict_data->conflict == (CCtsp_lpclique *) NULL) {
            fprintf (stderr, "Out of memory in conflict_callback\n");
            rval = 1; goto CLEANUP;
        }
        rval = CCtsp_array_to_lpclique (arr, cnt, conflict_data->conflict);
        if (rval) {
            fprintf (stderr, "CCtsp_array_to_lpclique failed\n");
            rval = 1; goto CLEANUP;
        }
        *stop = 1;
    }

    rval = 0;
  CLEANUP:
    if (rval) {
        CC_IFFREE (conflict_data->conflict, CCtsp_lpclique);
    }
    return rval;
}

#define HANDLE 1
#define TOOTH1 2
#define TOOTH2 4
#define TOOTH3 8
#define NREGIONS 16

static const int flipseq[NREGIONS] = {0, TOOTH1, TOOTH2, TOOTH3,
    TOOTH1 | TOOTH2, TOOTH1 | TOOTH3, TOOTH2 | TOOTH3,
    TOOTH1 | TOOTH2 | TOOTH3, HANDLE | TOOTH1, HANDLE | TOOTH2,
    HANDLE | TOOTH3, HANDLE | TOOTH1 | TOOTH2, HANDLE | TOOTH1 | TOOTH3,
    HANDLE | TOOTH2 | TOOTH3};

static int found_comb_conflict (CCtsp_lpcut_in **cuts, int *cutcount,
        int nodecount, CCtsp_lpclique *conflict_cliques[4])
{
    int *mark = (int *) NULL;
    int i;
    int j;
    int k;
    int rval;
    int flip;
    int handle;
    int hbit;
    unsigned int hist[NREGIONS];
    unsigned int hist2[NREGIONS];

    mark = CC_SAFE_MALLOC (nodecount, int);
    if (mark == (int *) NULL) {
        fprintf (stderr, "Out of memory in found_comb_conflict\n");
        rval = 1; goto CLEANUP;
    }

    for (i=0; i<4; i++) {
        CC_FOREACH_NODE_IN_CLIQUE (j, *conflict_cliques[i], k) {
            mark[j] = 0;
        }
    }

    for (i=0; i<4; i++) {
        CC_FOREACH_NODE_IN_CLIQUE (j, *conflict_cliques[i], k) {
            mark[j] |= 1<<i;
        }
    }

    hist[0] = nodecount;
    for (i=1; i<NREGIONS; i++) hist[i] = 0;
    for (i=0; i<4; i++) {
        CC_FOREACH_NODE_IN_CLIQUE (j, *conflict_cliques[i], k) {
            hist[mark[j]]++;
            hist[0]--;
            mark[j] = 0;
        }
    }

    CC_FREE (mark, int);

    for (flip = 0; flip < NREGIONS; flip++) {
        for (handle = 0; handle < 4; handle++) {
            hbit = (1 << handle);
            for (i = 0; i < NREGIONS; i++) {
                j = i & ~(1 | hbit);
                if (i & 1)
                    j |= hbit;
                if (i & hbit)
                    j |= 1;
                j ^= flipseq[flip];
                hist2[i] = hist[j];
            }
            if (histok (hist2)) {
                rval = add_conflict_to_cuts (cuts, cutcount, conflict_cliques,
                        nodecount);
                if (rval) {
                    fprintf (stderr, "add_conflict_to_cliques failed\n");
                    goto CLEANUP;
                }
                rval = 0;
                goto CLEANUP;
            }
        }
    }
    rval = 0;

  CLEANUP:
    CC_IFFREE (mark, int);
    return rval;
}

static int histok (unsigned int *hist)
{
    if (!hist[TOOTH1])
        return 0;               /* tooth 1 has a cavity */
    if (!hist[TOOTH2])
        return 0;               /* tooth 2 has a cavity */
    if (!hist[TOOTH3])
        return 0;               /* tooth 3 has a cavity */
    if (!hist[TOOTH1 | HANDLE])
        return 0;               /* tooth 1 intersects the handle */
    if (!hist[TOOTH2 | HANDLE])
        return 0;               /* tooth 1 intersects the handle */
    if (!hist[TOOTH3 | HANDLE])
        return 0;               /* tooth 1 intersects the handle */
    if (hist[TOOTH1 | TOOTH2] ||/* Nothing else happens */
        hist[TOOTH1 | TOOTH3] ||
        hist[TOOTH2 | TOOTH3] ||
        hist[TOOTH1 | TOOTH2 | TOOTH3] ||
        hist[TOOTH1 | TOOTH2 | HANDLE] ||
        hist[TOOTH1 | TOOTH3 | HANDLE] ||
        hist[TOOTH2 | TOOTH3 | HANDLE] ||
        hist[TOOTH1 | TOOTH2 | TOOTH3 | HANDLE]) {
        return 0;
    }
    return 1;
}


static int add_conflict_to_cuts (CCtsp_lpcut_in **cuts, int *cutcount,
        CCtsp_lpclique *conflict_cliques[4], int ncount)
{
    CCtsp_lpcut_in *c = (CCtsp_lpcut_in *) NULL;
    int rval;
    int i;

    c = CC_SAFE_MALLOC (1, CCtsp_lpcut_in);
    if (c == (CCtsp_lpcut_in *) NULL) {
        fprintf (stderr, "Out of memory in add_conflict_to_cuts\n");
        rval = 1; goto CLEANUP;
    }
    CCtsp_init_lpcut_in (c);
    c->cliquecount = 4;
    c->cliques = CC_SAFE_MALLOC (4, CCtsp_lpclique);
    if (c->cliques == (CCtsp_lpclique *) NULL) {
        fprintf (stderr, "Out of memory in add_conflict_to_cuts\n");
        rval = 1; goto CLEANUP;
    }

    for (i=0; i<4; i++) {
        c->cliques[i].nodes = (CCtsp_segment *) NULL;
        c->cliques[i].segcount = 0;
    }

    for (i=0; i<4; i++) {
        rval = CCtsp_copy_lpclique (conflict_cliques[i], &c->cliques[i]);
        if (rval) {
            fprintf (stderr, "CCtsp_copy_lpclique failed\n");
            goto CLEANUP;
        }
    }

    c->rhs = 10;
    c->sense = 'G';
    c->branch = 0;
    rval = CCtsp_construct_skeleton (c, ncount);
    if (rval) {
        fprintf (stderr, "CCtsp_construct_skeleton failed\n");
        goto CLEANUP;
    }
        
    c->next = *cuts;
    *cuts = c;
    (*cutcount)++;

    rval = 0;

  CLEANUP:
    if (rval) {
        if (c != (CCtsp_lpcut_in *) NULL) {
            if (c->cliques != (CCtsp_lpclique *) NULL) {
                for (i=0; i<4; i++) {
                    CCtsp_free_lpclique (&c->cliques[i]);
                }
                CC_FREE (c->cliques, CCtsp_lpclique);
            }
            CC_FREE (c, CCtsp_lpcut_in);
        }
    }
    return rval;
}
