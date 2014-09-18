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
/*                    Interface to the Cutters                              */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: February 17, 1997                                                 */
/*                                                                          */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCtsp_connect_cuts (CCtsp_lpcut_in **cuts, int *cutcount,           */
/*      int ncount, int ecount, int *elist, double *x)                      */
/*    FINDS violated subtour inequalities via connectivity.                 */
/*     -cuts will return any new cuts found (they will be added to the      */
/*      head of the linked list)                                            */
/*     -cutcount will return the number of new cuts added                   */
/*     -ncount is the number of nodes                                       */
/*     -ecount is the number of edges                                       */
/*     -elist contains the LP edges in node node format                     */
/*     -x is an LP solution                                                 */
/*                                                                          */
/*  int CCtsp_segment_cuts (CCtsp_lpcut_in **cuts, int *cutcount,           */
/*      int ncount, int ecount, int *elist, double *x)                      */
/*    FINDS violated subtour inequalities via linsub.                       */
/*                                                                          */
/*  int CCtsp_exact_subtours (CCtsp_lpcut_in **cuts, int *cutcount,         */
/*      int ncount, int ecount, int *elist, double *x)                      */
/*    FINDS violated subtour inequalities via a mincut algorithm.           */
/*                                                                          */
/*  int CCtsp_tighten_lp (CCtsp_lpcuts *cuts, CCtsp_tighten_info *stats,    */
/*      CCtsp_lpcut_in **cutsout, int *cutcount, int ncount,                */
/*      int ecount, int *elist, double *x, double testtol,                  */
/*      int maxcuts, double *viol, CCrandstate *rstate)                     */
/*    CALLS tighten for each cut in the cuts.                               */
/*     -stats contains some running statistics of tighten                   */
/*     -cutsout returns the tightened cuts that are violated (they are      */
/*      added to the tail of the linked list)                               */
/*     -cutcount is the number of cuts in cutsout                           */
/*     -testtol is a tolerance for calling tighten (call only when the      */
/*      cut has slack value within testtol)                                 */
/*     -maxcuts is a bound on the number of cuts to be returned             */
/*                                                                          */
/*  int CCtsp_double_decker_lp (CCtsp_lpcuts *cuts,                         */
/*      CCtsp_tighten_info *stats, CCtsp_lpcut_in **cutsout,                */
/*      int *cutcount, int ncount, int ecount, int *elist, double *x,       */
/*      double testtol, int maxcuts, double *viol, CCrandstate *rstate)     */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCtsp_cliquetree_lp (CCtsp_lpcuts *cuts,                            */
/*      CCtsp_tighten_info *stats, CCtsp_lpcut_in **cutsout,                */
/*      int *cutcount, int ncount, int ecount, int *elist, double *x,       */
/*      double testtol, int maxcuts, double *viol, CCrandstate *rstate)     */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCtsp_star_lp (CCtsp_lpcuts *cuts, CCtsp_tighten_info *stats,       */
/*      CCtsp_lpcut_in **cutsout, int *cutcount, int ncount,                */
/*      int ecount, int *elist, double *x, double testtol,                  */
/*      int maxcuts, double *viol, CCrandstate *rstate)                     */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCtsp_handling_lp (CCtsp_lpcuts *cuts, CCtsp_tighten_info *stats,   */
/*      CCtsp_lpcut_in **cutsout, int *cutcount, int ncount,                */
/*      int ecount, int *elist, double *x, double testtol,                  */
/*      int maxcuts, double *viol, CCrandstate *rstate)                     */
/*    CALLS CCtsp_comb_handling for each comb in cuts.                      */
/*     -agruments as in CCtsp_tighten_lp.                                   */
/*                                                                          */
/*  int CCtsp_handling_lp (CCtsp_lpcuts *cuts, CCtsp_tighten_info *stats,   */
/*      CCtsp_lpcut_in **cutsout, int *cutcount, int ncount,                */
/*      int ecount, int *elist, double *x, double testtol,                  */
/*      int maxcuts, double *viol, CCrandstate *rstate)                     */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCtsp_teething_lp (CCtsp_lpcuts *cuts, CCtsp_tighten_info *stats,   */
/*      CCtsp_lpcut_in **cutsout, int *cutcount, int ncount,                */
/*      int ecount, int *elist, double *x, double testtol,                  */
/*      int maxcuts, double *viol, CCrandstate *rstate)                     */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCtsp_file_cuts (char *cutfile, CCtsp_lpcut_in **cuts,              */
/*      int *cutcount, int ncount, int *tour)                               */
/*    READS a set of cuts from a file; the format of the cuts can be        */
/*     found by examining the code                                          */
/*     -cutfile is an asci file with a list of cuts                         */
/*     -cuts will return any new cuts found (they will be added to the      */
/*      tail of the linked list)                                            */
/*     -cutcount with return the number of new cuts added                   */
/*     -ncount is the number of nodes                                       */
/*     -tour the permutation tour (used to map the incoming nodes)          */
/*                                                                          */
/*  int CCtsp_file_cuts_write (const char *cutfile, CCtsp_lpcuts *cuts,     */
/*      int *tour)                                                          */
/*    WRITES a set of cuts in a text file that can be read by               */
/*     tsp_file_cuts                                                        */
/*     -cutfile is the name of the file to be written                       */
/*     -cuts is the set of cuts to be written                               */
/*     -tour is a permutation tour (used to map the outgoing nodes)         */
/*                                                                          */
/*  int CCtsp_test_pure_comb (int ncount, CCtsp_lpcut_in *c, int *yes_no,   */
/*      int *handle)                                                        */
/*    TEST if the cut is a comb (without flipped teeth or intersections)    */
/*     -ncount is the number of nodes in the TSP                            */
/*     -yes_no will be set to either 0 or 1, with 1 meaning yes             */
/*     -handle with return the index of the handle if the cut is a comb     */
/*      (handle can be NULL)                                                */
/*                                                                          */
/*  int CCtsp_test_pseudocomb (int ncount, CCtsp_lpcut_in *c, int handle,   */
/*      int *yes_no)                                                        */
/*    TEST if the cut is a pseudocomb.                                      */
/*     -handle gives the index of the handle of the pseudocomb              */
/*                                                                          */
/*  int CCtsp_test_teeth_disjoint (int ncount, CCtsp_lpcut_in *c,           */
/*      int handle, int *yes_no)                                            */
/*    TEST if the cliques other than handle are pairwise disjoint.          */
/*     -yes_no is 1 if disjoint and 0 otherwise.                            */
/*                                                                          */
/*  int CCtsp_find_pure_handle (int ncount, CCtsp_lpcut_in *c,              */
/*      int *handle)                                                        */
/*    FINDS a clique that is c's handle if c is a comb; the search          */
/*     assumes that the teeth are disjoint, so if the comb has              */
/*     extra intersections then a tooth may be returned.                    */
/*     -handle returns the potential handle (it will return -1 if no        */
/*      clique is a potential handle)                                       */
/*                                                                          */
/*  int CCtsp_buildcut_begin (CCtsp_cutinfo *cuts, int init_cliquecount)    */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCtsp_buildcut_addclique (CCtsp_cutinfo *cuts, *arr, int size)      */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCtsp_buildcut_finish (CCtsp_cutinfo *cuts, int rhs)                */
/*    MISSING                                                               */
/*                                                                          */
/*  void CCtsp_buildcut_abort (CCtsp_cutinfo *cuts)                         */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCtsp_truncate_cutlist (CCtsp_lpcut_in **cuts, int ncount,          */
/*      int ecount, int *elist, double *x, int maxcuts,                     */
/*      CCrandstate *rstate)                                                */
/*    RETURNS the maxcuts most violated cuts in the linked list, the        */
/*     remaining cuts are freed.                                            */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "macrorus.h"
#include "util.h"
#include "tsp.h"
#include "cut.h"
#include "combs.h"
#include "verify.h"

#ifdef CCtsp_USE_DOMINO_CUTS
int DPseparator (int nnodes, int nedges, int* edges, double* weigh, int* nIneq,
    int** nDominoes, int*** nAset, int*** nBset, int** nHandle, int**** Aset,
    int **** Bset, int*** Handle, const char *boss_name);
#endif

#define X_FLUFF (1e-10)
#undef  DUMP_BUILDCUT

typedef struct exactsub_param {
    int             nodecount;
    int             cutcount;
    CCtsp_lpcut_in *cuts;
} exactsub_param;


static int
    add_segment (double val, int a, int b, void *pass_param),
    add_exact (double val, int count, int *cutarray, void *pass_param),
    work_on_combs_in_lp (CCtsp_lpcuts *cuts, CCtsp_tighten_info *stats,
        CCtsp_lpcut_in **cutsout, int *cutcount, int ncount, int ecount,
        int *elist, double *x, double testtol, int maxcuts, int caller,
        double *viol, CCrandstate *rstate),
    grab_nonzero_x (int ecount, int *elist, double *x, int *new_ecount,
        int **new_elist, double **new_x, double tol),
    comb_to_domino (CCtsp_lpcut_in *c, CCtsp_lpcut_in **d, int ncount,
        int *marks, int *in_marker);


int CCtsp_connect_cuts (CCtsp_lpcut_in **cuts, int *cutcount, int ncount,
        int ecount, int *elist, double *x)
{
    int rval = 0;
    int i, k, ncomp;
    CCtsp_lpcut_in *c     = (CCtsp_lpcut_in *) NULL;
    int *comps      = (int *) NULL;
    int *compscount = (int *) NULL;

    *cutcount = 0;
    rval = CCcut_connect_components (ncount, ecount, elist, x, &ncomp,
                                     &compscount, &comps);
    if (rval) {
        fprintf (stderr, "CCcut_connect_components failed\n"); goto CLEANUP;
    }

    for (i = 0, k = 0; i < ncomp - 1; k += compscount[i], i++) {
        rval = CCtsp_array_to_subtour (&c, comps + k, compscount[i], ncount);
        if (rval) {
            fprintf (stderr, "CCtsp_array_to_subtour failed\n");
            rval = 1; goto CLEANUP;
        }
        c->next = *cuts;
        *cuts = c;
        (*cutcount)++;
    }

CLEANUP:

    CC_IFFREE (comps, int);
    CC_IFFREE (compscount, int);

    return rval;
}

int CCtsp_segment_cuts (CCtsp_lpcut_in **cuts, int *cutcount, int ncount,
        int ecount, int *elist, double *x)
{
    int rval;
    exactsub_param p;
    int i;
    int *endmark = (int *) NULL;

    *cutcount = 0;

    p.nodecount = ncount;
    p.cutcount  = 0;
    p.cuts      = *cuts;

    endmark = CC_SAFE_MALLOC (ncount, int);
    if (endmark == (int *) NULL) {
        fprintf (stderr, "Out of memory in CCtsp_segment_cuts\n");
        rval = 1; goto CLEANUP;
    }

    for (i=0; i<ncount; i++) {
        endmark[i] = 0;
    }
    for (i=0; i<ecount; i++) {
        if (x[i] >= 0.999999) {
            endmark[elist[2*i]]++;
            endmark[elist[2*i+1]]++;
        }
    }
    for (i=0; i<ncount; i++) {
        if (endmark[i] == 2) {
            endmark[i] = CC_LINSUB_NO_END;
        } else {
            endmark[i] = CC_LINSUB_BOTH_END;
        }
    }

    rval = CCcut_linsub (ncount, ecount, endmark, elist, x, 2.0 - 0.0001,
                         (void *) &p, add_segment);
    if (rval) {
        fprintf (stderr, "CCcut_linsub failed\n");
        goto CLEANUP;
    }

    *cutcount = p.cutcount;
    *cuts = p.cuts;

    rval = 0;

  CLEANUP:
    CC_IFFREE (endmark, int);
    return rval;
}

int CCtsp_shrink_subtours (CCtsp_lpcut_in **cuts, int *cutcount, int ncount,
            int ecount, int *elist, double *x)
{
    int rval;
    exactsub_param p;

    *cutcount = 0;
/*
    rval = CCtsp_connect_cuts (cuts, cutcount, ncount, ecount, elist, x);
    if (rval) {
        fprintf (stderr, "CCtsp_connect_cuts failed\n"); goto CLEANUP;
    }

    if (*cutcount > 0) {
        rval = 0; goto CLEANUP;
    }
*/

    p.nodecount = ncount;
    p.cutcount  = 0;
    p.cuts      = *cuts;

    rval = CCcut_shrink_cuts (ncount, ecount, elist, x, 2.0 - 0.0001,
                       add_exact, (void *) &p);
    if (rval) {
        fprintf (stderr, "CCcut_violated_cuts failed\n"); goto CLEANUP;
    }

    *cutcount = p.cutcount;
    *cuts = p.cuts;

CLEANUP:

    return rval;
}


int CCtsp_exact_subtours (CCtsp_lpcut_in **cuts, int *cutcount, int ncount,
            int ecount, int *elist, double *x)
{
    int rval;
    exactsub_param p;

    *cutcount = 0;
    rval = CCtsp_connect_cuts (cuts, cutcount, ncount, ecount, elist, x);
    if (rval) {
        fprintf (stderr, "CCtsp_connect_cuts failed\n"); goto CLEANUP;
    }

    if (*cutcount > 0) {
        rval = 0; goto CLEANUP;
    }

    p.nodecount = ncount;
    p.cutcount  = 0;
    p.cuts      = *cuts;

    rval = CCcut_violated_cuts (ncount, ecount, elist, x, 2.0 - 0.0001,
                       add_exact, (void *) &p);
    if (rval) {
        fprintf (stderr, "CCcut_violated_cuts failed\n"); goto CLEANUP;
    }

    *cutcount = p.cutcount;
    *cuts = p.cuts;

#if 0
  - this is just to check the values of the exact cuts
    if (*cutcount) {
        CCtsp_lpgraph lg;
        CCtsp_lpcut_in *c;
        double t;

        CCtsp_init_lpgraph_struct (&lg);

        rval = CCtsp_build_lpgraph (&lg, ncount, ecount, elist, (int *) NULL);
        if (rval) {
            fprintf (stderr, "CCtsp_build_lpgraph failed\n"); goto CLEANUP;
        }
        rval = CCtsp_build_lpadj (&lg, 0, ecount);
        if (rval) {
            CCtsp_free_lpgraph (&lg);
            fprintf (stderr, "CCtsp_build_lpadj failed\n"); goto CLEANUP;
        }
        for (c = p.cuts; c; c = c->next) {
            t = CCtsp_cutprice (&lg, c, x);
            printf ("[%f] ", 2.0 + t); fflush (stdout);
        }
        printf ("\n"); fflush (stdout);
        CCtsp_free_lpgraph (&lg);
    }
#endif
    rval = 0;

CLEANUP:

    return rval;
}

static int add_segment (double val, int a, int b, void *pass_param)
{
    int rval = 0;
    exactsub_param *p = (exactsub_param *) pass_param;
    CCtsp_lpcut_in *c = (CCtsp_lpcut_in *) NULL;

    if (val > 2.0) {
        printf ("Warning: Cut of value %f in add_segment\n", val);
        fflush (stdout);
        goto CLEANUP;
    }

    rval = CCtsp_segment_to_subtour (&c, a, b, p->nodecount);
    if (rval) {
        fprintf (stderr, "CCtsp_segment_to_subtour failed\n");
        rval = 1; goto CLEANUP;
    }
    c->next = p->cuts;
    p->cuts = c;
    p->cutcount++;

CLEANUP:

    return rval;
}

static int add_exact (double val, int count, int *cutarray, void *pass_param)
{
    int rval = 0;
    exactsub_param *p = (exactsub_param *) pass_param;
    CCtsp_lpcut_in *c = (CCtsp_lpcut_in *) NULL;

    if (count >= p->nodecount) goto CLEANUP;

    if (val > 2.0) {
        printf ("Warning: Cut of value %f in add_exact\n", val);
        fflush (stdout);
        goto CLEANUP;
    }

    rval = CCtsp_array_to_subtour (&c, cutarray, count, p->nodecount);
    if (rval) {
        fprintf (stderr, "CCtsp_array_to_subtour failed\n");
        rval = 1; goto CLEANUP;
    }
    c->next = p->cuts;
    p->cuts = c;
    p->cutcount++;

CLEANUP:

    return rval;
}

int CCtsp_tighten_lp (CCtsp_lpcuts *cuts, CCtsp_tighten_info *stats,
        CCtsp_lpcut_in **cutsout, int *cutcount, int ncount, int ecount,
        int *elist, double *x, double testtol, int maxcuts,
        double *viol, CCrandstate *rstate)
{
    CCtsp_lpcut_in new, old;
    CCtsp_lpcut_in *c;
    int i;
    int rval = 0;
    double improve;
    CCtsp_lpgraph lg;
    double *newx = (double *) NULL;
    int *newelist = (int *) NULL;
    int newecount;
    CCtsp_lpcut_in **clist = (CCtsp_lpcut_in **) NULL;
    double *vlist = (double *) NULL;
    double maxviol = 0.0;
    int clistsize = 0, vlistsize = 0;
    int count = 0;
    int *perm = (int *) NULL;
    double *cutval = (double *) NULL;

    *cutcount = 0;
    if (!cuts || !cuts->cutcount) return 0;


    rval = grab_nonzero_x (ecount, elist, x, &newecount, &newelist, &newx,
                           X_FLUFF);
    if (rval) {
        fprintf (stderr, "grab_nonzero_x failed\n"); goto CLEANUP;
    }

    cutval = CC_SAFE_MALLOC (cuts->cutcount, double);
    if (!cutval) {
        fprintf (stderr, "out of memory in CCtsp_tighten_lp\n");
        rval = 1; goto CLEANUP;
    }
    rval = CCtsp_price_cuts (cuts, ncount, newecount, newelist, newx, cutval);
    if (rval) {
        fprintf (stderr, "CCtsp_price_cuts failed\n"); goto CLEANUP;
    }

    CCtsp_init_lpgraph_struct (&lg);

    rval = CCtsp_build_lpgraph (&lg, ncount, newecount, newelist,
                                (int *) NULL);
    if (rval) {
        fprintf (stderr, "CCtsp_build_lpgraph failed\n"); goto CLEANUP;
    }
    CC_FREE (newelist, int);
    rval = CCtsp_build_lpadj (&lg, 0, newecount);
    if (rval) {
        fprintf (stderr, "CCtsp_build_lpadj failed\n"); goto CLEANUP;
    }

    for (i = 0; i < cuts->cutcount; i++) {
        if (cutval[i] < testtol && !cuts->cuts[i].branch
            && (cuts->cuts[i].cliquecount > 1 || cutval[i] < 0.1*testtol)
            /* && cuts->cuts[i].age < 3 */) {
            rval = CCtsp_lpcut_to_lpcut_in (cuts, &(cuts->cuts[i]), &old);
            if (rval) {
                fprintf (stderr, "CCtsp_lpcut_to_lpcut_in failed\n");
                goto CLEANUP;
            }
            rval = CCtsp_tighten_lpcut_in (&lg, &old, newx, &new, stats,
                                           &improve);
            if (rval) {
                fprintf (stderr, "CCtsp_tighten_lpcut failed\n");
                goto CLEANUP;
            }
            CCtsp_free_lpcut_in (&old);

            if (improve - cutval[i] > CCtsp_MIN_VIOL) {
                c = CC_SAFE_MALLOC (1, CCtsp_lpcut_in);
                CCcheck_NULL (c, "out of memory in CCtsp_tighte_lp");
                *c = new;

                if (count >= clistsize) {
                    rval = CCutil_reallocrus_scale ((void **) &clist,
                                &clistsize, count + 1, 1.3,
                                sizeof (CCtsp_lpcut_in *));
                    if (rval) {
                        fprintf (stderr, "CCutil_reallocrus_scale failed\n");
                        rval = 1; goto CLEANUP;
                    }
                }
                if (count >= vlistsize) {
                    rval = CCutil_reallocrus_scale ((void **) &vlist,
                                &vlistsize, count + 1, 1.3, sizeof (double));
                    if (rval) {
                        fprintf (stderr, "CCutil_reallocrus_scale failed\n");
                        rval = 1; goto CLEANUP;
                    }
                }
                clist[count] = c;
                vlist[count] = cutval[i] - improve;
                count++;
            } else {
                CCtsp_free_lpcut_in (&new);
            }
        }
    }

    if (count) {
        perm = CC_SAFE_MALLOC (count, int);
        if (!perm) {
            fprintf (stderr, "out of memory in CCtsp_tighten_lp\n");
            rval = 1; goto CLEANUP;
        }
        for (i = 0; i < count; i++) {
            perm[i] = i;
        }
        if (count > maxcuts) {
            CCutil_rselect (perm, 0, count - 1, maxcuts, vlist, rstate);
            for (i = maxcuts; i < count; i++) {
                CCtsp_free_lpcut_in (clist[perm[i]]);
            }
            count = maxcuts;
        }
        for (i = 0; i < count; i++) {
            if (vlist[perm[i]] < maxviol)
                maxviol = vlist[perm[i]];
            clist[perm[i]]->next = *cutsout;
            *cutsout = clist[perm[i]];
        }
    }

    *cutcount = count;
    if (viol) *viol = -maxviol;

CLEANUP:

    CC_IFFREE (newelist, int);
    CC_IFFREE (newx, double);
    CC_IFFREE (clist, CCtsp_lpcut_in *);
    CC_IFFREE (vlist, double);
    CC_IFFREE (perm, int);
    CC_IFFREE (cutval, double);
    CCtsp_free_lpgraph (&lg);
    return rval;
}

#define CALL_TEETHING   1
#define CALL_DDECKER    2
#define CALL_CLIQUETREE 3
#define CALL_STAR       4
#define CALL_HANDLING   5

int CCtsp_double_decker_lp (CCtsp_lpcuts *cuts, CCtsp_tighten_info *stats,
        CCtsp_lpcut_in **cutsout, int *cutcount, int ncount, int ecount,
        int *elist, double *x, double testtol, int maxcuts,
        double *viol, CCrandstate *rstate)
{
    int rval = 0;

    rval = work_on_combs_in_lp (cuts, stats, cutsout, cutcount, ncount, ecount,
                elist, x, testtol, maxcuts, CALL_DDECKER, viol, rstate);
    if (rval) {
        fprintf (stderr, "work_on_combs_in_lp failed\n");
        goto CLEANUP;
    }

CLEANUP:

    return rval;
}

int CCtsp_teething_lp (CCtsp_lpcuts *cuts, CCtsp_tighten_info *stats,
        CCtsp_lpcut_in **cutsout, int *cutcount, int ncount, int ecount,
        int *elist, double *x, double testtol, int maxcuts,
        double *viol, CCrandstate *rstate)
{
    int rval = 0;

    rval = work_on_combs_in_lp (cuts, stats, cutsout, cutcount, ncount, ecount,
                elist, x, testtol, maxcuts, CALL_TEETHING, viol, rstate);
    if (rval) {
        fprintf (stderr, "work_on_combs_in_lp failed\n");
        goto CLEANUP;
    }

CLEANUP:

    return rval;
}

int CCtsp_cliquetree_lp (CCtsp_lpcuts *cuts, CCtsp_tighten_info *stats,
        CCtsp_lpcut_in **cutsout, int *cutcount, int ncount, int ecount,
        int *elist, double *x, double testtol, int maxcuts,
        double *viol, CCrandstate *rstate)
{
    int rval = 0;

    rval = work_on_combs_in_lp (cuts, stats, cutsout, cutcount, ncount, ecount,
                elist, x, testtol, maxcuts, CALL_CLIQUETREE, viol, rstate);
    if (rval) {
        fprintf (stderr, "work_on_combs_in_lp failed\n");
        goto CLEANUP;
    }

CLEANUP:

    return rval;
}

int CCtsp_star_lp (CCtsp_lpcuts *cuts, CCtsp_tighten_info *stats,
        CCtsp_lpcut_in **cutsout, int *cutcount, int ncount, int ecount,
        int *elist, double *x, double testtol, int maxcuts,
        double *viol, CCrandstate *rstate)
{
    int rval = 0;

    rval = work_on_combs_in_lp (cuts, stats, cutsout, cutcount, ncount, ecount,
                elist, x, testtol, maxcuts, CALL_STAR, viol, rstate);
    if (rval) {
        fprintf (stderr, "work_on_combs_in_lp failed\n");
        goto CLEANUP;
    }

CLEANUP:

    return rval;
}

int CCtsp_handling_lp (CCtsp_lpcuts *cuts, CCtsp_tighten_info *stats,
        CCtsp_lpcut_in **cutsout, int *cutcount, int ncount, int ecount,
        int *elist, double *x, double testtol, int maxcuts,
        double *viol, CCrandstate *rstate)
{
    int rval = 0;

    rval = work_on_combs_in_lp (cuts, stats, cutsout, cutcount, ncount, ecount,
                elist, x, testtol, maxcuts, CALL_HANDLING, viol, rstate);
    if (rval) {
        fprintf (stderr, "work_on_combs_in_lp failed\n");
        goto CLEANUP;
    }

CLEANUP:

    return rval;
}

static int work_on_combs_in_lp (CCtsp_lpcuts *cuts, CCtsp_tighten_info *stats,
        CCtsp_lpcut_in **cutsout, int *cutcount, int ncount, int ecount,
        int *elist, double *x, double testtol, int maxcuts, int caller,
        double *viol, CCrandstate *rstate)
{
    CCtsp_lpcut_in new, old;
    CCtsp_lpcut_in *c  = (CCtsp_lpcut_in *) NULL;
    CCtsp_lpcut_in *dd = (CCtsp_lpcut_in *) NULL;
    CCtsp_lpcut_in *ddnext;
    int i, test;
    int rval = 0;
    double improve, newslack, dslack;
    CCtsp_lpgraph lg;
    CC_GCgraph gg;
    double *newx = (double *) NULL;
    int *newelist = (int *) NULL;
    int newecount;
    CCtsp_lpcut_in **clist = (CCtsp_lpcut_in **) NULL;
    double *vlist = (double *) NULL;
    double maxviol = 0.0;
    int clistsize = 0, vlistsize = 0;
    int count = 0;
    int *perm = (int *) NULL;
    double *cutval = (double *) NULL;

    *cutcount = 0;
    if (!cuts || !cuts->cutcount) return 0;

    CCtsp_init_lpgraph_struct (&lg);
    CCcombs_GC_init_graph (&gg);

    rval = grab_nonzero_x (ecount, elist, x, &newecount, &newelist, &newx,
                           X_FLUFF);
    if (rval) {
        fprintf (stderr, "grab_nonzero_x failed\n"); goto CLEANUP;
    }

    cutval = CC_SAFE_MALLOC (cuts->cutcount, double);
    if (!cutval) {
        fprintf (stderr, "out of memory in CCtsp_tighten_lp\n");
        rval = 1; goto CLEANUP;
    }
    rval = CCtsp_price_cuts (cuts, ncount, newecount, newelist, newx, cutval);
    if (rval) {
        fprintf (stderr, "CCtsp_price_cuts failed\n"); goto CLEANUP;
    }

    rval = CCtsp_build_lpgraph (&lg, ncount, newecount, newelist,
                                (int *) NULL);
    if (rval) {
        fprintf (stderr, "CCtsp_build_lpgraph failed\n"); goto CLEANUP;
    }
    if (caller == CALL_DDECKER || caller == CALL_CLIQUETREE ||
        caller == CALL_STAR    || caller == CALL_HANDLING) {
        rval = CCcombs_GC_build_graph (&gg, ncount, newecount, newelist, newx);
        if (rval) {
            fprintf (stderr, "CCcombs_GC_build_graph failed\n"); goto CLEANUP;
        }
    }

    CC_FREE (newelist, int);
    rval = CCtsp_build_lpadj (&lg, 0, newecount);
    if (rval) {
        fprintf (stderr, "CCtsp_build_lpadj failed\n"); goto CLEANUP;
    }

    for (i = 0; i < cuts->cutcount; i++) {
        if (cuts->cuts[i].branch || cuts->cuts[i].cliquecount % 2 ||
            cuts->cuts[i].cliquecount < 4 || cutval[i] >= testtol) {
            continue;
        }
        rval = CCtsp_lpcut_to_lpcut_in (cuts, &(cuts->cuts[i]), &old);
        if (rval) {
            fprintf (stderr, "CCtsp_lpcut_to_lpcut_in failed\n"); goto CLEANUP;
        }
        rval = CCtsp_test_pure_comb (ncount, &old, &test, (int *) NULL);
        if (rval) {
            fprintf (stderr, "CCtsp_test_pure_comb failed\n");
            CCtsp_free_lpcut_in (&old);
            goto CLEANUP;
        }
        if (test == 1) {
            switch (caller) {
            case CALL_TEETHING:
                rval = CCtsp_teething (&lg, newx, &old, &dd);
                if (rval) {
                    fprintf (stderr, "CCtsp_teething failed\n"); goto CLEANUP;
                }
                break;
            case CALL_DDECKER:
                rval = CCtsp_comb_to_double_decker (&lg, &gg, newx, &old, &dd);
                if (rval) {
                    fprintf (stderr, "CCtsp_comb_to_double_decker failed\n");
                    goto CLEANUP;
                }
                break;
            case CALL_CLIQUETREE:
                rval = CCtsp_comb_to_cliquetree (&lg, &gg, newx, &old, &dd);
                if (rval) {
                    fprintf (stderr, "CCtsp_comb_to_cliquetree failed\n");
                    goto CLEANUP;
                }
                break;
            case CALL_STAR:
                rval = CCtsp_comb_to_star (&lg, &gg, newx, &old, &dd);
                if (rval) {
                    fprintf (stderr, "CCtsp_comb_to_star failed\n");
                    goto CLEANUP;
                }
                break;
            case CALL_HANDLING:
                rval = CCtsp_comb_handling (&lg, &gg, newx, &old, &dd);
                if (rval) {
                    fprintf (stderr, "CCtsp_comb_handling failed\n");
                    goto CLEANUP;
                }
                break;
            default:
                fprintf (stderr, "unknown caller in work_on_combs_in_lp\n");
                rval = 1; goto CLEANUP;
            }

            CCtsp_free_lpcut_in (&old);
            while (dd) {
                ddnext = dd->next;
                dslack = CCtsp_cutprice (&lg, dd, newx);
                if (dslack >= 1.0) {
                    CCtsp_free_lpcut_in (dd);
                    CC_FREE (dd, CCtsp_lpcut_in);
                } else {
                    rval = CCtsp_tighten_lpcut_in (&lg, dd, newx, &new,
                                                   stats, &improve);
                    if (rval) {
                        fprintf (stderr, "CCtsp_tighten_lpcut failed\n");
                        goto CLEANUP;
                    }
                    CCtsp_free_lpcut_in (dd);
                    CC_FREE (dd, CCtsp_lpcut_in);

                    newslack = dslack - 2.0*improve;
                    if (-newslack > CCtsp_MIN_VIOL) {
                        c = CC_SAFE_MALLOC (1, CCtsp_lpcut_in);
                        if (!c) {
                            fprintf (stderr,
                               "out of memory in work_on_combs_in_lp\n");
                            CCtsp_free_lpcut_in (&new);
                            rval = 1; goto CLEANUP;
                        }
                        *c = new;

                        if (count >= clistsize) {
                            rval = CCutil_reallocrus_scale ((void **) &clist,
                                    &clistsize, count + 1, 1.3,
                                    sizeof (CCtsp_lpcut_in *));
                            if (rval) {
                                fprintf (stderr,
                                    "CCutil_reallocrus_scale failed\n");
                                rval = 1; goto CLEANUP;
                            }
                        }
                        if (count >= vlistsize) {
                            rval = CCutil_reallocrus_scale ((void **) &vlist,
                                     &vlistsize, count + 1, 1.3,
                                     sizeof (double));
                            if (rval) {
                                fprintf (stderr,
                                    "CCutil_reallocrus_scale failed\n");
                                rval = 1; goto CLEANUP;
                            }
                        }
                        clist[count] = c;
                        vlist[count] = newslack;
                        count++;
                    } else {
                        CCtsp_free_lpcut_in (&new);
                    }
                }
                dd = ddnext;
            }
        } else {
            CCtsp_free_lpcut_in (&old);
        }
    }

    if (count) {
        perm = CC_SAFE_MALLOC (count, int);
        if (!perm) {
            fprintf (stderr, "out of memory in work_on_combs_in_lp\n");
            rval = 1; goto CLEANUP;
        }
        for (i = 0; i < count; i++) {
            perm[i] = i;
        }
        if (count > maxcuts) {
            CCutil_rselect (perm, 0, count - 1, maxcuts, vlist, rstate);
            for (i = maxcuts; i < count; i++) {
                CCtsp_free_lpcut_in (clist[perm[i]]);
            }
            count = maxcuts;
        }
        for (i = 0; i < count; i++) {
            if (vlist[perm[i]] < maxviol)
                maxviol = vlist[perm[i]];
            clist[perm[i]]->next = *cutsout;
            *cutsout = clist[perm[i]];
        }
    }

    *cutcount = count;
    if (viol) *viol = -maxviol;

CLEANUP:

    CC_IFFREE (newelist, int);
    CC_IFFREE (newx, double);
    CC_IFFREE (clist, CCtsp_lpcut_in *);
    CC_IFFREE (vlist, double);
    CC_IFFREE (perm, int);
    CC_IFFREE (cutval, double);
    CCtsp_free_lpgraph (&lg);
    CCcombs_GC_free_graph (&gg);
    if (dd) {
        CCtsp_free_lpcut_in (dd);
    }
    return rval;
}

int CCtsp_domino_trial (CCtsp_lpcut_in **cuts, int *cutcount, int ncount,
        int ecount, int *elist, double *x, CCrandstate *rstate)
{
    int rval = 0;
    int i;
    CCtsp_lpcut_in *c, *d, *cnext;
    CCtsp_lpcut_in *gcuts = (CCtsp_lpcut_in *) NULL;
    int gcutcount = 0;
    int *marks = (int *) NULL;
    int marker = 0;

    *cutcount = 0;

    rval = CCtsp_exactblossom (&gcuts, &gcutcount, ncount, ecount, elist,
                               x, rstate);
    CCcheck_rval (rval, "CCtsp_exactblossom failed");

    marks = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (marks, "out of memory in CCtsp_domino_trial");
    for (i = 0; i < ncount; i++) {
        marks[i] = 0;
    }

    for (c = gcuts; c; c = cnext) {
        cnext = c->next;
        rval = comb_to_domino (c, &d, ncount, marks, &marker); 
        CCcheck_rval (rval, "comb_to_domino failed");
        if (d) {
            d->next = *cuts;
            *cuts = d;
            (*cutcount)++;
        }
        CCtsp_free_lpcut_in (c);
        CC_FREE (c, CCtsp_lpcut_in);
    }

CLEANUP:

    CC_IFFREE (marks, int);
    return rval;
}

static int comb_to_domino (CCtsp_lpcut_in *c, CCtsp_lpcut_in **d, int ncount,
        int *marks, int *in_marker)
{
    int rval = 0;
    int yes_no, ihand, i, j, k, tmp;
    CCtsp_lpclique *hand;
    CCtsp_lpcut_in *new = (CCtsp_lpcut_in *) NULL;
    int marker = *in_marker;
    int incount, outcount, imin, omin;
    int *inside = (int *) NULL;
    int *outside = (int *) NULL;
    int *tmpside;

/*
    printf ("comb_to_domino ...\n"); fflush (stdout);
*/

    *d = (CCtsp_lpcut_in *) NULL;
    CCtsp_init_lpcut_in (new);

    if (!c) goto CLEANUP;

    rval = CCtsp_test_pure_comb (ncount, c, &yes_no, &ihand);
    CCcheck_rval (rval, "CCtsp_test_pure_comb failed");

    if (!yes_no) goto CLEANUP;

    hand = &c->cliques[ihand];

    new = CC_SAFE_MALLOC (1, CCtsp_lpcut_in);
    CCcheck_NULL (new, "out of memory in comb_to_domino");

    CCtsp_init_lpcut_in (new);

    new->cliquecount = 1;
    new->cliques = CC_SAFE_MALLOC (1, CCtsp_lpclique);
    CCcheck_NULL (new->cliques, "out of memory in comb_to_domino");

    rval = CCtsp_copy_lpclique (hand, &new->cliques[0]);
    CCcheck_rval (rval, "CCtsp_copy_lpclique failed");

    new->dominos = CC_SAFE_MALLOC (c->cliquecount - 1, CCtsp_lpdomino);
    CCcheck_NULL (new->dominos, "out of memory in comb_to_domino");
    new->dominocount = c->cliquecount -1;
    for (i = 0; i < new->dominocount; i++) {
        CCtsp_init_lpdomino (&new->dominos[i]);
    }

    marker++;
    CCtsp_mark_clique  (hand, marks, marker);

    for (i = 0, k = 0; i < c->cliquecount; i++) {
        if (i != ihand) {
            incount = outcount = 0;
            CC_FOREACH_NODE_IN_CLIQUE (j, c->cliques[i], tmp) {
                if (marks[j] == marker) {
                    incount++;
                } else {
                    outcount++;
                }
            }
            if (!incount|| !outcount) {
                fprintf (stderr, "error in converted domino\n");
                rval = 1;  goto CLEANUP;
            }
            inside = CC_SAFE_MALLOC (incount, int);
            CCcheck_NULL (inside, "out of memory in comb_to_domino");
            outside = CC_SAFE_MALLOC (outcount, int);
            CCcheck_NULL (inside, "out of memory in comb_to_domino");

            incount = outcount = 0;
            imin = omin = ncount;
            CC_FOREACH_NODE_IN_CLIQUE (j, c->cliques[i], tmp) {
                if (marks[j] == marker) {
                    inside[incount] = j;
                    incount++;
                    if (j < imin) imin = j; 
                } else {
                    outside[outcount] = j;
                    outcount++;
                    if (j < omin) omin = j;
                }
            }
            if (omin < imin) {  /* to get a canonical domino */
                CC_SWAP (incount, outcount, tmp);
                CC_SWAP (inside, outside, tmpside);
            }
            rval = CCtsp_array_to_lpclique (inside, incount,
                                            &new->dominos[k].sets[0]);
            CCcheck_rval (rval, "CCtsp_array_to_lpclique failed");
            rval = CCtsp_array_to_lpclique (outside, outcount,
                                            &new->dominos[k].sets[1]);
            CCcheck_rval (rval, "CCtsp_array_to_lpclique failed");
            k++;
            CC_FREE (inside, int);
            CC_FREE (outside, int);
        }
    }

    new->rhs = 3*new->dominocount + 1;
    new->sense = 'G';
    
    *d = new;

CLEANUP:

/*
    printf ("DONE comb_to_domino: %d\n", rval); fflush (stdout);
*/

    *in_marker = marker;
    if (rval) {
        if (new) CCtsp_free_lpcut_in (new);
    }
    CC_IFFREE (inside, int);
    CC_IFFREE (outside, int);
    return rval;
}


int CCtsp_file_cuts (char *cutfile, CCtsp_lpcut_in **cuts, int *cutcount,
        int ncount, int *tour)
{
    FILE *in = (FILE *) NULL;
    int *inv = (int *) NULL;
    CCtsp_lpcut_in *c;
    CCtsp_lpcut_in **clast;
    int i, j, k;
    int ncliques, size;
    int *icliq = (int *) NULL;
    int rval = 0;

    *cutcount = 0;

    in = fopen (cutfile, "r");
    if  (in == (FILE *) NULL) {
        fprintf (stderr, "unable to open %s for reading\n", cutfile);
        return 0;
    }

    inv = CC_SAFE_MALLOC (ncount, int);
    if (!inv) {
        fprintf (stderr, "out of memory in CCtsp_file_cuts\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < ncount; i++) {
        inv[tour[i]] = i;
    }

    clast = cuts;
    while ((*clast) != (CCtsp_lpcut_in *) NULL) {
        clast = &((*clast)->next);
    }
    
    while (fscanf (in, "%d", &ncliques) != EOF) {
        c = CC_SAFE_MALLOC (1, CCtsp_lpcut_in);
        CCcheck_NULL (c, "out of memory in CCtsp_file_cuts");
        CCtsp_init_lpcut_in (c);

        c->cliquecount = ncliques;
        c->cliques = CC_SAFE_MALLOC (ncliques, CCtsp_lpclique);
        if (!c->cliques) {
            fprintf (stderr, "out of memory in CCtsp_file_cuts\n");
            rval = 1; goto CLEANUP;
        }
        for (i = 0; i < ncliques; i++) {
            fscanf (in, "%d", &size);
            icliq = CC_SAFE_MALLOC (size, int);
            if (!icliq) {
                fprintf (stderr, "out of memory in CCtsp_file_cuts\n");
                rval = 1; goto CLEANUP;
            }
            for (j = 0; j < size; j++) {
                fscanf (in, "%d", &k);
                icliq[j] = inv[k];
            }
            rval = CCtsp_array_to_lpclique (icliq, size, &(c->cliques[i]));
            if (rval) {
                fprintf (stderr, "CCtsp_array_to_lpclique failed\n");
                goto CLEANUP;
            }
            CC_FREE (icliq, int);
        }
        fscanf (in, "%d", &(c->rhs));
        c->sense = 'G';
        c->branch = 0;
        rval = CCtsp_construct_skeleton (c, ncount);
        if (rval) {
            fprintf (stderr, "CCtsp_construct_skeleton failed\n");
            goto CLEANUP;
        }

        (*clast) = c;
        c->next = (CCtsp_lpcut_in *) NULL;
        clast = &(c->next);
        (*cutcount)++;
#if 1
        rval = CCverify_cut (c, CC_TYPE_ALL, &i);
        if (rval) {
            fprintf (stderr, "Invalid file cut\n");
            goto CLEANUP;
        } else {
            printf ("File cut type %d\n", i);
        }
#endif
    }

CLEANUP:

    CC_IFFREE (inv, int);
    fclose (in);
    return  rval;
}

int CCtsp_file_cuts_write (const char *cutfile, CCtsp_lpcuts *cuts, int *tour)
{
    FILE *out = (FILE *) NULL;
    int i, j, k, p;
    int cutcount = cuts->cutcount;
    CCtsp_lpcut *c;
    CCtsp_lpclique *cl;
    int isize;

    out = fopen (cutfile, "w");
    if  (out == (FILE *) NULL) {
        fprintf (stderr, "unable to open %s for writing\n", cutfile);
        return 1;
    }

    for (i = 0; i < cutcount; i++) {
        c = &cuts->cuts[i];
        if (!c->branch) {
            fprintf (out, "%d\n", c->cliquecount);
            for (j = 0; j < c->cliquecount; j++) {
                cl = &cuts->cliques[c->cliques[j]];
                for (k = 0, isize = 0; k < cl->segcount; k++) {
                    isize += (cl->nodes[k].hi - cl->nodes[k].lo + 1);
                }
                fprintf (out, "%d  ", isize);
                CC_FOREACH_NODE_IN_CLIQUE (p, *cl, k) {
                    fprintf (out, "%d ", tour[p]);
                }
                fprintf (out, "\n");
            }
            fprintf (out, "%d\n", c->rhs);
        }
    }

    fclose (out);
    return 0;
}

int CCtsp_buildcut_begin (CCtsp_cutinfo *cuts, int init_cliquecount)
{
    cuts->current = CC_SAFE_MALLOC (1, CCtsp_lpcut_in);
    if (!cuts->current) return -1;
    CCtsp_init_lpcut_in (cuts->current);
    cuts->current->cliques = CC_SAFE_MALLOC (init_cliquecount, CCtsp_lpclique);
    if (!cuts->current->cliques) {
        CC_FREE (cuts->current, CCtsp_lpcut_in);
        return -1;
    }
    return 0;
}

int CCtsp_buildcut_addclique (CCtsp_cutinfo *cuts, int *arr, int size)
{
    int i;
    int *newarr = (int *) NULL;
    int newsize;
    int rval;
    CCtsp_lpcut_in *c = cuts->current;

    if (!c) {
        fprintf (stderr, "Trying to add to nonexistent clique\n");
        return -1;
    }

    rval = CCcut_SRK_expand (&cuts->expand, arr, size, &newarr, &newsize);
    if (rval) {
        fprintf (stderr, "CCcut_SRK_expand failed\n");
        CCtsp_buildcut_abort (cuts);
        return rval;
    }

    rval = CCutil_reallocrus_count ((void **) &(c->cliques), c->cliquecount+1,
                             sizeof (c->cliques[0]));
    if (rval) {
        fprintf (stderr, "couldn't realloc cliques\n");
        CC_IFFREE (newarr, int);
        CCtsp_buildcut_abort (cuts);
        return rval;
    }
    
    i = c->cliquecount;

    rval = CCtsp_array_to_lpclique (newarr, newsize, &(c->cliques[i]));
    if (rval) {
        fprintf (stderr, "CCtsp_array_to_lpclique failed\n");
        CC_IFFREE (newarr, int);
        CCtsp_buildcut_abort (cuts);
        return rval;
    }
    c->cliquecount++;
    CC_IFFREE (newarr, int);
    return 0;
}

void CCtsp_buildcut_abort (CCtsp_cutinfo *cuts)
{
    CCtsp_free_lpcut_in (cuts->current);
    CC_IFFREE (cuts->current, CCtsp_lpcut_in);
}

int CCtsp_buildcut_finish (CCtsp_cutinfo *cuts, int rhs)
{
    CCtsp_lpcut_in *c = cuts->current;
    int rval;

#ifdef DUMP_BUILDCUT
    {
        int i, j, tmp;
        printf ("new buildcut (%d):", c->cliquecount);
        for (i=0; i<c->cliquecount; i++) {
            printf (" (");
            CC_FOREACH_NODE_IN_CLIQUE (j, c->cliques[i], tmp) {
                printf ("%d ",j);
            }
            printf (")");
        }
        printf (" >= %d\n", rhs);
        fflush (stdout);
    }
#endif

    c->rhs = rhs;
    c->sense = 'G';
    c->branch = 0;

    rval = CCtsp_construct_skeleton (c,
            CCcut_SRK_original_ncount (&cuts->expand));
    if (rval) {
        fprintf (stderr, "CCtsp_construct_skeleton failed\n");
        goto CLEANUP;
    }

    c->next = *cuts->clist;
    (*cuts->clist) = c;
    cuts->current = (CCtsp_lpcut_in *) NULL;
    (*cuts->cutcount)++;

    rval = 0;
 CLEANUP:
    if (rval) {
        CCtsp_free_lpcut_in (c);
    }
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

int CCtsp_test_pure_comb (int ncount, CCtsp_lpcut_in *c, int *yes_no,
        int *handle)
{
    int rval = 0;
    int i, marked, ihandle;
    int *marks = (int *) NULL;

    *yes_no = 0;
    if (handle) *handle = -1;

    if (c->cliquecount < 4 || c->cliquecount % 2 ||
        c->sense != 'G') {
        goto CLEANUP;
    }

    rval = CCtsp_find_pure_handle (ncount, c, &ihandle);
    if (rval) {
        fprintf (stderr, "CCtsp_find_pure_handle failed\n");
        goto CLEANUP;
    }
    if (ihandle == -1) goto CLEANUP;

    marks = CC_SAFE_MALLOC (ncount, int);
    if (!marks) {
        fprintf (stderr, "out of memory in CCtsp_test_pure_comb\n");
        rval = 1; goto CLEANUP;
    }
    CCtsp_mark_cut (c, marks, 0);

    CCtsp_mark_clique (&c->cliques[ihandle], marks, 1);
    for (i = 0; i < c->cliquecount; i++) {
        if (i != ihandle) {
            CCtsp_is_clique_marked (&c->cliques[i], marks, 1, &marked);
            if (!marked) goto CLEANUP;
            CCtsp_is_clique_marked (&c->cliques[i], marks, 0, &marked);
            if (!marked) goto CLEANUP;
        }
    }
    CCtsp_mark_clique (&c->cliques[ihandle], marks, 0);

    for (i = 0; i < c->cliquecount; i++) {
        if (i != ihandle) {
            CCtsp_is_clique_marked (&c->cliques[i], marks, 1, &marked);
            if (marked) goto CLEANUP;
            CCtsp_mark_clique (&c->cliques[i], marks, 1);
        }
    }

    *yes_no = 1;
    if (handle) *handle = ihandle;

CLEANUP:

    CC_IFFREE (marks, int);
    return rval;
}

int CCtsp_test_pseudocomb (int ncount, CCtsp_lpcut_in *c, int handle,
        int *yes_no)
{
    int rval = 0;
    int i, k, marked;
    int *ends = (int *) NULL;
    int *marks = (int *) NULL;

    *yes_no = 0;
    if (c->cliquecount <= 1 || c->cliquecount % 2 || c->sense != 'G') {
        printf ("bad cliquecount or sense in pseudocomb\n"); fflush (stdout);
        goto CLEANUP;
    }

    marks = CC_SAFE_MALLOC (ncount, int);
    if (!marks) {
        fprintf (stderr, "out of memory in CCtsp_test_pseudocomb\n");
        rval = 1; goto CLEANUP;
    }
    CCtsp_mark_cut (c, marks, 0);

    /* Teeth intersect H and are not contained in H */

    CCtsp_mark_clique (&c->cliques[handle], marks, 1);
    for (i = 0; i < c->cliquecount; i++) {
        if (i != handle) {
            CCtsp_is_clique_marked (&c->cliques[i], marks, 1, &marked);
            if (!marked) goto CLEANUP;
            CCtsp_is_clique_marked (&c->cliques[i], marks, 0, &marked);
            if (!marked) goto CLEANUP;
        }
    }
    CCtsp_mark_clique (&c->cliques[0], marks, 0);

    /* Big teeth are pairwise disjoint */

    for (i = 0; i < c->cliquecount; i++) {
        if (i != handle) {
            CCtsp_clique_count (&c->cliques[i], &k);
            if (k >= 3) {
                CCtsp_is_clique_marked (&c->cliques[i], marks, 1, &marked);
                if (marked) goto CLEANUP;
                CCtsp_mark_clique (&c->cliques[i], marks, 1);
            }
        }
    }
    for (i = 1; i < c->cliquecount; i++) {
        CCtsp_mark_clique (&c->cliques[i], marks, 0);
    }

    /* No small tooth is contained in a big tooth */

    for (i = 0; i < c->cliquecount; i++) {
        if (i != handle) {
            CCtsp_clique_count (&c->cliques[i], &k);
            if (k >= 3) {
                CCtsp_mark_clique (&c->cliques[i], marks, i + 1);
            }
        }
    }
    for (i = 0; i < c->cliquecount; i++) {
        if (i != handle) {
            CCtsp_clique_count (&c->cliques[i], &k);
            if (k < 3) {
                rval = CCtsp_clique_to_array (&c->cliques[i], &ends, &k);
                if (rval) {
                    fprintf (stderr, "CCtsp_clique_to_array failed\n");
                    goto CLEANUP;
                }
                if (ends[0] != 0 && ends[0] == ends[1]) goto CLEANUP;
                CC_IFFREE (ends, int);
            }
        }
    }


    *yes_no = 1;


CLEANUP:

    CC_IFFREE (marks, int);
    CC_IFFREE (ends, int);
    return rval;
}

int CCtsp_test_teeth_disjoint (int ncount, CCtsp_lpcut_in *c, int handle,
        int *yes_no)
{
    int rval = 0;
    int i, marked;
    int *marks = (int *) NULL;

    *yes_no = 0;

    marks = CC_SAFE_MALLOC (ncount, int);
    if (!marks) {
        fprintf (stderr, "out of memory in CCtsp_teeth_disjoint\n");
        rval = 1; goto CLEANUP;
    }
    CCtsp_mark_cut (c, marks, 0);

    for (i = 0; i < c->cliquecount; i++) {
        if (i != handle) {
            CCtsp_is_clique_marked (&c->cliques[i], marks, 1, &marked);
            if (marked) goto CLEANUP;
            CCtsp_mark_clique (&c->cliques[i], marks, 1);
        }
    }

    *yes_no = 1;

CLEANUP:

    CC_IFFREE (marks, int);
    return rval;
}

int CCtsp_find_pure_handle (int ncount, CCtsp_lpcut_in *c, int *handle)
{
    int rval = 0;
    int *marks = (int *) NULL;
    int i, test;

    *handle = -1;
    if (c->cliquecount % 2 || c->cliquecount < 4) goto CLEANUP;

    marks = CC_SAFE_MALLOC (ncount, int);
    if (!marks) {
        fprintf (stderr, "out of memory in CCtsp_pure_find_handle\n");
        rval = 1; goto CLEANUP;
    }
    CCtsp_mark_cut (c, marks, 0);

    CCtsp_mark_clique (&c->cliques[0], marks, 1);
    CCtsp_is_clique_marked (&c->cliques[1], marks, 1, &test);
    if (test) {
        CCtsp_is_clique_marked (&c->cliques[2], marks, 1, &test);
        if (test) {
            *handle = 0; goto CLEANUP;
        } else {
            *handle = 1; goto CLEANUP;
        }
    } else {
        for (i = 2; i < c->cliquecount; i++) {
            CCtsp_is_clique_marked (&c->cliques[i], marks, 1, &test);
            if (test) {
                *handle = i;
                goto CLEANUP;
            }
        }
    }

CLEANUP:

    CC_IFFREE (marks, int);
    return rval;
}

int CCtsp_truncate_cutlist (CCtsp_lpcut_in **cuts, int ncount, int ecount,
        int *elist, double *x, int maxcuts, CCrandstate *rstate)
{
    int i;
    int count = 0;
    int rval = 0;
    CCtsp_lpcut_in *c, *cnext;
    CCtsp_lpcut_in **clist = (CCtsp_lpcut_in **) NULL;
    double *vlist = (double *) NULL;
    int *perm = (int *) NULL;
    CCtsp_lpgraph lg;

    CCtsp_init_lpgraph_struct (&lg);

    if (maxcuts <= 0) {
        for (c = *cuts; c; c = cnext) {
            cnext = c->next;
            CCtsp_free_lpcut_in (c);
        }
        *cuts = (CCtsp_lpcut_in *) NULL;
        goto CLEANUP;
    }

    for (c = *cuts; c; c = c->next) {
        count++;
    }

    if (count > maxcuts) {
        rval = CCtsp_build_lpgraph (&lg, ncount, ecount, elist, (int *) NULL);
        if (rval) {
            fprintf (stderr, "CCtsp_build_lpgraph failed\n"); goto CLEANUP;
        }
        rval = CCtsp_build_lpadj (&lg, 0, ecount);
        if (rval) {
            fprintf (stderr, "CCtsp_build_lpadj failed\n"); goto CLEANUP;
        }

        vlist = CC_SAFE_MALLOC (count, double);
        clist = CC_SAFE_MALLOC (count, CCtsp_lpcut_in *);
        perm  = CC_SAFE_MALLOC (count, int);
        if (!vlist || !clist || !perm) {
            fprintf (stderr, "out of memory in CCtsp_tighten_lp\n");
            rval = 1; goto CLEANUP;
        }
        for (i = 0, c = *cuts; c; c = c->next, i++) {
            clist[i] = c;
            vlist[i] = CCtsp_cutprice (&lg, c, x);
            perm[i] = i;
        }

        CCutil_rselect (perm, 0, count - 1, maxcuts, vlist, rstate);
        for (i = maxcuts; i < count; i++) {
            CCtsp_free_lpcut_in (clist[perm[i]]);
        }

        *cuts = (CCtsp_lpcut_in *) NULL;
        for (i = 0; i < maxcuts; i++) {
            clist[perm[i]]->next = *cuts;
            *cuts = clist[perm[i]];
        }
    }

CLEANUP:

    CC_IFFREE (clist, CCtsp_lpcut_in *);
    CC_IFFREE (vlist, double);
    CC_IFFREE (perm, int);
    CCtsp_free_lpgraph (&lg);

    return rval;
}

#ifdef CCtsp_USE_DOMINO_CUTS

int DPfindBadEdge(int nnodes, int nedges, int* edges, double *weight);
int DPfindBadEdgeK(int nnodes, int nedges, int* edges, int k);

static int shrink_to_planar_graph (int ncount, int ecount, int *elist,
        double *elen, int *p_oncount, int *p_oecount, int **p_oelist,
        double **p_oelen, CC_SRKexpinfo *expand, int quickshrink,
        int rand_minor, CCrandstate *rstate);
static int try_greedy_cut (int ncount, int ecount, int *elist, double *elen,
      int scount, int *slist, double *cval, int **clist, int *ccount);
static int grab_sorted_edges (CC_SRKgraph *G, int *oncount, int *oecount,
        int **polist, double **polen, CC_SRKexpinfo *expand);
static int build_dp_cut_expand (CCtsp_lpcut_in **cut, int ndomino, int *Acount,
        int **A, int *Bcount, int **B, int handlecount, int *handle,
        CC_SRKexpinfo *expand, int ncount, int *comb);
static int remove_light_edges (int ecount, int *elist, double *elen,
        int *new_ecount, int **new_elist, double **new_elen, double tol);
static void free_raw_dominos (int nIneq, int *nHandle, int **Handle,
    int *nDominoes, int **nAset, int ***Aset, int **nBset, int ***Bset);
static int check_raw_domino (int ncount, int hcount, int *hand, int dcount, 
        int *Acount, int **A, int *Bcount, int **B, int *valid, int *comb);
static void  print_raw_domino (int hcount, int *hand, int dcount, 
        int *Acount, int **A, int *Bcount, int **B);

int CCtsp_new_domino (CCtsp_lpcut_in **cuts, int *cutcount, int ncount,
        int ecount, int *elist, double *x, const char *dombossname)
{
    int rval = 0;
    int i, valid, comb, nIneq = 0, nCombs = 0;
    int *nDominoes = (int *) NULL;
    int **nBset = (int **) NULL;
    int *nHandle = (int *) NULL;
    int ***Aset = (int ***) NULL;
    int ***Bset = (int ***) NULL;
    int **nAset = (int **) NULL;
    int **Handle = (int **) NULL;
    CCtsp_lpcut_in *c;

    printf ("Call DPseparator ...\n"); fflush (stdout);
    if (dombossname) {
        printf ("Use Domino Boss: %s\n", dombossname); fflush (stdout);
    }

    *cutcount = 0;

    rval = DPseparator (ncount, ecount, elist, x, &nIneq, &nDominoes, &nAset,
                        &nBset, &nHandle, &Aset, &Bset, &Handle,
                        dombossname);
    CCcheck_rval (rval, "DPseparator failed");

    for (i = 0; i < nIneq; i++) {
        rval = check_raw_domino (ncount, nHandle[i], Handle[i], nDominoes[i], 
                                 nAset[i], Aset[i], nBset[i], Bset[i], &valid,
                                 &comb);
        CCcheck_rval (rval, "check_raw_domino failed");
        nCombs += comb;

        rval = CCtsp_build_dp_cut (&c, nDominoes[i], nAset[i], Aset[i],
                                   nBset[i], Bset[i], nHandle[i], Handle[i]);
        CCcheck_rval (rval, "CCtsp_build_dp_cut failed");
        c->next = *cuts;
        *cuts = c;
        (*cutcount)++;
    }

    printf ("Found %d domino cuts (%d combs) ...\n", nIneq, nCombs);
    fflush (stdout);

CLEANUP:

    free_raw_dominos (nIneq, nHandle, Handle, nDominoes, nAset, Aset, nBset,
                      Bset);
    return rval;
}

int CCtsp_shrink_domino (CCtsp_lpcut_in **cuts, int *cutcount, int ncount,
        int ecount, int *elist, double *x, int quickshrink, int rand_minor,
        CCrandstate *rstate, const char *dombossname)
{
    int rval = 0;
    int i, comb, nIneq = 0, nCombs = 0;
    int *nDominoes = (int *) NULL;
    int **nBset = (int **) NULL;
    int *nHandle = (int *) NULL;
    int ***Aset = (int ***) NULL;
    int ***Bset = (int ***) NULL;
    int **nAset = (int **) NULL;
    int **Handle = (int **) NULL;
    CCtsp_lpcut_in *c;
    CC_SRKexpinfo expand;
    int new_ecount, planar_ecount, planar_ncount;
    int *new_elist = (int *) NULL;
    int *planar_elist = (int *) NULL;
    double *new_x = (double *) NULL;
    double *planar_x = (double *) NULL;

    CCcut_SRK_init_expinfo (&expand);

    printf ("Call DPseparator with shrinking ...\n"); fflush (stdout);
    if (dombossname) {
        printf ("Use Domino Boss: %s\n", dombossname); fflush (stdout);
    }

    *cutcount = 0;

    rval = remove_light_edges (ecount, elist, x, &new_ecount,
                               &new_elist, &new_x, 0.00001);
    CCcheck_rval (rval, "remove_light_edges failed");

    rval = shrink_to_planar_graph (ncount, new_ecount, new_elist, new_x,
                 &planar_ncount, &planar_ecount, &planar_elist, &planar_x,
                 &expand, quickshrink, rand_minor, rstate); 
    CCcheck_rval (rval, "shrink_to_planar_graph failed");

    if (planar_ncount > 6) {
        rval = DPseparator (planar_ncount, planar_ecount, planar_elist,
                 planar_x, &nIneq, &nDominoes, &nAset, &nBset, &nHandle,
                 &Aset, &Bset, &Handle, dombossname);
        CCcheck_rval (rval, "DPseparator failed");

        for (i = 0; i < nIneq; i++) {
            rval = build_dp_cut_expand (&c, nDominoes[i], nAset[i], Aset[i],
                         nBset[i], Bset[i], nHandle[i], Handle[i], &expand,
                         ncount, &comb);
            CCcheck_rval (rval, "_build_dp_cut_expand failed");
            nCombs += comb;
            c->next = *cuts;
            *cuts = c;
            (*cutcount)++;
        }
    } else {
        printf ("Skipping shrunk graph with %d nodes\n", planar_ncount);
        fflush (stdout); 
    }

    printf ("Found %d shunk domino cuts (%d combs) ...\n", nIneq, nCombs);
    fflush (stdout);


CLEANUP:

    CCcut_SRK_free_expinfo (&expand);
    CC_IFFREE (new_elist, int);
    CC_IFFREE (new_x, double);
    CC_IFFREE (planar_elist, int);
    CC_IFFREE (planar_x, double);

    free_raw_dominos (nIneq, nHandle, Handle, nDominoes, nAset, Aset, nBset,
                     Bset);
    return rval;
}

#define PLANAR_SHRINK_TOL  0.1
#define PLANAR_CUT_TOL     0.025

static int shrink_to_planar_graph (int ncount, int ecount, int *elist,
        double *elen, int *p_oncount, int *p_oecount, int **p_oelist,
        double **p_oelen, CC_SRKexpinfo *expand, int quickshrink, 
        int rand_minor, CCrandstate *rstate)
{
    int rval = 0;
    int bad_edge, bad_count;
    int oncount = 0, oecount = 0;
    int *oelist = (int *) NULL;
    double *oelen = (double *) NULL;
    int bad_ends[2];
    double bad_val;
    int *bad_list = (int *) NULL;
    CC_SRKgraph G;

    CCcut_SRK_init_graph (&G);
    rval = CCcut_SRK_buildgraph (&G, ncount, ecount, elist, elen);
    CCcheck_rval (rval, "CCcut_SRK_buildgraph failed");

    rval = grab_sorted_edges (&G, &oncount, &oecount, &oelist, &oelen,
                              expand);
    CCcheck_rval (rval, "grab_sorted_edges failed");
    if (rand_minor <= 1) {
        bad_edge = DPfindBadEdge (oncount, oecount, oelist, oelen);
    } else {
        bad_edge = DPfindBadEdgeK (oncount, oecount, oelist, rand_minor);
    }

    while (bad_edge < oecount && oelen[bad_edge] > PLANAR_SHRINK_TOL) {
        CCcut_SRK_free_expinfo (expand);

        bad_ends[0] = oelist[2*bad_edge];
        bad_ends[1] = oelist[2*bad_edge+1];

        if (quickshrink <= 1) {
            rval = try_greedy_cut (oncount, oecount, oelist, oelen, 2,
                                  bad_ends, &bad_val, &bad_list, &bad_count);
            CCcheck_rval (rval, "try_greedy_cut failed");
            printf ("Greedy Returned Cut Val = %f, Cnt = %d\n", bad_val,
                     bad_count);
            fflush (stdout);

            if (bad_val >= 2.0 + PLANAR_CUT_TOL) {
                CC_IFFREE (bad_list, int);

                rval = CCcut_mincut_containing_set (oncount, oecount, oelist,
                        oelen, 2, bad_ends, &bad_val, &bad_list, &bad_count,
                        quickshrink, rstate);
                CCcheck_rval (rval, "CCcut_mincut_containing_set failed")
            }

            if (bad_count < 2) {
                fprintf (stderr, "illegal containing set\n");
                rval = 1;  goto CLEANUP;
            }

            rval = CCcut_SRK_identify_set (&G, bad_count, bad_list);
            CCcheck_rval (rval, "CCcut_SRK_identify_set failed");

            CC_IFFREE (bad_list, int);
        } else {
            rval = CCcut_SRK_identify_set (&G, 2, bad_ends);
            CCcheck_rval (rval, "CCcut_SRK_identify_set failed");
        }


        CC_IFFREE (oelist, int);
        CC_IFFREE (oelen, double);

        rval = grab_sorted_edges (&G, &oncount, &oecount, &oelist, &oelen,
                                  expand);
        CCcheck_rval (rval, "CCcut_SRK_grab_edges failed");

        if (rand_minor <= 1) {
            bad_edge = DPfindBadEdge (oncount, oecount, oelist, oelen);
        } else {
            bad_edge = DPfindBadEdgeK (oncount, oecount, oelist, rand_minor);
        }
        printf ("Shrunk ecount = %d, ncount = %d, bad edge = %d\n",
                oecount, oncount, bad_edge);
        fflush (stdout);
    }

    *p_oncount = oncount;
    *p_oecount = oecount;
    *p_oelist  = oelist;
    *p_oelen = oelen;

CLEANUP:

    CCcut_SRK_free_graph (&G);
    CC_IFFREE (bad_list, int);
    return rval;
}

static int try_greedy_cut (int ncount, int ecount, int *elist, double *elen,
      int scount, int *slist, double *cval, int **clist, int *ccount)
{
    int i, rval = 0;
    CC_GCgraph G;

    CCcombs_GC_init_graph (&G);

    rval = CCcombs_GC_build_graph (&G, ncount, ecount, elist, elen);
    CCcheck_rval (rval, "CCcombs_GC_build_graph");
    for (i = 0; i < scount; i++) G.nodelist[slist[i]].mark = 2;

    *ccount = scount;
    *clist = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (*clist, "out of memory for clist");

    for (i = 0; i < scount; i++) (*clist)[i] = slist[i];

    printf ("Call greedy with count = %d\n", *ccount); fflush (stdout);

    rval = CCcombs_greedy_cut (&G, ccount, *clist, 2, 0, 0, 0, (int *) NULL,
                               cval);
    CCcheck_rval (rval, "CCcombs_greedy_cut failed");

    printf ("Greedy returns with with count = %d\n", *ccount); fflush (stdout);
    


CLEANUP:

    CCcombs_GC_free_graph (&G);
    return rval;
}

static int grab_sorted_edges (CC_SRKgraph *G, int *oncount, int *oecount,
        int **polist, double **polen, CC_SRKexpinfo *expand)
{
    int rval = 0;
    int i, ncount, ecount;
    int *elist = (int *) NULL;
    int *perm = (int *) NULL;
    int *olist = (int *) NULL;
    double *olen = (double *) NULL;
    double *elen = (double *) NULL;

    rval = CCcut_SRK_grab_edges (G, &ncount, &ecount, &elist, &elen, expand);
    CCcheck_rval (rval, "CCcut_SRK_graph_edges failed");

    perm = CC_SAFE_MALLOC (ecount, int);
    CCcheck_NULL (perm, "out of memory for perm");
    for (i = 0; i < ecount; i++) {
        perm[i] = i;
        elen[i] = -elen[i];
    }

    CCutil_double_perm_quicksort (perm, elen, ecount);


    olen = CC_SAFE_MALLOC (ecount, double);
    CCcheck_NULL (olen, "out of memory for olen");
    olist = CC_SAFE_MALLOC (2*ecount, int);
    CCcheck_NULL (olist, "out of memory for olist");

    for (i = 0; i < ecount; i++) {
        olen[i] = -elen[perm[i]];
        olist[2*i]   = elist[2*perm[i]];
        olist[2*i+1] = elist[2*perm[i]+1];
    }

    *oncount = ncount;
    *oecount = ecount;
    *polist = olist;
    *polen = olen;
    

CLEANUP:

    CC_IFFREE (perm, int);
    CC_IFFREE (elist, int);
    CC_IFFREE (elen, double);

    return rval;

}

static int build_dp_cut_expand (CCtsp_lpcut_in **cut, int ndomino, int *Acount,
        int **A, int *Bcount, int **B, int handlecount, int *handle,
        CC_SRKexpinfo *expand, int ncount, int *comb)
{
    int rval = 0;
    int *new_Acount = (int *) NULL;
    int *new_Bcount = (int *) NULL;
    int **new_A = (int **) NULL;
    int **new_B = (int **) NULL;
    int new_handlecount = 0;
    int *new_handle = (int *) NULL;
    int i, valid;

    new_Acount = CC_SAFE_MALLOC (ndomino, int);
    CCcheck_NULL (new_Acount, "out of memory for new_Acount");
    new_Bcount = CC_SAFE_MALLOC (ndomino, int);
    CCcheck_NULL (new_Bcount, "out of memory for new_Bcount");
    new_A = CC_SAFE_MALLOC (ndomino, int *);
    CCcheck_NULL (new_A, "out of memory for new_A");
    new_B = CC_SAFE_MALLOC (ndomino, int *);
    CCcheck_NULL (new_B, "out of memory for new_B");

    for (i = 0; i < ndomino; i++) {
        rval = CCcut_SRK_expand (expand, A[i], Acount[i], &(new_A[i]),
                                 &(new_Acount[i]));
        CCcheck_rval (rval, "CCcut_SRK_expand failed");
        rval = CCcut_SRK_expand (expand, B[i], Bcount[i], &(new_B[i]),
                                 &(new_Bcount[i]));
        CCcheck_rval (rval, "CCcut_SRK_expand failed");
    }

    rval = CCcut_SRK_expand (expand, handle, handlecount, &new_handle,
                             &new_handlecount);
    CCcheck_rval (rval, "CCcut_SRK_expand failed");

    rval = check_raw_domino (ncount, new_handlecount, new_handle, ndomino, 
                             new_Acount, new_A, new_Bcount, new_B, &valid,
                             comb);
    CCcheck_rval (rval, "check_raw_domino failed");
    if (!valid) {
        fprintf (stderr, "Expanded domino is not valid\n");
        rval = 1;  goto CLEANUP;
    }

    rval = CCtsp_build_dp_cut (cut, ndomino, new_Acount, new_A, new_Bcount,
              new_B, new_handlecount, new_handle);
    CCcheck_rval (rval, "CCtsp_build_dp_cut failed");
   
CLEANUP:

    CC_IFFREE (new_Acount, int);
    CC_IFFREE (new_Bcount, int);
    if (new_A) {
        for (i = 0; i < ndomino; i++) {
            CC_IFFREE (new_A[i], int);
        }
        CC_IFFREE (new_A, int *);
    }
    if (new_B) {
        for (i = 0; i < ndomino; i++) {
            CC_IFFREE (new_B[i], int);
        }
        CC_IFFREE (new_B, int *);
    }
    CC_IFFREE (new_handle, int);

    return rval;
}

static int remove_light_edges (int ecount, int *elist, double *elen,
        int *new_ecount, int **new_elist, double **new_elen, double tol)
{
    int rval = 0;
    int i, pecount = 0;
    int *pelist = (int *) NULL;
    double *pelen = (double *) NULL;

    for (i = 0; i < ecount; i++) {
        if (elen[i] >= tol) pecount++;
    }

    if (pecount > 0) {
        pelist = CC_SAFE_MALLOC (2*pecount, int);
        CCcheck_NULL (pelist, "out of memory for pelist");
        pelen = CC_SAFE_MALLOC (pecount, double);
        CCcheck_NULL (pelist, "out of memory for pelen");
    }

    pecount = 0;
    for (i = 0; i < ecount; i++) {
        if (elen[i] >= tol) {
            pelist[2*pecount] = elist[2*i];
            pelist[2*pecount + 1] = elist[2*i + 1];
            pelen[pecount] = elen[i];
            pecount++;
        }
    }

CLEANUP:

    *new_ecount = pecount;
    *new_elist = pelist;
    *new_elen = pelen;

    return rval;
}

static void free_raw_dominos (int nIneq, int *nHandle, int **Handle, 
        int *nDominoes, int **nAset, int ***Aset, int **nBset, int ***Bset)
{
    int i, j;

    if (Aset) {
        for (i = 0; i < nIneq; i++) {
            for (j = 0; j < nDominoes[i]; j++) {
                CC_IFFREE (Aset[i][j], int);
            }
            CC_IFFREE (Aset[i], int *);
        }
        CC_IFFREE (Aset, int **);
    }
    if (Bset) {
        for (i = 0; i < nIneq; i++) {
            for (j = 0; j < nDominoes[i]; j++) {
                CC_IFFREE (Bset[i][j], int);
            }
            CC_IFFREE (Bset[i], int *);
        }
        CC_IFFREE (Bset, int **);
    }
    if (nAset) {
        for (i = 0; i < nIneq; i++) {
            CC_IFFREE (nAset[i], int);
        }
        CC_IFFREE (nAset, int *);
    }
    if (nBset) {
        for (i = 0; i < nIneq; i++) {
            CC_IFFREE (nBset[i], int);
        }
        CC_IFFREE (nBset, int *);
    }
    if (Handle) {
        for (i = 0; i < nIneq; i++) {
            CC_IFFREE (Handle[i], int); 
        }
        CC_IFFREE (Handle, int *);
    }
    CC_IFFREE (nHandle, int);
    CC_IFFREE (nDominoes, int);
}

static int check_raw_domino (int ncount, int hcount, int *hand, int dcount, 
        int *Acount, int **A, int *Bcount, int **B, int *valid, int *comb)
{
    int rval = 0;
    int *marks = (int *) NULL;
    int i, j, marker = 0;
    int ain, aout, bin, bout;

    *valid = 0;
    *comb = 0;

    marks = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (marks, "out of memory for marks");
    for (i = 0; i < ncount; i++) marks[i] = 0;

    if (hcount <= 0 || hcount >= ncount) goto CLEANUP;
    marker++;
    for (i = 0; i < hcount; i++) {
        if (hand[i] < 0 || hand[i] >= ncount) goto CLEANUP;
        if (marks[hand[i]] == marker) goto CLEANUP;
        marks[hand[i]] = marker;
    }
    for (i = 0; i < dcount; i++) {
        if (Acount[i] + Bcount[i] >= ncount) goto CLEANUP;
        if (Acount[i] <= 0 || Bcount[i] <= 0) goto CLEANUP;
        marker++;
        for (j = 0; j < Acount[i]; j++) {
            if (marks[A[i][j]] == marker) goto CLEANUP;
            marks[A[i][j]] = marker;
        }
        for (j = 0; j < Bcount[i]; j++) {
            if (marks[B[i][j]] == marker) goto CLEANUP;
            marks[B[i][j]] = marker;
        }
    }
    *valid = 1;

     /* Test for comb */

     /* First check in dominos intersect handle correctly */

    marker++;
    for (i = 0; i < hcount; i++) marks[hand[i]] = marker;
    for (i = 0; i < dcount; i++) {
        ain = aout = 0;
        for (j = 0; j < Acount[i]; j++) {
            if (marks[A[i][j]] == marker) ain++;
            else                          aout++;
        }
        if (ain && aout) goto CLEANUP;

        bin = bout = 0;
        for (j = 0; j < Bcount[i]; j++) {
            if (marks[B[i][j]] == marker) bin++;
            else                          bout++;
        }
        if (bin && bout) goto CLEANUP;
        if ((aout && bout) || (ain && bin)) goto CLEANUP;
    }

    /* Now check that teeth are disjoint */

    marker++;
    for (i = 0; i < dcount; i++) {
        for (j = 0; j < Acount[i]; j++) {
            if (marks[A[i][j]] == marker) goto CLEANUP;
            marks[A[i][j]] = marker;
        }
        for (j = 0; j < Bcount[i]; j++) {
            if (marks[B[i][j]] == marker) goto CLEANUP;
            marks[B[i][j]] = marker;
        }
    }

    *comb = 1;

CLEANUP:

    CC_IFFREE (marks, int);
    return rval;
}

static void  print_raw_domino (int hcount, int *hand, int dcount, 
        int *Acount, int **A, int *Bcount, int **B)
{
    int i, j;

    printf ("Raw Domino\n");
    printf ("    Handle:");
    for (i = 0; i < hcount; i++)  {
        printf (" %d", hand[i]);
        if (i % 15 == 14) printf ("\n           ");
    }
    printf ("\n");
    for (i = 0; i < dcount; i++) {
        printf ("    Domino %d:", i);
        for (j = 0; j < Acount[i]; j++) printf (" %d", A[i][j]); 
        printf ("  |  ");
        for (j = 0; j < Bcount[i]; j++) printf (" %d", B[i][j]); 
        printf ("\n");
    }
    fflush (stdout);
}

#endif
