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
/*  int CCchunk_lift (CCchunk_graph *c, CCchunk_fault *fault,               */
/*      CCchunk_lift_timer *timer, CCchunk_cut_callback *callback)          */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCchunk_ineq_to_lpcut_in (int ncount, int ecount, int *elist,       */
/*      int *ecoef, int rhs, CCtsp_lpcut_in *c)                             */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCchunk_ineq_to_cut (int ncount, int ecount, int *elist,            */
/*      int *ecoef, int rhs, int outside, CCchunk_cut_callback *callback)   */
/*    MISSING                                                               */
/*                                                                          */
/****************************************************************************/

/* THIS MAKES THE ASSUMPTION THAT THERE ARE NO PATHS OF FIXED==1 EDGES OF
   LENGTH > 1 IN THE CHUNK */

#include "machdefs.h"
#include "util.h"
#include "localcut.h"
#include "macrorus.h"

#define EPS 0.01                /* disregard cuts violated by less */

#define KILL_EDGES              /* use edge_dead to fix edges to 0 */

#undef  DEBUG
#undef  DUMPCHUNKS
#undef  DUMPCHUNKSOLS
#undef  DUMP_FACETS

#undef  FIND_ALL_CUTS

#undef  EXTRA_VERIFY

#define MAXEFFORT 400           /* the effort limit for the oracle */

#define ADD_NONFACET            /* when decompose fails, add the cut anyway */

#define LIFT_OVERFLOW (3)

#define MULTLIMIT(x) (((x) == 0) ? (INT_MAX/2) : ((INT_MAX/2)/CC_OURABS(x)))

#define TRISIZE(n) (TRIMAT((n-1),(n-1)))
#define TRIMAT(i,j) (((i)>(j)) ? TRIMAT2((i),(j)) : TRIMAT2((j),(i)))
#define TRIMAT2(i,j) ((i)*((i)-1)/2+(j))

typedef struct paths_info {
    int *deg_sum;
    int *edge_sum;
    int *edge_dead;
    int *pathend;
    int pathcount;
    int ncount;
} paths_info;

typedef struct stripnode {
    int label;
    int done;
    int nadj;
    int *adj;
    int *coef;
} stripnode;


static int
    fix_paths (CCchunk_graph *c),
   *member_dup (int *omem),
    count_unfixed (CCchunk_graph *c),
    load_initial_sols (CCchunk_intmat *mat, CCchunk_graph *ch, CCchunk_ineq *a, int nsols,
        int *sols, paths_info *pinfo),
    strengthen_equality (CCchunk_graph *ch, CCchunk_ineq *a, CCchunk_intmat *mat,
        paths_info *pinfo, CCchunk_lift_timer *timer),
    strengthen_edges (CCchunk_graph *ch, CCchunk_ineq *a, CCchunk_intmat *mat,
        paths_info *pinfo, CCchunk_lift_timer *timer),
    strengthen_work (CCchunk_graph *ch, CCchunk_ineq *a, CCchunk_intmat *mat, CCchunk_ineq *v,
        paths_info *pinfo, CCchunk_lift_timer *timer),
    decompose (CCchunk_graph *ch, CCchunk_ineq *a, CCchunk_intmat *mat, paths_info *pinfo,
        CCchunk_lift_timer *timer, CCchunk_cut_callback *callback, int *finished),
    find_nontight (CCchunk_graph *ch, CCchunk_ineq *a, int *y),
    adjust (int ecount, CCchunk_ineq *a, CCchunk_ineq *c, int *y),
    slack (int ecount, CCchunk_ineq *a, int *y),
    find_orthovec (CCchunk_graph *ch, CCchunk_intmat *mat, CCchunk_ineq *a, CCchunk_ineq *ortho),
    tilt (CCchunk_graph *ch, CCchunk_ineq *a, CCchunk_ineq *cin, int *yin, CCchunk_ineq *c, int *y,
        CCchunk_lift_timer *timer),
    intmat_newpath (CCchunk_intmat *mat, CCchunk_graph *c, int *x),
    paths_build (paths_info *p, int ncount),
    paths_copy (paths_info *p1, paths_info *p2),
#ifdef EXTRA_VERIFY
    verify_cut (CCchunk_graph *ch, CCchunk_ineq *a, CCutil_timer *timer),
#endif
    liberate_fixed (CCchunk_graph *ch, CCchunk_ineq *a, paths_info *pinfo,
        CCchunk_lift_timer *timer, CCchunk_cut_callback *callback, int *finished),
    collect_facet (CCchunk_graph *ch, CCchunk_ineq *a, CCchunk_cut_callback *callback,
        int *finished),
    strip_graph (CCchunk_graph *ch, int nnodes, stripnode *nodes,
        CCchunk_cut_callback *callback, int rhs, int *finished),
    strip_graph_nomembers (int nnodes, stripnode *nodes,
        CCchunk_cut_callback *callback, int rhs, int *finished),
    strip_maxclique (int nnodes, stripnode *nodes, int *maxclique,
        int *work1),
    build_graph (CCchunk_graph *ch, CCchunk_ineq *a, stripnode **pnodes,
                 int **palladj, int **pallcoef);

static void
    liberate_equality (CCchunk_graph *ch, CCchunk_ineq *a, CCchunk_lift_timer *timer),
    flip (int ecount, CCchunk_ineq *a),
    copy_ineq (int ecount, CCchunk_ineq *a, CCchunk_ineq *out),
    copy_tour (int ecount, int *a, int *out),
    scale_down (int ecount, CCchunk_ineq *c),
    paths_init (paths_info *p),
    paths_free (paths_info *p),
    paths_newpath (paths_info *pinfo, CCchunk_graph *c, CCchunk_ineq *a, int *x),
    report_facet (CCchunk_graph *ch, CCchunk_ineq *a),
    subtract_clique (stripnode *nodes, int maxsize, int *maxclique, int label);

static double
    chunk_slack (CCchunk_graph *ch, CCchunk_ineq *a);

static CCchunk_graph
   *chunk_dup (CCchunk_graph *c),
   *chunk_complete (CCchunk_graph *c, CCchunk_ineq *a, CCchunk_ineq *anew);


int CCchunk_lift (CCchunk_graph *c, CCchunk_fault *fault,
        CCchunk_lift_timer *timer, CCchunk_cut_callback *callback)
{
    int rval = 0;
    CCchunk_intmat mat;
    CCchunk_graph *cnew = (CCchunk_graph *) NULL;
    paths_info pinfo;
    CCchunk_ineq a;
    int i;
    int finished = 0;

    CCutil_start_timer (&timer->all);
    a.coef = (int *) NULL;
    
#ifdef DUMPCHUNKS
    printf ("CCchunk_lift\n");
    printf ("%d %d\n", c->ncount, c->ecount);
    for (i=0; i<c->ecount; i++) {
        printf ("%d %d %.16f %d\n",c->end0[i], c->end1[i], c->weight[i],
                c->fixed[i]);
    }
    for (i=0; i<c->ncount; i++) {
        if (i>0) printf (" ");
        printf ("%d",c->equality[i]);
    }
    printf ("\n");
    for (i=0; i<c->ecount; i++) {
        printf ("%d ", fault->a.coef[i]);
    }
    printf ("<= %d\n", fault->a.rhs);
#ifdef DUMPCHUNKSOLS
    printf ("%d\n", fault->nsols);
    for (i=0; i<fault->nsols; i++) {
        int j;
        for (j=0; j<c->ecount; j++) {
            printf ("%d ", fault->sols[i*c->ecount+j]);
        }
        printf ("\n");
    }
#else /* DUMPCHUNKSOLS */
    printf ("0\n");
#endif /* DUMPCHUNKSOLS */
    printf ("viol %.6f\n", -chunk_slack (c, &fault->a));
    fflush (stdout);
#endif /* DUMPCHUNKS */

/*    printf ("Slack %.6f\n", chunk_slack (c, &fault->a));*/

    CCchunk_intmat_init (&mat);
    paths_init (&pinfo);

    if (chunk_slack (c, &fault->a) >= -EPS) {
#ifdef DEBUG
        printf ("No violation - nothing to do\n");
#endif
        rval = 0; goto CLEANUP;
    }

    a.coef = CC_SAFE_MALLOC (c->ecount, int);
    if (a.coef == (int *) NULL) {
        fprintf (stderr, "Out of memory in CCchunk_lift\n");
        rval = 1; goto CLEANUP;
    }
    for (i=0; i<c->ecount; i++) {
        a.coef[i] = fault->a.coef[i];
    }
    a.rhs = fault->a.rhs;

    rval = paths_build (&pinfo, c->ncount);
    if (rval) {
        fprintf (stderr, "paths_build failed\n");
        goto CLEANUP;
    }

    /* make a copy of the chunk since fix_paths, liberate_equality,
       and liberate_fixed modify the chunk */
    cnew = chunk_dup (c);
    if (!cnew) {
        fprintf (stderr, "chunk_dup failed\n");
        rval = -1;
        goto CLEANUP;
    }

    rval = fix_paths (cnew);
    if (rval) {
        fprintf (stderr, "fix_paths failed\n");
        rval = -1;
        goto CLEANUP;
    }

    rval = load_initial_sols (&mat, cnew, &a, fault->nsols,
                              fault->sols, &pinfo);
    if (rval) {
        fprintf (stderr, "load_initial_sols failed\n");
        goto CLEANUP;
    }

    liberate_equality (cnew, &a, timer);

#ifdef DEBUG
    printf ("Equality liberated to:");
    for (i=0; i<cnew->ecount; i++) {
        printf (" %d", a.coef[i]);
    }
    printf (" <= %d (viol %f)\n", a.rhs, -chunk_slack (cnew, &a));
    fflush (stdout);
#endif /* DEBUG */

    if (chunk_slack (cnew, &a) >= -EPS) {
#ifdef DEBUG
        printf ("violation vanished\n");
        fflush (stdout);
#endif
        rval = 0; goto CLEANUP;
    }

    rval = strengthen_edges (cnew, &a, &mat, &pinfo, timer);
    if (rval) {
        fprintf (stderr, "strengthen_edges failed\n");
        goto CLEANUP;
    }

#ifdef DEBUG
    printf ("Edge strengthened to:");
    for (i=0; i<cnew->ecount; i++) {
        printf (" %d", a.coef[i]);
    }
    printf (" <= %d (viol %f)\n", a.rhs, -chunk_slack (cnew, &a));
    fflush (stdout);
#endif

    rval = strengthen_equality (cnew, &a, &mat, &pinfo, timer);
    if (rval) {
        fprintf (stderr, "strengthen_equality failed\n");
        goto CLEANUP;
    }

#ifdef DEBUG
    printf ("Equality strengthened to:");
    for (i=0; i<cnew->ecount; i++) {
        printf (" %d", a.coef[i]);
    }
    printf (" <= %d (viol %f)\n", a.rhs, -chunk_slack (cnew, &a));
    fflush (stdout);
#endif

    if (chunk_slack (cnew, &a) >= -EPS) {
#ifdef DEBUG
        printf ("violation vanished\n");
        fflush (stdout);
#endif
        rval = 0; goto CLEANUP;
    }

    CCutil_start_timer (&timer->decompose);
    rval = decompose (cnew, &a, &mat, &pinfo, timer, callback, &finished);
    CCutil_stop_timer (&timer->decompose, 0);
    if (rval) {
        fprintf (stderr, "decompose failed\n");
        goto CLEANUP;
    }

    rval = 0;

  CLEANUP:
    CCutil_stop_timer (&timer->all, 0);
    CCchunk_intmat_free (&mat);
    paths_free (&pinfo);
    if (cnew) CCchunk_graph_free (cnew);
    CC_IFFREE (a.coef, int);
    return rval;
}

static int fix_paths (CCchunk_graph *c)
{
    int nfixed = 0;
    int *deg = (int *) NULL;
    int i;
    int ecount = c->ecount;
    int ncount = c->ncount;
    int *fixed = c->fixed;
    int *end0 = c->end0;
    int *end1 = c->end1;
    int rval;

    deg = CC_SAFE_MALLOC (ncount, int);
    if (deg == (int *) NULL) {
        fprintf (stderr, "Out of memory\n");
        rval = 1; goto CLEANUP;
    }

    for (i=0; i<ncount; i++) {
        deg[i] = 0;
    }

    for (i=0; i<ecount; i++) {
        if (fixed[i] == 1) {
            deg[end0[i]]++;
            deg[end1[i]]++;
        }
    }

    for (i=0; i<ecount; i++) {
        if (fixed[i] == 1 && deg[end0[i]] == 2 && deg[end1[i]] == 2) {
            fixed[i] = -1;
            deg[end0[i]]--;
            deg[end1[i]]--;
            nfixed++;
        }
    }

    for (i=0; i<ecount; i++) {
        if (fixed[i] == 1 && (deg[end0[i]] == 2 || deg[end1[i]] == 2)) {
            fixed[i] = -1;
            deg[end0[i]]--;
            deg[end1[i]]--;
            nfixed++;
        }
    }

    for (i=0; i<ncount; i++) {
        if (deg[i] == 2) {
            fprintf (stderr, "COULDN'T FIX PATHS\n");
            rval = 1; goto CLEANUP;
        }
    }

#ifdef DEBUG
    printf ("%d edges unfixed to break paths\n", nfixed);
#endif
    rval = 0;

  CLEANUP:
    CC_IFFREE (deg, int);
    return rval;
}

static CCchunk_graph *chunk_complete (CCchunk_graph *c, CCchunk_ineq *a, CCchunk_ineq *anew)
{
    int i, j, k;
    int ecount = TRISIZE(c->ncount);
    CCchunk_graph *cnew = (CCchunk_graph *) NULL;

    CC_IFFREE (anew->coef, int);

    cnew = CCchunk_graph_alloc (c->ncount, ecount);
    if (!cnew) {
        fprintf (stderr, "CCchunk_graph_alloc failed\n");
        goto CLEANUP;
    }

    anew->coef = CC_SAFE_MALLOC (ecount, int);
    if (!anew->coef) {
        fprintf (stderr, "Out of memory\n");
        goto CLEANUP;
    }

    for (i=0; i<cnew->ncount; i++) {
        cnew->equality[i] = c->equality[i];
        cnew->members[i] = member_dup (c->members[i]);
        if (!cnew->members[i]) {
            fprintf (stderr, "member_dup failed\n");
            goto CLEANUP;
        }
        for (j=0; j<i; j++) {
            k = TRIMAT (i,j);
            cnew->end0[k] = j;
            cnew->end1[k] = i;
            cnew->fixed[k] = 0;
            cnew->weight[k] = 0.0;
            anew->coef[k] = 0;
        }
    }

    for (i=0; i<c->ecount; i++) {
        k = TRIMAT (c->end0[i], c->end1[i]);
        cnew->fixed[k] = c->fixed[i];
        cnew->weight[k] = c->weight[i];
        anew->coef[k] = a->coef[i];
    }
    anew->rhs = a->rhs;

    return cnew;

  CLEANUP:
    if (cnew) CCchunk_graph_free (cnew);
    CC_IFFREE (a->coef, int);
    return (CCchunk_graph *) NULL;
}

static CCchunk_graph *chunk_dup (CCchunk_graph *c)
{
    int i;
    CCchunk_graph *cnew = (CCchunk_graph *) NULL;

    cnew = CCchunk_graph_alloc (c->ncount, c->ecount);
    if (!cnew) {
        fprintf (stderr, "CCchunk_graph_alloc failed\n");
        goto CLEANUP;
    }

    for (i=0; i<cnew->ncount; i++) {
        cnew->equality[i] = c->equality[i];
        cnew->members[i] = member_dup (c->members[i]);
        if (!cnew->members[i]) {
            fprintf (stderr, "member_dup failed\n");
            goto CLEANUP;
        }
    }

    for (i=0; i<c->ecount; i++) {
        cnew->end0[i] = c->end0[i];
        cnew->end1[i] = c->end1[i];
        cnew->fixed[i] = c->fixed[i];
        cnew->weight[i] = c->weight[i];
    }

    return cnew;

  CLEANUP:
    if (cnew) CCchunk_graph_free (cnew);
    return (CCchunk_graph *) NULL;
}

static int *member_dup (int *omem)
{
    int i;
    int *nmem = (int *) NULL;
    int cnt;

    for (i=0, cnt=0; omem[i] != -1; i++) cnt++;
    cnt++;

    nmem = CC_SAFE_MALLOC (cnt, int);
    if (!nmem) {
        fprintf (stderr, "Out of memory in member_dup\n");
        return (int *) NULL;
    }

    for (i=0; omem[i] != -1; i++) nmem[i] = omem[i];
    nmem[i] = -1;
    return nmem;
}

static int count_unfixed (CCchunk_graph *c)
{
    int i;
    int cnt = 0;

    for (i=0; i<c->ecount; i++) {
        if (c->fixed[i] == -1) cnt++;
    }
    return cnt;
}

static int load_initial_sols (CCchunk_intmat *mat, CCchunk_graph *ch, CCchunk_ineq *a,
                              int nsols, int *sols, paths_info *pinfo)
{
    int unfixed = count_unfixed (ch);
    int ecount = ch->ecount;
    int rval;
    int i;

    rval = CCchunk_intmat_build (mat, unfixed);
    if (rval) {
        fprintf (stderr, "CCchunk_intmat_build failed\n");
        goto CLEANUP;
    }

    for (i=0; i<nsols; i++) {
        if (slack (ecount, a, sols + i*ecount) == 0) {
#ifdef DEBUG
            int j;
            printf ("pilgrim:");
            for (j=0; j<ecount; j++) {
                printf (" %d", sols[i*ecount+j]);
            }
            printf ("\n");
#endif /* DEBUG */
            paths_newpath (pinfo, ch, a, sols + i*ecount);
            rval = intmat_newpath (mat, ch, sols + i*ecount);
            if (rval) {
                fprintf (stderr, "intmat_newpath failed\n");
                goto CLEANUP;
            }
        }
    }
    rval = 0;

  CLEANUP:
    return rval;
}

static void liberate_equality (CCchunk_graph *ch, CCchunk_ineq *a, CCchunk_lift_timer *timer)
{
    int ncount = ch->ncount;
    int ecount = ch->ecount;
    int m = 0;
    int i;

    CCutil_start_timer (&timer->liberate_equality);

    m = 0;
    for (i=0; i<ecount; i++) {
        if (a->coef[i] > 0) m += a->coef[i];
    }
    m -= a->rhs;
    if (m < 0) m = 0;

    for (i=0; i<ecount; i++) {
        if (ch->fixed[i] == -1) {
            if (ch->equality[ch->end0[i]]) {
                a->coef[i] += m;
            }
            if (ch->equality[ch->end1[i]]) {
                a->coef[i] += m;
            }
        } else if (ch->fixed[i] == 1) {
            /* the += 2*m below overcounts these stars */
            if (ch->equality[ch->end0[i]]) {
                a->rhs -= m;
            }
            if (ch->equality[ch->end1[i]]) {
                a->rhs -= m;
            }
        }
    }
    for (i=0; i<ncount; i++) {
        if (ch->equality[i]) {
            a->rhs += 2*m;
            ch->equality[i] = 0;
        }
    }
    CCutil_stop_timer (&timer->liberate_equality, 0);
}

static int strengthen_equality (CCchunk_graph *ch, CCchunk_ineq *a, CCchunk_intmat *mat,
        paths_info *pinfo, CCchunk_lift_timer *timer)
{
    CCchunk_ineq v;
    int ecount = ch->ecount;
    int ncount = ch->ncount;
    int i;
    int j;
    int *deg_sum = pinfo->deg_sum;
    int rval = 0;

    CCutil_start_timer (&timer->strengthen_equality);

    v.coef = (int *) NULL;

    v.coef = CC_SAFE_MALLOC (ecount, int);
    if (v.coef == (int *) NULL) {
        fprintf (stderr, "Out of memory in strengthen_equality\n");
        rval = -1;
        goto CLEANUP;
    }
    for (i=0; i<ecount; i++) {
        v.coef[i] = 0;
    }
    v.rhs = 0;

    for (i=0; i<ncount; i++) {
        if (deg_sum[i] == 2 * pinfo->pathcount) {
            v.rhs = -2;
            for (j=0; j<ecount; j++) {
                v.coef[j] = 0;
                if (ch->end0[j] == i || ch->end1[j] == i) {
                    if (ch->fixed[j] == 1) {
                        v.rhs += 1;
                    } else if (ch->fixed[j] == -1) {
                        v.coef[j] = -1;
                    }
                }
            }
            rval = strengthen_work (ch, a, mat, &v, pinfo, timer);
            if (rval) {
                fprintf (stderr, "strengthen_work failed\n");
                goto CLEANUP;
            }
        }
    }
    rval = 0;

  CLEANUP:
    CCutil_stop_timer (&timer->strengthen_equality, 0);
    CC_IFFREE (v.coef, int);
    return rval;
}

static int strengthen_edges (CCchunk_graph *ch, CCchunk_ineq *a, CCchunk_intmat *mat,
        paths_info *pinfo, CCchunk_lift_timer *timer)
{
    CCchunk_ineq v;
    int ecount = ch->ecount;
    int i;
    int j;
    int *edge_sum = pinfo->edge_sum;
    int rval = 0;

    CCutil_start_timer (&timer->strengthen_edges);

    v.coef = (int *) NULL;

    v.coef = CC_SAFE_MALLOC (ecount, int);
    if (v.coef == (int *) NULL) {
        fprintf (stderr, "Out of memory in strengthen_edges\n");
        rval = -1;
        goto CLEANUP;
    }
    for (i=0; i<ecount; i++) {
        v.coef[i] = 0;
    }
    v.rhs = 0;

    for (i=0; i<ecount; i++) {
        if (ch->fixed[i] == -1) {
            j = TRIMAT(ch->end0[i],ch->end1[i]);
            if (edge_sum[j] == pinfo->pathcount) {
                v.coef[i] = -1;
                v.rhs = -1;
                rval = strengthen_work (ch, a, mat, &v, pinfo, timer);
                if (rval) {
                    fprintf (stderr, "strengthen_work failed\n");
                    goto CLEANUP;
                }
            }
            if (edge_sum[j] == 0) {
                v.coef[i] = 1;
                v.rhs = 0;
                rval = strengthen_work (ch, a, mat, &v, pinfo, timer);
                if (rval) {
                    fprintf (stderr, "strengthen_work failed\n");
                    goto CLEANUP;
                }
            }
            v.coef[i] = 0;
        }
    }
    rval = 0;

  CLEANUP:
    CCutil_stop_timer (&timer->strengthen_edges, 0);
    CC_IFFREE (v.coef, int);
    return rval;
}

static int strengthen_work (CCchunk_graph *ch, CCchunk_ineq *a, CCchunk_intmat *mat, CCchunk_ineq *v,
        paths_info *pinfo, CCchunk_lift_timer *timer)
{
    int *z = (int *) NULL;
    int *z1 = (int *) NULL;
    CCchunk_ineq a1;
    int ecount = ch->ecount;
    int i;
    int rval = 0;

    CCutil_start_timer (&timer->strengthen_work);

    a1.coef = (int *) NULL;

    z = CC_SAFE_MALLOC (ecount, int);
    z1 = CC_SAFE_MALLOC (ecount, int);
    a1.coef = CC_SAFE_MALLOC (ecount, int);
    if (z == (int *) NULL ||
        z1 == (int *) NULL ||
        a1.coef == (int *) NULL) {
        fprintf (stderr, "Out of memory in strengthen_work\n");
        rval = -1;
        goto CLEANUP;
    }

    for (i=0; i<ecount; i++) z[i] = 0;

    rval = tilt (ch, a, v, z, &a1, z1, timer);
    if (rval) {
        fprintf (stderr, "tilt failed\n");
        goto CLEANUP;
    }

    paths_newpath (pinfo, ch, &a1, z1);
    rval = intmat_newpath (mat, ch, z1);
    if (rval) {
        fprintf (stderr, "intmat_newpath failed\n");
        goto CLEANUP;
    }

    copy_ineq (ecount, &a1, a);

    rval = 0;

  CLEANUP:
    CCutil_stop_timer (&timer->strengthen_work, 0);
    CC_IFFREE (z, int);
    CC_IFFREE (z1, int);
    CC_IFFREE (a1.coef, int);
    return rval;
}

static int decompose (CCchunk_graph *ch, CCchunk_ineq *a, CCchunk_intmat *mat,
        paths_info *pinfo, CCchunk_lift_timer *timer, CCchunk_cut_callback *callback,
        int *finished)
{
    int rval = 0;
    CCchunk_ineq c;
    CCchunk_ineq a1;
    CCchunk_ineq a2;
    CCchunk_ineq itmp;
    int *y = (int *) NULL;
    int *z1 = (int *) NULL;
    int *z2 = (int *) NULL;
    int *ztmp = (int *) NULL;
    int ecount = ch->ecount;
    int nadded = 0;
    double s1, s2, stmp;
    paths_info psave;

    c.coef = (int *) NULL;
    a1.coef = (int *) NULL;
    a2.coef = (int *) NULL;
    paths_init (&psave);

    c.coef = CC_SAFE_MALLOC (ecount, int);
    a1.coef = CC_SAFE_MALLOC (ecount, int);
    a2.coef = CC_SAFE_MALLOC (ecount, int);
    y = CC_SAFE_MALLOC (ecount, int);
    z1 = CC_SAFE_MALLOC (ecount, int);
    z2 = CC_SAFE_MALLOC (ecount, int);
    if (!c.coef || !a1.coef || !a2.coef || !y || !z1 || !z2) {
        fprintf (stderr, "Out of memory in decompose\n");
        rval = -1;
        goto CLEANUP;
    }

    for (;;) {
        rval = find_orthovec (ch, mat, a, &c);
        if (rval == CC_CHUNK_INTMAT_NOORTHO) {
            CCutil_suspend_timer (&timer->decompose);
            rval = liberate_fixed (ch, a, pinfo, timer, callback, finished);
            CCutil_resume_timer (&timer->decompose);
            if (rval) {
                fprintf (stderr, "liberate_fixed failed\n");
            }
            goto CLEANUP;
        } else if (rval) {
            fprintf (stderr, "find_orthovec failed\n");
            goto FAILED_ADD;
        }

        rval = find_nontight (ch, a, y);
        if (rval) {
            fprintf (stderr, "find_nontight failed\n");
            goto CLEANUP;
        }

        rval = adjust (ecount, a, &c, y);
        if (rval) {
            fprintf (stderr, "adjust failed\n");
            goto FAILED_ADD;
        }

#ifdef DEBUG
        printf ("TILTING BY ORTHOGONAL\n");
        fflush (stdout);
#endif

        flip (ecount, &c);

#ifdef DEBUG
        printf ("first tilt\n");
#endif
        rval = tilt (ch, a, &c, y, &a1, z1, timer);
        if (rval) {
            fprintf (stderr, "tilt failed\n");
            goto FAILED_ADD;
        }

        flip (ecount, &c);
        if (slack (ecount, a, z1) == 0) {
#ifdef DEBUG
            printf ("first slack 0\n");
#endif
            paths_newpath (pinfo, ch, a, z1);
            rval = intmat_newpath (mat, ch, z1);
            if (rval) {
                fprintf (stderr, "intmat_newpath failed\n");
                goto CLEANUP;
            }

            nadded++;
        } else {
#ifdef DEBUG
            printf ("second tilt\n");
#endif
            rval = tilt (ch, a, &c, y, &a2, z2, timer);
            if (rval) {
                fprintf (stderr, "tilt failed\n");
                goto FAILED_ADD;
            }
            if (slack (ecount, a, z2) == 0) {
#ifdef DEBUG
                printf ("second slack 0\n");
#endif
                paths_newpath (pinfo, ch, a, z2);
                rval = intmat_newpath (mat, ch, z2);
                if (rval) {
                    fprintf (stderr, "intmat_newpath failed\n");
                    goto CLEANUP;
                }
                nadded++;
            } else {
                s1 = chunk_slack (ch, &a1);
                s2 = chunk_slack (ch, &a2);
                if (s1 > s2) {
                    CC_SWAP (a1, a2, itmp);
                    CC_SWAP (s1, s2, stmp);
                    CC_SWAP (z1, z2, ztmp);
                }
                rval = paths_copy (pinfo, &psave);
                if (rval) {
                    fprintf (stderr, "paths_copy failed\n");
                    goto CLEANUP;
                }

                if (s1 < -EPS) {
                    paths_newpath (pinfo, ch, &a1, z1);
                    rval = intmat_newpath (mat, ch, z1);
                    if (rval) {
                        fprintf (stderr, "intmat_newpath failed\n");
                        goto CLEANUP;
                    }
                    nadded++;
#ifdef DEBUG
                    printf ("recursive decompose 1\n");
#endif
                    rval = decompose (ch, &a1, mat, pinfo, timer, callback,
                                      finished);
                    if (rval == CC_CHUNK_ORACLE_SEARCHLIMITEXCEEDED) {
                        fprintf (stderr, "decompose node limit exceeded\n");
                    } else if (rval) {
                        fprintf (stderr, "decompose failed\n");
                        goto CLEANUP;
                    }
                    CCchunk_intmat_dellastrows (mat, 1);
                    nadded--;
                }

                if (*finished == 0 && s2 < -EPS) {
                    paths_newpath (&psave, ch, &a2, z2);
                    rval = intmat_newpath (mat, ch, z2);
                    if (rval) {
                        fprintf (stderr, "intmat_newpath failed\n");
                        goto CLEANUP;
                    }
                    nadded++;
#ifdef DEBUG
                    printf ("recursive decompose 2\n");
#endif
                    rval = decompose (ch, &a2, mat, &psave, timer, callback,
                                      finished);
                    if (rval == CC_CHUNK_ORACLE_SEARCHLIMITEXCEEDED) {
                        fprintf (stderr, "decompose node limit exceeded\n");
                    } else if (rval) {
                        fprintf (stderr, "decompose failed\n");
                        goto CLEANUP;
                    }
                    CCchunk_intmat_dellastrows (mat, 1);
                    nadded--;
                }

                rval = 0;
                goto CLEANUP;
            }
        }
    }

 FAILED_ADD:
#ifdef ADD_NONFACET
    fprintf (stderr, "ADDING POSSIBLE NON-FACET\n");
    CCutil_suspend_timer (&timer->decompose);
    rval = liberate_fixed (ch, a, pinfo, timer, callback, finished);
    CCutil_resume_timer (&timer->decompose);
    if (rval) {
        fprintf (stderr, "liberate_fixed failed\n");
    }
#endif /* ADD_NONFACET */

  CLEANUP:
    if (nadded) CCchunk_intmat_dellastrows (mat, nadded);
    CC_IFFREE (c.coef, int);
    CC_IFFREE (a1.coef, int);
    CC_IFFREE (a2.coef, int);
    CC_IFFREE (y, int);
    CC_IFFREE (z1, int);
    CC_IFFREE (z2, int);
    paths_free (&psave);

    return rval;
}

static int find_nontight (CCchunk_graph *ch, CCchunk_ineq *a, int *y)
{
/* this assumes that there are no paths of fixed=1 edges in the chunk.
   Thus, setting only 1 nonfixed edge to 1 will not create a subtour,
   nor make a node have degree > 2. */

    int ecount = ch->ecount;
    int i;
    int rhs = a->rhs;

    for (i=0; i<ecount; i++) {
        if (ch->fixed[i] == 1) {
            y[i] = 1;
            rhs -= a->coef[i];
        } else {
            y[i] = 0;
        }
    }
    if (rhs != 0) {
        return 0;
    }

    for (i=0; i<ecount; i++) {
        if (ch->fixed[i] == -1 && a->coef[i] < 0) {
            y[i] = 1;
            return 0;
        }
    }
    fprintf (stderr, "WHOA, unable to find nontight\n");
    return -1;
}

static int adjust (int ecount, CCchunk_ineq *a, CCchunk_ineq *c, int *y)
{
    int lambda = slack (ecount, c, y);
    int mu = slack (ecount, a, y);
    int divv = CCutil_our_gcd (lambda, mu);
    int i;
    int a_lim, c_lim;

#ifdef DEBUG
    printf ("Adjust a (slack %d):", mu);
    for (i=0; i<ecount; i++) printf (" %d", a->coef[i]);
    printf (" <= %d\n", a->rhs);
    printf ("Adjust c (slack %d, div %d):", lambda, divv);
    for (i=0; i<ecount; i++) printf (" %d", c->coef[i]);
    printf (" <= %d\n", c->rhs);
    fflush (stdout);
#endif
    if (divv > 1) {
        lambda /= divv;
        mu /= divv;
    }

    a_lim = MULTLIMIT(lambda);
    c_lim = MULTLIMIT(mu);
    for (i=0; i<ecount; i++) {
        if (CC_OURABS(c->coef[i]) > c_lim || CC_OURABS(a->coef[i]) > a_lim) {
            fprintf (stderr, "overflow in adjust: lambda %d mu %d\n", lambda,
                     mu);
            return LIFT_OVERFLOW;
        }
        c->coef[i] = mu * c->coef[i] - lambda * a->coef[i];
    }
    if (CC_OURABS(c->rhs) > c_lim || CC_OURABS(a->rhs) > a_lim) {
        fprintf (stderr, "overflow in adjust: lambda %d mu %d\n", lambda,
                 mu);
        return LIFT_OVERFLOW;
    }
    c->rhs = mu * c->rhs - lambda * a->rhs;
    scale_down (ecount, c);
    return 0;
}

static int slack (int ecount, CCchunk_ineq *a, int *y)
{
    int i;
    int sl = a->rhs;

    for (i=0; i<ecount; i++) {
        sl -= a->coef[i] * y[i];
    }
    return sl;
}

static double chunk_slack (CCchunk_graph *ch, CCchunk_ineq *a)
{
    int i;
    double sl = a->rhs;

    for (i=0; i<ch->ecount; i++) {
        sl -= a->coef[i] * ch->weight[i];
    }
    return sl;
}

static void flip (int ecount, CCchunk_ineq *a)
{
    int i;

    for (i=0; i<ecount; i++) a->coef[i] = -a->coef[i];
    a->rhs = -a->rhs;
}

static void copy_ineq (int ecount, CCchunk_ineq *a, CCchunk_ineq *out)
{
    int i;

    for (i=0; i<ecount; i++) out->coef[i] = a->coef[i];
    out->rhs = a->rhs;
}

static void copy_tour (int ecount, int *a, int *out)
{
    int i;

    for (i=0; i<ecount; i++) out[i] = a[i];
}

static int find_orthovec (CCchunk_graph *ch, CCchunk_intmat *mat,
        CCchunk_ineq *a, CCchunk_ineq *ortho)
{
    int rval = 0;
    int *amat = (int *) NULL;
    int *orthomat = (int *) NULL;
    int unfixed = count_unfixed (ch);
    int ecount = ch->ecount;
    int i, j;

    amat = CC_SAFE_MALLOC (unfixed, int);
    orthomat = CC_SAFE_MALLOC (unfixed, int);
    if (amat == (int *) NULL || orthomat == (int *) NULL) {
        fprintf (stderr, "Out of memory in find_orthovec\n");
        goto CLEANUP;
    }

    for (i=0, j=0; i<ecount; i++) {
        if (ch->fixed[i] == -1) {
            amat[j++] = a->coef[i];
        }
    }

    rval = CCchunk_intmat_ortho (mat, orthomat, &ortho->rhs, amat);
    if (rval && rval != CC_CHUNK_INTMAT_NOORTHO) {
        fprintf (stderr, "CCchunk_intmat_ortho failed\n");
    }
    if (rval == 0) {
        for (i=0, j=0; i<ecount; i++) {
            if (ch->fixed[i] == -1) {
                ortho->coef[i] = orthomat[j++];
            } else {
                ortho->coef[i] = 0;
            }
        }
    }
  CLEANUP:
    CC_IFFREE (orthomat, int);
    CC_IFFREE (amat, int);
    return rval;
}

static int tilt (CCchunk_graph *ch, CCchunk_ineq *a, CCchunk_ineq *cin,
        int *yin, CCchunk_ineq *c, int *y, CCchunk_lift_timer *timer)
{
    int ecount = ch->ecount;
    int mu;
    int *z = (int *) NULL;
    int rval = 0;

#ifdef DEBUG
    {
        int i;
        printf ("Tilting:");
        for (i=0; i<ecount; i++) printf (" %d", a->coef[i]);
        printf (" <= %d\n", a->rhs);
        printf ("by:");
        for (i=0; i<ecount; i++) printf (" %d", cin->coef[i]);
        printf (" <= %d\n", cin->rhs);
        fflush (stdout);
    }
#endif

    z = CC_SAFE_MALLOC (ecount, int);
    if (!z) {
        fprintf (stderr, "Out of memory in tilt\n");
        return -1;
    }

    copy_tour (ecount, yin, y);

    copy_ineq (ecount, cin, c);

    for (;;) {
        CCutil_start_timer (&timer->tilt_oracle);
        rval = CCchunk_oracle (ch, c, z, (int *) NULL, 1, MAXEFFORT, &timer->oracle);
        CCutil_stop_timer (&timer->tilt_oracle, 0);
        if (rval == CC_CHUNK_ORACLE_INFEASIBLE) {
            rval = 0;
            goto CLEANUP;
        } else if (rval == CC_CHUNK_ORACLE_SEARCHLIMITEXCEEDED) {
            fprintf (stderr, "CCchunk_oracle node limit exceeded\n");
            goto CLEANUP;
        } else if (rval) {
            fprintf (stderr, "CCchunk_oracle failed\n");
            goto CLEANUP;
        }
#ifdef DEBUG
        {
            int i;
            printf ("tsp oracle found:");
            for (i=0; i<ecount; i++) {
                printf (" %d", z[i]);
            }
            printf ("\n");
            fflush (stdout);
        }
#endif
        mu = slack (ecount, a, z);
        if (mu < 0) {
            fprintf (stderr, "ERROR - tilt found point with slack %d\n", mu);
            rval = -1;
            goto CLEANUP;
        } else if (mu == 0) {
#ifdef DEBUG
            printf ("Didn't need to tilt\n");
#endif
            copy_ineq (ecount, a, c);
            copy_tour (ecount, z, y);
            rval = 0;
            goto CLEANUP;
        } else {
            rval = adjust (ecount, a, c, z);
            if (rval) {
                fprintf (stderr, "adjust failed\n");
                goto CLEANUP;
            }
            copy_tour (ecount, z, y);
        }
    }


  CLEANUP:
#ifdef DEBUG
    if (rval == 0) {
        int i;

        printf ("Tilted to:");
        for (i=0; i<ecount; i++) printf (" %d", c->coef[i]);
        printf (" <= %d (viol %f)\n", c->rhs, -chunk_slack (ch, c));
        printf ("Tilt found:");
        for (i=0; i<ecount; i++) {
            printf (" %d", y[i]);
        }
        printf ("\n");
    }
#endif

    CC_IFFREE (z, int);
    return rval;
}

static void scale_down (int ecount, CCchunk_ineq *c)
{
    int i;
    int divv;

    divv = c->rhs;
    for (i = 0; i < ecount; i++) {
        divv = CCutil_our_gcd (divv, c->coef[i]);
    }

    if (divv > 1) {
        for (i = 0; i < ecount; i++) {
            c->coef[i] /= divv;
        }
        c->rhs /= divv;
    }
}

static int intmat_newpath (CCchunk_intmat *mat, CCchunk_graph *c, int *x)
{
    int rval;
    int *xmat = (int *) NULL;
    int ecount = c->ecount;
    int *fixed = c->fixed;
    int i;
    int j;

    xmat = (int *) CC_SAFE_MALLOC (ecount, int);
    if (xmat == (int *) NULL) {
        fprintf (stderr, "Out of memory in intmat_newpath\n");
        rval = 1; goto CLEANUP;
    }

    for (i=0, j=0; i<ecount; i++) {
        if (fixed[i] == -1) {
            xmat[j++] = x[i];
        }
    }

    rval = CCchunk_intmat_addrow (mat, xmat);
    if (rval) {
        fprintf (stderr, "CCchunk_intmat_addrow failed\n");
        goto CLEANUP;
    }

    rval = 0;

  CLEANUP:
    CC_IFFREE (xmat, int);
    return rval;
}

static void paths_init (paths_info *p)
{
    p->deg_sum   = (int *) NULL;
    p->edge_sum  = (int *) NULL;
    p->edge_dead = (int *) NULL;
    p->pathend   = (int *) NULL;
}

static void paths_free (paths_info *p)
{
    CC_IFFREE (p->deg_sum, int);
    CC_IFFREE (p->edge_sum, int);
    CC_IFFREE (p->edge_dead, int);
    CC_IFFREE (p->pathend, int);
}

static int paths_build (paths_info *p, int ncount)
{
    int i;
    int ecount;

    paths_free (p);

    p->ncount = ncount;
    p->pathcount = 0;

    ecount = TRISIZE(ncount);

    p->deg_sum   = CC_SAFE_MALLOC (ncount, int);
    p->edge_sum  = CC_SAFE_MALLOC (ecount, int);
    p->edge_dead = CC_SAFE_MALLOC (ecount, int);
    p->pathend   = CC_SAFE_MALLOC (ncount, int);
    if (p->deg_sum   == (int *) NULL ||
        p->edge_sum  == (int *) NULL ||
        p->edge_dead == (int *) NULL ||
        p->pathend   == (int *) NULL) {
        paths_free (p);
        fprintf (stderr, "Out of memory in paths_build\n");
        return 1;
    }
    for (i=0; i<p->ncount; i++) {
        p->deg_sum[i] = 0;
        p->pathend[i] = 0;
    }
    for (i=0; i<ecount; i++) {
        p->edge_sum[i] = 0;
        p->edge_dead[i] = 0;
    }

    return 0;
}

static int paths_copy (paths_info *p1, paths_info *p2)
{
    int ncount = p1->ncount;
    int ecount = TRISIZE(ncount);
    int i;
    int rval;

    paths_free (p2);

    rval = paths_build (p2, ncount);
    if (rval) {
        fprintf (stderr, "paths_build failed\n");
        goto CLEANUP;
    }

    p2->pathcount = p1->pathcount;

    for (i=0; i<ncount; i++) {
        p2->deg_sum[i] = p1->deg_sum[i];
        p2->pathend[i] = p1->pathend[i];
    }
    for (i=0; i<ecount; i++) {
        p2->edge_sum[i] = p1->edge_sum[i];
        p2->edge_dead[i] = p1->edge_dead[i];
    }

    rval = 0;

  CLEANUP:
    if (rval) {
        paths_free (p2);
    }
    return rval;
}

static void paths_newpath (paths_info *pinfo, CCchunk_graph *c,
        CC_UNUSED CCchunk_ineq *a, int *x)
{
    int i;
    int j;
    int ncount     = pinfo->ncount;
    int *deg_sum   = pinfo->deg_sum;
    int *edge_sum  = pinfo->edge_sum;
    int *edge_dead = pinfo->edge_dead;
    int *pathend   = pinfo->pathend;
    int ecount     = c->ecount;
    int n1, n2;
    int n1end, n2end;

    for (i=0; i<ecount; i++) {
        n1 = c->end0[i];
        n2 = c->end1[i];
        j = TRIMAT(n1,n2);
        if (x[i]) {
            deg_sum[n1]++;
            deg_sum[n2]++;
            edge_sum[j]++;
        }
    }
    for (i=0; i<ncount; i++) {
        pathend[i] = i;
    }

    for (i=0; i<ecount; i++) {
        if (x[i]) {
            n1 = c->end0[i];
            n2 = c->end1[i];
            n1end = pathend[n1];
            n2end = pathend[n2];
            pathend[n1] = -1;
            pathend[n2] = -1;
            pathend[n1end] = n2end;
            pathend[n2end] = n1end;
        }
    }

#ifdef KILL_EDGES
    for (i=0; i<ncount; i++) {
        if (pathend[i] != -1) {
            for (j=i+1; j<ncount; j++) {
                if (pathend[j] != -1 &&
                    j != pathend[i]) {
#ifdef DEBUG
                    if (edge_dead[TRIMAT(i,j)] == 0) {
                        printf ("edge %d-%d is dead\n", i, j);
                    }
#endif
                    edge_dead[TRIMAT(i,j)] = 1;
                }
            }
        }
    }
#endif

    pinfo->pathcount++;
}

static void report_facet (CCchunk_graph *ch, CCchunk_ineq *a)
{
    int i;
    int ecount = ch->ecount;

    printf ("FACET: ");
    for (i=0; i<ecount; i++) {
        printf ("%d ", a->coef[i]);
    }
    printf ("<= %d (viol %f)\n", a->rhs, -chunk_slack(ch, a));
    fflush (stdout);
}

#ifdef EXTRA_VERIFY
static int verify_cut (CCchunk_graph *ch, CCchunk_ineq *a, CCutil_timer *timer)
{
    int i;
    int ecount = ch->ecount;
    int ncount = ch->ncount;
    int k, j;
    int rval;

    for (i=0; i<ecount; i++) {
        if (a->coef[i] < 0) {
            fprintf (stderr, "WHOA, SUPPOSED FACET HAS NEGATIVE COEFFICIENT\n");
            return -1;
        }
    }
    if (a->rhs <= 0) {
        fprintf (stderr, "WHOA, SUPPOSED FACET HAS NONPOSITIVE RHS\n");
        return -1;
    }

    for (i=0, k=0; i<ncount; i++) {
        if (ch->equality[i]) {
            fprintf (stderr, "Verifying cut, equality %d, expected 0\n",
                     ch->equality[i]);
            return -1;
        }
        for (j=0; j<i; j++) {
            if (ch->end0[k] != j || ch->end1[k] != i) {
                fprintf (stderr, "Verifying cut, ends (%d,%d) expected (%d,%d)\n",
                         ch->end0[k], ch->end1[k],i,j);
                return -1;
            }
            if (ch->fixed[k] != -1) {
                fprintf (stderr, "Verifying cut, fixed %d, expected -1\n",
                         ch->fixed[k]);
                return -1;
            }
            k++;
        }
    }
    if (k != ecount) {
        fprintf (stderr, "Verifying cut, ecount %d, expected %d\n", ecount, k);
        return -1;
    }

    CCutil_start_timer (timer);
    rval = CCchunk_verify (ch, a);
    CCutil_stop_timer (timer, 0);
    if (rval) {
        fprintf (stderr, "CCchunk_verify failed\n");
        return rval;
    }

    return 0;
}
#endif /* EXTRA_VERIFY */

static int liberate_fixed (CCchunk_graph *ch, CCchunk_ineq *a, paths_info *pinfo,
        CCchunk_lift_timer *timer, CCchunk_cut_callback *callback, int *finished)
{
    int ecount;
    int rval = 0;
    int i;
    int v0, v1;
    CCchunk_ineq anew;
    CCchunk_graph *cnew = (CCchunk_graph *) NULL;
    int *xsol = (int *) NULL;
    int *fixed;
    int *edge_dead = pinfo->edge_dead;
    int j;

    CCutil_start_timer (&timer->liberate_fixed);

    anew.coef = (int *) NULL;

    cnew = chunk_complete (ch, a, &anew);
    if (cnew == (CCchunk_graph *) NULL) {
        fprintf (stderr, "chunk_complete failed\n");
        rval = -1;
        goto CLEANUP;
    }

    ecount = cnew->ecount;
    fixed = cnew->fixed;

#ifdef KILL_EDGES
    for (i=0; i<ecount; i++) {
        if (anew.coef[i] == 0 && fixed[i] == -1) {
#ifdef DEBUG
            printf ("new edge %d-%d fixed at 0\n", cnew->end0[i],
                    cnew->end1[i]);
#endif
            fixed[i] = 0;
            edge_dead[TRIMAT(cnew->end0[i],cnew->end1[i])] = 1;
        }
    }
#endif

    xsol = CC_SAFE_MALLOC (ecount, int);
    if (xsol == (int *) NULL) {
        fprintf (stderr, "Out of memory in liberate_fixed\n");
        rval = -1;
        goto CLEANUP;
    }

    for (i=0; i<ecount; i++) {
        /* first liberate edges fixed at 1, since they might generate */
        /* additional dead edges */
        if (fixed[i] == 1) {
            v1 = anew.rhs;
            fixed[i] = 0;
            CCutil_start_timer (&timer->liberate_oracle);
            rval = CCchunk_oracle (cnew, &anew, xsol, &v0, 0, MAXEFFORT,
                               &timer->oracle);
            CCutil_stop_timer (&timer->liberate_oracle, 0);
            if (rval) {
/*
                fprintf (stderr, "CCchunk_oracle failed in liberate_fixed\n");
*/
                goto CLEANUP;
            }
#ifdef DEBUG
            printf ("Edge %d (%d-%d) was fixed at 1.  coef %d->%d, rhs %d->%d, tour:\n",
                    i, cnew->end0[i], cnew->end1[i], anew.coef[i],
                    anew.coef[i] + v0 - v1, anew.rhs, v0);
            for (j=0; j<ecount; j++) {
                printf ("%d ", xsol[j]);
            }
            printf ("\n");
            fflush (stdout);
#endif /* DEBUG */
            anew.coef[i] += v0 - v1;
            anew.rhs = v0;
#ifdef KILL_EDGES
            if (anew.coef[i] == 0) {
                fixed[i] = 0;
                xsol[i] = 0;
                edge_dead[TRIMAT(cnew->end0[i],cnew->end1[i])] = 1;
#ifdef DEBUG
                printf ("surprise - fixed edge %d-%d killed\n",
                        cnew->end0[i], cnew->end1[i]);
#endif
            } else {
                fixed[i] = -1;
            }
#else /* KILL_EDGES */
            fixed[i] = -1;
#endif /* KILL_EDGES */
            paths_newpath (pinfo, cnew, &anew, xsol);
            scale_down (ecount, &anew);
        }
    }

    for (i=0; i<ecount; i++) {
        /* now, the non-dead edges fixed at 0 */
        j = TRIMAT(cnew->end0[i],cnew->end1[i]);
#ifdef KILL_EDGES
        if (fixed[i] == 0 && edge_dead[j] == 0) {
#else
        if (fixed[i] == 0) {
#endif
            v0 = anew.rhs;
            fixed[i] = 1;
            CCutil_start_timer (&timer->liberate_oracle);
            rval = CCchunk_oracle (cnew, &anew, xsol, &v1, 0, MAXEFFORT,
                               &timer->oracle);
            CCutil_stop_timer (&timer->liberate_oracle, 0);
            if (rval) {
/*
                fprintf (stderr, "CCchunk_oracle failed in liberate_fixed\n");
*/
                goto CLEANUP;
            }
#ifdef DEBUG
            printf ("Edge %d (%d-%d) was fixed at 0.  coef %d->%d, rhs %d->%d, tour:\n",
                    i, cnew->end0[i], cnew->end1[i], anew.coef[i],
                    anew.coef[i] + v0 - v1, anew.rhs, v0);
            for (j=0; j<ecount; j++) {
                printf ("%d ", xsol[j]);
            }
            printf ("\n");
            fflush (stdout);
#endif /* DEBUG */
            anew.coef[i] += v0 - v1;
            anew.rhs = v0;
#ifdef KILL_EDGES
            if (anew.coef[i] == 0) {
                fixed[i] = 0;
                xsol[i] = 0;
                edge_dead[TRIMAT(cnew->end0[i],cnew->end1[i])] = 1;
#ifdef DEBUG
                printf ("zero edge %d-%d killed\n", cnew->end0[i],
                        cnew->end1[i]);
#endif
            } else {
                fixed[i] = -1;
            }
#else /* KILL_EDGES */
            fixed[i] = -1;
#endif /* KILL_EDGES */
            paths_newpath (pinfo, cnew, &anew, xsol);
            scale_down (ecount, &anew);
        }
#if 0
    }
#endif
    }

#ifdef KILL_EDGES
    for (i=0; i<ecount; i++) {
        j = TRIMAT(cnew->end0[i],cnew->end1[i]);
        if (fixed[i] == 0 && edge_dead[j]) {
            fixed[i] = -1;
        }
    }
#endif


#ifdef EXTRA_VERIFY
    rval = verify_cut (cnew, &anew, &timer->verify_oracle);
    if (rval) {
        fprintf (stderr, "verify_cut failed\n");
        goto CLEANUP;
    }
#endif

    CCutil_suspend_timer (&timer->liberate_fixed);
    rval = collect_facet (cnew, &anew, callback, finished);
    CCutil_resume_timer (&timer->liberate_fixed);
    if (rval) {
        fprintf (stderr, "collect_facet failed\n");
        goto CLEANUP;
    }
    rval = 0;

  CLEANUP:
    CCutil_stop_timer (&timer->liberate_fixed, 0);
    CC_IFFREE (anew.coef, int);
    CCchunk_graph_free (cnew);
    CC_IFFREE (xsol, int);
    return 0;
}

static int collect_facet (CCchunk_graph *ch, CCchunk_ineq *a,
        CCchunk_cut_callback *callback, int *finished)
{
    stripnode *nodes = (stripnode *) NULL;
    int *alladj = (int *) NULL;
    int *allcoef = (int *) NULL;
    int rval = 0;

#ifdef DEBUG
    report_facet (ch, a);
#endif

    rval = build_graph (ch, a, &nodes, &alladj, &allcoef);
    if (rval) {
        fprintf (stderr, "build_graph failed\n");
        goto CLEANUP;
    }

    rval = strip_graph (ch, ch->ncount, nodes, callback, -2*a->rhs, finished);
    if (rval) {
        fprintf (stderr, "strip_graph failed\n");
        goto CLEANUP;
    }

#ifdef DUMP_FACETS
    printf (" (viol %.6f)\n", -chunk_slack (ch, a));
#endif

    rval = 0;

  CLEANUP:
    CC_IFFREE (allcoef, int);
    CC_IFFREE (alladj, int);
    CC_IFFREE (nodes, stripnode);

    return rval;
}

static int lpcut_begin (void *u_data)
{
    CCtsp_lpcut_in *c = (CCtsp_lpcut_in *) u_data;

    c->cliques = CC_SAFE_MALLOC (10, CCtsp_lpclique);
    if (c->cliques == (CCtsp_lpclique *) NULL) {
        return -1;
    }
    return 0;
}

static int lpcut_add_clique (int *arr, int size, void *u_data)
{
    int rval;
    CCtsp_lpcut_in *c = (CCtsp_lpcut_in *) u_data;

    rval = CCutil_reallocrus_count ((void **) &(c->cliques), c->cliquecount+1,
                                    sizeof (c->cliques[0]));
    if (rval) {
        fprintf (stderr, "couldn't realloc cliques\n");
        return rval;
    }
    rval = CCtsp_array_to_lpclique (arr, size, &(c->cliques[c->cliquecount]));
    if (rval) {
        fprintf (stderr, "CCtsp_array_to_lpclique failed\n");
        return rval;
    }
    c->cliquecount++;
    return 0;
}

static int lpcut_abort (CC_UNUSED void *u_data)
{
    return -1;
}

static int lpcut_finish (int rhs, int *finished, void *u_data)
{
    CCtsp_lpcut_in *c = (CCtsp_lpcut_in *) u_data;
    int rval;

    c->rhs = rhs;
    c->sense = 'G';
    c->branch = 0;
    
    rval = CCtsp_construct_skeleton (c, CCtsp_max_node (c));
    if (rval) {
        fprintf (stderr, "CCtsp_construct_skeleton failed\n");
        return rval;
    }

    *finished = 0;
    
    return 0;
}
    
/* ineq is (ecoef * x) >= rhs */
int CCchunk_ineq_to_lpcut_in (int ncount, int ecount, int *elist, int *ecoef,
        int rhs, CCtsp_lpcut_in *c)
{
    int rval;
    
    CCchunk_cut_callback ccb;
    
    CCtsp_init_lpcut_in (c);

    ccb.begin_cut = lpcut_begin;
    ccb.add_clique = lpcut_add_clique;
    ccb.abort_cut = lpcut_abort;
    ccb.finish_cut = lpcut_finish;
    ccb.u_data = (void *) c;

    rval = CCchunk_ineq_to_cut (ncount, ecount, elist, ecoef, rhs, 0, &ccb);
    if (rval) {
        fprintf (stderr, "CCchunk_ineq_to_cut failed\n");
        CCtsp_free_lpcut_in (c);
    }
    return rval;
}

/* ineq is (ecoef * x) >= rhs */
int CCchunk_ineq_to_cut (int ncount, int ecount, int *elist, int *ecoef,
        int rhs, int outside, CCchunk_cut_callback *callback)
{
    int rval;
    stripnode *nodes = (stripnode *) NULL;
    int *alladj = (int *) NULL;
    int *allcoef = (int *) NULL;
    int i;
    int j;
    int k;
    int c;
    int minc;
    int gotmin;
    int xdiv;
    int finished = 0;

    nodes = CC_SAFE_MALLOC (ncount, stripnode);
    alladj = CC_SAFE_MALLOC (ncount*ncount, int);
    allcoef = CC_SAFE_MALLOC (ncount*ncount, int);
    if (nodes == (stripnode *) NULL ||
        alladj == (int *) NULL ||
        allcoef == (int *) NULL) {
        fprintf (stderr, "Out of memory in CCchunk_ineq_to_cut\n");
        rval = 1; goto CLEANUP;
    }
    
    for (i=0; i<ncount; i++) {
        nodes[i].adj = alladj + i*ncount;
        nodes[i].coef = allcoef + i*ncount;
        for (j=0; j<ncount; j++) {
            nodes[i].adj[j] = j;
            nodes[i].coef[j] = 0;
        }
        nodes[i].nadj = ncount;
    }
    
    for (i=0; i<ecount; i++) {
        if (ecoef[i]) {
            j = elist[2*i];
            k = elist[2*i+1];
            nodes[j].coef[k] = -ecoef[i];
            nodes[k].coef[j] = -ecoef[i];
        }
    }
    rhs = -rhs;

    for (i=0; i<ncount; i++) {
        if (i != outside) {
            c = -nodes[i].coef[outside];
            for (j=0; j<ncount; j++) {
                nodes[i].coef[j] += c;
                nodes[j].coef[i] += c;
            }
            rhs += 2*c;
        }
    }

    gotmin = 0;
    minc = 0;
    for (i=0; i<ncount; i++) {
        if (i != outside) {
            for (j=0; j<ncount; j++) {
                if (j != outside && i != j) {
                    if (!gotmin || nodes[i].coef[j] < minc) {
                        minc = nodes[i].coef[j];
                        gotmin = 1;
                    }
                }
            }
        }
    }

    for (i=0; i<ncount; i++) {
        if (i != outside) {
            for (j=0; j<ncount; j++) {
                if (j != outside) {
                    nodes[i].coef[j] -= minc;
                }
            }
        }
    }
    rhs -= (ncount-2) * minc;

    for (i=0; i<ncount; i++) {
        nodes[i].coef[i] = 0;
    }

    for (i=0; i<ncount; i++) {
        nodes[i].nadj = 0;
        for (j=0; j<ncount; j++) {
            if (nodes[i].coef[j]) {
                nodes[i].adj[nodes[i].nadj] = j;
                nodes[i].coef[nodes[i].nadj] = nodes[i].coef[j];
                nodes[i].nadj++;
            }
        }
    }

    xdiv = rhs;
    for (i=0; i<ncount; i++) {
        for (j=0; j<nodes[i].nadj; j++) {
            xdiv = CCutil_our_gcd (xdiv, nodes[i].coef[j]);
        }
    }

    rhs /= xdiv;
    for (i=0; i<ncount; i++) {
        for (j=0; j<nodes[i].nadj; j++) {
            nodes[i].coef[j] /= xdiv;
        }
    }

    rval = strip_graph_nomembers (ncount, nodes, callback, -2*rhs, &finished);
    if (rval) {
        fprintf (stderr, "strip_graph_nomembers failed\n");
        goto CLEANUP;
    }

    rval = 0;

  CLEANUP:
    CC_IFFREE (nodes, stripnode);
    CC_IFFREE (alladj, int);
    CC_IFFREE (allcoef, int);
    return rval;
}

static int strip_graph (CCchunk_graph *ch, int nnodes, stripnode *nodes,
        CCchunk_cut_callback *callback, int rhs, int *finished)
{
    int deg;
    int *maxclique = (int *) NULL;
    int *work1 = (int *) NULL;
    int cliquesize;
    int i;
    int rval = 0;
    int *p = (int *) NULL;

    rval = (*callback->begin_cut) (callback->u_data);
    if (rval) {
        fprintf (stderr, "callback begin_cut failed\n");
        goto CLEANUP;
    }

    maxclique = CC_SAFE_MALLOC (nnodes, int);
    if (!maxclique) goto CLEANUP;
    work1 = CC_SAFE_MALLOC (nnodes, int);
    if (!work1) goto CLEANUP;

#ifdef DUMP_FACETS
    printf ("FACET: ");
#endif

    while (1) {
        deg = strip_maxclique (nnodes, nodes, maxclique, work1);
        if (deg <= 1) {
            CC_FREE (work1, int);
            CC_FREE (maxclique, int);
#ifdef DUMP_FACETS
            printf (" >= %d", rhs);
#endif
            rval = (*callback->finish_cut) (rhs, finished, callback->u_data);
            if (rval) {
                fprintf (stderr, "callback finish_cut failed\n");
                goto CLEANUP;
            }
            return 0;
        }
#ifdef  DUMP_FACETS
        printf ("(");
        for (i=0; i<deg; i++) {
            printf ("%d",maxclique[i]);
            if (i+1<deg) printf (" ");
        }
        printf (")");
#endif
        for (i=0, cliquesize=0; i<deg; i++) {
            for (p = ch->members[maxclique[i]]; *p != -1; p++) {
                work1[cliquesize++] = *p;
            }
        }
        rhs += 2*deg;
        rval = (*callback->add_clique) (work1, cliquesize, callback->u_data);
        if (rval) {
            fprintf (stderr, "callback add_clique failed\n");
            goto CLEANUP;
        }
    }

  CLEANUP:
    CC_IFFREE (work1, int);
    CC_IFFREE (maxclique, int);
    (*callback->abort_cut) (callback->u_data);
    return rval;
}

static int strip_graph_nomembers (int nnodes, stripnode *nodes,
        CCchunk_cut_callback *callback, int rhs, int *finished)
{
    int deg;
    int *maxclique = (int *) NULL;
    int *work1 = (int *) NULL;
    int rval = 0;

    rval = (*callback->begin_cut) (callback->u_data);
    if (rval) {
        fprintf (stderr, "callback begin_cut failed\n");
        goto FAILURE;
    }

    maxclique = CC_SAFE_MALLOC (nnodes, int);
    if (!maxclique) goto FAILURE;
    work1 = CC_SAFE_MALLOC (nnodes, int);
    if (!work1) goto FAILURE;

    while (1) {
        deg = strip_maxclique (nnodes, nodes, maxclique, work1);
        if (deg <= 1) {
            CC_FREE (work1, int);
            CC_FREE (maxclique, int);
            rval = (*callback->finish_cut) (rhs, finished, callback->u_data);
            if (rval) {
                fprintf (stderr, "callback finish_cut failed\n");
                goto FAILURE;
            }
            return 0;
        }
        rhs += 2*deg;
        rval = (*callback->add_clique) (maxclique, deg, callback->u_data);
        if (rval) {
            fprintf (stderr, "callback add_clique failed\n");
            goto FAILURE;
        }
    }

  FAILURE:
    CC_IFFREE (work1, int);
    CC_IFFREE (maxclique, int);
    (*callback->abort_cut) (callback->u_data);
    return 1;
}

static int strip_maxclique (int nnodes, stripnode *nodes, int *maxclique,
                            int *work1)
{
    int simpsize = 0;
    int nonsize = 0;
    int i, ia, j, ja, k;
    int label = 0;
    int mindeg;
    int size;

    for (i=0; i<nnodes; i++) {
        nodes[i].label = 0;
        nodes[i].done = 0;
    }
    for (i=0; i<nnodes; i++) {
        if (!nodes[i].done) {
            mindeg = nodes[i].nadj;
            size = 1;
            label++;
            nodes[i].label = label;
            for (ia=0; ia<nodes[i].nadj; ia++) {
                j = nodes[i].adj[ia];
                nodes[j].label = label;
            }
            for (ia=0; ia<nodes[i].nadj; ia++) {
                j = nodes[i].adj[ia];
                if (nodes[j].label == label) {
                    size++;
                    if (nodes[j].nadj < mindeg) mindeg = nodes[j].nadj;
                    label++;
                    nodes[j].label = label;
                    for (ja=0; ja<nodes[j].nadj; ja++) {
                        k = nodes[j].adj[ja];
                        if (nodes[k].label == label-1) {
                            nodes[k].label = label;
                        }
                    }
                }
            }
            if (mindeg == size-1 && size > simpsize) {
                simpsize = size;
                maxclique[0] = i;
                for (ia=0, k=1; ia<nodes[i].nadj; ia++) {
                    j = nodes[i].adj[ia];
                    if (nodes[j].label == label) {
                        maxclique[k++] = j;
                    }
                }
                if (size != k) {
                    fprintf (stderr, "Size mismatch 1\n");
                    simpsize = k;
                }
            } else if (mindeg > size-1 && size > nonsize) {
                nonsize = size;
                work1[0] = i;
                for (ia=0, k=1; ia<nodes[i].nadj; ia++) {
                    j = nodes[i].adj[ia];
                    if (nodes[j].label == label) {
                        work1[k++] = j;
                    }
                }
                if (size != k) {
                    fprintf (stderr, "Size mismatch 2\n");
                    nonsize = k;
                }
                for (ia=0; ia<nodes[i].nadj; ia++) {
                    j = nodes[i].adj[ia];
                    if (nodes[j].label == label) {
                        nodes[j].done = 1;
                    }
                }
            }
        }
    }
    if (simpsize <= 1 && nonsize > 1) {
        for (i=0; i<nonsize; i++) {
            maxclique[i] = work1[i];
        }
        simpsize = nonsize;
    }

    label++;
    subtract_clique (nodes, simpsize, maxclique, label);

    return simpsize;
}

static void subtract_clique (stripnode *nodes, int maxsize, int *maxclique,
                             int label)
{
    int i, j;
    int k;

    for (i=0; i<maxsize; i++) {
        nodes[maxclique[i]].label = label;
    }

    for (i=0; i<maxsize; i++) {
        j = maxclique[i];
        k = 0;
        while (k < nodes[j].nadj) {
            if (nodes[nodes[j].adj[k]].label == label) {
                nodes[j].coef[k]--;
            }
            if (nodes[j].coef[k]) {
                k++;
            } else {
                nodes[j].nadj--;
                nodes[j].adj[k] = nodes[j].adj[nodes[j].nadj];
                nodes[j].coef[k] = nodes[j].coef[nodes[j].nadj];
            }
        }
    }
}



static int build_graph (CCchunk_graph *ch, CCchunk_ineq *a, stripnode **pnodes,
                        int **palladj, int **pallcoef)
{
    int i, j, k;
    int ncoef;
    stripnode *nodes = (stripnode *) NULL;
    int *alladj = (int *) NULL;
    int *allcoef = (int *) NULL;

    nodes = CC_SAFE_MALLOC (ch->ncount, stripnode);
    if (!nodes) {
        return -1;
    }
    for (i=0; i<ch->ncount; i++) {
        nodes[i].nadj = 0;
    }

    ncoef = 0;
    for (i=0; i<ch->ecount; i++) {
        if (a->coef[i]) {
            nodes[ch->end0[i]].nadj++;
            nodes[ch->end1[i]].nadj++;
            ncoef++;
        }
    }
    alladj = CC_SAFE_MALLOC (2*ncoef, int);
    if (!alladj) {
        CC_FREE (nodes, stripnode);
        return -1;
    }
    allcoef = CC_SAFE_MALLOC (2*ncoef, int);
    if (!allcoef) {
        CC_FREE (alladj, int);
        CC_FREE (nodes, stripnode);
        return -1;
    }
    ncoef = 0;
    for (i=0; i<ch->ncount; i++) {
        nodes[i].adj = alladj + ncoef;
        nodes[i].coef = allcoef + ncoef;
        ncoef += nodes[i].nadj;
        nodes[i].nadj = 0;
    }
    for (i=0; i<ch->ecount; i++) {
        if (a->coef[i]) {
            j = ch->end0[i];
            k = ch->end1[i];
            nodes[j].adj[nodes[j].nadj] = k;
            nodes[j].coef[nodes[j].nadj] = a->coef[i];
            nodes[j].nadj++;
            nodes[k].adj[nodes[k].nadj] = j;
            nodes[k].coef[nodes[k].nadj] = a->coef[i];
            nodes[k].nadj++;
        }
    }

    *pnodes = nodes;
    *palladj = alladj;
    *pallcoef = allcoef;

    return 0;
}

static void print_vector (int ecount, int v[])
{
    int i;

    for (i=0; i<ecount; i++) {
        printf ("%d ", v[i]);
    }
    printf ("\n");
    fflush (stdout);
}

static void print_ineq (int ecount, CCchunk_ineq *a)
{
    int i;

    for (i=0; i<ecount; i++) {
        printf ("%d ", a->coef[i]);
    }
    printf ("<= %d\n", a->rhs);
    fflush (stdout);
}
