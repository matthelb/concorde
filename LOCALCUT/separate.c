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
/*  int CCchunk_lpinit (CCchunklp **lp_p, const char *lp_name,              */
/*      int lp_nrows, double *xstar)                                        */
/*    Creates the *lp_p object (see 'typedef struct CCchunklp' below)       */
/*    and initializes it to have name lp_name, N empty, and the initial     */
/*    (empty) lp from above.                                                */
/*                                                                          */
/*  void CCchunk_lpfree (CCchunklp **lp_p)                                  */
/*    Destroys *lp_p, freeing all its parts.                                */
/*                                                                          */
/*  int CCchunk_lpaddcol (CCchunklp *lp, double *x)                         */
/*    Adds a new entry to N, and updates the lp to reflect the new x^j.     */
/*    Note that it is possible some rows of the LP have been relaxed;       */
/*    however, x is still expected to have lp_nrows entries; those for      */
/*    the relaxed rows are skipped when adding x.                           */
/*                                                                          */
/*  int CCchunk_lprelaxrow (CCchunklp *lp, int del_row)                     */
/*    Deletes row del_row+1 from lp -- that is, changes lp so that          */
/*    coordinate del_row of x^* is ignored.  This makes easier for the      */
/*    lp to be feasible.                                                    */
/*                                                                          */
/*  int CCchunk_lpsolve (CCchunklp *lp, int *lpstatus_p, double *c,         */
/*      double *alpha_p)                                                    */
/*    "Solves" lp: *lp_status_p and (c,*alpha_p) will be set as             */
/*    described above, depending on the settings of SEPARATE_OPTIMIZE,      */
/*    SEPARATE_OPTIMIZE, SEPARATE_NORML1, and SEPARATE_EPSILON.             */
/*                                                                          */
/*  int CCchunk_lpbasis (CCchunklp *lp, int ncols, int *basis)              */
/*    Sets basis[i] to 1 if column i is in the currentbasis, 0              */
/*    otherwise.                                                            */
/*                                                                          */
/****************************************************************************/

/*

The purpose of the following is to solve feasibility problems
of the form:

(0)   exists?       lambda_j
      subject to    sum_{j \in N} lambda_j     = 1
                    sum_{j \in N} lambda_j x^j = x^*
                    lambda_j >= 0 (j \in N)

where {x^j \in {0,1}^m: j \in N} is a finite set of incidence
vectors of tours indexed on some edge-set.

To solve this, there are several linear programs that can be used,
controlled by the #defines below.  In all cases, if the LP used is
feasible, then CCchunk_lpsolve sets *lpstatus_p = CC_CHUNK_LPFEASIBLE.  If
the LP used is infeasible, then CCchunk_lpsolve sets *lpstatus_p =
CC_CHUNK_LPINFEASIBLE and (c,*alpha_p) will satisfy

          c' x^j <= (*alpha_p)  for j \in N
          c' x^* >  (*alpha_p)

If SEPARATE_OPTIMIZE is defined, when *lpstatus_p == CC_CHUNK_LPINFEASIBLE,
(c,*alpha_p) will satisfy the stronger

          c' x^* >  (*alpha_p) + SEPARATE_EPSILON

and will maximize the slack of this constraint subject to the other
constraints.  In some cases, when N is empty, this maximization is
unbounded.  In this case, c' == 1 and (*alpha_p) == -1 is returned.

The exact meaning of feasible, and further properties of (c,*alpha_p)
depend on the variant used.

----------

If !defined(SEPARATE_OPTIMIZE), then the lp solved is:

(1)   minimize      0
      subject to    sum_{j \in N} lambda_j     = 1
                    sum_{j \in N} lambda_j x^j = x^*
                    lambda_j >= 0 (j \in N)

When the lp is infeasible, we also have

    c'x^* - (*alpha_p) = (sum_{j in \N, lambda_j < 0} -lambda_j) 
                         + (sum_{i=1 to m} |delta_i|),

for some vectors lambda and delta such that

        sum_{j \in N} lambda_j     +                delta_0     = 1
        sum_{j \in N} lambda_j x^j + sum_{i=1 to m} delta_i e_i = x^*

----------

If defined(SEPARATE_OPTIMIZE) && !defined(SEPARATE_STRICTAFFINE) &&
!defined(SEPARATE_NORML1), then the lp solved is:

(2)   minimize      e'u + e'v + s + t
      subject to    sum_{j \in N} lambda_j     + s - t = 1
                    sum_{j \in N} lambda_j x^j + u - v = x^*
                    lambda_j >= 0 (j \in N)
                    u >= 0, v >= 0, s >= 0, t >= 0

(u, v are vectors, s, t are scalars).  The problem is considered
feasible if e'u + e'v + s + t <= SEPARATE_EPSILON.  Otherwise
(c,*alpha_p) satisfy:

                    c' x^* > (*alpha_p) + SEPARATE_EPSILON
                    c' x^j <= (*alpha_p) for j \in N
                    -1 <= c' <= 1
                    -1 <= (*alpha_p) <= 1

----------

If defined(SEPARATE_OPTIMIZE) && defined(SEPARATE_STRICTAFFINE) &&
!defined(SEPARATE_NORML1), then the lp solved is:

(3)   minimize      e'u + e'v
      subject to    sum_{j \in N} lambda_j             = 1
                    sum_{j \in N} lambda_j x^j + u - v = x^*
                    lambda_j >= 0 (j \in N)
                    u >= 0, v >= 0

The problem is considered feasible if e'u + e'v <= SEPARATE_EPSILON.
Otherwise (c,*alpha_p) satisfy:

                    c' x^* > (*alpha_p) + SEPARATE_EPSILON
                    c' x^j <= (*alpha_p) for j \in N
                    -1 <= c' <= 1

----------

If defined(SEPARATE_OPTIMIZE) && !defined(SEPARATE_STRICTAFFINE) &&
defined(SEPARATE_NORML1), then the lp solved is:

(4)   minimize      -t
      subject to    sum_{j \in N} lambda_j     + s - t = 0
                    sum_{j \in N} lambda_j x^j + w - t x^* = 0
                    lambda_j >= 0 (j \in N)
                    -1 <= w <= 1
                    -1 <= s <= 1
                     0 <= t <= 2/SEPARATE_EPSILON

(w is vector, s, t are scalar).  The problem is considered feasible if
1/t <= SEPARATE_EPSILON.  Otherwise (c,*alpha_p) satisfy:

                    c' x^* > (*alpha_p) + SEPARATE_EPSILON
                    c' x^j <= (*alpha_p) for j \in N
                    ||c||_1 + |(*alpha_p)| = 1  

Note that SEPARATE_NORML1 does not work with SEPARATE_EPSILON == 0.0

----------

If defined(SEPARATE_OPTIMIZE) && defined(SEPARATE_STRICTAFFINE) &&
defined(SEPARATE_NORML1), then the lp solved is:

(5)   minimize      -t
      subject to    sum_{j \in N} lambda_j         - t = 0
                    sum_{j \in N} lambda_j x^j + w - t x^* = 0
                    lambda_j >= 0 (j \in N)
                    -1 <= w <= 1
                     0 <= t <= 2/SEPARATE_EPSILON

(w is vector, t is scalar).  The problem is considered feasible if
1/t <= SEPARATE_EPSILON.  Otherwise (c,*alpha_p) satisfy:

                    c' x^* > (*alpha_p) + SEPARATE_EPSILON
                    c' x^j <= (*alpha_p) for j \in N
                    ||c||_1 = 1  

Note that SEPARATE_NORML1 does not work with SEPARATE_EPSILON == 0.0

----------

The following 5 functions deal with the above problem class:

  int CCchunk_lpinit (CCchunklp **lp_p, const char *lp_name, int lp_nrows,
                      double *xstar)

     Creates the *lp_p object (see 'typedef struct CCchunklp' below)
     and initializes it to have name lp_name, N empty, and the initial
     (empty) lp from above.

  void CCchunk_lpfree (CCchunklp **lp_p)

     Destroys *lp_p, freeing all its parts.

  int CCchunk_lpaddcol (CCchunklp *lp, double *x)

     Adds a new entry to N, and updates the lp to reflect the new x^j.
     Note that it is possible some rows of the LP have been relaxed;
     however, x is still expected to have lp_nrows entries; those for
     the relaxed rows are skipped when adding x.

  int CCchunk_lprelaxrow (CCchunklp *lp, int del_row)

     Deletes row del_row+1 from lp -- that is, changes lp so that
     coordinate del_row of x^* is ignored.  This makes easier for the
     lp to be feasible.

  int CCchunk_lpsolve (CCchunklp *lp, int *lpstatus_p, double *c,
                       double *alpha_p)

     "Solves" lp: *lp_status_p and (c,*alpha_p) will be set as
     described above, depending on the settings of SEPARATE_OPTIMIZE,
     SEPARATE_OPTIMIZE, SEPARATE_NORML1, and SEPARATE_EPSILON.

*/

#include "machdefs.h"
#include "util.h"
#include "lp.h"
#include "localcut.h"

/* if SEPARATE_OPTIMIZE is defined, and *lpstatus_p == CC_CHUNK_LPINFEASIBLE,
   the returned c and *alpha_p maximize sum_{i=0 to m}... */
#define SEPARATE_OPTIMIZE

/* if SEPARATE_STRICTAFFINE is defined, then the affine hull constraint's
   violation will be constrained to zero, rather than included in the
   optimization.*/
#define SEPARATE_STRICTAFFINE

/* if SEPARATE_NORML1 is defined, then feasibility will be based on
   getting close in the L_1 norm.  Otherwise feasibility is defined
   based on getting close in the L_infinity norm. */
#define SEPARATE_NORML1

/* SEPARATE_EPSILON defines "close" for feasibility.  Only used if
   SEPARATE_OPTIMIZE is defined. */
#define SEPARATE_EPSILON 0.01


#ifndef SEPARATE_OPTIMIZE
#ifdef SEPARATE_STRICTAFFINE
#undef SEPARATE_STRICTAFFINE
#endif
#ifdef SEPARATE_NORML1
#undef SEPARATE_NORML1
#endif
#endif

int CCchunk_lpinit (CCchunklp **lp_p, const char *lp_name, int lp_nrows,
        double *xstar)
{
    int  rval = 0;
    int  i;
    char  sense = 'E';
    double  rhs;
#ifdef SEPARATE_OPTIMIZE
    double obj[1];
    int cmatbeg[1];
    int cmatind[1];
    double cmatval[1];
    double lb[1];
    double ub[1];
#endif

    *lp_p = (CCchunklp *) NULL;

    *lp_p = CC_SAFE_MALLOC (1, CCchunklp);
    if ( *lp_p == (CCchunklp *) NULL ) {
        rval = 1; fprintf (stderr, "lp_p allocation failed\n");
        return rval;
    }

    (*lp_p)->lp        = (CClp *)    NULL;
    (*lp_p)->active    = (int *)     NULL;
    (*lp_p)->cmatind   = (int *)     NULL;
    (*lp_p)->cmatval   = (double *)  NULL;
    (*lp_p)->pi        = (double *)  NULL;
    (*lp_p)->nrows     = 0;
    (*lp_p)->extracols = 0;

    rval = CClp_init (&((*lp_p)->lp));
    if ( rval ) {
        fprintf (stderr, "CClp_init failed\n");
        goto CLEANUP;
    }

    rval = CClp_create ((*lp_p)->lp, lp_name);
    if ( rval ) {
        fprintf (stderr, "CClp_create failed\n");
        goto CLEANUP;
    }

    rval = CClp_tune_small ((*lp_p)->lp);
    if (rval) {
        fprintf (stderr, "CClp_tune_small failed\n");
        goto CLEANUP;
    }

    rval = CClp_disable_presolve ((*lp_p)->lp);
    if (rval) {
        fprintf (stderr, "CClp_disable_presolve failed\n");
        goto CLEANUP;
    }

    /* Number of edges in chunk.  Will be used later, together with
       lp->active, for reporting separating hyperplanes.  The calling
       routine expects hyperplane normals that are indexed on all
       edges, even if some are no longer active */

    if ( lp_nrows <= 0 ) {
       rval = 1;  fprintf (stderr, "lp_nrows <= 0\n");
       goto CLEANUP;
    }

    (*lp_p)->nrows = lp_nrows;

    (*lp_p)->active  = CC_SAFE_MALLOC (lp_nrows, int);
    (*lp_p)->cmatind = CC_SAFE_MALLOC (lp_nrows+1, int);
    (*lp_p)->cmatval = CC_SAFE_MALLOC (lp_nrows+1, double);
    (*lp_p)->pi      = CC_SAFE_MALLOC (lp_nrows+1, double);
    if ( (*lp_p)->active  == (int *)    NULL ||
         (*lp_p)->cmatind == (int *)    NULL ||
         (*lp_p)->cmatval == (double *) NULL ||
         (*lp_p)->pi      == (double *) NULL   ) {
        rval = 1; fprintf (stderr, "Not enough memory\n");
        goto CLEANUP;
    }

#ifdef SEPARATE_NORML1
    rhs = 0.0;
#else
    rhs = 1.0;
#endif

    rval = CClp_new_row ((*lp_p)->lp, sense, rhs);
    if ( rval ) {
        fprintf (stderr, "CClp_new_row failed\n");
        goto CLEANUP;
    }

    for (i = 0; i < lp_nrows; i++) {
        (*lp_p)->active[i] = 1; /* yes <==> not 0 */
#ifdef SEPARATE_NORML1
        rhs = 0.0;
#else
        rhs = xstar[i];
#endif
        rval = CClp_new_row ((*lp_p)->lp, sense, rhs);
        if ( rval ) {
            fprintf (stderr, "CClp_new_row failed\n");
            goto CLEANUP;
        }
    }

#ifdef SEPARATE_NORML1
    obj[0] = -1.0;
    cmatbeg[0] = 0;
    lb[0] = 0.0;
    if (SEPARATE_EPSILON == 0.0) {
        fprintf (stderr, "SEPARATE_NORML1 does not work with SEPARATE_EPSILON == 0.0\n");
        rval = 1; goto CLEANUP;
    }
    ub[0] = 2.0/SEPARATE_EPSILON;
    (*lp_p)->cmatval[0] = -1.0;
    (*lp_p)->cmatind[0] = 0;
    for (i=0; i<lp_nrows; i++) {
        (*lp_p)->cmatind[i+1] = i+1;
        (*lp_p)->cmatval[i+1] = -xstar[i];
    }
    rval = CClp_addcols ((*lp_p)->lp, 1, lp_nrows+1, obj, cmatbeg, (*lp_p)->cmatind,
                         (*lp_p)->cmatval, lb, ub);
    if (rval) {
        fprintf (stderr, "CClp_addcols failed\n");
        goto CLEANUP;
    }
    (*lp_p)->extracols++;
#endif /* SEPARATE_NORML1 */
    
#ifdef SEPARATE_OPTIMIZE
    for (i=0; i<=lp_nrows; i++) {

#ifdef SEPARATE_STRICTAFFINE
        if (i == 0) continue;
#endif /* SEPARATE_STRICTAFFINE */

#ifdef SEPARATE_NORML1

        obj[0] = 0.0;
        cmatbeg[0] = 0;
        cmatind[0] = i;
        cmatval[0] = 1.0;
        lb[0] = -1.0;
        ub[0] = 1.0;
        rval = CClp_addcols ((*lp_p)->lp, 1, 1, obj, cmatbeg, cmatind,
                             cmatval, lb, ub);
        if (rval) {
            fprintf (stderr, "CClp_addcols failed\n");
            goto CLEANUP;
        }
        (*lp_p)->extracols++;

#else /* SEPARATE_NORML1 */

        obj[0] = 1.0;
        cmatbeg[0] = 0;
        cmatind[0] = i;
        cmatval[0] = 1.0;
        lb[0] = 0.0;
        ub[0] = 1e10;
        rval = CClp_addcols ((*lp_p)->lp, 1, 1, obj, cmatbeg, cmatind,
                             cmatval, lb, ub);
        if (rval) {
            fprintf (stderr, "CClp_addcols failed\n");
            goto CLEANUP;
        }
        (*lp_p)->extracols++;

        obj[0] = 1.0;
        cmatbeg[0] = 0;
        cmatind[0] = i;
        cmatval[0] = -1.0;
        lb[0] = 0.0;
        ub[0] = 1e10;
        rval = CClp_addcols ((*lp_p)->lp, 1, 1, obj, cmatbeg, cmatind,
                             cmatval, lb, ub);
        if (rval) {
            fprintf (stderr, "CClp_addcols failed\n");
            goto CLEANUP;
        }
        (*lp_p)->extracols++;

#endif /* SEPARATE_NORML1 */
    }        
#endif /* SEPARATE_OPTIMIZE */

    /* The first entry in every column is the same */

    (*lp_p)->cmatind[0] = 0;
    (*lp_p)->cmatval[0] = 1.0;

    rval = 0;

CLEANUP:

    if ( rval )  CCchunk_lpfree (lp_p);
    return rval;

} /* END CCchunk_lpinit */



void CCchunk_lpfree (CCchunklp **lp_p)
{
    if ( *lp_p != (CCchunklp *) NULL ) {
        CC_IFFREE ((*lp_p)->pi, double);
        CC_IFFREE ((*lp_p)->cmatval, double);
        CC_IFFREE ((*lp_p)->cmatind, int);
        CC_IFFREE ((*lp_p)->active, int);
        CClp_free (&((*lp_p)->lp));
        CC_IFFREE (*lp_p, CCchunklp);
    }
} /* END CCchunk_lpfree */



int CCchunk_lpaddcol (CCchunklp *lp, double *x)
{
    int  rval = 0;
    int  i, k = 1, nzs = 1, cmatbeg = 0;
    double  obj = 0.0;
    double  lb = 0.0;
    double  ub = 1e10;

    for (i = 0; i < lp->nrows; i++) {
        if ( lp->active[i] ) {
            if ( x[i] ) {
               lp->cmatind[nzs] = k;
               lp->cmatval[nzs] = x[i];
               nzs++;
            }
            k++;
        }
    }

    rval = CClp_addcols (lp->lp, 1 /* no. cols */,
                       nzs /* no. nzs */, &obj, &cmatbeg, lp->cmatind,
                       lp->cmatval, &lb, &ub);
    if ( rval ) {
        fprintf (stderr, "CClp_addcols failed\n");
        goto CLEANUP;
    }

CLEANUP:

    return rval;

} /* END CCchunk_lpaddcol */



int CCchunk_lprelaxrow (CCchunklp *lp, int del_row)
{
    int  rval = 0;
    int  i, k = 0;

#ifdef SEPARATE_STRICTAFFINE
    if (CClp_ncols (lp->lp) == lp->extracols) {
/*        fprintf (stderr, "Asking for trouble, relaxing row from initial problem\n");*/
    }
#endif
        
    if ( del_row >= lp->nrows  ||  del_row < 0 ) {
        rval = 1; fprintf (stderr, "Illegal row index %d\n", del_row);
        return rval;
    }

    if ( lp->active[del_row] == 0 ) {
        fprintf (stderr, "Row %d already inactive\n", del_row);
        return rval;
    }

    for (i = 0; i <= del_row; i++) {
        if ( lp->active[i] )  k++;
    }

    rval = CClp_delete_row (lp->lp, k);
    if ( rval ) {
        fprintf (stderr, "CClp_delete_row failed\n");
    } else {
        lp->active[del_row] = 0;
    }

    return rval;

} /* END CCchunk_lprelaxrow */



int CCchunk_lpsolve (CCchunklp *lp, int *lpstatus_p, double *c,
                     double *alpha_p)
{
    int  rval = 0;
    int  i, k = 0;
#ifdef SEPARATE_OPTIMIZE
    double obj = 0.0;
#endif

#ifdef SEPARATE_STRICTAFFINE
    if (CClp_ncols (lp->lp) == lp->extracols) {
        /* the strict affine version of optimize is infeasible in this
           special case */
        *alpha_p = -1;
        for (i = 0; i < lp->nrows; i++) {
            if (lp->active[i]) c[i] = 1.0;
            else c[i] = 0.0;
        }
        *lpstatus_p = CC_CHUNK_LPINFEASIBLE;
        return 0;
    }
#endif

    *lpstatus_p = 0;
    *alpha_p = 0.0;
    for (i = 0; i < lp->nrows; i++)  c[i] = 0.0;

    rval = CClp_opt (lp->lp, CClp_METHOD_PRIMAL);
#if 0
    CClp_dump_lp (lp->lp, "chunk.sav");
#endif
#ifdef SEPARATE_OPTIMIZE
    if ( rval ) {
        fprintf (stderr, "CClp_opt failed, rval %d\n", rval);
        return 1;
    }
    rval = CClp_objval (lp->lp, &obj);
    if (rval) {
        fprintf (stderr, "CClp_objval failed, rval %d\n", rval);
        return 1;
    }
#ifdef SEPARATE_NORML1
    if (obj == 0.0) obj = 2*SEPARATE_EPSILON;
    else obj = -1.0/obj;
#endif /* SEPARATE_NORML1 */
    if ( obj <= SEPARATE_EPSILON ) {
        *lpstatus_p = CC_CHUNK_LPFEASIBLE;
        return 0;
    }
#else /* SEPARATE_OPTIMIZE */
    if (rval != 0 && rval != 2) {
        fprintf (stderr, "CClp_opt failed\n");
        return rval;
    }
    if (rval == 0) {
        *lpstatus_p = CC_CHUNK_LPFEASIBLE;
        return 0;
    }
#endif /* SEPARATE_OPTIMIZE */
    rval = CClp_pi (lp->lp, lp->pi);
    if (rval) {
        fprintf (stderr, "CClp_pi failed\n");
        return rval;
    }
    k = 0;
#ifdef SEPARATE_NORML1
    *alpha_p = - lp->pi[k++] * obj;
#else
    *alpha_p = - lp->pi[k++];
#endif
    for (i = 0; i < lp->nrows; i++) {
#ifdef SEPARATE_NORML1
        if ( lp->active[i] )  c[i] = lp->pi[k++] * obj;
        else                  c[i] = 0.0;
#else
        if ( lp->active[i] )  c[i] = lp->pi[k++];
        else                  c[i] = 0.0;
#endif
    }
    *lpstatus_p = CC_CHUNK_LPINFEASIBLE;
    return 0;
} /* END CCchunk_lpsolve */


int CCchunk_lpbasis (CCchunklp *lp, int ncols, int *basis)
{
    int rval;
    CClp_info *info = (CClp_info *) NULL;
    int i;

    rval = CClp_get_info (lp->lp, &info);
    if (rval) {
        fprintf (stderr, "CClp_get_info failed\n");
        goto CLEANUP;
    }

    for (i=0; i<ncols; i++) {
        basis[i] = CClp_is_col_active (info, lp->extracols + i);
    }

    rval = 0;

  CLEANUP:
    CClp_free_info (&info);
    return rval;
}

#if 0

/* ONLY HERE FOR USE IN RAYS; USES DIRECT CPLEX LIBRARY CALLS INSTEAD OF
   THE LP INTERFACE, BECAUSE THE LP INTERFACE DOESN'T SUPPORT DIRECT
   MODIFICATION OF COEFFICIENTS */
/* #include <cplex.h> */
typedef struct CPXENV *CPXENVptr;
typedef struct CPXLP *CPXLPptr;
#ifdef __stdcall
#define CPXPUBLIC __stdcall
#else
#define CPXPUBLIC
#endif

int CPXPUBLIC CPXchgcoef (CPXENVptr, CPXLPptr, int, int, double);

typedef struct CClp_parameters {
    int    scrind;
    int    simdisplay;
    int    advind;
    int    dpriind;
    int    ppriind;
    double epper;
    double epopt;
    double eprhs;
    int perind;
    int preind;
    int aggind;
} CClp_parameters;

struct CClp {
    CPXENVptr cplex_env;
    CPXLPptr  cplex_lp;
    CClp_parameters cplex_params;
};


int CCchunk_lprhs (CCchunklp *lp, double *xstar);
int CCchunk_lprhs (CCchunklp *lp, double *xstar)
{
    int i;
    int j;
    int rval;
#ifdef SEPARATE_NORML1
    int col = 0;
    double mul = -1.0;
#else
    int col = -1;
    double mul = 1.0;
#endif
    j = 1;
    for (i=0; i<lp->nrows; i++) {
        if (lp->active[i]) {
            rval = CPXchgcoef (lp->lp->cplex_env, lp->lp->cplex_lp, j, col,
                               mul * xstar[i]);
            if (rval) {
                fprintf (stderr, "CPXchgcoef failed\n");
                return rval;
            }
            j++;
        }
    }
    return 0;
}

void CCchunk_lpdump (CCchunklp *lp, char *fname);
void CCchunk_lpdump (CCchunklp *lp, char *fname)
{
    CClp_dump_lp (lp->lp, fname);
}

#endif /* RAY SHOOTING */
