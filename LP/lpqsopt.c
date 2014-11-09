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
/*      Interface Routines to The beta-version of the QSopt LP Solver       */
/*                                                                          */
/*  NOTE: Use this code in place of lp_none.c to access the QSopt beta      */
/*   libraries. You will also need to link the QSopt  library via the       */
/*   makefile.                                                              */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: May 17, 2000 (bico)                                               */
/*        January 13, 2002 (bico)                                           */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CClp_init (CClp **lp)                                               */
/*    INITIALIZES the LP.                                                   */
/*                                                                          */
/*  int CClp_force_perturb (CClp *lp)                                       */
/*    Forces a perturbation in the LP simplex solves                        */
/*                                                                          */
/*  int CClp_tune_small (CClp *lp)                                          */
/*    SETS solver options for tiny problems.                                */
/*                                                                          */
/*  int CClp_disable_presolve (CClp *lp)                                    */
/*    DISABLES the solvers presolve.                                        */
/*                                                                          */
/*  void CClp_free (CClp **lp)                                              */
/*    FREES the LP, both the allocation in CClp_init () and the             */
/*    allocation in CClp_loadlp.                                            */
/*                                                                          */
/*  void CClp_freelp (CClp **lp)                                            */
/*    FREES the LP loaded in CClp_loadlp (or CClp_create), but does         */
/*    not free the data allocated by CClp_init.                             */
/*                                                                          */
/*  int CClp_loadlp (CClp *lp, const char *name, int ncols, int nrows,      */
/*      int objsense, double *obj, double *rhs, char *sense,                */
/*      int *matbeg, int *matcnt,                                           */
/*      int *matind, double *matval,                                        */
/*      double *lb, double *ub)                                             */
/*    LOADS the data into the LP.                                           */
/*      -name attaches a name to the LP (it can be used by the LP solver    */
/*       in io routines)                                                    */
/*      -ncols and nrows give the number of columns and rows in the LP      */
/*      -objsense should be 1 for minimize and -1 for maximize              */
/*      -obj and rhs are arrays giving the objective function and rhs       */
/*      -sense is an array specifying 'L', 'E', or 'G' for each of the      */
/*       rows                                                               */
/*      -matbeg, matcnt, matind, and matval give the coefficients of        */
/*       the constraint matrix in column by column order. matbeg gives      */
/*       gives the index of the start of each column; matcnt gives the      */
/*       number of coefficients in each column; matind gives the indices    */
/*       of the rows where the coefficients are located in the constraint   */
/*       matrix (so for column j, the indices are given in matcnt[j]        */
/*       locations starting at matind[matbeg[j]]; and matval gives the      */
/*       actual coefficients (organized like matind).                       */
/*      -lb and ub are arrays giving the upper and lower bounds of          */
/*       the variables.                                                     */
/*                                                                          */
/*  int CClp_create (CClp *lp, const char *name)                            */
/*    CREATES an empty lp.  This supports an alternative to CClp_loadlp     */
/*    for loading a problem.                                                */
/*      -name attaches a name to the LP (it can be used by the LP solver    */
/*       in io routines)                                                    */
/*                                                                          */
/*  int CClp_new_row (CClp *lp, char sense, double rhs)                     */
/*    ADDS a new empty row to the lp                                        */
/*      -sense is 'L', 'E', or 'G' for a <=, =, or >= constraint            */
/*      -rhs is the right-hand side of the row                              */
/*                                                                          */
/*  int CClp_change_sense (CClp *lp, int row, char sense)                   */
/*    CHANGES the sense of a row                                            */
/*      -row is the row number to change                                    */
/*      -sense is 'L', 'E', or 'G' to change to <=, =, or >=                */
/*                                                                          */
/*  int CClp_opt (CClp *lp, int method)                                     */
/*    CALLS designated LP solution method.                                  */
/*                                                                          */
/*  int CClp_limited_dualopt (CClp *lp, int lim, int *status,               */
/*      double *upperbound)                                                 */
/*    CALLS the dual simplex method with a limit on the number of pivots.   */
/*      -upperbound it is used to cutoff the dual simplex method (when      */
/*       the objective value reaches upperbound); it can be NULL            */
/*      -status returns the status of the optimization (it can be NULL)     */
/*                                                                          */
/*  int CClp_addrows (CClp *lp, int newrows, int newnz, double *rhs,        */
/*      char *sense, int *rmatbeg, int *rmatind, double *rmatval)           */
/*    ADDS the rows to the LP.                                              */
/*      -newrows is the number of rows to be added                          */
/*      -newnz is the number of nonzero coefficients in the new rows        */
/*      -rhs is an array of the rhs values for the new rows                 */
/*      -sense is 'L', 'E', or 'G' for each of the new rows                 */
/*      -rmatbeg, rmatind, and rmatval give the coefficients of the         */
/*       new rows in sparse format. The arrays can be freed after the       */
/*       call.                                                              */
/*                                                                          */
/*  int CClp_addcols (CClp *lp, int newcols, int newnz, double *obj,        */
/*      int *cmatbeg, int *cmatind, double *cmatval,                        */
/*      double *lb, double *ub)                                             */
/*    ADDS the columns to the LP.                                           */
/*                                                                          */
/*  int CClp_delete_row (CClp *lp, int i)                                   */
/*    DELETES row i of the LP.                                              */
/*                                                                          */
/*  int CClp_delete_set_of_rows (CClp *lp, int *delstat)                    */
/*    DELETES the rows corresponding to 1 entries in delstat.               */
/*      -delstat is a 0/1 array having an entry for each row                */
/*                                                                          */
/*  int CClp_delete_column (CClp *lp, int i)                                */
/*    DELETES column i from the LP.                                         */
/*                                                                          */
/*  int CClp_delete_set_of_columns (CClp *lp, int *delstat)                 */
/*    DELETES the columns corresponding to the 1 entries in delstat.        */
/*      -delstat is a 0/1 array having an entry for each column             */
/*                                                                          */
/*  int CClp_setbnd (CClp *lp, int col, char lower_or_upper, double bnd)    */
/*    SETS the bound on the variable index by col.                          */
/*      -lower_or_upper should be either 'L' or 'U'                         */
/*                                                                          */
/*  int CClp_get_warmstart (CClp *lp, CClp_warmstart **w)                   */
/*    SAVES information for efficiently resolving the current lp in w,      */
/*    for example, basis or norm information                                */
/*                                                                          */
/*  int CClp_load_warmstart (CClp *lp, CClp_warmstart *w)                   */
/*    RESTORES the warmstart information in w.                              */
/*                                                                          */
/*  int CClp_build_warmstart (CClp_warmstart **w, CClp_info *i)             */
/*    BUILDS some warmstart information from the row/column information     */
/*    in i.                                                                 */
/*                                                                          */
/*  void CClp_free_warmstart (CClp_warmstart **w)                           */
/*    FREES the memory used by w.                                           */
/*                                                                          */
/*  int CClp_sread_warmstart (CC_SFILE *f, CClp_warmstart **w)              */
/*    READS warmstart information from the f.                               */
/*                                                                          */
/*  int CClp_swrite_warmstart (CC_SFILE *f, CClp_warmstart *w)              */
/*    WRITES warmstart information from the f.                              */
/*                                                                          */
/*  int CClp_get_info (CClp *lp, CClp_info **i)                             */
/*    BUILDS information useful for efficiently answering questions         */
/*    about the status of rows and columns                                  */
/*                                                                          */
/*  int CClp_create_info (CClp_info **i, int rcount, int ccount)            */
/*    CREATES a structure for storing information about the status of       */
/*    rows and columns.                                                     */
/*                                                                          */
/*  int CClp_is_col_active (CClp_info *i, int c)                            */
/*    returns 1 if column e is active, 0 otherwise.                         */
/*    "active" means participating in the current solution (for example,    */
/*    it could mean basic or nonbasic at upper bound)                       */
/*                                                                          */
/*  int CClp_is_row_active (CClp_info *i, int r)                            */
/*    returns 1 if row e is active, 0 otherwise.                            */
/*    "active" means participating in the current solution (for example,    */
/*    it could mean the row's slack is non-basic)                           */
/*                                                                          */
/*  void CClp_set_col_active (CClp_info *i, int c)                          */
/*    marks column e as active (for eventual CClp_build_warmstart)          */
/*                                                                          */
/*  void CClp_set_col_inactive (CClp_info *i, int c)                        */
/*    marks column e as inactive (for eventual CClp_build_warmstart)        */
/*                                                                          */
/*  void CClp_set_col_upper (CClp_info *i, int c)                           */
/*    marks column e as active at upper bound (for eventual                 */
/*    CClp_build_warmstart)                                                 */
/*                                                                          */
/*  void CClp_set_row_active (CClp_info *i, int r)                          */
/*    marks row r as active (for eventual CClp_build_warmstart)             */
/*                                                                          */
/*  void CClp_set_row_inactive (CClp_info *i, int r)                        */
/*    marks row r as inactive (for eventual CClp_build_warmstart)           */
/*                                                                          */
/*  void CClp_free_info (CClp_info **i)                                     */
/*    FREES the memory used by i.                                           */
/*                                                                          */
/*  int CClp_x (CClp *lp, double *x)                                        */
/*    RETURNS the current LP solution.                                      */
/*      -x should be an array of length at least ncols                      */
/*                                                                          */
/*  int CClp_rc (CClp *lp, double *rc)                                      */
/*    RETURNS the current reduced costs.                                    */
/*      -rc should be an array of length at least ncols                     */
/*                                                                          */
/*  int CClp_pi (CClp *lp, double *pi)                                      */
/*    RETURNS the dual values on the constraints.                           */
/*      -pi should be an array of length at least nrows                     */
/*    NOTES: If the lp and dual lp are feasible, these pi values are        */
/*      the traditional dual solution.  If the dual is unbounded, these     */
/*      pi satisfy                                                          */
/*                                                                          */
/*       pi_i <= 0  for <= constraints                                      */
/*       pi_i >= 0  for >= constraints                                      */
/*                                                                          */
/*       pi'b - sum (pi'A_j * u_j: pi'A_j > 0)                              */
/*            - sum (pi'A_j * l_j: pi'A_j < 0) > 0                          */
/*                                                                          */
/*    where b is the rhs vector, u_j is the upper bound on variable x_j,    */
/*    l_j the lower bound, and A_j the constraint matrix column for x_j.    */
/*                                                                          */
/*  int CClp_objval (CClp *lp, double *obj)                                 */
/*    RETURNS the objective value of the lp.                                */
/*                                                                          */
/*  int CClp_nrows (CClp *lp)                                               */
/*    RETURNS the number of rows in the LP.                                 */
/*                                                                          */
/*  int CClp_ncols (CClp *lp)                                               */
/*    RETURNS the number of columns in the LP.                              */
/*                                                                          */
/*  int CClp_nnonzeros (CClp *lp)                                           */
/*    RETURNS the number of nonzeros in the LP.                             */
/*                                                                          */
/*  int CClp_status (CClp *lp, int *status)                                 */
/*    CHECKS whether the current lp is infeasible or whether an optimal     */
/*     solution has been found. It returns an error if the LP has not       */
/*     not been optimized.                                                  */
/*      -lp is the lp                                                       */
/*      -status returns 0 if the lp has an optimal solution and 1 if it     */
/*       is infeasible.                                                     */
/*                                                                          */
/*  int CClp_getweight (CClp *lp, int nrows, int *rmatbeg, int *rmatind,    */
/*      double *rmatval, double *weight)                                    */
/*    COMPUTES the duals of the steepest edge norms for the n rows          */
/*     specified in rmatbeg, rmatind, and rmatval.                          */
/*      -weight returns the array of weights; the array should be at        */
/*       least nrows long                                                   */
/*                                                                          */
/*  int CClp_dump_lp (CClp *lp, const char *fname)                          */
/*    WRITES the LP to file fname.                                          */
/*                                                                          */
/*  int CClp_getgoodlist (CClp *lp, int *goodlist, int *goodlen_p,          */
/*      double *downpen, double *uppen)                                     */
/*    RETURNS an array of the column indices corresponding to variables     */
/*     that move in both the up and down directions. This is a useful       */
/*     list of candidates for strong branching.                             */
/*      -goodlist, downpen and uppen should be arrays of length at          */
/*       least ncols.                                                       */
/*                                                                          */
/*  int CClp_strongbranch (CClp *lp, int *candidatelist, int ncand,         */
/*      double *downpen, double *uppen, int iterations,                     */
/*      double upperbound)                                                  */
/*    RETURNS estimates of the lp values obtained by setting each of the    */
/*      ncand variables listed in candidatelist to 0 and 1. The estimates   */
/*      are obtained by performing iterations pivots of dual simplex        */
/*      method. upperbound is used to cutoff the dual simplex method.       */
/*      downpen and uppen should never be > upperbound.                     */
/*       -downpen and uppen should be arrays of length at least ncand       */
/*                                                                          */
/****************************************************************************/



#include "machdefs.h"
#include "util.h"
#include "macrorus.h"
#include "lp.h"
#include "tsp.h"
#include <qsopt.h>

struct CClp {
    QSprob p;
};

struct CClp_warmstart {
    int nstruct;
    int nrows;
    char *cstat;
    char *rstat;
    double *dnorm;
};

struct CClp_info {
    int nstruct;
    int nrows;
    char *cstat;
    char *rstat;
};

#define SOLVER_WARMSTART_NAME "QSO"


static void
    init_warmstart (CClp_warmstart *w),
    lp_message (void);

static int
     primalopt (CClp *lp),
     dualopt (CClp *lp);



static void lp_message (void)
{
    fprintf (stderr, "need to link an lp solver to use this function\n");
}

int CClp_init (CClp **lp)
{
    int rval = 0;
 
    CClp_free (lp);

    (*lp) = CC_SAFE_MALLOC (1, CClp);
    if ((*lp) == (CClp *) NULL) {
        fprintf (stderr, "Out of memory in CClp_init\n");
        rval = 1; goto CLEANUP;
    }

    (*lp)->p = (QSprob) NULL;

CLEANUP:

    return rval;
}

int CClp_force_perturb (CClp *lp)
{
    int rval = 0;

    /* Not implemented in illuin */

    if (!lp) {
        fprintf (stderr, "CClp_force_perturb called without an lp\n");
        rval = 1; goto CLEANUP;
    }

CLEANUP:

    return rval;
}

int CClp_tune_small (CClp *lp)
{
    int rval = 0;

    if (!lp) {
        fprintf (stderr, "CClp_tune_small called without an LP\n");
        rval = 1; goto CLEANUP;
    }

    rval = QSset_param (lp->p, QS_PARAM_PRIMAL_PRICING, QS_PRICE_PDEVEX);
    if (rval) {
        fprintf (stderr, "QSset_param failed\n");
        rval = 1; goto CLEANUP;
    }

    rval = QSset_param (lp->p, QS_PARAM_DUAL_PRICING, QS_PRICE_DDANTZIG);
    if (rval) {
        fprintf (stderr, "QSset_param failed\n");
        rval = 1; goto CLEANUP;
    }

    rval = QSset_param (lp->p, QS_PARAM_SIMPLEX_DISPLAY, 0);
    if (rval) {
        fprintf (stderr, "QSset_param failed\n");
        rval = 1; goto CLEANUP;
    }

    rval = QSset_param (lp->p, QS_PARAM_SIMPLEX_SCALING, 0);
    CCcheck_rval (rval, "QSset_param failed");

CLEANUP:

    return rval;
}

int CClp_disable_presolve (CClp *lp)
{
    int rval = 0;

    if (!lp) {
        fprintf (stderr, "CClp_tune_small called without an LP\n");
        rval = 1; goto CLEANUP;
    }

CLEANUP:

    return rval;
}

void CClp_free (CClp **lp)
{
    if (*lp) {
        QSfree_prob ((*lp)->p);
        CC_FREE (*lp, CClp);
    }
}

void CClp_freelp (CClp **lp)
{
    if (*lp) {
        QSfree_prob ((*lp)->p);
        (*lp)->p = (QSprob) NULL;
    }
}

int CClp_loadlp (CClp *lp, const char *name, int ncols, int nrows,
        int objsense, double *obj, double *rhs, char *sense, int *matbeg,
        int *matcnt, int *matind, double *matval, double *lb, double *ub)

{
    int rval = 0;

    lp->p = QSload_prob (name, ncols, nrows, matcnt, matbeg, matind,
                matval, objsense, obj, rhs, sense, lb, ub, (const char **) NULL,
                (const char **) NULL);
    if (lp->p == (QSprob) NULL) {
        fprintf (stderr, "QSload_prob failed\n");
        rval = 1; goto CLEANUP;
    }

    rval = QSset_param (lp->p, QS_PARAM_DUAL_PRICING, QS_PRICE_DSTEEP);
    if (rval) {
        fprintf (stderr, "QSset_param failed\n");
        rval = 1; goto CLEANUP;
    }

    rval = QSset_param (lp->p, QS_PARAM_SIMPLEX_DISPLAY, 0);
    if (rval) {
        fprintf (stderr, "QSset_param failed\n");
        rval = 1; goto CLEANUP;
    }

    rval = QSset_param (lp->p, QS_PARAM_SIMPLEX_SCALING, 0);
    CCcheck_rval (rval, "QSset_param failed");

CLEANUP:

    return rval;
}

int CClp_create (CClp *lp, const char *name)
{
    int rval = 0;

    if (!lp) {
        fprintf (stderr, "CClp_create called without an lp structure\n");
        rval = 1; goto CLEANUP;
    }

    lp->p = QScreate_prob (name, QS_MIN);
    if (lp->p == (QSprob) NULL) {
        fprintf (stderr, "QScreate_prob failed\n");
        rval = 1; goto CLEANUP;
    }

    rval = QSset_param (lp->p, QS_PARAM_DUAL_PRICING, QS_PRICE_DSTEEP);
    if (rval) {
        fprintf (stderr, "QSset_param failed\n");
        rval = 1; goto CLEANUP;
    }

    rval = QSset_param (lp->p, QS_PARAM_SIMPLEX_SCALING, 0);
    CCcheck_rval (rval, "QSset_param failed");

CLEANUP:

    return rval;
}

int CClp_new_row (CClp *lp, char sense, double rhs)
{
    int rval = 0;

    rval = QSnew_row (lp->p, rhs, sense, (char *) NULL);
    if (rval) {
        fprintf (stderr, "QSnew_row failed\n");
        rval = 1; goto CLEANUP;
    }

CLEANUP:

    return rval;
}

int CClp_change_sense (CClp *lp, int row, char sense)
{
    int rval = 0;

    rval = QSchange_sense (lp->p, row, sense);
    if (rval) {
        fprintf (stderr, "QSchange_sense failed\n");
        goto CLEANUP;
    }

CLEANUP:

    return rval;
}

int CClp_opt (CClp *lp, int method)
{
    int rval = 0;

    switch (method) {
    case CClp_METHOD_PRIMAL:
        rval = primalopt (lp);
        break;
    case CClp_METHOD_DUAL:
        rval = dualopt (lp);
        break;
    case CClp_METHOD_BARRIER:
        fprintf (stderr, "qsopt does not yet have a barrier code\n");
        rval = 1; goto CLEANUP;
    default:
        rval = 1;
        fprintf (stderr, "Nonexistent method in CClp_opt\n");
        break;
    }

CLEANUP:

    return rval;
}

static int primalopt (CClp *lp)
{
    int rval = 0;
    int status; 
/*
    double szeit = CCutil_zeit ();
    static double trickit = 0.0;
    static int tiii = 0;
    int ncols, nrows;
    static double tiii_cols = 0.0;
    static double tiii_rows = 0.0;
*/

    rval = QSopt_primal (lp->p, &status);
    if (rval) {
        fprintf (stderr, "QSopt_primal failed\n"); goto CLEANUP;
    }

/*
    ncols = QSget_colcount (lp->p);
    nrows = QSget_rowcount (lp->p);

    tiii_cols += (double) ncols;
    tiii_rows += (double) nrows;

    trickit += (CCutil_zeit () - szeit);
    if (tiii++ % 1000 == 999) {
        printf ("I-LP: %f (%d, %f), rows = %f, cols = %f\n",
             trickit, tiii, trickit / (double) tiii, 
             tiii_rows / (double) tiii, tiii_cols / (double) tiii);
        fflush (stdout);
    }
*/

    if (status == QS_LP_ITER_LIMIT) {
        printf ("Primal LP Solver reached iteration limit\n"); fflush (stdout);
        rval = QSwrite_prob (lp->p, "piter.lp", "LP");
        if (rval) {
            fprintf (stderr, "QSwrite_prob failed\n"); goto CLEANUP;
        }
        printf ("Saved LP as piter.lp\n"); fflush (stdout);
        rval = QSwrite_basis (lp->p, (QSbas) NULL, "piter.bas");
        if (rval) {
            fprintf (stderr, "QSwrite_basis failed\n"); goto CLEANUP;
        }
        printf ("Saved LP as piter.bas\n"); fflush (stdout);
    } else if (status == QS_LP_TIME_LIMIT) {
        printf ("Primal LP Solver reached time limit\n"); fflush (stdout);
    }

    if (status == QS_LP_INFEASIBLE) {
        rval = 2; goto CLEANUP;
    } else if (status != QS_LP_OPTIMAL) {
        fprintf (stderr, "no optimal LP-solution exists: %d\n", status);
        rval = 1; goto CLEANUP;
    }

CLEANUP:

    return rval;
}

static int dualopt (CClp *lp)
{
    int rval = 0;
    int status; 

    rval = QSopt_dual (lp->p, &status);
    if (rval) {
        fprintf (stderr, "QSopt_dual failed\n"); goto CLEANUP;
    }

    if (status == QS_LP_ITER_LIMIT) {
        printf ("Dual LP Solver reached iteration limit\n"); fflush (stdout);
        rval = QSwrite_prob (lp->p, "iter.lp", "LP");
        if (rval) {
            fprintf (stderr, "QSwrite_prob failed\n");
            goto CLEANUP;
        }
        printf ("Saved LP as iter.lp\n"); fflush (stdout);
        rval = QSwrite_basis (lp->p, (QSbas) NULL, "iter.bas");
        if (rval) {
            fprintf (stderr, "QSwrite_basis failed\n"); goto CLEANUP;
        }
        printf ("Saved LP as iter.bas\n"); fflush (stdout);
    } else if (status == QS_LP_TIME_LIMIT) {
        printf ("Dual LP Solver reached time limit\n"); fflush (stdout);
    }

    if (status == QS_LP_INFEASIBLE) {
        rval = 2; goto CLEANUP;
    } else if (status != QS_LP_OPTIMAL) {
        fprintf (stderr, "no optimal LP-solution exists: %d\n", status);
        rval = 1; goto CLEANUP;
    }

CLEANUP:

    return rval;
}


int CClp_limited_dualopt (CClp *lp, int iterationlim, int *status,
        double *objupperlim)
{
    int rval = 0;
    int istatus, olditer = -1;
    double tt;

    if (objupperlim) tt = *objupperlim;   /* just to stop complier warning */

    /* Hack, till we understand what is happening */

    iterationlim = (iterationlim * 0.25) + 1;

    rval = QSget_param (lp->p, QS_PARAM_SIMPLEX_MAX_ITERATIONS, &olditer);
    if (rval) {
        fprintf (stderr, "QSget_param failed\n"); goto CLEANUP;
    }

    rval = QSset_param (lp->p, QS_PARAM_SIMPLEX_MAX_ITERATIONS, iterationlim);
    if (rval) {
        fprintf (stderr, "QSset_param failed\n"); goto CLEANUP;
    }

    rval = QSopt_dual (lp->p, &istatus);
    if (rval) {
        fprintf (stderr, "QSopt_dual failed\n"); goto CLEANUP;
    }

    if (istatus == QS_LP_INFEASIBLE) {
        if (status) *status = CClp_INFEASIBLE;
    } else if (istatus == QS_LP_UNBOUNDED) {
        if (status) *status = CClp_UNBOUNDED;
    } else if (istatus == QS_LP_UNSOLVED) {
        fprintf (stderr, "no optimal LP-solution exists\n");
        if (status) *status = CClp_FAILURE;
    } else {
        if (status) *status = CClp_SUCCESS;
    }

CLEANUP:

    if (olditer != -1) {
        rval = QSset_param (lp->p, QS_PARAM_SIMPLEX_MAX_ITERATIONS, olditer);
    }

    return rval;
}

int CClp_addrows (CClp *lp, int newrows, int newnz, double *rhs, char *sense,
        int *rmatbeg, int *rmatind, double *rmatval)
{
    int i, rval = 0;
    int *rmatcnt = (int *) NULL;

    rmatcnt = CC_SAFE_MALLOC (newrows, int);
    if (!rmatcnt) {
        fprintf (stderr, "out of memory in CClp_addrows\n");
        rval = 1; goto CLEANUP;
    }

    for (i = 0; i < newrows - 1; i++) {
        rmatcnt[i] = rmatbeg[i+1] - rmatbeg[i];
    }
    rmatcnt[newrows - 1] = newnz - rmatbeg[newrows - 1];

    rval = QSadd_rows (lp->p, newrows, rmatcnt, rmatbeg, rmatind, rmatval,
                        rhs, sense, (const char **) NULL);
    if (rval) {
        fprintf (stderr, "QSadd_rows failed\n"); goto CLEANUP;
    }

CLEANUP:

    CC_IFFREE (rmatcnt, int);
    return rval;
}

int CClp_addcols (CClp *lp, int newcols, int newnz, double *obj,
                  int *cmatbeg, int *cmatind, double *cmatval,
                  double *lb, double *ub)
{
    int i, rval = 0;
    int *cmatcnt = (int *) NULL;

    cmatcnt = CC_SAFE_MALLOC (newcols, int);
    if (!cmatcnt) {
        fprintf (stderr, "out of memory in CClp_addcols\n");
        rval = 1; goto CLEANUP;
    }

    for (i = 0; i < newcols - 1; i++) {
        cmatcnt[i] = cmatbeg[i+1] - cmatbeg[i];
    }
    cmatcnt[newcols - 1] = newnz - cmatbeg[newcols - 1];

    rval = QSadd_cols (lp->p, newcols, cmatcnt, cmatbeg, cmatind, cmatval,
                        obj, lb, ub, (const char **) NULL);
    if (rval) {
        fprintf (stderr, "QSadd_cols failed\n"); goto CLEANUP;
    }

CLEANUP:

    CC_IFFREE (cmatcnt, int);
    return rval;
}

int CClp_delete_row (CClp *lp, int i)
{
    int rval = 0;
    int dellist[1];

    dellist[0] = i;

    rval = QSdelete_rows (lp->p, 1, dellist);
    if (rval) {
        fprintf (stderr, "QSdelete_cols failed\n"); goto CLEANUP;
    }

CLEANUP:

    return rval;
}

int CClp_delete_set_of_rows (CClp *lp, int *delstat)
{
    int rval = 0;
    int *dellist = (int *) NULL;
    int delcnt = 0;
    int i, j, rcnt;

    rcnt = QSget_rowcount (lp->p);

    for (i = 0; i < rcnt; i++) {
        if (delstat[i]) delcnt++;
    }
    if (delcnt == 0) {
        fprintf (stderr, "delete_set_of_rows with no deleted rows\n");
        goto CLEANUP;
    }

    dellist = CC_SAFE_MALLOC (delcnt, int);
    if (dellist == (int *) NULL) {
        fprintf (stderr, "out of memory in delete_set_of_rows\n");
        return 1;
    }
    for (i = 0, j = 0; i < rcnt; i++) {
        if (delstat[i]) {
            dellist[j++] = i;
        }
    }

    rval = QSopt_pivotin_row (lp->p, delcnt, dellist);
    if (rval) {
        fprintf (stderr, "QSopt_pivotin_row failded, continuing anyway\n");
    }

    rval = QSdelete_rows (lp->p, delcnt, dellist);
    if (rval) {
        fprintf (stderr, "QSdelete_rows failed\n"); goto CLEANUP;
    }

CLEANUP:

    CC_IFFREE (dellist, int);
    return rval;
}

int CClp_delete_column (CClp *lp, int i)
{
    int rval = 0;
    int dellist[1];

    dellist[0] = i;
    
    rval = QSdelete_cols (lp->p, 1, dellist);
    if (rval) {
        fprintf (stderr, "QSdelete_cols failed\n"); goto CLEANUP;
    }

CLEANUP:

    return rval;
}

int CClp_delete_set_of_columns (CClp *lp, int *delstat)
{
    int rval = 0;
    int *dellist = (int *) NULL;
    int delcnt = 0;
    int i, j, ccnt;

    ccnt = QSget_colcount (lp->p);

    for (i = 0; i < ccnt; i++) {
        if (delstat[i]) delcnt++;
    }
    if (delcnt == 0) {
        fprintf (stderr, "delete_set_of_columns with no deleted columns\n");
        goto CLEANUP;
    }

    dellist = CC_SAFE_MALLOC (delcnt, int);
    if (dellist == (int *) NULL) {
        fprintf (stderr, "out of memory in delete_set_of_rows\n");
        return 1;
    }
    for (i = 0, j = 0; i < ccnt; i++) {
        if (delstat[i]) {
            dellist[j++] = i;
        }
    }

    rval = QSdelete_cols (lp->p, delcnt, dellist);
    if (rval) {
        fprintf (stderr, "QSdelete_cols failed\n"); goto CLEANUP;
    }

CLEANUP:

    return rval;
}

int CClp_setbnd (CClp *lp, int col, char lower_or_upper, double bnd)
{
    int rval = 0;
    int collist[1];
    char lu[1];
    double bounds[1];

    collist[0] = col;
    lu[0] = lower_or_upper;
    bounds[0] = bnd;

    rval = QSchange_bounds (lp->p, 1, collist, lu, bounds);
    if (rval) {
        fprintf (stderr, "QSchange_bounds failed\n"); goto CLEANUP;
    }

CLEANUP:

    return rval;
}

int CClp_get_warmstart (CClp *lp, CClp_warmstart **w)
{
    int rval = 0;
    int k;

    CClp_free_warmstart (w);

    (*w) = CC_SAFE_MALLOC (1, CClp_warmstart);
    if ((*w) == (CClp_warmstart *) NULL) {
        fprintf (stderr, "Out of memory in CClp_get_warmstart\n");
        rval = 1; goto CLEANUP;
    }
    init_warmstart (*w);

    k = QSget_colcount (lp->p);
    (*w)->nstruct = k;

    k = QSget_rowcount (lp->p);
    (*w)->nrows = k;

    (*w)->cstat = CC_SAFE_MALLOC ((*w)->nstruct, char);
    (*w)->rstat = CC_SAFE_MALLOC ((*w)->nrows, char);
    (*w)->dnorm = CC_SAFE_MALLOC ((*w)->nrows, double);

    if (!(*w)->cstat || !(*w)->rstat || !(*w)->dnorm) {
        fprintf (stderr, "out of memory in CClp_get_warmstart\n");
        rval = 1; goto CLEANUP;
    }

    if (QStest_row_norms (lp->p) == 0) {
        printf ("recomputing rownorms ...\n"); fflush (stdout);
        rval = QScompute_row_norms (lp->p);
        if (rval) {
            fprintf (stderr, "QScompute_row_norms failed\n");
            goto CLEANUP;
        }
    }

    rval = QSget_basis_and_row_norms_array (lp->p, (*w)->cstat, (*w)->rstat,
                                                    (*w)->dnorm);
    if (rval) {
        fprintf (stderr, "QSget_basis_and_row_norms_array failed\n");
        fprintf (stderr, "Trying to get basis\n");
        CC_IFFREE ((*w)->dnorm, double);
        rval = QSget_basis_array (lp->p, (*w)->cstat, (*w)->rstat);
        if (rval) {
            fprintf (stderr, "QSget_basis_array failed\n"); goto CLEANUP;
        }
    }

CLEANUP:

    return rval;
}

int CClp_load_warmstart (CClp *lp, CClp_warmstart *w)
{
    int rval = 0;

    if (w->cstat == (char *) NULL || w->rstat == (char *) NULL) {
        fprintf (stderr, "WARNING: no basis in call to load_warmstart\n");
        goto CLEANUP;
    }

    if (w->dnorm != (double *) NULL) {
        rval = QSload_basis_and_row_norms_array (lp->p, w->cstat, w->rstat,
                                                 w->dnorm);
        if (rval) {
            fprintf (stderr, "QSload_basis_and_row_norms_array failed");
            goto CLEANUP;
        }
    } else {
        rval = QSload_basis_array (lp->p, w->cstat, w->rstat);
        if (rval) {
            fprintf (stderr, "QSload_basis_array failed");
            goto CLEANUP;
        }
    }

CLEANUP:

    return rval;
}

int CClp_build_warmstart (CClp_warmstart **w, CClp_info *i)
{
    int rval = 0;
    int j;

    CClp_free_warmstart (w);

    (*w) = CC_SAFE_MALLOC (1, CClp_warmstart);
    if ((*w) == (CClp_warmstart *) NULL) {
        fprintf (stderr, "Out of memory in CClp_get_warmstart\n");
        rval = 1; goto CLEANUP;
    }
    init_warmstart (*w);

    (*w)->nstruct = i->nstruct;
    if ((*w)->nstruct == 0) {
        fprintf (stderr, "No columns in CClp_info\n");
        rval = 1; goto CLEANUP;
    }
    (*w)->nrows = i->nrows;
    if ((*w)->nrows == 0) {
        fprintf (stderr, "No rows in CClp_info\n");
        rval = 1; goto CLEANUP;
    }

    (*w)->cstat = CC_SAFE_MALLOC ((*w)->nstruct, char);
    (*w)->rstat = CC_SAFE_MALLOC ((*w)->nrows, char);
    if (!(*w)->cstat || !(*w)->rstat) {
        fprintf (stderr, "out of memory in CClp_get_warmstart\n");
        rval = 1; goto CLEANUP;
    }

    for (j = 0; j < (*w)->nstruct; j++) {
        (*w)->cstat[j] = i->cstat[j];
    }
    for (j = 0; j < (*w)->nrows; j++) {
        (*w)->rstat[j] = i->rstat[j];
    }

CLEANUP:

    if (rval) {
        CClp_free_warmstart (w);
    }

    return rval;
}

static void init_warmstart (CClp_warmstart *w)
{
    if (w) {
        w->nstruct = 0;
        w->nrows   = 0;
        w->cstat   = (char *) NULL;
        w->rstat   = (char *) NULL;
        w->dnorm   = (double *) NULL;
    }
}

void CClp_free_warmstart (CClp_warmstart **w)
{
    if (w && (*w)) {
        CC_IFFREE ((*w)->cstat, char);
        CC_IFFREE ((*w)->rstat, char);
        CC_IFFREE ((*w)->dnorm, double);
        CC_FREE (*w, CClp_warmstart);

    }
}

int CClp_sread_warmstart (CC_SFILE *f, CClp_warmstart **w)
{
    int rval = 0;
    char name[5];
    int i, k, nstruct, nrows, has_dnorms;

    CClp_free_warmstart (w);

    for (i=0; i<4; i++) {
        if (CCutil_sread_char (f, &name[i])) {
            rval = 1;  goto CLEANUP;
        }
    }
    name[4] = '\0';

    if (strncmp (name, SOLVER_WARMSTART_NAME, 4) &&
        strncmp (name, "CPL5", 4)) {
        fprintf (stderr, "warmstart for another solver (%s) ignored\n", name);
        goto CLEANUP;
    }

    if (CCutil_sread_int (f, &nstruct)) goto CLEANUP;
    if (CCutil_sread_int (f, &nrows)) goto CLEANUP;

    (*w) = CC_SAFE_MALLOC (1, CClp_warmstart);
    if ((*w) == (CClp_warmstart *) NULL) {
        fprintf (stderr, "Out of memory in CClp_sread_warmstart\n");
        rval = 1; goto CLEANUP;
    }
    init_warmstart (*w);

    (*w)->cstat = CC_SAFE_MALLOC (nstruct, char);
    (*w)->rstat = CC_SAFE_MALLOC (nrows, char);
    if ((*w)->cstat == (char *) NULL ||
        (*w)->rstat == (char *) NULL) {
        fprintf (stderr, "out of memory in CClp_sread_warmstart\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < nstruct; i++) {
        if (CCutil_sread_bits (f, &k, 2)) {
            rval = 1; goto CLEANUP;
        }
        (*w)->cstat[i] = '0' + k;
    }
    for (i = 0; i < nrows; i++) {
        if (CCutil_sread_bits (f, &k, 1)) {
            rval = 1; goto CLEANUP;
        }
        (*w)->rstat[i] = '0' + k;
    }

    if (CCutil_sread_int (f, &has_dnorms)) { rval = 1; goto CLEANUP; }

    if (has_dnorms) {
        (*w)->dnorm = CC_SAFE_MALLOC (nrows, double);
        if ((*w)->dnorm == (double *) NULL) {
            fprintf (stderr, "out of memory in CClp_sread_warmstart\n");
            rval = 1; goto CLEANUP;
        }
        for (i = 0; i < nrows; i++) {
            if (CCutil_sread_double (f, &(((*w)->dnorm)[i]))) {
                rval = 1; goto CLEANUP;
            }
        }
    }

    (*w)->nstruct = nstruct;
    (*w)->nrows   = nrows;

CLEANUP:

    if (rval) CClp_free_warmstart (w);
    return rval;
}

int CClp_swrite_warmstart (CC_SFILE *f, CClp_warmstart *w)
{
    int i, rval = 0;
    const char *name = SOLVER_WARMSTART_NAME;

    for (i = 0; i < 4; i++) {
        if (CCutil_swrite_char (f, name[i])) { rval = 1; goto CLEANUP; }
    }

    if (CCutil_swrite_int (f, w->nstruct)) { rval = 1; goto CLEANUP; }
    if (CCutil_swrite_int (f, w->nrows)) { rval = 1; goto CLEANUP; }

    for (i = 0; i < w->nstruct; i++) {
        if (CCutil_swrite_bits (f, w->cstat[i], 2)) { rval = 1; goto CLEANUP; }
    }

    for (i = 0; i < w->nrows; i++) {
        if (CCutil_swrite_bits (f, w->rstat[i], 1)) { rval = 1; goto CLEANUP; }
    }

    if (w->dnorm == (double *) NULL) {
        if (CCutil_swrite_int (f, 0)) { rval = 1; goto CLEANUP; }
    } else {
        if (CCutil_swrite_int (f, 1)) { rval = 1; goto CLEANUP; }
        for (i = 0; i < w->nrows; i++) {
            if (CCutil_swrite_double (f, w->dnorm[i])) {
                rval = 1; goto CLEANUP;
            }
        }
    }

CLEANUP:

    return rval;
}

int CClp_get_info (CClp *lp, CClp_info **i)
{
    int rval = 0;
    int k;

    CClp_free_info (i);

    (*i) = CC_SAFE_MALLOC (1, CClp_info);
    if ((*i) == (CClp_info *) NULL) {
        fprintf (stderr, "Out of memory in CClp_get_info\n");
        rval = 1; goto CLEANUP;
    }

    (*i)->nstruct = 0;
    (*i)->nrows  = 0;
    (*i)->cstat = (char *) NULL;
    (*i)->rstat = (char *) NULL;


    k = QSget_colcount (lp->p);
    (*i)->nstruct = k;

    k = QSget_rowcount (lp->p);
    (*i)->nrows = k;


    (*i)->cstat = CC_SAFE_MALLOC ((*i)->nstruct, char);
    (*i)->rstat = CC_SAFE_MALLOC ((*i)->nrows, char);
    if (!(*i)->cstat || !(*i)->rstat) {
        fprintf (stderr, "out of memory in CClp_get_info\n");
        rval = 1; goto CLEANUP;
    }

    rval = QSget_basis_array (lp->p, (*i)->cstat, (*i)->rstat);
    if (rval) {
        fprintf (stderr, "QSget_basis_array failed\n");
        goto CLEANUP;
    }

CLEANUP:

    if (rval) {
        CClp_free_info (i);
    }
    return rval;
}

int CClp_create_info (CClp_info **i, int rcount, int ccount)
{
    int rval = 0;
    int j;

    CClp_free_info (i);

    (*i) = CC_SAFE_MALLOC (1, CClp_info);
    if ((*i) == (CClp_info *) NULL) {
        fprintf (stderr, "Out of memory in CClp_create_info\n");
        rval = 1; goto CLEANUP;
    }

    (*i)->cstat = (char *) NULL;
    (*i)->rstat = (char *) NULL;

    (*i)->nstruct = ccount;
    if (ccount == 0) {
        fprintf (stderr, "No columns in CClp_create_info\n");
        rval = 1; goto CLEANUP;
    }
    (*i)->nrows = rcount;
    if (rcount == 0) {
        fprintf (stderr, "No rows in CClp_create_info\n");
        rval = 1; goto CLEANUP;
    }

    (*i)->cstat = CC_SAFE_MALLOC ((*i)->nstruct, char);
    (*i)->rstat = CC_SAFE_MALLOC ((*i)->nrows, char);
    if (!(*i)->cstat || !(*i)->rstat) {
        fprintf (stderr, "out of memory in CClp_create_info\n");
        rval = 1; goto CLEANUP;
    }

    for (j = 0; j < ccount; j++) {
        (*i)->cstat[j] = 0;
    }
    for (j = 0; j < rcount; j++) {
        (*i)->rstat[j] = 0;
    }

CLEANUP:

    if (rval) {
        CClp_free_info (i);
    }

    return rval;
}

int CClp_is_col_active (CClp_info *i, int c)
{
    if (c < 0 || c >= i->nstruct) return 0;
    return i->cstat[c] == QS_COL_BSTAT_BASIC ||
           i->cstat[c] == QS_COL_BSTAT_UPPER;
}

int CClp_is_row_active (CClp_info *i, int r)
{
    if (r < 0 || r >= i->nrows) return 0;
    return i->rstat[r] == QS_ROW_BSTAT_LOWER;
}

void CClp_set_col_active (CClp_info *i, int c)
{
    if (c >= 0 && c < i->nstruct) i->cstat[c] = QS_COL_BSTAT_BASIC; 
}

void CClp_set_col_inactive (CClp_info *i, int c)
{
    if (c >= 0 && c < i->nstruct) i->cstat[c] = QS_COL_BSTAT_LOWER; 
}

void CClp_set_col_upper (CClp_info *i, int c)
{
    if (c >= 0 && c < i->nstruct) i->cstat[c] = QS_COL_BSTAT_UPPER; 
}

void CClp_set_row_active (CClp_info *i, int r)
{
    if (r >= 0 && r < i->nrows) i->rstat[r] = QS_ROW_BSTAT_LOWER;
}

void CClp_set_row_inactive (CClp_info *i, int r)
{
    if (r >= 0 && r < i->nrows) i->rstat[r] = QS_ROW_BSTAT_BASIC;
}

void CClp_free_info (CClp_info **i)
{
    if (i && (*i)) {
        CC_IFFREE ((*i)->cstat, char);
        CC_IFFREE ((*i)->rstat, char);
        CC_FREE (*i, CClp_info);
    }
}

int CClp_x (CClp *lp, double *x)
{
    int rval = 0;

    rval = QSget_x_array (lp->p, x);
    if (rval) {
        fprintf (stderr, "QSget_x_array\n"); goto CLEANUP;
    }

CLEANUP:

    return rval;
}

int CClp_rc (CClp *lp, double *rc)
{
    int rval = 0;

    rval = QSget_rc_array (lp->p, rc);
    if (rval) {
        fprintf (stderr, "QSget_rc_array failed"); goto CLEANUP;
    }

CLEANUP:

    return rval;
}

int CClp_pi (CClp *lp, double *pi)
{
    int rval = 0;
    int status;

    rval = QSget_status (lp->p, &status);
    if (rval) {
        fprintf (stderr, "QSget_status failed"); goto CLEANUP;
    }

    if (status == QS_LP_INFEASIBLE) {
        rval = QSget_infeas_array (lp->p, pi);
        if (rval) {
            fprintf (stderr, "QSget_infeas_array failed"); goto CLEANUP;
        }
    } else {
        rval = QSget_pi_array (lp->p, pi);
        if (rval) {
            fprintf (stderr, "QSget_pi_array failed"); goto CLEANUP;
        }
    }

CLEANUP:

    return rval;
}

int CClp_objval (CClp *lp, double *obj)
{
    int rval = 0;

    rval = QSget_objval (lp->p, obj);
    if (rval) {
        fprintf (stderr, "QSget_objval failed");
        goto CLEANUP;
    }

CLEANUP:

    return rval;
}

int CClp_nrows (CClp *lp)
{
    int k;

    k = QSget_rowcount (lp->p);
    if (k == 0) {
        fprintf (stderr, "QSget_rowcount failed - continue anyway\n");
    }

    return k;
}

int CClp_ncols (CClp *lp)
{
    int k;

    k = QSget_colcount (lp->p);
    if (k == 0) {
        fprintf (stderr, "QSget_colcount failed - continue anyway\n");
    }

    return k;
}

int CClp_nnonzeros (CClp *lp)
{
    int k;

    k = QSget_nzcount (lp->p);
    if (k == 0) {
        fprintf (stderr, "QSget_nzcount failed - continue anyway\n");
    }

    return k;
}

int CClp_status (CClp *lp, int *status)
{
    printf ("CClp_status ...\n"); fflush (stdout);

    if (lp || status) {
        lp_message (); return 1;
    } else {
        return 0;
    }
}

int CClp_getweight (CClp *lp, int nrows, int *rmatbeg, int *rmatind,
                    double *rmatval, double *weight)
{
    printf ("CClp_getweight ...\n"); fflush (stdout);

    if (lp || nrows || rmatbeg || rmatind || rmatval || weight) {
        lp_message (); return 1;
    } else {
        return 0;
    }
}

int CClp_dump_lp (CClp *lp, const char *fname)
{
    int rval = 0;

    rval = QSwrite_prob (lp->p, fname, "LP");
    if (rval) {
        fprintf (stderr, "QSwrite_prob failed\n");
        goto CLEANUP;
    }

CLEANUP:

    return rval;
}

int CClp_getgoodlist (CClp *lp, int *goodlist, int *goodlen_p,
        double *downpen, double *uppen)
{
    /* Replace by something looking for non-degenerate pivots */

    int rval = 0;
    int i, j, ncols = 0;
    double *x = (double *) NULL;

    *goodlen_p = 0;

    ncols = QSget_colcount (lp->p);

    x = CC_SAFE_MALLOC (ncols, double);
    if (!x) {
        fprintf (stderr, "out of memory in CClp_getgoodlist\n");
        rval = 1; goto CLEANUP;
    }

    rval = QSget_x_array (lp->p, x);
    if (rval) {
        fprintf (stderr, "QSget_x_array\n"); goto CLEANUP;
    }

    for (i = 0, j = 0; i < ncols; i++) {
        if (x[i] >= CCtsp_INTTOL && x[i] <= 1.0 - CCtsp_INTTOL) {
            goodlist[j] = i;
            downpen[j]  = x[i];
            uppen[j]    = 1.0 - x[i];
            j++;
        }
    }

    *goodlen_p = j;
    
CLEANUP:

    CC_IFFREE (x, double);
    return rval;
}

int CClp_strongbranch (CClp *lp, int *candidatelist, int ncand,
        double *downpen, double *uppen, int iterations, double upperbound)
{
    int rval = 0;

    /* Hack until we understand what is happening */

    iterations = (iterations * 0.25) + 1;

    rval = QSopt_strongbranch (lp->p, ncand, candidatelist, (double *) NULL,
                               downpen, uppen, iterations, upperbound);
    if (rval) {
        fprintf (stderr, "QSopt_strongbranch failed\n"); goto CLEANUP;
    }


CLEANUP:
    
    return rval;
}
