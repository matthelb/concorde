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
/*              Interface Routines to an LP Solver                          */
/*                                                                          */
/*                                                                          */
/*  NOTE: These are just dummy routines. To make use of the LP-based        */
/*   portion of concorde you will need to write an interface between these  */
/*   routines and an LP solver.  This file documents the interface          */
/*   routines that are needed, and provides dummy routines so that the      */
/*   concorde can be compiled without an LP solver.                         */
/*                                                                          */
/*   To use the CPLEX 5.0 library, the code in lpcplex5.c can be used       */
/*   instead of this file.  To use the CPLEX 4.0 library, the code in       */
/*   lpcplex4.c can be used instead of this file.                           */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: April 17, 1997                                                    */
/*        June 19, 1997 (bico, REB)                                         */
/*        January 13, 2002 (bico)                                           */
/*                                                                          */
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
/*                                                                          */
/****************************************************************************/



#include "machdefs.h"
#include "util.h"
#include "macrorus.h"
#include "lp.h"

struct CClp {
    int dummy;
};

struct CClp_warmstart {
    int dummy;
};

struct CClp_info {
    int dummy;
};


static void
    lp_message (void);



static void lp_message (void)
{
    fprintf (stderr, "need to link an lp solver to use this function\n");
}

int CClp_init (CClp **lp)
{
    if (lp && *lp) {
        lp_message (); return 1;
    } else {
        return 0;
    }
}

int CClp_force_perturb (CClp *lp)
{
    if (lp) {
        lp_message (); return 1;
    } else {
        return 0;
    }
}

int CClp_tune_small (CClp *lp)
{
    if (lp) {
        lp_message (); return 1;
    } else {
        return 0;
    }
}

int CClp_disable_presolve (CClp *lp)
{
    if (lp) {
        lp_message (); return 1;
    } else {
        return 0;
    }
}

void CClp_free (CClp **lp)
{
    if (lp && (*lp)) lp_message ();
}

void CClp_freelp (CClp **lp)
{
    if (lp && (*lp)) lp_message ();
}

int CClp_loadlp (CClp *lp, const char *name, int ncols, int nrows,
        int objsense, double *obj, double *rhs, char *sense, int *matbeg,
        int *matcnt, int *matind, double *matval, double *lb, double *ub)

{
    if (lp || name || ncols || nrows || objsense || obj || rhs || sense
           || matbeg || matcnt || matind || matval || lb || ub) {
        lp_message (); return 1;
    } else {
        return 0;
    }
}

int CClp_create (CClp *lp, const char *name)
{
    if (lp || name) {
        lp_message (); return 1;
    } else {
        return 0;
    }
}

int CClp_new_row (CClp *lp, char sense, double rhs)
{
    if (lp || sense || rhs) {
        lp_message (); return 1;
    } else {
        return 0;
    }
}

int CClp_change_sense (CClp *lp, int row, char sense)
{
    if (lp || row || sense) {
        lp_message (); return 1;
    } else {
        return 0;
    }
}

int CClp_opt (CClp *lp, int method)
{
    if (lp || method) {
        lp_message (); return 1;
    } else {
        return 0;
    }
}

int CClp_limited_dualopt (CClp *lp, int iterationlim, int *status,
        double *objupperlim)
{
    if (lp || iterationlim || status || objupperlim) {
        lp_message (); return 1;
    } else {
        return 0;
    }
}

int CClp_addrows (CClp *lp, int newrows, int newnz, double *rhs, char *sense,
                  int *rmatbeg, int *rmatind, double *rmatval)
{
    if (lp || newrows || newnz || rhs || sense || rmatbeg || rmatind
           || rmatval) {
        lp_message (); return 1;
    } else {
        return 0;
    }
}

int CClp_addcols (CClp *lp, int newcols, int newnz, double *obj,
                  int *cmatbeg, int *cmatind, double *cmatval,
                  double *lb, double *ub)
{
    if (lp || newcols || newnz || obj || cmatbeg || cmatind || cmatval
           || lb || ub) {
        lp_message (); return 1;
    } else {
        return 0;
    }
}

int CClp_delete_row (CClp *lp, int i)
{
    if (lp || i) {
        lp_message (); return 1;
    } else {
        return 0;
    }
}

int CClp_delete_set_of_rows (CClp *lp, int *delstat)
{
    if (lp || delstat) {
        lp_message (); return 1;
    } else {
        return 0;
    }
}

int CClp_delete_column (CClp *lp, int i)
{
    if (lp || i) {
        lp_message (); return 1;
    } else {
        return 0;
    }
}

int CClp_delete_set_of_columns (CClp *lp, int *delstat)
{
    if (lp || delstat) {
        lp_message (); return 1;
    } else {
        return 0;
    }
}

int CClp_setbnd (CClp *lp, int col, char lower_or_upper, double bnd)
{
    if (lp || col || lower_or_upper || bnd) {
        lp_message (); return 1;
    } else {
        return 0;
    }
}

int CClp_get_warmstart (CClp *lp, CClp_warmstart **w)
{
    if (lp || w) {
        lp_message (); return 1;
    } else {
        return 0;
    }
}

int CClp_load_warmstart (CClp *lp, CClp_warmstart *w)
{
    if (lp || w) {
        lp_message (); return 1;
    } else {
        return 0;
    }
}

int CClp_build_warmstart (CClp_warmstart **w, CClp_info *i)
{
    if (w || i) {
        lp_message (); return 1;
    } else {
        return 0;
    }
}

void CClp_free_warmstart (CClp_warmstart **w)
{
    if (w) lp_message ();
}

int CClp_sread_warmstart (CC_SFILE *f, CClp_warmstart **w)
{
    if (f || w) {
        lp_message (); return 1;
    } else {
        return 0;
    }
}

int CClp_swrite_warmstart (CC_SFILE *f, CClp_warmstart *w)
{
    if (f || w) {
        lp_message (); return 1;
    } else {
        return 0;
    }
}

int CClp_get_info (CClp *lp, CClp_info **i)
{
    if (lp || i) {
        lp_message (); return 1;
    } else {
        return 0;
    }
}

int CClp_create_info (CClp_info **i, int rcount, int ccount)
{
    if (i || rcount || ccount) {
        lp_message (); return 1;
    } else {
        return 0;
    }
}

int CClp_is_col_active (CClp_info *i, int c)
{
    if (i || c) {
        lp_message (); return 1;
    } else {
        return 0;
    }
}

int CClp_is_row_active (CClp_info *i, int r)
{
    if (i || r) {
        lp_message (); return 1;
    } else {
        return 0;
    }
}

void CClp_set_col_active (CClp_info *i, int c)
{
    if (i || c) lp_message ();
}

void CClp_set_col_inactive (CClp_info *i, int c)
{
    if (i || c) lp_message ();
}

void CClp_set_col_upper (CClp_info *i, int c)
{
    if (i || c) lp_message ();
}

void CClp_set_row_active (CClp_info *i, int r)
{
    if (i || r) lp_message ();
}

void CClp_set_row_inactive (CClp_info *i, int r)
{
    if (i || r) lp_message ();
}

void CClp_free_info (CClp_info **i)
{
    if (i) lp_message ();
}

int CClp_x (CClp *lp, double *x)
{
    if (lp || x) {
        lp_message (); return 1;
    } else {
        return 0;
    }
}

int CClp_rc (CClp *lp, double *rc)
{
    if (lp || rc) {
        lp_message (); return 1;
    } else {
        return 0;
    }
}

int CClp_pi (CClp *lp, double *pi)
{
    if (lp || pi) {
        lp_message (); return 1;
    } else {
        return 0;
    }
}

int CClp_objval (CClp *lp, double *obj)
{
    if (lp || obj) {
        lp_message (); return 1;
    } else {
        return 0;
    }
}

int CClp_nrows (CClp *lp)
{
    if (lp) {
        lp_message (); return 1;
    } else {
        return 0;
    }
}

int CClp_ncols (CClp *lp)
{
    if (lp) {
        lp_message (); return 1;
    } else {
        return 0;
    }
}

int CClp_nnonzeros (CClp *lp)
{
    if (lp) {
        lp_message (); return 1;
    } else {
        return 0;
    }
}

int CClp_status (CClp *lp, int *status)
{
    if (lp || status) {
        lp_message (); return 1;
    } else {
        return 0;
    }
}

int CClp_getweight (CClp *lp, int nrows, int *rmatbeg, int *rmatind,
                    double *rmatval, double *weight)
{
    if (lp || nrows || rmatbeg || rmatind || rmatval || weight) {
        lp_message (); return 1;
    } else {
        return 0;
    }
}

int CClp_dump_lp (CClp *lp, const char *fname)
{
    if (lp || fname ) {
        lp_message (); return 1;
    } else {
        return 0;
    }
}

int CClp_getgoodlist (CClp *lp, int *goodlist, int *goodlen_p,
        double *downpen, double *uppen)
{
    if (lp || goodlist || goodlen_p || downpen || uppen) {
        lp_message (); return 1;
    } else {
        return 0;
    }
}

int CClp_strongbranch (CClp *lp, int *candidatelist, int ncand,
        double *downpen, double *uppen, int iterations, double upperbound)
{
    if (lp || candidatelist || ncand || downpen || uppen || iterations
           || upperbound) {
        lp_message (); return 1;
    } else {
        return 0;
    }
}
