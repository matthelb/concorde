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

#ifndef __LP_H
#define __LP_H

#include "util.h"

#define  CClp_METHOD_DUAL    1
#define  CClp_METHOD_PRIMAL  2
#define  CClp_METHOD_BARRIER 3

#define  CClp_SUCCESS        0
#define  CClp_FAILURE        1
#define  CClp_UNBOUNDED      2
#define  CClp_INFEASIBLE     3
#define  CClp_UNKNOWN        4

struct CClp;
typedef struct CClp CClp;

struct CClp_warmstart;
typedef struct CClp_warmstart CClp_warmstart;

struct CClp_info;
typedef struct CClp_info CClp_info;


int
    CClp_init (CClp **lp),
    CClp_force_perturb (CClp *lp),
    CClp_tune_small (CClp *lp),
    CClp_disable_presolve (CClp *lp),
    CClp_loadlp (CClp *lp, const char *name, int ncols, int nrows,
        int objsense, double *obj, double *rhs, char *sense, int *matbeg,
        int *matcnt, int *matind, double *matval, double *lb, double *ub),
    CClp_create (CClp *lp, const char *name),
    CClp_new_row (CClp *lp, char sense, double rhs),
    CClp_change_sense (CClp *lp, int row, char sense),
    CClp_opt (CClp *lp, int method),
    CClp_limited_dualopt (CClp *lp, int lim, int *status, double *upperbound),
    CClp_addrows (CClp *lp, int newrows, int newnz, double *rhs, char *sense,
            int *rmatbeg, int *rmatind, double *rmatval),
    CClp_addcols (CClp *lp, int newcols, int newnz, double *obj,
            int *cmatbeg, int *cmatind, double *cmatval, double *lb,
            double *ub),
    CClp_delete_row (CClp *lp, int i),
    CClp_delete_set_of_rows (CClp *lp, int *delstat),
    CClp_delete_column (CClp *lp, int i),
    CClp_delete_set_of_columns (CClp *lp, int *delstat),
    CClp_setbnd (CClp *lp, int col, char lower_or_upper, double bnd),
    CClp_get_warmstart (CClp *lp, CClp_warmstart **w),
    CClp_load_warmstart (CClp *lp, CClp_warmstart *w),
    CClp_build_warmstart (CClp_warmstart **w, CClp_info *i),
    CClp_sread_warmstart (CC_SFILE *f, CClp_warmstart **w),
    CClp_swrite_warmstart (CC_SFILE *f, CClp_warmstart *w),
    CClp_get_info (CClp *lp, CClp_info **i),
    CClp_create_info (CClp_info **i, int rcount, int ccount),
    CClp_is_col_active (CClp_info *i, int c),
    CClp_is_row_active (CClp_info *i, int c),
    CClp_x (CClp *lp, double *x),
    CClp_rc (CClp *lp, double *rc),
    CClp_pi (CClp *lp, double *pi),
    CClp_objval (CClp *lp, double *obj),
    CClp_nrows (CClp *lp),
    CClp_ncols (CClp *lp),
    CClp_nnonzeros (CClp *lp),
    CClp_status (CClp *lp, int *status),
    CClp_getweight (CClp *lp, int nrows, int *rmatbeg, int *rmatind,
            double *rmatval, double *weight),
    CClp_dump_lp (CClp *lp, const char *fname),
    CClp_getgoodlist (CClp *lp, int *goodlist, int *goodlen_p,
            double *downpen, double *uppen),
    CClp_strongbranch (CClp *lp, int *candidatelist, int ncand,
            double *downpen, double *uppen, int iterations,
            double upperbound);

void
    CClp_free (CClp **lp),
    CClp_freelp (CClp **lp),
    CClp_free_warmstart (CClp_warmstart **w),
    CClp_set_col_active (CClp_info *i, int c),
    CClp_set_col_inactive (CClp_info *i, int c),
    CClp_set_col_upper (CClp_info *i, int c),
    CClp_set_row_active (CClp_info *i, int r),
    CClp_set_row_inactive (CClp_info *i, int r),
    CClp_free_info (CClp_info **i),
    CClp_pivotout (CClp *lp, int j),
    CClp_pivotin (CClp *lp, int i);



#endif  /* __LP_H */
