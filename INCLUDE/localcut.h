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

#ifndef __LOCALCUT_H_
#define __LOCALCUT_H_

#include "tsp.h"
#include "util.h"
#include "lp.h"

typedef struct CCchunk_flag {
    unsigned dummy : 1;
    unsigned permute : 1;
    unsigned weighted : 1;
    unsigned spheres : 1;
    unsigned uncivilized : 1;
    unsigned noshrink : 1;
    unsigned nolift : 2;
    unsigned maxchunksize : 8;
    unsigned spheresize : 8;
} CCchunk_flag;

typedef struct CCchunk_find_timer {
    CCutil_timer shrink;
    CCutil_timer locate;
    CCutil_timer all;
} CCchunk_find_timer;

typedef struct CCchunk_oracle_timer {
    CCutil_timer bnbtsp;
    CCutil_timer tinytsp;
    CCutil_timer all;
} CCchunk_oracle_timer;

typedef struct CCchunk_separate_timer {
    CCutil_timer lpsolver;
    CCchunk_oracle_timer oracle;
    CCutil_timer all;
} CCchunk_separate_timer;

typedef struct CCchunk_lift_timer {
    CCutil_timer all;
    CCutil_timer liberate_equality;
    CCutil_timer liberate_fixed;
    CCutil_timer strengthen_edges;
    CCutil_timer strengthen_equality;
    CCutil_timer strengthen_work;
    CCutil_timer decompose;
    CCutil_timer tilt_oracle;
    CCutil_timer liberate_oracle;
    CCutil_timer verify_oracle;
    CCchunk_oracle_timer oracle;
} CCchunk_lift_timer;

typedef struct CCchunk_localcut_timer {
    CCchunk_find_timer     find;
    CCchunk_separate_timer separate;
    CCchunk_lift_timer     lift;
    CCutil_timer   all;
} CCchunk_localcut_timer;


int
    CCchunk_localcuts (CCtsp_lpcut_in **clist, int *cutcount, int ncount,
        int ecount, int *elist, double *x, double eps, CCchunk_flag flags,
        CCchunk_localcut_timer *timer, int silent, CCrandstate *rstate);

void
    CCchunk_init_separate_timer (CCchunk_separate_timer *timer),
    CCchunk_init_find_timer (CCchunk_find_timer *timer),
    CCchunk_init_lift_timer (CCchunk_lift_timer *timer),
    CCchunk_init_oracle_timer (CCchunk_oracle_timer *timer),
    CCchunk_init_localcut_timer (CCchunk_localcut_timer *timer),
    CCchunk_print_separate_timer (CCchunk_separate_timer *timer),
    CCchunk_print_find_timer (CCchunk_find_timer *timer),
    CCchunk_print_lift_timer (CCchunk_lift_timer *timer),
    CCchunk_print_oracle_timer (CCchunk_oracle_timer *timer),
    CCchunk_print_localcut_timer (CCchunk_localcut_timer *timer);

typedef struct CCchunklp {
   CClp       *lp;
   int        *active;
   int        *cmatind;
   double     *cmatval;
   double     *pi;
   int        nrows;
   int        extracols;
} CCchunklp;

#define CC_CHUNK_LPFEASIBLE 0
#define CC_CHUNK_LPINFEASIBLE 1

int
    CCchunk_lpinit     (CCchunklp **lp_p, const char *lp_name, int lp_nrows,
                        double *xstar),
    CCchunk_lpaddcol   (CCchunklp *lp, double *x),
    CCchunk_lprelaxrow (CCchunklp *lp, int del_row),
    CCchunk_lpsolve    (CCchunklp *lp, int *lpstatus_p, double *c,
                        double *alpha_p),
    CCchunk_lpbasis    (CCchunklp *lp, int ncols, int *basis);

void
    CCchunk_lpfree     (CCchunklp **lp_p);


typedef struct CCchunk_graph {
    int                 ncount;
    int                 ecount;
    int                *end0;
    int                *end1;
    int                *fixed;
    double             *weight;
    int                *equality;
    int               **members;
}   CCchunk_graph;

typedef struct CCchunk_ineq {
    int *coef;
    int rhs;
} CCchunk_ineq;

typedef struct CCchunk_fault {
    CCchunk_ineq a;
    int nsols;
    int *sols;
} CCchunk_fault;

typedef struct CCchunk_chunk_callback {
    int (*func) (CCchunk_graph *chunk, int *faulty, void *u_data);
    void *u_data;
} CCchunk_chunk_callback;

typedef struct CCchunk_fault_callback {
    int (*func) (CCchunk_graph *chunk, CCchunk_fault *fault, int *finished,
        void *u_data);
    void *u_data;
} CCchunk_fault_callback;

typedef struct CCchunk_cut_callback {
    int (*begin_cut) (void *u_data);
    int (*add_clique) (int *arr, int size, void *u_data);
    int (*abort_cut) (void *u_data);
    int (*finish_cut) (int rhs, int *finished, void *u_data);
    void *u_data;
} CCchunk_cut_callback;

#define CC_CHUNK_ORACLE_ERROR (-1)
#define CC_CHUNK_ORACLE_SEARCHLIMITEXCEEDED (1)
#define CC_CHUNK_ORACLE_INFEASIBLE (2)


int
    CCchunk_finder (int ncount, int ecount, int *elist, double *elen,
        double eps, CCchunk_flag flags, CCchunk_find_timer *timer,
        CCchunk_chunk_callback *callback, CCrandstate *rstate),
    CCchunk_separate (CCchunk_graph *chunk, CCchunk_separate_timer *timer,
        CCchunk_fault_callback *callback),
    CCchunk_lift (CCchunk_graph *chunk, CCchunk_fault *fault,
        CCchunk_lift_timer *timer,
        CCchunk_cut_callback *callback),
    CCchunk_ineq_to_lpcut_in (int nnodes, int ecount, int *elist, int *ecoef,
        int rhs, CCtsp_lpcut_in *c),
    CCchunk_ineq_to_cut (int nnodes, int ecount, int *elist, int *ecoef,
        int rhs, int outside, CCchunk_cut_callback *callback),
    CCchunk_oracle (CCchunk_graph *ch, CCchunk_ineq *c, int *xsol, int *objval,
        int rhsvalid, int effort_limit, CCchunk_oracle_timer *timer),
    CCchunk_verify (CCchunk_graph *ch, CCchunk_ineq *c);

CCchunk_graph
   *CCchunk_graph_alloc (int ncount, int ecount);

void
    CCchunk_graph_free (CCchunk_graph *c);


#define CC_CHUNK_INTMAT_NOORTHO (1)
#define CC_CHUNK_INTMAT_MEMORY (-1)
#define CC_CHUNK_INTMAT_OVERFLOW_M (-2)
#define CC_CHUNK_INTMAT_OVERFLOW_A (-3)
#define CC_CHUNK_INTMAT_OVERFLOW_S (-4)
#define CC_CHUNK_INTMAT_ERROR (-5)

typedef long int CCmatval;
#define CC_MATVAL_MAX ((((CCmatval) 1) << (sizeof (CCmatval) * 8 - 2))|((((CCmatval) 1) << (sizeof (CCmatval) * 8 - 2))-1))

typedef struct CCchunk_intmat {
    int  *matrix;
    CCmatval  *factor;
    int  *csize;
    int  *rperm;
    int  *cperm;
    CCmatval *x;
    CCmatval *best_x;
    int  nrows;
    int  ncols;
    int  rowspace;
} CCchunk_intmat;

/* Exported functions */


int
    CCchunk_intmat_build (CCchunk_intmat *mat_p, int ncols),
    CCchunk_intmat_addrow (CCchunk_intmat *mat_p, int *row),
    CCchunk_intmat_ortho (CCchunk_intmat *mat_p, int *ortho,
        int *pcol_p, int *taboo);

void
    CCchunk_intmat_init (CCchunk_intmat *mat_p),
    CCchunk_intmat_free (CCchunk_intmat *mat_p),
    CCchunk_intmat_dellastrows (CCchunk_intmat *mat_p, int ndel);


#endif /* __LOCALCUT_H_ */
