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
/*              Interface Routines to The CPLEX 4.0 LP Solver               */
/*                                                                          */
/*  NOTE: Use this code in place of lp_none.c to access the Cplex 4.0       */
/*   library. You will also need to link the cplex library via the          */
/*   makefile.                                                              */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: April 17, 1997                                                    */
/*        June 19, 1997 (bico, REB)                                         */
/*        September 14, 1997 (REB)                                          */
/*                                                                          */
/*  See the start of the file lp_none.c for a description of the            */
/*  functions exported by this file.                                        */
/*                                                                          */
/*                                                                          */
/*  NOTES:                                                                  */
/*                                                                          */
/****************************************************************************/


#include "machdefs.h"
#include "util.h"
#include "macrorus.h"
#include "lp.h"
#include <cplex.h>

typedef struct CCcplex40_data {
    char   *probname;
    double *obj;
    double *rhs;
    char   *sense;
    int    *matbeg;
    int    *matcnt;
    int    *matind;
    double *matval;
    double *lb;
    double *ub;
    int    triv_row_ind;
    int    triv_col_ind;
} CCcplex40_data;

struct CClp {
    CPXENVptr       cplex_env;
    CPXLPptr        cplex_lp;
    CCcplex40_data *cplex40_data;
};

struct CClp_warmstart {
    int      rcount;
    int      ccount;
    int     *rstat;
    int     *cstat;
    double  *dnorm;
    int     *bhead;
    int      dlen;
};

struct CClp_info {
    int  rcount;
    int  ccount;
    int *rstat;
    int *cstat;
};

#define SOLVER_WARMSTART_NAME "CPL4"


static int
    init_cplex40_data (CCcplex40_data **data_p, int ncols, int nrows,
        const char *name, double *obj, double *rhs, char *sense, int *matbeg,
        int *matcnt, int *matind, double *matval, double *lb, double *ub,
        int *cspace_p, int *rspace_p, int *nzspace_p),
    reallocprob (CClp *lp, int cneeded, int rneeded, int nzneeded),
    getspace (CClp *lp, int *cspace_p, int *rspace_p, int *nzspace_p),
    getsurplus (CClp *lp, int *csurplus_p, int *rsurplus_p, int *nzsurplus_p),
    primalopt (CClp *lp),
    dualopt (CClp *lp),
    baropt (CClp *lp),
    getfarkasmultipliers (CClp *lp, double *y);

static void
    free_cplex40_data (CCcplex40_data **data_p);


#undef CC_CPLEX_DISPLAY

int CClp_init (CClp **lp)
{
    int rval = 0;

    CClp_free (lp);

    (*lp) = CC_SAFE_MALLOC (1, CClp);
    if ((*lp) == (CClp *) NULL) {
        fprintf (stderr, "Out of memory in CClp_init\n");
        rval = 1; goto CLEANUP;
    }

    (*lp)->cplex_env    = (CPXENVptr) NULL;
    (*lp)->cplex_lp     = (CPXLPptr) NULL;
    (*lp)->cplex40_data = (CCcplex40_data *) NULL;

    (*lp)->cplex_env = CPXopenCPLEX (&rval);
    if (rval) {
        fprintf (stderr, "CPXopenCPLEX failed\n"); goto CLEANUP;
    }

#ifdef CC_CPLEX_DISPLAY
    /* the documentation doesn't say what the return value means */
    rval = CPXsetintparam ((*lp)->cplex_env, CPX_PARAM_SCRIND, 1);
    if (rval) {
        fprintf (stderr, "CPXsetintparam CPX_PARAM_SCRIND failed\n");
        goto CLEANUP;
    }

    rval = CPXsetintparam ((*lp)->cplex_env, CPX_PARAM_SIMDISPLAY, 1);
    if (rval) {
        fprintf (stderr, "CPXsetintparam CPX_PARAM_SIMDISPLAY failed\n");
        goto CLEANUP;
    }
#endif

    rval = CPXsetintparam ((*lp)->cplex_env, CPX_PARAM_ADVIND, 1);
    if (rval) {
        fprintf (stderr, "CPXsetintparam CPX_PARAM_ADVIND failed\n");
        goto CLEANUP;
    }
    rval = CPXsetintparam ((*lp)->cplex_env, CPX_PARAM_DPRIIND,
                           CPX_DPRIIND_STEEP);
    if (rval) {
        fprintf (stderr, "CPXsetintparam CPX_PARAM_DPRIIND failed\n");
        goto CLEANUP;
    }
    rval = CPXsetintparam ((*lp)->cplex_env, CPX_PARAM_PPRIIND,
                           CPX_PPRIIND_STEEP);
    if (rval) {
        fprintf (stderr, "CPXsetintparam CPX_PARAM_PPRIIND failed\n");
        goto CLEANUP;
    }

    /* REB, 14 October 1997:  The following three parameter settings 
       help fl3795 a bunch, and are probably not a bad idea in general */

    rval = CPXsetdblparam ((*lp)->cplex_env, CPX_PARAM_EPPER, 1.0E-6);
    if (rval) {
        fprintf (stderr, "CPXsetdblparam CPX_PARAM_EPPER failed\n");
        goto CLEANUP;
    }
    rval = CPXsetdblparam ((*lp)->cplex_env, CPX_PARAM_EPOPT, 1.0E-9);
    if (rval) {
        fprintf (stderr, "CPXsetdblparam CPX_PARAM_EPOPT failed\n");
        goto CLEANUP;
    }
    rval = CPXsetdblparam ((*lp)->cplex_env, CPX_PARAM_EPRHS, 1.0E-9);
    if (rval) {
        fprintf (stderr, "CPXsetdblparam CPX_PARAM_EPRHS failed\n");
        goto CLEANUP;
    }

CLEANUP:

    if (rval) {
        CClp_free (lp);
    }

    return rval;
}

int CClp_force_perturb (CClp *lp)
{
    int rval = 0;

    rval = CPXsetintparam (lp->cplex_env, CPX_PARAM_PERIND, 1);
    if (rval) {
        fprintf (stderr, "CPXsetintparam CPX_PARAM_PERIND failed\n");
    }

    return rval;
}

int CClp_tune_small (CClp *lp)
{
    int rval;

    rval = CPXsetintparam (lp->cplex_env, CPX_PARAM_DPRIIND,
                           CPX_DPRIIND_FULL);
    if (rval) {
        fprintf (stderr, "CPXsetintparam CPX_PARAM_DPRIIND failed\n");
        return rval;
    }
    rval = CPXsetintparam (lp->cplex_env, CPX_PARAM_PPRIIND,
                           CPX_PPRIIND_AUTO);
    if (rval) {
        fprintf (stderr, "CPXsetintparam CPX_PARAM_PPRIIND failed\n");
        return rval;
    }
    return 0;
}

int CClp_disable_presolve (CClp *lp)
{
    if (CPXsetintparam (lp->cplex_env, CPX_PARAM_PREIND, CPX_OFF)) {
        fprintf (stderr, "CPXsetintparam CPX_PARAM_PREIND failed\n");
        return 1;
    }
    if (CPXsetintparam (lp->cplex_env, CPX_PARAM_AGGIND, 0)) {
        fprintf (stderr, "CPXsetintparam CPX_PARAM_AGGIND failed\n");
        return 1;
    }
    return 0;
}

void CClp_free (CClp **lp)
{
    if (*lp) {
        if ((*lp)->cplex_env) {
            if ((*lp)->cplex_lp) {
                CPXfreeprob ((*lp)->cplex_env, &((*lp)->cplex_lp));
            }
            CPXcloseCPLEX (&((*lp)->cplex_env));
        }
        if ((*lp)->cplex40_data) {
            free_cplex40_data (&((*lp)->cplex40_data));
        }
        CC_FREE (*lp, CClp);
    }
}

void CClp_freelp (CClp **lp)
{
    if (*lp) {
        if ((*lp)->cplex_lp) {
            CPXfreeprob ((*lp)->cplex_env, &((*lp)->cplex_lp));
        }
        if ((*lp)->cplex40_data) {
            free_cplex40_data (&((*lp)->cplex40_data));
        }
    }
}

int CClp_loadlp (CClp *lp, const char *name, int ncols, int nrows,
        int objsense, double *obj, double *rhs, char *sense, int *matbeg,
        int *matcnt, int *matind, double *matval, double *lb, double *ub)
{
    int rval = 0;
    CCcplex40_data *data = (CCcplex40_data *) NULL;
    int cspace, rspace, nzspace;

    rval = init_cplex40_data (&lp->cplex40_data, ncols, nrows, name, obj,
                              rhs, sense, matbeg, matcnt, matind, matval,
                              lb, ub, &cspace, &rspace, &nzspace);
    if (rval) {
        fprintf (stderr, "init_cplex40_data failed\n"); return 1;
    }
    data = lp->cplex40_data;
    lp->cplex_lp = CPXloadlp (lp->cplex_env, data->probname,
                              ncols, nrows, objsense,
                              data->obj, data->rhs,
                              data->sense, data->matbeg,
                              data->matcnt, data->matind,
                              data->matval, data->lb,
                              data->ub, (double *) NULL, cspace,
                              rspace, nzspace);
    if (!lp->cplex_lp) {
       fprintf (stderr, "CPXloadlp failed\n");  return 1;
    }

    return 0;
}

int CClp_create (CClp *lp, const char *name)
{
    CCcplex40_data *data = (CCcplex40_data *) NULL;

    int ncols = 1;
    int nrows = 1;
    double obj[1], rhs[1];
    char sense[1];
    int matbeg[1], matcnt[1], matind[1];
    double matval[1];
    double lb[1], ub[1];
    int cspace, rspace, nzspace;

    obj[0] = 0.0;
    rhs[0] = 0.0;
    sense[0] = 'E';
    matbeg[0] = 0;
    matcnt[0] = 0;
    lb[0] = 0.0;
    ub[0] = 0.0;

    if ( init_cplex40_data (&lp->cplex40_data, ncols, nrows, name, obj,
                            rhs, sense, matbeg, matcnt, matind, matval,
                            lb, ub, &cspace, &rspace, &nzspace) ) {
        fprintf (stderr, "init_cplex40_data failed\n"); return 1;
    }
    data = lp->cplex40_data;
    lp->cplex_lp = CPXloadlp (lp->cplex_env, data->probname,
                              ncols, nrows, 1,
                              data->obj, data->rhs,
                              data->sense, data->matbeg,
                              data->matcnt, data->matind,
                              data->matval, data->lb,
                              data->ub, (double *) NULL, cspace,
                              rspace, nzspace);

    if (!lp->cplex_lp) {
        fprintf (stderr, "CPXloadlp failed\n");  return 1;
    }

    /* The CPLEX 4.0 CPXloadlp() doesn't accept problems with
       nrows=0 or ncols=0, so we start with nrows=ncols=1, and
       then delete the bogus row/column the first time that real
       rows/columns are added.  That's the reason for the 'triv'
       indicators defined below.

       It does NOT work to immediately delete the initial row and
       column after loading, even though CPXdelrows() and CPXdelcols()
       would allow that, since there are various function calls
       in Cplex that assume nrows and ncols are both positive */

    data->triv_row_ind = 1;
    data->triv_col_ind = 1;

    return 0;
}

int CClp_new_row (CClp *lp, char sense, double rhs)
{
    double localrhs[1];
    char   localsense[1];
    int    rmatbeg[1];
    int    rmatind[1];
    double rmatval[1];

    localrhs[0] = rhs;
    localsense[0] = sense;
    rmatbeg[0] = 0;

    if ( CClp_addrows (lp, 1, 0, localrhs, localsense,
                       rmatbeg, rmatind, rmatval) ) {
        fprintf (stderr, "CClp_addrows failed\n");  return 1;
    }
    return 0;
}

int CClp_change_sense (CClp *lp, int row, char sense)
{
    int xindex[1];
    char asense[1];
    int rval;

    xindex[0] = row;
    asense[0] = sense;
    rval = CPXchgsense (lp->cplex_env, lp->cplex_lp, 1, xindex, asense);
    if (rval) {
        fprintf (stderr, "CPXchgsense failed\n");
        return rval;
    }
    return 0;
}

int CClp_opt (CClp *lp, int method)
{
    int  rval = 0;

    switch (method) {
        case CClp_METHOD_PRIMAL:
            rval = primalopt (lp);
            break;
        case CClp_METHOD_DUAL:
            rval = dualopt (lp);
            break;
        case CClp_METHOD_BARRIER:
            rval = baropt (lp);
            break;
        default:
            rval = 1;
            fprintf (stderr, "Nonexistent method in CClp_opt\n");
            break;
    }
    return rval;
}

static int primalopt (CClp *lp)
{
    int rval;
    int solstat;
#ifdef  CC_CPLEX_WRITE_PRIMAL
    static int  probcnt = 0;
    char probname[100];

    sprintf (probname, "prim%d.sav", probcnt);
    probcnt++;
    printf ("Writing %s\n", probname);
    CPXsavwrite (lp->cplex_env, lp->cplex_lp, probname);
#endif

    rval = CPXoptimize (lp->cplex_env, lp->cplex_lp);
    if (rval) {
        if (rval == CPXERR_PRESLV_INForUNBD) {
            int old, oldagg;
            printf ("Cplex presolve failed, switch to simplex\n");
            fflush (stdout);
            if (CPXgetintparam (lp->cplex_env, CPX_PARAM_PREIND, &old)) {
                fprintf (stderr, "CPXgetintparam CPX_PARAM_PREIND failed\n");
                return 1;
            }
            if (CPXsetintparam (lp->cplex_env, CPX_PARAM_PREIND, CPX_OFF)) {
                fprintf (stderr, "CPXsetintparam CPX_PARAM_PREIND failed\n");
                return 1;
            }
            if (CPXgetintparam (lp->cplex_env, CPX_PARAM_AGGIND, &oldagg)) {
                fprintf (stderr, "CPXgetintparam CPX_PARAM_AGGIND failed\n");
                return 1;
            }
            if (CPXsetintparam (lp->cplex_env, CPX_PARAM_AGGIND, 0)) {
                fprintf (stderr, "CPXsetintparam CPX_PARAM_AGGIND failed\n");
                return 1;
            }
            rval = CPXoptimize (lp->cplex_env, lp->cplex_lp);
            if (rval) {
                fprintf (stderr, "CPXoptimize failed\n");
                return 1;
            }
            if (CPXsetintparam (lp->cplex_env, CPX_PARAM_PREIND, old)) {
                fprintf (stderr, "CPXsetintparam CPX_PARAM_PREIND failed\n");
                return 1;
            }
            if (CPXsetintparam (lp->cplex_env, CPX_PARAM_AGGIND, oldagg)) {
                fprintf (stderr, "CPXsetintparam CPX_PARAM_AGGIND failed\n");
                return 1;
            }
        } else {
            fprintf (stderr, "CPXoptimize failed\n");
            return 1;
        }
    }
    solstat = CPXgetstat (lp->cplex_env, lp->cplex_lp);
    if (solstat == CPX_INFEASIBLE) {
        return 2;
    } else if (solstat != CPX_OPTIMAL && solstat != CPX_OPTIMAL_INFEAS) {
        fprintf (stderr, "Cplex optimization status %d\n", solstat);
        return 1;
    }
    return 0;
}

static int dualopt (CClp *lp)
{
    int rval;
    int solstat;
#ifdef  CC_CPLEX_WRITE_DUAL
    static int  probcnt = 0;
    char probname[100];

    sprintf (probname, "dual%d.sav", probcnt);
    probcnt++;
    printf ("Writing %s\n", probname);
    CPXsavwrite (lp->cplex_env, lp->cplex_lp, probname);
#endif

    rval = CPXdualopt (lp->cplex_env, lp->cplex_lp);
    if (rval) {
        if (rval == CPXERR_PRESLV_INForUNBD) {
            int old, oldagg;
            printf ("Cplex presolve failed, switch to simplex\n");
            fflush (stdout);
            if (CPXgetintparam (lp->cplex_env, CPX_PARAM_PREIND, &old)) {
                fprintf (stderr, "CPXgetintparam CPX_PARAM_PREIND failed\n");
                return 1;
            }
            if (CPXsetintparam (lp->cplex_env, CPX_PARAM_PREIND, CPX_OFF)) {
                fprintf (stderr, "CPXsetintparam CPX_PARAM_PREIND failed\n");
                return 1;
            }
            if (CPXgetintparam (lp->cplex_env, CPX_PARAM_AGGIND, &oldagg)) {
                fprintf (stderr, "CPXgetintparam CPX_PARAM_AGGIND failed\n");
                return 1;
            }
            if (CPXsetintparam (lp->cplex_env, CPX_PARAM_AGGIND, 0)) {
                fprintf (stderr, "CPXsetintparam CPX_PARAM_AGGIND failed\n");
                return 1;
            }
            rval = CPXdualopt (lp->cplex_env, lp->cplex_lp);
            if (rval) {
                fprintf (stderr, "CPXdualopt failed\n");
                return 1;
            }
            if (CPXsetintparam (lp->cplex_env, CPX_PARAM_PREIND, old)) {
                fprintf (stderr, "CPXsetintparam CPX_PARAM_PREIND failed\n");
                return 1;
            }
            if (CPXsetintparam (lp->cplex_env, CPX_PARAM_AGGIND, oldagg)) {
                fprintf (stderr, "CPXsetintparam CPX_PARAM_AGGIND failed\n");
                return 1;
            }
        } else {
            fprintf (stderr, "CPXdualopt failed\n");
            return 1;
        }
    }
    solstat = CPXgetstat (lp->cplex_env, lp->cplex_lp);
    if (solstat == CPX_UNBOUNDED) {
        return 2;
    } else if (solstat != CPX_OPTIMAL && solstat != CPX_OPTIMAL_INFEAS) {
        fprintf (stderr, "Cplex optimization status %d\n", solstat);
        if (solstat == CPX_IT_LIM_FEAS) {
            int itlim;
            rval = CPXgetintparam (lp->cplex_env, CPX_PARAM_ITLIM, &itlim);
            if (!rval) {
                printf ("cplex iteration limit: %d\n", itlim);
                fflush (stdout);
            }
        }
        return 1;
    }
    return 0;
}

static int baropt (CClp *lp)
{
    int rval;
    int solstat;
#ifdef CC_CPLEX_WRITE_BARRIER
    static int  probcnt = 0;
    char probname[100];

    sprintf (probname, "barrier%d.sav", probcnt);
    probcnt++;
    printf ("Writing %s\n", probname);
    CPXsavwrite (lp->cplex_env, lp->cplex_lp, probname);
#endif

    rval = CPXbaropt (lp->cplex_env, lp->cplex_lp);
    if (rval) {
        printf ("CPXbaropt failed, calling CPXdualopt\n");
        return dualopt (lp);
    }
    solstat = CPXgetstat (lp->cplex_env, lp->cplex_lp);
    if (solstat != CPX_OPTIMAL) {
        printf ("CPXbaropt returned non-optimal solution, calling CPXdualopt\n");
        return dualopt (lp);
    }
    return 0;
}

#define CC_MAX_REFACTORFREQ  150

int CClp_limited_dualopt (CClp *lp, int iterationlim, int *status,
                          double *objupperlim)
{
    int rval = 0;
    int sval = 0;
    int solstat;

    int got_iterationlim = 0;
    int got_presolveind  = 0;
    int got_aggregateind = 0;
    int got_objupperlim  = 0;
    int got_refactorfreq = 0;
    int got_perind       = 0;

    int    old_iterationlim;
    int    old_presolveind;
    int    old_aggregateind;
    double old_objupperlim;
    int    old_refactorfreq;
    int    old_perind;

    /* REB, 14 October 1997:  If perturbation has been turned on elsewhere,
       we need to turn it off for branch selection */

    rval = CPXgetintparam (lp->cplex_env, CPX_PARAM_PERIND, &old_perind);
    if (rval) {
         fprintf (stderr, "CPXgetintparam CPX_PARAM_PERIND failed\n");
         goto CLEANUP;
    }
    got_perind = 1;

    rval = CPXgetintparam (lp->cplex_env, CPX_PARAM_ITLIM, &old_iterationlim);
    if (rval) {
        fprintf (stderr, "CPXgetintparam CPX_PARAM_ITLIM failed\n");
        goto CLEANUP;
    }
    got_iterationlim = 1;

    rval = CPXgetdblparam (lp->cplex_env, CPX_PARAM_OBJULIM, &old_objupperlim);
    if (rval) {
        fprintf (stderr, "CPXgetdblparam CPX_PARAM_OBJULIM failed\n");
        goto CLEANUP;
    }
    got_objupperlim = 1;

    rval = CPXgetintparam (lp->cplex_env, CPX_PARAM_PREIND, &old_presolveind);
    if (rval) {
        fprintf (stderr, "CPXgetintparam CPX_PARAM_PREIND failed\n");
        goto CLEANUP;
    }
    got_presolveind = 1;

    rval = CPXgetintparam (lp->cplex_env, CPX_PARAM_AGGIND, &old_aggregateind);
    if (rval) {
        fprintf (stderr, "CPXgetintparam CPX_PARAM_AGGIND failed\n");
        goto CLEANUP;
    }
    got_aggregateind = 1;

    rval = CPXgetintparam (lp->cplex_env, CPX_PARAM_REINV, &old_refactorfreq);
    if (rval) {
        fprintf (stderr, "CPXgetintparam CPX_PARAM_REINV failed\n");
        goto CLEANUP;
    }
    got_refactorfreq = 1;

    rval = CPXsetintparam (lp->cplex_env, CPX_PARAM_PERIND, 0);
    if (rval) {
         fprintf (stderr, "CPXsetintparam CPX_PARAM_PERIND failed\n");
         goto CLEANUP;
    }

    rval = CPXsetintparam (lp->cplex_env, CPX_PARAM_ITLIM, iterationlim);
    if (rval) {
        fprintf (stderr, "CPXsetintparam CPX_PARAM_ITLIM failed\n");
        goto CLEANUP;
    }

    if ( iterationlim < CC_MAX_REFACTORFREQ ) {
        rval = CPXsetintparam (lp->cplex_env, CPX_PARAM_REINV, iterationlim+1);
        if (rval) {
            fprintf (stderr, "CPXsetintparam CPX_PARAM_REINV failed\n");
            goto CLEANUP;
        }
    }

    if (objupperlim) {
        rval = CPXsetdblparam (lp->cplex_env, CPX_PARAM_OBJULIM, *objupperlim);
        if (rval) {
            fprintf (stderr, "CPXsetdblparam CPX_PARAM_OBJULIM failed\n");
            goto CLEANUP;
        }
    }

    rval = CPXdualopt (lp->cplex_env, lp->cplex_lp);
    if (rval) {
        if (rval == CPXERR_PRESLV_INForUNBD) {
            printf ("Cplex presolve failed, force simplex\n");
            fflush (stdout);

            if (CPXsetintparam (lp->cplex_env, CPX_PARAM_PREIND, CPX_OFF)) {
                fprintf (stderr, "CPXsetintparam CPX_PARAM_PREIND failed\n");
                goto CLEANUP;
            }
            if (CPXsetintparam (lp->cplex_env, CPX_PARAM_AGGIND, 0)) {
                fprintf (stderr, "CPXsetintparam CPX_PARAM_AGGIND failed\n");
                goto CLEANUP;
            }
            rval = CPXdualopt (lp->cplex_env, lp->cplex_lp);
            if (rval) {
                fprintf (stderr, "CPXdualopt failed\n"); goto CLEANUP;
            }
        } else {
            fprintf (stderr, "CPXdualopt failed\n"); goto CLEANUP;
        }
    }
    solstat = CPXgetstat (lp->cplex_env, lp->cplex_lp);
    if (solstat==CPX_IT_LIM_INFEAS) {
        rval = CPXsetdblparam (lp->cplex_env, CPX_PARAM_OBJULIM, -1.0E75);
        if (rval) {
            fprintf (stderr, "CPXsetdblparam CPX_PARAM_OBJULIM failed\n");
            goto CLEANUP;
        }
        /* We could be even more aggressive here and make the iteration
           limit infinite, but that approach seems contrary to the
           intent of this function.  Hence, the repeat test for
           CPX_IT_LIM_INFEAS below -- REB, 1 July 97 */
        rval = CPXdualopt (lp->cplex_env, lp->cplex_lp);
        if (rval) {
            fprintf (stderr, "CPXdualopt failed\n"); goto CLEANUP;
        }
    }

    if (solstat == CPX_UNBOUNDED) {
        printf ("Infeasible in CPXdualopt\n"); fflush (stdout);
        if (status) *status = CClp_INFEASIBLE;
    } else if (solstat == CPX_IT_LIM_INFEAS) {
        printf ("LP infeasible after the limited number of iterations\n");
        fflush (stdout);
        if (status) *status = CClp_UNKNOWN;
    } else if (solstat != CPX_OPTIMAL && solstat != CPX_OPTIMAL_INFEAS &&
               solstat != CPX_IT_LIM_FEAS && solstat != CPX_OBJ_LIM) {
        fprintf (stderr, "Cplex optimization status %d\n", solstat);
        if (status) *status = CClp_FAILURE;
    } else {
        if (status) *status = CClp_SUCCESS;
    }

CLEANUP:

    if (got_iterationlim == 1) {
        sval = CPXsetintparam (lp->cplex_env, CPX_PARAM_ITLIM,
                               old_iterationlim);
        if (sval) {
            fprintf (stderr, "CPXsetintparam CPX_PARAM_ITLIM failed\n");
            rval = 1;
        }
    }

    if (got_objupperlim == 1) {
        sval = CPXsetdblparam (lp->cplex_env, CPX_PARAM_OBJULIM,
                               old_objupperlim);
        if (sval) {
            fprintf (stderr, "CPXsetdblparam CPX_PARAM_OBJULIM failed\n");
            rval = 1;
        }
    }

    if (got_presolveind == 1) {
        sval = CPXsetintparam (lp->cplex_env, CPX_PARAM_PREIND,
                               old_presolveind);
        if (sval) {
            fprintf (stderr, "CPXsetintparam CPX_PARAM_PREIND failed\n");
            rval = 1;
        }
    }

    if (got_aggregateind == 1) {
        sval = CPXsetintparam (lp->cplex_env, CPX_PARAM_AGGIND,
                               old_aggregateind);
        if (sval) {
            fprintf (stderr, "CPXsetintparam CPX_PARAM_AGGIND failed\n");
            rval = 1;
        }
    }

    if (got_refactorfreq == 1) {
        sval = CPXsetintparam (lp->cplex_env, CPX_PARAM_REINV,
                               old_refactorfreq);
        if (sval) {
            fprintf (stderr, "CPXsetintparam CPX_PARAM_REINV failed\n");
            rval = 1;
        }
    }

    if (got_perind == 1) {
        sval = CPXsetintparam (lp->cplex_env, CPX_PARAM_PERIND, old_perind);
        if (sval) {
            fprintf (stderr, "CPXsetintparam CPX_PARAM_PERIND failed\n");
            rval = 1;;
        }
    }

    if (rval && status) {
        *status = CClp_FAILURE;
    }

    return rval;
}

int CClp_addrows (CClp *lp, int newrows, int newnz, double *rhs, char *sense,
                  int *rmatbeg, int *rmatind, double *rmatval)
{
    int rval = 0;

    rval = CPXfaddrows (lp->cplex_env, lp->cplex_lp, 0, newrows, newnz,
                        rhs, sense, rmatbeg, rmatind, rmatval,
                        (char **) NULL, (char **) NULL);
    if (rval) {
        rval = reallocprob (lp, 0, newrows, newnz);
        if (rval) {
            fprintf (stderr, "reallocprob failed\n");
            return rval;
        }
        rval = CPXfaddrows (lp->cplex_env, lp->cplex_lp, 0, newrows, newnz,
                            rhs, sense, rmatbeg, rmatind, rmatval,
                            (char **) NULL, (char **) NULL);
        if (rval) {
            fprintf (stderr, "CPXfaddrows failed\n");
            return rval;
        }
    }
    if (lp->cplex40_data &&
        lp->cplex40_data->triv_row_ind &&
        CPXgetnumrows (lp->cplex_env, lp->cplex_lp) > 1) {
        rval = CPXdelrows (lp->cplex_env, lp->cplex_lp, 0, 0);
        if (rval) {
           fprintf (stderr, "CPXdelrows failed\n");  return rval;
        }
        lp->cplex40_data->triv_row_ind = 0;
    }
    return rval;
}

int CClp_addcols (CClp *lp, int newcols, int newnz, double *obj,
                  int *cmatbeg, int *cmatind, double *cmatval,
                  double *lb, double *ub)
{
    int rval = 0;

    rval = CPXaddcols (lp->cplex_env, lp->cplex_lp, newcols, newnz, obj,
                   cmatbeg, cmatind, cmatval, lb, ub, (char **) NULL);
    if (rval) {
        rval = reallocprob (lp, newcols, 0, newnz);
        if (rval) {
            fprintf (stderr, "reallocprob failed\n");
            return rval;
        }
        rval = CPXaddcols (lp->cplex_env, lp->cplex_lp, newcols, newnz, obj,
                           cmatbeg, cmatind, cmatval, lb, ub, (char **) NULL);
        if (rval) {
            fprintf (stderr, "CPXaddcols failed\n");
            return rval;
        }
    }
    if (lp->cplex40_data &&
        lp->cplex40_data->triv_col_ind &&
        CPXgetnumcols (lp->cplex_env, lp->cplex_lp) > 1 ) {
        rval = CPXdelcols (lp->cplex_env, lp->cplex_lp, 0, 0);
        if (rval) {
           fprintf (stderr, "CPXdelcols failed\n");  return rval;
        }
        lp->cplex40_data->triv_col_ind = 0;
    }
    return rval;
}

int CClp_delete_row (CClp *lp, int i)
{
    int rval = 0;
    int locali[1];

    locali[0] = i;
    if (CPXpivotin (lp->cplex_env, lp->cplex_lp, locali, 1)) {
        fprintf (stderr, "CPXpivotin failed, continuing anyway\n");
    }
    rval = CPXdelrows (lp->cplex_env, lp->cplex_lp, i, i);
    if (rval) fprintf (stderr, "CPXdelrows failed\n");
    return rval;
}

int CClp_delete_set_of_rows (CClp *lp, int *delstat)
{
    int rval = 0;
    int *dellist = (int *) NULL;
    int delcnt = 0;
    int i;
    int j;
    int rcnt = CPXgetnumrows (lp->cplex_env, lp->cplex_lp);

    for (i=0; i<rcnt; i++) {
        if (delstat[i]) delcnt++;
    }
    if (delcnt == 0) {
        fprintf (stderr, "delete_set_of_rows with no deleted rows\n");
        return 0;
    }
    dellist = CC_SAFE_MALLOC (delcnt, int);
    if (dellist == (int *) NULL) {
        fprintf (stderr, "Out of memory in delete_set_of_rows\n");
        return 1;
    }
    for (i=0, j=0; i<rcnt; i++) {
        if (delstat[i]) {
            dellist[j++] = i;
        }
    }
    if (j != delcnt) {
        fprintf (stderr, "Lost some deleted rows\n");
        CC_FREE (dellist, int);
        return 1;
    }

    if (CPXpivotin (lp->cplex_env, lp->cplex_lp, dellist, delcnt)) {
        fprintf (stderr, "CPXpivotin failed, continuing anyway\n");
    }
    CC_FREE (dellist, int);
    
    rval = CPXdelsetrows (lp->cplex_env, lp->cplex_lp, delstat);
    if (rval) fprintf (stderr, "CPXdelsetrows failed\n");
    return rval;
}

int CClp_delete_column (CClp *lp, int i)
{
    int rval = 0;

    rval = CPXdelcols (lp->cplex_env, lp->cplex_lp, i, i);
    if (rval) fprintf (stderr, "CPXdelcols failed\n");
    return rval;
}

int CClp_delete_set_of_columns (CClp *lp, int *delstat)
{
    int rval = 0;

    rval = CPXdelsetcols (lp->cplex_env, lp->cplex_lp, delstat);
    if (rval) fprintf (stderr, "CPXdelsetcols failed\n");
    return rval;
}

int CClp_setbnd (CClp *lp, int col, char lower_or_upper, double bnd)
{
    int cindex[1];
    double bd[1];
    char lu[1];
    int rval;

    cindex[0] = col;
    lu[0] = lower_or_upper;
    bd[0] = bnd;

    rval = CPXchgbds (lp->cplex_env, lp->cplex_lp, 1, cindex, lu, bd);
    if (rval) {
        fprintf (stderr, "Couldn't set bnd on variable %d in cplex\n", col);
        return rval;
    }
    return 0;
}

int CClp_get_warmstart (CClp *lp, CClp_warmstart **w)
{
    int rval = 0;

    CClp_free_warmstart (w);

    (*w) = CC_SAFE_MALLOC (1, CClp_warmstart);
    if ((*w) == (CClp_warmstart *) NULL) {
        fprintf (stderr, "Out of memory in CClp_get_warmstart\n");
        rval = 1; goto CLEANUP;
    }

    (*w)->ccount = 0;
    (*w)->rcount = 0;
    (*w)->cstat  = (int *) NULL;
    (*w)->rstat  = (int *) NULL;
    (*w)->dnorm  = (double *) NULL;
    (*w)->bhead  = (int *) NULL;
    (*w)->dlen   = 0;

    (*w)->ccount = CPXgetnumcols (lp->cplex_env, lp->cplex_lp);
    if ((*w)->ccount == 0) {
        fprintf (stderr, "No columns in LP\n");
        rval = 1; goto CLEANUP;
    }
    (*w)->rcount = CPXgetnumrows (lp->cplex_env, lp->cplex_lp);
    if ((*w)->rcount == 0) {
        fprintf (stderr ,"No rows in LP\n");
        rval = 1; goto CLEANUP;
    }

    (*w)->cstat = CC_SAFE_MALLOC ((*w)->ccount, int);
    (*w)->rstat = CC_SAFE_MALLOC ((*w)->rcount, int);
    (*w)->dnorm = CC_SAFE_MALLOC ((*w)->rcount, double);
    (*w)->bhead = CC_SAFE_MALLOC ((*w)->rcount, int);
    if (!(*w)->cstat || !(*w)->rstat || !(*w)->dnorm || !(*w)->bhead) {
        fprintf (stderr, "out of memory in CClp_get_warmstart\n");
        rval = 1; goto CLEANUP;
    }

    rval = CPXgetbase (lp->cplex_env, lp->cplex_lp, (*w)->cstat, (*w)->rstat);
    if (rval) {
        fprintf (stderr, "CPXgetbase failed\n"); goto CLEANUP;
    }
    rval = CPXgetdnorms (lp->cplex_env, lp->cplex_lp, (*w)->dnorm, (*w)->bhead,
                         &(*w)->dlen);
    if (rval) {
        fprintf (stderr, "CPXgetdnorms failed\n");
        CC_IFFREE ((*w)->dnorm, double);
        CC_IFFREE ((*w)->bhead, int);
        (*w)->dlen = 0;
    }

    return 0;

CLEANUP:

    CClp_free_warmstart (w);
    return rval;
}

int CClp_load_warmstart (CClp *lp, CClp_warmstart *w)
{
    int rval = 0;

    if (w->cstat && w->rstat) {
        rval = CPXloadbase (lp->cplex_env, lp->cplex_lp, w->cstat, w->rstat);
        if (rval) {
            fprintf (stderr, "CPXloadbase failed\n"); goto CLEANUP;
        }
        if (w->dnorm && w->bhead) {
            rval = CPXloaddnorms (lp->cplex_env, lp->cplex_lp, w->dnorm,
                                  w->bhead, w->dlen);
            if (rval) {
                fprintf (stderr, "CPXloaddnorms failed\n"); goto CLEANUP;
            }
        }
    } else {
        printf ("WARNING: No basis in call to load_warmstart\n");
        fflush (stdout);
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

    (*w)->ccount = 0;
    (*w)->rcount = 0;
    (*w)->cstat = (int *) NULL;
    (*w)->rstat = (int *) NULL;
    (*w)->dnorm = (double *) NULL;
    (*w)->bhead = (int *) NULL;

    (*w)->ccount = i->ccount;
    if ((*w)->ccount == 0) {
        fprintf (stderr, "No columns in CClp_info\n");
        rval = 1; goto CLEANUP;
    }
    (*w)->rcount = i->rcount;
    if ((*w)->rcount == 0) {
        fprintf (stderr, "No rows in CClp_info\n");
        rval = 1; goto CLEANUP;
    }

    (*w)->cstat = CC_SAFE_MALLOC ((*w)->ccount, int);
    (*w)->rstat = CC_SAFE_MALLOC ((*w)->rcount, int);
    if (!(*w)->cstat || !(*w)->rstat) {
        fprintf (stderr, "out of memory in CClp_get_warmstart\n");
        rval = 1; goto CLEANUP;
    }

    for (j=0; j<(*w)->ccount; j++) {
        (*w)->cstat[j] = i->cstat[j];
    }
    for (j=0; j<(*w)->rcount; j++) {
        (*w)->rstat[j] = i->rstat[j];
    }

    return 0;

  CLEANUP:
    CClp_free_warmstart (w);
    return rval;
}

void CClp_free_warmstart (CClp_warmstart **w)
{
    if ((*w) != (CClp_warmstart *) NULL) {
        CC_IFFREE ((*w)->cstat, int);
        CC_IFFREE ((*w)->rstat, int);
        CC_IFFREE ((*w)->dnorm, double);
        CC_IFFREE ((*w)->bhead, int);
        CC_FREE (*w, CClp_warmstart);
    }
}

int CClp_sread_warmstart (CC_SFILE *f, CClp_warmstart **w)
{
    char name[5];
    int i;
    int ccount;
    int rcount;
    int has_dnorms;

    CClp_free_warmstart (w);

    for (i=0; i<4; i++) {
        if (CCutil_sread_char (f, &name[i])) goto CLEANUP;
    }
    name[4] = '\0';

    if (strncmp (name, SOLVER_WARMSTART_NAME, 4)) {
        fprintf (stderr, "warmstart for another solver (%s) ignored\n", name);
        return 0;
    }

    if (CCutil_sread_int (f, &ccount)) goto CLEANUP;
    if (CCutil_sread_int (f, &rcount)) goto CLEANUP;

    (*w) = CC_SAFE_MALLOC (1, CClp_warmstart);
    if ((*w) == (CClp_warmstart *) NULL) {
        fprintf (stderr, "Out of memory in CClp_sread_warmstart\n");
        goto CLEANUP;
    }

    (*w)->ccount = 0;
    (*w)->rcount = 0;
    (*w)->cstat  = (int *) NULL;
    (*w)->rstat  = (int *) NULL;
    (*w)->dnorm  = (double *) NULL;
    (*w)->bhead  = (int *) NULL;

    (*w)->cstat = CC_SAFE_MALLOC (ccount, int);
    (*w)->rstat = CC_SAFE_MALLOC (rcount, int);
    if ((*w)->cstat == (int *) NULL ||
        (*w)->rstat == (int *) NULL) {
        fprintf (stderr, "out of memory in CClp_sread_warmstart\n");
        goto CLEANUP;
    }
    for (i = 0; i < ccount; i++) {
        if (CCutil_sread_bits (f, &(((*w)->cstat)[i]), 2)) goto CLEANUP;
    }
    for (i = 0; i < rcount; i++) {
        if (CCutil_sread_bits (f, &(((*w)->rstat)[i]), 1)) goto CLEANUP;
    }

    if (CCutil_sread_int (f, &has_dnorms)) goto CLEANUP;

    if (has_dnorms) {
        (*w)->dnorm = CC_SAFE_MALLOC (rcount, double);
        (*w)->bhead = CC_SAFE_MALLOC (rcount, int);
        if ((*w)->dnorm == (double *) NULL ||
            (*w)->bhead == (int *) NULL) {
            fprintf (stderr, "out of memory in CClp_sread_warmstart\n");
            goto CLEANUP;
        }
        for (i = 0; i < rcount; i++) {
            if (CCutil_sread_double (f, &(((*w)->dnorm)[i]))) goto CLEANUP;
        }
        for (i = 0; i < rcount; i++) {
            if (CCutil_sread_int (f, &(((*w)->bhead)[i]))) goto CLEANUP;
        }
        if (CCutil_sread_int (f, &((*w)->dlen))) goto CLEANUP;
    }

    (*w)->ccount = ccount;
    (*w)->rcount = rcount;
    
    return 0;

CLEANUP:

    CClp_free_warmstart (w);
    return 1;
}

int CClp_swrite_warmstart (CC_SFILE *f, CClp_warmstart *w)
{
    int i;
    const char *name = SOLVER_WARMSTART_NAME;

    for (i=0; i<4; i++) {
        if (CCutil_swrite_char (f, name[i])) return 1;
    }

    if (CCutil_swrite_int (f, w->ccount)) return 1;
    if (CCutil_swrite_int (f, w->rcount)) return 1;

    for (i = 0; i < w->ccount; i++) {
        if (CCutil_swrite_bits (f, w->cstat[i], 2)) return 1;
    }

    for (i = 0; i < w->rcount; i++) {
        if (CCutil_swrite_bits (f, w->rstat[i], 1)) return 1;
    }

    if (w->dnorm == (double *) NULL) {
        if (CCutil_swrite_int (f, 0)) return 1;
    } else {
        if (CCutil_swrite_int (f, 1)) return 1;
        for (i = 0; i < w->rcount; i++) {
            if (CCutil_swrite_double (f, w->dnorm[i])) return 1;
        }
        for (i = 0; i < w->rcount; i++) {
            if (CCutil_swrite_int (f, w->bhead[i])) return 1;
        }
        if (CCutil_swrite_int (f, w->dlen)) return 1;
    }

    return 0;
}

int CClp_get_info (CClp *lp, CClp_info **i)
{
    int rval = 0;

    CClp_free_info (i);

    (*i) = CC_SAFE_MALLOC (1, CClp_info);
    if ((*i) == (CClp_info *) NULL) {
        fprintf (stderr, "Out of memory in CClp_get_info\n");
        rval = 1; goto CLEANUP;
    }

    (*i)->ccount = 0;
    (*i)->rcount = 0;
    (*i)->cstat = (int *) NULL;
    (*i)->rstat = (int *) NULL;

    (*i)->ccount = CPXgetnumcols (lp->cplex_env, lp->cplex_lp);
    if ((*i)->ccount == 0) {
        fprintf (stderr, "No columns in CClp_get_info\n");
        rval = 1; goto CLEANUP;
    }
    (*i)->rcount = CPXgetnumrows (lp->cplex_env, lp->cplex_lp);
    if ((*i)->rcount == 0) {
        fprintf (stderr ,"No rows in CClp_get_info\n");
        rval = 1; goto CLEANUP;
    }

    (*i)->cstat = CC_SAFE_MALLOC ((*i)->ccount, int);
    (*i)->rstat = CC_SAFE_MALLOC ((*i)->rcount, int);
    if (!(*i)->cstat || !(*i)->rstat) {
        fprintf (stderr, "out of memory in CClp_get_info\n");
        rval = 1; goto CLEANUP;
    }

    rval = CPXgetbase (lp->cplex_env, lp->cplex_lp, (*i)->cstat, (*i)->rstat);
    if (rval) {
        fprintf (stderr, "CPXgetbase failed\n"); goto CLEANUP;
    }

    return 0;

CLEANUP:

    CClp_free_info (i);
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

    (*i)->ccount = 0;
    (*i)->rcount = 0;
    (*i)->cstat = (int *) NULL;
    (*i)->rstat = (int *) NULL;

    (*i)->ccount = ccount;
    if (ccount == 0) {
        fprintf (stderr, "No columns in CClp_create_info\n");
        rval = 1; goto CLEANUP;
    }
    (*i)->rcount = rcount;
    if (rcount == 0) {
        fprintf (stderr, "No rows in CClp_create_info\n");
        rval = 1; goto CLEANUP;
    }

    (*i)->cstat = CC_SAFE_MALLOC ((*i)->ccount, int);
    (*i)->rstat = CC_SAFE_MALLOC ((*i)->rcount, int);
    if (!(*i)->cstat || !(*i)->rstat) {
        fprintf (stderr, "out of memory in CClp_create_info\n");
        rval = 1; goto CLEANUP;
    }

    for (j=0; j<ccount; j++) {
        (*i)->cstat[j] = 0;
    }
    for (j=0; j<rcount; j++) {
        (*i)->rstat[j] = 0;
    }

    return 0;

  CLEANUP:
    CClp_free_info (i);
    return rval;
}

int CClp_is_col_active (CClp_info *i, int c)
{
    if (c < 0 || c >= i->ccount) return 0;
    return i->cstat[c] == 1 || i->cstat[c] == 2;
}

int CClp_is_row_active (CClp_info *i, int r)
{
    if (r < 0 || r >= i->rcount) return 0;
    return i->rstat[r] == 0;
}

void CClp_set_col_active (CClp_info *i, int c)
{
    if (c >= 0 && c < i->ccount) i->cstat[c] = 1;
}

void CClp_set_col_inactive (CClp_info *i, int c)
{
    if (c >= 0 && c < i->ccount) i->cstat[c] = 0;
}

void CClp_set_col_upper (CClp_info *i, int c)
{
    if (c >= 0 && c < i->ccount) i->cstat[c] = 2;
}

void CClp_set_row_active (CClp_info *i, int r)
{
    if (r >= 0 && r < i->rcount) i->rstat[r] = 0;
}

void CClp_set_row_inactive (CClp_info *i, int r)
{
    if (r >= 0 && r < i->rcount) i->rstat[r] = 1;
}

void CClp_free_info (CClp_info **i)
{
    if ((*i) != (CClp_info *) NULL) {
        CC_IFFREE ((*i)->cstat, int);
        CC_IFFREE ((*i)->rstat, int);
        CC_FREE (*i, CClp_info);
    }
}

int CClp_x (CClp *lp, double *x)
{
    int rval = 0;
    int ncols;

    ncols = CPXgetnumcols (lp->cplex_env, lp->cplex_lp);
    if (ncols == 0) {
        fprintf (stderr, "No columns in LP\n");
        return 1;
    }
    rval = CPXgetx (lp->cplex_env, lp->cplex_lp, x, 0, ncols - 1);
    if (rval) {
        fprintf (stderr, "CPXgetx failed\n");
        return rval;
    }
    return 0;
}

int CClp_rc (CClp *lp, double *rc)
{
    int rval = 0;
    int ncols;

    ncols = CPXgetnumcols (lp->cplex_env, lp->cplex_lp);
    if (ncols == 0) {
        fprintf (stderr, "No columns in LP\n"); return 1;
    }
    rval = CPXgetdj (lp->cplex_env, lp->cplex_lp, rc, 0, ncols - 1);
    if (rval) {
        fprintf (stderr, "CPXgetdj failed\n"); return rval;
    }
    return 0;
}

int CClp_pi (CClp *lp, double *pi)
{
    int rval = 0;
    int nrows;

    if ( CPXgetmethod (lp->cplex_env, lp->cplex_lp) == CPXALG_DUAL  &&
         CPXgetstat (lp->cplex_env, lp->cplex_lp )  == CPX_UNBOUNDED  ) {
        rval = getfarkasmultipliers (lp, pi);
        if (rval) {
            fprintf (stderr, "getfarkasmultipliers failed\n"); return rval;
        }
        return 0;
    }

    nrows = CPXgetnumrows (lp->cplex_env, lp->cplex_lp);
    if (nrows == 0) {
        fprintf (stderr, "No rows in LP\n"); return 1;
    }
    rval = CPXgetpi (lp->cplex_env, lp->cplex_lp, pi, 0, nrows - 1);
    if (rval) {
        fprintf (stderr, "CPXgetpi failed\n"); return rval;
    }
    return 0;
}

int CClp_objval (CClp *lp, double *obj)
{
    int rval;

    rval = CPXgetobjval (lp->cplex_env, lp->cplex_lp, obj);
    if (rval) {
        fprintf (stderr, "CPXgetobjval failed\n");
        return rval;
    }
    return 0;
}

int CClp_nrows (CClp *lp)
{
    return CPXgetnumrows (lp->cplex_env, lp->cplex_lp);
}

int CClp_ncols (CClp *lp)
{
    return CPXgetnumcols (lp->cplex_env, lp->cplex_lp);
}

int CClp_nnonzeros (CClp *lp)
{
    return CPXgetnumnz (lp->cplex_env, lp->cplex_lp);
}

int CClp_status (CClp *lp, int *status)
{
    int solmethod, solstat;

    solmethod = CPXgetmethod (lp->cplex_env, lp->cplex_lp);
    if (solmethod == CPXALG_PRIMAL || solmethod == CPXALG_DUAL) {
        solstat = CPXgetstat (lp->cplex_env, lp->cplex_lp);
        if (solstat == CPX_OPTIMAL || solstat == CPX_OPTIMAL_INFEAS) {
            *status = 0;
            return 0;
        } else if (solstat == CPX_UNBOUNDED && solmethod == CPXALG_DUAL) {
            *status = 1;
            return 0;
        } else {
            fprintf (stderr, "lp in an unknown state: %d %d\n",
                           solmethod, solstat);
            *status = -1;
            return 1;
        }
    } else {
        fprintf (stderr, "lp not solved by usual methods: %d\n", solmethod);
        *status = -2;
        return 1;
    }
}

int CClp_getweight (CClp *lp, int nrows, int *rmatbeg, int *rmatind,
                    double *rmatval, double *weight)
{
    int rval = 0;

    rval = CPXgetweight (lp->cplex_env, lp->cplex_lp, nrows,
                         rmatbeg, rmatind, rmatval, weight, CPX_DPRIIND_STEEP);
    if (rval) {
        fprintf (stderr, "CPXgetweight failed\n");
    }
    return rval;
}

int CClp_dump_lp (CClp *lp, const char *fname)
{
    int rval = 0;
    char nambuf[32];

    /* We copy the name since CPXsavwrite doesn't declare fname as const */
    strncpy (nambuf, fname, sizeof (nambuf));
    nambuf[sizeof(nambuf)-1] = '\0';

    rval = CPXsavwrite (lp->cplex_env, lp->cplex_lp, nambuf);
    if (rval) {
        fprintf (stderr, "CPXsavwrite failed\n");
    }
    return rval;
}

#define OURCPLEXZERO    (1.0E-10)
#define OURCPLEX_INTTOL (0.0001)

int CClp_getgoodlist (CClp *lp, int *goodlist, int *goodlen_p,
                      double *downpen, double *uppen)
{
    if ( CPXmdualopt (lp->cplex_env, lp->cplex_lp, goodlist, goodlen_p,
                      downpen, uppen) ) {
        fprintf (stderr, "CPXmdualopt failed\n");
        return 1;
    } else {
        return 0;
    }
}

int CClp_strongbranch (CClp *lp, int *candidatelist, int ncand,
                       double *downpen, double *uppen, int iterations,
                       double upperbound)
{
    double oldupperbound;
    int rval = 0;
    int sval = 0;
    int old_perind;
    int i;

    rval = CPXgetdblparam (lp->cplex_env, CPX_PARAM_OBJULIM,
                           &oldupperbound);
    if (rval) {
        fprintf (stderr, "CPXgetdblparam failed\n"); return rval;
    }
    rval = CPXsetdblparam (lp->cplex_env, CPX_PARAM_OBJULIM,
                           upperbound);
    if (rval) {
        fprintf (stderr, "CPXsetdblparam failed\n"); return rval;
    }

    /* REB, 14 October 1997:  If perturbation has been turned on elsewhere,
       we need to turn it off for branch selection */

    rval = CPXgetintparam (lp->cplex_env, CPX_PARAM_PERIND, &old_perind);
    if (rval) {
         fprintf (stderr, "CPXgetintparam failed\n"); return rval;
    }

    rval = CPXsetintparam (lp->cplex_env, CPX_PARAM_PERIND, 0);
    if (rval) {
         fprintf (stderr, "CPXsetintparam failed\n"); return rval;
    }

    rval = CPXstrongbranch (lp->cplex_env, lp->cplex_lp, candidatelist,
                            ncand, downpen, uppen, iterations);
    if (rval) {
        fprintf (stderr, "CPXstrongbranch failed\n");
        sval = CPXsetdblparam (lp->cplex_env, CPX_PARAM_OBJULIM,
                               oldupperbound);
        if (sval) {
            fprintf (stderr,
                     "CPXsetdblparam failed with return code %d\n",
                     sval);
        }
        return rval;
    }

    rval = CPXsetintparam (lp->cplex_env, CPX_PARAM_PERIND, old_perind);
    if (rval) {
         fprintf (stderr, "CPXsetintparam failed\n"); return rval;
    }

    rval = CPXsetdblparam (lp->cplex_env, CPX_PARAM_OBJULIM,
                           oldupperbound);
    if (rval) {
        fprintf (stderr, "CPXsetdblparam failed\n"); return rval;
    }

    for (i=0; i<ncand; i++) {
        if (downpen[i] > upperbound) downpen[i] = upperbound;
        if (uppen[i] > upperbound) uppen[i] = upperbound;
    }

    return 0;
}

static int getfarkasmultipliers (CClp *lp, double *y)
{
    int  rval = 0;

    int  i = 0, nrows, idiv, jdiv;
    double  val, lb, ub;
    int  *bhead =  (int *) NULL;
    char  *sense = (char *) NULL;

    if ( lp->cplex_env == (struct cpxenv *) NULL ||
         lp->cplex_lp  == (struct cpxlp *)  NULL) {
        rval = 1;  fprintf (stderr, "env object or lp object is NULL\n");
        goto CLEANUP;
    }

    if ( CPXgetmethod (lp->cplex_env, lp->cplex_lp) != CPXALG_DUAL  ||
         CPXgetstat (lp->cplex_env, lp->cplex_lp )  != CPX_UNBOUNDED  ) {
        rval = 1;  fprintf (stderr, "Incorrect solution type\n");
        goto CLEANUP;
    }

    if ( CPXgetijdiv (lp->cplex_env, lp->cplex_lp, &idiv, &jdiv) ) {
        rval = 1;  fprintf (stderr, "CPXgetijdiv failed\n");
        goto CLEANUP;
    }

    if ( (jdiv == -1  &&  idiv == -1) ||
         (jdiv != -1  &&  idiv != -1)   ) {
        rval = 1;  fprintf (stderr, "CPLEX returned illegal indices\n");
        goto CLEANUP;
    }

    nrows = CPXgetnumrows (lp->cplex_env, lp->cplex_lp);
    if ( nrows == 0 ) {
        rval = 1;  fprintf (stderr, "lp->cplex_lp has no rows\n");
        goto CLEANUP;
    }

    bhead = CC_SAFE_MALLOC (nrows, int);
    sense = CC_SAFE_MALLOC (nrows, char);
    if ( bhead == (int *) NULL ||
         sense == (char *) NULL   ) {
        rval = -1;  fprintf (stderr, "Out of memory\n");
        goto CLEANUP;
    }

    if ( CPXgetbhead (lp->cplex_env, lp->cplex_lp, bhead, NULL) ) {
        rval = 1;  fprintf (stderr, "CPXgetbhead failed\n");
        goto CLEANUP;
    }

    if ( CPXgetsense (lp->cplex_env, lp->cplex_lp, sense, 0, nrows-1) ) {
        rval = 1;  fprintf (stderr, "CPXgetsense failed\n");
        goto CLEANUP;
    }

    if ( jdiv >= 0 ) {
        for (i = 0; i < nrows; i++) {
            if ( bhead[i] == jdiv )  break;
        }
        if ( i == nrows ) {
            rval = 1;  fprintf (stderr, "Basis index not found\n");
            goto CLEANUP;
        }
        if ( CPXgetx (lp->cplex_env, lp->cplex_lp, &val, jdiv, jdiv) ) {
            rval = 1;  fprintf (stderr, "CPXgetx failed\n");
            goto CLEANUP;
        }
        if ( CPXgetlb (lp->cplex_env, lp->cplex_lp, &lb, jdiv, jdiv) ) {
            rval = 1;  fprintf (stderr, "CPXgetlb failed\n");
            goto CLEANUP;
        }
        if ( CPXgetub (lp->cplex_env, lp->cplex_lp, &ub, jdiv, jdiv) ) {
            rval = 1;  fprintf (stderr, "CPXgetub failed\n");
            goto CLEANUP;
        }
    } else {
        for (i = 0; i < nrows; i++) {
            if ( bhead[i] ==  -idiv-1 )  break;
        }
        if ( i == nrows ) {
            rval = 1;  fprintf (stderr, "Basis index not found\n");
            goto CLEANUP;
        }
        if ( CPXgetslack (lp->cplex_env, lp->cplex_lp, &val, idiv, idiv) ) {
            rval = 1;  fprintf (stderr, "CPXgetslack failed\n");
            goto CLEANUP;
        }
        lb = 0.0;
        if ( sense[idiv] == 'E' )  ub = 0.0;
        else                       ub = INFBOUND;
        if ( sense[idiv] == 'G' )  val *= -1.0;
    }

    if ( CPXbinvrow (lp->cplex_env, lp->cplex_lp, i, y) ) {
        rval = 1;  fprintf (stderr, "CPXbinvrow failed\n");
        goto CLEANUP;
    }

    if ( val < lb ) {
        for (i = 0; i < nrows; i++)  y[i] *= -1.0;
    }

    for (i = 0; i < nrows; i++) {
        if ( sense[i] == 'L'  &&  y[i] > 0.0 )  y[i] = 0.0;
        if ( sense[i] == 'G'  &&  y[i] < 0.0 )  y[i] = 0.0;
    }

CLEANUP:

    CC_IFFREE (bhead, int);
    CC_IFFREE (sense, char);

    return rval;
}

#define CC_INITIAL_RSPACE 100
#define CC_INITIAL_CSPACE 100
#define CC_INITIAL_NZSPACE 2000

#define CC_SPACE_MULTIPLE 1.10

static int init_cplex40_data (CCcplex40_data **data_p, int ncols, int nrows,
        const char *name, double *obj, double *rhs, char *sense, int *matbeg,
        int *matcnt, int *matind, double *matval, double *lb, double *ub,
        int *cspace_p, int *rspace_p, int *nzspace_p)
{
    int rval = 0;
    int i, j, k;

    if (!name || !obj || !rhs || !sense || !matbeg ||
        !matcnt || !matind || !matval || !lb || !ub) {
        rval = 1;
        fprintf (stderr, "Null pointer in CCinit_lpdata\n"); goto CLEANUP;
    }

    *cspace_p = ncols;
    *rspace_p = nrows;

    /* Doing the loop backwards avoids handles the case ncols=0, and
       when ncols>0, picks a good initial value for *nzspace_p */

    *nzspace_p = 0;
    for (j = ncols - 1; j >= 0; j--) {
       if (matbeg[j]+matcnt[j] > *nzspace_p) {
           *nzspace_p = matbeg[j]+matcnt[j];
       }
    }
    (*nzspace_p)++;

    if ( *cspace_p < CC_INITIAL_CSPACE ) {
        *cspace_p = CC_INITIAL_CSPACE;
    }
    if ( *rspace_p < CC_INITIAL_RSPACE ) {
        *rspace_p = CC_INITIAL_RSPACE;
    }
    if ( *nzspace_p < CC_INITIAL_NZSPACE ) {
        *nzspace_p = CC_INITIAL_NZSPACE;
    }

    *cspace_p  *= CC_SPACE_MULTIPLE;
    *rspace_p  *= CC_SPACE_MULTIPLE;
    *nzspace_p *= CC_SPACE_MULTIPLE;

    *data_p = CC_SAFE_MALLOC (1, CCcplex40_data);
    if (!(*data_p)) {
        rval = 1;
        fprintf (stderr, "CClp_init_cplex40_data memory allocation failed\n");
        goto CLEANUP;
    }

    (*data_p)->triv_row_ind = 0;
    (*data_p)->triv_col_ind = 0;

    (*data_p)->probname = CC_SAFE_MALLOC (strlen (name)+1, char);
    (*data_p)->obj      = CC_SAFE_MALLOC (*cspace_p, double);
    (*data_p)->rhs      = CC_SAFE_MALLOC (*rspace_p, double);
    (*data_p)->sense    = CC_SAFE_MALLOC (*rspace_p, char);
    (*data_p)->matbeg   = CC_SAFE_MALLOC (*cspace_p, int);
    (*data_p)->matcnt   = CC_SAFE_MALLOC (*cspace_p, int);
    (*data_p)->matind   = CC_SAFE_MALLOC (*nzspace_p, int);
    (*data_p)->matval   = CC_SAFE_MALLOC (*nzspace_p, double);
    (*data_p)->lb       = CC_SAFE_MALLOC (*cspace_p, double);
    (*data_p)->ub       = CC_SAFE_MALLOC (*cspace_p, double);
    if ( !(*data_p)->probname ||
         !(*data_p)->rhs      ||
         !(*data_p)->sense    ||
         !(*data_p)->matbeg   ||
         !(*data_p)->matcnt   ||
         !(*data_p)->matind   ||
         !(*data_p)->matval   ||
         !(*data_p)->lb       ||
         !(*data_p)->ub         ) {
        rval = 1;
        fprintf (stderr, "init_cplex40_data memory allocation failed\n");
        goto CLEANUP;
    }

    strcpy ((*data_p)->probname, name);
    for (j = 0; j < ncols; j++) {
        (*data_p)->obj[j]    = obj[j];
        (*data_p)->matbeg[j] = matbeg[j];
        (*data_p)->matcnt[j] = matcnt[j];
        (*data_p)->lb[j]     = lb[j];
        (*data_p)->ub[j]     = ub[j];
        for (k = matbeg[j]; k < matbeg[j]+matcnt[j]; k++) {
            (*data_p)->matind[k] = matind[k];
            (*data_p)->matval[k] = matval[k];
        }
    }

    for (i = 0; i < nrows; i++) {
        (*data_p)->rhs[i]   = rhs[i];
        (*data_p)->sense[i] = sense[i];
    }

CLEANUP:

    if (rval) {
        free_cplex40_data (data_p);
        *cspace_p  = 0;
        *rspace_p  = 0;
        *nzspace_p = 0;
    }
    return rval;
}

static void free_cplex40_data (CCcplex40_data **data_p)
{
   if (*data_p) {
      CC_IFFREE ((*data_p)->probname, char);
      CC_IFFREE ((*data_p)->obj, double);
      CC_IFFREE ((*data_p)->rhs, double);
      CC_IFFREE ((*data_p)->sense, char);
      CC_IFFREE ((*data_p)->matbeg, int);
      CC_IFFREE ((*data_p)->matcnt, int);
      CC_IFFREE ((*data_p)->matind, int);
      CC_IFFREE ((*data_p)->matval, double);
      CC_IFFREE ((*data_p)->lb, double);
      CC_IFFREE ((*data_p)->ub, double);
      CC_IFFREE (*data_p, CCcplex40_data);
   }
}

static int reallocprob (CClp *lp, int cneeded, int rneeded, int nzneeded)
{
    int cspace = 0, rspace = 0, nzspace = 0;
    int csurplus = 0, rsurplus = 0, nzsurplus = 0;
    CCcplex40_data *data = lp->cplex40_data;

    if (getspace (lp, &cspace, &rspace, &nzspace)) {
        fprintf (stderr, "getspace failed\n");
        return 1;
    }
    if (getsurplus (lp, &csurplus, &rsurplus, &nzsurplus)) {
        fprintf (stderr, "getsurplus failed\n");
        return 1;
    }

    if (cneeded > csurplus) {
        cspace = (int) ((double) CC_SPACE_MULTIPLE*(cspace-csurplus+cneeded));
    }
    if (rneeded > rsurplus) {
        rspace = (int) ((double) CC_SPACE_MULTIPLE*(rspace-rsurplus+rneeded));
    }
    if (nzneeded > nzsurplus) {
        nzspace = (int) ((double) CC_SPACE_MULTIPLE*
                         (nzspace-nzsurplus+nzneeded));
    }

    if (CPXreallocprob (lp->cplex_env, lp->cplex_lp, &data->obj,
                        &data->rhs, &data->sense, &data->matbeg,
                        &data->matcnt, &data->matind, &data->matval,
                        &data->lb, &data->ub, (double **) NULL,
                        (char ***) NULL, (char **) NULL, (char ***) NULL,
                        (char **) NULL, (char **) NULL, cspace, rspace,
                        nzspace, 0, 0)) {
        fprintf (stderr, "CPXreallocprob failed\n");
        return 1;
    }

    printf ("Cplex_data reallocated:  cspace %d  rspace %d  nzspace %d\n",
            cspace, rspace, nzspace); fflush (stdout);
    return 0;
}

static int getspace (CClp *lp, int *cspace_p, int *rspace_p, int *nzspace_p)
{
    *cspace_p = CPXgetcolspace (lp->cplex_env, lp->cplex_lp);
    if (!(*cspace_p)) {
        fprintf (stderr, "CPXgetcolspace failed\n");
        return 1;
    }
    *rspace_p = CPXgetrowspace (lp->cplex_env, lp->cplex_lp);
    if (!(*rspace_p)) {
        fprintf (stderr, "CPXgetrowspace failed\n");
        return 1;
    }
    *nzspace_p = CPXgetnzspace (lp->cplex_env, lp->cplex_lp);
    if (!(*nzspace_p)) {
        fprintf (stderr, "CPXgetnzspace failed\n");
        return 1;
    }
    return 0;
}

static int getsurplus (CClp *lp, int *csurplus_p, int *rsurplus_p,
                       int *nzsurplus_p)
{
    if (CPXgetspace (lp->cplex_env, lp->cplex_lp, csurplus_p, rsurplus_p,
                     nzsurplus_p, (unsigned *) NULL, (unsigned *) NULL)) {
        fprintf (stderr, "CPXgetspace failed\n");
        return 1;
    }
    return 0;
}
