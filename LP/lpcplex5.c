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
/*              Interface Routines to The CPLEX 5.0-6.5 LP Solvers          */
/*                                                                          */
/*  NOTE: Use this code in place of lp_none.c to access the Cplex 5.0-6.5   */
/*   libraries. You will also need to link the cplex library via the        */
/*   makefile.                                                              */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: April 17, 1997                                                    */
/*        June 19, 1997 (bico, REB)                                         */
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

#undef  CC_CPLEX_DISPLAY

#define CC_ONE_ENV

#ifdef CC_ONE_ENV
static CPXENVptr CClp_cplex_env = (CPXENVptr) NULL;
static int CClp_env_count = 0;

typedef struct CClp_parameters {
    int    scrind;
    int    simdisplay;
    int    fastmip;
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
#endif

struct CClp {
    CPXENVptr cplex_env;
    CPXLPptr  cplex_lp;
#ifdef CC_ONE_ENV
    CClp_parameters cplex_params;
#endif
};

struct CClp_warmstart {
    int      rcount;
    int      ccount;
    int     *rstat;
    int     *cstat;
    double  *dnorm;
};

struct CClp_info {
    int  rcount;
    int  ccount;
    int *rstat;
    int *cstat;
};

#define SOLVER_WARMSTART_NAME "CPL5"


static int
    primalopt (CClp *lp),
    dualopt (CClp *lp),
    baropt (CClp *lp),
    getfarkasmultipliers (CClp *lp, double *y);

#ifdef CC_ONE_ENV
static int
    set_parameters (CPXENVptr cplex_env, CClp_parameters *params);
#endif


int CClp_init (CClp **lp)
{
    int rval = 0;

    CClp_free (lp);

    (*lp) = CC_SAFE_MALLOC (1, CClp);
    if ((*lp) == (CClp *) NULL) {
        fprintf (stderr, "Out of memory in CClp_init\n");
        rval = 1; goto CLEANUP;
    }

    (*lp)->cplex_env = (CPXENVptr) NULL;
    (*lp)->cplex_lp = (CPXLPptr) NULL;

#ifdef CC_ONE_ENV
    if (CClp_cplex_env == (CPXENVptr) NULL) {
        CClp_cplex_env = CPXopenCPLEXdevelop (&rval);
        if (rval) {
            fprintf (stderr, "CPXopenCPLEXdevelop failed, return code %d\n", rval);
            goto CLEANUP;
        }
        CClp_env_count = 0;
    }
    (*lp)->cplex_env = CClp_cplex_env;
    CClp_env_count++;
#else /* CC_ONE_ENV */
    (*lp)->cplex_env = CPXopenCPLEXdevelop (&rval);
    if (rval) {
        fprintf (stderr, "CPXopenCPLEXdevelop failed, return code %d\n", rval);
        goto CLEANUP;
    }
#endif /* CC_ONE_ENV */

#ifdef CC_ONE_ENV
#ifdef CC_CPLEX_DISPLAY
    (*lp)->cplex_params.scrind = 1;
    (*lp)->cplex_params.simdisplay = 1;
#else
    (*lp)->cplex_params.scrind = 0;
    (*lp)->cplex_params.simdisplay = 0;
#endif
    (*lp)->cplex_params.fastmip = 1;
    (*lp)->cplex_params.advind = 1;
    (*lp)->cplex_params.dpriind = CPX_DPRIIND_STEEP;
    (*lp)->cplex_params.ppriind = CPX_PPRIIND_STEEP;
    (*lp)->cplex_params.epper = 1.0E-6;
    (*lp)->cplex_params.epopt = 1.0E-9;
    (*lp)->cplex_params.eprhs = 1.0E-9;

    (*lp)->cplex_params.perind = 0;
    (*lp)->cplex_params.preind = 1;
    (*lp)->cplex_params.aggind = 1;

#else /* CC_ONE_ENV */

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

    rval = CPXsetintparam ((*lp)->cplex_env, CPX_PARAM_FASTMIP, 1);
    if (rval) {
        fprintf (stderr, "CPXsetintparam CPX_PARAM_FASTMIP failed\n");
        goto CLEANUP;
    }

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
#endif /* CC_ONE_ENV */

CLEANUP:

    if (rval) {
        CClp_free (lp);
    }

    return rval;
}

int CClp_force_perturb (CClp *lp)
{
    int rval = 0;

#ifdef CC_ONE_ENV
    lp->cplex_params.perind = 1;
#else /* CC_ONE_ENV */
    rval = CPXsetintparam (lp->cplex_env, CPX_PARAM_PERIND, 1);
    if (rval) {
        fprintf (stderr, "CPXsetintparam CPX_PARAM_PERIND failed\n");
    }
#endif /* CC_ONE_ENV */
    return rval;
}

int CClp_tune_small (CClp *lp)
{
#ifdef CC_ONE_ENV
    lp->cplex_params.dpriind = CPX_DPRIIND_FULL;
    lp->cplex_params.ppriind = CPX_PPRIIND_AUTO;
#else /* CC_ONE_ENV */
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
#endif /* CC_ONE_ENV */
    return 0;
}

int CClp_disable_presolve (CClp *lp)
{
#ifdef CC_ONE_ENV
    lp->cplex_params.preind = CPX_OFF;
    lp->cplex_params.aggind = 0;
#else /* CC_ONE_ENV */
    if (CPXsetintparam (lp->cplex_env, CPX_PARAM_PREIND, CPX_OFF)) {
        fprintf (stderr, "CPXsetintparam CPX_PARAM_PREIND failed\n");
        return 1;
    }
    if (CPXsetintparam (lp->cplex_env, CPX_PARAM_AGGIND, 0)) {
        fprintf (stderr, "CPXsetintparam CPX_PARAM_AGGIND failed\n");
        return 1;
    }
#endif /* CC_ONE_ENV */
    return 0;
}

void CClp_free (CClp **lp)
{
    if (*lp) {
        if ((*lp)->cplex_env) {
            if ((*lp)->cplex_lp) {
                CPXfreeprob ((*lp)->cplex_env, &((*lp)->cplex_lp));
            }
#ifdef CC_ONE_ENV
            (*lp)->cplex_env = (CPXENVptr) NULL;
            CClp_env_count--;
            if (CClp_env_count == 0) {
                CPXcloseCPLEX (&CClp_cplex_env);
                CClp_cplex_env = (CPXENVptr) NULL;
            }
#else
            CPXcloseCPLEX (&((*lp)->cplex_env));
#endif
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
    }
}

int CClp_loadlp (CClp *lp, const char *name, int ncols, int nrows,
        int objsense, double *obj, double *rhs, char *sense, int *matbeg,
        int *matcnt, int *matind, double *matval, double *lb, double *ub)
{
    int rval = 0;
    char nambuf[32];

    /* We copy the name since CPXcreateprob doesn't declare name as const */
    strncpy (nambuf, name, sizeof (nambuf));
    nambuf[sizeof(nambuf)-1] = '\0';

    lp->cplex_lp = CPXcreateprob (lp->cplex_env, &rval, nambuf);
    if (!lp->cplex_lp || rval) {
       fprintf (stderr, "CPXcreateprob failed, return code %d\n", rval);
       return 1;
    }

    rval = CPXcopylp (lp->cplex_env, lp->cplex_lp, ncols, nrows,
                      objsense, obj, rhs, sense, matbeg, matcnt,
                      matind, matval, lb, ub, (double *) NULL);
    if (rval) {
       fprintf (stderr, "CPXcopylp failed, return code %d\n", rval);
       return 1;
    }

    return 0;
}

int CClp_create (CClp *lp, const char *name)
{
    int rval;
    char nambuf[32];

    /* We copy the name since CPXcreateprob doesn't declare name as const */
    strncpy (nambuf, name, sizeof (nambuf));
    nambuf[sizeof(nambuf)-1] = '\0';

    lp->cplex_lp = CPXcreateprob (lp->cplex_env, &rval, nambuf);
    if (!lp->cplex_lp || rval) {
       fprintf (stderr, "CPXcreateprob failed, return code %d\n", rval);
       return 1;
    }
    return 0;
}

int CClp_new_row (CClp *lp, char sense, double rhs)
{
    int rval;
    char asense[1];
    double arhs[1];

    asense[0] = sense;
    arhs[0] = rhs;

    rval = CPXnewrows (lp->cplex_env, lp->cplex_lp, 1, arhs, asense,
                       (double *) NULL, (char **) NULL);
    if (rval) {
        fprintf (stderr, "CPXnewrows failed\n");
        return rval;
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

#ifdef CC_ONE_ENV
    if (set_parameters (lp->cplex_env, &lp->cplex_params)) {
        fprintf (stderr, "Unable to set optimization parameters\n");
    }
#endif

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

    rval = CPXprimopt (lp->cplex_env, lp->cplex_lp);
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
            rval = CPXprimopt (lp->cplex_env, lp->cplex_lp);
            if (rval) {
                fprintf (stderr, "CPXprimopt failed, return code %d\n", rval);
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
            fprintf (stderr, "CPXprimopt failed, return code %d\n", rval);
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
                fprintf (stderr, "CPXdualopt failed, return code %d\n", rval);
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
            fprintf (stderr, "CPXdualopt failed, return code %d\n", rval);
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
        printf ("CPXbaropt failed, return code %d, calling CPXdualopt\n",
                rval);
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

#ifdef CC_ONE_ENV
    if (set_parameters (lp->cplex_env, &lp->cplex_params)) {
        fprintf (stderr, "Unable to set optimization parameters\n");
    }
#endif

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
                fprintf (stderr, "CPXdualopt failed, return code %d\n", rval);
                goto CLEANUP;
            }
        } else {
            fprintf (stderr, "CPXdualopt failed, return code %d\n", rval);
            goto CLEANUP;
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
            fprintf (stderr, "CPXdualopt failed, return code %d\n", rval);
            goto CLEANUP;
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

    rval = CPXaddrows (lp->cplex_env, lp->cplex_lp, 0, newrows, newnz,
                       rhs, sense, rmatbeg, rmatind, rmatval,
                       (char **) NULL, (char **) NULL);
    if (rval) fprintf (stderr, "CPXaddrows failed\n");
    return rval;
}

int CClp_addcols (CClp *lp, int newcols, int newnz, double *obj,
                  int *cmatbeg, int *cmatind, double *cmatval,
                  double *lb, double *ub)
{
    int rval = 0;

    rval = CPXaddcols (lp->cplex_env, lp->cplex_lp, newcols, newnz, obj,
                   cmatbeg, cmatind, cmatval, lb, ub, (char **) NULL);
    if (rval) fprintf (stderr, "CPXaddcols failed\n");
    return rval;
}

int CClp_delete_row (CClp *lp, int i)
{
    int rval = 0;
    int locali[1];

#ifdef CC_ONE_ENV
    if (set_parameters (lp->cplex_env, &lp->cplex_params)) {
        fprintf (stderr, "Unable to set optimization parameters\n");
    }
#endif

    locali[0] = i;
    if (CPXpivotin (lp->cplex_env, lp->cplex_lp, locali, 1)) {
/*        fprintf (stderr, "CPXpivotin failed, continuing anyway\n");*/
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

#ifdef CC_ONE_ENV
    if (set_parameters (lp->cplex_env, &lp->cplex_params)) {
        fprintf (stderr, "Unable to set optimization parameters\n");
    }
#endif

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
    int locali[1];
    char lu[1];
    double bd[1];

#ifdef CC_ONE_ENV
    if (set_parameters (lp->cplex_env, &lp->cplex_params)) {
        fprintf (stderr, "Unable to set optimization parameters\n");
    }
#endif

    locali[0] = i;
    lu[0] = 'B';
    bd[0] = 0.0;

    if (CPXchgbds (lp->cplex_env, lp->cplex_lp, 1, locali, lu, bd)) {
        fprintf (stderr, "CPXchgbds failed, continuing anyway\n");
    }

    if (CPXdualopt (lp->cplex_env, lp->cplex_lp)) {
        fprintf (stderr, "CPXdualopt failed, continuing anyway\n");
    }

    if (CPXpivotout (lp->cplex_env, lp->cplex_lp, locali, 1)) {
        fprintf (stderr, "CPXpivotout failed, continuing anyway\n");
    }

    rval = CPXdelcols (lp->cplex_env, lp->cplex_lp, i, i);
    if (rval) fprintf (stderr, "CPXdelcols failed\n");
    return rval;
}

int CClp_delete_set_of_columns (CClp *lp, int *delstat)
{
    int rval = 0;
    int *dellist = (int *) NULL;
    char *lu = (char *) NULL;
    double *bd = (double *) NULL;
    int delcnt = 0;
    int i;
    int j;
    int ccnt = CPXgetnumcols (lp->cplex_env, lp->cplex_lp);

#ifdef CC_ONE_ENV
    if (set_parameters (lp->cplex_env, &lp->cplex_params)) {
        fprintf (stderr, "Unable to set optimization parameters\n");
    }
#endif

    for (i=0; i<ccnt; i++) {
        if (delstat[i]) delcnt++;
    }
    if (delcnt == 0) {
        fprintf (stderr, "delete_set_of_columns with no deleted columns\n");
        return 0;
    }
    dellist = CC_SAFE_MALLOC (delcnt, int);
    lu = CC_SAFE_MALLOC (delcnt, char);
    bd = CC_SAFE_MALLOC (delcnt, double);
    if (dellist == (int *) NULL ||
        lu == (char *) NULL ||
        bd == (double *) NULL) {
        fprintf (stderr, "Out of memory in delete_set_of_columns\n");
        CC_IFFREE (dellist, int);
        CC_IFFREE (lu, char);
        CC_IFFREE (bd, double);
        return 1;
    }
    for (i=0, j=0; i<ccnt; i++) {
        if (delstat[i]) {
            lu[j] = 'B';
            bd[j] = 0.0;
            dellist[j++] = i;
        }
    }
    if (j != delcnt) {
        fprintf (stderr, "Lost some deleted columns\n");
        CC_FREE (dellist, int);
        CC_FREE (lu, char);
        CC_FREE (bd, double);
        return 1;
    }

    if (CPXchgbds (lp->cplex_env, lp->cplex_lp, delcnt, dellist, lu, bd)) {
        fprintf (stderr, "CPXchgbds failed, stumbling on anyway\n");
    }
    
    if (CPXdualopt (lp->cplex_env, lp->cplex_lp)) {
        fprintf (stderr, "CPXdualopt failed, continuing anyway\n");
    }

    if (CPXpivotout (lp->cplex_env, lp->cplex_lp, dellist, delcnt)) {
        fprintf (stderr, "CPXpivotout failed, continuing anyway\n");
    }

    CC_FREE (dellist, int);
    CC_FREE (lu, char);
    CC_FREE (bd, double);
    
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
    if (!(*w)->cstat || !(*w)->rstat || !(*w)->dnorm) {
        fprintf (stderr, "out of memory in CClp_get_warmstart\n");
        rval = 1; goto CLEANUP;
    }

    rval = CPXgetbasednorms (lp->cplex_env, lp->cplex_lp, (*w)->cstat,
                             (*w)->rstat, (*w)->dnorm);
    if (rval) {
        fprintf (stderr, "CPXgetbasednorms failed, trying to get basis\n");
        CC_IFFREE ((*w)->dnorm, double);
        rval = CPXgetbase (lp->cplex_env, lp->cplex_lp, (*w)->cstat,
                           (*w)->rstat);
        if (rval) {
            fprintf (stderr, "CPXgetbase failed\n"); goto CLEANUP;
        }
    }

    return 0;

CLEANUP:

    CClp_free_warmstart (w);
    return rval;
}

int CClp_load_warmstart (CClp *lp, CClp_warmstart *w)
{
    int rval = 0;

    if (w->cstat && w->rstat && w->dnorm) {
        rval = CPXcopybasednorms (lp->cplex_env, lp->cplex_lp, w->cstat,
                                  w->rstat, w->dnorm);
        if (rval) {
            fprintf (stderr, "CPXcopybasednorms failed\n"); goto CLEANUP;
        }
    } else if (w->cstat && w->rstat) {
        rval = CPXloadbase (lp->cplex_env, lp->cplex_lp, w->cstat, w->rstat);
        if (rval) {
            fprintf (stderr, "CPXloadbase failed\n"); goto CLEANUP;
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
    if (*w != (CClp_warmstart *) NULL) {
        CC_IFFREE ((*w)->cstat, int);
        CC_IFFREE ((*w)->rstat, int);
        CC_IFFREE ((*w)->dnorm, double);
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

    (*w)->cstat = CC_SAFE_MALLOC (ccount, int);
    (*w)->rstat = CC_SAFE_MALLOC (rcount, int);
    if ((*w)->cstat == (int *) NULL ||
        (*w)->rstat == (int *) NULL) {
        fprintf (stderr, "out of memory in CClp_sread_warmstart\n");
        goto CLEANUP;
    }
    for (i = 0; i < ccount; i++) {
        if (CCutil_sread_bits (f, &(((*w)->cstat)[i]), 2))
            goto CLEANUP;
    }
    for (i = 0; i < rcount; i++) {
        if (CCutil_sread_bits (f, &(((*w)->rstat)[i]), 1))
            goto CLEANUP;
    }

    if (CCutil_sread_int (f, &has_dnorms)) goto CLEANUP;

    if (has_dnorms) {
        (*w)->dnorm = CC_SAFE_MALLOC (rcount, double);
        if ((*w)->dnorm == (double *) NULL) {
            fprintf (stderr, "out of memory in CClp_sread_warmstart\n");
            goto CLEANUP;
        }
        for (i = 0; i < rcount; i++) {
            if (CCutil_sread_double (f, &(((*w)->dnorm)[i]))) goto CLEANUP;
        }
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

#ifdef CC_ONE_ENV
    if (set_parameters (lp->cplex_env, &lp->cplex_params)) {
        fprintf (stderr, "Unable to set optimization parameters\n");
    }
#endif

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

#ifdef CC_ONE_ENV
    if (set_parameters (lp->cplex_env, &lp->cplex_params)) {
        fprintf (stderr, "Unable to set optimization parameters\n");
    }
#endif

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
    int  rval = 0;
    int  ncols, i, k;
    int  *cstat = (int *) NULL;
    double *x = (double *) NULL;

    /* Call CPXdualopt and verify optimality */

#ifdef CC_ONE_ENV
    if (set_parameters (lp->cplex_env, &lp->cplex_params)) {
        fprintf (stderr, "Unable to set optimization parameters\n");
    }
#endif

    rval = CPXdualopt (lp->cplex_env, lp->cplex_lp);
    if (rval) {
        fprintf (stderr, "CPXdualopt failed, return code %d\n", rval);
        rval = 1; goto CLEANUP;
    }

    ncols = CPXgetnumcols (lp->cplex_env, lp->cplex_lp);
    if ( ncols == 0 ) {
        fprintf (stderr, "No columns in LP\n");
        rval = 1; goto CLEANUP;
    }

    x = CC_SAFE_MALLOC (ncols, double);
    if (x == (double *) NULL) {
        fprintf (stderr, "out of memory in branch_getgoodlist\n");
        rval = 1; goto CLEANUP;
    }
    if (CPXgetx (lp->cplex_env, lp->cplex_lp, x, 0, ncols-1)) {
        fprintf (stderr, "CPXgetx failed\n");
        rval = 1; goto CLEANUP;
    }

    cstat = CC_SAFE_MALLOC (ncols, int);
    if ( cstat == (int *) NULL ) {
        fprintf (stderr, "Out of memory\n");
        rval = 1; goto CLEANUP;
    }

    /* Get basis */

    if ( CPXgetbase (lp->cplex_env, lp->cplex_lp, cstat, (int *) NULL) ) {
        fprintf (stderr, "CPXgetbase failed\n");
        rval = 1; goto CLEANUP;
    }

    /* Make initial goodlist and goodlen */

    *goodlen_p = 0;
    for (i = 0; i < ncols; i++) {
       if ( cstat[i] == 1 ) {
          goodlist[(*goodlen_p)++] = i;
       }
    }

    /* Call CPXmdleave */

    if ( CPXmdleave (lp->cplex_env, lp->cplex_lp, goodlist, *goodlen_p,
                     downpen, uppen)) {
       fprintf (stderr, "CPXmdleave failed\n");
       rval = 1; goto CLEANUP;
    }

    /* Keep only the nondegenerate ones */

    k = *goodlen_p;
    *goodlen_p = 0;
    for (i = 0; i < k; i++) {
       if ( CC_OURABS (downpen[i]) > OURCPLEXZERO   &&
            CC_OURABS (uppen[i])   > OURCPLEXZERO   &&
            x[goodlist[i]] >= OURCPLEX_INTTOL      &&
            x[goodlist[i]] <= 1.0 - OURCPLEX_INTTOL  ) {
          goodlist[*goodlen_p]  = goodlist[i];
          downpen[*goodlen_p]   = x[goodlist[i]] * downpen[i];
          uppen[*goodlen_p]     = (1.0 - x[goodlist[i]]) * uppen[i];
          (*goodlen_p)++;
       }
    }

    if (*goodlen_p == 0) {
        /* All edges have degenerate pivots */
        for (i = 0; i < k; i++) {
            if (x[goodlist[i]] >= OURCPLEX_INTTOL &&
                x[goodlist[i]] <= 1.0 - OURCPLEX_INTTOL) {
                goodlist[*goodlen_p] = goodlist[i];
                downpen[*goodlen_p]  = x[goodlist[i]];
                uppen[*goodlen_p]    = 1.0 - x[goodlist[i]];
                (*goodlen_p)++;
            }
        }
    }

    rval = 0;
    
CLEANUP:

    CC_IFFREE (cstat, int);
    CC_IFFREE (x, double);
    return rval;
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

#ifdef CC_ONE_ENV
    if (set_parameters (lp->cplex_env, &lp->cplex_params)) {
        fprintf (stderr, "Unable to set optimization parameters\n");
    }
#endif

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
        fprintf (stderr, "CPXstrongbranch failed, return code %d\n", rval);
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

#ifdef CC_ONE_ENV
    if (set_parameters (lp->cplex_env, &lp->cplex_params)) {
        fprintf (stderr, "Unable to set optimization parameters\n");
    }
#endif

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

#ifdef CC_ONE_ENV

static int set_parameters (CPXENVptr cplex_env, CClp_parameters *params)
{
    int rval;
    
    /* the documentation doesn't say what the return value means */
    rval = CPXsetintparam (cplex_env, CPX_PARAM_SCRIND, params->scrind);
    if (rval) {
        fprintf (stderr, "CPXsetintparam CPX_PARAM_SCRIND failed\n");
        goto CLEANUP;
    }

    rval = CPXsetintparam (cplex_env, CPX_PARAM_SIMDISPLAY,
                           params->simdisplay);
    if (rval) {
        fprintf (stderr, "CPXsetintparam CPX_PARAM_SIMDISPLAY failed\n");
        goto CLEANUP;
    }

    rval = CPXsetintparam (cplex_env, CPX_PARAM_FASTMIP, params->fastmip);
    if (rval) {
        fprintf (stderr, "CPXsetintparam CPX_PARAM_FASTMIP failed\n");
        goto CLEANUP;
    }
    rval = CPXsetintparam (cplex_env, CPX_PARAM_ADVIND, params->advind);
    if (rval) {
        fprintf (stderr, "CPXsetintparam CPX_PARAM_ADVIND failed\n");
        goto CLEANUP;
    }
    rval = CPXsetintparam (cplex_env, CPX_PARAM_DPRIIND, params->dpriind);
    if (rval) {
        fprintf (stderr, "CPXsetintparam CPX_PARAM_DPRIIND failed\n");
        goto CLEANUP;
    }
    rval = CPXsetintparam (cplex_env, CPX_PARAM_PPRIIND, params->ppriind);
    if (rval) {
        fprintf (stderr, "CPXsetintparam CPX_PARAM_PPRIIND failed\n");
        goto CLEANUP;
    }

    rval = CPXsetdblparam (cplex_env, CPX_PARAM_EPPER, params->epper);
    if (rval) {
        fprintf (stderr, "CPXsetdblparam CPX_PARAM_EPPER failed\n");
        goto CLEANUP;
    }
    rval = CPXsetdblparam (cplex_env, CPX_PARAM_EPOPT, params->epopt);
    if (rval) {
        fprintf (stderr, "CPXsetdblparam CPX_PARAM_EPOPT failed\n");
        goto CLEANUP;
    }
    rval = CPXsetdblparam (cplex_env, CPX_PARAM_EPRHS, params->eprhs);
    if (rval) {
        fprintf (stderr, "CPXsetdblparam CPX_PARAM_EPRHS failed\n");
        goto CLEANUP;
    }
    
    rval = CPXsetintparam (cplex_env, CPX_PARAM_PERIND, params->perind);
    if (rval) {
        fprintf (stderr, "CPXsetintparam CPX_PARAM_PERIND failed\n");
        goto CLEANUP;
    }

    rval = CPXsetintparam (cplex_env, CPX_PARAM_PREIND, params->preind);
    if (rval) {
        fprintf (stderr, "CPXsetintparam CPX_PARAM_PREIND failed\n");
        goto CLEANUP;
    }

    rval = CPXsetintparam (cplex_env, CPX_PARAM_AGGIND, params->aggind);
    if (rval) {
        fprintf (stderr, "CPXsetintparam CPX_PARAM_AGGIND failed\n");
        goto CLEANUP;
    }

    rval = 0;
  CLEANUP:
    return rval;
}
#endif /* CC_ONE_ENV */
