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
/*  int CCchunk_oracle (CCchunk_graph *ch, CCchunk_ineq *c, int *xsol,      */
/*      int *objval, int rhsvalid, int effort_limit,                        */
/*      CCchunk_oracle_timer *timer)                                        */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCchunk_verify (CCchunk_graph *ch, CCchunk_ineq *c)                 */
/*    MISSING                                                               */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "localcut.h"
#include "tinytsp.h"

#undef  DEBUG
#undef  DUMP_HARDMSPS
#undef  DUMP_FAILEDMSPS


int CCchunk_oracle (CCchunk_graph *ch, CCchunk_ineq *c, int *xsol, int *objval,
                int rhsvalid, int effort_limit, CCchunk_oracle_timer *timer)
{
    int ncount = ch->ncount + 1;
    int ecount = ch->ecount + ch->ncount;
    int *elist = (int *) NULL;
    int *weight = (int *) NULL;
    int *lbound = (int *) NULL;
    int *ubound = (int *) NULL;
    int *tmpxsol = (int *) NULL;
    int i, k;
    int rval = 0;
    double drhs, dobjval;

/*
    printf ("+"); fflush (stdout);
*/

    CCutil_start_timer (&timer->all);

#ifdef DEBUG
    printf ("CCchunk_oracle:");
    for (i=0; i<ch->ecount; i++) {
        printf (" %d", c->coef[i]);
    }
    printf (" <= %d\n", c->rhs);
    fflush (stdout);
#endif
    elist = CC_SAFE_MALLOC (ecount*2, int);
    weight = CC_SAFE_MALLOC (ecount, int);
    lbound = CC_SAFE_MALLOC (ecount, int);
    ubound = CC_SAFE_MALLOC (ecount, int);
    if (!elist || !weight || !lbound || !ubound) {
        fprintf (stderr, "Out of memory in CCchunk_oracle\n");
        rval = CC_CHUNK_ORACLE_ERROR;
        goto CLEANUP;
    }
    if (xsol) {
        tmpxsol = CC_SAFE_MALLOC (ecount, int);
        if (!tmpxsol) {
            fprintf (stderr, "Out of memory in CCchunk_oracle\n");
            rval = CC_CHUNK_ORACLE_ERROR;
            goto CLEANUP;
        }
    }
    for (i=0; i<ch->ecount; i++) {
        elist[2*i] = ch->end0[i] + 1;
        elist[2*i+1] = ch->end1[i] + 1;
        weight[i] = c->coef[i];
        if (ch->fixed[i] == 0) {
            lbound[i] = 0;
            ubound[i] = 0;
        } else if (ch->fixed[i] == 1) {
            lbound[i] = 1;
            ubound[i] = 1;
        } else {
            lbound[i] = 0;
            ubound[i] = 1;
        }
        if (tmpxsol) tmpxsol[i] = 0;
    }
    for (i=0; i<ch->ncount; i++) {
        k = ch->ecount + i;
        elist[2*k] = 0;
        elist[2*k+1] = i+1;
        weight[k] = 0;
        if (ch->equality[i]) {
            lbound[k] = 0;
            ubound[k] = 0;
        } else {
            lbound[k] = 0;
            ubound[k] = 2;
        }
        if (tmpxsol) tmpxsol[k] = 0;
    }

    drhs = c->rhs;
    CCutil_start_timer (&timer->bnbtsp);
    rval = CCtiny_bnb_msp (ncount, ecount, elist, weight, 0 /* depot */,
                     lbound, ubound, (rhsvalid ? &drhs : (double *) NULL),
                     CC_TINYTSP_MAXIMIZE, &dobjval, tmpxsol,
                     20 * effort_limit /* searchlimit */);
    CCutil_stop_timer (&timer->bnbtsp, 0);
    if (rval == CC_TINYTSP_INFEASIBLE) {
        rval = CC_CHUNK_ORACLE_INFEASIBLE;
        goto CLEANUP;
    } else if (rval == 0 && rhsvalid && dobjval <= drhs) {
        fprintf (stderr, "bnbtsp obj %.0f c->rhs %.0f, should have reported infeas\n",
                 dobjval, drhs);
        rval = CC_CHUNK_ORACLE_INFEASIBLE;
        goto CLEANUP;
    } else if (rval == 0) {
        if (objval) *objval = dobjval;
        if (xsol) {
            for (i=0; i<ch->ecount; i++) {
                xsol[i] = tmpxsol[i];
            }
        }
        goto CLEANUP;
    } else if (rval != CC_TINYTSP_SEARCHLIMITEXCEEDED) {
        fprintf (stderr, "CCtiny_bnb_msp rval %d\n", rval);
        rval = CC_CHUNK_ORACLE_ERROR;
#ifdef DUMP_FAILEDMSPS
        printf ("%d %d\n",ncount, ecount);
        for (i=0; i<ecount; i++) {
            printf ("%d %d %d %d %d\n",elist[2*i],elist[2*i+1],lbound[i],
                    ubound[i],weight[i]);
        }
        if (rhsvalid) printf ("%d\n", c->rhs);
        printf ("\n");
        fflush (stdout);
#endif
        goto CLEANUP;
    }

    drhs = c->rhs;
    CCutil_start_timer (&timer->tinytsp);
    rval = CCtiny_bnc_msp (ncount, ecount, elist, weight, 0 /* depot */,
                    lbound, ubound, (rhsvalid ? &drhs : (double *) NULL),
                    CC_TINYTSP_MAXIMIZE, &dobjval, tmpxsol, 1 /* checkresult */,
                    1 * effort_limit /* searchlimit */);
    CCutil_stop_timer (&timer->tinytsp, 0);
    if (rval == CC_TINYTSP_INFEASIBLE) {
        rval = CC_CHUNK_ORACLE_INFEASIBLE;
        goto CLEANUP;
    } else if (rval == 0 && rhsvalid && dobjval <= drhs) {
        fprintf (stderr, "CCtiny_bnc_msp obj %.0f c->rhs %.0f, should have reported infeas\n",
                 dobjval, drhs);
        rval = CC_CHUNK_ORACLE_INFEASIBLE;
        goto CLEANUP;
    } else if (rval == 0) {
        if (objval) *objval = dobjval;
        if (xsol) {
            for (i=0; i<ch->ecount; i++) {
                xsol[i] = tmpxsol[i];
            }
        }
        goto CLEANUP;
    } else if (rval != CC_TINYTSP_SEARCHLIMITEXCEEDED) {
        fprintf (stderr, "CCtiny_bnc_msp rval %d\n", rval);
        rval = CC_CHUNK_ORACLE_ERROR;
#ifdef DUMP_FAILEDMSPS
        printf ("%d %d\n",ncount, ecount);
        for (i=0; i<ecount; i++) {
            printf ("%d %d %d %d %d\n",elist[2*i],elist[2*i+1],lbound[i],
                    ubound[i],weight[i]);
        }
        if (rhsvalid) printf ("%d\n", c->rhs);
        printf ("\n");
        fflush (stdout);
#endif
        goto CLEANUP;
    }

    rval = CC_CHUNK_ORACLE_SEARCHLIMITEXCEEDED;
#ifdef DUMP_HARDMSPS
    printf ("TSPORACLE Search Limit Exceeded");
    if (rhsvalid) printf ("TSP objlimit %d", c->rhs);
    printf (":\n");
    printf ("%d %d\n",ncount, ecount);
    for (i=0; i<ecount; i++) {
        printf ("%d %d %d %d %d\n",elist[2*i],elist[2*i+1],lbound[i],
                ubound[i],weight[i]);
    }
    if (rhsvalid) printf ("%d\n", c->rhs);
    printf ("\n");
    fflush (stdout);
#endif

    printf ("-"); fflush (stdout);

  CLEANUP:
    CC_IFFREE (elist, int);
    CC_IFFREE (weight, int);
    CC_IFFREE (lbound, int);
    CC_IFFREE (ubound, int);
    CC_IFFREE (tmpxsol, int);
    CCutil_stop_timer (&timer->all, 0);
    return rval;
}

int CCchunk_verify (CCchunk_graph *ch, CCchunk_ineq *c)
{
    int ncount = ch->ncount + 1;
    int ecount = ch->ecount + ch->ncount;
    int *elist = (int *) NULL;
    int *weight = (int *) NULL;
    int *lbound = (int *) NULL;
    int *ubound = (int *) NULL;
    int i, k;
    int rval = 0;
    double drhs, dobjval;

#ifdef DEBUG
    printf ("CCchunk_verify:");
    for (i=0; i<ch->ecount; i++) {
        printf (" %d", c->coef[i]);
    }
    printf (" <= %d\n", c->rhs);
    fflush (stdout);
#endif
    elist = CC_SAFE_MALLOC (ecount*2, int);
    weight = CC_SAFE_MALLOC (ecount, int);
    lbound = CC_SAFE_MALLOC (ecount, int);
    ubound = CC_SAFE_MALLOC (ecount, int);
    if (!elist || !weight || !lbound || !ubound) {
        fprintf (stderr, "Out of memory in CCchunk_oracle\n");
        rval = CC_CHUNK_ORACLE_ERROR;
        goto CLEANUP;
    }

    for (i=0; i<ch->ecount; i++) {
        elist[2*i] = ch->end0[i] + 1;
        elist[2*i+1] = ch->end1[i] + 1;
        weight[i] = c->coef[i];
        lbound[i] = 0;
        ubound[i] = 1;
    }
    for (i=0; i<ch->ncount; i++) {
        k = ch->ecount + i;
        elist[2*k] = 0;
        elist[2*k+1] = i+1;
        weight[k] = 0;
        lbound[k] = 0;
        ubound[k] = 2;
    }

    drhs = c->rhs;
    rval = CCtiny_bnc_msp (ncount, ecount, elist, weight, 0 /* depot */,
                    lbound, ubound, &drhs, CC_TINYTSP_MAXIMIZE,
                    &dobjval, (int *) NULL /* xsol */, 1 /* checkresult */,
                    2000 /* searchlimit */);

    if (rval == CC_TINYTSP_INFEASIBLE) {
        rval = 0;
    } else if (rval == 0 && dobjval <= drhs) {
        fprintf (stderr, "CCtiny_bnc_msp obj %.0f c->rhs %.0f, should have reported infeas\n",
                 dobjval, drhs);
        rval = 0;
    } else if (rval == 0) {
        rval = -1;
    }

    if (rval) {
        printf ("TSPVERIFY Failed rval %d objlimit %d:\n", rval, c->rhs);
        printf ("%d %d\n",ncount, ecount);
        for (i=0; i<ecount; i++) {
            printf ("%d %d %d %d %d\n",elist[2*i],elist[2*i+1],lbound[i],
                    ubound[i],weight[i]);
        }
        printf ("\n");
    }

  CLEANUP:
    CC_IFFREE (elist, int);
    CC_IFFREE (weight, int);
    CC_IFFREE (lbound, int);
    CC_IFFREE (ubound, int);
    return rval;
}

