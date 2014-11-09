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

#include "machdefs.h"
#include "localcut.h"
#include "util.h"

static int
    begin_cut (void *u_data),
    add_clique (int *arr, int size, void *u_data),
    abort_cut (void *u_data),
    finish_cut (int rhs, int *finished, void *u_data);

static int begin_cut (CC_UNUSED void *u_data)
{
    return 0;
}

static int add_clique (int *arr, int size, CC_UNUSED void *u_data)
{
    int i;
    
    printf ("(");
    for (i=0; i<size; i++) {
        printf ("%d", arr[i]);
        if (i < size-1) printf (" ");
    }
    printf (") ");

    return 0;
}

static int abort_cut (CC_UNUSED void *u_data)
{
    printf ("ABORTED\n");
    fflush (stdout);
    return 0;
}

static int finish_cut (int rhs, int *finished, CC_UNUSED void *u_data)
{
    printf (">= %d\n", rhs);
    fflush (stdout);
    *finished = 0;
    return 0;
}

int main (CC_UNUSED int ac, CC_UNUSED char **av)
{
    int nodecount;
    int edgecount;
    int *elist = (int *) NULL;
    int *coef = (int *) NULL;
    int rhs;
    int i;
    int rval;
    CCchunk_cut_callback callback;
    CCtsp_lpcut_in c;

    CCtsp_init_lpcut_in (&c);
    
    scanf ("%d%d", &nodecount, &edgecount);

    elist = CC_SAFE_MALLOC (edgecount*2, int);
    coef = CC_SAFE_MALLOC (edgecount, int);
    if (elist == (int *) NULL ||
        coef == (int *) NULL) {
        fprintf (stderr, "Out of memory\n");
        rval = -1;
        goto CLEANUP;
    }

    for (i=0; i<edgecount; i++) {
        scanf ("%d%d%d",&elist[2*i], &elist[2*i+1], &coef[i]);
        coef[i] = coef[i];
    }

    scanf ("%d", &rhs);
    rhs = rhs;

    callback.begin_cut = begin_cut;
    callback.add_clique = add_clique;
    callback.abort_cut = abort_cut;
    callback.finish_cut = finish_cut;
    callback.u_data = (void *) NULL;
    
    for (i=0; i<nodecount; i++) {
        printf ("With %d outside:\n", i);
        rval = CCchunk_ineq_to_cut (nodecount, edgecount, elist, coef, rhs, i,
                            &callback);
        if (rval) {
            fprintf (stderr, "CCchunk_ineq_to_cut failed\n");
        }
        printf ("\n");
    }

    printf ("Without callback:\n");
    rval = CCchunk_ineq_to_lpcut_in (nodecount, edgecount, elist, coef, rhs, &c);
    if (rval) {
        fprintf (stderr, "CCchunk_ineq_to_lpcut_in failed\n");
    }
    CCtsp_print_lpcut_in (&c);
    
    rval = 0;

  CLEANUP:
    CC_IFFREE (elist, int);
    CC_IFFREE (coef, int);
    CCtsp_free_lpcut_in (&c);
    return rval;
}

