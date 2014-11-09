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

typedef struct CCchunk_cut_callback_data {
    int all_cuts;
    int cut_count;
    double cutval;
    CCchunk_graph *chunk;
} CCchunk_cut_callback_data;

typedef struct CCchunk_fault_callback_data {
    CCchunk_lift_timer *lift_t;
    int fault_count;
} CCchunk_fault_callback_data;


static int
    begin_CCchunk_cut_callback (void *u_data),
    add_clique_callback (int *arr, int size, void *u_data),
    abort_CCchunk_cut_callback (void *u_data),
    finish_CCchunk_cut_callback (int rhs, int *finished, void *u_data),
    found_CCchunk_fault_callback (CCchunk_graph *chunk, CCchunk_fault *fault,
        int *finished, void *u_data);



int main (CC_UNUSED int ac, CC_UNUSED char **av)
{
    CCchunk_graph *c = (CCchunk_graph *) NULL;
    CCchunk_fault f;
    int rval;
    int i;
    int j;
    int ncount, ecount;
    CCchunk_cut_callback ccb;
    CCchunk_fault_callback fcb;
    CCchunk_lift_timer lift_t;
    CCchunk_separate_timer separate_t;
    CCchunk_cut_callback_data ccd;
    CCchunk_fault_callback_data fcd;

    f.a.coef = (int *) NULL;
    f.nsols = 0;
    f.sols = (int *) NULL;

    ccd.cut_count = 0;
    ccd.all_cuts = 1;
    
    ccb.begin_cut = begin_CCchunk_cut_callback;
    ccb.add_clique = add_clique_callback;
    ccb.abort_cut = abort_CCchunk_cut_callback;
    ccb.finish_cut = finish_CCchunk_cut_callback;
    ccb.u_data = (void *) &ccd;

    fcd.fault_count = 0;
    fcd.lift_t = &lift_t;

    fcb.func = found_CCchunk_fault_callback;
    fcb.u_data = (void *) &fcd;

    CCchunk_init_separate_timer (&separate_t);
    CCchunk_init_lift_timer (&lift_t);
    
    scanf("%d%d",&ncount,&ecount);
    c = CCchunk_graph_alloc (ncount, ecount);
    if (!c) {
        fprintf (stderr, "Unable to allocate CCchunk_graph\n");
        rval = 1; goto CLEANUP;
    }

    for (i=0; i<ecount; i++) {
        scanf ("%d%d%lf%d", &c->end0[i],&c->end1[i], &c->weight[i],
               &c->fixed[i]);
    }

    for (i=0; i<c->ncount; i++) {
        scanf ("%d", &c->equality[i]);
    }
        
    for (i=0; i<c->ncount; i++) {
        c->members[i] = CC_SAFE_MALLOC (2, int);
        if (!c->members[i]) {
            fprintf (stderr, "Unable to allocate members\n");
            rval = 1; goto CLEANUP;
        }
        c->members[i][0] = i;
        c->members[i][1] = -1;
    }

    scanf (" ");
    if (!feof (stdin)) {
        f.a.coef = CC_SAFE_MALLOC (c->ecount, int);
        if (f.a.coef == (int *) NULL) {
            fprintf (stderr, "Out of memory\n");
            rval = 1; goto CLEANUP;
        }

        for (i=0; i<ecount; i++) {
            scanf ("%d", &f.a.coef[i]);
        }
        scanf (" <= %d", &f.a.rhs);
        scanf ("%d", &f.nsols);
        if (f.nsols > 0) {
            f.sols = CC_SAFE_MALLOC (c->ecount * f.nsols, int);
            if (f.sols == (int *) NULL) {
                fprintf (stderr, "Out of memory\n");
                rval = 1; goto CLEANUP;
            }
            for (i=0; i<f.nsols; i++) {
                for (j=0; j<c->ecount; j++) {
                    scanf ("%d", &f.sols[i*c->ecount+j]);
                }
            }
        }

        ccd.chunk = c;

        printf ("lifting...\n"); fflush (stdout);
        rval = CCchunk_lift (c, &f, &lift_t, &ccb);
        if (rval) {
            fprintf (stderr, "CCchunk_lift failed, return code %d\n", rval);
        }
        CCchunk_print_lift_timer (&lift_t);
    } else {
        printf ("separating...\n"); fflush (stdout);
        rval = CCchunk_separate (c, &separate_t, &fcb);
        if (rval) {
            fprintf (stderr, "chunk_fault failed, return code %d", rval);
            if (rval == -1) fprintf (stderr, " (memory)\n");
            else if (rval == -2) fprintf (stderr, " (determinant)\n");
            else if (rval == -3) fprintf (stderr, " (simplex iterations)\n");
            else fprintf (stderr, " (unknown)\n");
        }
        CCchunk_print_separate_timer (&separate_t);
        CCchunk_print_lift_timer (&lift_t);
    }

 CLEANUP:
    CCchunk_graph_free (c);
    CC_IFFREE (f.sols, int);
    CC_IFFREE (f.a.coef, int);
    return rval;
}

static int begin_CCchunk_cut_callback (void *u_data)
{
    CCchunk_cut_callback_data *ccd = (CCchunk_cut_callback_data *) u_data;

    ccd->cutval = 0.0;
    printf ("CUT FOUND:");
    return 0;
}

static int add_clique_callback (int *arr, int size, void *u_data)
{
    int i;
    CCchunk_cut_callback_data *ccd = (CCchunk_cut_callback_data *) u_data;
    CCchunk_graph *c = ccd->chunk;
    int *mark = (int *) NULL;
    double v;

    mark = CC_SAFE_MALLOC (c->ncount, int);
    if (mark == (int *) NULL) {
        fprintf (stderr, "Out of memory in add_clique_callback\n");
        return -1;
    }
    for (i=0; i<c->ncount; i++) mark[i] = 0;
    for (i=0; i<size; i++) {
        mark[arr[i]] = 1;
    }
    v = 0.0;
    for (i=0; i<c->ecount; i++) {
        if (mark[c->end0[i]] == 1 && mark[c->end1[i]] == 1) {
            v += c->weight[i];
        }
    }
    ccd->cutval += 2*(size - v);
    CC_FREE (mark, int);
    
    printf ("[");
    for (i=0; i<size; i++) {
        printf ("%d", arr[i]);
        if (i != size-1) printf (" ");
    }
    printf ("]");
    
    return 0;
}

static int abort_CCchunk_cut_callback (CC_UNUSED void *u_data)
{
    printf ("ABORT BUILDCUT\n");
    fflush (stdout);
    return 0;
}

static int finish_CCchunk_cut_callback (int rhs, int *finished, void *u_data)
{
    CCchunk_cut_callback_data *ccd = (CCchunk_cut_callback_data *) u_data;

    printf ("[ >= %d]  viol %.6f\n", rhs, rhs - ccd->cutval);
    ccd->cut_count++;
    *finished = !ccd->all_cuts;
    return 0;
}

static int found_CCchunk_fault_callback (CCchunk_graph *chunk, CCchunk_fault *fault,
        int *finished, void *u_data)
{
    CCchunk_fault_callback_data *fcd = (CCchunk_fault_callback_data *) u_data;
    CCchunk_cut_callback ccb;
    CCchunk_cut_callback_data ccd;
    int rval;
    int i;
    double s;

    ccd.all_cuts = 0;
    ccd.cut_count = 0;
    ccd.chunk = chunk;

    ccb.begin_cut = begin_CCchunk_cut_callback;
    ccb.add_clique = add_clique_callback;
    ccb.abort_cut = abort_CCchunk_cut_callback;
    ccb.finish_cut = finish_CCchunk_cut_callback;
    ccb.u_data = (void *) &ccd;

    s = fault->a.rhs;
    for (i=0; i<chunk->ecount; i++) {
        s -= fault->a.coef[i] * chunk->weight[i];
    }
    
    printf ("Found fault:");
    for (i=0; i<chunk->ecount; i++) {
        printf (" %d", fault->a.coef[i]);
    }
    printf (" <= %d\n", fault->a.rhs);
    printf ("viol %.6f, now lifting...\n", -s);
    fflush (stdout);

    rval = CCchunk_lift (chunk, fault, fcd->lift_t, &ccb);
    if (rval) return rval;

    *finished = (ccd.all_cuts == 0 && ccd.cut_count > 0);
    return 0;
}

