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
/*  int CCchunk_localcuts (CCtsp_lpcut_in **clist, int *cutcount,           */
/*      int ncount, int ecount, int *elist, double *x, double eps,          */
/*      CCchunk_flag flags, CCchunk_localcut_timer *timer, int silent,      */
/*      CCrandstate *rstate)                                                */
/*    MISSING                                                               */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "tsp.h"
#include "localcut.h"
#include "cut.h"

#undef  DUMPGRAPHS
#define OLD_VERSION

#define BLOCK_INTERIOR      /* block only the interior of a faulty chunk */
#define CHECK_INTERIOR      /* block chunk only if the interior is blocked */
#define OVERLAP_BOUND 2

typedef struct chunk_callback_data {
    CCtsp_cutinfo          *cuts;
    CCchunk_separate_timer *separate_t;
    CCchunk_lift_timer     *lift_t;
    int                    *blocked;
    int                     nolift;
} chunk_callback_data;

typedef struct fault_callback_data {
    CCtsp_cutinfo      *cuts;
    CCchunk_lift_timer *lift_t;
    int                 fault_count;
    int                *blocked;
    int                 nolift;
} fault_callback_data;

typedef struct cut_callback_data {
    CCtsp_cutinfo *cuts;
    int            all_cuts;
    int            cut_count;
} cut_callback_data;


static int
    found_chunk_callback (CCchunk_graph *chunk, int *faulty, void *u_data),
    found_fault_callback (CCchunk_graph *chunk, CCchunk_fault *fault,
        int *finished, void *u_data),
    found_fault_callback_nolift (CCchunk_graph *chunk, CCchunk_fault *fault,
        int *finished, void *u_data),
    begin_cut_callback (void *u_data),
    add_clique_callback (int *arr, int size, void *u_data),
    abort_cut_callback (void *u_data),
    finish_cut_callback (int rhs, int *finished, void *u_data),
    shrink_ones (int ncount, int ecount, int *elist, double *dlen,
        int *oncount, int *oecount, int **olist, double **olen,
        CC_SRKexpinfo *expand);

#ifdef DUMPGRAPHS
static void
    dumpgraph (char *fname, int ncount, int ecount, int *elist, double *x);
#endif



int CCchunk_localcuts (CCtsp_lpcut_in **clist, int *cutcount, int ncount,
        int ecount, int *elist, double *x, double eps, CCchunk_flag flags,
        CCchunk_localcut_timer *timer, int silent, CCrandstate *rstate)
{
    int rval;
    CCtsp_cutinfo cuts;
    chunk_callback_data ccd;
    CCchunk_chunk_callback ccb;
    int sncount;
    int secount;
    int *selist = (int *) NULL;
    double *sx = (double *) NULL;
    double st;
#ifndef OLD_VERSION
    int i;
    CCchunk_flag newflags;
#endif

    CCutil_start_timer (&timer->all);
    ccd.blocked = (int *) NULL;

    if (flags.nolift) {
        printf ("WARNING: nolift %d is turned on.  Will not lift faults\n",
                flags.nolift);
        fflush (stdout);
    }

#ifdef DUMPGRAPHS
    {
        static int cnt = 0;
        char nam[1024];

        sprintf (nam, "local.x.%c%d.%d",flags.spheres ? 's' : 'c',
                 flags.maxchunksize, cnt);
        dumpgraph (nam, ncount, ecount, elist, x);
        cnt++;
    }
#endif

    CCcut_SRK_init_expinfo (&cuts.expand);

    *cutcount = 0;
    cuts.clist = clist;
    cuts.cutcount = cutcount;

    ccd.cuts = &cuts;
    ccd.separate_t = &timer->separate;
    ccd.lift_t = &timer->lift;
    ccd.nolift = flags.nolift;

    ccb.func = found_chunk_callback;
    ccb.u_data = (void *) &ccd;
    
    if (!silent) {
        printf ("localcuts");
        if (flags.noshrink) printf (" noshrink");
        if (flags.uncivilized) printf (" uncivilized");
#ifdef OLD_VERSION
        if (flags.dummy) {
            printf (" dummy(%ud/%ud)", flags.spheresize, flags.maxchunksize);
        } else if (flags.permute) {
            printf (" permuted spheres(%ud/%ud)", flags.spheresize,
                    flags.maxchunksize);
        } else if (flags.weighted) {
            printf (" weighted spheres(%ud/%ud)", flags.spheresize,
                    flags.maxchunksize);
        } else if (flags.spheres) {
            printf (" spheres(%ud/%ud)", flags.spheresize, flags.maxchunksize);
        } else {
            printf (" classes(%ud)",flags.maxchunksize);
        }
#else /* OLD_VERSION */
        printf ("(%ud)",flags.maxchunksize);
#endif /* OLD_VERSION */
        printf ("...");
        fflush (stdout);
    }

    if (flags.noshrink) {
        CCutil_start_timer (&timer->find.shrink);
        rval = CCcut_SRK_trivial (ncount, &cuts.expand);
        CCutil_stop_timer (&timer->find.shrink, 0);
        if (rval) {
            fprintf (stderr, "CCcut_SRK_trivial failed\n");
            goto CLEANUP;
        }
        sncount = ncount;
        secount = ecount;
        selist = elist;
        sx = x;
    } else {
        CCutil_start_timer (&timer->find.shrink);
        rval = shrink_ones (ncount, ecount, elist, x, &sncount, &secount,
                            &selist, &sx, &cuts.expand);
        CCutil_stop_timer (&timer->find.shrink, 0);
        if (rval) {
            fprintf (stderr, "shrink_ones failed\n");
            goto CLEANUP;
        }
    }        

    if (!silent) {
        printf (" %d nodes %d edges\n", sncount, secount);
        fflush (stdout);
    }

#ifdef OLD_VERSION
    rval = CCchunk_finder (sncount, secount, selist, sx, eps, flags,
                           &timer->find, &ccb, rstate);
    if (rval) {
        fprintf (stderr, "CCchunk_finder failed\n");
        goto CLEANUP;
    }
#else  /* OLD_VERSION */
    ccd.blocked = CC_SAFE_MALLOC (sncount, int);
    if (ccd.blocked == (int *) NULL) {
        fprintf (stderr, "Out of memory in CCchunk_localcuts\n");
        rval = 1; goto CLEANUP;
    }
    for (i=0; i<sncount; i++) {
        ccd.blocked[i] = 0;
    }

    newflags = flags;
    for (i=8; i<=(int) flags.maxchunksize; i++) {
        newflags.dummy = 0;
        newflags.permute = 0;
        newflags.weighted = 0;
        newflags.spheres = 0;
        newflags.maxchunksize = i;
        newflags.spheresize = i-2;

        if (!silent) {
            printf ("%1.1d", i%10);
            fflush (stdout);
        }

#if 0
        rval = CCchunk_finder (sncount, secount, selist, sx, eps, newflags,
                               &timer->find, &ccb, rstate);
        if (rval) {
            fprintf (stderr, "CCchunk_finder failed\n");
            goto CLEANUP;
        }

#endif
        newflags.spheres = 1;
        
        rval = CCchunk_finder (sncount, secount, selist, sx, eps, newflags,
                               &timer->find, &ccb, rstate);
        if (rval) {
            fprintf (stderr, "CCchunk_finder failed\n");
            goto CLEANUP;
        }

        newflags.dummy = 1;

        rval = CCchunk_finder (sncount, secount, selist, sx, eps, newflags,
                               &timer->find, &ccb, rstate);
        if (rval) {
            fprintf (stderr, "CCchunk_finder failed\n");
            goto CLEANUP;
        }

        newflags.dummy = 0;
        newflags.permute = 1;

        rval = CCchunk_finder (sncount, secount, selist, sx, eps, newflags,
                               &timer->find, &ccb, rstate);
        if (rval) {
            fprintf (stderr, "CCchunk_finder failed\n");
            goto CLEANUP;
        }

        newflags.permute = 0;
        newflags.weighted = 1;
        
        rval = CCchunk_finder (sncount, secount, selist, sx, eps, newflags,
                               &timer->find, &ccb, rstate);
        if (rval) {
            fprintf (stderr, "CCchunk_finder failed\n");
            goto CLEANUP;
        }
    }

    if (!silent)  printf ("\n");
#endif /* OLD_VERSION */
    st =  CCutil_stop_timer (&timer->all, 0);
    if (!silent) {
        printf ("localcuts done in %.2f seconds\n", st);
        fflush (stdout);
    }

    rval = 0;

  CLEANUP:
    CCcut_SRK_free_expinfo (&cuts.expand);
    CC_IFFREE (ccd.blocked, int);
    if (selist != elist) CC_IFFREE (selist, int);
    if (sx != x) CC_IFFREE (sx, double);
    return rval;
}

static int found_chunk_callback (CCchunk_graph *chunk, int *faulty, void *u_data)
{
    chunk_callback_data *ccd = (chunk_callback_data *) u_data;
    fault_callback_data fcd;
    CCchunk_fault_callback fcb;
    int rval;
#ifndef OLD_VERSION
    int i;
    int *members;
    int *blocked = ccd->blocked;
    int j;
#endif

#ifndef OLD_VERSION
    for (i=0; i<chunk->ncount; i++) {
#ifdef CHECK_INTERIOR
        if (chunk->equality[i]) {
#endif
            members = chunk->members[i];
            for (j=0; members[j] != -1; j++) {
                if (blocked[members[j]] >= OVERLAP_BOUND) return 0;
            }
#ifdef CHECK_INTERIOR
        }
#endif
    }
#endif /* OLD_VERSION */

    fcd.cuts = ccd->cuts;
    fcd.lift_t = ccd->lift_t;
    fcd.fault_count = 0;
#ifndef OLD_VERSION
    fcd.blocked = blocked;
#endif
    fcd.nolift = ccd->nolift;

    if (ccd->nolift) {
        fcb.func = found_fault_callback_nolift;
    } else {
        fcb.func = found_fault_callback;
    }
    fcb.u_data = (void *) &fcd;
    
    rval = CCchunk_separate (chunk, ccd->separate_t, &fcb);
    if (rval) return rval;

    *faulty = (fcd.fault_count > 0);

    return 0;
}

static int found_fault_callback (CCchunk_graph *chunk, CCchunk_fault *fault,
        int *finished, void *u_data)
{
    fault_callback_data *fcd = (fault_callback_data *) u_data;
    CCchunk_cut_callback ccb;
    cut_callback_data ccd;
    int rval;
#ifndef OLD_VERSION
    int i;
    int *members;
    int *blocked = fcd->blocked;
    int j;
#endif

#ifndef OLD_VERSION
    for (i=0; i<chunk->ncount; i++) {
#ifdef BLOCK_INTERIOR
        if (chunk->equality[i]) {
#endif
            members = chunk->members[i];
            for (j=0; members[j] != -1; j++) {
                blocked[members[j]]++;
            }
#ifdef BLOCK_INTERIOR
        }
#endif
    }
#endif

    ccd.cuts = fcd->cuts;
    ccd.all_cuts = 0;
    ccd.cut_count = 0;

    ccb.begin_cut = begin_cut_callback;
    ccb.add_clique = add_clique_callback;
    ccb.abort_cut = abort_cut_callback;
    ccb.finish_cut = finish_cut_callback;
    ccb.u_data = (void *) &ccd;
    
    rval = CCchunk_lift (chunk, fault, fcd->lift_t, &ccb);
    if (rval) return rval;

    *finished = (ccd.all_cuts == 0 && ccd.cut_count > 0);
    return 0;
}

static int found_fault_callback_nolift (CCchunk_graph *chunk,
        CCchunk_fault *fault, int *finished, void *u_data)
{
    fault_callback_data *fcd = (fault_callback_data *) u_data;
    int i;
#ifndef OLD_VERSION
    int *members;
    int *blocked = fcd->blocked;
#endif
    double s;

#ifndef OLD_VERSION
    for (i=0; i<chunk->ncount; i++) {
#ifdef BLOCK_INTERIOR
        if (chunk->equality[i]) {
#endif
            int j;
            members = chunk->members[i];
            for (j=0; members[j] != -1; j++) {
                blocked[members[j]]++;
            }
#ifdef BLOCK_INTERIOR
        }
#endif
    }
#endif

    *finished = 1;
    CCutil_start_timer (&fcd->lift_t->all);
    CCutil_stop_timer (&fcd->lift_t->all, 0);
    if (fcd->nolift >= 2) {
        printf ("faulty chunk\n");
        printf ("%d %d\n", chunk->ncount, chunk->ecount);
        for (i=0; i<chunk->ecount; i++) {
            printf ("%d %d %.16f %d\n",chunk->end0[i], chunk->end1[i], chunk->weight[i],
                    chunk->fixed[i]);
        }
        for (i=0; i<chunk->ncount; i++) {
            if (i>0) printf (" ");
            printf ("%d",chunk->equality[i]);
        }
        printf ("\n");
        for (i=0; i<chunk->ecount; i++) {
            printf ("%d ", fault->a.coef[i]);
        }
        printf ("<= %d\n", fault->a.rhs);
        if (fcd->nolift >= 3) {
            printf ("%d\n", fault->nsols);
            for (i=0; i<fault->nsols; i++) {
                int j;
                for (j=0; j<chunk->ecount; j++) {
                    printf ("%d ", fault->sols[i*chunk->ecount+j]);
                }
                printf ("\n");
            }
        } else {
            printf ("0\n");
        }
        s = fault->a.rhs;
        for (i=0; i<chunk->ecount; i++) {
            s -= fault->a.coef[i] * chunk->weight[i];
        }
        printf ("viol %.6f\n", -s);
        fflush (stdout);
    }
    return 0;
}

static int begin_cut_callback (void *u_data)
{
    cut_callback_data *ccd = (cut_callback_data *) u_data;

    return CCtsp_buildcut_begin (ccd->cuts, 5);
}

static int add_clique_callback (int *arr, int size, void *u_data)
{
    cut_callback_data *ccd = (cut_callback_data *) u_data;

    return CCtsp_buildcut_addclique (ccd->cuts, arr, size);
}

static int abort_cut_callback (void *u_data)
{
    cut_callback_data *ccd = (cut_callback_data *) u_data;

    CCtsp_buildcut_abort (ccd->cuts);
    return 0;
}

static int finish_cut_callback (int rhs, int *finished, void *u_data)
{
    cut_callback_data *ccd = (cut_callback_data *) u_data;
    int rval;

    rval = CCtsp_buildcut_finish (ccd->cuts, rhs);
    if (rval) return rval;

    ccd->cut_count++;
    *finished = (ccd->all_cuts == 0 && ccd->cut_count > 0);

    return 0;
}

static int shrink_ones (int ncount, int ecount, int *elist, double *dlen,
        int *oncount, int *oecount, int **olist, double **olen,
        CC_SRKexpinfo *expand)
{
    CC_SRKgraph G;
    int k;

    if (CCcut_SRK_buildgraph (&G, ncount, ecount, elist, dlen)) {
        fprintf (stderr, "buildgraph failed in shrink_ones\n");
        return 1;
    }

    if (CCcut_SRK_defluff (&G)) {
        fprintf (stderr, "CCcut_SRK_defluff failed in shrink_ones\n");
        return 1;
    }

    CCcut_SRK_identify_paths_to_edges (&G, &k, 0);
    CCcut_SRK_identify_one_triangles (&G, &k, (CC_SRKnode *) NULL, 0.001, 2.0,
                                      0);

    if (CCcut_SRK_grab_edges (&G, oncount, oecount, olist, olen, expand)) {
        fprintf (stderr, "grab edges failed in shrink_ones\n");
        CCcut_SRK_free_graph (&G);
        return 1;
    }

    CCcut_SRK_free_graph (&G);
    return 0;
}

#ifdef DUMPGRAPHS
static void dumpgraph (char *fname, int ncount, int ecount, int *elist,
                       double *x)
{
    int i;
    FILE *f = fopen (fname, "w");

    if (!f) {
        perror (fname);
        fprintf (stderr, "Unable to open %s for output\n", fname);
        return;
    }

    fprintf (f, "%d %d\n", ncount, ecount);
    for (i=0; i<ecount; i++) {
        fprintf (f, "%d %d %.16f\n", elist[2*i], elist[2*i+1], x[i]);
    }
    fclose (f);
}
#endif
