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
#include "util.h"
#include "cut.h"


static int
    dump_segment (double cut_val, int cut_start, int cut_end,
                  void *u_data);


int main (int ac, char **av)
{
    FILE *f = (FILE *) NULL;
    int rval;
    int *elist = (int *) NULL;
    double *x = (double *) NULL;
    int *endmark = (int *) NULL;
    int *perm = (int *) NULL;
    int ecount;
    int ncount;
    int i;
    double maxval;
    double szeit;

    if (ac < 2) {
        fprintf (stderr, "Usage: %s edge_file\n", av[0]);
        rval = -1;
        goto CLEANUP;
    }
    if (ac > 2) {
        maxval = atof (av[2]);
    } else {
        maxval = 2.001;
    }

    f = fopen (av[1], "r");
    if (f == (FILE *) NULL) {
        perror (av[1]);
        fprintf (stderr, "Unable to open %s for input\n", av[1]);
        rval = -1;
        goto CLEANUP;
    }

    fscanf (f, "%d%d", &ncount, &ecount);
    elist = CC_SAFE_MALLOC (ecount*2, int);
    x = CC_SAFE_MALLOC (ecount, double);
    endmark = CC_SAFE_MALLOC (ncount, int);
    perm = CC_SAFE_MALLOC (ncount, int);
    if (elist == (int *) NULL ||
        x == (double *) NULL ||
        endmark == (int *) NULL ||
        perm == (int *) NULL) {
        fprintf (stderr, "Out of memory\n");
        rval = -1;
        goto CLEANUP;
    }

    for (i=0; i<ecount; i++) {
        fscanf (f, "%d%d%lf", &elist[2*i], &elist[2*i+1], &x[i]);
    }

    fclose (f);
    f = (FILE *) NULL;

    for (i=0; i<ncount; i++) {
        endmark[i] = CC_LINSUB_BOTH_END;
    }

    printf ("generating all CCtsp_segment cuts of weight <= %f\n", maxval);

    szeit = CCutil_zeit();

    rval = CCcut_linsub_allcuts (ncount, ecount, (int *) NULL, endmark,
            elist, x, maxval, (void *) NULL, dump_segment);
    if (rval) {
        fprintf (stderr, "CCcut_linsub_allcuts failed\n");
        goto CLEANUP;
    }

    printf ("done in %.2f seconds\n", CCutil_zeit() - szeit);
    fflush (stdout);

    printf ("generating all CCtsp_segment cuts (even to *3) of weight <= %f\n",
            maxval);

    for (i=0; i<ncount; i++) {
        endmark[i] = 0;
    }
    for (i=0; i<ncount; i+=2) {
        endmark[i] |= CC_LINSUB_LEFT_END;
    }
    for (i=0; i<ncount; i+=3) {
        endmark[i] |= CC_LINSUB_RIGHT_END;
    }

    szeit = CCutil_zeit();

    rval = CCcut_linsub_allcuts (ncount, ecount, (int *) NULL, endmark,
            elist, x, maxval, (void *) NULL, dump_segment);
    if (rval) {
        fprintf (stderr, "CCcut_linsub_allcuts failed\n");
        goto CLEANUP;
    }

    printf ("done in %.2f seconds\n", CCutil_zeit() - szeit);
    fflush (stdout);

    printf ("generating all CCtsp_segment cuts rotated by 5 (even to *3) of weight <= %f\n",
            maxval);

    for (i=0; i<ncount; i++) {
        endmark[i] = 0;
    }
    for (i=0; i<ncount; i+=2) {
        endmark[i] |= CC_LINSUB_LEFT_END;
    }
    for (i=0; i<ncount; i+=3) {
        endmark[i] |= CC_LINSUB_RIGHT_END;
    }

    for (i=0; i<ncount; i++) {
        perm[i] = (i+5) % ncount;
    }

    szeit = CCutil_zeit();

    rval = CCcut_linsub_allcuts (ncount, ecount, perm, endmark, elist, x,
                           maxval, (void *) NULL, dump_segment);
    if (rval) {
        fprintf (stderr, "CCcut_linsub_allcuts failed\n");
        goto CLEANUP;
    }

    printf ("done in %.2f seconds\n", CCutil_zeit() - szeit);
    fflush (stdout);

    rval = 0;

  CLEANUP:
    if (f != (FILE *) NULL) fclose (f);
    CC_IFFREE (elist, int);
    CC_IFFREE (x, double);
    return rval;
}

static int dump_segment (double cut_val, int cut_start, int cut_end,
                         CC_UNUSED void *u_data)
{
    printf ("%d %d %.6f\n", cut_start, cut_end, cut_val);
    return 0;
}
