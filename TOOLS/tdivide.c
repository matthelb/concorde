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
/*     A PROGRAM TO COMPUTE SUBDIVIDE A TSP INSTANCE FOR PARALLEL CODE      */
/*                                                                          */
/*                                TSP CODE                                  */
/*                                                                          */
/*  Written by:  Cook                                                       */
/*  Date: March 12, 2002                                                    */
/*                                                                          */
/*  Splits the instance into subproblems for lb, lkh, or linkern.  For      */
/*  lb, codes uses m depot nodes (m is set to 4*sqrt(nodes) by default).    */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "kdtree.h"

static char *tourfname = (char *) NULL;
static char *tspfname  = (char *) NULL;
static char *probname  = (char *) NULL;
static int simpletour = 0;
static int depotcount = 0;
static int binary_in = 0;
static int tsplib_in = 1;
static int seed = 0;
static int bound_part = 0;
static int lkh_part = 0;
static int tour_part = 0;
static int first_bucket = 0;
static int norm = CC_EUCLIDEAN;

int
    main (int ac, char **av);

static int
    create_tour_problems (CCdatagroup *dat, char *pname, int tcount,
        CCsubdiv *trac, int **slist),
    build_tour_subproblem (CCdatagroup *dat, char *pname, int id,
        int scount, int *slist),
    create_subproblems (int ncount, CCdatagroup *dat, int *invtour,
        char *pname, int ndepot, int tcount, CCsubdiv *trac, int **slist,
        CCrandstate *rstate, int norm),
    build_subproblem (CCdatagroup *dat, CCkdtree *kt, int *invtour,
        char *pname, int ndepot, int id, int scount, int *slist, char *hit,
        int norm, CCsubdiv *sbox),
    lkh_subproblem (int id, CCdatagroup *dat, int scount, int *sub, int start,
        CCsubdiv_lkh *plist, char *pname),
    geom_len (double x0, double y0, double x1, double y1, int *len),
    create_lkh (char *name, int ncount, CCdatagroup *dat, int *tour,
        int bucketsize, int first, double tourlen),
    find_full_tour (char *name, int scount, CCsubdiv *trac),
    find_part_tour (char *name, int scount, CCsubdiv *trac),
    parseargs (int ac, char **av);

static char
    *get_problabel (const char *probloc);

static void
    no_overlap (int tcount, CCsubdiv *trac, int *yesno),
    in_box (double x, double y, CCsubdiv *b, int *yesno, int i, int j,
        CCsubdiv *trac),
    usage (char *fname);


int main (int ac, char **av)
{
    int i, lap, ncount, pcount, rval = 0;
    int *tour = (int *) NULL;
    int *invtour = (int *) NULL;
    int **plist = (int **) NULL;
    double val;
    CCdatagroup dat;
    char *name = (char *) NULL;
    CCsubdiv *trac = (CCsubdiv *) NULL;
    CCrandstate rstate;

    CCutil_init_datagroup (&dat);
    seed = (int) CCutil_real_zeit ();

    rval = parseargs (ac, av);
    if (rval) return 1;

    if (!tspfname) {
        fprintf (stderr, "No TSPLIB or DAT file specified\n");
        goto CLEANUP;
    }


    CCutil_sprand (seed, &rstate);

    if (probname) name = get_problabel (probname);
    else          name = get_problabel (tspfname);

    printf ("Name: %s\n", name); fflush (stdout);

    if (tsplib_in) {
        rval = CCutil_gettsplib (tspfname, &ncount, &dat);
        CCcheck_rval (rval, "CCutil_gettsplib failed");
        CCutil_dat_getnorm (&dat, &norm);
    } else {
        rval = CCutil_getdata (tspfname, binary_in, norm, &ncount, &dat,
                               0, 0, &rstate);
        CCcheck_rval (rval, "CCutil_getdata failed");
    }

    if ((norm & CC_NORM_SIZE_BITS) != CC_D2_NORM_SIZE) {
        fprintf (stderr, "Only set up for 2D norms\n");
        rval = 1;  goto CLEANUP;
    }

    if (tourfname) {
        tour = CC_SAFE_MALLOC (ncount, int);
        CCcheck_NULL (tour, "out of memory in main");
        invtour = CC_SAFE_MALLOC (ncount, int);
        CCcheck_NULL (invtour, "out of memory in main");

        if (simpletour) {
            rval = CCutil_getcycle (ncount, tourfname, tour, 0);
            CCcheck_rval (rval, "CCutil_getcycle failed");
        } else {
            rval = CCutil_getcycle_tsplib (ncount, tourfname, tour);
            CCcheck_rval (rval, "CCutil_getcycle_tsplib failed");
        }

        CCutil_cycle_len (ncount, &dat, tour, &val);
        printf ("Tour Length: %.0f\n", val); fflush (stdout);

        for (i = 0; i < ncount; i++) {
           invtour[tour[i]] = i;
        }
    }

    if (lkh_part > 0) {
        if (!tour) {
            fprintf (stderr, "Need to specify a tour for LKH subproblems\n");
            rval = 1;  goto CLEANUP; 
        }
        rval = create_lkh (name, ncount, &dat, tour, lkh_part, first_bucket,
                           val);
        CCcheck_rval (rval, "create_lkh failed");
    } else if (bound_part > 0) {
        if ((norm & CC_NORM_BITS) != CC_KD_NORM_TYPE &&
                             norm != CC_GEOM) {
            fprintf (stderr, "Only set up for GEOM and KD-tree norms\n");
            rval = 1;  goto CLEANUP;
        }
        if (!tour) {
            fprintf (stderr, "Need to specify a tour for bound subproblems\n");
            rval = 1;  goto CLEANUP; 
        }
        rval = CCutil_karp_partition (ncount, &dat, bound_part, &pcount,
                                      &trac, &plist, &rstate);
        CCcheck_rval (rval, "CCutil_karp_partition failed");

        no_overlap (pcount, trac, &lap);
        if (lap == 1) {
            fprintf (stderr, "Error: regions overlap\n");
            rval = 1;  goto CLEANUP;
        }

        rval = create_subproblems (ncount, &dat, invtour, name, depotcount,
                                   pcount, trac, plist, &rstate, norm);
        CCcheck_rval (rval, "create_subproblems failed");

        rval = CCutil_write_subdivision_index (name, ncount, pcount, trac);
        CCcheck_rval (rval, "CCutil_write_subdivision_index failed");
    } else if (tour_part > 0)  {
        rval = CCutil_karp_partition (ncount, &dat, tour_part, &pcount,
                                      &trac, &plist, &rstate);
        CCcheck_rval (rval, "CCutil_karp_partition failed");

        rval = create_tour_problems (&dat, name, pcount, trac, plist);
        CCcheck_rval (rval, "create_tour_problems");

        rval = CCutil_write_subdivision_index (name, ncount, pcount, trac);
        CCcheck_rval (rval, "CCutil_write_subdivision_index failed");
    } else {
        fprintf (stderr, "No partition type specified\n");
        goto CLEANUP;
    }

CLEANUP:

    CC_IFFREE (tour, int);
    CC_IFFREE (invtour, int);
    CC_IFFREE (name, char);
    CC_IFFREE (trac, CCsubdiv);
    CCutil_freedatagroup (&dat);
    if (plist) {
        for (i = 0; i < pcount; i++) {
            CC_IFFREE (plist[i], int);
        }
        CC_FREE (plist, int *);
    }

    return rval;
}

static int create_tour_problems (CCdatagroup *dat, char *pname,
       int tcount, CCsubdiv *trac, int **slist)
{
    int i, rval = 0;

    for (i = 0; i < tcount; i++) {
        rval = build_tour_subproblem (dat, pname, trac[i].id, trac[i].cnt,
                                      slist[i]); 
        CCcheck_rval (rval, "build_tour_subproblem failed");
    }

CLEANUP:

    return rval;
}

static int build_tour_subproblem (CCdatagroup *dat, char *pname, int id,
        int scount, int *slist)
{
    int i, inorm, rval = 0;
    char buf[1024];
    FILE *out = (FILE *) NULL;
    CCdatagroup sdat;

    CCutil_init_datagroup (&sdat);
    CCutil_dat_getnorm (dat, &inorm);
    CCutil_dat_setnorm (&sdat, inorm);

    sdat.x = CC_SAFE_MALLOC (scount, double);
    CCcheck_NULL (sdat.x, "out of memory in build_tour_subproblem");
    sdat.y = CC_SAFE_MALLOC (scount, double);
    CCcheck_NULL (sdat.y, "out of memory in build_tour_subproblem");

    for (i = 0; i < scount; i++) {
        sdat.x[i] = dat->x[slist[i]];
        sdat.y[i] = dat->y[slist[i]];
    }

    sprintf (buf, "%s_%d.dat", pname, id);
    printf ("Create %s\n", buf); fflush (stdout);
    rval = CCutil_writedata (buf, 0, scount, &sdat);
    CCcheck_rval (rval, "CCutil_writedata failed");

    sprintf (buf, "%s_%d.nam", pname, id);
    printf ("Create %s\n", buf); fflush (stdout);

    out = fopen (buf, "w");
    if (out == (FILE *) NULL) {
        fprintf (stderr, "Could not open %s for output\n", buf);
        rval = 1;  goto CLEANUP;
    }
    fprintf (out, "%d\n", scount);
    for (i = 0; i < scount; i++) {
        fprintf (out, "%d\n", slist[i]);
    }

CLEANUP:

    if (out) fclose (out);
    CCutil_freedatagroup (&sdat);

    return rval;
}


static void no_overlap (int tcount, CCsubdiv *trac, int *yesno)
{
    int i, j;
    CCsubdiv *p, *q;
    double x0, x1, y0, y1;

    *yesno = 0;

    for (i = 0; i < tcount; i++) {
        for (j = i+1; j < tcount; j++) {
            p = &trac[i];
            q = &trac[j];

            x0 = p->xrange[0]; x1 = p->xrange[1];
            y0 = p->yrange[0]; y1 = p->yrange[1];
            in_box (x0, y0, q, yesno, i, j, trac); if (*yesno) return;
            in_box (x0, y1, q, yesno, i, j, trac); if (*yesno) return;
            in_box (x1, y0, q, yesno, i, j, trac); if (*yesno) return;
            in_box (x1, y1, q, yesno, i, j, trac); if (*yesno) return;

            x0 = q->xrange[0]; x1 = q->xrange[1];
            y0 = q->yrange[0]; y1 = q->yrange[1];
            in_box (x0, y0, p, yesno, i, j, trac); if (*yesno) return;
            in_box (x0, y1, p, yesno, i, j, trac); if (*yesno) return;
            in_box (x1, y0, p, yesno, i, j, trac); if (*yesno) return;
            in_box (x1, y1, p, yesno, i, j, trac); if (*yesno) return;
        }
    }
    printf ("Boxes do not overlap -- good\n"); fflush (stdout);
}

static void in_box (double x, double y, CCsubdiv *b, int *yesno, int i, int j,
       CCsubdiv *trac)
{
    if (x > b->xrange[0] && x < b->xrange[1] && 
        y > b->yrange[0] && y < b->yrange[1]) {
        fprintf (stderr, "Box %d (%f, %f, %f, %f) intersects\n",
            i, trac[i].xrange[0], trac[i].xrange[1], trac[i].yrange[0],
            trac[i].yrange[1]);
        fprintf (stderr, "Box %d (%f, %f, %f, %f)\n",
            j, trac[j].xrange[0], trac[j].xrange[1], trac[j].yrange[0],
            trac[j].yrange[1]);
        *yesno = 1;
    } else {
        *yesno = 0;
    }
}

static int create_subproblems (int ncount, CCdatagroup *dat, int *invtour,
        char *pname, int ndepot, int tcount, CCsubdiv *trac, int **slist,
        CCrandstate *rstate, int norm)
{
    CCkdtree kt;
    CCkdtree *p_kt = (CCkdtree *) NULL;
    int i, rval = 0;
    char *hit = (char *) NULL;

    if (norm != CC_GEOM) {
        rval = CCkdtree_build (&kt, ncount, dat, (double *) NULL, rstate);
        CCcheck_rval (rval, "CCkdtree_build failed");
        p_kt = &kt;
    }

    hit = CC_SAFE_MALLOC (ncount, char);
    CCcheck_NULL (hit, "out of memory in create_subproblems");
    for (i = 0; i < ncount; i++) hit[i] = 0;

    for (i = 0; i < tcount; i++) {
        rval = build_subproblem (dat, p_kt, invtour, pname, ndepot,
                           trac[i].id, trac[i].cnt, slist[i], hit, norm,
                           &trac[i]); 
        CCcheck_rval (rval, "build_subproblem failed");
    }

    for (i = 0; i < ncount; i++) {
        if (hit[i] == 0) {
            fprintf (stderr, "missed a node in partitions\n");
            rval = 1; goto CLEANUP;
        }
    }

CLEANUP:

    if (p_kt) CCkdtree_free (p_kt);
    CC_IFFREE (hit, char);
    return rval;
}

static int build_subproblem (CCdatagroup *dat, CCkdtree *kt, int *invtour,
        char *pname, int ndepot, int id, int scount, int *slist, char *hit,
        int norm, CCsubdiv *sbox)
{
    int rval = 0;
    int i, neighbor, len, t;
    int *tpos = (int *) NULL;
    int *perm = (int *) NULL;
    int *p_slist = (int *) NULL;
    CC_SFILE *out = (CC_SFILE *) NULL;
    char buf[CCutil_FILE_NAME_LEN];
    double xn, yn;

    if (ndepot == 0) {
        ndepot = (int) sqrt ((double) scount);
        ndepot = 4*(ndepot + 1);
    }


    for (i = 0; i < scount; i++) {
        xn = dat->x[slist[i]];
        yn = dat->y[slist[i]];

        if (xn < sbox->xrange[0] || xn > sbox->xrange[1] ||
            yn < sbox->yrange[0] || yn > sbox->yrange[1]) {
            fprintf (stderr, "Point (%f,%f) not in Box (%f, %f, %f, %f)\n",
                  xn, yn, sbox->xrange[0], sbox->xrange[1], sbox->yrange[0],
                  sbox->yrange[1]);
            rval = 1;  goto CLEANUP;

        }
    }

    sprintf (buf, "%s_%d.mas", pname, id);
    printf ("Create %s, with %d depots\n", buf, ndepot); fflush (stdout);

    for (i = 0; i < scount; i++) {
        if (hit[slist[i]]) {
            fprintf (stderr, "duplicate node in partitions");
            rval = 1;  goto CLEANUP;
        }
        hit[slist[i]] = 1;
    }

    /* Build the permutation from the exisiting tour */

    tpos = CC_SAFE_MALLOC (scount, int);
    CCcheck_NULL (tpos, "out of memory in build_subproblem");
    perm = CC_SAFE_MALLOC (scount, int);
    CCcheck_NULL (perm, "out of memory in build_subproblem");

    for (i = 0; i < scount; i++) {
        tpos[i] = invtour[slist[i]];
        perm[i] = i;
    }
    CCutil_int_perm_quicksort (perm, tpos, scount);
    CC_FREE (tpos, int);

    p_slist = CC_SAFE_MALLOC (scount, int);
    CCcheck_NULL (p_slist, "out of memory in build_subproblem");
    for (i = 0; i < scount; i++) {
        p_slist[i] = slist[perm[i]];
    }

    /* The header information */

    out = CCutil_sopen (buf, "w");
    if (out == (CC_SFILE *) NULL) {
        fprintf (stderr, "Could not open %s for output\n", buf);
        rval = 1;  goto CLEANUP;
    }

    rval = CCutil_swrite_int (out, scount);
    CCcheck_rval (rval, "CCutil_swrite_int failed");
    rval = CCutil_swrite_int (out, CC_MASTER_DAT);
    CCcheck_rval (rval, "CCutil_swrite_int failed");
    rval = CCutil_swrite_int (out, CC_SUBDIVISION);
    CCcheck_rval (rval, "CCutil_swrite_int failed");
    rval = CCutil_swrite_int (out, norm);
    CCcheck_rval (rval, "CCutil_swrite_int failed");
    rval = CCutil_swrite_int (out, ndepot);
    CCcheck_rval (rval, "CCutil_swrite_int failed");

    /* The depot costs */

    if (norm != CC_GEOM) {
        for (i = 0; i < scount; i++) {
            CCkdtree_delete (kt, slist[i]);
        }
    }

    for (i = 0; i < scount; i++) {
        if (norm == CC_GEOM) {
            xn = dat->x[p_slist[i]];
            yn = dat->y[p_slist[i]];

            len = CCutil_MAXINT;
            rval = geom_len (xn, yn, sbox->xrange[0], yn, &t);
            CCcheck_rval (rval, "geom_len failed");
            if (t < len) len = t;

            rval = geom_len (xn, yn, sbox->xrange[1], yn, &t);
            CCcheck_rval (rval, "geom_len failed");
            if (t < len) len = t;

            rval = geom_len (xn, yn, xn, sbox->yrange[0], &t);
            CCcheck_rval (rval, "geom_len failed");
            if (t < len) len = t;

            rval = geom_len (xn, yn, xn, sbox->yrange[1], &t);
            CCcheck_rval (rval, "geom_len failed");
            if (t < len) len = t;
        } else {
            CCkdtree_undelete (kt, p_slist[i]);
            neighbor = CCkdtree_node_nearest (kt, p_slist[i], dat,
                                             (double *) NULL);
            len = CCutil_dat_edgelen (p_slist[i], neighbor, dat);
            len = len / 2;    /* Since edge might get hit in two regions */
            CCkdtree_delete (kt, p_slist[i]);
        }
        rval = CCutil_swrite_int (out, len);
        CCcheck_rval (rval, "CCutil_swrite_int failed");
    }

    if (norm != CC_GEOM) {
        for (i = 0; i < scount; i++) {
            CCkdtree_undelete (kt, slist[i]);
        }
    }

    /* The permuted datagroup */

    for (i = 0; i < scount; i++) {
        rval = CCutil_swrite_double (out, dat->x[p_slist[i]]);
        CCcheck_rval (rval, "CCutil_swrite_double failed");
        rval = CCutil_swrite_double (out, dat->y[p_slist[i]]);
        CCcheck_rval (rval, "CCutil_swrite_double failed");
    }

    /* The permutation */

    for (i = 0; i < scount; i++) {
        rval = CCutil_swrite_int (out, perm[i]);
        CCcheck_rval (rval, "CCutil_swrite_int failed");
    }

    /* Original node names (not permuted) */

    for (i = 0; i < scount; i++) {
        rval = CCutil_swrite_int (out, slist[i]);
        CCcheck_rval (rval, "CCutil_swrite_int failed");
    }

CLEANUP:

    CC_IFFREE (tpos, int);
    CC_IFFREE (perm, int);
    CC_IFFREE (p_slist, int);
    if (out) CCutil_sclose (out);
    return rval;
}

static int geom_len (double x0, double y0, double x1, double y1, int *len)
{
    int rval = 0;
    CCdatagroup dat;

    CCutil_init_datagroup (&dat);
    CCutil_dat_setnorm (&dat, CC_GEOM);

    dat.x = CC_SAFE_MALLOC (2, double);
    CCcheck_NULL (dat.x, "out of memory in geom_len");
    dat.y = CC_SAFE_MALLOC (2, double);
    CCcheck_NULL (dat.y, "out of memory in geom_len");

    dat.x[0] = x0;
    dat.y[0] = y0;
    dat.x[1] = x1;
    dat.y[1] = y1;

    *len = CCutil_dat_edgelen (0, 1, &dat);

CLEANUP:

    CCutil_freedatagroup (&dat);
    return rval;
}

static int create_lkh (char *name, int ncount, CCdatagroup *dat, int *tour,
        int bucketsize, int first, double tourlen)
{
    int rval = 0;
    FILE *out = (FILE *) NULL;
    int *sub = (int *) NULL;
    int i, k, nsub, remain, extra, start;
    CCsubdiv_lkh *plist = (CCsubdiv_lkh *) NULL;

    if (first <= 0 || first >= ncount) first = bucketsize;

    nsub = 1 + ((ncount - first) / bucketsize);
    remain = ncount - first  - ((nsub-1) * bucketsize);
    extra = (remain > bucketsize/16 ? 1 : 0);

    if (first < bucketsize + remain) {
        sub = CC_SAFE_MALLOC (bucketsize + remain, int);
    } else {
        sub = CC_SAFE_MALLOC (first, int);
    }
    CCcheck_NULL (sub, "out of memory in create_lkh");

    plist = CC_SAFE_MALLOC (nsub + extra, CCsubdiv_lkh);
    CCcheck_NULL (plist, "out of memory in create_lkh");

    printf ("Create %d LKH subproblems\n", nsub + extra);
    fflush (stdout);

    start = 0;
    for (i = 0; i < first; i++) {
        sub[i] = tour[start+i];
    }
    rval = lkh_subproblem (0, dat, first, sub, start, plist, name);
    CCcheck_rval (rval, "lkh_subproblem failed");
    start += first;

    for (k = 1; k < nsub-1; k++) {
        for (i = 0; i < bucketsize; i++) {
            sub[i] = tour[start+i];
        }
        rval = lkh_subproblem (k, dat, bucketsize, sub, start, plist, name);
        CCcheck_rval (rval, "lkh_subproblem failed");
        start += bucketsize;
    }
    if (extra == 0) {
        for (i = 0; i < bucketsize+remain; i++) {
            sub[i] = tour[start+i];
        }
        rval = lkh_subproblem (k, dat, bucketsize+remain, sub, start, plist,
                               name);
        CCcheck_rval (rval, "lkh_subproblem failed");
        start += (bucketsize+remain);

    } else {
        for (i = 0; i < bucketsize; i++) {
            sub[i] = tour[start+i];
        }
        rval = lkh_subproblem (k, dat, bucketsize, sub, start, plist, name);
        CCcheck_rval (rval, "lkh_subproblem failed");
        start += bucketsize;
        k++;

        for (i = 0; i < remain; i++) {
            sub[i] = tour[start+i];
        }
        rval = lkh_subproblem (k, dat, remain, sub, start, plist, name);
        CCcheck_rval (rval, "lkh_subproblem failed");
    }

    printf ("Subproblems\n");
    for (i = 0; i < nsub + extra; i++) {
        printf ("%d %d %d %.0f %.0f\n", plist[i].id, plist[i].cnt,
                        plist[i].start, plist[i].origlen, plist[i].newlen);
        fflush (stdout);
    }

    rval = CCutil_write_subdivision_lkh_index (name, ncount, nsub+extra,
                                               plist, tourlen);
    CCcheck_rval (rval, "CCutil_write_subdivision_lkh_index failed");

CLEANUP:

    CC_IFFREE (sub, int);
    CC_IFFREE (plist, CCsubdiv_lkh);
    if (out) fclose (out);
    return rval;
}

static int lkh_subproblem (int id, CCdatagroup *dat, int scount, int *sub,
        int start, CCsubdiv_lkh *plist, char *pname)
{
    int i;
    int rval = 0;
    double total = 0;


    for (i = 1; i < scount; i++) {
        total += ((double) CCutil_dat_edgelen (sub[i-1], sub[i], dat));
    }

    plist[id].id = id; 
    plist[id].cnt = scount;
    plist[id].start = start;
    plist[id].origlen = total;
    plist[id].newlen = -1.0;

    rval =  build_tour_subproblem (dat, pname, id, scount, sub);
    CCcheck_rval (rval, "build_tour_subproblem failed");

CLEANUP:

    return rval;
}

static int parseargs (int ac, char **av)
{
    int c, inorm;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "bBf:k:j:l:m:o:StT:P:N:?", &boptind, &boptarg)) != EOF) { 
        switch (c) {
        case 'B':
            binary_in = 1;
            break;
        case 'b':
            binary_in = 2;
            break;
        case 'f':
            first_bucket = atoi(boptarg);
            break;
        case 'l':
            lkh_part = atoi(boptarg);
            break;
        case 'j':
            tour_part = atoi(boptarg);
            break;
        case 'k':
            bound_part = atoi(boptarg);
            break;
        case 'm':
            depotcount = atoi(boptarg);
            break;
        case 'P':
            probname = boptarg;
            break;
        case 't':
            simpletour = 1;
            break;
        case 'T':
            tourfname  = boptarg;
            break;
        case 'N':
            inorm = atoi (boptarg);
            switch (inorm) {
            case 0: norm = CC_MAXNORM; break;
            case 1: norm = CC_MANNORM; break;
            case 2: norm = CC_EUCLIDEAN; break;
            case 3: norm = CC_EUCLIDEAN_3D; break;
            case 17: norm = CC_GEOM; break;
            case 18: norm = CC_EUCLIDEAN_CEIL; break;
            default:
                usage (av[0]);
                return 1;
            }
            tsplib_in = 0;
            break;
        case CC_BIX_GETOPT_UNKNOWN:
        case '?':
        default:
            usage (av[0]);
            return 1;
        }
    }

    if (boptind < ac) {
        tspfname = av[boptind++];
    } else {
        fprintf (stderr, "Missing tspfile\n");
        usage (av[0]);
        return 1;
    }

    return 0;
}

static void usage (char *fname)
{
    fprintf (stderr, "Usage: %s [-flags below] tsp_or_dat file\n", fname);
    fprintf (stderr, "   -b    datfile in double binary format\n");
    fprintf (stderr, "   -B    datfile in double integer format\n");
    fprintf (stderr, "   -f #  size of first bucket for LKH subproblems\n");
    fprintf (stderr, "   -j #  create tour subproblems with bucketsize #\n");
    fprintf (stderr, "   -k #  create bound subproblems with bucketsize #\n");
    fprintf (stderr, "   -l #  create LKH subproblems with bucketsize #\n");
    fprintf (stderr, "   -m #  number of depots (default 4*sqrt(n)\n");
    fprintf (stderr, "   -T f  specify a tour (needed for LKH and bound)\n");
    fprintf (stderr, "   -t    tour file in concorde format (default TSPLIB)\n");
    fprintf (stderr, "   -P s  specify a name for the output master files\n");
    fprintf (stderr, "   -N #  norm (must specify if dat file is not a TSPLIB file)\n");
    fprintf (stderr, "         0=MAX, 1=L1, 2=L2, 17=GEOM, 18=JOHNSON\n");
}

static char *get_problabel (const char *probloc)
{
    const char *p;
    const char *problabel = probloc;
    char *probcopy = (char *) NULL;
    char *p2;

    p = CCutil_strrchr_c (problabel, ':');
    if (p != (const char *) NULL) problabel = p+1;
    p = CCutil_strrchr_c (problabel, '/');
    if (p != (const char *) NULL) problabel = p+1;
    probcopy = CCutil_strdup (problabel);
    if (probcopy == (char *) NULL) return (char *) NULL;
    p2 = CCutil_strchr (probcopy, '.');
    if (p2 != (char *) NULL) *p2 = '\0';
    return probcopy;
}

