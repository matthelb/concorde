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
/*       A PROGRAM TO PATCH TOGETHER A TOUR FOR A SUBDIVIDED PROBLEM        */
/*                                                                          */
/*                                TSP CODE                                  */
/*                                                                          */
/*  Written by:  Cook                                                       */
/*  Date: April 16, 2003                                                    */
/*                                                                          */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "edgegen.h"
#include "linkern.h"
#include "macrorus.h"

#define B_POINTS 4000
#define LK_POINTS 50000

static char *tourfname = (char *) NULL;
static char *fulltourfname = (char *) NULL;
static char *tspfname  = (char *) NULL;
static char *indexfname = (char *) NULL;
static int binary_in = 0;
static int tsplib_in = 1;
static int seed = 0;
static int dump_partition_tsp = 0;
static int border_crossings = 0;
static int border_optimization = 0;
static int cat_tour = 0;
static int innorm = CC_EUCLIDEAN;


int
    main (int ac, char **av);

static int
    border_opt (int ncount, CCdatagroup *dat, char *name, int pcount,
        CCsubdiv *trac, int *tour, CCrandstate *rstate),
    insert_crossings (int ncount, CCdatagroup *dat, char *name,
        int pcount, CCsubdiv *trac, int *ptour, int *ftour),
    improve_border (char *name, int ncount, int ip, int iq,
        CCdatagroup *dat, CCsubdiv *trac, int *tour, int dateline, int pole,
        double ylo, double yhi, CCrandstate *rstate),
    insert_closest_pair (char *name, int ncount, int ip, int iq,
        CCdatagroup *dat, CCsubdiv *trac, int *tour, int dateline, int pole,
        double ylo, double yhi, double *tdelta, int *hit),
    cat_full_tour (int ncount, CCdatagroup *dat, char *name, int pcount,
        CCsubdiv_lkh *trac, double origval),
    find_full_tour (int ncount, CCdatagroup *dat, char *name,
        int pcount, CCsubdiv *trac, int *ptour),
    find_closest_pair (char *name, int ip, int iq, int *p_p0, int *p_p1,
        int *p_q0, int *p_q1, CCdatagroup *dat, CCsubdiv *trac, int *p_bcount,
        int **p_btour),
    grab_border_points (int btype, int scount, int *slist0,
        CCdatagroup *dat, int *p_bcount, int **p_blist, double xhi,
        double xlow, int maxpoints),
    get_partial_tour_sum (char *name, int pcount, double *sumval,
        CCdatagroup *dat),
    get_part_data (char *name, int ind, int *p_count, int **p_namelist,
        int **p_tour, char *tname, int use_default),
    run_lk_subproblem (int scount, int *slist, int ncount, int *tour,
        CCdatagroup *dat, CCrandstate *rstate),
    find_part_tour (char *name, int scount, CCsubdiv *trac),
    parseargs (int ac, char **av);

static char
    *get_problabel (const char *probloc);

static void
    subdiv_adj (CCsubdiv *s, CCsubdiv *t, int *len),
    dateline_adj (CCsubdiv *s, CCsubdiv *t, double ylo, double yhi,
        int *len),
    pole_adj (CCsubdiv *s, CCsubdiv *t, double xlo, double xhi,
        int *pole),
    find_common_border (CCsubdiv *p, CCsubdiv *q, int *border,
        int *flip),
    usage (char *fname);

/*
int hackarray[100];
*/


int main (int ac, char **av)
{
    int ncount, pcount, rval = 0;
    int *ptour = (int *) NULL;
    int *ftour = (int *) NULL;
    CCdatagroup dat;
    char *name = (char *) NULL;
    CCsubdiv *trac = (CCsubdiv *) NULL;
    CCsubdiv_lkh *lkhtrac = (CCsubdiv_lkh *) NULL;
    CCrandstate rstate;
    double origlen;

/*
    for (i = 0; i < 100; i++) hackarray[i] = 0;
*/

    CCutil_init_datagroup (&dat);
    if (!seed) seed = (int) CCutil_real_zeit ();

    rval = parseargs (ac, av);
    if (rval) return 1;

    if (!dump_partition_tsp && !tspfname) {
        fprintf (stderr, "No TSPLIB or DAT file specified\n");
        usage (av[0]);
        goto CLEANUP;
    }

    if ((!dump_partition_tsp && !cat_tour && !border_optimization) &&
         !tourfname) {
        fprintf (stderr, "No path file specfied\n");
        usage (av[0]);
        goto CLEANUP;
    }

    if ((border_crossings || border_optimization) && !fulltourfname) {
        fprintf (stderr, "No full tour specfied\n");
        usage (av[0]);
        goto CLEANUP;
    }

    if (!indexfname) {
        fprintf (stderr, "No index file specfied\n");
        usage (av[0]);
        goto CLEANUP;
    }

    CCutil_sprand (seed, &rstate);

    if (cat_tour) {
        rval = CCutil_read_subdivision_lkh_index (indexfname, &name, &ncount,
                                                  &pcount, &lkhtrac, &origlen);
        CCcheck_rval (rval, "CCutil_read_subdivision_lkh_index failed");
    } else {
        rval = CCutil_read_subdivision_index (indexfname, &name, &ncount,
                                              &pcount, &trac);
        CCcheck_rval (rval, "CCutil_read_subdivision_index failed");
    }
    CCcheck_rval (rval, "CCutil_read_subdivision_index failed");
    printf ("Name: %s\n", name);
    printf ("ncount = %d, partitions = %d\n", ncount, pcount);
    fflush (stdout);

    if (dump_partition_tsp) {
        rval = find_part_tour (name, pcount, trac);
        CCcheck_rval (rval, "find_part_tour failed");
        goto CLEANUP;
    }

    if (tsplib_in) {
        rval = CCutil_gettsplib (tspfname, &ncount, &dat);
        CCcheck_rval (rval, "CCutil_gettsplib failed");
        CCutil_dat_getnorm (&dat, &innorm);
    } else {
        rval = CCutil_getdata (tspfname, binary_in, innorm, &ncount, &dat,
                               0, 0, &rstate);
        CCcheck_rval (rval, "CCutil_getdata failed");
    }

    if ((innorm & CC_NORM_SIZE_BITS) != CC_D2_NORM_SIZE) {
        fprintf (stderr, "Only set up for 2D norms\n");
        rval = 1;  goto CLEANUP;
    }

    if (tourfname) {
        ptour = CC_SAFE_MALLOC (pcount, int);
        CCcheck_NULL (ptour, "out of memory in main");

        rval = CCutil_getcycle (pcount, tourfname, ptour, 0);
        CCcheck_rval (rval, "CCutil_getcycle failed");
    }

    if (cat_tour) {
        rval = cat_full_tour (ncount, &dat, name, pcount, lkhtrac, origlen);
        CCcheck_rval (rval, "cat_full_tour failed");
    } else if (border_optimization) {
        ftour = CC_SAFE_MALLOC (ncount, int);
        CCcheck_NULL (ftour, "out of memory in main");

        rval = CCutil_getcycle (ncount, fulltourfname, ftour, 0);
        CCcheck_rval (rval, "CCutil_getcycle failed");

        rval = border_opt (ncount, &dat, name, pcount, trac, ftour, &rstate);
        CCcheck_rval (rval, "insert_crossings failed");
    } else if (border_crossings) {
        ftour = CC_SAFE_MALLOC (ncount, int);
        CCcheck_NULL (ftour, "out of memory in main");

        rval = CCutil_getcycle (ncount, fulltourfname, ftour, 0);
        CCcheck_rval (rval, "CCutil_getcycle failed");

        rval = insert_crossings (ncount, &dat, name, pcount, trac, ptour,
                                 ftour);
        CCcheck_rval (rval, "insert_crossings failed");
    } else {
        rval = find_full_tour (ncount, &dat, name, pcount, trac, ptour);
        CCcheck_rval (rval, "find_find_tour failed");
    }


CLEANUP:

    CC_IFFREE (ptour, int);
    CC_IFFREE (ftour, int);
    CC_IFFREE (name, char);
    CC_IFFREE (trac, CCsubdiv);
    CC_IFFREE (lkhtrac, CCsubdiv_lkh);
    CCutil_freedatagroup (&dat);

    return rval;
}

static int border_opt (int ncount, CCdatagroup *dat, char *name, int pcount,
        CCsubdiv *trac, int *tour, CCrandstate *rstate)
{
    int rval = 0;
    double yhi, ylo, xhi, xlo, val, newval;
    int i, j, norm, len, pole, cnt = 0;
    char buf[1024];

    printf ("Border optimization ...\n"); fflush (stdout);
    CCutil_cycle_len (ncount, dat, tour, &val);
    printf ("Starting tour length:       %16.0f\n", val);  fflush (stdout);

    CCutil_dat_getnorm (dat, &norm);

    xlo = trac[0].xrange[0];
    xhi = trac[0].xrange[1];
    ylo = trac[0].yrange[0];
    yhi = trac[0].yrange[1];
    for (i = 1; i < pcount; i++) {
        if (trac[i].xrange[0] < xlo) xlo = trac[i].xrange[0];
        if (trac[i].xrange[1] > xhi) xhi = trac[i].xrange[1];
        if (trac[i].yrange[0] < ylo) ylo = trac[i].yrange[0];
        if (trac[i].yrange[1] > yhi) yhi = trac[i].yrange[1];
    }

    for (i = 0; i < pcount; i++) {
        for (j = i+1; j < pcount; j++) {
            subdiv_adj (&trac[i], &trac[j], &len);
            if (len == 0) {
                rval = improve_border (name, ncount, i, j, dat, trac,
                                      tour, 0, 0, ylo, yhi, rstate);
                CCcheck_rval (rval, "improve_border failed");
                cnt++;
            } else if (norm == CC_GEOM) {
                dateline_adj (&trac[i], &trac[j], ylo, yhi, &len);
                if (len == 0) {
                    rval = improve_border (name, ncount, i, j, dat, trac,
                                     tour, 1, 0, ylo, yhi, rstate);
                    CCcheck_rval (rval, "improve_border failed");
                    cnt++;
                } else {
                    pole_adj (&trac[i], &trac[j], xlo, xhi, &pole);
                    if (pole) {
                        rval = improve_border (name, ncount, i, j, dat, trac,
                                    tour, 0, pole, ylo, yhi, rstate);
                           CCcheck_rval (rval, "insert_closest_pair failed");
                        cnt++;
                    }
                }
            }
        }
    }

    printf ("Re-optimized %d borders\n", cnt); fflush (stdout);

    CCutil_cycle_len (ncount, dat, tour, &newval);
    printf ("Starting tour length:       %16.0f\n", val); 
    printf ("Adjusted tour length:       %16.0f\n", newval);
    fflush (stdout);

    sprintf (buf, "%s_window.tour", name);
    rval = CCutil_writecycle (ncount, buf, tour, 0);
    CCcheck_rval (rval, "CCutil_writecycle failed");

CLEANUP:

    return rval;
}

static int insert_crossings (int ncount, CCdatagroup *dat, char *name,
        int pcount, CCsubdiv *trac, int *ptour, int *ftour)
{
    int rval = 0;
    int i, j, cnt = 0, hit;
    int norm, len, pole;
    int *inv = (int *) NULL;
    double yhi, ylo, xhi, xlo, val, newval, tdelta, delta = 0.0;
    char buf[1024];
 
    printf ("Insert border crossings in tour\n");  fflush (stdout);
    CCutil_cycle_len (ncount, dat, ftour, &val);
    printf ("Starting tour length:       %16.0f\n", val);  fflush (stdout);

    CCutil_dat_getnorm (dat, &norm);

    inv = CC_SAFE_MALLOC (pcount, int);
    CCcheck_NULL (inv, "out of memory in insert_crossings");

    for (i = 0; i < pcount; i++) {
        inv[ptour[i]] = i;
    }

    xlo = trac[0].xrange[0];
    xhi = trac[0].xrange[1];
    ylo = trac[0].yrange[0];
    yhi = trac[0].yrange[1];
    for (i = 1; i < pcount; i++) {
        if (trac[i].xrange[0] < xlo) xlo = trac[i].xrange[0];
        if (trac[i].xrange[1] > xhi) xhi = trac[i].xrange[1];
        if (trac[i].yrange[0] < ylo) ylo = trac[i].yrange[0];
        if (trac[i].yrange[1] > yhi) yhi = trac[i].yrange[1];
    }

    for (i = 0; i < pcount; i++) {
        for (j = i+1; j < pcount; j++) {
            if (inv[i] != inv[j]+1 && inv[i] != inv[j]-1) {
                subdiv_adj (&trac[i], &trac[j], &len);
                if (len == 0) {
                    rval = insert_closest_pair (name, ncount, i, j, dat, trac,
                                         ftour, 0, 0, ylo, yhi, &tdelta, &hit);
                    CCcheck_rval (rval, "insert_closest_pair failed");
                    if (hit) {
                        delta += tdelta; cnt++;
                    }
                } else if (norm == CC_GEOM) {
                    dateline_adj (&trac[i], &trac[j], ylo, yhi, &len);
                    if (len == 0) {
                        rval = insert_closest_pair (name, ncount, i, j, dat,
                                    trac, ftour, 1, 0, ylo, yhi, &tdelta, &hit);
                        CCcheck_rval (rval, "insert_closest_pair failed");
                        if (hit) {
                            delta += tdelta; cnt++;
                        }
                    } else {
                       pole_adj (&trac[i], &trac[j], xlo, xhi, &pole);
                       if (pole) {
                           rval = insert_closest_pair (name, ncount, i, j, dat,
                                            trac, ftour, 0, pole, ylo, yhi,
                                            &tdelta, &hit);
                           CCcheck_rval (rval, "insert_closest_pair failed");
                           if (hit) {
                               delta += tdelta; cnt++;
                           }
                        }
                    }
                }
            }
        }
    }
    printf ("\n");

    printf ("%d total inserts, with delta = %.0f\n", cnt, delta);
    fflush (stdout);

    CCutil_cycle_len (ncount, dat, ftour, &newval);
    printf ("Starting tour length:       %16.0f\n", val);  fflush (stdout);
    printf ("Adjusted tour length:       %16.0f\n", newval);  fflush (stdout);
    if (newval - val != delta) {
         printf ("Note: not the same as delta\n"); fflush (stdout);
    }

    sprintf (buf, "%s_insert.tour", name);
    rval = CCutil_writecycle (ncount, buf, ftour, 0);
    CCcheck_rval (rval, "CCutil_writecycle failed");

CLEANUP:

    CC_IFFREE (inv, int);
    return rval;
}

static int cat_full_tour (int ncount, CCdatagroup *dat, char *name, int pcount,
        CCsubdiv_lkh *lkhtrac, double origval)
{
    int rval = 0;
    int i, j, k = 0;
    int *tour = (int *) NULL;
    char *hits = (char *) NULL;
    int scount;
    int *snam = (int *) NULL;
    int *stour = (int *) NULL;
    char tnam[128], buf[1024];
    double newval;

    printf ("Cat full tour from paths\n"); fflush (stdout);

    sprintf (tnam, "lkpath");

    tour = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (tour, "out of memory in cat_full_tour");

    for (i = 0; i < pcount; i++) {
        rval = get_part_data (name, i, &scount, &snam, &stour, tnam, 1);
        CCcheck_rval (rval, "get_part_data failed");

        if (scount != lkhtrac[i].cnt) {
            fprintf (stderr, "subpath count does not match index\n");
            rval = 1;  goto CLEANUP;
        } 

        for (j = 0; j < scount; j++) {
            tour[k++] = snam[stour[j]];
        }

        CC_IFFREE (snam, int);
        CC_IFFREE (stour, int);
    }

    if (k != ncount) {
        fprintf (stderr, "the subpaths do not make a full tour\n");
        rval = 1;  goto CLEANUP;
    }

    hits = CC_SAFE_MALLOC (ncount, char);
    CCcheck_NULL (hits, "out of memory in cat_full_tour");
    for (i = 0; i < ncount; i++) hits[i] = 0;

    for (i = 0; i < ncount; i++) {
        if (hits[tour[i]]) {
            fprintf (stderr, "duplicate node in tour\n");
            rval = 1;  goto CLEANUP;
        } else {
            hits[tour[i]] = 1;
        }
    }

    CCutil_cycle_len (ncount, dat, tour, &newval);
    printf ("Starting tour length:       %16.0f\n", origval);  fflush (stdout);
    printf ("Adjusted tour length:       %16.0f\n", newval);  fflush (stdout);
    
    sprintf (buf, "%s_lksegment.tour", name);
    rval = CCutil_writecycle (ncount, buf, tour, 0);
    CCcheck_rval (rval, "CCutil_writecycle failed");


CLEANUP:

    CC_IFFREE (snam, int);
    CC_IFFREE (stour, int);
    CC_IFFREE (tour, int);
    CC_IFFREE (hits, char);

    return rval;
}

static int find_full_tour (int ncount, CCdatagroup *dat, char *name,
        int pcount, CCsubdiv *trac, int *ptour)
{
    int i, rval = 0;
    int n0, n1, m0, m1, bcount = 0;
    int *btour = (int *) NULL;
    int *bnamelist = (int *) NULL;
    char *hits = (char *) NULL;
    char buf[1024];
    char tnam[128];
    double bval, sumval;

    sprintf (tnam, "lkcyc");

    printf ("Build full tour from components\n"); fflush (stdout);

    printf ("Partition Path\n");
    for (i = 0; i < pcount; i++) {
        printf ("%d ", ptour[i]);
        if (i % 10 == 9) printf ("\n");
    }
    printf ("\n"); fflush (stdout);

    rval = get_partial_tour_sum (name, pcount, &sumval, dat);
    CCcheck_rval (rval, "get_partial_tour_sum failed");
    printf ("Sum of partial tours:  %.0f\n", sumval); fflush (stdout);

    rval = get_part_data (name, ptour[0], &bcount, &bnamelist, &btour, tnam, 0);
    CCcheck_rval (rval, "get_part_data failed");
    for (i = 0; i < bcount; i++) {
        btour[i] = bnamelist[btour[i]];
    }
    CC_IFFREE (bnamelist, int);

    for (i = 1; i < pcount; i++) {
        rval = find_closest_pair (name, ptour[i-1], ptour[i], &n0, &n1,
                                  &m0, &m1, dat, trac, &bcount, &btour); 
        CCcheck_rval (rval, "find_closest_pair failed");
        printf ("%d %d %d\n", n0, m0, CCutil_dat_edgelen (n0, m0, dat));
        printf ("%d %d %d\n", n1, m1, CCutil_dat_edgelen (n1, m1, dat));
    }

    if (bcount != ncount) {
        fprintf (stderr, "full tour is missing nodes\n");
        rval = 1;  goto CLEANUP;
    }
 
    hits = CC_SAFE_MALLOC (bcount, char);
    CCcheck_NULL (hits, "out of memory in find_full_tour");
    for (i = 0; i < bcount; i++) hits[i] = 0;
    for (i = 0; i < bcount; i++) {
        if (hits[btour[i]] != 0) {
            fprintf (stderr, "repeated node in full tour\n");
            rval = 1; goto CLEANUP;
        }
        hits[btour[i]] = 1;
    }

    sprintf (buf, "%s_full.tour", name);
    rval = CCutil_writecycle (bcount, buf, btour, 0);
    CCcheck_rval (rval, "CCutil_writecycle failed");

    CCutil_cycle_len (bcount, dat, btour, &bval);
    printf ("Length of full tour:       %16.0f\n", bval);  fflush (stdout);
    printf ("Delta from partial  tours: %16.0f\n", bval - sumval);

CLEANUP:

    CC_IFFREE (bnamelist, int);
    return rval;
}

static int improve_border (char *name, int ncount, int ip, int iq,
        CCdatagroup *dat, CCsubdiv *trac, int *tour, int dateline, int pole,
        double ylo, double yhi, CCrandstate *rstate)
{
    int rval = 0;
    CCsubdiv *p = &(trac[ip]);
    CCsubdiv *q = &(trac[iq]);
    CCsubdiv *tempsub;
    int flip = 0, border = 0;
    int i, tempi, pcount, qcount, pbcnt, qbcnt;
    int *pnamelist = (int *) NULL;
    int *qnamelist = (int *) NULL;
    int *pbord = (int *) NULL;
    int *qbord = (int *) NULL;
    int *ulist = (int *) NULL;

    printf ("Improve border %d-%d", ip, iq);
    if (dateline)  printf (" dateline\n");
    else if (pole) printf (" polar\n");
    else           printf ("\n");
    fflush (stdout);

    if (pole) {
        border = 0;
        flip = 0;
    } else if (dateline) {
        border = 1;
        if (       p->yrange[0] == ylo && q->yrange[1] == yhi) {
            flip = 0;
        } else if (p->yrange[1] == yhi && q->yrange[0] == ylo) {
            flip = 1;
        } else {
            fprintf (stderr, "the adjacent pair does not share dateline\n");
            rval = 1;  goto CLEANUP;
        }
    } else {
        find_common_border (p, q, &border, &flip);
        if (border == -1) {
            fprintf (stderr, "the adjacent pair does not share a border\n");
            rval = 1;  goto CLEANUP;
        }
    }

    if (flip) {
        CC_SWAP (p, q, tempsub);
        CC_SWAP (ip, iq, tempi);
    }

    rval = get_part_data (name, ip, &pcount, &pnamelist, (int **) NULL,
                         (char *) NULL, 0);
    CCcheck_rval (rval, "get_part_data failed");

    rval = get_part_data (name, iq, &qcount, &qnamelist, (int **) NULL,
                         (char *) NULL, 0);
    CCcheck_rval (rval, "get_part_data failed");

    if (pole) {
        rval = grab_border_points (pole-1, pcount, pnamelist, dat, &pbcnt,
                                   &pbord, ylo, yhi, LK_POINTS);
        CCcheck_rval (rval, "grab_border_points failed");

        rval = grab_border_points (pole-1, qcount, qnamelist, dat, &qbcnt,
                                   &qbord, ylo, yhi, LK_POINTS);
        CCcheck_rval (rval, "grab_border_points failed");
    } else if (border == 0) {
        rval = grab_border_points (0, pcount, pnamelist, dat, &pbcnt, &pbord,
                                   q->yrange[0], q->yrange[1], LK_POINTS);
        CCcheck_rval (rval, "grab_border_points failed");

        rval = grab_border_points (1, qcount, qnamelist, dat, &qbcnt, &qbord,
                                   p->yrange[0], p->yrange[1], LK_POINTS);
        CCcheck_rval (rval, "grab_border_points failed");
    } else {
        rval = grab_border_points (2, pcount, pnamelist, dat, &pbcnt, &pbord,
                                   q->xrange[0], q->xrange[1], LK_POINTS);
        CCcheck_rval (rval, "grab_border_points failed");

        rval = grab_border_points (3, qcount, qnamelist, dat, &qbcnt, &qbord,
                                   p->xrange[0], p->xrange[1], LK_POINTS);
        CCcheck_rval (rval, "grab_border_points failed");
    }

    ulist = CC_SAFE_MALLOC (pbcnt + qbcnt, int);
    CCcheck_NULL (ulist, "out of memory in improve_border");

    for (i = 0; i < pbcnt; i++) ulist[i] = pbord[i];
    for (i = 0; i < qbcnt; i++) ulist[pbcnt+i] = qbord[i];

    rval = run_lk_subproblem (pbcnt+qbcnt, ulist, ncount, tour, dat, 
                              rstate);
    CCcheck_rval (rval, "run_lk_subproblem failed");

CLEANUP:

    CC_IFFREE (pnamelist, int);
    CC_IFFREE (qnamelist, int);
    CC_IFFREE (ulist, int);
    CC_IFFREE (pbord, int);
    CC_IFFREE (qbord, int);

    return rval;
}

static int insert_closest_pair (char *name, int ncount, int ip, int iq,
        CCdatagroup *dat, CCsubdiv *trac, int *tour, int dateline, int pole,
        double ylo, double yhi, double *tdelta, int *hit)
{
    int rval = 0;
    CCsubdiv *p = &(trac[ip]);
    CCsubdiv *q = &(trac[iq]);
    CCsubdiv *tempsub;
    int flip = 0, border = 0;
    int i, j, k, tempi, pcount, qcount, pbcnt, qbcnt;
    int p0, p1, q0, q1, p0index, p1index, q0index, q1index;
    double t, mini;
    int *pnamelist = (int *) NULL;
    int *qnamelist = (int *) NULL;
    int *newtour = (int *) NULL;
    int *pbord = (int *) NULL;
    int *qbord = (int *) NULL;
    int *ptempi;
    char tnam[128];

    sprintf (tnam, "lkcyc");

    printf ("Insert border %d-%d", ip, iq);
    if (dateline)  printf (" dateline\n");
    else if (pole) printf (" polar\n");
    else           printf ("\n");
    fflush (stdout);

    if (pole) {
        border = 0;
        flip = 0;
    } else if (dateline) {
        border = 1;
        if (       p->yrange[0] == ylo && q->yrange[1] == yhi) {
            flip = 0;
        } else if (p->yrange[1] == yhi && q->yrange[0] == ylo) {
            flip = 1;
        } else {
            fprintf (stderr, "the adjacent pair does not share dateline\n");
            rval = 1;  goto CLEANUP;
        }
    } else {
        find_common_border (p, q, &border, &flip);
        if (border == -1) {
            fprintf (stderr, "the adjacent pair does not share a border\n");
            rval = 1;  goto CLEANUP;
        }
    }

    if (flip) {
        CC_SWAP (p, q, tempsub);
        CC_SWAP (ip, iq, tempi);
    }

    rval = get_part_data (name, ip, &pcount, &pnamelist, (int **) NULL, tnam,
                          0);
    CCcheck_rval (rval, "get_part_data failed");

    rval = get_part_data (name, iq, &qcount, &qnamelist, (int **) NULL, tnam,
                          0);
    CCcheck_rval (rval, "get_part_data failed");

    if (pole) {
        rval = grab_border_points (pole-1, pcount, pnamelist, dat, &pbcnt,
                                   &pbord, ylo, yhi, B_POINTS);
        CCcheck_rval (rval, "grab_border_points failed");

        rval = grab_border_points (pole-1, qcount, qnamelist, dat, &qbcnt,
                                   &qbord, ylo, yhi, B_POINTS);
        CCcheck_rval (rval, "grab_border_points failed");
    } else if (border == 0) {
        rval = grab_border_points (0, pcount, pnamelist, dat, &pbcnt, &pbord,
                                   q->yrange[0], q->yrange[1], B_POINTS);
        CCcheck_rval (rval, "grab_border_points failed");

        rval = grab_border_points (1, qcount, qnamelist, dat, &qbcnt, &qbord,
                                   p->yrange[0], p->yrange[1], B_POINTS);
        CCcheck_rval (rval, "grab_border_points failed");
    } else {
        rval = grab_border_points (2, pcount, pnamelist, dat, &pbcnt, &pbord,
                                   q->xrange[0], q->xrange[1], B_POINTS);
        CCcheck_rval (rval, "grab_border_points failed");

        rval = grab_border_points (3, qcount, qnamelist, dat, &qbcnt, &qbord,
                                   p->xrange[0], p->xrange[1], B_POINTS);
        CCcheck_rval (rval, "grab_border_points failed");
    }

    mini = 1e30;
    p0 = q0 = p0index = q0index = -1;
    for (i = 0; i < pbcnt; i++) {
        for (j = 0; j < qbcnt; j++) {
            t = CCutil_dat_edgelen (pbord[i], qbord[j], dat);
            if (t < mini) {
                mini = t;
                p0 = pbord[i];
                q0 = qbord[j];
            }
        }
    }

    if (flip) {
        CC_SWAP (p, q, tempsub);
        CC_SWAP (ip, iq, tempi);
        CC_SWAP (p0, q0, tempi);
        CC_SWAP (pcount, qcount, tempi);
        CC_SWAP (pnamelist, qnamelist, ptempi);
    }

    for (i = 0; i < ncount; i++) {
        if (tour[i] == p0) p0index = i;
        if (tour[i] == q0) q0index = i;
    }

    if (p0index == ncount-1) {
        p1 = tour[0]; p1index = 0;
    } else {
        p1 = tour[p0index+1]; p1index = p0index+1;
    }

    if (q0index == ncount-1) {
        q1 = tour[0]; q1index = 0;
    } else {
        q1 = tour[q0index+1]; q1index = q0index+1;
    }

    if (p0 == q1 || p1 == q0) {
        fprintf (stderr, "paired edges not disjoint: p0 = %d, p1 = %d, q0 = %d, q1 = %d\n", p0, p1, q0, q1);
        fprintf (stderr, "skipping this pair\n");
        *hit = 0;
        goto CLEANUP;
    }

    *tdelta = (double) (CCutil_dat_edgelen (q0, p0, dat) +
                        CCutil_dat_edgelen (q1, p1, dat) -
                        CCutil_dat_edgelen (p0, p1, dat) -
                        CCutil_dat_edgelen (q0, q1, dat)); 
    printf ("Swap (%d,%d), (%d,%d), delta = %.0f\n", p0, p1, q0, q1, *tdelta);
    fflush (stdout);
/*
    printf ("EDGE %d %d 1\n", q0, p0);
    printf ("EDGE %d %d 1\n", q1, p1);
*/

    newtour = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (newtour, "out of memory in insert_closest_pair");

    /* newtour is p1-->q0, p0-(rev)->q1 */

    k = 0;
    if (q0index < p1index) {
         for (i = p1index; i < ncount; i++) newtour[k++] = tour[i];
         for (i = 0; i <= q0index; i++) newtour[k++] = tour[i];
    } else {
         for (i = p1index; i <= q0index; i++) newtour[k++] = tour[i];
    }

    if (q1index > p0index) {
        for (i = p0index; i >= 0; i--) newtour[k++] = tour[i];
        for (i = ncount-1; i >= q1index; i--) newtour[k++] = tour[i];
    } else {
        for (i = p0index; i >= q1index; i--) newtour[k++] = tour[i];
    }

    if (k != ncount) {
        fprintf (stderr, "error in 2-swap\n");  
        rval = 1;  goto CLEANUP;
    }

    for (i = 0; i < ncount; i++) tour[i] = newtour[i];

    *hit = 1;

CLEANUP:

    CC_IFFREE (pnamelist, int);
    CC_IFFREE (qnamelist, int);
    CC_IFFREE (newtour, int);
    CC_IFFREE (pbord, int);
    CC_IFFREE (qbord, int);

    return rval;
}

static int find_closest_pair (char *name, int ip, int iq, int *p_p0, int *p_p1,
        int *p_q0, int *p_q1, CCdatagroup *dat, CCsubdiv *trac, int *p_bcount,
        int **p_btour)
{
    int rval = 0;
    int pcount, qcount, i, j, k, tempi, pbcnt, qbcnt;
    int p0 = -1, p1 = -1, q0 = -1, q1 = -1;
    int p0index, p1index, q0index, q1index;
    int *pnamelist = (int *) NULL;
    int *qnamelist = (int *) NULL;
    int *ptour = (int *) NULL;
    int *qtour = (int *) NULL;
    int *newtour = (int *) NULL;
    int *pbord = (int *) NULL;
    int *qbord = (int *) NULL;
    int flip = 0, border = 0;
    CCsubdiv *p = &(trac[ip]);
    CCsubdiv *q = &(trac[iq]);
    CCsubdiv *tempsub;
    double mini, t;
    int *btour = *p_btour;
    int bcount = *p_bcount;
    int *ptempi;
    char tnam[128];

    sprintf (tnam, "lkcyc");

    find_common_border (p, q, &border, &flip);
    if (border == -1) {
        fprintf (stderr, "the adjacent pair does not share a border\n");
        rval = 1;  goto CLEANUP;
    }

    if (flip) {
        CC_SWAP (p, q, tempsub);
        CC_SWAP (ip, iq, tempi);
    }

    rval = get_part_data (name, ip, &pcount, &pnamelist, &ptour, tnam, 0);
    CCcheck_rval (rval, "get_part_data failed");

    rval = get_part_data (name, iq, &qcount, &qnamelist, &qtour, tnam, 0);
    CCcheck_rval (rval, "get_part_data failed");

    if (border == 0) {
        rval = grab_border_points (0, pcount, pnamelist, dat, &pbcnt, &pbord,
                                   q->yrange[0], q->yrange[1], B_POINTS);
        CCcheck_rval (rval, "grab_border_points failed");

        rval = grab_border_points (1, qcount, qnamelist, dat, &qbcnt, &qbord,
                                   p->yrange[0], p->yrange[1], B_POINTS);
        CCcheck_rval (rval, "grab_border_points failed");
    } else {
        rval = grab_border_points (2, pcount, pnamelist, dat, &pbcnt, &pbord,
                                   q->xrange[0], q->xrange[1], B_POINTS);
        CCcheck_rval (rval, "grab_border_points failed");

        rval = grab_border_points (3, qcount, qnamelist, dat, &qbcnt, &qbord,
                                   p->xrange[0], p->xrange[1], B_POINTS);
        CCcheck_rval (rval, "grab_border_points failed");
    }

    mini = 1e30;
    for (i = 0; i < pbcnt; i++) {
        for (j = 0; j < qbcnt; j++) {
            t = CCutil_dat_edgelen (pbord[i], qbord[j], dat);
            if (t < mini) {
                mini = t;
                p0 = pbord[i];
                q0 = qbord[j];
            }
        }
    }

    if (p0 == -1) {
        fprintf (stderr, "could not find close pair\n");
        rval = 1;  goto CLEANUP;
    }

    if (flip) {
        CC_SWAP (p, q, tempsub);
        CC_SWAP (ip, iq, tempi);
        CC_SWAP (p0, q0, tempi);
        CC_SWAP (p1, q1, tempi);
        CC_SWAP (pcount, qcount, tempi);
        CC_SWAP (ptour, qtour, ptempi);
        CC_SWAP (pnamelist, qnamelist, ptempi);
    }

    newtour = CC_SAFE_MALLOC (bcount + qcount, int);
    CCcheck_NULL (newtour, "out of memory in find_closest_pair");

    p0index = -1;
    for (i = 0; i < bcount; i++) {
        if (btour[i] == p0) p0index = i; 
    }
    if (p0index == -1) {
        fprintf (stderr, "error in btour\n");
        rval = 1;  goto CLEANUP;
    }
    if (p0index < bcount-1) p1index = p0index+1;
    else                    p1index = 0;

    for (i = 0; i < qcount; i++) {
        qtour[i] = qnamelist[qtour[i]];
    }
    q0index = -1;
    for (i = 0; i < qcount; i++) {
        if (qtour[i] == q0) q0index = i; 
    }
    if (q0index == -1) {
        fprintf (stderr, "error in partial tour B\n");
        rval = 1;  goto CLEANUP;
    }
    if (q0index < qcount-1) q1index = q0index+1;
    else                    q1index = 0;

    p1 = btour[p1index];
    q1 = qtour[q1index];

    if (CCutil_dat_edgelen (p0, q0, dat) + CCutil_dat_edgelen (p1, q1, dat) >
        CCutil_dat_edgelen (p0, q1, dat) + CCutil_dat_edgelen (p1, q0, dat)) {
        CC_SWAP (p0, p1, tempi);
        CC_SWAP (p0index, p1index, tempi);
    } 


    /* Order the btour from p1 to p0 */

    k = 0;
    if (p1index == 0 && p0index == bcount-1) {
        for (i = 0; i < bcount; i++) newtour[k++] = btour[i];
    } else if (p1index == bcount-1 && p0index == 0) {
        for (i = bcount-1; i >= 0; i--) newtour[k++] = btour[i];
    } else if (p1index == p0index+1) {
        for (i = p1index; i < bcount; i++) newtour[k++] = btour[i];
        for (i = 0; i < p1index; i++) newtour[k++] = btour[i];
    } else if (p1index == p0index-1) {
        for (i = p1index; i >= 0; i--) newtour[k++] = btour[i];
        for (i = bcount-1; i > p1index; i--)  newtour[k++] = btour[i];
    } else {
        printf ("p-pair not adjacent, use alternative\n");
        for (i = p1index; i < bcount; i++) newtour[k++] = btour[i];
        for (i = 0; i < p1index; i++) newtour[k++] = btour[i];
        p0 = btour[k-1]; 
    }

 
    /* Order the qtour from q0 to q1 */

    if (q0index == 0 && q1index == qcount-1) {
        for (i = 0; i < qcount; i++) newtour[k++] = qtour[i];
    } else if (q0index == qcount-1 && q1index == 0) {
        for (i = qcount-1; i >= 0; i--) newtour[k++] = qtour[i];
    } else if (q0index == q1index+1) {
        for (i = q0index; i < qcount; i++) newtour[k++] = qtour[i];
        for (i = 0; i < q0index; i++) newtour[k++] = qtour[i];
    } else if (q0index == q1index-1) {
        for (i = q0index; i >= 0; i--) newtour[k++] = qtour[i];
        for (i = qcount-1; i > q0index; i--)  newtour[k++] = qtour[i];
    } else {
        fprintf (stderr, "error in partial tour C\n");
        rval = 1;  goto CLEANUP;
    }

    CC_IFFREE (btour, int);
    *p_btour = newtour;
    *p_bcount = k;

    *p_p0 = p0; *p_p1 = p1;
    *p_q0 = q0; *p_q1 = q1;

CLEANUP:

    CC_IFFREE (pnamelist, int);
    CC_IFFREE (qnamelist, int);
    CC_IFFREE (ptour, int);
    CC_IFFREE (qtour, int);
    CC_IFFREE (pbord, int);
    CC_IFFREE (qbord, int);
    if (rval) {
        CC_IFFREE (newtour, int);
    }
    return rval;
}

static void find_common_border (CCsubdiv *p, CCsubdiv *q, int *border,
        int *flip)
{
    if (p->xrange[0] == q->xrange[1]) {
        *border = 0;  *flip = 0;
    } else if (p->xrange[1] == q->xrange[0]) {
        *border = 0;  *flip = 1;
    } else if (p->yrange[0] == q->yrange[1]) {
        *border = 1;  *flip = 0;
    } else if (p->yrange[1] == q->yrange[0]) {
        *border = 1;  *flip = 1;
    } else {
        *border = -1; *flip = 0;
    }
}

static int grab_border_points (int btype, int scount, int *slist,
        CCdatagroup *dat, int *p_bcount, int **p_blist, double xlow,
        double xhi, int maxpoints)
{
    int rval = 0;
    int bcount = 0;
    int *blist = (int *) NULL;
    double *x = (double *) NULL;
    int *xnames = (int *) NULL;
    int *xperm = (int *) NULL;
    int i, k, xcount = 0;

    x = CC_SAFE_MALLOC (scount, double);
    CCcheck_NULL (x, "out of memory in grab_border_points");

    xnames = CC_SAFE_MALLOC (scount, int);
    CCcheck_NULL (xnames, "out of memory in grab_border_points");

    if (btype == 0 || btype == 1) {
        for (i = 0; i < scount; i++) {
            k = slist[i];
            if (dat->y[k] >= xlow && dat->y[k] <= xhi) {
                xnames[xcount] = k;
                if (btype == 0) x[xcount++] = dat->x[k];
                else            x[xcount++] = -dat->x[k];
            }
        }
    } else {
        for (i = 0; i < scount; i++) {
            k = slist[i];
            if (dat->x[k] >= xlow && dat->x[k] <= xhi) {
                xnames[xcount] = k;
                if (btype == 2) x[xcount++] = dat->y[k];
                else            x[xcount++] = -dat->y[k];
            }
        }
    }

/*
    printf ("xcount = %d, btype = %d\n", xcount, btype); fflush (stdout);
*/

    if (xcount == 0) {
        fprintf (stderr, "no points to grab\n");
        rval = 1;  goto CLEANUP;
    }

    xperm = CC_SAFE_MALLOC (xcount, int);
    CCcheck_NULL (xperm, "out of memory in grab_border_points");

    for (i = 0; i < xcount; i++) {
        xperm[i] = i;
    }

    CCutil_double_perm_quicksort (xperm, x, xcount);
    if (xcount < maxpoints) bcount = xcount;
    else                    bcount = maxpoints;

    blist = CC_SAFE_MALLOC (bcount, int);
    CCcheck_NULL (blist, "out of memory in grab_border_points");

    for (i = 0; i < bcount; i++) {
        blist[i] = xnames[xperm[i]];
    }

    *p_bcount = bcount;
    *p_blist = blist;

CLEANUP:

    CC_IFFREE (x, double);
    CC_IFFREE (xnames, int);
    CC_IFFREE (xperm, int);
    if (rval) {
        CC_IFFREE (blist, int);
    }
    return rval;
}

static int  get_partial_tour_sum (char *name, int pcount, double *sumval,
        CCdatagroup *dat)
{
    int *tour = (int *) NULL;
    int *namelist = (int *) NULL;
    int i, j, count;
    double val = 0.0;
    int rval = 0;
    char tnam[128];

    sprintf (tnam, "lkcyc");

    for (i = 0; i < pcount; i++) {
        rval = get_part_data (name, i, &count, &namelist, &tour, tnam, 0);
        CCcheck_rval (rval, "get_part_data failed");
        for (j = 1; j < count; j++) {
            val += CCutil_dat_edgelen (namelist[tour[j-1]], namelist[tour[j]],
                                       dat);
        }
        val += CCutil_dat_edgelen (namelist[tour[count-1]], namelist[tour[0]],
                                   dat);
        CC_IFFREE (tour, int);
        CC_IFFREE (namelist, int);
    }

    *sumval = val;

CLEANUP:

    CC_IFFREE (tour, int);
    CC_IFFREE (namelist, int);
    return rval;
}

static int get_part_data (char *name, int ind, int *p_count, int **p_namelist,
        int **p_tour, char *tname, int use_default)
{
    int rval = 0;
    char nambuf[1024];
    char tourbuf[1024];
    FILE *in = (FILE *) NULL;
    int i, count = 0;
    int *namelist = (int *) NULL;
    int *tour = (int *) NULL;

    sprintf (nambuf, "%s_%d.nam", name, ind);
    if (p_tour) {
        sprintf (tourbuf, "%s.%d", tname, ind);
    }

    in = fopen (nambuf, "r");
    if (!in) {
        fprintf (stderr, "could not open %s for reading\n", nambuf);
        rval = 1; goto CLEANUP;
    }
    fscanf (in, "%d", &count);

    if (p_namelist) {
        namelist = CC_SAFE_MALLOC (count, int);
        CCcheck_NULL (namelist, "out of memory in get_part_data");
        for (i = 0; i < count; i++) {
            fscanf (in, "%d", &(namelist[i]));
        }
    }
    fclose (in);
    in = (FILE *) NULL;

    if (p_tour) {
        tour = CC_SAFE_MALLOC (count, int);
        CCcheck_NULL (tour, "out of memory in get_part_data");
        rval = CCutil_getcycle (count, tourbuf, tour, 0);
        if (rval) {
            if (use_default) {
                printf ("Using default permutation tour\n"); fflush (stdout);
                for (i = 0; i < count; i++) tour[i] = i;
                rval = 0;
            } else {
                fprintf (stderr, "CCutil_getcycle failed\n");
                goto CLEANUP;
            }
        }
    }

/*
    if (hackarray[ind] == 0) {
        for (i = 1; i < count; i++) {
            printf ("%d %d 1\n", namelist[tour[i-1]], namelist[tour[i]]);
        }
        printf ("%d %d 1\n", namelist[tour[count-1]], namelist[tour[0]]);
        hackarray[ind] = 1;
    }
*/


    if (p_count) *p_count = count;
    if (p_tour) *p_tour = tour;
    if (p_namelist) *p_namelist = namelist;

CLEANUP:

    if (rval) {
        CC_IFFREE (tour, int);
        CC_IFFREE (namelist, int);
    }

    if (in) fclose (in);
    return rval;
}

static int run_lk_subproblem (int scount, int *slist, int ncount, int *tour,
        CCdatagroup *dat, CCrandstate *rstate)
{
    int rval = 0;
    CCdatagroup sdat;
    CCedgegengroup plan;
    int *goodlist = (int *) NULL;
    int *stour = (int *) NULL;
    int *fixed = (int *) NULL, *sfixed = (int *) NULL;
    int *newstour = (int *) NULL;
    int *new_tour = (int *) NULL;
    int i, j, goodcount, inorm;
    double sval;
    char *in_out = (char *) NULL;
    int *big_nr = (int *) NULL;   /* array of length scount */
    int *small_nr = (int *) NULL;   /* array of length scount */
    int *big_tournr = (int *) NULL;   /* array of length scount */
    int fixed_count;
    double *x, *y;

    CCutil_init_datagroup (&sdat);

    printf ("Improve subproblem with %d points\n", scount);  fflush (stdout);

    CCedgegen_init_edgegengroup (&plan);
    plan.quadnearest = 2;

    CCutil_dat_getnorm (dat, &inorm);
    CCutil_dat_setnorm (&sdat, inorm);

    small_nr = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (small_nr, "out of memory in run_lk_linkern");

    in_out = CC_SAFE_MALLOC (ncount, char);
    CCcheck_NULL (in_out, "out of memory in run_lk_linkern");

    for (i=0;i<ncount;i++) in_out[i] = 0;
    for (i=0;i<scount;i++) in_out[slist[i]] = 1;

    /** x, y and stour will be sorted like they appear in tour */

    x = CC_SAFE_MALLOC (scount, double);
    CCcheck_NULL (x, "out of memory in build_tour_subproblem");
    y = CC_SAFE_MALLOC (scount, double);
    CCcheck_NULL (y, "out of memory in build_tour_subproblem");

    stour = CC_SAFE_MALLOC (scount, int);
    CCcheck_NULL (stour, "out of memory in run_lk_subproblem");

    big_nr = CC_SAFE_MALLOC (scount, int);
    CCcheck_NULL (big_nr, "out of memory in run_lk_subproblem");

    big_tournr = CC_SAFE_MALLOC (scount, int);
    CCcheck_NULL (big_tournr, "out of memory in run_lk_subproblem");

  
    /* calc x and y, don't use the isolated points */
    /* also calc stour, small_nr, & big_nr */

    for (scount=0,i=0;i<ncount;i++) {
        if (in_out[tour[i]]==1) {
            if (in_out[tour[(i-1+ncount)%ncount]]==1 ||
	        in_out[tour[(i+1)%ncount]]==1) {
	        y[scount]=dat->y[tour[i]];
	        x[scount]=dat->x[tour[i]];
	        stour[scount] = scount;
	        small_nr[tour[i]]=scount;
	        big_nr[scount]=i;
	        big_tournr[scount]=tour[i];
	        scount++;
            } else {
	        in_out[tour[i]]=0;
	        small_nr[tour[i]]=-1;
            }
        }
    }

    if (scount < 10) {
        printf ("Skipping subproblem with %d points\n", scount);
        fflush (stdout);
        goto CLEANUP;
    }
  
    for (fixed_count=0,i=0;i<scount;i++) {
        if ( ((big_nr[i]+1)%ncount) != big_nr[(i+1)%scount] ) {
            fixed_count++;
        }
    }

    fixed = CC_SAFE_MALLOC (2*fixed_count, int);
    CCcheck_NULL (fixed, "out of memory in run_lk_subproblem");
    sfixed = CC_SAFE_MALLOC (scount, int);
    CCcheck_NULL (sfixed, "out of memory in run_lk_subproblem");

    for (i=0;i<scount;i++) 
        sfixed[i] = 0;
    /* set all fixed edges in array. */
    for (fixed_count=0,i=0;i<scount;i++) {
        if ( ((big_nr[i]+1)%ncount) != big_nr[(i+1)%scount] ) {
            fixed[2*fixed_count]=i;
            fixed[2*fixed_count+1]=(i+1)%scount;      
            fixed_count++;
            sfixed[i]=fixed_count;
            sfixed[(i+1)%scount]=fixed_count;
        }
    }

    printf ("Tour had %d fixed edges\n", fixed_count); fflush (stdout);
    printf ("scount = %d\n", scount);

    sdat.x = x;
    sdat.y = y;

    rval = CCedgegen_edges (&plan, scount, &sdat, (double *) NULL,
                            &goodcount, &goodlist, 0, rstate);
    CCcheck_rval (rval, "CCedgegen_edges failed");

    newstour = CC_SAFE_MALLOC (scount, int);
    CCcheck_NULL (newstour, "out of memory in run_lk_subproblem");

    rval = CClinkern_fixed (scount, &sdat, goodcount, goodlist,
                 scount, stour, newstour, &sval, fixed_count, fixed, 0, rstate);

    {
        char *hit;
        hit = CC_SAFE_MALLOC (scount, char);

        for (i = 0; i < scount; i++) hit[i] = 0;
        for (i = 0; i < scount; i++) {
            if (hit[newstour[i]]) {
                fprintf (stderr, "BAD NEWS LK %d\n", i);
                rval = 1;  goto CLEANUP;
            } else {
                hit[newstour[i]] = 1;
            }
        }
        free (hit);
    }

    /* output of linkern is a new tour in an array newstour */

    new_tour = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (new_tour, "out of memory in run_lk_subproblem");

    /* make the new_tour */  
    j = 0; /* this is the counter for the big tour */
    for (i=0;i<scount;i++) {
        int snode = newstour[i];
        int next_snode = newstour[(i+1)%scount];

        if (j >= ncount)  {
            fprintf (stderr, "Counting too many nodes in tour patch\n");
            rval = 1;  goto CLEANUP;
        }

        if (sfixed[snode]==0 || sfixed[snode]!=sfixed[next_snode]) {	
            /* Take single edge from newstour */
            /* printf ("1: %d-%d\n", j, big_tournr[next_snode]); */
            new_tour[j++] = big_tournr[next_snode];
        } else {
            int i2;
            /* Take set of edges from the big tour */
            int big_tour_snode = big_nr[snode];
            if (small_nr[tour[big_tour_snode]] != snode)  {
                printf("Something is wrong !!\n");
                rval = 1;  goto CLEANUP;
            }

            i2 = big_tour_snode;
            do { /** loop until we are at a node inside subset again */
  	        i2 = (i2+1)%ncount;
            } while (in_out[tour[i2]]==0);
            if (small_nr[tour[i2]] == next_snode) {
  	        /* go forward */
	        i2 = big_tour_snode;
	        do { /** loop until we are at a node inside subset again */
	            i2 = (i2+1)%ncount;
                    /* printf ("2: %d-%d\n", j, tour[i2]); */
	            new_tour[j++] = tour[i2];
	        } while (in_out[tour[i2]]==0);
	        if (small_nr[tour[i2]] != next_snode) {
	            printf("Something is wrong II !!\n");	  
                    rval = 1;  goto CLEANUP;
                }
            } else {
	         /* go backward */
	        i2 = big_tour_snode;
                do { /** loop until we are at a node inside subset again */
	            i2 = (i2-1+ncount)%ncount;
                    /* printf ("3: %d-%d\n", j, tour[i2]); */
	            new_tour[j++] = tour[i2];
	        } while (in_out[tour[i2]]==0);
	        if (small_nr[tour[i2]] != next_snode) {
	            printf("Something is wrong III !!\n");	  
                    rval = 1;  goto CLEANUP;
                }
            }
        }
    }

    for (i = 0; i < ncount; i++) tour[i] = new_tour[i];
 
    {
        char *hit;
        hit = CC_SAFE_MALLOC (ncount, char);

        for (i = 0; i < ncount; i++) hit[i] = 0;
        for (i = 0; i < ncount; i++) {
            if (hit[tour[i]]) {
                fprintf (stderr, "BAD NEWS %d - %d\n", i, tour[i]);
                exit (1);
            } else {
                hit[tour[i]] = 1;
            }
        }
        free (hit);
    }

    printf ("End of function\n"); fflush (stdout);

CLEANUP:

    CC_IFFREE (stour, int);
    CC_IFFREE (newstour, int);
    CC_IFFREE (new_tour, int);
    CC_IFFREE (goodlist, int);
    CC_IFFREE (fixed, int);
    CC_IFFREE (sfixed, int);
    CC_IFFREE (in_out, char);
    CC_IFFREE (small_nr, int);
    CC_IFFREE (big_nr, int);
    CC_IFFREE (big_tournr, int);
    CCutil_freedatagroup (&sdat);

    return rval;
}

static int find_part_tour (char *name, int scount, CCsubdiv *trac)
{
    int i, j, len, rval = 0;
    char buf[1024];
    FILE *out = (FILE *) NULL;

    sprintf (buf, "%s_part.edg", name);

    out = fopen (buf, "w");
    if (!out) {
        fprintf (stderr, "could not open %s for output\n", buf);
        rval = 1;  goto CLEANUP;
    }

    fprintf (out, "%d %d\n", scount, (scount * (scount - 1)) / 2);
    for (i = 0; i < scount; i++) {
        for (j = i+1; j < scount; j++) {
            subdiv_adj (&(trac[i]), &(trac[j]), &len);
            fprintf (out, "%d %d %d\n", i, j, len);
        }
    }

CLEANUP:

    if (out) fclose (out);
    return rval;
}

static void subdiv_adj (CCsubdiv *s, CCsubdiv *t, int *len)
{
    *len = 1;

    if ((s->xrange[0] == t->xrange[1]) || (s->xrange[1] == t->xrange[0])) {
        if ((t->yrange[0] <= s->yrange[1] && t->yrange[0] >= s->yrange[0]) ||
            (t->yrange[1] <= s->yrange[1] && t->yrange[1] >= s->yrange[0]) ||
            (s->yrange[1] <= t->yrange[1] && s->yrange[1] >= t->yrange[0])) {
                *len = 0;
        }
    } else if ((s->yrange[0]==t->yrange[1]) || (s->yrange[1]==t->yrange[0])) {
        if ((t->xrange[0] <= s->xrange[1] && t->xrange[0] >= s->xrange[0]) ||
            (t->xrange[1] <= s->xrange[1] && t->xrange[1] >= s->xrange[0]) ||
            (s->xrange[1] <= t->xrange[1] && s->xrange[1] >= t->xrange[0])) {
                *len = 0;
        }
    }
}

static void dateline_adj (CCsubdiv *s, CCsubdiv *t, double ylo, double yhi,
        int *len)
{
    *len = 1;

    if ((s->yrange[0] == ylo && t->yrange[1] == yhi) || 
        (s->yrange[1] == yhi && t->yrange[0]==ylo)) { 
        if ((t->xrange[0] <= s->xrange[1] && t->xrange[0] >= s->xrange[0]) ||
            (t->xrange[1] <= s->xrange[1] && t->xrange[1] >= s->xrange[0]) ||
            (s->xrange[1] <= t->xrange[1] && s->xrange[1] >= t->xrange[0])) {
                *len = 0;
        }
    }
}

static void pole_adj (CCsubdiv *s, CCsubdiv *t, double xlo, double xhi,
        int *pole)
{
    *pole = 0;

    if (       s->xrange[1] == xhi && t->xrange[1] == xhi) {
        *pole = 2;
    } else if (s->xrange[0] == xlo && t->xrange[0] == xlo) {
        *pole = 1;
    }
}

static int parseargs (int ac, char **av)
{
    int c, inorm;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "bcBD:I:P:rs:T:qwN:?", &boptind, &boptarg)) != EOF) { 
        switch (c) {
        case 's':
            seed = atoi(boptarg);
            break;
        case 'B':
            binary_in = 1;
            break;
        case 'b':
            binary_in = 2;
            break;
        case 'c':
            cat_tour = 1;
            break;
        case 'D':
            tspfname = boptarg;
            break;
        case 'I':
            indexfname = boptarg;
            break;
        case 'q':
            dump_partition_tsp = 1;
            break;
        case 'P':
            tourfname  = boptarg;
            break;
        case 'r':
            border_optimization = 1;
            break;
        case 'T':
            fulltourfname  = boptarg;
            break;
        case 'w':
            border_crossings = 1;
            break;
        case 'N':
            inorm = atoi (boptarg);
            switch (inorm) {
            case 0: innorm = CC_MAXNORM; break;
            case 1: innorm = CC_MANNORM; break;
            case 2: innorm = CC_EUCLIDEAN; break;
            case 3: innorm = CC_EUCLIDEAN_3D; break;
            case 17: innorm = CC_GEOM; break;
            case 18: innorm = CC_EUCLIDEAN_CEIL; break;
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
        fprintf (stderr, "extra items\n");
        usage (av[0]);
        return 1;
    }

    return 0;
}

static void usage (char *fname)
{
    fprintf (stderr, "Usage: %s [-flags below]\n", fname);
    fprintf (stderr, "   -b    datfile in double binary format\n");
    fprintf (stderr, "   -B    datfile in integer binary format\n");
    fprintf (stderr, "   -c    cat a tour (need -I -D)\n");
    fprintf (stderr, "   -D f  specify a tsp_or_dat file\n");
    fprintf (stderr, "   -I f  specify an index file\n");
    fprintf (stderr, "   -P f  specify a path file (for partitions)\n");
    fprintf (stderr, "   -r    use border reoptimization (need -I -T -D)\n");
    fprintf (stderr, "   -T f  specify a tour file (for full set)\n");
    fprintf (stderr, "   -w    insert border crossings (need -I -P -T -D)\n");
    fprintf (stderr, "   -q    dump partition tsp as edg file (need -I)\n");
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

