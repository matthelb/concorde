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
#include "bigguy.h"

#define EPSILON ((1.0/(1<<15))*(1.0/(1<<15)))


static CCbigguy
    bgrand (double *delta, int maxx, CCrandstate *rstate);

static int
    bgcheck (CC_SFILE *sfin, CC_SFILE *sfout, CCbigguy x, const char *name),
    dotest (int maxx, int maxy, int maxm, int ntrial, double *delta,
        CCbigguy *sum, CC_SFILE *sfin, CC_SFILE *sfout, CCrandstate *rstate);

static void
    sum_chk (CCbigguy *sum, CCbigguy x);

static int
    intrand (int maxi, CCrandstate *rstate);


#define MOREDELTA(sum,bgx,dbx,msg) { \
    if (fabs((dbx)-CCbigguy_bigguytod((bgx))) > 1e-13 * fabs((dbx))) { \
        fprintf (stderr, "Warning: (%s) %.16f != %.16f by %.16f\n", (msg), (dbx), CCbigguy_bigguytod((bgx)), fabs((dbx)-CCbigguy_bigguytod((bgx)))); \
    } \
   sum += fabs((dbx)-CCbigguy_bigguytod((bgx))); \
}
                                  
int main (int ac, char **av)
{
    int i;
    int x, y, d;
    int nerrs = 0;
    CC_SFILE *sfin = CCutil_sdopen (0, "r");
    CC_SFILE *sfout = CCutil_sdopen (1, "w");
    double delta;
    CCbigguy sum;
    CCrandstate rstate;

    CCutil_sprand (47, &rstate);

    if (ac > 1) {
        fprintf (stderr, "Usage: %s < std_file > new_file\n", av[0]);
        return 1;
    }

    nerrs += bgcheck (sfin, sfout, CCbigguy_MINBIGGUY, "CCbigguy_MINBIGGUY");
    nerrs += bgcheck (sfin, sfout, CCbigguy_MAXBIGGUY, "CCbigguy_MAXBIGGUY");
    nerrs += bgcheck (sfin, sfout, CCbigguy_ZERO, "CCbigguy_ZERO");
    nerrs += bgcheck (sfin, sfout, CCbigguy_ONE, "CCbigguy_ONE");

    delta = 0.0;

    sum = CCbigguy_ZERO;
    nerrs += bgcheck (sfin, sfout, sum, "sum_init");

    for (i=6; i<=30; i+=8) {
        for (x=(1<<i), y=(1<<i); y >= 1; y /= 16, x += 15*y) {
            for (d=y; d>0; d /= 6) {
                nerrs += dotest (x, y/d, d, 100, &delta, &sum, sfin, sfout,
                                 &rstate);
            }
            for (d=y; d>0; d /= 6) {
                nerrs += dotest (y, x/d, d, 100, &delta, &sum, sfin, sfout,
                                 &rstate);
            }
        }
    }

    CCutil_sclose (sfin);
    CCutil_sclose (sfout);

    fprintf (stderr, "%d total errors\n", nerrs);
    fprintf (stderr, "Total delta: %.20f\n", delta);

    return nerrs;
}

static int dotest (int maxx, int maxy, int maxm, int ntrial, double *delta,
        CCbigguy *sum, CC_SFILE *sfin, CC_SFILE *sfout, CCrandstate *rstate)
{
    int i;
    int j;
    CCbigguy x;
    double dx;
    CCbigguy y;
    double dy;
    CCbigguy a, b, c, d, e, f, g;
    int m;
    int nerrs = 0;

    for (i=0; i<ntrial; i++) {
        x = bgrand (delta, maxx, rstate);
        dx = CCbigguy_bigguytod(x);
        y = bgrand (delta, maxy, rstate);
        dy = CCbigguy_bigguytod(y);
        m = intrand (maxm, rstate);

        a = x;
        CCbigguy_add (&a, y);
        MOREDELTA (*delta, a, dx + dy, "add");

        b = x;
        CCbigguy_sub (&b, y);
        MOREDELTA (*delta, b, dx - dy, "sub");

        c = x;
        CCbigguy_addmult (&c, y, CCbigguy_cmp (x, y));
        MOREDELTA (*delta, c, (dx < dy) ? dx-dy : dx+dy, "addmult cmp");

        d = x;
        CCbigguy_addmult (&d, y, m);
        MOREDELTA (*delta, d, dx + m*dy, "addmult");

        e = CCbigguy_ceil (x);
        MOREDELTA (*delta, e, ceil(dx), "ceil");

        j = CCutil_lprand (rstate);

        f = CCbigguy_itobigguy (j);
        MOREDELTA (*delta, f, (double) j, "itobigguy");

        g = CCbigguy_itobigguy (-j);
        MOREDELTA (*delta, g, (double) -j, "itobigguy neg");

        sum_chk (sum, x);
        sum_chk (sum, y);
        sum_chk (sum, a);
        sum_chk (sum, b);
        sum_chk (sum, c);
        sum_chk (sum, d);
        sum_chk (sum, e);
        sum_chk (sum, f);
        sum_chk (sum, g);

        if (i < 3) {
            char nambuf[5];
            sprintf (nambuf, "x%d", i);
            nerrs += bgcheck (sfin, sfout, x, nambuf);
            nambuf[0] = 'y';
            nerrs += bgcheck (sfin, sfout, y, nambuf);
            nambuf[0] = 'a';
            nerrs += bgcheck (sfin, sfout, a, nambuf);
            nambuf[0] = 'b';
            nerrs += bgcheck (sfin, sfout, b, nambuf);
            nambuf[0] = 'c';
            nerrs += bgcheck (sfin, sfout, c, nambuf);
            nambuf[0] = 'd';
            nerrs += bgcheck (sfin, sfout, d, nambuf);
            nambuf[0] = 'e';
            nerrs += bgcheck (sfin, sfout, e, nambuf);
            nambuf[0] = 'f';
            nerrs += bgcheck (sfin, sfout, f, nambuf);
            nambuf[0] = 'g';
            nerrs += bgcheck (sfin, sfout, g, nambuf);
        }
    }
    nerrs += bgcheck (sfin, sfout, *sum, "end sum");

    return nerrs;
}

static void sum_chk (CCbigguy *sum, CCbigguy x)
{
    CCbigguy y;

    if (CCbigguy_cmp (x, CCbigguy_ZERO) > 0 &&
        CCbigguy_cmp (*sum, CCbigguy_ZERO) > 0) {
        y = CCbigguy_MAXBIGGUY;
        CCbigguy_sub (&y, x);
        if (CCbigguy_cmp (*sum, y) >= 0) {
            CCbigguy_sub (sum, CCbigguy_MAXBIGGUY);
        }
    } else if (CCbigguy_cmp (x, CCbigguy_ZERO) < 0 &&
        CCbigguy_cmp (*sum, CCbigguy_ZERO) < 0) {
        y = CCbigguy_MINBIGGUY;
        CCbigguy_sub (&y, x);
        if (CCbigguy_cmp (*sum, y) <= 0) {
            CCbigguy_sub (sum, CCbigguy_MINBIGGUY);
        }
    }
    CCbigguy_add (sum, x);
}

static int bgcheck (CC_SFILE *sfin, CC_SFILE *sfout, CCbigguy x,
        const char *name)
{
    CCbigguy v;

    CCbigguy_sread (sfin, &v);
    CCbigguy_swrite (sfout, x);
    if (CCbigguy_cmp(v,x)) {
        fprintf (stderr, "%s misread (%.15f != %.15f)\n", name,
                 CCbigguy_bigguytod(v), CCbigguy_bigguytod(x));
        return 1;
    } else {
        return 0;
    }
}

static CCbigguy bgrand (double *delta, int maxv, CCrandstate *rstate)
{
    CCbigguy x, y;
    double d16 = 1.0 / (1 << 16);
    double d32 = d16 * d16;
    double dx;

    dx = d32 * (CCutil_lprand (rstate) & 0xffff);
    x = CCbigguy_dtobigguy (dx);
    MOREDELTA(*delta, x, dx, "bgrand 1");

    dx = d16 * (CCutil_lprand (rstate) & 0xffff);
    y = CCbigguy_dtobigguy(dx);
    MOREDELTA(*delta, y, dx, "bgrand 2");
    CCbigguy_add (&x, y);

    dx = 1.0 * (CCutil_lprand (rstate) % maxv);
    y = CCbigguy_dtobigguy(dx);
    MOREDELTA(*delta, y, dx, "bgrand 3");
    CCbigguy_add (&x, y);

    if (CCutil_lprand (rstate) & 1) {
        return x;
    } else {
        dx = CCbigguy_bigguytod(x);
        y = CCbigguy_ZERO;
        CCbigguy_sub (&y, x);
        MOREDELTA(*delta, y, -dx, "bgrand 4");
        return y;
    }
}

static int intrand (int maxi, CCrandstate *rstate)
{
    int s = CCutil_lprand (rstate) % maxi;

    if (CCutil_lprand (rstate) % 1) s = -s;

    return s;
}
