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
/*                KD-TREE BASED GEOMETRIC 2-OPT ROUTINE                     */
/*                                                                          */
/*  (Based on Jon Bentley's paper "Fast algorithms for geometric            */
/*   traveling salesman problems", ORSA Journal on Computing.)              */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook (mainly from Jon       */
/*               Bentley's paper)                                           */
/*  Date: February 16, 1995  (cofeb16)                                      */
/*        Modified: October 11, 1995 (Bico) - removed globals.              */
/*                                                                          */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCkdtree_twoopt_tour (CCkdtree *kt, int ncount, CCdatagroup *dat,   */
/*      int *incycle, int *outcycle, double *val,                           */
/*      int run_two_and_a_half_opt, int silent, CCrandstate *rstate)        */
/*    RETURNS a 2-opted cycle (well, approximately 2-opted)                 */
/*      -kt can be NULL.                                                    */
/*      -Does not use node weights.                                         */
/*      -incycle is the starting cycle.                                     */
/*      -If outcycle is not NULL, then it should point to an array of       */
/*       length at least ncount (allocated by the calling routine). The     */
/*       final tour will be placed in this array.                           */
/*      -The length of the tour is returned in val.                         */
/*      -If in_run_two_and_a_half_opt is nonzero,  then some limited        */
/*       3-swapping is performded.                                          */
/*      -silent (if nonzero then very little info will be printed)          */
/*                                                                          */
/*  int CCkdtree_3opt_tour (CCkdtree *kt, int ncount, CCdatagroup *dat,     */
/*      int *incycle, int *outcycle, double *val, silent,                   */
/*      CCrandstate *rstate)                                                */
/*    RETURNS an approximately 3-opted tour.                                */
/*      -kt can be NULL.                                                    */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/*  NOTES:                                                                  */
/*     Uses around 13n bytes of memory for an n node problem (plus the      */
/*  the memory for the CCkdtree). Returns 0 if successful, and nonzero if   */
/*  it failed (usually due to running out of memory.                        */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "kdtree.h"
#include "macrorus.h"

#define Edgelen(a, b)  CCutil_dat_edgelen ((a), (b), dat)
#define ADD_TO_ACTIVE_QUEUE(n, ip, q) {                                    \
    if (!(q)->active[(n)]) {                                               \
        (q)->active[(n)] = 1;                                              \
        (ip) = intptralloc (&(q)->intptr_world);                           \
        (ip)->this = (n);                                                  \
        (ip)->next = (intptr *) NULL;                                      \
        if ((q)->bottom)                                                   \
            (q)->bottom->next = (ip);                                      \
        else                                                               \
            (q)->queue = (ip);                                             \
        (q)->bottom = (ip);                                                \
    }                                                                      \
}

typedef struct intptr {
    int this;
    struct intptr *next;
} intptr;

typedef struct optqueue {
    char *active;
    intptr *queue;
    intptr *bottom;
    CCptrworld intptr_world;
} optqueue;

typedef struct flipper {
    int reversed;
    int cycle_size;
    int short_size;
    int *cyc;
    int *cyc_inv;
} flipper;

typedef struct twoopt_param {
    CCdatagroup *dat;
    optqueue  *activequeue;
    flipper   *flipstuff;
    int node_b;
} twoopt_param;

typedef struct twoopt_and_a_half_param {
    CCdatagroup *dat;
    optqueue  *activequeue;
    flipper   *flipstuff;
    int node_b;
} twoopt_and_a_half_param;

typedef struct threeopt_param {
    CCdatagroup *dat;
    optqueue  *activequeue;
    flipper   *flipstuff;
    int node_a;
    int node_b;
    int node_c;
    int ab_dist;
    int ac_dist;
    int cd_dist;
    CCkdtree *thetree;
} threeopt_param;

CC_PTRWORLD_ROUTINES (intptr, intptralloc, intptr_bulk_alloc, intptrfree)
CC_PTRWORLD_LEAKS_ROUTINE (intptr, intptr_check_leaks, this, int)

static void
    flipper_cycle (flipper *f, int *x),
    flipper_finish (flipper *f),
    flipper_flip (flipper *f, int x, int y),
    randcycle (int ncount, int *cyc, CCrandstate *rstate),
    twoopt_free_world (CCptrworld *intptr_world);
static int
    run_two_opt (int ncount, int silent, int *tour,
          int run_two_and_a_half_swap, int *neighbor, CCkdtree *thetree,
          CCdatagroup *dat, optqueue *activequeue),
    run_3_opt (int ncount, int silent, int *tour, int *neighbor,
           CCkdtree *thetree, CCdatagroup *dat, optqueue *activequeue),
    swap_from_node (int a, int run_two_an_a_half_swap, int *neighbor,
                    CCkdtree *thetree, CCdatagroup *dat, optqueue *activequeue,
                    flipper *f),
    try_two_swap (int a, int c, void *pass_param),
    try_two_and_a_half_swap (int a, int c, void *pass_param),
    three_swap_from_node (int a, int *neighbor, CCkdtree *thetree,
           CCdatagroup *dat, optqueue *activequeue, flipper *f),
    try_three_swap (int a, int c, void *pass_param),
    try_three_swap_2 (int a, int c, void *pass_param),
    flipper_init (flipper *f, int *incyc, int n),
    flipper_next (flipper *f, int x),
    flipper_prev (flipper *f, int x),
    flipper_sequence (flipper *f, int x, int y, int z);
static double
    cycle_length (int *cyc, int ncount, CCdatagroup *dat);


int CCkdtree_twoopt_tour (CCkdtree *kt, int ncount, CCdatagroup *dat,
        int *incycle, int *outcycle, double *val, int run_two_and_a_half_opt,
        int silent, CCrandstate *rstate)
{
    CCkdtree localkt;
    int i;
    int rval = 0;
    double szeit;
    int newtree = 0;
    intptr *ip;
    int *tour = (int *) NULL;
    int *neighbor = (int *) NULL;
    CCkdtree *thetree = (CCkdtree *) NULL;
    optqueue activequeue;

    *val = 0.0;
    CCptrworld_init (&activequeue.intptr_world);

    if (kt == (CCkdtree *) NULL) {
        if (CCkdtree_build (&localkt, ncount, dat, (double *) NULL, rstate)) {
            fprintf (stderr, "Unable to build CCkdtree\n");
            return 1;
        }
        thetree = &localkt;
        newtree = 1;
    } else {
        thetree = kt;
    }

    activequeue.active = (char *) NULL;
    activequeue.queue = (intptr *) NULL;
    activequeue.bottom = (intptr *) NULL;

    if (!silent) {
        printf ("Find 2-opt Tour starting with tour of length %.2f\n",
                   cycle_length (incycle, ncount, dat));
        fflush (stdout);
    }

    szeit = CCutil_zeit ();
    tour = CC_SAFE_MALLOC (ncount, int);
    if (!tour) {
        rval = 1;
        goto CLEANUP;
    }
    neighbor = CC_SAFE_MALLOC (ncount, int);
    if (!neighbor) {
        rval = 1;
        goto CLEANUP;
    }
    activequeue.active = CC_SAFE_MALLOC (ncount, char);
    if (!activequeue.active) {
        rval = 1;
        goto CLEANUP;
    }

    for (i = 0; i < ncount; i++)
        activequeue.active[i] = 0;

    randcycle (ncount, tour, rstate); /* init active_queue with random order */
    for (i = 0; i < ncount; i++) {
        ADD_TO_ACTIVE_QUEUE (tour[i], ip, &activequeue);
    }

    for (i = 0; i < ncount; i++) {
        tour[i] = incycle[i];
    }
    run_two_opt (ncount, silent, tour, run_two_and_a_half_opt,
                 neighbor, thetree, dat, &activequeue);

    *val = cycle_length (tour, ncount, dat);
    if (!silent) {
        printf ("Length of Two-opt Cycle: %.2f\n", *val);
        fflush (stdout);
    }
    if (outcycle != (int *) NULL) {
        for (i = 0; i < ncount; i++)
            outcycle[i] = tour[i];
    }
    if (!silent) {
        printf ("Running time for Two Opt: %.2f\n", CCutil_zeit () - szeit);
        fflush (stdout);
    }

CLEANUP:

    if (newtree)
        CCkdtree_free (&localkt);
    if (tour)
        CC_FREE (tour, int);
    if (neighbor)
        CC_FREE (neighbor, int);
    if (activequeue.active)
        CC_FREE (activequeue.active, char);
    twoopt_free_world (&activequeue.intptr_world);
    return rval;
}

static int run_two_opt (int ncount, int silent, int *tour,
        int run_two_and_a_half_opt, int *neighbor, CCkdtree *thetree,
        CCdatagroup *dat, optqueue *activequeue)
{
    int i, hit;
    int start;
    intptr *ip;
    flipper f;

    for (i = 0; i < ncount; i++) {
        CCkdtree_delete (thetree, i);
        neighbor[i] = CCkdtree_node_nearest (thetree, i, dat, (double *) NULL);
        CCkdtree_undelete (thetree, i);
    }

    flipper_init (&f, tour, ncount);
    hit = 0;
    while (activequeue->queue) {
        ip = activequeue->queue;
        start = ip->this;
        activequeue->queue = ip->next;
        if (ip == activequeue->bottom)
            activequeue->bottom = (intptr *) NULL;
        intptrfree (&activequeue->intptr_world, ip);
        activequeue->active[start] = 0;
        if (swap_from_node (start, run_two_and_a_half_opt, neighbor,
                            thetree, dat, activequeue, &f)) {
            hit++;
            if (!silent && hit % 1000 == 0) {
                printf (".");
                fflush (stdout);
                if (hit % 50000 == 0) {
                    flipper_cycle (&f, tour);
                    printf ("\nCurrent length: %.2f\n",
                            cycle_length (tour, ncount, dat));
                    fflush (stdout);
                }
            }
        }
    }

    if (!silent) {
        printf ("\nMade %d swaps\n", hit);
        fflush (stdout);
    }

    flipper_cycle (&f, tour);
    flipper_finish (&f);
    return 0;
}

static void twoopt_free_world (CCptrworld *intptr_world)
{
    int total, onlist;

    if (intptr_check_leaks (intptr_world, &total, &onlist)) {
        fprintf (stderr, "WARNING: %d outstanding intptrs\n", total - onlist);
    }
    CCptrworld_delete (intptr_world);
}


int CCkdtree_3opt_tour (CCkdtree *kt, int ncount, CCdatagroup *dat,
        int *incycle, int *outcycle, double *val, int silent,
        CCrandstate *rstate)
{
    CCkdtree localkt;
    int i;
    int rval = 0;
    double szeit;
    int newtree = 0;
    intptr *ip;
    int *tour = (int *) NULL;
    int *neighbor = (int *) NULL;
    CCkdtree *thetree = (CCkdtree *) NULL;
    optqueue activequeue;

    *val = 0.0;
    CCptrworld_init (&activequeue.intptr_world);

    if (kt == (CCkdtree *) NULL) {
        if (CCkdtree_build (&localkt, ncount, dat, (double *) NULL, rstate)) {
            fprintf (stderr, "Unable to build CCkdtree\n");
            return 1;
        }
        thetree = &localkt;
        newtree = 1;
    } else {
        thetree = kt;
    }

    activequeue.active = (char *) NULL;
    activequeue.queue = (intptr *) NULL;
    activequeue.bottom = (intptr *) NULL;

    if (!silent) {
        printf ("Find 3-opt Tour starting with tour of length %.2f\n",
                       cycle_length (incycle, ncount, dat));
        fflush (stdout);
    }

    szeit = CCutil_zeit ();
    tour = CC_SAFE_MALLOC (ncount, int);
    if (!tour) {
        rval = 1;
        goto CLEANUP;
    }
    activequeue.active = CC_SAFE_MALLOC (ncount, char);
    if (!activequeue.active) {
        rval = 1;
        goto CLEANUP;
    }
    neighbor = CC_SAFE_MALLOC (ncount, int);
    if (!neighbor) {
        rval = 1;
        goto CLEANUP;
    }

    for (i = 0; i < ncount; i++)
        activequeue.active[i] = 0;

    randcycle (ncount, tour, rstate); /* init active_queue with random order */
    for (i = 0; i < ncount; i++) {
        ADD_TO_ACTIVE_QUEUE (tour[i], ip, &activequeue);
    }

    for (i = 0; i < ncount; i++) {
        tour[i] = incycle[i];
    }
    run_3_opt (ncount, silent, tour, neighbor, thetree, dat,
               &activequeue);

    *val = cycle_length (tour, ncount, dat);
    if (!silent) {
        printf ("Length of 3-opt Cycle: %.2f\n", *val);
        fflush (stdout);
    }
    if (outcycle != (int *) NULL) {
        for (i = 0; i < ncount; i++)
            outcycle[i] = tour[i];
    }
    if (!silent) {
        printf ("Running time for 3-Opt: %.2f\n", CCutil_zeit () - szeit);
        fflush (stdout);
    }

CLEANUP:

    if (newtree)
        CCkdtree_free (&localkt);
    if (tour)
        CC_FREE (tour, int);
    if (activequeue.active)
        CC_FREE (activequeue.active, char);
    if (neighbor)
        CC_FREE (neighbor, int);
    twoopt_free_world (&activequeue.intptr_world);
    return rval;
}

static int run_3_opt (int ncount, int silent, int *tour,
                      int *neighbor, CCkdtree *thetree, CCdatagroup *dat,
                      optqueue *activequeue)
{
    int i, hit;
    int start;
    intptr *ip;
    flipper f;

    for (i = 0; i < ncount; i++) {
        CCkdtree_delete (thetree, i);
        neighbor[i] = CCkdtree_node_nearest (thetree, i, dat, (double *) NULL);
        CCkdtree_undelete (thetree, i);
    }

    flipper_init (&f, tour, ncount);
    hit = 0;
    while (activequeue->queue) {
        ip = activequeue->queue;
        start = ip->this;
        activequeue->queue = ip->next;
        if (ip == activequeue->bottom)
            activequeue->bottom = (intptr *) NULL;
        intptrfree (&activequeue->intptr_world, ip);
        activequeue->active[start] = 0;
        if (three_swap_from_node (start, neighbor, thetree, dat, activequeue,
                                  &f)) {
            hit++;
            if (!silent && hit % 1000 == 0) {
                printf (".");
                fflush (stdout);
                if (hit % 50000 == 0) {
                    flipper_cycle (&f, tour);
                    printf ("\nCurrent length: %.2f\n",
                            cycle_length (tour, ncount, dat));
                    fflush (stdout);
                }
            }
        }
    }

    if (!silent) {
        printf ("\nMade %d swaps\n", hit);
        fflush (stdout);
    }

    flipper_cycle (&f, tour);
    flipper_finish (&f);
    return 0;
}

static int swap_from_node (int a, int run_two_and_a_half_opt,
                           int *neighbor, CCkdtree *thetree, CCdatagroup *dat,
                           optqueue *activequeue, flipper *f)
{
    twoopt_param param2;
    twoopt_param param2h;
    int b;

    b = flipper_next (f, a);
    if (neighbor[a] ==b) {
        return 0;
    } else {
        int dist = Edgelen (a, b);
        if (run_two_and_a_half_opt) {
            param2h.dat = dat;
            param2h.activequeue = activequeue;
            param2h.node_b = b;
            param2h.flipstuff = f;
            return CCkdtree_fixed_radius_nearest (thetree, dat,
                    (double *) NULL, a, (double) dist, try_two_and_a_half_swap,
                    (void *) (&param2h));
        } else {
            param2.dat = dat;
            param2.activequeue = activequeue;
            param2.node_b = b;
            param2.flipstuff = f;
            return CCkdtree_fixed_radius_nearest (thetree, dat,
                    (double *) NULL, a, (double) dist, try_two_swap,
                    (void *) (&param2));
        }
    }
}

static int try_two_swap (int a, int c, void *pass_param)
{
    int d;
    int b = ((twoopt_param *) pass_param)->node_b;
    CCdatagroup *dat = ((twoopt_param *) pass_param)->dat;
    optqueue *activequeue = ((twoopt_param *) pass_param)->activequeue;
    flipper *flipstuff = ((twoopt_param *) pass_param)->flipstuff;
    intptr *ip;

    d = flipper_next (flipstuff, c);

    if (d == a)
        return 0;

    if (Edgelen (a, b) + Edgelen (c, d) >
        Edgelen (a, c) + Edgelen (b, d)) {
        flipper_flip (flipstuff, d, a);
        ADD_TO_ACTIVE_QUEUE (a, ip, activequeue);
        ADD_TO_ACTIVE_QUEUE (b, ip, activequeue);
        ADD_TO_ACTIVE_QUEUE (c, ip, activequeue);
        ADD_TO_ACTIVE_QUEUE (d, ip, activequeue);
        return 1;
    } else
        return 0;
}

static int try_two_and_a_half_swap (int a, int c, void *pass_param)
{
    int d, e;
    int b = ((twoopt_and_a_half_param *) pass_param)->node_b;
    CCdatagroup *dat = ((twoopt_param *) pass_param)->dat;
    optqueue *activequeue = ((twoopt_param *) pass_param)->activequeue;
    flipper *flipstuff = ((twoopt_param *) pass_param)->flipstuff;
    intptr *ip;

    d = flipper_next (flipstuff, c);

    if (d == a)
        return 0;

    if (Edgelen (a, b) + Edgelen (c, d) >
        Edgelen (a, c) + Edgelen (b, d)) {
        flipper_flip (flipstuff, d, a);
        ADD_TO_ACTIVE_QUEUE (a, ip, activequeue);
        ADD_TO_ACTIVE_QUEUE (b, ip, activequeue);
        ADD_TO_ACTIVE_QUEUE (c, ip, activequeue);
        ADD_TO_ACTIVE_QUEUE (d, ip, activequeue);
        return 1;
    } else {
        e = flipper_prev (flipstuff, c);
        if (e != b &&
            Edgelen (a, b) + Edgelen (c, d) + Edgelen (c, e) >
            Edgelen (a, c) + Edgelen (b, c) + Edgelen (d, e)) {
            flipper_flip (flipstuff, b, e);
            flipper_flip (flipstuff, e, c);
            ADD_TO_ACTIVE_QUEUE (a, ip, activequeue);
            ADD_TO_ACTIVE_QUEUE (b, ip, activequeue);
            ADD_TO_ACTIVE_QUEUE (c, ip, activequeue);
            ADD_TO_ACTIVE_QUEUE (d, ip, activequeue);
            ADD_TO_ACTIVE_QUEUE (e, ip, activequeue);
            return 1;
        } else {
            return 0;
        }
    }
}

static int three_swap_from_node (int a, int *neighbor, CCkdtree *thetree,
        CCdatagroup *dat, optqueue *activequeue, flipper *f)
{
    threeopt_param param3;
    int b, ab;

    b = flipper_next (f, a);

    if (neighbor[a] == b) {
        return 0;
    } else {
        ab = Edgelen (a, b);
        param3.dat = dat;
        param3.activequeue = activequeue;
        param3.node_a = a;
        param3.node_b = b;
        param3.thetree = thetree;
        param3.ab_dist = ab;
        param3.flipstuff = f;
        return CCkdtree_fixed_radius_nearest (thetree, dat, (double *) NULL,
                     a, (double) ab, try_three_swap, (void *) (&param3));
    }
}

static int try_three_swap (int a, int c, void *pass_param)
{
    int d, b, ac, cd;
    threeopt_param *p = (threeopt_param *) pass_param;
    CCdatagroup *dat;
    optqueue *activequeue;
    intptr *ip;

    dat = p->dat;
    activequeue = p->activequeue;
    b = p->node_b;
    d = flipper_next (p->flipstuff, c);

    ac = Edgelen (a, c);
    cd = Edgelen (c, d);
    if (p->ab_dist + cd > ac + Edgelen (b, d)) {
        flipper_flip (p->flipstuff, d, a);
        ADD_TO_ACTIVE_QUEUE (a, ip, activequeue);
        ADD_TO_ACTIVE_QUEUE (b, ip, activequeue);
        ADD_TO_ACTIVE_QUEUE (c, ip, activequeue);
        ADD_TO_ACTIVE_QUEUE (d, ip, activequeue);
        return 1;
    } else {
        int e = flipper_prev (p->flipstuff, c);
        if (e != b && d != a &&
            p->ab_dist + cd + Edgelen (c, e) >
            ac + Edgelen (b, c) + Edgelen (d, e)) {
            flipper_flip (p->flipstuff, b, e);
            flipper_flip (p->flipstuff, e, c);
            ADD_TO_ACTIVE_QUEUE (a, ip, activequeue);
            ADD_TO_ACTIVE_QUEUE (b, ip, activequeue);
            ADD_TO_ACTIVE_QUEUE (c, ip, activequeue);
            ADD_TO_ACTIVE_QUEUE (d, ip, activequeue);
            ADD_TO_ACTIVE_QUEUE (e, ip, activequeue);
            return 1;
        } else {
            p->node_c = c;
            p->ac_dist = ac;
            p->cd_dist = cd;

            return CCkdtree_fixed_radius_nearest (p->thetree, dat,
              (double *) NULL, d, (double) (p->ab_dist + cd - ac),
              try_three_swap_2, pass_param);
        }
    }
}

static int try_three_swap_2 (int d, int e, void *pass_param)
{
    threeopt_param *p = (threeopt_param *) pass_param;
    CCdatagroup *dat;
    optqueue *activequeue;
    int a, b, c, f;
    intptr *ip;

    a = p->node_a;
    b = p->node_b;
    c = p->node_c;
    dat = p->dat;
    activequeue = p->activequeue;

    if (e == c)
        return 0;
    else if (flipper_sequence (p->flipstuff, b, e, c))
        f = flipper_next (p->flipstuff, e);
    else
        f = flipper_prev (p->flipstuff, e);

    if (p->ab_dist + p->cd_dist + Edgelen (e, f) >
        p->ac_dist + Edgelen (b, f) + Edgelen (d, e)) {
        flipper_flip (p->flipstuff, d, a);
        flipper_flip (p->flipstuff, b, e);
        ADD_TO_ACTIVE_QUEUE (a, ip, activequeue);
        ADD_TO_ACTIVE_QUEUE (b, ip, activequeue);
        ADD_TO_ACTIVE_QUEUE (c, ip, activequeue);
        ADD_TO_ACTIVE_QUEUE (d, ip, activequeue);
        ADD_TO_ACTIVE_QUEUE (e, ip, activequeue);
        ADD_TO_ACTIVE_QUEUE (f, ip, activequeue);
        return 1;
    } else {
        return 0;
    }
}

static double cycle_length (int *cyc, int ncount, CCdatagroup *dat)
{
    int i;
    double len;

    for (len = 0.0, i = 1; i < ncount; i++)
        len += CCutil_dat_edgelen (cyc[i - 1], cyc[i], dat);
    len += CCutil_dat_edgelen (cyc[0], cyc[ncount - 1], dat);

    return len;
}

static void randcycle (int ncount, int *cyc, CCrandstate *rstate)
{
    int i, k, temp;

    for (i = 0; i < ncount; i++)
        cyc[i] = i;

    for (i = ncount; i > 1; i--) {
        k = CCutil_lprand (rstate) % i;
        CC_SWAP (cyc[i - 1], cyc[k], temp);
    }
}


/**********   modified version of flipper3.c **********/

/*
   flipper routines:

   flipper_init (flipper *f, int n)
      initializes flipper3 to an initial * cycle from 0 to n-1.
   flipper_cycle (flipper *f, int *p)
      places the current * cycle in p.
   flipper_finish (flipper *f, )
      frees up space * allocated by flipper.
   flipper_next (flipper *f, int x)
      returns the * successor to x in the current cycle.
   flipper_prev (flipper *f, int x)
      returns the predecessor of x in the current cycle.
   flipper_flip (flipper *f, int * x, int y)
      flips the portion of the cycle from x to y (inclusive).
   flipper_sequence (flipper *f, int x, int y, int z)
      returns 1 if x, y, z occur as an increasing
      subsequence of the cycle, 0 otherwise.

   The basic implementation - uses a raw cycle and its inverse, but flips
   the short side.
*/


static int flipper_init (flipper *f, int *incyc, int n)
{
    int i;

    f->cyc_inv = (int *) NULL;
    f->cyc = CC_SAFE_MALLOC (n, int);
    if (!f->cyc)
        return 1;
    f->cyc_inv = CC_SAFE_MALLOC (n, int);
    if (!f->cyc_inv) {
        CC_FREE (f->cyc, int);
        return 1;
    }

    for (i = 0; i < n; i++) {
        f->cyc[i] = incyc[i];
        f->cyc_inv[incyc[i]] = i;
    }

    f->cycle_size = n;
    f->short_size = n / 2;
    f->reversed = 0;
    return 0;
}

static void flipper_cycle (flipper *f, int *x)
{
    int *p;

    if (f->reversed) {
        p = f->cyc + f->cycle_size;
        while (p > f->cyc) {
            *x++ = *--p;
        }
    } else {
        p = f->cyc + f->cycle_size;
        x += f->cycle_size;
        while (p > f->cyc) {
            *--x = *--p;
        }
    }
}

static void flipper_finish (flipper *f)
{
    if (f->cyc)
        CC_FREE (f->cyc, int);
    if (f->cyc_inv)
        CC_FREE (f->cyc_inv, int);
    f->cycle_size = 0;
    f->short_size = 0;
    f->reversed = 0;
}

static int flipper_next (flipper *f, int x)
{
    int y;

    if (f->reversed) {
        y = f->cyc_inv[x] - 1;
        return (y >= 0) ? f->cyc[y] : f->cyc[f->cycle_size - 1];
    } else {
        y = f->cyc_inv[x] + 1;
        return (y < f->cycle_size) ? f->cyc[y] : f->cyc[0];
    }
}

static int flipper_prev (flipper *f, int x)
{
    int y;

    if (f->reversed) {
        y = f->cyc_inv[x] + 1;
        return (y < f->cycle_size) ? f->cyc[y] : f->cyc[0];
    } else {
        y = f->cyc_inv[x] - 1;
        return (y >= 0) ? f->cyc[y] : f->cyc[f->cycle_size - 1];
    }
}

static void flipper_flip (flipper *f, int x, int y)
{
    int xloc = f->cyc_inv[x];
    int yloc = f->cyc_inv[y];
    int temp;
    int gap;

    if (f->reversed) {
        CC_SWAP (xloc, yloc, temp);
    }
    gap = yloc - xloc;
    if (gap < 0)
        gap += f->cycle_size;
    if (gap > f->short_size) {
        CC_SWAP (xloc, yloc, temp);
        f->reversed ^= 1;
        xloc++;
        if (xloc >= f->cycle_size)
            xloc = 0;
        yloc--;
        if (yloc < 0)
            yloc = f->cycle_size - 1;
        gap = f->cycle_size - gap - 2;
    }

    if (xloc > yloc) {
        gap++;
        gap /= 2;
        for (; gap; gap--) {
            x = f->cyc[xloc];
            y = f->cyc[yloc];
            f->cyc[xloc] = y;
            f->cyc[yloc] = x;
            f->cyc_inv[x] = yloc--;
            f->cyc_inv[y] = xloc++;
            if (xloc >= f->cycle_size)
                xloc = 0;
            if (yloc < 0)
                yloc = f->cycle_size - 1;
        }
    } else {
        gap++;
        gap /= 2;
        for (; gap; gap--) {
            x = f->cyc[xloc];
            y = f->cyc[yloc];
            f->cyc[xloc] = y;
            f->cyc[yloc] = x;
            f->cyc_inv[x] = yloc--;
            f->cyc_inv[y] = xloc++;
        }
    }

/*
    while (gap > 0) {
        CC_SWAP (f->cyc[xloc], f->cyc[yloc], temp);
        CC_SWAP (f->cyc_inv[cyc[xloc]], f->cyc_inv[cyc[yloc]], temp);
        xloc++;
        if (xloc >= f->cycle_size)
            xloc = 0;
        yloc--;
        if (yloc < 0)
            yloc = f->cycle_size - 1;
        gap -= 2;
    }
*/
}

static int flipper_sequence (flipper *f, int x, int y, int z)
{
    int xloc = f->cyc_inv[x];
    int yloc = f->cyc_inv[y];
    int zloc = f->cyc_inv[z];

    if (f->reversed) {
        if (xloc >= yloc)
            return yloc >= zloc || zloc >= xloc;
        else
            return yloc >= zloc && zloc >= xloc;
    } else {
        if (xloc <= yloc)
            return yloc <= zloc || zloc <= xloc;
        else
            return yloc <= zloc && zloc <= xloc;
    }
}


