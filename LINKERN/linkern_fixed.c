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
/*           CHAINED LIN-KERNIGHAN  -  Modified for Fixed Edges             */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: March 22, 1995                                                    */
/*        May 1, 1998 (bico)                                                */
/*        May 6, 2003 (bico)                                                */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CClinkern_fixed (int ncount, CCdatagroup *dat, int ecount,          */
/*      int *elist, int nkicks, int *incycle, int *outcycle,                */
/*      double *val, int fcount, int *flist, int silent,                    */
/*      CCrandstate *rstate)                                                */
/*    RUNS Chained Lin-Kernighan, with fixed edges given in node-node       */
/*      form flist                                                          */
/*    -ncount (the number of nodes int the graph)                           */
/*    -dat (coordinate dat)                                                 */
/*    -ecount (the number of good edges - should not be 0)                  */
/*    -elist (the good edges in end1 end2 format)                           */
/*    -nkicks (the number of 4-swap kicks)                                  */
/*    -incycle  (a starting cycle in node node node format)                 */
/*    -outcycle (returns the cycle - can be NULL)                           */
/*    -silent (if nonzero, then very little info will be printed)           */
/*                                                                          */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "linkern.h"
#include "kdtree.h"
#include "util.h"
#include "macrorus.h"

#define MAXDEPTH       25   /* Shouldn't really be less than 2.             */
#define KICK_MAXDEPTH  50
#define LONG_KICKER  
#define ACCEPT_TIES 

#define USE_LESS_OR_EQUAL 
#define LATE_DEPTH 10      /* Should be less than MAXDEPTH                 */

#define MARK_LEVEL 10       /* Number of tour neighbors after 4-swap kick   */
#define BACKTRACK   4
#define MAX_BACK   12       /* Upper bound on the XXX_count entries         */
static const int backtrack_count[BACKTRACK] = {4, 3, 3, 2};
static const int weird_backtrack_count[3] = {4, 3, 3};

#define BIGINT 2000000000
#define Edgelen(n1, n2, D)  dist (n1, n2, D)
/*
#define Edgelen(n1, n2, D)  CCutil_dat_edgelen (n1, n2, D->dat) 
*/

#define FLIP(aprev, a, b, bnext, f, x) {                                   \
    CClinkern_flipper_flip ((x),(a), (b));                                 \
    (f)->stack[(f)->counter].first = (a);                                  \
    (f)->stack[(f)->counter++].last = (b);                                 \
}

#define UNFLIP(aprev, a, b, bnext, f, x) {                                 \
    CClinkern_flipper_flip ((x), (b), (a));                                \
    (f)->counter--;                                                        \
}

#define MARK(xn, xQ, xF, xD, xG, xW)  turn ((xn), (xQ), (xW))

#define markedge_add(n1, n2, E)    E->add_edges[n1 ^ n2] = 1
#define markedge_del(n1, n2, E)    E->del_edges[n1 ^ n2] = 1
#define unmarkedge_add(n1, n2, E)  E->add_edges[n1 ^ n2] = 0
#define unmarkedge_del(n1, n2, E)  E->del_edges[n1 ^ n2] = 0
#define is_it_added(n1, n2, E)     E->add_edges[n1 ^ n2]
#define is_it_deleted(n1, n2, E)   E->del_edges[n1 ^ n2]

#define Fixededge(n1, n2) (G->fixlist && fixed_edge (n1, n2, G->fixlist))

typedef struct edge {
    int other;
    int weight;
} edge;

typedef struct edgelook {
    struct edgelook *next;
    int other;
    int diff;
    int over;
    int seq;
    int side;
    int mm;
} edgelook;

typedef struct intptr {
    int this;
    struct intptr *next;
} intptr;

typedef struct flippair {
    int firstprev;
    int first;
    int last;
    int lastnext;
} flippair;

typedef struct flipstack {
    flippair *stack;
    int counter;
    int max;
} flipstack;

typedef struct graph {
    edge **goodlist;
    edge *edgespace;
    int  *degree;
    int  *weirdmark;
    int   weirdmagic;
    int   ncount;
    intptr **fixlist;
    CCrandstate *rstate;
} graph;

typedef struct distobj {
    CCdatagroup *dat;
    int       *cacheval;
    int       *cacheind;
    int        cacheM;
} distobj;

typedef struct adddel {
    char *add_edges;
    char *del_edges;
} adddel;

typedef struct aqueue {
    char *active;
    intptr *active_queue;
    intptr *bottom_active_queue;
    CCdheap *h;
} aqueue;

static void
   lin_kernighan (graph *G, distobj *D, adddel *E, aqueue *Q, CClk_flipper *F,
       double *val, int *win_cycle, flipstack *w, flipstack *fstack,
       CCptrworld *intptr_world, CCptrworld *edgelook_world),
   look_ahead_noback (graph *G, distobj *D, adddel *E, CClk_flipper *F,
       int first, int last, int gain, edgelook *winner),
   turn (int n, aqueue *Q, CCptrworld *intptr_world),
   kickturn (int n, aqueue *Q, distobj *D, graph *G, CClk_flipper *F,
        CCptrworld *intptr_world),
   bigturn (graph *G, int n, int tonext, aqueue *Q, CClk_flipper *F,
        distobj *D, CCptrworld *intptr_world),
   first_kicker (graph *G, distobj *D, CClk_flipper *F, int *t1, int *t2),
   find_walk_four (graph *G, distobj *D, CClk_flipper *F, int *t1,
        int *t2, int *t3, int *t4, int *t5, int *t6, int *t7, int *t8),
   randcycle (int ncount, int *cyc, CCrandstate *rstate),
   insertedge (graph *G, int n1, int n2, int w),
   initgraph (graph *G),
   freegraph (graph *G),
   init_adddel (adddel *E),
   free_adddel (adddel *E),
   init_aqueue (aqueue *Q),
   free_aqueue (aqueue *Q, CCptrworld *intptr_world),
   add_to_active_queue (int n, aqueue *Q, CCptrworld *intptr_world),
   init_distobj (distobj *D),
   free_distobj (distobj *D),
   linkern_free_world (CCptrworld *intptr_world, CCptrworld *edgelook_world),
   free_flipstack (flipstack *f);

static int
   buildgraph (graph *G, int ncount, int ecount, int *elist, distobj *D),
   repeated_lin_kernighan (graph *G, distobj *D, int *cyc,
       int repeatcount, double *val, int silent, CCptrworld *intptr_world,
       CCptrworld *edgelook_world),
   weird_second_step (graph *G, distobj *D, adddel *E, aqueue *Q,
       CClk_flipper *F, int gain, int t1, int t2, flipstack *fstack,
       CCptrworld *intptr_world, CCptrworld *edgelook_world),
   step (graph *G, distobj *D, adddel *E, aqueue *Q, CClk_flipper *F,
       int level, int gain, int *Gstar, int first, int last, flipstack *fstack,
       CCptrworld *intptr_world, CCptrworld *edgelook_world),
   step_noback (graph *G, distobj *D, adddel *E, aqueue *Q, CClk_flipper *F,
       int level, int gain, int *Gstar, int first, int last,
       flipstack *fstack, CCptrworld *intptr_world),
   random_four_swap (graph *G, distobj *D, aqueue *Q, CClk_flipper *F,
       int *delta, flipstack *win, flipstack *fstack, CCptrworld *intptr_world),
   build_adddel (adddel *E, int ncount),
   build_aqueue (aqueue *Q, int ncount, CCptrworld *intptr_world),
   pop_from_active_queue (aqueue *Q, CCptrworld *intptr_world),
   build_distobj (distobj *D, int ncount, CCdatagroup *dat),
   dist (int i, int j, distobj *D), 
   init_flipstack (flipstack *f, int total, int single),
   check_cycle (graph *G, int *path, int infcount),
   fixed_edge (int n1, int n2, intptr **fixlist),
   init_fixededges (graph *G, int fcount, int *flist,
        CCptrworld *intptr_world);

static double
   cycle_length (int ncount, int *cyc, distobj *D),
   improve_tour (graph *G, distobj *D, adddel *E, aqueue *Q, CClk_flipper *F,
       int start, flipstack *fstack, CCptrworld *intptr_world,
       CCptrworld *edgelook_world);

static edgelook
   *look_ahead (graph *G, distobj *D, adddel *E, CClk_flipper *F, int first,
       int last, int gain, int level, CCptrworld *edgelook_world),
   *weird_look_ahead  (graph *G, distobj *D, CClk_flipper *F, int gain, int t1,
       int t2, CCptrworld *edgelook_world),
   *weird_look_ahead2 (graph *G, distobj *D, CClk_flipper *F, int gain, int t2,
       int t3, int t4, CCptrworld *edgelook_world),
   *weird_look_ahead3 (graph *G, distobj *D, CClk_flipper *F, int gain, int t2,
       int t3, int t6, CCptrworld *edgelook_world);


CC_PTRWORLD_ROUTINES(intptr, intptralloc, intptr_bulkalloc, intptrfree)
CC_PTRWORLD_LISTFREE_ROUTINE(intptr, intptr_listfree, intptrfree)
CC_PTRWORLD_LEAKS_ROUTINE(intptr, intptr_check_leaks, this, int)

CC_PTRWORLD_ROUTINES(edgelook, edgelookalloc, edgelook_bulkalloc, edgelookfree)
CC_PTRWORLD_LISTFREE_ROUTINE(edgelook, edgelook_listfree, edgelookfree)
CC_PTRWORLD_LEAKS_ROUTINE(edgelook, edgelook_check_leaks, diff, int)


int CClinkern_fixed (int ncount, CCdatagroup *dat, int ecount, int *elist,
        int nkicks, int *incycle, int *outcycle, double *val, int fcount,
        int *flist, int silent, CCrandstate *rstate)
{
    int rval = 0;
    int i;
    int *tcyc = (int *) NULL;
    double startzeit;
    graph G;
    distobj D;
    CCptrworld intptr_world;
    CCptrworld edgelook_world;

    if (silent == 0) {
        printf ("linkern ...\n"); fflush (stdout);
    }
    startzeit = CCutil_zeit ();

    initgraph (&G);
    init_distobj (&D);
    CCptrworld_init (&intptr_world);
    CCptrworld_init (&edgelook_world);
    G.rstate = rstate;

    if (ncount < 10 && nkicks > 0) {
        printf ("Less than 10 nodes, setting nkicks to 0\n");  
        fflush (stdout);
        nkicks = 0;
    }

    /* These bulkalloc's allocate sufficient objects that the individual
     * allocs will not fail, and thus do not need to be tested */
    rval = intptr_bulkalloc (&intptr_world, 3*ncount);
    if (rval) {
        fprintf (stderr, "Unable to allocate initial intptrs\n");
        goto CLEANUP;
    }

    rval = edgelook_bulkalloc (&edgelook_world, MAX_BACK * (BACKTRACK + 3));
    if (rval) {
        fprintf (stderr, "Unable to allocate initial edgelooks\n");
        goto CLEANUP;
    }

    tcyc = CC_SAFE_MALLOC (ncount, int);
    CCcheck_NULL (tcyc, "out of memory in CClinkern_fixed");

    rval = build_distobj (&D, ncount, dat);
    if (rval) goto CLEANUP;
    
    rval = buildgraph (&G, ncount, ecount, elist, &D);
    if (rval) {
        fprintf (stderr, "buildgraph failed\n"); goto CLEANUP;
    }

    rval = init_fixededges (&G, fcount, flist,  &intptr_world);
    CCcheck_rval (rval, "init_fixededges failed");

    if (incycle) {
        for (i = 0; i < ncount; i++) tcyc[i] = incycle[i];
    } else {
        fprintf (stderr, "with fixed linkern a cycle must be given");
        rval = 1;  goto CLEANUP;
    }

    *val = cycle_length (ncount, tcyc, &D);
    if (silent == 0) {
        printf ("Starting Cycle: %.0f\n", *val); fflush (stdout);
    }

    rval = repeated_lin_kernighan (&G, &D, tcyc, nkicks,
                 val, silent, &intptr_world, &edgelook_world);
    CCcheck_rval (rval, "repeated_lin_kernighan failed");

    if (silent == 0) {
        printf ("Best path length: %.0f\n", *val);
        printf ("Lin-Kernighan Running Time: %.2f\n",
                  CCutil_zeit () - startzeit);
        fflush (stdout);
    }

    rval = check_cycle (&G, tcyc, fcount);
    CCcheck_rval (rval, "ERROR: Fixed edge was lost");

    if (outcycle) {
        for (i = 0; i < ncount; i++) outcycle[i] = tcyc[i];
    }

    if (G.fixlist) {
        for (i = 0; i < ncount; i++) {
            intptr_listfree (&intptr_world, G.fixlist[i]);
        }
    }

CLEANUP:

    CC_IFFREE (tcyc, int);
    freegraph (&G);
    free_distobj (&D);
    linkern_free_world (&intptr_world, &edgelook_world);

    return rval;
}

static int repeated_lin_kernighan (graph *G, distobj *D, int *cyc,
        int count, double *val, int silent,
        CCptrworld *intptr_world, CCptrworld *edgelook_world)
{
    int rval    = 0;
    int round   = 0;
    int newtree = 0;
    int hit, delta;
    int *win_cycle = (int *) NULL;
    CCkdtree kdt;
    flipstack winstack, fstack;
    double t, best = *val;
    double szeit = CCutil_zeit ();
    int ncount = G->ncount;
    adddel E;
    CClk_flipper F;
    aqueue Q;

    init_aqueue (&Q);
    init_adddel (&E);
    rval = build_aqueue (&Q, ncount, intptr_world);
    if (rval) {
        fprintf (stderr, "build_aqueue failed\n"); goto CLEANUP;
    }
    rval = build_adddel (&E, ncount);
    if (rval) {
        fprintf (stderr, "build_adddel failed\n"); goto CLEANUP;
    }

    hit = 2 * (MAXDEPTH + 7 + KICK_MAXDEPTH);
    rval = init_flipstack (&fstack, hit, 0);
    if (rval) {
        fprintf (stderr, "init_flipstack failed\n"); goto CLEANUP;
    }
    rval = init_flipstack (&winstack, 500 + ncount / 50, hit);
    if (rval) {
        fprintf (stderr, "init_flipstack failed\n"); goto CLEANUP;
    }

    win_cycle = CC_SAFE_MALLOC (ncount, int);
    if (win_cycle == (int *) NULL) {
        fprintf (stderr, "out of memory in repeated_lin_kernighan\n");
        rval = 1; goto CLEANUP;
    }
    win_cycle[0] = -1;

    CClinkern_flipper_init (&F, ncount, cyc);
    fstack.counter = 0;
    winstack.counter = 0;
    win_cycle[0] = -1;

    {
        int *tcyc = (int *) NULL;
        int i;

        tcyc = CC_SAFE_MALLOC (ncount, int);
        if (tcyc == (int *) NULL) {
            fprintf (stderr, "out of memory in repeated_lin_kernighan\n");
            rval = 1; goto CLEANUP;
        }
        /* init active_queue with random order */
        randcycle (ncount, tcyc, G->rstate);
        for (i = 0; i < ncount; i++) {
            add_to_active_queue (tcyc[i], &Q, intptr_world);
        }
        CC_IFFREE (tcyc, int);
    }

    lin_kernighan (G, D, &E, &Q, &F, &best, win_cycle, &winstack, &fstack,
                   intptr_world, edgelook_world);

    winstack.counter = 0;
    win_cycle[0] = -1;

    if (silent == 0) {
        if (count > 0) {
            printf ("%4d Steps   Best: %.0f   %.2f seconds\n", round, best,
                                CCutil_zeit () - szeit);
        } else {
            printf ("LK Cycle: %.0f\n", best);
        }
        fflush (stdout);
    }

    while (round < count) {
        hit = 0;
        fstack.counter = 0;

        rval = random_four_swap (G, D, &Q, &F, &delta,
                                 &winstack, &fstack, intptr_world);
        CCcheck_rval (rval, "random_four_swap failed");

        fstack.counter = 0;
        t = best + delta;
        lin_kernighan (G, D, &E, &Q, &F, &t, win_cycle, &winstack, &fstack,
                       intptr_world, edgelook_world);

#ifdef ACCEPT_TIES
        if (t <= best) {
#else
        if (t < best) {
#endif /* ACCEPT_TIES */
            winstack.counter = 0;
            win_cycle[0] = -1;
            if (t < best) {
                best = t;
                hit++;
            }
        } else {
            if (win_cycle[0] == -1) {
                while (winstack.counter) {
                    winstack.counter--;
                    CClinkern_flipper_flip (&F,
                                      winstack.stack[winstack.counter].last, 
                                      winstack.stack[winstack.counter].first);
                }
            } else {
                CClinkern_flipper_finish (&F);
                CClinkern_flipper_init (&F, ncount, win_cycle);
                while (winstack.counter) {
                    winstack.counter--;
                    CClinkern_flipper_flip (&F,
                                      winstack.stack[winstack.counter].last, 
                                      winstack.stack[winstack.counter].first);
                }
                win_cycle[0] = -1;
            }
        }

        round++;
        if (silent == 0 && (hit || (round % 1000 == 999))) {
            printf ("%4d Steps   Best: %.0f   %.2f seconds\n",
                               round, best, CCutil_zeit () - szeit);
            fflush (stdout);
        }
    }
    if (silent == 0 && round > 0) {
        printf ("%4d Total Steps.\n", round); fflush (stdout);
    }

    CClinkern_flipper_cycle (&F, cyc);
    CClinkern_flipper_finish (&F);

    t = cycle_length (ncount, cyc, D);
    if (t != best) {
        printf ("WARNING: LK incremental counter was off by %.0f\n", t-best);
        fflush (stdout);
        best = t;
    }
    *val = best;

CLEANUP:

    free_aqueue (&Q, intptr_world);
    free_adddel (&E);
    free_flipstack (&fstack);
    free_flipstack (&winstack);
    CC_IFFREE (win_cycle, int);
    if (newtree) CCkdtree_free (&kdt);
    return rval;
}

static void lin_kernighan (graph *G, distobj *D, adddel *E, aqueue *Q,
        CClk_flipper *F, double *val, int *win_cycle, flipstack *win,
        flipstack *fstack, CCptrworld *intptr_world,
        CCptrworld *edgelook_world)
{
    int start, i;
    double delta, totalwin = 0.0;

    while (1) {
        start = pop_from_active_queue (Q, intptr_world);
        if (start == -1) break;

        delta = improve_tour (G, D, E, Q, F, start, fstack, intptr_world,
                              edgelook_world);
        if (delta > 0.0) {
            totalwin += delta;
            if (win->counter < win->max) {
                for (i = 0; i < fstack->counter; i++) {
                    win->stack[win->counter].first = fstack->stack[i].first;
                    win->stack[win->counter].last  = fstack->stack[i].last;
                    win->stack[win->counter].firstprev =
                             fstack->stack[i].firstprev;
                    win->stack[win->counter].lastnext =
                             fstack->stack[i].lastnext;
                    win->counter++;
                }
            } else if (win_cycle[0] == -1) {
                for (i = 0; i < fstack->counter; i++) {
                    win->stack[win->counter].first = fstack->stack[i].first;
                    win->stack[win->counter].last  = fstack->stack[i].last;
                    win->counter++;
                }
                CClinkern_flipper_cycle (F, win_cycle);
            }
            fstack->counter = 0;
        }
    }

    if (win_cycle[0] == -1) {
        for (i = 0; i < fstack->counter; i++) {
            win->stack[win->counter].first = fstack->stack[i].first;
            win->stack[win->counter].last  = fstack->stack[i].last;
            win->counter++;
        }
    }
    (*val) -= totalwin;
}

static double improve_tour (graph *G, distobj *D, adddel *E, aqueue *Q,
        CClk_flipper *F, int t1, flipstack *fstack,
        CCptrworld *intptr_world, CCptrworld *edgelook_world)
{
    int t2 = CClinkern_flipper_next (F, t1);
    int win, gain, Gstar = 0;

    if (Fixededge (t1, t2)) return 0.0;

    gain = Edgelen (t1, t2, D);
    markedge_del (t1, t2, E);

    win = step (G, D, E, Q, F, 0, gain, &Gstar, t1, t2, fstack, intptr_world,
               edgelook_world);
    if (win == 0) {
        Gstar = weird_second_step (G, D, E, Q, F, gain, t1, t2, fstack,
                                   intptr_world, edgelook_world); 
    }
    unmarkedge_del (t1, t2, E);

    if (Gstar) {
        MARK (t1, Q, F, D, G, intptr_world);
        MARK(t2, Q, F, D, G, intptr_world);
    }
    return (double) Gstar;
}

static int step (graph *G, distobj *D, adddel *E, aqueue *Q, CClk_flipper *F,
        int level, int gain, int *Gstar, int first, int last,
        flipstack *fstack, CCptrworld *intptr_world, 
        CCptrworld *edgelook_world)
{
    int val, this, newlast, hit = 0, oldG = gain;
    edgelook *list, *e;

    if (level >= BACKTRACK) {
        return step_noback (G, D, E, Q, F, level, gain, Gstar, first, last,
                            fstack, intptr_world);
    }

    list = look_ahead (G, D, E, F, first, last, gain, level, edgelook_world);
    for (e = list; e; e = e->next) {
        this = e->other;
        newlast = e->over;
   
        if (Fixededge (this, newlast)) continue;

        gain = oldG - e->diff;
        val = gain - Edgelen (newlast, first, D);
        if (val > *Gstar) {
            *Gstar = val;
            hit++;
        }

        FLIP (first, last, newlast, this, fstack, F);

        if (level < MAXDEPTH) {
            markedge_add (last, this, E);
            markedge_del (this, newlast, E);
            hit += step (G, D, E, Q, F, level + 1, gain, Gstar, first,
                         newlast, fstack, intptr_world, edgelook_world);
            unmarkedge_add (last, this, E);
            unmarkedge_del (this, newlast, E);
        }

        if (!hit) {
            UNFLIP (first, last, newlast, this, fstack, F);
        } else {
            MARK (this, Q, F, D, G, intptr_world);
            MARK (newlast, Q, F, D, G, intptr_world);
            edgelook_listfree (edgelook_world, list);
            return 1;
        }
    }
    edgelook_listfree (edgelook_world, list);
    return 0;
}

static int step_noback (graph *G, distobj *D, adddel *E, aqueue *Q,
        CClk_flipper *F, int level, int gain, int *Gstar, int first, int last,
        flipstack *fstack, CCptrworld *intptr_world)
{
    edgelook e;

    look_ahead_noback (G, D, E, F, first, last, gain - *Gstar - level, &e);

    if (e.diff < BIGINT) {
        if (e.mm) {
            int hit = 0;
            int this = e.other;
            int newfirst = e.over;
            int val;

            if (Fixededge (this, newfirst)) return 0;

            gain -= e.diff;
            val = gain - Edgelen (newfirst, last, D);
            if (val > *Gstar) {
                *Gstar = val;
                hit++;
            }
            FLIP (this, newfirst, first, last, fstack, F);

            if (level < MAXDEPTH) {
                markedge_add (first, this, E);
                markedge_del (this, newfirst, E);
                hit += step_noback (G, D, E, Q, F, level + 1, gain, Gstar,
                                    newfirst, last, fstack, intptr_world);
                unmarkedge_add (first, this, E);
                unmarkedge_del (this, newfirst, E);
            }

            if (!hit) {
                UNFLIP (this, newfirst, first, last, fstack, F);
                return 0;
            } else {
                MARK (this, Q, F, D, G, intptr_world);
                MARK (newfirst, Q, F, D, G, intptr_world);
                return 1;
            }
        } else {
            int hit = 0;
            int this = e.other;
            int newlast = e.over;
            int val;

            if (Fixededge (this, newlast)) return 0;

            gain -= e.diff;
            val = gain - Edgelen (newlast, first, D);
            if (val > *Gstar) {
                *Gstar = val;
                hit++;
            }

            FLIP (first, last, newlast, this, fstack, F);

            if (level < MAXDEPTH) {
                markedge_add (last, this, E);
                markedge_del (this, newlast, E);
                hit += step_noback (G, D, E, Q, F, level + 1, gain, Gstar,
                                    first, newlast, fstack, intptr_world);
                unmarkedge_add (last, this, E);
                unmarkedge_del (this, newlast, E);
            }

            if (!hit) {
                UNFLIP (first, last, newlast, this, fstack, F);
                return 0;
            } else {
                MARK (this, Q, F, D, G, intptr_world);
                MARK (newlast, Q, F, D, G, intptr_world);
                return 1;
            }
        }
    } else {
        return 0;
    }
}

static int weird_second_step (graph *G, distobj *D, adddel *E, aqueue *Q,
        CClk_flipper *F, int len_t1_t2, int t1, int t2, flipstack *fstack,
        CCptrworld *intptr_world, CCptrworld *edgelook_world)
{
    int t3, t4, t5, t6, t7, t8;
    int oldG, gain, tG, Gstar = 0, val, hit;
    int t3prev, t4next;
    edgelook *e, *f, *h, *list, *list2, *list3;

    list = weird_look_ahead (G, D, F, len_t1_t2, t1, t2, edgelook_world);
    for (h = list; h; h = h->next) {
        t3 = h->other;
        t4 = h->over;

        if (Fixededge (t3, t4)) continue;

        oldG = len_t1_t2 - h->diff;
  
        t3prev = CClinkern_flipper_prev (F, t3);
        t4next = CClinkern_flipper_next (F, t4);
  
        markedge_add (t2, t3, E);
        markedge_del (t3, t4, E);
        G->weirdmagic++;
        G->weirdmark[t1] = G->weirdmagic;
        G->weirdmark[t2] = G->weirdmagic;
        G->weirdmark[t3] = G->weirdmagic;
        G->weirdmark[t4next] = G->weirdmagic;

        list2 = weird_look_ahead2 (G, D, F, oldG, t2, t3, t4, edgelook_world);
        for (e = list2; e; e = e->next) {
            t5 = e->other;
            t6 = e->over;

            if (Fixededge (t5, t6)) continue;

            markedge_add (t4, t5, E);
            if (e->seq) {
                if (!e->side) {
                    gain = oldG - e->diff;
                    val = gain - Edgelen (t6, t1, D);
                    if (val > Gstar)
                        Gstar = val;
                    FLIP (t1, t2, t6, t5, fstack, F);
                    FLIP (t2, t5, t3, t4, fstack, F);

                    markedge_del (t5, t6, E);
                    hit = step (G, D, E, Q, F, 2, gain, &Gstar, t1, t6, 
                                fstack, intptr_world, edgelook_world);
                    unmarkedge_del (t5, t6, E);

                    if (!hit && Gstar) 
                        hit = 1;

                    if (!hit) {
                        UNFLIP (t2, t5, t3, t4, fstack, F);
                        UNFLIP (t1, t2, t6, t5, fstack, F);
                    } else {
                        unmarkedge_add (t2, t3, E);
                        unmarkedge_del (t3, t4, E);
                        unmarkedge_add (t4, t5, E);
                        MARK (t3, Q, F, D, G, intptr_world);
                        MARK (t4, Q, F, D, G, intptr_world);
                        MARK (t5, Q, F, D, G, intptr_world);
                        MARK (t6, Q, F, D, G, intptr_world);
                        edgelook_listfree (edgelook_world, list);
                        edgelook_listfree (edgelook_world, list2);
                        return Gstar;
                    }
                } else {   
                    gain = oldG - e->diff;
                    val = gain - Edgelen (t6, t1, D);
                    if (val > Gstar)
                        Gstar = val;
                    FLIP (t1, t2, t3, t4, fstack, F);
                    FLIP (t6, t5, t2, t4, fstack, F);
                    FLIP (t1, t3, t6, t2, fstack, F);

                    markedge_del (t5, t6, E);
                    hit = step (G, D, E, Q, F, 2, gain, &Gstar, t1, t6,
                                fstack, intptr_world, edgelook_world);
                    unmarkedge_del (t5, t6, E);

                    if (!hit && Gstar) 
                        hit = 1;

                    if (!hit) {
                        UNFLIP (t1, t3, t6, t2, fstack, F);
                        UNFLIP (t6, t5, t2, t4, fstack, F);
                        UNFLIP (t1, t2, t3, t4, fstack, F); 
                    } else {
                        unmarkedge_add (t2, t3, E);
                        unmarkedge_del (t3, t4, E);
                        unmarkedge_add (t4, t5, E);
                        MARK (t3, Q, F, D, G, intptr_world);
                        MARK (t4, Q, F, D, G, intptr_world);
                        MARK (t5, Q, F, D, G, intptr_world);
                        MARK (t6, Q, F, D, G, intptr_world);
                        edgelook_listfree (edgelook_world, list);
                        edgelook_listfree (edgelook_world, list2);
                        return Gstar;
                    }
                }
            } else {
                tG = oldG - e->diff;
                markedge_del (t5, t6, E);
                list3 = weird_look_ahead3 (G, D, F, tG, t2, t3, t6,
                                           edgelook_world);
                for (f = list3; f; f = f->next) {
                    t7 = f->other;
                    t8 = f->over;
                    if (Fixededge (t7, t8)) continue;

                    gain = tG - f->diff;
                    if (!f->side) {
                        val = gain - Edgelen (t8, t1, D);
                        if (val > Gstar)
                            Gstar = val;
                        FLIP (t1, t2, t8, t7, fstack, F);
                        FLIP (t2, t7, t3, t4, fstack, F);
                        FLIP (t7, t4, t6, t5, fstack, F);

                        markedge_add (t6, t7, E);
                        markedge_del (t7, t8, E);
                        hit = step (G, D, E, Q, F, 3, gain, &Gstar, t1, t8,
                                    fstack, intptr_world, edgelook_world);
                        unmarkedge_del (t6, t7, E);
                        unmarkedge_del (t7, t8, E);

                        if (!hit && Gstar) 
                            hit = 1;

                        if (!hit) {
                            UNFLIP (t7, t4, t6, t5, fstack, F);
                            UNFLIP (t2, t7, t3, t4, fstack, F);
                            UNFLIP (t1, t2, t8, t7, fstack, F);
                        } else {
                            unmarkedge_add (t2, t3, E);
                            unmarkedge_del (t3, t4, E);
                            unmarkedge_add (t4, t5, E);
                            unmarkedge_del (t5, t6, E);
                            MARK (t3, Q, F, D, G, intptr_world);
                            MARK (t4, Q, F, D, G, intptr_world);
                            MARK (t5, Q, F, D, G, intptr_world);
                            MARK (t6, Q, F, D, G, intptr_world);
                            MARK (t7, Q, F, D, G, intptr_world);
                            MARK (t8, Q, F, D, G, intptr_world);
                            edgelook_listfree (edgelook_world, list);
                            edgelook_listfree (edgelook_world, list2);
                            edgelook_listfree (edgelook_world, list3);
                            return Gstar;
                        }
                    } else {
                        val = gain - Edgelen (t8, t1, D);
                        if (val > Gstar)
                            Gstar = val;
                        FLIP (t1, t2, t6, t5, fstack, F);
                        FLIP (t1, t6, t8, t7, fstack, F);
                        FLIP (t3, t4, t2, t5, fstack, F);

                        markedge_add (t6, t7, E);
                        markedge_del (t7, t8, E);
                        hit = step (G, D, E, Q, F, 3, gain, &Gstar, t1, t8,
                                    fstack, intptr_world, edgelook_world);
                        unmarkedge_add (t6, t7, E);
                        unmarkedge_del (t7, t8, E);

                        if (!hit && Gstar) 
                            hit = 1;

                        if (!hit) {
                            UNFLIP (t3, t4, t2, t5, fstack, F);
                            UNFLIP (t1, t6, t8, t7, fstack, F);
                            UNFLIP (t1, t2, t6, t5, fstack, F);
                        } else {
                            unmarkedge_add (t2, t3, E);
                            unmarkedge_del (t3, t4, E);
                            unmarkedge_add (t4, t5, E);
                            unmarkedge_del (t5, t6, E);
                            MARK (t3, Q, F, D, G, intptr_world);
                            MARK (t4, Q, F, D, G, intptr_world);
                            MARK (t5, Q, F, D, G, intptr_world);
                            MARK (t6, Q, F, D, G, intptr_world);
                            MARK (t7, Q, F, D, G, intptr_world);
                            MARK (t8, Q, F, D, G, intptr_world);
                            edgelook_listfree (edgelook_world, list);
                            edgelook_listfree (edgelook_world, list2);
                            edgelook_listfree (edgelook_world, list3);
                            return Gstar;
                        }
                    }
                }
                edgelook_listfree (edgelook_world, list3);
                unmarkedge_del (t5, t6, E);
            }
            unmarkedge_add (t4, t5, E);
        }
        edgelook_listfree (edgelook_world, list2);
        unmarkedge_add (t2, t3, E);
        unmarkedge_del (t3, t4, E);
    }
    edgelook_listfree (edgelook_world, list);
    return 0;
}

static edgelook *look_ahead (graph *G, distobj *D, adddel *E, CClk_flipper *F,
       int first, int last, int gain, int level, CCptrworld *edgelook_world)
{
    edgelook *list = (edgelook *) NULL, *el;
    int i, val;
    int this, prev;
    int lastnext = CClinkern_flipper_next (F, last);
    int other[MAX_BACK], save[MAX_BACK];
    int value[MAX_BACK + 1];
    int k, ahead = backtrack_count[level];
    edge **goodlist = G->goodlist;

    for (i = 0; i < ahead; i++) {
        value[i] = BIGINT;
    }
    value[ahead] = -BIGINT;

#ifdef USE_LESS_OR_EQUAL
    for (i = 0; goodlist[last][i].weight <= gain; i++) {
#else
    for (i = 0; goodlist[last][i].weight < gain; i++) {
#endif
        this = goodlist[last][i].other;
        if (!is_it_deleted (last, this, E) && this != first &&
                                              this != lastnext) {
            prev = CClinkern_flipper_prev (F, this);
            if (!is_it_added (this, prev, E)) {
                val = goodlist[last][i].weight - Edgelen (this, prev, D);
                if (val < value[0]) {
                    for (k = 0; value[k+1] > val; k++) {
                        value[k] = value[k+1];
                        other[k] = other[k+1];
                        save[k] = save[k+1];
                    }
                    value[k] = val;
                    other[k] = this;
                    save[k] = prev;
                }
            }
        }
    }


    for (i = 0; i < ahead; i++) {
        if (value[i] < BIGINT) {
            el = edgelookalloc (edgelook_world);
            el->diff = value[i];
            el->other = other[i];
            el->over = save[i];
            el->next = list;
            list = el;
        }
    }

    return list;
}

static void look_ahead_noback (graph *G, distobj *D, adddel *E, CClk_flipper *F,
        int first, int last, int gain, edgelook *winner)
{
    int val;
    int this, prev;
    int lastnext = CClinkern_flipper_next (F, last);
    int i, next;
    edge **goodlist = G->goodlist;

    winner->diff = BIGINT;
    for (i = 0; goodlist[last][i].weight < gain; i++) {
        this = goodlist[last][i].other;
        if (!is_it_deleted (last, this, E) && this != first &&
                                              this != lastnext) {
            prev = CClinkern_flipper_prev (F, this);
            if (!is_it_added (this, prev, E)) {
                val = goodlist[last][i].weight - Edgelen (this, prev, D);
                if (val < winner->diff) {
                    winner->diff = val;
                    winner->other = this;
                    winner->over = prev;
                    winner->mm = 0;
                }
            }
        }
    }
    {
        int firstprev = CClinkern_flipper_prev (F, first);

        for (i = 0; goodlist[first][i].weight < gain; i++) {
            this = goodlist[first][i].other;
            if (!is_it_deleted (first, this, E) && this != last && 
                                                   this != firstprev) {
                next = CClinkern_flipper_next (F, this);
                if (!is_it_added (this, next, E)) {
                    val = goodlist[first][i].weight - Edgelen (this, next, D);
                    if (val < winner->diff) {
                        winner->diff = val;
                        winner->other = this;
                        winner->over = next;
                        winner->mm = 1;
                    }
                }
            }
        }
    }
}

static edgelook *weird_look_ahead (graph *G, distobj *D, CClk_flipper *F,
        int gain, int t1, int t2, CCptrworld *edgelook_world)
{
    edgelook *list, *el;
    int i, this, next;
    int other[MAX_BACK], save[MAX_BACK];
    int value[MAX_BACK + 1];
    int k, val, ahead;
    edge **goodlist = G->goodlist;

    list = (edgelook *) NULL;
    ahead = weird_backtrack_count[0];
    for (i = 0; i < ahead; i++)
        value[i] = BIGINT;
    value[ahead] = -BIGINT;

#ifdef USE_LESS_OR_EQUAL
    for (i = 0; goodlist[t2][i].weight <= gain; i++) {
#else
    for (i = 0; goodlist[t2][i].weight < gain; i++) {
#endif
        this = goodlist[t2][i].other;
        if (this != t1) {
            next = CClinkern_flipper_next (F, this);
            val = goodlist[t2][i].weight - Edgelen (this, next, D);
            if (val < value[0]) {
                for (k = 0; value[k+1] > val; k++) {
                    value[k] = value[k+1];
                    other[k] = other[k+1];
                    save[k] = save[k+1];
                }
                value[k] = val;
                other[k] = this;
            save[k] = next;
            }
        }
    }
    for (i = 0; i < ahead; i++) {
        if (value[i] < BIGINT) {
            el = edgelookalloc (edgelook_world);
            el->diff = value[i];
            el->other = other[i];
            el->over = save[i];
            el->next = list;
            list = el;
        }
    }
    return list;
}

static edgelook *weird_look_ahead2 (graph *G, distobj *D, CClk_flipper *F,
       int gain, int t2, int t3, int t4, CCptrworld *edgelook_world)
{
    edgelook *list = (edgelook *) NULL;
    edgelook *el;
    int i, t5, t6;
    int other[MAX_BACK], save[MAX_BACK], seq[MAX_BACK], side[MAX_BACK];
    int value[MAX_BACK + 1];
    int k, val;
    int ahead = weird_backtrack_count[1];
    edge **goodlist = G->goodlist;
    int  *weirdmark = G->weirdmark;
    int  weirdmagic = G->weirdmagic;

    for (i = 0; i < ahead; i++)
        value[i] = BIGINT;
    value[ahead] = -BIGINT;

#ifdef USE_LESS_OR_EQUAL
    for (i = 0; goodlist[t4][i].weight <= gain; i++) {
#else
    for (i = 0; goodlist[t4][i].weight < gain; i++) {
#endif
        t5 = goodlist[t4][i].other;
        if (weirdmark[t5] != weirdmagic) {
            if (CClinkern_flipper_sequence (F, t2, t5, t3)) {
                t6 = CClinkern_flipper_prev (F, t5);
                val = goodlist[t4][i].weight - Edgelen (t5, t6, D);
                if (val < value[0]) {
                    for (k = 0; value[k+1] > val; k++) {
                        value[k] = value[k+1];
                        other[k] = other[k+1];
                        save[k] = save[k+1];
                        seq[k] = seq[k+1];
                        side[k] = side[k+1];
                    }
                    value[k] = val;
                    other[k] = t5;
                    save[k] = t6;
                    seq[k] = 1;
                    side[k] = 0;
                }
                t6 = CClinkern_flipper_next (F, t5);
                val = goodlist[t4][i].weight - Edgelen (t5, t6, D);
                if (val < value[0]) {
                    for (k = 0; value[k+1] > val; k++) {
                        value[k] = value[k+1];
                        other[k] = other[k+1];
                        save[k] = save[k+1];
                        seq[k] = seq[k+1];
                        side[k] = side[k+1];
                    }
                    value[k] = val;
                    other[k] = t5;
                    save[k] = t6;
                    seq[k] = 1;
                    side[k] = 1;
                }
            } else {
                t6 = CClinkern_flipper_prev (F, t5);
                val = goodlist[t4][i].weight - Edgelen (t5, t6, D);
                if (val < value[0]) {
                    for (k = 0; value[k+1] > val; k++) {
                        value[k] = value[k+1];
                        other[k] = other[k+1];
                        save[k] = save[k+1];
                        seq[k] = seq[k+1];
                        side[k] = side[k+1];
                    }
                    value[k] = val;
                    other[k] = t5;
                    save[k] = t6;
                    seq[k] = 0;
                    side[k] = 0;
                }
            }
        }
    }

    for (i = 0; i < ahead; i++) {
        if (value[i] < BIGINT) {
            el = edgelookalloc (edgelook_world);
            el->diff = value[i];
            el->other = other[i];
            el->over = save[i];
            el->seq = seq[i];
            el->side = side[i];
            el->next = list;
            list = el;
        }
    }
    return list;
}

static edgelook *weird_look_ahead3 (graph *G, distobj *D, CClk_flipper *F,
        int gain, int t2, int t3, int t6, CCptrworld *edgelook_world)
{
    edgelook *list = (edgelook *) NULL;
    edgelook *el;
    int i, t7, t8;
    int other[MAX_BACK], save[MAX_BACK], side[MAX_BACK];
    int value[MAX_BACK + 1];
    int k, val;
    int ahead = weird_backtrack_count[2];
    edge **goodlist = G->goodlist;
    int  *weirdmark = G->weirdmark;
    int  weirdmagic = G->weirdmagic;

    for (i = 0; i < ahead; i++)
        value[i] = BIGINT;
    value[ahead] = -BIGINT;

#ifdef USE_LESS_OR_EQUAL
    for (i = 0; goodlist[t6][i].weight <= gain; i++) {
#else 
    for (i = 0; goodlist[t6][i].weight < gain; i++) {
#endif
        t7 = goodlist[t6][i].other;   /* Need t7 != t2, t3, t2next, t3prev */
        if (weirdmark[t7] != weirdmagic &&
                   CClinkern_flipper_sequence (F, t2, t7, t3)) {
            t8 = CClinkern_flipper_prev (F, t7);
            val = goodlist[t6][i].weight - Edgelen (t7, t8, D);
            if (val < value[0]) {
                for (k = 0; value[k+1] > val; k++) {
                    value[k] = value[k+1];
                    other[k] = other[k+1];
                    save[k] = save[k+1];
                    side[k] = side[k+1];
                }
                value[k] = val;
                other[k] = t7;
                save[k] = t8;
                side[k] = 0;
            }
            t8 = CClinkern_flipper_next (F, t7);
            val = goodlist[t6][i].weight - Edgelen (t7, t8, D);
            if (val < value[0]) {
                for (k = 0; value[k+1] > val; k++) {
                    value[k] = value[k+1];
                    other[k] = other[k+1];
                    save[k] = save[k+1];
                    side[k] = side[k+1];
                }
                value[k] = val;
                other[k] = t7;
                save[k] = t8;
                side[k] = 1;
            }
        }
    }

    for (i = 0; i < ahead; i++) {
        if (value[i] < BIGINT) {
            el = edgelookalloc (edgelook_world);
            el->diff = value[i];
            el->other = other[i];
            el->over = save[i];
            el->side = side[i];
            el->next = list;
            list = el;
        }
    }
    return list;
}

static double cycle_length (int ncount, int *cyc, distobj *D)
{
    int i;
    double val = 0.0;
    
    for (i = 1; i < ncount; i++) {
        val += (double) Edgelen (cyc[i - 1], cyc[i], D);
    }
    val += (double) Edgelen (cyc[0], cyc[ncount - 1], D);

    return val;
}

static int random_four_swap (graph *G, distobj *D, aqueue *Q, CClk_flipper *F,
       int *delta, flipstack *win, flipstack *fstack, CCptrworld *intptr_world)
{
    int t1, t2, t3, t4, t5, t6, t7, t8, temp;

    find_walk_four (G, D, F, &t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8);

    if (!CClinkern_flipper_sequence (F, t1, t3, t5)) {
        CC_SWAP (t3, t5, temp);
        CC_SWAP (t4, t6, temp);
    }
    if (!CClinkern_flipper_sequence (F, t1, t5, t7)) {
        CC_SWAP (t5, t7, temp);
        CC_SWAP (t6, t8, temp);
        if (!CClinkern_flipper_sequence (F, t1, t3, t5)) {
            CC_SWAP (t3, t5, temp);
            CC_SWAP (t4, t6, temp);
        }
    }
    FLIP (t1, t2, t5, t6, fstack, F);
    FLIP (t4, t3, t7, t8, fstack, F);
    FLIP (t1, t5, t6, t2, fstack, F);

    if (win->counter < win->max) {
        win->stack[win->counter].first = t2;
        win->stack[win->counter].last = t5;
        win->counter++;
    }
    if (win->counter < win->max) {
        win->stack[win->counter].first = t3;
        win->stack[win->counter].last = t7;
        win->counter++;
    }
    if (win->counter < win->max) {
        win->stack[win->counter].first = t5;
        win->stack[win->counter].last = t6;
        win->counter++;
    }

    bigturn (G ,t1, 0, Q, F, D, intptr_world);
    bigturn (G, t2, 1, Q, F, D, intptr_world);
    bigturn (G, t3, 0, Q, F, D, intptr_world);
    bigturn (G, t4, 1, Q, F, D, intptr_world);
    bigturn (G, t5, 0, Q, F, D, intptr_world);
    bigturn (G, t6, 1, Q, F, D, intptr_world);
    bigturn (G, t7, 0, Q, F, D, intptr_world);
    bigturn (G, t8, 1, Q, F, D, intptr_world);

    *delta = 
           Edgelen (t1, t6, D) + Edgelen (t2, t5, D) +
           Edgelen (t3, t8, D) + Edgelen (t4, t7, D) -
           Edgelen (t1, t2, D) - Edgelen (t3, t4, D) -
           Edgelen (t5, t6, D) - Edgelen (t7, t8, D);
    return 0;
}

#define HUNT_PORTION_LONG 0.001 

static void first_kicker (graph *G, distobj *D, CClk_flipper *F, int *t1,
        int *t2)
{
#ifdef LONG_KICKER
    int longcount = (int) ((double) G->ncount * HUNT_PORTION_LONG) + 10;
    int i, best, try1, len, next, prev, nextl, prevl;
    int ncount = G->ncount;
    edge **goodlist = G->goodlist;

    try1 = CCutil_lprand (G->rstate) % ncount;
    next = CClinkern_flipper_next (F, try1);
    prev = CClinkern_flipper_prev (F, try1);
    nextl = Edgelen (try1, next, D);
    prevl = Edgelen (try1, prev, D);
    if (nextl >= prevl) {
        *t1 = try1;
        *t2 = next;
        best = nextl - goodlist[*t1][0].weight;
    } else {
        *t1 = prev;
        *t2 = try1;
        best = prevl - goodlist[*t1][0].weight;
    }

    for (i = 0; i < longcount; i++) {
        try1 = CCutil_lprand (G->rstate) % ncount;
        next = CClinkern_flipper_next (F, try1);
        prev = CClinkern_flipper_prev (F, try1);
        nextl = Edgelen (try1, next, D);
        prevl = Edgelen (try1, prev, D);
        if (nextl >= prevl) {
            len = nextl - goodlist[try1][0].weight;
            if (len > best) {
                *t1 = try1;
                *t2 = next;
            }
        } else {
            len = prevl - goodlist[try1][0].weight;
            if (len > best) {
                *t1 = prev;
                *t2 = try1;
            }
        }
    }
#else   /* LONG_KICKER */
    *t1 = CCutil_lprand (G->rstate) % G->ncount;
    *t2 = CClinkern_flipper_next (F, *t1);
#endif  /* LONG_KICKER */
}


#define WALK_STEPS 50

static void find_walk_four (graph *G, distobj *D, CClk_flipper *F, int *t1,
        int *t2, int *t3, int *t4, int *t5, int *t6, int *t7, int *t8)
{
    int s1, s2, s3, s4, s5, s6, s7, s8;
    int old, n, i, j;

/*
    s1 = CCutil_lprand (G->rstate) % G->ncount;
    s2 = CClinkern_flipper_next (F, s1);
*/

    do {
        first_kicker (G, D, F, &s1, &s2);
    } while (Fixededge (s1, s2));

    do {
        old = -1;
        n = s2;

        for (i = 0;  i < WALK_STEPS; i++) {
            j = CCutil_lprand (G->rstate) % (G->degree[n]);
            if (old != G->goodlist[n][j].other) {
                old = n;
                n = G->goodlist[n][j].other;
            }
        }
        s3 = n;
        s4 = CClinkern_flipper_next (F, s3);

        n = s4;
        for (i = 0; i < WALK_STEPS; i++) {
            j = CCutil_lprand (G->rstate) % (G->degree[n]);
            if (old != G->goodlist[n][j].other) {
                old = n;
                n = G->goodlist[n][j].other;
            }
        }
        s5 = n;
        s6 = CClinkern_flipper_next (F, s5);

        n = s6;
        for (i = 0; i < WALK_STEPS; i++) {
            j = CCutil_lprand (G->rstate) % (G->degree[n]);
            if (old != G->goodlist[n][j].other) {
                old = n;
                n = G->goodlist[n][j].other;
            }
        }
        s7 = n;
        s8 = CClinkern_flipper_next (F, s7);
    } while (s1 == s3 || s1 == s4 || s1 == s5 || s1 == s6 || s1 == s7 ||
             s1 == s8 ||
             s2 == s3 || s2 == s4 || s2 == s5 || s2 == s6 || s2 == s7 ||
             s2 == s8 ||
             s3 == s5 || s3 == s6 || s3 == s7 || s3 == s8 ||
             s4 == s5 || s4 == s6 || s4 == s7 || s4 == s8 ||
             s5 == s7 || s5 == s8 ||
             s6 == s7 || s6 == s8 ||
             Fixededge (s1, s2) || Fixededge (s3, s4) ||
             Fixededge (s5, s6) || Fixededge (s7, s8));

    *t1 = s1;  *t2 = s2;  *t3 = s3;  *t4 = s4;
    *t5 = s5;  *t6 = s6;  *t7 = s7;  *t8 = s8;
}

static void turn (int n, aqueue *Q, CCptrworld *intptr_world)
{
    add_to_active_queue (n, Q, intptr_world);
}

static void kickturn (int n, aqueue *Q, CC_UNUSED distobj *D,
        CC_UNUSED graph *G, CClk_flipper *F, CCptrworld *intptr_world)
{
    add_to_active_queue (n, Q, intptr_world);
    {
        int k;
        k = CClinkern_flipper_next (F, n);
        add_to_active_queue (k, Q, intptr_world);
        k = CClinkern_flipper_next (F, k);
        add_to_active_queue (k, Q, intptr_world);
        k = CClinkern_flipper_prev (F, n);
        add_to_active_queue (k, Q, intptr_world);
        k = CClinkern_flipper_prev (F, k);
        add_to_active_queue (k, Q, intptr_world);
    }
}

static void bigturn (graph *G, int n, int tonext, aqueue *Q, CClk_flipper *F,
        CC_UNUSED distobj *D, CCptrworld *intptr_world)
{
    int i, k;

    add_to_active_queue (n, Q, intptr_world);
    if (tonext) {
        for (i = 0, k = n; i < MARK_LEVEL; i++) {
            k = CClinkern_flipper_next (F, k);
            add_to_active_queue (k, Q, intptr_world);
        }
    } else {
        for (i = 0, k = n; i < MARK_LEVEL; i++) {
            k = CClinkern_flipper_prev (F, k);
            add_to_active_queue (k, Q, intptr_world);
        }
    }

    for (i = 0; i < G->degree[n]; i++) {
        add_to_active_queue (G->goodlist[n][i].other, Q, intptr_world);
    }
}

static void randcycle (int ncount, int *cyc, CCrandstate *rstate)
{
    int i, k, temp;

    for (i = 0; i < ncount; i++) cyc[i] = i;
    for (i = ncount; i > 1; i--) {
        k = CCutil_lprand (rstate) % i;
        CC_SWAP (cyc[i - 1], cyc[k], temp);
    }
}

static int check_cycle (graph *G, int *path, int infcount)
{
    int i, fcount = 0;

    for (i = 1; i < G->ncount; i++) {
        if (Fixededge (path[i-1], path[i])) fcount++;;
    }
    if (Fixededge (path[0], path[G->ncount-1])) fcount++;

    if (infcount == fcount) return 0;
    else                    return 1;
}

static void initgraph (graph *G)
{
    G->goodlist   = (edge **) NULL;
    G->edgespace  = (edge *) NULL;    
    G->degree     = (int *) NULL;
    G->weirdmark  = (int *) NULL;
    G->weirdmagic = 0;
    G->ncount     = 0;
    G->fixlist    = (intptr **) NULL;
}

static void freegraph (graph *G)
{
    if (G) {
        CC_IFFREE (G->goodlist, edge *);
        CC_IFFREE (G->edgespace, edge);
        CC_IFFREE (G->degree, int);
        CC_IFFREE (G->weirdmark, int);
        G->weirdmagic = 0;
        G->ncount     = 0;
        CC_IFFREE (G->fixlist, intptr *);
    }
}

static int buildgraph (graph *G, int ncount, int ecount, int *elist,
        distobj *D)
{
    int rval = 0;
    int n1, n2, w, i;
    edge *p;

    G->goodlist  = CC_SAFE_MALLOC (ncount, edge *);
    G->degree    = CC_SAFE_MALLOC (ncount, int);
    G->weirdmark = CC_SAFE_MALLOC (ncount, int);
    G->edgespace = CC_SAFE_MALLOC ((2 * ecount) + ncount, edge);
    if (G->goodlist == (edge **) NULL || G->degree == (int *) NULL ||
        G->edgespace == (edge *) NULL)  {
        fprintf (stderr, "out of memory in buildgraph\n");
        rval = 1; goto CLEANUP;
    }

    for (i = 0; i < ncount; i++) {
        G->degree[i] = 1;
        G->weirdmark[i] = 0;
    }
    for (i = ecount - 1; i >= 0; i--) {
        G->degree[elist[2 * i]]++;
        G->degree[elist[(2 * i) + 1]]++;
    }

    for (i = 0, p = G->edgespace; i < ncount; i++) {
        G->goodlist[i] = p;
        p += (G->degree[i]);
        G->goodlist[i][G->degree[i] - 1].weight = BIGINT;
        G->degree[i] = 0;
    }

    for (i = ecount - 1; i >= 0; i--) {
        n1 = elist[2 * i];
        n2 = elist[(2 * i) + 1];
        w = Edgelen (n1, n2, D);
        insertedge (G, n1, n2, w);
        insertedge (G, n2, n1, w);
    }
    G->ncount     = ncount;
    G->weirdmagic = 0;

CLEANUP:

    if (rval) freegraph (G);
    return rval;
}

static void insertedge (graph *G, int n1, int n2, int w)
{
    int i;
    edge *e = G->goodlist[n1];

    for (i = G->degree[n1] - 1; i >= 0 && e[i].weight >= w; i--) {
        e[i + 1].weight = e[i].weight;
        e[i + 1].other  = e[i].other;
    } 
    e[i + 1].weight = w;
    e[i + 1].other  = n2;
    G->degree[n1]++;
}

static void linkern_free_world (CCptrworld *intptr_world,
        CCptrworld *edgelook_world)
{
    int total, onlist;

    if (intptr_check_leaks (intptr_world, &total, &onlist)) {
        fprintf (stderr, "WARNING: %d outstanding intptrs\n", total-onlist);
    }
    if (edgelook_check_leaks (edgelook_world, &total, &onlist)) {
        fprintf (stderr, "WARNING: %d outstanding edgelooks\n", total-onlist);
    }
    CCptrworld_delete (intptr_world);
    CCptrworld_delete (edgelook_world);
}

static int init_flipstack (flipstack *f, int total, int single)
{
    f->counter = 0;
    f->max     = 0;
    f->stack   = (flippair *) NULL;

    f->stack = CC_SAFE_MALLOC (total + single, flippair);
    if (f->stack == (flippair *) NULL) {
        fprintf (stderr, "out of memory in init_flipstack\n"); return 1;
    }
    f->max = total;

    return 0;
}

static void free_flipstack (flipstack *f)
{
    f->counter = 0;
    f->max     = 0;
    CC_IFFREE (f->stack, flippair);
}

static void init_adddel (adddel *E)
{
    E->add_edges = (char *) NULL;
    E->del_edges = (char *) NULL;
}

static void free_adddel (adddel *E)
{
    if (E) {
        CC_IFFREE (E->add_edges, char);
        CC_IFFREE (E->del_edges, char);
    }
}

static int build_adddel (adddel *E, int ncount)
{
    int rval = 0;
    int i, M;

    i = 0;
    while ((1 << i) < ncount)
        i++;
    M = (1 << i);

    E->add_edges = CC_SAFE_MALLOC (M, char);
    E->del_edges = CC_SAFE_MALLOC (M, char);
    if (E->add_edges == (char *) NULL || E->del_edges == (char *) NULL) {
        fprintf (stderr, "out of memory in build_adddel\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < M; i++) {
        E->add_edges[i] = 0;
        E->del_edges[i] = 0;
    }

CLEANUP:

    if (rval) {
        free_adddel (E);
    }
    return rval;
}

static void init_aqueue (aqueue *Q)
{
    Q->active = (char *) NULL;
    Q->active_queue = (intptr *) NULL;
    Q->bottom_active_queue = (intptr *) NULL;
    Q->h = (CCdheap *) NULL;
}

static void free_aqueue (aqueue *Q, CCptrworld *intptr_world)
{
    if (Q) {
        CC_IFFREE (Q->active, char);
        intptr_listfree (intptr_world, Q->active_queue);
        Q->active_queue = (intptr *) NULL;
        Q->bottom_active_queue = (intptr *) NULL;
        if (Q->h) {
            CCutil_dheap_free (Q->h);
            Q->h = (CCdheap *) NULL;
        }
    }
}

static int build_aqueue (aqueue *Q, int ncount, CCptrworld *intptr_world)
{
    int rval = 0;
    int i;

    init_aqueue (Q);

    Q->active = CC_SAFE_MALLOC (ncount, char);
    if (Q->active == (char *) NULL) {
        fprintf (stderr, "out of memory in build_aqueue\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < ncount; i++) Q->active[i] = 0;

CLEANUP:

    if (rval) {
        free_aqueue (Q, intptr_world);
    }
    return rval;
}

static void add_to_active_queue (int n, aqueue *Q, CCptrworld *intptr_world)
{
    intptr *ip;

    /* intptralloc will not fail - the initial supply reserved with
     * intptr_bulkalloc is large enough */

    if (Q->active[n] == 0) { 
        Q->active[n] = 1;
        ip = intptralloc (intptr_world);
        ip->this = n;
        ip->next = (intptr *) NULL;
        if (Q->bottom_active_queue) {
            Q->bottom_active_queue->next = ip;
        } else {
            Q->active_queue = ip;
        }
        Q->bottom_active_queue = ip;
    }
}

static int pop_from_active_queue (aqueue *Q, CCptrworld *intptr_world)
{
    intptr *ip;
    int n = -1;

    if (Q->active_queue != (intptr *) NULL) {
        ip = Q->active_queue;
        n = ip->this;
        Q->active_queue = ip->next;
        if (ip == Q->bottom_active_queue) {
            Q->bottom_active_queue = (intptr *) NULL;
        }
        intptrfree (intptr_world, ip);
        Q->active[n] = 0;
    }
    return n;
}

static void init_distobj (distobj *D)
{
    D->dat = (CCdatagroup *) NULL;
    D->cacheind  = (int *) NULL;
    D->cacheval  = (int *) NULL;
    D->cacheM = 0;
}

static void free_distobj (distobj *D)
{
    if (D) {
         D->dat = (CCdatagroup *) NULL;
         CC_IFFREE (D->cacheind, int);
         CC_IFFREE (D->cacheval, int);
         D->cacheM = 0;
    }
}

static int build_distobj (distobj *D, int ncount, CCdatagroup *dat)
{
    int rval = 0;
    int i;

    init_distobj (D);
    D->dat = dat;

#ifndef BENTLEY_CACHE
    i = 0;
    while ((1 << i) < (ncount << 2))
        i++;
    D->cacheM = (1 << i);  
#else
    i = 0;
    while ((1 << i) < ncount)
        i++;
    D->cacheM = (1 << i);
#endif

    D->cacheind = CC_SAFE_MALLOC (D->cacheM, int);
    D->cacheval = CC_SAFE_MALLOC (D->cacheM, int);
    if (D->cacheind == (int *) NULL || D->cacheval == (int *) NULL) {
        fprintf (stderr, "out of memory in build_distobj\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < D->cacheM; i++) {
        D->cacheind[i] = -1;
    }

#ifndef BENTLEY_CACHE
    D->cacheM--;
#endif

CLEANUP:

    if (rval) {
        free_distobj (D);
    }
    return rval; 
}


static int dist (int i, int j, distobj *D)   /* As in Bentley's kdtree paper */
{
    int ind;

    if (i > j) {
        int temp;
        CC_SWAP (i, j, temp);
    }

#ifndef BENTLEY_CACHE
    ind = (((i << 8) + i + j) & (D->cacheM));
#else
    ind = i ^ j;
#endif

    if (D->cacheind[ind] != i) {
        D->cacheind[ind] = i;
        D->cacheval[ind] = CCutil_dat_edgelen (i, j, D->dat);
    }
    return D->cacheval[ind];
}

static int init_fixededges (graph *G, int fcount, int *flist,
        CCptrworld *intptr_world)
{
    int i, k, n1, n2, temp;
    intptr *ip;
    int ncount = G->ncount;
    int rval = 0;

    G->fixlist = CC_SAFE_MALLOC (ncount, intptr *);
    CCcheck_NULL (G->fixlist, "out of memory in init_fixededges");

    for (i = 0; i < ncount; i++) {
        G->fixlist[i] = (intptr *) NULL;
    }

    for (i = 0, k = 0; i < fcount; i++) {
        n1 = flist[k++];
        n2 = flist[k++];
        if (n1 > n2) {
            CC_SWAP (n1, n2, temp);
        }
        ip = intptralloc (intptr_world);
        ip->this = n2;
        ip->next = G->fixlist[n1];
        G->fixlist[n1] = ip;
    }

CLEANUP:

    return rval;
}

static int fixed_edge (int n1, int n2, intptr **fixlist)
{
    intptr *ip;

    if (n1 > n2) {
        for (ip = fixlist[n2]; ip; ip = ip->next)
            if (ip->this == n1)
                return 1;
        return 0;
    } else {
        for (ip = fixlist[n1]; ip; ip = ip->next)
            if (ip->this == n2)
                return 1;
        return 0;
    }
}

