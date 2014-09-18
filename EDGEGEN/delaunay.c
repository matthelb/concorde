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
/*           EDGESET OF DELAUNAY TRIANGULATION (L2-NORM)                    */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*    *** ONLY MINOR MODIFICATIONS TO STEVE FORTUNE'S SWEEP2 CODE ***       */
/*               *** SWEEP2 IS AVAILABLE FROM NETLIB ***                    */
/*                                                                          */
/*            *** SWEEP2 COMES WITH THE FOLLOWING NOTICE ***                */
/*                                                                          */
/* The author of this software is Steven Fortune.  Copyright (c) 1994 by    */
/* AT&T Bell Laboratories. Permission to use, copy, modify, and distribute  */
/* this software for any purpose without fee is hereby granted, provided    */
/* that this entire notice is included in all copies of any software which  */
/* is or includes a copy or modification of this software and in all copies */
/* of the supporting documentation for such software. THIS SOFTWARE IS      */
/* BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED WARRANTY.  IN     */
/* PARTICULAR, NEITHER THE AUTHORS NOR AT&T MAKE ANY REPRESENTATION OR      */
/* WARRANTY OF ANY KIND CONCERNING THE MERCHANTABILITY OF THIS SOFTWARE OR  */
/* ITS FITNESS FOR ANY PARTICULAR PURPOSE.                                  */
/*                                                                          */
/*                    *** END OF SWEEP2 NOTICE ***                          */
/*                                                                          */
/*                                                                          */
/*  Modifications by:  Applegate, Bixby, Chvatal, and Cook                  */
/*  Date: March 12, 1996                                                    */
/*                                                                          */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCedgegen_delaunay (int ncount, CCdatagroup *dat, int wantlist,     */
/*      int *ecount, int **elist)                                           */
/*    RETURNS the edgeset of the (Euclidean-norm) Delaunay triangulation    */
/*    of the point set for norms of size CC_D2_NORM_SIZE (2-coordinates).   */
/*      -ncount is the number of nodes)                                     */
/*      -dat contains the info to generate edge lengths                     */
/*      -wantlist is 1 if you want the function to return the edges.        */
/*      -ecount returns the number of edges if wantlist is 1                */
/*      -elist returns the edges in end1 end2 format if wantlist is 1       */
/*    NOTES:                                                                */
/*      To get the triangles (rather than the edges), modify the            */
/*      function out_triple as in the comments.                             */
/*                                                                          */
/*      The code is NOT threadsafe and does not clean up its memory.        */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "delaunay.h"
#include "macrorus.h"

#define DELETED -2
#define le 0
#define re 1

struct Freenode {
    struct Freenode *nextfree;
};

struct Freelist {
    struct Freenode *head;
    int nodesize;
};

struct Point {
    double x, y;
};

/* structure used both for sites and for vertices */

struct Site {
    struct Point coord;
    int sitenbr;
    int refcnt;
};


struct Edge {
    double a, b, c;
    struct Site *ep[2];
    struct Site *reg[2];
    int edgenbr;
};

struct Halfedge {
    struct Halfedge *ELleft, *ELright;
    struct Edge *ELedge;
    int ELrefcnt;
    char ELpm;
    struct Site *vertex;
    double ystar;
    struct Halfedge *PQnext;
};

typedef struct intptr {
    int this;
    struct intptr *next;
} intptr;

typedef struct delaunaydat {
    int nvertices;
    int PQmin;
    int PQcount;
    int PQhashsize;
    int ELhashsize;
    int siteidx;
    struct Halfedge *PQhash;
    double xmin, xmax, ymin, ymax, deltax, deltay;
    int nedges;
    struct Freelist efl;
    struct Freelist hfl;
    struct Halfedge *ELleftend, *ELrightend;
    struct Halfedge **ELhash;
    struct Site *sites;
    int nsites;
    int sqrt_nsites;
    struct Freelist sfl;
    struct Site *bottomsite;
    intptr **table;
    CCptrworld intptr_world;
    int total_alloc;
} delaunaydat;


CC_PTRWORLD_LIST_ROUTINES (intptr, int, intptralloc, intptr_bulk_alloc,
        intptrfree, intptr_listadd, intptr_listfree)
CC_PTRWORLD_LEAKS_ROUTINE (intptr, intptr_check_leaks, this, int)


static int
    put_in_table (delaunaydat *dd, int i, int j, int *added);



static void
    PQinsert (delaunaydat *dd, struct Halfedge *he, struct Site *v,
        double offset),
    PQdelete (delaunaydat *dd, struct Halfedge *he),
    PQinitialize (delaunaydat *dd),
    ELinitialize (delaunaydat *dd),
    ELinsert (struct Halfedge *lb, struct Halfedge *new),
    ELdelete (struct Halfedge *he),
    freeinit (struct Freelist *fl, int size),
    makefree (struct Freenode *curr, struct Freelist *fl),
    geominit (delaunaydat *dd),
    endpoint (delaunaydat *dd, struct Edge *e, int lr, struct Site *s),
    deref (delaunaydat *dd, struct Site *v),
    ref (struct Site *v),
    makevertex (delaunaydat *dd, struct Site *v);

static char
    *getfree (delaunaydat *dd, struct Freelist *fl),
    *vor_myalloc (delaunaydat *dd, unsigned n);

static int
    out_triple (delaunaydat *dd, struct Site *s1, struct Site *s2,
        struct Site *s3, int *ntotal),
    scomp (const void *s1_, const void *s2_),
    set_up_sites (delaunaydat *dd, int ncount, double *x, double *y),
    PQbucket (delaunaydat *dd, struct Halfedge *he),
    PQempty (delaunaydat *dd),
    right_of (struct Halfedge *el, struct Point *p);

static double
    dist (struct Site *s, struct Site *t);

static struct Site
    *nextsite (delaunaydat *dd),
    *leftreg (delaunaydat *dd, struct Halfedge *he),
    *rightreg (delaunaydat *dd, struct Halfedge *he),
    *intersect (delaunaydat *dd, struct Halfedge *el1, struct Halfedge *el2);

static struct Edge
    *bisect (delaunaydat *dd, struct Site *s1, struct Site *s2);

static struct Point
    PQ_min (delaunaydat *dd);

static struct Halfedge
    *PQextractmin (delaunaydat *dd),
    *HEcreate (delaunaydat *dd, struct Edge *e, int pm),
    *ELgethash (delaunaydat *dd, int b),
    *ELleftbnd (delaunaydat *dd, struct Point *p),
    *ELright (struct Halfedge *he),
    *ELleft (struct Halfedge *he);



int CCedgegen_delaunay (int ncount, CCdatagroup *dat, int wantlist,
        int *ecount, int **elist)
{
    struct Site *newsite, *bot, *top, *temp, *p;
    struct Site *v;
    struct Point newintstar;
    int i, pm;
    struct Halfedge *lbnd, *rbnd, *llbnd, *rrbnd, *bisector;
    struct Edge *e;
    int ntotal = 0;
    int rval = 0;
    intptr *ip, *ipnext;
    int onlist, total;
    int norm;
    delaunaydat dd;

    CCptrworld_init (&dd.intptr_world);
    dd.total_alloc = 0;
    
    dd.siteidx = 0;

    if (wantlist) {
        *ecount = 0;
        *elist = (int *) NULL;
    }

    CCutil_dat_getnorm (dat, &norm);
    if ((norm & CC_NORM_SIZE_BITS) != CC_D2_NORM_SIZE) {
        printf ("Cannot compute Delaunay triangulation with norm %d\n",
                 norm);
        fflush (stdout);
        return 0;
    }

    dd.table = CC_SAFE_MALLOC (ncount, intptr *);
    if (!dd.table) {
        rval = 1;
        goto CLEANUP;
    }
    for (i = 0; i < ncount; i++)
        dd.table[i] = (intptr *) NULL;

    freeinit (&dd.sfl, sizeof *dd.sites);
    if (set_up_sites (&dd, ncount, dat->x, dat->y)) {
        fprintf (stderr, "set_up_sites failed\n");
        rval = 1;
        goto CLEANUP;
    }
    geominit (&dd);

    PQinitialize (&dd);
    dd.bottomsite = nextsite (&dd);
    ELinitialize (&dd);

    newsite = nextsite (&dd);
    while (1) {
        if (!PQempty (&dd))
            newintstar = PQ_min (&dd);

        if (newsite != (struct Site *) NULL
            && (PQempty (&dd)
                || newsite->coord.y < newintstar.y
                || (newsite->coord.y == newintstar.y
                    && newsite->coord.x < newintstar.x))) {
            /* new site is smallest */
            lbnd = ELleftbnd (&dd, &(newsite->coord));
            rbnd = ELright (lbnd);
            bot = rightreg (&dd, lbnd);
            e = bisect (&dd, bot, newsite);
            bisector = HEcreate (&dd, e, le);
            ELinsert (lbnd, bisector);
            if ((p = intersect (&dd, lbnd, bisector)) != (struct Site *) NULL) {
                PQdelete (&dd, lbnd);
                PQinsert (&dd, lbnd, p, dist (p, newsite));
            }
            lbnd = bisector;
            bisector = HEcreate (&dd, e, re);
            ELinsert (lbnd, bisector);
            if ((p = intersect (&dd, bisector, rbnd)) != (struct Site *) NULL) {
                PQinsert (&dd, bisector, p, dist (p, newsite));
            }
            newsite = nextsite (&dd);
        } else if (!PQempty (&dd)) {
            /* intersection is smallest */
            lbnd = PQextractmin (&dd);
            llbnd = ELleft (lbnd);
            rbnd = ELright (lbnd);
            rrbnd = ELright (rbnd);
            bot = leftreg (&dd, lbnd);
            top = rightreg (&dd, rbnd);
            if (wantlist)
                out_triple (&dd, bot, top, rightreg (&dd, lbnd), &ntotal);
            v = lbnd->vertex;
            makevertex (&dd, v);
            endpoint (&dd, lbnd->ELedge, lbnd->ELpm, v);
            endpoint (&dd, rbnd->ELedge, rbnd->ELpm, v);
            ELdelete (lbnd);
            PQdelete (&dd, rbnd);
            ELdelete (rbnd);
            pm = le;
            if (bot->coord.y > top->coord.y) {
                temp = bot;
                bot = top;
                top = temp;
                pm = re;
            }
            e = bisect (&dd, bot, top);
            bisector = HEcreate (&dd, e, pm);
            ELinsert (llbnd, bisector);
            endpoint (&dd, e, re - pm, v);
            deref (&dd, v);
            if ((p = intersect (&dd, llbnd, bisector)) != (struct Site *) NULL) {
                PQdelete (&dd, llbnd);
                PQinsert (&dd, llbnd, p, dist (p, bot));
            }
            if ((p = intersect (&dd, bisector, rrbnd)) != (struct Site *) NULL) {
                PQinsert (&dd, bisector, p, dist (p, bot));
            }
        } else
            break;
    }

    if (wantlist) {
        int j = 0;
        *elist = CC_SAFE_MALLOC (2 * ntotal, int);
        if (!(*elist)) {
            rval = 1;
            goto CLEANUP;
        }
        *ecount = ntotal;
        for (i = 0; i < ncount; i++) {
            for (ip = dd.table[i]; ip; ip = ipnext) {
                ipnext =  ip->next;
                (*elist)[j++] = i;
                (*elist)[j++] = ip->this;
                intptrfree (&dd.intptr_world, ip);
            }
            dd.table[i] = (intptr *) NULL;
        }
    } else {
        for (i = 0; i < ncount; i++) {
            intptr_listfree (&dd.intptr_world, dd.table[i]);
            dd.table[i] = (intptr *) NULL;
        }
    }

    if (intptr_check_leaks (&dd.intptr_world, &total, &onlist)) {
        fprintf (stderr, "WARNING: %d outstanding intptrs in delaunay\n",
                 total - onlist);
    }

CLEANUP:

    CCptrworld_delete (&dd.intptr_world);
    if (dd.table)
        CC_FREE (dd.table, intptr *);
    if (dd.sites)
        CC_FREE (dd.sites, struct Site);

    return rval;
}

static int put_in_table (delaunaydat *dd, int i, int j, int *added)
{
    intptr *ip;

    if (j < i) {
        int temp;
        CC_SWAP(i, j, temp);
    }

    for (ip = dd->table[i]; ip; ip = ip->next) {
        if (ip->this == j) {
            *added = 0;
            return 0;
        }
    }
    if (intptr_listadd (&dd->table[i], j, &dd->intptr_world)) {
        *added = 0;
        return 1;
    }
    *added = 1;
    return 0;
}

/* sort sites on y, then x, coord */

static int scomp (const void *s1_, const void *s2_)
{
    const struct Point *s1,*s2;
    s1=(const struct Point*) s1_;
    s2=(const struct Point*) s2_;

    if (s1->y < s2->y)
        return (-1);
    if (s1->y > s2->y)
        return (1);
    if (s1->x < s2->x)
        return (-1);
    if (s1->x > s2->x)
        return (1);
    return (0);
}

static int set_up_sites (delaunaydat *dd, int ncount, double *x, double *y)
{
    int i;

    dd->nsites = ncount;
    dd->sites = CC_SAFE_MALLOC (ncount, struct Site);
    if (!dd->sites) {
        fprintf (stderr, "out of memory in set_up_sites\n");
        fflush (stdout);
        return 1;
    }

    for (i = 0; i < ncount; i++) {
        dd->sites[i].coord.x = x[i];
        dd->sites[i].coord.y = y[i];
        dd->sites[i].sitenbr = i;
        dd->sites[i].refcnt = 0;
    }
    qsort ((void *)dd->sites, (size_t)dd->nsites, sizeof *dd->sites, scomp);
    dd->xmin = dd->sites[0].coord.x;
    dd->xmax = dd->sites[0].coord.x;
    for (i = 1; i < dd->nsites; i += 1) {
        if (dd->sites[i].coord.x < dd->xmin)
            dd->xmin = dd->sites[i].coord.x;
        if (dd->sites[i].coord.x > dd->xmax)
            dd->xmax = dd->sites[i].coord.x;
    }
    dd->ymin = dd->sites[0].coord.y;
    dd->ymax = dd->sites[dd->nsites - 1].coord.y;

    return 0;
}

static int out_triple (delaunaydat *dd, struct Site *s1, struct Site *s2,
        struct Site *s3, int *ntotal)
{
    int added;
/*
    to just print the triangles, use:

    printf ("%d %d %d\n", s1->sitenbr, s2->sitenbr, s3->sitenbr);
*/

    if (put_in_table (dd, s1->sitenbr, s2->sitenbr, &added)) {
        fprintf (stderr, "put_in_table failed\n");
        return 1;
    }
    *ntotal += added;
    if (put_in_table (dd, s2->sitenbr, s3->sitenbr, &added)) {
        fprintf (stderr, "put_in_table failed\n");
        return 1;
    }
    *ntotal += added;
    if (put_in_table (dd, s3->sitenbr, s1->sitenbr, &added)) {
        fprintf (stderr, "put_in_table failed\n");
        return 1;
    }
    *ntotal += added;
    return 0;
}

static struct Site *nextsite (delaunaydat *dd)
{
    struct Site *s;

    if (dd->siteidx < dd->nsites) {
        s = &dd->sites[dd->siteidx];
        dd->siteidx += 1;
        return (s);
    } else
        return ((struct Site *) NULL);
}


/******************* edgelist.c  ***********************/

static void ELinitialize (delaunaydat *dd)
{
    int i;

    freeinit (&dd->hfl, sizeof **dd->ELhash);
    dd->ELhashsize = 2 * dd->sqrt_nsites;
    dd->ELhash = (struct Halfedge **) vor_myalloc (dd, (unsigned) (sizeof (*dd->ELhash) * dd->ELhashsize));
    for (i = 0; i < dd->ELhashsize; i += 1)
        dd->ELhash[i] = (struct Halfedge *) NULL;
    dd->ELleftend = HEcreate (dd, (struct Edge *) NULL, 0);
    dd->ELrightend = HEcreate (dd, (struct Edge *) NULL, 0);
    dd->ELleftend->ELleft = (struct Halfedge *) NULL;
    dd->ELleftend->ELright = dd->ELrightend;
    dd->ELrightend->ELleft = dd->ELleftend;
    dd->ELrightend->ELright = (struct Halfedge *) NULL;
    dd->ELhash[0] = dd->ELleftend;
    dd->ELhash[dd->ELhashsize - 1] = dd->ELrightend;
}

static struct Halfedge *HEcreate (delaunaydat *dd, struct Edge *e, int pm)
{
    struct Halfedge *answer;

    answer = (struct Halfedge *) getfree (dd, &dd->hfl);
    answer->ELedge = e;
    answer->ELpm = pm;
    answer->PQnext = (struct Halfedge *) NULL;
    answer->vertex = (struct Site *) NULL;
    answer->ELrefcnt = 0;
    return (answer);
}

static void ELinsert (struct Halfedge *lb, struct Halfedge *new)
{
    new->ELleft = lb;
    new->ELright = lb->ELright;
    (lb->ELright)->ELleft = new;
    lb->ELright = new;
}

/* Get entry from hash table, pruning any deleted nodes */

static struct Halfedge *ELgethash (delaunaydat *dd, int b)
{
    struct Halfedge *he;

    if (b < 0 || b >= dd->ELhashsize)
        return ((struct Halfedge *) NULL);
    he = dd->ELhash[b];
    if (he == (struct Halfedge *) NULL ||
        he->ELedge != (struct Edge *) (size_t) DELETED)
        return (he);

    /* Hash table points to deleted half edge.  Patch as necessary. */
    dd->ELhash[b] = (struct Halfedge *) NULL;
    if ((he->ELrefcnt -= 1) == 0)
        makefree ((struct Freenode *) he, &dd->hfl);
    return ((struct Halfedge *) NULL);
}

static struct Halfedge *ELleftbnd (delaunaydat *dd, struct Point *p)
{
    int i, bucket;
    struct Halfedge *he;

    /* Use hash table to get close to desired halfedge */
    bucket = (p->x - dd->xmin) / dd->deltax * dd->ELhashsize;
    if (bucket < 0)
        bucket = 0;
    if (bucket >= dd->ELhashsize)
        bucket = dd->ELhashsize - 1;
    he = ELgethash (dd, bucket);
    if (he == (struct Halfedge *) NULL) {
        for (i = 1; 1; i += 1) {
            if ((he = ELgethash (dd, bucket - i)) != (struct Halfedge *) NULL)
                break;
            if ((he = ELgethash (dd, bucket + i)) != (struct Halfedge *) NULL)
                break;
        }
    }
    /* Now search linear list of halfedges for the corect one */
    if (he == dd->ELleftend || (he != dd->ELrightend && right_of (he, p))) {
        do {
            he = he->ELright;
        } while (he != dd->ELrightend && right_of (he, p));
        he = he->ELleft;
    } else
        do {
            he = he->ELleft;
        } while (he != dd->ELleftend && !right_of (he, p));

    /* Update hash table and reference counts */
    if (bucket > 0 && bucket < dd->ELhashsize - 1) {
        if (dd->ELhash[bucket] != (struct Halfedge *) NULL)
            dd->ELhash[bucket]->ELrefcnt -= 1;
        dd->ELhash[bucket] = he;
        dd->ELhash[bucket]->ELrefcnt += 1;
    };
    return (he);
}


/* This delete routine can't reclaim node, since pointers from hash table may
 * be present.   */

static void ELdelete (struct Halfedge *he)
{
    (he->ELleft)->ELright = he->ELright;
    (he->ELright)->ELleft = he->ELleft;
    he->ELedge = (struct Edge *) (size_t) DELETED;
}


static struct Halfedge *ELright (struct Halfedge *he)
{
    return (he->ELright);
}

static struct Halfedge *ELleft (struct Halfedge *he)
{
    return (he->ELleft);
}


static struct Site *leftreg (delaunaydat *dd, struct Halfedge *he)
{
    if (he->ELedge == (struct Edge *) NULL)
        return (dd->bottomsite);
    return (he->ELpm == le ?
            he->ELedge->reg[le] : he->ELedge->reg[re]);
}

static struct Site *rightreg (delaunaydat *dd, struct Halfedge *he)
{
    if (he->ELedge == (struct Edge *) NULL)
        return (dd->bottomsite);
    return (he->ELpm == le ?
            he->ELedge->reg[re] : he->ELedge->reg[le]);
}


/***************** geometry.c ********************/

static void geominit (delaunaydat *dd)
{
    struct Edge e;
    double sn;

    freeinit (&dd->efl, sizeof e);
    dd->nvertices = 0;
    dd->nedges = 0;
    sn = dd->nsites + 4;
    dd->sqrt_nsites = sqrt (sn);
    dd->deltay = dd->ymax - dd->ymin;
    dd->deltax = dd->xmax - dd->xmin;
}

static struct Edge *bisect (delaunaydat *dd, struct Site *s1, struct Site *s2)
{
    double dx, dy, adx, ady;
    struct Edge *newedge;

    newedge = (struct Edge *) getfree (dd, &dd->efl);

    newedge->reg[0] = s1;
    newedge->reg[1] = s2;
    ref (s1);
    ref (s2);
    newedge->ep[0] = (struct Site *) NULL;
    newedge->ep[1] = (struct Site *) NULL;

    dx = s2->coord.x - s1->coord.x;
    dy = s2->coord.y - s1->coord.y;
    adx = dx > 0 ? dx : -dx;
    ady = dy > 0 ? dy : -dy;
    newedge->c = s1->coord.x * dx + s1->coord.y * dy +
                        (dx * dx + dy * dy) * 0.5;
    if (adx > ady) {
        newedge->a = 1.0;
        newedge->b = dy / dx;
        newedge->c /= dx;
    } else {
        newedge->b = 1.0;
        newedge->a = dx / dy;
        newedge->c /= dy;
    };

    newedge->edgenbr = dd->nedges;
    dd->nedges += 1;
    return (newedge);
}

static struct Site *intersect (delaunaydat *dd, struct Halfedge *el1,
        struct Halfedge *el2)
{
    struct Edge *e1, *e2, *e;
    struct Halfedge *el;
    double d, xint, yint;
    int right_of_site;
    struct Site *v;

    e1 = el1->ELedge;
    e2 = el2->ELedge;
    if (e1 == (struct Edge *) NULL || e2 == (struct Edge *) NULL)
        return ((struct Site *) NULL);
    if (e1->reg[1] == e2->reg[1])
        return ((struct Site *) NULL);

    d = e1->a * e2->b - e1->b * e2->a;
    if (-1.0e-10 < d && d < 1.0e-10)
        return ((struct Site *) NULL);

    xint = (e1->c * e2->b - e2->c * e1->b) / d;
    yint = (e2->c * e1->a - e1->c * e2->a) / d;

    if ((e1->reg[1]->coord.y < e2->reg[1]->coord.y) ||
        (e1->reg[1]->coord.y == e2->reg[1]->coord.y &&
         e1->reg[1]->coord.x < e2->reg[1]->coord.x)) {
        el = el1;
        e = e1;
    } else {
        el = el2;
        e = e2;
    };
    right_of_site = xint >= e->reg[1]->coord.x;
    if ((right_of_site && el->ELpm == le) ||
        (!right_of_site && el->ELpm == re))
        return ((struct Site *) NULL);

    v = (struct Site *) getfree (dd, &dd->sfl);
    v->refcnt = 0;
    v->coord.x = xint;
    v->coord.y = yint;
    return (v);
}

/* returns 1 if p is to right of halfedge e */

static int right_of (struct Halfedge *el, struct Point *p)
{
    struct Edge *e;
    struct Site *topsite;
    int right_of_site, above, fast;
    double dxp, dyp, dxs, t1, t2, t3, yl;

    e = el->ELedge;
    topsite = e->reg[1];
    right_of_site = p->x > topsite->coord.x;
    if (right_of_site && el->ELpm == le)
        return (1);
    if (!right_of_site && el->ELpm == re)
        return (0);

    if (e->a == 1.0) {
        dyp = p->y - topsite->coord.y;
        dxp = p->x - topsite->coord.x;
        fast = 0;
        if ((!right_of_site & (e->b < 0.0)) |
            (right_of_site & (e->b >= 0.0))) {
            above = dyp >= e->b * dxp;
            fast = above;
        } else {
            above = p->x + p->y * e->b > e->c;
            if (e->b < 0.0)
                above = !above;
            if (!above)
                fast = 1;
        };
        if (!fast) {
            dxs = topsite->coord.x - (e->reg[0])->coord.x;
            above = e->b * (dxp * dxp - dyp * dyp) <
                dxs * dyp * (1.0 + 2.0 * dxp / dxs + e->b * e->b);
            if (e->b < 0.0)
                above = !above;
        };
    } else {                    /* e->b==1.0 */
        yl = e->c - e->a * p->x;
        t1 = p->y - yl;
        t2 = p->x - topsite->coord.x;
        t3 = yl - topsite->coord.y;
        above = t1 * t1 > t2 * t2 + t3 * t3;
    };
    return (el->ELpm == le ? above : !above);
}


static void endpoint (delaunaydat *dd, struct Edge *e, int lr, struct Site *s)
{
    e->ep[lr] = s;
    ref (s);
    if (e->ep[re - lr] == (struct Site *) NULL)
        return;
    deref (dd, e->reg[le]);
    deref (dd, e->reg[re]);
    makefree ((struct Freenode *) e, &dd->efl);
}


static double dist (struct Site *s, struct Site *t)
{
    double dx, dy;

    dx = s->coord.x - t->coord.x;
    dy = s->coord.y - t->coord.y;
    return (sqrt (dx * dx + dy * dy));
}


static void makevertex (delaunaydat *dd, struct Site *v)
{
    v->sitenbr = dd->nvertices;
    dd->nvertices += 1;
}


static void deref (delaunaydat *dd, struct Site *v)
{
    v->refcnt -= 1;
    if (v->refcnt == 0)
        makefree ((struct Freenode *) v, &dd->sfl);
}

static void ref (struct Site *v)
{
    v->refcnt += 1;
}


/*********************** heap.c *******************/

static void PQinsert (delaunaydat *dd, struct Halfedge *he, struct Site *v,
        double offset)
{
    struct Halfedge *last, *next;

    he->vertex = v;
    ref (v);
    he->ystar = v->coord.y + offset;
    last = &dd->PQhash[PQbucket (dd, he)];
    while ((next = last->PQnext) != (struct Halfedge *) NULL &&
           (he->ystar > next->ystar ||
        (he->ystar == next->ystar && v->coord.x > next->vertex->coord.x))) {
        last = next;
    }
    he->PQnext = last->PQnext;
    last->PQnext = he;
    dd->PQcount += 1;
}

static void PQdelete (delaunaydat *dd, struct Halfedge *he)
{
    struct Halfedge *last;

    if (he->vertex != (struct Site *) NULL) {
        last = &dd->PQhash[PQbucket (dd, he)];
        while (last->PQnext != he)
            last = last->PQnext;
        last->PQnext = he->PQnext;
        dd->PQcount -= 1;
        deref (dd, he->vertex);
        he->vertex = (struct Site *) NULL;
    }
}

static int PQbucket (delaunaydat *dd, struct Halfedge *he)
{
    int bucket;

    bucket = (he->ystar - dd->ymin) / dd->deltay * dd->PQhashsize;
    if (bucket < 0)
        bucket = 0;
    if (bucket >= dd->PQhashsize)
        bucket = dd->PQhashsize - 1;
    if (bucket < dd->PQmin)
        dd->PQmin = bucket;
    return (bucket);
}



static int PQempty (delaunaydat *dd)
{
    return (dd->PQcount == 0);
}


static struct Point PQ_min (delaunaydat *dd)
{
    struct Point answer;

    while (dd->PQhash[dd->PQmin].PQnext == (struct Halfedge *) NULL) {
        dd->PQmin += 1;
    }
    answer.x = dd->PQhash[dd->PQmin].PQnext->vertex->coord.x;
    answer.y = dd->PQhash[dd->PQmin].PQnext->ystar;
    return (answer);
}

static struct Halfedge *PQextractmin (delaunaydat *dd)
{
    struct Halfedge *curr;

    curr = dd->PQhash[dd->PQmin].PQnext;
    dd->PQhash[dd->PQmin].PQnext = curr->PQnext;
    dd->PQcount -= 1;
    return (curr);
}


static void PQinitialize (delaunaydat *dd)
{
    int i;

    dd->PQcount = 0;
    dd->PQmin = 0;
    dd->PQhashsize = 4 * dd->sqrt_nsites;
    dd->PQhash = (struct Halfedge *) vor_myalloc (dd, (unsigned) (dd->PQhashsize * sizeof (*dd->PQhash)));
    for (i = 0; i < dd->PQhashsize; i += 1)
        dd->PQhash[i].PQnext = (struct Halfedge *) NULL;
}


/***************** memory.c *****************/

static void freeinit (struct Freelist *fl, int size)
{
    fl->head = (struct Freenode *) NULL;
    fl->nodesize = size;
}

static char *getfree (delaunaydat *dd, struct Freelist *fl)
{
    int i;
    struct Freenode *t;

    if (fl->head == (struct Freenode *) NULL) {
        t = (struct Freenode *) vor_myalloc (dd, (unsigned) (dd->sqrt_nsites * fl->nodesize));
        for (i = 0; i < dd->sqrt_nsites; i += 1)
            makefree ((struct Freenode *) ((char *) t + i * fl->nodesize), fl);
    };
    t = fl->head;
    fl->head = (fl->head)->nextfree;
    return ((char *) t);
}

static void makefree (struct Freenode *curr, struct Freelist *fl)
{
    curr->nextfree = fl->head;
    fl->head = curr;
}

static char *vor_myalloc (delaunaydat *dd, unsigned n)
{
    char *t;

    if ((t = malloc (n)) == (char *) 0) {
        fprintf (stderr, "Out of memory processing %d (%d bytes in use)\n",
                 dd->siteidx, dd->total_alloc);
        exit (1);
    };
    dd->total_alloc += n;
    return (t);
}
