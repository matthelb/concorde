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
/*      TOUR MAINTANENCE ROUTINES FOR LIN-KERNIGHAN - Two-level Tree        */
/*                                                                          */
/*                             TSP CODE                                     */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: May 2, 1995                                                       */
/*        May 1, 1998 (bico)                                                */
/*                                                                          */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CClinkern_flipper_init (CClk_flipper *f, int ncount, int *cyc)      */
/*    initializes flipper to an initial cycle given in cyc.                 */
/*    returns 0 on success, nonzero on failure.                             */
/*                                                                          */
/*  void CClinkern_flipper_cycle (CClk_flipper *F, int *x)                  */
/*    places the current cycle in x.                                        */
/*                                                                          */
/*  void CClinkern_flipper_finish (CClk_flipper *F)                         */
/*    frees up space allocated by CClinkern_flipper_init.                   */
/*    every CClinkern_flipper_init should lead to a                         */
/*    CClinkern_flipper_finish call.                                        */
/*                                                                          */
/*  int CClinkern_flipper_next (CClk_flipper *f, int x)                     */
/*    returns the successor to x in the current cycle.                      */
/*                                                                          */
/*  int CClinkern_flipper_prev (CClk_flipper *f, int x)                     */
/*    returns the predecessor of x in the current cycle.                    */
/*                                                                          */
/*  void CClinkern_flipper_flip (CClk_flipper *F, int x, int y)             */
/*    flips the portion of the cycle from x to y (inclusive).               */
/*                                                                          */
/*  int CClinkern_flipper_sequence (CClk_flipper *f, int * x, int y,        */
/*      int z)                                                              */
/*    returns 1 if xyz occur as an increasing subsequence of the cycle,     */
/*    returns 0 otherwise.                                                  */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* NOTES:                                                                   */
/*       This is desribed in the paper "Data structures for traveling       */
/*   salesman" by Fredman, Johnson, McGeoch, and Ostheimer.                 */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "linkern.h"

/****************************************************************************/
/*                                                                          */
/* TWO-LEVEL TREES:                                                         */
/*                                                                          */
/*     1. Uses the "groupsize" approach described in the paper.             */
/*                                                                          */
/****************************************************************************/

#define GROUPSIZE_FACTOR 0.50
#define SEGMENT_SPLIT_CUTOFF 0.30 


static void
    same_segment_flip (CClk_flipper *F, CClk_childnode *a, CClk_childnode *b),
    consecutive_segment_flip (CClk_flipper *F, CClk_parentnode *a,
        CClk_parentnode *b),
    segment_split (CClk_flipper *F, CClk_parentnode *p, CClk_childnode *aprev,
        CClk_childnode *a, int left_or_right),
    init_flipper (CClk_flipper *Fl),
    free_flipper (CClk_flipper *Fl);

static int
    build_flipper (CClk_flipper *Fl, int ncount);


#define SAME_SEGMENT(a, b)                                                  \
     (a->parent == b->parent &&                                             \
      ((!((F->reversed)^(a->parent->rev)) && a->id <= b->id) ||             \
       (((F->reversed)^(a->parent->rev)) && a->id >= b->id)))


int CClinkern_flipper_init (CClk_flipper *F, int ncount, int *cyc)
{
    int i, j, cind, remain;
    int rval = 0;
    CClk_childnode *c, *cprev;
    CClk_parentnode *p;

    init_flipper (F);
    rval = build_flipper (F, ncount); 
    if (rval) {
        fprintf (stderr, "build_flipper failed\n"); goto CLEANUP;
    }

    remain = ncount;
    i = 0;
    j = 2 * F->groupsize;
    while (remain >= j) {
        F->parents[i].size = F->groupsize;  
        remain -= F->groupsize;
        i++;
    }
    if (remain > F->groupsize) {
        F->parents[i].size = remain / 2;
        remain -= (remain / 2);
        i++;
    }
    F->parents[i].size = remain;
    i++;

    if (i != F->nsegments) {
        fprintf (stderr, "seg count is wrong\n");
        rval = 1; goto CLEANUP;
    }

    c = &(F->children[cyc[ncount - 1]]);
    for (i=0,p=F->parents,cind=0; i < F->nsegments; p++,i++) {
        p->id = i;
        p->rev = 0;
        p->ends[0] = &(F->children[cyc[cind]]);
        for (j = p->size; j > 0; j--) {
            cprev = c;
            c = &(F->children[cyc[cind]]);
            c->id = cind;
            c->name = cyc[cind];
            c->parent = p;
            c->adj[0] = cprev;
            cprev->adj[1] = c;
            cind++;
        }
        p->ends[1] = c;
        p->adj[0] = p - 1;
        p->adj[1] = p + 1;
    }
    F->parents[0].adj[0] = &(F->parents[F->nsegments - 1]);
    F->parents[F->nsegments - 1].adj[1] = &(F->parents[0]);

CLEANUP:
   
    if (rval) {
        free_flipper (F);
    }
    return rval;
}

void CClinkern_flipper_cycle (CClk_flipper *F, int *x)
{
    CClk_childnode *c, *start;
    int k = 0;

    start = &(F->children[0]);
    c = start->adj[!((F->reversed)^(start->parent->rev))];

    x[k++] = start->name;
    while (c != start) {
        x[k++] = c->name;
        c = c->adj[!((F->reversed)^(c->parent->rev))];
    }
}

void CClinkern_flipper_finish (CClk_flipper *F)
{
    free_flipper (F);
}

int CClinkern_flipper_next (CClk_flipper *F, int x)
{
    return
      F->children[x].adj[!((F->reversed)^(F->children[x].parent->rev))]->name;
}

int CClinkern_flipper_prev (CClk_flipper *F, int x)
{
    return
      F->children[x].adj[(F->reversed)^(F->children[x].parent->rev)]->name;
}

void CClinkern_flipper_flip (CClk_flipper *F, int x, int y)
{
    CClk_childnode *xc = &(F->children[x]);
    CClk_childnode *yc = &(F->children[y]);

    if (SAME_SEGMENT (xc, yc)) {
        if (xc != yc) {
            same_segment_flip (F, xc, yc);
        }
    } else {
        int xdir = ((F->reversed)^(xc->parent->rev));
        int ydir = ((F->reversed)^(yc->parent->rev));
        CClk_childnode *xprev = xc->adj[xdir];
        CClk_childnode *ynext = yc->adj[!ydir];
        if (SAME_SEGMENT (ynext, xprev)) {
            if (ynext != xprev) {
                same_segment_flip (F, ynext, xprev);
            }
            (F->reversed) ^= 1;
        } else {
            int side;
            if (xc->parent->ends[xdir] == xc &&
                yc->parent->ends[!ydir] == yc) {
                if (F->reversed)
                    side = xc->parent->id - yc->parent->id;
                else
                    side = yc->parent->id - xc->parent->id;
                if (side < 0) 
                    side += F->nsegments;
                if (side < F->nsegments / 2) {
                    consecutive_segment_flip (F, xc->parent, yc->parent);
                } else {
                    consecutive_segment_flip (F,yc->parent->adj[!F->reversed],
                                                xc->parent->adj[F->reversed]);
                    (F->reversed) ^= 1;
                }
            } else {
                if (xprev->parent == xc->parent) {
                    segment_split (F, xc->parent, xprev, xc, 0);
                    if (SAME_SEGMENT (xc, yc)) {
                        if (xc != yc)
                            same_segment_flip (F, xc, yc);
                        return;
                    } else if (SAME_SEGMENT (ynext, xprev)) {
                        if (ynext != xprev) {
                            same_segment_flip (F, ynext, xprev);
                        }
                        (F->reversed) ^= 1;
                        return;
                    }
                }
                if (ynext->parent == yc->parent) {
                    segment_split (F, yc->parent, yc, ynext, 0);
                    if (SAME_SEGMENT (xc, yc)) {
                        if (xc != yc)
                            same_segment_flip (F, xc, yc);
                        return;
                    } else if (SAME_SEGMENT (ynext, xprev)) {
                        if (ynext != xprev) {
                            same_segment_flip (F, ynext, xprev);
                        }
                        (F->reversed) ^= 1;
                        return;
                    }
                }
                if (F->reversed)
                    side = xc->parent->id - yc->parent->id;
                else
                    side = yc->parent->id - xc->parent->id;
                if (side < 0)
                    side += F->nsegments;
                if (side < F->nsegments / 2) {
                    consecutive_segment_flip (F, xc->parent, yc->parent);
                } else {
                    consecutive_segment_flip(F, yc->parent->adj[!F->reversed],
                                                xc->parent->adj[F->reversed]);
                    (F->reversed) ^= 1;
                }

            }
        }
    }
}

static void same_segment_flip (CClk_flipper *F, CClk_childnode *a,
        CClk_childnode *b)
{
    CClk_parentnode *parent = a->parent;
    int dir = ((F->reversed)^(parent->rev)); 
    CClk_childnode *aprev = a->adj[dir];
    CClk_childnode *bnext = b->adj[!dir];
    CClk_childnode *c, *cnext;

    if ((dir && a->id - b->id > F->split_cutoff) ||
       (!dir && b->id - a->id > F->split_cutoff)) {
        if (aprev->parent == parent) 
            segment_split (F, parent, aprev, a, 1);
        if (bnext->parent == parent)
            segment_split (F, parent, b, bnext, 2);
        aprev->adj[!((F->reversed)^(aprev->parent->rev))] = b;
        bnext->adj[(F->reversed)^(bnext->parent->rev)] = a;
        a->adj[dir] = bnext;
        b->adj[!dir] = aprev;
        parent->rev ^= 1;
        return;
    }

    if (dir) {
        int id = a->id;
        aprev->adj[!((F->reversed)^(aprev->parent->rev))] = b;
        bnext->adj[(F->reversed)^(bnext->parent->rev)] = a;
        cnext = b->adj[1];
        b->adj[1] = aprev;
        b->adj[0] = cnext;
        b->id = id--;
        c = cnext;
        while (c != a) {
            cnext = c->adj[1];
            c->adj[1] = c->adj[0];
            c->adj[0] = cnext;
            c->id = id--;
            c = cnext;
        }
        a->adj[1] = a->adj[0];
        a->adj[0] = bnext;
        a->id = id;
        if (parent->ends[1] == a)
            parent->ends[1] = b;
        if (parent->ends[0] == b)
            parent->ends[0] = a;
    } else {
        int id = a->id;
        aprev->adj[!((F->reversed)^(aprev->parent->rev))] = b;
        bnext->adj[(F->reversed)^(bnext->parent->rev)] = a;
        c = b->adj[0];
        b->adj[0] = aprev;
        b->adj[1] = c;
        b->id = id++;
        while (c != a) {
            cnext = c->adj[0];
            c->adj[0] = c->adj[1];
            c->adj[1] = cnext;
            c->id = id++;
            c = cnext;
        }
        a->adj[0] = a->adj[1];
        a->adj[1] = bnext;
        a->id = id;
        if (parent->ends[0] == a)
            parent->ends[0] = b;
        if (parent->ends[1] == b)
            parent->ends[1] = a;
    }
}

static void consecutive_segment_flip (CClk_flipper *F, CClk_parentnode *a,
        CClk_parentnode *b)
{
    CClk_parentnode *aprev = a->adj[F->reversed];
    CClk_parentnode *bnext = b->adj[!F->reversed];
    CClk_parentnode *c, *cnext;
    CClk_childnode *achild = a->ends[(F->reversed)^(a->rev)];
    CClk_childnode *bchild = b->ends[!((F->reversed)^(b->rev))];
    CClk_childnode *childprev, *childnext;
    int id = a->id;

    if (F->reversed) {
        childprev = achild->adj[!a->rev];
        childnext = bchild->adj[b->rev];
        childprev->adj[childprev->parent->rev] = bchild;
        childnext->adj[!childnext->parent->rev] = achild;
        bchild->adj[b->rev] = childprev;
        achild->adj[!a->rev] = childnext;
  
        aprev->adj[0] = b;
        bnext->adj[1] = a;
        c = b->adj[1];
        b->adj[1] = aprev;
        b->adj[0] = c;
        b->id = id--;
        b->rev ^= 1;
        while (c != a) {
            cnext = c->adj[1];
            c->adj[1] = c->adj[0];
            c->adj[0] = cnext;
            c->id = id--;
            c->rev ^= 1;
            c = cnext;
        }
        a->adj[1] = a->adj[0];
        a->adj[0] = bnext;
        a->id = id;
        a->rev ^= 1;
    } else {
        childprev = achild->adj[a->rev];
        childnext = bchild->adj[!b->rev];
        childprev->adj[!childprev->parent->rev] = bchild;
        childnext->adj[childnext->parent->rev] = achild;
        bchild->adj[!b->rev] = childprev;
        achild->adj[a->rev] = childnext;

        aprev->adj[1] = b;
        bnext->adj[0] = a;
        c = b->adj[0];
        b->adj[0] = aprev;
        b->adj[1] = c;
        b->id = id++;
        b->rev ^= 1;
        while (c != a) {
            cnext = c->adj[0];
            c->adj[0] = c->adj[1];
            c->adj[1] = cnext;
            c->id = id++;
            c->rev ^= 1;
            c = cnext;
        }
        a->adj[0] = a->adj[1];
        a->adj[1] = bnext;
        a->id = id;
        a->rev ^= 1;
    }
}

/* split between a and aprev */

static void segment_split (CClk_flipper *F, CClk_parentnode *p,
        CClk_childnode *aprev, CClk_childnode *a, int left_or_right)
{
    int side;
    int dir = ((F->reversed)^(p->rev));
    int id;
    CClk_parentnode *pnext;
    CClk_childnode *b, *bnext;

    if (dir) side = p->ends[1]->id - aprev->id + 1;
    else     side = aprev->id - p->ends[0]->id + 1;

    if ((left_or_right == 0 && side <= p->size / 2) || left_or_right == 1) {
        pnext = p->adj[F->reversed];
        pnext->size += side;
        p->size -= side;
        if (pnext->rev == p->rev) {
            b = pnext->ends[!dir];
            id = b->id;
            if (dir) {
                do {
                    b = b->adj[0];
                    b->id = --id;
                    b->parent = pnext;
                } while (b != aprev);
            } else {
                do {
                    b = b->adj[1];
                    b->id = ++id;
                    b->parent = pnext;
                } while (b != aprev);
            }
            pnext->ends[!dir] = aprev;
            p->ends[dir] = a;
        } else {
            b = pnext->ends[dir];
            id = b->id;
            if (!dir) {
                bnext = b->adj[0];
                do {
                    b = bnext;
                    b->id = --id;
                    b->parent = pnext;
                    bnext = b->adj[1];
                    b->adj[1] = b->adj[0];
                    b->adj[0] = bnext;
                } while (b != aprev);
            } else {
                bnext = b->adj[1];
                do {
                    b = bnext;
                    b->id = ++id;
                    b->parent = pnext;
                    bnext = b->adj[0];
                    b->adj[0] = b->adj[1];
                    b->adj[1] = bnext;
                } while (b != aprev);
            }
            pnext->ends[dir] = aprev;
            p->ends[dir] = a;
        }
    } else {
        pnext = p->adj[!F->reversed];
        pnext->size += (p->size - side);
        p->size = side;
        if (pnext->rev == p->rev) {
            b = pnext->ends[dir];
            id = b->id;
            if (dir) {
                do {
                    b = b->adj[1];
                    b->id = ++id;
                    b->parent = pnext;
                } while (b != a);
            } else {
                do {
                    b = b->adj[0];
                    b->id = --id;
                    b->parent = pnext;
                } while (b != a);
            }
            pnext->ends[dir] = a;
            p->ends[!dir] = aprev;
        } else {
            b = pnext->ends[!dir];
            id = b->id;
            if (!dir) {
                bnext = b->adj[1];
                do {
                    b = bnext;
                    b->id = ++id;
                    b->parent = pnext;
                    bnext = b->adj[0];
                    b->adj[0] = b->adj[1];
                    b->adj[1] = bnext;
                } while (b != a);
            } else {
                bnext = b->adj[0];
                do {
                    b = bnext;
                    b->id = --id;
                    b->parent = pnext;
                    bnext = b->adj[1];
                    b->adj[1] = b->adj[0];
                    b->adj[0] = bnext;
                } while (b != a);
            }
            pnext->ends[!dir] = a;
            p->ends[!dir] = aprev;
        }
    }
}

int CClinkern_flipper_sequence (CClk_flipper *F, int x, int y, int z)
{
    CClk_childnode *a = &(F->children[x]);
    CClk_childnode *b = &(F->children[y]);
    CClk_childnode *c = &(F->children[z]);
    CClk_parentnode *pa = a->parent;
    CClk_parentnode *pb = b->parent;
    CClk_parentnode *pc = c->parent;

    if (pa == pb) {
        if (pa == pc) {
            if ((F->reversed)^(pa->rev)) {
                if (a->id >= b->id) {
                    return (b->id >= c->id || c->id >= a->id); 
                } else {
                    return (b->id >= c->id && c->id >= a->id);
                }
            } else {
                if (a->id <= b->id) {
                    return (b->id <= c->id || c->id <= a->id); 
                } else {
                    return (b->id <= c->id && c->id <= a->id);
                }
            }
        } else {
            if ((F->reversed)^(pa->rev)) {
                return (a->id >= b->id);
            } else {
                return (a->id <= b->id);
            }
        }
    } else if (pa == pc) {
            if ((F->reversed)^(pa->rev)) {
                return (a->id <= c->id);
            } else {
                return (a->id >= c->id);
            }
    } else if (pb == pc) {
            if ((F->reversed)^(pb->rev)) {
                return (b->id >= c->id);
            } else {
                return (b->id <= c->id);
            }
    } else {
        if (F->reversed) {
            if (pa->id >= pb->id) {
                return (pb->id >= pc->id || pc->id >= pa->id); 
            } else {
                return (pb->id >= pc->id && pc->id >= pa->id);
            }
        } else {
            if (pa->id <= pb->id) {
                return (pb->id <= pc->id || pc->id <= pa->id); 
            } else {
                return (pb->id <= pc->id && pc->id <= pa->id);
            }
        }
    }
}


static void init_flipper (CClk_flipper *Fl)
{
    Fl->parents = (CClk_parentnode *) NULL;
    Fl->children = (CClk_childnode *) NULL;
    Fl->reversed = 0;
    Fl->nsegments = 0;
    Fl->groupsize = 100;
    Fl->split_cutoff = 100;
}

static void free_flipper (CClk_flipper *Fl)
{
    if (Fl) {
        CC_IFFREE (Fl->parents, CClk_parentnode);
        CC_IFFREE (Fl->children, CClk_childnode);
        Fl->reversed = 0;
        Fl->nsegments = 0;
        Fl->groupsize = 0;
        Fl->split_cutoff = 0;
    }
}

static int build_flipper (CClk_flipper *Fl, int ncount)
{
    int rval = 0;

    Fl->reversed = 0;
    Fl->groupsize = (int) (sqrt ((double) ncount) * GROUPSIZE_FACTOR);
    Fl->nsegments =  (ncount + Fl->groupsize - 1) / Fl->groupsize;
    Fl->split_cutoff = Fl->groupsize * SEGMENT_SPLIT_CUTOFF;

    Fl->parents = CC_SAFE_MALLOC (Fl->nsegments, CClk_parentnode);
    Fl->children = CC_SAFE_MALLOC (ncount + 1, CClk_childnode); 
                 /* The +1 will stop a purify burp later */
    if (Fl->parents  == (CClk_parentnode *) NULL ||
        Fl->children == (CClk_childnode *) NULL) {
        fprintf (stderr, "out of memory in build_flipper\n");
        rval = 1; goto CLEANUP;
    }

CLEANUP:

    if (rval) {
        free_flipper (Fl);
    }
    return rval;
}


