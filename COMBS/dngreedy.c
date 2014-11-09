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
/*                      GREEDY CUT HEURISTICS                               */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: September 4, 1997                                                 */
/*                                                                          */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCcombs_greedy_cut (CC_GCgraph *g, int *setsize, int *set,          */
/*      int mark_fixed, int forced_moves, int bad_moves,                    */
/*      int fixed_moves, int *moves_done, double *cut_val)                  */
/*    RETURNS a set of nodes of having a small coboundary.                  */
/*     -g is the graph to grow in, including marks                          */
/*     -setsize is the size of the set (modified by greedy)                 */
/*     -set is the set of vertices (modified by greedy). This should        */
/*      point to an array of size g->ncount (to handle any size set).       */
/*     -mark_fixed, nodes with this mark are fixed (forbidden from          */
/*      entering or leaving the set)                                        */
/*     -forced_moves, make at least this many moves, even if they make      */
/*      the cut worse.  If a move is forced (ie, would not have happened    */
/*      otherwise) then the result is fixed.                                */
/*     -bad_moves, make at least this many moves that make the cut worse.   */
/*      The result of the bad move is fixed.                                */
/*     -fixed_moves, make at least this many moves, even if they make       */
/*      the cut worse.  Fix the result of each of these moves.              */
/*     -moves_done, if non-NULL, should point to an int which will be       */
/*      set to 1 if all of the forced moves, bad moves, and fixed moves     */
/*      were made, and 0 if the cut ran out of available nodes before       */
/*      making the required moves.                                          */
/*     -cut_val is the set's cut value (set by greedy)                      */
/*    NOTES:                                                                */
/*     -assumes that status is 0 for every node in g, and restores          */
/*      this condition on successful exit                                   */
/*     -modifies the node fields qhandle, flow, setloc, setdeg, status      */
/*                                                                          */
/*  int CCcombs_GC_build_graph (CC_GCgraph *G, int ncount, int ecount,      */
/*      int *elist, double *x)                                              */
/*    BUILDS the graph corresponding to the edge lists, using x for the     */
/*     weights on the edges.                                                */
/*                                                                          */
/*  void CCcombs_GC_init_graph (CC_GCgraph *G)                              */
/*    INITIALIZES the fields of the graph struct to NULL.                   */
/*                                                                          */
/*  void CCcombs_GC_free_graph (CC_GCgraph *G)                              */
/*    FREES the fields of the graph struct.                                 */
/*                                                                          */
/*    NOTES: differences from Denis Naddef's greedy grower:                 */
/*     1.  Denis only adds nodes to the step.  Greedy can delete nodes.     */
/*     2.  Denis selects any improving move.  Greedy selects the most       */
/*         improving move.                                                  */
/*     3.  Denis has several methods of selecting a worsening move.         */
/*         Greedy selects the least worsening move.                         */
/*     4.  Denis continues until a cut < target is found, or until          */
/*         badmoves worsening moves have been made.  The best cut seen      */
/*         during the process is returned.  Greedy continues until          */
/*         at least forced_moves moves have been made, and at least         */
/*         bad_moves worsening moves have been made, and until there are    */
/*         no more improving moves available.  Greedy returns the final     */
/*         cut.                                                             */
/*     5.  Greedy supports fixed nodes.                                     */
/*                                                                          */
/****************************************************************************/

/*
 * qhandle == data for priority queue
 * flow(v) == x(S:v)
 * setloc(v) == -1 ==> v not in S
 * setloc(v) >= 0 ==> set[setloc(v)] = v
 * setdeg(v) = |(S:v)| - used to free a node when it ceases to be in the
 *     neighborhood of S
 * status: 1 = initialized (in neighborhood)
 *         2 = last move tied delete
 *         4 = fixed
 */

#include "machdefs.h"
#include "util.h"
#include "combs.h"

#define GCUT_EPS (1e-6)

#define STAT_INITIALIZED 1
#define STAT_LASTDELTIED 2
#define STAT_FIXED       4

typedef struct greedy_data {
    CCpriority q;
    double cut_val;
    int *set;
    int mark_fixed;
    int setsize;
} greedy_data;


static int
    init_greedy (CC_GCgraph *g, greedy_data *gd),
    update_node (greedy_data *gd, CC_GCnode *n),
    add_node (greedy_data *gd, CC_GCgraph *g, CC_GCnode *n, int do_update),
    del_node (greedy_data *gd, CC_GCgraph *g, CC_GCnode *n, int do_update);

static void
    init_node (greedy_data *gd, CC_GCnode *n),
    clear_greedy (CC_GCgraph *g, greedy_data *gd);


int CCcombs_greedy_cut (CC_GCgraph *g, int *setsize, int *set, int mark_fixed,
        int forced_moves, int bad_moves, int fixed_moves, int *moves_done,
        double *cut_val)
{
    greedy_data gd;
    CC_GCnode *n;
    double delta;
    int rval;

    rval = CCutil_priority_init (&gd.q, *setsize);
    if (rval) {
        fprintf (stderr, "CCutil_priority_init failed\n"); goto CLEANUP;
    }
    gd.set = set;
    gd.setsize = *setsize;
    gd.mark_fixed = mark_fixed;
    if (moves_done != (int *) NULL) *moves_done = 0;

    rval = init_greedy (g, &gd);
    if (rval) {
        fprintf (stderr, "init_greedy failed\n"); goto CLEANUP;
    }

    while ((n = (CC_GCnode *) CCutil_priority_findmin (&gd.q, &delta)) !=
           (CC_GCnode *) NULL) {
        if (delta > 0.0) {
            if (forced_moves <= 0 && bad_moves <= 0 && fixed_moves <= 0) {
                break;
            }
            n->status |= STAT_FIXED;
            bad_moves--;
        }
        if (fixed_moves > 0) n->status |= STAT_FIXED;
        forced_moves--;
        fixed_moves--;
        if (n->setloc >= 0) {
            rval = del_node (&gd, g, n, 1);
            if (rval) {
                fprintf (stderr, "del_node failed\n"); goto CLEANUP;
            }
        } else {
            rval = add_node (&gd, g, n, 1);
            if (rval) {
                fprintf (stderr, "add_node failed\n"); goto CLEANUP;
            }
        }
    }

    *setsize = gd.setsize;
    *cut_val = gd.cut_val;
    if (moves_done != (int *) NULL && forced_moves <= 0 &&
        bad_moves <= 0 && fixed_moves <= 0) {
        *moves_done = 1;
    }

    rval = 0;

  CLEANUP:
    clear_greedy (g, &gd);
    CCutil_priority_free (&gd.q);
    return rval;
}

static int init_greedy (CC_GCgraph *g, greedy_data *gd)
{
    int i;
    int j;
    int in_setsize = gd->setsize;
    CC_GCnode *nodelist = g->nodelist;
    CC_GCnode *n;
    int *set = gd->set;
    int rval;

    gd->setsize = 0;
    gd->cut_val = 0.0;

    for (i = 0; i < in_setsize; i++) {
        init_node (gd, &nodelist[set[i]]);
    }

    for (i = 0; i < in_setsize; i++) {
        rval = add_node (gd, g, &nodelist[set[i]], 0);
        if (rval) {
            fprintf (stderr, "add_node failed\n");
            return rval;
        }
    }

    for (i = 0; i < in_setsize; i++) {
        n = &nodelist[set[i]];
        rval = update_node (gd, n);
        if (rval) {
            fprintf (stderr, "update_node failed\n");
            return rval;
        }
        for (j = 0; j < n->deg; j++) {
            rval = update_node (gd, &nodelist[n->adj[j].to]);
            if (rval) {
                fprintf (stderr, "update_node failed\n");
                return rval;
            }
        }
    }
    return 0;
}

static void init_node (greedy_data *gd, CC_GCnode *n)
{
    n->setdeg = 0;
    n->flow = 0.0;
    n->setloc = -1;
    n->qhandle = -1;
    n->status = STAT_INITIALIZED;
    if (n->mark == gd->mark_fixed) {
        n->status |= STAT_FIXED;
    }
}

static void clear_greedy (CC_GCgraph *g, greedy_data *gd)
{
    int i;
    int j;
    CC_GCnode *n;
    CC_GCnode *nodelist = g->nodelist;
    int setsize = gd->setsize;
    int *set = gd->set;

    for (i = 0; i < setsize; i++) {
        n = &nodelist[set[i]];
        n->status = 0;
        for (j = 0; j < n->deg; j++) {
            nodelist[n->adj[j].to].status = 0;
        }
    }
}

static int update_node (greedy_data *gd, CC_GCnode *n)
{
    double delta;

    if (n->setdeg == 0 && n->setloc == -1) {
        if (n->qhandle >= 0) {
            CCutil_priority_delete (&(gd->q), n->qhandle);
            n->qhandle = -1;
        }
        n->status = 0;
    } else if (n->status & STAT_FIXED) {
        if (n->qhandle >= 0) {
            CCutil_priority_delete (&(gd->q), n->qhandle);
            n->qhandle = -1;
        }
    } else if ((n->status & STAT_LASTDELTIED) &&
               n->setloc < 0 &&
               n->flow >= 1.0 - GCUT_EPS &&
               n->flow <= 1.0 + GCUT_EPS) {
        /* tied adds undoing a tied delete are disabled */
        if (n->qhandle >= 0) {
            CCutil_priority_delete (&(gd->q), n->qhandle);
            n->qhandle = -1;
        }
    } else {
        if (n->setloc < 0) delta = 1.0 - n->flow;
        else delta = n->flow - 1.0;

        /* improving moves and tied adds have higher priority than
           tied deletes and unimproving moves.  delta is shifted so
           values < 0.0 are "good" and values > 0.0 are "bad" */

        if (delta > GCUT_EPS || (delta >= -GCUT_EPS && n->setloc >= 0)) {
            delta += 0.5;
        } else {
            delta -= 0.5;
        }

        if (n->qhandle < 0) {
            n->qhandle = CCutil_priority_insert (&(gd->q), (void *) n, delta);
            if (n->qhandle < 0) {
                fprintf (stderr, "CCutil_priority_insert failed\n");
                return 1;
            }
        } else {
            CCutil_priority_changekey (&(gd->q), n->qhandle, delta);
        }
    }
    return 0;
}

static int add_node (greedy_data *gd, CC_GCgraph *g, CC_GCnode *n,
        int do_update)
{
    CC_GCnode *nodelist = g->nodelist;
    CC_GCedge *adj = n->adj;
    int deg = n->deg;
    int i;
    int m;
    double w;
    int rval;

    n->setloc = gd->setsize;
    gd->set[gd->setsize] = (int) (n - nodelist);
    gd->setsize++;

    for (i=0; i<deg; i++) {
        m = adj[i].to;
        w = adj[i].weight;
        if ((nodelist[m].status & STAT_INITIALIZED) == 0) {
            init_node (gd, &nodelist[m]);
        }
        if (nodelist[m].setloc < 0) {
            gd->cut_val += w;
        } else {
            gd->cut_val -= w;
        }
        nodelist[m].flow += w;
        nodelist[m].setdeg++;
        if (do_update) {
            rval = update_node (gd, &nodelist[m]);
            if (rval) {
                fprintf (stderr, "update_node failed\n");
                return rval;
            }
        }
    }
    if (do_update) {
        n->status &= ~STAT_LASTDELTIED;
        rval = update_node (gd, n);
        if (rval) {
            fprintf (stderr, "update_node failed\n");
            return rval;
        }
    }
    return 0;
}

static int del_node (greedy_data *gd, CC_GCgraph *g, CC_GCnode *n,
        int do_update)
{
    CC_GCnode *nodelist = g->nodelist;
    CC_GCedge *adj = n->adj;
    int deg = n->deg;
    int i;
    int m;
    double w;
    int rval;

    gd->setsize--;
    m = gd->set[gd->setsize];
    gd->set[n->setloc] = m;
    nodelist[m].setloc = n->setloc;
    n->setloc = -1;

    for (i=0; i<deg; i++) {
        m = adj[i].to;
        w = adj[i].weight;
        if ((nodelist[m].status & STAT_INITIALIZED) == 0) {
            init_node (gd, &nodelist[m]);
        }
        if (nodelist[m].setloc < 0) {
            gd->cut_val -= w;
        } else {
            gd->cut_val += w;
        }
        nodelist[m].flow -= w;
        nodelist[m].setdeg--;
        if (do_update) {
            rval = update_node (gd, &nodelist[m]);
            if (rval) {
                fprintf (stderr, "update_node failed\n");
                return rval;
            }
        }
    }
    if (do_update) {
        if (1.0 - GCUT_EPS <= n->flow && n->flow <= 1.0 + GCUT_EPS) {
            n->status |= STAT_LASTDELTIED;
        } else {
            n->status &= ~STAT_LASTDELTIED;
        }
        rval = update_node (gd, n);
        if (rval) {
            fprintf (stderr, "update_node failed\n");
            return rval;
        }
    }
    return 0;
}

int CCcombs_GC_build_graph (CC_GCgraph *G, int ncount, int ecount, int *elist,
        double *x)
{
    int rval = 0;
    CC_GCnode *n;
    CC_GCedge *e;
    int i, j;

    CCcombs_GC_init_graph (G);

    if (ncount) {
        G->nodelist = CC_SAFE_MALLOC (ncount, CC_GCnode);
        if (G->nodelist == (CC_GCnode *) NULL) {
            fprintf (stderr, "out of memory in CCcombs_GC_build_graph\n");
            rval = 1; goto CLEANUP;
        }
    }
    if (ecount) {
        G->edgespace = CC_SAFE_MALLOC (2 * ecount, CC_GCedge);
        if (G->edgespace == (CC_GCedge *) NULL) {
            fprintf (stderr, "out of memory in CCcombs_GC_build_graph\n");
            CC_IFFREE (G->nodelist, CC_GCnode);
            rval = 1; goto CLEANUP;
        }
    }

    n = G->nodelist;
    for (i = 0; i < ncount; i++) {
        n[i].deg     = 0;
        n[i].adj     = (CC_GCedge *) NULL;
        n[i].mark    = 0;
        n[i].qhandle = 0;
        n[i].flow    = 0.0;
        n[i].setloc  = 0;
        n[i].setdeg  = 0;
        n[i].status  = 0;
    }

    for (i = 0; i < ecount; i++) {
        n[elist[2*i]].deg++;
        n[elist[2*i+1]].deg++;
    }

    e = G->edgespace;
    for (i = 0; i < ncount; i++) {
        n[i].adj = e;
        e += n[i].deg;
        n[i].deg = 0;
    }
    for (i = 0; i < ecount; i++) {
        j = elist[2*i];
        n[j].adj[n[j].deg].to     =  elist[2*i+1];
        n[j].adj[n[j].deg].weight = x[i];
        n[j].deg++;
        j = elist[2*i+1];
        n[j].adj[n[j].deg].to     =  elist[2*i];
        n[j].adj[n[j].deg].weight = x[i];
        n[j].deg++;
    }

    G->ncount = ncount;
    G->ecount = ecount;

CLEANUP:

    return rval;
}

void CCcombs_GC_init_graph (CC_GCgraph *G)
{
    if (G) {
        G->nodelist  = (CC_GCnode *) NULL;
        G->edgespace = (CC_GCedge *) NULL;
        G->ncount    = 0;
        G->ecount    = 0;
    }
}

void CCcombs_GC_free_graph (CC_GCgraph *G)
{
    if (G) {
        CC_IFFREE (G->nodelist, CC_GCnode);
        CC_IFFREE (G->edgespace, CC_GCedge);
    }
}
