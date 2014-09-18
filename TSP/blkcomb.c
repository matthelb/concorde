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
/*               Find Combs Based on Block Decomposition                    */
/*                                                                          */
/*                              TSP CODE                                    */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: October 2, 1997                                                   */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCtsp_block_combs (CCtsp_lpcut_in **cuts, int *cutcount,            */
/*      int ncount, int ecount, int *elist, double *x, int silent)          */
/*    -ncount is the number of nodes                                        */
/*    -ecount is the number of edges                                        */
/*    -elist is the edge list in node node format                           */
/*    -x is an lp solution vector                                           */
/*    -silent turns off output if set to a nonzero value                    */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "cut.h"
#include "tsp.h"
#include "combs.h"

#define DOUBLE_HAND_CUTOFF 10
#define CUT_BATCH          1000   

#define ONE_TRIANGLE       1
#define TRIANGLE_SQUARE    2
#define ONE_SQUARE         3
#define TIGHT_TRIANGLE     4
#define TIGHT_SQUARE       5

typedef struct intptr {
    int this;
    struct intptr *next;
} intptr;

typedef struct node {
    struct edge *adj;
    int degree;
    int magiclabel;
    int cutmark;
} node;

typedef struct edge {
    double x;
    int to;
} edge;

typedef struct graph {
    node *nodelist;
    edge *adjspace;
    int   ncount;
    int   ecount;
    int   magicnum;
} graph;



static int
    block_combs_work (CCtsp_lpgraph *L, double *origx, int ncount,
        int ecount, int *elist, double *x, CC_SRKexpinfo *expand,
        CCtsp_lpcut_in **cuts, int *cutcount, double *maxviol, double *czeit),
    greedy_doubleblock_comb (CC_GCgraph *G, CCtsp_lpgraph *L, graph *H,
        double *x, int bcnt1, int *block1, int bcnt2, int *block2, int *set,
        CCtsp_lpcut_in **newcut, double *viol, CC_SRKexpinfo *expand,
        int marker),
    greedy_block_comb (CC_GCgraph *G, CCtsp_lpgraph *L, graph *H,
        double *x, int bcnt, int *block, int *set, CCtsp_lpcut_in **newcut,
        double *viol, CC_SRKexpinfo *expand, int marker),
    buildgraph (graph *G, int ncount, int ecount, int *elist, double *x),
    shrink_ones (CC_SRKgraph *S, int ncount, int ecount, int *elist,
        double *dlen, int *oncount, int *oecount, int **olist, double **olen,
        CC_SRKexpinfo *expand),
    shrink_extra (CC_SRKgraph *G, int *oncount, int *oecount, int **olist,
        double **olen, CC_SRKexpinfo *expand, int level),
    shrink_olaf_new (CC_SRKgraph *G, int ncount, int ecount, int *elist,
        double *dlen, int *oncount, int *oecount, int **olist, double **olen,
        CC_SRKexpinfo *expand, int level),
    shrink_olaf (CC_SRKgraph *G, int *oncount, int *oecount, int **olist,
        double **olen, CC_SRKexpinfo *expand, int level);

static void
    initgraph (graph *G),
    freegraph (graph *G),
    blockcomb_free_world (CCptrworld *intptr_world);


CC_PTRWORLD_LIST_ROUTINES (intptr, int, intptralloc, intptr_bulkalloc,
        intptrfree, intptr_listadd, intptr_listfree)
CC_PTRWORLD_LEAKS_ROUTINE (intptr, intptr_check_leaks, this, int)

int CCtsp_block_combs (CCtsp_lpcut_in **cuts, int *cutcount, int ncount,
        int ecount, int *elist, double *x, int silent)
{
    int rval = 0;
    int oncount = 0;
    int oecount = 0;
    int kcount  = 0;
    int level, oldcount;
    int *oelist = (int *) NULL;
    double szeit = CCutil_zeit ();
    double czeit, viol;
    double *ox = (double *) NULL;
    CCtsp_lpgraph L;
    CC_SRKgraph S;
    CC_SRKexpinfo expand;

    CCtsp_init_lpgraph_struct (&L);
    CCcut_SRK_init_graph (&S);
    CCcut_SRK_init_expinfo (&expand);

    if (cutcount) *cutcount = 0;

    if (!silent) {
        printf ("Block Combs\n"); fflush (stdout);
    }

    rval = CCtsp_build_lpgraph (&L, ncount, ecount, elist, (int *) NULL);
    if (rval) {
        fprintf (stderr, "CCtsp_build_lpgraph failed\n"); goto CLEANUP;
    }
    rval = CCtsp_build_lpadj (&L, 0, ecount);
    if (rval) {
        fprintf (stderr, "CCtsp_build_lpadj failed\n"); goto CLEANUP;
    }

    rval = shrink_ones (&S, ncount, ecount, elist, x, &oncount, &oecount,
                        &oelist, &ox, &expand);
    if (rval) {
        fprintf (stderr, "shrink_ones failed\n"); goto CLEANUP;
    }

    rval = block_combs_work (&L, x, oncount, oecount, oelist, ox, &expand,
                             cuts, &kcount, &viol, &czeit);
    if (rval) {
        fprintf (stderr, "block_combs_work failed\n"); goto CLEANUP;
    }
    if (cutcount) *cutcount += kcount;
    if (!silent) {
        printf ("  Type 0: %d cuts  %.4f max violation  %.2f seconds\n",
                             kcount, viol, czeit);
        fflush (stdout);
    }

    CC_IFFREE (oelist, int);
    CC_IFFREE (ox, double);
    CCcut_SRK_free_expinfo (&expand);

    oldcount = oncount;
    level = 1;
    while (kcount < CUT_BATCH && level <= 5) {
        CCcut_SRK_increment_marker (&S);
        rval = shrink_extra (&S, &oncount, &oecount, &oelist, &ox, &expand,
                             level);
        if (rval) {
            fprintf (stderr, "shrink_extra failed\n"); goto CLEANUP;
        }

        if (oncount < (.9 * oldcount)) {
            rval = block_combs_work (&L, x, oncount, oecount, oelist, ox,
                                     &expand, cuts, &kcount, &viol, &czeit);
            if (rval) {
                fprintf (stderr, "block_combs_work failed\n"); goto CLEANUP;
            }
            if (cutcount) *cutcount += kcount;
            if (!silent) {
                printf ("  Type %d: %d cuts  %.4f max violation  %.2f seconds\n",
                                 level, kcount, viol, czeit);
                fflush (stdout);
            }
            oldcount = oncount;
        }
        CC_IFFREE (oelist, int);
        CC_IFFREE (ox, double);
        CCcut_SRK_free_expinfo (&expand);
        level++;
    }
    CCcut_SRK_free_graph (&S);

    level = 1;
    while (kcount < CUT_BATCH && level <= 4) {
        rval = shrink_olaf_new (&S, ncount, ecount, elist, x, &oncount,
                        &oecount, &oelist, &ox, &expand, level);
        if (rval) {
            fprintf (stderr, "shrink_olaf_new failed\n"); goto CLEANUP;
        }
        rval = block_combs_work (&L, x, oncount, oecount, oelist, ox,
                                 &expand, cuts, &kcount, &viol, &czeit);
        if (rval) {
            fprintf (stderr, "block_combs_work failed\n"); goto CLEANUP;
        }
        if (cutcount) *cutcount += kcount;
        if (!silent) {
            printf ("  Olaf %d: %d cuts  %.4f max violation  %.2f seconds\n",
                                 level, kcount, viol, czeit);
            fflush (stdout);
        }

        CC_IFFREE (oelist, int);
        CC_IFFREE (ox, double);
        CCcut_SRK_free_expinfo (&expand);
        CCcut_SRK_free_graph (&S);
        level++;
    }

    if (!silent) {
        printf ("  Total Time in block_combs: %.2f\n", CCutil_zeit () - szeit);
        fflush (stdout);
    }


CLEANUP:

    CCtsp_free_lpgraph (&L);
    CCcut_SRK_free_graph (&S);
    CCcut_SRK_free_expinfo (&expand);

    CC_IFFREE (oelist, int);
    CC_IFFREE (ox, double);

    return rval;
}

static int block_combs_work (CCtsp_lpgraph *L, double *origx, int ncount,
        int ecount, int *elist, double *x, CC_SRKexpinfo *expand,
        CCtsp_lpcut_in **cuts, int *cutcount, double *maxviol, double *czeit)
{
    int rval = 0;
    int kcut = 0;
    int nblocks = 0;
    int marker = 0;
    int i, j;
    int ncutnodes = 0;
    int *cutnodes = (int *) NULL;
    int *set = (int *) NULL;
    int *blockcnt = (int *) NULL;
    int **blocks = (int **) NULL;
    intptr **cutadj = (intptr **) NULL;
    CC_GCgraph G;
    graph H;
    double szeit, viol;
    CCptrworld intptr_world;

    CCcombs_GC_init_graph (&G);
    initgraph (&H);
    CCptrworld_init (&intptr_world);
    if (cutcount) *cutcount = 0;
    if (maxviol) *maxviol = 0.0;

    szeit = CCutil_zeit ();

    if (ncount == 0 || ecount == 0) goto CLEANUP;

    rval = CCcombs_find_blocks (ncount, ecount, elist, x,
             &nblocks, &blockcnt, &blocks, &ncutnodes, &cutnodes);
    if (rval) {
        fprintf (stderr, "CCcombs_find_blocks failed\n"); goto CLEANUP;
    }

    rval = CCcombs_GC_build_graph (&G, ncount, ecount, elist, x);
    if (rval) {
        fprintf (stderr, "CCcombs_GC_build_graph failed\n"); goto CLEANUP;
    }

    set = CC_SAFE_MALLOC (ncount, int);
    if (set == (int *) NULL) {
        fprintf (stderr, "out of memory in CCtsp_block_combs\n");
        rval = 1; goto CLEANUP;
    }

    rval = buildgraph (&H, ncount, ecount, elist, x);
    if (rval) {
        fprintf (stderr, "buildgraph failed\n"); goto CLEANUP;
    }

    for (i = 0; i < ncutnodes; i++) {
        H.nodelist[cutnodes[i]].cutmark = 1;
    }

    for (i = 0; i < nblocks; i++) {
        CCtsp_lpcut_in *newcut = (CCtsp_lpcut_in *) NULL;
        marker++;
        rval = greedy_block_comb (&G, L, &H, origx, blockcnt[i], blocks[i], set,
                                  &newcut, &viol, expand, marker);
        if (rval) {
            fprintf (stderr, "greedy_block_comb failed\n"); goto CLEANUP;
        }
        if (newcut) {
            if (viol > CCtsp_MIN_VIOL) {
                newcut->next = *cuts;
                *cuts = newcut;
                kcut++;
            } else {
                CCtsp_free_lpcut_in (newcut);
                CC_IFFREE (newcut, CCtsp_lpcut_in);
            }
            if (maxviol) {
                if (viol > *maxviol) *maxviol = viol;
            }
        }
    }


    /* Now try to build combs with union of two blocks as the handle */

    cutadj = CC_SAFE_MALLOC (H.ncount, intptr *);
    if (cutadj == (intptr **) NULL) {
        fprintf (stderr, "out of memory in CCtsp_block_combs\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < H.ncount; i++) {
        cutadj[i] = (intptr *) NULL;
    }
    for (i = 0; i < nblocks; i++) {
        for (j = 0; j < blockcnt[i]; j++) {
            if (H.nodelist[blocks[i][j]].cutmark == 1) {
                rval = intptr_listadd (&cutadj[blocks[i][j]], i,
                                       &intptr_world);
                if (rval) goto CLEANUP;
            }
        }
    }
    for (i = 0; i < ncutnodes; i++) {
        intptr *ip, *jp;
        int ct = 0;
        for (ip = cutadj[cutnodes[i]]; ip; ip = ip->next) {
            ct++;
        }
        if (ct > 1) {
            H.nodelist[cutnodes[i]].cutmark = 0;
            if (ct < DOUBLE_HAND_CUTOFF) {
                for (ip = cutadj[cutnodes[i]]; ip; ip = ip->next) {
                    for (jp = ip->next; jp; jp = jp->next) {
                        CCtsp_lpcut_in *newcut = (CCtsp_lpcut_in *) NULL;
                        marker++;
                        rval = greedy_doubleblock_comb (&G, L, &H, origx,
                                 blockcnt[ip->this], blocks[ip->this],
                                 blockcnt[jp->this], blocks[jp->this],
                                 set, &newcut, &viol, expand, marker);
                        if (rval) {
                            fprintf (stderr, "doubleblock_comb failed\n");
                            goto CLEANUP;
                        }
                        if (newcut) {
                            if (viol > CCtsp_MIN_VIOL) {
                                newcut->next = *cuts;
                                *cuts = newcut;
                                kcut++;
                            } else {
                                CCtsp_free_lpcut_in (newcut);
                                CC_IFFREE (newcut, CCtsp_lpcut_in);
                            }
                            if (maxviol) {
                                if (viol > *maxviol) *maxviol = viol;
                            }
                        }
                    }
                }
            } else {
                for (ip = cutadj[cutnodes[i]]; ip && ip->next; ip = ip->next) {
                    CCtsp_lpcut_in *newcut = (CCtsp_lpcut_in *) NULL;
                    marker++;
                    rval = greedy_doubleblock_comb (&G, L, &H, origx,
                             blockcnt[ip->this], blocks[ip->this],
                             blockcnt[ip->next->this], blocks[ip->next->this],
                             set, &newcut, &viol, expand, marker);
                    if (rval) {
                        fprintf (stderr, "greedy_doubleblock_comb failed\n");
                        goto CLEANUP;
                    }
                    if (newcut) {
                        if (viol > CCtsp_MIN_VIOL) {
                            newcut->next = *cuts;
                            *cuts = newcut;
                            kcut++;
                        } else {
                            CCtsp_free_lpcut_in (newcut);
                            CC_IFFREE (newcut, CCtsp_lpcut_in);
                        }
                        if (maxviol) {
                            if (viol > *maxviol) *maxviol = viol;
                        }
                    }
                }
            }
            H.nodelist[cutnodes[i]].cutmark = 1;
        }
    }
    if (cutcount) *cutcount = kcut;


CLEANUP:

    if (czeit) *czeit = CCutil_zeit () - szeit;

    CCcombs_GC_free_graph (&G);
    CC_IFFREE (cutnodes, int);
    for (i = 0; i < nblocks; i++) {
        CC_IFFREE (blocks[i], int);
    }
    CC_IFFREE (blocks, int *);
    CC_IFFREE (blockcnt, int);
    CC_IFFREE (set, int);
    if (cutadj) {
        for (i = 0; i < H.ncount; i++) {
            intptr_listfree (&intptr_world, cutadj[i]);
        }
        CC_FREE (cutadj, intptr *);
    }

    blockcomb_free_world (&intptr_world);
    freegraph (&H);

    return rval;
}

static int greedy_doubleblock_comb (CC_GCgraph *G, CCtsp_lpgraph *L, graph *H,
        double *x, int bcnt1, int *block1, int bcnt2, int *block2, int *set,
        CCtsp_lpcut_in **newcut, double *viol, CC_SRKexpinfo *expand,
        int marker)
{
    int i, ct;
    int rval = 0;
    int *dblock = (int *) NULL;

    ct = bcnt1 + bcnt2;
    if (ct) {
        dblock = CC_SAFE_MALLOC (ct, int);
        if (dblock == (int *) NULL) {
            fprintf (stderr, "out of memory in block_combs\n");
            rval = 0; goto CLEANUP;
        }
        ct = 0;
        H->magicnum++;
        for (i = 0; i < bcnt1; i++) {
            dblock[ct++] = block1[i];
            H->nodelist[block1[i]].magiclabel = H->magicnum;
        }
        for (i = 0; i < bcnt2; i++) {
            if (H->nodelist[block2[i]].magiclabel != H->magicnum) {
                dblock[ct++] = block2[i];
            }
        }
        rval = greedy_block_comb (G, L, H, x, ct, dblock, set,
                                  newcut, viol, expand, marker);
        if (rval) {
            fprintf (stderr, "greedy_block_comb failed\n");
            goto CLEANUP;
        }
    }

CLEANUP:

    CC_IFFREE (dblock, int);
    return rval;
}

static int greedy_block_comb (CC_GCgraph *G, CCtsp_lpgraph *L, graph *H,
        double *x, int bcnt, int *block, int *set, CCtsp_lpcut_in **newcut,
        double *viol, CC_SRKexpinfo *expand, int marker)
{
    int rval = 0;
    CCtsp_lpclique handle;
    CCtsp_lpclique *teeth = (CCtsp_lpclique *) NULL;
    CCtsp_lpclique **bigteeth = (CCtsp_lpclique **) NULL;
    int tcnt = 0;
    int cutcnt = 0;
    int i, j, setsize;
    double cutval, slack;
    int *wset = (int *) NULL;

    if (viol) *viol = 0.0;
    CCtsp_init_lpclique (&handle);

    if (bcnt < 3) return 0;

    wset = CC_SAFE_MALLOC (L->ncount, int);
    if (wset == (int *) NULL) {
        fprintf (stderr, "out of memory in greedy_block_comb\n");
        rval = 1; goto CLEANUP;
    }

    rval = CCtsp_shrunk_set_to_lpclique (bcnt, block, wset, expand, &handle);
    if (rval) {
        fprintf (stderr, "CCtsp_shrunk_set_to_lpclique failed\n"); goto CLEANUP;
    }
    for (i = 0; i < bcnt; i++) {
        G->nodelist[block[i]].mark = marker;
    }

    for (i = 0; i < bcnt; i++) {
        if (H->nodelist[block[i]].cutmark == 1) {
            cutcnt++;
        }
    }


    if (cutcnt) {
        teeth    = CC_SAFE_MALLOC (cutcnt, CCtsp_lpclique);
        bigteeth = CC_SAFE_MALLOC (cutcnt + 1, CCtsp_lpclique *);
        if (!teeth || !bigteeth) {
            fprintf (stderr, "out of memory in greedy_block_comb\n");
            rval = 1; goto CLEANUP;
        }
        for (i = 0; i < bcnt; i++) {
            if (H->nodelist[block[i]].cutmark == 1) {
                int wsetsize = 0;
                set[0] = block[i];
                setsize = 1;

                rval = CCcombs_greedy_cut (G, &setsize, set, marker, 2, 0, 2,
                                           (int *) NULL, &cutval);
                if (rval) {
                    fprintf (stderr, "CCcombs_greedy_cut failed\n");
                    goto CLEANUP;
                }

                if (setsize > 1) {
                    for (j = 0; j < setsize; j++) {
                        wsetsize += (expand->memindex[set[j] + 1] - 
                                     expand->memindex[set[j]]);
                    }
               
                    if (wsetsize >= 3) {
                        rval = CCtsp_shrunk_set_to_lpclique (setsize, set, wset,
                                           expand, &teeth[tcnt]);
                        if (rval) {
                            fprintf (stderr, "shrunk_set_to_lpclique failed\n");
                            goto CLEANUP;
                        }
                        bigteeth[tcnt + 1] = &teeth[tcnt];
                        tcnt++;
                        for (j = 0; j < setsize; j++) {
                            G->nodelist[set[j]].mark = marker;
                        }
                    }
                }
            }
        }
    }

    rval = CCtsp_teething_list (L, x, &handle, tcnt, bigteeth, newcut);
    if (rval) {
        fprintf (stderr, "CCtsp_teething_list failed\n"); goto CLEANUP;
    }  

    if (*newcut) {
        slack = CCtsp_cutprice (L, *newcut, x);
        if (viol) *viol = -slack;
    }

CLEANUP:

    CCtsp_free_lpclique (&handle);
    for (i = 0; i < tcnt; i++) {
         CCtsp_free_lpclique (&teeth[i]);
    }
    CC_IFFREE (teeth, CCtsp_lpclique);
    CC_IFFREE (bigteeth, CCtsp_lpclique *);
    CC_IFFREE (wset, int);

    return rval;
}

static int buildgraph (graph *G, int ncount, int ecount, int *elist, double *x)
{
    int rval = 0;
    node *n;
    edge *e;
    int i, j;

    if (ncount) {
        G->nodelist = CC_SAFE_MALLOC (ncount, node);
        if (G->nodelist == (node *) NULL) {
            fprintf (stderr, "out of memory in buildgraph\n");
            rval = 1; goto CLEANUP;
        }
    }
    if (ecount) {
        G->adjspace = CC_SAFE_MALLOC (2 * ecount, edge);
        if (G->adjspace == (edge *) NULL) {
            fprintf (stderr, "out of memory in buildgraph\n");
            CC_IFFREE (G->nodelist, node);
            rval = 1; goto CLEANUP;
        }
    }

    n = G->nodelist;
    for (i = 0; i < ncount; i++) {
        n[i].degree = 0;
        n[i].adj    = (edge *) NULL;
        n[i].magiclabel = 0;
        n[i].cutmark = 0;
    }

    for (i = 0; i < ecount; i++) {
        n[elist[2*i]].degree++;
        n[elist[2*i+1]].degree++;
    }

    e = G->adjspace;
    for (i = 0; i < ncount; i++) {
        n[i].adj = e;
        e += n[i].degree;
        n[i].degree = 0;
    }
    for (i = 0; i < ecount; i++) {
        j = elist[2*i];
        n[j].adj[n[j].degree].to =  elist[2*i+1];
        n[j].adj[n[j].degree].x  = x[i];
        n[j].degree++;
        j = elist[2*i+1];
        n[j].adj[n[j].degree].to =  elist[2*i];
        n[j].adj[n[j].degree].x  = x[i];
        n[j].degree++;
    }

    G->ncount = ncount;
    G->ecount = ecount;
    G->magicnum = 0;

CLEANUP:

    return rval;
}

static void initgraph (graph *G)
{
    if (G) {
        G->nodelist = (node *) NULL;
        G->adjspace = (edge *) NULL;
        G->ncount   = 0;
        G->ecount   = 0;
        G->magicnum = 0;
    }
}

static void freegraph (graph *G)
{
    if (G) {
        CC_IFFREE (G->nodelist, node);
        CC_IFFREE (G->adjspace, edge);
    }
}

static void blockcomb_free_world (CCptrworld *intptr_world)
{
    int total, onlist;

    if (intptr_check_leaks (intptr_world, &total, &onlist)) {
        fprintf (stderr, "WARNING: %d outstanding intptrs\n",
                 total - onlist);
    }
    CCptrworld_delete (intptr_world);
}

static int shrink_ones (CC_SRKgraph *G, int ncount, int ecount, int *elist,
        double *dlen, int *oncount, int *oecount, int **olist, double **olen,
        CC_SRKexpinfo *expand)
{
    int k;

    if (CCcut_SRK_buildgraph (G, ncount, ecount, elist, dlen)) {
        fprintf (stderr, "buildgraph failed in shrink_ones\n");
        return 1;
    }
    CCcut_SRK_increment_marker (G);

    if (CCcut_SRK_defluff (G)) {
        fprintf (stderr, "CCcut_SRK_defluff failed in shrink_ones\n");
        return 1;
    }

    CCcut_SRK_identify_paths_to_edges (G, &k, 0);
    if (CCcut_SRK_grab_edges (G, oncount, oecount, olist, olen, expand)) {
        fprintf (stderr, "grab edges failed in shrink_ones\n");
        return 1;
    }

    return 0;
}

static int shrink_extra (CC_SRKgraph *G, int *oncount, int *oecount,
        int **olist, double **olen, CC_SRKexpinfo *expand, int level)
{
    int rval = 0;

    switch (level) {
    case ONE_TRIANGLE:
        CCcut_SRK_identify_one_triangles (G, (int *) NULL, (CC_SRKnode *) NULL,
                                          0.001, 2.0, 1);
        break;
    case TIGHT_TRIANGLE:
        CCcut_SRK_identify_tight_triangles (G, (int *) NULL, 2.0, 1);
        break;
    case TIGHT_SQUARE:
        CCcut_SRK_identify_tight_squares (G, (int *) NULL, 2.0, 1);
        break;
    case TRIANGLE_SQUARE:
        CCcut_SRK_identify_triangle_square (G, (int *) NULL, 0.001, 1);
        break;
    case ONE_SQUARE:
        CCcut_SRK_identify_one_square (G, (int *) NULL, 0.001, 1.0, 1);
        break;
    }

    rval = CCcut_SRK_grab_edges (G, oncount, oecount, olist, olen, expand);
    if (rval) {
        fprintf (stderr, "grab edges failed\n");
        goto CLEANUP;
    }

CLEANUP:

    return rval;
}

static int shrink_olaf_new (CC_SRKgraph *G, int ncount, int ecount, int *elist,
        double *dlen, int *oncount, int *oecount, int **olist, double **olen,
        CC_SRKexpinfo *expand, int level)
{
    int rval = 0;

    if (CCcut_SRK_buildgraph (G, ncount, ecount, elist, dlen)) {
        fprintf (stderr, "buildgraph failed in shrink_ones\n");
        return 1;
    }
    CCcut_SRK_increment_marker (G);

    rval = CCcut_SRK_defluff (G);
    if (rval) {
        fprintf (stderr, "CCcut_SRK_defluff failed in shrink_ones\n");
        goto CLEANUP;
    }

    rval = shrink_olaf (G, oncount, oecount, olist, olen, expand, level);
    if (rval) {
        fprintf (stderr, "shrink_olaf failed\n");
        goto CLEANUP;
    }

CLEANUP:

    return rval;
}

static int shrink_olaf (CC_SRKgraph *G, int *oncount, int *oecount,
        int **olist, double **olen, CC_SRKexpinfo *expand, int level)
{
    int rval = 0;
    int k;

    CCcut_SRK_increment_marker (G);

    switch (level) {
    case 1:
        CCcut_SRK_identify_paths_to_edges (G, &k, 0);
        CCcut_SRK_identify_one_triangles (G, (int *) NULL, (CC_SRKnode *) NULL,
                                          0.001, 2.0, 1);
        CCcut_SRK_identify_tight_squares (G, (int *) NULL, 2.0, 1);
        break;
    case 2:
        CCcut_SRK_identify_paths_to_edges (G, &k, 0);
        CCcut_SRK_identify_one_square (G, (int *) NULL, 0.001, 1.0, 1);
        CCcut_SRK_identify_tight_squares (G, (int *) NULL, 2.0, 1);
        break;
    case 3:
        CCcut_SRK_identify_paths_to_edges (G, &k, 0);
        CCcut_SRK_identify_one_triangles (G, (int *) NULL, (CC_SRKnode *) NULL,
                                          0.001, 1.5, 1);
        CCcut_SRK_identify_tight_squares (G, (int *) NULL, 2.0, 1);
        break;
    case 4:
        CCcut_SRK_identify_paths_to_edges (G, &k, 0);
        CCcut_SRK_identify_one_square (G, (int *) NULL, 0.001, 1.0, 1);
        CCcut_SRK_identify_one_triangles (G, (int *) NULL, (CC_SRKnode *) NULL,
                                          0.001, 1.5, 1);
        break;
    }

    rval = CCcut_SRK_grab_edges (G, oncount, oecount, olist, olen, expand);
    if (rval) {
        fprintf (stderr, "grab edges failed\n");
        goto CLEANUP;
    }

CLEANUP:

    return rval;
}
