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
/*                    Compute the Blocks of a Graph                         */
/*                                                                          */
/*                              TSP CODE                                    */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: September 29, 1997                                                */
/*                                                                          */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCcombs_find_blocks (int ncount, int ecount, int *elist,            */
/*      double *x, int *nblocks, int **blockcnts, int ***blocks,            */
/*      int *ncutnodes, int **cutnodes)                                     */
/*    RETURNS the 2-connected components of the graph specified by the      */
/*     edgeset                                                              */
/*     -ncount is the number of nodes                                       */
/*     -ecount is the number of edges                                       */
/*     -elist is the edge list in node node format                          */
/*     -x is an vector of length ecount (it can be NULL); if it is not      */
/*      NULL, then the blocks components will be for the graph consisting   */
/*      of the edges e with epsilon < x_e < 1 - epsilon                     */
/*     -nblocks will return the number of blocks (it can be NULL)           */
/*     -blockcnts with return the size of the blocks (it can be NULL)       */
/*     -blocks will return the members of the blocks (it can be NULL)       */
/*     -ncutnodes with return the number of cutnodes (it can be NULL)       */
/*     -cutnodes will return the cutnodes (it can be NULL)                  */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "combs.h"

typedef struct node {
    int *adj;
    int lowpt;
    int number;
    int degree;
    int magiclabel;
    int mark;
} node;

typedef struct graph {
    node *nodelist;
    int  *adjspace;
    int ncount;
    int ecount;
    int magicnum;
    CCptrworld intptr_world;
    CCptrworld intptrptr_world;
} graph;

typedef struct shortedge {
    int ends[2];
} shortedge;

typedef struct intptr {
    int this;
    struct intptr *next;
} intptr;

typedef struct intptrptr {
    intptr *this;
    struct intptrptr *next;
} intptrptr;



static int
    findblocks (graph *G, int v, int u, int *lastnumber,
        shortedge **top_pointer, intptr **clist, intptrptr **blist),
    build_graph (graph *G, int ncount, int ecount, int *elist,
        double *x),
    add_node_to_list (graph *G, int n, intptr **list);

static void
    init_graph (graph *G),
    free_graph (graph *G),
    block_free_world (CCptrworld *intptr_world, CCptrworld *intptrptr_world);


CC_PTRWORLD_LIST_ROUTINES (intptr, int, intptralloc, intptr_bulk_alloc,
        intptrfree, intptr_listadd, intptr_listfree)
CC_PTRWORLD_LEAKS_ROUTINE (intptr, intptr_check_leaks, this, int)

CC_PTRWORLD_LIST_ROUTINES (intptrptr, intptr *, intptrptralloc,
        intptrptr_bulk_alloc, intptrptrfree, intptrptr_listadd,
        intptrptr_listfree)
CC_PTRWORLD_LEAKS_ROUTINE (intptrptr, intptrptr_check_leaks, this, intptr *)


int CCcombs_find_blocks (int ncount, int ecount, int *elist, double *x,
        int *nblocks, int **blockcnts, int ***blocks, int *ncutnodes,
        int **cutnodes)
{
    graph G;
    shortedge *edgestack = (shortedge *) NULL;
    shortedge *top;
    int i, lastnumber;
    int rval = 0;
    intptr *clist = (intptr *) NULL;
    intptr **clistptr;
    intptrptr *blist = (intptrptr *) NULL;
    intptrptr **blistptr;

    if (nblocks) *nblocks = 0;
    if (blockcnts) *blockcnts = (int *) NULL;
    if (blocks)  *blocks = (int **) NULL;
    if (ncutnodes) *ncutnodes = 0;
    if (cutnodes) *cutnodes = (int *) NULL;
    init_graph (&G);

    if (cutnodes && !ncutnodes) {
        fprintf (stderr, "must use ncutnodes if using cutnodes\n");
        rval = 1; goto CLEANUP;
    }
    if (blocks && (!nblocks || !blockcnts)) {
        fprintf (stderr, "must use nblocks and blockcnts if using blocks\n");
        rval = 1; goto CLEANUP;
    }

    rval = build_graph (&G, ncount, ecount, elist, x);
    if (rval) {
        fprintf (stderr, "build_graph failed\n");
        goto CLEANUP;
    }

    if (G.ecount > 0) {
        edgestack = CC_SAFE_MALLOC (G.ecount, shortedge);
        if (!edgestack) {
            fprintf (stderr, "out of memory in CCcombs_find_blocks\n");
            rval = 1; goto CLEANUP;
        }
    } else {
        goto CLEANUP;
    }

    lastnumber = 1;

    if (cutnodes) {
        clistptr = &clist;
    } else {
        clistptr = (intptr **) NULL;
    }

    if (blocks) {
        blistptr = &blist;
    } else {
        blistptr = (intptrptr **) NULL;
    }

    for (i = 0; i < G.ncount; i++) {
        if (!G.nodelist[i].number) {
            top = edgestack;
            rval = findblocks (&G, i, -1, &lastnumber, &top, clistptr,
                               blistptr);
            if (rval) {
                fprintf (stderr, "findblocks failed\n"); goto CLEANUP;
            }
        }
    }

    if (cutnodes) {
        intptr *ip;
        int k = 0;

        if (clist) {
            for (ip = clist->next; ip; ip = ip->next) {
                k++;
            }
            if (k) {
                *cutnodes = CC_SAFE_MALLOC (k, int);
                if (*cutnodes == (int *) NULL) {
                    fprintf (stderr, "out of memory in CCcombs_find_blocks\n");
                    rval = 1; goto CLEANUP;
                }
                k = 0;
                for (ip = clist->next; ip; ip = ip->next) {
                    (*cutnodes)[k++] = ip->this;
                }
                *ncutnodes = k;
            }
        }
    }

    if (blocks) {
        intptr *ip;
        intptrptr *ipp;
        int *b = (int *) NULL;
        int k = 0;
        int bnt = 0;

        for (ipp = blist; ipp; ipp = ipp->next) {
            bnt++;
        }
        if (bnt) {
            *blocks = CC_SAFE_MALLOC (bnt, int *);
            *blockcnts = CC_SAFE_MALLOC (bnt, int);
            if (*blocks == (int **) NULL) {
                fprintf (stderr, "out of memory in CCcombs_find_blocks\n");
                if (cutnodes) {
                    CC_IFFREE (*cutnodes, int);
                    *ncutnodes = 0;
                }
                rval = 1; goto CLEANUP;
            }
            bnt = 0;
            for (ipp = blist; ipp; ipp = ipp->next) {
                k = 0;
                for (ip = ipp->this; ip; ip = ip->next) {
                    k++;
                }
                b = CC_SAFE_MALLOC (k, int);
                if (b == (int *) NULL) {
                    fprintf (stderr, "out of memory in CCcombs_find_blocks\n");
                    for (k = 0; k < bnt; k++) {
                        CC_IFFREE ((*blocks)[k], int);
                    }
                    CC_IFFREE (*blocks, int *);
                    CC_IFFREE (*blockcnts, int);
                    if (cutnodes) {
                        CC_IFFREE (*cutnodes, int);
                        *ncutnodes = 0;
                    }
                    rval = 1; goto CLEANUP;
                }
                k = 0;
                for (ip = ipp->this; ip; ip = ip->next) {
                    b[k++] = ip->this;
                }
                (*blocks)[bnt] = b;
                (*blockcnts)[bnt++] = k;
            }
            *nblocks = bnt;
        }
    }
             

CLEANUP:

    free_graph (&G);
    CC_IFFREE (edgestack, shortedge);
    intptr_listfree (&G.intptr_world, clist); 
    if (blist) {
        intptrptr *ipp;
        for (ipp = blist; ipp; ipp = ipp->next) {
            intptr_listfree (&G.intptr_world, ipp->this);
        }
        intptrptr_listfree (&G.intptrptr_world, blist);
    }
    block_free_world (&G.intptr_world, &G.intptrptr_world);
    return rval;
}

static int findblocks (graph *G, int v, int u, int *lastnumber,
       shortedge **top_pointer, intptr **clist, intptrptr **blist)
{
    shortedge *top = *top_pointer;
    node *nlist = G->nodelist;
    int i, w;
    int rval = 0;
    intptr *b;

    /* Based on R. Tarjan, "Depth-first search and linear graph algorithms" */
    /*    SIAM Journal on Computing 1 (1972) 146-160.                       */
    /*    - Finds edgesets of the blocks.                                   */

    nlist[v].number = nlist[v].lowpt = (*lastnumber)++;
    for (i = 0; i < nlist[v].degree; i++) {
        w = nlist[v].adj[i];
        if (!nlist[w].number) {
            top->ends[0] = v;
            top->ends[1] = w;
            top++;
            rval = findblocks (G, w, v, lastnumber, &top, clist, blist);
            if (rval) goto CLEANUP;
            if (nlist[w].lowpt < nlist[v].lowpt)
                nlist[v].lowpt = nlist[w].lowpt;
            if (nlist[w].lowpt >= nlist[v].number) {   /* v is a cutnode */
                rval = intptr_listadd (clist, v, &G->intptr_world);
                if (rval) goto CLEANUP;
                G->magicnum++;
                b = (intptr *) NULL;
                while (G->nodelist[((top - 1))->ends[0]].number >=
                                               nlist[w].number) {
                    top--;
                    if (blist) {
                        rval = add_node_to_list (G, top->ends[0], &b);
                        if (rval) goto CLEANUP;
                        rval = add_node_to_list (G, top->ends[1], &b);
                        if (rval) goto CLEANUP;
                    }
                }
                top--;
                if (blist) {
                    rval = add_node_to_list (G, top->ends[0], &b);
                    if (rval) goto CLEANUP;
                    rval = add_node_to_list (G, top->ends[1], &b);
                    if (rval) goto CLEANUP;
                    rval = intptrptr_listadd (blist, b, &G->intptrptr_world);
                    if (rval) goto CLEANUP;
                }
            }
        } else if (nlist[w].number < nlist[v].number && w != u) {
            top->ends[0] = v;
            top->ends[1] = w;
            top++;
            if (nlist[w].number < nlist[v].lowpt)
                nlist[v].lowpt = nlist[w].number;
        }
    }
    *top_pointer = top;

CLEANUP:

    return rval;
}

static int build_graph (graph *G, int ncount, int ecount, int *elist,
                 double *x)
{
    int rval = 0;
    int i;
    int *p;
    node *nodelist;

    G->nodelist = (node *) NULL;
    G->adjspace = (int *) NULL;
    G->ncount = ncount;
    if (x) {
        G->ecount = 0;
        for (i = 0; i < ecount; i++) {
            if (x[i] > CCcombs_BLOCK_ZERO_EPSILON &&
                x[i] < 1.0 - CCcombs_BLOCK_ZERO_EPSILON)
                G->ecount++;
        }
    } else {
       G->ecount = ecount;
    }

    if (G->ncount > 0) {
        G->nodelist = CC_SAFE_MALLOC (G->ncount, node);
        if (!G->nodelist) {
            fprintf (stderr, "out of memory in build_graph\n");
            rval = 1; goto CLEANUP;
        }
    }

    if (G->ecount > 0) {
        G->adjspace = CC_SAFE_MALLOC (2*G->ecount, int);
        if (!G->adjspace) {
            fprintf (stderr, "out of memory in build_graph\n");
            rval = 1; goto CLEANUP;
        }
    }


    nodelist = G->nodelist;
    G->magicnum = 0;

    for (i = 0; i < ncount; i++) {
        nodelist[i].degree = 0;
        nodelist[i].mark = 0;
        nodelist[i].magiclabel = 0;
        nodelist[i].lowpt = 0;
        nodelist[i].number = 0;
    }

    if (x) {
        for (i = 0; i < ecount; i++) {
            if (x[i] > CCcombs_BLOCK_ZERO_EPSILON &&
                x[i] < 1.0 - CCcombs_BLOCK_ZERO_EPSILON) {
                nodelist[elist[2*i]].degree++;
                nodelist[elist[2*i+1]].degree++;
            }
        }
    } else {
        for (i = 0; i < ecount; i++) {
            nodelist[elist[2*i]].degree++;
            nodelist[elist[2*i+1]].degree++;
        }
    }

    p = G->adjspace;
    for (i = 0; i < ncount; i++) {
        nodelist[i].adj = p;
        p += nodelist[i].degree;
        nodelist[i].degree = 0;
    }

    if (x) {
        for (i = 0; i < ecount; i++) {
            if (x[i] > CCcombs_BLOCK_ZERO_EPSILON &&
                x[i] < 1.0 - CCcombs_BLOCK_ZERO_EPSILON) {
                nodelist[elist[2*i]].adj[nodelist[elist[2*i]].degree++]
                                                = elist[2*i+1];
                nodelist[elist[2*i+1]].adj[nodelist[elist[2*i+1]].degree++]
                                                = elist[2*i];
            }
        }
    } else {
        for (i = 0; i < ecount; i++) {
            nodelist[elist[2*i]].adj[nodelist[elist[2*i]].degree++]
                                            = elist[2*i+1];
            nodelist[elist[2*i+1]].adj[nodelist[elist[2*i+1]].degree++]
                                            = elist[2*i];
        }
    }

CLEANUP:

    if (rval) {
        CC_IFFREE (G->nodelist, node);
        CC_IFFREE (G->adjspace, int);
    }
    return rval;
}

static void init_graph (graph *G)
{
    G->nodelist = (node *) NULL;
    G->adjspace = (int *) NULL;
    G->ncount   = 0;
    G->ecount   = 0;
    G->magicnum = 0;
    CCptrworld_init (&G->intptr_world);
    CCptrworld_init (&G->intptrptr_world);
}

static void free_graph (graph *G)
{
    if (G) {
        CC_IFFREE (G->nodelist, node);
        CC_IFFREE (G->adjspace, int);
    }
}

static int add_node_to_list (graph *G, int n, intptr **list)
{
    if (list && G->nodelist[n].magiclabel != G->magicnum) {
        G->nodelist[n].magiclabel = G->magicnum;
        return intptr_listadd (list, n, &G->intptr_world);
    } else {
        return 0;
    }
}

static void block_free_world (CCptrworld *intptr_world,
        CCptrworld *intptrptr_world)
{
    int total, onlist;

    if (intptr_check_leaks (intptr_world, &total, &onlist)) {
        fprintf (stderr, "WARNING: %d outstanding intptrs\n",
                 total - onlist);
    }
    CCptrworld_delete (intptr_world);

    if (intptrptr_check_leaks (intptrptr_world, &total, &onlist)) {
        fprintf (stderr, "WARNING: %d outstanding intptrptrs\n",
                 total - onlist);
    }
    CCptrworld_delete (intptrptr_world);
}
