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
/*               ROUTINES TO HASH CLIQUES AND DOMINOS                       */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: October 10, 1995                                                  */
/*  Date: January 28, 2003 (bico)                                           */
/*                                                                          */
/*                                                                          */
/*  EXPORTED FUNCTIONS:                                                     */
/*                                                                          */
/*  int CCtsp_init_cliquehash (CCtsp_lpcuts *cuts, int size)                */
/*    initializes the clique hash storage in cuts.                          */
/*    int size is an estimate of the number of cliques                      */
/*                                                                          */
/*  int CCtsp_register_clique (CCtsp_lpcuts *cuts, CCtsp_lpclique *c)       */
/*    returns an integer index for c, adding c to cuts if necessary         */
/*    -1 ==> failure                                                        */
/*                                                                          */
/*  void CCtsp_free_cliquehash (CCtsp_lpcuts *cuts)                         */
/*    frees the clique hashtable space                                      */
/*                                                                          */
/*  void CCtsp_unregister_clique (CCtsp_lpcuts *cuts, int c)                */
/*    deletes a reference to clique c, and deletes the clique if no         */
/*    references remain                                                     */
/*                                                                          */
/*  void CCtsp_clique_eq (CCtsp_lpclique *c, CCtsp_lpclique *d,             */
/*      int *yes_no)                                                        */
/*    MISSING                                                               */
/*                                                                          */
/*  unsigned int CCtsp_hashclique (CCtsp_lpclique *c)                       */
/*    RETURNS a hash value for the clique.                                  */
/*                                                                          */
/*  int CCtsp_init_dominohash (CCtsp_lpcuts *cuts, int size)                */
/*  void CCtsp_free_dominohash (CCtsp_lpcuts *cuts)                         */
/*  unsigned int CCtsp_hashdomino (CCtsp_lpdomino *d)                       */
/*  void CCtsp_domino_eq (CCtsp_lpdomino *c, CCtsp_lpdomino *d,             */
/*      int *yes_no)                                                        */
/*  int CCtsp_register_domino (CCtsp_lpcuts *cuts, CCtsp_lpdomino *c)       */
/*  void CCtsp_unregister_domino (CCtsp_lpcuts *cuts, int c)                */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "tsp.h"



int CCtsp_init_cliquehash (CCtsp_lpcuts *cuts, int size)
{
    int i;

    cuts->cliquehashsize = (int) CCutil_nextprime ((unsigned int) size);
    cuts->cliquehash = CC_SAFE_MALLOC (cuts->cliquehashsize, int);
    if (!cuts->cliquehash) {
        cuts->cliquehashsize = 0;
        return 1;
    }
    for (i=0; i<cuts->cliquehashsize; i++) {
        cuts->cliquehash[i] = -1;
    }
    cuts->cliquefree = -1;
    return 0;
}

void CCtsp_free_cliquehash (CCtsp_lpcuts *cuts)
{
    CC_IFFREE (cuts->cliquehash, int);
    cuts->cliquehashsize = 0;
}

unsigned int CCtsp_hashclique (CCtsp_lpclique *c)
{
    unsigned int x = 0;
    int i;

    for (i=0; i<c->segcount; i++) {
        x = x * 65537 + c->nodes[i].lo * 4099 + c->nodes[i].hi;
    }
    return x;
}

void CCtsp_clique_eq (CCtsp_lpclique *c, CCtsp_lpclique *d, int *yes_no)
{
    int i;

    *yes_no = 0;
    if (c->segcount != d->segcount) return;
    for (i=0; i<c->segcount; i++) {
        if (c->nodes[i].lo != d->nodes[i].lo) return;
        if (c->nodes[i].hi != d->nodes[i].hi) return;
    }
    *yes_no = 1;
}

int CCtsp_register_clique (CCtsp_lpcuts *cuts, CCtsp_lpclique *c)
{
    int x = CCtsp_hashclique (c) % cuts->cliquehashsize;
    int y = cuts->cliquehash[x];
    CCtsp_segment *new = (CCtsp_segment *) NULL;
    int i;
    int test;

    while (y != -1) {
        CCtsp_clique_eq (c, &cuts->cliques[y], &test);
        if (test) {
            cuts->cliques[y].refcount++;
            return y;
        }
        y = cuts->cliques[y].hashnext;
    }

    new = CC_SAFE_MALLOC (c->segcount, CCtsp_segment);
    if (!new) {
        fprintf (stderr, "out of memory in CCtsp_register_clique\n");
        return -1;
    }

    if (cuts->cliquefree != -1) {
        y = cuts->cliquefree;
        cuts->cliquefree = cuts->cliques[y].hashnext;
    } else {
        if (cuts->cliqueend >= cuts->cliquespace) {
            if (CCutil_reallocrus_scale ((void **) &cuts->cliques,
                    &cuts->cliquespace, cuts->cliqueend + 1, 1.3,
                    sizeof (CCtsp_lpclique))) {
                CC_FREE (new, CCtsp_segment);
                return -1;
            }
        }
        y = cuts->cliqueend++;
    }
    cuts->cliques[y].segcount = c->segcount;
    for (i=0; i<c->segcount; i++) {
        new[i] = c->nodes[i];
    }
    cuts->cliques[y].nodes = new;
    cuts->cliques[y].refcount = 1;
    cuts->cliques[y].hashnext = cuts->cliquehash[x];
    cuts->cliquehash[x] = y;

    return y;
}

void CCtsp_unregister_clique (CCtsp_lpcuts *cuts, int c)
{
    int x, y, yprev;

    cuts->cliques[c].refcount--;
    if (cuts->cliques[c].refcount) return;
    x = CCtsp_hashclique (&cuts->cliques[c]) % cuts->cliquehashsize;
    y = cuts->cliquehash[x];
    if (y == c) {
        cuts->cliquehash[x] = cuts->cliques[c].hashnext;
    } else {
        yprev = y;
        y = cuts->cliques[y].hashnext;
        while (y != c && y != -1) {
            yprev = y;
            y = cuts->cliques[y].hashnext;
        }
        if (y == -1) {
            fprintf (stderr, "Couldn't find clique to delete from hash\n");
            return;
        }
        cuts->cliques[yprev].hashnext = cuts->cliques[c].hashnext;
    }
    CC_FREE (cuts->cliques[c].nodes, CCtsp_segment);
    cuts->cliques[c].segcount = -1;
    cuts->cliques[c].hashnext = cuts->cliquefree;
    cuts->cliquefree = c;
}

/**********  New material for dominos **********/

int CCtsp_init_dominohash (CCtsp_lpcuts *cuts, int size)
{
    int i;

    cuts->dominohashsize = (int) CCutil_nextprime ((unsigned int) size);
    cuts->dominohash = CC_SAFE_MALLOC (cuts->dominohashsize, int);
    if (!cuts->dominohash) {
        cuts->dominohashsize = 0;
        return 1;
    }
    for (i=0; i<cuts->dominohashsize; i++) {
        cuts->dominohash[i] = -1;
    }
    cuts->dominofree = -1;
    return 0;
}

void CCtsp_free_dominohash (CCtsp_lpcuts *cuts)
{
    CC_IFFREE (cuts->dominohash, int);
    cuts->dominohashsize = 0;
}

unsigned int CCtsp_hashdomino (CCtsp_lpdomino *d)
{
    unsigned int x = 0;
    int i, k;
    CCtsp_lpclique *c;

    for (k = 0; k < 2; k++) {
        c = &(d->sets[k]);
        for (i=0; i<c->segcount; i++) {
            x = x * 65537 + c->nodes[i].lo * 4099 + c->nodes[i].hi;
        }
    }
    return x;
}

void CCtsp_domino_eq (CCtsp_lpdomino *c, CCtsp_lpdomino *d, int *yes_no)
{
    int k;

    for (k = 0; k < 2; k++) {
        CCtsp_clique_eq (&(c->sets[k]), &(d->sets[k]), yes_no);
        if (*yes_no == 0) return;
    }
}

int CCtsp_register_domino (CCtsp_lpcuts *cuts, CCtsp_lpdomino *c)
{
    int x = CCtsp_hashdomino (c) % cuts->dominohashsize;
    int y = cuts->dominohash[x];
    CCtsp_segment *new[2];
    int i, k;
    int test;

    for (k = 0; k < 2; k++) {
        new[k] = (CCtsp_segment *) NULL;
    }

    while (y != -1) {
        CCtsp_domino_eq (c, &cuts->dominos[y], &test);
        if (test) {
            cuts->dominos[y].refcount++;
            return y;
        }
        y = cuts->dominos[y].hashnext;
    }

    for (k = 0; k < 2; k++) {
        new[k] = CC_SAFE_MALLOC (c->sets[k].segcount, CCtsp_segment);
        if (!new[k]) {
            fprintf (stderr, "out of memory in CCtsp_register_domino\n");
            if (k == 1) { CC_FREE (new[0], CCtsp_segment); }
            return -1;
        }
    }

    if (cuts->dominofree != -1) {
        y = cuts->dominofree;
        cuts->dominofree = cuts->dominos[y].hashnext;
    } else {
        if (cuts->dominoend >= cuts->dominospace) {
            if (CCutil_reallocrus_scale ((void **) &cuts->dominos,
                    &cuts->dominospace, cuts->dominoend + 1, 1.3,
                    sizeof (CCtsp_lpdomino))) {
                CC_FREE (new[0], CCtsp_segment);
                CC_FREE (new[1], CCtsp_segment);
                return -1;
            }
        }
        y = cuts->dominoend++;
    }

    for (k = 0; k < 2; k++) {
        cuts->dominos[y].sets[k].segcount = c->sets[k].segcount;
        for (i = 0; i < c->sets[k].segcount; i++) {
            new[k][i] = c->sets[k].nodes[i];
        }
        cuts->dominos[y].sets[k].nodes = new[k];
        cuts->dominos[y].sets[k].hashnext = 0;   /* Not used */
        cuts->dominos[y].sets[k].refcount = 0;   /* Not used */
    }
    cuts->dominos[y].refcount = 1;
    cuts->dominos[y].hashnext = cuts->dominohash[x];
    cuts->dominohash[x] = y;

    return y;
}

void CCtsp_unregister_domino (CCtsp_lpcuts *cuts, int c)
{
    int k, x, y, yprev;

    cuts->dominos[c].refcount--;
    if (cuts->dominos[c].refcount) return;
    x = CCtsp_hashdomino (&cuts->dominos[c]) % cuts->dominohashsize;
    y = cuts->dominohash[x];
    if (y == c) {
        cuts->dominohash[x] = cuts->dominos[c].hashnext;
    } else {
        yprev = y;
        y = cuts->dominos[y].hashnext;
        while (y != c && y != -1) {
            yprev = y;
            y = cuts->dominos[y].hashnext;
        }
        if (y == -1) {
            fprintf (stderr, "Couldn't find domino to delete from hash\n");
            return;
        }
        cuts->dominos[yprev].hashnext = cuts->dominos[c].hashnext;
    }
    for (k = 0; k < 2; k++) {
        CC_FREE (cuts->dominos[c].sets[k].nodes, CCtsp_segment);
        cuts->dominos[c].sets[k].segcount = -1;
    }
    cuts->dominos[c].hashnext = cuts->dominofree;
    cuts->dominofree = c;
}
