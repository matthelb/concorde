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
/*            ROUTINES TO MAP NODE PAIRS TO EDGES                           */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: September 27, 1995                                                */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCutil_edgehash_init (CCutil_edgehash *h, int size)                 */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCutil_edgehash_add (CCutil_edgehash *h, int end1, int end2,        */
/*      int val)                                                            */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCutil_edgehash_set (CCutil_edgehash *h, int end1, int end2,        */
/*      int val)                                                            */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCutil_edgehash_del (CCutil_edgehash *h, int end1, int end2)        */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCutil_edgehash_find (CCutil_edgehash *h, int end1, int end2,       */
/*      int *val)                                                           */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCutil_edgehash_getall (CCutil_edgehash *h, int *ecount,            */
/*      int **elist, int **elen);                                           */
/*    MISSING                                                               */
/*                                                                          */
/*  void CCutil_edgehash_delall (CCutil_edgehash *h)                        */
/*    MISSING                                                               */
/*                                                                          */
/*  void CCutil_edgehash_free (CCutil_edgehash *h)                          */
/*    MISSING                                                               */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "macrorus.h"

CC_PTRWORLD_ROUTINES (CCutil_edgeinf, edgeinf_alloc, edgeinf_bulkalloc,
        edgeinf_free)
CC_PTRWORLD_LISTFREE_ROUTINE (CCutil_edgeinf, edgeinf_listfree, edgeinf_free)

int CCutil_edgehash_init (CCutil_edgehash *h, int size)
{
    unsigned int i;

    h->size = CCutil_nextprime ((unsigned int) size);
    h->mult = (int) sqrt ((double) h->size);
    h->table = CC_SAFE_MALLOC ((int) h->size, CCutil_edgeinf *);
    CCptrworld_init (&h->edgeinf_world);
    
    if (!h->table) {
        h->size = 0;
        return 1;
    }
    for (i=0; i<h->size; i++) {
        h->table[i] = (CCutil_edgeinf *) NULL;
    }
    return 0;
}

int CCutil_edgehash_add (CCutil_edgehash *h, int end1, int end2, int val)
{
    int t;
    unsigned int loc;
    CCutil_edgeinf *e;

    if (h->size == 0) return 1;
    e = edgeinf_alloc(&h->edgeinf_world);
    if (!e) return 1;

    if (end1 > end2) CC_SWAP (end1, end2, t);
    loc = (end1 * h->mult + end2) % h->size;
    e->ends[0] = end1;
    e->ends[1] = end2;
    e->val = val;
    e->next = h->table[loc];
    h->table[loc] = e;
    return 0;
}

int CCutil_edgehash_set (CCutil_edgehash *h, int end1, int end2, int val)
{
    int t;
    unsigned int loc;
    CCutil_edgeinf *e;

    if (h->size == 0) return 1;
    if (end1 > end2) CC_SWAP (end1, end2, t);
    loc = (end1 * h->mult + end2) % h->size;
    for (e = h->table[loc]; e; e = e->next) {
        if (e->ends[0] == end1 && e->ends[1] == end2) {
            e->val = val;
            return 0;
        }
    }
    e = edgeinf_alloc(&h->edgeinf_world);
    if (!e) return 1;

    e->ends[0] = end1;
    e->ends[1] = end2;
    e->val = val;
    e->next = h->table[loc];
    h->table[loc] = e;
    return 0;
}

int CCutil_edgehash_del (CCutil_edgehash *h, int end1, int end2)
{
    int t;
    CCutil_edgeinf **prev;
    CCutil_edgeinf *p;

    if (end1 > end2) CC_SWAP (end1, end2, t);
    if (h->size == 0) return 1;

    prev  = &h->table[(end1 * h->mult + end2) % h->size];
    p = *prev;
    while (p) {
        if (p->ends[0] == end1 && p->ends[1] == end2) {
            *prev = p->next;
            edgeinf_free (&h->edgeinf_world, p);
            return 0;
        }
        prev = &p->next;
        p = *prev;
    }
    return 1;
}

int CCutil_edgehash_getall (CCutil_edgehash *h, int *ecount, int **elist,
        int **elen)
{
    unsigned int i;
    int k = 0;
    CCutil_edgeinf *p;

    for (i = 0; i < h->size; i++) {
        for (p = h->table[i]; p; p = p->next) k++;
    } 
   
    if (k > 0) {
        *elist = CC_SAFE_MALLOC (2*k, int);
        *elen  = CC_SAFE_MALLOC (k, int);
        if (!(*elist) || !(*elen)) {
            fprintf (stderr, "out of memory in CCutil_edgehash_getall\n");
            CC_IFFREE (*elist, int);
            CC_IFFREE (*elen, int);
            return 1;
        }
        *ecount = k;
        k = 0;
        for (i = 0; i < h->size; i++) {
            for (p = h->table[i]; p; p = p->next) {
                (*elist)[2*k]   = p->ends[0];
                (*elist)[2*k+1] = p->ends[1];
                (*elen)[k++]    = p->val;
            }
        } 
    } else {
        *elist = (int *) NULL;
        *elen = (int *) NULL;
        *ecount = 0;
    }
    
    return 0;
}

void CCutil_edgehash_delall (CCutil_edgehash *h)
{
    unsigned int i;

    for (i=0; i<h->size; i++) {
        if (h->table[i]) {
            edgeinf_listfree (&h->edgeinf_world, h->table[i]);
            h->table[i] = (CCutil_edgeinf *) NULL;
        }
    }
}

int CCutil_edgehash_find (CCutil_edgehash *h, int end1, int end2, int *val)
{
    int t;
    CCutil_edgeinf *p;

    *val = 0;
    if (h->size == 0) return -1;
    if (end1 > end2) CC_SWAP (end1, end2, t);

    for (p = h->table[(end1 * h->mult + end2) % h->size]; p; p = p->next) {
        if (p->ends[0] == end1 && p->ends[1] == end2) {
            *val = p->val;
            return 0;
        }
    }
    return -1;
}

void CCutil_edgehash_free (CCutil_edgehash *h)
{
    CCutil_edgehash_delall (h);
    CC_FREE (h->table, CCutil_edgeinf *);
    CCptrworld_delete (&h->edgeinf_world);
    h->size = 0;
}
