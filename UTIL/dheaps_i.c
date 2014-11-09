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
/*                           DHEAP ROUTINES                                 */
/*                                                                          */
/*                                                                          */
/*                              TSP CODE                                    */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: February 9, 1995                                                  */
/*  Reference: R.E. Tarjan, Data Structures and Network Algorithms          */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCutil_dheap_init (CCdheap *h, int k)                               */
/*        -h should point to a CCdheap struct.                              */
/*        -k the max number of elements in the dheap.                       */
/*                                                                          */
/*  void CCutil_dheap_free (CCdheap *h)                                     */
/*    -frees the spaces allocated by CCutil_dheap_init                      */
/*                                                                          */
/*  int CCutil_dheap_resize (CCdheap *h, int newsize)                       */
/*    -REALLOCs h so it can contain newsize elements.                       */
/*    -returns -1 if it can't resize the heap.                              */
/*                                                                          */
/*  int CCutil_dheap_findmin (CCdheap *h)                                   */
/*    -returns the index of the element with min value h->key[i]            */
/*    -returns -1 if no elements in heap.                                   */
/*                                                                          */
/*  int CCutil_dheap_insert (CCdheap *h, int i)                             */
/*    -inserts the element with index i (so its key should be loaded        */
/*     beforehand in h->key[i]).                                            */
/*                                                                          */
/*  void CCutil_dheap_delete (CCdheap *h, int i)                            */
/*    -deletes the element with index i.                                    */
/*                                                                          */
/*  int CCutil_dheap_deletemin (CCdheap *h)                                 */
/*    -returns the min element in the heap, and deletes the min element     */
/*    -returns -1 if no elements in heap.                                   */
/*                                                                          */
/*  void CCutil_dheap_changekey (CCdheap *h, int i, double newkey)          */
/*    -changes the key of the element with index i to newkey.               */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/*  NOTES:                                                                  */
/*      A k-element heap will malloc 16k bytes of memory. If memory is      */
/*  tight, using integer keys (instead of doubles), brings it down to       */
/*  12k bytes, and if arbitrary deletions are not required, with a little   */
/*  rewriting, the h->loc field can be eliminated, bring the space down     */
/*  to 8k bytes.                                                            */
/*      These routines work with indices into the h->key array, so in       */
/*  some cases, you will need to maintain a separate names array to know    */
/*  what element belongs to index i. For an example, see the k_nearest      */
/*  code in kdnear.c.                                                       */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"

#define HEAP_D 3
#define HEAP_UP(x) (((x)-1)/HEAP_D)
#define HEAP_DOWN(x) (((x)*HEAP_D)+1)


static void
    dheap_siftup (CCdheap *h, int i, int x),
    dheap_siftdown (CCdheap *h, int i, int x);

static int
    dheap_minchild (int x, CCdheap *h);



int CCutil_dheap_init (CCdheap *h, int k)
{
    h->loc = (int *) NULL;
    h->key = (double *) NULL;
    h->entry = CC_SAFE_MALLOC (k, int);
    if (!h->entry)
        return 1;
    h->loc = CC_SAFE_MALLOC (k, int);
    if (!h->loc) {
        CC_FREE (h->entry, int);
        return 1;
    }
    h->key = CC_SAFE_MALLOC (k, double);
    if (!h->key) {
        CC_FREE (h->entry, int);
        CC_FREE (h->loc, int);
        return 1;
    }
    h->total_space = k;
    h->size = 0;
    return 0;
}

void CCutil_dheap_free (CCdheap *h)
{
    CC_IFFREE (h->entry, int);
    CC_IFFREE (h->loc, int);
    CC_IFFREE (h->key, double);
}

int CCutil_dheap_resize (CCdheap *h, int newsize)
{
    if (newsize < h->size || newsize < h->total_space) return 0;
    if (CCutil_reallocrus_count ((void **) &(h->key), newsize,
                                 sizeof (double))) {
        return -1;
    }
    if (CCutil_reallocrus_count ((void **) &(h->entry), newsize,
                                 sizeof (int))) {
        return -1;
    }
    if (CCutil_reallocrus_count ((void **) &(h->loc), newsize, sizeof (int))) {
        return -1;
    }
    h->total_space = newsize;

    return 0;
}

int CCutil_dheap_findmin (CCdheap *h)
{
    if (h->size == 0)
        return -1;
    else
        return h->entry[0];
}

int CCutil_dheap_insert (CCdheap *h, int i)
{
    if (h->size >= h->total_space) {
        fprintf (stderr, "Error - heap already full\n");
        return 1;
    }
    h->size++;
    dheap_siftup (h, i, h->size - 1);
    return 0;
}

void CCutil_dheap_delete (CCdheap *h, int i)
{
    int j;

    h->size--;
    j = h->entry[h->size];
    h->entry[h->size] = -1;

    if (j != i) {
        if (h->key[j] <= h->key[i]) {
            dheap_siftup (h, j, h->loc[i]);
        } else {
            dheap_siftdown (h, j, h->loc[i]);
        }
    }
}

int  CCutil_dheap_deletemin (CCdheap *h)
{
    int i;

    if (h->size == 0)
        return -1;
    else {
        i = h->entry[0];
        CCutil_dheap_delete (h, i);
        return i;
    }
}

void CCutil_dheap_changekey (CCdheap *h, int i, double newkey)
{
    if (newkey < h->key[i]) {
        h->key[i] = newkey;
        dheap_siftup (h, i, h->loc[i]);
    } else if (newkey > h->key[i]) {
        h->key[i] = newkey;
        dheap_siftdown (h, i, h->loc[i]);
    }
}

static void dheap_siftup (CCdheap *h, int i, int x)
{
    int p;

    p = HEAP_UP (x);
    while (x && h->key[h->entry[p]] > h->key[i]) {
        h->entry[x] = h->entry[p];
        h->loc[h->entry[p]] = x;
        x = p;
        p = HEAP_UP (p);
    }
    h->entry[x] = i;
    h->loc[i] = x;
}

static void dheap_siftdown (CCdheap *h, int i, int x)
{
    int c;

    c = dheap_minchild (x, h);

    while (c >= 0 && h->key[h->entry[c]] < h->key[i]) {
        h->entry[x] = h->entry[c];
        h->loc[h->entry[c]] = x;
        x = c;
        c = dheap_minchild (c, h);
    }
    h->entry[x] = i;
    h->loc[i] = x;
}

static int dheap_minchild (int x, CCdheap *h)
{
    int c = HEAP_DOWN (x);
    int cend;
    double minval;
    int minloc;

    if (c >= h->size)
        return -1;
    minval = h->key[h->entry[c]];
    minloc = c;
    cend = c + HEAP_D;
    if (h->size < cend)
        cend = h->size;
    for (c++; c < cend; c++) {
        if (h->key[h->entry[c]] < minval) {
            minval = h->key[h->entry[c]];
            minloc = c;
        }
    }
    return minloc;
}
