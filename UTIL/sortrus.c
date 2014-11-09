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
/*                         SORTING ROUTINES                                 */
/*                                                                          */
/*                             TSP CODE                                     */
/*                                                                          */
/*   Written by:  Applegate, Bixby, Chvatal, and Cook                       */
/*   DATE:  February 24, 1994                                               */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  char *CCutil_linked_radixsort (char *data, char *datanext,              */
/*      char *dataval, int valsize)                                         */
/*    USAGE:                                                                */
/*      head = (bar *) CCutil_linked_radixsort ((char *) head,              */
/*         (char *) &(head->next), (char *) &(head->val), sizeof (int));    */
/*    Then head is the start of the linked list in increasing order of      */
/*    val, with next as the field that links the bars.                      */
/*    WARNING: DOES NOT HANDLE NEGATIVE NUMBERS PROPERLY.                   */
/*                                                                          */
/*  void CCutil_int_array_quicksort (int *len, int n)                       */
/*    len - the array to be sorted                                          */
/*    n - the number of elements in len                                     */
/*    Uses quicksort to put len in increasing order.                        */
/*                                                                          */
/*  void CCutil_int_perm_quicksort (int *perm, int *len, int n)             */
/*    perm - must be allocated and initialized by the calling routine,      */
/*           it will be arranged in increasing order of len.                */
/*    n - the number of elements in perm and len.                           */
/*                                                                          */
/*  void CCutil_double_perm_quicksort (int *perm, double *len, int n)       */
/*    perm - must be allocated and initialized by the calling routine,      */
/*           it will be arranged in increasing order of len.                */
/*    n - the number of elements in perm and len.                           */
/*                                                                          */
/*  void CCutil_rselect (int *arr, int l, int r, int m,                     */
/*      double *coord, CCrandstate *rstate)                                 */
/*    arr - permutation that will be rearranged                             */
/*    l,r - specify the range of arr that we are interested in              */
/*    m - is the index into l,r that is the break point for the perm        */
/*    coord - gives the keys that determine the ordering                    */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "macrorus.h"

#define BITS_PER_PASS (8)

#define NBINS (1<<BITS_PER_PASS)


static void
    select_split (int *arr, int n, double v, int *start, int *end,
           double *coord),
    select_sort (int *arr, int n, double *coord),
    select_sort_dsample (double *samp, int n);



char *CCutil_linked_radixsort (char *data, char *datanext, char *dataval,
        int valsize)
{
    size_t nextoff = datanext - data;
    size_t valoff = dataval - data;
    int i;
    char *head[NBINS];
    char **tail[NBINS];
    char *p;
    char **last;
    int j;
    int v;


    for (j = valsize - 1; j >= 0; j--) {
        for (i = 0; i < NBINS; i++) {
            head[i] = (char *) NULL;
            tail[i] = &head[i];
        }
        for (p = data; p; p = *(char **) (p + nextoff)) {
            v = (unsigned char) p[valoff + j];
            *tail[v] = p;
            tail[v] = (char **) (p + nextoff);
        }
        last = &data;
        for (i = 0; i < NBINS; i++) {
            if (head[i]) {
                *last = head[i];
                last = tail[i];
            }
        }
        *last = (char *) NULL;
    }
    return data;
}

void CCutil_int_array_quicksort (int *len, int n)
{
    int i, j, temp, t;

    if (n <= 1)
        return;

    CC_SWAP (len[0], len[(n - 1)/2], temp);

    i = 0;
    j = n;
    t = len[0];

    while (1) {
        do i++; while (i < n && len[i] < t);
        do j--; while (len[j] > t);
        if (j < i) break;
        CC_SWAP (len[i], len[j], temp);
    }
    CC_SWAP (len[0], len[j], temp);

    CCutil_int_array_quicksort (len, j);
    CCutil_int_array_quicksort (len + i, n - i);
}

void CCutil_int_perm_quicksort (int *perm, int *len, int n)
{
    int i, j, temp, t;

    if (n <= 1)
        return;

    CC_SWAP (perm[0], perm[(n - 1)/2], temp);

    i = 0;
    j = n;
    t = len[perm[0]];

    while (1) {
        do i++; while (i < n && len[perm[i]] < t);
        do j--; while (len[perm[j]] > t);
        if (j < i) break;
        CC_SWAP (perm[i], perm[j], temp);
    }
    CC_SWAP (perm[0], perm[j], temp);

    CCutil_int_perm_quicksort (perm, len, j);
    CCutil_int_perm_quicksort (perm + i, len, n - i);
}


void CCutil_double_perm_quicksort (int *perm, double *len, int n)
{
    int i, j, temp;
    double t;

    if (n <= 1)
        return;

    CC_SWAP (perm[0], perm[(n - 1)/2], temp);

    i = 0;
    j = n;
    t = len[perm[0]];

    while (1) {
        do i++; while (i < n && len[perm[i]] < t);
        do j--; while (len[perm[j]] > t);
        if (j < i) break;
        CC_SWAP (perm[i], perm[j], temp);
    }
    CC_SWAP (perm[0], perm[j], temp);

    CCutil_double_perm_quicksort (perm, len, j);
    CCutil_double_perm_quicksort (perm + i, len, n - i);
}


/**********  Median - Select Routines **********/

/* NSAMPLES should be odd */
#define NSAMPLES 3
#define SORTSIZE 20


void CCutil_rselect (int *arr, int l, int r, int m, double *coord,
        CCrandstate *rstate)
{
    double samplevals[NSAMPLES];
    int i;
    int st, en;
    int n;

    arr += l;
    n = r - l + 1;
    m -= l;

    while (n > SORTSIZE) {
        for (i = 0; i < NSAMPLES; i++) {
            samplevals[i] = coord[arr[CCutil_lprand (rstate) % n]];
        }
        select_sort_dsample (samplevals, NSAMPLES);
        select_split (arr, n, samplevals[(NSAMPLES - 1) / 2], &st, &en, coord);
        if (st > m) {
            n = st;
        } else if (en <= m) {
            arr += en;
            n -= en;
            m -= en;
        } else {
            return;
        }
    }

    select_sort (arr, n, coord);
    return;
}

static void select_split (int *arr, int n, double v, int *start, int *end,
                          double *coord)
{
    int i, j, k;
    int t;

    i = 0;
    j = k = n;

    while (i < j) {
        if (coord[arr[i]] < v) {
            i++;
        } else if (coord[arr[i]] == v) {
            j--;
            CC_SWAP (arr[i], arr[j], t);
        } else {
            j--;
            k--;
            t = arr[i];
            arr[i] = arr[j];
            arr[j] = arr[k];
            arr[k] = t;
        }
    }
    *start = j;
    *end = k;
    return;
}

static void select_sort (int *arr, int n, double *coord)
{
    int i, j;
    int t;

    for (i = 1; i < n; i++) {
        t = arr[i];
        for (j = i; j > 0 && coord[arr[j - 1]] > coord[t]; j--) {
            arr[j] = arr[j - 1];
        }
        arr[j] = t;
    }
}

static void select_sort_dsample (double *samp, int n)
{
    int i, j;
    double t;

    for (i = 1; i < n; i++) {
        t = samp[i];
        for (j = i; j > 0 && samp[j - 1] > t; j--) {
            samp[j] = samp[j - 1];
        }
        samp[j] = t;
    }
}
