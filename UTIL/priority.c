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
/*                      PRIORITY QUEUE ROUTINES                             */
/*                                                                          */
/*                                                                          */
/*                              TSP CODE                                    */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: March 3, 1997                                                     */
/*  Reference: R.E. Tarjan, Data Structures and Network Algorithms          */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCutil_priority_init (CCpriority *pri, int k)                       */
/*    -h should point to a CCpriority struct.                               */
/*    -k an initial allocation for the priority queue.                      */
/*                                                                          */
/*  void CCutil_priority_free (CCpriority *pri)                             */
/*    -frees the spaces allocated for the CCpriority queue.                 */
/*                                                                          */
/*  void *CCutil_priority_findmin (CCpriority *pri, double *keyval)         */
/*    -returns the entry with least key value.                              */
/*    -returns NULL if no entries in heap.                                  */
/*    -if (keyval != NULL), *keyval will be the minimum key value.          */
/*                                                                          */
/*  int CCutil_priority_insert (CCpriority *pri, void *data,                */
/*      double keyval)                                                      */
/*    -adds (data, keyval) to h.                                            */
/*    -returns a handle (>= 0) to use when deleting or changing the entry   */
/*    -returns -1 if out of memory.                                         */
/*                                                                          */
/*  void CCutil_priority_delete (CCpriority *pri, int handle)               */
/*    -deletes an entry from the queue.  handle is the value returned by    */
/*     CCutil_priority_insert.                                              */
/*                                                                          */
/*  void *CCutil_priority_deletemin (CCpriority *pri, double *keyval)       */
/*    -like CCutil_priority_findmin, but also deletes the entry.            */
/*                                                                          */
/*  void CCutil_priority_changekey (CCpriority *pri, int handle,            */
/*      double newkey)                                                      */
/*    -changes the key of an entry in the queue.  handle is the value       */
/*     returned by CCutil_priority_insert.                                  */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/*  NOTES:                                                                  */
/*      These priority queue routines use the CCdheap routines to maintain  */
/*  the priority queue.                                                     */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"

int CCutil_priority_init (CCpriority *pri, int k)
{
    int i;
    int list;

    pri->space = k;
    pri->pri_info = CC_SAFE_MALLOC (k, union CCpri_data);
    if (!pri->pri_info) {
        return -1;
    }
    if (CCutil_dheap_init (&pri->heap, k)) {
        CC_FREE (pri->pri_info, union CCpri_data);
        return -1;
    }

    list = -1;
    for (i=k-1; i>=0; i--) {
        pri->pri_info[i].next = list;
        list = i;
    }
    pri->freelist = list;
    return 0;
}

void CCutil_priority_free (CCpriority *pri)
{
    CCutil_dheap_free (&pri->heap);
    CC_FREE (pri->pri_info, union CCpri_data);
    pri->space = 0;
}

void *CCutil_priority_findmin (CCpriority *pri, double *keyval)
{
    int handle = CCutil_dheap_findmin (&pri->heap);

    if (handle < 0) {
        return (void *) NULL;
    } else {
        if (keyval) *keyval = pri->heap.key[handle];
        return pri->pri_info[handle].data;
    }
}

int CCutil_priority_insert (CCpriority *pri, void *data, double keyval)
{
    int newsize;
    int i;
    int list;

    if (pri->freelist == -1) {
        newsize = (int) (1.3 * pri->space);
        if (newsize < pri->space + 1000) newsize = pri->space + 1000;
        if (CCutil_dheap_resize (&pri->heap, newsize)) {
            return -1;
        }
        if (CCutil_reallocrus_count ((void **) &pri->pri_info, newsize,
                              sizeof (union CCpri_data))) {
            return -1;
        }
        list = -1;
        for (i=newsize-1; i>=pri->space; i--) {
            pri->pri_info[i].next = list;
            list = i;
        }
        pri->space = newsize;
        pri->freelist = list;
    }

    i = pri->freelist;
    pri->freelist = pri->pri_info[i].next;
    pri->pri_info[i].data = data;
    pri->heap.key[i] = keyval;
    if (CCutil_dheap_insert (&pri->heap, i)) {
        fprintf (stderr, "Unexpected failure from CCutil_dheap_insert\n");
        return -1;
    }
    return i;
}

void CCutil_priority_delete (CCpriority *pri, int handle)
{
    CCutil_dheap_delete (&pri->heap, handle);
    pri->pri_info[handle].next = pri->freelist;
    pri->freelist = handle;
}

void *CCutil_priority_deletemin (CCpriority *pri, double *keyval)
{
    int handle = CCutil_dheap_deletemin (&pri->heap);
    void *data;

    if (handle < 0) {
        return (void *) NULL;
    } else {
        if (keyval) *keyval = pri->heap.key[handle];
        data = pri->pri_info[handle].data;
        pri->pri_info[handle].next = pri->freelist;
        pri->freelist = handle;
        return data;
    }
}

void CCutil_priority_changekey (CCpriority *pri, int handle, double newkey)
{
    CCutil_dheap_changekey (&pri->heap, handle, newkey);
}
