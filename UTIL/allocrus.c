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
/*                   MEMORY ALLOCATION MACROS                               */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: February 2, 1995 (cofeb16)                                        */
/*                                                                          */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  void *CCutil_allocrus (size_t size)                                     */
/*    RETURNS a pointer to an allocated block of "size" memory.             */
/*                                                                          */
/*  void CCutil_freerus (void *ptr)                                         */
/*    FREES ptr.                                                            */
/*                                                                          */
/*  void *CCutil_reallocrus (void *ptr, size_t size)                        */
/*    REALLOCS ptr to size bytes.                                           */
/*                                                                          */
/*  int CCutil_reallocrus_scale (void **pptr, int *pnnum, int count,        */
/*      double scale, size_t size)                                          */
/*    void **pptr (a reference to the pointer to the allocated space)       */
/*    int *pnnum (a reference to the number of objects in the               */
/*                allocated space)                                          */
/*    int count (a minimum value for the new nnum)                          */
/*    double scale (a scale factor to apply to nnum)                        */
/*    int size (the size of objects to be realloced)                        */
/*    RETURNS 0 if *pptr was successfully changed to point to at            */
/*            least max(*pnnum*scale, *pnnum+1000, count) objects.          */
/*            *pnnum is changed to the new object count.                    */
/*            Otherwise, prints an error message, leaves *pptr and          */
/*            *pnnum alone, and returns nonzero.                            */
/*                                                                          */
/*  int CCutil_reallocrus_count (void **pptr, int count,                    */
/*      size_t size)                                                        */
/*    void **pptr (a reference to the pointer to the allocated space)       */
/*    int count (number of objects to be realloced)                         */
/*    int size (the size of the objects to be realloced)                    */
/*    RETURNS 0 is successful, and 1 if the realloc failed.                 */
/*                                                                          */
/*  CCbigchunkptr *CCutil_bigchunkalloc (void)                              */
/*         RETURNS a CCbigchunkptr with the "this_one" field loaded with a  */
/*                 a pointer to a bigchunk of memory.                       */
/*    NOTES:                                                                */
/*       The idea is to use bigchunks (the size of a bigchunk is defined    */
/*       by CC_BIGCHUNK in util.h) to supply local routines with memory     */
/*       for ptrs, so the memory can be shared with other                   */
/*       local routines.                                                    */
/*                                                                          */
/*  CCutil_bigchunkfree (CCbigchunkptr *bp)                                 */
/*    ACTION: Frees a CCbigchunkptr.                                        */
/*                                                                          */
/*  void CCptrworld_init (CCptrworld *world)                                */
/*     initialize a CCptrworld with 1 reference                             */
/*                                                                          */
/*  void CCptrworld_add (CCptrworld *world)                                 */
/*     add a reference to a CCptrworld                                      */
/*                                                                          */
/*  void CCptrworld_delete (CCptrworld *world)                              */
/*     delete a reference to a ptrworld, and free if no more references     */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"

typedef struct CCbigchunk {
    char space[CC_BIGCHUNK];
    CCbigchunkptr ptr;
} CCbigchunk;
    
void *CCutil_allocrus (size_t size)
{
    void *mem = (void *) NULL;

    if (size == 0) {
        fprintf (stderr, "Warning: 0 bytes allocated\n");
    }

    mem = (void *) malloc (size);
    if (mem == (void *) NULL) {
        fprintf (stderr, "Out of memory. Asked for %d bytes\n", (int) size);
    }
    return mem;
}

void CCutil_freerus (void *p)
{
    if (!p) {
        fprintf (stderr, "Warning: null pointer freed\n");
        return;
    }

    free (p);
}

void *CCutil_reallocrus (void *ptr, size_t size)
{
    void *newptr;

    if (!ptr) {
        return CCutil_allocrus (size);
    } else {
        newptr = (void *) realloc (ptr, size);
        if (!newptr) {
            fprintf (stderr, "Out of memory.  Tried to grow to %d bytes\n",
                     (int) size);
        }
        return newptr;
    }
}

int CCutil_reallocrus_scale (void **pptr, int *pnnum, int count, double scale,
        size_t size)
{
    int newsize = (int) (((double) *pnnum) * scale);
    void *p;

    if (newsize < *pnnum+1000) newsize = *pnnum+1000;
    if (newsize < count) newsize = count;
    p = CCutil_reallocrus (*pptr, newsize * size);
    if (!p) {
        return 1;
    } else {
        *pptr = p;
        *pnnum = newsize;
        return 0;
    }
}

int CCutil_reallocrus_count (void **pptr, int count, size_t size)
{
    void *p = CCutil_reallocrus (*pptr, count * size);

    if (!p) {
        return 1;
    } else {
        *pptr = p;
        return 0;
    }
}


CCbigchunkptr *CCutil_bigchunkalloc (void)
{
    CCbigchunk *p = CC_SAFE_MALLOC (1, CCbigchunk);

    if (p == (CCbigchunk *) NULL) {
        fprintf (stderr, "Out of memory in CCutil_bigchunkalloc\n");
        return (CCbigchunkptr *) NULL;
    }
    p->ptr.this_chunk = p;
    p->ptr.this_one = (void *) p->space;
    return &(p->ptr);
}

void CCutil_bigchunkfree (CCbigchunkptr *bp)
{
    /* This copy is necessary since CC_FREE zeros its first argument */
    CCbigchunk *p = bp->this_chunk;
    
    CC_FREE (p, CCbigchunk);
}

void CCptrworld_init (CCptrworld *world)
{
    world->refcount = 1;
    world->freelist = (void *) NULL;
    world->chunklist = (CCbigchunkptr *) NULL;
}

void CCptrworld_add (CCptrworld *world)
{
    world->refcount++;
}

void CCptrworld_delete (CCptrworld *world)
{
    world->refcount--;
    if (world->refcount <= 0) {
        CCbigchunkptr *bp, *bpnext;

        for (bp = world->chunklist ; bp; bp = bpnext) {
            bpnext = bp->next;
            CCutil_bigchunkfree (bp);
        }
        world->chunklist = (CCbigchunkptr *) NULL;
        world->freelist = (void *) NULL;
        world->refcount = 0;
    }
}
