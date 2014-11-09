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
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  CCpq_set_INIT(CCpq_set s)                                               */
/*    initializes s to the empty set                                        */
/*                                                                          */
/*  CCpq_set_ISEMPTY(CCpq_set s)                                            */
/*    is an expression which tests if s is empty                            */
/*                                                                          */
/*  CCpq_set_SIZE(CCpq_set s)                                               */
/*    is the size of s                                                      */
/*                                                                          */
/*  CCpq_set_PTR_TO(CCpq_elem e, CCpq_node *q)                              */
/*    is the field of e which points to q                                   */
/*                                                                          */
/*  CCpq_set_PTR_REPLACE(CCpq_elem e, CCpq_node *q, CCpq_node *r)           */
/*    replaces q with r in e                                                */
/*                                                                          */
/*  CCpq_set_PTR_AWAY(CCpq_elem e, CCpq_node *q)                            */
/*    is the field of e which doesn't point to q                            */
/*                                                                          */
/*  CCpq_set_ADD_WORK(CCpq_node *x, CCpq_set s, ELEM_FIELD efield,          */
/*      DIRECTION dir)                                                      */
/*    adds x to s at the dir end, using the efield field to link things     */
/*    together. It is intended to be an internal macro used by              */
/*    CCpq_set_ADD_LEFT and CCpq_set_ADD_RIGHT.                             */
/*                                                                          */
/*  CCpq_set_ADD_LEFT(CCpq_node *x, CCpq_set s, ELEM_FIELD efield)          */
/*    adds x to s at the left end.                                          */
/*                                                                          */
/*  CCpq_set_ADD_RIGHT(CCpq_node *x, CCpq_set s, ELEM_FIELD efield)         */
/*    adds x to s at the right end.                                         */
/*                                                                          */
/*  CCpq_set_ADD(CCpq_node *x, CCpq_set s, ELEM_FIELD efield)               */
/*    adds x to s                                                           */
/*                                                                          */
/*  CCpq_set_DELETE(CCpq_node *x, CCpq_set s, ELEM_FIELD efield)            */
/*    deletes x from s                                                      */
/*                                                                          */
/*  CCpq_set_DELETE2(CCpq_node *x, SET_FIELD sfield,                        */
/*      ELEM_FIELD efield)                                                  */
/*    deletes x from the sfield set of its parent.  If x is an endmost      */
/*    child in this set, then DELETE2 uses the parent pointer to find       */
/*    the set field.  Otherwise, the parent pointer is not used.            */
/*                                                                          */
/*  CCpq_set_LEFT_ELEM(CCpq_set s)                                          */
/*    is the left element of s, or NULL if s is empty                       */
/*                                                                          */
/*  CCpq_set_RIGHT_ELEM(CCpq_set s)                                         */
/*    is the right element of s, or NULL if s is empty                      */
/*                                                                          */
/*  CCpq_set_FOREACH(CCpq_set s, CCpq_node *x, ELEM_FIELD efield,           */
/*      CCpq_node *xprev, CCpq_node *xnext)                                 */
/*    iterates x over elements of s using temporary variables xprev and     */
/*    xnext.                                                                */
/*                                                                          */
/*  CCpq_set_FOREACH_FROM(CCpq_node *x, ELEM_FIELD efield,                  */
/*      CCpq_node *xprev, CCpq_node *xnext)                                 */
/*    iterates x over elements of s starting at x and going away from       */
/*    xprev.  xnext is a temporary variable used in the loop, which also    */
/*    changes xprev and x                                                   */
/*                                                                          */
/*  CCpq_set_FOREACH_DEL(CCpq_set s, CCpq_node *x, ELEM_FIELD efield,       */
/*      CCpq_node *xprev, CCpq_node *xnext)                                 */
/*    iterates x over elements of s using temporary variables xprev and     */
/*    xnext.  x may be deleted from the set in the body of the loop (but    */
/*    delete any other element from s at your own risk).                    */
/*                                                                          */
/*  CCpq_set_FOREACH_ADJ(CCpq_node *x, ELEM_FIELD efield, CCpq_node *z,     */
/*      int itemp)                                                          */
/*    iterates z over the immediate neighbors of x using itemp as a         */
/*    temporary variable. It is just used to save code replication when     */
/*    you want to do something for each of the two neighbors.               */
/*                                                                          */
/****************************************************************************/

#ifndef __PQSETS_H
#define __PQSETS_H

typedef struct CCpq_set {
    int size;
    struct CCpq_node *left;
    struct CCpq_node *right;
} CCpq_set;

typedef struct CCpq_elem {
    struct CCpq_node *ptr1;
    struct CCpq_node *ptr2;
} CCpq_elem;



/* CCpq_set_INIT(CCpq_set s) initializes s to the empty set */

#define CCpq_set_INIT(s) {                                                \
    (s).size = 0;                                                         \
    (s).left = (CCpq_node *) NULL;                                        \
    (s).right = (CCpq_node *) NULL;                                       \
}

/* In the comments below, something of type ELEM_FIELD is the name of a field
   in a CCpq_node of type CCpq_elem, and something of type DIRECTION is either
   left or right. */

/* CCpq_set_ISEMPTY(CCpq_set s) is an expression which tests if s is empty  */

#define CCpq_set_ISEMPTY(s)       ((s).left == (CCpq_node *) NULL)

/* CCpq_set_SIZE(CCpq_set s) is the size of s */

#define CCpq_set_SIZE(s)  ((s).size)

/* CCpq_set_PTR_TO(CCpq_elem e, CCpq_node *q) is the field of e which points
   to q */

#define CCpq_set_PTR_TO(e,q)      (((e).ptr1 == (q)) ? ((e).ptr1) : ((e).ptr2))

/* CCpq_set_PTR_REPLACE(CCpq_elem e, CCpq_node *q, CCpq_node *r) replaces q
   with r in e */

#define CCpq_set_PTR_REPLACE(e,q,r) {                                     \
        if ((e).ptr1 == (q)) {                                            \
                (e).ptr1 = (r);                                           \
        } else {                                                          \
                (e).ptr2 = (r);                                           \
        }                                                                 \
}

/* CCpq_set_PTR_AWAY(pq_elem e, CCpq_node *q) is the field of e which doesn't
   point to q */

#define CCpq_set_PTR_AWAY(e,q)    (((e).ptr1 == (q)) ? ((e).ptr2) : ((e).ptr1))

/* CCpq_set_ADD_WORK(CCpq_node *x, CCpq_set s, ELEM_FIELD efield, DIRECTION
   dir) adds x to s at the dir end, using the efield field to link things
   together. It is intended to be an internal macro used by CCpq_set_ADD_LEFT
   and CCpq_set_ADD_RIGHT. */

#define CCpq_set_ADD_WORK(x,s,efield,dir) {                               \
    (x)->efield.ptr1 = (s).dir;                                           \
    (x)->efield.ptr2 = (CCpq_node *) NULL;                                \
    if ((s).dir) {                                                        \
        CCpq_set_PTR_REPLACE((s).dir->efield,(CCpq_node *) NULL, (x));    \
        (s).dir = (x);                                                    \
    } else {                                                              \
        (s).left = (s).right = (x);                                       \
    }                                                                     \
    (s).size++;                                                           \
}

/* CCpq_set_ADD_LEFT(CCpq_node *x, CCpq_set s, ELEM_FIELD efield) adds x to
   s at the left end. */

#define CCpq_set_ADD_LEFT(x,s,efield) CCpq_set_ADD_WORK(x,s,efield,left)

/* CCpq_set_ADD_RIGHT(CCpq_node *x, CCpq_set s, ELEM_FIELD efield) adds x to
   s at the right end. */

#define CCpq_set_ADD_RIGHT(x,s,efield) CCpq_set_ADD_WORK(x,s,efield,right)

/* CCpq_set_ADD(CCpq_node *x, CCpq_set s, ELEM_FIELD efield) adds x to s */

#define CCpq_set_ADD(x,s,efield) CCpq_set_ADD_LEFT(x,s,efield)

/* CCpq_set_DELETE(CCpq_node *x, CCpq_set s, ELEM_FIELD efield) deletes x
   from s */

#define CCpq_set_DELETE(x,s,efield) {                                     \
    if (CCpq_set_ISEMPTY(s)) {                                            \
        fprintf (stderr, "Error - attempt to delete from empty set\n");   \
    }                                                                     \
    if ((x)->efield.ptr1) {                                               \
        CCpq_set_PTR_REPLACE((x)->efield.ptr1->efield,(x),                \
                             (x)->efield.ptr2);                           \
    } else {                                                              \
        if ((s).left == (x)) {                                            \
                (s).left = (x)->efield.ptr2;                              \
        } else {                                                          \
                (s).right = (x)->efield.ptr2;                             \
        }                                                                 \
    }                                                                     \
    if ((x)->efield.ptr2) {                                               \
        CCpq_set_PTR_REPLACE((x)->efield.ptr2->efield,(x),                \
                             (x)->efield.ptr1);                           \
    } else {                                                              \
        if ((s).right == (x)) {                                           \
                (s).right = (x)->efield.ptr1;                             \
        } else {                                                          \
                (s).left = (x)->efield.ptr1;                              \
        }                                                                 \
    }                                                                     \
    (s).size--;                                                           \
}

/* CCpq_set_DELETE2(CCpq_node *x, SET_FIELD sfield, ELEM_FIELD efield) deletes
   x from the sfield set of its parent.  If x is an endmost child in this
   set, then DELETE2 uses the parent pointer to find the set field.
   Otherwise, the parent pointer is not used.  If the parent pointer is
   NULL, then the parent update is ignored. */

#define CCpq_set_DELETE2(x,sfield,efield) {                               \
    if ((x)->efield.ptr1) {                                               \
        CCpq_set_PTR_REPLACE((x)->efield.ptr1->efield,(x),                \
                             (x)->efield.ptr2);                           \
    } else if ((x)->parent) {                                             \
        if ((x)->parent->sfield.left == (x)) {                            \
                (x)->parent->sfield.left = (x)->efield.ptr2;              \
        } else {                                                          \
                (x)->parent->sfield.right = (x)->efield.ptr2;             \
        }                                                                 \
    }                                                                     \
    if ((x)->efield.ptr2) {                                               \
        CCpq_set_PTR_REPLACE((x)->efield.ptr2->efield,(x),                \
                             (x)->efield.ptr1);                           \
    } else if ((x)->parent) {                                             \
        if ((x)->parent->sfield.right == (x)) {                           \
                (x)->parent->sfield.right = (x)->efield.ptr1;             \
        } else {                                                          \
                (x)->parent->sfield.left = (x)->efield.ptr1;              \
        }                                                                 \
    }                                                                     \
    if ((x)->parent) {                                                    \
        (x)->parent->sfield.size--;                                       \
    }                                                                     \
}

/* CCpq_set_LEFT_ELEM(CCpq_set s) is the left element of s, or NULL if s is
   empty */

#define CCpq_set_LEFT_ELEM(s) ((s).left)

/* CCpq_set_RIGHT_ELEM(CCpq_set s) is the right element of s, or NULL if s is
   empty */

#define CCpq_set_RIGHT_ELEM(s) ((s).right)

/* CCpq_set_FOREACH(CCpq_set s, CCpq_node *x, ELEM_FIELD efield,
   CCpq_node *xprev, CCpq_node *xnext) iterates x over elements of s using
   temporary variables xprev and xnext. */

#define CCpq_set_FOREACH(s,x,efield,xprev,xnext)                          \
        for ((xprev) = (CCpq_node *) NULL,                                \
             (x) = (s).left;                                              \
            (x);                                                          \
             (xnext) = CCpq_set_PTR_AWAY((x)->efield,xprev),              \
             (xprev) = (x),                                               \
             (x) = (xnext)                                                \
)

/* CCpq_set_FOREACH_FROM(CCpq_node *x, ELEM_FIELD efield, CCpq_node *xprev,
   CCpq_node *xnext) iterates x over elements of s starting at x and going
   away from xprev.  xnext is a temporary variable used in the loop, which
   also changes xprev and x */

#define CCpq_set_FOREACH_FROM(x,efield,xprev,xnext)                       \
        for (;                                                            \
            (x);                                                          \
             (xnext) = CCpq_set_PTR_AWAY((x)->efield,xprev),              \
             (xprev) = (x),                                               \
             (x) = (xnext)                                                \
)

/* CCpq_set_FOREACH_DEL(CCpq_set s, CCpq_node *x, ELEM_FIELD efield,
   CCpq_node *xprev, CCpq_node *xnext) iterates x over
   elements of s using temporary variables xprev and xnext.  x may be deleted
   from the set in the body of the loop (but delete any other element from s
   at your own risk). */

#define CCpq_set_FOREACH_DEL(s,x,efield,xprev,xnext)                      \
        for ((xprev) = (CCpq_node *) NULL,                                \
             (x) = (s).left,                                              \
             (xnext) = ((x) ? CCpq_set_PTR_AWAY((x)->efield,              \
                                                (CCpq_node *) NULL)       \
                            : (CCpq_node *) NULL);                        \
            (x);                                                          \
             (xprev) = (((xprev) ? (((xprev)->efield.ptr1 == (x))         \
                                 || ((xprev)->efield.ptr2 == (x)))        \
                                 : ((s).left == (x)))                     \
                        ? (x) : (xprev)),                                 \
             (x) = (xnext),                                               \
             (xnext) = ((x) ? CCpq_set_PTR_AWAY((x)->efield, xprev)       \
                            : NULL)                                       \
)

/* CCpq_set_FOREACH_ADJ(CCpq_node *x, ELEM_FIELD efield, CCpq_node *z,
   int itemp) iterates z over the immediate neighbors of x using itemp as
   a temporary variable. It is just used to save code replication when you
   want to do something for each of the two neighbors. */

#define CCpq_set_FOREACH_ADJ(x,efield,z,itemp)                            \
        for (z = x->efield.ptr1, itemp = 0;                               \
             itemp < 2;                                                   \
             z = x->efield.ptr2, itemp++)

#endif  /* __PQSETS_H */
