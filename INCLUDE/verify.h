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

#ifndef __VERIFY_H
#define __VERIFY_H

#include "tsp.h"

#define CC_TYPE_SUBTOUR 1
#define CC_TYPE_COMB 2
#define CC_TYPE_STAR 4
#define CC_TYPE_BIPARTITION 8
#define CC_TYPE_OTHER 16
#define CC_TYPE_ALL    (CC_TYPE_SUBTOUR     | CC_TYPE_COMB  | CC_TYPE_STAR | \
                        CC_TYPE_BIPARTITION | CC_TYPE_OTHER)
#define CC_TYPE_SIMPLE (CC_TYPE_SUBTOUR     | CC_TYPE_COMB  | CC_TYPE_STAR | \
                        CC_TYPE_BIPARTITION)

typedef struct CCverify_cutclass {
    int type;
    int nhandles;
    int nfamilies;
    int *cliques;
    int *inverted;
    int *family_start;
} CCverify_cutclass;


int
    CCverify_cut (CCtsp_lpcut_in *cut, int check_types, int *type),
    CCverify_classify (CCtsp_lpcut_in *cut, int check_types,
        CCverify_cutclass *class);

void
    CCverify_initcutclass (CCverify_cutclass *class),
    CCverify_freecutclass (CCverify_cutclass *class);


#endif  /* __VERIFY_H */
