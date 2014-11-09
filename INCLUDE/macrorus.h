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
/*  CC_SWAP(a,b,t)                                                          */
/*    swaps a and b, using t as temporary space.  a, b, and t should all    */
/*    be the same type.                                                     */
/*                                                                          */
/*  CC_OURABS(a)                                                            */
/*    returns the absolute value of a.                                      */
/*                                                                          */
/****************************************************************************/

#ifndef  __MACRORUS_H
#define  __MACRORUS_H

#define CC_SWAP(a,b,t) (((t)=(a)),((a)=(b)),((b)=(t)))

#define CC_OURABS(a) (((a) >= 0) ? (a) : -(a))

#endif  /* __MACRORUS_H */
