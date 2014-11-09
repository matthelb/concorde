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
/*                   CODE TO READ ASCI INTS FAST                            */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*                                                                          */
/*  Written by:  Dave                                                       */
/*  Date: Septmember 1994 (Bonn)  (cofeb16)                                 */
/*                                                                          */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCutil_readint (FILE *f)                                            */
/*      - Returns the next int in the file f.                               */
/*    This is much faster that scanf. It is useful for big files and        */
/*    and for profiling.                                                    */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"

int CCutil_readint (FILE *f)
{
    int v = 0;
    int c;

    while (( c = getc(f)) != EOF && !((c >= '0' && c <= '9') || c == '-'));
    if (c == '-') {
        v = 0;
        while ((c = getc(f)) != EOF && c >= '0' && c <= '9') {
            v = v * 10 + c - '0';
        }
        return -v;
    } else {
        v = c - '0';
        while ((c = getc(f)) != EOF && c >= '0' && c <= '9') {
            v = v * 10 + c - '0';
        }
        return v;
    }
}
