#!/usr/bin/nawk -f
#
#   This file is part of CONCORDE
#
#   (c) Copyright 1995--1999 by David Applegate, Robert Bixby,
#   Vasek Chvatal, and William Cook
#
#   Permission is granted for academic research use.  For other uses,
#   contact the authors for licensing options.
#
#   Use at your own risk.  We make no guarantees about the
#   correctness or usefulness of this code.
#

BEGIN{
    if (N == "") N=10;
    for (i=0; i<N; i++) {
        x[i] = int(rand() * 1000);
        y[i] = int(rand() * 1000);
    }
    print N, N*(N-1)/2;
    for (i=0; i<N; i++) {
        for (j=i+1; j<N; j++) {
            print i, j, int(sqrt((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j])));
        }
    }
}

    
