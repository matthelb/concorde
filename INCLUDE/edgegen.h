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
/****************************************************************************/
/*                                                                          */
/*                      PROTOTYPES FOR FILES IN EDGEGEN                     */
/*                                                                          */
/****************************************************************************/
/****************************************************************************/

#ifndef __EDGEGEN_H
#define __EDGEGEN_H

#include "util.h"


/****************************************************************************/
/*                                                                          */
/*                             edgegen.c                                    */
/*                                                                          */
/****************************************************************************/

typedef struct CCedgegengroup {
    struct {
        int count;
        int quadnearest;
        int nearest;
        int nearest_start;
        int greedy_start;
        int boruvka_start;
        int qboruvka_start;
        int random_start;
        int nkicks;
    } linkern;

    struct {
        int twoopt_count;
        int twoopt5_count;
        int threeopt_count;
        int greedy;
        int boruvka;
        int qboruvka;
        int nearest_count;
        int random_count;
    } tour;

    struct {
        int wantit;
        int basic;
        int priced;
    } f2match;

    struct {
        int number;
        int basic;
        int priced;
    } f2match_nearest;

    int    random;
    int    nearest;
    int    quadnearest;
    int    want_tree;
    int    nearest_twomatch_count;
    int    delaunay;
    int    mlinkern;
} CCedgegengroup;



int
    CCedgegen_read (char *egname, CCedgegengroup *plan),
    CCedgegen_edges (CCedgegengroup *plan, int ncount, CCdatagroup *dat,
        double *wcoord, int *ecount, int **elist, int silent,
        CCrandstate *rstate);
void
    CCedgegen_init_edgegengroup (CCedgegengroup *plan);




/****************************************************************************/
/*                                                                          */
/*                             xnear.c                                      */
/*                                                                          */
/****************************************************************************/

typedef struct CCxnear {
    struct CCdatagroup dat;
    double            *w;
    int               *nodenames;
    int               *invnames;
} CCxnear;



int
    CCedgegen_x_k_nearest (int ncount, int num, CCdatagroup *dat,
        double *wcoord, int wantlist, int *ecount, int **elist, int silent),
    CCedgegen_x_quadrant_k_nearest (int ncount, int num, CCdatagroup *dat,
        double *wcoord, int wantlist, int *ecount, int **elist, int silent),
    CCedgegen_x_node_k_nearest (CCxnear *xn, int n, int nearnum, int ncount,
        int *list),
    CCedgegen_x_node_quadrant_k_nearest (CCxnear *xn, int n, int nearnum,
        int ncount, int *list),
    CCedgegen_x_node_nearest (CCxnear *xn, int ncount, int ni, char *marks),
    CCedgegen_x_nearest_neighbor_tour (int ncount, int start, CCdatagroup *dat,
        int *outcycle, double *val),
    CCedgegen_x_greedy_tour (int ncount, CCdatagroup *dat, int *outcycle,
        double *val, int ecount, int *elist, int silent),
    CCedgegen_x_qboruvka_tour (int ncount, CCdatagroup *dat, int *outcycle,
        double *val, int ecount, int *elist, int silent),
    CCedgegen_junk_k_nearest (int ncount, int num, CCdatagroup *dat,
        double *wcoord, int wantlist, int *ecount, int **elist, int silent),
    CCedgegen_junk_node_k_nearest (CCdatagroup *dat, double *wcoord, int n,
        int nearnum, int ncount, int *list),
    CCedgegen_junk_node_nearest (CCdatagroup *dat, double *wcoord, int ncount,
        int n, char *marks),
    CCedgegen_junk_nearest_neighbor_tour (int ncount, int start,
        CCdatagroup *dat, int *outcycle, double *val, int silent),
    CCedgegen_junk_greedy_tour (int ncount, CCdatagroup *dat, int *outcycle,
        double *val, int ecount, int *elist, int silent),
    CCedgegen_junk_qboruvka_tour (int ncount, CCdatagroup *dat, int *outcycle,
        double *val, int ecount, int *elist, int silent),
    CCedgegen_xnear_build (int ncount, CCdatagroup *dat, double *wcoord,
        CCxnear *xn);

void
    CCedgegen_xnear_free (CCxnear *xn);


#endif  /* __EDGEGEN_H */
