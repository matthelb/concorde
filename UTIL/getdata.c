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
/*                   SOME DATA READING ROUTINES                             */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: March 2, 1995                                                     */
/*  Changes: 17.7.96 (Bico)                                                 */
/*                                                                          */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCutil_getdata (char *datname, int binary_in, int innorm,           */
/*      int *ncount, CCdatagroup *dat, int gridsize, int allow_dups,        */
/*      CCrandstate *rstate)                                                */
/*    RETURNS the data to generate edge lengths in the dat structure.       */
/*            The calling routine should be sure that dat points to         */
/*            a structure. If datname is NULL then random entries will be   */
/*            generated.                                                    */
/*     -datname is the name of the datfile or the matrix file, if NULL      */
/*      then random data will be generated, according to the norm type.     */
/*      For D2 and D3 norms, the coordinates will be uniform between 0      */
/*      and ncount -1 (GEOGRAPHIC and GEOM norms have x between -90 and 90  */
/*      and y between -180 and 180). (For D2, the points will be distinct.) */
/*      For MATRIX norms, the entries will be                               */
/*      uniform between 0 and MATRAND_SCALE * ncount - 1 (currently         */
/*      10*ncount - 1. For CRYSTAL norms, a random matrix and bounds in     */
/*      range of the TSPLIB problems is generated - the wavelength is       */
/*      chosen to be 1.0, 1.35, or 1.70 depending on the ncount (but the    */
/*      problem will not be very close to hitting ncount.).                 */
/*     -binary_in should be 1 if the datname file is in binary integers,    */
/*      or 2 if the datname file is in binary doubles.                      */
/*     -innorm is the norm.                                                 */
/*     -ncount will return the number of nodes. If datname is NULL, then    */
/*      ncount should be passed in with the number of nodes to be used in   */
/*      the random problem generation.                                      */
/*     -dat will contain the info to call the edgelen function.             */
/*     -gridsize specifies the size of the square grid random points are    */
/*      drawn from.                                                         */
/*     -allow_dups indicates whether or not to allow duplicate points in    */
/*      the random point set.                                               */
/*                                                                          */
/*  int CCutil_putmaster (char *mastername, int ncount, CCdatagroup *dat,   */
/*      int *perm)                                                          */
/*    WRITES the dat information and the permutation into a binary file.    */
/*           This is used in the TSP, where the dat file has usually        */
/*           been permuted to put the nodes into tour order.                */
/*     -mastername is the name of the file (cannot be NULL)                 */
/*     -ncount is the number of nodes                                       */
/*     -dat contains the edgelen info (e.g. x,y coordinates), it can be     */
/*      NULL                                                                */
/*     -perm contains a permutation of 0 to ncount - 1 (so a tour in        */
/*      node node node format)                                              */
/*                                                                          */
/*  int CCutil_writemaster (CC_SFILE *out, int ncount, CCdatagroup *dat,    */
/*      int *perm)                                                          */
/*    WRITES the dat information and the permutation into a binary file.    */
/*           This is used in the TSP, where the dat file has usually        */
/*           been permuted to put the nodes into tour order.                */
/*     -f is the CC_SFILE to write into                                     */
/*     -ncount is the number of nodes                                       */
/*     -dat contains the edgelen info (e.g. x,y coordinates), it can be     */
/*      NULL                                                                */
/*     -perm contains a permutation of 0 to ncount - 1 (so a tour in        */
/*      node node node format)                                              */
/*                                                                          */
/*  int CCutil_getmaster (char *mastername, int *ncount, CCdatagroup *dat,  */
/*      int **perm)                                                         */
/*    RETURNS the dat information and the permutation from a binary file    */
/*            (written by a call to CCutil_writemaster). Used by the TSP    */
/*            code.                                                         */
/*     -mastername is the name of the file (cannot be NULL)                 */
/*     -ncount returns the number of nodes                                  */
/*     -dat returns the edgelen info (e.g. x,y coordinates), or NULL        */
/*     -perm returns a permutation of 0 to ncount - 1 (so a tour in         */
/*      node node node format)                                              */
/*                                                                          */
/*  int CCutil_readmaster (CC_SFILE *in, int *ncount, CCdatagroup *dat,     */
/*      int **perm)                                                         */
/*    RETURNS the dat information and the permutation from a binary file    */
/*            (written by a call to CCutil_writemaster). Used by the TSP    */
/*            code.                                                         */
/*     -f is the CC_SFILE to read from                                      */
/*     -ncount returns the number of nodes                                  */
/*     -dat returns the edgelen info (e.g. x,y coordinates), or NULL        */
/*     -perm returns a permutation of 0 to ncount - 1 (so a tour in         */
/*      node node node format)                                              */
/*                                                                          */
/*  int CCutil_getnodeweights (char *weightname, int ncount,                */
/*      int weight_limit, double **wcoord, CCrandstate *rstate)             */
/*    RETURNS a list of nonnegative nodesweights (they will be translated   */
/*            if the weightname file has negative entries). If weightname   */
/*            is NULL then random entries (from 0 to ncount - 1 will be     */
/*            generated.                                                    */
/*     -weightname is the name of a file containing nodeweights. It will    */
/*      only be read if inusenodeweights is 1. If NULL, random entries      */
/*      will be generated (between 0 and weight_limit - 1)                  */
/*     -wcoord will contain the nonnegative node weights (it will be        */
/*      here in CCutil_getnodeweights)                                      */
/*                                                                          */
/*  int CCutil_gettsplib (char *datname, int *ncount, CCdatagroup *dat)     */
/*    READS an xxx.tsp TSPLIB file, and returns the dat structure to        */
/*            generate edge lengths.                                        */
/*     -datname should be the name of a TSPLIB xxx.tsp file.                */
/*     -ncount returns the number of nodes.                                 */
/*     -dat returns the data.                                               */
/*                                                                          */
/*  int CCutil_writetsplib (const char *fname, int ncount,                  */
/*      CCdatagroup *dat)                                                   */
/*    WRITES a TSPLIB file for the instance specified by dat.               */
/*                                                                          */
/*  int CCutil_getedgelist (int ncount, char *fname, int *ecount,           */
/*      int **elist, int **elen, int binary_in)                             */
/*    READS an edgelist in end1 end2 length format.                         */
/*     -fname name of the file                                              */
/*     -ecount returns the number of edges                                  */
/*     -elist returns the edges in end1 end2 format (it will be allocated   */
/*      by CCutil_getedgelist)                                              */
/*     -elen returns the length of the edges in len len len format          */
/*                                                                          */
/*  int CCutil_getedgelist_n (int *ncount, char *fname, int *ecount,        */
/*      int **elist, int **elen, int binary_in)                             */
/*    READS an edgelist in end1 end2 length format.                         */
/*    Like CCutil_getedgelist (), but it also returns ncount.               */
/*                                                                          */
/*  int CCutil_genedgelist (int ncount, int ecount, int **elist,            */
/*      int **elen, CCdatagroup *dat, int maxlen, CCrandstate *rstate)      */
/*    GENERATES a random graph with ncount nodes and ecount edges, with     */
/*     with edgelengths either determined by dat or random between 0 and    */
/*     maxlen-1.                                                            */
/*     -dat specifies the function for computing the edge lengths (it can   */
/*      be NULL.                                                            */
/*     -maxlen gives the range for random edge lengths (it is used if dat   */
/*      is NULL.                                                            */
/*                                                                          */
/*  int CCutil_getcycle_tsplib (int ncount, char *cyclename, int *incycle)  */
/*    READS a cycle in TSPLIB TOUR format and returns the cycle in node     */
/*            node format in the array outcycle.                            */
/*     -outcycle should be allocated by the calling routine (and should     */
/*      be at least ncount long)                                            */
/*                                                                          */
/*  int CCutil_getcycle_edgelist (int ncount, char *cyclename,              */
/*      int *incycle, int binary_in)                                        */
/*    READS a cycle in end1 end2 length format, and returns the cycle in    */
/*            node node format in the array outcycle.                       */
/*     -outcycle should be allocated by the calling routine (and should     */
/*      be at least ncount long)                                            */
/*                                                                          */
/*  int CCutil_getcycle (int ncount, char *cyclename, int *incycle,         */
/*      int binary_in)                                                      */
/*    READS a cycle in node node format, and returns the cycle in node      */
/*            node format in the array outcycle.                            */
/*     -incycle should be allocated by the calling routine                  */
/*                                                                          */
/*  void CCutil_cycle_len (int ncount, CCdatagroup *dat, int *cycle,        */
/*      double *len)                                                        */
/*    COMPUTES the length of a cycle (in permutation format).               */
/*                                                                          */
/*  int CCutil_getedges_double (int *ncount, char *fname, int *ecount,      */
/*      int **elist, double **elen, int binary_in)                          */
/*    READS the edgelist in end1 end2 length format, where the length       */
/*            is a double (used for x-vectors)                              */
/*     -returns the number of nodes                                         */
/*     -fname is the name of the file to read.                              */
/*     -ecount returns the number of edges.                                 */
/*     -elist returns the list of edges in end1 end2 format (it will be     */
/*      allocated by CCutil_getedges_double).                               */
/*     -elen returns the lengths of the edges.                              */
/*     -binary_in indicates whether the file should be read in binary       */
/*      or in ascii (1 is binary, 0 is ascii)                               */
/*                                                                          */
/*  int CCutil_writeedges (int ncount, char *outedgename, int ecount,       */
/*      int *elist, CCdatagroup *dat, int binary_out)                       */
/*    WRITES the edgelist in end1 end2 length format.                       */
/*     -ncount the number of nodes                                          */
/*     -outedgename is the name of the file to write to.                    */
/*     -ecount is the number of edges.                                      */
/*     -elist is the list of edges in end1 end2 format.                     */
/*     -dat contains the data to compute edgelengths.                       */
/*     -binary_in indicates whether the file should be written in binary    */
/*      or in ascii (1 is binary, 0 is ascii)                               */
/*                                                                          */
/*  int CCutil_writecycle_edgelist (int ncount, char *outedgename,          */
/*      int *cycle, CCdatagroup *dat, int binary_out)                       */
/*    WRITES the cycle in edgelist format.                                  */
/*                                                                          */
/*  int CCutil_writecycle (int ncount, char *outcyclename, int *cycle,      */
/*      int binary_out)                                                     */
/*    WRITES the cycle in node node node format.                            */
/*                                                                          */
/*  int CCutil_writeedges_int (int ncount, char *outedgename, int ecount,   */
/*      int *elist, int *elen, int binary_out)                              */
/*    WRITES the edgelist in end1 end2 length format.                       */
/*    Like CCutil_writeedges, but lengths are specified in elen.            */
/*                                                                          */
/*  int CCutil_writeedges_double (int ncount, char *outedgename,            */
/*      int ecount, int *elist, double *elen, int binary_out)               */
/*    WRITES the edgelist in end1 end2 length format, where the length      */
/*            is a double (used for x-vectors)                              */
/*     -ncount the number of nodes                                          */
/*     -outedgename is the name of the file to write to.                    */
/*     -ecount is the number of edges.                                      */
/*     -elist is the list of edges in end1 end2 format.                     */
/*     -elen are the lengths of the edges.                                  */
/*     -binary_out indicates whether the file should be written in binary   */
/*        or in ascii (1 is binary, 0 is ascii)                             */
/*                                                                          */
/*  int CCutil_datagroup_perm (int ncount, CCdatagroup *dat, int *perm)     */
/*    reorders the nodes to match the order given in perm.                  */
/*                                                                          */
/*  int CCutil_copy_datagroup (int ncount, CCdatagroup *indat,              */
/*      CCdatagroup *outdat)                                                */
/*    COPIES indat to outdat.                                               */
/*                                                                          */
/*  int CCutil_tri2dat (int ncount, int *elen, CCdatagroup *dat)            */
/*     convert an upper triangular matrix of edge lengths to a CCdatagroup  */
/*                                                                          */
/*  int CCutil_graph2dat_matrix (int ncount, int ecount, int *elist,        */
/*      int *elen, int defaultlen, CCdatagroup *dat)                        */
/*     convert a matrix of edge lengths to a CCdatagroup                    */
/*                                                                          */
/*  int CCutil_graph2dat_sparse (int ncount, int ecount, int *elist,        */
/*      int *elen, int defaultlen, CCdatagroup *dat)                        */
/*    BUILDS a datagroup from the edge list using the CC_SPARSE norm.       */
/*     -dat will be filled in by the function.                              */
/*     -defaultlen gives the value for the non-edges (if it is <= 0, then   */
/*      a default value of (max_edgelen + 1) * nnodes if fiven).            */
/*                                                                          */
/*  int CCutil_get_sparse_dat_edges (int ncount, CCdatagroup *dat,          */
/*      int *ecount, int **elist, int **elen)                               */
/*    GRABS the enitre set of edges in a CC_SPARSE datagroup.  It returns   */
/*      an error if the datagroup norm is not CC_SPARSE.                    */
/*                                                                          */
/*  int CCutil_sparse_strip_edges (CCdatagroup *dat, int in_ecount,         */
/*      int *in_elist, int *in_elen, int *ecount, int **elist, int **elen)  */
/*    REMOVES the artifical edges from in_elist.                            */
/*     -ecount, elist, and elen return the stripped edgelist.               */
/*                                                                          */
/*  int CCutil_sparse_real_tour (int ncount, CCdatagroup *dat, int *cyc,    */
/*      int *yesno)                                                         */
/*    CHECKS whether each edge of cyc is in the CC_SPARSE datagroup.        */
/*     -cyc gives a tour in node node node format.                          */
/*     -yesno will be set to 1 if it is a tour, and 0 otherwise.            */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/*  NOTES:                                                                  */
/*     Functions prototyped in util.h. Functions return 0 when they         */
/*    succeed and nonzero when they fail (usually do to bad filenames or    */
/*    not enough memory).                                                   */
/*     The TSPLIB reader works for all problems in TSPLIB_1.2, but does     */
/*    not include all of the options listed in Reinelt's orginal TSPLIB     */
/*    paper. It returns a failure on linhp318.tsp, since there is no        */
/*    place for fixed edges in our edge length dat structure.               */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "macrorus.h"

#define MATRAND_SCALE 10  /* Range of edge lengths: [0, MATRAND_SCALE * n) */
#define SPARSE_ECOUNT  5  /* # of edges in random graph:  SPARSE_ECOUNT * n */

static void
    make_weights_nonnegative (int ncount, double *wcoord);

static int
    read_crystal (char *datname, int binary_in, int *ncount, CCdatagroup *dat,
        CCrandstate *rstate),
    read_d2 (char *datname, int binary_in, int *ncount,
        CCdatagroup *dat, CCrandstate *rstate, int innorm, int gridsize,
        int allow_dups),
    write_d2 (char *datname, int binary_out, int ncount, CCdatagroup *dat),
    build_d2 (int *ncount, CCdatagroup *dat, CCrandstate *rstate,
        int innorm, int gridsize, int allow_dups),
    read_d2_binary (char *datname, int *ncount, CCdatagroup *dat, int btype),
    read_d2_text (char *datname, int *ncount, CCdatagroup *dat),
    write_d2_binary (char *datname, int ncount, CCdatagroup *dat),
    write_d2_text (char *datname, int ncount, CCdatagroup *dat),
    read_d3 (char *datname, int binary_in, int *ncount,
        CCdatagroup *dat, CCrandstate *rstate, int gridsize),
    write_d3 (char *datname, int binary_out, int ncount, CCdatagroup *dat),
    build_d3 (int *ncount, CCdatagroup *dat, CCrandstate *rstate,
        int gridsize),
    read_d3_binary (char *datname, int *ncount, CCdatagroup *dat),
    read_d3_text (char *datname, int *ncount, CCdatagroup *dat),
    write_d3_text (char *datname, int ncount, CCdatagroup *dat),
    read_matrix (char *datname, int binary_in, int *ncount,
        CCdatagroup *dat, CCrandstate *rstate),
    write_matrix (char *datname, int binary_out, int ncount, CCdatagroup *dat),
    build_matrix (int *ncount, CCdatagroup *dat, CCrandstate *rstate),
    read_matrix_binary (char *datname, int *ncount, CCdatagroup *dat),
    read_matrix_text (char *datname, int *ncount, CCdatagroup *dat),
    write_matrix_text (char *datname, int ncount, CCdatagroup *dat),
    read_dsjrand (char *datname, int binary_in, int *ncount,
        CCdatagroup *dat),
    read_sparse (char *datname, int binary_in, int *ncount,
        CCdatagroup *dat, CCrandstate *rstate),
    write_sparse (char *datname, int binary_out, int ncount, CCdatagroup *dat),
    build_sparse (int *ncount, CCdatagroup *dat, CCrandstate *rstate),
    read_sparse_text (char *datname, int *ncount, CCdatagroup *dat,
        int binary_in),
    write_sparse_text (char *datname, int ncount, CCdatagroup *dat),
    read_user (char *datname, int binary_in, int *ncount,
        CCdata_user *userdat, CCrandstate *rstate, int gridsize,
        int allow_dups),
    build_user (int *ncount, CCdata_user *userdat, CCrandstate *rstate,
                int gridsize, int allow_dups),
    read_user_binary (char *datname, int *ncount, CCdata_user *userdat),
    read_user_text (char *datname, int *ncount, CCdata_user *userdat),
    permute_user (int ncount, CCdata_user *userdat, int *perm),
    writemaster_user (CC_SFILE *out, int ncount, CCdata_user *userdat),
    readmaster_user (CC_SFILE *in, int ncount, CCdata_user *userdat),
    read_rhdata (char *datname, int innorm, int binary_in, int *ncount,
        CCdata_rhvector *dat),
    read_rhdata_text (char *datname, int innorm, int *ncount,
        CCdata_rhvector *dat),
    writemaster_rhdata (CC_SFILE *out, int ncount, CCdata_rhvector *dat),
    readmaster_rhdata (CC_SFILE *in, int ncount, CCdata_rhvector *dat),
    permute_rhdata (int ncount, CCdata_rhvector *dat, int *perm),
    build_sparse_dat (int ncount, int ecount, int *elist, int *elen,
        CCdatagroup *dat, int defaultlen),
    writedat_sparse (CC_SFILE *out, int ncount, CCdatagroup *dat),
    readdat_sparse (CC_SFILE *in, int ncount, CCdatagroup *dat),
    permute_sparse (int ncount, CCdatagroup *dat, int *perm),
    copy_sparse (int ncount, CCdatagroup *indat, CCdatagroup *outdat);


int CCutil_getdata (char *datname, int binary_in, int innorm,
        int *ncount, CCdatagroup *dat, int gridsize, int allow_dups,
        CCrandstate *rstate)
{
    CCutil_init_datagroup (dat);
    CCutil_dat_setnorm (dat, innorm);

    if (datname == (char *) NULL && *ncount == 0) {
        fprintf (stderr, "CCutil_getdata needs a datfile or a nodecount\n");
        return 1;
    }

    if (innorm == CC_CRYSTAL) {
        return read_crystal (datname, binary_in, ncount, dat, rstate);
    } else if ((innorm & CC_NORM_SIZE_BITS) == CC_D2_NORM_SIZE) {
        return read_d2 (datname, binary_in, ncount, dat, rstate, innorm,
                        gridsize, allow_dups);
    } else if ((innorm & CC_NORM_SIZE_BITS) == CC_D3_NORM_SIZE) {
        return read_d3 (datname, binary_in, ncount, dat, rstate, gridsize);
    } else if ((innorm & CC_NORM_SIZE_BITS) == CC_MATRIX_NORM_SIZE){
        return read_matrix (datname, binary_in, ncount, dat, rstate);
    } else if (innorm == CC_DSJRANDNORM) {
        return read_dsjrand (datname, binary_in, ncount, dat);
    } else if (innorm == CC_USER) {
        return read_user (datname, binary_in, ncount, &dat->userdat, rstate,
                          gridsize, allow_dups);
    } else if (innorm == CC_SPARSE) {
        return read_sparse (datname, binary_in, ncount, dat, rstate);
    } else if (innorm == CC_RHMAP1 || innorm == CC_RHMAP2 ||
               innorm == CC_RHMAP3 || innorm == CC_RHMAP4 ||
               innorm == CC_RHMAP5) {
        return read_rhdata (datname, innorm, binary_in, ncount,
                            &dat->rhdat);
    } else {
        fprintf (stderr, "Unknown norm %d\n", innorm);
        return 1;
    }
}

int CCutil_writedata (char *datname, int binary_out, int ncount,
        CCdatagroup *dat)
{
    int norm;

    CCutil_dat_getnorm (dat, &norm);

    if ((norm & CC_NORM_SIZE_BITS) == CC_D2_NORM_SIZE) {
        return write_d2 (datname, binary_out, ncount, dat);
    } else if ((norm & CC_NORM_SIZE_BITS) == CC_D3_NORM_SIZE) {
        return write_d3 (datname, binary_out, ncount, dat);
    } else if ((norm & CC_NORM_SIZE_BITS) == CC_MATRIX_NORM_SIZE){
        return write_matrix (datname, binary_out, ncount, dat);
    } else if (norm == CC_SPARSE) {
        return write_sparse (datname, binary_out, ncount, dat);
    } else {
        fprintf (stderr, "Output of this norm not yet implemented\n");
        return 1;
    }
}

int CCutil_putmaster (char *mastername, int ncount, CCdatagroup *dat,
        int *perm)
{
    CC_SFILE *out = (CC_SFILE *) NULL;

    if (mastername == (char *) NULL) {
        fprintf (stderr, "CCutil_writemaster needs a filename\n");
        return 1;
    }

    if (!dat) {
        fprintf (stderr, "Cannot put a master without a datagroup\n");
        return 1;
    }
    if (!perm) {
        fprintf (stderr, "Cannot put a master without a permutation\n");
        return 1;
    }

    out = CCutil_sopen (mastername, "w");
    if (out == (CC_SFILE *) NULL) {
        fprintf (stderr, "Unable to open %s for output\n", mastername);
        return 1;
    }

    if (CCutil_writemaster (out, ncount, dat, perm)) {
        fprintf (stderr, "CCutil_writemaster failed\n");
        CCutil_sclose (out);
        return 1;
    }

    CCutil_sclose (out);
    return 0;
}

int CCutil_writemaster (CC_SFILE *out, int ncount, CCdatagroup *dat, int *perm)
{
    int i, j;
    int norm;
    int ndepot = dat->ndepot;

    if (!dat) {
        fprintf (stderr, "Cannot write a master without a datagroup\n");
        return 1;
    }
    if (!perm) {
        fprintf (stderr, "Cannot write a master without a permutation\n");
        return 1;
    }

    CCutil_dat_getnorm (dat, &norm);

    if (ndepot > 0) {
        ncount = ncount - ndepot;
    }
    
    if (CCutil_swrite_int (out, ncount)) {
        return 1;
    }

    if (CCutil_swrite_int (out, CC_MASTER_DAT)) {
        return 1;
    }

    if (ndepot > 0) {
        if (CCutil_swrite_int (out, CC_SUBDIVISION)) {
            return 1;
        }
    }

    if (CCutil_swrite_int (out, norm)) {
        return 1;
    }

    if (ndepot > 0) {
        if (CCutil_swrite_int (out, ndepot)) {
            return 1;
        }
        for (i = 0; i < ncount; i++) {
            if (CCutil_swrite_int (out, dat->depotcost[i])) {
                return 1;
            }
        }
    }

    if ((norm & CC_NORM_SIZE_BITS) == CC_D2_NORM_SIZE) {
        for (i = 0; i < ncount; i++) {
            if (CCutil_swrite_double (out, dat->x[i])) {
                return 1;
            }
            if (CCutil_swrite_double (out, dat->y[i])) {
                return 1;
            }
        }
    } else if ((norm & CC_NORM_SIZE_BITS) == CC_D3_NORM_SIZE) {
        for (i = 0; i < ncount; i++) {
            if (CCutil_swrite_double (out, dat->x[i])) {
                return 1;
            }
            if (CCutil_swrite_double (out, dat->y[i])) {
                return 1;
            }
            if (CCutil_swrite_double (out, dat->z[i])) {
                return 1;
            }
        }
    } else if ((norm & CC_NORM_SIZE_BITS) == CC_MATRIX_NORM_SIZE){
        /* Matrix is the lower triangle plus the diagonal */
        for (i = 0; i < ncount; i++) {
            for (j = 0; j <= i; j++) {
                if (CCutil_swrite_int (out, dat->adj[i][j])) {
                    return 1;
                }
            }
        }
    } else if (norm == CC_DSJRANDNORM) {
        for (i = 0; i < ncount; i++) {
            if (CCutil_swrite_double (out, dat->x[i])) {
                return 1;
            }
        }
    } else if (norm == CC_SPARSE) {
        if (writedat_sparse (out, ncount, dat)) {
            fprintf (stderr, "writedat_sparse failed\n");
            return 1;
        }
    } else if (norm == CC_USER) {
        if (writemaster_user (out, ncount, &dat->userdat)) {
            fprintf (stderr, "writemaster_user failed\n");
            return 1;
        }
    } else if (norm == CC_RHMAP1 || norm == CC_RHMAP2 || norm == CC_RHMAP3 ||
               norm == CC_RHMAP4 || norm == CC_RHMAP5) {
        if (writemaster_rhdata (out, ncount, &dat->rhdat)) {
            fprintf (stderr, "writemaster_rhdata failed\n");
            return 1;
        }
    } else {
        fprintf (stderr, "unknown norm: %d\n", norm);
        return 1;
    }

    for (i = 0; i < ncount; i++) {
        if (perm[i] < 0 || perm[i] >= ncount) {
            fprintf (stderr, "permutation in wrong format\n");
            return 1;
        }
        if (CCutil_swrite_int (out, perm[i])) {
            return 1;
        }
    }

    if (ndepot > 0) {
        for (i = 0; i < ncount; i++) {
            if (CCutil_swrite_int (out, dat->orig_names[i])) {
                return 1;
            }
        }
    }

    return 0;
}

int CCutil_getmaster (char *mastername, int *ncount, CCdatagroup *dat,
        int **perm)
{
    CC_SFILE *in = (CC_SFILE *) NULL;

    if (mastername == (char *) NULL) {
        fprintf (stderr, "CCutil_getmaster needs a filename\n");
        return 1;
    }

    in = CCutil_sopen (mastername, "r");
    if (in == (CC_SFILE *) NULL) {
        fprintf (stderr, "Unable to open %s for input\n", mastername);
        return 1;
    }

    if (CCutil_readmaster (in, ncount, dat, perm)) {
        fprintf (stderr, "CCutil_readmaster failed\n");
        CCutil_sclose (in);
        return 1;
    }

    CCutil_sclose (in);
    return 0;
}

int CCutil_readmaster (CC_SFILE *in, int *ncount, CCdatagroup *dat, int **perm)
{
    int i, j;
    int havedat = 0;
    int norm, subdivision = 0, ndepot = 0;
    int rval = 0;

    *ncount = 0;
    CCutil_init_datagroup (dat);
    *perm = (int *) NULL;

    rval = CCutil_sread_int (in, ncount);
    CCcheck_rval (rval, "CCutil_sread_int failed");

    rval = CCutil_sread_int (in, &havedat);
    CCcheck_rval (rval, "CCutil_sread_int failed");

    if (havedat != CC_MASTER_DAT) {
        fprintf (stderr, "masterfile does not have a dat section\n");
        rval = 1; goto CLEANUP;
    }

    rval = CCutil_sread_int (in, &norm);
    CCcheck_rval (rval, "CCutil_sread_int failed");

    if (norm == CC_SUBDIVISION) {
        subdivision = 1;
        rval = CCutil_sread_int (in, &norm);
        CCcheck_rval (rval, "CCutil_sread_int failed");

        rval = CCutil_sread_int (in, &ndepot);
        CCcheck_rval (rval, "CCutil_sread_int failed");

        dat->ndepot = ndepot;

        dat->depotcost = CC_SAFE_MALLOC ((*ncount) + ndepot, int);
        CCcheck_NULL (dat->depotcost, "out of memory in CCutil_getmaster");
        for (i = 0; i < *ncount; i++) {
            rval = CCutil_sread_int (in, &dat->depotcost[i]);
            CCcheck_rval (rval, "CCutil_sread_int failed");
        }
        for (i = *ncount; i < *ncount + ndepot; i++) {
            dat->depotcost[i] = 0;
        }
    }

    if (CCutil_dat_setnorm (dat, norm)) {
        goto CLEANUP;
    }

    if ((norm & CC_NORM_SIZE_BITS) == CC_D2_NORM_SIZE) {
        dat->x = CC_SAFE_MALLOC (*ncount + ndepot, double);
        dat->y = CC_SAFE_MALLOC (*ncount + ndepot, double);

        CCcheck_NULL (dat->x, "out of memory in CCutil_getmaster");
        CCcheck_NULL (dat->y, "out of memory in CCutil_getmaster");

        for (i = 0; i < *ncount; i++) {
            rval = CCutil_sread_double (in, &(dat->x[i]));
            CCcheck_rval (rval, "CCutil_sread_double failed");
            rval = CCutil_sread_double (in, &(dat->y[i]));
            CCcheck_rval (rval, "CCutil_sread_double failed");
        }
        for (i = 0; i < ndepot; i++) {
            dat->x[*ncount + i] = 0.0;
            dat->y[*ncount + i] = 0.0;
        }
        if (norm == CC_EUCTOROIDAL) {
            dat->gridsize = 1000000;
            printf ("ASSUMING GRIDSIZE 1000000\n");
        }
    } else if ((norm & CC_NORM_SIZE_BITS) == CC_D3_NORM_SIZE) {
        dat->x = CC_SAFE_MALLOC (*ncount + ndepot, double);
        dat->y = CC_SAFE_MALLOC (*ncount + ndepot, double);
        dat->z = CC_SAFE_MALLOC (*ncount + ndepot, double);

        CCcheck_NULL (dat->x, "out of memory in CCutil_getmaster");
        CCcheck_NULL (dat->y, "out of memory in CCutil_getmaster");
        CCcheck_NULL (dat->z, "out of memory in CCutil_getmaster");

        for (i = 0; i < *ncount; i++) {
            rval = CCutil_sread_double (in, &(dat->x[i]));
            CCcheck_rval (rval, "CCutil_sread_double failed");
            rval = CCutil_sread_double (in, &(dat->y[i]));
            CCcheck_rval (rval, "CCutil_sread_double failed");
            rval = CCutil_sread_double (in, &(dat->z[i]));
            CCcheck_rval (rval, "CCutil_sread_double failed");
        }
        for (i = 0; i < ndepot; i++) {
            dat->x[*ncount + i] = 0.0;
            dat->y[*ncount + i] = 0.0;
            dat->z[*ncount + i] = 0.0;
        }
    } else if ((norm & CC_NORM_SIZE_BITS) == CC_MATRIX_NORM_SIZE){
        /* Matrix is the lower triangle plus the diagonal */

        dat->adj = CC_SAFE_MALLOC (*ncount + ndepot, int *);
        CCcheck_NULL (dat->adj, "out of memory in CCutil_getmaster");
        dat->adjspace =
            CC_SAFE_MALLOC ((*ncount+ndepot) * (*ncount+ndepot+1) / 2, int);
        CCcheck_NULL (dat->adjspace, "out of memory in CCutil_getmaster");

        for (i = 0, j = 0; i < *ncount; i++) {
            dat->adj[i] = dat->adjspace + j;
            j += (i+1);
        }
        for (i = 0; i < *ncount; i++) {
            for (j = 0; j <= i; j++) {
                rval = CCutil_sread_int (in, &(dat->adj[i][j]));
                CCcheck_rval (rval, "CCutil_sread_int failed");
            }
        }

        /* The initialization of the depot entries is not really needed;   */
        /* we add it in for now to help in debugging the subdivision code. */

        for (i = *ncount; i < *ncount + ndepot; i++) {
            for (j = 0; j < i; j++) {
                dat->adj[i][j] = dat->depotcost[j];
            }
            dat->adj[i][i] = 0;
       }
    } else if (norm == CC_DSJRANDNORM) {
        if (ndepot > 0) {
            fprintf (stderr, "cannot have depots with DSJRANDOM norm\n");
            rval = 1; goto CLEANUP;
        }
        dat->x = CC_SAFE_MALLOC (*ncount, double);
        CCcheck_NULL (dat->x, "out of memory in CCutil_getmaster");
        for (i = 0; i < *ncount; i++) {
            rval = CCutil_sread_double (in, &(dat->x[i]));
            CCcheck_rval (rval, "CCutil_sread_double failed");
        }
    } else if (norm == CC_SPARSE) {
        if (ndepot > 0) {
            fprintf (stderr, "cannot have depots with SPARSE norm\n");
            rval = 1; goto CLEANUP;
        }
        rval = readdat_sparse (in, *ncount, dat);
        CCcheck_rval (rval, "readmaster_sparse failed");
    } else if (norm == CC_USER) {
        if (ndepot > 0) {
            fprintf (stderr, "cannot have depots with USER norm\n");
            rval = 1; goto CLEANUP;
        }
        rval = readmaster_user (in, *ncount, &dat->userdat);
        CCcheck_rval (rval, "readmaster_user failed");
    } else if (norm == CC_RHMAP1 || norm == CC_RHMAP2 ||
               norm == CC_RHMAP3 || norm == CC_RHMAP4 ||
               norm == CC_RHMAP5) {
        if (ndepot > 0) {
            fprintf (stderr, "cannot have depots with RH norm\n");
            rval = 1; goto CLEANUP;
        }
        rval = readmaster_rhdata (in, *ncount, &dat->rhdat);
        CCcheck_rval (rval, "readmaster_rhdata failed");
    } else {
        fprintf (stderr, "unknown norm: %d\n", norm);
        rval = 1; goto CLEANUP;
    }

    *perm = CC_SAFE_MALLOC (*ncount + ndepot, int);
    CCcheck_NULL (*perm, "out of memory in CCutil_getmaster");

    for (i = 0; i < *ncount; i++) {
        rval = CCutil_sread_int (in, &((*perm)[i]));
        CCcheck_rval (rval, "CCutil_sread_int failed");
    }
    for (i = *ncount; i < *ncount + ndepot; i++) {
        (*perm)[i] = i;
    }

    if (subdivision == 1) {
        /* Read the original names (not permuted) */
        dat->orig_names = CC_SAFE_MALLOC (*ncount, int);
        CCcheck_NULL (dat->orig_names,
                     "out of memory in CCutil_getmaster");
        for (i = 0; i < *ncount; i++) {
            rval = CCutil_sread_int (in, &(dat->orig_names[i]));
            CCcheck_rval (rval, "CCutil_sread_int failed");
        }
    }

    dat->orig_ncount = *ncount;
    *ncount = dat->orig_ncount + ndepot;

CLEANUP:

    if (rval) {
        CC_IFFREE (*perm, int);
        CCutil_freedatagroup (dat);
    }

    return rval;
}

int CCutil_getnodeweights (char *weightname, int ncount, int weight_limit,
        double **wcoord, CCrandstate *rstate)
{
    int i, k;

    *wcoord = CC_SAFE_MALLOC (ncount, double);
    if (!(*wcoord))
        return 1;
    if (weightname != (char *) NULL) {
        FILE *weightin = fopen (weightname, "r");
        if (weightin == (FILE *) NULL) {
            perror (weightname);
            fprintf (stderr, "Unable to open %s for input\n", weightname);
            CC_FREE (*wcoord, double);
            return 1;
        }
        fscanf (weightin, "%d", &k);
        if (k != ncount) {
            fprintf (stderr, "Weight file does not match node file\n");
            fclose (weightin);
            CC_FREE (*wcoord, double);
            return 1;
        }
        for (i = 0; i < k; i++) {
            fscanf (weightin, "%lf", &((*wcoord)[i]));
        }
        make_weights_nonnegative (ncount, *wcoord);
        fclose (weightin);
    } else {
        for (i = 0; i < ncount; i++)
            (*wcoord)[i] = (double) (CCutil_lprand (rstate) % (weight_limit));
    }

    return 0;
}

static void make_weights_nonnegative (int ncount, double *wcoord)
{
    int i;
    double minx;

    minx = wcoord[0];

    for (i = 1; i < ncount; i++) {
        if (wcoord[i] < minx)
            minx = wcoord[i];
    }
    if (minx < 0.0) {
        printf ("****WARNING**** Adjusting node weights by %f\n", minx);
        for (i = 0; i < ncount; i++)
            wcoord[i] -= minx;
    }
}

#define MATRIX_LOWER_DIAG_ROW  0
#define MATRIX_UPPER_ROW       1
#define MATRIX_UPPER_DIAG_ROW  2
#define MATRIX_FULL_MATRIX     3

int CCutil_gettsplib(char *datname, int *ncount, CCdatagroup *dat)
{
    char buf[256], key[256], field[256];
    char *p;
    FILE *in;
    int matrixform = MATRIX_LOWER_DIAG_ROW;
    int norm = -1;

    CCutil_init_datagroup (dat);
    *ncount = -1;

    if ((in = fopen (datname, "r")) == (FILE *) NULL) {
        perror (datname);
        fprintf (stderr, "Unable to open %s for input\n", datname);
        return 1;
    }

    while (fgets (buf, 254, in) != (char *) NULL) {
        p = buf;
        while (*p != '\0') {
            if (*p == ':')
                *p = ' ';
            p++;
        }
        p = buf;
        if (sscanf (p, "%s", key) != EOF) {
            p += strlen (key);
            while (*p == ' ')
                p++;
            if (!strcmp (key, "NAME")) {
                printf ("Problem Name: %s", p);
            } else if (!strcmp (key, "TYPE")) {
                printf ("Problem Type: %s", p);
                if (sscanf (p, "%s", field) == EOF || strcmp (field, "TSP")) {
                    fprintf (stderr, "Not a TSP problem\n");
                    return 1;
                }
            } else if (!strcmp (key, "COMMENT")) {
                printf ("%s", p);
            } else if (!strcmp (key, "DIMENSION")) {
                if (sscanf (p, "%s", field) == EOF) {
                    fprintf (stderr, "ERROR in DIMENSION line\n");
                    return 1;
                }
                *ncount = atoi (field);
                printf ("Number of Nodes: %d\n", *ncount);
            } else if (!strcmp (key, "EDGE_WEIGHT_TYPE")) {
                if (sscanf (p, "%s", field) == EOF) {
                    fprintf (stderr, "ERROR in EDGE_WEIGHT_TYPE line\n");
                    return 1;
                }
                if (!strcmp (field, "EXPLICIT")) {
                    norm = CC_MATRIXNORM;
                    printf ("Explicit Lengths (CC_MATRIXNORM)\n");
                } else if (!strcmp (field, "EUC_2D")) {
                    norm = CC_EUCLIDEAN;
                    printf ("Rounded Euclidean Norm (CC_EUCLIDEAN)\n");
                } else if (!strcmp (field, "EUC_3D")) {
                    norm = CC_EUCLIDEAN_3D;
                    printf ("Rounded Euclidean 3D Norm (CC_EUCLIDEAN_3D)\n");
                } else if (!strcmp (field, "MAX_2D")) {
                    norm = CC_MAXNORM;
                    printf ("Max Norm (CC_MAXNORM)\n");
                } else if (!strcmp (field, "MAN_2D")) {
                    norm = CC_MANNORM;
                    printf ("Max Norm (CC_MAXNORM)\n");
                } else if (!strcmp (field, "GEO")) {
                    norm = CC_GEOGRAPHIC;
                    printf ("Geographical Norm (CC_GEOGRAPHIC)\n");
                } else if (!strcmp (field, "GEOM")) {
                    norm = CC_GEOM;
                    printf ("Geographical Norm in Meters (CC_GEOM)\n");
                } else if (!strcmp (field, "ATT")) {
                    norm = CC_ATT;
                    printf ("ATT Norm (CC_ATT)\n");
                } else if (!strcmp (field, "CEIL_2D")) {
                    norm = CC_EUCLIDEAN_CEIL;
                    printf ("Rounded Up Euclidean Norm (CC_EUCLIDEAN_CEIL)\n");
                } else if (!strcmp (field, "DSJRAND")) {
                    norm = CC_DSJRANDNORM;
                    printf ("David Johnson Random Norm (CC_DSJRANDNORM)\n");
                } else {
                    fprintf (stderr, "ERROR: Not set up for norm %s\n", field);
                    return 1;
                }
                if (CCutil_dat_setnorm (dat, norm)) {
                    fprintf (stderr, "ERROR: Couldn't set norm %d\n", norm);
                    return 1;
                }
            } else if (!strcmp (key, "EDGE_WEIGHT_FORMAT")) {
                if (sscanf (p, "%s", field) == EOF) {
                    fprintf (stderr, "ERROR in EDGE_WEIGHT_FORMAT line\n");
                    return 1;
                }
                if (!strcmp (field, "LOWER_DIAG_ROW")) {
                    matrixform = MATRIX_LOWER_DIAG_ROW;
                } else if (!strcmp (field, "UPPER_ROW")) {
                    matrixform = MATRIX_UPPER_ROW;
                } else if (!strcmp (field, "UPPER_DIAG_ROW")) {
                    matrixform = MATRIX_UPPER_DIAG_ROW;
                } else if (!strcmp (field, "FULL_MATRIX")) {
                    matrixform = MATRIX_FULL_MATRIX;
                } else if (strcmp (field, "FUNCTION")) {
                    fprintf (stderr, "Cannot handle format: %s\n", field);
                    return 1;
                }
            } else if (!strcmp (key, "NODE_COORD_SECTION")) {
                int i;
                if (*ncount <= 0) {
                    fprintf (stderr, "ERROR: Dimension not specified\n");
                    return 1;
                }
                if (dat->x != (double *) NULL) {
                    fprintf (stderr, "ERROR: A second NODE_COORD_SECTION?\n");
                    CCutil_freedatagroup (dat);
                    return 1;
                }
                if ((norm & CC_NORM_SIZE_BITS) == CC_D2_NORM_SIZE) {
                    dat->x = CC_SAFE_MALLOC (*ncount, double);
                    if (!dat->x) {
                        CCutil_freedatagroup (dat);
                        return 1;
                    }
                    dat->y = CC_SAFE_MALLOC (*ncount, double);
                    if (!dat->y) {
                        CCutil_freedatagroup (dat);
                        return 1;
                    }
                    for (i = 0; i < *ncount; i++) {
                        fscanf (in, "%*d %lf %lf", &(dat->x[i]), &(dat->y[i]));
                    }
                } else if ((norm & CC_NORM_SIZE_BITS) == CC_D3_NORM_SIZE) {
                    dat->x = CC_SAFE_MALLOC (*ncount, double);
                    if (!dat->x) {
                        CCutil_freedatagroup (dat);
                        return 1;
                    }
                    dat->y = CC_SAFE_MALLOC (*ncount, double);
                    if (!dat->y) {
                        CCutil_freedatagroup (dat);
                        return 1;
                    }
                    dat->z = CC_SAFE_MALLOC (*ncount, double);
                    if (!dat->z) {
                        CCutil_freedatagroup (dat);
                        return 1;
                    }
                    for (i = 0; i < *ncount; i++) {
                        fscanf (in, "%*d %lf %lf %lf",
                               &(dat->x[i]), &(dat->y[i]), &(dat->z[i]));
                    }
                } else {
                    fprintf (stderr, "ERROR: Node coordinates with norm %d?\n",
                                 norm);
                    return 1;
                }
            } else if (!strcmp (key, "EDGE_WEIGHT_SECTION")) {
                int i, j;
                if (*ncount <= 0) {
                    fprintf (stderr, "ERROR: Dimension not specified\n");
                    return 1;
                }
                if (dat->adj != (int **) NULL) {
                    fprintf (stderr, "ERROR: A second NODE_COORD_SECTION?\n");
                    CCutil_freedatagroup (dat);
                    return 1;
                }
                if ((norm & CC_NORM_SIZE_BITS) == CC_MATRIX_NORM_SIZE) {
                    dat->adj = CC_SAFE_MALLOC (*ncount, int *);
                    dat->adjspace = CC_SAFE_MALLOC ((*ncount)*(*ncount+1)/2,
                                                    int);
                    if (dat->adj == (int **) NULL ||
                        dat->adjspace == (int *) NULL) {
                        CCutil_freedatagroup (dat);
                        return 1;
                    }
                    for (i = 0, j = 0; i < *ncount; i++) {
                        dat->adj[i] = dat->adjspace + j;
                        j += (i+1);
                    }
                    if (matrixform == MATRIX_LOWER_DIAG_ROW) {
                        for (i = 0; i < *ncount; i++) {
                            for (j = 0; j <= i; j++)
                                fscanf (in, "%d", &(dat->adj[i][j]));
                        }
                    } else if (matrixform == MATRIX_UPPER_ROW ||
                               matrixform == MATRIX_UPPER_DIAG_ROW ||
                               matrixform == MATRIX_FULL_MATRIX) {
                        int **tempadj = (int **) NULL;
                        int *tempadjspace = (int *) NULL;
                        tempadj = CC_SAFE_MALLOC (*ncount, int *);
                        tempadjspace = CC_SAFE_MALLOC ((*ncount) * (*ncount),
                                                       int);
                        if (tempadj == (int **) NULL ||
                            tempadjspace == (int *) NULL) {
                            CC_IFFREE (tempadj, int *);
                            CC_IFFREE (tempadjspace, int);
                            CCutil_freedatagroup (dat);
                            return 1;
                        }
                        for (i = 0; i < *ncount; i++) {
                            tempadj[i] = tempadjspace + i * (*ncount);
                            if (matrixform == MATRIX_UPPER_ROW) {
                                tempadj[i][i] = 0;
                                for (j = i + 1; j < *ncount; j++)
                                    fscanf (in, "%d", &(tempadj[i][j]));
                            } else if (matrixform == MATRIX_UPPER_DIAG_ROW) {
                                for (j = i; j < *ncount; j++)
                                    fscanf (in, "%d", &(tempadj[i][j]));
                            } else {
                                for (j = 0; j < *ncount; j++)
                                    fscanf (in, "%d", &(tempadj[i][j]));
                            }
                        }
                        for (i = 0; i < *ncount; i++) {
                            for (j = 0; j <= i; j++)
                                dat->adj[i][j] = tempadj[j][i];
                        }
                        CC_FREE (tempadjspace, int);
                        CC_FREE (tempadj, int *);
                    }
                } else {
                    fprintf (stderr, "ERROR: Matrix with norm %d?\n",
                             norm);
                    return 1;
                }
            } else if (!strcmp (key, "FIXED_EDGES_SECTION")) {
                fprintf (stderr, "ERROR: Not set up for fixed edges\n");
                return 1;
            }
        }
    }
    fclose (in);

    if (dat->x == (double *) NULL && dat->adj == (int **) NULL) {
        fprintf (stderr, "ERROR: Didn't find the data\n");
        return 1;
    } else {
        return 0;
    }
}

int CCutil_writetsplib (const char *fname, int ncount, CCdatagroup *dat)
{
    int rval = 0;
    int i, j, k;
    FILE *out = (FILE *) NULL;

    if ((out = fopen (fname, "w")) == (FILE *) NULL) {
        perror (fname);
        fprintf (stderr, "Unable to open %s for input\n", fname);
        rval = 1; goto CLEANUP;
    }

    fprintf (out, "NAME: concorde%d\n", ncount);
    fprintf (out, "TYPE: TSP\n");
    fprintf (out, "COMMENT: Generated by CCutil_writetsplib\n");
    fprintf (out, "DIMENSION: %d\n", ncount);

    fprintf (out, "EDGE_WEIGHT_TYPE: ");
    switch (dat->norm) {
    case CC_MAXNORM:
        fprintf (out, "MAX_2D\n");
        break;
    case CC_MANNORM:
        fprintf (out, "MAN_2D\n");
        break;
    case CC_EUCLIDEAN_CEIL:
        fprintf (out, "CEIL_2D\n");
        break;
    case CC_EUCLIDEAN:
        fprintf (out, "EUC_2D\n");
        break;
    case CC_EUCLIDEAN_3D:
        fprintf (out, "EUC_3D\n");
        break;
    case CC_USER:
        fprintf (out, "USER\n");
        fprintf (stderr, "Warning: Norm not a supported in TSPLIB\n");
        rval = 1; goto CLEANUP;
    case CC_ATT:
        fprintf (out, "ATT\n");
        break;
    case CC_GEOGRAPHIC:
        fprintf (out, "GEO\n");
        break;
    case CC_GEOM:
        fprintf (out, "GEOM\n");
        break;
    case CC_MATRIXNORM:
        fprintf (out, "EXPLICIT\n");
        break;
    case CC_DSJRANDNORM:
        fprintf (out, "DSJ_RANDOM\n");
        fprintf (stderr, "Warning: Norm not a supported in TSPLIB\n");
        rval = 1; goto CLEANUP;
    case CC_CRYSTAL:
        fprintf (out, "CRYSTAL\n");
        fprintf (stderr, "Warning: Crystal is not supported in TSPLIB\n");
        rval = 1; goto CLEANUP;
    case CC_RHMAP1:
    case CC_RHMAP2:
    case CC_RHMAP3:
    case CC_RHMAP4:
    case CC_RHMAP5:
        fprintf (out, "RHMAPx\n");
        fprintf (stderr, "Warning: Norm not a supported in TSPLIB\n");
        rval = 1; goto CLEANUP;
    default:
        fprintf (stderr, "unknown NORM\n");
        rval = 1; goto CLEANUP;
    }

    switch (dat->norm & CC_NORM_SIZE_BITS) {
    case CC_D2_NORM_SIZE:
        fprintf (out, "NODE_COORD_SECTION\n");
        for (i = 0; i < ncount; i++) {
            fprintf (out, "%d %f %f\n", i+1, dat->x[i], dat->y[i]);
        }
        break;
    case CC_D3_NORM_SIZE:
        fprintf (out, "NODE_COORD_SECTION\n");
        for (i = 0; i < ncount; i++) {
            fprintf (out, "%d %f %f %f\n", i+1, dat->x[i], dat->y[i], dat->z[i]);
        }
        break;
    case CC_MATRIX_NORM_SIZE:
        fprintf (out, "EDGE_WEIGHT_FORMAT: LOWER_DIAG_ROW\n");
        fprintf (out, "EDGE_WEIGHT_SECTION\n");
        for (i = 0, k = 0; i < ncount; i++) {
            for (j = 0; j <= i; j++) {
                fprintf (out, "%2d ", dat->adj[i][j]);
/*
                if ((k++ % 5) == 4) {
                    fprintf (out, "\n");
                }
*/
            }
            fprintf (out, "\n"); 
        }
        if ((k % 5) != 4) {
            fprintf (out, "\n");
        }
        break;
    default:
        fprintf (stderr, "unknown NORM_SIZE\n");
        rval = 1; goto CLEANUP;
    }

CLEANUP:
    
    if (out != (FILE *) NULL) fclose (out);
    return rval;
}

int CCutil_getcycle_tsplib (int ncount, char *cyclename, int *outcycle)
{
    int rval = 0;
    FILE *in = (FILE *) NULL;
    char buf[256], key[256], field[256];
    char *p;
    int icount = 0;

    in = fopen (cyclename, "r");
    if (in == (FILE *) NULL) {
        perror (cyclename);
        fprintf (stderr, "Unable to open %s for input\n", cyclename);
        rval = 1; goto CLEANUP;
    }

    while (fgets (buf, 254, in) != (char *) NULL) {
        p = buf;
        while (*p != '\0') {
            if (*p == ':')
                *p = ' ';
            p++;
        }
        p = buf;
        if (sscanf (p, "%s", key) != EOF) {
            p += strlen (key);
            while (*p == ' ')
                p++;
            if (!strcmp (key, "NAME")) {
                printf ("Name: %s", p);
            } else if (!strcmp (key, "TYPE")) {
                printf ("Problem Type: %s", p);
                if (sscanf (p, "%s", field) == EOF || strcmp (field, "TOUR")) {
                    fprintf (stderr, "Not a TOUR File\n");
                    rval = 1; goto CLEANUP;
                }
            } else if (!strcmp (key, "COMMENT")) {
                printf ("%s", p);
            } else if (!strcmp (key, "DIMENSION")) {
                if (sscanf (p, "%s", field) == EOF) {
                    fprintf (stderr, "ERROR in DIMENSION line\n");
                    rval = 1; goto CLEANUP;
                }
                icount = atoi (field);
                if (icount != ncount) {
                    fprintf (stderr, "Number of nodes does not agree\n");
                    rval = 1; goto CLEANUP;
                }
            } else if (!strcmp (key, "TOUR_SECTION")) {
                int i, k;
                if (icount <= 0) {
                    fprintf (stderr, "ERROR: Dimension not specified\n");
                    rval = 1; goto CLEANUP;
                }
                for (i = 0; i < icount; i++) {
                    fscanf (in, "%d", &k);
                    if (k < 1 && k > ncount) {
                        fprintf (stderr, "ERROR: Bad format in TSPLIB tour\n");
                        rval = 1; goto CLEANUP;
                    }
                    outcycle[i] = k - 1;
                }
                fscanf (in, "%d", &k);
                if (k != -1) {
                    fprintf (stderr, "Warning: tour not -1 terminated\n");
                }
            } else if (!strcmp (key, "EOF")) {
                goto CLEANUP;
            } else {
                fprintf (stderr, "ERROR %s: Bad section in TSPLIB tour\n", key);
                rval = 1; goto CLEANUP;
            }
        }
    }


CLEANUP:

    if (in) fclose (in);

    return rval;
}

int CCutil_datagroup_perm (int ncount, CCdatagroup *dat, int *perm)
{
    int i, j;

    if (dat->norm == CC_USER) {
        if (permute_user (ncount, &dat->userdat, perm)) {
            fprintf (stderr, "permute_user failed\n");
            return 1;
        }
        return 0;
    } else if (dat->norm == CC_RHMAP1 || dat->norm == CC_RHMAP2 ||
               dat->norm == CC_RHMAP3 || dat->norm == CC_RHMAP4 ||
               dat->norm == CC_RHMAP5) {
        if (permute_rhdata (ncount, &dat->rhdat, perm)) {
            fprintf (stderr, "permute_rhdata failed\n");
            return 1;
        }
    } else if (dat->norm == CC_SPARSE) {
        if (permute_sparse (ncount, dat, perm)) {
            fprintf (stderr, "permute_sparse failed\n");
            return 1;
        }
        return 0;
    }
    
    if (dat->x != (double *) NULL) {
        double *tempx;

        tempx = CC_SAFE_MALLOC (ncount, double);
        if (!tempx)
            return 1;
        for (i = 0; i < ncount; i++) {
            tempx[i] = dat->x[perm[i]];
        }
        CC_FREE (dat->x, double);
        dat->x = tempx;
    }
    if (dat->y != (double *) NULL) {
        double *tempy;

        tempy = CC_SAFE_MALLOC (ncount, double);
        if (!tempy)
            return 1;
        for (i = 0; i < ncount; i++) {
            tempy[i] = dat->y[perm[i]];
        }
        CC_FREE (dat->y, double);
        dat->y = tempy;
    }
    if (dat->z != (double *) NULL) {
        double *tempz;

        tempz = CC_SAFE_MALLOC (ncount, double);
        if (!tempz)
            return 1;
        for (i = 0; i < ncount; i++) {
            tempz[i] = dat->z[perm[i]];
        }
        CC_FREE (dat->z, double);
        dat->z = tempz;
    }
    if (dat->adj != (int **) NULL) {
        int **tempadj = (int **) NULL;
        int *tempadjspace = (int *) NULL;

        tempadj = CC_SAFE_MALLOC (ncount, int *);
        tempadjspace = CC_SAFE_MALLOC (ncount * (ncount+1) / 2, int);
        
        if (tempadj == (int **) NULL ||
            tempadjspace == (int *) NULL) {
            CC_IFFREE (tempadj, int *);
            CC_IFFREE (tempadjspace, int);
            return 1;
        }

        for (i = 0, j = 0; i < ncount; i++) {
            tempadj[i] = tempadjspace + j;
            j += (i+1);
        }
        for (i = 0, j = 0; i < ncount; i++) {
            for (j = 0; j <= i; j++) {
                if (perm[i] <  perm[j])
                    tempadj[i][j] = dat->adj[perm[j]][perm[i]];
                else
                    tempadj[i][j] = dat->adj[perm[i]][perm[j]];
            }
        }
        CC_FREE (dat->adj, int *);
        CC_FREE (dat->adjspace, int);
        dat->adj = tempadj;
        dat->adjspace = tempadjspace;
    }
    return 0;
}

int CCutil_copy_datagroup (int ncount, CCdatagroup *indat, CCdatagroup *outdat)
{
    int i, j;
    int rval = 0;

    CCutil_init_datagroup (outdat);
    CCutil_dat_setnorm (outdat, indat->norm);

    if (indat->norm == CC_USER) {
        fprintf (stderr, "CCutil_copy_datagroup not set up for user norm\n");
        return 1;
    }

    if (indat->norm == CC_RHMAP1 || indat->norm == CC_RHMAP2 ||
        indat->norm == CC_RHMAP3 || indat->norm == CC_RHMAP4 ||
        indat->norm == CC_RHMAP5) {
        fprintf (stderr, "CCutil_copy_datagroup not set up for rh vectors\n");
        return 1;
    }

    if (indat->norm == CC_SPARSE) {
        rval = copy_sparse (ncount, indat, outdat);
        if (rval) {
            fprintf (stderr, "copy_sparse failed\n"); goto CLEANUP;
        }
        goto CLEANUP;
    }

    if (indat->x != (double *) NULL) {
        double *tempx;

        tempx = CC_SAFE_MALLOC (ncount, double);
        if (!tempx) {
            fprintf (stderr, "out of memory in CCutil_copy_datagroup\n");
            rval = 1; goto CLEANUP;
        }
        for (i = 0; i < ncount; i++) {
            tempx[i] = indat->x[i];
        }
        outdat->x = tempx;
    }
    if (indat->y != (double *) NULL) {
        double *tempy;

        tempy = CC_SAFE_MALLOC (ncount, double);
        if (!tempy) {
            fprintf (stderr, "out of memory in CCutil_copy_datagroup\n");
            rval = 1; goto CLEANUP;
        }
        for (i = 0; i < ncount; i++) {
            tempy[i] = indat->y[i];
        }
        outdat->y = tempy;
    }
    if (indat->z != (double *) NULL) {
        double *tempz;

        tempz = CC_SAFE_MALLOC (ncount, double);
        if (!tempz) {
            fprintf (stderr, "out of memory in CCutil_copy_datagroup\n");
            rval = 1; goto CLEANUP;
        }
        for (i = 0; i < ncount; i++) {
            tempz[i] = indat->z[i];
        }
        outdat->z = tempz;
    }
    if (indat->adj != (int **) NULL) {
        int **tempadj     = (int **) NULL;
        int *tempadjspace = (int *) NULL;

        tempadj      = CC_SAFE_MALLOC (ncount, int *);
        tempadjspace = CC_SAFE_MALLOC (ncount * (ncount+1) / 2, int);
        
        if (!tempadj ||!tempadjspace) {
            fprintf (stderr, "out of memory in CCutil_copy_datagroup\n");
            CC_IFFREE (tempadj, int *);
            CC_IFFREE (tempadjspace, int);
            rval = 1; goto CLEANUP;
        }

        for (i = 0, j = 0; i < ncount; i++) {
            tempadj[i] = tempadjspace + j;
            j += (i+1);
        }
        for (i = 0, j = 0; i < ncount; i++) {
            for (j = 0; j <= i; j++) {
                    tempadj[i][j] = indat->adj[i][j];
            }
        }
        outdat->adj = tempadj;
        outdat->adjspace = tempadjspace;
    }

CLEANUP:

    if (rval) {
        CCutil_freedatagroup (outdat);
    }
    return rval;
}

int CCutil_getedgelist (int ncount, char *fname, int *ecount, int **elist,
        int **elen, int binary_in)
{
    int k;

    if (CCutil_getedgelist_n (&k, fname, ecount, elist, elen, binary_in)) {
        fprintf (stderr, "CCutil_getedgelist_n failed\n");
        return 1;
    }

    if (k != ncount) {
        fprintf (stderr, "Edge file does not match problem\n");
        return 1;
    }

    return 0;
}

int CCutil_getedgelist_n (int *ncount, char *fname, int *ecount, int **elist,
        int **elen, int binary_in)
{
    FILE *f_in = (FILE *) NULL;
    CC_SFILE *s_in = (CC_SFILE *) NULL;
    int i;
    int rval;

    *elist = (int *) NULL;
    *elen = (int *) NULL;

    if (binary_in) {
        if ((s_in = CCutil_sopen (fname, "r")) == (CC_SFILE *) NULL) {
            fprintf (stderr, "Unable to open %s for input\n", fname);
            rval = 1; goto CLEANUP;
        }
        if (CCutil_sread_int (s_in, ncount) ||
            CCutil_sread_int (s_in, ecount)) {
            fprintf (stderr, "CCutil_sread_int failed\n");
            rval = 1; goto CLEANUP;
        }
    } else {
        if ((f_in = fopen (fname, "r")) == (FILE *) NULL) {
            perror (fname);
            fprintf (stderr, "Unable to open %s for input\n", fname);
            rval = 1; goto CLEANUP;
        }
        *ncount = CCutil_readint (f_in);
        *ecount = CCutil_readint (f_in);
    }
    
    *elist = CC_SAFE_MALLOC(2 * (*ecount), int);
    if (!(*elist)) {
        fprintf (stderr, "out of memory in CCutil_getedgelist_binary_n\n");
        rval = 1; goto CLEANUP;
    }
    *elen = CC_SAFE_MALLOC(*ecount, int);
    if (!(*elen)) {
        fprintf (stderr, "out of memory in CCutil_getedgelist_binary_in\n");
        rval = 1; goto CLEANUP;
    }

    for (i = 0; i < *ecount; i++) {
        if (binary_in) {
            if (CCutil_sread_int (s_in, &((*elist)[2*i])) ||
                CCutil_sread_int (s_in, &((*elist)[2*i+1])) ||
                CCutil_sread_int (s_in, &((*elen)[i]))) {
                fprintf (stderr, "CCutil_sread_int failed\n");
                rval = 1; goto CLEANUP;
            }
        } else {
            (*elist)[2*i] = CCutil_readint (f_in);
            (*elist)[2*i+1] = CCutil_readint (f_in);
            (*elen)[i] = CCutil_readint (f_in);
        }
    }
    rval = 0;

CLEANUP:
    CCutil_sclose (s_in);
    if (f_in != (FILE *) NULL) fclose (f_in);
    if (rval) {
        CC_IFFREE (*elist, int);
        CC_IFFREE (*elen, int);
    }
    return rval;
}

int CCutil_genedgelist (int ncount, int ecount, int **elist, int **elen,
        CCdatagroup *dat, int maxlen, CCrandstate *rstate)
{
    int i, head, tail, temp;
    int rval = 0;
    CCutil_edgehash eh;
    int have_eh = 0;
    int val = 0;

    *elist = (int *) NULL;
    *elen = (int *) NULL;

    if (ecount > (ncount * (ncount - 1)) / 2) {
        fprintf (stderr, "Cannot generate %d edges in a %d node graph\n",
                 ecount, ncount);
        rval = 1; goto CLEANUP;
    }

    *elist = CC_SAFE_MALLOC(2 * (ecount), int);
    *elen  = CC_SAFE_MALLOC(ecount, int);
    if (!(*elist) || !(*elen)) {
        fprintf (stderr, "out of memory in CCutil_genedgelist\n");
        rval = 1; goto CLEANUP;
    }

    if (ecount > (ncount * (ncount - 1)) / 20) {
        int eleft = (ncount * (ncount - 1)) / 2;
        int j;
        for (i=0; i<ncount; i++) {
            for (j=i+1; j<ncount; j++) {
                if (CCutil_lprand (rstate) % eleft < ecount) {
                    ecount--;
                    (*elist)[2*ecount] = i;
                    (*elist)[2*ecount+1] = j;
                    if (dat) {
                        (*elen)[ecount] = CCutil_dat_edgelen (i, j, dat);
                    } else {
                        (*elen)[ecount] = CCutil_lprand (rstate) % maxlen;
                    }
                }
                eleft--;
            }
        }
    } else {
        rval = CCutil_edgehash_init (&eh, (int) (ecount * 1.5));
        if (rval) {
            fprintf (stderr, "CCutil_edgehash_init failed\n");
            goto CLEANUP;
        }
        have_eh = 1;
        
        for (i = 0; i < ecount; i++) {
            do {
                head = CCutil_lprand (rstate) % ncount;
                tail = CCutil_lprand (rstate) % ncount;
                if (head > tail) {
                    CC_SWAP (head, tail, temp); 
                }
            } while ((head == tail) || CCutil_edgehash_find (&eh, head, tail, &val) == 0);
            rval = CCutil_edgehash_add (&eh, head, tail, 1);
            if (rval) {
                fprintf (stderr, "CCutil_edgehash_add failed\n");
                goto CLEANUP;
            }
            (*elist)[2*i]   = head;
            (*elist)[2*i+1] = tail;
            if (dat) {
                (*elen)[i] = CCutil_dat_edgelen (head, tail, dat);
            } else {
                (*elen)[i] = CCutil_lprand (rstate) % maxlen;
            }
        }
    }

CLEANUP:

    if (rval) {
        CC_IFFREE (*elist, int);
        CC_IFFREE (*elen, int);
    }
    if (have_eh) {
        CCutil_edgehash_free (&eh);
    }
    return rval;
}

int CCutil_getcycle_edgelist (int ncount, char *cyclename, int *outcycle,
        int binary_in)
{
    FILE *cycfin = (FILE *) NULL;
    CC_SFILE *cycsin = (CC_SFILE *) NULL;
    int *elist = (int *) NULL;
    int i, k, istour, rval = 0;

    if (binary_in) {
        cycsin = CCutil_sopen (cyclename, "r");
        if (cycsin == (CC_SFILE *) NULL) {
            fprintf (stderr, "Unable to open %s for input\n", cyclename);
            rval = 1; goto CLEANUP;
        }
    } else {
        cycfin = fopen (cyclename, "r");
        if (cycfin == (FILE *) NULL) {
            perror (cyclename);
            fprintf (stderr, "Unable to open %s for input\n", cyclename);
            rval = 1; goto CLEANUP;
        }
    }

    elist = CC_SAFE_MALLOC (2 * ncount, int);
    if (!elist) {
        fprintf (stderr, "out of memory in CCutil_getcycle_edgelist\n");
        rval = 1; goto CLEANUP;
    }

    if (binary_in) {
        if (CCutil_sread_int (cycsin, &i) ||
            CCutil_sread_int (cycsin, &k)) {
            fprintf (stderr, "CCutil_sread_int failed\n");
            rval = 1; goto CLEANUP;
        }
    } else {
        fscanf (cycfin, "%d %d", &i, &k);
    }
    if (i != ncount || k != ncount) {
        fprintf (stderr, "file is not a cycle-edge file for this problem\n");
        rval = 1; goto CLEANUP;
    }

    for (i = 0; i < ncount; i++) {
        if (binary_in) {
            if (CCutil_sread_int (cycsin, &(elist[2*i])) ||
                CCutil_sread_int (cycsin, &(elist[2*i+1]))) {
                fprintf (stderr, "CCutil_sread_int failed\n");
                rval = 1; goto CLEANUP;
            }
        } else {
            fscanf (cycfin, "%d %d %*d", &(elist[2 * i]),
                    &(elist[(2 * i) + 1]));
        }
    }

    rval = CCutil_edge_to_cycle (ncount, elist, &istour, outcycle);
    if (rval) {
        fprintf (stderr, "CCutil_edge_to_cycle failed\n"); goto CLEANUP;
    }
    if (istour == 0) {
        fprintf (stderr, "Edge-file is not a tour\n");
        rval = 1; goto CLEANUP;
    }
    rval = 0;

CLEANUP:
    CC_IFFREE (elist, int);
    if (cycfin != (FILE *) NULL) fclose (cycfin);
    CCutil_sclose (cycsin);
    return rval;
}

int CCutil_getcycle (int ncount, char *cyclename, int *outcycle, int binary_in)
{
    CC_SFILE *cycsin = (CC_SFILE *) NULL;
    FILE *cycfin = (FILE *) NULL;
    int i;
    int rval;

    if (binary_in) {
        cycsin = CCutil_sopen (cyclename, "r");
        if (cycsin == (CC_SFILE *) NULL) {
            fprintf (stderr, "Unable to open %s for input\n", cyclename);
            rval = 1; goto CLEANUP;
        }
        if (CCutil_sread_int (cycsin, &i)) {
            fprintf (stderr, "CCutil_sread_int failed\n");
            rval = 1; goto CLEANUP;
        }
    } else {
        cycfin = fopen (cyclename, "r");
        if (cycfin == (FILE *) NULL) {
            perror (cyclename);
            fprintf (stderr, "Unable to open %s for input\n", cyclename);
            rval = 1; goto CLEANUP;
        }
        i = CCutil_readint (cycfin);
    }

    if (i != ncount) {
        fprintf (stderr, "Cycle files has wrong number of nodes\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < ncount; i++) {
        if (binary_in) {
            if (CCutil_sread_int (cycsin, &outcycle[i])) {
                fprintf (stderr, "CCutil_sread_int failed\n");
                rval = 1; goto CLEANUP;
            }
        } else {
            outcycle[i] = CCutil_readint (cycfin);
        }
    }
    rval = 0;

CLEANUP:
    CCutil_sclose (cycsin);
    if (cycfin != (FILE *) NULL) fclose (cycfin);
    return rval;
}

void CCutil_cycle_len (int ncount, CCdatagroup *dat, int *cycle, double *len)
{
    int i;

    *len = 0.0;
    for (i = 1; i < ncount; i++) {
        (*len) += (double) CCutil_dat_edgelen (cycle[i-1], cycle[i], dat);
    } 
    (*len) += (double) CCutil_dat_edgelen (cycle[ncount-1], cycle[0], dat);
}

int CCutil_getedges_double (int *ncount, char *fname, int *ecount, int **elist,
                double **elen, int binary_in)
{
    int i;
    CC_SFILE *s_in = (CC_SFILE *) NULL;
    FILE *f_in = (FILE *) NULL;
    int rval;

    *elist = (int *) NULL;
    *elen = (double *) NULL;

    if (binary_in) {
        s_in = CCutil_sopen (fname, "r");
        if (s_in == (CC_SFILE *) NULL) {
            fprintf (stderr, "Unable to open %s for input\n", fname);
            rval = 1; goto CLEANUP;
        }
        if (CCutil_sread_int (s_in, ncount) ||
            CCutil_sread_int (s_in, ecount)) {
            fprintf (stderr, "CCutil_sread_int failed\n");
            rval = 1; goto CLEANUP;
        }
    } else {
        f_in = fopen (fname, "r");
        if (f_in == (FILE *) NULL) {
            perror (fname);
            fprintf (stderr, "Unable to open %s for input\n", fname);
            return 1;
        }

        *ncount = CCutil_readint (f_in);
        *ecount = CCutil_readint (f_in);
    }

    *elist = CC_SAFE_MALLOC(2 * (*ecount), int);
    if (!(*elist)) {
        rval = 1; goto CLEANUP;
    }
    *elen = CC_SAFE_MALLOC(*ecount, double);
    if (!(*elen)) {
        rval = 1; goto CLEANUP;
    }

    for (i = 0; i < *ecount; i++) {
        if (binary_in) {
            if (CCutil_sread_int (s_in, &((*elist)[2*i])) ||
                CCutil_sread_int (s_in, &((*elist)[2*i+1])) ||
                CCutil_sread_double (s_in, &((*elen)[i]))) {
                fprintf (stderr, "CCutil_sread_int/double failed\n");
                rval = 1; goto CLEANUP;
            }
        } else {
            (*elist)[2*i] = CCutil_readint (f_in);
            (*elist)[2*i+1] = CCutil_readint (f_in);
            if (fscanf (f_in, "%lf", &((*elen)[i])) != 1) {
                fprintf (stderr, "input file is in the wrong format\n");
                rval = 1; goto CLEANUP;
            }
        }
    }
    rval = 0;
CLEANUP:
    CCutil_sclose (s_in);
    if (f_in != (FILE *) NULL) fclose (f_in);
    if (rval) {
        CC_IFFREE ((*elist), int);
        CC_IFFREE ((*elen), double);
    }
    return rval;
}

int CCutil_writeedges (int ncount, char *outedgename, int ecount, int *elist,
        CCdatagroup *dat, int binary_out)
{
    FILE *fout = (FILE *) NULL;
    CC_SFILE *sout = (CC_SFILE *) NULL;
    int i;
    int l;
    int rval;

    if (binary_out) {
        sout = CCutil_sopen (outedgename, "w");
        if (sout == (CC_SFILE *) NULL) {
            fprintf (stderr, "Unable to open %s for output\n", outedgename);
            rval = 1; goto CLEANUP;
        }
        if (CCutil_swrite_int (sout, ncount)) {
            fprintf (stderr, "CCutil_swrite_int failed\n");
            rval = 1; goto CLEANUP;
        }
        if (CCutil_swrite_int (sout, ecount)) {
            fprintf (stderr, "CCutil_swrite_int failed\n");
            rval = 1; goto CLEANUP;
        }
        for (i = 0; i < ecount; i++) {
            if (CCutil_swrite_int (sout, elist[2*i])) {
                fprintf (stderr, "CCutil_swrite_int failed\n");
                rval = 1; goto CLEANUP;
            }
            if (CCutil_swrite_int (sout, elist[2*i+1])) {
                fprintf (stderr, "CCutil_swrite_int failed\n");
                rval = 1; goto CLEANUP;
            }
            l = CCutil_dat_edgelen (elist[2 * i], elist[(2 * i) + 1], dat);
            if (CCutil_swrite_int (sout, l)) {
                fprintf (stderr, "CCutil_swrite_int failed\n");
                rval = 1; goto CLEANUP;
            }
            
        }
    } else {
        fout = fopen (outedgename, "w");
        if (fout == (FILE *) NULL) {
            perror (outedgename);
            fprintf (stderr, "Unable to open %s for output\n", outedgename);
            rval = 1; goto CLEANUP;
        }
        fprintf (fout, "%d %d\n", ncount, ecount);
        for (i = 0; i < ecount; i++) {
            fprintf (fout, "%d %d %d\n", elist[2 * i], elist[(2 * i) + 1],
                   CCutil_dat_edgelen (elist[2 * i], elist[(2 * i) + 1], dat));
        }
    }
    rval = 0;

CLEANUP:
    CCutil_sclose (sout);
    if (fout != (FILE *) NULL) fclose (fout);

    return rval;
}

int CCutil_writeedges_int (int ncount, char *outedgename, int ecount,
        int *elist, int *elen, int binary_out)
{
    FILE *fout = (FILE *) NULL;
    CC_SFILE *sout = (CC_SFILE *) NULL;
    int i;
    int rval;

    if (binary_out) {
        sout = CCutil_sopen (outedgename, "w");
        if (sout == (CC_SFILE *) NULL) {
            fprintf (stderr, "Unable to open %s for output\n", outedgename);
            rval = 1; goto CLEANUP;
        }

        if (CCutil_swrite_int (sout, ncount)) {
            fprintf (stderr, "CCutil_swrite_int failed\n");
            rval = 1; goto CLEANUP;
        }
        if (CCutil_swrite_int (sout, ecount)) {
            fprintf (stderr, "CCutil_swrite_int failed\n");
            rval = 1; goto CLEANUP;
        }
        for (i = 0; i < ecount; i++) {
            if (CCutil_swrite_int (sout, elist[2*i])) {
                fprintf (stderr, "CCutil_swrite_int failed\n");
                rval = 1; goto CLEANUP;
            }
            if (CCutil_swrite_int (sout, elist[2*i+1])) {
                fprintf (stderr, "CCutil_swrite_int failed\n");
                rval = 1; goto CLEANUP;
            }
            if (CCutil_swrite_int (sout, elen[i])) {
                fprintf (stderr, "CCutil_swrite_int failed\n");
                rval = 1; goto CLEANUP;
            }
        }
    } else {
        fout = fopen (outedgename, "w");
        if (fout == (FILE *) NULL) {
            perror (outedgename);
            fprintf (stderr, "Unable to open %s for output\n", outedgename);
            rval = 1; goto CLEANUP;
        }
        fprintf (fout, "%d %d\n", ncount, ecount);
        for (i = 0; i < ecount; i++) {
            fprintf (fout, "%d %d %d\n", elist[2 * i], elist[(2 * i) + 1],
                     elen[i]);
        }
    }
    rval = 0;

CLEANUP:
    CCutil_sclose (sout);
    if (fout != (FILE *) NULL) fclose (fout);

    return rval;
}

int CCutil_writecycle_edgelist (int ncount, char *outedgename, int *cycle,
        CCdatagroup *dat, int binary_out)
{
    FILE *fout = (FILE *) NULL;
    CC_SFILE *sout = (CC_SFILE *) NULL;
    int i;
    int l;
    int rval;

    if (binary_out) {
        sout = CCutil_sopen (outedgename, "w");
        if (sout == (CC_SFILE *) NULL) {
            fprintf (stderr, "Unable to open %s for output\n", outedgename);
            rval = 1; goto CLEANUP;
        }
        if (CCutil_swrite_int (sout, ncount)) {
            fprintf (stderr, "CCutil_swrite_int failed\n");
            rval = 1; goto CLEANUP;
        }
        if (CCutil_swrite_int (sout, ncount)) {
            fprintf (stderr, "CCutil_swrite_int failed\n");
            rval = 1; goto CLEANUP;
        }
        for (i = 1; i < ncount; i++) {
            if (CCutil_swrite_int (sout, cycle[i-1])) {
                fprintf (stderr, "CCutil_swrite_int failed\n");
                rval = 1; goto CLEANUP;
            }
            if (CCutil_swrite_int (sout, cycle[i])) {
                fprintf (stderr, "CCutil_swrite_int failed\n");
                rval = 1; goto CLEANUP;
            }
            l = CCutil_dat_edgelen (cycle[i - 1], cycle[i], dat);
            if (CCutil_swrite_int (sout, l)) {
                fprintf (stderr, "CCutil_swrite_int failed\n");
                rval = 1; goto CLEANUP;
            }
        }
        if (CCutil_swrite_int (sout, cycle[ncount - 1])) {
            fprintf (stderr, "CCutil_swrite_int failed\n");
            rval = 1; goto CLEANUP;
        }
        if (CCutil_swrite_int (sout, cycle[0])) {
            fprintf (stderr, "CCutil_swrite_int failed\n");
            rval = 1; goto CLEANUP;
        }
        l = CCutil_dat_edgelen (cycle[ncount - 1], cycle[0], dat);
        if (CCutil_swrite_int (sout, l)) {
            fprintf (stderr, "CCutil_swrite_int failed\n");
            rval = 1; goto CLEANUP;
        }
    } else {
        fout = fopen (outedgename, "w");
        if (fout == (FILE *) NULL) {
            perror (outedgename);
            fprintf (stderr, "Unable to open %s for output\n", outedgename);
            rval = 1; goto CLEANUP;
        }
        fprintf (fout, "%d %d\n", ncount, ncount);
        for (i = 1; i < ncount; i++) {
            fprintf (fout, "%d %d %d\n", cycle[i - 1], cycle[i],
                     CCutil_dat_edgelen (cycle[i - 1], cycle[i], dat));
        }
        fprintf (fout, "%d %d %d\n", cycle[ncount - 1], cycle[0],
                 CCutil_dat_edgelen (cycle[ncount - 1], cycle[0], dat));
    }
    rval = 0;

CLEANUP:
    CCutil_sclose (sout);
    if (fout != (FILE *) NULL) fclose (fout);
    
    return rval;
}

int CCutil_writecycle (int ncount, char *outcyclename, int *cycle,
        int binary_out)
{
    FILE *cycfout = (FILE *) NULL;
    CC_SFILE *cycsout = (CC_SFILE *) NULL;
    int i;
    int rval;

    if (binary_out) {
        cycsout = CCutil_sopen (outcyclename, "w");
        if (cycsout == (CC_SFILE *) NULL) {
            fprintf (stderr, "Unable to open %s for output\n", outcyclename);
            rval = 1; goto CLEANUP;
        }
        if (CCutil_swrite_int (cycsout, ncount)) {
            fprintf (stderr, "CCutil_swrite_int failed\n");
            rval = 1; goto CLEANUP;
        }
        for (i = 0; i < ncount; i++) {
            if (CCutil_swrite_int (cycsout, cycle[i])) {
                fprintf (stderr, "CCutil_swrite_int failed\n");
                rval = 1; goto CLEANUP;
            }
        }
    } else {
        cycfout = fopen (outcyclename, "w");
        if (cycfout == (FILE *) NULL) {
            perror (outcyclename);
            fprintf (stderr, "Unable to open %s for output\n", outcyclename);
            rval = 1; goto CLEANUP;
        }
        fprintf (cycfout, "%d\n", ncount);
        for (i = 0; i < ncount; i++) {
            fprintf (cycfout, "%d ", cycle[i]);
            if (i % 10 == 9)
                fprintf (cycfout, "\n");
        }
        if (i % 10)
            fprintf (cycfout, "\n");
    }
    rval = 0;

CLEANUP:
    CCutil_sclose (cycsout);
    if (cycfout != (FILE *) NULL) fclose (cycfout);

    return rval;
}

int CCutil_writeedges_double (int ncount, char *outedgename, int ecount,
        int *elist, double *elen, int binary_out)
{
    int i;

    if (binary_out) {
        CC_SFILE *out = CCutil_sopen (outedgename, "w");

        if (out == (CC_SFILE *) NULL) {
            fprintf (stderr, "Unable to open %s for output\n", outedgename);
            return 1;
        }

        if (CCutil_swrite_int (out, ncount)) {
            CCutil_sclose (out);
            return 1;
        }
        if (CCutil_swrite_int (out, ecount)) {
            CCutil_sclose (out);
            return 1;
        }

        for (i = 0; i < ecount; i++) {
            if (CCutil_swrite_int (out, elist[2 * i])) {
                CCutil_sclose (out);
                return 1;
            }
            if (CCutil_swrite_int (out, elist[(2 * i) + 1])) {
                CCutil_sclose (out);
                return 1;
            }
            if (CCutil_swrite_double (out, elen[i])) {
                CCutil_sclose (out);
                return 1;
            }
        }
        CCutil_sclose (out);
    } else {
        FILE *out = fopen (outedgename, "w");

        if (out == (FILE *) NULL) {
            perror (outedgename);
            fprintf (stderr, "Unable to open %s for output\n", outedgename);
            return 1;
        }

        fprintf (out, "%d %d\n", ncount, ecount);
        for (i = 0; i < ecount; i++) {
            fprintf (out, "%d %d %f\n", elist[2 * i], elist[(2 * i) + 1],
                                        elen[i]);
        }
        fclose (out);
    }

    return 0;
}

int CCutil_tri2dat (int ncount, int *elen, CCdatagroup *dat)
{
    int i;
    int j;
    int k;
    int rval;
    
    CCutil_init_datagroup (dat);
    dat->adj = CC_SAFE_MALLOC (ncount, int *);
    dat->adjspace = CC_SAFE_MALLOC (ncount * (ncount+1) / 2, int);
    if (dat->adj == (int **) NULL ||
        dat->adjspace == (int *) NULL) {
        fprintf (stderr, "Our of memory in CCutil_tri2dat\n");
        rval = 1;
        goto CLEANUP;
    }

    for (i=0, j=0; i<ncount; i++) {
        dat->adj[i] = dat->adjspace + j;
        j += (i+1);
    }
    for (i=0, k=0; i<ncount; i++) {
        for (j=0; j<i; j++) {
            dat->adj[i][j] = elen[k++];
        }
        dat->adj[i][i] = 0;
    }
    if (CCutil_dat_setnorm (dat, CC_MATRIXNORM)) {
        fprintf (stderr, "CCutil_dat_setnorm failed\n");
        rval = 1;
        goto CLEANUP;
    }
    rval = 0;

  CLEANUP:
    if (rval) {
        CCutil_freedatagroup (dat);
    }
    return rval;
}

int CCutil_graph2dat_matrix (int ncount, int ecount, int *elist, int *elen,
        int defaultlen, CCdatagroup *dat)
{
    int i;
    int j;
    int k;
    int rval;
    
    CCutil_init_datagroup (dat);
    dat->adj = CC_SAFE_MALLOC (ncount, int *);
    dat->adjspace = CC_SAFE_MALLOC (ncount * (ncount+1) / 2, int);
    if (dat->adj == (int **) NULL ||
        dat->adjspace == (int *) NULL) {
        fprintf (stderr, "Our of memory in CCutil_graph2dat\n");
        rval = 1;
        goto CLEANUP;
    }

    for (i=0, j=0; i<ncount; i++) {
        dat->adj[i] = dat->adjspace + j;
        j += (i+1);
    }
    for (i=0, k=0; i<ncount; i++) {
        for (j=0; j<i; j++) {
            dat->adj[i][j] = defaultlen;
        }
        dat->adj[i][i] = 0;
    }
    for (i=0; i<ecount; i++) {
        j = elist[2*i];
        k = elist[2*i+1];
        if (j < k) dat->adj[k][j] = elen[i];
        else dat->adj[j][k] = elen[i];
    }
    if (CCutil_dat_setnorm (dat, CC_MATRIXNORM)) {
        fprintf (stderr, "CCutil_dat_setnorm failed\n");
        rval = 1;
        goto CLEANUP;
    }
    rval = 0;

CLEANUP:
    if (rval) {
        CCutil_freedatagroup (dat);
    }
    return rval;
}

int CCutil_graph2dat_sparse (int ncount, int ecount, int *elist, int *elen,
        int defaultlen, CCdatagroup *dat)
{
    int rval = 0;

    CCutil_init_datagroup (dat);
    rval = CCutil_dat_setnorm (dat, CC_SPARSE);
    if (rval) {
        fprintf (stderr, "CCutil_dat_setnorm failed\n"); goto CLEANUP;
    }
    rval = build_sparse_dat (ncount, ecount, elist, elen, dat, defaultlen);
    if (rval) {
        fprintf (stderr, "build_sparse_dat failed\n"); goto CLEANUP;
    }

CLEANUP:

    if (rval) {
        CCutil_freedatagroup (dat);
    }
    return rval;
}

static int build_sparse_dat (int ncount, int ecount, int *elist, int *elen,
        CCdatagroup *dat, int defaultlen)
{
    int i, k, head, tail, temp;

    dat->adjspace   = CC_SAFE_MALLOC (ecount, int);
    dat->lenspace   = CC_SAFE_MALLOC (ecount, int);
    dat->adj        = CC_SAFE_MALLOC (ncount, int *);
    dat->len        = CC_SAFE_MALLOC (ncount, int *);
    dat->degree     = CC_SAFE_MALLOC (ncount, int);

    if (!dat->adj || !dat->adjspace || !dat->len || !dat->lenspace ||
        !dat->degree) {
        fprintf (stderr, "out of memory in build_sparse_dat\n");
        CCutil_freedatagroup (dat);
        return 1;
    }

    for (i = 0; i < ncount; i++) {
        dat->degree[i] = 0;
    }

    for (i = 0; i < ecount; i++) {
        head = elist[2*i];
        tail = elist[2*i+1];
        if (head > tail) {
            CC_SWAP (head, tail, temp);
        }
        dat->degree[head]++;
    }

    for (i = 0, k = 0; i < ncount; i++) {
        dat->adj[i] = &(dat->adjspace[k]);
        dat->len[i] = &(dat->lenspace[k]);
        k += dat->degree[i];
        dat->degree[i] = 0;
    }

    for (i = 0; i < ecount; i++) {
        head = elist[2*i];
        tail = elist[2*i+1];
        if (head > tail) {
            CC_SWAP (head, tail, temp);
        }
        dat->adj[head][dat->degree[head]] = tail;
        dat->len[head][dat->degree[head]] = elen[i];
        dat->degree[head]++;
    }


    if (defaultlen <= 0) {
        double v;

        dat->default_len = 0;
        for (i = 0; i < ecount; i++) {
            if (elen[i] > dat->default_len) {
                dat->default_len = elen[i];
            }
        }
        v =  (double) (dat->default_len + 1) *  (double) ncount;
        if (256 * v > (double) CCutil_MAXINT) {
            printf ("WARNING: Large edge lengths in sparse graph\n");
            fflush (stdout);
            dat->default_len = CCutil_MAXINT / 256;
        } else {
            dat->default_len = (dat->default_len + 1) * ncount;
        }
        printf ("Default Edge Length: %d\n", dat->default_len); fflush (stdout);
        fflush (stdout);
    } else {
        dat->default_len = defaultlen;
    }
    dat->sparse_ecount = ecount;

    return 0;
}

int CCutil_get_sparse_dat_edges (int ncount, CCdatagroup *dat, int *ecount,
    int **elist, int **elen)
{
    int i, j, k, norm, rval = 0;

    CCutil_dat_getnorm (dat, &norm);
    if (norm != CC_SPARSE) {
        fprintf (stderr, "CCutil_get_sparse_dat_edges called with norm %d\n",
                          norm);
        rval = 1; goto CLEANUP;
    }
 
    *ecount = dat->sparse_ecount;
    *elist = CC_SAFE_MALLOC (2*(*ecount), int);
    *elen  = CC_SAFE_MALLOC (*ecount, int);

    if (!*elist || !(*elen)) {
        fprintf (stderr, "out of memory in CCutil_get_sparse_dat_edges\n");
        rval = 1; goto CLEANUP;
    }

    for (i = 0, k = 0; i < ncount; i++) {
        for (j = 0; j < dat->degree[i]; j++) {
            (*elist)[2*k] = i;
            (*elist)[2*k + 1] = dat->adj[i][j];
            (*elen)[k] = dat->len[i][j];
            k++;
        }
    }

CLEANUP:

    if (rval) {
        CC_IFFREE (*elist, int);
        CC_IFFREE (*elen, int);
        *ecount = 0;
    }

    return rval;
}

int CCutil_sparse_strip_edges (CCdatagroup *dat, int in_ecount,
    int *in_elist, int *in_elen, int *ecount, int **elist, int **elen)
{
    int i, count, norm, rval = 0;

    *ecount = 0;
    *elist = (int *) NULL;
    *elen  = (int *) NULL;

    CCutil_dat_getnorm (dat, &norm);
    if (norm != CC_SPARSE) {
        fprintf (stderr, "CCutil_sparse_strip_edges called with norm %d\n",
                 norm);
        rval = 1; goto CLEANUP;
    }

    for (i = 0, count = 0; i < in_ecount; i++) {
        if (in_elen[i] != dat->default_len) count++;
    }

    if (count > 0) {
        *ecount = count;
        *elist = CC_SAFE_MALLOC (2*count, int);
        *elen  = CC_SAFE_MALLOC (count, int);
        if (!(*elist) || !(*elen)) {
            fprintf (stderr, "out of memory in CCutil_sparse_strip_edges\n");
            rval = 1; goto CLEANUP;
        }

        for (i = 0, count = 0; i < in_ecount; i++) {
            if (in_elen[i] != dat->default_len) {
                (*elist)[2*count] = in_elist[2*i];
                (*elist)[2*count+1] = in_elist[2*i+1];
                (*elen)[count++] = in_elen[i];
            }
        }
    }

CLEANUP:
 
    if (rval) {
        *ecount = 0;
        CC_IFFREE (*elist, int);
        CC_IFFREE (*elen, int);
    }
    return rval;
}



int CCutil_sparse_real_tour (int ncount, CCdatagroup *dat, int *cyc,
    int *yesno)
{
    int norm, i, rval = 0;

    *yesno = 0;

    CCutil_dat_getnorm (dat, &norm);
    if (norm != CC_SPARSE) {
        fprintf (stderr, "CCutil_sparse_real_tour called with norm %d\n", norm);
        rval = 1; goto CLEANUP;
    }

    for (i = 1; i < ncount; i++) {
        if (CCutil_dat_edgelen (cyc[i-1], cyc[i], dat) == dat->default_len) {
            goto CLEANUP; 
        }
    }
    if (CCutil_dat_edgelen (cyc[ncount-1], cyc[0], dat) == dat->default_len) {
        goto CLEANUP; 
    }

    *yesno = 1;

CLEANUP:

    return rval;
}

static int writedat_sparse (CC_SFILE *out, int ncount, CCdatagroup *dat)
{
    int ecount = dat->sparse_ecount;
    int i, j;

    if (CCutil_swrite_int (out, ecount)) {
        fprintf (stderr, "CCutil_swrite_int failed\n"); return 1;
    }
    
    for (i = 0; i < ncount; i++) {
        if (CCutil_swrite_int (out, dat->degree[i])) {
            fprintf (stderr, "CCutil_swrite_int failed\n"); return 1;
        }
        for (j = 0; j < dat->degree[i]; j++) {
            if (CCutil_swrite_int (out, dat->adj[i][j])) {
                fprintf (stderr, "CCutil_swrite_int failed\n"); return 1;
            }
            if (CCutil_swrite_int (out, dat->len[i][j])) {
                fprintf (stderr, "CCutil_swrite_int failed\n"); return 1;
            }
        }
    }
    return 0;
}

static int readdat_sparse (CC_SFILE *in, int ncount, CCdatagroup *dat)
{
    int i, j, k, ecount;
    int rval = 0;

    rval = CCutil_sread_int (in, &ecount);
    if (rval) {
        fprintf (stderr, "CCutil_sread_int failed\n"); return rval;
    }

    dat->adjspace = CC_SAFE_MALLOC (ecount, int);
    dat->lenspace = CC_SAFE_MALLOC (ecount, int);
    dat->adj      = CC_SAFE_MALLOC (ncount, int *);
    dat->len      = CC_SAFE_MALLOC (ncount, int *);
    dat->degree   = CC_SAFE_MALLOC (ncount, int);
   
    if (!dat->adjspace || !dat->lenspace || !dat->adj || !dat->len ||
        !dat->degree) {
        fprintf (stderr, "out of memory in readdat_sparse\n");
        rval = 1; goto CLEANUP;
    }

    for (i = 0, k = 0; i < ncount; i++) {
        dat->adj[i] = &(dat->adjspace[k]);
        dat->len[i] = &(dat->lenspace[k]);

        rval = CCutil_sread_int (in, &(dat->degree[i]));
        if (rval) {
            fprintf (stderr, "CCutil_sread_int failed\n"); goto CLEANUP;
        }
        k += dat->degree[i];

        for (j = 0; j < dat->degree[i]; j++) {
            rval = CCutil_sread_int (in, &(dat->adj[i][j]));
            if (rval) {
                fprintf (stderr, "CCutil_sread_int failed\n"); goto CLEANUP;
            }
            rval = CCutil_sread_int (in, &(dat->len[i][j]));
            if (rval) {
                fprintf (stderr, "CCutil_sread_int failed\n"); goto CLEANUP;
            }
        }
    }
    dat->sparse_ecount = ecount;

CLEANUP:
    if (rval) {
        CCutil_freedatagroup (dat);
    }
    return rval;
}

static int permute_sparse (int ncount, CCdatagroup *dat, int *perm)
{
    int i, j, k, ecount = dat->sparse_ecount;
    int rval = 0;
    int **tadj     = (int **) NULL;
    int *tadjspace = (int *) NULL;
    int **tlen     = (int **) NULL;
    int *tlenspace = (int *) NULL;
    int *tdegree   = (int *) NULL; 
    int *inv       = (int *) NULL;

    tadjspace = CC_SAFE_MALLOC (ecount, int);
    tlenspace = CC_SAFE_MALLOC (ecount, int);
    tadj      = CC_SAFE_MALLOC (ncount, int *);
    tlen      = CC_SAFE_MALLOC (ncount, int *);
    tdegree   = CC_SAFE_MALLOC (ncount, int);
    inv       = CC_SAFE_MALLOC (ncount, int);

    if (!tdegree || !tadj || !tadjspace || !tlen || !tlenspace || !inv) {
        fprintf (stderr, "out of memory in permute_sparse\n");
        rval = 1; goto CLEANUP;
    }

    for (i = 0; i < ncount; i++) {
        tdegree[i] = 0;
        inv[perm[i]] = i;
    }
    for (i = 0; i < ncount; i++) {
        for (j = 0; j < dat->degree[i]; j++) {
            if (inv[i] < inv[dat->adj[i][j]]) {
                tdegree[inv[i]]++;
            } else {
                tdegree[inv[dat->adj[i][j]]]++;
            }
        }
    }

    for (i = 0, k = 0; i < ncount; i++) {
        tadj[i] = &(tadjspace[k]);
        tlen[i] = &(tlenspace[k]);
        k += tdegree[i];
        tdegree[i] = 0;
    }
   
    for (i = 0; i < ncount; i++) {
        for (j = 0; j < dat->degree[i]; j++) {
            k = inv[dat->adj[i][j]];
            if (inv[i] < k) {
                tadj[inv[i]][tdegree[inv[i]]] = k;
                tlen[inv[i]][tdegree[inv[i]]] = dat->len[i][j];
                tdegree[inv[i]]++;
            } else {
                tadj[k][tdegree[k]] = inv[i];
                tlen[k][tdegree[k]] = dat->len[i][j];
                tdegree[k]++;
            }
        }
    }

    CC_IFFREE (dat->degree, int);
    dat->degree = tdegree;
    CC_IFFREE (dat->adj, int *);
    dat->adj = tadj;
    CC_IFFREE (dat->adjspace, int);
    dat->adjspace = tadjspace;
    CC_IFFREE (dat->len, int *);
    dat->len = tlen;
    CC_IFFREE (dat->lenspace, int);
    dat->lenspace = tlenspace;

CLEANUP:

    if (rval) {
        CC_IFFREE (tdegree, int);
        CC_IFFREE (tadj, int *);
        CC_IFFREE (tadjspace, int);
        CC_IFFREE (tlen, int *);
        CC_IFFREE (tlenspace, int);
    }
    CC_IFFREE (inv, int);
    return rval;
}

static int copy_sparse (int ncount, CCdatagroup *indat, CCdatagroup *outdat)
{
    int i, k, ecount = indat->sparse_ecount;
    int rval = 0;

    CCutil_init_datagroup (outdat);
    CCutil_dat_setnorm (outdat, CC_SPARSE);

    outdat->adjspace = CC_SAFE_MALLOC (ecount, int);
    outdat->lenspace = CC_SAFE_MALLOC (ecount, int);
    outdat->adj      = CC_SAFE_MALLOC (ncount, int *);
    outdat->len      = CC_SAFE_MALLOC (ncount, int *);
    outdat->degree   = CC_SAFE_MALLOC (ncount, int);

    if (!outdat->degree || !outdat->adj || !outdat->adjspace || !outdat->len ||
        !outdat->lenspace) {
        fprintf (stderr, "out of memory in permute_sparse\n");
        rval = 1; goto CLEANUP;
    }

    for (i = 0, k = 0; i < ncount; i++) {
        outdat->degree[i] = indat->degree[i];
        outdat->adj[i] = &(outdat->adjspace[k]);
        outdat->len[i] = &(outdat->lenspace[k]);
        k += outdat->degree[i];
    }
    for (i = 0; i < ecount; i++) {
        outdat->adjspace[i] = indat->adjspace[i];
        outdat->lenspace[i] = indat->lenspace[i];
    }
    outdat->sparse_ecount = ecount;
    outdat->default_len = indat->default_len;

CLEANUP:

    if (rval) {
        CCutil_freedatagroup (outdat);
    }
    return rval;
}

/* The crystal functions are based on David Shallcross's fortran code */

#define CRYSTAL_SCALE 10000

typedef struct three_d {
    double phi;
    double chi;
    double twoth;
} three_d;



static void
    cry_quicksort (three_d *x, int l, int u),
    cryswap (three_d *x, three_d *y);
static int
    point_compare (three_d *x, three_d *y),
    crygenpts (double orient[3][3], double lambda, int bounds[3][2],
        three_d **p_crypoints, int *p_ncrypoints),
    cryangles (int h, int k, int l, double orient[3][3], double lambda,
            double omega, double *phi, double *chi, double *twoth);
static double
    userint (double f);



#define CRY_MIN_ORIENT -0.2    /* For random problems */
#define CRY_MAX_ORIENT  0.2

#define CRY_MIN_BOUND  10
#define CRY_MAX_BOUND  20

#define CRY_1_00_CUTOFF 15000
#define CRY_1_35_CUTOFF 10000

static int read_crystal (char *datname, int binary_in, int *ncount,
        CCdatagroup *dat, CCrandstate *rstate)
{
    FILE *datin;
    double lambda;
    double orient[3][3];
    int bounds[3][2];
    int i, j;
    three_d *crypoints = (three_d *) NULL;
    int ncrypoints = 0;

    if (!datname) {
        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                orient[i][j] = CRY_MIN_ORIENT +
                   ((CRY_MAX_ORIENT - CRY_MIN_ORIENT) *
                    ((double) (CCutil_lprand (rstate) % 1000) / 1000.0));
            }
        }
        for (i = 0; i < 3; i++) {
            bounds[i][1] = CRY_MIN_BOUND +
               (CCutil_lprand (rstate) % (CRY_MAX_BOUND - CRY_MIN_BOUND + 1));
            bounds[i][0] = -bounds[i][1];
        }

        if (*ncount > CRY_1_00_CUTOFF)
            lambda = 1.0;
        else if (*ncount > CRY_1_35_CUTOFF)
            lambda = 1.35;
        else
            lambda = 1.70;

        printf ("Random crystal problem\n");
        printf ("Note that the number of nodes will not match the request\n");
        printf ("Orient:\n");
        for (i = 0; i < 3; i++)
            printf (" %.4f  %.4f  %.4f\n", orient[i][0], orient[i][1],
                                                         orient[i][2]);
        printf ("Bounds:\n");
        for (i = 0; i < 3; i++)
            printf (" %d %d ", bounds[i][0], bounds[i][1]);
        printf ("\nWavelength:\n");
            printf (" %.2f\n", lambda);
        fflush (stdout);
    } else if (binary_in) {
        fprintf (stderr, "CRYSTAL norms do not support binary input\n");
        return 1;
    } else {
        datin = fopen (datname, "r");
        if (datin == (FILE *) NULL) {
            perror (datname);
            fprintf (stderr, "Unable to open %s for input\n", datname);
            return 1;
        }

        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                fscanf (datin, "%lf", &orient[i][j]);
            }
        }
        fscanf (datin, "%lf", &lambda);
        for (i = 0; i < 3; i++) {
            for (j = 0; j < 2; j++) {
                fscanf (datin, "%d", &(bounds[i][j]));
            }
        }
        fclose (datin);
    }

    if (crygenpts (orient, lambda, bounds, &crypoints, &ncrypoints)) {
        fprintf (stderr, "crygenpts failed\n");
        if (crypoints)
            CC_FREE (crypoints, three_d);
        return 1;
    }
    cry_quicksort (crypoints, 0, ncrypoints -1);
    printf ("Number of crystal points: %d\n", ncrypoints);

    dat->x = CC_SAFE_MALLOC (ncrypoints, double);
    if (!dat->x) {
        CC_FREE (crypoints, three_d);
        return 1;
    }
    dat->y = CC_SAFE_MALLOC (ncrypoints, double);
    if (!dat->y) {
        CC_FREE (crypoints, three_d);
        CCutil_freedatagroup (dat);
        return 1;
    }
    dat->z = CC_SAFE_MALLOC (ncrypoints, double);
    if (!dat->z) {
        CC_FREE (crypoints, three_d);
        CCutil_freedatagroup (dat);
        return 1;
    }

    for (i = 0; i < ncrypoints; i++) {
        dat->x[i] = userint (crypoints[i].chi * CRYSTAL_SCALE),
        dat->y[i] = userint (crypoints[i].phi * CRYSTAL_SCALE),
        dat->z[i] = userint (crypoints[i].twoth * CRYSTAL_SCALE);
    }
    *ncount = ncrypoints;

    if (crypoints)
        CC_FREE (crypoints, three_d);

    return 0;
}

static int point_compare (three_d *x, three_d *y)
{
    int del;

    del = (int) userint (x->chi * CRYSTAL_SCALE) -
          (int) userint (y->chi * CRYSTAL_SCALE);
    if (del)
        return del;
    del = (int) userint (x->phi * CRYSTAL_SCALE) -
          (int) userint (y->phi * CRYSTAL_SCALE);
    if (del)
        return del;
    del = (int) userint (x->twoth * CRYSTAL_SCALE) -
          (int) userint (y->twoth * CRYSTAL_SCALE);
    return del;
}

static double userint (double f)
{
    double t;
    int u;

    if (f > 0.0) {
        u = (int) f;
        t = (double) u;
        return (f - t < 0.5 ? t : t + 1.0);
    } else {
        u = (int) f;
        t = (double) u;
        return (t - f < 0.5 ? t : t - 1.0);
    }
}


/* genpts builds points, the array of accessible points of the lattice */

static int crygenpts (double orient[3][3], double lambda, int bounds[3][2],
        three_d **p_crypoints, int *p_ncrypoints)
{
    int h, k, l;
    double p = 0.0, c = 0.0, t = 0.0;
    double minp = 360.0, minc = 180.0, mint = 155.0;
    double maxp = 0.0, maxc = -180.0, maxt = -55.0;
    three_d *crypoints;
    int ncrypoints;

    *p_crypoints = (three_d *) NULL;
    *p_ncrypoints = 0;

    ncrypoints = 0;

    crypoints = CC_SAFE_MALLOC (100000, three_d);
    if (!crypoints)
        return 1;

    for (h = bounds[0][0]; h <= bounds[0][1]; h++) {
        for (k = bounds[1][0]; k <= bounds[1][1]; k++) {
            for (l = bounds[2][0]; l <= bounds[2][1]; l++) {
                if (cryangles (h, k, l, orient, lambda, 0.0, &p, &c, &t)) {
                    crypoints[ncrypoints].phi = p;
                    crypoints[ncrypoints].chi = c;
                    crypoints[ncrypoints].twoth = t;
                    ncrypoints++;
                    if (p < minp)
                        minp = p;
                    if (p > maxp)
                        maxp = p;
                    if (c < minc)
                        minc = c;
                    if (c > maxc)
                        maxc = c;
                    if (t < mint)
                        mint = t;
                    if (t > maxt)
                        maxt = t;
                }
            }
        }
    }

    *p_crypoints = crypoints;
    *p_ncrypoints = ncrypoints;
    return 0;
}

/*
   angles computes positioning information for the detector given the
   miller indices.  (From Matt Small, April 5, 1984)
*/

static int cryangles (int h, int k, int l, double orient[3][3], double lambda,
                  double omega, double *phi, double *chi, double *twoth)
{
    /* The original had pi = 3.14159265368979; */
    const double pi = 3.14159265358979;
    double x, y, z, d, dum, cosomg, sinchi, sinphi, cosphi, q, r;
    double rh = h, rk = k, rl = l;
    double t1, t2, t3;

    omega = omega / 180.0 * pi;
    cosomg = cos (omega);

    x = orient[0][0] * rh + orient[0][1] * rk + orient[0][2] * rl;
    y = orient[1][0] * rh + orient[1][1] * rk + orient[1][2] * rl;
    z = orient[2][0] * rh + orient[2][1] * rk + orient[2][2] * rl;
    d = sqrt (x * x + y * y + z * z);
    dum = lambda * d / 2.0;
    if (dum < .000001 || dum >= 1.0) {
        return 0;
    }
    *twoth = asin (dum) * 2.0;
    t1 = x / d;
    t2 = y / d;
    t3 = z / d;
    if ((t3 < 0.0 ? -t3 : t3) >= (cosomg < 0.0 ? -cosomg : cosomg)) {
        return 0;
    }
    sinchi = -t3 / cosomg;
    *chi = asin (sinchi);
    q = sin (omega);
    r = cosomg * cos (*chi);
    cosphi = (q * t1 + r * t2) / (t1 * t1 + t2 * t2);
    sinphi = (r * t1 - q * t2) / (t1 * t1 + t2 * t2);
    if (sinphi <= -.7)
        *phi = -acos (cosphi);
    else if (sinphi >= .7)
        *phi = acos (cosphi);
    else if (cosphi <= 0.0)
        *phi = pi - asin (sinphi);
    else
        *phi = asin (sinphi);
    *phi *= 180.0 / pi;
    *chi *= 180.0 / pi;
    *twoth *= 180.0 / pi;
    *phi /= 1.25;   /* FOR VARYING MORTOR SPEEDS */
    *chi /= 1.5;
    *twoth /= 1.15;
    return 1;
}

static void cry_quicksort (three_d *x, int l, int u)
{
    int i, j;
    three_d *t;

    if (l >= u)
        return;

    cryswap (x + l, x + ((l+u)/2));

    i = l;
    j = u + 1;
    t = x + l;

    while (1) {
        do i++; while (i <= u && (point_compare (x + i, t) < 0));
        do j--; while (point_compare (x + j, t) > 0);
        if (j < i) break;
        cryswap (x + i, x + j);
    }
    cryswap (x + l, x + j);
    cry_quicksort (x, l, j - 1);
    cry_quicksort (x, i, u);
}

static void cryswap (three_d *x, three_d *y)
{
    three_d t;

    t.phi = x->phi;
    t.chi = x->chi;
    t.twoth = x->twoth;

    x->phi = y->phi;
    x->chi = y->chi;
    x->twoth = y->twoth;

    y->phi = t.phi;
    y->chi = t.chi;
    y->twoth = t.twoth;
}

static int read_d2 (char *datname, int binary_in, int *ncount,
        CCdatagroup *dat, CCrandstate *rstate, int innorm, int gridsize,
        int allow_dups)
{
    if (innorm == CC_EUCTOROIDAL) {
        dat->gridsize = gridsize;
        printf ("Using TOROIDAL gridsize %.0f\n", dat->gridsize);
    }
    if (datname == (char *) NULL) {
        return build_d2 (ncount, dat, rstate, innorm, gridsize, allow_dups);
    } else if (binary_in) {
        return read_d2_binary (datname, ncount, dat, binary_in);
    } else {
        return read_d2_text (datname, ncount, dat);
    }
}

static int write_d2 (char *datname, int binary_out, int ncount,
        CCdatagroup *dat)
{
    if (binary_out) {
        return write_d2_binary (datname, ncount, dat);
    } else {
        return write_d2_text (datname, ncount, dat);
    }
}

static int build_d2 (int *ncount, CCdatagroup *dat, CCrandstate *rstate,
        int innorm, int gridsize, int allow_dups)
{
    int i;
    int j;
    
    printf ("Random %d point set\n", *ncount);
    fflush (stdout);
    dat->x = CC_SAFE_MALLOC (*ncount, double);
    if (!dat->x)
        return 1;
    dat->y = CC_SAFE_MALLOC (*ncount, double);
    if (!dat->y) {
        CCutil_freedatagroup (dat);
        return 1;
    }
    if (innorm == CC_GEOGRAPHIC || innorm == CC_GEOM) {
        for (i = 0; i < *ncount; i++) {
            dat->x[i] = (double) (CCutil_lprand (rstate) % 180) - 90.0;
            dat->y[i] = (double) (CCutil_lprand (rstate) % 360) - 180.0;
        }
    } else {
        int **hit = (int **) NULL;
        int *hitcount = (int *) NULL;
        int winner, x, y;

        if (!allow_dups) {
            hit      = CC_SAFE_MALLOC (*ncount, int *);
            hitcount = CC_SAFE_MALLOC (*ncount, int);
            if (!hit || !hitcount) {
                fprintf (stderr, "out of memory in CCutil_getdata\n");
                CC_IFFREE (hit, int *);
                CC_IFFREE (hitcount, int);
                CCutil_freedatagroup (dat);
                return 1;
            }
        
            for (i = 0; i < *ncount; i++) {
                hit[i] = (int *) NULL;
                hitcount[i] = 0;
            }
        }
        
        for (i = 0; i < *ncount; i++) {
            winner = 0;
            do {
                x = CCutil_lprand (rstate) % gridsize;
                y = CCutil_lprand (rstate) % gridsize;
                if (!allow_dups) {
                    for (j = 0; j < hitcount[x]; j++) {
                        if (hit[x][j] == y) break;
                    }
                    if (j == hitcount[x]) {
                        if (CCutil_reallocrus_count ((void **) &(hit[x]),
                                                     hitcount[x] + 1, sizeof (int))) {
                            fprintf (stderr, "CCutil_reallocrus_count failed\n");
                            for (i = 0; i < *ncount; i++) {
                                CC_IFFREE (hit[i], int);
                            }
                            CC_IFFREE (hit, int *);
                            CC_IFFREE (hitcount, int);
                            CCutil_freedatagroup (dat);
                            return 1;
                        }
                        hit[x][hitcount[x]] = y;
                        hitcount[x]++;
                        winner = 1;
                    }
                    if (!winner) {
                        printf ("X"); fflush (stdout);
                    }
                }
            } while (!allow_dups && !winner);
            dat->x[i] = (double) x;
            dat->y[i] = (double) y;
        }
        if (!allow_dups) {
            for (i = 0; i < *ncount; i++) {
                CC_IFFREE (hit[i], int);
            }
        }
        CC_IFFREE (hit, int *);
        CC_IFFREE (hitcount, int);
    }
   
/*
    printf ("%d\n", *ncount);
    for (i = 0; i < *ncount; i++) {
        printf ("%d %d\n", (int) dat->x[i], (int) dat->y[i]);

    }
*/
    return 0;
}

static int read_d2_binary (char *datname, int *ncount, CCdatagroup *dat,
        int btype)
{
    int i, xi, yi, tval, rval = 0;
    CC_SFILE *f = CCutil_sopen (datname, "r");

    if (f == (CC_SFILE *) NULL) {
        fprintf (stderr, "could not open %s for reading\n", datname);
        rval = 1;  goto CLEANUP;
    }

    rval = CCutil_sread_int (f, ncount);
    CCcheck_rval (rval, "CCutil_sread_int failed");

    printf ("nnodes = %d\n", *ncount); fflush (stdout);

    dat->x = CC_SAFE_MALLOC (*ncount, double);
    CCcheck_NULL (dat->x, "out of memory in read_d2_binary");
    dat->y = CC_SAFE_MALLOC (*ncount, double);
    CCcheck_NULL (dat->y, "out of memory in read_d2_binary");
  
    if (btype == 2) {
        for (i = 0; i < *ncount; i++) {
            rval = CCutil_sread_double (f, &(dat->x[i]));
            CCcheck_rval (rval, "CCutil_sread_double failed");
            rval = CCutil_sread_double (f, &(dat->y[i]));
            CCcheck_rval (rval, "CCutil_sread_double failed");
        }
    } else {
        for (i = 0; i < *ncount; i++) {
            rval = CCutil_sread_int (f, &xi);
            CCcheck_rval (rval, "CCutil_sread_int failed");
            rval = CCutil_sread_int (f, &yi);
            CCcheck_rval (rval, "CCutil_sread_int failed");
            dat->x[i] = (double) xi;
            dat->y[i] = (double) yi;
        }
    }

CLEANUP:

    if (f) {
        tval = CCutil_sclose (f);
        if (tval) {
            fprintf (stderr, "CCutil_sclose failed\n");
            if (rval == 0)  rval = tval;
        }
    }
    if (rval) {
        CCutil_freedatagroup (dat);
    }
    return rval;
}

static int read_d2_text (char *datname, int *ncount, CCdatagroup *dat)
{
    int i;
    FILE *datin = fopen (datname, "r");
    
    if (datin == (FILE *) NULL) {
        perror (datname);
        fprintf (stderr, "Unable to open %s for input\n", datname);
        return 1;
    }
    fscanf (datin, "%d", ncount);
    printf ("nnodes = %d\n", *ncount);
    dat->x = CC_SAFE_MALLOC (*ncount, double);
    if (!dat->x) {
        fclose (datin);
        return 1;
    }
    dat->y = CC_SAFE_MALLOC (*ncount, double);
    if (!dat->y) {
        fclose (datin);
        CCutil_freedatagroup (dat);
        return 1;
    }
    for (i = 0; i < *ncount; i++) {
        fscanf (datin, "%lf%lf", &(dat->x[i]), &(dat->y[i]));
    }
    fclose (datin);
    return 0;
}

static int write_d2_binary (char *datname, int ncount, CCdatagroup *dat)
{
    int i, tval, rval = 0;
    CC_SFILE *f = CCutil_sopen (datname, "w");

    if (f == (CC_SFILE *) NULL) {
        fprintf (stderr, "could not open %s for output\n", datname);
        rval = 1;  goto CLEANUP;
    }

    printf ("Creating binary-double-dat file, ncount %d\n", ncount);
    fflush (stdout);

    rval = CCutil_swrite_int (f, ncount);
    CCcheck_rval (rval, "CCutil_swrite_int failed");

    for (i = 0; i < ncount; i++) {
        rval = CCutil_swrite_double (f, dat->x[i]);
        CCcheck_rval (rval, "CCutil_sread_double failed");
        rval = CCutil_swrite_double (f, dat->y[i]);
        CCcheck_rval (rval, "CCutil_sread_double failed");
    }

CLEANUP:

    if (f) {
        tval = CCutil_sclose (f);
        if (tval) {
            fprintf (stderr, "CCutil_sclose failed\n");
            if (rval == 0)  rval = tval;
        }
    }
    return rval;
}
static int write_d2_text (char *datname, int ncount, CCdatagroup *dat)
{
    int i;
    FILE *datout = fopen (datname, "w");
    
    if (datout == (FILE *) NULL) {
        perror (datname);
        fprintf (stderr, "Unable to open %s for output\n", datname);
        return 1;
    }
    fprintf (datout, "%d\n", ncount);
    for (i = 0; i < ncount; i++) {
/*
        fprintf (datout, "%.12g %.12g\n", dat->x[i], dat->y[i]);
*/
        fprintf (datout, "%lf %lf\n", dat->x[i], dat->y[i]);
    }
    fclose (datout);
    return 0;
}

static int read_d3 (char *datname, int binary_in, int *ncount,
        CCdatagroup *dat, CCrandstate *rstate, int gridsize)
{
    if (datname == (char *) NULL) {
        return build_d3 (ncount, dat, rstate, gridsize);
    } else if (binary_in) {
        return read_d3_binary (datname, ncount, dat);
    } else {
        return read_d3_text (datname, ncount, dat);
    }
}

static int write_d3 (char *datname, int binary_out, int ncount,
        CCdatagroup *dat)
{
    if (binary_out) {
        fprintf (stderr, "Binary output of this norm is not yet implemented\n");
        return 1;
        
    } else {
        return write_d3_text (datname, ncount, dat);
    }
}

static int build_d3 (int *ncount, CCdatagroup *dat, CCrandstate *rstate,
        int gridsize)
{
    int i;
    
    printf ("Random %d point set\n", *ncount);
    fflush (stdout);
    dat->x = CC_SAFE_MALLOC (*ncount, double);
    if (!dat->x)
        return 1;
    dat->y = CC_SAFE_MALLOC (*ncount, double);
    if (!dat->y) {
        CCutil_freedatagroup (dat);
        return 1;
    }
    dat->z = CC_SAFE_MALLOC (*ncount, double);
    if (!dat->z) {
        CCutil_freedatagroup (dat);
        return 1;
    }
    for (i = 0; i < *ncount; i++) {
        dat->x[i] = (double) (CCutil_lprand (rstate) % gridsize);
        dat->y[i] = (double) (CCutil_lprand (rstate) % gridsize);
        dat->z[i] = (double) (CCutil_lprand (rstate) % gridsize);
    }
    return 0;
}

static int read_d3_binary (char *datname, int *ncount, CCdatagroup *dat)
{
    int i;
    int xi;
    int yi;
    int zi;
    CC_SFILE *f = CCutil_sopen (datname, "r");

    if (f == (CC_SFILE *) NULL)
        return 1;
    if (CCutil_sread_int (f, ncount)) {
        CCutil_sclose (f);
        return 1;
    }
    printf ("nnodes = %d\n", *ncount);
    fflush (stdout);
    dat->x = CC_SAFE_MALLOC (*ncount, double);
    if (!dat->x) {
        if (CCutil_sclose (f))
            fprintf (stderr, "Could not close file\n");
        return 1;
    }
    dat->y = CC_SAFE_MALLOC (*ncount, double);
    if (!dat->y) {
        if (CCutil_sclose (f))
            fprintf (stderr, "Could not close file\n");
        CCutil_freedatagroup (dat);
        return 1;
    }
    dat->z = CC_SAFE_MALLOC (*ncount, double);
    if (!dat->z) {
        if (CCutil_sclose (f))
            fprintf (stderr, "Could not close file\n");
        CCutil_freedatagroup (dat);
        return 1;
    }
    for (i = 0; i < *ncount; i++) {
        if (CCutil_sread_int (f, &xi)) {
            CCutil_sclose (f);
            CCutil_freedatagroup (dat);
            return 1;
        }
        if (CCutil_sread_int (f, &yi)) {
            CCutil_sclose (f);
            CCutil_freedatagroup (dat);
            return 1;
        }
        if (CCutil_sread_int (f, &zi)) {
            CCutil_sclose (f);
            CCutil_freedatagroup (dat);
            return 1;
        }
        dat->x[i] = (double) xi;
        dat->y[i] = (double) yi;
        dat->z[i] = (double) zi;
    }
    if (CCutil_sclose (f)) {
        CCutil_freedatagroup (dat);
        return 1;
    }
    return 0;
}

static int read_d3_text (char *datname, int *ncount, CCdatagroup *dat)
{
    int i;
    FILE *datin = fopen (datname, "r");
    
    if (datin == (FILE *) NULL) {
        perror (datname);
        fprintf (stderr, "Unable to open %s for input\n", datname);
        return 1;
    }
    fscanf (datin, "%d", ncount);
    printf ("nnodes = %d\n", *ncount);
    dat->x = CC_SAFE_MALLOC (*ncount, double);
    if (!dat->x) {
        fclose (datin);
        return 1;
    }
    dat->y = CC_SAFE_MALLOC (*ncount, double);
    if (!dat->y) {
        fclose (datin);
        CCutil_freedatagroup (dat);
        return 1;
    }
    dat->z = CC_SAFE_MALLOC (*ncount, double);
    if (!dat->z) {
        fclose (datin);
        CCutil_freedatagroup (dat);
        return 1;
    }
    for (i = 0; i < *ncount; i++) {
        fscanf (datin, "%lf%lf%lf",
                &(dat->x[i]), &(dat->y[i]), &(dat->z[i]));
    }
    fclose (datin);
    return 0;
}

static int write_d3_text (char *datname, int ncount, CCdatagroup *dat)
{
    int i;
    FILE *datout = fopen (datname, "w");
    
    if (datout == (FILE *) NULL) {
        perror (datname);
        fprintf (stderr, "Unable to open %s for output\n", datname);
        return 1;
    }
    fprintf (datout, "%d\n", ncount);
    for (i = 0; i < ncount; i++) {
        fprintf (datout, "%.12g %.12g %.12g\n", dat->x[i], dat->y[i],
                 dat->z[i]);
    }
    fclose (datout);
    return 0;
}

static int read_matrix (char *datname, int binary_in, int *ncount,
        CCdatagroup *dat, CCrandstate *rstate)
{
    /* Matrix is the lower triangle plus the diagonal */
    if (datname == (char *) NULL) {
        return build_matrix (ncount, dat, rstate);
    } else if (binary_in) {
        return read_matrix_binary (datname, ncount, dat);
    } else {
        return read_matrix_text (datname, ncount, dat);
    }
}

static int write_matrix (char *datname, int binary_out, int ncount,
        CCdatagroup *dat)
{
    /* Matrix is the lower triangle plus the diagonal */
    if (binary_out) {
        fprintf (stderr, "Binary output of this norm not yet implemented\n");
        return 1;
    } else {
        return write_matrix_text (datname, ncount, dat);
    }
}

static int build_matrix (int *ncount, CCdatagroup *dat, CCrandstate *rstate)
{
    int i;
    int j;
    
    printf ("Complete graph with %d nodes and random edge lengths\n",
            *ncount);
    fflush (stdout);
    dat->adj = CC_SAFE_MALLOC (*ncount, int *);
    dat->adjspace = CC_SAFE_MALLOC ((*ncount) * (*ncount+1) / 2,
                                    int);
    if (dat->adj == (int **) NULL ||
        dat->adjspace == (int *) NULL) {
        return 1;
    }
    for (i = 0, j = 0; i < *ncount; i++) {
        dat->adj[i] = dat->adjspace + j;
        j += (i+1);
    }
    for (i = 0; i < *ncount; i++) {
        for (j = 0; j < i; j++)
            dat->adj[i][j] = CCutil_lprand (rstate) %
                (MATRAND_SCALE * (*ncount));
        dat->adj[i][i] = 0;
    }
    return 0;
}

static int read_matrix_binary (char *datname, int *ncount, CCdatagroup *dat)
{
    int i;
    int j;
    CC_SFILE *f = CCutil_sopen (datname, "r");
    
    if (f == (CC_SFILE *) NULL)
        return 1;
    if (CCutil_sread_int (f, ncount)) {
        CCutil_sclose (f);
        return 1;
    }
    printf ("nnodes = %d\n", *ncount);
    fflush (stdout);
    dat->adj = CC_SAFE_MALLOC (*ncount, int *);
    dat->adjspace = CC_SAFE_MALLOC ((*ncount) * (*ncount+1) / 2,
                                    int);
    if (dat->adj == (int **) NULL ||
        dat->adjspace == (int *) NULL) {
        CC_IFFREE (dat->adj, int *);
        CC_IFFREE (dat->adjspace, int);
        if (CCutil_sclose (f))
            fprintf (stderr, "Could not close file\n");
        return 1;
    }
    for (i = 0, j = 0; i < *ncount; i++) {
        dat->adj[i] = dat->adjspace + j;
        j += (i+1);
    }
    for (i = 0; i < *ncount; i++) {
        for (j = 0; j <= i; j++) {
            if (CCutil_sread_int (f, &(dat->adj[i][j]))) {
                CCutil_sclose (f);
                CCutil_freedatagroup (dat);
                return 1;
            }
        }
    }
    if (CCutil_sclose (f)) {
        CCutil_freedatagroup (dat);
        return 1;
    }
    return 0;
}

static int read_matrix_text (char *datname, int *ncount, CCdatagroup *dat)
{
    int i;
    int j;
    FILE *datin = fopen (datname, "r");
    
    if (datin == (FILE *) NULL) {
        perror (datname);
        fprintf (stderr, "Unable to open %s for input\n", datname);
        return 1;
    }
    *ncount = CCutil_readint (datin);
    printf ("nnodes = %d\n", *ncount);
    dat->adj = CC_SAFE_MALLOC (*ncount, int *);
    dat->adjspace = CC_SAFE_MALLOC ((*ncount) * (*ncount+1) / 2,
                                    int);
    if (dat->adj == (int **) NULL ||
        dat->adjspace == (int *) NULL) {
        CC_IFFREE (dat->adj, int *);
        CC_IFFREE (dat->adjspace, int);
        fclose (datin);
        return 1;
    }
    for (i = 0, j = 0; i < *ncount; i++) {
        dat->adj[i] = dat->adjspace + j;
        j += (i+1);
    }
    for (i = 0; i < *ncount; i++) {
        for (j = 0; j <= i; j++) {
            dat->adj[i][j] = CCutil_readint (datin);
        }
    }
    
    fclose (datin);
    return 0;
}

static int write_matrix_text (char *datname, int ncount, CCdatagroup *dat)
{
    int i;
    int j;
    FILE *datout = fopen (datname, "w");
    
    if (datout == (FILE *) NULL) {
        perror (datname);
        fprintf (stderr, "Unable to open %s for output\n", datname);
        return 1;
    }
    printf ("%d\n", ncount);
    for (i = 0; i < ncount; i++) {
        for (j = 0; j <= i; j++) {
            printf ("%d ", dat->adj[i][j]);
            if (j%10 == 9) printf ("\n");
        }
        printf ("\n");
    }
    
    fclose (datout);
    return 0;
}

static int read_dsjrand (char *datname, int binary_in, int *ncount,
        CCdatagroup *dat)
{
    int i;
    
    dat->x = CC_SAFE_MALLOC (*ncount, double);
    if (!dat->x)
        return 1;
    for (i = 0; i < *ncount; i++)
        dat->x[i] = 0x12345672*(i+1) + 1;

    if (datname == (char *) NULL) {
        CCutil_dsjrand_init (dat, 1000000, 1);
    } else if (binary_in) {
        fprintf (stderr, "DSJRANDNORM doesn't support binary input\n");
        CC_IFFREE (dat->x, double);
        return 1;
    } else {
        int seed;
        int maxdist;
        FILE *datin = fopen (datname, "r");
        
        if (datin == (FILE *) NULL) {
            perror (datname);
            fprintf (stderr, "Unable to open %s for input\n", datname);
            return 1;
        }
        fscanf (datin, "%d", &seed);
        fscanf (datin, "%d", &maxdist);
        fclose (datin);
        CCutil_dsjrand_init (dat, 1000000, 1);
    }
    return 0;
}

static int read_sparse (char *datname, int binary_in, int *ncount,
        CCdatagroup *dat, CCrandstate *rstate)
{
    if (datname == (char *) NULL) {
        return build_sparse (ncount, dat, rstate);
    } else {
        return read_sparse_text (datname, ncount, dat, binary_in);
    }
}

static int write_sparse (char *datname, int binary_out, int ncount,
        CCdatagroup *dat)
{
    if (binary_out) {
        fprintf (stderr, "Binary output of sparse data not yet implemented\n");
        return 1;
    } else {
        return write_sparse_text (datname, ncount, dat);
    }
}

static int build_sparse (int *ncount, CCdatagroup *dat, CCrandstate *rstate)
{
    int ecount;
    int *elist = (int *) NULL;
    int *elen = (int *) NULL;
    int rval = 0;
    
    ecount = SPARSE_ECOUNT * (*ncount);
    printf ("Random graph with %d nodes, %d edges, and random edge lengths\n",
            *ncount, ecount);
    fflush (stdout);
    rval = CCutil_genedgelist (*ncount, ecount, &elist, &elen,
                     (CCdatagroup *) NULL, MATRAND_SCALE * (*ncount), rstate);
    if (rval) {
        fprintf (stderr, "CCutil_genedgelist failed\n"); goto CLEANUP;
    }
    rval = build_sparse_dat (*ncount, ecount, elist, elen, dat, 0);
    if (rval) {
        fprintf (stderr, "build_sparse_dat failed\n"); goto CLEANUP;
    }

 CLEANUP:
    CC_IFFREE (elist, int);
    CC_IFFREE (elen, int);
    return rval;
}

static int read_sparse_text (char *datname, int *ncount, CCdatagroup *dat,
        int binary_in)
{
    int ecount;
    int *elist = (int *) NULL;
    int *elen = (int *) NULL;
    int rval = 0;

    rval = CCutil_getedgelist_n (ncount, datname, &ecount, &elist,
                                 &elen, binary_in);
    if (rval) {
        fprintf (stderr, "CCutil_getedgelist failed\n"); goto CLEANUP;
    }
    printf ("Have a sparse graph with %d nodes and %d edges\n",
            *ncount, ecount);
    fflush (stdout);
    rval = build_sparse_dat (*ncount, ecount, elist, elen, dat, 0);
    if (rval) {
        fprintf (stderr, "build_sparse_dat failed\n"); goto CLEANUP;
    }
    
 CLEANUP:
    CC_IFFREE (elist, int);
    CC_IFFREE (elen, int);
    return rval;
}

static int write_sparse_text (char *datname, int ncount, CCdatagroup *dat)
{
    int i;
    int j;
    int ecount;
    FILE *datout = fopen (datname, "w");

    if (datout == (FILE *) NULL) {
        perror (datname);
        fprintf (stderr, "Unable to open %s for output\n", datname);
        return 1;
    }

    ecount=0;
    for (i=0; i<ncount; i++) {
        ecount += dat->degree[i];
    }
    fprintf (datout, "%d %d\n", ncount, ecount);
    for (i=0; i<ncount; i++) {
        for (j=0; j<dat->degree[i]; j++) {
            fprintf (datout, "%d %d %d\n", i, dat->adj[i][j], dat->len[i][j]);
        }
    }

    fclose (datout);
    return 0;
}

static int read_user (char *datname, int binary_in, int *ncount,
        CCdata_user *userdat, CCrandstate *rstate, int gridsize,
        int allow_dups)
{
    if (datname == (char *) NULL) {
        return build_user (ncount, userdat, rstate, gridsize, allow_dups);
    } else if (binary_in) {
        return read_user_binary (datname, ncount, userdat);
    } else {
        return read_user_text (datname, ncount, userdat);
    }
}

/****************************************************************************/
/*  build_user builds random data for the user-defined norm for computing   */
/*  edge lengths.  On input, *ncount is set to the number of nodes          */
/*  desired.  If necessary, the number of nodes generated can be different  */
/*  from the number of nodes desired (in which case, on return, *ncount     */
/*  should be set to the number of actual nodes generated).  On return,     */
/*  build_user should have set userdat to contain the information           */
/*  necessary for computing edge lengths.                                   */
/*                                                                          */
/*  build_user should return 0 if successful, nonzero if unsuccessful.      */
/*  If building a random data set is not supported for the user norm,       */
/*  build_user should output an error message and return 1.                 */
/****************************************************************************/
static int build_user (int *ncount, CCdata_user *userdat, CCrandstate *rstate,
        int gridsize, CC_UNUSED int allow_dups)
{
    int i;
    
    userdat->x = (double *) NULL;
    userdat->y = (double *) NULL;

    printf ("Random %d point set\n", *ncount); fflush (stdout);
    userdat->x = CC_SAFE_MALLOC (*ncount, double);
    userdat->y = CC_SAFE_MALLOC (*ncount, double);
    if (userdat->x == (double *) NULL ||
        userdat->y == (double *) NULL) {
        goto FAILURE;
    }
    for (i=0; i<*ncount; i++) {
        userdat->x[i] = (double) (CCutil_lprand (rstate) % gridsize);
        userdat->y[i] = (double) (CCutil_lprand (rstate) % gridsize);
    }
    return 0;

 FAILURE:
    CC_IFFREE (userdat->x, double);
    CC_IFFREE (userdat->y, double);
    return 1;
}

/****************************************************************************/
/*  read_user_binary reads a binary version of the data for the user-       */
/*  defined norm for computing edge lengths.  On return, read_user_binary   */
/*  should set *ncount to be the number of nodes read, and userdat to       */
/*  contain the information necessary for computing edge lengths.           */
/*                                                                          */
/*  read_user_binary should return 0 if successful, nonzero if              */
/*  unsuccessful.  If reading binary data is not supported for the user     */
/*  norm, read_user_binary should output an error message and return 1.     */
/****************************************************************************/
static int read_user_binary (char *datname, int *ncount, CCdata_user *userdat)
{
    int i;
    CC_SFILE *f = (CC_SFILE *) NULL;
    int xi;
    int yi;
    
    userdat->x = (double *) NULL;
    userdat->y = (double *) NULL;

    f = CCutil_sopen (datname, "r");
    if (f == (CC_SFILE *) NULL) goto FAILURE;
    if (CCutil_sread_int (f, ncount)) goto FAILURE;
    printf ("nnodes = %d\n", *ncount); fflush (stdout);
    userdat->x = CC_SAFE_MALLOC (*ncount, double);
    userdat->y = CC_SAFE_MALLOC (*ncount, double);
    if (userdat->x == (double *) NULL ||
        userdat->y == (double *) NULL) {
        goto FAILURE;
    }
    for (i=0; i<*ncount; i++) {
        if (CCutil_sread_int (f, &xi) ||
            CCutil_sread_int (f, &yi)) {
            goto FAILURE;
        }
        userdat->x[i] = (double) xi;
        userdat->y[i] = (double) yi;
    }
    if (CCutil_sclose (f)) goto FAILURE;
    f = (CC_SFILE *) NULL;
    return 0;

 FAILURE:
    CC_IFFREE (userdat->x, double);
    CC_IFFREE (userdat->y, double);
    if (f != (CC_SFILE *) NULL) CCutil_sclose (f);
    return 1;
}

/****************************************************************************/
/*  read_user_text reads a text version of the data for the user-defined    */
/*  norm for computing edge lengths.  On return, read_user_text should      */
/*  set *ncount to be the number of nodes, and userdat to contain the       */
/*  information necessary for computing edge lengths.                       */
/*                                                                          */
/*  read_user_text should return 0 if successful, nonzero if unsuccessful.  */
/****************************************************************************/
static int read_user_text (char *datname, int *ncount, CCdata_user *userdat)
{
    int i;
    FILE *datin = (FILE *) NULL;
    
    userdat->x = (double *) NULL;
    userdat->y = (double *) NULL;

    datin = fopen (datname, "r");
    if (datin == (FILE *) NULL) {
        perror (datname);
        fprintf (stderr, "Unable to open %s for input\n", datname);
        goto FAILURE;
    }
    fscanf (datin, "%d", ncount);
    printf ("nnodes = %d\n", *ncount); fflush (stdout);
    userdat->x = CC_SAFE_MALLOC (*ncount, double);
    userdat->y = CC_SAFE_MALLOC (*ncount, double);
    if (userdat->x == (double *) NULL ||
        userdat->y == (double *) NULL) {
        goto FAILURE;
    }
    for (i=0; i<*ncount; i++) {
        fscanf (datin, "%lf%lf", &(userdat->x[i]), &(userdat->y[i]));
    }
    fclose (datin);
    datin = (FILE *) NULL;
    return 0;
    
 FAILURE:
    CC_IFFREE (userdat->x, double);
    CC_IFFREE (userdat->y, double);
    if (datin != (FILE *) NULL) fclose (datin);
    return 1;
}

/****************************************************************************/
/*  writemaster_user writes the data in userdat in binary form to out.      */
/*                                                                          */
/*  writemaster_user should return 0 if successful, nonzero if unsuccessful*/
/****************************************************************************/
static int writemaster_user (CC_SFILE *out, int ncount, CCdata_user *userdat)
{
    int i;
    
    for (i=0; i<ncount; i++) {
        if (CCutil_swrite_double (out, userdat->x[i]) ||
            CCutil_swrite_double (out, userdat->y[i])) {
            return 1;
        }
    }
    return 0;
}

/****************************************************************************/
/*  readmaster_user reads the data from in into userdat.  The format read   */
/*  should be identical to that written by writemaster_user.                */
/*                                                                          */
/*  readmaster_user should return 0 if successful, nonzero if unsuccessful.*/
/****************************************************************************/
static int readmaster_user (CC_SFILE *in, int ncount, CCdata_user *userdat)
{
    int i;

    userdat->x = (double *) NULL;
    userdat->y = (double *) NULL;

    userdat->x = CC_SAFE_MALLOC (ncount, double);
    userdat->y = CC_SAFE_MALLOC (ncount, double);
    if (userdat->x == (double *) NULL ||
        userdat->y == (double *) NULL) {
        fprintf (stderr, "out of memory in readmaster_user\n");
        goto FAILURE;
    }
    for (i=0; i<ncount; i++) {
        if (CCutil_sread_double (in, &userdat->x[i]) ||
            CCutil_sread_double (in, &userdat->y[i])) {
            goto FAILURE;
        }
    }
    return 0;
 FAILURE:
    CC_IFFREE (userdat->x, double);
    CC_IFFREE (userdat->y, double);
    return 1;
}

/****************************************************************************/
/*  permute_user applies the permutation in perm to the data in userdat.    */
/*                                                                          */
/*  permute_user should return 0 if successful, nonzero if unsuccessful.    */
/****************************************************************************/
static int permute_user (int ncount, CCdata_user *userdat, int *perm)
{
    double *temp = (double *) NULL;
    int i;

    temp = CC_SAFE_MALLOC (ncount, double);
    if (temp == (double *) NULL) return 1;
    for (i=0; i<ncount; i++) {
        temp[i] = userdat->x[perm[i]];
    }
    CC_IFFREE (userdat->x, double);
    userdat->x = temp;

    temp = CC_SAFE_MALLOC (ncount, double);
    if (temp == (double *) NULL) return 1;
    for (i=0; i<ncount; i++) {
        temp[i] = userdat->y[perm[i]];
    }
    CC_IFFREE (userdat->y, double);
    userdat->y = temp;

    return 0;
}

static int read_rhdata (char *datname, int innorm, int binary_in, int *ncount,
        CCdata_rhvector *dat)
{
    if (datname == (char *) NULL) {
        fprintf (stderr, "generating random rhdata not supported\n");
        return 1;
    } else if (binary_in) {
        fprintf (stderr, "reading binary rhdata not supported\n");
        return 1;
    } else {
        return read_rhdata_text (datname, innorm, ncount, dat);
    }
}

static int read_rhdata_text (char *datname, int innorm, int *ncount,
        CCdata_rhvector *dat)
{
    int i;
    int j;
    FILE *datin = (FILE *) NULL;
    int rhlength;
    char *v;
    int xj;
    int nvectors;
    int n0;
    int n1;
    
    dat->space = (char *) NULL;
    dat->vectors = (char **) NULL;

    datin = fopen (datname, "r");
    if (datin == (FILE *) NULL) {
        perror (datname);
        fprintf (stderr, "Unable to open %s for input\n", datname);
        goto FAILURE;
    }
    fscanf (datin, "%d%d", &nvectors, &rhlength);
    *ncount = nvectors+1;
    printf ("nnodes = %d, rhlength %d, ", *ncount, rhlength);
    switch (innorm) {
    case CC_RHMAP1: printf ("rh norm 1\n"); break;
    case CC_RHMAP2: printf ("rh norm 2\n"); break;
    case CC_RHMAP3: printf ("rh norm 3\n"); break;
    case CC_RHMAP4: printf ("rh norm 4\n"); break;
    case CC_RHMAP5: printf ("rh norm 5\n"); break;
    default:        printf ("UNKNOWN rh norm\n"); break;
    }        
    fflush (stdout);
    dat->space = CC_SAFE_MALLOC (nvectors * rhlength, char);
    dat->vectors = CC_SAFE_MALLOC (nvectors+1, char *);
    if (dat->space == (char *) NULL ||
        dat->vectors == (char **) NULL) {
        goto FAILURE;
    }
    fscanf (datin, "%d", &dat->dist_00);
    fscanf (datin, "%d", &dat->dist_01);
    fscanf (datin, "%d", &dat->dist_02);
    fscanf (datin, "%d", &dat->dist_12);
    fscanf (datin, "%d", &dat->dist_22);

    dat->vectors[0] = (char *) NULL;
    for (i=0; i<nvectors; i++) {
        fscanf (datin, "%*d %*s "); /* vector labels */
        v = dat->space + i * rhlength;
        dat->vectors[i+1] = v;
        for (j=0; j<rhlength; j++) {
            xj = getc(datin);
            if (xj < '0' || xj > '2') {
                fprintf (stderr, "Syntax error in text rhvector %d data - entry (octal %o) %c\n", i, xj, xj);
                goto FAILURE;
            }
            v[j] = xj - '0';
        }
    }
    dat->rhlength = rhlength;
    fclose (datin);
    datin = (FILE *) NULL;

    n0 = 0;
    n1 = 0;
    for (i=0; i<*ncount; i++) {
        v = dat->vectors[i];
        if (v != (char *) NULL) {
            for (j=0; j<rhlength; j++) {
                if (v[j] == 0) n0++;
                else if (v[j] == 1) n1++;
            }
        }
    }
    if (n0 + n1 > 0) {
        dat->p = ((double) n1) / ((double) (n0 + n1));
    } else {
        dat->p = 0.0;
    }
            
    return 0;
 FAILURE:
    CC_IFFREE (dat->vectors, char *);
    CC_IFFREE (dat->space, char);
    if (datin != (FILE *) NULL) fclose (datin);
    return 1;
}

static int writemaster_rhdata (CC_SFILE *out, int ncount, CCdata_rhvector *dat)
{
    int i;
    int j;
    int rhlength = dat->rhlength;
    char *v;

    CCutil_swrite_int (out, rhlength);
    CCutil_swrite_int (out, dat->dist_00);
    CCutil_swrite_int (out, dat->dist_01);
    CCutil_swrite_int (out, dat->dist_02);
    CCutil_swrite_int (out, dat->dist_12);
    CCutil_swrite_int (out, dat->dist_22);
    for (i=0; i<ncount; i++) {
        v = dat->vectors[i];
        if (v == (char *) NULL) {
            if (CCutil_swrite_bits (out, 3, 2)) {
                return 1;
            }
        } else {
            for (j=0; j<rhlength; j++) {
                if (CCutil_swrite_bits (out, v[j], 2)) {
                    return 1;
                }
            }
        }
    }
    if (CCutil_sflush (out)) return 1;
    return 0;
}

static int readmaster_rhdata (CC_SFILE *in, int ncount, CCdata_rhvector *dat)
{
    int i;
    int j;
    int rhlength;
    char *v;
    int xj;
    int n0;
    int n1;
    
    dat->space = (char *) NULL;
    dat->vectors = (char **) NULL;

    if (CCutil_sread_int (in, &rhlength)) goto FAILURE;
    dat->rhlength = rhlength;
    
    dat->space = CC_SAFE_MALLOC (ncount * rhlength, char);
    dat->vectors = CC_SAFE_MALLOC (ncount, char *);
    if (dat->space == (char *) NULL ||
        dat->vectors == (char **) NULL) goto FAILURE;

    CCutil_sread_int (in, &dat->dist_00);
    CCutil_sread_int (in, &dat->dist_01);
    CCutil_sread_int (in, &dat->dist_02);
    CCutil_sread_int (in, &dat->dist_12);
    CCutil_sread_int (in, &dat->dist_22);

    v = dat->space;
    for (i=0; i<ncount; i++) {
        if (CCutil_sread_bits (in, &xj, 2)) {
            goto FAILURE;
        }
        if (xj == 3) {
            dat->vectors[i] = (char *) NULL;
        } else {
            dat->vectors[i] = v;
            v[0] = xj;
            for (j=1; j<rhlength; j++) {
                if (CCutil_sread_bits (in, &xj, 2)) {
                    goto FAILURE;
                }
                v[j] = xj;
            }
            v += rhlength;
        }
    }
    if (CCutil_sflush (in)) return 1;

    n0 = 0;
    n1 = 0;
    for (i=0; i<ncount; i++) {
        v = dat->vectors[i];
        if (v != (char *) NULL) {
            for (j=0; j<rhlength; j++) {
                if (v[j] == 0) n0++;
                else if (v[j] == 1) n1++;
            }
        }
    }
    if (n0 + n1 > 0) {
        dat->p = ((double) n1) / ((double) (n0 + n1));
    } else {
        dat->p = 0.0;
    }
            
    return 0;
 FAILURE:
    CC_IFFREE (dat->vectors, char *);
    CC_IFFREE (dat->space, char);
    return 1;
}

static int permute_rhdata (int ncount, CCdata_rhvector *dat, int *perm)
{
    char **temp = (char **) NULL;
    int i;

    temp = CC_SAFE_MALLOC (ncount, char *);
    if (temp == (char **) NULL) return 1;
    for (i=0; i<ncount; i++) {
        temp[i] = dat->vectors[perm[i]];
    }
    CC_IFFREE (dat->vectors, char *);
    dat->vectors = temp;

    return 0;
}
