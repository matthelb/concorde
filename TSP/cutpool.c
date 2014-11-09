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
/*                  STORING AND SEARCHING THE CUTPOOL                       */
/*                                                                          */
/*                            TSP CODE                                      */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: March 19, 1997                                                    */
/*        May 27, 1997 (bico)                                               */
/*        January 30, 2003 (bico)                                           */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCtsp_init_cutpool (int *ncount, char *poolfilename,                */
/*      CCtsp_lpcuts **pool)                                                */
/*        -ncount is a pointer to the number of nodes in the problem        */
/*        -poolfilename is a file containing an cutpool (it can be NULL)    */
/*        -CCtsp_lpcuts will return the pool                                */
/*        NOTES: poolfilename must be non-NULL or ncount must be            */
/*        non-NULL and *ncount nonzero.  If ncount is non-NULL but          */
/*        *ncount == zero, then *ncount will be set to the number of        */
/*        nodes in the cutpool in poolfilename                              */
/*    NOTES:                                                                */
/*        This version does not use the compressed set references.  Notes   */
/*    on the representation are given in "Chapter 4: The Linear             */
/*    Programming Problems".                                                */
/*                                                                          */
/*  int CCtsp_search_cutpool (CCtsp_lpcuts *pool, CCtsp_lpcut_in **cuts,    */
/*      int *cutcount, double *maxviol, int ncount, int ecount,             */
/*      int *elist, double *x, int nthreads, CCrandstate *rstate)           */
/*    RETURNS an array of cuts having x(delta(C)) < rhs(C)                  */
/*     -pool points to a cutpool (or cuts of an lp)                         */
/*     -cuts will return the array of cuts                                  */
/*     -cutcount with return the length of the array                        */
/*     -ncount is the number of nodes in the problem                        */
/*     -ecount is the number of edges in elist                              */
/*     -elist is a list of edges in end end format                          */
/*     -x is an ecount-long array of weights                                */
/*     -nthreads is the number of threads to use.  0 ==> sequential code    */
/*      threads are only used if CC_POSIXTHREADS is defined                 */
/*                                                                          */
/*  int CCtsp_search_remotepool (char *remotehost,                          */
/*      unsigned short remoteport, CCtsp_lpcut_in **cuts,                   */
/*      int *cutcount, double *maxviol, int ncount, int ecount,             */
/*      int *elist, double *x)                                              */
/*    RETURNS an array of cuts having x(delta(C)) < rhs(C) from a remote    */
/*      cutpool                                                             */
/*     -remotehost is the host with the cuts                                */
/*     -remoteport is the port on which to contact the remote server        */
/*      (if remoteport == 0, use CCtsp_CUT_PORT)                            */
/*     -cuts will return the array of cuts                                  */
/*     -cutcount with return the length of the array                        */
/*     -ncount is the number of nodes in the problem                        */
/*     -ecount is the number of edges in elist                              */
/*     -elist is a list of edges in end end format                          */
/*     -x is an ecount-long array of weights                                */
/*                                                                          */
/*  int CCtsp_search_cutpool_cliques (CCtsp_lpcuts *pool,                   */
/*      CCtsp_lpclique **cliques, int *cliquecount, int ncount,             */
/*      int ecount, int *elist, double *x, double maxdelta,                 */
/*      int maxcliques, double **cliquevals, CCrandstate *rstate)           */
/*    RETURNS an array of cliques having x(delta(C)) < maxdelta             */
/*     -pool points to a cutpool (or cuts of an lp)                         */
/*     -cliques will return the array of cliques                            */
/*     -cliquecount with return the length of the array                     */
/*     -ncount is the number of nodes in the problem                        */
/*     -ecount is the number of edges in elist                              */
/*     -elist is a list of edges in end end format                          */
/*     -x is an ecount-long array of weights                                */
/*     -maxdelta is a bound on x(delta(C))                                  */
/*     -maxcliques is an upperbound on the number of cliques that should    */
/*      be returned                                                         */
/*     -cliquevals will return the values of x(delta(C)) for the cliques    */
/*      (this parameter can be NULL)                                        */
/*                                                                          */
/*  int CCtsp_add_cut_to_cutlist (CCtsp_lpcuts *cuts, CCtsp_lpcut *c)       */
/*    NONE                                                                  */
/*                                                                          */
/*  void CCtsp_delete_cut_from_cutlist (CCtsp_lpcuts *cuts, int ind)        */
/*    NONE                                                                  */
/*                                                                          */
/*  void CCtsp_free_cutpool (CCtsp_lpcuts **pool)                           */
/*    FREES the pool of cuts.                                               */
/*                                                                          */
/*  int CCtsp_write_cutpool (int ncount, const char *poolfilename,          */
/*      CCtsp_lpcuts *pool)                                                 */
/*    WRITES pool to poolfilename.                                          */
/*                                                                          */
/*  int CCtsp_branch_cutpool_cliques (CCtsp_lpcuts *pool,                   */
/*      CCtsp_lpclique **cliques, int *cliquecount, int ncount,             */
/*      int ecount, int *elist, double *x, int nwant,                       */
/*      double **cliquevals, int silent)                                    */
/*    RETURNS an array of cliques having x(delta(C)) as close to 3.0 as     */
/*     possible.                                                            */
/*     -the parmeters are like those used by search_cutpool_cliques,        */
/*      where nwant is the number of cliques we would like to have in       */
/*      the array.                                                          */
/*                                                                          */
/*  int CCtsp_price_cuts (CCtsp_lpcuts *pool, int ncount, int ecount,       */
/*      int *elist, double *x, double *cutval)                              */
/*    COMPUTES the slack on each cut in the pool                            */
/*     -ecount, elist, and x give an x-vector                               */
/*     -cutval returns the array of slack values (it should be passed in    */
/*      as an array of length at least pool->cutcount)                      */
/*                                                                          */
/*  int CCtsp_price_cuts_threaded (CCtsp_lpcuts *pool, int ncount,          */
/*      int ecount, int *elist, double *x, double *cutval,                  */
/*      int numthreads)                                                     */
/*    COMPUTES the slack on each cut in the pool in parallel                */
/*     -ecount, elist, and x give an x-vector                               */
/*     -cutval returns the array of slack values (it should be passed in    */
/*      as an array of length at least pool->cutcount)                      */
/*     -nthreads is the number of parallel threads to use.                  */
/*                                                                          */
/*  int CCtsp_get_clique_prices (CCtsp_lpcuts *pool, int **p_cliquenums,    */
/*      double **p_cliquevals, double mindelta, double maxdelta,            */
/*      int *p_cliquecount, int ncount, int ecount, int *elist,             */
/*      double *x)                                                          */
/*    RETURNS the id's and x(delta(C)) for cliques in the pool.             */
/*     -the parameters pool, ncount, ecount, elist, and x are like those    */
/*      used by search_cutpool_cliques.                                     */
/*     -mindelta and maxdelta are bounds on x(delta(C))                     */
/*     -cliquenums and cliquevals return id's and x(delta(C)) for           */
/*      cliques with x(delta(C)) between mindelta and maxdelta              */
/*     -cliquecount returns the number of cliques in cliquenums/vals        */
/*     -use CCtsp_get_clique to retrive specific cliques.                   */
/*                                                                          */
/*  int CCtsp_get_clique (CCtsp_lpcuts *pool, int cliquenum,                */
/*      CCtsp_lpclique **p_clique)                                          */
/*    RETURNS p_clique, a pointer to the clique numbered clqiuenum in       */
/*      the pool.  Note that this clique is not a copy, and thus should     */
/*      not be freed.                                                       */
/*                                                                          */
/*  int CCtsp_display_cutpool (CCtsp_lpcuts *pool)                          */
/*    DISPLAYS the contents of a cutpool.                                   */
/*                                                                          */
/*  int CCtsp_add_to_cutpool (CCtsp_lpcuts *pool, CCtsp_lpcuts *cuts,       */
/*      CCtsp_lpcut *c)                                                     */
/*     -pool is the pool to add the cut to                                  */
/*     -cuts is the lpcuts the cut is from                                  */
/*     -c is the cut                                                        */
/*    ADDS a cut to a pool                                                  */
/*                                                                          */
/*  int CCtsp_add_to_cutpool_lpcut_in (CCtsp_lpcuts *pool,                  */
/*      CCtsp_lpcut_in *c)                                                  */
/*     -pool is the pool to add the cut to                                  */
/*     -c is the cut                                                        */
/*    ADDS a cut to a pool                                                  */
/*                                                                          */
/*  void CCtsp_free_lpcut_in (CCtsp_lpcut_in *c)                            */
/*    FREES the fields in the CCtsp_lpcut pointed to by c.                  */
/*                                                                          */
/*  void CCtsp_free_lpclique (CCtsp_lpclique *c)                            */
/*    FREES the fields in the CCtsp_lpclique pointed to by c.               */
/*                                                                          */
/*  void CCtsp_free_lpdomino (CCtsp_lpdomino *c)                            */
/*    FREES the fields in the CCtsp_lpdomino pointed to by c.               */
/*                                                                          */
/*  int CCtsp_read_cuts (CC_SFILE *f, int *ncount, CCtsp_lpcuts *cuts,      */
/*      int readmods, int buildhash)                                        */
/*    READS the cuts from f into cuts.                                      */
/*    -readmods indicates whether or not the file contains sparser mods     */
/*                                                                          */
/*  int CCtsp_read_lpcut_in (CC_SFILE *f, CCtsp_lpcut_in *c, int ncount)    */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCtsp_read_lpclique (CC_SFILE *f, CCtsp_lpclique *c, int ncount)    */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCtsp_read_lpdomino (CC_SFILE *f, CCtsp_lpdomino *d, int ncount)    */
/*    READS the lpdomino from the file stream.                              */
/*                                                                          */
/*  int CCtsp_send_newcuts (int ncount, CCtsp_lpcuts *pool,                 */
/*      char *remotehost, unsigned short remoteport)                        */
/*    SENDS the new cuts from pool to the remote host                       */
/*                                                                          */
/*  int CCtsp_write_cuts (CC_SFILE *f, int ncount, CCtsp_lpcuts *cuts,      */
/*      int writemods)                                                      */
/*    WRITES the cuts from cuts to f.                                       */
/*    -writemods indicates whether or not the file should contain           */
/*     sparser mods                                                         */
/*                                                                          */
/*  int CCtsp_write_lpcut_in (CC_SFILE *f, CCtsp_lpcut_in *c, int ncount)   */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCtsp_write_lpcut (CC_SFILE *f, CCtsp_lpcuts *cuts,                 */
/*      CCtsp_lpcut *c, int ncount)                                         */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCtsp_write_lpclique (CC_SFILE *f, CCtsp_lpclique *c, int ncount)   */
/*    MISSING                                                               */
/*                                                                          */
/*  int CCtsp_copy_cuts (CC_SFILE *f, CC_SFILE *t, int copymods)            */
/*    COPIES the cuts from f to t.                                          */
/*                                                                          */
/*  int CCtsp_register_cliques (CCtsp_lpcuts *cuts, CCtsp_lpcut_in *c,      */
/*      CCtsp_lpcut *new)                                                   */
/*    BUILDS the references to the cliques in c into the cut strucure       */
/*    pointed to by cuts and creates an array of the indices of the         */
/*    the cliques in CCtsp_lpcut new                                        */
/*     -cuts is the structure holding the set of cuts                       */
/*     -c describes the cut to be added to the structure                    */
/*     -new returns the array of clique indices                             */
/*                                                                          */
/*  void CCtsp_unregister_cliques (CCtsp_lpcuts *cuts, CCtsp_lpcut *c)      */
/*    REMOVES the references to the cliques in cut c (and deletes the       */
/*     cliques if they have no more references) and frees the array         */
/*     of clique indices in c                                               */
/*     -cuts is the structure holding the set of cuts                       */
/*     -c is the cut containing the cliques to be removed                   */
/*                                                                          */
/*  int CCtsp_register_dominos (CCtsp_lpcuts *cuts, CCtsp_lpcut_in *c,      */
/*      CCtsp_lpcut *new)                                                   */
/*    BUILDS the references to the dominos in c into the cut strucure       */
/*    pointed to by cuts and creates an array of the indices of the         */
/*    the dominos in CCtsp_lpcut new                                        */
/*                                                                          */
/*  void CCtsp_unregister_dominos (CCtsp_lpcuts *cuts, CCtsp_lpcut *c)      */
/*    REMOVES the references to the dominos in cut c (and deletes the       */
/*     dominos if they have no more references) and frees the array         */
/*     of domino indices in c                                               */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "macrorus.h"
#include "tsp.h"

#define ZERO_EPSILON 0.0000000001
#define POOL_MAXCUTS 500
#define POOL_MINVIOL 0.001

#define PROB_CUTS_VERSION 2   /* Version 1 is pre-dominos */

typedef struct pooledge {
    double x;
    int to;
} pooledge;

typedef struct poolnode {
    struct pooledge *adj;
    int mark;
    int deg;
} poolnode;


static int
    init_empty_cutpool_hash (int ncount, CCtsp_lpcuts *pool),
    cut_eq (void *v_cut1, void *v_cut2, void *u_data),
    read_cutpool (int *ncount, char *poolfilename, CCtsp_lpcuts *pool),
    register_lpcuts (CCtsp_lpcuts *pool),
    price_cliques (CCtsp_lpclique *cliques, int ncount, int ecount, int *elist,
            double *x, double *cval, int cend),
    make_pricing_graph (int ncount, int ecount, int *elist, double *x,
            poolnode **p_nlist, pooledge **p_espace);

static unsigned int
    cut_hash (void *v_cut, void *u_data);

static void
#ifdef CC_POSIXTHREADS
   *price_cliques_thread (void *args),
   *price_cuts_thread (void *args),
    rebalance_load (int nthreads, double *workload, double *worktime,
        double *work),
#endif
    price_cuts (CCtsp_lpcut *cuts, int cutcount, double *cval,
        double *cutval),
    sort_cliques (CCtsp_lpcut *c),
    sort_dominos (CCtsp_lpcut *c);

static double
    price_clique (poolnode *nlist, CCtsp_lpclique *c, int marker);



int CCtsp_init_cutpool (int *ncount, char *poolfilename, CCtsp_lpcuts **pool)
{
    int rval = 0;
    CCtsp_lpcuts *p = (CCtsp_lpcuts *) NULL;

    p = CC_SAFE_MALLOC (1, CCtsp_lpcuts);
    if (!p) {
        fprintf (stderr, "out of memory in CCtsp_init_cutpool\n");
        return 1;
    }
    *pool = p;

    p->cutcount    = 0;
    p->savecount   = 0;
    p->cuts        = (CCtsp_lpcut *) NULL;
    p->cutspace    = 0;
    p->cliqueend   = 0;
    p->cliques     = (CCtsp_lpclique *) NULL;
    p->cliquespace = 0;
    p->cliquehash  = (int *) NULL;
    p->dominoend   = 0;
    p->dominos     = (CCtsp_lpdomino *) NULL;
    p->dominospace = 0;
    p->dominohash  = (int *) NULL;
    p->cuthash     = (CCgenhash *) NULL;
    p->workloads   = (double *) NULL;

    if (poolfilename == (char *) NULL) {
        if (ncount == (int *) NULL || *ncount <= 0) {
            fprintf (stderr, "Neither poolfilename nor ncount\n");
            rval = 1; goto CLEANUP;
        }
        rval = init_empty_cutpool_hash (*ncount, p);
        if (rval) {
            fprintf (stderr, "init_empty_cutpool_hash failed\n"); goto CLEANUP;
        }
    } else {
        rval = read_cutpool (ncount, poolfilename, p);
        if (rval) {
            fprintf (stderr, "read_cutpool failed\n"); goto CLEANUP;
        }
    }

CLEANUP:
    return rval;
}

void CCtsp_free_cutpool (CCtsp_lpcuts **pool)
{
    int i, k;

    if (*pool) {
        if ((*pool)->cuts) {
            for (i = 0; i < (*pool)->cutcount; i++) {
                CC_IFFREE ((*pool)->cuts[i].cliques, int);
                CCtsp_free_skeleton (&(*pool)->cuts[i].skel);
            }
            CC_FREE ((*pool)->cuts, CCtsp_lpcut);
        }
        if ((*pool)->cliques) {
            for (i=0; i < (*pool)->cliqueend; i++) {
                CC_IFFREE ((*pool)->cliques[i].nodes, CCtsp_segment);
            }
            CC_FREE ((*pool)->cliques, CCtsp_lpclique);
        }
        if ((*pool)->dominos) {
            for (i=0; i < (*pool)->dominoend; i++) {
                for (k = 0; k < 2; k++) {
                    CC_IFFREE ((*pool)->dominos[i].sets[k].nodes,
                                CCtsp_segment);
                }
            }
            CC_FREE ((*pool)->dominos, CCtsp_lpdomino);
        }

        CCtsp_free_cliquehash (*pool);
        CCtsp_free_dominohash (*pool);

        if ((*pool)->cuthash) {
           CCutil_genhash_free ((*pool)->cuthash, NULL);
           CC_FREE ((*pool)->cuthash, CCgenhash);
        }
        CC_IFFREE ((*pool)->workloads, double);
        CC_FREE (*pool, CCtsp_lpcuts);
    }
}

static int init_empty_cutpool_hash (int ncount, CCtsp_lpcuts *pool)
{
    int rval = 0;

    rval = CCtsp_init_cliquehash (pool, 10 * ncount);
    CCcheck_rval (rval, "CCtsp_init_cliquehash failed");

    rval = CCtsp_init_dominohash (pool, 10 * ncount);
    CCcheck_rval (rval, "CCtsp_init_dominohash failed");

    pool->cuthash = CC_SAFE_MALLOC (1, CCgenhash);
    CCcheck_NULL (pool->cuthash, "out of memory in init_empty_cutpool_hash");

    rval = CCutil_genhash_init (pool->cuthash, 10 * ncount, cut_eq,
                         cut_hash, (void *) pool, 1.0, 0.6);
    CCcheck_rval (rval, "CCutil_genhash_init failed");

CLEANUP:

    return rval;
}

static int cut_eq (void *v_cut1, void *v_cut2, void *u_data)
{
    CCtsp_lpcuts *pool = (CCtsp_lpcuts *) u_data;
    CCtsp_lpcut *cut1 = pool->cuts + (long) v_cut1;
    CCtsp_lpcut *cut2 = pool->cuts + (long) v_cut2;
    int i;

    if (cut1->cliquecount != cut2->cliquecount) return 1;
    if (cut1->dominocount != cut2->dominocount) return 1;
    if (cut1->rhs != cut2->rhs) return 1;
    if (cut1->sense != cut2->sense) return 1;
    for (i=0; i<cut1->cliquecount; i++) {
        if (cut1->cliques[i] != cut2->cliques[i]) return 1;
    }
    for (i=0; i<cut1->dominocount; i++) {
        if (cut1->dominos[i] != cut2->dominos[i]) return 1;
    }
#if 0
    {
        int diff = 0;
        CCtsp_compare_skeletons (&cut1->skel, &cut2->skel, &diff);
        if (diff) {
          printf ("SURPRISE - cuts look equal but have different skeletons\n");
          return 1;
        }
    }
#endif
    return 0;
}

static unsigned int cut_hash (void *v_cut, void *u_data)
{
    CCtsp_lpcuts *pool = (CCtsp_lpcuts *) u_data;
    CCtsp_lpcut *cut = pool->cuts + (long) v_cut;
    unsigned int x = ((unsigned int) cut->rhs) * 257 +
                     ((unsigned int) cut->sense);
    int i;

    for (i=0; i<cut->cliquecount; i++) {
        x = x * 4099 + cut->cliques[i];
    }
    for (i=0; i<cut->dominocount; i++) {
        x = x * 4099 + cut->dominos[i];
    }
    return x;
}

static int read_cutpool (int *ncount, char *poolfilename, CCtsp_lpcuts *pool)
{
    CC_SFILE *in = (CC_SFILE *) NULL;
    int n;
    int rval = 0;

    if (poolfilename == (char *) NULL) {
        fprintf (stderr, "pool file name is not set\n");
        rval = 1; goto CLEANUP;
    }

    in = CCutil_sopen (poolfilename, "r");
    if (!in) {
        fprintf (stderr, "CCutil_sopen failed\n");
        rval = 1; goto CLEANUP;
    }

    rval = CCtsp_read_cuts (in, &n, pool, 0, 1);
    if (rval < 0) {
        fprintf (stderr, "CCtsp_read_cuts failed\n"); goto CLEANUP;
    }

    if (ncount != (int *) NULL && *ncount > 0 && n != *ncount) {
        fprintf (stderr, "cutpool %s does not have the correct ncount\n",
                            poolfilename);
        rval = 1; goto CLEANUP;
    }

    if (ncount != (int *) NULL) *ncount = n;

    rval = CCutil_sclose (in);
    in = (CC_SFILE *) NULL;
    if (rval) {
        fprintf (stderr, "CCutil_sclose failed\n"); goto CLEANUP;
    }

    rval = 0;
    
  CLEANUP:
    if (in != (CC_SFILE *) NULL) {
        CCutil_sclose (in);
    }
    
    return rval;
}

int CCtsp_read_cuts (CC_SFILE *f, int *ncount, CCtsp_lpcuts *cuts,
        int readmods, int buildhash)
{
    int i, j, k, ncliq, ndom, rhs, cliq, dom, nmod, n, mult;
    char sense;
    int cliqcount = 0, cutcount = 0, dominocount = 0;
    CCtsp_lpclique c;
    CCtsp_lpdomino d;
    CCtsp_lpcut u;
    int *hits = (int *) NULL;
    int *domhits = (int *) NULL;
    char version;
    int rval;
    int nbits, cbits, dbits;

    CCtsp_init_lpclique (&c);
    
    if (CCutil_sread_char (f, &version)) goto FAILURE;

    if (version != 1 && version != 2) {
        fprintf (stderr, "Unknown cuts version %d\n", (unsigned) version);
        goto FAILURE;
    } else {
        if (CCutil_sread_int (f, ncount)) goto FAILURE;
        nbits = CCutil_sbits (*ncount);

        if (buildhash) {
            rval = init_empty_cutpool_hash (*ncount, cuts);
            if (rval) {
                fprintf (stderr, "init_empty_cutpool_hash failed\n");
                goto FAILURE;
            }
        }

        if (CCutil_sread_int (f, &cliqcount)) goto FAILURE;
        for (i = 0; i < cliqcount; i++) {
            rval = CCtsp_read_lpclique (f, &c, *ncount);
            if (rval) {
                fprintf (stderr, "CCtsp_read_lpclique failed\n");
                goto FAILURE;
            }
            k = CCtsp_register_clique (cuts, &c);
            if (k == -1) {
                fprintf (stderr, "CCtsp_register_clique failed\n");
                goto FAILURE;
            }
            if (k != i) {
                fprintf (stderr, "clique registration number is out of seq\n");
                goto FAILURE;
            }
            CCtsp_free_lpclique (&c);
        }

        if (version == 1) {
            dominocount = 0;
        } else {
            if (CCutil_sread_int (f, &dominocount)) goto FAILURE;
        }
        for (i = 0; i < dominocount; i++) {
            rval = CCtsp_read_lpdomino (f, &d, *ncount);
            if (rval) {
                fprintf (stderr, "CCtsp_read_lpdomino failed\n");
                goto FAILURE;
            }
            k = CCtsp_register_domino (cuts, &d);
            if (k == -1) {
                fprintf (stderr, "CCtsp_register_domino failed\n");
                goto FAILURE;
            }
            if (k != i) {
                fprintf (stderr, "domino registration number is out of seq\n");
                goto FAILURE;
            }
            CCtsp_free_lpdomino (&d);
        }

        if (cliqcount) {
            hits = CC_SAFE_MALLOC (cliqcount, int);
            if (!hits) {
                fprintf (stderr, "out of memory in CCtsp_read_cuts\n");
                goto FAILURE;
            }
            for (i = 0; i < cliqcount; i++) hits[i] = 0;
        }
 
        if (dominocount) {
            domhits = CC_SAFE_MALLOC (dominocount, int);
            if (!domhits) {
                fprintf (stderr, "out of memory in CCtsp_read_cuts\n");
                goto FAILURE;
            }
            for (i = 0; i < dominocount; i++) domhits[i] = 0;
        }

        cbits = CCutil_sbits (cliqcount);
        dbits = CCutil_sbits (dominocount);
        if (CCutil_sread_int (f, &cutcount)) goto FAILURE;

        for (i=0; i<cutcount; i++) {
            CCtsp_init_lpcut (&u);
            if (CCutil_sread_int (f, &ncliq)) goto FAILURE;
            if (version == 1) {
                ndom = 0;
            } else {
                if (CCutil_sread_int (f, &ndom)) goto FAILURE;
            }
            if (CCutil_sread_int (f, &rhs)) goto FAILURE;
            if (CCutil_sread_char (f, &sense)) goto FAILURE;
            u.cliquecount = ncliq;
            u.dominocount = ndom;
            u.rhs = rhs;
            u.sense = sense;
            u.branch = 0;
            u.age = 0;
            u.cliques = CC_SAFE_MALLOC (ncliq, int);
            if (!u.cliques) {
                fprintf (stderr, "out of memory in CCtsp_read_cuts\n");
                goto FAILURE;
            }
            for (j = 0; j < ncliq; j++) {
                if (CCutil_sread_bits (f, &cliq, cbits))
                    goto FAILURE;
                u.cliques[j] = cliq;
                if (hits[cliq])
                    cuts->cliques[cliq].refcount++;
                else
                    hits[cliq] = 1;
            }
            if (ndom) {
                u.dominos = CC_SAFE_MALLOC (ndom, int);
                if (!u.dominos) {
                    fprintf (stderr, "out of memory in CCtsp_read_cuts\n");
                    goto FAILURE;
                }
                for (j = 0; j < ndom; j++) {
                    if (CCutil_sread_bits (f, &dom, dbits)) {
                        goto FAILURE;
                    }
                    u.dominos[j] = dom;
                    if (domhits[dom]) {
                        cuts->dominos[dom].refcount++;
                    } else {
                        domhits[dom] = 1;
                    }
                }
            } else {
                u.dominos = (int *) NULL;
            }
            if (readmods) {
                if (CCutil_sread_int (f, &nmod)) goto FAILURE;
                u.modcount = nmod;
                if (nmod) {
                    u.mods = CC_SAFE_MALLOC (nmod, CCtsp_sparser);
                    if (!u.mods) {
                        fprintf (stderr, "out of memory in CCtsp_read_cuts\n");
                        CC_FREE (u.cliques, int);
                        goto FAILURE;
                    }
                } else {
                    u.mods = (CCtsp_sparser *) NULL;
                }
                for (j = 0; j < nmod; j++) {
                    if (CCutil_sread_bits (f, &n, nbits))
                        goto FAILURE;
                    if (CCutil_sread_int (f, &mult))
                        goto FAILURE;
                    u.mods[j].node = (unsigned int) n;
                    u.mods[j].mult = (unsigned int) mult;
                }
            } else {
                u.mods = (CCtsp_sparser *) NULL;
                u.modcount = 0;
            }

            rval = CCtsp_read_skeleton (f, &u.skel, *ncount);
            if (rval) {
                fprintf (stderr, "CCtsp_read_skeleton failed\n");
                goto FAILURE;
            }

            k = CCtsp_add_cut_to_cutlist (cuts, &u);
            if (k == -1) {
                fprintf (stderr, "CCtsp_add_cut_to_cutlist failed\n");
                goto FAILURE;
            }
            if (k != i) {
                fprintf (stderr, "cut location is out of seq\n");
                goto FAILURE;
            }
        }

        CC_IFFREE (hits, int);
        CC_IFFREE (domhits, int);

        if (buildhash) {
            rval = register_lpcuts (cuts);
            if (rval) {
                fprintf (stderr, "register_lpcuts failed\n");
                goto FAILURE;
            }
        }
    }

    return 0;

FAILURE:

    CC_IFFREE (hits, int);
    CC_IFFREE (domhits, int);
    if (cutcount) {
        for (i = 0; i < cutcount; i++) {
            CCtsp_delete_cut_from_cutlist (cuts, i);
        }
    } else {
        for (i = 0; i < cliqcount; i++) {
            CCtsp_unregister_clique (cuts, i);
       }
        for (i = 0; i < dominocount; i++) {
            CCtsp_unregister_domino (cuts, i);
       }
    }
    return -1;
}

int CCtsp_read_lpcut_in (CC_SFILE *f, CCtsp_lpcut_in *c, int ncount)
{
    int rval = 0;
    int ncliques, ndominos;
    int i;
    int rhs;
    char sense;

    CCtsp_init_lpcut_in (c);
    
    rval = CCutil_sread_int (f, &ncliques);
    CCcheck_rval (rval, "CCutil_sread_int failed");

    rval = CCutil_sread_int (f, &ndominos);
    CCcheck_rval (rval, "CCutil_sread_int failed");

    c->cliquecount = ncliques;
    c->cliques = CC_SAFE_MALLOC (ncliques, CCtsp_lpclique);
    CCcheck_NULL (c->cliques, "out of memory in CCtsp_read_lpcut_in");
    for (i=0; i<ncliques; i++) {
        CCtsp_init_lpclique (&c->cliques[i]);
    }
    for (i = 0; i < ncliques; i++) {
        rval = CCtsp_read_lpclique (f, &c->cliques[i], ncount);
        CCcheck_rval (rval, "CCtsp_read_lpclique failed");
    }

    c->dominocount = ndominos;
    if (ndominos > 0) {
        c->dominos = CC_SAFE_MALLOC (ndominos, CCtsp_lpdomino);
        CCcheck_NULL (c->dominos, "out of memory in CCtsp_read_lpcut_in");
        for (i=0; i<ndominos; i++) {
            CCtsp_init_lpdomino (&c->dominos[i]);
        }
        for (i = 0; i < ndominos; i++) {
            rval = CCtsp_read_lpdomino (f, &c->dominos[i], ncount);
            CCcheck_rval (rval, "CCtsp_read_lpdomino failed");
        }
    }

    rval = CCutil_sread_int (f, &rhs);
    CCcheck_rval (rval, "CCutil_sread_int failed");
    c->rhs = rhs;

    rval = CCutil_sread_char (f, &sense);
    CCcheck_rval (rval, "CCutil_sread_int failed");
    c->sense = sense;

    c->branch = 0;

    rval = CCtsp_read_skeleton (f, &c->skel, ncount);
    if (rval) {
        fprintf (stderr, "CCtsp_read_skeleton failed\n");
        goto CLEANUP;
    }
                  
    rval = 0;

CLEANUP:

    if (rval) {
        CCtsp_free_lpcut_in (c);
    }
    return rval;
}

int CCtsp_read_lpclique (CC_SFILE *f, CCtsp_lpclique *c, int ncount)
{
    int rval;
    int size;
    int i;
    int lo;
    int hi;
    int nbits = CCutil_sbits (ncount);

    c->nodes = (CCtsp_segment *) NULL;
    
    rval = CCutil_sread_bits (f, &size, nbits);
    if (rval) {
        fprintf (stderr, "CCutil_sread_int failed\n");
        goto CLEANUP;
    }
    c->nodes = CC_SAFE_MALLOC (size, CCtsp_segment);
    if (c->nodes == (CCtsp_segment *) NULL) {
        fprintf (stderr, "Out of memory in CCtsp_read_lpclique\n");
        rval = 1; goto CLEANUP;
    }
    c->segcount = size;

    for (i=0; i < size; i++) {
        rval = CCutil_sread_bits (f, &lo, nbits);
        if (rval) {
            fprintf (stderr, "CCutil_sread_int failed\n");
            goto CLEANUP;
        }
        rval = CCutil_sread_bits (f, &hi, nbits);
        if (rval) {
            fprintf (stderr, "CCutil_sread_int failed\n");
            goto CLEANUP;
        }
        c->nodes[i].lo = lo;
        c->nodes[i].hi = hi;
    }
    rval = 0;
    
CLEANUP:

    if (rval) {
        CCtsp_free_lpclique (c);
    }
    return rval;
}

int CCtsp_read_lpdomino (CC_SFILE *f, CCtsp_lpdomino *d, int ncount)
{
    int rval = 0;
    int k;

    CCtsp_init_lpdomino (d);
    
    for (k = 0; k < 2; k++) {
        rval = CCtsp_read_lpclique (f, &(d->sets[k]), ncount);
        CCcheck_rval (rval, "CCtsp_read_lpclique failed");
    }

CLEANUP:

    if (rval) {
        CCtsp_free_lpdomino (d);
    }
    return rval;
}

int CCtsp_write_cuts (CC_SFILE *f, int ncount, CCtsp_lpcuts *cuts,
        int writemods)
{
    int *marks = (int *) NULL;
    int *dmarks = (int *) NULL;
    int cend = cuts->cliqueend;
    int dend = cuts->dominoend;
    int i, j;
    int rval;
    int cnt = 0, dcnt = 0;
    int cbits, dbits, nbits;

    rval = CCutil_swrite_char (f, PROB_CUTS_VERSION);
    CCcheck_rval (rval, "CCutil_swrite_char failed");
    rval = CCutil_swrite_int (f, ncount);
    CCcheck_rval (rval, "CCutil_swrite_int failed");
    nbits = CCutil_sbits (ncount);
    
    if (cend > 0) {
        marks = CC_SAFE_MALLOC (cend, int);
        CCcheck_NULL (marks, "out of memory in CCtsp_write_cuts");
        for (i = 0; i < cend; i++) marks[i] = 0;

        for (i = 0; i < cuts->cutcount; i++) {
            for (j = 0; j < cuts->cuts[i].cliquecount; j++) {
                marks[cuts->cuts[i].cliques[j]]++;
            }
        }
        for (i = 0; i < cend; i++) {
            if (marks[i]) {
                if (marks[i] != cuts->cliques[i].refcount) {
                    fprintf (stderr, "ERROR in refcount for clique %d\n", i);
                    rval = 1;  goto CLEANUP;
                }
                marks[i] = cnt+1;
                cnt++;
            }
        }
    }

    cbits = CCutil_sbits (cnt);
    rval = CCutil_swrite_int (f, cnt);
    if (rval) goto CLEANUP;
        
    for (i = 0; i < cend; i++) {
        if (marks[i]) {
            rval = CCtsp_write_lpclique (f, &cuts->cliques[i], ncount);
            CCcheck_rval (rval, "CCtsp_write_lpclique failed");
        }
    }


    if (dend > 0) {
        dmarks = CC_SAFE_MALLOC (dend, int);
        CCcheck_NULL (dmarks, "out of memory in CCtsp_write_cuts");
        for (i = 0; i < dend; i++) dmarks[i] = 0;

        for (i = 0; i < cuts->cutcount; i++) {
            for (j = 0; j < cuts->cuts[i].dominocount; j++) {
                dmarks[cuts->cuts[i].dominos[j]]++;
            }
        }
        for (i = 0; i < dend; i++) {
            if (dmarks[i]) {
                if (dmarks[i] != cuts->dominos[i].refcount) {
                    fprintf (stderr, "ERROR in ref for domino %d (%d, %d)\n",
                                      i, dmarks[i], cuts->dominos[i].refcount);
                    rval = 1;  goto CLEANUP;
                }
                dmarks[i] = dcnt+1;
                dcnt++;
            }
        }
    }

    dbits = CCutil_sbits (dcnt);
    rval = CCutil_swrite_int (f, dcnt);
    if (rval) goto CLEANUP;
        
    for (i = 0; i < dend; i++) {
        if (dmarks[i]) {
            rval = CCtsp_write_lpdomino (f, &cuts->dominos[i], ncount);
            CCcheck_rval (rval, "CCtsp_write_lpdomino failed");
        }
    }


    rval = CCutil_swrite_int (f, cuts->cutcount);
    if (rval) goto CLEANUP;

    for (i = 0; i < cuts->cutcount; i++) {
        rval = CCutil_swrite_int (f, cuts->cuts[i].cliquecount);
        if (rval) goto CLEANUP;
        rval = CCutil_swrite_int (f, cuts->cuts[i].dominocount);
        if (rval) goto CLEANUP;
        rval = CCutil_swrite_int (f, cuts->cuts[i].rhs);
        if (rval) goto CLEANUP;
        rval = CCutil_swrite_char (f, cuts->cuts[i].sense);
        if (rval) goto CLEANUP;

        for (j = 0; j < cuts->cuts[i].cliquecount; j++) {
            rval = CCutil_swrite_bits (f,marks[cuts->cuts[i].cliques[j]] - 1,
                                       cbits);
            CCcheck_rval (rval, "CCutil_swrite_bits failed");
        }
        for (j = 0; j < cuts->cuts[i].dominocount; j++) {
            rval = CCutil_swrite_bits (f,dmarks[cuts->cuts[i].dominos[j]] - 1,
                                       dbits);
            CCcheck_rval (rval, "CCutil_swrite_bits failed");
        }
        if (writemods) {
            rval = CCutil_swrite_int (f, cuts->cuts[i].modcount);
            if (rval) goto CLEANUP;
            for (j = 0; j < cuts->cuts[i].modcount; j++) {
                rval = CCutil_swrite_bits (f, cuts->cuts[i].mods[j].node,
                                           nbits);
                if (rval) goto CLEANUP;
                rval = CCutil_swrite_int (f, cuts->cuts[i].mods[j].mult);
                if (rval) goto CLEANUP;
            }
        }
        rval = CCtsp_write_skeleton (f, &cuts->cuts[i].skel, ncount);
        CCcheck_rval (rval, "CCtsp_write_skeleton failed");
    }

    rval = 0;

CLEANUP:

    CC_IFFREE (marks, int);
    CC_IFFREE (dmarks, int);
    return rval;
}

int CCtsp_send_newcuts (int ncount, CCtsp_lpcuts *pool, char *remotehost,
        unsigned short remoteport)
{
    CC_SFILE *f = (CC_SFILE *) NULL;
    int i;
    int rval = 0;

#ifdef CC_NETREADY
    f = CCutil_snet_open (remotehost, remoteport);
    if (f == (CC_SFILE *) NULL) {
        fprintf (stderr, "CCutil_snet_open failed\n");
        rval = 1; goto CLEANUP;
    }
#endif /* CC_NETREADY */

    rval = CCutil_swrite_char (f, CCtsp_POOL_PUTCUTS);
    CCcheck_rval (rval, "CCutil_swrite_char failed");

    rval = CCutil_swrite_int (f, ncount);
    CCcheck_rval (rval, "CCutil_swrite_int failed");

    rval = CCutil_swrite_int (f, pool->cutcount - pool->savecount);
    CCcheck_rval (rval, "CCutil_swrite_int failed");

    for (i=pool->savecount; i<pool->cutcount; i++) {
        rval = CCtsp_write_lpcut (f, pool, &pool->cuts[i], ncount);
        CCcheck_rval (rval, "CCtsp_write_lpcut failed");
    }

    rval = CCutil_sclose (f);
    f = (CC_SFILE *) NULL;
    CCcheck_rval (rval, "CCutil_sclose failed");

    pool->savecount = pool->cutcount;

CLEANUP:

    if (f != (CC_SFILE *) NULL) {
        CCutil_sclose (f);
    }
    return rval;
}

int CCtsp_write_lpcut_in (CC_SFILE *f, CCtsp_lpcut_in *c, int ncount)
{
    int rval = 0;
    int i;

    rval = CCutil_swrite_int (f, c->cliquecount);
    CCcheck_rval (rval, "CCutil_swrite_int failed");

    rval = CCutil_swrite_int (f, c->dominocount);
    CCcheck_rval (rval, "CCutil_swrite_int failed");

    for (i = 0; i < c->cliquecount; i++) {
        rval = CCtsp_write_lpclique (f, &c->cliques[i], ncount);
        CCcheck_rval (rval, "CCtsp_write_lpclique failed");
    }

    for (i = 0; i < c->dominocount; i++) {
        rval = CCtsp_write_lpdomino (f, &c->dominos[i], ncount);
        CCcheck_rval (rval, "CCtsp_write_lpdomino failed");
    }

    rval = CCutil_swrite_int (f, c->rhs);
    CCcheck_rval (rval, "CCutil_swrite_int failed");

    rval = CCutil_swrite_char (f, c->sense);
    CCcheck_rval (rval, "CCutil_swrite_char failed");

    rval = CCtsp_write_skeleton (f, &c->skel, ncount);
    CCcheck_rval (rval, "CCtsp_write_skeleton failed");
    
CLEANUP:

    return rval;
}

int CCtsp_write_lpcut (CC_SFILE *f, CCtsp_lpcuts *cuts, CCtsp_lpcut *c,
        int ncount)
{
    int rval = 0;
    int i;

    rval = CCutil_swrite_int (f, c->cliquecount);
    CCcheck_rval (rval, "CCutil_swrite_int failed");

    rval = CCutil_swrite_int (f, c->dominocount);
    CCcheck_rval (rval, "CCutil_swrite_int failed");

    for (i = 0; i < c->cliquecount; i++) {
        rval = CCtsp_write_lpclique (f, &cuts->cliques[c->cliques[i]], ncount);
        CCcheck_rval (rval, "CCtsp_write_lpclique failed");
    }

    for (i = 0; i < c->dominocount; i++) {
        rval = CCtsp_write_lpdomino (f, &cuts->dominos[c->dominos[i]], ncount);
        CCcheck_rval (rval, "CCtsp_write_lpdomino failed");
    }

    rval = CCutil_swrite_int (f, c->rhs);
    CCcheck_rval (rval, "CCutil_swrite_int failed");

    rval = CCutil_swrite_char (f, c->sense);
    CCcheck_rval (rval, "CCutil_swrite_char failed");

    rval = CCtsp_write_skeleton (f, &c->skel, ncount);
    CCcheck_rval (rval, "CCtsp_write_skeleton failed");
    
CLEANUP:

    return rval;
}

int CCtsp_write_lpclique (CC_SFILE *f, CCtsp_lpclique *c, int ncount)
{
    int rval = 0;
    int i;
    int nbits = CCutil_sbits (ncount);

    rval = CCutil_swrite_bits (f, c->segcount, nbits);
    CCcheck_rval (rval, "CCutil_swrite_bits failed");

    for (i=0; i < c->segcount; i++) {
        rval = CCutil_swrite_bits (f, c->nodes[i].lo, nbits);
        CCcheck_rval (rval, "CCutil_swrite_bits failed");

        rval = CCutil_swrite_bits (f, c->nodes[i].hi, nbits);
        CCcheck_rval (rval, "CCutil_swrite_bits failed");
    }

CLEANUP:

    return rval;
}

int CCtsp_write_lpdomino (CC_SFILE *f, CCtsp_lpdomino *c, int ncount)
{
    int rval = 0;
    int k;

    for (k = 0; k < 2; k++) {
        rval = CCtsp_write_lpclique (f, &(c->sets[k]), ncount);
        CCcheck_rval (rval, "CCutil_write_lpcplique failed");
    }

CLEANUP:

    return rval;
}

int CCtsp_write_cutpool (int ncount, const char *poolfilename,
        CCtsp_lpcuts *pool)
{
    CC_SFILE *out;
    int rval = 0;

    if (!poolfilename) {
        fprintf (stderr, "pool file name not set\n");
        return 1;
    }

    out = CCutil_sopen (poolfilename, "w");
    if (!out) {
        fprintf (stderr, "CCutil_sopen failed\n");
        return 1;
    }

    rval = CCtsp_write_cuts (out, ncount, pool, 0);
    if (rval) {
        fprintf (stderr, "CCtsp_write_cuts failed\n");
        CCutil_sclose (out);
        return 1;
    }

    CCutil_sclose (out);
    return 0;
}

int CCtsp_copy_cuts (CC_SFILE *f, CC_SFILE *t, int copymods)
{
    int rval;
    int cliqcount, dominocount;
    int cutcount;
    int cnt, ndom;
    int cbits, dbits;
    int nbits;
    char version;
    char ch;
    int ncount;
    CCtsp_lpclique c;
    CCtsp_lpdomino d;
    CCtsp_skeleton skel;
    int i;
    int j;
    int k;

    CCtsp_init_lpclique (&c);
    CCtsp_init_lpdomino (&d);
    CCtsp_init_skeleton (&skel);

    rval = CCutil_sread_char (f, &version);
    if (rval) goto CLEANUP;
    rval = CCutil_swrite_char (t, PROB_CUTS_VERSION);
    if (rval) goto CLEANUP;

    if (version != 1 && version != 2) {
        fprintf (stderr, "Unknown cuts version %d\n", (unsigned) version);
        rval = 1;  goto CLEANUP;
    }

    rval = CCutil_sread_int (f, &ncount);
    if (rval) goto CLEANUP;
    rval = CCutil_swrite_int (t, ncount);
    if (rval) goto CLEANUP;
    nbits = CCutil_sbits (ncount);
    
    rval = CCutil_sread_int (f, &cliqcount);
    if (rval) goto CLEANUP;
    rval = CCutil_swrite_int (t, cliqcount);
    if (rval) goto CLEANUP;
        
    for (i = 0; i < cliqcount; i++) {
        rval = CCtsp_read_lpclique (f, &c, ncount);
        if (rval) goto CLEANUP;
        rval = CCtsp_write_lpclique (t, &c, ncount);
        if (rval) goto CLEANUP;
        CCtsp_free_lpclique (&c);
    }

    if (version == 1) {
        dominocount = 0;
    } else {
        rval = CCutil_sread_int (f, &dominocount);
        CCcheck_rval (rval, "CCutil_sread_int failed");
    }
    rval = CCutil_swrite_int (t, dominocount);
    CCcheck_rval (rval, "CCutil_swrite_int failed");

    for (i = 0; i < dominocount; i++) {
        rval = CCtsp_read_lpdomino (f, &d, ncount);
        CCcheck_rval (rval, "CCtsp_read_lpdomino failed");
        rval = CCtsp_write_lpdomino (t, &d, ncount);
        CCcheck_rval (rval, "CCtsp_write_lpdomino failed");
    }

    cbits = CCutil_sbits (cliqcount);
    dbits = CCutil_sbits (dominocount);

    rval = CCutil_sread_int (f, &cutcount);
    if (rval) goto CLEANUP;
    rval = CCutil_swrite_int (t, cutcount);
    if (rval) goto CLEANUP;

    for (i = 0; i < cutcount; i++) {
        rval = CCutil_sread_int (f, &cnt);
        if (rval) goto CLEANUP;
        rval = CCutil_swrite_int (t, cnt);
        if (rval) goto CLEANUP;

        if (version == 1) {
            ndom = 0;
        } else {
            rval = CCutil_sread_int (f, &ndom);
            CCcheck_rval (rval, "CCutil_sread_int failed");
        }
        rval = CCutil_swrite_int (t, ndom);
        if (rval) goto CLEANUP;

        rval = CCutil_sread_int (f, &j); /* rhs */
        if (rval) goto CLEANUP;
        rval = CCutil_swrite_int (t, j);
        if (rval) goto CLEANUP;
        rval = CCutil_sread_char (f, &ch); /* sense */
        if (rval) goto CLEANUP;
        rval = CCutil_swrite_char (t, ch);
        if (rval) goto CLEANUP;

        for (j = 0; j < cnt; j++) {
            rval = CCutil_sread_bits (f, &k, cbits);
            if (rval) goto CLEANUP;
            rval = CCutil_swrite_bits (t, k, cbits);
            if (rval) goto CLEANUP;
        }

        for (j = 0; j < ndom; j++) {
            rval = CCutil_sread_bits (f, &k, dbits);
            if (rval) goto CLEANUP;
            rval = CCutil_swrite_bits (t, k, dbits);
            if (rval) goto CLEANUP;
        }

        if (copymods) {
            rval = CCutil_sread_int (f, &cnt);
            if (rval) goto CLEANUP;
            rval = CCutil_swrite_int (t, cnt);
            if (rval) goto CLEANUP;
            for (j = 0; j < cnt; j++) {
                rval = CCutil_sread_bits (f, &k, nbits);
                if (rval) goto CLEANUP;
                rval = CCutil_swrite_bits (t, k, nbits);
                if (rval) goto CLEANUP;
                rval = CCutil_sread_int (f, &k);/* mult */
                if (rval) goto CLEANUP;
                rval = CCutil_swrite_int (t, k);
                if (rval) goto CLEANUP;
            }
        }

        rval = CCtsp_read_skeleton (f, &skel, ncount);
        if (rval) goto CLEANUP;
        rval = CCtsp_write_skeleton (t, &skel, ncount);
        if (rval) goto CLEANUP;
        CCtsp_free_skeleton (&skel);
    }

CLEANUP:

    CCtsp_free_lpclique (&c);
    CCtsp_free_lpdomino (&d);
    CCtsp_free_skeleton (&skel);
    return rval;
}

int CCtsp_search_cutpool (CCtsp_lpcuts *pool, CCtsp_lpcut_in **cuts,
        int *cutcount, double *maxviol, int ncount, int ecount, int *elist,
        double *x, int nthreads, CCrandstate *rstate)
{
    int rval = 0;
    double *cval = (double *) NULL;
    int *ind = (int *) NULL;
    int i;
    CCtsp_lpcut_in *newc;
    double lmaxviol;

/*
    printf ("CCtsp_search_cutpool (%d)\n", pool->cutcount);
    fflush (stdout);
*/

    for (i = 0; i < pool->cutcount; i++) {
        if (pool->cuts[i].dominocount != 0) {
            fprintf (stderr, "POOL yipes %d\n", pool->cuts[i].dominocount);
            rval = 1; goto CLEANUP;
        }
    }

    *cutcount = 0;
    *maxviol = 0.0;
    *cuts = (CCtsp_lpcut_in *) NULL;

    if (pool->cutcount == 0) return 0;

    cval = CC_SAFE_MALLOC (pool->cutcount, double);
    if (!cval) {
        fprintf (stderr, "out of memory in CCtsp_search_cutpool\n");
        rval = 1; goto CLEANUP;
    }

    if (nthreads > 0) {
        rval = CCtsp_price_cuts_threaded (pool, ncount, ecount, elist, x, cval,
                                          nthreads);
        if (rval) {
            fprintf (stderr, "CCtsp_price_cuts_threaded failed\n");
            goto CLEANUP;
        }
    } else {
        rval = CCtsp_price_cuts (pool, ncount, ecount, elist, x, cval);
        if (rval) {
            fprintf (stderr, "CCtsp_price_cuts failed\n");
            goto CLEANUP;
        }
    }

    ind = CC_SAFE_MALLOC (pool->cutcount, int);
    if (!ind) {
        fprintf (stderr, "out of memory in CCtsp_search_cutpool\n");
        rval = 1; goto CLEANUP;
    }

    for (i = 0; i < pool->cutcount; i++) {
        ind[i] = i;
    }

    CCutil_rselect (ind, 0, pool->cutcount - 1, POOL_MAXCUTS, cval, rstate);

    lmaxviol = 0.0;
    for (i = 0; i < pool->cutcount && i < POOL_MAXCUTS; i++) {
        if (cval[ind[i]] < lmaxviol) lmaxviol = cval[ind[i]];
        if (cval[ind[i]] < -POOL_MINVIOL) {
            newc = CC_SAFE_MALLOC (1, CCtsp_lpcut_in);
            CCcheck_NULL (newc, "out of memory in CCtsp_search_cutpool");
            rval = CCtsp_lpcut_to_lpcut_in (pool, &pool->cuts[ind[i]], newc);
            if (rval) {
                fprintf (stderr, "CCtsp_lpcut_to_lpcut_in failed\n");
                CC_FREE (newc, CCtsp_lpcut_in);
                goto CLEANUP;
            }
            newc->next = *cuts;
            *cuts = newc;
            (*cutcount)++;
        }
    }
    *maxviol = -lmaxviol;
    rval = 0;

CLEANUP:

    CC_IFFREE (cval, double);
    CC_IFFREE (ind, int);
    return rval;
}

int CCtsp_search_remotepool (char *remotehost, unsigned short remoteport,
        CCtsp_lpcut_in **cuts, int *cutcount, double *maxviol, int ncount,
        int ecount, int *elist, double *x)
{
    CC_SFILE *f = (CC_SFILE *) NULL;
    int rval = 0;
    int i;
    CCtsp_lpcut_in *newc = (CCtsp_lpcut_in *) NULL;
    int cnt;

    *cutcount = 0;
    *cuts = (CCtsp_lpcut_in *) NULL;

#ifdef CC_NETREADY
    f = CCutil_snet_open (remotehost, remoteport);
    if (f == (CC_SFILE *) NULL) {
        fprintf (stderr, "CCutil_snet_open failed\n");
        rval = 1; goto CLEANUP;
    }
#endif /* CC_NETREADY */

    rval = CCutil_swrite_char (f, CCtsp_POOL_GETCUTS);
    if (rval) {
        fprintf (stderr, "CCutil_swrite_char failed\n");
        goto CLEANUP;
    }

    rval = CCutil_swrite_int (f, ncount);
    if (rval) {
        fprintf (stderr, "CCutil_swrite_int failed\n");
        goto CLEANUP;
    }

    rval = CCutil_swrite_int (f, ecount);
    if (rval) {
        fprintf (stderr, "CCutil_swrite_int failed\n");
        goto CLEANUP;
    }

    for (i=0; i<ecount; i++) {
        rval = CCutil_swrite_int (f, elist[2*i]);
        if (rval) {
            fprintf (stderr, "CCutil_swrite_int failed\n");
            goto CLEANUP;
        }
        rval = CCutil_swrite_int (f, elist[2*i+1]);
        if (rval) {
            fprintf (stderr, "CCutil_swrite_int failed\n");
            goto CLEANUP;
        }
        rval = CCutil_swrite_double (f, x[i]);
        if (rval) {
            fprintf (stderr, "CCutil_swrite_double failed\n");
            goto CLEANUP;
        }
    }

    rval = CCutil_sread_int (f, &cnt);
    if (rval) {
        fprintf (stderr, "CCutil_sread_int failed\n");
        goto CLEANUP;
    }

    rval = CCutil_sread_double (f, maxviol);
    if (rval) {
        fprintf (stderr, "CCutil_sread_double failed\n");
        goto CLEANUP;
    }

    for (i=0; i<cnt; i++) {
        newc = CC_SAFE_MALLOC (1, CCtsp_lpcut_in);
        CCcheck_NULL (newc, "out of memory in CCtsp_search_cutpool");
        rval = CCtsp_read_lpcut_in (f, newc, ncount);
        if (rval) {
            fprintf (stderr, "read_lpcut_in failed\n");
            goto CLEANUP;
        }
        newc->next = *cuts;
        *cuts = newc;
        (*cutcount)++;
        newc = (CCtsp_lpcut_in *) NULL;
    }

    rval = CCutil_sclose (f);
    f = (CC_SFILE *) NULL;
    if (rval) {
        fprintf (stderr, "CCutil_sclose failed\n");
        goto CLEANUP;
    }
    
    rval = 0;

CLEANUP:

    if (f != (CC_SFILE *) NULL) {
        CCutil_sclose (f);
    }
    CC_IFFREE (newc, CCtsp_lpcut_in);
    if (rval) {
        fprintf (stderr, "Failure in CCtsp_search_remotepool, continuing anyway\n");
        rval = 0;
    }
    
    return rval;
}

int CCtsp_search_cutpool_cliques (CCtsp_lpcuts *pool, CCtsp_lpclique **cliques,
        int *cliquecount, int ncount, int ecount, int *elist, double *x,
        double maxdelta, int maxcliques, double **cliquevals,
        CCrandstate *rstate)
{
    int rval = 0;
    double *cval = (double *) NULL;
    int *ind = (int *) NULL;
    double upperdelta, lowerdelta;
    int i, k;
    int ccount = 0;

    *cliquecount = 0;
    *cliques = (CCtsp_lpclique *) NULL;
    if (cliquevals) {
        *cliquevals = (double *) NULL;
    }

    if (pool->cutcount == 0) return 0;

    cval = CC_SAFE_MALLOC (pool->cliqueend, double);
    if (!cval) {
        fprintf (stderr, "out of memory in CCtsp_search_cutpool_cliques\n");
        rval = 1; goto CLEANUP;
    }

    rval = price_cliques (pool->cliques, ncount, ecount, elist, x, cval,
                          pool->cliqueend);
    if (rval) {
        fprintf (stderr, "price_cliques failed\n");
        goto CLEANUP;
    }

    ind = CC_SAFE_MALLOC (pool->cliqueend, int);
    if (!ind) {
        fprintf (stderr, "out of memory in CCtsp_search_cutpool_cliques\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < pool->cliqueend; i++) {
        ind[i] = i;
    }

    CCutil_rselect (ind, 0, pool->cliqueend - 1, maxcliques, cval, rstate);

    upperdelta = -1.0;
    lowerdelta = maxdelta;
    for (i = 0; i < pool->cliqueend && i < maxcliques; i++) {
        if (cval[ind[i]] < maxdelta) {
            ccount++;
            if (cval[ind[i]] < lowerdelta) lowerdelta = cval[ind[i]];
            if (cval[ind[i]] > upperdelta) upperdelta = cval[ind[i]];
        }
    }

    if (ccount == 0) {
        printf ("Found no nearly tight cliques\n"); fflush (stdout);
        goto CLEANUP;
    }

    *cliques = CC_SAFE_MALLOC (ccount, CCtsp_lpclique);
    if (!(*cliques)) {
        fprintf (stderr, "out of memory in CCtsp_search_cutpool_cliques\n");
        rval = 1; goto CLEANUP;
    }
    if (cliquevals) {
        *cliquevals = CC_SAFE_MALLOC (ccount, double);
        if (!(*cliquevals)) {
            fprintf (stderr, "out of memory in CCtsp_search_cutpool_cliques\n");
            CC_FREE (*cliques, CCtsp_lpclique);
            rval = 1; goto CLEANUP;
        }
    }

    ccount = 0;
    for (i = 0; i < pool->cliqueend && i < maxcliques; i++) {
        if (cval[ind[i]] < maxdelta) {
            rval = CCtsp_copy_lpclique (&(pool->cliques[ind[i]]),
                                        &((*cliques)[ccount]));
            if (rval) {
                fprintf (stderr, "CCtsp_copy_lpclique failed\n");
                for (k = 0; k < ccount; k++) {
                    CC_FREE ((*cliques)[k].nodes, CCtsp_segment);
                }
                CC_FREE (*cliques, CCtsp_lpclique);
                if (cliquevals) {
                    CC_FREE (*cliquevals, double);
                }
                goto CLEANUP;
            }
            if (cliquevals) {
                (*cliquevals)[ccount] = cval[ind[i]];
            }
            ccount++;
        }
    }
    *cliquecount = ccount;

    printf ("%d nearly tight cliques found, range (%.3f, %.3f)\n",
              *cliquecount, lowerdelta, upperdelta);
    fflush (stdout);

CLEANUP:

    CC_IFFREE (cval, double);
    CC_IFFREE (ind, int);
    return rval;
}

#define BRANCH_CLIQUE_GOAL 3.00
#define BRANCH_CLIQUE_TOL  0.99

int CCtsp_branch_cutpool_cliques (CCtsp_lpcuts *pool, CCtsp_lpclique **cliques,
        int *cliquecount, int ncount, int ecount, int *elist, double *x,
        int nwant, double **cliquevals, int silent)
{
    int rval = 0;
    double *cval = (double *) NULL;
    double upperdelta, lowerdelta;
    double t;
    int i, k;
    int ccount = 0;
    int *blist =   (int *) NULL;
    double *bval = (double *) NULL;

    *cliquecount = 0;
    *cliques = (CCtsp_lpclique *) NULL;
    if (cliquevals) {
        *cliquevals = (double *) NULL;
    }

    if (pool->cutcount == 0 || nwant <= 0) return 0;

    blist = CC_SAFE_MALLOC (nwant + 1, int);
    bval  = CC_SAFE_MALLOC (nwant + 1, double);
    cval = CC_SAFE_MALLOC (pool->cliqueend, double);
    if (!blist || !bval || !cval) {
        fprintf (stderr, "out of memory in CCtsp_search_cutpool_cliques\n");
        rval = 1; goto CLEANUP;
    }

    rval = price_cliques (pool->cliques, ncount, ecount, elist, x, cval,
                          pool->cliqueend);
    if (rval) {
        fprintf (stderr, "price_cliques failed\n");
        goto CLEANUP;
    }

    for (i = 0; i < nwant; i++) {
        blist[i] = -1;
        bval[i]  = CCtsp_LP_MAXDOUBLE;
    }
    blist[nwant] = -1;
    bval[i]      = -1.0;

    for (i = 0; i < pool->cliqueend; i++) {
        t = CC_OURABS (BRANCH_CLIQUE_GOAL - cval[i]);
        if (t < bval[0] && t < BRANCH_CLIQUE_TOL) {
            for (k = 0; t < bval[k+1]; k++) {
                blist[k] = blist[k+1];
                bval[k]  =  bval[k+1];
            }
            blist[k] = i;
            bval[k]  = t;
        }
    }

    upperdelta = -1.0;
    lowerdelta = CCtsp_LP_MAXDOUBLE;
    for (i = 0; i < nwant; i++) {
        if (blist[i] != -1) {
            if (upperdelta < cval[blist[i]]) {
                upperdelta = cval[blist[i]];
            }
            if (lowerdelta > cval[blist[i]]) {
                lowerdelta = cval[blist[i]];
            }
            ccount++;
        }
    }

    if (ccount == 0) {
        printf ("Found no nearly tight cliques\n"); fflush (stdout);
        goto CLEANUP;
    }

    *cliques = CC_SAFE_MALLOC (ccount, CCtsp_lpclique);
    if (!(*cliques)) {
        fprintf (stderr, "out of memory in CCtsp_search_cutpool_cliques\n");
        rval = 1; goto CLEANUP;
    }
    if (cliquevals) {
        *cliquevals = CC_SAFE_MALLOC (ccount, double);
        if (!(*cliquevals)) {
            fprintf (stderr, "out of memory in CCtsp_search_cutpool_cliques\n");
            CC_FREE (*cliques, CCtsp_lpclique);
            rval = 1; goto CLEANUP;
        }
    }

    ccount = 0;
    for (i = nwant - 1; i >= 0; i--) {
        if (blist[i] != -1) {
            rval = CCtsp_copy_lpclique (&(pool->cliques[blist[i]]),
                                      &((*cliques)[ccount]));
            if (rval) {
                fprintf (stderr, "CCtsp_copy_lpclique failed\n");
                for (k = 0; k < ccount; k++) {
                    CC_FREE ((*cliques)[k].nodes, CCtsp_segment);
                }
                CC_FREE (*cliques, CCtsp_lpclique);
                if (cliquevals) {
                    CC_FREE (*cliquevals, double);
                }
                goto CLEANUP;
            }
            if (cliquevals) {
                (*cliquevals)[ccount] = cval[blist[i]];
            }
            ccount++;
        }
    }
    *cliquecount = ccount;

    if (!silent) {
        printf ("%d candidate branching cliques, range (%.3f, %.3f)\n",
                  *cliquecount, lowerdelta, upperdelta);
        fflush (stdout);
    }


CLEANUP:

    CC_IFFREE (blist, int);
    CC_IFFREE (bval, double);
    CC_IFFREE (cval, double);
    return rval;
}

int CCtsp_get_clique_prices (CCtsp_lpcuts *pool, int **p_cliquenums,
        double **p_cliquevals, double mindelta, double maxdelta,
        int *p_cliquecount, int ncount, int ecount, int *elist, double *x)
{
    double *cval = (double *) NULL;
    int i;
    int rval;
    int cliquecount = 0;
    int *cliquenums = (int *) NULL;
    double *cliquevals = (double *) NULL;

    if (p_cliquevals) *p_cliquevals = (double *) NULL;
    *p_cliquenums = (int *) NULL;
    *p_cliquecount = 0;

    if (pool->cutcount == 0) return 0;

    cval = CC_SAFE_MALLOC (pool->cliqueend, double);
    if (cval == (double *) NULL) {
        fprintf (stderr, "Out of memory in CCtsp_get_clique_prices\n");
        rval = 1; goto CLEANUP;
    }

    rval = price_cliques (pool->cliques, ncount, ecount, elist, x, cval,
                          pool->cliqueend);
    if (rval) {
        fprintf (stderr, "price_cliques failed\n");
        goto CLEANUP;
    }

    cliquecount = 0;
    for (i=0; i<pool->cliqueend; i++) {
        if (pool->cliques[i].segcount > 0 &&
            cval[i] >= mindelta && cval[i] <= maxdelta) {
            cliquecount++;
        }
    }

    if (cliquecount == 0) {
        rval = 0; goto CLEANUP;
    }
    
    cliquenums = CC_SAFE_MALLOC (cliquecount, int);
    if (cliquenums == (int *) NULL) {
        fprintf (stderr, "Out of memory in CCtsp_get_clique_prices\n");
        rval = 1; goto CLEANUP;
    }

    if (p_cliquevals) {
        cliquevals = CC_SAFE_MALLOC (cliquecount, double);
        if (cliquevals == (double *) NULL) {
            fprintf (stderr, "Out of memory in CCtsp_get_clique_prices\n");
            rval = 1; goto CLEANUP;
        }
    }

    cliquecount = 0;
    for (i=0; i<pool->cliqueend; i++) {
        if (pool->cliques[i].segcount > 0 &&
            cval[i] >= mindelta && cval[i] <= maxdelta) {
            cliquenums[cliquecount] = i;
            if (cliquevals) cliquevals[cliquecount] = cval[i];
            cliquecount++;
        }
    }

    rval = 0;

  CLEANUP:

    if (rval) {
        CC_IFFREE (cliquenums, int);
        CC_IFFREE (cliquevals, double);
        cliquecount = 0;
    }
    CC_IFFREE (cval, double);
    *p_cliquenums = cliquenums;
    if (p_cliquevals) *p_cliquevals = cliquevals;
    *p_cliquecount = cliquecount;
    return rval;
}

int CCtsp_get_clique (CCtsp_lpcuts *pool, int cliquenum,
        CCtsp_lpclique **p_clique)
{
    if (cliquenum < 0 || cliquenum >= pool->cliqueend ||
        pool->cliques[cliquenum].segcount <= 0) {
        fprintf (stderr, "Illegal cliquenum in CCtsp_get_clique\n");
        return -1;
    }
    *p_clique = &pool->cliques[cliquenum];
    return 0;
}

int CCtsp_add_to_cutpool (CCtsp_lpcuts *pool, CCtsp_lpcuts *cuts,
        CCtsp_lpcut *c)
{
    int rval = 0;
    CCtsp_lpcut_in cin;

    if (!c || c->cliquecount <= 1)
        return 0;

    CCtsp_init_lpcut_in (&cin);

    rval = CCtsp_lpcut_to_lpcut_in (cuts, c, &cin);
    if (rval) {
        fprintf (stderr, "CCtsp_lpcut_to_lpcut_in failed\n");
        goto CLEANUP;
    }

    rval = CCtsp_add_to_cutpool_lpcut_in (pool, &cin);
    if (rval) {
        fprintf (stderr, "CCtsp_add_to_cutpool_lpcut_in failed\n");
        goto CLEANUP;
    }

CLEANUP:

    CCtsp_free_lpcut_in (&cin);
    return rval;
}

int CCtsp_add_to_dominopool (CCtsp_lpcuts *pool, CCtsp_lpcuts *cuts,
        CCtsp_lpcut *c)
{
    int rval = 0;
    CCtsp_lpcut_in cin;

/*
    printf ("CCtsp_add_to_dominopool (%d, %d)\n", c->cliquecount,
                                                  c->dominocount);
    fflush (stdout);
*/

    if (!pool) {
        fprintf (stderr, "NO domino pool!\n");
        rval = 1;  goto CLEANUP;
    }

    if (!c) return 0;

    CCtsp_init_lpcut_in (&cin);

    rval = CCtsp_lpcut_to_lpcut_in (cuts, c, &cin);
    if (rval) {
        fprintf (stderr, "CCtsp_lpcut_to_lpcut_in failed\n");
        goto CLEANUP;
    }

    rval = CCtsp_add_to_cutpool_lpcut_in (pool, &cin);
    if (rval) {
        fprintf (stderr, "CCtsp_add_to_cutpool_lpcut_in failed\n");
        goto CLEANUP;
    }

CLEANUP:

    CCtsp_free_lpcut_in (&cin);
    return rval;
}

int CCtsp_add_to_cutpool_lpcut_in (CCtsp_lpcuts *pool, CCtsp_lpcut_in *cut)
{
    int rval = 0;
    CCtsp_lpcut new;
    int cutloc;
    unsigned int hval;

    if (!pool) goto CLEANUP;

    CCtsp_init_lpcut (&new); 

    new.rhs      = cut->rhs;
    new.branch   = cut->branch;
    new.sense    = cut->sense;

    rval = CCtsp_register_cliques (pool, cut, &new);
    CCcheck_rval (rval, "CCtsp_register_cliques failed");

    rval = CCtsp_register_dominos (pool, cut, &new);
    CCcheck_rval (rval, "CCtsp_register_dominos failed");

    rval = CCtsp_copy_skeleton (&cut->skel, &new.skel);
    if (rval) {
        fprintf (stderr, "CCtsp_copy_skeleton failed\n");
        CCtsp_unregister_cliques (pool, &new);
        CCtsp_unregister_dominos (pool, &new);
        goto CLEANUP;
    }

    sort_cliques (&new);
    sort_dominos (&new);

    cutloc = CCtsp_add_cut_to_cutlist (pool, &new);
    if (cutloc < 0) {
        fprintf (stderr, "CCtsp_add_cut_to_cutlist failed\n");
        CCtsp_unregister_cliques (pool, &new);
        CCtsp_unregister_dominos (pool, &new);
        rval = cutloc;
        goto CLEANUP;
    }

    hval = CCutil_genhash_hash (pool->cuthash, (void *) ((long) cutloc));
    if (CCutil_genhash_lookup_h (pool->cuthash, hval,
                                 (void *) ((long) cutloc))) {
        /* cut was already in pool */
        CCtsp_delete_cut_from_cutlist (pool, cutloc);
        goto CLEANUP;
    }

    rval = CCutil_genhash_insert_h (pool->cuthash, hval,
            (void *) ((long) cutloc),  (void *) ((long) 1));
    if (rval) {
        fprintf (stderr, "CCutil_genhash_insert_h failed\n");
        CCtsp_delete_cut_from_cutlist (pool, cutloc);
        goto CLEANUP; 
    }

CLEANUP:

    return rval;
}

static void sort_cliques (CCtsp_lpcut *c)
{
    CCutil_int_array_quicksort (c->cliques, c->cliquecount);
}

static void sort_dominos (CCtsp_lpcut *c)
{
    CCutil_int_array_quicksort (c->dominos, c->dominocount);
}


static int register_lpcuts (CCtsp_lpcuts *pool)
{
    int i;
    unsigned int hval;
    int rval = 0;
    int ndup = 0;

    for (i=0; i<pool->cutcount; i++) {
        sort_cliques (&pool->cuts[i]);
        sort_dominos (&pool->cuts[i]);
        hval = CCutil_genhash_hash (pool->cuthash, (void *) ((long) i));
        if (CCutil_genhash_lookup_h (pool->cuthash, hval,
                                     (void *) ((long) i))) {
            ndup++;
        } else {
            rval = CCutil_genhash_insert_h (pool->cuthash, hval,
                    (void *) ((long) i), (void *) ((long) 1));
            if (rval) {
                fprintf (stderr, "CCutil_genhash_insert_h failed\n");
                return rval;
            }
        }
    }
    if (ndup) {
        printf ("%d duplicates detected in pool\n", ndup);
        fflush (stdout);
    }
    return 0;
}

int CCtsp_display_cutpool (CCtsp_lpcuts *pool)
{
    int i;
    CCtsp_lpcut_in c;

    for (i = 0; i < pool->cutcount; i++) {
        if (CCtsp_lpcut_to_lpcut_in (pool, &(pool->cuts[i]), &c)) {
            fprintf (stderr, "CCtsp_lpcut_to_lpcut_in failed\n");
            return 1;
        }
        CCtsp_print_lpcut_in (&c);
        CCtsp_free_lpcut_in (&c);
    }

    return 0;
}

int CCtsp_price_cuts (CCtsp_lpcuts *pool, int ncount, int ecount,
        int *elist, double *x, double *cutval)
{
    double *cval = (double *) NULL;
    int rval = 0;

    cval = CC_SAFE_MALLOC (pool->cliqueend, double);
    CCcheck_NULL (cval, "out of memory in CCtsp_price_cuts");

    rval = price_cliques (pool->cliques, ncount, ecount, elist, x, cval,
                          pool->cliqueend);
    CCcheck_rval (rval, "price_cliques failed");

    price_cuts (pool->cuts, pool->cutcount, cval, cutval);

CLEANUP:

    CC_IFFREE (cval, double);
    return rval;
}

#ifdef CC_POSIXTHREADS

typedef struct priceclique_args {
    CCtsp_lpclique *cliques;
    int ncount;
    int ecount;
    int *elist;
    double *x;
    double *cval;
    int cend;
    int rval;
    double real_zeit;
} priceclique_args;

typedef struct pricecut_args {
    CCtsp_lpcut *cuts;
    int cutcount;
    double *cval;
    double *cutval;
    int rval;
    double real_zeit;
} pricecut_args;

static void *price_cliques_thread (void *args)
{
    priceclique_args *clargs = (priceclique_args *) args;
    double szeit = CCutil_real_zeit();

    clargs->rval = price_cliques (clargs->cliques, clargs->ncount, 
        clargs->ecount, clargs->elist, clargs->x, clargs->cval, clargs->cend);

    clargs->real_zeit = CCutil_real_zeit() - szeit;
    
    return clargs;
}

static void *price_cuts_thread (void *args)
{
    pricecut_args *cuargs = (pricecut_args *) args;
    double szeit = CCutil_real_zeit();

    price_cuts (cuargs->cuts, cuargs->cutcount, cuargs->cval, cuargs->cutval); 

    cuargs->real_zeit = CCutil_real_zeit() - szeit;
    cuargs->rval = 0;
    
    return cuargs;
}

/* rebalance_load makes the (unrealistic) assumption that the worktime
 * in each thread was uniformly distributed over the workload for that
 * thread.
 */
static void rebalance_load (int nthreads, double *workload, double *worktime,
        double *work)
{
    double goal;
    int i;
    int j;
    double iload;
    double itime;
    double jload;
    double jtime;

    goal = 0.0;
    for (i=0; i<nthreads; i++) {
        /* try to cope with roundoff in real times */
        if (worktime[i] == 0.0) worktime[i] = 1.0;
        goal += worktime[i];
    }
    goal /= nthreads;

    i = 0;
    j = 0;
    jload = workload[j];
    jtime = worktime[j];
    iload = 0.0;
    itime = goal;
    while (i < nthreads) {
        if (jtime == 0.0 || itime >= jtime) {
            itime -= jtime;
            iload += jload;
            j++;
            if (j < nthreads) {
                jload = workload[j];
                jtime = worktime[j];
            } else {
                jload = 0.0;
                jtime = 1e10;
            }
        } else {
            iload += itime / jtime * jload;
            jload -= itime / jtime * jload;
            jtime -= itime;
            work[i] = iload;
            i++;
            iload = 0.0;
            itime = goal;
        }
    }
    
    goal = 0.0;
    for (i=0; i<nthreads; i++) {
        goal += work[i];
    }
    if (goal == 0.0) goal = 1.0;
    for (i=0; i<nthreads; i++) {
        workload[i] = work[i] / goal;
    }
}

#endif /* CC_POSIXTHREADS */

int CCtsp_price_cuts_threaded (CCtsp_lpcuts *pool, int ncount, int ecount,
        int *elist, double *x, double *cutval, CC_UNUSED int nthreads)
{
#ifndef CC_POSIXTHREADS
    return CCtsp_price_cuts (pool, ncount, ecount, elist, x, cutval);
#else /* CC_POSIXTHREADS */
    double *cval = (double *) NULL;
    priceclique_args *clargs = (priceclique_args *) NULL;
    priceclique_args *clrval;
    pricecut_args *cuargs = (pricecut_args *) NULL;
    pricecut_args *curval;
    void *thr_rval;
    pthread_t *thread_id = (pthread_t *) NULL;
    pthread_attr_t attr;
    int cstart;
    int cend;
    int i;
    double *cliqueworkload;
    double *cutworkload;
    double *worktimes;
    double *balancework;
    int rval = 0;

    if (pool->workloads != (double *) NULL &&
        pool->workloads[0] != (double) nthreads) {
        CC_FREE (pool->workloads, double);
    }
    
    if (pool->workloads == (double *) NULL) {
        pool->workloads = CC_SAFE_MALLOC (nthreads*4+1, double);
        if (pool->workloads == (double *) NULL) {
            fprintf (stderr, "Out of memory in CCtsp_price_cuts_threaded\n");
            rval = 1; goto CLEANUP;
        }
        pool->workloads[0] = (double) nthreads;
        cliqueworkload = pool->workloads + 1;
        cutworkload    = pool->workloads + 1 + nthreads;
        for (i=0; i<nthreads; i++) {
            cliqueworkload[i] = 1.0 / ((double) nthreads);
            cutworkload[i]    = 1.0 / ((double) nthreads);
        }
    }

    /* nthreads, two workloads, and worktimes are all packed in workloads */
    cliqueworkload = pool->workloads + 1;
    cutworkload    = pool->workloads + 1 + nthreads;
    worktimes      = pool->workloads + 1 + nthreads + nthreads;
    balancework    = pool->workloads + 1 + nthreads + nthreads + nthreads;
    
    cval      = CC_SAFE_MALLOC (pool->cliqueend, double);
    clargs    = CC_SAFE_MALLOC (nthreads, priceclique_args);
    cuargs    = CC_SAFE_MALLOC (nthreads, pricecut_args);
    thread_id = CC_SAFE_MALLOC (nthreads, pthread_t);
    if (cval      == (double *) NULL ||
        clargs    == (priceclique_args *) NULL ||
        cuargs    == (pricecut_args *) NULL ||
        thread_id == (pthread_t *) NULL) {
        fprintf (stderr, "out of memory in CCtsp_price_cuts\n");
        rval = 1; goto CLEANUP;
    }

    rval = pthread_attr_init (&attr);
    if (rval) {
        fprintf (stderr, "pthread_attr_init failed, rval %d\n", rval);
        rval = 1; goto CLEANUP;
    }

#ifdef PTHREAD_CREATE_JOINABLE
    rval = pthread_attr_setdetachstate (&attr, PTHREAD_CREATE_JOINABLE);
    if (rval) {
        fprintf (stderr, "pthread_attr_setdetachstate failed, rval %d\n",
                 rval);
        rval = 1; goto CLEANUP;
    }
#else
#ifdef PTHREAD_CREATE_UNDETACHED
    rval = pthread_attr_setdetachstate (&attr, PTHREAD_CREATE_UNDETACHED);
    if (rval) {
        fprintf (stderr, "pthread_attr_setdetachstate failed, rval %d\n",
                 rval);
        rval = 1; goto CLEANUP;
    }
#endif
#endif
    
    cstart = 0;
    for (i=0; i<nthreads; i++) {
        if (i == nthreads-1) {
            cend = pool->cliqueend;
        } else {
            cend = cstart + (int) (cliqueworkload[i] * pool->cliqueend);
        }
        clargs[i].cliques = pool->cliques + cstart;
        clargs[i].ncount  = ncount;
        clargs[i].ecount  = ecount;
        clargs[i].elist   = elist;
        clargs[i].x       = x;
        clargs[i].cval    = cval + cstart;
        clargs[i].cend    = cend - cstart;
        clargs[i].rval    = 0;
        cstart = cend;
        rval = pthread_create (&thread_id[i], &attr, price_cliques_thread,
                               &clargs[i]);
        if (rval) {
            perror ("pthread_create");
            fprintf (stderr, "pthread_create failed, rval %d\n", rval);
            goto CLEANUP;
        }
    }
    for (i=0; i<nthreads; i++) {
        rval = pthread_join (thread_id[i], &thr_rval);
        if (rval) {
            fprintf (stderr, "pthread_join failed\n");
            goto CLEANUP;
        }
        clrval = (priceclique_args *) thr_rval;
        if (clrval->rval) {
            fprintf (stderr, "pricing clique thread failed\n");
            rval = clrval->rval; goto CLEANUP;
        }
        worktimes[i] = clrval->real_zeit;
    }

    printf ("\nThread clique:");
    for (i=0; i<nthreads; i++) {
        printf (" %.0f", worktimes[i]);
    }
    fflush (stdout);
    
    rebalance_load (nthreads, cliqueworkload, worktimes, balancework);

    cstart = 0;
    for (i=0; i<nthreads; i++) {
        if (i == nthreads-1) {
            cend = pool->cutcount;
        } else {
            cend = cstart + cutworkload[i] * pool->cutcount;
        }
        cuargs[i].cuts     = pool->cuts + cstart;
        cuargs[i].cutcount = cend - cstart;
        cuargs[i].cval     = cval;
        cuargs[i].cutval   = cutval + cstart;
        cuargs[i].rval     = 0;
        cstart = cend;
        rval = pthread_create (&thread_id[i], &attr, price_cuts_thread,
                               &cuargs[i]);
        if (rval) {
            fprintf (stderr, "pthread_create failed\n");
            goto CLEANUP;
        }
    }
    for (i=0; i<nthreads; i++) {
        rval = pthread_join (thread_id[i], &thr_rval);
        if (rval) {
            fprintf (stderr, "pthread_join failed\n");
            goto CLEANUP;
        }
        curval = (pricecut_args *) thr_rval;
        if (curval->rval) {
            fprintf (stderr, "pricing cut thread failed\n");
            rval = curval->rval; goto CLEANUP;
        }
        worktimes[i] = curval->real_zeit;
    }
    
    printf (" cut:");
    for (i=0; i<nthreads; i++) {
        printf (" %.0f", worktimes[i]);
    }
    printf ("\n");
    fflush (stdout);
    
    rebalance_load (nthreads, cutworkload, worktimes, balancework);

    rval = 0;

 CLEANUP:
    CC_IFFREE (cval, double);
    CC_IFFREE (clargs, priceclique_args);
    CC_IFFREE (cuargs, pricecut_args);
    CC_IFFREE (thread_id, pthread_t);
    return rval;
#endif /* CC_POSIXTHREADS */
}

static void price_cuts (CCtsp_lpcut *cuts, int cutcount, double *cval,
        double *cutval)
{
    int i, j;
    CCtsp_lpcut *c;
    double v;
    
    for (i = 0, c = cuts; i < cutcount; i++, c++) {
        if (c->dominocount == 0) {
            v = (double) -(c->rhs);
            for (j  = 0; j < c->cliquecount; j++)  {
                v += cval[c->cliques[j]];
            }
            cutval[i] = v;
        } else {
            cutval[i] = 1000.0;  /* For now, do not price dominos */
            /* We can easily add domino pricing, but the parallel */
            /* code will be not be so easy, since dominos are     */
            /* shared among the threads.                          */  
        }
    }
}

static int price_cliques (CCtsp_lpclique *cliques, int ncount, int ecount,
        int *elist, double *x, double *cval, int cend)
{
    poolnode *nlist = (poolnode *) NULL;
    pooledge *espace = (pooledge *) NULL;
    int marker = 0;
    int i;
    int rval = 0;

    rval = make_pricing_graph (ncount, ecount, elist, x, &nlist, &espace);
    if (rval) {
        fprintf (stderr, "make_pricing_graph failed\n");
        goto CLEANUP;
    }
    for (i = 0; i < cend; i++) {
        if (cliques[i].segcount > 0) {
            marker++;
            cval[i] = price_clique (nlist, &(cliques[i]), marker);
        } else {
            cval[i] = -1.0;
        }
    }

CLEANUP:

    CC_IFFREE (nlist, poolnode);
    CC_IFFREE (espace, pooledge);
    return rval;
}

static int make_pricing_graph (int ncount, int ecount, int *elist, double *x,
    poolnode **p_nlist, pooledge **p_espace)
{
    poolnode *nlist = (poolnode *) NULL;
    pooledge *espace = (pooledge *) NULL;
    pooledge *p;
    int i;
    int count;
    int a, b;
    int rval;

    *p_nlist = (poolnode *) NULL;
    *p_espace = (pooledge *) NULL;
    
    nlist =  CC_SAFE_MALLOC (ncount, poolnode);
    if (nlist == (poolnode *) NULL) {
        fprintf (stderr, "out of memory in make_pricing_graph\n");
        rval = 1; goto CLEANUP;
    }

    for (i = 0; i < ncount; i++) {
        nlist[i].mark = 0;
        nlist[i].deg = 0;
    }

    count = 0;
    for (i = 0; i < ecount; i++) {
        if (x[i] >= ZERO_EPSILON) {
            nlist[elist[2*i]].deg++;
            nlist[elist[2*i+1]].deg++;
            count++;
        }
    }

    espace = CC_SAFE_MALLOC (2*count, pooledge);
    if (espace == (pooledge *) NULL) {
        fprintf (stderr, "out of memory in price_cliques\n");
        rval = 1; goto CLEANUP;
    }

    p = espace;
    for (i = 0; i < ncount; i++) {
        nlist[i].adj = p;
        p += nlist[i].deg;
        nlist[i].deg = 0;
    }
    for (i = 0; i < ecount; i++) {
        if (x[i] >= ZERO_EPSILON) {
            a = elist[2*i];
            b = elist[2*i+1];
            nlist[a].adj[nlist[a].deg].x = x[i];
            nlist[a].adj[nlist[a].deg++].to = b;
            nlist[b].adj[nlist[b].deg].x = x[i];
            nlist[b].adj[nlist[b].deg++].to = a;
        }
    }

    *p_nlist = nlist;
    *p_espace = espace;
    
    rval = 0;
 CLEANUP:
    if (rval) {
        CC_IFFREE (nlist, poolnode);
        CC_IFFREE (espace, pooledge);
    }
    return rval;
}

static double price_clique (poolnode *nlist, CCtsp_lpclique *c, int marker)
{
    double val = 0.0;
    poolnode *n;
    int tmp, j, k;

    CC_FOREACH_NODE_IN_CLIQUE (j, *c, tmp) {
        nlist[j].mark = marker;
    }
    CC_FOREACH_NODE_IN_CLIQUE (j, *c, tmp) {
        n = &(nlist[j]);
        for (k = 0; k < n->deg; k++) {
            if (nlist[n->adj[k].to].mark != marker) {
                val += n->adj[k].x;
            }
        }
    }
    return val;
}

void CCtsp_free_lpcut_in (CCtsp_lpcut_in *c)
{
    int i;

    if (c != (CCtsp_lpcut_in *) NULL) {
        if (c->cliques != (CCtsp_lpclique *) NULL) {
            for (i = 0; i < c->cliquecount; i++) {
                CCtsp_free_lpclique (&c->cliques[i]);
            }
            CC_IFFREE (c->cliques, CCtsp_lpclique);
        }
        if (c->dominos != (CCtsp_lpdomino *) NULL) {
            for (i = 0; i < c->dominocount; i++) {
                CCtsp_free_lpdomino (&c->dominos[i]);
            }
            CC_IFFREE (c->dominos, CCtsp_lpdomino);
        }
        CCtsp_free_skeleton (&c->skel);
    }
}

void CCtsp_free_lpclique (CCtsp_lpclique *c)
{
    if (c) {
        CC_IFFREE (c->nodes, CCtsp_segment);
        c->segcount = 0;
        c->hashnext = 0;
        c->refcount = 0;
    }
}

void CCtsp_free_lpdomino (CCtsp_lpdomino *c)
{
    if (c) {
        CCtsp_free_lpclique (&(c->sets[0]));
        CCtsp_free_lpclique (&(c->sets[1]));
        c->hashnext = 0;
        c->refcount = 0;
    }
}

int CCtsp_register_cliques (CCtsp_lpcuts *cuts, CCtsp_lpcut_in *c,
        CCtsp_lpcut *new)
{
    int i, j;

    new->cliques = CC_SAFE_MALLOC (c->cliquecount, int);
    if (!new->cliques) return 1;
    new->cliquecount = c->cliquecount;
    for (i = 0; i < c->cliquecount; i++) {
        new->cliques[i] = CCtsp_register_clique (cuts, &c->cliques[i]);
        if (new->cliques[i] == -1) {
            for (j=0; j<i; j++) {
                CCtsp_unregister_clique (cuts, new->cliques[j]);
            }
            CC_FREE (new->cliques, int);
            return 1;
        }
    }
    return 0;
}

void CCtsp_unregister_cliques (CCtsp_lpcuts *cuts, CCtsp_lpcut *c)
{
    int i;

    for (i = 0; i < c->cliquecount; i++) {
        CCtsp_unregister_clique (cuts, c->cliques[i]);
    }
    CC_FREE (c->cliques, int);
    c->cliquecount = 0;
}

int CCtsp_register_dominos (CCtsp_lpcuts *cuts, CCtsp_lpcut_in *c,
        CCtsp_lpcut *new)
{
    int i, j;

    if (c->dominocount > 0 ) {
        new->dominos = CC_SAFE_MALLOC (c->dominocount, int);
        if (!new->dominos) return 1;
        new->dominocount = c->dominocount;
        for (i = 0; i < c->dominocount; i++) {
            new->dominos[i] = CCtsp_register_domino (cuts, &c->dominos[i]);
            if (new->dominos[i] == -1) {
                for (j=0; j<i; j++) {
                    CCtsp_unregister_domino (cuts, new->dominos[j]);
                }
                CC_FREE (new->dominos, int);
                return 1;
            }
        }
    }
    return 0;
}

void CCtsp_unregister_dominos (CCtsp_lpcuts *cuts, CCtsp_lpcut *c)
{
    int i;

    if (c->dominocount > 0 ) {
        for (i = 0; i < c->dominocount; i++) {
            CCtsp_unregister_domino (cuts, c->dominos[i]);
        }
        CC_FREE (c->dominos, int);
        c->dominocount = 0;
    }
}

int CCtsp_add_cut_to_cutlist (CCtsp_lpcuts *cuts, CCtsp_lpcut *c)
{
    if (cuts->cutcount >= cuts->cutspace) {
        if (CCutil_reallocrus_scale ((void **) &cuts->cuts, &cuts->cutspace,
                cuts->cutcount + 1, 1.3, sizeof (CCtsp_lpcut))) {
            return -1;
        }
    }
    cuts->cuts[cuts->cutcount] = *c;
    return cuts->cutcount++;
}

void CCtsp_delete_cut_from_cutlist (CCtsp_lpcuts *cuts, int ind)
{
    int i;

/*
    printf ("CCtsp_delete_cut_from_cutlist (%d)\n",
              cuts->cuts[ind].dominocount);
    fflush (stdout);
*/

    CCtsp_unregister_cliques (cuts, &cuts->cuts[ind]);
    CCtsp_unregister_dominos (cuts, &cuts->cuts[ind]);
    CC_IFFREE (cuts->cuts[ind].mods, CCtsp_sparser);
    CCtsp_free_skeleton (&cuts->cuts[ind].skel);
    for (i = ind+1; i < cuts->cutcount; i++) {
        cuts->cuts[i-1] = cuts->cuts[i];
    }
    cuts->cutcount--;
}
