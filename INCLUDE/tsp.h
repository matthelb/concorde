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
/*                      PROTOTYPES FOR FILES IN TSP                         */
/*                                                                          */
/****************************************************************************/
/****************************************************************************/


#ifndef __TSP_H
#define __TSP_H

#include "util.h"
#include "edgegen.h"
#include "bigguy.h"
#include "lp.h"
#include "cut.h"
#include "kdtree.h"
#include "combs.h"

/*************** Tolerances for the LP and Cutting routines *****************/

#define CCtsp_MIN_VIOL (0.002)    /* min violation for cut to be added to lp */
#define CCtsp_CUTS_DELTA          /* define to set tolerances on ub-lb */
#define CCtsp_CUTS_NEXT_TOL (0.01)         /* to try next level  */
#define CCtsp_CUTS_NEXT_ROUND (0.001)      /* if improve is less, stop round */
#define CCtsp_TENTATIVE_CUTS_NEXT_TOL (0.1)    
#define CCtsp_TENTATIVE_CUTS_NEXT_ROUND (0.01)
#define CCtsp_PRICE_RCTHRESH  (-0.00001)   /* to add a bad edge */
#define CCtsp_PRICE_MAXPENALTY (0.10)      /* penalty permitted in addbad */
#define CCtsp_PHASE1_RCTHRESH (-0.000000001)
#define CCtsp_PHASE1_MAXPENALTY (0.00000001)
#define CCtsp_EDGE_LIFE (1000000) /* 1000000 */      /* 200 */  /* Large for subtour runs */
#define CCtsp_CUT_LIFE  (10)             /* 10 */
#define CCtsp_DUAL_DUST (0.01)           /* 0.001  */
#define CCtsp_EDGE_DUST (0.000001)       /* 0.0001 */

#define CCtsp_CUT_BATCH (250)     /* number of new cuts before lp optimize */
#define CCtsp_STORE_BATCH (250) /* 50 */    /* number of new cuts before lp addrows  */
#define CCtsp_INTTOL (0.0001)     /* used to check if lp soln is integral  */

/************************** Branching Strategies  ***************************/

#define CCtsp_BRANCH_MIDDLE 1
#define CCtsp_BRANCH_STRONG 2

/****************************************************************************/

/************************** Default Communication Ports *********************/

#define CCtsp_HOST_PORT   ((unsigned short) 24846)
#define CCtsp_PROB_PORT   ((unsigned short) 24847)
#define CCtsp_CUT_PORT    ((unsigned short) 24868)
#define CCtsp_DOMINO_PORT ((unsigned short) 24869)

/****************************************************************************/

/************************ Experimental Cutting Planes ***********************/

#undef  CCtsp_USE_DOMINO_CUTS

/****************************************************************************/

#define CCtsp_LP_MAXDOUBLE  1e30

#define CCtsp_COMBRHS(c) (3*(c)->cliquecount - 2)
#define CCtsp_DOMINORHS(c) (3*(c)->dominocount + 1)

typedef struct CCtsp_lpnode {
    int                 deg;
    int                 mark;
    struct CCtsp_lpadj *adj;
} CCtsp_lpnode;

typedef struct CCtsp_lpedge {
    int       ends[2];   /* ends[0] should always be < ends[1] */
    int       fixed;
    int       branch;    /* < 0 means set to 0 and > 0 means set to 1 */
    int       len;
    int       age;
    int       coef;      /* should be maintained at zero */
    int       coefnext;  /* should be maintained at -2 */
} CCtsp_lpedge;

typedef struct CCtsp_lpadj {
    int       to;
    int       edge;
} CCtsp_lpadj;

typedef struct CCtsp_lpgraph {
    int              ncount;
    int              espace;
    int              ecount;
    int              nodemarker;
    CCtsp_lpnode    *nodes;
    CCtsp_lpedge    *edges;
    CCtsp_lpadj     *adjspace;
    int              adjstart;
    int              adjend;
} CCtsp_lpgraph;

typedef struct CCtsp_predge {
    int        ends[2];
    int        len;
    double     rc;
} CCtsp_predge;

typedef struct CCtsp_pricegroup {
    int                    ncount;
    int                    espace;
    int                    ecount;
    CCtsp_lpnode          *nodes;
    CCtsp_predge          *edges;
    int                    cliquecount;
    struct CCtsp_lpclique *cliques; /* just a copy of the pointer */
    CCtsp_lpgraph         *graph;   /* pointer to the copy in a CCtsp_lp */
    CCtsp_lpadj           *adjspace;
    double                *node_pi;
    double                *clique_pi;
    double                 penalty;
} CCtsp_pricegroup;

typedef struct CCtsp_extraedge {
    int       ends[2];
} CCtsp_extraedge;

typedef struct CCtsp_sparser {
    unsigned int node : 24;
    unsigned int mult : 8;
} CCtsp_sparser;

typedef struct CCtsp_segment {
    int lo;
    int hi;
} CCtsp_segment;

typedef struct CCtsp_lpclique {
    int                   segcount;
    struct CCtsp_segment *nodes;
    int                   hashnext;
    int                   refcount;
} CCtsp_lpclique;

typedef struct CCtsp_lpdomino {
    CCtsp_lpclique        sets[2];
    int                   hashnext;
    int                   refcount;
} CCtsp_lpdomino;

#define CC_FOREACH_NODE_IN_CLIQUE(i,c,tmp) \
    for(tmp=0;tmp<(c).segcount;tmp++) \
        for(i=(c).nodes[tmp].lo;i<=(c).nodes[tmp].hi;i++)

typedef struct CCtsp_skeleton {
    int  atomcount;
    int *atoms;
} CCtsp_skeleton;

#define CCtsp_NEWCUT_AGE (-1)

typedef struct CCtsp_lpcut {
    int                   cliquecount;
    int                   dominocount;
    int                   modcount;
    int                   age;
    int                   rhs;
    char                  sense;
    char                  branch;
    int                  *cliques;
    int                  *dominos;
    struct CCtsp_sparser *mods;
    CCtsp_skeleton        skel;
} CCtsp_lpcut;

typedef struct CCtsp_lpcut_in {
    int                    cliquecount;
    int                    dominocount;
    int                    rhs;
    char                   sense;
    char                   branch;
    CCtsp_lpclique        *cliques;
    CCtsp_lpdomino        *dominos;
    CCtsp_skeleton         skel;
    struct CCtsp_lpcut_in *next;
    struct CCtsp_lpcut_in *prev;
} CCtsp_lpcut_in;

typedef struct CCtsp_lp_result {
    double         ub;
    double         lb;
    int            ecount;
    int           *elist;
    double        *x;
    double        *rc;
} CCtsp_lp_result;

typedef struct CCtsp_lpcuts {
    int             cutcount;
    int             savecount;
    int             cliqueend;
    int             cutspace;
    int             cliquespace;
    int             cliquehashsize;
    int             cliquefree;
    int            *cliquehash;
    CCtsp_lpcut    *cuts;
    CCtsp_lpclique *cliques;
    CCgenhash      *cuthash;
    char           *tempcuthash;
    int             tempcuthashsize;
    int             dominoend;
    int             dominospace;
    int             dominohashsize;
    int             dominofree;
    int            *dominohash;
    CCtsp_lpdomino *dominos;
    double         *workloads;
} CCtsp_lpcuts;

typedef struct CCtsp_bigdual {
    int           cutcount;
    CCbigguy     *node_pi;
    CCbigguy     *cut_pi;
} CCtsp_bigdual;

typedef struct CCtsp_tighten_info {
    int    ncall;
    int    nfail;
    int    nadd;
    int    nadd_tied;
    int    ndel;
    int    ndel_tied;
    double add_delta;
    double del_delta;
    double time;
} CCtsp_tighten_info;

typedef struct CCtsp_branchobj {
    int             depth;
    int             rhs;
    int             ends[2];
    char            sense;
    CCtsp_lpclique *clique;
} CCtsp_branchobj;

typedef struct CCtsp_cutnode {
#define CCtsp_CUT_INNODELIST(t) ((t)&4)
#define CCtsp_CUT_ROOT 0
#define CCtsp_CUT_PNODE 1
#define CCtsp_CUT_QNODE 2
#define CCtsp_CUT_LEAF 4
#define CCtsp_CUT_EXTERN 5
    int             type;
    struct CCtsp_cutnode *child;
    struct CCtsp_cutnode *sibling;
    struct CCtsp_cutnode *next;
} CCtsp_cutnode;

typedef struct CCtsp_cuttree {
    int      nodecount;
    int      extern_node;
    CCtsp_cutnode *nodelist;
    CCtsp_cutnode *root;
    CCptrworld cutnode_world;
} CCtsp_cuttree;

/* nodes are reordered to match compression tour */

typedef struct CCtsp_genadj {
    int                     deg;
    struct CCtsp_genadjobj *list;
} CCtsp_genadj;

typedef struct CCtsp_genadjobj {
    int end;
    int len;
} CCtsp_genadjobj;

typedef struct CCtsp_edgegenerator {
    double                    *node_piest;
    struct CCdatagroup        *dg;
    int                       *supply;
    CCkdtree                  *kdtree;
    CCxnear                   *xnear;
    struct CCtsp_xnorm_pricer *xprice;
    CCtsp_genadjobj           *adjobjspace;
    CCtsp_genadj              *adj;
    int                        ncount;
    int                        nneighbors;
    int                        start;
    int                        current;
    int                        supplyhead;
    int                        supplycount;
} CCtsp_edgegenerator;

typedef struct CCtsp_xnorm_pricer_val {
    double                         val;
    struct CCtsp_xnorm_pricer_val *next;
    struct CCtsp_xnorm_pricer_val *prev;
    int                            index;
} CCtsp_xnorm_pricer_val;

typedef struct CCtsp_xnorm_pricer {
    CCdatagroup            *dat;
    double                 *pi;
    int                    *order;
    CCtsp_xnorm_pricer_val *xminuspi_space;
    CCtsp_xnorm_pricer_val *xminuspi;
    int                    *invxminuspi;
    int                     ncount;
} CCtsp_xnorm_pricer;

typedef struct CCtsp_statistics {
    CCutil_timer       cutting_loop;
    CCutil_timer       cutting_inner_loop;
    CCutil_timer       cuts_filecut;
    CCutil_timer       cuts_filecut_opt;
    CCutil_timer       cuts_cutpool;
    CCutil_timer       cuts_cutpool_opt;
    CCutil_timer       cuts_connect;
    CCutil_timer       cuts_connect_opt;
    CCutil_timer       cuts_segment;
    CCutil_timer       cuts_segment_opt;
    CCutil_timer       cuts_remotepool;
    CCutil_timer       cuts_remotepool_opt;
    CCutil_timer       cuts_blockcomb;
    CCutil_timer       cuts_blockcomb_opt;
    CCutil_timer       cuts_growcomb;
    CCutil_timer       cuts_growcomb_opt;
    CCutil_timer       cuts_exactsubtour;
    CCutil_timer       cuts_exactsubtour_opt;
    CCutil_timer       cuts_tighten_lp;
    CCutil_timer       cuts_tighten_lp_opt;
    CCutil_timer       cuts_tighten_lp_close;
    CCutil_timer       cuts_tighten_lp_close_opt;
    CCutil_timer       cuts_decker_lp;
    CCutil_timer       cuts_decker_lp_opt;
    CCutil_timer       cuts_decker_lp_close;
    CCutil_timer       cuts_decker_lp_close_opt;
    CCutil_timer       cuts_star_lp;
    CCutil_timer       cuts_star_lp_opt;
    CCutil_timer       cuts_handling_lp;
    CCutil_timer       cuts_handling_lp_opt;
    CCutil_timer       cuts_cliquetree_lp;
    CCutil_timer       cuts_cliquetree_lp_opt;
    CCutil_timer       cuts_teething_lp;
    CCutil_timer       cuts_teething_lp_opt;
    CCutil_timer       cuts_fastblossom;
    CCutil_timer       cuts_fastblossom_opt;
    CCutil_timer       cuts_ghfastblossom;
    CCutil_timer       cuts_ghfastblossom_opt;
    CCutil_timer       cuts_exactblossom;
    CCutil_timer       cuts_exactblossom_opt;
    CCutil_timer       cuts_tighten_pool;
    CCutil_timer       cuts_tighten_pool_opt;
    CCutil_timer       cuts_decker_pool;
    CCutil_timer       cuts_decker_pool_opt;
    CCutil_timer       cuts_star_pool;
    CCutil_timer       cuts_star_pool_opt;
    CCutil_timer       cuts_handling_pool;
    CCutil_timer       cuts_handling_pool_opt;
    CCutil_timer       cuts_teething_pool;
    CCutil_timer       cuts_teething_pool_opt;
    CCutil_timer       cuts_consecutiveones;
    CCutil_timer       cuts_consecutiveones_opt;
    CCutil_timer       cuts_necklace;
    CCutil_timer       cuts_necklace_opt;
    CCutil_timer       cuts_localcut;
    CCutil_timer       cuts_localcut_opt;

    CCutil_timer       cuts_extraconnect;
    CCutil_timer       cuts_extraconnect_opt;

    CCutil_timer       sparse_edge_check;
    CCutil_timer       full_edge_check;

    CCutil_timer       addcuts;
    CCutil_timer       addcuts_opt;
    CCutil_timer       agecuts;
    CCutil_timer       agecuts_opt;
    CCutil_timer       ageedges;
    CCutil_timer       ageedges_opt;
    CCutil_timer       addbad;
    CCutil_timer       addbad_opt;
    CCutil_timer       strongbranch;
    CCutil_timer       strongbranch_opt;
    CCutil_timer       linkern;

    CCutil_timer       misc;
    CCutil_timer       misc_opt;
    CCutil_timer       total;
    int                problem_cnt;

    CCtsp_tighten_info tighten_stats;
    CCtsp_tighten_info extra_tighten_stats;
} CCtsp_statistics;
    
typedef struct CCtsp_lp {
    CCtsp_lpgraph              graph;
    CCtsp_lpcuts               cuts;
    CCtsp_lpcuts              *pool; 
    CCtsp_lpcuts              *remotepool;
    CCtsp_lpcuts              *dominopool;
    CClp                      *lp;
    int                       *perm;
    CCdatagroup               *dat;
    int                        fullcount;
    CCtsp_genadj              *fulladj;
    CCtsp_genadjobj           *fulladjspace;
    int                        nfixededges;
    int                       *fixededges;
    struct CCtsp_qsparsegroup *sparsifier;
    int                        edge_life;
    int                        cut_life;
    char                      *problabel;
    char                      *probloc;
    int                        id;
    int                        parent_id;
    int                        root;
    double                     upperbound;
    double                     lowerbound;
    CCbigguy                   exact_lowerbound;
    CCtsp_bigdual             *exact_dual;
    int                        infeasible;
    int                        full_edges_valid;
    CClp_warmstart            *warmstart;
    CCtsp_lpcut_in             cutqueue;    /* dummy entry for doubly-linked
                                               list */
    CCtsp_lp_result            result;
    int                        branchdepth;
    CCtsp_branchobj           *branchhistory;
    CCtsp_cuttree              tightcuts;
    CCtsp_statistics           stats;
} CCtsp_lp;

typedef struct CCtsp_lprow {
    int           rowcnt;
    int           nzcnt;
    char         *sense;
    double       *rhs;
    int          *begin;      /* offset into the array for start of row */
    int           indexspace;
    int          *indices;    /* the column indices of the row entries  */
    int           entryspace;
    double       *entries;    /* the matrix entries                     */
} CCtsp_lprow;

typedef struct CCtsp_cutselect {
    int    cutpool;
    int    remotepool;
    char  *remotehost;
    unsigned short remoteport;
    int    domboss;
    char  *dombosshost;
    int    connect;
    int    segments;
    int    exactsubtour;
    int    blockcombs;
    int    growcombs;
    int    prclique;
    int    tighten_lp;
    int    teething_lp;
    int    cliquetree_lp;
    int    tighten_pool;
    int    decker_lp;
    int    decker_pool;
    int    star_lp;
    int    star_pool;
    int    handling_lp;
    int    handling_pool;
    int    maxchunksize;
    int    filecuts;
    char  *filecutname;
    int    teething_pool;
    int    fastblossom;
    int    ghfastblossom;
    int    exactblossom;
    int    consecutiveones;
    int    dominos;
    int    shrunk_dominos;
    int    necklace;
    int    usetighten;     /* set to 1 to tighten before cuts are added */
    int    extra_connect;  /* set to 1 to force a connected solution */
    double nexttol;
    double roundtol;
    int    fastcuts;       /* set to 0 to stop the update of tols */
} CCtsp_cutselect;



/****************************************************************************/
/*                                                                          */
/*                            bcontrol.c                                    */
/*                                                                          */
/****************************************************************************/

#define CCtsp_BBTASK_BRANCH    'b'
#define CCtsp_BBREQ_BRANCHDONE 'B'
#define CCtsp_BBTASK_CUT       'c'
#define CCtsp_BBREQ_CUTDONE    'C'
#define CCtsp_BBREQ_DEADNODE   'D'
#define CCtsp_BBREQ_HELLO      'H'
#define CCtsp_BBREQ_NOBRANCH   'N'
#define CCtsp_BBREQ_TASK       'T'
#define CCtsp_BBREQ_TOUR       'U'
#define CCtsp_BBTASK_WAIT      'w'
#define CCtsp_BBTASK_EXIT      'x'
#define CCtsp_BBREQ_EXIT       'X'

#define CCtsp_BBTASK_TENTATIVE_CUT       'i'
#define CCtsp_BBREQ_TENTATIVE_CUTDONE    'I'
#define CCtsp_BBTASK_TENTATIVE_BRANCH    'j'
#define CCtsp_BBREQ_TENTATIVE_BRANCHDONE 'J'


int
    CCtsp_bfs_brancher (char *probloc, int id, double lowerbound,
        CCtsp_cutselect *sel, CCtsp_cutselect *tsel, double *upbound,
        int *bbcount, int usecliques, CCdatagroup *mydat, int *ptour,
        CCtsp_lpcuts *pool, int ncount, int *besttour, unsigned short hostport,
        double *branchzeit, int save_proof, int tentative_branch_num,
        int longedge_branching, double *timebound, int *hit_timebound,
        int silent, CCrandstate *rstate),
    CCtsp_bfs_restart (char *probloc, char *restart_name, CCtsp_cutselect *sel,
        CCtsp_cutselect *tsel, double *upbound, int *bbcount, int usecliques,
        CCdatagroup *dat, int *ptour, CCtsp_lpcuts *pool, int ncount,
        int *besttour, unsigned short hostport, double *branchzeit,
        int save_proof, int tentative_branch_num, int longedge_branching,
        double *timebound, int *hit_timebound, int silent,
        CCrandstate *rstate),
#ifdef CC_NETREADY
    CCtsp_grunt (char *hostname, unsigned short hostport, char *poolfname,
        char *cutbossname, char *probloc, int silent, 
        CCrandstate *rstate),
#endif /* CC_NETREADY */
    CCtsp_easy_dfs_brancher (CCtsp_lp *lp, CCtsp_cutselect *sel, int depth,
        double *upbound, int *bbcount, int usecliques, int *besttour,
        int longedge_branching, int simple_branching, int silent,
        CCrandstate *rstate),
    CCtsp_do_interactive_branch (CCtsp_lp *lp, int silent, CCrandstate *rstate);



/****************************************************************************/
/*                                                                          */
/*                            blkcomb.c                                     */
/*                                                                          */
/****************************************************************************/


int
    CCtsp_block_combs (CCtsp_lpcut_in **cuts, int *cutcount, int ncount,
        int ecount, int *elist, double *x, int silent);



/****************************************************************************/
/*                                                                          */
/*                            blossom.c                                     */
/*                                                                          */
/****************************************************************************/


int
    CCtsp_fastblossom (CCtsp_lpcut_in **cuts, int *cutcount, int ncount,
        int ecount, int *elist, double *x),
    CCtsp_ghfastblossom (CCtsp_lpcut_in **cuts, int *cutcount, int ncount,
        int ecount, int *elist, double *x),
    CCtsp_exactblossom (CCtsp_lpcut_in **cuts, int *cutcount, int ncount,
        int ecount, int *elist, double *x, CCrandstate *rstate);



/****************************************************************************/
/*                                                                          */
/*                            branch.c                                      */
/*                                                                          */
/****************************************************************************/


int
    CCtsp_find_branch (CCtsp_lp *lp, int nwant, int *ngot,
        CCtsp_branchobj **bobj, double *val, int **cyc, int usecliques,
        int longedge_branching, int silent),
    CCtsp_find_fast_branch (CCtsp_lp *lp, int *ngot, CCtsp_branchobj **bobj,
        double *val, int **cyc, int usecliques, int longedge_branching,
        int silent),
    CCtsp_find_branch_edge (CCtsp_lp *lp, int *n0, int *n1, double *val,
        int **cyc, int branchtype, int silent),
    CCtsp_check_integral (CCtsp_lp *lp, double *val, int **cyc, int *yesno,
        int silent),
    CCtsp_find_branch_cliques (CCtsp_lp *lp, int nwant, int longedge_branching,
        int *ngot, CCtsp_lpclique **bcliques, double **bval, int silent),
    CCtsp_execute_branch (CCtsp_lp *lp, CCtsp_branchobj *b,
        int silent, CCrandstate *rstate),
    CCtsp_execute_unbranch (CCtsp_lp *lp, CClp_warmstart *warmstart,
        int silent, CCrandstate *rstate),
    CCtsp_add_branchhistory_to_lp (CCtsp_lp *lp),
    CCtsp_bb_find_branch (char *probname, int probnum, int ncount,
        CCdatagroup *dat, int *ptour, double *upperbound, CCtsp_lpcuts *pool,
        int nwant, int *ngot, CCtsp_branchobj **b, int usecliques,
        int longedge_branching, int *prune, int *foundtour, int *besttour,
        int silent, CCrandstate *rstate),
    CCtsp_splitprob (CCtsp_lp *lp, CCtsp_branchobj *b, int child0, int child1,
        int silent, CCrandstate *rstate),
    CCtsp_bb_splitprob (char *probname, int probnum, int ncount,
        CCdatagroup *dat, int *ptour, double initial_ub, CCtsp_lpcuts *pool,
        CCtsp_branchobj *b, int child0, int child1, double *val0, double *val1,
        int *prune0, int *prune1, int silent, CCrandstate *rstate),
    CCtsp_dumptour (int ncount, CCdatagroup *dat, int *perm, char *probname,
        int *tour, char *fname, int writeedges, int silent);

void
    CCtsp_init_branchobj (CCtsp_branchobj *b),
    CCtsp_free_branchobj (CCtsp_branchobj *b),
    CCtsp_print_branchhistory (CCtsp_lp *lp);


/****************************************************************************/
/*                                                                          */
/*                             cliqhash.c                                   */
/*                                                                          */
/****************************************************************************/


int
    CCtsp_init_cliquehash (CCtsp_lpcuts *cuts, int size),
    CCtsp_register_clique (CCtsp_lpcuts *cuts, CCtsp_lpclique *c);

unsigned int
    CCtsp_hashclique (CCtsp_lpclique *c);

void
    CCtsp_free_cliquehash (CCtsp_lpcuts *cuts),
    CCtsp_unregister_clique (CCtsp_lpcuts *cuts, int c),
    CCtsp_clique_eq (CCtsp_lpclique *c, CCtsp_lpclique *d, int *yes_no);

int
    CCtsp_init_dominohash (CCtsp_lpcuts *cuts, int size),
    CCtsp_register_domino (CCtsp_lpcuts *cuts, CCtsp_lpdomino *c);

unsigned int
    CCtsp_hashdomino (CCtsp_lpdomino *d);

void
    CCtsp_free_dominohash (CCtsp_lpcuts *cuts),
    CCtsp_domino_eq (CCtsp_lpdomino *c, CCtsp_lpdomino *d, int *yes_no),
    CCtsp_unregister_domino (CCtsp_lpcuts *cuts, int c);



/****************************************************************************/
/*                                                                          */
/*                           cliqwork.c                                     */
/*                                                                          */
/****************************************************************************/

typedef struct CCtsp_cutinfo {
    CC_SRKexpinfo    expand;
    CCtsp_lpcut_in **clist;
    CCtsp_lpcut_in  *current;
    int             *cutcount;
} CCtsp_cutinfo;


int
    CCtsp_clique_to_array (CCtsp_lpclique *c, int **ar, int *count),
    CCtsp_clique_delta (CCtsp_lpgraph *g, double *x, CCtsp_lpclique *c,
        double *delta),
    CCtsp_copy_lpcut_in (CCtsp_lpcut_in *c, CCtsp_lpcut_in *new),
    CCtsp_segment_to_subtour (CCtsp_lpcut_in **cut, int a, int b, int ncount),
    CCtsp_array_to_subtour (CCtsp_lpcut_in **cut, int *ar, int acount,
        int ncount),
    CCtsp_array_to_lpclique (int *ar, int acount, CCtsp_lpclique *cliq),
    CCtsp_seglist_to_lpclique (int nseg, int *list, CCtsp_lpclique *cliq),
    CCtsp_shrunk_set_to_lpclique (int cnt, int *set, int *wset,
        CC_SRKexpinfo *expand, CCtsp_lpclique *cliq),
    CCtsp_add_nodes_to_lpclique (CCtsp_lpclique *cin, CCtsp_lpclique *cout,
         int addcount, int *adda),
    CCtsp_delete_nodes_from_lpclique (CCtsp_lpclique *cin,
         CCtsp_lpclique *cout, int delcount, int *del),
    CCtsp_lpcut_to_lpcut_in (CCtsp_lpcuts *cuts, CCtsp_lpcut *c,
        CCtsp_lpcut_in *new),
    CCtsp_copy_lpclique (CCtsp_lpclique *c, CCtsp_lpclique *new),
    CCtsp_copy_lpdomino (CCtsp_lpdomino *c, CCtsp_lpdomino *new),
    CCtsp_create_lpcliques (CCtsp_lpcut_in *c, int cliquecount),
    CCtsp_max_node (CCtsp_lpcut_in *c),
    CCtsp_build_dp_cut (CCtsp_lpcut_in **cut, int ndomino, int *Acount,
        int **A, int *Bcount, int **B, int handlecount, int *handle);

void
    CCtsp_mark_clique (CCtsp_lpclique *c, int *marks, int marker),
    CCtsp_mark_domino (CCtsp_lpdomino *c, int *marks, int marker),
    CCtsp_mark_clique_and_neighbors (CCtsp_lpgraph *g, CCtsp_lpclique *c,
        int *marks, int marker),
    CCtsp_mark_domino_and_neighbors (CCtsp_lpgraph *g, CCtsp_lpdomino *c,
        int *marks, int marker),
    CCtsp_mark_clique_and_neighbors_double (CCtsp_lpgraph *g,
        CCtsp_lpclique *c, double *marks, double marker),
    CCtsp_mark_cut (CCtsp_lpcut_in *c, int *marks, int marker),
    CCtsp_mark_cut_and_neighbors (CCtsp_lpgraph *g, CCtsp_lpcut_in *c,
        int *marks, int marker),
    CCtsp_is_clique_marked (CCtsp_lpclique *c, int *marks, int marker,
        int *yes_no),
    CCtsp_clique_count (CCtsp_lpclique *c, int *count),
    CCtsp_clique_marked_count (CCtsp_lpclique *c, int *marks, int marker,
         int *count),
    CCtsp_init_lpcut_in (CCtsp_lpcut_in *c),
    CCtsp_init_lpcut (CCtsp_lpcut *c),
    CCtsp_init_lpclique (CCtsp_lpclique *c),
    CCtsp_init_lpdomino (CCtsp_lpdomino *c),
    CCtsp_print_lpcut_in (CCtsp_lpcut_in *c),
    CCtsp_print_lpclique (CCtsp_lpclique *c),
    CCtsp_print_lpdomino (CCtsp_lpdomino *d),
    CCtsp_lpclique_compare (CCtsp_lpclique *a, CCtsp_lpclique *b, int *diff);



/****************************************************************************/
/*                                                                          */
/*                            control.c                                     */
/*                                                                          */
/****************************************************************************/


int
    CCtsp_cutting_multiple_loop (CCtsp_lp *lp, CCtsp_cutselect *sel,
        int savelp, int maxlocal, int update_tol, int silent,
        CCrandstate *rstate),
    CCtsp_cutting_loop (CCtsp_lp *lp, CCtsp_cutselect *sel, int savelp,
        int silent, CCrandstate *rstate),
    CCtsp_subtour_loop (CCtsp_lp *lp, int silent, CCrandstate *rstate),
    CCtsp_blossom_loop (CCtsp_lp *lp, int silent, CCrandstate *rstate),
    CCtsp_subtour_and_blossom_loop (CCtsp_lp *lp, int silent,
        CCrandstate *rstate),
    CCtsp_pricing_loop (CCtsp_lp *lp, double *bnd, int silent,
        CCrandstate *rstate),
    CCtsp_call_x_heuristic (CCtsp_lp *lp, double *val, int *outcyc,
        int silent, CCrandstate *rstate),
    CCtsp_bb_cutting (char *probname, int probnum, int prob_newnum, int ncount,
        CCdatagroup *dat, int *ptour, double *upbound, CCtsp_lpcuts *pool,
        CCtsp_cutselect *sel, double *val, int *prune, int *foundtour,
        int *besttour, int level, int silent, CCrandstate *rstate),
    CCtsp_cutselect_set_tols (CCtsp_cutselect *s, CCtsp_lp *lp, int level,
        int silent);

void
    CCtsp_init_cutselect (CCtsp_cutselect *s),
    CCtsp_cutselect_dominos (CCtsp_cutselect *s, int domsel),
    CCtsp_cutselect_tighten (CCtsp_cutselect *s, int tighten),
    CCtsp_cutselect_chunksize (CCtsp_cutselect *s, int chunksize),
    CCtsp_cutselect_filecuts (CCtsp_cutselect *s, char *fname),
    CCtsp_cutselect_remotepool (CCtsp_cutselect *s, char *cutbossname),
    CCtsp_cutselect_domboss (CCtsp_cutselect *s, char *dombossname),
    CCtsp_init_tentative_cutselect (CCtsp_cutselect *s),
    CCtsp_init_simple_cutselect (CCtsp_cutselect *s),
    CCtsp_init_fast_cutselect (CCtsp_cutselect *s);


/****************************************************************************/
/*                                                                          */
/*                             cutcall.c                                    */
/*                                                                          */
/****************************************************************************/


int
    CCtsp_connect_cuts (CCtsp_lpcut_in **cuts, int *cutcount, int ncount,
        int ecount, int *elist, double *x),
    CCtsp_segment_cuts (CCtsp_lpcut_in **cuts, int *cutcount, int ncount,
        int ecount, int *elist, double *x),
    CCtsp_shrink_subtours (CCtsp_lpcut_in **cuts, int *cutcount, int ncount,
        int ecount, int *elist, double *x),
    CCtsp_exact_subtours (CCtsp_lpcut_in **cuts, int *cutcount, int ncount,
        int ecount, int *elist, double *x),
    CCtsp_tighten_lp (CCtsp_lpcuts *cuts, CCtsp_tighten_info *stats,
        CCtsp_lpcut_in **cutsout, int *cutcount, int ncount, int ecount,
        int *elist, double *x, double testtol, int maxcuts,
        double *viol, CCrandstate *rstate),
    CCtsp_double_decker_lp (CCtsp_lpcuts *cuts, CCtsp_tighten_info *stats,
        CCtsp_lpcut_in **cutsout, int *cutcount, int ncount, int ecount,
        int *elist, double *x, double testtol, int maxcuts,
        double *viol, CCrandstate *rstate),
    CCtsp_cliquetree_lp (CCtsp_lpcuts *cuts, CCtsp_tighten_info *stats,
        CCtsp_lpcut_in **cutsout, int *cutcount, int ncount, int ecount,
        int *elist, double *x, double testtol, int maxcuts,
        double *viol, CCrandstate *rstate),
    CCtsp_star_lp (CCtsp_lpcuts *cuts, CCtsp_tighten_info *stats,
        CCtsp_lpcut_in **cutsout, int *cutcount, int ncount, int ecount,
        int *elist, double *x, double testtol, int maxcuts,
        double *viol, CCrandstate *rstate),
    CCtsp_handling_lp (CCtsp_lpcuts *cuts, CCtsp_tighten_info *stats,
        CCtsp_lpcut_in **cutsout, int *cutcount, int ncount, int ecount,
        int *elist, double *x, double testtol, int maxcuts,
        double *viol, CCrandstate *rstate),
    CCtsp_teething_lp (CCtsp_lpcuts *cuts, CCtsp_tighten_info *stats,
        CCtsp_lpcut_in **cutsout, int *cutcount, int ncount, int ecount,
        int *elist, double *x, double testtol, int maxcuts,
        double *viol, CCrandstate *rstate),
    CCtsp_domino_trial (CCtsp_lpcut_in **cuts, int *cutcount, int ncount,
        int ecount, int *elist, double *x, CCrandstate *rstate),
    CCtsp_file_cuts (char *cutfile, CCtsp_lpcut_in **cuts, int *cutcount,
        int ncount, int *tour),
    CCtsp_file_cuts_write (const char *cutfile, CCtsp_lpcuts *cuts, int *tour),
    CCtsp_test_pure_comb (int ncount, CCtsp_lpcut_in *c, int *yes_no,
        int *handle),
    CCtsp_test_pseudocomb (int ncount, CCtsp_lpcut_in *c, int handle,
        int *yes_no),
    CCtsp_test_teeth_disjoint (int ncount, CCtsp_lpcut_in *c, int handle,
        int *yes_no),
    CCtsp_find_pure_handle (int ncount, CCtsp_lpcut_in *c, int *handle),
    CCtsp_truncate_cutlist (CCtsp_lpcut_in **cuts, int ncount, int ecount,
        int *elist, double *x, int maxcuts, CCrandstate *rstate),
    CCtsp_buildcut_begin (CCtsp_cutinfo *cuts, int init_cliquecount),
    CCtsp_buildcut_addclique (CCtsp_cutinfo *cuts, int *arr, int size),
    CCtsp_buildcut_finish (CCtsp_cutinfo *cuts, int rhs),
    CCtsp_new_domino (CCtsp_lpcut_in **cuts, int *cutcount, int ncount,
        int ecount, int *elist, double *x, const char *bossname),
    CCtsp_shrink_domino (CCtsp_lpcut_in **cuts, int *cutcount, int ncount,
        int ecount, int *elist, double *x, int quickshrink, int rand_minor,
        CCrandstate *rstate, const char *bossname);

void
    CCtsp_buildcut_abort (CCtsp_cutinfo *cuts);



/****************************************************************************/
/*                                                                          */
/*                            cutpool.c                                     */
/*                                                                          */
/****************************************************************************/

#define CCtsp_POOL_GETCUTS     'G'
#define CCtsp_POOL_PUTCUTS     'P'
#define CCtsp_POOL_SAVECUTS    'S'
#define CCtsp_POOL_EXIT        'X'


int
    CCtsp_init_cutpool (int *ncount, char *poolfilename, CCtsp_lpcuts **pool),
    CCtsp_write_cutpool (int ncount, const char *poolfilename,
        CCtsp_lpcuts  *pool),
    CCtsp_search_cutpool (CCtsp_lpcuts *pool, CCtsp_lpcut_in **cuts,
        int *cutcount, double *maxviol, int ncount, int ecount, int *elist,
        double *x, int nthreads, CCrandstate *rstate),
    CCtsp_search_remotepool (char *remotehost, unsigned short remoteport,
        CCtsp_lpcut_in **cuts, int *cutcount, double *maxviol, int ncount,
        int ecount, int *elist, double *x),
    CCtsp_read_cuts (CC_SFILE *f, int *ncount, CCtsp_lpcuts *cuts,
        int readmods, int buildhash),
    CCtsp_read_lpcut_in (CC_SFILE *f, CCtsp_lpcut_in *c, int ncount),
    CCtsp_read_lpclique (CC_SFILE *f, CCtsp_lpclique *c, int ncount),
    CCtsp_read_lpdomino (CC_SFILE *f, CCtsp_lpdomino *d, int ncount),
    CCtsp_write_cuts (CC_SFILE *f, int ncount, CCtsp_lpcuts *cuts,
        int writemods),
    CCtsp_send_newcuts (int ncount, CCtsp_lpcuts *pool, char *remotehost,
        unsigned short remoteport),
    CCtsp_write_lpcut_in (CC_SFILE *f, CCtsp_lpcut_in *c, int ncount),
    CCtsp_write_lpcut (CC_SFILE *f, CCtsp_lpcuts *cuts, CCtsp_lpcut *c,
        int ncount),
    CCtsp_write_lpclique (CC_SFILE *f, CCtsp_lpclique *c, int ncount),
    CCtsp_write_lpdomino (CC_SFILE *f, CCtsp_lpdomino *c, int ncount),
    CCtsp_copy_cuts (CC_SFILE *f, CC_SFILE *t, int copymods),
    CCtsp_search_cutpool_cliques (CCtsp_lpcuts *pool, CCtsp_lpclique **cliques,
        int *cliquecount, int ncount, int ecount, int *elist, double *x,
        double maxdelta, int maxcliques, double **cliquevals,
        CCrandstate *rstate),
    CCtsp_branch_cutpool_cliques (CCtsp_lpcuts *pool, CCtsp_lpclique **cliques,
        int *cliquecount, int ncount, int ecount, int *elist, double *x,
        int nwant, double **cliquevals, int silent),
    CCtsp_get_clique_prices (CCtsp_lpcuts *pool, int **p_cliquenums,
        double **p_cliquevals, double mindelta, double maxdelta,
        int *p_cliquecount, int ncount, int ecount, int *elist, double *x),
    CCtsp_get_clique (CCtsp_lpcuts *pool, int cliquenum,
        CCtsp_lpclique **p_clique),
    CCtsp_add_to_cutpool (CCtsp_lpcuts *pool, CCtsp_lpcuts *cuts,
        CCtsp_lpcut *c),
    CCtsp_add_to_dominopool (CCtsp_lpcuts *pool, CCtsp_lpcuts *cuts,
        CCtsp_lpcut *c),
    CCtsp_add_to_cutpool_lpcut_in (CCtsp_lpcuts *pool, CCtsp_lpcut_in *cut),
    CCtsp_display_cutpool (CCtsp_lpcuts *pool),
    CCtsp_price_cuts (CCtsp_lpcuts *pool, int ncount, int ecount, int *elist,
        double *x, double *cutval),
    CCtsp_price_cuts_threaded (CCtsp_lpcuts *pool, int ncount, int ecount,
        int *elist, double *x, double *cutval, int numthreads),
    CCtsp_register_cliques (CCtsp_lpcuts *cuts, CCtsp_lpcut_in *c,
        CCtsp_lpcut *new),
    CCtsp_register_dominos (CCtsp_lpcuts *cuts, CCtsp_lpcut_in *c,
        CCtsp_lpcut *new),
    CCtsp_add_cut_to_cutlist (CCtsp_lpcuts *cuts, CCtsp_lpcut *c);

void
    CCtsp_free_cutpool (CCtsp_lpcuts **pool),
    CCtsp_free_lpcut_in (CCtsp_lpcut_in *c),
    CCtsp_free_lpclique (CCtsp_lpclique *c),
    CCtsp_free_lpdomino (CCtsp_lpdomino *c),
    CCtsp_unregister_cliques (CCtsp_lpcuts *cuts, CCtsp_lpcut *c),
    CCtsp_unregister_dominos (CCtsp_lpcuts *cuts, CCtsp_lpcut *c),
    CCtsp_delete_cut_from_cutlist (CCtsp_lpcuts *cuts, int ind);


/****************************************************************************/
/*                                                                          */
/*                            ddecker.c                                     */
/*                                                                          */
/****************************************************************************/


int
    CCtsp_test_pure_double_decker (CCtsp_lpcut_in *c, int *yes_no,
        int *handle1, int *handle2),
    CCtsp_comb_to_double_decker (CCtsp_lpgraph *g, CC_GCgraph *h,
        double *x, CCtsp_lpcut_in *c, CCtsp_lpcut_in **d),
    CCtsp_comb_to_star (CCtsp_lpgraph *g, CC_GCgraph *h, double *x,
        CCtsp_lpcut_in *c, CCtsp_lpcut_in **d),
    CCtsp_test_pure_simple_cliquetree (int ncount, CCtsp_lpcut_in *c,
       int *yes_no),
    CCtsp_comb_to_cliquetree (CCtsp_lpgraph *g, CC_GCgraph *h, double *x,
        CCtsp_lpcut_in *c, CCtsp_lpcut_in **d),
    CCtsp_comb_handling (CCtsp_lpgraph *g, CC_GCgraph *h, double *x,
        CCtsp_lpcut_in *c, CCtsp_lpcut_in **d);



/****************************************************************************/
/*                                                                          */
/*                            ex_price.c                                    */
/*                                                                          */
/****************************************************************************/


int
    CCtsp_exact_price (CCtsp_lp *lp, CCbigguy *bound, int complete_price,
        int phase1, int silent),
    CCtsp_edge_elimination (CCtsp_lp *lp, int eliminate_sparse, int silent),
    CCtsp_exact_dual (CCtsp_lp *lp),
    CCtsp_verify_infeasible_lp (CCtsp_lp *lp, int *yesno, int silent),
    CCtsp_verify_lp_prune (CCtsp_lp *lp, int *yesno, int silent);

void
    CCtsp_free_bigdual (CCtsp_bigdual **d);


/****************************************************************************/
/*                                                                          */
/*                             generate.c                                   */
/*                                                                          */
/****************************************************************************/


#define CCtsp_PRICE_COMPLETE_GRAPH -1
#define CCtsp_GEN_PRICE_EPSILON 0.0001 /* 0.0000001 */
#define CCtsp_GEN_USE_ADJ 50           /* Cutoff for using explicit adj list */


void
    CCtsp_free_edgegenerator (CCtsp_edgegenerator *eg);

int
    CCtsp_init_edgegenerator (CCtsp_edgegenerator *eg, int ncount,
        CCdatagroup *dg, CCtsp_genadj *adj, int nneighbors,
        int silent, CCrandstate *rstate),
    CCtsp_reset_edgegenerator (CCtsp_edgegenerator *eg, double *node_piest,
        int silent),
    CCtsp_generate_edges (CCtsp_edgegenerator *eg, int nwant, int *pngot,
        int *elist, int *elen, int *finished, int silent, CCrandstate *rstate),
    CCtsp_edgelist_to_genadj (int ncount, int ecount, int *elist, int *elen,
        CCtsp_genadj **adj, CCtsp_genadjobj **adjobjspace);



/****************************************************************************/
/*                                                                          */
/*                            growcomb.c                                    */
/*                                                                          */
/****************************************************************************/


int
    CCtsp_edge_comb_grower (CCtsp_lpcut_in **cuts, int *cutcount, int ncount,
        int ecount, int *elist, double *x, CCtsp_tighten_info *stats);



/****************************************************************************/
/*                                                                          */
/*                            prclique.c                                    */
/*                                                                          */
/****************************************************************************/


int
    CCtsp_pr_cliquetree (CCtsp_lpcut_in **cuts, int *cutcount, int ncount,
        int ecount, int *elist, double *x, CCtsp_tighten_info *stats);



/****************************************************************************/
/*                                                                          */
/*                             prob_io.c                                    */
/*                                                                          */
/****************************************************************************/

#define CCtsp_PROB_FILE_NAME_LEN 128

#define CCtsp_Pdelete    'D'
#define CCtsp_Pread      'R'
#define CCtsp_Pwrite     'W'
#define CCtsp_Pmaster    'M'
#define CCtsp_Pexit      'X'
#define CCtsp_Pcuts      'c'
#define CCtsp_Pdual      'd'
#define CCtsp_Pedges     'e'
#define CCtsp_Pfixed     'f'
#define CCtsp_Pfull      'g'
#define CCtsp_Pheader    'h'
#define CCtsp_Phistory   'i'
#define CCtsp_Ptour      't'
#define CCtsp_Pwarmstart 'w'

typedef struct CCtsp_PROB_FILE {
    CC_SFILE *f;
    int type;
    char name[CCtsp_PROB_FILE_NAME_LEN];
    int id;
    int parent;
    double ub;
    double lb;
    CCbigguy exactlb;
    int nnodes;
    int child0;
    int child1;
    int real;       /* Set to 1 when we know this is a real child */
    int processed;
    int infeasible;
    struct {
        int dat;
        int edge;
        int fulladj;
        int cut;
        int tour;
        int basis;  /* obsolete - replaced by warmstart */
        int norms;  /* obsolete - replaced by warmstart */
        int fix;
        int exactdual;
        int history;
        int warmstart;
    } offsets;
} CCtsp_PROB_FILE;


CCtsp_PROB_FILE
    *CCtsp_prob_read (char *f, int n),
    *CCtsp_prob_read_name (char *f),
    *CCtsp_prob_write (char *f, int n),
    *CCtsp_prob_write_name (char *fname);

int
    CCtsp_prob_file_delete (char *f, int n),
    CCtsp_prob_getname (CCtsp_PROB_FILE *p, char *name),
    CCtsp_prob_getid (CCtsp_PROB_FILE *p, int *id),
    CCtsp_prob_getparent (CCtsp_PROB_FILE *p, int *parent),
    CCtsp_prob_getub (CCtsp_PROB_FILE *p, double *ub),
    CCtsp_prob_getlb (CCtsp_PROB_FILE *p, double *lb),
    CCtsp_prob_getexactlb (CCtsp_PROB_FILE *p, CCbigguy *lb),
    CCtsp_prob_getnnodes (CCtsp_PROB_FILE *p, int *nnodes),
    CCtsp_prob_getchildren (CCtsp_PROB_FILE *p, int *child0, int *child1),
    CCtsp_prob_getreal (CCtsp_PROB_FILE *p, int *real),
    CCtsp_prob_getprocessed (CCtsp_PROB_FILE *p, int *processed),
    CCtsp_prob_getinfeasible (CCtsp_PROB_FILE *p, int *infeasible),
    CCtsp_prob_gettour (CCtsp_PROB_FILE *p, int ncount, int **tour, int silent),
    CCtsp_prob_getedges (CCtsp_PROB_FILE *p, int ncount, int *nedges,
        int **elist, int **elen, int silent),
    CCtsp_prob_getcuts (CCtsp_PROB_FILE *p, int *ncount, CCtsp_lpcuts *cuts,
        int silent),
    CCtsp_prob_getwarmstart (CCtsp_PROB_FILE *p, CClp_warmstart **w,
        int silent),
    CCtsp_prob_getfulladj (CCtsp_PROB_FILE *p, int ncount, int *fullcount,
        CCtsp_genadj **adj, CCtsp_genadjobj **adjspace, int silent),
    CCtsp_prob_getfixed (CCtsp_PROB_FILE *p, int ncount, int *ecount,
        int **elist, int silent),
    CCtsp_prob_getexactdual (CCtsp_PROB_FILE *p, int ncount,
        CCtsp_bigdual **d, int silent),
    CCtsp_prob_gethistory (CCtsp_PROB_FILE *p, int *depth,
        CCtsp_branchobj **history, int silent),
    CCtsp_prob_rclose (CCtsp_PROB_FILE *p),
    CCtsp_prob_putname (CCtsp_PROB_FILE *p, char *name),
    CCtsp_prob_putid (CCtsp_PROB_FILE *p, int id),
    CCtsp_prob_putparent (CCtsp_PROB_FILE *p, int parent),
    CCtsp_prob_putub (CCtsp_PROB_FILE *p, double ub),
    CCtsp_prob_putlb (CCtsp_PROB_FILE *p, double lb),
    CCtsp_prob_putexactlb (CCtsp_PROB_FILE *p, CCbigguy lb),
    CCtsp_prob_putnnodes (CCtsp_PROB_FILE *p, int nnodes),
    CCtsp_prob_putchildren (CCtsp_PROB_FILE *p, int child0, int child1),
    CCtsp_prob_putreal (CCtsp_PROB_FILE *p, int real),
    CCtsp_prob_putprocessed (CCtsp_PROB_FILE *p, int processed),
    CCtsp_prob_putinfeasible (CCtsp_PROB_FILE *p, int infeasible),
    CCtsp_prob_puttour (CCtsp_PROB_FILE *p, int ncount, int *tour),
    CCtsp_prob_putedges (CCtsp_PROB_FILE *p, int ncount, int nedges,
        int *elist, int *elen),
    CCtsp_prob_putcuts (CCtsp_PROB_FILE *p, int ncount, CCtsp_lpcuts *cuts),
    CCtsp_prob_putwarmstart (CCtsp_PROB_FILE *p, CClp_warmstart *w),
    CCtsp_prob_putfulladj (CCtsp_PROB_FILE *p, int ncount, int fullcount,
        CCtsp_genadj *adj),
    CCtsp_prob_putfixed (CCtsp_PROB_FILE *p, int ncount, int ecount,
        int *elist),
    CCtsp_prob_putexactdual (CCtsp_PROB_FILE *p, CCtsp_bigdual *d, int ncount),
    CCtsp_prob_puthistory (CCtsp_PROB_FILE *p, int depth,
        CCtsp_branchobj *history),
    CCtsp_prob_wclose (CCtsp_PROB_FILE *p),
    CCtsp_prob_copy_section (CCtsp_PROB_FILE *f, CCtsp_PROB_FILE *t,
        char section, int silent);

char
   *CCtsp_problabel (const char *probloc);

#ifdef CC_NETREADY
CCtsp_PROB_FILE
   *CCtsp_prob_read_remote (char *hname, char *pname, int n),
   *CCtsp_prob_write_remote (char *hname, char *pname, int n),
   *CCtsp_prob_server (CC_SFILE *s);

int
    CCtsp_prob_delete_remote (char *hname, char *pname, int n);
#endif /* CC_NETREADY */




/****************************************************************************/
/*                                                                          */
/*                             qsparse.c                                    */
/*                                                                          */
/****************************************************************************/

typedef struct CCtsp_qsparsegroup {
    CCdheap *add_queue;   /* An empty heap will be maintained */
    CCdheap *sub_queue;   /* An empty heap will be maintained */
    int *count_m1;        /* The array will be maintained at 0 */
    int *count_non0;      /* The array will be maintained at 0 */
    int *count_1;         /* The array will be maintained at 0 */
    int *on_add_queue;    /* The array will be maintained at 0 */
    int *on_sub_queue;    /* The array will be maintained at 0 */
    int *mults;           /* The array will be maintained at 0 */
} CCtsp_qsparsegroup;


void
    CCtsp_free_qsparsify (CCtsp_qsparsegroup **pqs);

int
    CCtsp_qsparsify (CCtsp_qsparsegroup **pqs, struct CCtsp_lpgraph *g,
        int *pnzlist, int *scount, struct CCtsp_sparser **slist,
        int *savedcount);


/****************************************************************************/
/*                                                                          */
/*                           skeleton.c                                     */
/*                                                                          */
/****************************************************************************/


int
    CCtsp_copy_skeleton (CCtsp_skeleton *old, CCtsp_skeleton *new),
    CCtsp_construct_skeleton (CCtsp_lpcut_in *c, int nodecount),
    CCtsp_read_skeleton (CC_SFILE *f, CCtsp_skeleton *skel, int ncount),
    CCtsp_write_skeleton (CC_SFILE *f, CCtsp_skeleton *skel, int ncount);

void
    CCtsp_init_skeleton (CCtsp_skeleton *skel),
    CCtsp_free_skeleton (CCtsp_skeleton *skel),
    CCtsp_compare_skeletons (CCtsp_skeleton *a, CCtsp_skeleton *b, int *diff);



/****************************************************************************/
/*                                                                          */
/*                           teething.c                                     */
/*                                                                          */
/****************************************************************************/


int
    CCtsp_teething (CCtsp_lpgraph *g, double *x, CCtsp_lpcut_in *cut,
        CCtsp_lpcut_in **newcut),
    CCtsp_teething_list (CCtsp_lpgraph *g, double *x, CCtsp_lpclique *handle,
        int nbig, CCtsp_lpclique **bigteeth, CCtsp_lpcut_in **newcut);



/****************************************************************************/
/*                                                                          */
/*                           tighten.c                                      */
/*                                                                          */
/****************************************************************************/


int
    CCtsp_tighten_lpcut_in (CCtsp_lpgraph *g, CCtsp_lpcut_in *c, double *x,
        CCtsp_lpcut_in *d, CCtsp_tighten_info *stats, double *pimprove),
    CCtsp_tighten_lpcut (CCtsp_lpgraph *g, CCtsp_lpclique *cliques,
        CCtsp_lpcut *c, double *x, CCtsp_lpcut_in *d,
        CCtsp_tighten_info *stats, double *pimprove);

void
    CCtsp_init_tighten_info (CCtsp_tighten_info *stats),
    CCtsp_print_tighten_info (CCtsp_tighten_info *stats);


/****************************************************************************/
/*                                                                          */
/*                            tsp_lp.c                                      */
/*                                                                          */
/****************************************************************************/


int
    CCtsp_bb_init_lp (CCtsp_lp **lp, char *probname, int probnum, int ncount,
        CCdatagroup *dat, int *ptour, double initial_ub, CCtsp_lpcuts *pool,
        int silent, CCrandstate *rstate),
    CCtsp_init_lp (CCtsp_lp **lp, char *probname, int probnum,
        char *probfilename, int ncount, CCdatagroup *dat, int ecount,
        int *elist, int *elen, int excount, int *exlist, int *exlen,
        int exvalid, int *ptour, double initial_ub, CCtsp_lpcuts *pool,
        CCtsp_lpcuts *dominopool, int silent, CCrandstate *rstate),
    CCtsp_build_lpgraph (CCtsp_lpgraph *g, int ncount, int ecount,
        int *elist, int *elen),
    CCtsp_build_lpadj (CCtsp_lpgraph *g, int estart, int eend),
    CCtsp_find_edge (CCtsp_lpgraph *g, int from, int to),
    CCtsp_inspect_full_edges (CCtsp_lp *lp),
    CCtsp_resparsify_lp (CCtsp_lp *lp, int silent),
    CCtsp_lpcut_nzlist (CCtsp_lpgraph *g, CCtsp_lpcut *c,
        CCtsp_lpclique *cliques, CCtsp_lpdomino *dominos, int do_mods),
    CCtsp_update_result (CCtsp_lp *lp),
    CCtsp_get_lp_result (CCtsp_lp *lp, double *lb, double *ub, int *ecount,
        int **elist, double **x, double **rc, double **node_pi,
        double **cut_pi),
    CCtsp_lpcut_in_nzlist (CCtsp_lpgraph *g, CCtsp_lpcut_in *c),
    CCtsp_process_cuts (CCtsp_lp *lp, int *pnadded, int tighten,
        int silent, CCrandstate *rstate),
    CCtsp_infeas_recover (CCtsp_lp *lp, int silent, CCrandstate *rstate),
    CCtsp_add_cut (CCtsp_lp *lp, CCtsp_lpcut_in *d, CCtsp_lprow *cr),
    CCtsp_add_nzlist_to_lp (CCtsp_lp *lp, int nzlist, int rhs, char sense,
        CCtsp_lprow *cr),
    CCtsp_addbad_variables (CCtsp_lp *lp, CCtsp_edgegenerator *eg,
        double *ppenalty, int *pnadded, double rcthresh,
        double maxpenalty, int phase1, int *feasible, int silent,
        CCrandstate *rstate),
    CCtsp_eliminate_variables (CCtsp_lp *lp, int eliminate_sparse, int silent),
    CCtsp_add_vars_to_lp (CCtsp_lp *lp, CCtsp_predge *prlist, int n),
    CCtsp_add_multiple_rows (CCtsp_lp *lp, CCtsp_lprow *cr),
    CCtsp_delete_cut (CCtsp_lp *lp, int i),
    CCtsp_reduced_cost_nearest (CCtsp_lp *lp, int k, int *ecount, int **elist,
        double **elen, int sparse),
    CCtsp_write_probfile_sav (CCtsp_lp *lp),
    CCtsp_write_probfile_id (CCtsp_lp *lp),
    CCtsp_write_probroot_id (char *probloc, CCtsp_lp *lp),
    CCtsp_write_probleaf_id (CCtsp_lp *lp),
    CCtsp_read_probfile (CCtsp_lp *lp, char *fname, char *probloc,
        int *ncount, int silent),
    CCtsp_read_probfile_id (CCtsp_lp *lp, char *fname, int id, int *ncount,
        int silent),
    CCtsp_dump_rc_nearest (CCtsp_lp *lp, int k, char *fname, int sparse),
    CCtsp_dump_x (CCtsp_lp *lp, char *fname),
    CCtsp_depot_valid (CCtsp_lp *lp, int ndepot, int *yesno);

double
    CCtsp_cutprice (CCtsp_lpgraph *g, CCtsp_lpcut_in *c, double *x);

void
    CCtsp_init_tsp_lpcuts_struct (CCtsp_lpcuts *c),
    CCtsp_init_tsp_lp_struct (CCtsp_lp *lp),
    CCtsp_free_tsp_lp_struct (CCtsp_lp **lp),
    CCtsp_init_lpgraph_struct (CCtsp_lpgraph *g),
    CCtsp_free_lpgraph (CCtsp_lpgraph *g),
    CCtsp_init_statistics (CCtsp_statistics *stats),
    CCtsp_output_statistics (CCtsp_statistics *stats),
    CCtsp_add_cuts_to_queue (CCtsp_lp *lp, CCtsp_lpcut_in **c),
    CCtsp_init_lprow (CCtsp_lprow *cr),
    CCtsp_free_lprow (CCtsp_lprow *cr);


/****************************************************************************/
/*                                                                          */
/*                            tsp_lp.c                                      */
/*                                                                          */
/****************************************************************************/

int
    CCtsp_solve_sparse (int ncount, int ecount, int *elist, int *elen,
        int *in_tour, int *out_tour, double *in_val, double *optval,
        int *success, int *foundtour, char *name, double *timebound,
        int *hit_timebound, int silent, CCrandstate *rstate),
    CCtsp_solve_dat (int ncount, CCdatagroup *indat, int *in_tour,
        int *out_tour, double *in_val, double *optval, int *success,
        int *foundtour, char *name, double *timebound, int *hit_timebound,
        int silent, CCrandstate *rstate);



/****************************************************************************/
/*                                                                          */
/*                             xtour.c                                      */
/*                                                                          */
/****************************************************************************/


int
    CCtsp_x_greedy_tour (CCdatagroup *dat, int ncount, int ecount, int *elist,
        double *x, int *cyc, double *val, int silent),
    CCtsp_x_greedy_tour_lk (CCdatagroup *dat, int ncount, int ecount,
        int *elist, double *x, int *cyc, double *val, int silent,
        CCrandstate *rstate);


/****************************************************************************/
/*                                                                          */
/*                           domboss.c                                      */
/*                                                                          */
/****************************************************************************/

#define CCtsp_DOMINO_WORK        'A'
#define CCtsp_DOMINO_GRAPH       'G'
#define CCtsp_DOMINO_NO          'N'
#define CCtsp_DOMINO_RECEIVE     'R'
#define CCtsp_DOMINO_SEND        'S'
#define CCtsp_DOMINO_WAIT        'W'
#define CCtsp_DOMINO_YES         'Y'
#define CCtsp_DOMINO_EXIT        'X'

#endif  /* __TSP_H */
