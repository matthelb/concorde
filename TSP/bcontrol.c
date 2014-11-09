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
/*                  THE CONTROLLER FOR BRANCHING RUNS                       */
/*                                                                          */
/*                             TSP CODE                                     */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: July 21, 1997                                                     */
/*                                                                          */
/*    EXPORTED FUNCTIONS:                                                   */
/*                                                                          */
/*  int CCtsp_bfs_brancher (char *probloc, int id, double lowerbound,       */
/*      CCtsp_cutselect *sel, CCtsp_cutselect *tsel,                        */
/*      double *upbound, int *bbcount, int usecliques,                      */
/*      CCdatagroup *mydat, int *ptour, CCtsp_lpcuts *pool,                 */
/*      int ncount, int *besttour, unsigned short hostport,                 */
/*      double *branchzeit, int save_proof,                                 */
/*      int tentative_branch_num, int longedge_branching,                   */
/*      double *timebound, int *hit_timebound, int silent,                  */
/*      CCrandstate *rstate)                                                */
/*    PERFORMS a best-first (best-bound) branch and cut search for          */
/*    an optimal tour.                                                      */
/*     -tentative_branch_num specifies the number of trial children         */
/*      created (this should be set to 0 to run standard branching)         */
/*     -timebound can specify an upperbound on the time to spend in the     */
/*      bfs search (it can be NULL; does not work with netgrunts)           */
/*     -hit_timebound will return 1 if timebound is reached (it can be      */
/*      NULL)                                                               */
/*     -if longedge_branching is nonzero, longedge branching will be used   */
/*      instead of the default branching                                    */
/*                                                                          */
/*  int CCtsp_bfs_restart (char *probloc, char *restart_name,               */
/*      CCtsp_cutselect *sel, CCtsp_cutselect *tsel,                        */
/*      double *upbound, int *bbcount, int usecliques,                      */
/*      CCdatagroup *dat, int *ptour, CCtsp_lpcuts *pool,                   */
/*      int ncount, int *besttour,  unsigned short hostport,                */
/*      double *branchzeit, int save_proof,                                 */
/*      int tentative_branch_num, int longedge_branching,                   */
/*      double *timebound, int *hit_timebound, int silent,                  */
/*      CCrandstate *rstate)                                                */
/*    CONTINUES a best-first (best-bound) branch and cut search for         */
/*    an optimal tour from the restart data in restart_name.                */
/*                                                                          */
/*  int CCtsp_grunt (char *hostname, unsigned short hostport,               */
/*      char *poolfname, char *cutbossname, char *probloc,                  */
/*      int silent, CCrandstate *rstate)                                    */
/*    RUNS a grunt served by hostname at port hostport                      */
/*    If probloc is NULL, the probloc will be received from the boss.       */
/*    Only exists if CC_NETREADY is defined                                 */
/*                                                                          */
/*  int CCtsp_easy_dfs_brancher (CCtsp_lp *lp, CCtsp_cutselect *sel,        */
/*      int depth, double *upbound, int *bbcount, int usecliques,           */
/*      int *besttour, int longedge_branching, int simple_branching,        */
/*      int silent, CCrandstate *rstate)                                    */
/*    PERFORMS a depth-first branch and cut search for an optimal tour.     */
/*    NOTES: this will be very inefficient if upbound is not a good bound.  */
/*                                                                          */
/*  int CCtsp_do_interactive_branch (CCtsp_lp *lp, int silent,              */
/*      CCrandstate *rstate)                                                */
/*    SPLITS a problem into two subproblems, prompting the user for         */
/*    information about what split to use.                                  */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"
#include "edgegen.h"
#include "tsp.h"

#define MAX_TENTATIVE_CHILDREN (50)

#define TSP_TENTATIVE_BRANCH_WEIGHT (10.0)
#define TSP_TENATIVE_BRANCH_VAL(v0,v1)                                     \
    (((v0) < (v1) ? (TSP_TENTATIVE_BRANCH_WEIGHT * (v0) + (v1))            \
                  : (TSP_TENTATIVE_BRANCH_WEIGHT * (v1) + (v0)))           \
                    / (TSP_TENTATIVE_BRANCH_WEIGHT + 1.0))


typedef struct tsp_bbnode {
    int id;
    int status;
    int workstatus;
    int numtentative;
    struct tsp_bbnode *prev;
    struct tsp_bbnode *next;
    struct tsp_bbnode *parent;
    struct tsp_bbnode *child0;
    struct tsp_bbnode *child1;
    struct tsp_tnode *tentative_nodes;
    struct tsp_tnode *tparent;
    double lowerbound;
    double cputime;
} tsp_bbnode;

typedef struct tsp_tnode {
    tsp_bbnode *parent;
    tsp_bbnode *child0;
    tsp_bbnode *child1;
} tsp_tnode;

#define BB_NEEDS_CUTTING             (1)
#define BB_NEEDS_TENTATIVE_CUTTING   (2)
#define BB_NEEDS_BRANCHING           (3)
#define BB_NEEDS_TENTATIVE_BRANCHING (4)
#define BB_DONE                      (5)
#define BB_IDLE              (1)
#define BB_WORKING           (2)
#define BB_PRUNED            (3)

#define TASK_WAIT_SECONDS  (60)

typedef struct tsp_bbtask {
#define TASK_EXIT             (1)
#define TASK_WAIT             (2)
#define TASK_TENTATIVE_CUT    (3)
#define TASK_CUT              (4)
#define TASK_BRANCH           (5)
#define TASK_TENTATIVE_BRANCH (6)
    int type;
    int id;
    int new_id;
    int numtentative;
    int tentative[2 * MAX_TENTATIVE_CHILDREN];
    int child0;
    int child1;
} tsp_bbtask;

typedef struct tsp_bbinfo {
    char *hostname;
    unsigned short hostport;
    char *problabel;
    char *probloc;
    int ncount;
    CCdatagroup *dat;
    int *ptour;
    double *upbound;
    CCtsp_lpcuts *pool;
    tsp_bbnode *bblist;
    tsp_bbnode *bbroot;
    CCtsp_cutselect *sel;
    CCtsp_cutselect *tsel;
    int usecliques;
    int *besttour;
    int *bbcount;
    double *branchzeit;
    int max_id;
    int taskcount;
    int changed;
    int finished;
    int gruntcount;
    int cutcount;
    int tcutcount;
    int branchcount;
    int save_proof;
    int tentative_branch_num;
    int longedge_branching;
    int silent;
    CCptrworld bbnode_world;
} tsp_bbinfo;

typedef struct tsp_treport {
    int id0;
    int id1;
    int prune0;
    int prune1;
    double val0;
    double val1;
} tsp_treport;

CC_PTRWORLD_ROUTINES (tsp_bbnode, tsp_bbnode_alloc, tsp_bbnode_bulk_alloc,
        tsp_bbnode_free)


static void
    init_bbnode (tsp_bbnode *bbnode),
    insert_bbnode (tsp_bbnode **firstbbnode, tsp_bbnode *bbnode),
    delete_bbnode (tsp_bbnode **firstbbnode, tsp_bbnode *bbnode),
    bblist_info (tsp_bbnode *bblist, int *cutavail, int *tcutavail,
        int *branchavail, int *active, double *lowerbound),
    collect_active_nodes (tsp_bbnode *b, tsp_bbnode **p_list,
        int *max_id),
    free_tree (tsp_bbnode **bbnode, CCptrworld *bbnode_world);

static int
    bfs_process (tsp_bbinfo *info, CCrandstate *rstate, double *tbound,
         int *hit_limit),
#ifdef CC_NETREADY
    net_process (tsp_bbinfo *info),
    process_connection (CC_SPORT *p, tsp_bbinfo *info),
    boss_send_task (CC_SFILE *f, tsp_bbinfo *info),
    grunt_receive_task (tsp_bbinfo *info, tsp_bbtask *task),
    boss_receive_tour (CC_SFILE *f, tsp_bbinfo *info),
    grunt_send_tour (tsp_bbinfo *info),
    boss_receive_cutnode (CC_SFILE *f, tsp_bbinfo *info),
    boss_receive_tentative_cutnode (CC_SFILE *f, tsp_bbinfo *info),
    grunt_send_cutnode (tsp_bbinfo *info, int id, int new_id, int prune,
        double val, double cputime),
    boss_receive_nobranch (CC_SFILE *f, tsp_bbinfo *info),
    grunt_send_nobranch (tsp_bbinfo *info, int id, double cputime),
    boss_receive_branch (CC_SFILE *f, tsp_bbinfo *info),
    boss_receive_tentative_branch (CC_SFILE *f, tsp_bbinfo *info),
    grunt_send_branch (tsp_bbinfo *info, int id, int child0, int child1,
        double val0, double val1, int prune0, int prune1, double cputime),
    grunt_send_tentative_branch (tsp_bbinfo *info, int id, int num,
        tsp_treport *children, double cputime),
    boss_receive_hello (CC_SFILE *f, tsp_bbinfo *info),
    grunt_send_hello (tsp_bbinfo *info),
    boss_receive_deadnode (CC_SFILE *f, tsp_bbinfo *info),
#endif /* CC_NETREADY */
    get_task (tsp_bbinfo *info, tsp_bbtask *task, int verbose),
    do_task (tsp_bbinfo *info, tsp_bbtask *task, CCrandstate *rstate),
    report_tour (tsp_bbinfo *info),
    report_cut (tsp_bbinfo *info, int id, int new_id, int prune, double val,
        double cputime, int standalone),
    report_tentative_cut (tsp_bbinfo *info, int id, int new_id,
       int prune, double val, double cputime),
    update_tentative_bbnode (tsp_bbinfo *info, tsp_bbnode *b),
    report_nobranch (tsp_bbinfo *info, int id, double cputime),
    report_branch (tsp_bbinfo *info, int id, int child0, int child1,
        double val0, double val1, int prune0, int prune1, double cputime),
    report_tentative_branch (tsp_bbinfo *info, int id, int num,
        tsp_treport *children, double cputime),
    write_restart (char *probname, tsp_bbnode *rootbbnode, double upbound,
        int ncount, int bbcount, double branchzeit),
    write_bbtree (FILE *f, tsp_bbnode *b),
    write_tentative_nodes (FILE *f, int count, tsp_tnode *list),
    read_restart (char *restart_name, char **p_probname,
        tsp_bbnode **p_rootbbnode, double *p_upbound, int *p_ncount,
        int *p_bbcount, double *p_branchzeit, CCptrworld *bbnode_world),
    read_bbtree (FILE *f, tsp_bbnode **p_b, CCptrworld *bbnode_world),
    read_tentative_nodes (FILE *f, int count, tsp_tnode **list,
        tsp_bbnode *parent, CCptrworld *bbnode_world),
    add_children (tsp_bbnode **firstbbnode, tsp_bbnode *parent,
        int id0, int id1, double val0, double val1, int prune0, int prune1,
        CCptrworld *bbnode_world);

static tsp_bbnode
    *find_bbnode (tsp_bbnode *bblist, int id),
    *select_bbnode (tsp_bbnode *firstbbnode, int verbose, int silent);

#ifdef CC_NETREADY
static void
    boss_status (int dir, const char *desc, tsp_bbinfo *info);
#endif /* CC_NETREADY */



static void init_bbnode (tsp_bbnode *bbnode)
{
    bbnode->id           = 0;
    bbnode->lowerbound   = 0.0;
    bbnode->status       = BB_NEEDS_CUTTING;
    bbnode->workstatus   = BB_IDLE;
    bbnode->prev         = (tsp_bbnode *) NULL;
    bbnode->next         = (tsp_bbnode *) NULL;
    bbnode->parent       = (tsp_bbnode *) NULL;
    bbnode->child0       = (tsp_bbnode *) NULL;
    bbnode->child1       = (tsp_bbnode *) NULL;
    bbnode->numtentative = 0;
    bbnode->tentative_nodes  = (tsp_tnode *) NULL;
    bbnode->tparent      = (tsp_tnode *) NULL;
    bbnode->cputime      = 0.0;
}

int CCtsp_bfs_brancher (char *probloc, int id, double lowerbound,
        CCtsp_cutselect *sel, CCtsp_cutselect *tsel, double *upbound,
        int *bbcount, int usecliques, CCdatagroup *dat, int *ptour,
        CCtsp_lpcuts *pool, int ncount, int *besttour,
        unsigned short hostport, double *branchzeit, int save_proof,
        int tentative_branch_num, int longedge_branching,
        double *timebound, int *hit_timebound, int silent,
        CCrandstate *rstate)
{
    int rval = 0;
    tsp_bbnode *rootbbnode  = (tsp_bbnode *) NULL;
    tsp_bbinfo info;
    char *problabel = (char *) NULL;

    CCptrworld_init (&info.bbnode_world);
    
    if (tentative_branch_num > MAX_TENTATIVE_CHILDREN) {
        tentative_branch_num = MAX_TENTATIVE_CHILDREN;
    } else if (tentative_branch_num < 0) {
        tentative_branch_num = 0;
    }

    rootbbnode = tsp_bbnode_alloc (&info.bbnode_world);
    if (!rootbbnode) {
        fprintf (stderr, "Failed to allocate root node\n");
        rval = 1; goto CLEANUP;
    }

    init_bbnode (rootbbnode);
    rootbbnode->id = id;
    rootbbnode->lowerbound = lowerbound;
    rootbbnode->status = BB_NEEDS_BRANCHING;

    *bbcount    = 1;
    *branchzeit = 0.0;

    problabel = CCtsp_problabel (probloc);
    if (problabel == (char *) NULL) {
        fprintf (stderr, "CCtsp_problabel failed\n");
        rval = 1; goto CLEANUP;
    }
    
    info.hostname     = (char *) NULL;
    info.hostport     = hostport;
    info.problabel    = problabel;
    info.probloc      = probloc;
    info.ncount       = ncount;
    info.dat          = dat;
    info.ptour        = ptour;
    info.upbound      = upbound;
    info.pool         = pool;
    info.bblist       = (tsp_bbnode *) NULL;
    info.bbroot       = rootbbnode;
    info.sel          = sel;
    info.tsel         = tsel;
    info.usecliques   = usecliques;
    info.besttour     = besttour;
    info.bbcount      = bbcount;
    info.branchzeit   = branchzeit;
    info.max_id       = -1;
    info.taskcount    = 0;
    info.changed      = 0;
    info.finished     = 0;
    info.gruntcount   = 0;
    info.cutcount     = 0;
    info.tcutcount    = 0;
    info.branchcount  = 0;
    info.save_proof   = save_proof;
    info.tentative_branch_num = tentative_branch_num;
    info.longedge_branching = longedge_branching;
    info.silent = silent;

    if ((unsigned int) hostport == 0) {
        rval = bfs_process (&info, rstate, timebound, hit_timebound);
        if (rval) {
            fprintf (stderr, "bfs_process failed\n"); goto CLEANUP;
        }
        rval = 0;
    } else {
#ifdef CC_NETREADY
        rval = net_process (&info);
        if (rval) {
            fprintf (stderr, "net_process failed\n"); goto CLEANUP;
        }
        rval = 0;
#else
        fprintf (stderr, "Network code is disabled\n");
        rval = 1; goto CLEANUP;
#endif
    }
    
  CLEANUP:
    free_tree (&rootbbnode, &info.bbnode_world);
    CCptrworld_delete (&info.bbnode_world);
    CC_IFFREE (problabel, char);
    return rval;
}

int CCtsp_bfs_restart (char *probloc, char *restart_name, CCtsp_cutselect *sel,
        CCtsp_cutselect *tsel, double *upbound, int *bbcount, int usecliques,
        CCdatagroup *dat, int *ptour, CCtsp_lpcuts *pool, int ncount,
        int *besttour, unsigned short hostport, double *branchzeit,
        int save_proof, int tentative_branch_num, int longedge_branching,
        double *timebound, int *hit_timebound, int silent,
        CCrandstate *rstate)
{
    int rval = 0;
    tsp_bbnode *rootbbnode  = (tsp_bbnode *) NULL;
    char *problabel = (char *) NULL;
    double restart_upbound = 0.0;
    int restart_ncount = 0;
    tsp_bbinfo info;

    CCptrworld_init (&info.bbnode_world);
    
    if (tentative_branch_num > MAX_TENTATIVE_CHILDREN) {
        tentative_branch_num = MAX_TENTATIVE_CHILDREN;
    } else if (tentative_branch_num < 0) {
        tentative_branch_num = 0;
    }

    rval = read_restart (restart_name, &problabel, &rootbbnode, 
                         &restart_upbound, &restart_ncount, bbcount,
                         branchzeit, &info.bbnode_world);
    if (rval) {
        fprintf (stderr, "read_restart failed\n");
        goto CLEANUP;
    }
    if (ncount != restart_ncount) {
        fprintf (stderr, "wrong ncount in restart file\n");
        rval = 1; goto CLEANUP;
    }
    if (restart_upbound < *upbound) *upbound = restart_upbound;

    info.hostname    = (char *) NULL;
    info.hostport    = hostport;
    info.probloc     = probloc;
    info.problabel   = problabel;
    info.ncount      = ncount;
    info.dat         = dat;
    info.ptour       = ptour;
    info.upbound     = upbound;
    info.pool        = pool;
    info.bblist      = (tsp_bbnode *) NULL;
    info.bbroot      = rootbbnode;
    info.sel         = sel;
    info.tsel        = tsel;
    info.usecliques  = usecliques;
    info.besttour    = besttour;
    info.bbcount     = bbcount;
    info.branchzeit  = branchzeit;
    info.max_id      = -1;
    info.taskcount   = 0;
    info.changed     = 0;
    info.finished    = 0;
    info.gruntcount  = 0;
    info.cutcount    = 0;
    info.tcutcount   = 0;
    info.branchcount = 0;
    info.save_proof  = save_proof;
    info.tentative_branch_num = tentative_branch_num;
    info.longedge_branching = longedge_branching;
    info.silent = silent;
    
    if ((unsigned int) hostport == 0) {
        rval = bfs_process (&info, rstate, timebound, hit_timebound);
        if (rval) {
            fprintf (stderr, "bfs_process failed\n"); goto CLEANUP;
        }
        rval = 0;
    } else {
#ifdef CC_NETREADY
        rval = net_process (&info);
        if (rval) {
            fprintf (stderr, "net_process failed\n"); goto CLEANUP;
        }
        rval = 0;
#else /* CC_NETREADY */
        fprintf (stderr, "Network branching not enabled\n");
        rval = 1;
        goto CLEANUP;
#endif /* CC_NETREADY */
    }
    
  CLEANUP:
    CC_IFFREE (problabel, char);
    free_tree (&rootbbnode, &info.bbnode_world);
    CCptrworld_delete (&info.bbnode_world);
    return rval;
}

static int bfs_process (tsp_bbinfo *info, CCrandstate *rstate, double *tbound,
        int *hit_tbound)
{
    int rval = 0;
    double szeit = CCutil_zeit ();
    tsp_bbtask bbtask;

    if (hit_tbound) *hit_tbound = 0;

    collect_active_nodes (info->bbroot, &info->bblist, &info->max_id);
    /* convert from circular to linear list */
    info->bblist->prev->next = (tsp_bbnode *) NULL;
    info->bblist->prev = (tsp_bbnode *) NULL;

    rval = write_restart (info->problabel, info->bbroot, *info->upbound,
                          info->ncount, *info->bbcount, *info->branchzeit);
    if (rval) {
        fprintf (stderr, "write_restart failed\n");
        return rval;
    }

    info->finished = 0;
    info->changed  = 0;

    while (info->finished == 0) {
        if (tbound) {
            if (CCutil_zeit () - szeit > *tbound) {
                fprintf (stderr, "Hit time limit in bfs branching\n");
                if (hit_tbound) *hit_tbound = 1;
                 return 0;
            }
        }

        rval = get_task (info, &bbtask, 1);
        if (rval) {
            fprintf (stderr, "get_task failed\n");
            return rval;
        }

        rval = do_task (info, &bbtask, rstate);
        if (rval) {
            fprintf (stderr, "do_task failed\n");
            return rval;
        }
        
        if (info->changed) {
            rval = write_restart (info->problabel, info->bbroot,
                    *info->upbound, info->ncount, *info->bbcount,
                    *info->branchzeit);
            if (rval) {
                fprintf (stderr, "write_restart failed\n");
                return rval;
            }
            info->changed = 0;
        }
    }
    return 0;
}

static void bblist_info (tsp_bbnode *bblist, int *cutavail, int *tcutavail,
        int *branchavail, int *active, double *lowerbound)
{
    tsp_bbnode *b;
    int ccnt = 0;
    int bcnt = 0;
    int tcnt = 0;
    int acnt = 0;
    int i;
    double l = CCtsp_LP_MAXDOUBLE;

    for (b = bblist; b; b = b->next) {
        if (b->lowerbound < l) l = b->lowerbound;
        if (b->status != BB_DONE) acnt++;
        if (b->status == BB_NEEDS_CUTTING) {
            ccnt++;
        } else if (b->status == BB_NEEDS_BRANCHING) {
            if (b->tentative_nodes == (tsp_tnode *) NULL) {
                bcnt++;
            } else {
                for (i = 0; i < b->numtentative; i++) {
                    if (b->tentative_nodes[i].child0->status ==
                        BB_NEEDS_TENTATIVE_CUTTING) tcnt++;
                    if (b->tentative_nodes[i].child1->status ==
                        BB_NEEDS_TENTATIVE_CUTTING) tcnt++;
                }
            }
        } else {
            fprintf (stderr, "Hmm, bbnode status %d on active list\n",
                     b->status);
        }
    }
    if (cutavail != (int *) NULL) *cutavail = ccnt;
    if (tcutavail != (int *) NULL) *tcutavail = tcnt;
    if (branchavail != (int *) NULL) *branchavail = bcnt;
    if (active != (int *) NULL) *active = acnt;
    if (lowerbound != (double *) NULL) *lowerbound = l;
}
                         
#ifdef CC_NETREADY

static int net_process (tsp_bbinfo *info)
{
    int rval = 0;
    CC_SPORT *p = (CC_SPORT *) NULL;

    printf ("\nBEGINNING NET PROCESSING\n\n");
    fflush (stdout);

    collect_active_nodes (info->bbroot, &info->bblist, &info->max_id);
    /* convert from circular to linear list */
    info->bblist->prev->next = (tsp_bbnode *) NULL;
    info->bblist->prev = (tsp_bbnode *) NULL;

    rval = write_restart (info->problabel, info->bbroot, *info->upbound,
                          info->ncount, *info->bbcount, *info->branchzeit);
    if (rval) {
        fprintf (stderr, "write_restart failed\n");
        goto CLEANUP;
    }

    p = CCutil_snet_listen (info->hostport);
    if (p == (CC_SPORT *) NULL) {
        fprintf (stderr, "CCutil_snet_listen failed\n");
        rval = 1; goto CLEANUP;
    }
    
    info->finished = 0;
    info->changed  = 0;

    while (info->finished == 0) {
        rval = process_connection (p, info);
        if (rval) {
            fprintf (stderr, "process_connection failed\n");
        }

        if (info->changed) {
            rval = write_restart (info->problabel, info->bbroot,
                    *info->upbound, info->ncount, *info->bbcount,
                    *info->branchzeit);
            if (rval) {
                fprintf (stderr, "write_restart failed\n"); goto CLEANUP;
            }
            info->changed = 0;
        }
    }
    rval = 0;
    
  CLEANUP:
    if (p != (CC_SPORT *) NULL) {
        CCutil_snet_unlisten (p);
    }
    return rval;
}

static int process_connection (CC_SPORT *p, tsp_bbinfo *info)
{
    CC_SFILE *f = (CC_SFILE *) NULL;
    int rval;
    char request;

    f = CCutil_snet_receive (p);
    if (f == (CC_SFILE *) NULL) {
        fprintf (stderr, "CCutil_snet_receive failed\n");
        rval = 1; goto CLEANUP;
    }

    rval = CCutil_sread_char (f, &request);
    if (rval) goto CLEANUP;

    switch (request) {
      case CCtsp_BBREQ_TASK:
        rval = boss_send_task (f, info);
        if (rval) {
            fprintf (stderr, "boss_send_task failed\n");
            goto CLEANUP;
        }
        break;
      case CCtsp_BBREQ_TOUR:
        rval = boss_receive_tour (f, info);
        if (rval) {
            fprintf (stderr, "boss_receive_tour failed\n");
            goto CLEANUP;
        }
        boss_status (1, "TOUR:", info);
        break;
      case CCtsp_BBREQ_CUTDONE:
        rval = boss_receive_cutnode (f, info);
        if (rval) {
            fprintf (stderr, "boss_receive_cutnode failed\n");
            goto CLEANUP;
        }
        info->cutcount--;
        boss_status (1, "CUT:", info);
        break;
      case CCtsp_BBREQ_TENTATIVE_CUTDONE:
        rval = boss_receive_tentative_cutnode (f, info);
        if (rval) {
            fprintf (stderr, "boss_receive_tentative_cutnode failed\n");
            goto CLEANUP;
        }
        info->tcutcount--;
        boss_status (1, "TCUT:", info);
        break;
      case CCtsp_BBREQ_NOBRANCH:
        rval = boss_receive_nobranch (f, info);
        if (rval) {
            fprintf (stderr, "boss_receive_nobranch failed\n");
            goto CLEANUP;
        }
        info->branchcount--;
        boss_status (1, "NOBRAN:", info);
        break;
      case CCtsp_BBREQ_BRANCHDONE:
        rval = boss_receive_branch (f, info);
        if (rval) {
            fprintf (stderr, "boss_receive_branch failed\n");
            goto CLEANUP;
        }
        info->branchcount--;
        boss_status (1, "BRANCH:", info);
        break;
      case CCtsp_BBREQ_TENTATIVE_BRANCHDONE:
        rval = boss_receive_tentative_branch (f, info);
        if (rval) {
            fprintf (stderr, "boss_receive_tentative_branch failed\n");
            goto CLEANUP;
        }
        info->branchcount--;
        boss_status (1, "TBRANCH:", info);
        break;
      case CCtsp_BBREQ_HELLO:
        rval = boss_receive_hello (f, info);
        if (rval) {
            fprintf (stderr, "boss_receive_hello failed\n");
            goto CLEANUP;
        }
        info->gruntcount++;
        boss_status (1, "HELLO:", info);
        break;
      case CCtsp_BBREQ_DEADNODE:
        rval = boss_receive_deadnode (f, info);
        if (rval) {
            fprintf (stderr, "boss_receive_deadnode failed\n");
            goto CLEANUP;
        }
        boss_status (1, "DEAD:", info);
        break;
      case CCtsp_BBREQ_EXIT:
        info->finished = 1;
        boss_status (1, "EXIT:", info);
        break;
      default:
        fprintf (stderr, "Unknown host request %c\n", request);
        rval = 1; goto CLEANUP;
    }
    rval = CCutil_sclose (f);
    if (rval) {
        fprintf (stderr, "CCutil_sclose failed\n");
        goto CLEANUP;
    }

    f = (CC_SFILE *) NULL;
    rval = 0;
  CLEANUP:
    if (f != (CC_SFILE *) NULL) {
        CCutil_sclose (f);
    }
    return rval;
}

static void boss_status (int dir, const char *desc, tsp_bbinfo *info)
{
    int cutavail = 0;
    int tcutavail = 0;
    int active = 0;
    int branchavail = 0;
    double lb = CCtsp_LP_MAXDOUBLE;

    bblist_info (info->bblist, &cutavail, &tcutavail, &branchavail,
                 &active, &lb);
    
    printf ("%-1s %-8s  cut %2d/%3d  tcut %2d/%3d  bra %2d/%3d  wait %d  cpu %.2f\n",
            dir == 1 ? "R" : "S", desc, info->cutcount, cutavail,
            info->tcutcount, tcutavail, info->branchcount, branchavail,
            info->gruntcount - info->cutcount - info->tcutcount -
            info->branchcount, *info->branchzeit);
    if (dir == 1) {
        printf ("  Lower Bound: %f   Remaining Nodes: %d\n", lb, active);
    }
    fflush (stdout);
}

int CCtsp_grunt (char *hostname, unsigned short hostport, char *poolfname,
        char *cutbossname, char *probloc, int silent,
        CCrandstate *rstate)
{
    char probbuf[CCtsp_PROB_FILE_NAME_LEN];
    CCtsp_cutselect cutselect, tcutselect;
    tsp_bbinfo info;
    CCdatagroup dat;
    double upbound;
    tsp_bbtask bbtask;
    int rval;

    CCptrworld_init (&info.bbnode_world);
    CCtsp_init_cutselect (&cutselect);
    CCtsp_init_tentative_cutselect (&tcutselect);
    CCutil_init_datagroup (&dat);
    
    info.hostname = hostname;
    info.hostport = hostport;
    info.probloc = probbuf;
    info.problabel = (char *) NULL;
    info.ncount = 0;
    info.dat = &dat;
    info.ptour = (int *) NULL;
    info.upbound = &upbound;
    info.pool = (CCtsp_lpcuts *) NULL;
    info.bblist = (tsp_bbnode *) NULL;
    info.bbroot = (tsp_bbnode *) NULL;
    info.sel  = &cutselect;
    info.tsel = &tcutselect;
    info.usecliques = 0;
    info.besttour = (int *) NULL;
    info.bbcount = (int *) NULL;
    info.branchzeit = (double *) NULL;
    info.max_id = 0;
    info.taskcount = 0;
    info.changed = 0;
    info.finished = 0;
    info.longedge_branching = 0;
    info.silent = silent;

    rval = grunt_send_hello (&info);
    if (rval) {
        fprintf (stderr, "send_hello failed\n");
        goto CLEANUP;
    }

    if (probloc) info.probloc = probloc;
    
    info.besttour = CC_SAFE_MALLOC (info.ncount, int);
    if (info.besttour == (int *) NULL) {
        fprintf (stderr, "Out of memory in CCtsp_grunt\n");
        goto CLEANUP;
    }

    rval = CCtsp_init_cutpool (&info.ncount, poolfname, &info.pool);
    if (rval) {
        fprintf (stderr, "CCtsp_init_cutpool failed\n");
        goto CLEANUP;
    }

    if (cutbossname) {
        cutselect.remotepool = 1;
        cutselect.remotehost = cutbossname;
        cutselect.remoteport = CCtsp_CUT_PORT;
        tcutselect.remotepool = 1;
        tcutselect.remotehost = cutbossname;
        tcutselect.remoteport = CCtsp_CUT_PORT;
    }
    
    while (info.finished == 0) {
        rval = grunt_receive_task (&info, &bbtask);
        if (rval) {
            fprintf (stderr, "grunt_receive_task failed\n");
            goto CLEANUP;
        }

        rval = do_task (&info, &bbtask, rstate);
        if (rval) {
            fprintf (stderr, "do_task failed\n");
            goto CLEANUP;
        }
    }

    rval = 0;

  CLEANUP:
    CC_IFFREE (info.ptour, int);
    CC_IFFREE (info.besttour, int);
    CC_IFFREE (info.problabel, char);
    CCutil_freedatagroup (&dat);
    if (info.pool) {
        CCtsp_free_cutpool (&info.pool);
    }
    CCptrworld_delete (&info.bbnode_world);
    
    return rval;
}

static int boss_send_task (CC_SFILE *f, tsp_bbinfo *info)
{
    tsp_bbtask bbtask;
    int i;
    int rval = 0;
    
    rval = get_task (info, &bbtask, 0);
    if (rval) {
        fprintf (stderr, "get_task failed\n");
        return rval;
    }
    switch (bbtask.type) {
      case TASK_EXIT:
        rval = CCutil_swrite_char (f, CCtsp_BBTASK_EXIT);
        if (rval) return rval;
        info->gruntcount--;
        boss_status (2, "EXIT:", info);
        if (info->gruntcount == 0) info->finished = 1;
        break;
      case TASK_WAIT:
        rval = CCutil_swrite_char (f, CCtsp_BBTASK_WAIT);
        if (rval) return rval;
/*
        boss_status (2, "WAIT:", info);
*/
        break;
      case TASK_TENTATIVE_CUT:
        rval = CCutil_swrite_char (f, CCtsp_BBTASK_TENTATIVE_CUT);
        if (rval) return rval;
        rval = CCutil_swrite_int (f, bbtask.id);
        if (rval) return rval;
        rval = CCutil_swrite_int (f, bbtask.new_id);
        if (rval) return rval;
        rval = CCutil_swrite_double (f, *info->upbound);
        if (rval) return rval;
        info->tcutcount++;
        boss_status (2, "TCUT:", info);
        break;
      case TASK_CUT:
        rval = CCutil_swrite_char (f, CCtsp_BBTASK_CUT);
        if (rval) return rval;
        rval = CCutil_swrite_int (f, bbtask.id);
        if (rval) return rval;
        rval = CCutil_swrite_int (f, bbtask.new_id);
        if (rval) return rval;
        rval = CCutil_swrite_double (f, *info->upbound);
        if (rval) return rval;
        info->cutcount++;
        boss_status (2, "CUT:", info);
        break;
      case TASK_BRANCH:
        rval = CCutil_swrite_char (f, CCtsp_BBTASK_BRANCH);
        if (rval) return rval;
        rval = CCutil_swrite_int (f, bbtask.id);
        if (rval) return rval;
        rval = CCutil_swrite_int (f, bbtask.child0);
        if (rval) return rval;
        rval = CCutil_swrite_int (f, bbtask.child1);
        if (rval) return rval;
        rval = CCutil_swrite_double (f, *info->upbound);
        if (rval) return rval;
        info->branchcount++;
        boss_status (2, "BRANCH:", info);
        break;
      case TASK_TENTATIVE_BRANCH:
        rval = CCutil_swrite_char (f, CCtsp_BBTASK_TENTATIVE_BRANCH);
        if (rval) return rval;
        rval = CCutil_swrite_int (f, bbtask.id);
        if (rval) return rval;
        rval = CCutil_swrite_int (f, bbtask.numtentative);
        if (rval) return rval;
        for (i = 0; i < (2*bbtask.numtentative); i++) {
            rval = CCutil_swrite_int (f, bbtask.tentative[i]);
            if (rval) return rval;
        }
        rval = CCutil_swrite_double (f, *info->upbound);
        if (rval) return rval;
        info->branchcount++;
        boss_status (2, "TBRANCH:", info);
        break;
      default:
        fprintf (stderr, "Unknown bbtask type %d\n", bbtask.type);
        return 1;
    }
    return 0;
}

static int grunt_receive_task (tsp_bbinfo *info, tsp_bbtask *task)
{
    char tasktype;
    CC_SFILE *f = (CC_SFILE *) NULL;
    int rval = 0;
    int i;

    f = CCutil_snet_open (info->hostname, info->hostport);
    if (f == (CC_SFILE *) NULL) {
        fprintf (stderr, "CCutil_snet_open failed\n");
        rval = 1; goto CLEANUP;
    }

    rval = CCutil_swrite_char (f, CCtsp_BBREQ_TASK);
    if (rval) goto CLEANUP;

    rval = CCutil_sread_char (f, &tasktype);
    if (rval) goto CLEANUP;

    switch (tasktype) {
      case CCtsp_BBTASK_EXIT:
        task->type = TASK_EXIT;
        break;
      case CCtsp_BBTASK_WAIT:
        task->type = TASK_WAIT;
        break;
      case CCtsp_BBTASK_TENTATIVE_CUT:
        task->type = TASK_TENTATIVE_CUT;
        rval = CCutil_sread_int (f, &task->id);
        if (rval) goto CLEANUP;
        rval = CCutil_sread_int (f, &task->new_id);
        if (rval) goto CLEANUP;
        rval = CCutil_sread_double (f, info->upbound);
        if (rval) goto CLEANUP;
        break;
      case CCtsp_BBTASK_CUT:
        task->type = TASK_CUT;
        rval = CCutil_sread_int (f, &task->id);
        if (rval) goto CLEANUP;
        rval = CCutil_sread_int (f, &task->new_id);
        if (rval) goto CLEANUP;
        rval = CCutil_sread_double (f, info->upbound);
        if (rval) goto CLEANUP;
        break;
      case CCtsp_BBTASK_BRANCH:
        task->type = TASK_BRANCH;
        rval = CCutil_sread_int (f, &task->id);
        if (rval) goto CLEANUP;
        rval = CCutil_sread_int (f, &task->child0);
        if (rval) goto CLEANUP;
        rval = CCutil_sread_int (f, &task->child1);
        if (rval) goto CLEANUP;
        rval = CCutil_sread_double (f, info->upbound);
        if (rval) goto CLEANUP;
        break;
      case CCtsp_BBTASK_TENTATIVE_BRANCH:
        task->type = TASK_TENTATIVE_BRANCH;
        rval = CCutil_sread_int (f, &task->id);
        if (rval) goto CLEANUP;
        rval = CCutil_sread_int (f, &task->numtentative);
        if (rval) goto CLEANUP;
        for (i = 0; i < (2*task->numtentative); i++) {
            rval = CCutil_sread_int (f, &task->tentative[i]);
            if (rval) return rval;
        }
        rval = CCutil_sread_double (f, info->upbound);
        if (rval) goto CLEANUP;
        break;
      default:
        fprintf (stderr, "Unknown bbtask code %c\n", tasktype);
        rval = 1; goto CLEANUP;
    }

    rval = CCutil_sclose (f);
    if (rval) {
        fprintf (stderr, "CCutil_sclose failed\n");
        goto CLEANUP;
    }
    f = (CC_SFILE *) NULL;
    rval = 0;

  CLEANUP:
    if (f != (CC_SFILE *) NULL) {
        CCutil_sclose (f);
    }
    return rval;
}

static int boss_receive_tour (CC_SFILE *f, tsp_bbinfo *info)
{
    double val;
    int rval;
    int ncount = info->ncount;
    int *tour = (int *) NULL;
    int i;

    tour = CC_SAFE_MALLOC (ncount, int);
    if (tour == (int *) NULL) {
        fprintf (stderr, "Out of memory in receive_tour\n");
        rval = 1; goto CLEANUP;
    }
    
    rval = CCutil_sread_double (f, &val);
    if (rval) goto CLEANUP;

    for (i=0; i<ncount; i++) {
        rval = CCutil_sread_int (f, &tour[i]);
        if (rval) goto CLEANUP;
    }
    
    if (val < *info->upbound) {
        *info->upbound = val;
        for (i=0; i<ncount; i++) {
            info->besttour[i] = tour[i];
        }
        rval = report_tour (info);
        if (rval) {
            fprintf (stderr, "report_tour failed\n");
            goto CLEANUP;
        }
    }

    rval = 0;
  CLEANUP:
    CC_IFFREE (tour, int);
    return rval;
}

static int grunt_send_tour (tsp_bbinfo *info)
{
    CC_SFILE *f = (CC_SFILE *) NULL;
    int rval;
    int i;

    f = CCutil_snet_open (info->hostname, info->hostport);
    if (f == (CC_SFILE *) NULL) {
        fprintf (stderr, "CCutil_snet_open failed\n");
        rval = 1; goto CLEANUP;
    }

    rval = CCutil_swrite_char (f, CCtsp_BBREQ_TOUR);
    if (rval) goto CLEANUP;

    rval = CCutil_swrite_double (f, *info->upbound);
    if (rval) goto CLEANUP;

    for (i=0; i<info->ncount; i++) {
        rval = CCutil_swrite_int (f, info->besttour[i]);
        if (rval) goto CLEANUP;
    }
    
    rval = CCutil_sclose (f);
    if (rval) {
        fprintf (stderr, "CCutil_sclose failed\n");
        goto CLEANUP;
    }
    f = (CC_SFILE *) NULL;
    rval = 0;

  CLEANUP:
    if (f != (CC_SFILE *) NULL) {
        CCutil_sclose (f);
    }
    return rval;
}

static int boss_receive_cutnode (CC_SFILE *f, tsp_bbinfo *info)
{
    int rval;
    int id;
    int new_id;
    int prune;
    double val;
    double cputime;

    rval = CCutil_sread_int (f, &id);
    if (rval) return rval;
    rval = CCutil_sread_int (f, &new_id);
    if (rval) return rval;
    rval = CCutil_sread_int (f, &prune);
    if (rval) return rval;
    rval = CCutil_sread_double (f, &val);
    if (rval) return rval;
    rval = CCutil_sread_double (f, &cputime);
    if (rval) return rval;

    rval = report_cut (info, id, new_id, prune, val, cputime, 0);
    if (rval) {
        fprintf (stderr, "report_cut failed\n");
        return rval;
    }
    return 0;
}

static int boss_receive_tentative_cutnode (CC_SFILE *f, tsp_bbinfo *info)
{
    int rval = 0;
    int id;
    int new_id;
    int prune;
    double val;
    double cputime;

    rval = CCutil_sread_int (f, &id);
    if (rval) return rval;
    rval = CCutil_sread_int (f, &new_id);
    if (rval) return rval;
    rval = CCutil_sread_int (f, &prune);
    if (rval) return rval;
    rval = CCutil_sread_double (f, &val);
    if (rval) return rval;
    rval = CCutil_sread_double (f, &cputime);
    if (rval) return rval;

    rval = report_tentative_cut (info, id, new_id, prune, val, cputime);
    if (rval) {
        fprintf (stderr, "report_tentative_cut failed\n");
        return rval;
    }
    return 0;
}


static int grunt_send_cutnode (tsp_bbinfo *info, int id, int new_id, int prune,
        double val, double cputime)
{
    CC_SFILE *f = (CC_SFILE *) NULL;
    int rval;

    f = CCutil_snet_open (info->hostname, info->hostport);
    if (f == (CC_SFILE *) NULL) {
        fprintf (stderr, "CCutil_snet_open failed\n");
        rval = 1; goto CLEANUP;
    }

    rval = CCutil_swrite_char (f, CCtsp_BBREQ_CUTDONE);
    if (rval) goto CLEANUP;
    
    rval = CCutil_swrite_int (f, id);
    if (rval) goto CLEANUP;
    rval = CCutil_swrite_int (f, new_id);
    if (rval) goto CLEANUP;
    rval = CCutil_swrite_int (f, prune);
    if (rval) goto CLEANUP;
    rval = CCutil_swrite_double (f, val);
    if (rval) goto CLEANUP;
    rval = CCutil_swrite_double (f, cputime);
    if (rval) goto CLEANUP;

    rval = CCutil_sclose (f);
    if (rval) {
        fprintf (stderr, "CCutil_sclose failed\n");
        goto CLEANUP;
    }
    f = (CC_SFILE *) NULL;
    rval = 0;

  CLEANUP:
    if (f != (CC_SFILE *) NULL) {
        CCutil_sclose (f);
    }
    return rval;
}

static int grunt_send_tentative_cutnode (tsp_bbinfo *info, int id, int new_id,
        int prune, double val, double cputime)
{
    CC_SFILE *f = (CC_SFILE *) NULL;
    int rval;

    f = CCutil_snet_open (info->hostname, info->hostport);
    if (f == (CC_SFILE *) NULL) {
        fprintf (stderr, "CCutil_snet_open failed\n");
        rval = 1; goto CLEANUP;
    }

    rval = CCutil_swrite_char (f, CCtsp_BBREQ_TENTATIVE_CUTDONE);
    if (rval) goto CLEANUP;
    
    rval = CCutil_swrite_int (f, id);
    if (rval) goto CLEANUP;
    rval = CCutil_swrite_int (f, new_id);
    if (rval) goto CLEANUP;
    rval = CCutil_swrite_int (f, prune);
    if (rval) goto CLEANUP;
    rval = CCutil_swrite_double (f, val);
    if (rval) goto CLEANUP;
    rval = CCutil_swrite_double (f, cputime);
    if (rval) goto CLEANUP;

    rval = CCutil_sclose (f);
    if (rval) {
        fprintf (stderr, "CCutil_sclose failed\n");
        goto CLEANUP;
    }
    f = (CC_SFILE *) NULL;
    rval = 0;

  CLEANUP:
    if (f != (CC_SFILE *) NULL) {
        CCutil_sclose (f);
    }
    return rval;
}

static int boss_receive_nobranch (CC_SFILE *f, tsp_bbinfo *info)
{
    int id;
    int rval;
    double cputime;

    rval = CCutil_sread_int (f, &id);
    if (rval) return rval;
    rval = CCutil_sread_double (f, &cputime);
    if (rval) return rval;

    rval = report_nobranch (info, id, cputime);
    if (rval) {
        fprintf (stderr, "report_nobranch failed\n");
        return rval;
    }
    return 0;
}

static int grunt_send_nobranch (tsp_bbinfo *info, int id, double cputime)
{
    CC_SFILE *f = (CC_SFILE *) NULL;
    int rval;

    f = CCutil_snet_open (info->hostname, info->hostport);
    if (f == (CC_SFILE *) NULL) {
        fprintf (stderr, "CCutil_snet_open failed\n");
        rval = 1; goto CLEANUP;
    }

    rval = CCutil_swrite_char (f, CCtsp_BBREQ_NOBRANCH);
    if (rval) goto CLEANUP;
    
    rval = CCutil_swrite_int (f, id);
    if (rval) goto CLEANUP;
    rval = CCutil_swrite_double (f, cputime);
    if (rval) goto CLEANUP;

    rval = CCutil_sclose (f);
    if (rval) {
        fprintf (stderr, "CCutil_sclose failed\n");
        goto CLEANUP;
    }
    f = (CC_SFILE *) NULL;
    rval = 0;

  CLEANUP:
    if (f != (CC_SFILE *) NULL) {
        CCutil_sclose (f);
    }
    return rval;
}

static int boss_receive_branch (CC_SFILE *f, tsp_bbinfo *info)
{
    int id;
    int child0;
    int child1;
    double val0;
    double val1;
    int prune0;
    int prune1;
    int rval;
    double cputime;

    rval = CCutil_sread_int (f, &id);
    if (rval) return rval;
    rval = CCutil_sread_int (f, &child0);
    if (rval) return rval;
    rval = CCutil_sread_double (f, &val0);
    if (rval) return rval;
    rval = CCutil_sread_int (f, &prune0);
    if (rval) return rval;
    rval = CCutil_sread_int (f, &child1);
    if (rval) return rval;
    rval = CCutil_sread_double (f, &val1);
    if (rval) return rval;
    rval = CCutil_sread_int (f, &prune1);
    if (rval) return rval;
    rval = CCutil_sread_double (f, &cputime);
    if (rval) return rval;

    rval = report_branch (info, id, child0, child1, val0, val1, prune0,
                          prune1, cputime);
    if (rval) {
        fprintf (stderr, "report_branch failed\n");
        return rval;
    }
    return 0;
}

static int boss_receive_tentative_branch (CC_SFILE *f, tsp_bbinfo *info)
{
    int rval = 0;
    int i, id, num;
    tsp_treport *children = (tsp_treport *) NULL;
    double cputime;

    rval = CCutil_sread_int (f, &id);
    if (rval) goto CLEANUP;
    rval = CCutil_sread_int (f, &num);
    if (rval) goto CLEANUP;
    if (num <= 0) {
        fprintf (stderr, "received %d tentative children\n", num);
        rval = 1; goto CLEANUP;
    }
    children = CC_SAFE_MALLOC (num, tsp_treport);
    if (children == (tsp_treport *) NULL) {
        fprintf (stderr, "out of memory in boss_receive_tentative_branch\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < num; i++) {
        rval = CCutil_sread_int (f, &children[i].id0);
        if (rval) goto CLEANUP;
        rval = CCutil_sread_double (f, &children[i].val0);
        if (rval) goto CLEANUP;
        rval = CCutil_sread_int (f, &children[i].prune0);
        if (rval) goto CLEANUP;
        rval = CCutil_sread_int (f, &children[i].id1);
        if (rval) goto CLEANUP;
        rval = CCutil_sread_double (f, &children[i].val1);
        if (rval) goto CLEANUP;
        rval = CCutil_sread_int (f, &children[i].prune1);
        if (rval) goto CLEANUP;
    }
    rval = CCutil_sread_double (f, &cputime);
    if (rval) goto CLEANUP;

    rval = report_tentative_branch (info, id, num, children, cputime);
    if (rval) {
        fprintf (stderr, "report_tentative_branch failed\n");
        return rval;
    }

CLEANUP:

    CC_IFFREE (children, tsp_treport);
    return rval;
}

static int grunt_send_branch (tsp_bbinfo *info, int id, int child0, int child1,
        double val0, double val1, int prune0, int prune1, double cputime)
{
    CC_SFILE *f = (CC_SFILE *) NULL;
    int rval;

    f = CCutil_snet_open (info->hostname, info->hostport);
    if (f == (CC_SFILE *) NULL) {
        fprintf (stderr, "CCutil_snet_open failed\n");
        rval = 1; goto CLEANUP;
    }

    rval = CCutil_swrite_char (f, CCtsp_BBREQ_BRANCHDONE);
    if (rval) goto CLEANUP;
    
    rval = CCutil_swrite_int (f, id);
    if (rval) goto CLEANUP;
    rval = CCutil_swrite_int (f, child0);
    if (rval) goto CLEANUP;
    rval = CCutil_swrite_double (f, val0);
    if (rval) goto CLEANUP;
    rval = CCutil_swrite_int (f, prune0);
    if (rval) goto CLEANUP;
    rval = CCutil_swrite_int (f, child1);
    if (rval) goto CLEANUP;
    rval = CCutil_swrite_double (f, val1);
    if (rval) goto CLEANUP;
    rval = CCutil_swrite_int (f, prune1);
    if (rval) goto CLEANUP;
    rval = CCutil_swrite_double (f, cputime);
    if (rval) goto CLEANUP;

    rval = CCutil_sclose (f);
    if (rval) {
        fprintf (stderr, "CCutil_sclose failed\n");
        goto CLEANUP;
    }
    f = (CC_SFILE *) NULL;
    rval = 0;

  CLEANUP:
    if (f != (CC_SFILE *) NULL) {
        CCutil_sclose (f);
    }
    return rval;
}

static int grunt_send_tentative_branch (tsp_bbinfo *info, int id, int num,
        tsp_treport *children, double cputime)
{
    CC_SFILE *f = (CC_SFILE *) NULL;
    int rval = 0;
    int i;

    f = CCutil_snet_open (info->hostname, info->hostport);
    if (f == (CC_SFILE *) NULL) {
        fprintf (stderr, "CCutil_snet_open failed\n");
        rval = 1; goto CLEANUP;
    }

    rval = CCutil_swrite_char (f, CCtsp_BBREQ_TENTATIVE_BRANCHDONE);
    if (rval) goto CLEANUP;
    
    rval = CCutil_swrite_int (f, id);
    if (rval) goto CLEANUP;
    rval = CCutil_swrite_int (f, num);
    if (rval) goto CLEANUP;

    for (i = 0; i < num; i++) {
        rval = CCutil_swrite_int (f, children[i].id0);
        if (rval) goto CLEANUP;
        rval = CCutil_swrite_double (f, children[i].val0);
        if (rval) goto CLEANUP;
        rval = CCutil_swrite_int (f, children[i].prune0);
        if (rval) goto CLEANUP;
        rval = CCutil_swrite_int (f, children[i].id1);
        if (rval) goto CLEANUP;
        rval = CCutil_swrite_double (f, children[i].val1);
        if (rval) goto CLEANUP;
        rval = CCutil_swrite_int (f, children[i].prune1);
        if (rval) goto CLEANUP;
    }

    rval = CCutil_swrite_double (f, cputime);
    if (rval) goto CLEANUP;

    rval = CCutil_sclose (f);
    if (rval) {
        fprintf (stderr, "CCutil_sclose failed\n");
        goto CLEANUP;
    }
    f = (CC_SFILE *) NULL;
    rval = 0;

  CLEANUP:
    if (f != (CC_SFILE *) NULL) {
        CCutil_sclose (f);
    }
    return rval;
}

static int grunt_send_hello (tsp_bbinfo *info)
{
    CC_SFILE *f = (CC_SFILE *) NULL;
    CCtsp_cutselect *s;
    int i;
    int rval;

    f = CCutil_snet_open (info->hostname, info->hostport);
    if (f == (CC_SFILE *) NULL) {
        fprintf (stderr, "CCutil_snet_open failed\n");
        rval = 1; goto CLEANUP;
    }

    rval = CCutil_swrite_char (f, CCtsp_BBREQ_HELLO);
    if (rval) goto CLEANUP;
    rval = CCutil_sread_string (f, info->probloc, CCtsp_PROB_FILE_NAME_LEN);
    if (rval) goto CLEANUP;

    CC_IFFREE (info->problabel, char);
    info->problabel = CCtsp_problabel (info->probloc);
    if (info->problabel == (char *) NULL) {
        fprintf (stderr, "CCtsp_problabel failed\n");
        rval = 1; goto CLEANUP;
    }
    
    for (i = 0; i <= 1; i++) {
        if (i == 0) s = info->sel;
        else        s = info->tsel;

        rval = CCutil_sread_int (f, &s->cutpool);
        if (rval) goto CLEANUP;
        rval = CCutil_sread_int (f, &s->connect);
        if (rval) goto CLEANUP;
        rval = CCutil_sread_int (f, &s->segments);
        if (rval) goto CLEANUP;
        rval = CCutil_sread_int (f, &s->exactsubtour);
        if (rval) goto CLEANUP;
        rval = CCutil_sread_int (f, &s->tighten_lp);
        if (rval) goto CLEANUP;
        rval = CCutil_sread_int (f, &s->teething_lp);
        if (rval) goto CLEANUP;
        rval = CCutil_sread_int (f, &s->tighten_pool);
        if (rval) goto CLEANUP;
        rval = CCutil_sread_int (f, &s->decker_lp);
        if (rval) goto CLEANUP;
        rval = CCutil_sread_int (f, &s->decker_pool);
        if (rval) goto CLEANUP;
        rval = CCutil_sread_int (f, &s->maxchunksize);
        if (rval) goto CLEANUP;
        rval = CCutil_sread_int (f, &s->teething_pool);
        if (rval) goto CLEANUP;
        rval = CCutil_sread_int (f, &s->exactblossom);
        if (rval) goto CLEANUP;
        rval = CCutil_sread_int (f, &s->consecutiveones);
        if (rval) goto CLEANUP;
        rval = CCutil_sread_int (f, &s->necklace);
        if (rval) goto CLEANUP;
        rval = CCutil_sread_int (f, &s->usetighten);
        if (rval) goto CLEANUP;
        rval = CCutil_sread_int (f, &s->extra_connect);
        if (rval) goto CLEANUP;
    }

    rval = CCutil_sread_int (f, &info->usecliques);
    if (rval) goto CLEANUP;
    rval = CCutil_sread_int (f, &info->longedge_branching);
    if (rval) goto CLEANUP;
    rval = CCutil_readmaster (f, &info->ncount, info->dat, &info->ptour);
    if (rval) goto CLEANUP;

    rval = CCutil_sclose (f);
    if (rval) {
        fprintf (stderr, "CCutil_sclose failed\n");
        goto CLEANUP;
    }
    f = (CC_SFILE *) NULL;
    rval = 0;

  CLEANUP:
    if (f != (CC_SFILE *) NULL) {
        CCutil_sclose (f);
    }
    return rval;
}

static int boss_receive_hello (CC_SFILE *f, tsp_bbinfo *info)
{
    int rval;
    int i;
    CCtsp_cutselect *s;
    
    rval = CCutil_swrite_string (f, info->probloc);
    if (rval) return rval;

    for (i = 0; i <= 1; i++) {
        if (i == 0) s = info->sel;
        else        s = info->tsel;

        rval = CCutil_swrite_int (f, s->cutpool);
        if (rval) return rval;
        rval = CCutil_swrite_int (f, s->connect);
        if (rval) return rval;
        rval = CCutil_swrite_int (f, s->segments);
        if (rval) return rval;
        rval = CCutil_swrite_int (f, s->exactsubtour);
        if (rval) return rval;
        rval = CCutil_swrite_int (f, s->tighten_lp);
        if (rval) return rval;
        rval = CCutil_swrite_int (f, s->teething_lp);
        if (rval) return rval;
        rval = CCutil_swrite_int (f, s->tighten_pool);
        if (rval) return rval;
        rval = CCutil_swrite_int (f, s->decker_lp);
        if (rval) return rval;
        rval = CCutil_swrite_int (f, s->decker_pool);
        if (rval) return rval;
        rval = CCutil_swrite_int (f, s->maxchunksize);
        if (rval) return rval;
        rval = CCutil_swrite_int (f, s->teething_pool);
        if (rval) return rval;
        rval = CCutil_swrite_int (f, s->exactblossom);
        if (rval) return rval;
        rval = CCutil_swrite_int (f, s->consecutiveones);
        if (rval) return rval;
        rval = CCutil_swrite_int (f, s->necklace);
        if (rval) return rval;
        rval = CCutil_swrite_int (f, s->usetighten);
        if (rval) return rval;
        rval = CCutil_swrite_int (f, s->extra_connect);
        if (rval) return rval;
    }

    rval = CCutil_swrite_int (f, info->usecliques);
    if (rval) return rval;
    rval = CCutil_swrite_int (f, info->longedge_branching);
    if (rval) return rval;
    rval = CCutil_writemaster (f, info->ncount, info->dat, info->ptour);
    if (rval) return rval;
    
    return 0;
}

static int boss_receive_deadnode (CC_SFILE *f, tsp_bbinfo *info)
{
    int id;
    int rval;
    tsp_bbnode *bbnode;

    rval = CCutil_sread_int (f, &id);
    if (rval) return rval;

    if (id == -1) {
        info->gruntcount--;
        return 0;
    }
    
    bbnode = find_bbnode (info->bblist, id);
    if (bbnode == (tsp_bbnode *) NULL) {
        printf ("BBnode %d no longer active\n", id);
        return 0;
    } else if (bbnode->workstatus != BB_WORKING) {
        printf ("BBnode %d is not working\n", id);
        return 0;
    }

    if (bbnode->status == BB_NEEDS_CUTTING) {
        info->cutcount--;
    } else if (bbnode->status == BB_NEEDS_BRANCHING) {
        info->branchcount--;
    } else if (bbnode->status == BB_NEEDS_TENTATIVE_CUTTING) {
        info->tcutcount--;
    }
    bbnode->workstatus = BB_IDLE;
    info->gruntcount--;

    return 0;
}

#endif /* CC_NETREADY */

static int get_task (tsp_bbinfo *info, tsp_bbtask *task, int verbose)
{
    tsp_bbnode *b = (tsp_bbnode *) NULL;
    tsp_bbnode *bblist = info->bblist;
    int i;
    
    if (bblist == (tsp_bbnode *) NULL) {
        task->type = TASK_EXIT;
        printf ("\n"); fflush (stdout);
        return 0;
    } else {
        b = select_bbnode (bblist, verbose, info->silent);
        if (b == (tsp_bbnode *) NULL) {
            task->type = TASK_WAIT;
            return 0;
        } else {
            b->workstatus = BB_WORKING;
            task->id = b->id;
            if (b->status == BB_NEEDS_CUTTING) {
                task->type = TASK_CUT;
                task->new_id = ++(info->max_id);
            } else if (b->status == BB_NEEDS_TENTATIVE_CUTTING) {
                task->type = TASK_TENTATIVE_CUT;
                task->new_id = ++(info->max_id);
            } else if (b->status == BB_NEEDS_BRANCHING) {
                if (info->tentative_branch_num > 0) {
                    task->type = TASK_TENTATIVE_BRANCH;
                    task->numtentative = info->tentative_branch_num;
                    for (i = 0; i < (2*task->numtentative); i++) {
                        task->tentative[i] = ++(info->max_id);
                    }
                } else {
                    task->type = TASK_BRANCH;
                    task->child0 = ++(info->max_id);
                    task->child1 = ++(info->max_id);
                }
            } else {
                fprintf (stderr, "Bogus bbnode status %d\n", b->status);
                return 1;
            }
            return 0;
        }
    }
}
    
static int do_task (tsp_bbinfo *info, tsp_bbtask *task, CCrandstate *rstate)
{
    int i, j, foundtour, ngot;
    int rval;
    double val, val0, val1;
    int id0, id1, prune, prune0, prune1;
    CCtsp_branchobj *b = (CCtsp_branchobj *) NULL;
    tsp_treport *trp = (tsp_treport *) NULL;
    double szeit;
    double cputime;
    
    switch (task->type) {
      case TASK_EXIT:
        printf ("Task %d: Exit\n", info->taskcount++);
        fflush (stdout);
        info->finished = 1;
        break;
      case TASK_WAIT:
        printf ("Task %d: Wait\n", info->taskcount++);
        fflush (stdout);
#ifdef HAVE_SLEEP
        sleep (TASK_WAIT_SECONDS);
#else
        {
            int j;
            for (j=0; j<10000000; j++)
                ;
        }
#endif /* HAVE_SLEEP */
        break;
      case TASK_TENTATIVE_CUT:
        printf ("Task %d: Tentative Cutting on node %d\n",
                 info->taskcount++, task->id);
        fflush (stdout);
        szeit = CCutil_zeit();
        rval = CCtsp_bb_cutting (info->probloc, task->id, task->new_id, 
                info->ncount, info->dat, info->ptour, info->upbound,
                info->pool, info->tsel, &val, &prune, &foundtour,
                info->besttour, 0, info->silent, rstate);
        if (rval) {
            fprintf (stderr, "CCtsp_bb_cutting failed\n"); goto CLEANUP;
        }
        cputime = CCutil_zeit() - szeit;
        if (foundtour) {
            rval = report_tour (info);
            if (rval) {
                fprintf (stderr, "report_tour failed\n"); goto CLEANUP;
            }
        }
        rval = report_tentative_cut (info, task->id, task->new_id, prune, val,
                                     cputime);
        if (rval) {
            fprintf (stderr, "report_tentative_cut failed\n"); goto CLEANUP;
        }
        break;
      case TASK_CUT:
        printf ("Task %d: Cutting on node %d\n", info->taskcount++, task->id);
        fflush (stdout);
        szeit = CCutil_zeit();
        rval = CCtsp_bb_cutting (info->probloc, task->id, task->new_id, 
                info->ncount, info->dat, info->ptour, info->upbound,
                info->pool, info->sel, &val, &prune, &foundtour,
                info->besttour, 1 /* 1 for normal, -1 for fast FAST HACK */,
                info->silent, rstate);
        if (rval) {
            fprintf (stderr, "CCtsp_bb_cutting failed\n"); goto CLEANUP;
        }
        cputime = CCutil_zeit() - szeit;
        if (foundtour) {
            rval = report_tour (info);
            if (rval) {
                fprintf (stderr, "report_tour failed\n"); goto CLEANUP;
            }
        }
        rval = report_cut (info, task->id, task->new_id, prune, val, cputime,
                           1);
        if (rval) {
            fprintf (stderr, "report_cut failed\n"); goto CLEANUP;
        }
        break;
      case TASK_BRANCH:
        printf ("Task %d: Branching on node %d\n", info->taskcount++,
                task->id);
        fflush (stdout);
        szeit = CCutil_zeit();
        rval = CCtsp_bb_find_branch (info->probloc, task->id, info->ncount,
                info->dat, info->ptour, info->upbound, info->pool, 1, &ngot,
                &b, info->usecliques, info->longedge_branching, &prune,
                &foundtour, info->besttour, info->silent, rstate);
        if (rval) {
            fprintf (stderr, "CCtsp_bb_find_branch failed\n"); goto CLEANUP;
        }
        if (prune) {
            cputime = CCutil_zeit() - szeit;
            rval = report_nobranch (info, task->id, cputime);
            if (rval) {
                fprintf (stderr, "report_nobranch failed\n"); goto CLEANUP;
            }
        } else if (foundtour) {
            cputime = CCutil_zeit() - szeit;
            rval = report_tour (info);
            if (rval) {
                fprintf (stderr, "report_tour failed\n"); goto CLEANUP;
            }
            rval = report_nobranch (info, task->id, cputime);
            if (rval) {
                fprintf (stderr, "report_nobranch failed\n"); goto CLEANUP;
            }
        } else {
            if (!info->silent) {
                printf ("Found Branch - split problem into children\n");
                fflush (stdout);
            }
            id0 = task->child0;
            id1 = task->child1;
            rval = CCtsp_bb_splitprob (info->probloc, task->id, info->ncount,
                    info->dat, info->ptour, *info->upbound, info->pool,
                    &b[0], id0, id1, &val0, &val1, &prune0, &prune1,
                    info->silent, rstate);
            CCtsp_free_branchobj (&b[0]);
            CC_IFFREE (b, CCtsp_branchobj);
            if (rval) {
                fprintf (stderr, "CCtsp_bb_splitprob failed\n"); goto CLEANUP;
            }
            cputime = CCutil_zeit() - szeit;
            rval = report_branch (info, task->id, id0, id1, val0, val1,
                                  prune0, prune1, cputime);
            if (rval) {
                fprintf (stderr, "report_branch failed\n"); goto CLEANUP;
            }
        }
        break;
      case TASK_TENTATIVE_BRANCH:
        printf ("Task %d: Tentative branching on node %d\n", info->taskcount++,
                task->id);
        fflush (stdout);
        szeit = CCutil_zeit();
        rval = CCtsp_bb_find_branch (info->probloc, task->id, info->ncount,
                info->dat, info->ptour, info->upbound, info->pool, 
                task->numtentative, &ngot, &b, 
                info->usecliques, info->longedge_branching, &prune, &foundtour,
                info->besttour, info->silent, rstate);
        if (rval) {
            fprintf (stderr, "CCtsp_bb_find_branch failed\n"); goto CLEANUP;
        }
        if (prune) {
            cputime = CCutil_zeit() - szeit;
            rval = report_nobranch (info, task->id, cputime);
            if (rval) {
                fprintf (stderr, "report_nobranch failed\n"); goto CLEANUP;
            }
        } else if (foundtour) {
            cputime = CCutil_zeit() - szeit;
            rval = report_tour (info);
            if (rval) {
                fprintf (stderr, "report_tour failed\n"); goto CLEANUP;
            }
            rval = report_nobranch (info, task->id, cputime);
            if (rval) {
                fprintf (stderr, "report_nobranch failed\n"); goto CLEANUP;
            }
        } else {
            if (!info->silent) {
                printf ("Found Tentative Branch for %d - create tchildren\n",
                             task->id);
                fflush (stdout);

                for (i = 0; i < ngot; i++) {
                    printf ("    Tbranchobj %d: ", i); 
                    if (b[i].ends[0] != -1) {
                        printf ("Edge (%d,%d)", b[i].ends[0], b[i].ends[1]);
                    } else {
                        printf ("Clique ");
                        for (j = 0; j < b[i].clique->segcount; j++) {
                            printf ("%d->%d ", b[i].clique->nodes[j].lo,
                                               b[i].clique->nodes[j].hi);
                        }
                    }
                    printf ("\n");
                    fflush (stdout);
                }
            }

            trp = CC_SAFE_MALLOC (ngot, tsp_treport);
            if (trp == (tsp_treport *) NULL) {
                fprintf (stderr, "out of memory in do_task\n");
                rval = 1; goto CLEANUP;
            }
            for (i = 0; i < ngot; i++) {
                trp[i].id0 = task->tentative[2*i];
                trp[i].id1 = task->tentative[2*i+1];
                rval = CCtsp_bb_splitprob (info->probloc, task->id,
                    info->ncount,
                    info->dat, info->ptour, *info->upbound, info->pool,
                    &b[i], trp[i].id0, trp[i].id1, &trp[i].val0, &trp[i].val1,
                    &trp[i].prune0, &trp[i].prune1, info->silent, rstate);
                if (rval) {
                    fprintf (stderr, "CCtsp_bb_splitprob failed\n");
                    
                    for (j = 0; j < i; j++) {
                        CCtsp_free_branchobj (&b[j]);
                    }
                    CC_IFFREE (b, CCtsp_branchobj);
                    goto CLEANUP;
                }
            }
            cputime = CCutil_zeit() - szeit;
            rval = report_tentative_branch (info, task->id, ngot, trp,
                                            cputime);
            if (rval) {
                fprintf (stderr, "report_tentative_branch failed\n");
                goto CLEANUP;
            }
        }
        break;
      default:
        fprintf (stderr, "BOGUS TASK: %d\n", task->type);
        rval = 1; goto CLEANUP;
    }
    
    rval = 0;
    
 CLEANUP:
    if (b != (CCtsp_branchobj *) NULL) {
        CCtsp_free_branchobj (&b[0]);
    }
    CC_IFFREE (b, CCtsp_branchobj);
    CC_IFFREE (trp, tsp_treport);
    return rval;
}

static int report_tour (tsp_bbinfo *info)
{
    int rval;

    printf ("TOUR FOUND - upperbound is %.2f\n", *info->upbound);
    fflush (stdout);

#ifdef CC_NETREADY
    if (info->hostname != (char *) NULL) {
        rval = grunt_send_tour (info);
        if (rval) {
            fprintf (stderr, "grunt_send_tour failed\n");
            return rval;
        }
        return 0;
    }
#endif /* CC_NETREADY */
    
    rval = CCtsp_dumptour (info->ncount, info->dat, info->ptour,
               info->problabel, info->besttour, (char *) NULL, 0, info->silent);
    if (rval) {
        fprintf (stderr, "CCtsp_dumptour failed\n");
        return rval;
    }
    info->changed = 1;

    return 0;
}

static int report_cut (tsp_bbinfo *info, int id, int new_id, int prune,
        double val, double cputime, int standalone)
{
    tsp_bbnode *bbnode;
    int rval;

    if (standalone && info->hostname == (char *) NULL) {
        char buf[1024];
        printf ("Writing Pool: %d cuts\n", info->pool->cutcount);
        fflush (stdout);
        sprintf (buf, "%s.pul", info->problabel);
        rval = CCtsp_write_cutpool (info->ncount, buf, info->pool);
        if (rval) {
            fprintf (stderr, "CCtsp_write_cutpool failed\n"); return rval;
        }
    }

#ifdef CC_NETREADY
    if (info->hostname != (char *) NULL) {
        rval = grunt_send_cutnode (info, id, new_id, prune, val, cputime);
        if (rval) {
            fprintf (stderr, "grunt_send_cutnode failed\n");
            return rval;
        }
        return 0;
    }
#endif /* CC_NETREADY */
    
    bbnode = find_bbnode (info->bblist, id);
    if (bbnode == (tsp_bbnode *) NULL) {
        printf ("BBnode %d no longer active\n", id);
        return 0;
    }

    if (bbnode->status != BB_NEEDS_CUTTING) {
        printf ("BBnode %d does not need cutting\n", id);
        return 0;
    }

    bbnode->id = new_id;
    bbnode->lowerbound = val;
    bbnode->cputime   += cputime;
    *info->branchzeit += cputime;

    if (prune) {
        printf ("BBnode %d (now %d) can be pruned: upperbound %.2f (%.2f seconds)\n",
                id, new_id, *info->upbound, cputime);
        fflush (stdout);
        delete_bbnode (&info->bblist, bbnode);
        info->changed = 1;
        rval = CCtsp_prob_file_delete (info->probloc, id);
        if (rval) return rval;
        if (info->save_proof == 0) {
            rval = CCtsp_prob_file_delete (info->probloc, new_id);
            if (rval) return rval;
        }
    } else {
        printf ("BBnode %d (now %d) done cutting: lowerbound %.2f (%.2f seconds)\n",
                id, new_id, val, cputime);
        fflush (stdout);
        bbnode->status     = BB_NEEDS_BRANCHING;
        bbnode->workstatus = BB_IDLE;
        info->changed = 1;
        rval = CCtsp_prob_file_delete (info->probloc, id);
        if (rval) return rval;
    }
    return 0;
}

static int report_tentative_cut (tsp_bbinfo *info, int id, int new_id,
       int prune, double val, double cputime)
{
    tsp_bbnode *bbnode;
    int rval;

#ifdef CC_NETREADY
    if (info->hostname != (char *) NULL) {
        rval = grunt_send_tentative_cutnode (info, id, new_id, prune, val, 
                                             cputime);
        if (rval) {
            fprintf (stderr, "grunt_send_tentative_cutnode failed\n");
            return rval;
        }
        return 0;
    }
#endif /* CC_NETREADY */
    
    bbnode = find_bbnode (info->bblist, id);
    if (bbnode == (tsp_bbnode *) NULL) {
        printf ("BBnode %d no longer active\n", id);
        return 0;
    }
    if (bbnode->tparent == (tsp_tnode *) NULL) {
        printf ("BBnode %d is not a tentative bbnode\n", id);
        return 0;
    }
    if (bbnode->status != BB_NEEDS_TENTATIVE_CUTTING) {
        printf ("BBnode %d does not need tentative cutting\n", id);
        return 0;
    }
    bbnode->id = new_id;
    bbnode->lowerbound = val;
    bbnode->cputime   += cputime;
    *info->branchzeit += cputime;

    printf ("Tnode %d (now %d) done tentative cutting: lowerbound %.2f (%.2f seconds)\n",
            id, new_id, val, cputime);
    fflush (stdout);
    bbnode->status = BB_DONE;
    if (prune) {
        bbnode->workstatus = BB_PRUNED;
    } else {
        bbnode->workstatus = BB_IDLE;
    }
    update_tentative_bbnode (info, bbnode->tparent->parent);
    info->changed = 1;
    rval = CCtsp_prob_file_delete (info->probloc, id);
    if (rval) return rval;

    return 0;
}

static int update_tentative_bbnode (tsp_bbinfo *info, tsp_bbnode *b)
{
    tsp_tnode *t;
    int rval = 0;
    int i;

    for (i = 0; i < b->numtentative; i++) {
        if (b->tentative_nodes[i].child0->status != BB_DONE) break;
        if (b->tentative_nodes[i].child1->status != BB_DONE) break;
    }

    if (i == b->numtentative) {
        tsp_tnode *tbest = (tsp_tnode *) NULL;
        double bestval = -CCtsp_LP_MAXDOUBLE;
        double val;

        printf ("Tentative Branching on BBnode %d\n", b->id);
        fflush (stdout);
        for (i = 0; i < b->numtentative; i++) {
            t = &(b->tentative_nodes[i]);
            printf ("    Tbranch %d:  %9.2f %9.2f\n", i,
                              t->child0->lowerbound, t->child1->lowerbound);
            fflush (stdout);
            val = TSP_TENATIVE_BRANCH_VAL (t->child0->lowerbound, 
                                           t->child1->lowerbound);
            if (val > bestval) {
                bestval = val;
                tbest = t;
            }
        }
        if (tbest == (tsp_tnode *) NULL) {
            fprintf (stderr, "error in update_tentative_bbnode\n");
            return 1;
        }

        b->child0 = tbest->child0;
        b->child1 = tbest->child1;
        if (b->child0->workstatus == BB_PRUNED) {
            b->child0->status = BB_DONE;
            if (info->save_proof == 0) {
                rval = CCtsp_prob_file_delete (info->probloc, b->child0->id);
                if (rval) return rval;
            }
        } else {
            b->child0->status = BB_NEEDS_CUTTING;
            insert_bbnode (&info->bblist, b->child0);
        }
        if (b->child1->workstatus == BB_PRUNED) {
            b->child1->status = BB_DONE;
            if (info->save_proof == 0) {
                rval = CCtsp_prob_file_delete (info->probloc, b->child1->id);
                if (rval) return rval;
            }
        } else {
            b->child1->status = BB_NEEDS_CUTTING;
            insert_bbnode (&info->bblist, b->child1);
        }
        b->child0->workstatus = BB_IDLE;
        b->child1->workstatus = BB_IDLE;
        b->status = BB_DONE;
        b->workstatus = BB_IDLE;

        for (i = 0; i < b->numtentative; i++) {
            t = &(b->tentative_nodes[i]);
            if (t != tbest) {
                rval = CCtsp_prob_file_delete (info->probloc, t->child0->id);
                if (rval) return rval;
                tsp_bbnode_free (&info->bbnode_world, t->child0);
                rval = CCtsp_prob_file_delete (info->probloc, t->child1->id);
                if (rval) return rval;
                tsp_bbnode_free (&info->bbnode_world, t->child1);
            }
        }
        CC_IFFREE (b->tentative_nodes, tsp_tnode);
        b->numtentative = 0;

        delete_bbnode (&info->bblist, b);
        if (info->save_proof == 0 || b->id != 0) { /* don't delete the root */
            rval = CCtsp_prob_file_delete (info->probloc, b->id);
            if (rval) return rval;
        }
        printf ("BBnode %d split into %d (%.2f%s) %d (%.2f%s)\n", b->id,
            b->child0->id, b->child0->lowerbound,
            (b->child0->status == BB_DONE) ? "X" : "",
            b->child1->id, b->child1->lowerbound,
            (b->child1->status == BB_DONE) ? "X" : "");
        fflush (stdout);
        *(info->bbcount) += 2;
    }
    return 0;
}

static int report_nobranch (tsp_bbinfo *info, int id, double cputime)
{
    tsp_bbnode *bbnode;
    int rval = 0;

#ifdef CC_NETREADY
    if (info->hostname != (char *) NULL) {
        rval = grunt_send_nobranch (info, id, cputime);
        if (rval) {
            fprintf (stderr, "grunt_send_nobranch failed\n");
            return rval;
        }
        return 0;
    }
#endif /* CC_NETREADY */

    bbnode = find_bbnode (info->bblist, id);
    if (bbnode == (tsp_bbnode *) NULL) {
        printf ("BBnode %d no longer active\n", id);
        return 0;
    } else if (bbnode->status != BB_NEEDS_BRANCHING) {
        printf ("BBnode %d does not need branching\n", id);
        return 0;
    }

    bbnode->cputime   += cputime;
    *info->branchzeit += cputime;

    printf ("BBnode %d is pruned - no branching (%.2f seconds)\n", id,
            cputime);
    fflush (stdout);
    delete_bbnode (&info->bblist, bbnode);
    info->changed = 1;
    if (info->save_proof == 0 || id != 0) { /* don't delete the root */
        rval = CCtsp_prob_file_delete (info->probloc, id);
        if (rval) return rval;
    }
    return 0;
}

static int report_branch (tsp_bbinfo *info, int id, int child0, int child1,
        double val0, double val1, int prune0, int prune1, double cputime)
{
    tsp_bbnode *bbnode;
    int rval;

#ifdef CC_NETREADY
    if (info->hostname != (char *) NULL) {
        rval = grunt_send_branch (info, id, child0, child1, val0, val1,
                prune0, prune1, cputime);
        if (rval) {
            fprintf (stderr, "grunt_send_branch failed\n");
            return rval;
        }
        return 0;
    }
#endif /* CC_NETREADY */

    bbnode = find_bbnode (info->bblist, id);
    if (bbnode == (tsp_bbnode *) NULL) {
        printf ("BBnode %d no longer active\n", id);
        return 0;
    } else if (bbnode->status != BB_NEEDS_BRANCHING) {
        printf ("BBnode %d does not need branching\n", id);
        return 0;
    }

    bbnode->cputime   += cputime;
    *info->branchzeit += cputime;
    
    printf ("BBnode %d split into %d (%.2f%s) %d (%.2f%s) (%.2f seconds)\n",
            id, child0, val0, prune0 ? "X" : "", child1, val1,
            prune1 ? "X" : "", cputime);
    fflush (stdout);
    
    rval = add_children (&info->bblist, bbnode, child0, child1, val0, val1,
                         prune0, prune1, &info->bbnode_world);
    if (rval) {
        fprintf (stderr, "add_children failed\n"); return rval;
    }
    info->changed = 1;
    *(info->bbcount) += 2;
    delete_bbnode (&info->bblist, bbnode);
    if (info->save_proof == 0 || id != 0) { /* don't delete the root */
        rval = CCtsp_prob_file_delete (info->probloc, id);
        if (rval) return rval;
    }
    if (info->save_proof == 0) {
        if (prune0) {
            rval = CCtsp_prob_file_delete (info->probloc, child0);
            if (rval) return rval;
        }
        if (prune1) {
            rval = CCtsp_prob_file_delete (info->probloc, child1);
            if (rval) return rval;
        }
    }
    
    return 0;
}

static int report_tentative_branch (tsp_bbinfo *info, int id, int num,
        tsp_treport *children, double cputime)
{
    tsp_bbnode *bbnode, *child;
    int i;
    int rval = 0;

#ifdef CC_NETREADY
    if (info->hostname != (char *) NULL) {
        rval = grunt_send_tentative_branch (info, id, num, children, cputime);
        if (rval) {
            fprintf (stderr, "grunt_send_tentative_branch failed\n");
            return rval;
        }
        return 0;
    }
#endif /* CC_NETREADY */

    bbnode = find_bbnode (info->bblist, id);
    if (bbnode == (tsp_bbnode *) NULL) {
        printf ("BBnode %d no longer active\n", id);
        goto CLEANUP;
    } else if (bbnode->status != BB_NEEDS_BRANCHING) {
        printf ("BBnode %d does not need branching\n", id);
        goto CLEANUP;
    }

    bbnode->cputime   += cputime;
    *info->branchzeit += cputime;
    bbnode->workstatus = BB_IDLE;
    
    printf ("BBnode %d tentative split with %d trials (%.2f seconds)\n",
            id, num, cputime);
    fflush (stdout);

    bbnode->numtentative = num;
    bbnode->tentative_nodes = CC_SAFE_MALLOC (num, tsp_tnode);
    if (bbnode->tentative_nodes == (tsp_tnode *) NULL) {
        fprintf (stderr, "out of memory in report_tentative_branch\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < num; i++) {
        printf ("  %3d %9.2f   %3d %9.2f\n",
                 children[i].id0, children[i].val0,
                 children[i].id1, children[i].val1);
        fflush (stdout);
        bbnode->tentative_nodes[i].parent = bbnode;
        child = tsp_bbnode_alloc (&info->bbnode_world);
        if (!child) {
            fprintf (stderr, "Failed to allocate child 0\n");
            rval = 1; goto CLEANUP;
        }
    
        init_bbnode (child);
        child->id = children[i].id0;
        child->lowerbound = children[i].val0;
        child->parent = bbnode;
        child->tparent = &(bbnode->tentative_nodes[i]);  
        bbnode->tentative_nodes[i].child0 = child;

        if (children[i].val0 == CCtsp_LP_MAXDOUBLE) {
            child->status = BB_DONE;
            child->workstatus = BB_PRUNED;
        } else if (children[i].prune0) {
            child->status = BB_DONE;
            child->workstatus = BB_PRUNED;
        } else {
            child->status = BB_NEEDS_TENTATIVE_CUTTING;
            child->workstatus = BB_IDLE;
        }

        child = tsp_bbnode_alloc (&info->bbnode_world);
        if (!child) {
            fprintf (stderr, "Failed to allocate child 0\n");
            rval = 1; goto CLEANUP;
        }
    
        init_bbnode (child);
        child->id = children[i].id1;
        child->lowerbound = children[i].val1;
        child->parent = bbnode;
        child->tparent = &(bbnode->tentative_nodes[i]);  
        bbnode->tentative_nodes[i].child1 = child;

        if (children[i].val1 == CCtsp_LP_MAXDOUBLE) {
            child->status = BB_DONE;
            child->workstatus = BB_PRUNED;
        } else if (children[i].prune1) {
            child->status = BB_DONE;
            child->workstatus = BB_PRUNED;
        } else {
            child->status = BB_NEEDS_TENTATIVE_CUTTING;
            child->workstatus = BB_IDLE;
        }
    }
    info->changed = 1;
    rval =  update_tentative_bbnode (info, bbnode);
    if (rval) {
        fprintf (stderr, "update_tentative_bbnode failed\n"); goto CLEANUP;
    }

CLEANUP:

    return rval;
}

static tsp_bbnode *find_bbnode (tsp_bbnode *bblist, int id)
{
    tsp_tnode *t;
    int i;

    while (bblist != (tsp_bbnode *) NULL) {
        if (bblist->id == id) return bblist;
        for (i = 0; i < bblist->numtentative; i++) {
            t = &(bblist->tentative_nodes[i]);
            if (t->child0->id == id) return t->child0;
            if (t->child1->id == id) return t->child1;
        }
        bblist = bblist->next;
    }
    return (tsp_bbnode *) NULL;
}

static void collect_active_nodes (tsp_bbnode *b, tsp_bbnode **p_list,
        int *max_id)
{
    int i;

    if (b->id > *max_id) *max_id = b->id;
    for (i = 0; i < b->numtentative; i++) {
        if (b->tentative_nodes[i].child0->id > *max_id) {
            *max_id = b->tentative_nodes[i].child0->id;
        }
        if (b->tentative_nodes[i].child1->id > *max_id) {
            *max_id = b->tentative_nodes[i].child1->id;
        }
    }
    
    if (b->status != BB_DONE) {
        if ((*p_list) == (tsp_bbnode *) NULL) {
            b->next = b;
            b->prev = b;
            (*p_list) = b;
        } else {
            b->next = (*p_list);
            b->prev = (*p_list)->prev;
            b->prev->next = b;
            b->next->prev = b;
        }
    }
    if (b->child0) collect_active_nodes (b->child0, p_list, max_id);
    if (b->child1) collect_active_nodes (b->child1, p_list, max_id);
}

static int write_restart (char *problabel, tsp_bbnode *rootbbnode,
        double upbound, int ncount, int bbcount, double branchzeit)
{
    FILE *f = (FILE *) NULL;
    char *restart_name = (char *) NULL;
    char *new_name = (char *) NULL;
    char *back_name = (char *) NULL;
    int rval = 0;
    size_t len;

    len = strlen(problabel) + 6;
    restart_name = CC_SAFE_MALLOC (len, char);
    new_name = CC_SAFE_MALLOC (len, char);
    back_name = CC_SAFE_MALLOC (len, char);
    if (restart_name == (char *) NULL ||
        new_name == (char *) NULL ||
        back_name == (char *) NULL) {
        fprintf (stderr, "Out of memory in write_restart\n");
        rval = 1; goto CLEANUP;
    }

    strcpy (restart_name, problabel);
    strcat (restart_name, ".res");
    strcpy (new_name, "N");
    strcat (new_name, restart_name);
    strcpy (back_name, "O");
    strcat (back_name, restart_name);
    
    f = fopen (new_name, "w");
    if (f == (FILE*) NULL) {
        perror (new_name);
        fprintf (stderr, "Unable to open %s for output in write_restart\n",
                 new_name);
        rval = 1; goto CLEANUP;
    }

    rval = fprintf (f, "%s %d %.0f %d %.2f\n", problabel, ncount, upbound,
            bbcount, branchzeit);
    if (rval <= 0) {
        perror (new_name);
        fprintf (stderr, "fprintf to %s failed\n", new_name);
        rval = 1; goto CLEANUP;
    }
    rval = write_bbtree (f, rootbbnode);
    if (rval) {
        fprintf (stderr, "write_bbtree failed\n"); goto CLEANUP;
    }
    
    rval = fclose (f);
    if (rval) {
        perror (new_name);
        fprintf (stderr, "fclose %s failed\n", new_name);
        rval = 1; goto CLEANUP;
    }
    f = (FILE *) NULL;

    rename (restart_name, back_name);
    rval = rename (new_name, restart_name);
    if (rval) {
        perror (restart_name);
        fprintf (stderr, "rename %s to %s failed\n", new_name, restart_name);
        rval = 1; goto CLEANUP;
    }
    
    rval = 0;

  CLEANUP:
    if (f != (FILE *) NULL) {
        fclose (f);
    }
    CC_IFFREE (new_name, char);
    CC_IFFREE (restart_name, char);
    CC_IFFREE (back_name, char);
    return rval;
}

static int write_bbtree (FILE *f, tsp_bbnode *b)
{
    int rval = 0;

    rval = fprintf (f, "%1d %4d %4d %4d %4d %13.2f %9.2f\n", b->status, b->id,
                    b->child0 ? b->child0->id : -1,
                    b->child1 ? b->child1->id : -1,
                    b->numtentative,
                    b->lowerbound, b->cputime);
    if (rval <= 0) {
        perror ("restart_file");
        fprintf (stderr, "fprintf failed writing restart file\n");
        return 1;
    }
    if (b->tentative_nodes) {
        rval = write_tentative_nodes (f, b->numtentative, b->tentative_nodes);
        if (rval) return rval;
    }
    if (b->child0) {
        rval = write_bbtree (f, b->child0);
        if (rval) return rval;
    }
    if (b->child1) {
        rval = write_bbtree (f, b->child1);
        if (rval) return rval;
    }
    return 0;
}

static int write_tentative_nodes (FILE *f, int count, tsp_tnode *list)
{
    int rval = 0;
    int i;
    tsp_tnode *s;

    for (i = 0; i < count; i++) {
        s = &list[i];
        if (s->child0 == (tsp_bbnode *) NULL ||
            s->child1 == (tsp_bbnode *) NULL) {
            fprintf (stderr, "tnode has NULL bbnodes\n");
            return 1;
        }
        rval = fprintf (f, "  %1d %4d %13.2f %9.2f    %1d %4d %13.2f %9.2f\n",
               s->child0->status, s->child0->id, s->child0->lowerbound,
               s->child0->cputime,
               s->child1->status, s->child1->id, s->child1->lowerbound,
               s->child1->cputime);
        if (rval <= 0) {
            perror ("restart_file");
            fprintf (stderr, "fprintf failed writing restart file\n");
            return 1;
        }
    } 
    return 0;
}

static int read_restart (char *restart_name, char **p_problabel,
        tsp_bbnode **p_rootbbnode, double *p_upbound, int *p_ncount,
        int *p_bbcount, double *p_branchzeit, CCptrworld *bbnode_world)
{
    FILE *f = (FILE *) NULL;
    char *problabel = (char *) NULL;
    int rval = 0;

    f = fopen (restart_name, "r");
    if (f == (FILE*) NULL) {
        perror (restart_name);
        fprintf (stderr, "Unable to open %s for input in read_restart\n",
                 restart_name);
        rval = 1; goto CLEANUP;
    }

    problabel = CC_SAFE_MALLOC (CCtsp_PROB_FILE_NAME_LEN, char);
    if (problabel == (char *) NULL) {
        fprintf (stderr, "Out of memory in read_restart\n");
        rval = 1; goto CLEANUP;
    }

    CCutil_readstr (f, problabel, CCtsp_PROB_FILE_NAME_LEN);
    
    rval = fscanf (f, "%d%lf%d%lf\n", p_ncount, p_upbound, p_bbcount,
            p_branchzeit);
    if (rval <= 0) {
        perror (restart_name);
        fprintf (stderr, "fscanf from %s failed\n", restart_name);
        rval = 1; goto CLEANUP;
    }
    rval = read_bbtree (f, p_rootbbnode, bbnode_world);
    if (rval) {
        fprintf (stderr, "read_bbtree failed\n"); goto CLEANUP;
    }
    
    rval = fclose (f);
    if (rval) {
        perror (restart_name);
        fprintf (stderr, "fclose %s failed\n", restart_name);
        rval = 1; goto CLEANUP;
    }
    f = (FILE *) NULL;

    *p_problabel = problabel;
    
    rval = 0;

  CLEANUP:
    if (f != (FILE *) NULL) {
        fclose (f);
    }
    if (rval) {
        CC_IFFREE (problabel, char);
        free_tree (p_rootbbnode, bbnode_world);
    }
    return rval;
}

static int read_bbtree (FILE *f, tsp_bbnode **p_b, CCptrworld *bbnode_world)
{
    int rval = 0;
    int child0, child1;
    tsp_bbnode *b = (tsp_bbnode *) NULL;

    b = tsp_bbnode_alloc (bbnode_world);
    if (b == (tsp_bbnode *) NULL) {
        fprintf (stderr, "tsp_bbnode_alloc failed\n");
        rval = 1; goto CLEANUP;
    }
    init_bbnode (b);
    
    rval = fscanf (f, "%d%d%d%d%d%lf%lf\n", &b->status, &b->id,
                    &child0, &child1, &b->numtentative, &b->lowerbound,
                    &b->cputime);
    if (rval <= 0) {
        perror ("restart_file");
        fprintf (stderr, "fscanf failed reading restart file\n");
        rval = 1; goto CLEANUP;
    }
    b->workstatus = BB_IDLE;

    if (b->numtentative > 0) {
        rval = read_tentative_nodes (f, b->numtentative, &b->tentative_nodes,
                                     b, bbnode_world);
        if (rval) {
            fprintf (stderr, "read_tentative_nodes failed\n");
            goto CLEANUP;
        }
    }

    if (child0 != -1) {
        rval = read_bbtree (f, &(b->child0), bbnode_world);
        if (rval) goto CLEANUP;
        if (b->child0->id != child0) {
            fprintf (stderr, "syntax error in restart file\n");
            rval = 1; goto CLEANUP;
        }
        b->child0->parent = b;
    }
    if (child1 != -1) {
        rval = read_bbtree (f, &(b->child1), bbnode_world);
        if (rval) goto CLEANUP;
        if (b->child1->id != child1) {
            fprintf (stderr, "syntax error in restart file\n");
            rval = 1; goto CLEANUP;
        }
        b->child1->parent = b;
    }

    *p_b = b;
    rval = 0;

  CLEANUP:
    if (rval) {
        tsp_bbnode_free (bbnode_world, b);
    }
    return rval;
}

static int read_tentative_nodes (FILE *f, int count, tsp_tnode **list,
        tsp_bbnode *parent, CCptrworld *bbnode_world)
{
    int rval = 0;
    int i;
    int obtained = 0;
    tsp_tnode *s;
    tsp_bbnode *child0 = (tsp_bbnode *) NULL;
    tsp_bbnode *child1 = (tsp_bbnode *) NULL;

    *list = CC_SAFE_MALLOC (count, tsp_tnode);
    if (*list == (tsp_tnode *) NULL) {
        fprintf (stderr, "out of memory in read_tentative_nodes\n");
        rval = 1; goto CLEANUP;
    }

    for (obtained = 0; obtained < count; obtained++) {
        s = &((*list)[obtained]);
        child0 = tsp_bbnode_alloc (bbnode_world);
        if (child0 == (tsp_bbnode *) NULL) {
            fprintf (stderr, "tsp_bbnode_alloc failed\n");
            rval = 1; goto CLEANUP;
        }
        init_bbnode (child0);
        child1 = tsp_bbnode_alloc (bbnode_world);
        if (child1 == (tsp_bbnode *) NULL) {
            fprintf (stderr, "tsp_bbnode_alloc failed\n");
            tsp_bbnode_free (bbnode_world, child0);
            rval = 1; goto CLEANUP;
        }
        init_bbnode (child1);

        rval = fscanf (f, "%d%d%lf%lf %d%d%lf%lf\n",
                    &child0->status, &child0->id,
                    &child0->lowerbound, &child0->cputime,
                    &child1->status, &child1->id,
                    &child1->lowerbound, &child1->cputime);

            /*  BICO 22 APRIL 01 - Causes a problem       */
            /*  since workstatus is used to indicate that */
            /*  a tnode is pruned.  Should add a          */
            /*  pruned field to restart.                  */

        if (rval <= 0) {
            perror ("restart_file");
            fprintf (stderr, "fscanf failed reading tentative line\n");
            rval = 1; goto CLEANUP;
        }
        child0->tparent = s;
        child1->tparent = s;
        s->child0 = child0;
        s->child1 = child1;
        s->parent = parent;
    }
    rval = 0;

CLEANUP:

    if (rval) {
        for (i = 0; i < obtained; i++) {
            tsp_bbnode_free (bbnode_world, (*list)[i].child0);
            tsp_bbnode_free (bbnode_world, (*list)[i].child1);
        }
        CC_IFFREE (*list, tsp_tnode);
    }
    return rval;
}

static int add_children (tsp_bbnode **firstbbnode, tsp_bbnode *parent,
        int id0, int id1, double val0, double val1, int prune0, int prune1,
        CCptrworld *bbnode_world)
{
    int rval = 0;
    tsp_bbnode *child = (tsp_bbnode *) NULL;

    child = tsp_bbnode_alloc (bbnode_world);
    if (!child) {
        fprintf (stderr, "Failed to allocate child 0\n");
        rval = 1; goto CLEANUP;
    }
    
    init_bbnode (child);
    child->id = id0;
    child->lowerbound = val0;
    child->parent = parent;
    parent->child0 = child;

    if (val0 == CCtsp_LP_MAXDOUBLE) {
        printf ("Child 0 is infeasible\n"); fflush (stdout);
        child->status = BB_DONE;
    } else if (prune0) {
        printf ("Child 0 is pruned\n"); fflush (stdout);
        child->status = BB_DONE;
    } else {
        child->status = BB_NEEDS_CUTTING;
        insert_bbnode (firstbbnode, child);
    }

    child = tsp_bbnode_alloc (bbnode_world);
    if (!child) {
        fprintf (stderr, "Failed to allocate child 0\n");
        rval = 1; goto CLEANUP;
    }

    init_bbnode (child);
    child->id = id1;
    child->lowerbound = val1;
    child->parent = parent;
    parent->child1 = child;

    if (val1 == CCtsp_LP_MAXDOUBLE) {
        printf ("Child 1 is infeasible\n"); fflush (stdout);
        child->status = BB_DONE;
    } else if (prune1) {
        printf ("Child 1 is pruned\n"); fflush (stdout);
        child->status = BB_DONE;
    } else {
        insert_bbnode (firstbbnode, child);
        child->status = BB_NEEDS_CUTTING;
    }

CLEANUP:

    return rval;
}

static void free_tree (tsp_bbnode **bbnode, CCptrworld *bbnode_world)
{
    if (!(*bbnode))  return;
    free_tree (&((*bbnode)->child0), bbnode_world);
    free_tree (&((*bbnode)->child1), bbnode_world);
    tsp_bbnode_free (bbnode_world, *bbnode);
    *bbnode = (tsp_bbnode *) NULL;
}

static tsp_bbnode *select_bbnode (tsp_bbnode *firstbbnode, int verbose,
        int silent)
{
    double bestbound  = CCtsp_LP_MAXDOUBLE;
    double lowerbound = CCtsp_LP_MAXDOUBLE;
    tsp_bbnode *bestbbnode = (tsp_bbnode *) NULL;
    tsp_bbnode *b;
    tsp_tnode *t;
    int cutavail = 0;
    int branchavail = 0;
    int tcutavail = 0;
    int active = 0;
    int i;

    if (verbose) {
        bblist_info (firstbbnode, &cutavail, &tcutavail, &branchavail,
                     &active, &lowerbound);

        printf ("LOWER BOUND: %f   ACTIVE NODES: %d\n\n", lowerbound,
                active);
        fflush (stdout);
    }

    if (firstbbnode) {
        /* Find the best bbnode */
        for (b = firstbbnode; b; b = b->next) {
            if (b->workstatus == BB_IDLE && b->lowerbound < bestbound) {
                if (b->tentative_nodes != (tsp_tnode *) NULL) {
                    for (i = 0; i < b->numtentative; i++) {
                        t = &(b->tentative_nodes[i]);
                        if (t->child0->workstatus == BB_IDLE &&
                            t->child0->status == BB_NEEDS_TENTATIVE_CUTTING) {
                            bestbound = b->lowerbound;
                            bestbbnode = t->child0;
                            break;
                        } else if (t->child1->workstatus == BB_IDLE &&
                            t->child1->status == BB_NEEDS_TENTATIVE_CUTTING) {
                            bestbound = b->lowerbound;
                            bestbbnode = t->child1;
                            break;
                        }
                    }
                } else {
                    bestbound = b->lowerbound;
                    bestbbnode = b;
                }
            }
        }
    }

    if (verbose && !silent) {
        if (!bestbbnode) {
            printf ("No idle bbnodes\n"); fflush (stdout);
        } else {
            printf ("Selected bbnode:  id %d  lowerbound %.2f\n",
                    bestbbnode->id, bestbound);
            fflush (stdout);
            if (active > 1) {
                printf ("Remaining active bbnodes:\n");
                for (b = firstbbnode; b; b = b->next) {
                    if (b->id != bestbbnode->id) {
                        printf ("  id %d  lowerbound %.2f\n",
                                b->id, b->lowerbound);
                    }
                }
                fflush (stdout);
            }
        }
    }

    return bestbbnode;
}

static void insert_bbnode (tsp_bbnode **firstbbnode, tsp_bbnode *bbnode)
{
    if (!bbnode) return;

    bbnode->prev = (tsp_bbnode *) NULL;
    bbnode->next = *firstbbnode;
    if (*firstbbnode)  (*firstbbnode)->prev = bbnode;
    *firstbbnode = bbnode;
}

static void delete_bbnode (tsp_bbnode **firstbbnode, tsp_bbnode *bbnode)
{
    tsp_bbnode *next, *prev;

    if (!bbnode) return;

    bbnode->status = BB_DONE;
    bbnode->workstatus = BB_IDLE;
    if (*firstbbnode == bbnode) {
        *firstbbnode = bbnode->next;
    }
    prev = bbnode->prev;
    next = bbnode->next;
    if (prev) prev->next = next;
    if (next) next->prev = prev;
}

int CCtsp_easy_dfs_brancher (CCtsp_lp *lp, CCtsp_cutselect *sel, int depth,
        double *upbound, int *bbcount, int usecliques, int *besttour,
        int longedge_branching, int simple_branching, int silent,
        CCrandstate *rstate)
{
    int rval = 0;
    int ngot, prune, i;
    int *cyc = (int *) NULL;
    double val, bnd, szeit, st;
    double oldbound = lp->lowerbound;
    CCtsp_branchobj *b = (CCtsp_branchobj *) NULL;

    if (!lp->full_edges_valid) {
        fprintf (stderr, "CCtsp_easy_dfs_brancher needs valid extra edges\n");
        rval = 1; goto CLEANUP;
    }

    printf ("Node %d, Depth %d: ", *bbcount, depth); fflush (stdout);
    (*bbcount)++;
    if (!silent) {
        printf ("\n");
        CCtsp_print_branchhistory (lp);
    }

    rval = CCtsp_pricing_loop (lp, &bnd, silent, rstate);
    if (rval) {
        fprintf (stderr, "CCtsp_pricing_loop failed\n");  goto CLEANUP;
    }
    lp->lowerbound = bnd;
    lp->upperbound = *upbound;
    if (silent) {
        printf ("%.2f -> ", bnd); fflush (stdout);
    }

    if (lp->lowerbound >= lp->upperbound - 0.9) {
        rval = CCtsp_verify_lp_prune (lp, &prune, silent);
        if (rval) {
            fprintf (stderr, "CCtsp_verify_lp_prune failed\n"); goto CLEANUP;
        }
        if (prune) {
            printf ("PRUNE SEARCH: upperbound = %f\n", *upbound);
            fflush (stdout);
            rval = 0; goto CLEANUP;
        } else {
            fprintf (stderr, "exact pricing could not prune the search\n");
            rval = 1; goto CLEANUP;
        }
    }

    szeit = CCutil_zeit ();
    rval = CCtsp_cutting_loop (lp, sel, 0, silent, rstate);
    if (rval == 2) {
        rval = CCtsp_verify_infeasible_lp (lp, &prune, silent);
        if (rval) {
            fprintf (stderr ,"CCtsp_verify_infeasible_lp failed\n");
            goto CLEANUP;
        }
        if (prune) {
            printf ("PRUNE SEARCH - infeasible LP\n"); fflush (stdout);
            rval = 0; goto CLEANUP;
        } else {
            fprintf (stderr, "exact pricing did not verify an infeasible LP\n");
            rval = 1; goto CLEANUP;
        }
    } else if (rval) {
        fprintf (stderr, "CCtsp_cutting_loop failed\n"); goto CLEANUP;
    }

    if (silent) {
        printf ("%.2f (%.2f seconds)\n", lp->lowerbound, CCutil_zeit () - szeit);
        fflush (stdout);
    }

    if (!simple_branching && lp->lowerbound < lp->upperbound - 0.9) {
        CCutil_start_timer (&lp->stats.linkern);
        rval = CCtsp_call_x_heuristic (lp, &val, besttour, silent, rstate);
        if (rval) {
            fprintf (stderr, "CCtsp_call_x_heuristic failed\n");
            goto CLEANUP;
        }
        if (silent) {
            CCutil_stop_timer (&lp->stats.linkern, 0);
        } else {
            CCutil_stop_timer (&lp->stats.linkern, 1);
        }
        if (val < lp->upperbound) {
            printf ("New upperbound from x-heuristic: %.2f\n", val);
            lp->upperbound = val;
            *upbound = val;
            rval = CCtsp_dumptour (lp->graph.ncount, lp->dat, lp->perm,
                     lp->problabel, besttour, (char *) NULL, 0, silent);
            if (rval) {
                fprintf (stderr, "CCtsp_dumptour failed\n"); goto CLEANUP;
            }
        }
    }

    if (lp->lowerbound >= lp->upperbound - 0.9) {
        rval = CCtsp_verify_lp_prune (lp, &prune, silent);
        if (rval) {
            fprintf (stderr, "CCtsp_verify_lp_prune failed\n"); goto CLEANUP;
        }
        if (prune) {
            printf ("PRUNE SEARCH: upperbound = %f\n", *upbound);
            fflush (stdout);
            rval = 0; goto CLEANUP;
        } else {
            fprintf (stderr, "exact pricing could not prune the search\n");
            rval = 1; goto CLEANUP;
        }
    }

    oldbound = lp->lowerbound;
    if (!silent) {
        printf ("Find branch object ...\n"); fflush (stdout);
    }

    szeit = CCutil_zeit ();
    if (simple_branching) {
        rval = CCtsp_find_fast_branch (lp, &ngot, &b, &val, &cyc, usecliques,
                                  longedge_branching, silent);
        if (rval) {
            fprintf (stderr, "CCtsp_find_fast_branch failed\n");
            goto CLEANUP;
        }
    } else {
        rval = CCtsp_find_branch (lp, 1, &ngot, &b, &val, &cyc, usecliques, 
                                  longedge_branching, silent);
        if (rval) {
            fprintf (stderr, "CCtsp_find_branch failed\n");
            goto CLEANUP;
        }
    }
    st = CCutil_zeit () - szeit;

    if (ngot == 0) {
        printf ("TOUR FOUND: %.2f\n", val); fflush (stdout);
        if (val < *upbound) {
            *upbound = val;
            lp->upperbound = val;
            for (i = 0; i < lp->graph.ncount; i++) {
                besttour[i] = cyc[i];
            }
            rval = CCtsp_dumptour (lp->graph.ncount, lp->dat, lp->perm,
                         lp->problabel, besttour, (char *) NULL, 0, silent);
            if (rval) {
                fprintf (stderr, "CCtsp_dumptour failed\n"); goto CLEANUP;
            }
        }
        CC_IFFREE (cyc, int);
        rval = CCtsp_verify_lp_prune (lp, &prune, silent);
        if (rval) {
            fprintf (stderr, "CCtsp_verify_lp_prune failed\n");
            goto CLEANUP;
        }
        if (prune) {
            printf ("with new tour, the node can be pruned\n");
            fflush (stdout);
            rval = 0; goto CLEANUP;
        } else {
            fprintf (stderr, "could not verify the pruning\n");
            rval = 1; goto CLEANUP;
        }
    }


    /**** Down-Side Branch ****/


    if (b[0].ends[0] != -1) {
        printf ("Branch: set edge (%d, %d) to 0 (%.2f seconds)\n",
                     b[0].ends[0], b[0].ends[1], st);
        fflush (stdout);
        b[0].rhs = 0;
    } else {
        printf ("Branch: set clique <= 2 (%.2f seconds)\n", st);
        fflush (stdout);
        b[0].rhs = 2; b[0].sense = 'L';
    }
    rval = CCtsp_execute_branch (lp, &b[0], silent, rstate);
    if (rval == 2) {
        fprintf (stderr, "branched lp was infeasible\n");
        rval = CCtsp_verify_infeasible_lp (lp, &prune, silent);
        if (rval) {
            fprintf (stderr ,"CCtsp_verify_infeasible_lp failed\n");
            goto CLEANUP;
        }
        if (prune) {
            printf ("PRUNE SIDE - infeasible LP\n"); fflush (stdout);
            rval = 0; 
        } else {
            fprintf (stderr, "exact pricing did not verify an infeasible LP\n");
            rval = 1; goto CLEANUP;
        }
    } else if (rval) {
        fprintf (stderr, "CCtsp_execute_branch failed\n");
        rval = 1; goto CLEANUP;
    } else {
        rval = CCtsp_easy_dfs_brancher (lp, sel, depth + 1, upbound, bbcount,
                      usecliques, besttour, longedge_branching, simple_branching,
                      silent, rstate);
        if (rval) {
            fprintf (stderr, "CCtsp_easy_dfs_brancher failed\n"); goto CLEANUP;
        }
    }
    rval = CCtsp_execute_unbranch (lp, (CClp_warmstart *) NULL, silent, rstate);
    if (rval == 2) {
        fprintf (stderr, "branched lp was infeasible\n");
        rval = CCtsp_verify_infeasible_lp (lp, &prune, silent);
        if (rval) {
            fprintf (stderr ,"CCtsp_verify_infeasible_lp failed\n");
            goto CLEANUP;
        }
        if (prune) {
            printf ("PRUNE BOTH SIDES - infeasible unbranched LP\n");
            fflush (stdout);
            rval = 0; 

            lp->lowerbound = oldbound;
            CCtsp_free_branchobj (&b[0]);
            CC_IFFREE (b, CCtsp_branchobj);

            goto CLEANUP;
        } else {
            fprintf (stderr, "exact pricing did not verify an infeasible LP\n");
            rval = 1; goto CLEANUP;
        }
    } else if (rval) {
        fprintf (stderr, "CCtsp_execute_unbranch failed\n");
        goto CLEANUP; 
    }
    lp->lowerbound = oldbound;


    /**** Up-Side Branch ****/


    if (b[0].ends[0] != -1) {
        printf ("Branch: set edge (%d, %d) to 1 (depth %d)\n",
                     b[0].ends[0], b[0].ends[1], depth);
        b[0].rhs = 1;
    } else {
        printf ("Branch: set clique >= 4 (depth %d)\n", depth);
        b[0].rhs   = 4; b[0].sense = 'G';
    }
    fflush (stdout);
    rval = CCtsp_execute_branch (lp, &b[0], silent, rstate);
    if (rval == 2) {
        fprintf (stderr, "branched lp was infeasible\n");
        rval = CCtsp_verify_infeasible_lp (lp, &prune, silent);
        if (rval) {
            fprintf (stderr ,"CCtsp_verify_infeasible_lp failed\n");
            goto CLEANUP;
        }
        if (prune) {
            printf ("PRUNE SIDE - infeasible LP\n"); fflush (stdout);
            rval = 0; 
        } else {
            fprintf (stderr, "exact pricing did not verify an infeasible LP\n");
            rval = 1; goto CLEANUP;
        }
        rval = 0;
    } else if (rval) {
        fprintf (stderr, "CCtsp_execute_branch failed\n");
        rval = 1; goto CLEANUP;
    } else {
        rval = CCtsp_easy_dfs_brancher (lp, sel, depth + 1, upbound, bbcount,
                      usecliques, besttour, longedge_branching, simple_branching,
                      silent, rstate);
        if (rval) {
            fprintf (stderr, "CCtsp_easy_dfs_brancher failed\n"); goto CLEANUP;
        }
    }
    rval = CCtsp_execute_unbranch (lp, (CClp_warmstart *) NULL, silent, rstate);
    if (rval == 2) {
        fprintf (stderr, "unbranched lp was infeasible\n");
        rval = CCtsp_verify_infeasible_lp (lp, &prune, silent);
        if (rval) {
            fprintf (stderr ,"CCtsp_verify_infeasible_lp failed\n");
            goto CLEANUP;
        }
        if (prune) {
            printf ("NOTE - infeasible unbranched LP\n");
            fflush (stdout);
            rval = 0; 
            goto CLEANUP;
        } else {
            fprintf (stderr, "exact pricing did not verify an infeasible LP\n");
            rval = 1; goto CLEANUP;
        }
    } else if (rval) {
        fprintf (stderr, "CCtsp_execute_unbranch failed\n");
        goto CLEANUP; 
    }
    lp->lowerbound = oldbound;

    CCtsp_free_branchobj (&b[0]);
    CC_IFFREE (b, CCtsp_branchobj);

CLEANUP:

    return rval;
}

int CCtsp_do_interactive_branch (CCtsp_lp *lp, int silent, CCrandstate *rstate)
{
    int bend0, bend1, ch0, ch1, tbran, nseg, i;
    CCtsp_branchobj b;
    CCtsp_lpclique *c = (CCtsp_lpclique *) NULL;
    int *slist  = (int *) NULL;
    int rval = 0;

    CCtsp_init_branchobj (&b);

    printf ("Enter the (integer) id's for the two child nodes: ");
    fflush (stdout);
    scanf ("%d %d", &ch0, &ch1);

    printf ("Enter 0 if edge-branch, 1 if clique-branch (internal),\n");
    printf ("      2 if clique-branch (original): ");
    fflush (stdout);
    scanf ("%d", &tbran);

    if (tbran == 0) {
        printf ("Enter ends of branching edge (use neg if original): ");
        fflush (stdout);
        scanf ("%d %d", &bend0, &bend1);
        if (bend0 < 0) {
            if (bend1 >= 0) {
                fprintf (stderr, "both ends must be from the same order\n");
                rval = 1; goto CLEANUP;
            }
            for (i = 0; i < lp->graph.ncount; i++) {
                if (lp->perm[i] == -bend0) bend0 = i;
                if (lp->perm[i] == -bend1) bend1 = i;
            }
            printf ("Current Names of the Ends: %d %d\n", bend0, bend1);
            fflush (stdout);
        }
        b.ends[0] = bend0;
        b.ends[1] = bend1;
        b.rhs     = 1;
    } else if (tbran == 1) {
        printf ("Enter the number of segments in clique: ");
        fflush (stdout);
        scanf ("%d", &nseg);
        slist = CC_SAFE_MALLOC (2*nseg, int);
        if (!slist) {
            fprintf (stderr, "out of memory\n");
            rval = 1; goto CLEANUP;
        }
        printf ("Enter the ends of the segments: ");
        fflush (stdout);
        for (i = 0; i < nseg; i++) {
            scanf ("%d %d", &slist[2*i], &slist[2*i+1]);
        }
        c = CC_SAFE_MALLOC (1, CCtsp_lpclique);
        if (!c) {
            fprintf (stderr, "out of memory\n");
            CC_IFFREE (slist, int);
            rval = 1; goto CLEANUP;
        }
        rval = CCtsp_seglist_to_lpclique (nseg, slist, c);
        if (rval) {
            fprintf (stderr, "CCtsp_seglist_to_lpclique failed\n");
            goto CLEANUP;
        }
        CC_IFFREE (slist, int);
        b.clique = c;
        b.rhs    = 4;
        b.sense  = 'G';
        CCtsp_print_lpclique (b.clique);
    } else {
        printf ("Enter the number of nodes in clique: ");
        fflush (stdout);
        scanf ("%d", &nseg);
        slist = CC_SAFE_MALLOC (nseg, int);
        if (!slist) {
            fprintf (stderr, "out of memory\n");
            rval = 1; goto CLEANUP;
        }
        printf ("Enter the nodes in the clique: ");
        fflush (stdout);
        for (i = 0; i < nseg; i++) {
            scanf ("%d", &bend0);
            for (bend1 = 0; bend1 < lp->graph.ncount; bend1++) {
                if (lp->perm[bend1] == bend0) {
                    slist[i] = bend1;
                    break;
                }
            }
        }
        c = CC_SAFE_MALLOC (1, CCtsp_lpclique);
        if (!c) {
            fprintf (stderr, "out of memory\n");
            CC_IFFREE (slist, int);
            rval = 1; goto CLEANUP;
        }
        rval = CCtsp_array_to_lpclique (slist, nseg, c);
        if (rval) {
            fprintf (stderr, "CCtsp_seglist_to_lpclique failed\n");
            goto CLEANUP;
        }
        CC_IFFREE (slist, int);
        b.clique = c;
        b.rhs    = 4;
        b.sense  = 'G';
        CCtsp_print_lpclique (b.clique);
    }

    rval = CCtsp_splitprob (lp, &b, ch0, ch1, silent, rstate);
    if (rval) {
        fprintf (stderr, "CCtsp_splitprob failed\n");
        goto CLEANUP;
    }

    CCtsp_free_branchobj (&b);

CLEANUP:

    return rval;
}
