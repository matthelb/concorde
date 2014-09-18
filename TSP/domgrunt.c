#include <stdio.h>
#include "machdefs.h"
#include "util.h"
#include "tsp.h"

#define DOM_TASK_WAIT_SECONDS  (5)

typedef struct domgraph {
    int ncount;
    int ecount;
    int gid;
} domgraph;

typedef struct dominfo {
    int count;
} dominfo;

static int process_subproblem (int id, domgraph *G, dominfo *D);
static int receive_graph (CC_SFILE *s, domgraph *G);
static int send_dominos (CC_SFILE *s, dominfo *D);

static void init_domgraph (domgraph *G);
static void free_domgraph (domgraph *G);
static void init_dominfo (dominfo *D);
static void free_dominfo (dominfo *D);

static int wolf = 1;

int main (int ac, char **av)
{
    char *bosshost = (char *) NULL;
    double rtime = 0.0;
    int id = -1;
    int rval = 0;
    CC_SFILE *s = (CC_SFILE *) NULL;
    double szeit;
    char task, need, ready;
    domgraph G;
    dominfo D;

    init_domgraph (&G);
    init_dominfo (&D);

    if (ac != 2) {
        fprintf (stderr, "Usage: %s boss\n", av[0]);
        rval = 1; goto CLEANUP;
    }

    CCutil_printlabel ();

    bosshost = av[1];

    while (1) {
        s = CCutil_snet_open (bosshost, CCtsp_DOMINO_PORT);
        if (!s) {
            fprintf (stderr, "CCutil_snet_open failed\n");
            rval = 1;  goto CLEANUP;
        }

        rval = CCutil_swrite_char (s, CCtsp_DOMINO_RECEIVE);
        CCcheck_rval (rval, "CCutil_swrite_char failed (RECEIVE)");

        rval = CCutil_sread_char (s, &ready);
        CCcheck_rval (rval, "CCutil_sread_char failed (ready)");

        if (ready == CCtsp_DOMINO_YES) {
            rval = CCutil_swrite_int (s, G.gid);
            CCcheck_rval (rval, "CCutil_swrite_int failed (gid");
            rval = CCutil_swrite_int (s, id);
            CCcheck_rval (rval, "CCutil_swrite_int failed (id)");

            if (id != -1) {
                rval = CCutil_swrite_double (s, rtime);
                CCcheck_rval (rval, "CCutil_swrite_double failed (rtime)");

                rval = CCutil_sread_char (s, &need);
                CCcheck_rval (rval, "CCutil_sread_char failed (need)");

                if (need == CCtsp_DOMINO_YES) {
                    rval = send_dominos (s, &D);
                    CCcheck_rval (rval, "send_dominos failed");
                }
                free_dominfo (&D);
                id = -1;
            }

            rval = CCutil_sread_char (s, &task);
            CCcheck_rval (rval, "CCutil_sread_char failed (task)");

            switch (task) {
            case CCtsp_DOMINO_WAIT:
                sleep (DOM_TASK_WAIT_SECONDS);
                break;
            case CCtsp_DOMINO_GRAPH:
                rval = receive_graph (s, &G);
                CCcheck_rval (rval, "receive_graph failed");
            case CCtsp_DOMINO_WORK:
                rval = CCutil_sread_int (s, &id);
                CCcheck_rval (rval, "CCutil_sread_int failed (id)");
                break;
            case CCtsp_DOMINO_EXIT:
                printf ("Shutting down the domino grunt\n");
                fflush (stdout);
                goto CLEANUP; 
            }
  
            CCutil_sclose (s);
            s = (CC_SFILE *) NULL;

            if (id != -1) {
                printf ("PROCESSING node %d, graph %d\n", id, G.gid);
                fflush (stdout);

                szeit = CCutil_zeit ();

                rval = process_subproblem (id, &G, &D);
                CCcheck_rval (rval, "process_subproblem failed");

                rtime = CCutil_zeit () - szeit;
            }
        } else {
            free_dominfo (&D);
            id = -1;
            sleep (DOM_TASK_WAIT_SECONDS);
        }
    }

CLEANUP:

    if (s != (CC_SFILE *) NULL) {
        CCutil_sclose (s);
    }
    free_domgraph (&G);
    free_dominfo (&D);

    return rval;
}

static int process_subproblem (int id, domgraph *G, dominfo *D)
{
    int rval = 0;

    if (!G) {
        fprintf (stderr, "no graph\n");
        rval = 1;  goto CLEANUP;
    }

    if (!D) {
        fprintf (stderr, "no dominfo\n");
        rval = 1;  goto CLEANUP;
    }

    printf ("process %d (ncount %d, ecount %d)\n", id, G->ncount, G->ecount);
    fflush (stdout);

    D->count = wolf++;

CLEANUP:

    return rval;
}

static int receive_graph (CC_SFILE *s, domgraph *G)
{
    int rval = 0;

    if (!G) {
        fprintf (stderr, "no graph structure\n");
        rval = 1;  goto CLEANUP;
    }

    free_domgraph (G);

    rval = CCutil_sread_int (s, &G->gid);
    CCcheck_rval (rval, "CCutil_sread_int failed (gid)");

    rval = CCutil_sread_int (s, &G->ncount);
    CCcheck_rval (rval, "CCutil_sread_int failed (ncount)");

    rval = CCutil_sread_int (s, &G->ecount);
    CCcheck_rval (rval, "CCutil_sread_int failed (ecount)");

    printf ("Graph id %d, ncount = %d, ecount = %d\n",
               G->gid, G->ncount, G->ecount);
    fflush (stdout);


CLEANUP:

    return rval;
}

static int send_dominos (CC_SFILE *s, dominfo *D)
{
    int rval = 0;

    if (!D) {
        fprintf (stderr, "no dominfo to send\n");
        rval = 1;  goto CLEANUP;
    }

    rval = CCutil_swrite_int (s, D->count);
    CCcheck_rval (rval, "CCutil_swrite_int failed (count)");

CLEANUP:

    return rval;
}

static void init_domgraph (domgraph *G)
{
    if (G) {
        G->ncount = 0;
        G->ecount = 0;
        G->gid = -1;
    }
}

static void free_domgraph (domgraph *G)
{
    if (G) {
        G->ncount = 0;
        G->ecount = 0;
        G->gid = -1;
    }
}

static void init_dominfo (dominfo *D)
{
    if (D) {
        D->count = 0;
    }
}

static void free_dominfo (dominfo *D)
{
    if (D) {
        D->count = 0;
    }
}
