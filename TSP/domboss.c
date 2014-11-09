#include <stdio.h>
#include "tsp.h"
#include "util.h"

#define DOM_STAT_OPEN 0
#define DOM_STAT_DONE 1
#define DOM_STAT_WORK 2

typedef struct bossgraph {
    int ncount;
    int ecount;
    int gid;
} bossgraph;

typedef struct domlist {
    int count;
    int space;
    int *list;
} domlist;

static int receive_graph (CC_SFILE *s, bossgraph *G);
static int send_graph (CC_SFILE *s, bossgraph *G);
static int send_dominos (CC_SFILE *s, domlist *D);
static int receive_dominos (CC_SFILE *s, domlist *D);

static void init_bossgraph (bossgraph *G);
static void free_bossgraph (bossgraph *G);
static void init_domlist (domlist *D);
static void free_domlist (domlist *D);

int main (int ac, char **av)
{
    int i, id, gid, nremain, curloc, try;
    int rval = 0;
    CC_SPORT *lport = (CC_SPORT *) NULL;
    CC_SFILE *s;
    double rtime, cumtime = 0.0;
    int *status = (int *) NULL;
    char request;
    bossgraph G;
    domlist D;

    init_bossgraph (&G);
    init_domlist (&D);

    if (ac != 1) {
        fprintf (stderr, "Usage: %s (no arguments)\n", av[0]);
        rval = 1; goto CLEANUP;
    }

    printf ("BEGINNING DOMINO NET PROCESSING\n\n");
    fflush (stdout);

    lport = CCutil_snet_listen (CCtsp_DOMINO_PORT);
    if (lport == (CC_SPORT *) NULL) {                                           
        fprintf (stderr, "CCutil_snet_listen failed\n");
        rval = 1; goto CLEANUP;
    }

    while (1) {
        request = 0;
        do {
            s = CCutil_snet_receive (lport);
            if (!s) {
                fprintf (stderr, "CCutil_snet_receive failed, ignoring\n");
                continue;
            }

            if (CCutil_sread_char (s, &request)) {
                fprintf (stderr, "CCutil_sread_char failed, abort con\n");
                CCutil_sclose (s);
                continue;
            }

            switch (request) {
            case CCtsp_DOMINO_RECEIVE:
                rval = CCutil_swrite_char (s, CCtsp_DOMINO_NO);
                CCcheck_rval (rval, "CCutil_swrite_char failed (NO)");
                CCutil_sclose (s);
                break;
            case CCtsp_DOMINO_EXIT:
                printf ("Shutting down the domino boss\n"); fflush (stdout);
                CCutil_sclose (s);
                rval = 1;  goto CLEANUP;
            case CCtsp_DOMINO_GRAPH:
                break;
            case CCtsp_DOMINO_SEND:
                fprintf (stderr, "No graph, cannot send dominos\n");
                CCutil_sclose (s);
                rval = 1;  goto CLEANUP;
            default:
                fprintf (stderr, "Invalid request %c\n", request);
            }
        } while (request != CCtsp_DOMINO_GRAPH);

        rval = receive_graph (s, &G);
        CCcheck_rval (rval, "receive_graph failed");

        status = CC_SAFE_MALLOC (G.ncount, int);
        CCcheck_NULL (status, "out of memory for status array");

        for (i = 0; i < G.ncount; i++) status[i] = DOM_STAT_OPEN;

        nremain = G.ncount;
        curloc = 0;

        while (nremain) {
            request = 0;
            do {
                s = CCutil_snet_receive (lport);
                if (!s) {
                    fprintf (stderr, "CCutil_snet_receive failed, ignoring\n");
                    continue;
                }

                if (CCutil_sread_char (s, &request)) {
                    fprintf (stderr, "CCutil_sread_char failed, abort con\n");
                    CCutil_sclose (s);
                    continue;
                }

                switch (request) {
                case CCtsp_DOMINO_RECEIVE:
                    rval = CCutil_swrite_char (s, CCtsp_DOMINO_YES);
                    CCcheck_rval (rval, "CCutil_swrite_char failed (YES)");
                    break;
                case CCtsp_DOMINO_EXIT:
                    printf ("Shutting down the domino boss\n"); fflush (stdout);
                    CCutil_sclose (s);
                    rval = 1;  goto CLEANUP;
                case CCtsp_DOMINO_GRAPH:
                    fprintf (stderr, "Cannot receive new graph\n");
                    CCutil_sclose (s);
                    rval = 1;  goto CLEANUP;
                case CCtsp_DOMINO_SEND:
                    CCutil_sclose (s);
                    break;
                default:
                    fprintf (stderr, "Invalid request %c\n", request);
                }
            } while (request != CCtsp_DOMINO_RECEIVE);

            rval = CCutil_sread_int (s, &gid);
            if (rval) {
                fprintf (stderr, "CCutil_sread_int failed, abort con\n");
                rval = 0;
                goto CLOSE_CONN;
            }
            rval = CCutil_sread_int (s, &id);
            if (rval) {
                fprintf (stderr, "CCutil_sread_int failed, abort con\n");
                rval = 0;
                goto CLOSE_CONN;
            }
            if (id != -1) { 
                rval = CCutil_sread_double (s, &rtime);
                if (rval) {
                    fprintf (stderr, "CCutil_sread_double failed, abort con\n");
                    rval = 0;
                    goto CLOSE_CONN;
                }

                if (gid != G.gid) {
                    printf ("Finished node with graph %d, ignoring\n", gid);
                    fflush (stdout);

                    rval = CCutil_swrite_char (s, CCtsp_DOMINO_NO);
                    CCcheck_rval (rval, "CCutil_swrite_char failed (NO)");
                } else if (status[id] == DOM_STAT_DONE) {
                    rval = CCutil_swrite_char (s, CCtsp_DOMINO_NO);
                    CCcheck_rval (rval, "CCutil_swrite_char failed (NO)");

                    printf ("Finished completed node %d, ignoring\n", id);
                    fflush (stdout);
                } else {
                    rval = CCutil_swrite_char (s, CCtsp_DOMINO_YES);
                    CCcheck_rval (rval, "CCutil_swrite_char failed (YES)");

                    rval = receive_dominos (s, &D);
                    CCcheck_rval (rval, "receive_dominos failed");
        
                    status[id] = DOM_STAT_DONE;
                    cumtime += rtime;
                    nremain--;

                    printf ("DONE %3d:  %4.0f sec, %7.0f total, %d remaining ",
                            id, rtime, cumtime, nremain);
                    fflush (stdout);
                }
            }
            
            try = 0;
            while (status[curloc % G.ncount] == DOM_STAT_DONE && 
                   try < G.ncount) {
                curloc++;
                try++;
            }

            if (try < G.ncount) {
                if (gid != G.gid) {
                    rval = CCutil_swrite_char (s, CCtsp_DOMINO_GRAPH);
                    CCcheck_rval (rval, "CCutil_swrite_char failed (GRAPH)");
                    rval = send_graph (s, &G);
                    CCcheck_rval (rval, "send_graph failed");
                } else {
                    rval = CCutil_swrite_char (s, CCtsp_DOMINO_WORK);
                    CCcheck_rval (rval, "CCutil_swrite_char failed");
                }
        
                rval = CCutil_swrite_int (s, curloc % G.ncount);
                CCcheck_rval (rval, "CCutil_swrite_int failed");

                status[curloc % G.ncount] = DOM_STAT_WORK;

                printf ("  New = %d\n", curloc % G.ncount); fflush (stdout);

                curloc++;
            } else {
                printf ("\n"); fflush (stdout);
                rval = CCutil_swrite_char (s, CCtsp_DOMINO_WAIT);
                CCcheck_rval (rval, "CCutil_swrite_char failed (WAIT)");
            }
                
        CLOSE_CONN:

            CCutil_sclose (s);
        }

        printf ("\nFINISHED Graph %d: %.2f seconds\n", G.gid, cumtime);
        fflush (stdout);

        request = 0;
        do {
            s = CCutil_snet_receive (lport);
            if (!s) {
                fprintf (stderr, "CCutil_snet_receive failed, ignoring\n");
                continue;
            }

            if (CCutil_sread_char (s, &request)) {
                fprintf (stderr, "CCutil_sread_char failed, abort con\n");
                CCutil_sclose (s);
                continue;
            }

            switch (request) {
            case CCtsp_DOMINO_RECEIVE:
                rval = CCutil_swrite_char (s, CCtsp_DOMINO_NO);
                CCcheck_rval (rval, "CCutil_swrite_char failed (NO)");
                CCutil_sclose (s);
                break;
            case CCtsp_DOMINO_EXIT:
                printf ("Shutting down the domino boss\n"); fflush (stdout);
                CCutil_sclose (s);
                rval = 1;  goto CLEANUP;
            case CCtsp_DOMINO_GRAPH:
                fprintf (stderr, "Cannot receive new graph\n");
                CCutil_sclose (s);
                rval = 1;  goto CLEANUP;
            case CCtsp_DOMINO_SEND:
                break;
            default:
                fprintf (stderr, "Invalid request %c\n", request);
            }
        } while (request != CCtsp_DOMINO_SEND);

        printf ("Send the dominos\n"); fflush (stdout);

        rval = send_dominos (s, &D);
        CCcheck_rval (rval, "send_dominos failed");

        CCutil_sclose (s);

        free_domlist (&D);
        free_bossgraph (&G);
    }

CLEANUP:

    CCutil_snet_unlisten (lport);
    CC_IFFREE (status, int);

    return rval;
}

static int receive_graph (CC_SFILE *s, bossgraph *G)
{
    int rval = 0;

    if (!G) {
        fprintf (stderr, "no graph structure\n");
        rval = 1;  goto CLEANUP;
    }

    rval = CCutil_sread_int (s, &G->gid);
    CCcheck_rval (rval, "CCutil_sread_int failed (gid)");

    rval = CCutil_sread_int (s, &G->ncount);
    CCcheck_rval (rval, "CCutil_sread_int failed (ncount)");

    rval = CCutil_sread_int (s, &G->ecount);
    CCcheck_rval (rval, "CCutil_sread_int failed (ecount)");

    printf ("New graph %d (%d nodes, %d edges)\n", G->gid, G->ncount,
                                                   G->ecount);
    fflush (stdout);

CLEANUP:

    return rval;
}

static int send_graph (CC_SFILE *s, bossgraph *G)
{
    int rval = 0;

    if (!G) {
        fprintf (stderr, "no graph to send\n");
        rval = 1;  goto CLEANUP;
    }

    rval = CCutil_swrite_int (s, G->gid);
    CCcheck_rval (rval, "CCutil_swrite_int failed (graph id)");

    rval = CCutil_swrite_int (s, G->ncount);
    CCcheck_rval (rval, "CCutil_swrite_int failed (graph ncount)");

    rval = CCutil_swrite_int (s, G->ecount);
    CCcheck_rval (rval, "CCutil_swrite_int failed (graph ecount)");

CLEANUP:

    return rval;
}

static int send_dominos (CC_SFILE *s, domlist *D)
{
    int rval = 0;
    int i;

    rval = CCutil_swrite_int (s, D->count);
    CCcheck_rval (rval, "CCutil_swrite_int failed (count)");

    for (i = 0; i < D->count; i++) {
        rval = CCutil_swrite_int (s, D->list[i]);
        CCcheck_rval (rval, "CCutil_swrite_int failed (list)");
    }

    printf ("Dom List: %d\n", D->count);
    for (i = 0; i < D->count; i++) {
        printf ("%d ", D->list[i]);
    }
    printf ("\n"); fflush (stdout);

CLEANUP:

    return rval;
}

static int receive_dominos (CC_SFILE *s, domlist *D)
{
    int rval = 0;
    int count;

    rval = CCutil_sread_int (s, &count);
    CCcheck_rval (rval, "CCutil_sread_int failed (count)");

    if (D->count >= D->space) {
        if (D->space == 0) {
            D->list = CC_SAFE_MALLOC (100, int);
            CCcheck_NULL (D->list, "out of memory for domino list");
            D->space = 100;
        } else {
            rval =  CCutil_reallocrus_scale ((void **) &D->list, &D->space,
                        D->count + 1, 1.3, sizeof (int));
            CCcheck_rval (rval, "CCutil_reallocrus_scale failed");
        }
    }

    D->list[D->count] = count;
    D->count++;

CLEANUP:

    return rval;
}

static void init_bossgraph (bossgraph *G)
{
    if (G) {
        G->ncount = 0;
        G->ecount = 0;
        G->gid = -1;
    }
}

static void free_bossgraph (bossgraph *G)
{
    if (G) {
        G->ncount = 0;
        G->ecount = 0;
        G->gid = 0;
    }
}

static void init_domlist (domlist *D)
{
    if (D) {
        D->count = 0;
        D->list = (int *) NULL;
    }
}

static void free_domlist (domlist *D)
{
    if (D) {
        D->count = 0;
        D->space = 0;
        CC_IFFREE (D->list, int);
    }
}
