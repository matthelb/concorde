#include <stdio.h>
#include "machdefs.h"
#include "util.h"
#include "tsp.h"

static int
    receive_dominos (CC_SFILE *s),
    parseargs (int ac, char **av);

static void
    usage (char *f);

static int graphid = 1;
static int sendgraph = 0;
static int stopboss = 0;
static int getdominos = 0;
static char *hostname = (char *) NULL;

int main (int ac, char **av)
{
    CC_SFILE *s = (CC_SFILE *) NULL;
    int rval = 0;

    rval = parseargs (ac, av);
    if (rval) return 0;

    CCutil_printlabel ();

    s = CCutil_snet_open (hostname, CCtsp_DOMINO_PORT);
    if (!s) {
        fprintf (stderr, "CCutil_snet_open failed\n");
        rval = 1;  goto CLEANUP;
    }

    if (stopboss) {
        printf ("stop the boss\n");  fflush (stdout);
        rval = CCutil_swrite_char (s, CCtsp_DOMINO_EXIT);
        CCcheck_rval (rval, "CCutil_swrite_char failed (EXIT)");
    } else if (sendgraph) {
        printf ("Send graph id = %d\n", graphid);  fflush (stdout);
        rval = CCutil_swrite_char (s, CCtsp_DOMINO_GRAPH);
        CCcheck_rval (rval, "CCutil_swrite_char failed (GRAPH)");
        rval = CCutil_swrite_int (s, graphid);
        CCcheck_rval (rval, "CCutil_swrite_int failed");
        rval = CCutil_swrite_int (s, 2*graphid);
        CCcheck_rval (rval, "CCutil_swrite_int failed");
        rval = CCutil_swrite_int (s, 3*graphid);
        CCcheck_rval (rval, "CCutil_swrite_int failed");
    } else if (getdominos) {
        printf ("get the domino list\n");  fflush (stdout);
        rval = CCutil_swrite_char (s, CCtsp_DOMINO_SEND);
        CCcheck_rval (rval, "CCutil_swrite_char failed (SEND)");
        rval = receive_dominos (s);
        CCcheck_rval (rval, "receive_dominos failed");
    } else {
        printf ("Nothing to do.\n"); fflush (stdout);
    }

CLEANUP:

    if (s) CCutil_sclose (s);
    return rval;
}

static int parseargs (int ac, char **av)
{
    int c;                                                                  
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "xg:d", &boptind, &boptarg)) != EOF) {
        switch (c) {                                                        
        case 'x':
            stopboss = 1;
            break;
        case 'g':                                                           
            sendgraph = 1;
            graphid = atoi(boptarg);
            break;                                                          
        case 'd':              
            getdominos = 1;
            break;
        default:
            usage (av[0]);
            return 1;
        }
    }

    if (boptind >= ac) {
        usage (av[0]);
        return 1;
    }
    hostname = av[boptind++];
    if (boptind != ac) {
        usage (av[0]);
        return 1;
    }

    return 0;
}

static int receive_dominos (CC_SFILE *s)
{
    int rval = 0;
    int i, count;
    int *list = (int *) NULL;

    rval = CCutil_sread_int (s, &count);
    CCcheck_rval (rval, "CCutil_sread_int failed (count)");

    list = CC_SAFE_MALLOC (count, int);
    CCcheck_NULL (list, "out memory for list");

    for (i = 0; i < count; i++) {
        rval = CCutil_sread_int (s, &list[i]);
        CCcheck_rval (rval, "CCutil_sread_int failed (list)");
    }

    printf ("Dom List: %d\n", count);
    for (i = 0; i < count; i++) {
        printf ("%d ", list[i]);
    }
    printf ("\n"); fflush (stdout);

CLEANUP:

    CC_IFFREE (list, int);
    return rval;
}

static void usage (char *f)
{
    fprintf (stderr, "Usage: %s [-see below-] hostname\n", f);
    fprintf (stderr, "   -d   tells the boss to send dominos\n");
    fprintf (stderr, "   -g # sends graph is id #\n");
    fprintf (stderr, "   -x   tells the boss to exit\n");
}
