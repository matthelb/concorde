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
/*               A PROGRAM TO ADD LENGTHS TO A FILE OF EDGES                */
/*                                                                          */
/*                                TSP CODE                                  */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: October 5, 1999                                                   */
/*                                                                          */
/*  SEE short decsribtion in usage ().                                      */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"

static char *tspfname  = (char *) NULL;
static char *edgefname = (char *) NULL;
static char *outfname  = (char *) NULL;
static int eformat = 0;
static int binary_out = 0;
static int seed = 0;



int
    main (int ac, char **av);

static int
    parseargs (int ac, char **av),
    get_edges (char *fname, int thirdfield, int ncount, int *ecount,
        int **elist);

static void
    usage (char *fname);


int main (int ac, char **av)
{
    int ncount, ecount, rval = 0;
    int *elist = (int *) NULL;
    CCdatagroup dat;
    CCrandstate rstate;

    CCutil_init_datagroup (&dat);

    seed = (int) CCutil_real_zeit ();

    rval = parseargs (ac, av);
    if (rval) return 1;

    if ((!edgefname || !tspfname)) {
        usage (av[0]);
        return 1;
    }

    CCutil_sprand (seed, &rstate);

    rval = CCutil_gettsplib (tspfname, &ncount, &dat);
    if (rval) {
        fprintf (stderr, "CCutil_gettsplib failed\n"); goto CLEANUP;
    }

    rval = get_edges (edgefname, eformat, ncount, &ecount, &elist);
    if (rval) {
        fprintf (stderr, "get_edges failed\n"); goto CLEANUP;
    }

    printf ("Number of Edges: %d\n", ecount); fflush (stdout);
    if (outfname) {
        rval = CCutil_writeedges (ncount, outfname, ecount, elist, &dat,
                                  binary_out);
        if (rval) {
            fprintf (stderr, "CCutil_writeedges failed\n");
        }
    }

CLEANUP:

    CC_IFFREE (elist, int);
    CCutil_freedatagroup (&dat);

    return rval;
}

static int parseargs (int ac, char **av)
{
    int c;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "beE:o:T:", &boptind, &boptarg)) != EOF) {
        switch (c) {
        case 'b':
            binary_out = 1;
            break;
        case 'e':
            eformat = 1;
            break;
        case 'E':
            edgefname = boptarg;
            break;
        case 'o':
            outfname = boptarg;
            break;
        case 'T':
            tspfname = boptarg;
            break;
        case CC_BIX_GETOPT_UNKNOWN:
        case '?':
        default:
            usage (av[0]);
            return 1;
        }
    }

    return 0;
}

static void usage (char *fname)
{
    fprintf (stderr, "Usage: %s [-flags below]\n", fname);
    fprintf (stderr, "   -b    write a binary file\n");
    fprintf (stderr, "   -e    the edge file has floating point third field (ignore it)\n");
    fprintf (stderr, "   -E f  file with list of edges (default is no length field\n");
    fprintf (stderr, "   -T f  TSPLIB file to specify lengths\n");
    fprintf (stderr, "   -o f  output file (for the edge list)\n");
    fprintf (stderr, " NOTE: -E and -T must be specified\n"); 
}


static int get_edges (char *fname, int thirdfield, int ncount, int *ecount,
        int **elist)
{
    FILE *f = (FILE *) NULL;
    int i, k;
    int rval = 0;

    f = fopen (fname, "r");
    if (f == (FILE *) NULL) {
        perror (fname);
        fprintf (stderr, "Unable to open %s for input\n", fname);
        return 1;
    }

    k = CCutil_readint (f);
    if (k != ncount) {
        fprintf (stderr, "TSP file and edge file do not match\n");
        rval = 1; goto CLEANUP;
    }

    *ecount = CCutil_readint (f);
    *elist = CC_SAFE_MALLOC(2 * (*ecount), int);
    if (!(*elist)) {
        fprintf (stderr, "out of memory n get_edges\n");
        rval = 1; goto CLEANUP;
    }

    for (i = 0; i < *ecount; i++) {
        (*elist)[2*i] = CCutil_readint (f);
        (*elist)[2*i+1] = CCutil_readint (f);
        if (thirdfield) {
            fscanf (f, "%*f");
        }
    }

CLEANUP:

    fclose (f);
    return rval;
}
