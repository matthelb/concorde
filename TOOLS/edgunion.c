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
/*                  A PROGRAM TO MERGE EDGE FILES                           */
/*                                                                          */
/*                             TSP CODE                                     */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: May 3, 1999                                                       */
/*                                                                          */
/*  SEE short decsribtion in usage ().                                      */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"

static int wantlen = 0;
static int nfiles = 0;
static char **filelist = (char **) NULL;
static char  *outfname = (char *) NULL;

int
    main (int ac, char **av);

static int
    parseargs (int ac, char **av);

static void
    usage (char *fname);


int main (int ac, char **av)
{
    int ncount, ecount, foundtour = 0, rval = 0;
    int *mytour;
    int *elist = (int *) NULL;
    int *elen  = (int *) NULL;
    double bestlen, *mylen;


    rval = parseargs (ac, av);
    if (rval) return 1;

    if (wantlen) {
        mytour = &foundtour;
        mylen  = &bestlen;
    } else {
        mytour = (int *) NULL;
        mylen  = (double *) NULL;
    }

    rval = CCutil_getedgelist_n (&ncount, filelist[0], &ecount, &elist,
                                 &elen, 0);
    if (rval) {
        fprintf (stderr, "CCutil_getedgelist_n failed\n"); goto CLEANUP;
    }
    CC_IFFREE (elist, int);
    CC_IFFREE (elen, int);

    rval = CCutil_edge_file_union (ncount, nfiles, filelist, &ecount, &elist,
                                   &elen, mytour, mylen);
    if (rval) {
        fprintf (stderr, "CCutil_edge_file_union failed\n"); goto CLEANUP;
    }
    printf ("Merged Edge List: %d edges\n", ecount); fflush (stdout);
      
    if (outfname) {
        rval = CCutil_writeedges_int (ncount, outfname, ecount, elist,
                                      elen, 0);
        if (rval) {
            fprintf (stderr, "CCutil_writeedges_int failed\n"); goto CLEANUP;
        }
    }

    if (wantlen) {
        if (foundtour == 1) {
            printf ("Best Tour:  %.0f\n", bestlen); fflush (stdout);
        } else {
            printf ("No tours\n"); fflush (stdout);
        }
    }


CLEANUP:

    CC_IFFREE (elist, int);
    CC_IFFREE (elen, int);

    return rval;
}

static int parseargs (int ac, char **av)
{
    int c;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "o:t", &boptind, &boptarg)) != EOF) {
        switch (c) {
        case 't':
            wantlen = 1;
            break;
        case 'o':
            outfname = boptarg;
            break;
        case CC_BIX_GETOPT_UNKNOWN:
        case '?':
        default:
            usage (av[0]);
            return 1;
        }
    }

    if (ac <= boptind) {
        usage (av[0]);
        return 1;
    }

    nfiles = ac - boptind;
    filelist = &(av[boptind]);
    return 0;
}

static void usage (char *fname)
{
    fprintf (stderr, "Usage: %s [-flags below] edge_files\n", fname);
    fprintf (stderr, "   -o f  write merged edge file\n");
    fprintf (stderr, "   -t    print best tour len (if any are tours)\n");
}
