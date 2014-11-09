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
/*         A PROGRAM TO COMPUTE CONVERT A TOUR FILE TO AN EDGE FILE         */
/*                 or CONVERT TO A TSPLIB TOUR FILE                         */
/*                                                                          */
/*                                TSP CODE                                  */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: January 24, 2001                                                  */
/*        May 7, 2002 (bico)                                                */
/*                                                                          */
/*  SEE short decsribtion in usage ().                                      */
/*                                                                          */
/****************************************************************************/

#include "machdefs.h"
#include "util.h"

static char *cycfname = (char *) NULL;
static char *tspfname = (char *) NULL;
static char *outfname = (char *) NULL;
static int seed = 0;
static int simpletour = 0;
static int tsplibtour = 0;


int
    main (int ac, char **av);

static int
    parseargs (int ac, char **av);

static char
    *get_problabel (const char *probloc);

static void
    usage (char *fname);


int main (int ac, char **av)
{
    int ncount, rval = 0;
    int *tour = (int *) NULL;
    double val;
    CCdatagroup dat;
    CCrandstate rstate;
    char *name = (char *) NULL;
    FILE *out = (FILE *) NULL;

    CCutil_init_datagroup (&dat);

    seed = (int) CCutil_real_zeit ();

    rval = parseargs (ac, av);
    if (rval) return 1;

    if (!tspfname) {
        fprintf (stderr, "No TSPLIB file specified for edge lengths\n");
        usage (av[0]);
        return 1;
    }

    CCutil_sprand (seed, &rstate);

    rval = CCutil_gettsplib (tspfname, &ncount, &dat);
    if (rval) {
        fprintf (stderr, "CCutil_gettsplib failed\n"); goto CLEANUP;
    }

    tour = CC_SAFE_MALLOC (ncount, int);
    if (!tour) {
        fprintf (stderr, "out of memory in main\n");
        rval = 1; goto CLEANUP;
    }

    if (simpletour) {
        rval = CCutil_getcycle (ncount, cycfname, tour, 0);
        if (rval) {
            fprintf (stderr, "CCutil_getcycle_tsplib failed\n");
            goto CLEANUP;
        }
    } else {
        rval = CCutil_getcycle_tsplib (ncount, cycfname, tour);
        if (rval) {
            fprintf (stderr, "CCutil_getcycle_tsplib failed\n");
            goto CLEANUP;
        }
    }

    CCutil_cycle_len (ncount, &dat, tour, &val);
    printf ("Tour Length: %.0f\n", val); fflush (stdout);

    if (outfname) {
        int i;
        out = fopen (outfname, "w");
        if (!out) {
            fprintf (stderr, "could not open %s for writing\n", outfname);
            rval = 1; goto CLEANUP;
        }
        
        if (tsplibtour) {
            name = get_problabel (tspfname);
            fprintf (out, "NAME: %s\n", name);
            fprintf (out, "COMMENT: Tour length %.0f\n", val);
            fprintf (out, "TYPE: TOUR\n");
            fprintf (out, "DIMENSION: %d\n", ncount);
            fprintf (out, "TOUR_SECTION\n");
            for (i = 0; i < ncount; i++) {
                fprintf (out, "%d\n", tour[i] + 1);
            }
            fprintf (out, "-1\n");
            fprintf (out, "EOF\n");
        } else {
            fprintf (out, "%d %d\n", ncount, ncount);
            for (i = 0; i < ncount - 1; i++) {
                fprintf (out, "%d %d %d\n", tour[i], tour[i+1],
                         CCutil_dat_edgelen (tour[i], tour[i+1], &dat));
            }
            fprintf (out, "%d %d %d\n", tour[0], tour[ncount-1],
                         CCutil_dat_edgelen (tour[0], tour[ncount-1], &dat));
        }
    }

CLEANUP:

    CC_IFFREE (tour, int);
    CC_IFFREE (name, char);
    CCutil_freedatagroup (&dat);
    if (out) fclose (out);

    return rval;
}

static int parseargs (int ac, char **av)
{
    int c;
    int boptind = 1;
    char *boptarg = (char *) NULL;

    while ((c = CCutil_bix_getopt (ac, av, "o:StT:", &boptind, &boptarg)) != EOF) {
        switch (c) {
        case 'o':
            outfname = boptarg;
            break;
        case 'S':
            tsplibtour = 1;
            break;
        case 't':
            simpletour = 1;
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

    if (boptind < ac) {
        cycfname = av[boptind++];
    } else {
        usage (av[0]);
        return 1;
    }

    return 0;
}

static char *get_problabel (const char *probloc)
{
    const char *p;
    const char *problabel = probloc;
    char *probcopy = (char *) NULL;
    char *p2;

    p = CCutil_strrchr_c (problabel, ':');
    if (p != (const char *) NULL) problabel = p+1;
    p = CCutil_strrchr_c (problabel, '/');
    if (p != (const char *) NULL) problabel = p+1;
    probcopy = CCutil_strdup (problabel);
    if (probcopy == (char *) NULL) return (char *) NULL;
    p2 = CCutil_strchr (probcopy, '.');
    if (p2 != (char *) NULL) *p2 = '\0';
    return probcopy;
}

static void usage (char *fname)
{
    fprintf (stderr, "Usage: %s [-flags below] -T TSPLIB_file tour_file\n", fname);
    fprintf (stderr, "   -t    tour file in concorde format (default TSPLIB)\n");
    fprintf (stderr, "   -o f  output file (for the edge list)\n");
    fprintf (stderr, "   -S    write a TSPLIB tour file with tour length\n");
}
