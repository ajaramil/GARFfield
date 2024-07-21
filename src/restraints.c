/* ----------------------------------------------------------------------
   GARFFIELD: Genetic-Algorithm Reactive Force FIELD parameter optimizer.

   Copyright (2012-2014) California Institute of Technology
   Andres Jaramillo-Botero (ajaramil@caltech.edu)
   http://www.wag.caltech.edu/home/ajaramil/GARFfield.html
------------------------------------------------------------------------- */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "restraints.h"
#include "mpi.h"

#ifdef DEBUG
# define DEBUG_PRINT(x) printf x
#else
# define DEBUG_PRINT(x)
#endif

#define MAX_LINE 512
#define MAX_TOKENS 512
#define MAX_TOKEN_LEN       512

extern int tokenize_string (char *, char ***);
extern void *scalloc (int, int, char *);
extern int debug_level;

/****************************************************************************
*  Function     : parse_restraints
*  Description  : Parses restraint file, if any
*  Parameters   : 
*  Effects      : No side effects, returns scalar value of dot product
****************************************************************************/

int parse_restraints (char *restraints_file, int rank) {

  FILE *fi;
  char **strtmp;
  char filen[MAX_LINE];
  char *line;
  int i, current_rst;

  line = (char *) scalloc (MAX_LINE, sizeof (char), "line");
  strtmp = (char **) scalloc (MAX_TOKENS, sizeof (char *), "strtmp");
  for (i = 0; i < MAX_TOKENS; i++)
    strtmp[i] = (char *) scalloc (MAX_TOKEN_LEN, sizeof (char), "tmpi");

  //rudimentary error checking
  sprintf(filen, "%s", restraints_file);
  fi = fopen (filen, "r");
  if (fi == NULL) {
    if (rank==0) fprintf (stderr, "ERROR: [%d] opening the restraints file %s...\n", rank, restraints_file);
    MPI_Abort (MPI_COMM_WORLD,1);
  }

  //count structures and restraints
  nrst = 0;
  while (fgets (line, MAX_LINE, fi) != NULL) 
    if (strstr (line, "RESTRAINT")) nrst++;
  if (rank == 0) 
    DEBUG_PRINT ((">> [%d] Found %d restraints in %s\n", rank, nrst, restraints_file));

  // Allocate space for restraints
  if (nrst > 0)
    rstrain = (restraint *) scalloc (nrst, sizeof (restraint), "rstrain");
  current_rst = 0;

  rewind(fi);
  i = -1;
  while (fgets (line, MAX_LINE, fi) != NULL) {
  // Restraint types: BOND, ANGLE, TORSION
// RESTRAINT Ca2_Si0 bond 1 2 2000 2000 1.2

    if (strstr (line, "RESTRAINT")) {
      tokenize_string (line, &strtmp);
      strcpy (rstrain[current_rst].name, strtmp[1]);
      rstrain[current_rst].atom[0] = atoi (strtmp[3]);
      rstrain[current_rst].atom[1] = atoi (strtmp[4]);

      // FORMAT: RESTRAINT structureName BOND atom1 atom2 Kstop Kstart distance
      if (strstr (strtmp[2], "BOND")) {
        rstrain[current_rst].nbody = 2;
        rstrain[current_rst].val = atof (strtmp[7]);
        rstrain[current_rst].f1 = atof (strtmp[5]);
        rstrain[current_rst].f2 = atof (strtmp[6]);
        if (rank == 0 && debug_level == 1) DEBUG_PRINT ((">> [%d] Bond restraint on %s: atoms %d %d %4.2f %4.2f %4.2f\n", rank, rstrain[current_rst].name,
          rstrain[current_rst].atom[0], rstrain[current_rst].atom[1],
          rstrain[current_rst].val, rstrain[current_rst].f1, rstrain[current_rst].f2));
      }
      else if (strstr (strtmp[2], "TS")) {
        rstrain[current_rst].nbody = 5;
        rstrain[current_rst].atom[2] = atoi (strtmp[5]);
        rstrain[current_rst].val = atof (strtmp[8]);
        rstrain[current_rst].f1 = atof (strtmp[6]);
        rstrain[current_rst].f2 = atof (strtmp[7]);
        if (rank == 0 && debug_level == 1) DEBUG_PRINT ((">> [%d] Transition State restraint on %s: atoms %d %d %d %4.2f %4.2f %4.2f\n", rank, rstrain[current_rst].name,
          rstrain[current_rst].atom[0], rstrain[current_rst].atom[1],
                rstrain[current_rst].atom[2], rstrain[current_rst].val, rstrain[current_rst].f1,
                rstrain[current_rst].f2));
      }
      // FORMAT: RESTRAINT structureName ANGLE atom1 atom2 atom3 Kstop Kstart angleValue
      else if (strstr (strtmp[2], "ANGLE")) {
        rstrain[current_rst].nbody = 3;
        rstrain[current_rst].atom[2] = atoi (strtmp[5]);
        rstrain[current_rst].val = atof (strtmp[8]);
        rstrain[current_rst].f1 = atof (strtmp[6]);
        rstrain[current_rst].f2 = atof (strtmp[7]);
        if (rank == 0 && debug_level == 1) DEBUG_PRINT ((">> [%d] Angle restraint on %s: atoms %d %d %d %4.2f %4.2f %4.2f\n", rank, rstrain[current_rst].name,
          rstrain[current_rst].atom[0], rstrain[current_rst].atom[1],
                rstrain[current_rst].atom[2], rstrain[current_rst].val, rstrain[current_rst].f1,
                rstrain[current_rst].f2));
      }
      // FORMAT: RESTRAINT structureName TORSION atom1 atom2 atom3 atom4 Kstop Kstart torsionValue
      else if (strstr (strtmp[2], "TORSION")) {
        rstrain[current_rst].nbody = 4;
        rstrain[current_rst].atom[2] = atoi (strtmp[5]);
        rstrain[current_rst].atom[3] = atoi (strtmp[6]);
        rstrain[current_rst].val = atof (strtmp[9]);
        rstrain[current_rst].f1 = atof (strtmp[7]);
        rstrain[current_rst].f2 = atof (strtmp[8]);
        if (rank == 0 && debug_level == 1) DEBUG_PRINT ((">> [%d] Torsion restraint on %s: atoms %d %d %d %d %4.2f %4.2f %4.2f\n",rank, rstrain[current_rst].name,
                rstrain[current_rst].atom[0], rstrain[current_rst].atom[1],
                rstrain[current_rst].atom[2], rstrain[current_rst].atom[3],
                rstrain[current_rst].val, rstrain[current_rst].f1, rstrain[current_rst].f2));
      }
      else {
        if (rank==0) printf ("ERROR: [%d] Unkown RESTRAINT in %s.\n\n", rank, rstrain[current_rst].name);
        MPI_Abort (MPI_COMM_WORLD, 1);
      }
      current_rst++;
    }
  }
  
  for (i = 0; i < MAX_TOKENS; i++)
    free(strtmp[i]);
  free (strtmp);
  free(line);
  fclose (fi);

  return current_rst;
}
