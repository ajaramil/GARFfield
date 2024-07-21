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

#include "cg_ffield.h"

extern int tokenize_string (char *, char ***);
extern void *scalloc (int, int, char *);

#ifdef DEBUG
# define DEBUG_PRINT(x) printf x
#else
# define DEBUG_PRINT(x) 
#endif

/****************************************************************************
*  Function     : Write_Force_Field_CG
*  Description  : Writes contents of the ffid_pointer into a ffield file 
                  with format of type Coarse-Grain, associated with user-selected 
                  force field (-F cg)
*  Parameters   : Name of force field file, ffield parameters structure, 
                  forcefield array, and MPI process rank
*  Effects      : Writes a new force field file per evaluation pass, per MPI
                  rank
****************************************************************************/

void Write_Force_Field_CG (char *ffield_file, cg_interaction * ffdata, double **ffid_pointer, int rank)
{
  int i, token;
  FILE *fo;

  /* open force field for writing */
  if ((fo = fopen (ffield_file, "w")) == NULL) {
    fprintf (stderr, "ERROR: [%d] Cannot write to %s...\n", rank,ffield_file);
    exit (1);
  }

  token = 0;
  fprintf (fo, "CG\n\n");
  fprintf (fo, "%3d       ! Number of CG atom types\n", ffdata->num_atom_types);

  /* print atom type section */
  for (i = 0; i < ffdata->num_atom_types; i++) 
    fprintf (fo, "%d %-2s\n", ffdata->cgatypes[i].type, ffdata->cgatypes[i].label);
  /* print pair type section */
  fprintf (fo, "\n%3d       ! Number of CG pairwise parameters\n", ffdata->num_pair_types);
  for (i = 0; i < ffdata->num_pair_types; i++) {
    fprintf (fo, "%d\t%d\t", ffdata->cgptypes[i].type1, ffdata->cgptypes[i].type2);
    fprintf (fo, "%12.6f", *ffid_pointer[token++]);
    fprintf (fo, "%12.6f", *ffid_pointer[token++]);
    fprintf (fo, "%12.6f\n", *ffid_pointer[token++]);
  }
  fclose (fo);
}

/****************************************************************************
*  Function     : Read_Force_Field_CG
*  Description  : Reads and parses contents of the user-provided ffield file
                  into the appropriate coarse-grain force field data structure format,
                  associated with user-selected force field (-F cg)
*  Parameters   : Name of force field file, ffield parameters structure,
                  forcefield array, and MPI process rank
*  Effects      : Reads force field parameters per GA evaluation pass, per MPI 
                  rank
****************************************************************************/

int Read_Force_Field_CG (char *ffield_file, cg_interaction * ffdata, double ***ffid_pointer,int rank)
{
  FILE *fi;
  char *s;
  char **tmp;
  int a, token, i, j;
  double val;

  /* open force field file */
  if ((fi = fopen (ffield_file, "r")) == NULL) {
    fprintf (stderr, "ERROR: [%d] opening the force field file! terminating...\n",rank);
    exit (1);
  }

  s = (char *) scalloc (MAX_LINE, sizeof (char), "s3");
  tmp = (char **) scalloc (MAX_TOKENS, sizeof (char *), "tmp3");
  for (i = 0; i < MAX_TOKENS; i++)
    tmp[i] = (char *) scalloc (MAX_TOKEN_LEN, sizeof (char), "tmp3i");

  /* reading first header comment */
  fgets (s, MAX_LINE, fi);
  /* reading empty line */
  fgets (s, MAX_LINE, fi);
  /* line 2 is number of cg atom types */
  fgets (s, MAX_LINE, fi);
  token = tokenize_string (s, &tmp);
  /* reading the number of cg atom types*/
  ffdata->num_atom_types = atoi (tmp[0]);
  if (ffdata->num_atom_types < 1) {
    fprintf (stderr, "WARNING: [%d] number of cg atoms in ffield file is 0!\n",rank);
    exit (1);
  }
  for (i = 0; i < ffdata->num_atom_types; i++)
    fgets (s, MAX_LINE, fi);

  /* reading empty line */
  fgets (s, MAX_LINE, fi);
  /* Number of cg pair types */
  fgets (s, MAX_LINE, fi);
  token = tokenize_string (s, &tmp);
  /* reading the number of cg pair types*/
  ffdata->num_pair_types = atoi (tmp[0]);
  if (ffdata->num_pair_types < 1) {
    fprintf (stderr, "WARNING: [%d] number of cg pair types in ffield file is 0!\n",rank);
    exit (1);
  }

  // allocate memory
  ffdata->cgatypes = (cg_atom_types *) scalloc (ffdata->num_atom_types, sizeof (cg_atom_types), "cg_atom_types");
  ffdata->cgptypes = (cg_pair_types *) scalloc (ffdata->num_pair_types, sizeof (cg_pair_types), "cg_pair_types");
  *ffid_pointer = (double **) scalloc (ffdata->num_pair_types * 5, sizeof (double *), "ffid_pointer");

  // return to start, skip globals and comments for general and sbp
  rewind (fi);
  for (i = 0; i < 3; i++)
    fgets (s, MAX_LINE, fi);

  for (i = 0; i < ffdata->num_atom_types; i++) {
    fgets (s, MAX_LINE, fi);
    token = tokenize_string (s, &tmp);
    ffdata->cgatypes[i].type = atoi(tmp[0]);
    for (j = 0; j < (int) strlen (tmp[1]); ++j)
      ffdata->cgatypes[i].label[j] = tmp[1][j];
  }

  fgets (s, MAX_LINE, fi);
  fgets (s, MAX_LINE, fi);
  a = 0;
  for (i = 0; i < ffdata->num_pair_types; i++) {
    fgets (s, MAX_LINE, fi);
    token = tokenize_string (s, &tmp);
    ffdata->cgptypes[i].type1 = atoi(tmp[0]);
    ffdata->cgptypes[i].type2 = atoi(tmp[1]);

    val = (double) atof (tmp[2]);
    ffdata->cgptypes[i].p1 = val;
    (*ffid_pointer)[a++] = &ffdata->cgptypes[i].p1;
    val = (double) atof (tmp[3]);
    ffdata->cgptypes[i].p2 = val;
    (*ffid_pointer)[a++] = &ffdata->cgptypes[i].p2;
    val = (double) atof (tmp[4]);
    ffdata->cgptypes[i].p3 = val;
    (*ffid_pointer)[a++] = &ffdata->cgptypes[i].p3;
  }

  ffdata->total_cg_parameters=a;
  if (rank == 0 && debug_level == 1) DEBUG_PRINT (("FFIELD: [%d] Read %d CG parameters\n", rank,ffdata->total_cg_parameters));

  // close file
  fclose (fi);

  /* deallocate helper storage */
  for (i = 0; i < MAX_TOKENS; i++)
    free (tmp[i]);
  free (tmp);
  free (s);

  if (rank == 0 && debug_level == 1) DEBUG_PRINT ((">> [%d] Finished reading force field\n",rank));

  return a;
}
