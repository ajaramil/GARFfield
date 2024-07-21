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

#include "eff_ffield.h"

extern int tokenize_string (char *, char ***);
void *scalloc (int, int, char *);

#ifdef DEBUG
# define DEBUG_PRINT(x) printf x
#else
# define DEBUG_PRINT(x) 
#endif

/****************************************************************************
*  Function     : Write_Force_Field_Eff
*  Description  : Writes contents of the ffid_pointer into a ffield file 
                  with format of type eFF, associated with user-selected 
                  force field (-F eff)
*  Parameters   : Name of force field file, ffield parameters structure, 
                  forcefield array, and MPI process rank
*  Effects      : Writes a new force field file per evaluation pass, per MPI
                  rank
****************************************************************************/

void Write_Force_Field_Eff (char *ffield_file, eff_interaction * ffdata, double **ffid_pointer, int rank)
{
  int i, token;
  FILE *fo;

  /* open force field for writing */
  if ((fo = fopen (ffield_file, "w")) == NULL) {
    fprintf (stderr, "ERROR: Cannot write to %s...\n", ffield_file);
    exit (1);
  }

  token = 0;
  fprintf (fo, "eFF-ECP\n");
  fprintf (fo, "%3d       ! Number of ECP parameters\n", ffdata->num_atom_types);

  for (i = 0; i < ffdata->num_atom_types; i++) {
//    fprintf (fo, " %-2d", ffdata->ecp[i].atype);
    fprintf (fo, " %-2s", ffdata->ecp[i].name);
    fprintf (fo, " %-2s", ffdata->ecp[i].ecptype);
    if (strcmp (ffdata->ecp[i].ecptype, "s") == 0) {
      fprintf (fo, "%12.6f", *ffid_pointer[token++]);
      fprintf (fo, "%12.6f", *ffid_pointer[token++]); 
      fprintf (fo, "%12.6f", *ffid_pointer[token++]); 
      fprintf (fo, "%12.6f\n", *ffid_pointer[token++]);
      token++;
      token++;
      token++;
      token++;
//      token++;
    }
    else if ((strcmp (ffdata->ecp[i].ecptype, "p") == 0) || (strcmp (ffdata->ecp[i].ecptype, "h") == 0)) {
      fprintf (fo, "%12.6f", *ffid_pointer[token++]);
      fprintf (fo, "%12.6f", *ffid_pointer[token++]);
      fprintf (fo, "%12.6f", *ffid_pointer[token++]);
      fprintf (fo, "%12.6f", *ffid_pointer[token++]);
      fprintf (fo, "%12.6f", *ffid_pointer[token++]);
      fprintf (fo, "%12.6f\n", *ffid_pointer[token++]);
      token++;
      token++;
//      token++;
    }
    else if (strcmp (ffdata->ecp[i].ecptype, "x") == 0) {
//      fprintf (fo, "%12.6f", *ffid_pointer[token++]);
      fprintf (fo, "%12.6f", *ffid_pointer[token++]);
      fprintf (fo, "%12.6f", *ffid_pointer[token++]);
      fprintf (fo, "%12.6f", *ffid_pointer[token++]);
      fprintf (fo, "%12.6f", *ffid_pointer[token++]);
      fprintf (fo, "%12.6f", *ffid_pointer[token++]);
      fprintf (fo, "%12.6f\n", *ffid_pointer[token++]);
      token++;
      token++;
    }
    else if (strcmp (ffdata->ecp[i].ecptype, "f") == 0) {
//      fprintf (fo, "%12.6f", *ffid_pointer[token++]);
      fprintf (fo, "%12.6f", *ffid_pointer[token++]);
      fprintf (fo, "%12.6f", *ffid_pointer[token++]);
      fprintf (fo, "%12.6f", *ffid_pointer[token++]);
      fprintf (fo, "%12.6f", *ffid_pointer[token++]);
      fprintf (fo, "%12.6f", *ffid_pointer[token++]);
      fprintf (fo, "%12.6f", *ffid_pointer[token++]);
      fprintf (fo, "%12.6f", *ffid_pointer[token++]);
      fprintf (fo, "%12.6f\n", *ffid_pointer[token++]);
    }
  }
  fclose (fo);
}

/****************************************************************************
*  Function     : Read_Force_Field_Eff
*  Description  : Reads and parses contents of the user-provided ffield file
                  into the appropriate eFF force field data structure format,
                  associated with user-selected force field (-F eff)
*  Parameters   : Name of force field file, ffield parameters structure,
                  forcefield array, and MPI process rank
*  Effects      : Reads force field parameters per GA evaluation pass, per MPI 
                  rank
****************************************************************************/

int Read_Force_Field_Eff (char *ffield_file, eff_interaction * ffdata, double ***ffid_pointer, int rank)
{
  FILE *fi;
  char *s;
  char **tmp;
  int a, token, i, j, n;
  double val;
  int total_ecp_parameters;

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
  /* line 2 is number of ecp parameters */
  fgets (s, MAX_LINE, fi);
  token = tokenize_string (s, &tmp);
  /* reading the number of ecp parameters */
  n = atoi (tmp[0]);

  if (n < 1) {
    fprintf (stderr, "WARNING: [%d] number of ECP parameters in ffield file is 0!\n",rank);
    exit (1);
  }

  ffdata->num_atom_types = n;
//  ffdata->gp.n = n;
  total_ecp_parameters = 0;	// ecpradius, p1-p3 for s type and p1-p5 for p type
  ffdata->ecp = (ecp_parameters *) scalloc (ffdata->num_atom_types, sizeof (ecp_parameters), "ecp");
  *ffid_pointer = (double **) scalloc (ffdata->num_atom_types * 8, sizeof (double *), "ffid_pointer");

  a = 0;
  for (i = 0; i < ffdata->num_atom_types; i++) {
    fgets (s, MAX_LINE, fi);
    token = tokenize_string (s, &tmp);
//    ffdata->ecp[i].atype = atoi (tmp[0]);
    for (j = 0; j < (int) strlen (tmp[0]); ++j)
      ffdata->ecp[i].name[j] = tmp[0][j];
    for (j = 0; j < (int) strlen (tmp[1]); ++j)
      ffdata->ecp[i].ecptype[j] = tmp[1][j];
    if (strcmp (ffdata->ecp[i].ecptype, "s") == 0) {
      val = (double) atof (tmp[2]);
      ffdata->ecp[i].ecpradius = val;
      (*ffid_pointer)[a++] = &ffdata->ecp[i].ecpradius; 
      val = (double) atof (tmp[3]);
      ffdata->ecp[i].p1 = val;
      (*ffid_pointer)[a++] = &ffdata->ecp[i].p1;
      val = (double) atof (tmp[4]);
      ffdata->ecp[i].p2 = val;
      (*ffid_pointer)[a++] = &ffdata->ecp[i].p2;
      val = (double) atof (tmp[5]);
      ffdata->ecp[i].p3 = val;
      (*ffid_pointer)[a++] = &ffdata->ecp[i].p3;
      val = 0.0;
      ffdata->ecp[i].p4 = val;
      (*ffid_pointer)[a++] = &ffdata->ecp[i].p4;
      ffdata->ecp[i].p5 = val;
      (*ffid_pointer)[a++] = &ffdata->ecp[i].p5;
      ffdata->ecp[i].p6 = val;
      (*ffid_pointer)[a++] = &ffdata->ecp[i].p6;
      ffdata->ecp[i].p7 = val;
      (*ffid_pointer)[a++] = &ffdata->ecp[i].p7;
//      ffdata->ecp[i].p8 = val;
//      (*ffid_pointer)[a++] = &ffdata->ecp[i].p8;
      total_ecp_parameters += 4;
    }
    else if ((strcmp (ffdata->ecp[i].ecptype, "p") == 0) || (strcmp (ffdata->ecp[i].ecptype, "h") == 0)) {
      val = (double) atof (tmp[2]);
      ffdata->ecp[i].ecpradius = val;
      (*ffid_pointer)[a++] = &ffdata->ecp[i].ecpradius;
      val = (double) atof (tmp[3]);
      ffdata->ecp[i].p1 = val;
      (*ffid_pointer)[a++] = &ffdata->ecp[i].p1;
      val = (double) atof (tmp[4]);
      ffdata->ecp[i].p2 = val;
      (*ffid_pointer)[a++] = &ffdata->ecp[i].p2;
      val = (double) atof (tmp[5]);
      ffdata->ecp[i].p3 = val;
      (*ffid_pointer)[a++] = &ffdata->ecp[i].p3;
      val = (double) atof (tmp[6]);
      ffdata->ecp[i].p4 = val;
      (*ffid_pointer)[a++] = &ffdata->ecp[i].p4;
      val = (double) atof (tmp[7]);
      ffdata->ecp[i].p5 = val;
      (*ffid_pointer)[a++] = &ffdata->ecp[i].p5;
      val = 0.0;
      ffdata->ecp[i].p6 = val;
      (*ffid_pointer)[a++] = &ffdata->ecp[i].p6;
      ffdata->ecp[i].p7 = val;
      (*ffid_pointer)[a++] = &ffdata->ecp[i].p7;
//      ffdata->ecp[i].p8 = val;
//      (*ffid_pointer)[a++] = &ffdata->ecp[i].p8;
      total_ecp_parameters += 6;
    }
    else if (strcmp (ffdata->ecp[i].ecptype, "x") == 0) {
      val = (double) atof (tmp[2]);
      ffdata->ecp[i].p1 = val;
      (*ffid_pointer)[a++] = &ffdata->ecp[i].p1;
      val = (double) atof (tmp[3]);
      ffdata->ecp[i].p2 = val;
      (*ffid_pointer)[a++] = &ffdata->ecp[i].p2;
      val = (double) atof (tmp[4]);
      ffdata->ecp[i].p3 = val;
      (*ffid_pointer)[a++] = &ffdata->ecp[i].p3;
      val = (double) atof (tmp[5]);
      ffdata->ecp[i].p4 = val;
      (*ffid_pointer)[a++] = &ffdata->ecp[i].p4;
      val = (double) atof (tmp[6]);
      ffdata->ecp[i].p5 = val;
      (*ffid_pointer)[a++] = &ffdata->ecp[i].p5;
      val = (double) atof (tmp[7]);
      ffdata->ecp[i].p6 = val;
      (*ffid_pointer)[a++] = &ffdata->ecp[i].p6;
      val = 0.0;
      ffdata->ecp[i].p7 = val;
      (*ffid_pointer)[a++] = &ffdata->ecp[i].p7;
      ffdata->ecp[i].p8 = val;
      (*ffid_pointer)[a++] = &ffdata->ecp[i].p8;
      total_ecp_parameters += 6;
    }
    else if (strcmp (ffdata->ecp[i].ecptype, "f") == 0) {
      val = (double) atof (tmp[2]);
      ffdata->ecp[i].p1 = val;
      (*ffid_pointer)[a++] = &ffdata->ecp[i].p1;
      val = (double) atof (tmp[3]);
      ffdata->ecp[i].p2 = val;
      (*ffid_pointer)[a++] = &ffdata->ecp[i].p2;
      val = (double) atof (tmp[4]);
      ffdata->ecp[i].p3 = val;
      (*ffid_pointer)[a++] = &ffdata->ecp[i].p3;
      val = (double) atof (tmp[5]);
      ffdata->ecp[i].p4 = val;
      (*ffid_pointer)[a++] = &ffdata->ecp[i].p4;
      val = (double) atof (tmp[6]);
      ffdata->ecp[i].p5 = val;
      (*ffid_pointer)[a++] = &ffdata->ecp[i].p5;
      val = (double) atof (tmp[7]);
      ffdata->ecp[i].p6 = val;
      (*ffid_pointer)[a++] = &ffdata->ecp[i].p6;
      val = (double) atof (tmp[8]);
      ffdata->ecp[i].p7 = val;
      (*ffid_pointer)[a++] = &ffdata->ecp[i].p7;
      val = (double) atof (tmp[9]);
      ffdata->ecp[i].p8 = val;
      (*ffid_pointer)[a++] = &ffdata->ecp[i].p8;

      total_ecp_parameters += 8;
    }
    else if (strcmp (ffdata->ecp[i].ecptype, "L") == 0) {
      val = (double) atof (tmp[2]);
      ffdata->ecp[i].p1 = val;
      (*ffid_pointer)[a++] = &ffdata->ecp[i].p1;
      val = (double) atof (tmp[3]);
      ffdata->ecp[i].p2 = val;
      (*ffid_pointer)[a++] = &ffdata->ecp[i].p2;

      val = (double) atof (tmp[4]);
      ffdata->ecp[i].p3 = val;
      (*ffid_pointer)[a++] = &ffdata->ecp[i].p3;
      val = (double) atof (tmp[5]);
      ffdata->ecp[i].p4 = val;
      (*ffid_pointer)[a++] = &ffdata->ecp[i].p4;
      val = (double) atof (tmp[6]);
      ffdata->ecp[i].p5 = val;
      (*ffid_pointer)[a++] = &ffdata->ecp[i].p5;
      val = (double) atof (tmp[7]);
      ffdata->ecp[i].p6 = val;
      (*ffid_pointer)[a++] = &ffdata->ecp[i].p6;
      val = (double) atof (tmp[8]);
      ffdata->ecp[i].p7 = val;
      (*ffid_pointer)[a++] = &ffdata->ecp[i].p7;
      val = (double) atof (tmp[9]);
      ffdata->ecp[i].p8 = val;
      (*ffid_pointer)[a++] = &ffdata->ecp[i].p8;

      total_ecp_parameters += 8;
    }
  }
  ffdata->gp.n = a;

#if DEBUG
  if (rank == 0 && debug_level == 1) {
    printf("eFF-ECP force field\n");
    for (i=0; i < a; i++) {
      printf("%f ", *(*ffid_pointer)[i]);
//  if (i%6!=0) printf("\n");
    }
    printf ("\n");
  }
#endif

  if (rank == 0 && debug_level == 1) DEBUG_PRINT (("FFIELD: [%d] Read %d (non-zero) eFF-ECP parameters\n", rank,total_ecp_parameters));

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
