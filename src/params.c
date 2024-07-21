/* ----------------------------------------------------------------------
    GARFFIELD: Genetic-Algorithm Reactive Force FIELD parameter optimizer.

    Copyright (2012-2014) California Institute of Technology
    Andres Jaramillo-Botero (ajaramil@caltech.edu)
    http://www.wag.caltech.edu/home/ajaramil/GARFfield.html
*------------------------------------------------------------------------- */

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "params.h"
#include "mpi.h"

#define MAX_LINE 512
#define MAX_TOKEN_LEN 512

enum {REAX, REAXC, EFF, PQEQ, PPR, CG, MORSE, COMB, TERSOFF, TERSOFF_MOD, ZHOU_EAM};

#ifdef DEBUG
# define DEBUG_PRINT(x) printf x
#else
# define DEBUG_PRINT(x)
#endif

/****************************************************************************
*  Function     : parse_params
*  Description  : Parses user parameter selection file, depending on the
                  force field defined by the user (with -F) and populates
                  the ffield parameters pointer (params_ptr)
*  Parameters   : parameter file, user forcefield selection, pointer to
                  ffield vector (parameters to be modified), and 
                  process rank
*  Effects      : Populates global array of parameter indices in the global
                  forcefield array (ffid_pointer) 
****************************************************************************/

void parse_params (char *params_file, int forcefield2optimize, int *params_ptr, int rank)
{

  char **strtmp;
  FILE *fp;
  int i, c, section, entry, parameter;
  char *line;
  char section_label[MAX_LINE];
  char filen[MAX_LINE];

  sprintf (filen, "%s", params_file);
  fp = fopen (filen, "r");
  if (fp == NULL) {
    if (rank == 0)
      printf ("ERROR: Cannot open params file %s!\n", filen);
    MPI_Abort (MPI_COMM_WORLD, 1);
  }

  line = (char *) scalloc (MAX_LINE, sizeof (char), "line");
  strtmp = (char **) scalloc (6, sizeof (char *), "strtmp");
  for (i = 0; i < 6; i++)
    strtmp[i] = (char *) scalloc (MAX_TOKEN_LEN, sizeof (char), "tmpi");

  i = 0;
  while (fgets (line, MAX_LINE, fp) != NULL) {
    c = tokenize_string (line, &strtmp);
    if (c != 0 && c > 5) {
      section = atoi (strtmp[0]);
      entry = atoi (strtmp[1]) - 1;
      parameter = atoi (strtmp[2]) - 1;

      switch (section) {
      case 1:			// general parameters
	if ((forcefield2optimize == REAX) || (forcefield2optimize == REAXC)) {
	  if (parameter <= ffdata->gp.n_global) {
	    params_ptr[i] = parameter;
	    strcpy (section_label, "GENERAL");
	  }
	  else {
	    if (rank == 0)
	      printf ("ERROR: %d %d does not exist !\n", section, parameter + 1);
	    MPI_Abort (MPI_COMM_WORLD, 1);
	  }
	}
	else if (forcefield2optimize == EFF) {
	  if (parameter <= effdata->gp.n) {
	    params_ptr[i] = parameter + entry * 8;
	    strcpy (section_label, "GENERAL");
	  }
	  else {
	    if (rank == 0)
	      printf ("ERROR: %d %d does not exist !\n", section, parameter + 1);
	    MPI_Abort (MPI_COMM_WORLD, 1);
	  }
	}
        else if (forcefield2optimize == PQEQ) {
          if (parameter <= pqeqffdata->total_pqeq_parameters) {
            params_ptr[i] = parameter + entry * 7;
            strcpy (section_label, "GENERAL");
          }
          else {
            if (rank == 0)
              printf ("ERROR: %d %d does not exist !\n", section, parameter + 1);
            MPI_Abort (MPI_COMM_WORLD, 1);
          }
        }
        else if (forcefield2optimize == PPR) {
          if (parameter <= pprffdata->total_ppr_parameters) {
            params_ptr[i] = parameter + entry * 7;
            strcpy (section_label, "GENERAL");
          }
          else {
            if (rank == 0)
              printf ("ERROR: %d %d does not exist !\n", section, parameter + 1);
            MPI_Abort (MPI_COMM_WORLD, 1);
          }
        }
	else if (forcefield2optimize == CG) {
	  if (parameter <= cgffdata->total_cg_parameters) {
	    params_ptr[i] = parameter + entry * 3;
	    strcpy (section_label, "GENERAL");
	  }
	  else {
	    if (rank == 0)
	      printf ("ERROR: %d %d does not exist !\n", section, parameter + 1);
	    MPI_Abort (MPI_COMM_WORLD, 1);
	  }
	}
	else if (forcefield2optimize == MORSE) {
	  if (parameter <= morseffdata->total_morse_parameters) {
	    params_ptr[i] = parameter + entry * 3;
	    strcpy (section_label, "GENERAL");
	  }
	  else {
	    if (rank == 0)
	      printf ("ERROR: %d %d does not exist !\n", section, parameter + 1);
	    MPI_Abort (MPI_COMM_WORLD, 1);
	  }
	}
        else if (forcefield2optimize == COMB) {
          if (parameter <= combffdata->num_entries*46) {
            params_ptr[i] = parameter + entry * 46;
            strcpy (section_label, "GENERAL");
          }
          else {
            if (rank == 0)
              printf ("ERROR: %d %d does not exist !\n", section, parameter + 1);
            MPI_Abort (MPI_COMM_WORLD, 1);
          }
        }    
        else if (forcefield2optimize == TERSOFF) {
          if (parameter <= tersoffffdata->num_entries*14) {
            params_ptr[i] = parameter + entry * 14;
            strcpy (section_label, "GENERAL");
          }
          else {
            if (rank == 0)
              printf ("ERROR: %d %d does not exist !\n", section, parameter + 1);
            MPI_Abort (MPI_COMM_WORLD, 1);
          }
        }
        else if (forcefield2optimize == TERSOFF_MOD) {
          if (parameter <= tersoff_modffdata->num_entries*17) {
            params_ptr[i] = parameter + entry * 17;
            strcpy (section_label, "GENERAL");
          }
          else {
            if (rank == 0)
              printf ("ERROR: %d %d does not exist !\n", section, parameter + 1);
            MPI_Abort (MPI_COMM_WORLD, 1);
          }
        }
      else if (forcefield2optimize == ZHOU_EAM) {
          if (parameter <= zhou_EAMffdata->num_entries*20) {
            params_ptr[i] = parameter + entry * 20;
            strcpy (section_label, "GENERAL");
          }
          else {
            if (rank == 0)
              printf ("ERROR: %d %d does not exist !\n", section, parameter + 1);
            MPI_Abort (MPI_COMM_WORLD, 1);
          }
        }
	break;
      case 2:			// single body parameters
	if (entry <= ffdata->num_atom_types) {
	  if (lg_reax) params_ptr[i] = gp_idx + entry * 34 + parameter;
          else params_ptr[i] = gp_idx + entry * 32 + parameter;
	  strcpy (section_label, "ATOMS");
	}
	else {
	  if (rank == 0)
	    printf ("ERROR: %d %d %d does not exist !\n", section, entry + 1, parameter + 1);
	  MPI_Abort (MPI_COMM_WORLD, 1);
	}
	break;
      case 3:			// two body parameters
	if (entry <= ffdata->num_bond_types) {
	  params_ptr[i] = gp_idx + sbp_idx + entry * 16 + parameter;
	  strcpy (section_label, "BONDS");
	}
	else {
	  if (rank == 0)
	    printf ("ERROR: %d %d %d does not exist !\n", section, entry + 1, parameter + 1);
	  MPI_Abort (MPI_COMM_WORLD, 1);
	}
	break;
      case 4:			// two body off diagonal parameters
	if (entry <= ffdata->num_off_diag_types) {
	  if (lg_reax) params_ptr[i] = gp_idx + sbp_idx + tbp_idx + entry * 7 + parameter;
	  else params_ptr[i] = gp_idx + sbp_idx + tbp_idx + entry * 6 + parameter;
          strcpy (section_label, "OFF-DIAGONALS");
	}
	else {
	  if (rank == 0)
	    printf ("ERROR: %d %d %d does not exist !\n", section, entry + 1, parameter + 1);
	  MPI_Abort (MPI_COMM_WORLD, 1);
	}
	break;
      case 5:			// three body parameters
	if (entry <= ffdata->num_angle_types) {
	  params_ptr[i] = gp_idx + sbp_idx + tbp_idx + tbodp_idx + entry * 7 + parameter;
	  strcpy (section_label, "ANGLES");
	}
	else {
	  if (rank == 0)
	    printf ("ERROR: %d %d %d does not exist !\n", section, entry + 1, parameter + 1);
	  MPI_Abort (MPI_COMM_WORLD, 1);
	}
	break;
      case 6:			// four body parameters
	if (entry <= ffdata->num_torsion_types) {
	  params_ptr[i] = gp_idx + sbp_idx + tbp_idx + tbodp_idx + thbp_idx + entry * 7 + parameter;
	  strcpy (section_label, "TORSIONS");
	}
	else {
	  if (rank == 0)
	    printf ("ERROR: %d %d %d does not exist !\n", section, entry + 1, parameter + 1);
	  MPI_Abort (MPI_COMM_WORLD, 1);
	}
	break;
      case 7:			// hydrogen body parameters
	if (entry <= ffdata->num_hbond_types) {
	  params_ptr[i] =
	    gp_idx + sbp_idx + tbp_idx + tbodp_idx + thbp_idx + fbp_idx + entry * 4 + parameter;
	  strcpy (section_label, "H-BONDS");
	}
	else {
	  if (rank == 0)
	    printf ("ERROR: %d %d %d does not exist !\n", section, entry + 1, parameter + 1);
	  MPI_Abort (MPI_COMM_WORLD, 1);
	}
	break;
      }
      if (section == 1) {
	if (rank == 0 && debug_level == 1)
	  DEBUG_PRINT (("--> %s: change parameter %d [with global index %d] within [%f:%f]\n",
			section_label, parameter + 1, params_ptr[i] + 1, Lower[i], Upper[i]));
      }
      else {
	if (rank == 0 && debug_level == 1)
	  DEBUG_PRINT (("--> %s: change entry %d, parameter %d [with global index %d] within [%f:%f]\n", section_label, entry + 1, parameter + 1, params_ptr[i] + 1, Lower[i], Upper[i]));
      }
      i++;

    }
    else if (rank == 0)
      printf ("ERROR: Check params file format\n");
  }

  for (i = 0; i < 6; i++)
    free (strtmp[i]);
  free (strtmp);
  free (line);
  fclose (fp);
}

/****************************************************************************
*  Function     : get_num_params
*  Description  : Returns the total number of parameters to be trained
*  Parameters   : User-defined parameter file and MPI process rank 
*  Effects      : Returns total number of parameters to be trained
****************************************************************************/

int get_num_params (char *params_file, int rank)
{
  char **strtmp;
  FILE *fp;
  int i, j, c;
  char *line;
  char filen[MAX_LINE];

  sprintf (filen, "%s", params_file);
  fp = fopen (filen, "r");
  if (fp == NULL) {
    if (rank == 0)
      printf ("ERROR: Cannot open params file %s!\n", filen);
    MPI_Abort (MPI_COMM_WORLD, 1);
  }

  line = (char *) scalloc (MAX_LINE, sizeof (char), "line");
  strtmp = (char **) scalloc (MAX_LINE, sizeof (char *), "strtmp");
  for (i = 0; i < 6; i++)
    strtmp[i] = (char *) scalloc (MAX_TOKEN_LEN, sizeof (char), "tmpi");

  // Count number of parameters to train from params file
  j = 0;
  while (fgets (line, MAX_LINE, fp) != NULL) {
    c = tokenize_string (line, &strtmp);
    if (c > 0 && c > 5)
      j++;
  }

  for (i = 0; i < 6; i++)
    free (strtmp[i]);
  free (strtmp);
  free (line);

  fclose (fp);
  return j;
}

/****************************************************************************
*  Function     : get_params
*  Description  : Populates the Lower and Upper arrays with the parameters
                  1) lower and upper ranges or 2) the mean and percent offset (if
                  the user selects to start with the existing ffield values)
*  Parameters   : User-defined parameters file, flag that defines if user has
                  selected 1) or 2)
*  Effects      : Populates Lower and Upper global arrays
****************************************************************************/

void get_params (char *params_file, int flag, int rank)
{
  char **strtmp;
  FILE *fp;
  int i, j, c;
  char filen[MAX_LINE];
  char *line;

  sprintf (filen, "%s", params_file);
  fp = fopen (filen, "r");
  if (fp == NULL) {
    if (rank == 0)
      printf ("ERROR: Cannot open params file %s!\n", filen);
    MPI_Abort (MPI_COMM_WORLD, 1);
  }

  line = (char *) scalloc (MAX_LINE, sizeof (char), "line");
  strtmp = (char **) scalloc (MAX_LINE, sizeof (char *), "strtmp");
  for (i = 0; i < 6; i++)
    strtmp[i] = (char *) scalloc (MAX_TOKEN_LEN, sizeof (char), "tmpi");

  // Parse params file
  i = 0;
  while (fgets (line, MAX_LINE, fp) != NULL) {
    c = tokenize_string (line, &strtmp);
    if (c > 0 && c > 5) {
      Lower[i] = atof (strtmp[4]);
      Upper[i] = atof (strtmp[5]);
      if (!flag) {
	if ((Upper[i] > 1.0) || (Upper[i] <= 0.0)) {
	  if (rank == 0)
	    printf (">> Incorrect percent for parameter %d (0-1] in params file\n", i);
	  MPI_Abort (MPI_COMM_WORLD, 1);
	}
      }
      i++;
    }
    else {
      if (rank == 0)
	printf ("ERROR: Check format of params file");
      MPI_Abort (MPI_COMM_WORLD, 1);
    }
  }

  // Free local memory
  for (j = 0; j < 6; j++)
    free (strtmp[j]);
  free (strtmp);
  free (line);
  fclose (fp);
}
