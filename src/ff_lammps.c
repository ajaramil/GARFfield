/* ----------------------------------------------------------------------
   GARFFIELD: Genetic-Algorithm Reactive Force FIELD parameter optimizer.

   Copyright (2012-2014) California Institute of Technology
   Andres Jaramillo-Botero (ajaramil@caltech.edu)
   http://www.wag.caltech.edu/home/ajaramil/GARFfield.html
------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "library.h"		/* this is the LAMMPS library include file */
#include "mpi.h"

#include "ff_lammps.h"
#include "restraints.h"
#include "reaxc_types.h"
#include "eff_types.h"
#include "cg_types.h"
#include "trainset.h"
#include "structures.h"

extern double calcerror (double *, int);

#ifdef DEBUG
# define DEBUG_PRINT(x) printf x
#else
# define DEBUG_PRINT(x)
#endif

/****************************************************************************
*  Function     : get_lammps_eng
*  Description  : Performs the evaluation of the fitness function by calling
                  an instance of LAMMPS per chromosome in the GA population.
                  This is executed concurrently for all user-defined training
                  cases and distributed among the number of MPI processes
                  available for computation.  The fitness function is currently
                  a weighted sum of all training cases (i.e. single objective),
                  including charges, geometries, energies, etc.  The control 
                  flow over the kind of LAMMPS calculations performed is 
                  determined by the force field type (i.e. -F).  Tallies all
                  values into the array energies.
*  Parameters   : Force field file name, force field type, MPI process rank,
                  and MPI communicator 
*  Effects      : Allocates and deallocates the LAMMPS environment for every
                  evaluation pass of user-defined training set
****************************************************************************/

double get_lammps_eng (char *forcefield, int forcefield2optimize, int rank, MPI_Comm comm)
{
  int i, j, k, n, id, id2, natoms, first, tokens;
  void *ptr = NULL;
  double *ff_energies, *energies, *cell_stress, *cell_press, *rms_force, rs1, rs2;
  double *xlo, *xhi, *ylo, *yhi, *zlo, *zhi, *xy, *xz, *yz;
  double *chrg, **atoms, *lmp_atoms, val, *force, *apos;
  int geom_flag, cell_flag, charge_flag, force_flag, stress_flag, energy_flag, atom_force_flag;
  double k1, k2, cutoff;
  char st1[40] = "cp ";
  char *cmd, *line;
  char **strtmp;
  double total_error, total_rms_force;
  FILE *fi;
  double *ff_tmp;

  if (forcefield2optimize == REAX) {	// use Fortran version of reaxFF
    strcat (st1, forcefield);
    strcat (st1, " ffield.reax");
    cmd = st1;			// cp ffield.new into ffield.reax
  }

  int XX = 1 << 6;
  int YY = 1 << 5;
  int ZZ = 1 << 4;
  int XY = 1 << 3;
  int XZ = 1 << 2;
  int YZ = 1 << 1;

  /* Lammps initialization */
  if (debug_level == 1 || debug_level == 2) 
    lammps_open (0, NULL, comm, &ptr);
  else 
    lammps_open (5, cmd_line, comm, &ptr);
// lammps_open_no_mpi (5,cmd_line,&ptr);

  energies = (double *) scalloc (num, sizeof (double), "lenergies");
//  ff_energies = (double *) scalloc (1, sizeof (double), "ff_energies");
  line = (char *) scalloc (MAX_LINE, sizeof (char), "lines");

  total_rms_force = 0;

  for (i = 0; i < num; i++) {
    if (rank == 0 && debug_level == 1)
      DEBUG_PRINT (("--> LAMMPS: Structure [%i/%i] from datafile %s\n", i + 1, num, data_files[i]));

//    lammps_open (5, cmd_line, comm, &ptr);
    geom_flag = charge_flag = cell_flag = force_flag = stress_flag = energy_flag = atom_force_flag = 0;

    if ((forcefield2optimize == EFF) && (eff_units == 0))
      strcpy (line, "units           electron");
    else if (forcefield2optimize == COMB)
      strcpy (line, "units           metal");
    else
      strcpy (line, "units           real");
    if (rank == 0 && debug_level == 1)
      DEBUG_PRINT (("%s\n", line));
    lammps_command (ptr, line);
    if (forcefield2optimize == COMB)
      strcpy (line, "newton on");
    if (rank == 0 && debug_level == 1)
      DEBUG_PRINT (("%s\n", line));
    lammps_command (ptr, line);

    if (cellflag[i])
      strcpy (line, "boundary p p p");
    else
      strcpy (line, "boundary f f f");
    if (rank == 0 && debug_level == 1)
      DEBUG_PRINT (("%s\n", line));
    lammps_command (ptr, line);

    if ((forcefield2optimize == REAX) || (forcefield2optimize == REAXC)
	|| (forcefield2optimize == COMB) || (forcefield2optimize == PQEQ) || (forcefield2optimize == PPR)) {
      if (pqeq_flag)
	strcpy (line, "atom_style	pqeq");
      else
	strcpy (line, "atom_style	charge");
    }
    else if (forcefield2optimize == EFF)
      strcpy (line, "atom_style      electron");
    else if (forcefield2optimize == CG) {
      strcpy (line, "atom_style     full");
      if (rank == 0 && debug_level == 1)
	DEBUG_PRINT (("%s\n", line));
      lammps_command (ptr, line);
      strcpy (line, "special_bonds     lj 0 1 1");
    }
    else if (forcefield2optimize == MORSE)
      strcpy (line, "atom_style     atomic");
    if (rank == 0 && debug_level == 1)
      DEBUG_PRINT (("%s\n", line));
    lammps_command (ptr, line);

    if (forcefield2optimize == CG) {
      strcpy (line, "bond_style      harmonic");
      if (rank == 0 && debug_level == 1)
	DEBUG_PRINT (("%s\n", line));
      lammps_command (ptr, line);
    }
    // apply only once if geometry has restrains
/*    for (n = 0; n < nrst; n++) {
      if (strcmp (fname[i],rstrain[n].name) == 0) {
        sprintf (line, "atom_modify map array");
        lammps_command (ptr, line); 
        break;
      }
    }
*/
    strcpy (line, "atom_modify map array");
    if (rank == 0 && debug_level == 1)
      DEBUG_PRINT (("%s\n", line));
    lammps_command (ptr, line);

    sprintf (line, "box tilt large");	// just in case we have a triclinic cell with large tilts
    if (rank == 0 && debug_level == 1)
      DEBUG_PRINT (("%s\n", line));
    lammps_command (ptr, line);

    if ((forcefield2optimize == REAX) || (forcefield2optimize == REAXC)
	|| (forcefield2optimize == MORSE) || (forcefield2optimize == COMB)
	|| (forcefield2optimize == TERSOFF) || (forcefield2optimize == TERSOFF_MOD)
	|| (forcefield2optimize == ZHOU_EAM) || (forcefield2optimize == PQEQ) || (forcefield2optimize == PPR)) {
      sprintf (line, "read_data data.%s.%d", fname[i], rank);
      if (rank == 0 && debug_level == 1)
	DEBUG_PRINT (("%s\n", line));
      lammps_command (ptr, line);
    }
    else if (forcefield2optimize == EFF || forcefield2optimize == CG) {
      sprintf (line, "read_data data.%s", fname[i]);
      if (rank == 0 && debug_level == 1)
	DEBUG_PRINT (("%s\n", line));
      lammps_command (ptr, line);

      double boxxlo = *((double *) lammps_extract_global (ptr, "boxxlo"));
      double boxxhi = *((double *) lammps_extract_global (ptr, "boxxhi"));
      double boxylo = *((double *) lammps_extract_global (ptr, "boxylo"));
      double boxyhi = *((double *) lammps_extract_global (ptr, "boxyhi"));
      double boxzlo = *((double *) lammps_extract_global (ptr, "boxzlo"));
      double boxzhi = *((double *) lammps_extract_global (ptr, "boxzhi"));

      double xprd = (boxxhi - boxxlo);
      double yprd = (boxyhi - boxylo);
      double zprd = (boxzhi - boxzlo);

      cutoff = xprd;
      if (yprd > cutoff)
	cutoff = yprd;
      if (zprd > cutoff)
	cutoff = zprd;
      // eFF uses minimum image convention with cutoff=1/2 shorted box length
      cutoff /= 2;		// 1/2 box max length
    }

    if (forcefield2optimize == REAXC) {
      if (pqeq_flag) {
	double boxxlo = *((double *) lammps_extract_global (ptr, "boxxlo"));
	double boxxhi = *((double *) lammps_extract_global (ptr, "boxxhi"));
	double boxylo = *((double *) lammps_extract_global (ptr, "boxylo"));
	double boxyhi = *((double *) lammps_extract_global (ptr, "boxyhi"));
	double boxzlo = *((double *) lammps_extract_global (ptr, "boxzlo"));
	double boxzhi = *((double *) lammps_extract_global (ptr, "boxzhi"));

	double xprd = (boxxhi - boxxlo);
	double yprd = (boxyhi - boxylo);
	double zprd = (boxzhi - boxzlo);

	cutoff = xprd;
	if (yprd > cutoff)
	  cutoff = yprd;
	if (zprd > cutoff)
	  cutoff = zprd;
	cutoff /= 2;
	// Use a fixed cutoff, instead of minimage
	sprintf (line,
//               "pair_style hybrid/overlay coul/pqeqgauss 0.0 10.0 reax/c NULL coulomb_off yes checkqeq no");
		 "pair_style hybrid/overlay coul/pqeqgauss 0 12.5 reax/c NULL lgvdw yes coulomb_off yes checkqeq no");
      }
      else
	strcpy (line, "pair_style reax/c NULL");
      if (lg_reax)
	strcat (line, " lgvdw yes");
      if (safezone_flag) {
	sprintf (st1, " safezone %2.2f", safezone_value);
	strcat (line, st1);
      }
      if (mincap_flag) {
	sprintf (st1, " mincap %d", mincap_value);
	strcat (line, st1);
      }
    }
    if (forcefield2optimize == PPR) {
      if (pqeq_flag) {
        double boxxlo = *((double *) lammps_extract_global (ptr, "boxxlo"));
        double boxxhi = *((double *) lammps_extract_global (ptr, "boxxhi"));
        double boxylo = *((double *) lammps_extract_global (ptr, "boxylo"));
        double boxyhi = *((double *) lammps_extract_global (ptr, "boxyhi"));
        double boxzlo = *((double *) lammps_extract_global (ptr, "boxzlo"));
        double boxzhi = *((double *) lammps_extract_global (ptr, "boxzhi"));

        double xprd = (boxxhi - boxxlo);
        double yprd = (boxyhi - boxylo);
        double zprd = (boxzhi - boxzlo);

        cutoff = xprd;
        if (yprd > cutoff)
          cutoff = yprd;
        if (zprd > cutoff)
          cutoff = zprd;
        cutoff /= 2;
        // Use a fixed cutoff, instead of minimage
        sprintf (line,
                 "pair_style      hybrid/overlay coul/pqeqgauss 0.00 12.50 ppr 0.0 12.5");
      }
      else
        strcpy (line, "pair_style ppr 0.0 12.5");
    }
    else if (forcefield2optimize == REAX)
      strcpy (line, "pair_style      reax 10.0 0 1 1.0e-5");
    else if (forcefield2optimize == EFF)
      sprintf (line, "pair_style      eff/cut %8.4f ecp %s", cutoff, atypes[i]);
    else if (forcefield2optimize == CG)
      strcpy (line, "pair_style      morse 10");
    else if (forcefield2optimize == MORSE)
      strcpy (line, "pair_style      morse 2.5");
    else if (forcefield2optimize == COMB)
      strcpy (line, "pair_style      comb");
    else if (forcefield2optimize == TERSOFF)
      strcpy (line, "pair_style      tersoff");
    else if (forcefield2optimize == TERSOFF_MOD)
      strcpy (line, "pair_style      tersoff/mod");
    else if (forcefield2optimize == ZHOU_EAM)
      strcpy (line, "pair_style      eam/alloy");

    if (rank == 0 && debug_level == 1)
      DEBUG_PRINT (("%s\n", line));
    lammps_command (ptr, line);

    // allocate helper storage
    strtmp = (char **) scalloc (MAX_TOKENS, sizeof (char *), "tmp3");
    for (k = 0; k < MAX_TOKENS; k++)
      strtmp[k] = (char *) scalloc (MAX_TOKEN_LEN, sizeof (char), "tmp3i");

    if (forcefield2optimize == REAXC) {
      if (pqeq_flag) {
	char pqeq_id[5];
	int polarizability;
	float chi, idem, Rcore, Qcore, Rshell, k2, k4;
	char qeqfile[25];

	sprintf (line, "pair_coeff * * reax/c %s %s", forcefield, atypes[i]);
	if (rank == 0 && debug_level == 1)
	  DEBUG_PRINT (("%s\n", line));
	lammps_command (ptr, line);

        sprintf (line, "pair_coeff * * coul/pqeqgauss 0 0 0 0 0 0 0 0");
        if (rank == 0 && debug_level == 1)
          DEBUG_PRINT (("%s\n", line));
        lammps_command (ptr, line);

        tokens = tokenize_string (atypes[i], &strtmp);
        /* open pqeq par file */
	sprintf (qeqfile, "%s.%d", pqeq_par_file, rank);
	if ((fi = fopen (qeqfile, "r")) == NULL) {
	  fprintf (stderr, "ERROR: opening %s! terminating...\n", qeqfile);
	  MPI_Abort (MPI_COMM_WORLD, 1);
	}
	while (fscanf (fi, "%s %d %f %f %f %f %f %f %f", pqeq_id, &polarizability, &chi, &idem, &Qcore, &Rcore, &Rshell, &k2, &k4) == 9) {
	  for (j = 0; j < tokens; j++) {
	    if (strcmp (pqeq_id, strtmp[j]) == 0) {
	      // pair_coeff type_i type_i  coul/pqeqgauss   Xo       Jo        Rc     P   Qc       Rc       K2       K4  
	      // pair_coeff    1        1   coul/pqeqgauss 8.74100  13.36400  0.66900  1 1.00000  0.66900 415.08000  0.00000  # O
//            sprintf (line, "pair_coeff %d %d coul/pqeqgauss %f %f %f %d %f %f %f %f", j + 1,
//                     j + 1, chi * 23.061, idem * 23.061, Rcore, polarizability, Qcore, Rshell, k2, k4);
	      sprintf (line, "pair_coeff %d %d coul/pqeqgauss %f %f %f %d %f %f %f %f", j + 1,
		       j + 1, chi, idem, Rcore, polarizability, Qcore, Rshell, k2, k4);
	      if (rank == 0 && debug_level == 1)
		DEBUG_PRINT (("%s\n", line));
	      lammps_command (ptr, line);
	    }
	  }
	}
	if (feof (fi))
	  fclose (fi);
	else {
	  fprintf (stderr, "ERROR: parsing %s (check parameters in file)! ...\n", qeqfile);
	  MPI_Abort (MPI_COMM_WORLD, 1);
	}
      }
      else {
	sprintf (line, "pair_coeff      * * %s %s", forcefield, atypes[i]);
	if (rank == 0 && debug_level == 1)
	  DEBUG_PRINT (("%s\n", line));
	lammps_command (ptr, line);
      }
    }
    if (forcefield2optimize == PPR) {
      if (pqeq_flag) {
        char pqeq_id[5];
        int polarizability;
        float chi, idem, Rcore, Qcore, Rshell, k2, k4;
        char qeqfile[25];

        sprintf (line, "pair_coeff * * ppr %s %s", forcefield, atypes[i]);
        if (rank == 0 && debug_level == 1)
          DEBUG_PRINT (("%s\n", line));
        lammps_command (ptr, line);

        sprintf (line, "pair_coeff * * coul/pqeqgauss 0 0 0 0 0 0 0 0");
        if (rank == 0 && debug_level == 1)
          DEBUG_PRINT (("%s\n", line));
        lammps_command (ptr, line);

        tokens = tokenize_string (atypes[i], &strtmp);
        /* open pqeq par file */
        sprintf (qeqfile, "%s.%d", pqeq_par_file, rank);
        if ((fi = fopen (qeqfile, "r")) == NULL) {
          fprintf (stderr, "ERROR: opening %s! terminating...\n", qeqfile);
          MPI_Abort (MPI_COMM_WORLD, 1);
        }
        while (fscanf (fi, "%s %d %f %f %f %f %f %f %f", pqeq_id, &polarizability, &chi, &idem, &Qcore, &Rcore, &Rshell, &k2, &k4) == 9) {
          for (j = 0; j < tokens; j++) {
            if (strcmp (pqeq_id, strtmp[j]) == 0) {
              // pair_coeff type_i type_i  coul/pqeqgauss   Xo       Jo        Rc     P   Qc       Rc       K2       K4
              sprintf (line, "pair_coeff %d %d coul/pqeqgauss %f %f %f %d %f %f %f %f", j + 1,
                       j + 1, chi, idem, Rcore, polarizability, Qcore, Rshell, k2, k4);
              if (rank == 0 && debug_level == 1)
                DEBUG_PRINT (("%s\n", line));
              lammps_command (ptr, line);
            }
          }
        }
        if (feof (fi))
          fclose (fi);
        else {
          fprintf (stderr, "ERROR: parsing %s (check parameters in file)! ...\n", qeqfile);
          MPI_Abort (MPI_COMM_WORLD, 1);
        }
      }
      // Loop over ppr ffield parameters
      int a1, a2;
      tokens = tokenize_string (atypes[i], &strtmp);
      for (j = 0; j < pprffdata->num_pair_types; j++) {
        a1 = a2 = 0;
        for (k = 0; k < tokens; k++)
          if (strcmp (pprffdata->pprptypes[j].type1, strtmp[k]) == 0)
            a1 = k + 1;
        for (k = 0; k < tokens; k++)
          if (strcmp (pprffdata->pprptypes[j].type2, strtmp[k]) == 0)
            a2 = k + 1;
        if ((a1 != 0) && (a2 != 0)) {
          if (a1 < a2) 
            sprintf (line, "pair_coeff    %d %d %f %f %f %f %f %f %f", a1, a2, pprffdata->pprptypes[j].p1, pprffdata->pprptypes[j].p2, pprffdata->pprptypes[j].p3, pprffdata->pprptypes[j].p4, pprffdata->pprptypes[j].p5, pprffdata->pprptypes[j].p6, pprffdata->pprptypes[j].p7);
          else 
            sprintf (line, "pair_coeff    %d %d %f %f %f %f %f %f %f", a2, a1, pprffdata->pprptypes[j].p1, pprffdata->pprptypes[j].p2, pprffdata->pprptypes[j].p3, pprffdata->pprptypes[j].p4, pprffdata->pprptypes[j].p5, pprffdata->pprptypes[j].p6, pprffdata->pprptypes[j].p7);
   
          if (rank == 0 && debug_level == 1)
            DEBUG_PRINT (("%s\n", line));
          lammps_command (ptr, line);
        }
      }
    }
    else if (forcefield2optimize == REAX) {
      system (cmd);		// copy ffield.new into ffield.reax for reax
      sprintf (line, "pair_coeff      * * ffield.reax %s", atypes[i]);
      if (rank == 0 && debug_level == 1)
	DEBUG_PRINT (("%s\n", line));
      lammps_command (ptr, line);
    }
    else if (forcefield2optimize == EFF) {
      strcpy (line, "pair_coeff      * *");
      if (rank == 0 && debug_level == 1)
	DEBUG_PRINT (("%s\n", line));
      lammps_command (ptr, line);

      tokens = tokenize_string (atypes[i], &strtmp);
      for (j = 1; j < tokens; j += 2) {
	for (k = 0; k < effdata->num_atom_types; k++) {
	  if (strcmp (effdata->ecp[k].name, strtmp[j]) == 0) {
	    if (strcmp (effdata->ecp[k].ecptype, "s") == 0) {
	      sprintf (line, "pair_coeff      %s s %f %f %f %f", strtmp[j - 1],
		       effdata->ecp[k].ecpradius, effdata->ecp[k].p1, effdata->ecp[k].p2, effdata->ecp[k].p3);
	    }
	    else if (strcmp (effdata->ecp[k].ecptype, "p") == 0) {
	      sprintf (line, "pair_coeff      %s p %f %f %f %f %f %f", strtmp[j - 1],
		       effdata->ecp[k].ecpradius, effdata->ecp[k].p1, effdata->ecp[k].p2, effdata->ecp[k].p3, effdata->ecp[k].p4, effdata->ecp[k].p5);
	    }
	    else if (strcmp (effdata->ecp[k].ecptype, "h") == 0) {
	      sprintf (line, "pair_coeff      %s h %f %f %f %f %f %f", strtmp[j - 1],
		       effdata->ecp[k].ecpradius, effdata->ecp[k].p1, effdata->ecp[k].p2, effdata->ecp[k].p3, effdata->ecp[k].p4, effdata->ecp[k].p5);
	    }
	    else if (strcmp (effdata->ecp[k].ecptype, "x") == 0) {
	      sprintf (line, "pair_coeff      %s x %f %f %f %f %f %f",
		       strtmp[j - 1], effdata->ecp[k].p1,
		       effdata->ecp[k].p2, effdata->ecp[k].p3, effdata->ecp[k].p4, effdata->ecp[k].p5, effdata->ecp[k].p6);
	    }
	    else if (strcmp (effdata->ecp[k].ecptype, "f") == 0) {
	      sprintf (line, "pair_coeff      %s f %f %f %f %f %f %f %f %f",
		       strtmp[j - 1], effdata->ecp[k].p1,
		       effdata->ecp[k].p2, effdata->ecp[k].p3, effdata->ecp[k].p4,
		       effdata->ecp[k].p5, effdata->ecp[k].p6, effdata->ecp[k].p7, effdata->ecp[k].p8);
	    }
	  }
	}
	if (rank == 0 && debug_level == 1)
	  DEBUG_PRINT (("%s\n", line));
	lammps_command (ptr, line);
      }
      if (cellflag[i]) {
	sprintf (line, "neigh_modify every 2 delay 10 check yes one 10000 page 100000");
	if (rank == 0 && debug_level == 1)
	  DEBUG_PRINT (("%s\n", line));
	lammps_command (ptr, line);
      }
    }
    else if (forcefield2optimize == CG) {
      for (j = 0; j < cgffdata->num_pair_types; j++) {
	sprintf (line, "pair_coeff    %d %d %f %f %f", cgffdata->cgptypes[j].type1,
		 cgffdata->cgptypes[j].type2, cgffdata->cgptypes[j].p1, cgffdata->cgptypes[j].p2, cgffdata->cgptypes[j].p3);
	if (rank == 0 && debug_level == 1)
	  DEBUG_PRINT (("%s\n", line));
	lammps_command (ptr, line);
      }
    }
    else if (forcefield2optimize == MORSE) {
      int a1, a2;
      tokens = tokenize_string (atypes[i], &strtmp);
      for (j = 0; j < morseffdata->num_pair_types; j++) {
	a1 = a2 = 0;
	for (k = 0; k < tokens; k++)
	  if (strcmp (morseffdata->morseptypes[j].type1, strtmp[k]) == 0)
	    a1 = k + 1;
	for (k = 0; k < tokens; k++)
	  if (strcmp (morseffdata->morseptypes[j].type2, strtmp[k]) == 0)
	    a2 = k + 1;
	if (a1 != 0 && a2 != 0) {
	  if (a1 < a2)
	    sprintf (line, "pair_coeff    %d %d %f %f %f 3.0", a1, a2,
		     morseffdata->morseptypes[j].p1, morseffdata->morseptypes[j].p2, morseffdata->morseptypes[j].p3);
	  else
	    sprintf (line, "pair_coeff    %d %d %f %f %f 3.0", a2, a1,
		     morseffdata->morseptypes[j].p1, morseffdata->morseptypes[j].p2, morseffdata->morseptypes[j].p3);
	  if (rank == 0 && debug_level == 1)
	    DEBUG_PRINT (("%s\n", line));
	  lammps_command (ptr, line);
	}
      }
    }
    else if (forcefield2optimize == COMB) {
      sprintf (line, "pair_coeff      * * %s %s", forcefield, atypes[i]);
      if (rank == 0 && debug_level == 1)
	DEBUG_PRINT (("%s\n", line));
      lammps_command (ptr, line);
    }
    else if (forcefield2optimize == TERSOFF) {
      sprintf (line, "pair_coeff      * * %s %s", forcefield, atypes[i]);
      if (rank == 0 && debug_level == 1)
	DEBUG_PRINT (("%s\n", line));
      lammps_command (ptr, line);
    }
    else if (forcefield2optimize == TERSOFF_MOD) {
      sprintf (line, "pair_coeff      * * %s %s", forcefield, atypes[i]);
      if (rank == 0 && debug_level == 1)
	DEBUG_PRINT (("%s\n", line));
      lammps_command (ptr, line);
    }
    else if (forcefield2optimize == ZHOU_EAM) {
      sprintf (line, "pair_coeff      * *  zhou.eam.%d  %s", rank, atypes[i]);	// will this work ?
      if (rank == 0 && debug_level == 1)
	DEBUG_PRINT (("%s\n", line));
      lammps_command (ptr, line);
    }

    /* deallocate helper storage */
    for (k = 0; k < MAX_TOKENS; k++)
      free (strtmp[k]);
    free (strtmp);

    if (forcefield2optimize == REAXC) {
      if (pqeq_flag) {
	// fix pqeq all pqeq 1 0.0 12.5 1.0e-6
	strcpy (line, "fix           pqeq all pqeq 1 0.0 12.5 1.0e-6");
	if (rank == 0 && debug_level == 1)
	  DEBUG_PRINT (("%s\n", line));
	lammps_command (ptr, line);
	strcpy (line, "run 10");
	if (rank == 0 && debug_level == 1)
	  DEBUG_PRINT (("%s\n", line));
	lammps_command (ptr, line);
      }
      else {
	strcpy (line, "fix             1 all qeq/reax 1 0.0 10.0 1.0e-6 reax/c");
	if (rank == 0 && debug_level == 1)
	  DEBUG_PRINT (("%s\n", line));
	lammps_command (ptr, line);
      }
    }

    if ((forcefield2optimize == REAX) || (forcefield2optimize == REAXC))
      sprintf (line, "timestep 0.2");
    else if (forcefield2optimize == EFF)
      sprintf (line, "timestep 0.005");
    else if (forcefield2optimize == CG)
      sprintf (line, "timestep 10.0");
    else if (forcefield2optimize == MORSE)
      sprintf (line, "timestep 1.0");
    else if (forcefield2optimize == COMB)
      sprintf (line, "timestep 0.0002");
    else if (forcefield2optimize == TERSOFF)
      sprintf (line, "timestep 0.001");
    else if (forcefield2optimize == ZHOU_EAM)
      sprintf (line, "timestep 0.001");
    else if (forcefield2optimize == TERSOFF_MOD)
      sprintf (line, "timestep 0.001");
    if (rank == 0 && debug_level == 1)
      DEBUG_PRINT (("%s\n", line));
    lammps_command (ptr, line);

    if (forcefield2optimize == REAXC) {
      if (pqeq_flag) {
	strcpy (line, "compute reax all pair reax/c");
	if (rank == 0 && debug_level == 1)
	  DEBUG_PRINT (("%s\n", line));
	lammps_command (ptr, line);
	strcpy (line, "compute pqeq all pair coul/pqeqgauss ecoul");
      }
      else
	strcpy (line, "compute reax all pair reax/c");
    }
    else if (forcefield2optimize == REAX)
      strcpy (line, "compute reax all pair reax");
    else if (forcefield2optimize == EFF) {
      strcpy (line, "compute eff all pair eff/cut");
      if (rank == 0 && debug_level == 1)
	DEBUG_PRINT (("%s\n", line));
      lammps_command (ptr, line);
      strcpy (line, "communicate     single vel yes");
    }
    else if (forcefield2optimize == CG)
      strcpy (line, "compute cg all pair morse");
    else if (forcefield2optimize == MORSE)
      strcpy (line, "compute emorse all pair morse");
    else if (forcefield2optimize == COMB)
      strcpy (line, "compute ecomb all pair comb");
    else if (forcefield2optimize == TERSOFF)
      strcpy (line, "compute etersoff all pair tersoff");
    else if (forcefield2optimize == ZHOU_EAM)
      strcpy (line, "compute ezhou_EAM all pair eam/alloy");
    else if (forcefield2optimize == TERSOFF_MOD)
      strcpy (line, "compute etersoff_mod all pair tersoff/mod");

    if (rank == 0 && debug_level == 1)
      DEBUG_PRINT (("%s\n", line));
    lammps_command (ptr, line);

    natoms = lammps_get_natoms (ptr);

    //If FORCE case, get forces from single point energy calculation
    for (n = 0; n < nfit[FORCE]; n++) {
      if (strcmp (tset->force[n].sname, fname[i]) != 0)
	continue;
      force_flag = 1;
      strcpy (line, "variable rms_force equal sqrt(fcm(all,x)^2+fcm(all,y)^2+fcm(all,z)^2)");
      if (rank == 0 && debug_level == 1)
	DEBUG_PRINT (("%s\n", line));
      lammps_command (ptr, line);
      strcpy (line, "run 0");
      if (rank == 0 && debug_level == 1)
	DEBUG_PRINT (("%s\n", line));
      lammps_command (ptr, line);
      rms_force = (double *) lammps_extract_variable (ptr, "rms_force", "all");
      tset->force[n].ff_val = *rms_force;
      strcpy (line, "variable rms_force delete");
      if (rank == 0 && debug_level == 1)
	DEBUG_PRINT (("%s\n", line));
      lammps_command (ptr, line);
      // free (rms_force);
      if (rank == 0 && debug_level == 1)
	DEBUG_PRINT (("--> LAMMPS: RMS Force for %s = %4.2f\n", tset->force[n].sname, tset->force[n].ff_val));
    }

    // If CELL STRESS, calculate stress tensor and compute trace
    for (n = 0; n < nfit[STRESS]; n++) {
      if (strcmp (tset->stress[n].sname, fname[i]) != 0)
	continue;
      stress_flag = 1;
      strcpy (line, "compute cstress all pressure thermo_temp virial");
      if (rank == 0 && debug_level == 1)
	DEBUG_PRINT (("%s\n", line));
      lammps_command (ptr, line);
      strcpy (line, "run 0");
      if (rank == 0 && debug_level == 1)
	DEBUG_PRINT (("%s\n", line));
      lammps_command (ptr, line);
      cell_stress = (double *) lammps_extract_compute (ptr, "cstress", 0, 1);
      if (tset->stress[n].id & XX)
	tset->stress[n].ff_val += cell_stress[1] * cell_stress[1];
      if (tset->stress[n].id & YY)
	tset->stress[n].ff_val += cell_stress[2] * cell_stress[2];
      if (tset->stress[n].id & ZZ)
	tset->stress[n].ff_val += cell_stress[3] * cell_stress[3];
      if (tset->stress[n].id & XY)
	tset->stress[n].ff_val += cell_stress[4] * cell_stress[4];
      if (tset->stress[n].id & XZ)
	tset->stress[n].ff_val += cell_stress[5] * cell_stress[5];
      if (tset->stress[n].id & YZ)
	tset->stress[n].ff_val += cell_stress[6] * cell_stress[6];
      tset->stress[n].ff_val = sqrt (tset->stress[n].ff_val);
      //free (cell_stress);
      strcpy (line, "uncompute cstress");
      if (rank == 0 && debug_level == 1)
	DEBUG_PRINT (("%s\n", line));
      lammps_command (ptr, line);
      if (rank == 0 && debug_level == 1)
	DEBUG_PRINT (("--> LAMMPS: Cell RMS Stress for %s = %4.2f\n", tset->stress[n].sname, tset->stress[n].ff_val));
    }

    // If CELL PRESSURE, calculate scalar pressure
    for (n = 0; n < nfit[PRESSURE]; n++) {
      if (strcmp (tset->press[n].sname, fname[i]) != 0)
	continue;
//      press_flag = 1;
      strcpy (line, "compute press all pressure thermo_temp virial");
      if (rank == 0 && debug_level == 1)
	DEBUG_PRINT (("%s\n", line));
      lammps_command (ptr, line);
      strcpy (line, "run 0");
      if (rank == 0 && debug_level == 1)
	DEBUG_PRINT (("%s\n", line));
      lammps_command (ptr, line);
      if (tset->press[n].id != 0) {	// returns scalar pressure
	cell_press = (double *) lammps_extract_compute (ptr, "press", 0, 0);
	tset->press[n].ff_val = cell_press[1];
      }
      else {			// returns pressure tensor
	cell_press = (double *) lammps_extract_compute (ptr, "press", 0, 1);
	tset->press[n].ff_val = cell_press[tset->press[n].id];
      }
      strcpy (line, "uncompute press");
      if (rank == 0 && debug_level == 1)
	DEBUG_PRINT (("%s\n", line));
      lammps_command (ptr, line);
      if (rank == 0 && debug_level == 1)
	DEBUG_PRINT (("--> LAMMPS: Cell pressure for %s = %4.2f\n", tset->press[n].sname, tset->press[n].ff_val));
    }

    // Calculate geometry cases to set literal values of reference
    // Assumes exact geometries
    if (initial_write && calc_initial_geom) {
      first = 1;
      if (nfit[GEOMETRY] > 0)
	lmp_atoms = (double *) scalloc (3 * natoms, sizeof (double), "lmp_atoms");
      for (n = 0; n < nfit[GEOMETRY]; n++) {
	if (strcmp (tset->geom[n].sname, fname[i]) != 0)
	  continue;
	geom_flag = 1;

	if (first) {		// only do it once, if the structure appears several times in GEOMETRY
	  lammps_gather_atoms (ptr, "x", 1, 3, lmp_atoms);
	  first = 0;
	}

	atoms = (double **) scalloc (tset->geom[n].num_atom, sizeof (double *), "geomatoms");
	for (j = 0; j < tset->geom[n].num_atom; j++) {
	  atoms[j] = (double *) scalloc (3, sizeof (double), "atomj");
	  id = tset->geom[n].atom[j] - 1;
	  atoms[j][0] = lmp_atoms[id * 3];
	  atoms[j][1] = lmp_atoms[id * 3 + 1];
	  atoms[j][2] = lmp_atoms[id * 3 + 2];
	}

	if (rank == 0 && debug_level == 1)
	  DEBUG_PRINT (("--> LAMMPS: Calculating initial geometry for %s ", tset->geom[n].sname));

	if (tset->geom[n].num_atom == 2) {
	  val = calc_bond (atoms);
	  if (rank == 0 && debug_level == 1)
	    DEBUG_PRINT (("(Bond = %4.2f)\n", val));
	}
	else if (tset->geom[n].num_atom == 3) {
	  val = calc_angle (atoms);
	  if (rank == 0 && debug_level == 1)
	    DEBUG_PRINT (("(Angle = %4.2f)\n", val));
	}
	else {
	  val = calc_torsion (atoms);
	  if (rank == 0 && debug_level == 1)
	    DEBUG_PRINT (("(Torsion = %4.2f)\n", val));
	}
	tset->geom[n].lit = val;

	// Cleanup
	for (j = 0; j < tset->geom[n].num_atom; j++)
	  free (atoms[j]);
	free (atoms);
      }
      if (nfit[GEOMETRY] > 0)
	free (lmp_atoms);
    }
    // FORCES, STRESSES and GEOMETRIC pre-calculations do not require structure optimization
//    if (!force_flag && !stress_flag && !pqeq_flag) {
    if (!force_flag && !stress_flag) {

      // Impose any defined restrains before any minimization or MD
      int numrst = 0;
      strcpy (line, "run 0");
      if (rank == 0 && debug_level == 1)
	DEBUG_PRINT (("%s\n", line));
      lammps_command (ptr, line);

      for (n = 0; n < nrst; n++) {
	if (strcmp (fname[i], rstrain[n].name) != 0)
	  continue;
	// minimize molecule energy with restraints
	// using harmonic restrain
	k1 = rstrain[n].f1;
	k2 = rstrain[n].f2;
	if (rstrain[n].nbody == 2)
	  sprintf (line, "fix REST%i all restrain %s %d %d %.3f %.3f %.3f", numrst,
		   "bond", rstrain[n].atom[0], rstrain[n].atom[1], k1, k2, rstrain[n].val);
	else if (rstrain[n].nbody == 3)
	  sprintf (line, "fix REST%i all restrain %s %d %d %d %.3f %.3f %.3f", numrst,
		   "angle", rstrain[n].atom[0], rstrain[n].atom[1], rstrain[n].atom[2], k1, k2, rstrain[n].val);
	else if (rstrain[n].nbody == 4)
	  sprintf (line, "fix REST%i all restrain %s %d %d %d %d %.3f %.3f %.3f", numrst,
		   "dihedral", rstrain[n].atom[0], rstrain[n].atom[1], rstrain[n].atom[2], rstrain[n].atom[3], k1, k2, rstrain[n].val + 180.0);
	else if (rstrain[n].nbody == 5) {	// it's really 3 bodies, for a transition state
	  lmp_atoms = (double *) scalloc (3 * natoms, sizeof (double), "lmp_atoms");
	  lammps_gather_atoms (ptr, "x", 1, 3, lmp_atoms);
	  id = tset->geom[i].atom[rstrain[n].atom[0]] - 1;
	  id2 = tset->geom[i].atom[rstrain[n].atom[1]] - 1;
	  apos = (double *) scalloc (6, sizeof (double), "atoms");
	  apos[0] = lmp_atoms[id * 3];
	  apos[1] = lmp_atoms[id * 3 + 1];
	  apos[2] = lmp_atoms[id * 3 + 2];
	  apos[3] = lmp_atoms[id2 * 3];
	  apos[4] = lmp_atoms[id2 * 3 + 1];
	  apos[5] = lmp_atoms[id2 * 3 + 2];
	  rs1 =
	    sqrt ((apos[0] - apos[3]) * (apos[0] - apos[3]) + (apos[1] - apos[4]) * (apos[1] - apos[4]) + (apos[2] - apos[5]) * (apos[2] - apos[5]));
	  rs2 = rs1 / rstrain[n].val;
	  sprintf (line, "fix REST%ia all restrain %s %d %d %.3f %.3f %.3f", numrst, "bond", rstrain[n].atom[0], rstrain[n].atom[1], k1, k2, rs1);
	  if (rank == 0 && debug_level == 1)
	    DEBUG_PRINT (("%s\n", line));
	  sprintf (line, "fix REST%ib all restrain %s %d %d %.3f %.3f %.3f", numrst, "bond", rstrain[n].atom[0], rstrain[n].atom[1], k1, k2, rs2);
	  free (lmp_atoms);
	  free (apos);
	}
	if (rank == 0 && debug_level == 1)
	  DEBUG_PRINT (("%s\n", line));
	lammps_command (ptr, line);
	if (rstrain[n].nbody == 5) {
	  sprintf (line, "fix_modify REST%ia energy yes", numrst);
	  sprintf (line, "fix_modify REST%ib energy yes", numrst);
	}
	else
	  sprintf (line, "fix_modify REST%i energy yes", numrst);
	if (rank == 0 && debug_level == 1)
	  DEBUG_PRINT (("%s\n", line));
	lammps_command (ptr, line);
	numrst++;
      }

      // Use FIRE when system has restraints in REAXFF
      if (minimization) {
//	if ((nfit[CHARGE] == 0) && (!pqeq_flag)) {	// don't minimize structure if optimizing charges
//        if ((nfit[CHARGE] == 0)) {      // minimize structure if optimizing charges
	  switch (min_method) {
	  case 1:
	    strcpy (line, "min_style fire");
	    break;
	  case 2:
	    strcpy (line, "min_style sd");
	    break;
	  default:
	    strcpy (line, "min_style cg");
	  }
	  if (rank == 0 && debug_level == 1)
	    DEBUG_PRINT (("%s\n", line));
	  lammps_command (ptr, line);

	  // If CELL PARAMETERS case, relax the unit cell before optimizing atomic positions
	  first = 1;
	  for (n = 0; n < nfit[CELL]; n++) {
	    if (strcmp (tset->cell[n].sname, fname[i]) != 0)
	      continue;
	    cell_flag = 1;
	    if (first) {	// only once per cell entry
	      strcpy (line, "fix 2 all box/relax iso 0.0 vmax 0.001");
	      if (rank == 0 && debug_level == 1)
		DEBUG_PRINT (("%s\n", line));
	      lammps_command (ptr, line);
	      first = 0;
	    }
	  }

	  // Need to make sure minimization line search is not zero
	  if (forcefield2optimize == EFF) {
	    strcpy (line, "min_modify line forcezero");
            if (rank == 0 && debug_level == 1)
              DEBUG_PRINT (("%s\n", line));
	    lammps_command (ptr, line);
	  }
	  // Now optimize atomic positions
	  sprintf (line, "minimize 0.0 1.0e-4 %d %d", struc_min_steps, struc_min_steps * 5);
	  if (rank == 0 && debug_level == 1)
	    DEBUG_PRINT (("%s\n", line));
	  lammps_command (ptr, line);
//	}
      }
      else {
	if (cell_flag) {
	  printf ("Need to use minimization to relax box (i.e. remove -f option)");
	  MPI_Abort (MPI_COMM_WORLD, 1);
	}
/*
              strcpy(line, "variable rms_force equal sqrt(fcm(all,x)^2+fcm(all,y)^2+fcm(all,z)^2)");
              if (rank == 0 && debug_level == 1)
                  DEBUG_PRINT(("%s\n", line));
              lammps_command(ptr, line);
*/
	strcpy (line, "run 0");
	if (rank == 0 && debug_level == 1)
	  DEBUG_PRINT (("%s\n", line));
	lammps_command (ptr, line);
/*
              rms_force = (double *) lammps_extract_variable(ptr, "rms_force", "all");
              total_rms_force += *rms_force;
              strcpy(line, "variable rms_force delete");
              if (rank == 0 && debug_level == 1)
                  DEBUG_PRINT(("%s\n", line));
              lammps_command(ptr, line);
              if (rank == 0 && debug_level == 1)
                  DEBUG_PRINT(("--> LAMMPS: RMS Force for %s = %4.2f\n",
                             uname[i], *rms_force));
              free (rms_force);
*/
	lmp_atoms = (double *) scalloc (3 * natoms, sizeof (double), "lmp_atoms");
	lammps_gather_atoms (ptr, "f", 1, 3, lmp_atoms);
	force = (double *) scalloc (3, sizeof (double), "geomforce");
	for (j = 0; j < natoms; j++) {
	  force[0] += lmp_atoms[j * 3] * lmp_atoms[j * 3];
	  force[1] += lmp_atoms[j * 3 + 1] * lmp_atoms[j * 3 + 1];
	  force[2] += lmp_atoms[j * 3 + 2] * lmp_atoms[j * 3 + 2];
	}
	total_rms_force = sqrt (force[0] + force[1] + force[2]);
	free (lmp_atoms);
	free (force);
      }

      // If GEOMETRY case, get coords and calculate corresponding geometric property
      first = 1;
      if (nfit[GEOMETRY] > 0)
	lmp_atoms = (double *) scalloc (3 * natoms, sizeof (double), "lmp_atoms");
      for (n = 0; n < nfit[GEOMETRY]; n++) {
	if (strcmp (tset->geom[n].sname, fname[i]) != 0)
	  continue;
	geom_flag = 1;

	if (first) {		// only do it once, if the structure appears several times in GEOMETRY
	  lammps_gather_atoms (ptr, "x", 1, 3, lmp_atoms);
	  first = 0;
	}

	atoms = (double **) scalloc (tset->geom[n].num_atom, sizeof (double *), "geomatoms");
	for (j = 0; j < tset->geom[n].num_atom; j++) {
	  atoms[j] = (double *) scalloc (3, sizeof (double), "atomj");
	  id = tset->geom[n].atom[j] - 1;
	  atoms[j][0] = lmp_atoms[id * 3];
	  atoms[j][1] = lmp_atoms[id * 3 + 1];
	  atoms[j][2] = lmp_atoms[id * 3 + 2];
	}

	if (rank == 0 && debug_level == 1)
	  DEBUG_PRINT (("--> LAMMPS: geometric properties of %s ", tset->geom[n].sname));

	// Perform geometric or RMS force calculations
	if (tset->geom[n].num_atom == 2) {
	  void read_CGFF (char *, double *);
	  val = calc_bond (atoms);
	  if (rank == 0 && debug_level == 1)
	    DEBUG_PRINT (("(Bond = %4.2f)\n", val));
	}
	else if (tset->geom[n].num_atom == 3) {
	  val = calc_angle (atoms);
	  if (rank == 0 && debug_level == 1)
	    DEBUG_PRINT (("(Angle = %4.2f)\n", val));
	}
	else {
	  val = calc_torsion (atoms);
	  if (rank == 0 && debug_level == 1)
	    DEBUG_PRINT (("(Torsion = %4.2f)\n", val));
	}
	tset->geom[n].ff_val = val;

	// clean atom arrays
	for (j = 0; j < tset->geom[n].num_atom; j++)
	  free (atoms[j]);
	free (atoms);
      }
      if (nfit[GEOMETRY] > 0)
	free (lmp_atoms);

      // If CG force field, we'll do some dynamics after optimizing geometry
      if (forcefield2optimize == CG) {
	strcpy (line, "log log.lammps");
	if (rank == 0 && debug_level == 1)
	  DEBUG_PRINT (("%s\n", line));
	lammps_command (ptr, line);

	strcpy (line, "thermo_style    custom step etotal ke temp pe ebond eangle edihed eimp evdwl ecoul elong press vol");
	if (rank == 0 && debug_level == 1)
	  DEBUG_PRINT (("%s\n", line));
	lammps_command (ptr, line);

	strcpy (line, "dump 1 all custom 1 CG.nvt.lammpstrj id type xu yu zu");
	if (rank == 0 && debug_level == 1)
	  DEBUG_PRINT (("%s\n", line));
	lammps_command (ptr, line);

	strcpy (line, "dump_modify 1 flush yes");
	if (rank == 0 && debug_level == 1)
	  DEBUG_PRINT (("%s\n", line));
	lammps_command (ptr, line);

	strcpy (line, "fix 3 all nvt temp 298 298 100");
	if (rank == 0 && debug_level == 1)
	  DEBUG_PRINT (("%s\n", line));
	lammps_command (ptr, line);

	strcpy (line, "run 100");
	if (rank == 0 && debug_level == 1)
	  DEBUG_PRINT (("%s\n", line));
	lammps_command (ptr, line);

	// now the dynamics is done, hacks to get RDFs from VMD
	// first, vmdRDF.pl calls VMD to make and parse the RDFs
	if (rank == 0 && debug_level == 1)
	  DEBUG_PRINT (("Running ./vmdRDF.pl CG.nvt.lammpstrj\n"));
	system ("/usr/bin/perl ./vmdRDF.pl CG.nvt.lammpstrj");
	// next, read the RDFs written by the perl scripts
	if (rank == 0 && debug_level == 1)
	  DEBUG_PRINT (("Reading RDF files\n"));
	readRDFs ();
      }

      // If CHARGE case, NOTE: Charges are computed AFTER structure minimization
      // except when using -f option
      if (nfit[CHARGE] > 0)
	chrg = (double *) scalloc (natoms, sizeof (double), "charge");
      for (n = 0; n < nfit[CHARGE]; n++) {
	if (strcmp (tset->charge[n].sname, fname[i]) != 0)
	  continue;
	charge_flag = 1;
	strcpy (line, "run 0");
	if (rank == 0 && debug_level == 1)
	  DEBUG_PRINT (("%s\n", line));
	lammps_command (ptr, line);

	if (pqeq_flag) {
	  // 1 0.0 10.0 1.0e-6
	  strcpy (line, "fix pqeq all pqeq 1 0.0 10.0 1.0e-6");
	  if (rank == 0 && debug_level == 1)
	    DEBUG_PRINT (("%s\n", line));
	  lammps_command (ptr, line);
	  strcpy (line, "run 2000");
	  lammps_command (ptr, line);
	}

	lammps_gather_atoms (ptr, "q", 1, 1, chrg);
	tset->charge[n].ff_val = chrg[tset->charge[n].atomid - 1];
	if (rank == 0 && debug_level == 1)
	  DEBUG_PRINT (("--> LAMMPS: Charge for %s = %4.2f\n", tset->charge[n].sname, tset->charge[n].ff_val));
      }
      if (nfit[CHARGE] > 0)
	free (chrg);

      // If ATOM_FORCE case

      if (nfit[ATOM_FORCE] > 0)
	lmp_atoms = (double *) scalloc (3 * natoms, sizeof (double), "lmp_atoms");


      for (n = 0; n < nfit[ATOM_FORCE]; n++) {
	if (strcmp (tset->atom_force[n].sname, fname[i]) != 0)
	  continue;
	atom_force_flag = 1;

	lammps_gather_atoms (ptr, "f", 1, 3, lmp_atoms);
	tset->atom_force[n].ff_val_x = lmp_atoms[3 * (tset->atom_force[n].atomid - 1)];
	tset->atom_force[n].ff_val_y = lmp_atoms[3 * (tset->atom_force[n].atomid - 1) + 1];
	tset->atom_force[n].ff_val_z = lmp_atoms[3 * (tset->atom_force[n].atomid - 1) + 2];

	if (rank == 0 && debug_level == 1)
	  DEBUG_PRINT (("--> LAMMPS: ATOM_FORCE for %s = %4.2f %4.2f %4.2f\n",
			tset->atom_force[n].sname, tset->atom_force[n].ff_val_x, tset->atom_force[n].ff_val_y, tset->atom_force[n].ff_val_z));
      }
      if (nfit[ATOM_FORCE] > 0)
	free (lmp_atoms);

      // Now complete the CELL case by calculating cell parameters on relaxed box with optimized atomic positions
      first = 1;
      for (n = 0; n < nfit[CELL]; n++) {
	if (strcmp (tset->cell[n].sname, fname[i]) != 0)
	  continue;
	//cell_flag = 1;
	if (first) {
	  xlo = (double *) lammps_extract_global (ptr, "boxxlo");
	  xhi = (double *) lammps_extract_global (ptr, "boxxhi");
	  ylo = (double *) lammps_extract_global (ptr, "boxylo");
	  yhi = (double *) lammps_extract_global (ptr, "boxyhi");
	  zlo = (double *) lammps_extract_global (ptr, "boxzlo");
	  zhi = (double *) lammps_extract_global (ptr, "boxzhi");
	  xy = (double *) lammps_extract_global (ptr, "boxxy");
	  xz = (double *) lammps_extract_global (ptr, "boxxz");
	  yz = (double *) lammps_extract_global (ptr, "boxyz");
	  first = 0;
	  strcpy (line, "unfix 2");
	  if (rank == 0 && debug_level == 1)
	    DEBUG_PRINT (("%s\n", line));
	  lammps_command (ptr, line);
	}

	if (strcmp (tset->cell[n].type, "a") == 0)
	  tset->cell[n].ff_val = xhi[0] - xlo[0];	// a
	else if (strcmp (tset->cell[n].type, "b") == 0)
	  tset->cell[n].ff_val = sqrt ((yhi[0] - ylo[0]) * (yhi[0] - ylo[0]) + xy[0] * xy[0]);	// b
	else if (strcmp (tset->cell[n].type, "c") == 0)
	  tset->cell[n].ff_val = sqrt (xz[0] * xz[0] + yz[0] * yz[0] + (zhi[0] - zlo[0]) * (zhi[0] - zlo[0]));	//c
	if (rank == 0 && debug_level == 1)
	  DEBUG_PRINT (("--> LAMMPS: Cell parameter %s = %4.2f for %s\n", tset->cell[n].type, tset->cell[n].ff_val, tset->cell[n].sname));
/*                free (xlo); free (ylo); free (zlo);
                free (xhi); free (yhi); free (zhi);
                free (xy); free (xz); free (yz);
*/ }
      if (final_write) {
	sprintf (line, "dump 1 all xyz 1 %s.last.xyz", fname[i]);
	if (rank == 0 && debug_level == 1)
	  DEBUG_PRINT (("%s\n", line));
	lammps_command (ptr, line);
	sprintf (line, "dump_modify 1 element %s", atypes[i]);
	if (rank == 0 && debug_level == 1)
	  DEBUG_PRINT (("%s\n", line));
	lammps_command (ptr, line);
	strcpy (line, "run 0");
	if (rank == 0 && debug_level == 1)
	  DEBUG_PRINT (("%s\n", line));
	lammps_command (ptr, line);
	strcpy (line, "undump 1");
	if (rank == 0 && debug_level == 1)
	  DEBUG_PRINT (("%s\n", line));
	lammps_command (ptr, line);
      }
    }

    // Extract all energy components
    if ((forcefield2optimize == REAX) || (forcefield2optimize == REAXC)) {
      if (pqeq_flag) {
	ff_energies = (double *) lammps_extract_compute (ptr, "reax", 0, 0);
	ff_tmp = (double *) lammps_extract_compute (ptr, "pqeq", 0, 0);
	*ff_energies = *ff_energies + *ff_tmp;
      }
      else
	ff_energies = (double *) lammps_extract_compute (ptr, "reax", 0, 0);
    }
    else if (forcefield2optimize == EFF)
      ff_energies = (double *) lammps_extract_compute (ptr, "eff", 0, 0);
    else if (forcefield2optimize == CG)
      ff_energies = (double *) lammps_extract_compute (ptr, "cg", 0, 0);
    else if (forcefield2optimize == MORSE)
      ff_energies = (double *) lammps_extract_compute (ptr, "emorse", 0, 0);
    else if (forcefield2optimize == COMB)
      ff_energies = (double *) lammps_extract_compute (ptr, "ecomb", 0, 0);
    else if (forcefield2optimize == TERSOFF)
      ff_energies = (double *) lammps_extract_compute (ptr, "etersoff", 0, 0);
    else if (forcefield2optimize == ZHOU_EAM)
      ff_energies = (double *) lammps_extract_compute (ptr, "ezhou_EAM", 0, 0);
    else if (forcefield2optimize == TERSOFF_MOD)
      ff_energies = (double *) lammps_extract_compute (ptr, "etersoff_mod", 0, 0);

    // Just use the total energy
    energies[i] = *ff_energies;
//        free (ff_energies);

    if (rank == 0 && debug_level == 1)
      DEBUG_PRINT (("--> LAMMPS: (%s) Total energy = %4.2f\n", fname[i], energies[i]));
//      ff_energies = NULL;

    // Clear the LAMMPS memory
    if (rank == 0 && debug_level == 1)
      DEBUG_PRINT (("clear\n"));
    lammps_command (ptr, "clear");

//    lammps_close (ptr);
  }

  total_error = calcerror (energies, rank);
  calc_initial_geom = 0;

  // Cleanup
  free (line);
  free (energies);

  //  free (ff_energies);
  lammps_close (ptr);

  if (minimization)
    return total_error;
  else
    return total_rms_force;
}


/****************************************************************************
*  Function     : calc_bond
*  Description  : Calculates the bond distance between two atoms
*  Parameters   : The atom coordinates array
*  Effects      : Returns scalar value of bond distance
****************************************************************************/

double calc_bond (double **atom)
{
  int i;
  double bond;

  bond = 0.0;
  for (i = 0; i < 3; i++)
    bond += pow ((atom[0][i] - atom[1][i]), 2);

  bond = sqrt (bond);

  return bond;
}

/****************************************************************************
*  Function     : calc_angle
*  Description  : Calculates the angle between three atoms
*  Parameters   : The atom coordinates array
*  Effects      : Returns scalar value of angle
****************************************************************************/

double calc_angle (double **atom)
{
  double theta, c, s;
  double delx1, delx2, dely1, dely2, delz1, delz2;
  double r1, r2, rsq1, rsq2;

  theta = 0.0;
  // 1st bond

  delx1 = atom[0][0] - atom[1][0];
  dely1 = atom[0][1] - atom[1][1];
  delz1 = atom[0][2] - atom[1][2];
//  domain->minimum_image(delx1,dely1,delz1);

  rsq1 = delx1 * delx1 + dely1 * dely1 + delz1 * delz1;
  r1 = sqrt (rsq1);

  // 2nd bond

  delx2 = atom[2][0] - atom[1][0];
  dely2 = atom[2][1] - atom[1][1];
  delz2 = atom[2][2] - atom[1][2];
//  domain->minimum_image(delx2,dely2,delz2);

  rsq2 = delx2 * delx2 + dely2 * dely2 + delz2 * delz2;
  r2 = sqrt (rsq2);

  // angle (cos and sin)

  c = delx1 * delx2 + dely1 * dely2 + delz1 * delz2;
  c /= r1 * r2;

  if (c > 1.0)
    c = 1.0;
  if (c < -1.0)
    c = -1.0;

  s = sqrt (1.0 - c * c);
  if (s < SMALL)
    s = SMALL;
  s = 1.0 / s;

//  theta = atan2 (s,c);
  theta = acos (c);
  return theta * 180 / PI;
}

/****************************************************************************
*  Function     : calc_torsion
*  Description  : Calculates the torsion angle between four atoms
*  Parameters   : The atom coordinates array
*  Effects      : Returns scalar value of torsion angle
****************************************************************************/

double calc_torsion (double **atom)
{
  double theta, costh;

  // C1-C2-C3-C4
  // v[0] = R(C2-C1)    
  // v[1] = R(C3-C2)   
  // v[2] = R(C4-C3)  
  // n    = v[1]xv[2]
  double v[3][3];		//vectors of the three bonds
  double n[3];
  int i, j, k;
  double ip00, ip11, ip22, ip01, ip02, ip12, ip01x2;
  double x, y;			//x: component of v[0] on projected (perpendicular to v[1]) v[2]
  //y: component of v[0] on the normal vector of C2-C3-C4 plane

  ip00 = ip11 = ip22 = ip01 = ip02 = ip12 = ip01x2 = 0.0;

  for (k = 0; k < 3; k++) {
    v[0][k] = (atom[1][k]) - (atom[0][k]);
    v[1][k] = (atom[2][k]) - (atom[1][k]);
    v[2][k] = (atom[3][k]) - (atom[2][k]);
    ip00 += (v[0][k] * v[0][k]);
    ip11 += (v[1][k] * v[1][k]);
    ip22 += (v[2][k] * v[2][k]);
    ip01 += (v[0][k] * v[1][k]);
    ip02 += (v[0][k] * v[2][k]);
    ip12 += (v[1][k] * v[2][k]);
  }

  i = 1;
  j = 2;
  n[0] = v[i][1] * v[j][2] - v[i][2] * v[j][1];
  n[1] = v[i][2] * v[j][0] - v[i][0] * v[j][2];
  n[2] = v[i][0] * v[j][1] - v[i][1] * v[j][0];

  for (k = 0; k < 3; k++)
    ip01x2 += (v[0][k] * n[k]);
  //x=-v[0]*v[2]+(v[0]*v[1])(v[1]*v[2])= v[0]*(-v[2]+(v[1]*v[2])v[1])
  x = -ip02 / sqrt (ip00 * ip22) + (ip01 * ip12) / (ip11 * sqrt (ip00 * ip22));
  //y=v[0]*(v[1]xv[2])
  y = -ip01x2 / sqrt (ip00 * ip11 * ip22);
  //printf("%f %f %f %f %f %f %f \n",ip00,ip11,ip22,ip01,ip02,ip12,ip01x2);
  //printf("%f %f \n",x,y);

  costh = x / sqrt (x * x + y * y);
  theta = acos (costh);
  if (y < 0)
    theta *= -1;

  return theta * 180 / PI;

}
