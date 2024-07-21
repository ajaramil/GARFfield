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

#include "tool_box.h"
#include "restraints.h"
#include "structures.h"
#include "mpi.h"

enum {REAX, REAXC, EFF, PQEQ, PPR, CG, MORSE, COMB, TERSOFF, TERSOFF_MOD, ZHOU_EAM };

#ifdef DEBUG
# define DEBUG_PRINT(x) printf x
#else
# define DEBUG_PRINT(x)
#endif

/****************************************************************************
*  Function     : dotProduct
*  Description  : Calculates dot product between two vectors a and b of size length
*  Parameters   : Array a and b, and array size in length
*  Effects      : No side effects, returns scalar value of dot product
****************************************************************************/

double dotProduct (double *a, double *b, int length)
{
  double runningSum = 0;
  int index;
  for (index = 0; index < length; index++)
    runningSum += a[index] * b[index];
  return runningSum;
}

/****************************************************************************
*  Function     : geo2data
*  Description  : Maps the biograf (bgf) or xyz format structure files in the user-defined
                  input geo file onto lammps data format files
*  Parameters   : bgf/xyz structure file, data file pointer to array pointer, atom types 
                  pointer to array pointer, unit cell type flag (periodic or finite)
                  array pointer, data file names array pointer, reax force field
                  data structure, total number of structure files in bgf/xyz file,
                  force field to optimize (-F reax or reaxc), and MPI process rank
*  Effects      : Writes out all the molecular structure data files required
                  by the LAMMPS-calculated fitness function, counts/returns the total
                  number of structures in the input bgf/xyz file, allocates and populates
                  several return arrays (dfile, atypes, cellfalg, fname)
****************************************************************************/

void geo2data (char *sfile, char ***dfile, char ***atypes, int **cellflag, char ***fname, ffieldtype * fftype, int *num, int forcefield2optimize,
	       int rank)
{
  FILE *fi, *fo;
  char massstr[MAX_LINE], type_id[MAX_LINE], typestr[MAX_LINE], type_idstr[MAX_LINE];
  char *fstr, *stok, ch;
  double bounds[6], cell[9];
  char filen[MAX_LINE];
  char *line, *fout, *aname, *header, *junk, *hjunk;

  int i, j, c, k, atype_id, n, has_cell, current_rst, ntypes;

  if ((forcefield2optimize == REAX) || (forcefield2optimize == REAXC))
    ntypes = fftype->rxff->num_atom_types;
  else if (forcefield2optimize == PQEQ)
    ntypes = fftype->pqeqff->num_atom_types;
  else if (forcefield2optimize == PPR)
    ntypes = fftype->pprff->num_atom_types;
  else if (forcefield2optimize == MORSE)
    ntypes = fftype->morseff->num_atom_types;
  else if (forcefield2optimize == COMB)
    ntypes = fftype->combff->num_atypes;
  else if (forcefield2optimize == TERSOFF)
    ntypes = fftype->tersoffff->num_atypes;
  else if (forcefield2optimize == ZHOU_EAM)
    ntypes = fftype->zhou_EAMff->num_atypes;
  else if (forcefield2optimize == TERSOFF_MOD)
    ntypes = fftype->tersoff_modff->num_atypes;

  double x, y, z, amass;
  char **strtmp;
  int file_format = 0;		// 0=.bgf or 1=.xyz (extended)

  // Allocate temporals
  hjunk = (char *) scalloc (MAX_LINE, sizeof (char), "line");
  junk = (char *) scalloc (MAX_LINE, sizeof (char), "line");
  header = (char *) scalloc (MAX_LINE, sizeof (char), "line");
  aname = (char *) scalloc (MAX_LINE, sizeof (char), "line");
  line = (char *) scalloc (MAX_LINE, sizeof (char), "line");
  fout = (char *) scalloc (MAX_LINE, sizeof (char), "fout");
  strtmp = (char **) scalloc (MAX_TOKENS, sizeof (char *), "strtmp");
  for (i = 0; i < MAX_TOKENS; i++)
    strtmp[i] = (char *) scalloc (MAX_TOKEN_LEN, sizeof (char), "tmpi");
  double **Rot;
  Rot = (double **) scalloc (3, sizeof (double), "structure:rot");
  for (i = 0; i < 3; i++)
    Rot[i] = (double *) scalloc (3, sizeof (double), "structure:rot");

  //rudimentary error checking
  sprintf (filen, "%s.%d", sfile, rank);
  fi = fopen (filen, "r");
  if (fi == NULL) {
    if (rank == 0)
      fprintf (stderr, "ERROR: opening the geo file %s terminating...\n", sfile);
    MPI_Abort (MPI_COMM_WORLD, 1);
  }

  //count structures and restraints
  nrst = 0;
  while (fgets (line, MAX_LINE, fi) != NULL) {
    if (strstr (line, "BIOGRF") || strstr (line, "XTLGRF") || strstr (line, "XYZ"))
      (*num)++;
    if (strstr (line, "RESTRAINT"))
      nrst++;
  }
  if (rank == 0) {
    DEBUG_PRINT ((">> Found %d structures in %s\n", *num, sfile));
    DEBUG_PRINT ((">> Found %d restraints in %s\n", nrst, sfile));
  }

  if ((*num) == 0) {
    if (rank == 0)
      printf ("ERROR: No valid structures read from %s!\n", sfile);
    MPI_Abort (MPI_COMM_WORLD, 1);
  }
  fclose (fi);

  //initialize arrays and variables

  *dfile = (char **) scalloc ((*num), sizeof (char *), "dfile");
  *atypes = (char **) scalloc ((*num), sizeof (char *), "atypes");
  *fname = (char **) scalloc ((*num), sizeof (char *), "fnames");
  *cellflag = (int *) scalloc ((*num), sizeof (int), "cellflag");
  for (i = 0; i < (*num); i++) {
    (*dfile)[i] = (char *) scalloc (MAX_LINE, sizeof (char), "dfilei");
    (*fname)[i] = (char *) scalloc (MAX_LINE, sizeof (char), "fnamei");
    (*atypes)[i] = (char *) scalloc (MAX_LINE, sizeof (char), "atypei");
    (*cellflag)[i] = 0;
  }

  // restraints
  if (nrst > 0)
    rstrain = (restraint *) scalloc (nrst, sizeof (restraint), "rstrain");
  //  for (i = 0; i < nrst; i++)
  //    rstrain[i] = (restraint *) scalloc (1, sizeof (restraint), "rsti");
  current_rst = 0;

  //get data from geo file
  sprintf (filen, "%s.%d", sfile, rank);
  if ((fi = fopen (filen, "r")) == NULL) {
    if (rank == 0)
      fprintf (stderr, "ERROR: opening structure file %s -- terminating !\n", sfile);
    MPI_Abort (MPI_COMM_WORLD, 1);
  }
  i = -1;
  while (fgets (line, MAX_LINE, fi) != NULL) {
    sscanf (line, "%c", &ch);	// to check atom number line in xyz file
    if (strstr (line, "BIOGRF") || strstr (line, "XTLGRF") || strstr (line, "XYZ")) {
      if (strstr (line, "XYZ")) {
	sscanf (line, "XYZ %s", hjunk);
	file_format = 1;	// xyz
      }
      if (i > -1 && fo != NULL) {	//write previous structure data
	fclose (fo);
	write_fdata (header, massstr, bounds, has_cell, cell, n, j, i, rank, file_format);
	(*cellflag)[i] = has_cell;
	sprintf ((*dfile)[i], "data.%s", header);
	strcpy ((*atypes)[i], type_idstr);
	strcpy ((*fname)[i], header);
      }
      i++;
      if (rank == 0 && debug_level == 1)
	DEBUG_PRINT ((">> Reading structure %d of %d\n", (i + 1), *num));
      j = n = 0;
      has_cell = 0;
      bounds[0] = bounds[2] = bounds[4] = LARGE;
      bounds[1] = bounds[3] = bounds[5] = SMALL;
      cell[0] = cell[1] = cell[2] = cell[3] = cell[4] = cell[5] = cell[6] = cell[7] = cell[8] = 0.0;

      strcpy (typestr, "");
      strcpy (type_idstr, "");
      strcpy (massstr, "");
      strcpy (header, "created by reax_lammps");
      sprintf (fout, "_tmp.dat.%d", rank);
      if ((fo = fopen (fout, "w")) == NULL) {
	if (rank == 0)
	  printf ("ERROR: Cannot create LAMMPS data file %s!\n", fout);
	MPI_Abort (MPI_COMM_WORLD, 1);
      }
      fprintf (fo, "Atoms\n\n");
    }
    else if (i > -1 && (strstr (line, "CRYSTX"))) {
      k = sscanf (line, "%s %lf %lf %lf %lf %lf %lf", junk, &cell[0], &cell[1], &cell[2], &cell[3], &cell[4], &cell[5]);
      if (k == 7 && cell[0] > 0.0 && cell[1] > 0.0 && cell[2] > 0.0 && cell[3] > 0.0 && cell[4] > 0.0 && cell[5] > 0.0)
	has_cell = 1;

      init_cell (cell, i, rank, file_format, Rot);
    }
    else if (i > -1 && (strstr (line, "Lattice"))) {
      stok = strtok (line, "Lattice=\"");
      k =
	sscanf (stok, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", &cell[0], &cell[1], &cell[2], &cell[3], &cell[4], &cell[5], &cell[6], &cell[7],
		&cell[8]);
      if (k == 9)
	has_cell = 1;

      init_cell (cell, i, rank, file_format, Rot);
    }
    else if (i > -1 && sscanf (line, "DESCRP %s", junk)) {
      strcpy (header, junk);
    }
    else if (i > -1 && strstr (line, "RESTRAINT")) {
      tokenize_string (line, &strtmp);
      rstrain[current_rst].atom[0] = atoi (strtmp[2]);
      rstrain[current_rst].atom[1] = atoi (strtmp[3]);
      strcpy (rstrain[current_rst].name, header);

      if (strstr (strtmp[0], "BOND")) {
	rstrain[current_rst].nbody = 2;
	rstrain[current_rst].val = atof (strtmp[4]);
	rstrain[current_rst].f1 = atof (strtmp[5]);
	rstrain[current_rst].f2 = atof (strtmp[6]);
	if (rank == 0 && debug_level == 1)
	  DEBUG_PRINT ((">> Bond restraint on %s: atoms %d %d %4.2f %4.2f %4.2f\n", header,
			rstrain[current_rst].atom[0], rstrain[current_rst].atom[1], rstrain[current_rst].val, rstrain[current_rst].f1,
			rstrain[current_rst].f2));
      }
      else if (strstr (strtmp[0], "TS")) {
	rstrain[current_rst].nbody = 5;
	rstrain[current_rst].atom[2] = atoi (strtmp[4]);
	rstrain[current_rst].val = atof (strtmp[5]);
	rstrain[current_rst].f1 = atof (strtmp[6]);
	rstrain[current_rst].f2 = atof (strtmp[7]);
	if (rank == 0 && debug_level == 1)
	  DEBUG_PRINT ((">> Transition State restraint on %s: atoms %d %d %d %4.2f %4.2f %4.2f\n",
			rstrain[current_rst].name, rstrain[current_rst].atom[0],
			rstrain[current_rst].atom[1], rstrain[current_rst].atom[2], rstrain[current_rst].val, rstrain[current_rst].f1,
			rstrain[current_rst].f2));
      }
      else if (strstr (strtmp[0], "ANGLE")) {
	rstrain[current_rst].nbody = 3;
	rstrain[current_rst].atom[2] = atoi (strtmp[4]);
	rstrain[current_rst].val = atof (strtmp[5]);
	rstrain[current_rst].f1 = atof (strtmp[6]);
	rstrain[current_rst].f2 = atof (strtmp[7]);
	if (rank == 0 && debug_level == 1)
	  DEBUG_PRINT ((">> Angle restraint on %s: atoms %d %d %d %4.2f %4.2f %4.2f\n", header,
			rstrain[current_rst].atom[0], rstrain[current_rst].atom[1],
			rstrain[current_rst].atom[2], rstrain[current_rst].val, rstrain[current_rst].f1, rstrain[current_rst].f2));
      }
      else if (strstr (strtmp[0], "TORSION")) {
	rstrain[current_rst].nbody = 4;
	rstrain[current_rst].atom[2] = atoi (strtmp[4]);
	rstrain[current_rst].atom[3] = atoi (strtmp[5]);
	rstrain[current_rst].val = atof (strtmp[6]);
	rstrain[current_rst].f1 = atof (strtmp[7]);
	rstrain[current_rst].f2 = atof (strtmp[8]);
	if (rank == 0 && debug_level == 1)
	  DEBUG_PRINT ((">> Torsion restraint on %s: atoms %d %d %d %d %4.2f %4.2f %4.2f\n", header,
			rstrain[current_rst].atom[0], rstrain[current_rst].atom[1],
			rstrain[current_rst].atom[2], rstrain[current_rst].atom[3], rstrain[current_rst].val, rstrain[current_rst].f1,
			rstrain[current_rst].f2));
      }
      else {
	if (rank == 0)
	  printf ("ERROR: Unkown RESTRAINT in %s.\n\n", rstrain[current_rst].name);
	MPI_Abort (MPI_COMM_WORLD, 1);
      }
      current_rst++;
    }
    else if (i > -1 && (strstr (line, "HETATM"))) {
      j++;			//number of atoms
      sscanf (line, "%s %s %s %lf %lf %lf %s", junk, junk, aname, &x, &y, &z, type_id);
      if (x < bounds[0])
	bounds[0] = x;
      if (x > bounds[1])
	bounds[1] = x;
      if (y < bounds[2])
	bounds[2] = y;
      if (y > bounds[3])
	bounds[3] = y;
      if (z < bounds[4])
	bounds[4] = z;
      if (z > bounds[5])
	bounds[5] = z;

      atype_id = get_atom_type_id (type_id, fftype, forcefield2optimize, ntypes, rank);
      if ((forcefield2optimize == REAX) || (forcefield2optimize == REAXC))
	amass = fftype->rxff->sbp[atype_id].mass;
      else if (forcefield2optimize == PQEQ)
        amass = fftype->pqeqff->pqeqatypes[atype_id].mass;
      else if (forcefield2optimize == PPR)
        amass = fftype->pprff->ppratypes[atype_id].mass;
      else if (forcefield2optimize == MORSE)
	amass = fftype->morseff->morseatypes[atype_id].mass;
      else if (forcefield2optimize == COMB)
	amass = fftype->combff->masses[atype_id];
      else if (forcefield2optimize == TERSOFF)
	amass = fftype->tersoffff->masses[atype_id];
      else if (forcefield2optimize == ZHOU_EAM)
	amass = fftype->zhou_EAMff->masses[atype_id];
      else if (forcefield2optimize == TERSOFF_MOD)
	amass = fftype->tersoff_modff->masses[atype_id];

      // AJB: this causes unwanted behavior (e.g. silicates, Ca2_O4_Si0)
      //      sprintf (junk, " %s ", type_id);
      sprintf (junk, " %s", type_id);
      //      fstr = strstr (typestr, junk);
      fstr = strstr (typestr, type_id);
      if (!fstr) {		//recorded new type
	n++;			//number of atom types
	c = n;
	strcat (typestr, junk);
	sprintf (junk, "%d %lf\n", c, amass);
	strcat (massstr, junk);
	if (forcefield2optimize == 0)
	  sprintf (junk, "%d ", (atype_id + 1));
	else
	  sprintf (junk, "%s ", type_id);
	strcat (type_idstr, junk);
      }
      // AJB: related to same problem above (e.g. silicates, Ca2_O4_Si)
      else {
	c = (fstr - typestr + 1) / 3 + 1;
	//        if (strlen(fstr)%2 != 0)
	//          c = (fstr-typestr+1)/3;
	//        printf("atomtype=%d fstr=%d [%s] typestr=%d [%s] c=%d\n", atype_id, fstr, fstr, typestr, typestr, c);
      }
      if (forcefield2optimize == REAX || forcefield2optimize == REAXC || forcefield2optimize == COMB) {
	if (!pqeq_flag)
	  fprintf (fo, "%d %d %lf %lf %lf %lf\n", j, c, 0.00, x, y, z);
	else
	  fprintf (fo, "%d %d %lf %lf %lf %lf 0.0 0.0 0.0\n", j, c, -1.00, x, y, z);
      }
      else if (forcefield2optimize == MORSE || (forcefield2optimize == TERSOFF) || (forcefield2optimize == PPR) 
	       || (forcefield2optimize == TERSOFF_MOD) || (forcefield2optimize == ZHOU_EAM))
	fprintf (fo, "%d %d %lf %lf %lf\n", j, c, x, y, z);
    }
    else if (i > -1 && isdigit (ch)) {
      //      sscanf (line, "%d", &j);        // num atoms in xyz format file
      //      file_format = 1;
      strcpy (header, hjunk);
      continue;
    }
    else if (i > -1 && file_format && (strcmp (line, "\n") != 0)) {
      j++;
      sscanf (line, "%s %lf %lf %lf", type_id, &x, &y, &z);
      if (x < bounds[0])
	bounds[0] = x;
      if (x > bounds[1])
	bounds[1] = x;
      if (y < bounds[2])
	bounds[2] = y;
      if (y > bounds[3])
	bounds[3] = y;
      if (z < bounds[4])
	bounds[4] = z;
      if (z > bounds[5])

	atype_id = get_atom_type_id (type_id, fftype, forcefield2optimize, ntypes, rank);
      // Customize if not using reax force field
      if ((forcefield2optimize == REAX) || (forcefield2optimize == REAXC))
	amass = fftype->rxff->sbp[atype_id].mass;
      else if (forcefield2optimize == PQEQ)
        amass = fftype->pqeqff->pqeqatypes[atype_id].mass;
      else if (forcefield2optimize == PPR)
        amass = fftype->pprff->ppratypes[atype_id].mass;
      else if (forcefield2optimize == MORSE)
	amass = fftype->morseff->morseatypes[atype_id].mass;
      else if (forcefield2optimize == COMB)
	amass = fftype->combff->masses[atype_id];
      else if (forcefield2optimize == TERSOFF)
	amass = fftype->tersoffff->masses[atype_id];
      else if (forcefield2optimize == ZHOU_EAM)
	amass = fftype->zhou_EAMff->masses[atype_id];
      else if (forcefield2optimize == TERSOFF_MOD)
	amass = fftype->tersoff_modff->masses[atype_id];

      sprintf (junk, " %s", type_id);
      //      fstr = strstr (typestr, junk);
      fstr = strstr (typestr, type_id);
      if (!fstr) {		//recorded new type
	n++;			//number of atom types
	c = n;
	strcat (typestr, junk);
	sprintf (junk, "%d %lf\n", c, amass);
	strcat (massstr, junk);
	if (forcefield2optimize == 0)
	  sprintf (junk, "%d ", (atype_id + 1));
	else
	  sprintf (junk, "%s ", type_id);
	strcat (type_idstr, junk);
      }
      else
	c = ((fstr - typestr + 1) / 3 + 1);
      // atom_type charge in LAMMPS
      // for triclinic cells, need to apply orientation changes to coordinates
      double rotx, roty, rotz;
      rotx = x * Rot[0][0] + y * Rot[1][0] + z * Rot[2][0];
      roty = x * Rot[0][1] + y * Rot[1][1] + z * Rot[2][1];
      rotz = x * Rot[0][2] + y * Rot[1][2] + z * Rot[2][2];
      if (forcefield2optimize == REAX || forcefield2optimize == REAXC || forcefield2optimize == COMB) {
	if (!pqeq_flag)
	  fprintf (fo, "%d %d %lf %lf %lf %lf\n", j, c, 0.00, rotx, roty, rotz);
	else
	  fprintf (fo, "%d %d %lf %lf %lf %lf 0.0 0.0 0.0\n", j, c, -1.00, rotx, roty, rotz);
      }
      else if (forcefield2optimize == MORSE || (forcefield2optimize == PPR) || (forcefield2optimize == TERSOFF)
	       || (forcefield2optimize == TERSOFF_MOD) || (forcefield2optimize == ZHOU_EAM))
	fprintf (fo, "%d %d %lf %lf %lf\n", j, c, rotx, roty, rotz);
    }
  }

  fclose (fi);
  fclose (fo);
  write_fdata (header, massstr, bounds, has_cell, cell, n, j, i, rank, file_format);
  (*cellflag)[i] = has_cell;
  sprintf ((*dfile)[i], "data.%s", header);
  strcpy ((*atypes)[i], type_idstr);
  strcpy ((*fname)[i], header);

  for (i = 0; i < MAX_TOKENS; i++)
    free (strtmp[i]);
  free (strtmp);
  free (line);
  free (fout);
  free (aname);
  free (header);
  free (junk);
  free (hjunk);
  for (i = 0; i < 3; i++)
    free (Rot[i]);
  free (Rot);

  // deallocate memory
  /*  for (i = 0; i < (*num); i++) {
     free ((*dfile)[i]);
     free ((*fname)[i]);
     free ((*atypes)[i]);
     }
     free (*dfile);
     free (*atypes);
     free (*fname);
     free (*cellflag);
   */
}

void init_cell (double *cell, int i, int rank, int file_format, double **Rot)
{
  double sq_c, cos_alpha, cos_gamma, sin_gamma, cos_beta, sin_beta;
  double a, b, c, A, B, C, xy, xz, yz;
  double alpha, beta, gamma;

  alpha = beta = gamma = 0.0;
  i = i + 0;

  if (file_format == 0) {	// bgf
    A = cell[0];
    B = cell[1];
    C = cell[2];
    alpha = cell[3];
    beta = cell[4];
    gamma = cell[5];

    sq_c = C * C;
    cos_alpha = cos (alpha * PI / 180.0);
    cos_gamma = cos (gamma * PI / 180.0);
    sin_gamma = sin (gamma * PI / 180.0);
    cos_beta = cos (beta * PI / 180.0);
    sin_beta = sin (beta * PI / 180.0);

    a = A;
    xy = B * cos_gamma;
    xz = C * cos_beta;
    b = B * sin_gamma;
    yz = C * (cos_alpha - cos_gamma * cos_beta) / sin_gamma;
    // Updated to match lammps description by Byung-Hyung Kim (Uppsala University, 6/7/2016)
    c = sqrt (C * C - xz * xz - yz * yz);

  }
  else {			// xyz
    double *vnorma, *vnormb, *vnormc;
    double prod;

    vnorma = (double *) scalloc (3, sizeof (double), "normv");;
    vnormb = (double *) scalloc (3, sizeof (double), "normv");;
    vnormc = (double *) scalloc (3, sizeof (double), "normv");;

    A = sqrt (cell[0] * cell[0] + cell[1] * cell[1] + cell[2] * cell[2]);
    B = sqrt (cell[3] * cell[3] + cell[4] * cell[4] + cell[5] * cell[5]);
    C = sqrt (cell[6] * cell[6] + cell[7] * cell[7] + cell[8] * cell[8]);

    // alpha in degrees
    vnormb[0] = cell[3] / B;
    vnormb[1] = cell[4] / B;
    vnormb[2] = cell[5] / B;

    vnormc[0] = cell[6] / C;
    vnormc[1] = cell[7] / C;
    vnormc[2] = cell[8] / C;

    prod = dotProduct (vnormb, vnormc, 3);

    if (fabs (prod) < 1.)
      alpha = acos (prod) * 180.0 / PI;
    else {
      if (rank == 0)
	printf ("Error: scalar product alpha %f\n", prod);
      MPI_Abort (MPI_COMM_WORLD, 1);
    }

    // beta in degrees
    vnorma[0] = cell[0] / A;
    vnorma[1] = cell[1] / A;
    vnorma[2] = cell[2] / A;

    vnormc[0] = cell[6] / C;
    vnormc[1] = cell[7] / C;
    vnormc[2] = cell[8] / C;

    prod = dotProduct (vnorma, vnormc, 3);

    if (fabs (prod) < 1.)
      beta = acos (prod) * 180.0 / PI;
    else {
      if (rank == 0)
	printf ("Error: scalar product beta %f\n", prod);
      MPI_Abort (MPI_COMM_WORLD, 1);
    }

    // gamma in degrees
    vnorma[0] = cell[0] / A;
    vnorma[1] = cell[1] / A;
    vnorma[2] = cell[2] / A;

    vnormb[0] = cell[3] / B;
    vnormb[1] = cell[4] / B;
    vnormb[2] = cell[5] / B;

    prod = dotProduct (vnorma, vnormb, 3);

    if (fabs (prod) < 1.)
      gamma = acos (prod) * 180.0 / PI;
    else {
      if (rank == 0)
	printf ("Error: scalar product gamma %f\n", prod);
      MPI_Abort (MPI_COMM_WORLD, 1);
    }
    sq_c = C * C;
    cos_alpha = cos (alpha * PI / 180.0);
    cos_gamma = cos (gamma * PI / 180.0);
    sin_gamma = sin (gamma * PI / 180.0);
    cos_beta = cos (beta * PI / 180.0);
    sin_beta = sin (beta * PI / 180.0);

    a = A;
    xy = B * cos_gamma;
    xz = C * cos_beta;
    b = B * sin_gamma;
    yz = C * (cos_alpha - cos_gamma * cos_beta) / sin_gamma;
    c = sqrt (C * C - xz * xz - yz * yz);
  }

  // unit cell spatial transformations, rotation/inversion
  // Updated to match lammps' triclinic description 

  int m, n;
  double **aa, **bb, **Apre;
  aa = (double **) scalloc (3, sizeof (double), "structure:aa");
  for (m = 0; m < 3; m++)
    aa[m] = (double *) scalloc (3, sizeof (double), "structure:aa");
  bb = (double **) scalloc (3, sizeof (double), "structure:bb");
  for (m = 0; m < 3; m++)
    bb[m] = (double *) scalloc (3, sizeof (double), "structure:bb");
  Apre = (double **) scalloc (3, sizeof (double), "structure:apre");
  for (m = 0; m < 3; m++)
    Apre[m] = (double *) scalloc (3, sizeof (double), "structure:apre");

  // populate with cell vectors
  // contributed by Byung-Hyung Kim (Uppsala University, 6/7/2016)
  for (m = 0; m < 3; m++) {
    for (n = 0; n < 3; n++) {
      aa[m][n] = cell[m * 3 + n];
    }
  }

  double determinant;
  for (m = 0; m < 3; m++)
    determinant = determinant + (aa[0][m] * (aa[1][(m + 1) % 3] * aa[2][(m + 2) % 3] - aa[1][(m + 2) % 3] * aa[2][(m + 1) % 3]));
  for (n = 0; n < 3; n++) {
    for (m = 0; m < 3; m++)
      bb[n][m] =
	((aa[(m + 1) % 3][(n + 1) % 3] * aa[(m + 2) % 3][(n + 2) % 3]) - (aa[(m + 1) % 3][(n + 2) % 3] * aa[(m + 2) % 3][(n + 1) % 3])) / determinant;
  }

  Apre[0][0] = a;
  Apre[0][1] = 0.0;
  Apre[0][2] = 0.0;
  Apre[1][0] = xy;
  Apre[1][1] = b;
  Apre[1][2] = 0.0;
  Apre[2][0] = xz;
  Apre[2][1] = yz;
  Apre[2][2] = c;

  // Construct rotation matrix
  Rot[0][0] = bb[0][0] * Apre[0][0] + bb[0][1] * Apre[1][0] + bb[0][2] * Apre[2][0];
  Rot[0][1] = bb[0][0] * Apre[0][1] + bb[0][1] * Apre[1][1] + bb[0][2] * Apre[2][1];
  Rot[0][2] = bb[0][0] * Apre[0][2] + bb[0][1] * Apre[1][2] + bb[0][2] * Apre[2][2];
  Rot[1][0] = bb[1][0] * Apre[0][0] + bb[1][1] * Apre[1][0] + bb[1][2] * Apre[2][0];
  Rot[1][1] = bb[1][0] * Apre[0][1] + bb[1][1] * Apre[1][1] + bb[1][2] * Apre[2][1];
  Rot[1][2] = bb[1][0] * Apre[0][2] + bb[1][1] * Apre[1][2] + bb[1][2] * Apre[2][2];
  Rot[2][0] = bb[2][0] * Apre[0][0] + bb[2][1] * Apre[1][0] + bb[2][2] * Apre[2][0];
  Rot[2][1] = bb[2][0] * Apre[0][1] + bb[2][1] * Apre[1][1] + bb[2][2] * Apre[2][1];
  Rot[2][2] = bb[2][0] * Apre[0][2] + bb[2][1] * Apre[1][2] + bb[2][2] * Apre[2][2];

  // deallocate temporals
  free (Apre);
  free (aa);
  free (bb);
}

void write_fdata (char *header, char *masses, double *bounds, int cell_flag, double *cell, int ntype, int natom, int i, int rank, int file_format)
{

  FILE *fi, *fo;
  char *fout, *line;
  double sq_c, cos_alpha, cos_gamma, sin_gamma, cos_beta, sin_beta;
  double a, b, c, A, B, C, xy, xz, yz;
  double alpha, beta, gamma;
  double xprd, yprd, xprdinv, yprdinv;

  alpha = beta = gamma = 0.0;
  i = i + 0;

  // Allocate string temporals
  line = (char *) scalloc (MAX_LINE, sizeof (char), "line");
  fout = (char *) scalloc (MAX_LINE, sizeof (char), "fout");

  sprintf (fout, "data.%s.%d", header, rank);
  if ((fo = fopen (fout, "w")) == NULL) {
    if (rank == 0)
      fprintf (stderr, "ERROR: Writing to %s! terminating...\n", fout);
    MPI_Abort (MPI_COMM_WORLD, 1);
  }
  sprintf (fout, "_tmp.dat.%d", rank);
  if ((fi = fopen (fout, "r")) == NULL) {
    if (rank == 0)
      printf ("ERROR: Cannot read from temporary file file %s!\n", fout);
    MPI_Abort (MPI_COMM_WORLD, 1);
  }

  fprintf (fo, "%s\n", header);
  fprintf (fo, "%d atoms\n%d atom types\n\n", natom, ntype);

  if (!cell_flag) {
    fprintf (fo, "%lf %lf xlo xhi\n%lf %lf ylo yhi\n%lf %lf zlo zhi\n\n", bounds[0] - 100, bounds[1] + 100, bounds[2] - 100, bounds[3] + 100,
	     bounds[4] - 100, bounds[5] + 100);
    //    fprintf (fo, "%lf %lf xlo xhi\n%lf %lf ylo yhi\n%lf %lf zlo zhi\n\n", bounds[0],
    //             bounds[1], bounds[2], bounds[3], bounds[4], bounds[5]);
  }
  else {
    if (file_format == 0) {	// bgf
      A = cell[0];
      B = cell[1];
      C = cell[2];
      alpha = cell[3];
      beta = cell[4];
      gamma = cell[5];

      sq_c = C * C;
      cos_alpha = cos (alpha * PI / 180.0);
      cos_gamma = cos (gamma * PI / 180.0);
      sin_gamma = sin (gamma * PI / 180.0);
      cos_beta = cos (beta * PI / 180.0);
      sin_beta = sin (beta * PI / 180.0);

      a = A;
      xy = B * cos_gamma;
      xz = C * cos_beta;
      b = B * sin_gamma;
      yz = C * (cos_alpha - cos_gamma * cos_beta) / sin_gamma;
      // Updated to match lammps' description (Byung-Hyung Kim, Uppsala University, 6/7/2016)
      c = sqrt (C * C - xz * xz - yz * yz);
    }
    else {			// xyz
      double *vnorma, *vnormb, *vnormc;
      double prod;

      vnorma = (double *) scalloc (3, sizeof (double), "normv");;
      vnormb = (double *) scalloc (3, sizeof (double), "normv");;
      vnormc = (double *) scalloc (3, sizeof (double), "normv");;

      A = sqrt (cell[0] * cell[0] + cell[1] * cell[1] + cell[2] * cell[2]);
      B = sqrt (cell[3] * cell[3] + cell[4] * cell[4] + cell[5] * cell[5]);
      C = sqrt (cell[6] * cell[6] + cell[7] * cell[7] + cell[8] * cell[8]);

      // alpha in degrees
      vnormb[0] = cell[3] / B;
      vnormb[1] = cell[4] / B;
      vnormb[2] = cell[5] / B;

      vnormc[0] = cell[6] / C;
      vnormc[1] = cell[7] / C;
      vnormc[2] = cell[8] / C;

      prod = dotProduct (vnormb, vnormc, 3);

      if (fabs (prod) < 1.)
	alpha = acos (prod) * 180.0 / PI;
      else {
	if (rank == 0)
	  printf ("Error: scalar product alpha %f\n", prod);
	MPI_Abort (MPI_COMM_WORLD, 1);
      }

      // beta in degrees
      vnorma[0] = cell[0] / A;
      vnorma[1] = cell[1] / A;
      vnorma[2] = cell[2] / A;

      vnormc[0] = cell[6] / C;
      vnormc[1] = cell[7] / C;
      vnormc[2] = cell[8] / C;

      prod = dotProduct (vnorma, vnormc, 3);

      if (fabs (prod) < 1.)
	beta = acos (prod) * 180.0 / PI;
      else {
	if (rank == 0)
	  printf ("Error: scalar product beta %f\n", prod);
	MPI_Abort (MPI_COMM_WORLD, 1);
      }

      // gamma in degrees
      vnorma[0] = cell[0] / A;
      vnorma[1] = cell[1] / A;
      vnorma[2] = cell[2] / A;

      vnormb[0] = cell[3] / B;
      vnormb[1] = cell[4] / B;
      vnormb[2] = cell[5] / B;

      prod = dotProduct (vnorma, vnormb, 3);

      if (fabs (prod) < 1.)
	gamma = acos (prod) * 180.0 / PI;
      else {
	if (rank == 0)
	  printf ("Error: scalar product gamma %f\n", prod);
	MPI_Abort (MPI_COMM_WORLD, 1);
      }
      sq_c = C * C;
      cos_alpha = cos (alpha * PI / 180.0);
      cos_gamma = cos (gamma * PI / 180.0);
      sin_gamma = sin (gamma * PI / 180.0);
      cos_beta = cos (beta * PI / 180.0);
      sin_beta = sin (beta * PI / 180.0);

      a = A;
      xy = B * cos_gamma;
      xz = C * cos_beta;
      b = B * sin_gamma;
      yz = C * (cos_alpha - cos_gamma * cos_beta) / sin_gamma;
      c = sqrt (C * C - xz * xz - yz * yz);

      free (vnorma);
      free (vnormb);
      free (vnormc);
    }

    // Write cell definition to lammps' data structure file
    fprintf (fo, "%lf %lf xlo xhi\n%lf %lf ylo yhi\n%lf %lf zlo zhi\n", 0.0, a, 0.0, b, 0.0, c);
    if (fabs (xy) > 1e-8 || fabs (xz) > 1e-8 || fabs (yz) > 1e-8) {
      xprd = a;
      yprd = b;
      xprdinv = 1.0 / xprd;
      yprdinv = 1.0 / yprd;
      // Adjust tilts factors when larger than half the box length, BK (2016-6-7)
      if (xy*xprdinv > 0.5) xy-=xprd;
      else if (xy*xprdinv < -0.5) xy+=xprd;
      if (xz*xprdinv > 0.5) xz-=xprd;
      else if (xz*xprdinv < -0.5) xz+=xprd;
      if (yz*yprdinv > 0.5) yz-=yprd;
      else if (yz*yprdinv < -0.5) yz-=yprd;
      fprintf (fo, "%lf %lf %lf xy xz yz\n", xy, xz, yz);     // iff triclinic
    }
  }

  fprintf (fo, "\nMasses\n\n%s\n", masses);
  while (fgets (line, MAX_LINE, fi) != NULL) {
    fputs (line, fo);
  }
  free (line);
  free (fout);
  fclose (fo);
  fclose (fi);
}

/****************************************************************************
*  Function     : get_atom_type_id
*  Description  : Uses atom type name to index into the force field data structure
                  (ffdata) that relates atom id with atom name (single body parameters in reaxFF)
*  Parameters   : Atom type LAMMPS id, reax/pqeq/morse/comb force field data structure (ffdata), and
                  MPI process rank
*  Effects      : Returns integer value corresponding to the LAMMPS atom type id
****************************************************************************/

int get_atom_type_id (char *type_id, ffieldtype * ffdata, int forcefield2optimize, int n, int rank)
{
  int i, id;

  id = -1;
  for (i = 0; i < n; i++) {
    if ((forcefield2optimize == REAX) || (forcefield2optimize == REAXC)) {
      if (strcmp (ffdata->rxff->sbp[i].name, type_id) == 0) {
	id = i;
	break;
      }
    }
    else if (forcefield2optimize == PQEQ) {
      if (strcmp (ffdata->pqeqff->pqeqatypes[i].label, type_id) == 0) {
        id = i;
        break;
      }
    }
    else if (forcefield2optimize == PPR) {
      if (strcmp (ffdata->pprff->ppratypes[i].label, type_id) == 0) {
        id = i;
        break;
      }
    }
    else if (forcefield2optimize == MORSE) {
      if (strcmp (ffdata->morseff->morseatypes[i].label, type_id) == 0) {
	id = i;
	break;
      }
    }
    else if (forcefield2optimize == COMB) {
      if (strcmp (ffdata->combff->atypes[i], type_id) == 0) {
	id = i;
	break;
      }
    }
    else if (forcefield2optimize == TERSOFF) {
      if (strcmp (ffdata->tersoffff->atypes[i], type_id) == 0) {
	id = i;
	break;
      }
    }
    else if (forcefield2optimize == ZHOU_EAM) {
      if (strcmp (ffdata->zhou_EAMff->atypes[i], type_id) == 0) {
	id = i;
	break;
      }
    }
    else if (forcefield2optimize == TERSOFF_MOD) {
      if (strcmp (ffdata->tersoff_modff->atypes[i], type_id) == 0) {
	id = i;
	break;
      }
    }
  }

  if (id == -1) {
    if (rank == 0)
      printf ("ERROR: Atom type %s not found in forcefield\n", type_id);
    MPI_Abort (MPI_COMM_WORLD, 1);
  }

  return id;
}
