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

#include "reaxc_ffield.h"
#include "mpi.h"

#ifdef DEBUG
# define DEBUG_PRINT(x) printf x
#else
# define DEBUG_PRINT(x)
#endif

int s (int a)
{
  return a ? a * s (a - 1) : 1;
}

int binomial (int n, int k)
{
  return s (n) / s (k) / s (n - k);
}

/****************************************************************************
*  Function     : Write_Force_Field_Reax
*  Description  : Writes contents of the ffid_pointer into a ffield file 
                  with format of type ReaxFF, associated with user-selected 
                  force field (-F reax)
*  Parameters   : Name of force field file, ffield parameters structure, 
                  forcefield array, and MPI process rank
*  Effects      : Writes a new force field file per evaluation pass, per MPI
                  rank
****************************************************************************/

void Write_Force_Field_Reax (char *ffield_file, reax_interaction * ffdata, double **ffid_pointer,
			     int rank)
{
  int i, j, k, l, token;
  FILE *fo;

  /* open force field for writing */
  if ((fo = fopen (ffield_file, "w")) == NULL) {
    fprintf (stderr, "ERROR: Cannot write to %s...\n", ffield_file);
    MPI_Abort (MPI_COMM_WORLD, 1);
  }

  token = 0;
  fprintf (fo, "Reactive MD-force field\n");
  fprintf (fo, "%3d       ! Number of general parameters\n", ffdata->gp.n_global);
  for (i = 0; i < ffdata->gp.n_global; i++)
    fprintf (fo, "%10.4f !Comment here\n", *ffid_pointer[token++]);

  //1-body data
  fprintf (fo, "%3d    !Nr of atoms; cov.r; valency;a.m;Rvdw;Evdw;gammaEEM;cov.r2;\n",
	   ffdata->num_atom_types);
  fprintf (fo, "            alfa;gammavdW;valency;Eunder;Eover;chiEEM;etaEEM;n.u.\n");
  fprintf (fo, "            cov r3;Elp;Heat inc.;n.u.;n.u.;n.u.;n.u.\n");
  fprintf (fo, "            ov/un;val1;n.u.;val3,vval4\n");
  for (i = 0; i < ffdata->num_atom_types; i++) {
    fprintf (fo, " %-2s", ffdata->sbp[i].name);

    fprintf (fo, "%9.4f", *ffid_pointer[token++]);
    fprintf (fo, "%9.4f", *ffid_pointer[token++]);
    fprintf (fo, "%9.4f", *ffid_pointer[token++]);
    fprintf (fo, "%9.4f", *ffid_pointer[token++]);
    fprintf (fo, "%9.4f", *ffid_pointer[token++]);
    fprintf (fo, "%9.4f", *ffid_pointer[token++]);
    fprintf (fo, "%9.4f", *ffid_pointer[token++]);
    fprintf (fo, "%9.4f\n", *ffid_pointer[token++]);

    fprintf (fo, "   %9.4f", *ffid_pointer[token++]);
    fprintf (fo, "%9.4f", *ffid_pointer[token++]);
    fprintf (fo, "%9.4f", *ffid_pointer[token++]);
    fprintf (fo, "%9.4f", *ffid_pointer[token++]);
    fprintf (fo, "%9.4f", *ffid_pointer[token++]);	// not used
    fprintf (fo, "%9.4f", *ffid_pointer[token++]);
    fprintf (fo, "%9.4f", *ffid_pointer[token++]);
    fprintf (fo, "%9.4f\n", *ffid_pointer[token++]);

    fprintf (fo, "   %9.4f", *ffid_pointer[token++]);
    fprintf (fo, "%9.4f", *ffid_pointer[token++]);
    fprintf (fo, "%9.4f", *ffid_pointer[token++]);	// not used
    fprintf (fo, "%9.4f", *ffid_pointer[token++]);
    fprintf (fo, "%9.4f", *ffid_pointer[token++]);
    fprintf (fo, "%9.4f", *ffid_pointer[token++]);
    fprintf (fo, "%9.4f", *ffid_pointer[token++]);	// not used
    fprintf (fo, "%9.4f\n", *ffid_pointer[token++]);	// not used

    fprintf (fo, "   %9.4f", *ffid_pointer[token++]);
    fprintf (fo, "%9.4f", *ffid_pointer[token++]);
    fprintf (fo, "%9.4f", *ffid_pointer[token++]);	// not used
    fprintf (fo, "%9.4f", *ffid_pointer[token++]);
    fprintf (fo, "%9.4f", *ffid_pointer[token++]);
    fprintf (fo, "%9.4f", *ffid_pointer[token++]);
    fprintf (fo, "%9.4f", *ffid_pointer[token++]);
    fprintf (fo, "%9.4f\n", *ffid_pointer[token++]);

    if (lg_reax) {
      fprintf (fo, "   %9.4f", *ffid_pointer[token++]);
      fprintf (fo, "%9.4f\n", *ffid_pointer[token++]);
    }

  }

  //2-body data
  fprintf (fo, "%3d      ! Nr of bonds; Edis1;LPpen;n.u.;pbe1;pbo5;13corr;pbo6\n",
	   ffdata->num_bond_types);
  fprintf (fo, "                         pbe2;pbo3;pbo4;n.u.;pbo1;pbo2;ovcorr\n");
  for (i = 0; i < ffdata->num_atom_types; i++) {
    for (j = 0; j < ffdata->num_atom_types; j++) {
      if (ffdata->tbp[i][j].used) {
	token = ffdata->tbp[i][j].ffid_ps;
	fprintf (fo, "%3d%3d", ffdata->tbp[i][j].i, ffdata->tbp[i][j].j);

	fprintf (fo, "%9.4f", *ffid_pointer[token++]);
	fprintf (fo, "%9.4f", *ffid_pointer[token++]);
	fprintf (fo, "%9.4f", *ffid_pointer[token++]);
	fprintf (fo, "%9.4f", *ffid_pointer[token++]);
	fprintf (fo, "%9.4f", *ffid_pointer[token++]);
	fprintf (fo, "%9.4f", *ffid_pointer[token++]);
	fprintf (fo, "%9.4f", *ffid_pointer[token++]);
	fprintf (fo, "%9.4f\n", *ffid_pointer[token++]);

	fprintf (fo, "      %9.4f", *ffid_pointer[token++]);
	fprintf (fo, "%9.4f", *ffid_pointer[token++]);
	fprintf (fo, "%9.4f", *ffid_pointer[token++]);
	fprintf (fo, "%9.4f", *ffid_pointer[token++]);
	fprintf (fo, "%9.4f", *ffid_pointer[token++]);
	fprintf (fo, "%9.4f", *ffid_pointer[token++]);
	fprintf (fo, "%9.4f", *ffid_pointer[token++]);
	fprintf (fo, "%9.4f\n", *ffid_pointer[token++]);

      }
    }
  }
  //off-diagonal
  fprintf (fo, "%3d    ! Nr of off-diagonal terms; Ediss;Ro;gamma;rsigma;rpi;rpi2\n",
	   ffdata->num_off_diag_types);
  for (i = 0; i < ffdata->num_atom_types; i++) {
    for (j = 0; j < ffdata->num_atom_types; j++) {
      if (ffdata->tbp[i][j].used_offdiag) {
	token = ffdata->tbp[i][j].od_ffid_ps;
	fprintf (fo, "%3d%3d", ffdata->tbp[i][j].i, ffdata->tbp[i][j].j);
	fprintf (fo, "%9.4f", *ffid_pointer[token++]);
	fprintf (fo, "%9.4f", *ffid_pointer[token++]);
	fprintf (fo, "%9.4f", *ffid_pointer[token++]);
	fprintf (fo, "%9.4f", *ffid_pointer[token++]);
	fprintf (fo, "%9.4f", *ffid_pointer[token++]);
	if (lg_reax)
	  fprintf (fo, "%9.4f", *ffid_pointer[token++]);
	fprintf (fo, "%9.4f\n", *ffid_pointer[token++]);
      }
    }
  }

  //3-body
  fprintf (fo, "%3d    ! Nr of angles;at1;at2;at3;Thetao,o;ka;kb;pv1;pv2\n",
	   ffdata->num_angle_types);
  for (i = 0; i < ffdata->num_atom_types; i++) {
    for (j = 0; j < ffdata->num_atom_types; j++) {
      for (k = 0; k < ffdata->num_atom_types; k++) {
	if (ffdata->thbp[i][j][k].used) {
	  token = ffdata->thbp[i][j][k].ffid_ps;
	  fprintf (fo, "%3d%3d%3d", ffdata->thbp[i][j][k].i, ffdata->thbp[i][j][k].j,
		   ffdata->thbp[i][j][k].k);

	  fprintf (fo, "%9.4f", *ffid_pointer[token++]);
	  fprintf (fo, "%9.4f", *ffid_pointer[token++]);
	  fprintf (fo, "%9.4f", *ffid_pointer[token++]);
	  fprintf (fo, "%9.4f", *ffid_pointer[token++]);
	  fprintf (fo, "%9.4f", *ffid_pointer[token++]);
	  fprintf (fo, "%9.4f", *ffid_pointer[token++]);
	  fprintf (fo, "%9.4f\n", *ffid_pointer[token++]);
	}
      }
    }
  }

  //4-body
// 4/15/2015: thanks to Byung-Hyun Kim (Uppsala University) for reporting a bug when
// the number of atoms is equal to the number of torsions in the force field file.
// The bug appeared when using the compact torsion representation, and it was fixed
// by looping over num_atom_types+1 for i,j,k, and l
  fprintf (fo, "%3d    ! Nr of torsions;at1;at2;at3;at4;;V1;V2;V3;V2(BO);vconj;n.u;n\n",
	   ffdata->num_torsion_types);
  for (i = 0; i < ffdata->num_atom_types+1; i++) {
    for (j = 0; j < ffdata->num_atom_types+1; j++) {
      for (k = 0; k < ffdata->num_atom_types+1; k++) {
	for (l = 0; l < ffdata->num_atom_types+1; l++) {
	  if (ffdata->fbp[i][j][k][l].used) {
            token = ffdata->fbp[i][j][k][l].ffid_ps;
	    fprintf (fo, "%3d%3d%3d%3d", i, j, k, l);

	    fprintf (fo, "%9.4f", *ffid_pointer[token++]);
	    fprintf (fo, "%9.4f", *ffid_pointer[token++]);
	    fprintf (fo, "%9.4f", *ffid_pointer[token++]);
	    fprintf (fo, "%9.4f", *ffid_pointer[token++]);
	    fprintf (fo, "%9.4f", *ffid_pointer[token++]);
	    fprintf (fo, "%9.4f", *ffid_pointer[token++]);
	    fprintf (fo, "%9.4f\n", *ffid_pointer[token++]);
	  }
	}
      }
    }
  }

  //h-bond
  fprintf (fo, "%3d    ! Nr of hydrogen bonds;at1;at2;at3;Rhb;Dehb;vhb1\n",
	   ffdata->num_hbond_types);
  for (i = 0; i < ffdata->num_atom_types; i++) {
    for (j = 0; j < ffdata->num_atom_types; j++) {
      for (k = 0; k < ffdata->num_atom_types; k++) {
	if (ffdata->hbp[i][j][k].used) {
	  token = ffdata->hbp[i][j][k].ffid_ps;
	  fprintf (fo, "%3d%3d%3d", ffdata->hbp[i][j][k].i, ffdata->hbp[i][j][k].j,
		   ffdata->hbp[i][j][k].k);

	  fprintf (fo, "%9.4f", *ffid_pointer[token++]);
	  fprintf (fo, "%9.4f", *ffid_pointer[token++]);
	  fprintf (fo, "%9.4f", *ffid_pointer[token++]);
	  fprintf (fo, "%9.4f\n", *ffid_pointer[token++]);

	}
      }
    }
  }
  fclose (fo);
}

/****************************************************************************
*  Function     : Read_Force_Field_Reax
*  Description  : Reads and parses contents of the user-provided ffield file
                  into the appropriate ReaxFF force field data structure format,
                  associated with user-selected force field (-F reax)
*  Parameters   : Name of force field file, ffield parameters structure,
                  forcefield array, and MPI process rank
*  Effects      : Reads force field parameters per GA evaluation pass, per MPI 
                  rank
****************************************************************************/

int Read_Force_Field_Reax (char *ffield_file, reax_interaction * reax, double ***ffid_pointer,
			   int rank)
{
  FILE *fi;
  char *s;
  char **strtmp;
  int a, token, i, j, k, m, n, cnt;
  double val;
  int sblines, sb_params, off_diag_params;
  char file[MAX_LINE];

#if (DEBUG)
  int idx;
#endif

  /* open force field file */
  sprintf (file, "%s", ffield_file);

  if ((fi = fopen (file, "r")) == NULL) {
    fprintf (stderr, "ERROR: opening the force field file! terminating...\n");
    MPI_Abort (MPI_COMM_WORLD, 1);
  }

  s = (char *) scalloc (MAX_LINE, sizeof (char), "s3");
  strtmp = (char **) scalloc (MAX_TOKENS, sizeof (char *), "tmp3");
  for (i = 0; i < MAX_TOKENS; i++)
    strtmp[i] = (char *) scalloc (MAX_TOKEN_LEN, sizeof (char), "tmp3i");

  /* reading first header comment */
  fgets (s, MAX_LINE, fi);
  /* line 2 is number of global parameters */
  fgets (s, MAX_LINE, fi);
  token = tokenize_string (s, &strtmp);
  /* reading the number of global parameters */
  n = atoi (strtmp[0]);
  if (n < 1) {
    fprintf (stderr, "WARNING: number of globals in ffield file is 0!\n");
    MPI_Abort (MPI_COMM_WORLD, 1);
  }
  reax->gp.n_global = n;
  if (initial_write)
    reax->gp.l = (double *) scalloc (n, sizeof (double), "reax_gp.l");
  // get general types, and skip
  for (i = 0; i < n; i++)
    fgets (s, MAX_LINE, fi);
  // get number of atom types, and skip
  fgets (s, MAX_LINE, fi);
  token = tokenize_string (s, &strtmp);
  reax->num_atom_types = atoi (strtmp[0]);
  // if low-gradient used, then atom entries have 5 lines, else 4
  // and off-diagonal parameters are 7, vs 6 if no lg
  if (lg_reax) {
    sblines = 5;
    sb_params = 34;
    off_diag_params = 7;
  }
  else {
    sblines = 4;
    sb_params = 32;
    off_diag_params = 6;
  }
  for (i = 0; i < reax->num_atom_types * sblines + 3; i++)
    fgets (s, MAX_LINE, fi);
  // get number of bond types, and skip
  fgets (s, MAX_LINE, fi);
  token = tokenize_string (s, &strtmp);
  reax->num_bond_types = atoi (strtmp[0]);
  for (i = 0; i < reax->num_bond_types * 2 + 1; i++)
    fgets (s, MAX_LINE, fi);
  // get number of off-diagonals types, and skip
  fgets (s, MAX_LINE, fi);
  token = tokenize_string (s, &strtmp);
  reax->num_off_diag_types = atoi (strtmp[0]);
  for (i = 0; i < reax->num_off_diag_types; i++)
    fgets (s, MAX_LINE, fi);
  // get number of angle types, and skip
  fgets (s, MAX_LINE, fi);
  token = tokenize_string (s, &strtmp);
  reax->num_angle_types = atoi (strtmp[0]);
  for (i = 0; i < reax->num_angle_types; i++)
    fgets (s, MAX_LINE, fi);
  // get number of torsion types, and skip
  fgets (s, MAX_LINE, fi);
  token = tokenize_string (s, &strtmp);
  reax->num_torsion_types = atoi (strtmp[0]);
  for (i = 0; i < reax->num_torsion_types; i++)
    fgets (s, MAX_LINE, fi);
  // get number of hbond types, and skip
  fgets (s, MAX_LINE, fi);
  token = tokenize_string (s, &strtmp);
  reax->num_hbond_types = atoi (strtmp[0]);
  for (i = 0; i < reax->num_hbond_types; i++)
    fgets (s, MAX_LINE, fi);

  if (rank == 0 && debug_level == 1)
    DEBUG_PRINT (("FFIELD: %d GENERAL, %d ATOM, %d BOND, %d OFF-DIAG, %d ANGLES, %d TORSION, and %d H-BONDS\n", reax->gp.n_global, reax->num_atom_types, reax->num_bond_types, reax->num_off_diag_types, reax->num_angle_types, reax->num_torsion_types, reax->num_hbond_types));

  // return to start, skip globals and comments for general and sbp
  rewind (fi);

  /* Allocating structures in reax_interaction */
  if (initial_write) {
    reax->sbp = (single_body_parameters *)
      scalloc (reax->num_atom_types, sizeof (single_body_parameters), "sbp");
    reax->tbp = (two_body_parameters **)
      scalloc (reax->num_atom_types, sizeof (two_body_parameters *), "tbp");
    reax->thbp = (three_body_header ***)
      scalloc (reax->num_atom_types, sizeof (three_body_header **), "thbp");
    reax->hbp = (hbond_parameters ***)
      scalloc (reax->num_atom_types, sizeof (hbond_parameters **), "hbp");
    reax->fbp = (four_body_header ****)
      scalloc (reax->num_atom_types+1, sizeof (four_body_header ***), "fbp");

    for (i = 0; i < reax->num_atom_types; i++) {
      reax->tbp[i] = (two_body_parameters *)
	scalloc (reax->num_atom_types, sizeof (two_body_parameters), "tbp[i]");
      reax->thbp[i] = (three_body_header **)
	scalloc (reax->num_atom_types, sizeof (three_body_header *), "thbp[i]");
      reax->hbp[i] = (hbond_parameters **)
	scalloc (reax->num_atom_types, sizeof (hbond_parameters *), "hbp[i]");

      for (j = 0; j < reax->num_atom_types; j++) {
	reax->thbp[i][j] = (three_body_header *)
	  scalloc (reax->num_atom_types, sizeof (three_body_header), "thbp[i,j]");
	reax->hbp[i][j] = (hbond_parameters *)
	  scalloc (reax->num_atom_types, sizeof (hbond_parameters), "hbp[i,j]");
      }
    }

    for (i = 0; i < reax->num_atom_types+1; i++) {
      reax->fbp[i] = (four_body_header ***)
        scalloc (reax->num_atom_types+1, sizeof (four_body_header **), "fbp[i]");

      for (j = 0; j < reax->num_atom_types+1; j++) {
        reax->fbp[i][j] = (four_body_header **)
          scalloc (reax->num_atom_types+1, sizeof (four_body_header *), "fbp[i,j]");

        for (k = 0; k < reax->num_atom_types+1; k++) {
          reax->fbp[i][j][k] = (four_body_header *)
            scalloc (reax->num_atom_types+1, sizeof (four_body_header), "fbp[i,j,k]");
        }
      }
    }
  }
  // Define number of entries for each ffield section
  gp_idx = reax->gp.n_global;
  sbp_idx = reax->num_atom_types * sb_params;
  tbp_idx = reax->num_bond_types * 16;
  tbodp_idx = reax->num_off_diag_types * off_diag_params;
  thbp_idx = reax->num_angle_types * 7;
  fbp_idx = reax->num_torsion_types * 7;
  hbp_idx = reax->num_hbond_types * 4;

  // vdWaals type: 1: Shielded Morse, no inner-wall 
  //               2: inner wall, no shielding  
  //               3: inner wall+shielding
  reax->gp.vdw_type = 0;

  a = gp_idx + sbp_idx + tbp_idx + tbodp_idx + thbp_idx + fbp_idx + hbp_idx;
  if (initial_write)
    *ffid_pointer = (double **) scalloc (a, sizeof (double *), "ffid_pointer");

  a = 0;

#if (DEBUG)
  if (rank == 0 && debug_level == 1) {
    if (reax->gp.n_global > 0)
      printf ("1: GENERAL PARAMETERS [%d]\n", gp_idx);
    else {
      printf ("NO GENERAL PARAMETERS FOUND IN FFIELD !\n");
      MPI_Abort (MPI_COMM_WORLD, 1);
    }
  }
#endif

  // reading general parameters
  /* reading first header comment */
  fgets (s, MAX_LINE, fi);
  /* line 2 is number of global parameters */
  fgets (s, MAX_LINE, fi);
  /* see reax_types.h for mapping between l[i] and the lambdas used in ff */
  for (i = 0; i < reax->gp.n_global; i++) {
    fgets (s, MAX_LINE, fi);
    token = tokenize_string (s, &strtmp);
    val = (double) atof (strtmp[0]);
    reax->gp.l[i] = val;
    (*ffid_pointer)[a++] = &reax->gp.l[i];
    if (rank == 0 && debug_level == 1)
      DEBUG_PRINT (("GENERAL: [%-2d]%7.3f\n", a, reax->gp.l[i]));
  }

  // skip atom header lines
  for (i = 0; i < 4; i++)
    fgets (s, MAX_LINE, fi);

  /* reading single atom parameters */
  /* there are 4 lines of each single atom parameters in ff files. these
     parameters later determine some of the pair and triplet parameters using
     combination rules. */

#if (DEBUG)
  if (rank == 0 && debug_level == 1) {
    if (reax->num_atom_types > 0)
      printf ("2: ATOM PARAMETERS [%d]\n", sbp_idx);
    else
      printf ("NO ATOM PARAMETERS FOUND IN FFIELD !\n");
  }
#endif

  for (i = 0; i < reax->num_atom_types; i++) {
    /* line one */
    fgets (s, MAX_LINE, fi);
    token = tokenize_string (s, &strtmp);
    for (j = 0; j < (int) strlen (strtmp[0]); ++j)
      reax->sbp[i].name[j] = strtmp[0][j];

    val = atof (strtmp[1]);
    reax->sbp[i].r_s = val;
    (*ffid_pointer)[a++] = &reax->sbp[i].r_s;
    val = atof (strtmp[2]);
    reax->sbp[i].valency = val;
    (*ffid_pointer)[a++] = &reax->sbp[i].valency;
    val = atof (strtmp[3]);
    reax->sbp[i].mass = val;
    (*ffid_pointer)[a++] = &reax->sbp[i].mass;
    val = atof (strtmp[4]);
    reax->sbp[i].r_vdw = val;
    (*ffid_pointer)[a++] = &reax->sbp[i].r_vdw;
    val = atof (strtmp[5]);
    reax->sbp[i].epsilon = val;
    (*ffid_pointer)[a++] = &reax->sbp[i].epsilon;
    val = atof (strtmp[6]);
    reax->sbp[i].gamma = val;
    (*ffid_pointer)[a++] = &reax->sbp[i].gamma;
    val = atof (strtmp[7]);
    reax->sbp[i].r_pi = val;
    (*ffid_pointer)[a++] = &reax->sbp[i].r_pi;
    val = atof (strtmp[8]);
    reax->sbp[i].valency_e = val;
    (*ffid_pointer)[a++] = &reax->sbp[i].valency_e;
    reax->sbp[i].nlp_opt = 0.5 * (reax->sbp[i].valency_e - reax->sbp[i].valency);

    /* line two */
    fgets (s, MAX_LINE, fi);
    token = tokenize_string (s, &strtmp);

    val = atof (strtmp[0]);
    reax->sbp[i].alpha = val;
    (*ffid_pointer)[a++] = &reax->sbp[i].alpha;
    val = atof (strtmp[1]);
    reax->sbp[i].gamma_w = val;
    (*ffid_pointer)[a++] = &reax->sbp[i].gamma_w;
    val = atof (strtmp[2]);
    reax->sbp[i].valency_boc = val;
    (*ffid_pointer)[a++] = &reax->sbp[i].valency_boc;
    val = atof (strtmp[3]);
    reax->sbp[i].p_ovun5 = val;
    (*ffid_pointer)[a++] = &reax->sbp[i].p_ovun5;
    val = atof (strtmp[4]);	// not used
    reax->sbp[i].p_ovun6 = val;
    (*ffid_pointer)[a++] = &reax->sbp[i].p_ovun6;
    val = atof (strtmp[5]);
    reax->sbp[i].chi = val;
    (*ffid_pointer)[a++] = &reax->sbp[i].chi;
    val = atof (strtmp[6]);
    reax->sbp[i].eta = val;
    (*ffid_pointer)[a++] = &reax->sbp[i].eta;
    val = atof (strtmp[7]);
    reax->sbp[i].p_hbond = (int) val;
    (*ffid_pointer)[a++] = (double *) &reax->sbp[i].p_hbond;
//    reax->sbp[i].p_hbond = val;
//    (*ffid_pointer)[a++] = &reax->sbp[i].p_hbond;

    /* line 3 */
    fgets (s, MAX_LINE, fi);
    token = tokenize_string (s, &strtmp);

    val = atof (strtmp[0]);
    reax->sbp[i].r_pi_pi = val;
    (*ffid_pointer)[a++] = &reax->sbp[i].r_pi_pi;
    val = atof (strtmp[1]);
    reax->sbp[i].p_lp2 = val;
    (*ffid_pointer)[a++] = &reax->sbp[i].p_lp2;
    val = atof (strtmp[2]);	// not used
    reax->sbp[i].heat = val;
    (*ffid_pointer)[a++] = &reax->sbp[i].heat;
    val = atof (strtmp[3]);
    reax->sbp[i].p_boc3 = val;
    (*ffid_pointer)[a++] = &reax->sbp[i].p_boc3;
    val = atof (strtmp[4]);
    reax->sbp[i].p_boc4 = val;
    (*ffid_pointer)[a++] = &reax->sbp[i].p_boc4;
    val = atof (strtmp[5]);
    reax->sbp[i].p_boc5 = val;
    (*ffid_pointer)[a++] = &reax->sbp[i].p_boc5;
    val = atof (strtmp[6]);	// not used
    reax->sbp[i].sbp_23 = val;
    (*ffid_pointer)[a++] = &reax->sbp[i].sbp_23;
    val = atof (strtmp[7]);	// not used
    reax->sbp[i].sbp_24 = val;
    (*ffid_pointer)[a++] = &reax->sbp[i].sbp_24;

    /* line 4  */
    fgets (s, MAX_LINE, fi);
    token = tokenize_string (s, &strtmp);

    val = atof (strtmp[0]);
    reax->sbp[i].p_ovun2 = val;
    (*ffid_pointer)[a++] = &reax->sbp[i].p_ovun2;
    val = atof (strtmp[1]);
    reax->sbp[i].p_val3 = val;
    (*ffid_pointer)[a++] = &reax->sbp[i].p_val3;
    val = atof (strtmp[2]);	// not used
    reax->sbp[i].sbp_27 = val;
    (*ffid_pointer)[a++] = &reax->sbp[i].sbp_27;
    val = atof (strtmp[3]);
    reax->sbp[i].valency_val = val;
    (*ffid_pointer)[a++] = &reax->sbp[i].valency_val;
    val = atof (strtmp[4]);
    reax->sbp[i].p_val5 = val;
    (*ffid_pointer)[a++] = &reax->sbp[i].p_val5;
    val = atof (strtmp[5]);
    reax->sbp[i].rcore2 = val;
    (*ffid_pointer)[a++] = &reax->sbp[i].rcore2;
    val = atof (strtmp[6]);
    reax->sbp[i].ecore2 = val;
    (*ffid_pointer)[a++] = &reax->sbp[i].ecore2;
    val = atof (strtmp[7]);
    reax->sbp[i].acore2 = val;
    (*ffid_pointer)[a++] = &reax->sbp[i].acore2;

    /* line 5 - iff sb_params */
    if (lg_reax) {
      fgets (s, MAX_LINE, fi);
      token = tokenize_string (s, &strtmp);

      val = atof (strtmp[0]);
      reax->sbp[i].lg1 = val;
      (*ffid_pointer)[a++] = &reax->sbp[i].lg1;
      val = atof (strtmp[1]);
      reax->sbp[i].lg2 = val;
      (*ffid_pointer)[a++] = &reax->sbp[i].lg2;
    }

#if (DEBUG)
    if (rank == 0 && debug_level == 1) {
      for (idx = 0; idx < sb_params; idx++) {
	if (idx == 0)
	  printf ("--%2d--%3s\n", i + 1, reax->sbp[i].name);
	else if ((idx - 1) % 8 == 0)
	  printf (" ");
	printf ("[%-4d]%7.3f ", idx + (i * sb_params) + gp_idx + 1,
		*(*ffid_pointer)[idx + (i * sb_params) + gp_idx]);
	if (idx == 7 || idx == 15 || idx == 23 || idx == 31 || idx == 33)
	  printf ("\n");
      }
    }
    printf ("\n");
#endif

    if (reax->sbp[i].rcore2 > 0.01 && reax->sbp[i].acore2 > 0.01) {	// Inner-wall
      if (reax->sbp[i].gamma_w > 0.5) {	// Shielding vdWaals
	if (reax->gp.vdw_type != 0 && reax->gp.vdw_type != 3)
	  fprintf (stderr, "WARNING: inconsistent vdWaals-parameters\n"
		   "Force field parameters for element %s\n"
		   "indicate inner wall+shielding, but earlier\n"
		   "atoms indicate different vdWaals-method.\n"
		   "This may cause division-by-zero errors.\n"
		   "Keeping vdWaals-setting for earlier atoms.\n", reax->sbp[i].name);
	else {
	  reax->gp.vdw_type = 3;
	  fprintf (stderr, "----> vdWaals type for element %s: Shielding+inner-wall\n",
		   reax->sbp[i].name);
	}
      }
      else {			// No shielding vdWaals parameters present
	if (reax->gp.vdw_type != 0 && reax->gp.vdw_type != 2)
	  fprintf (stderr, "WARNING: inconsistent vdWaals-parameters\n"
		   "Force field parameters for element %s\n"
		   "indicate inner wall without shielding, but earlier\n"
		   "atoms indicate different vdWaals-method.\n"
		   "This may cause division-by-zero errors.\n"
		   "Keeping vdWaals-setting for earlier atoms.\n", reax->sbp[i].name);
	else {
	  reax->gp.vdw_type = 2;
	  fprintf (stderr, "----> vdWaals type for element %s: No Shielding,inner-wall\n",
		   reax->sbp[i].name);
	}
      }
    }
    else {			// No Inner wall parameters present
      if (reax->sbp[i].gamma_w > 0.5) {	// Shielding vdWaals
	if (reax->gp.vdw_type != 0 && reax->gp.vdw_type != 1)
	  fprintf (stderr, "WARNING: inconsistent vdWaals-parameters\n"
		   "Force field parameters for element %s\n"
		   "indicate  shielding without inner wall, but earlier\n"
		   "atoms indicate different vdWaals-method.\n"
		   "This may cause division-by-zero errors.\n"
		   "Keeping vdWaals-setting for earlier atoms.\n", reax->sbp[i].name);
	else {
	  reax->gp.vdw_type = 1;
	  fprintf (stderr, "----> vdWaals type for element %s: Shielding,no inner-wall\n",
		   reax->sbp[i].name);
	}
      }
      else {
	fprintf (stderr, "ERROR: inconsistent vdWaals-parameters\n"
		 "No shielding or inner-wall set for element %s\n", reax->sbp[i].name);
      }
    }
  }

  if (rank == 0 && debug_level == 1)
    DEBUG_PRINT ((">> Done with ATOM interaction (%i types and %i parameters)\n",
		  reax->num_atom_types, a));
  fprintf (stderr, "----> vdWaals type: %d\n", reax->gp.vdw_type);

  /* Equate vval3 to valf for first-row elements (25/10/2004) */
  for (i = 0; i < reax->num_atom_types; i++)
    if (reax->sbp[i].mass < 21 && reax->sbp[i].valency_val != reax->sbp[i].valency_boc) {
      fprintf (stderr, "\nWARNING: changed valency_val to valency_boc for %s\n", reax->sbp[i].name);
      reax->sbp[i].valency_val = reax->sbp[i].valency_boc;
    }

  /* next line is number of two body combination and some comments */
  fgets (s, MAX_LINE, fi);
  /* a line of comments */
  fgets (s, MAX_LINE, fi);

#if (DEBUG)
  if (rank == 0 && debug_level == 1) {
    if (reax->num_bond_types > 0)
      printf ("2: 2-BODY PARAMETERS [%d]\n", tbp_idx);
    else
      printf ("NO 2-BODY PARAMETERS FOUND IN FFIELD !\n");
  }
#endif

  for (i = 0; i < reax->num_bond_types; i++) {
    /* line 1 */
    fgets (s, MAX_LINE, fi);
    token = tokenize_string (s, &strtmp);

    j = atoi (strtmp[0]) - 1;
    k = atoi (strtmp[1]) - 1;

    //if (j < reax->num_atom_types && k < reax->num_atom_types) {
    reax->tbp[j][k].used = 1;
    reax->tbp[j][k].ffid_ps = a;
    reax->tbp[j][k].i = j + 1;
    reax->tbp[j][k].j = k + 1;

    val = atof (strtmp[2]);
    reax->tbp[j][k].De_s = val;
    reax->tbp[k][j].De_s = val;
    (*ffid_pointer)[a++] = &reax->tbp[j][k].De_s;
    val = atof (strtmp[3]);
    reax->tbp[j][k].De_p = val;
    reax->tbp[k][j].De_p = val;
    (*ffid_pointer)[a++] = &reax->tbp[j][k].De_p;
    val = atof (strtmp[4]);
    reax->tbp[j][k].De_pp = val;
    reax->tbp[k][j].De_pp = val;
    (*ffid_pointer)[a++] = &reax->tbp[j][k].De_pp;
    val = atof (strtmp[5]);
    reax->tbp[j][k].p_be1 = val;
    reax->tbp[k][j].p_be1 = val;
    (*ffid_pointer)[a++] = &reax->tbp[j][k].p_be1;
    val = atof (strtmp[6]);
    reax->tbp[j][k].p_bo5 = val;
    reax->tbp[k][j].p_bo5 = val;
    (*ffid_pointer)[a++] = &reax->tbp[j][k].p_bo5;
    val = atof (strtmp[7]);
    reax->tbp[j][k].v13cor = val;
    reax->tbp[k][j].v13cor = val;
    (*ffid_pointer)[a++] = &reax->tbp[j][k].v13cor;
    val = atof (strtmp[8]);
    reax->tbp[j][k].p_bo6 = val;
    reax->tbp[k][j].p_bo6 = val;
    (*ffid_pointer)[a++] = &reax->tbp[j][k].p_bo6;
    val = atof (strtmp[9]);
    reax->tbp[j][k].p_ovun1 = val;
    reax->tbp[k][j].p_ovun1 = val;
    (*ffid_pointer)[a++] = &reax->tbp[j][k].p_ovun1;

    /* line 2 */
    fgets (s, MAX_LINE, fi);
    token = tokenize_string (s, &strtmp);

    val = atof (strtmp[0]);
    reax->tbp[j][k].p_be2 = val;
    reax->tbp[k][j].p_be2 = val;
    (*ffid_pointer)[a++] = &reax->tbp[j][k].p_be2;
    val = atof (strtmp[1]);
    reax->tbp[j][k].p_bo3 = val;
    reax->tbp[k][j].p_bo3 = val;
    (*ffid_pointer)[a++] = &reax->tbp[j][k].p_bo3;
    val = atof (strtmp[2]);
    reax->tbp[j][k].p_bo4 = val;
    reax->tbp[k][j].p_bo4 = val;
    (*ffid_pointer)[a++] = &reax->tbp[j][k].p_bo4;
    val = atof (strtmp[3]);
    reax->tbp[j][k].tbp_12 = val;
    reax->tbp[k][j].tbp_12 = val;
    (*ffid_pointer)[a++] = &reax->tbp[j][k].tbp_12;
    val = atof (strtmp[4]);
    reax->tbp[j][k].p_bo1 = val;
    reax->tbp[k][j].p_bo1 = val;
    (*ffid_pointer)[a++] = &reax->tbp[j][k].p_bo1;
    val = atof (strtmp[5]);
    reax->tbp[j][k].p_bo2 = val;
    reax->tbp[k][j].p_bo2 = val;
    (*ffid_pointer)[a++] = &reax->tbp[j][k].p_bo2;
    val = atof (strtmp[6]);
    reax->tbp[j][k].tbp_15 = val;
    reax->tbp[k][j].tbp_15 = val;
    (*ffid_pointer)[a++] = &reax->tbp[j][k].tbp_15;
    val = atof (strtmp[7]);
    reax->tbp[j][k].tbp_16 = val;
    reax->tbp[k][j].tbp_16 = val;
    (*ffid_pointer)[a++] = &reax->tbp[j][k].tbp_16;

#if (DEBUG)
    if (rank == 0 && debug_level == 1) {
      for (idx = 0; idx < 16; idx++) {
	if (idx == 0)
	  printf ("%3d: %3d-%-2d\n", i + 1, reax->tbp[j][k].i, reax->tbp[j][k].j);
	else if ((idx - 1) % 8 == 0)
	  printf ("  ");
	printf ("[%4d]%7.3f ", idx + (i * 16) + gp_idx + sbp_idx + 1,
		*(*ffid_pointer)[idx + (i * 16) + gp_idx + sbp_idx]);
	if (idx == 7 || idx == 15)
	  printf ("\n");
      }
    }
#endif
  }

  if (rank == 0 && debug_level == 1)
    DEBUG_PRINT ((">> Done with 2-BODY interactions [%i types and %i parameters]\n",
		  reax->num_bond_types, a));

  /* calculating combination rules and filling up remaining fields. */
  /*
     for (i=0; i < reax->num_atom_types; i++)
     for (j=i; j < reax->num_atom_types; j++) {

     reax->tbp[i][j].r_s = 0.5 * (reax->sbp[i].r_s + reax->sbp[j].r_s);
     reax->tbp[j][i].r_s = 0.5 * (reax->sbp[j].r_s + reax->sbp[i].r_s);

     reax->tbp[i][j].r_p = 0.5 * (reax->sbp[i].r_pi + reax->sbp[j].r_pi);
     reax->tbp[j][i].r_p = 0.5 * (reax->sbp[j].r_pi + reax->sbp[i].r_pi);

     reax->tbp[i][j].r_pp = 0.5 * (reax->sbp[i].r_pi_pi + reax->sbp[j].r_pi_pi);
     reax->tbp[j][i].r_pp = 0.5 * (reax->sbp[j].r_pi_pi + reax->sbp[i].r_pi_pi);

     reax->tbp[i][j].p_boc3 = sqrt(reax->sbp[i].b_o_132 * reax->sbp[j].b_o_132);
     reax->tbp[j][i].p_boc3 = sqrt(reax->sbp[j].b_o_132 * reax->sbp[i].b_o_132);

     reax->tbp[i][j].p_boc4 = sqrt(reax->sbp[i].b_o_131 * reax->sbp[j].b_o_131);
     reax->tbp[j][i].p_boc4 = sqrt(reax->sbp[j].b_o_131 * reax->sbp[i].b_o_131);

     reax->tbp[i][j].p_boc5 = sqrt(reax->sbp[i].b_o_133 * reax->sbp[j].b_o_133);
     reax->tbp[j][i].p_boc5 = sqrt(reax->sbp[j].b_o_133 * reax->sbp[i].b_o_133);

     reax->tbp[i][j].D = sqrt(reax->sbp[i].epsilon * reax->sbp[j].epsilon);
     reax->tbp[j][i].D = sqrt(reax->sbp[j].epsilon * reax->sbp[i].epsilon);

     reax->tbp[i][j].alpha = sqrt(reax->sbp[i].alpha * reax->sbp[j].alpha);
     reax->tbp[j][i].alpha = sqrt(reax->sbp[j].alpha * reax->sbp[i].alpha);

     reax->tbp[i][j].r_vdW = 2.0 * sqrt(reax->sbp[i].r_vdw * reax->sbp[j].r_vdw);
     reax->tbp[j][i].r_vdW = 2.0 * sqrt(reax->sbp[j].r_vdw * reax->sbp[i].r_vdw);

     reax->tbp[i][j].gamma_w = sqrt(reax->sbp[i].gamma_w * reax->sbp[j].gamma_w);
     reax->tbp[j][i].gamma_w = sqrt(reax->sbp[j].gamma_w * reax->sbp[i].gamma_w);

     reax->tbp[i][j].gamma = pow(reax->sbp[i].gamma * reax->sbp[j].gamma,-1.5);
     reax->tbp[j][i].gamma = pow(reax->sbp[j].gamma * reax->sbp[i].gamma,-1.5);

     // additions for additional vdWaals interaction types - inner core
     reax->tbp[i][j].rcore = reax->tbp[j][i].rcore = sqrt( reax->sbp[i].rcore2 * reax->sbp[j].rcore2 );
     reax->tbp[i][j].ecore = reax->tbp[j][i].ecore = sqrt( reax->sbp[i].ecore2 * reax->sbp[j].ecore2 );
     reax->tbp[i][j].acore = reax->tbp[j][i].acore = sqrt( reax->sbp[i].acore2 * reax->sbp[j].acore2 );
     }

   */
  /* next line is number of two body offdiagonal combinations and comments */
  /* these are two body offdiagonal terms that are different from the
     combination rules defined above */

  fgets (s, MAX_LINE, fi);

#if (DEBUG)
  if (rank == 0 && debug_level == 1) {
    if (reax->num_off_diag_types > 0)
      printf ("3: 2-BODY OFF-DIAGONAL PARAMETERS [%d]\n", tbodp_idx);
    else
      printf ("NO 2-BODY OFF-DIAGONAL PARAMETERS FOUND IN FFIELD !\n");
  }
#endif

  for (i = 0; i < reax->num_off_diag_types; i++) {
    fgets (s, MAX_LINE, fi);
    token = tokenize_string (s, &strtmp);
    j = atoi (strtmp[0]) - 1;
    k = atoi (strtmp[1]) - 1;

    if (j < reax->num_atom_types && k < reax->num_atom_types) {
      reax->tbp[j][k].used_offdiag = 1;

      reax->tbp[j][k].od_ffid_ps = a;
      reax->tbp[j][k].i = j + 1;
      reax->tbp[j][k].j = k + 1;

      val = atof (strtmp[2]);
      //if (val > 0.0) {
      reax->tbp[j][k].D = val;
      reax->tbp[k][j].D = val;
      (*ffid_pointer)[a++] = &reax->tbp[j][k].D;
      //}

      val = atof (strtmp[3]);
      //if (val > 0.0) {
      reax->tbp[j][k].r_vdW = val;
      reax->tbp[k][j].r_vdW = val;
      (*ffid_pointer)[a++] = &reax->tbp[j][k].r_vdW;
      //}

      val = atof (strtmp[4]);
      //if (val > 0.0) {
      reax->tbp[j][k].alpha = val;
      reax->tbp[k][j].alpha = val;
      (*ffid_pointer)[a++] = &reax->tbp[j][k].alpha;
      //}

      val = atof (strtmp[5]);
      //if (val > 0.0) {
      reax->tbp[j][k].r_s = val;
      reax->tbp[k][j].r_s = val;
      (*ffid_pointer)[a++] = &reax->tbp[j][k].r_s;
      //}

      val = atof (strtmp[6]);
      //if (val > 0.0) {
      reax->tbp[j][k].r_p = val;
      reax->tbp[k][j].r_p = val;
      (*ffid_pointer)[a++] = &reax->tbp[j][k].r_p;
      //}

      val = atof (strtmp[7]);
      //if (val > 0.0) {
      reax->tbp[j][k].r_pp = val;
      reax->tbp[k][j].r_pp = val;
      (*ffid_pointer)[a++] = &reax->tbp[j][k].r_pp;
      //}

      if (lg_reax) {
	val = atof (strtmp[8]);
	//if (val > 0.0) {
	reax->tbp[j][k].lg_cij = val;
	reax->tbp[k][j].lg_cij = val;
	(*ffid_pointer)[a++] = &reax->tbp[j][k].lg_cij;
      }
    }
    else {
      printf ("Atom not defined in Off-diagonal entry\n");
      MPI_Abort (MPI_COMM_WORLD, 1);
    }

#if (DEBUG)
    if (rank == 0 && debug_level == 1) {
      for (idx = 0; idx < off_diag_params; idx++) {
	if (idx == 0)
	  printf ("%3d: %3d-%-2d\n", i + 1, reax->tbp[j][k].i, reax->tbp[j][k].j);
	else if ((idx - 1) % 8 == 0)
	  printf ("  ");
	printf ("[%4d]%7.3f ", idx + (i * off_diag_params) + gp_idx + sbp_idx + tbp_idx + 1,
		*(*ffid_pointer)[idx + (i * off_diag_params) + gp_idx + sbp_idx + tbp_idx]);
	if (idx == off_diag_params - 1)
	  printf ("\n");
      }
    }
#endif

  }

  if (rank == 0 && debug_level == 1)
    DEBUG_PRINT ((">> Done with Off-diag interactions [%i types and %i parameters]\n",
		  reax->num_off_diag_types, a));

  /* 3-body parameters - 
     supports multi-well potentials (upto MAX_3BODY_PARAM in mytypes.h) */
  /* clear entries first */
  for (i = 0; i < reax->num_atom_types; ++i)
    for (j = 0; j < reax->num_atom_types; ++j)
      for (k = 0; k < reax->num_atom_types; ++k)
	reax->thbp[i][j][k].cnt = 0;

  /* next line is number of 3-body params and some comments */
  fgets (s, MAX_LINE, fi);

#if (DEBUG)
  if (rank == 0 && debug_level == 1) {
    if (reax->num_angle_types > 0)
      printf ("4: 3-BODY PARAMETERS [%d]\n", thbp_idx);
    else
      printf ("NO 3-BODY PARAMETERS FOUND IN FFIELD !\n");
  }
#endif

  for (i = 0; i < reax->num_angle_types; i++) {
    fgets (s, MAX_LINE, fi);
    token = tokenize_string (s, &strtmp);

    j = atoi (strtmp[0]) - 1;
    k = atoi (strtmp[1]) - 1;
    m = atoi (strtmp[2]) - 1;

    //if (j < reax->num_atom_types && k < reax->num_atom_types && 
    //m < reax->num_atom_types) {
    reax->thbp[j][k][m].used = 1;
    reax->thbp[j][k][m].ffid_ps = a;
    reax->thbp[j][k][m].i = j + 1;
    reax->thbp[j][k][m].j = k + 1;
    reax->thbp[j][k][m].k = m + 1;
    cnt = reax->thbp[j][k][m].cnt;
    reax->thbp[j][k][m].cnt++;
    reax->thbp[m][k][j].cnt++;

    val = atof (strtmp[3]);
    reax->thbp[j][k][m].prm[cnt].theta_00 = val;
    reax->thbp[m][k][j].prm[cnt].theta_00 = val;
    (*ffid_pointer)[a++] = &reax->thbp[j][k][m].prm[cnt].theta_00;

    val = atof (strtmp[4]);
    reax->thbp[j][k][m].prm[cnt].p_val1 = val;
    reax->thbp[m][k][j].prm[cnt].p_val1 = val;
    (*ffid_pointer)[a++] = &reax->thbp[j][k][m].prm[cnt].p_val1;

    val = atof (strtmp[5]);
    reax->thbp[j][k][m].prm[cnt].p_val2 = val;
    reax->thbp[m][k][j].prm[cnt].p_val2 = val;
    (*ffid_pointer)[a++] = &reax->thbp[j][k][m].prm[cnt].p_val2;

    val = atof (strtmp[6]);
    reax->thbp[j][k][m].prm[cnt].p_coa1 = val;
    reax->thbp[m][k][j].prm[cnt].p_coa1 = val;
    (*ffid_pointer)[a++] = &reax->thbp[j][k][m].prm[cnt].p_coa1;

    val = atof (strtmp[7]);
    reax->thbp[j][k][m].prm[cnt].p_val7 = val;
    reax->thbp[m][k][j].prm[cnt].p_val7 = val;
    (*ffid_pointer)[a++] = &reax->thbp[j][k][m].prm[cnt].p_val7;

    val = atof (strtmp[8]);
    reax->thbp[j][k][m].prm[cnt].p_pen1 = val;
    reax->thbp[m][k][j].prm[cnt].p_pen1 = val;
    (*ffid_pointer)[a++] = &reax->thbp[j][k][m].prm[cnt].p_pen1;

    val = atof (strtmp[9]);
    reax->thbp[j][k][m].prm[cnt].p_val4 = val;
    reax->thbp[m][k][j].prm[cnt].p_val4 = val;
    (*ffid_pointer)[a++] = &reax->thbp[j][k][m].prm[cnt].p_val4;

#if (DEBUG)
    if (rank == 0 && debug_level == 1) {
      for (idx = 0; idx < 7; idx++) {
	if (idx == 0)
	  printf ("%3d: %3d %2d %2d\n", i + 1, reax->thbp[j][k][m].i,
		  reax->thbp[j][k][m].j, reax->thbp[j][k][m].k);
	printf ("[%4d]%7.3f ", idx + (i * 7) + gp_idx + sbp_idx + tbp_idx + tbodp_idx + 1,
		*(*ffid_pointer)[idx + (i * 7) + gp_idx + sbp_idx + tbp_idx + tbodp_idx]);
	if (idx == 6)
	  printf ("\n");
      }
    }
#endif
  }

  if (rank == 0 && debug_level == 1)
    DEBUG_PRINT ((">> Done with three body interaction [%i types and %i parameters]\n",
		  reax->num_angle_types, a));

  /* 4-body parameters are entered in compact form. i.e. 0-X-Y-0
     correspond to any type of pair of atoms in 1 and 4
     position. However, explicit X-Y-Z-W takes precedence over the
     default description.
     supports multi-well potentials (upto MAX_4BODY_PARAM in mytypes.h)
     IMPORTANT: for now, directions on how to read multi-entries from ffield 
     is not clear */

  /* next line is number of 4-body params and some comments */
  fgets (s, MAX_LINE, fi);

#if (DEBUG)
  if (rank == 0 && debug_level == 1) {
    if (reax->num_torsion_types > 0)
      printf ("5: 4-BODY PARAMETERS [%d]\n", fbp_idx);
    else
      printf ("NO 4-BODY PARAMETERS FOUND IN FFIELD !\n");
  }
#endif

  for (i = 0; i < reax->num_torsion_types; i++) {
    fgets (s, MAX_LINE, fi);
    token = tokenize_string (s, &strtmp);

    j = atoi (strtmp[0]);
    k = atoi (strtmp[1]);
    m = atoi (strtmp[2]);
    n = atoi (strtmp[3]);

    reax->fbp[j][k][m][n].used = 1;
    reax->fbp[j][k][m][n].ffid_ps = a;

    val = atof (strtmp[4]);
    reax->fbp[j][k][m][n].prm[0].V1 = val;
    reax->fbp[n][m][k][j].prm[0].V1 = val;
    (*ffid_pointer)[a++] = &reax->fbp[j][k][m][n].prm[0].V1;

    val = atof (strtmp[5]);
    reax->fbp[j][k][m][n].prm[0].V2 = val;
    reax->fbp[n][m][k][j].prm[0].V2 = val;
    (*ffid_pointer)[a++] = &reax->fbp[j][k][m][n].prm[0].V2;

    val = atof (strtmp[6]);
    reax->fbp[j][k][m][n].prm[0].V3 = val;
    reax->fbp[n][m][k][j].prm[0].V3 = val;
    (*ffid_pointer)[a++] = &reax->fbp[j][k][m][n].prm[0].V3;

    val = atof (strtmp[7]);
    reax->fbp[j][k][m][n].prm[0].p_tor1 = val;
    reax->fbp[n][m][k][j].prm[0].p_tor1 = val;
    (*ffid_pointer)[a++] = &reax->fbp[j][k][m][n].prm[0].p_tor1;

    val = atof (strtmp[8]);
    reax->fbp[j][k][m][n].prm[0].p_cot1 = val;
    reax->fbp[n][m][k][j].prm[0].p_cot1 = val;
    (*ffid_pointer)[a++] = &reax->fbp[j][k][m][n].prm[0].p_cot1;

    val = atof (strtmp[9]);
    reax->fbp[j][k][m][n].prm[0].fbp_6 = val;
    reax->fbp[n][m][k][j].prm[0].fbp_6 = val;
    (*ffid_pointer)[a++] = &reax->fbp[j][k][m][n].prm[0].fbp_6;

    val = atof (strtmp[10]);
    reax->fbp[j][k][m][n].prm[0].fbp_7 = val;
    reax->fbp[n][m][k][j].prm[0].fbp_7 = val;
    (*ffid_pointer)[a++] = &reax->fbp[j][k][m][n].prm[0].fbp_7;

#if (DEBUG)
    if (rank == 0 && debug_level == 1) {
      for (idx = 0; idx < 7; idx++) {
	if (idx == 0)
	  printf ("%3d: %3d %2d %2d %2d: ", i + 1, j, k, m, n);
	printf ("[%4d]%7.3f ",
		idx + (i * 7) + gp_idx + sbp_idx + tbp_idx + tbodp_idx + thbp_idx + 1,
		*(*ffid_pointer)[idx + (i * 7) + gp_idx + sbp_idx + tbp_idx + tbodp_idx +
				 thbp_idx]);
	if (idx == 6)
	  printf ("\n");
      }
    }
#endif
  }

  if (rank == 0 && debug_level == 1)
    DEBUG_PRINT ((">> Done with 4-BODY interactions [%i types and %i parameters]\n",
		  reax->num_torsion_types, a));

  /* next line is number of hydrogen bond params and some comments */
  fgets (s, MAX_LINE, fi);

#if (DEBUG)
  if (rank == 0 && debug_level == 1) {
    if (reax->num_hbond_types > 0)
      printf ("6: H-BOND PARAMETERS [%d]\n", hbp_idx);
    else
      printf ("NO H-BOND PARAMETERS FOUND IN FFIELD !\n");
  }
#endif

  for (i = 0; i < reax->num_hbond_types; i++) {
    fgets (s, MAX_LINE, fi);
    token = tokenize_string (s, &strtmp);

    j = atoi (strtmp[0]) - 1;
    k = atoi (strtmp[1]) - 1;
    m = atoi (strtmp[2]) - 1;

    //if( j < reax->num_atom_types && m < reax->num_atom_types ) {
    reax->hbp[j][k][m].used = 1;
    reax->hbp[j][k][m].ffid_ps = a;
    reax->hbp[j][k][m].i = j + 1;
    reax->hbp[j][k][m].j = k + 1;
    reax->hbp[j][k][m].k = m + 1;
    val = atof (strtmp[3]);
    reax->hbp[j][k][m].r0_hb = val;
    (*ffid_pointer)[a++] = &reax->hbp[j][k][m].r0_hb;

    val = atof (strtmp[4]);
    reax->hbp[j][k][m].p_hb1 = val;
    (*ffid_pointer)[a++] = &reax->hbp[j][k][m].p_hb1;

    val = atof (strtmp[5]);
    reax->hbp[j][k][m].p_hb2 = val;
    (*ffid_pointer)[a++] = &reax->hbp[j][k][m].p_hb2;

    val = atof (strtmp[6]);
    reax->hbp[j][k][m].p_hb3 = val;
    (*ffid_pointer)[a++] = &reax->hbp[j][k][m].p_hb3;

#if (DEBUG)
    if (rank == 0 && debug_level == 1) {
      for (idx = 0; idx < 4; idx++) {
	if (idx == 0)
	  printf ("%3d: %3d %2d %2d: ", i + 1, reax->hbp[j][k][m].i, reax->hbp[j][k][m].j,
		  reax->hbp[j][k][m].k);
	printf ("[%4d]%7.3f ",
		idx + (i * 4) + gp_idx + sbp_idx + tbp_idx + tbodp_idx + thbp_idx + fbp_idx + 1,
		*(*ffid_pointer)[idx + (i * 4) + gp_idx + sbp_idx + tbp_idx + tbodp_idx + thbp_idx +
				 fbp_idx]);
	if (idx == 3)
	  printf ("\n");
      }
    }
#endif
  }

  if (rank == 0 && debug_level == 1)
    DEBUG_PRINT ((">> Done with H-BOND interaction [%i types and %i parameters]\n",
		  reax->num_hbond_types, a));

  // close file
  fclose (fi);

  /* deallocate helper storage */
  for (i = 0; i < MAX_TOKENS; i++)
    free (strtmp[i]);
  free (strtmp);
  free (s);

  if (rank == 0 && debug_level == 1)
    DEBUG_PRINT ((">> Finished reading force field\n"));

  return a;
}
