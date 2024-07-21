/* ----------------------------------------------------------------------
   GARFFIELD: Genetic-Algorithm Reactive Force FIELD parameter optimizer.

   Copyright (2012-2014) California Institute of Technology
   Andres Jaramillo-Botero (ajaramil@caltech.edu)
   http://www.wag.caltech.edu/home/ajaramil/GARFfield.html
------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

//#include "reaxc_types.h"
#include "tool_box.h"
#include "trainset.h"
#include "mpi.h"

#ifdef DEBUG
# define DEBUG_PRINT(x) printf x
#else
# define DEBUG_PRINT(x) 
#endif

//#define num_objectives 12

/****************************************************************************
*  Function     : parse_tset
*  Description  : Reads and parses the user-defined training set (tset_file) 
                  and populates the data structure tset. It also calculates
                  the number of training cases per fitness objective (i.e.
                  charges, geometry, energies, etc.).
*  Parameters   : User-defined training set file (tset_file), training set
                  data structure, number of training cases per fitness objective
                  (nfit), and MPI process rank
*  Effects      : Populates the training set data structure tset, and the
                  nfit array of number of training cases per fitness objective
                  (nfit)
****************************************************************************/

void parse_tset (char *tset_file, struct tset_data *tset, int **nfit, int rank)
{

  FILE *fp;
  int i, j, c, n, valid;
  char **strtmp;
  char file[MAX_LINE];
  char *line;

  line = (char*) scalloc(MAX_LINE, sizeof(char), "line");
  strtmp = (char **) scalloc (MAX_TOKENS, sizeof (char *), "strtmp");
  for (i = 0; i < MAX_TOKENS; i++)
    strtmp[i] = (char *) scalloc (MAX_TOKEN_LEN, sizeof (char), "tmpi");

  //check for charge
  i = get_num_entries (tset_file, "CHARGE", "ENDCHARGE", 3, rank);
  (*nfit)[CHARGE] = i;
  if (!i) {
    if (rank == 0) printf ("WARNING: No valid charge training set data read from %s!\n", tset_file);}
  else {
    if (rank == 0 && debug_level == 1) DEBUG_PRINT ((">> Found %i charge records in %s\n", i, tset_file));
    tset->charge = (charge_tset_item *) scalloc (i, sizeof (charge_tset_item), "tset_charge");
  }

  //check for cell parameters
  i = get_num_entries (tset_file, "CELL PARAMETERS", "ENDCELL PARAMETERS", 3, rank);
  (*nfit)[CELL] = i;
  if (!i) {
    if (rank == 0) printf ("WARNING: No valid cell training set data read from %s!\n", tset_file);}
  else {
    if (rank == 0 && debug_level == 1) DEBUG_PRINT ((">> Found %i cell records in %s\n", i, tset_file));
    tset->cell = (cell_tset_item *) scalloc (i, sizeof (cell_tset_item), "tset_cell");
  }

  //check for stress parameters
  i = get_num_entries (tset_file, "STRESS", "ENDSTRESS", 2, rank);
  (*nfit)[STRESS] = i;
  if (!i) { 
    if (rank == 0) printf ("WARNING: No valid cell stress training set data read from %s!\n", tset_file);}
  else {
    if (rank == 0 && debug_level == 1) DEBUG_PRINT ((">> Found %i cell stress records in %s\n", i, tset_file));
    tset->stress = (stress_tset_item *) scalloc (i, sizeof (stress_tset_item), "tset_stress");
  }

  //check for cell pressure
  i = get_num_entries (tset_file, "PRESSURE", "ENDPRESS", 3, rank);
  (*nfit)[PRESSURE] = i;
  if (!i) {
    if (rank == 0) printf ("WARNING: No valid cell pressure training set data read from %s!\n", tset_file);}
  else {
    if (rank == 0 && debug_level == 1) DEBUG_PRINT ((">> Found %i cell pressure records in %s\n", i, tset_file));
    tset->press = (press_tset_item *) scalloc (i, sizeof (press_tset_item), "tset_pressure");
  }

  //check for frequencies
  //not yet implemented

  //check for heat of formations
  //not yet implemented

  //check for geometry, changed min parameters to allow calc_initial_geom
  i = get_num_entries (tset_file, "GEOMETRY", "ENDGEOMETRY", 3, rank);
  (*nfit)[GEOMETRY] = i;
  if (!i) {
    if (rank == 0) printf ("WARNING: No valid geometry training set data read from %s!\n", tset_file);}
  else {
    if (rank == 0 && debug_level == 1) DEBUG_PRINT ((">> Found %i geometry records in %s\n", i, tset_file));
    tset->geom = (geom_tset_item *) scalloc (i, sizeof (geom_tset_item), "tset_geom");
  }

  //check for structure (RDFs)
  i = get_num_entries (tset_file, "STRUCTURE", "ENDSTRUCTURE", 4, rank);
  (*nfit)[STRUCTURE] = i;
  if (!i) {
    if (rank == 0) printf ("WARNING: No valid structure training set data read from %s!\n",tset_file);}
  else {
    if (rank == 0 && debug_level == 1) DEBUG_PRINT ((">> Found %i structure records in %s\n", i, tset_file));
    tset->struc = (struc_tset_item *) scalloc (i, sizeof (struc_tset_item), "tset_struc");
  }

  // check rms force
  i = get_num_entries (tset_file, "FORCE", "ENDFORCE", 1, rank);
  (*nfit)[FORCE] = i;
  if (!i) {
    if (rank == 0) printf ("WARNING: No valid rms force training set data read from %s!\n", tset_file);}
  else {
    if (rank == 0 && debug_level == 1) DEBUG_PRINT ((">> Found %i force records in %s\n", i, tset_file));
    tset->force = (force_tset_item *) scalloc (i, sizeof (force_tset_item), "tset_force");
  }

  // check atomic force
  i = get_num_entries (tset_file, "ATOM_FORCE", "ENDATOM_FORCE", 5, rank);
  (*nfit)[ATOM_FORCE] = i;
  if (!i) {
    if (rank == 0) printf ("WARNING: No valid atomic force training set data read from %s!\n", tset_file);}
  else {
    if (rank == 0 && debug_level == 1) DEBUG_PRINT ((">> Found %i force records in %s\n", i, tset_file));
    tset->atom_force = (atom_force_tset_item *) scalloc (i, sizeof (atom_force_tset_item), "tset_atom_force");
  }


  //check for energies
  i = get_num_entries (tset_file, "ENERGY", "ENDENERGY", 4, rank);
  (*nfit)[ENERGY] = i;
  if (!i) {
    if (rank == 0) printf ("WARNING: No valid energy training set data read from %s!\n", tset_file);}
  else {
    if (rank == 0 && debug_level == 1) DEBUG_PRINT ((">> Found %i energy records in %s\n", i, tset_file));
    tset->eng = (eng_tset_item *) scalloc (i, sizeof (eng_tset_item), "tset_eng");
  }

  if ((*nfit)[ENERGY] == 0 && (*nfit)[STRUCTURE] == 0 && (*nfit)[CHARGE] == 0 && (*nfit)[CELL] == 0
      && (*nfit)[FORCE] == 0 && (*nfit)[STRESS] == 0 && (*nfit)[GEOMETRY] == 0) {
    if (rank == 0) printf ("ERROR: No valid data in training set read from %s !\n", tset_file);
    exit (1);
  }

  //read energies
  sprintf(file, "%s", tset_file);
  fp = fopen (file, "r");
  if (fp == NULL) {
    if (rank == 0) printf ("ERROR: Cannot open training set file %s!\n", tset_file);
    exit (1);
  }
  valid = 0;
  int first=1;
  calc_initial_geom=0;
  char key[20];
  // Initialize secweight
  for (i = 0; i < num_objectives; i++) secweight[i] = 0.0;
  // Read in secweight
  while (fgets (line, MAX_LINE, fp) != NULL) {
    if (strncmp (line, "CHARGE", 6) == 0) {
      sscanf (line,"%s %lf",key,&secweight[CHARGE]);
      if (secweight[CHARGE]==0.0) secweight[CHARGE]=1.0;
      valid = CHARGE;
      i = -1;
    }
    else if (strncmp (line, "ENDCHARGE", 9) == 0)
      valid = 0;
    else if (strncmp (line, "FORCE", 5) == 0) {
      sscanf (line,"%s %lf",key,&secweight[FORCE]);
      if (secweight[FORCE]==0.0) secweight[FORCE]=1.0;
      valid = FORCE;
      i = -1;
    }
    else if (strncmp (line, "ENDFORCE", 8) == 0)
      valid = 0;
      
         else if (strncmp (line, "ATOM_FORCE", 10) == 0) {
      sscanf (line,"%s %lf",key,&secweight[ATOM_FORCE]);
      if (secweight[ATOM_FORCE]==0.0) secweight[ATOM_FORCE]=1.0;
      valid = ATOM_FORCE;
      i = -1;
    }
    else if (strncmp (line, "ENDATOM_FORCE", 13) == 0)
      valid = 0;
      
      
    else if (strncmp (line, "STRESS", 6) == 0) {
      sscanf (line,"%s %lf",key,&secweight[STRESS]);
      if (secweight[STRESS]==0.0) secweight[STRESS]=1.0;
      valid = STRESS;
      i = -1;
    }
    else if (strncmp (line, "ENDSTRESS", 9) == 0)
      valid = 0;
    else if (strncmp (line, "PRESSURE", 8) == 0) {
      sscanf (line,"%s %lf",key,&secweight[PRESSURE]);
      if (secweight[PRESSURE]==0.0) secweight[PRESSURE]=1.0;
      valid = PRESSURE;
      i = -1;
    }
    else if (strncmp (line, "ENDPRESS", 8) == 0)
      valid = 0;
    else if (strncmp (line, "CELL PARAMETERS", 15) == 0) {
      sscanf (line,"%s %lf",key,&secweight[CELL]);
      if (secweight[CELL]==0.0) secweight[CELL]=1.0;
      valid = CELL;
      i = -1;
    }
    else if (strncmp (line, "ENDCELL PARAMETERS", 18) == 0)
      valid = 0;
    else if (strncmp (line, "FREQUENCIES", 11) == 0) {
      sscanf (line,"%s %lf",key,&secweight[FREQUENCY]);
      if (secweight[FREQUENCY]==0.0) secweight[FREQUENCY]=1.0;
      valid = FREQUENCY;
      i = -1;
    }
    else if (strncmp (line, "ENDFREQUENCIES", 14) == 0)
      valid = 0;
    else if (strncmp (line, "HEATFO", 6) == 0) {
      sscanf (line,"%s %lf",key,&secweight[HEATFORM]);
      if (secweight[HEATFORM]==0.0) secweight[HEATFORM]=1.0;
      valid = HEATFORM;
      i = -1;
    }
    else if (strncmp (line, "ENDHEATFO", 9) == 0)
      valid = 0;
    else if (strncmp (line, "GEOMETRY", 8) == 0) {
      sscanf (line,"%s %lf",key,&secweight[GEOMETRY]);
      if (secweight[GEOMETRY]==0.0) secweight[GEOMETRY]=1.0;
      valid = GEOMETRY;
      i = -1;
    }
    else if (strncmp (line, "ENDGEOMETRY", 11) == 0)
      valid = 0;
    else if (strncmp (line, "STRUCTURE", 9) == 0) {
      sscanf (line,"%s %lf",key,&secweight[STRUCTURE]);
      if (secweight[STRUCTURE]==0.0) secweight[STRUCTURE]=1.0;
      valid = STRUCTURE;
      i = -1;
    }
    else if (strncmp (line, "ENDSTRUCTURE", 12) == 0)
      valid = 0;
    else if (strncmp (line, "ENERGY", 6) == 0) {
      sscanf (line,"%s %lf",key,&secweight[ENERGY]);
      if (secweight[ENERGY]==0.0) secweight[ENERGY]=1.0;
      valid = ENERGY;
      i = -1;
    }
    else if (strncmp (line, "ENDENERGY", 9) == 0)
      valid = 0;
    else if (valid && line[0]!='#') {		// if first char in line is #, it's a comment, else parse it
      c = tokenize_string (line, &strtmp);
      if (valid == CHARGE && c < 4)
	continue;
      else if (valid == FORCE && c < 2)
	continue;
	 else if (valid == ATOM_FORCE && c < 6)
	continue;
      else if (valid == STRESS && c < 3)
	continue;
      else if (valid == CELL && c < 4)
	continue;
      else if (valid == FREQUENCY && c < 4)
	continue;
      else if (valid == HEATFORM && c < 3)
	continue;
      else if (valid == GEOMETRY && c < 4)
	continue;
      else if (valid == STRUCTURE && c < 5)
        continue;
      else if (valid == ENERGY && c < 5)
	continue;
      i++;
      if (valid == CHARGE) {
	strcpy (tset->charge[i].line, line);
	strcpy (tset->charge[i].sname, strtmp[0]);
	tset->charge[i].weight = atof (strtmp[1]);
	tset->charge[i].atomid = atoi (strtmp[2]);
	tset->charge[i].charge = atof (strtmp[3]);
        if (rank == 0 && debug_level == 1) DEBUG_PRINT (("CHARGE: %s %4.2f %d %4.2f\n",tset->charge[i].sname,
          tset->charge[i].weight, tset->charge[i].atomid, tset->charge[i].charge));
      }
      else if (valid == CELL) {
	strcpy (tset->cell[i].line, line);
	strcpy (tset->cell[i].sname, strtmp[0]);
	tset->cell[i].weight = atof (strtmp[1]);
	strcpy (tset->cell[i].type, strtmp[2]);
	tset->cell[i].lit = atof (strtmp[3]);
        if (rank == 0 && debug_level == 1) DEBUG_PRINT (("CELL: %s %4.2f %s %4.2f\n",tset->cell[i].sname,
          tset->cell[i].weight, tset->cell[i].type, tset->cell[i].lit));
      }
      else if (valid == FORCE) {
	strcpy (tset->force[i].line, line);
	strcpy (tset->force[i].sname, strtmp[0]);
	tset->force[i].weight = atof (strtmp[1]);
        if (rank == 0 && debug_level == 1) DEBUG_PRINT (("FORCE: %s %4.2f\n", tset->force[i].sname, tset->force[i].weight));
      }
      
            else if (valid == ATOM_FORCE) {
	strcpy (tset->atom_force[i].line, line);
	strcpy (tset->atom_force[i].sname, strtmp[0]);
	tset->atom_force[i].weight = atof (strtmp[1]);
	tset->atom_force[i].atomid = atof (strtmp[2]);
        tset->atom_force[i].fx = atof (strtmp[3]);	
   	tset->atom_force[i].fy = atof (strtmp[4]);
	tset->atom_force[i].fz = atof (strtmp[5]);
        if (rank == 0 && debug_level == 1) DEBUG_PRINT (("ATOM_FORCE: %s %4.2f %d %4.2f %4.2f %4.2f\n",
         tset->atom_force[i].sname, tset->atom_force[i].weight, tset->atom_force[i].atomid, 
         tset->atom_force[i].fx, tset->atom_force[i].fy, tset->atom_force[i].fz ));
      }
      
      else if (valid == STRESS) {
	strcpy (tset->stress[i].line, line);
	strcpy (tset->stress[i].sname, strtmp[0]);
	tset->stress[i].weight = atof (strtmp[1]);
	tset->stress[i].id = get_stress_component_id (strtmp, c);
//        tset->stress[i].lit = atof (strtmp[3]);
        if (rank == 0 && debug_level == 1) DEBUG_PRINT (("STRESS: %s %4.2f %d\n", tset->stress[i].sname,
          tset->stress[i].weight, tset->stress[i].id));
      }
      else if (valid == PRESSURE) {
        strcpy (tset->press[i].line, line);
        strcpy (tset->press[i].sname, strtmp[0]);
        tset->press[i].weight = atof (strtmp[1]);
        tset->press[i].id = get_press_component_id (strtmp[2], c);
        tset->press[i].lit = atof (strtmp[3]);
        if (rank == 0 && debug_level == 1) DEBUG_PRINT (("PRESS: %s %4.2f %d %4.2f\n", tset->press[i].sname,
          tset->press[i].weight, tset->press[i].id, tset->press[i].lit));
      }
      else if (valid == GEOMETRY) {
	strcpy (tset->geom[i].line, line);
	strcpy (tset->geom[i].sname, strtmp[0]);
	tset->geom[i].weight = atof (strtmp[1]);

        // if last value in entry is real, then QM values are included, else need to pre-calculate
        if (first) {
          calc_initial_geom=check_int_or_float (strtmp[c-1], calc_initial_geom);
          if (calc_initial_geom) 
            printf(">> Literal geometric values not provided, will compute assuming geometries are exact\n");
          first=0;
        }

	if (calc_initial_geom) {
          tset->geom[i].num_atom = n = c - 2;
          tset->geom[i].lit = 0.0;      // when no QM values provided, initialize to zero
        } else {
          tset->geom[i].num_atom = n = c - 3;
          tset->geom[i].lit = atof (strtmp[c - 1]);
        }
        if (rank == 0 && debug_level == 1) 
          DEBUG_PRINT (("GEOMETRY: %s %4.2f [%d]",tset->geom[i].sname, tset->geom[i].weight,
          tset->geom[i].num_atom));
	for (j = 0; j < n; j++) {
	  tset->geom[i].atom[j] = atoi (strtmp[j + 2]);
          if (rank == 0 && debug_level == 1) 
            DEBUG_PRINT ((" %d ", tset->geom[i].atom[j]));
	}
        if (rank == 0 && debug_level == 1) 
          DEBUG_PRINT (("%f\n", tset->geom[i].lit));
      }
      else if (valid == STRUCTURE) {
        strcpy (tset->struc[i].line, line);
        for (j = 0; j < (int) strlen (strtmp[0]); ++j)
          tset->struc[i].sname[j] = strtmp[0][j];
//        strcpy (tset->struc[i].sname, strtmp[0]);
        tset->struc[i].weight = atof (strtmp[1]);
        tset->struc[i].dis = atof (strtmp[2]);
        tset->struc[i].lit = atof (strtmp[3]);
        if (rank == 0 && debug_level == 1) DEBUG_PRINT (("STRUCTURE: %s %4.2f %4.2f = %4.2f\n", tset->struc[i].sname, tset->struc[i].weight,
          tset->struc[i].dis, tset->struc[i].lit));
      }
      else if (valid == ENERGY) {
	n = (int) ((c - 2) / 3);
	strcpy (tset->eng[i].line, line);
	tset->eng[i].weight = atof (strtmp[0]);
	tset->eng[i].eng = atof (strtmp[c - 1]);
	tset->eng[i].n = n;
	tset->eng[i].sname = (char **) scalloc (n, sizeof (char *), "tset_engi_sname");
	tset->eng[i].op = (char **) scalloc (n, sizeof (char *), "tset_engi_op");
	tset->eng[i].factor = (double *) scalloc (n, sizeof (double), "tset_engi_factor");
        if (rank == 0 && debug_level == 1) DEBUG_PRINT (("ENERGY: %4.2f", tset->eng[i].weight));
	for (j = 0; j < n; j++) {
	  tset->eng[i].op[j] = (char *) scalloc (2, sizeof (char), "tset_engi_opj");
	  tset->eng[i].sname[j] = (char *) scalloc (MAX_LINE, sizeof (char), "tset_engi_snamej");
	  strcpy (tset->eng[i].op[j], strtmp[j * 3 + 1]);
	  strcpy (tset->eng[i].sname[j], strtmp[j * 3 + 2]);
	  tset->eng[i].factor[j] = atof (strtmp[j * 3 + 3]);
          if (tset->eng[i].factor[j] == 0.0) {
            printf ("Please check your training set file Energy section for a zero factor\n");
            MPI_Abort (MPI_COMM_WORLD, 1);
          }
          if (rank == 0 && debug_level == 1) DEBUG_PRINT ((" %s %s / %d", tset->eng[i].op[j], tset->eng[i].sname[j],
            (int) tset->eng[i].factor[j]));
	}
        if (rank == 0 && debug_level == 1) DEBUG_PRINT ((" = %4.2f\n", tset->eng[i].eng));
      }
    }
  }

  // If -W option set global section weights randomly
  if (random_wflag) {
    srand(time(NULL));
    for (i=0; i < num_objectives; i++) {
      int r = rand() % (100 - 1) + 1;
        secweight[i] = (double) r;
    }
  }

  // Correct section weights if no cases are included in training set section
  for (i = 0; i < num_objectives; i++) 
    if ((*nfit)[i]==0) secweight[i]=0.0;

  // Normalize (num_objectives) objective weighs
  // i.e. CHARGE, CELL,FREQUENCY,HEATFORM,GEOMETRY,STRUCTURE,ENERGY,FORCE,STRESS
  for (i = 0; i < num_objectives; i++)
    tset->sumsecweights += secweight[i];

  tset->sumentryweights = (double *) scalloc (num_objectives, sizeof (double), "tset_sumentryweights");   
  for (i = 0; i < (*nfit)[CHARGE]; i++) 
    tset->sumentryweights[CHARGE]+=tset->charge[i].weight;
  for (i = 0; i < (*nfit)[CELL]; i++) 
    tset->sumentryweights[CELL]+=tset->cell[i].weight;
/*  for (i = 0; i < (*nfit)[FREQUENCY]; i++) 
    tset->sumentryweights[FREQUENCY]+=tset->freq[i].weight;
  for (i = 0; i < (*nfit)[HEATFORM]; i++) 
    tset->sumentryweights[HEATFORM]+=tset->heat[i].weight;
*/
  for (i = 0; i < (*nfit)[GEOMETRY]; i++) 
      tset->sumentryweights[GEOMETRY]+=tset->geom[i].weight;
  for (i = 0; i < (*nfit)[STRUCTURE]; i++) 
    tset->sumentryweights[STRUCTURE]+=tset->struc[i].weight;
  for (i = 0; i < (*nfit)[ENERGY]; i++) 
    tset->sumentryweights[ENERGY]+=tset->eng[i].weight;
  for (i = 0; i < (*nfit)[FORCE]; i++) 
    tset->sumentryweights[FORCE]+=tset->force[i].weight;
  for (i = 0; i < (*nfit)[ATOM_FORCE]; i++)
    tset->sumentryweights[ATOM_FORCE]+=tset->atom_force[i].weight;
  for (i = 0; i < (*nfit)[STRESS]; i++) 
    tset->sumentryweights[STRESS]+=tset->stress[i].weight;
  for (i = 0; i < (*nfit)[PRESSURE]; i++)
    tset->sumentryweights[PRESSURE]+=tset->press[i].weight;

  if (normalizederrors) {
    for (i = 0; i < (*nfit)[CHARGE]; i++)
      tset->charge[i].nweight = (secweight[CHARGE]/tset->sumsecweights)*(tset->charge[i].weight/tset->sumentryweights[CHARGE]);
    for (i = 0; i < (*nfit)[CELL]; i++)
      tset->cell[i].nweight = (secweight[CELL]/tset->sumsecweights)*(tset->cell[i].weight/tset->sumentryweights[CELL]);
/*    for (i = 0; i < (*nfit)[FREQUENCY]; i++)
      tset->freq[i].nweight = (secweight[FREQUENCY]/tset->sumsecweights)*(tset->freq[i].weight/tset->sumentryweights[FREQUENCY]);
    for (i = 0; i < (*nfit)[HEATFORM]; i++)
      tset->heat[i].nweight = (secweight[HEATFORM]/tset->sumsecweights)*(tset->heat[i].weight/tset->sumentryweights[HEATFORM]);
*/    for (i = 0; i < (*nfit)[GEOMETRY]; i++)
      tset->geom[i].nweight = (secweight[GEOMETRY]/tset->sumsecweights)*(tset->geom[i].weight/tset->sumentryweights[GEOMETRY]);
    for (i = 0; i < (*nfit)[STRUCTURE]; i++)
      tset->struc[i].nweight = (secweight[STRUCTURE]/tset->sumsecweights)*(tset->struc[i].weight/tset->sumentryweights[STRUCTURE]);
    for (i = 0; i < (*nfit)[ENERGY]; i++)
      tset->eng[i].nweight = (secweight[ENERGY]/tset->sumsecweights)*(tset->eng[i].weight/tset->sumentryweights[ENERGY]);
    for (i = 0; i < (*nfit)[FORCE]; i++)
      tset->force[i].nweight = (secweight[FORCE]/tset->sumsecweights)*(tset->force[i].weight/tset->sumentryweights[FORCE]);
    for (i = 0; i < (*nfit)[ATOM_FORCE]; i++)
      tset->atom_force[i].nweight = (secweight[ATOM_FORCE]/tset->sumsecweights)*(tset->atom_force[i].weight/tset->sumentryweights[ATOM_FORCE]);
    for (i = 0; i < (*nfit)[STRESS]; i++)
      tset->stress[i].nweight = (secweight[STRESS]/tset->sumsecweights)*(tset->stress[i].weight/tset->sumentryweights[STRESS]);
    for (i = 0; i < (*nfit)[PRESSURE]; i++)
      tset->press[i].nweight = (secweight[PRESSURE]/tset->sumsecweights)*(tset->press[i].weight/tset->sumentryweights[PRESSURE]);
  }

  // Cleanup
  fclose (fp);
  for (i = 0; i < MAX_TOKENS; i++)
    free (strtmp[i]);
  free (strtmp);
  free (line);
}


/****************************************************************************
*  Function     : get_numb_entries
*  Description  : Counts the number of entries per line in the training set file
                  for a particular fitness function objective (i.e. charges, 
                  geometries, energies, etc.)
*  Parameters   : Training set file name (file_name), fitness function objective
                  start label, fitness function objective end label, minimum number 
                  of parameters per entry type, and MPI process rank 
*  Effects      : Returns integer value corresponding to the number of entries
                  for the particular fitness function objective 
****************************************************************************/

int get_num_entries (char *file_name, char *start_label, char *end_label, int min, int rank)
{
  int i, c, sl, el, valid;
  char **strtmp;
  char filen[MAX_LINE];
  FILE *fpg;
  char *line;

  sl = strlen (start_label);
  el = strlen (end_label);

  line = (char*) scalloc(MAX_LINE, sizeof(char), "line2");
  strtmp = (char **) scalloc (MAX_TOKENS, sizeof (char *), "tmp2");
  for (i = 0; i < MAX_TOKENS; i++)
    strtmp[i] = (char *) scalloc (MAX_TOKEN_LEN, sizeof (char), "tmp2i");

  sprintf(filen, "%s", file_name);
  fpg = fopen (filen, "r");
  if (fpg == NULL) {
    if (rank == 0) printf ("ERROR: Cannot open training set file %s!\n", file_name);
    exit (1);
  }
  valid = i = 0;
  while (fgets (line, MAX_LINE, fpg) != NULL) {
    if (strncmp (line, start_label, sl) == 0)
      valid = 1;
    else if (strncmp (line, end_label, el) == 0)
      valid = 0;
    else if (valid && line[0]!='#') {
      c = tokenize_string (line, &strtmp);
      if (c > min)
	i++;
    }
  }
  fclose (fpg);
  for (c = 0; c < MAX_TOKENS; c++)
    free (strtmp[c]);
  free (strtmp);
  free (line);

  return i;
}

/****************************************************************************
*  Function     : get_stress_component_id
*  Description  : Returns a particular stress component int id from the stress tensor,
                  as specified by the user in the input argument vals (principals
                  and products of stress, i.e. xx, yy, zz, xy, yz, xz) 
*  Parameters   : Returns the integer stress component id
*  Effects      : None
****************************************************************************/

int get_stress_component_id (char **vals, int n)
{
  int stress_id = 0;
  int i;

  int XX = 1 << 6;
  int YY = 1 << 5;
  int ZZ = 1 << 4;
  int XY = 1 << 3;
  int XZ = 1 << 2;
  int YZ = 1 << 1;

  for (i = 2; i < n; i++) {
    if (strcmp (vals[i], "xx") == 0) {
      stress_id |= XX;
    }
    else if (strcmp (vals[i], "yy") == 0) {
      stress_id |= YY;
    }
    else if (strcmp (vals[i], "zz") == 0) {
      stress_id |= ZZ;
    }
    else if (strcmp (vals[i], "xy") == 0 || strcmp (vals[i], "yx") == 0) {
      stress_id |= XY;
    }
    else if (strcmp (vals[i], "xz") == 0 || strcmp (vals[i], "zx") == 0) {
      stress_id |= XZ;
    }
    else if (strcmp (vals[i], "yz") == 0 || strcmp (vals[i], "zy") == 0) {
      stress_id |= YZ;
    }
  }

  return stress_id;
}

/****************************************************************************
   Function	: get_press_component_id
*  Description  : Returns a particular pressure component int id from the tensor,
                  as specified by the user in the input argument vals (principals
                  and products of press, i.e. xx, yy, zz, xy, yz, xz) 
*  Parameters   : Returns the integer pressure component id
*  Effects      : None
****************************************************************************/

int get_press_component_id (char *vals, int n)
{
  int press_id = 0;
  // xx, yy, zz, xy, xz, yz
  
  if (strcmp (vals, "xx") == 0)
    press_id = 1;
  else if (strcmp (vals, "yy") == 0) 
    press_id = 2;
  else if (strcmp (vals, "zz") == 0) 
    press_id = 3;
  else if (strcmp (vals, "xy") == 0 || strcmp (vals, "yx") == 0) 
    press_id = 4;
  else if (strcmp (vals, "xz") == 0 || strcmp (vals, "zx") == 0) 
    press_id = 5;
  else if (strcmp (vals, "yz") == 0 || strcmp (vals, "zy") == 0) 
    press_id = 6;
// press_id = 0 -> scalar 'press'  

  return press_id;
}

/****************************************************************************
 *  Function     : check_int_or_float
 *  Description  : Parses a string to determine if it's an integer or a real
 *  Parameters   : Returns a 1 if string is an integer, 0 otherwise
 *  Effects      : None
 *****************************************************************************/

int check_int_or_float (char *string, int result) 
{
  int idx, slen;
  int numbersfound = 0;
  int dotsfound = 0;

  slen=strlen(string);
  for (idx=0; (idx<slen) && (string[idx]); idx++){
    numbersfound += ((string[idx] >= '0') && (string[idx] <= '9'));
    dotsfound += (string[idx] == '.');
  }
  
  if (numbersfound == slen) 
    result=1;	// integer
  else if ((dotsfound == 1) && (numbersfound == slen-1))
    result=0;	// float

  return result;
}
