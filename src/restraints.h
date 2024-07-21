/* ----------------------------------------------------------------------
   GARFFIELD: Genetic-Algorithm Reactive Force FIELD parameter optimizer.

   Copyright (2012-2014) California Institute of Technology
   Andres Jaramillo-Botero (ajaramil@caltech.edu)
   http://www.wag.caltech.edu/home/ajaramil/GARFfield.html
------------------------------------------------------------------------- */
#ifndef RESTRAIN_H
#define RESTRAIN_H

typedef struct
{
  int nbody;			// # bodies in the restraint
  int atom[4];			// atom ids
  char name[100];		// nametag of geo file as set by DESCRP line in bgf data or XYZ line in xyz data
  double f1, f2;		// force1 and force2
  double val;			// restrained value
} restraint;

restraint *rstrain;
extern char *restraints_file;
int nrst;
int parse_restraints (char *, int);
#endif
