/* ----------------------------------------------------------------------
   GARFFIELD: Genetic-Algorithm Reactive Force FIELD parameter optimizer.

   Copyright (2012-2014) California Institute of Technology
   Andres Jaramillo-Botero (ajaramil@caltech.edu)
   http://www.wag.caltech.edu/home/ajaramil/GARFfield.html
------------------------------------------------------------------------- */

#ifndef TERSOFF_MOD_FFIELD_H
#define TERSOFF_MOD_FFIELD_H

#include "tersoff_mod_types.h"

void Write_Force_Field_TERSOFF_MOD (char *, tersoff_mod_interaction *, double **, int);
int Read_Force_Field_TERSOFF_MOD (char *, tersoff_mod_interaction *, double ***, int);

extern int tokenize_string (char *, char ***);
extern int debug_level;
extern int initial_write;
extern void *scalloc (int, int, char *);

#endif
