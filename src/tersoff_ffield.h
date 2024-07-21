/* ----------------------------------------------------------------------
   GARFFIELD: Genetic-Algorithm Reactive Force FIELD parameter optimizer.

   Copyright (2012-2014) California Institute of Technology
   Andres Jaramillo-Botero (ajaramil@caltech.edu)
   http://www.wag.caltech.edu/home/ajaramil/GARFfield.html
------------------------------------------------------------------------- */
/* ----------------------------------------------------------------------
    Code derived and contributed by Anurag Chaudhry, 2015
------------------------------------------------------------------------ */

#ifndef TERSOFF_FFIELD_H
#define TERSOFF_FFIELD_H

#include "tersoff_types.h"

void Write_Force_Field_TERSOFF (char *, tersoff_interaction *, double **, int);
int Read_Force_Field_TERSOFF (char *, tersoff_interaction *, double ***, int);

extern int tokenize_string (char *, char ***);
extern int debug_level;
extern int initial_write;
extern void *scalloc (int, int, char *);

#endif
