/* ----------------------------------------------------------------------
   GARFFIELD: Genetic-Algorithm Reactive Force FIELD parameter optimizer.

   Copyright (2012-2014) California Institute of Technology
   Andres Jaramillo-Botero (ajaramil@caltech.edu)
   http://www.wag.caltech.edu/home/ajaramil/GARFfield.html
------------------------------------------------------------------------- */

#ifndef PQEQ_FFIELD_H
#define PQEQ_FFIELD_H

#include "pqeq_types.h"

extern int debug_level;

void Write_Force_Field_PQeq (char *, pqeq_interaction *, double **, int);
int Read_Force_Field_PQeq (char *, pqeq_interaction *, double ***, int);

#endif
