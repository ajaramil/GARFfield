/* ----------------------------------------------------------------------
   GARFFIELD: Genetic-Algorithm Reactive Force FIELD parameter optimizer.

   Copyright (2012-2014) California Institute of Technology
   Andres Jaramillo-Botero (ajaramil@caltech.edu)
   http://www.wag.caltech.edu/home/ajaramil/GARFfield.html
------------------------------------------------------------------------- */

#ifndef EFF_FFIELD_H
#define EFF_FFIELD_H

#include "eff_types.h"

extern int debug_level;

void Write_Force_Field_Eff (char *, eff_interaction *, double **, int);
int Read_Force_Field_Eff (char *, eff_interaction *, double ***, int);

#endif
