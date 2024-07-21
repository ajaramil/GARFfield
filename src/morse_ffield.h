/* ----------------------------------------------------------------------
   GARFFIELD: Genetic-Algorithm Reactive Force FIELD parameter optimizer.

   Copyright (2012-2014) California Institute of Technology
   Andres Jaramillo-Botero (ajaramil@caltech.edu)
   http://www.wag.caltech.edu/home/ajaramil/GARFfield.html
------------------------------------------------------------------------- */

#ifndef MORSE_FFIELD_H
#define MORSE_FFIELD_H

#include "morse_types.h"

extern int debug_level;

void Write_Force_Field_Morse (char *, morse_interaction *, double **, int);
int Read_Force_Field_Morse (char *, morse_interaction *, double ***, int);

#endif
