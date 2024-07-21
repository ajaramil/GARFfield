/* ----------------------------------------------------------------------
   GARFFIELD: Genetic-Algorithm Reactive Force FIELD parameter optimizer.

   Copyright (2012-2014) California Institute of Technology
   Andres Jaramillo-Botero (ajaramil@caltech.edu)
   http://www.wag.caltech.edu/home/ajaramil/GARFfield.html
------------------------------------------------------------------------- */

#ifndef REAXC_FFIELD_H
#define REAXC_FFIELD_H

#include "reaxc_types.h"
#include "tool_box.h"

void Write_Force_Field_Reax (char *, reax_interaction *, double **, int);
int Read_Force_Field_Reax (char *, reax_interaction *, double ***, int);

extern int tokenize_string (char *, char ***);
extern int debug_level;
extern int lg_reax;
extern int initial_write;

int gp_idx, sbp_idx, tbp_idx, tbodp_idx, thbp_idx, fbp_idx, hbp_idx;

#endif
