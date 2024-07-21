/* ----------------------------------------------------------------------
   GARFFIELD: Genetic-Algorithm Reactive Force FIELD parameter optimizer.

   Copyright (2012-2014) California Institute of Technology
   Andres Jaramillo-Botero (ajaramil@caltech.edu)
   http://www.wag.caltech.edu/home/ajaramil/GARFfield.html
------------------------------------------------------------------------- */
/* ----------------------------------------------------------------------
    Code derived and contributed by Anurag Chaudhry, 2015
 ------------------------------------------------------------------------ */

#ifndef ZHOU_EAM_FFIELD_H
#define ZHOU_EAM_FFIELD_H

#include "zhou_EAM_types.h"

double V(double,double,double,double,double,double,double,double );
double rho(double,double,double,double,double );
double F(double,double,double,double,double,double,double,double,double,double,double,double,double );
void Write_Force_Field_ZHOU_EAM (char *, zhou_EAM_interaction *, double **, int);
int Read_Force_Field_ZHOU_EAM (char *, zhou_EAM_interaction *, double ***, int);

extern int tokenize_string (char *, char ***);
extern int debug_level;
extern int initial_write, final_write;
extern void *scalloc (int, int, char *);

#endif
