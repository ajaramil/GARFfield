/* ----------------------------------------------------------------------
    GARFFIELD: Genetic-Algorithm Reactive Force FIELD parameter optimizer.

    Copyright (2012) California Institute of Technology
    Andres Jaramillo-Botero (ajaramil@caltech.edu) 
    http://www.wag.caltech.edu/home/ajaramil/GARFfield.html
*------------------------------------------------------------------------- */

#ifndef PARAMS_H
#define PARAMS_H

void parse_params (char *, int, int *, int);
int get_num_params (char *, int);
void get_params (char *, int, int);

#include "eff_ffield.h"
#include "reaxc_ffield.h"
#include "pqeq_ffield.h"
#include "ppr_ffield.h"
#include "morse_ffield.h"
#include "comb_ffield.h"
#include "tersoff_ffield.h"
#include "tersoff_mod_ffield.h"
#include "zhou_EAM_ffield.h"
#include "cg_ffield.h"

extern double *Lower, *Upper;
extern reax_interaction *ffdata;
extern eff_interaction *effdata;
extern pqeq_interaction *pqeqffdata;
extern ppr_interaction *pprffdata;
extern cg_interaction *cgffdata;
extern morse_interaction *morseffdata;
extern comb_interaction *combffdata;
extern tersoff_interaction *tersoffffdata;
extern tersoff_mod_interaction *tersoff_modffdata;
extern zhou_EAM_interaction *zhou_EAMffdata;
extern int lg_reax;

extern int gp_idx, sbp_idx, tbp_idx, tbodp_idx, thbp_idx, fbp_idx, hbp_idx;

#endif
