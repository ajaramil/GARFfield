/* ----------------------------------------------------------------------
   GARFFIELD: Genetic-Algorithm Reactive Force FIELD parameter optimizer.

   Copyright (2012-2014) California Institute of Technology
   Andres Jaramillo-Botero (ajaramil@caltech.edu)
   http://www.wag.caltech.edu/home/ajaramil/GARFfield.html
------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "library.h"		/* this is the LAMMPS library include file */
#include "mpi.h"

#include "ff_lammps.h"
#include "errors.h"
#include "restraints.h"
#include "reaxc_types.h"
#include "eff_types.h"
#include "cg_types.h"
#include "trainset.h"
#include "structures.h"

#ifdef DEBUG
# define DEBUG_PRINT(x) printf x
#else
# define DEBUG_PRINT(x)
#endif

/****************************************************************************
*  Function     : calcerror
*  Description  : Computes the total error for all training cases accumulated
                  during the fitness function evaluation (get_lammps_eng). The
                  error function may be selected from the command line (-e) to
                  be PME (Percent Mean Error), RMSE (Root Mean Square Error),
                  NRMSE (max/min normalized RMSE), and SIWE (Square inverse
                  Weighted Error).
     
0. SIWE = sum((QM-FF)^2/wi)		// square inverse weight error
1. MPE = sum(wi'*(100*|QM-FF|/QM)/n)	// normalized mean percent error
2. RMSE = sum(sqrt(wi'*(QM-FF)^2)/n)	// root mean square error
3. NRMSE = RMSE/(FFmax-FFmin)		// normalized RMSE
4. MSE = sum(wi(QM-FF)^2)		// mean square error
where wi'=[sec_wi/sum(sec_wi)]*wi/sum(wi)

*  Parameters   : Energies array and MPI process rank
*  Effects      : Returns the total weighted error, based on calculated energies
                  for the entire training set, per GA chromosome entry
****************************************************************************/

double calcerror(double *energies, int rank)
{
    int i, j, id;
    double eval, partial_error, Emax, Emin, Emax_minus_Emin;
    FILE *fp;
    double press;
    char filename[MAX_LINE];
    double total_error, diff, force_ff, force_qm, acc_percent_diff, percent_diff, acc_error;

    /* open error report file */
    if (initial_write)
	sprintf(filename, "trainset.err.initial");
    else if (final_write)
	sprintf(filename, "trainset.err.best");
    else sprintf(filename, "trainset.err.%d", rank);
    fp = fopen(filename, "w");

    fprintf(fp, "\n== Error evaluation == \n");
    // clear error accumulator and initialize vars
    partial_error = acc_percent_diff = acc_error = press = 0.0;
    for (i = 0; i < 11; i++)
	error[i] = 0.0;

/*
0. SIWE = sum((QM-FF)^2/wi)             // square inverse weight error
1. MPE = sum(wi'*(100*|QM-FF|/QM)/n)   // normalized mean percent error
2. RMSE = sum(sqrt(wi'*(QM-FF)^2)/n)    // root mean square error
3. NRMSE = RMSE/(FFmax-FFmin)           // normalized RMSE
4. MSE = sum(wi(QM-FF)^2)               // mean square error
where wi'=[sec_wi/sum(sec_wi)]*wi/sum(wi)
*/
    char errorname[MAX_LINE];
    if (errorfunction == 0) {
	strcpy(errorname, "SIWE");
	fprintf(fp, "W. Error (SIWE) = sum((QM-FF)^2/wi)\n");
    } else if (errorfunction == 1) {
	strcpy(errorname, "MPE");
	fprintf(fp,
		"W. Error (MPE) = sum(wi'*(100*|QM-FF|/QM)/n) where wi'=[sec_wi/sum(sec_wi)]*wi/sum(wi)\n");
    } else if (errorfunction == 2) {
	strcpy(errorname, "RMSE");
	fprintf(fp,
		"W. Error (RMSE) = sum(sqrt(wi'*(QM-FF)^2)/n) where wi'=[sec_wi/sum(sec_wi)]*wi/sum(wi)\n");
    } else if (errorfunction == 3) {
	strcpy(errorname, "NRMSE");
	fprintf(fp,
		"W. Error (NRMSE) = sum(sqrt(wi'*(QM-FF)^2)/n)/(FFmax-FFmin) where wi'=[sec_wi/sum(sec_wi)]*wi/sum(wi)\n");
    } else if (errorfunction == 4) {
	strcpy(errorname, "MSE");
	fprintf(fp, "W. Error (MSE) = sum(wi*(QM-FF)^2)\n");
    }
    if (!minimization) fprintf (fp, "Using RMS force as fitness function\n");

    // CHARGE ERROR CALCULATION
    for (i = 0; i < nfit[CHARGE]; i++) {
	if (i == 0) {
	    if (normalizederrors)
		fprintf(fp,
			"CHARGES %lf\nWeight\t[nWeight]\tStructure\tQM\t\tForce Field\tW. Error\tAcc. W. Error\t%% Difference\n",
			secweight[CHARGE]);
	    else
		fprintf(fp,
			"CHARGES\nWeight\tStructure\tQM\tForce Field\tW. Error\tAcc. W. Error\t%% Difference\n");
	}
	diff = tset->charge[i].charge - tset->charge[i].ff_val;
	percent_diff = 100 * fabs(diff / tset->charge[i].charge);
	acc_percent_diff += percent_diff;

	if (errorfunction == 0)	// SIWE (old)
	    partial_error = pow(diff / tset->charge[i].weight, 2);
	else if (errorfunction == 1)	// MPE (default)
	    partial_error = tset->charge[i].nweight * percent_diff;
	else if (errorfunction == 2 || errorfunction == 3)	// RMSE or NRMSE
	    partial_error = tset->charge[i].nweight * pow(diff, 2);
	else if (errorfunction == 4)	// MSE
	    partial_error = pow(tset->charge[i].weight * diff, 2);

	error[CHARGE] += partial_error;
	acc_error += partial_error;

	if (errorfunction == 3) {
	    if (i == 0) {
		Emax = tset->charge[i].ff_val;
		Emin = tset->charge[i].ff_val;
	    } else if (tset->charge[i].ff_val > Emax)
		Emax = tset->charge[i].ff_val;
	    else if (tset->charge[i].ff_val < Emin)
		Emin = tset->charge[i].ff_val;
	}
	// Write to trainset.err
	if (normalizederrors)
	    fprintf(fp, "%-4.3f[%4.3f]\t%-10s\t%-10.4f%-10.4f%-10.4f[%-10.4f]\t[%-10.4f]\n",
		    tset->charge[i].weight, tset->charge[i].nweight, tset->charge[i].sname,
		    tset->charge[i].charge, tset->charge[i].ff_val, partial_error, acc_error,
		    percent_diff);
	else
	    fprintf(fp, "%-4.3f\t%-10s\t%-10.4f%-10.4f%-10.4f[%-10.4f]\t[%-10.4f]\n",
		    tset->charge[i].weight, tset->charge[i].sname,
		    tset->charge[i].charge, tset->charge[i].ff_val, partial_error, acc_error,
		    percent_diff);
    }
    // Print objective function section summary
    if (nfit[CHARGE] > 0) {
	if (errorfunction == 0)	// SIWE
	    fprintf(fp,
		    "[%s] Total Inverse Weighted Error [CELL] = %-10.4f (|Acc. %% diff| = %10.4f)\n",
		    errorname, error[CHARGE], acc_percent_diff);
	else if (errorfunction == 1 || errorfunction == 4) {	// MPE or MSE
	    fprintf(fp, "[%s] Total Mean (Percent or Squared) Error [CHARGE] = %4.2f/%d = ",
		    errorname, error[CHARGE], nfit[CHARGE]);
	    error[CHARGE] /= nfit[CHARGE];
	    fprintf(fp, "%-10.4f (|Acc. %% diff| = %10.4f)\n", error[CHARGE], acc_percent_diff);
	} else if (errorfunction == 2) {	// RMSE
	    fprintf(fp, "[%s] Total Root Mean Square Error [CHARGE] = sqrt(%4.2f/%d) = ", errorname,
		    error[CHARGE], nfit[CHARGE]);
	    error[CHARGE] = sqrt(error[CHARGE] / nfit[CHARGE]);
	    fprintf(fp, "%-10.4f (|Acc. %% diff| = %10.4f)\n", error[CHARGE], acc_percent_diff);
	} else if (errorfunction == 3) {	// NRMSE = RMSE/(max-min)
	    if (Emax != Emin)
		Emax_minus_Emin = Emax - Emin;
	    else
		Emax_minus_Emin = 1.0;
	    fprintf(fp, "[%s] Total normalized RMSE [CHARGE] = sqrt(%4.2f/%d)/%4.2f = ", errorname,
		    error[CHARGE], nfit[CHARGE], Emax_minus_Emin);
	    error[CHARGE] = sqrt(error[CHARGE] / nfit[CHARGE]) / Emax_minus_Emin;
	    fprintf(fp, "%-10.4f (|Acc. %% diff| = %10.4f)\n", error[CHARGE], acc_percent_diff);
	}
    }
    
  
      // ATOM_FORCE ERROR CALCULATION
    for (i = 0; i < nfit[ATOM_FORCE]; i++) {
	if (i == 0) {
	    if (normalizederrors)
		fprintf(fp,
			"ATOM_FORCES %lf\nWeight\t[nWeight]\tStructure\tQM\t\tForce Field\tW. Error\tAcc. W. Error\t%% Difference\n",
			secweight[ATOM_FORCE]);
	    else
		fprintf(fp,
			"ATOM_FORCES\nWeight\tStructure\tQM\tForce Field\tW. Error\tAcc. W. Error\t%% Difference\n");
	}
	diff = sqrt((tset->atom_force[i].fx - tset->atom_force[i].ff_val_x)*(tset->atom_force[i].fx - tset->atom_force[i].ff_val_x) + 
	 (tset->atom_force[i].fy - tset->atom_force[i].ff_val_y)*(tset->atom_force[i].fy - tset->atom_force[i].ff_val_y) +  
	(tset->atom_force[i].fz - tset->atom_force[i].ff_val_z)*(tset->atom_force[i].fz - tset->atom_force[i].ff_val_z)) ;
    force_ff =  sqrt( (tset->atom_force[i].ff_val_x*tset->atom_force[i].ff_val_x + tset->atom_force[i].ff_val_y*tset->atom_force[i].ff_val_y +  
	tset->atom_force[i].ff_val_z*tset->atom_force[i].ff_val_z ) );       
    force_qm =  sqrt( (tset->atom_force[i].fx*tset->atom_force[i].fx + tset->atom_force[i].fy*tset->atom_force[i].fy+
	tset->atom_force[i].fz*tset->atom_force[i].fz ) );      
	percent_diff = 100 * fabs( diff / force_qm );
	acc_percent_diff += percent_diff;
   
	if (errorfunction == 0)	// SIWE (old)
	    partial_error = pow(diff / tset->atom_force[i].weight, 2);
	else if (errorfunction == 1)	// MPE (default)
	    partial_error = tset->atom_force[i].nweight * percent_diff;
	else if (errorfunction == 2 || errorfunction == 3)	// RMSE or NRMSE
	    partial_error = tset->atom_force[i].nweight * pow(diff, 2);
	else if (errorfunction == 4)	// MSE
	    partial_error = pow(tset->atom_force[i].weight * diff, 2);
	error[ATOM_FORCE] += partial_error;
	acc_error += partial_error;
	if (errorfunction == 3) {
	    if (i == 0) {
		Emax = force_ff;
		Emin = force_ff;
	    } else if (force_ff > Emax)
		Emax = force_ff;
	    else if (force_ff < Emin)
		Emin = force_ff;
        }
	// Write to trainset.err
	if (normalizederrors)
	    fprintf(fp, "%-4.3f[%4.3f]\t%-10s\t%-10.4f%-10.4f%-10.4f[%-10.4f]\t[%-10.4f]\n",
		    tset->atom_force[i].weight, tset->atom_force[i].nweight, tset->atom_force[i].sname,
		    force_qm, force_ff, partial_error, acc_error,
		    percent_diff);
	else
	    fprintf(fp, "%-4.3f\t%-10s\t%-10.4f%-10.4f%-10.4f[%-10.4f]\t[%-10.4f]\n",
		    tset->atom_force[i].weight, tset->atom_force[i].sname,
		    force_qm, force_ff, partial_error, acc_error,
		    percent_diff);
    }
 
    // Print objective function section summary
    if (nfit[ATOM_FORCE] > 0) {
	if (errorfunction == 0)	// SIWE
	    fprintf(fp,
		    "[%s] Total Inverse Weighted Error [CELL] = %-10.4f (|Acc. %% diff| = %10.4f)\n",
		    errorname, error[ATOM_FORCE], acc_percent_diff);
	else if (errorfunction == 1 || errorfunction == 4) {	// MPE or MSE
	    fprintf(fp, "[%s] Total Mean (Percent or Squared) Error [ATOM_FORCE] = %4.2f/%d = ",
		    errorname, error[ATOM_FORCE], nfit[ATOM_FORCE]);
	    error[ATOM_FORCE] /= nfit[ATOM_FORCE];
	    fprintf(fp, "%-10.4f (|Acc. %% diff| = %10.4f)\n", error[ATOM_FORCE], acc_percent_diff);
	} else if (errorfunction == 2) {	// RMSE
	    fprintf(fp, "[%s] Total Root Mean Square Error [ATOM_FORCE] = sqrt(%4.2f/%d) = ", errorname,
		    error[ATOM_FORCE], nfit[ATOM_FORCE]);
	    error[ATOM_FORCE] = sqrt(error[ATOM_FORCE] / nfit[ATOM_FORCE]);
	    fprintf(fp, "%-10.4f (|Acc. %% diff| = %10.4f)\n", error[ATOM_FORCE], acc_percent_diff);
	} else if (errorfunction == 3) {	// NRMSE = RMSE/(max-min)
	    if (Emax != Emin)
		Emax_minus_Emin = Emax - Emin;
	    else
		Emax_minus_Emin = 1.0;
	    fprintf(fp, "[%s] Total normalized RMSE [ATOM_FORCE] = sqrt(%4.2f/%d)/%4.2f = ", errorname,
		    error[ATOM_FORCE], nfit[ATOM_FORCE], Emax_minus_Emin);
	    error[ATOM_FORCE] = sqrt(error[ATOM_FORCE] / nfit[ATOM_FORCE]) / Emax_minus_Emin;
	    fprintf(fp, "%-10.4f (|Acc. %% diff| = %10.4f)\n", error[ATOM_FORCE], acc_percent_diff);
	}
    }
    
 
    
    
    // CELL PARAMETERS ERROR CALCULATION
    for (i = 0; i < nfit[CELL]; i++) {
	if (i == 0) {
	    if (normalizederrors)
		fprintf(fp,
			"CELL PARAMETERS %lf\nWeight\t[nWeight]\tStructure\tQM\t\tForce Field\tW. Error\tAcc. W. Error\t%% Difference\n",
			secweight[CELL]);
	    else
		fprintf(fp,
			"CELL PARAMETERS\nWeight\tStructure\tQM\tForce Field\tW. Error\tAcc. W. Error\t[%% difference]\n");
	}

	diff = tset->cell[i].lit - tset->cell[i].ff_val;
	percent_diff = 100 * fabs(diff / tset->cell[i].lit);
	acc_percent_diff += percent_diff;

	if (errorfunction == 0)	// SIWE
	    partial_error = pow(diff / tset->cell[i].weight, 2);
	else if (errorfunction == 1)	// MPE (default)
	    partial_error = tset->cell[i].nweight * percent_diff;
	else if (errorfunction == 2 || errorfunction == 3)	// RMSE or NRMSE
	    partial_error = tset->cell[i].nweight * pow(diff, 2);
	else if (errorfunction == 4)	// MSE
	    partial_error = pow(tset->cell[i].weight * diff, 2);

	error[CELL] += partial_error;
	acc_error += partial_error;

	if (errorfunction == 3) {
	    if (i == 0) {
		Emax = tset->cell[i].ff_val;
		Emin = tset->cell[i].ff_val;
	    } else if (tset->cell[i].ff_val > Emax)
		Emax = tset->cell[i].ff_val;
	    else if (tset->cell[i].ff_val < Emin)
		Emin = tset->cell[i].ff_val;
	}

	if (normalizederrors)
	    fprintf(fp, "%-4.3f\t[%4.3f]\t%-10s\t%-10.4f\t%-10.4f\t%-10.4f\t[%-10.4f]\t[%-10.4f]\n",
		    tset->cell[i].weight, tset->cell[i].nweight, tset->cell[i].sname,
		    tset->cell[i].lit, tset->cell[i].ff_val, partial_error, acc_error,
		    percent_diff);
	else
	    fprintf(fp, "%-4.3f\t%-10s\t%-10.4f\t%-10.4f\t%-10.4f\t[%-10.4f]\t[%-10.4f]\n",
		    tset->cell[i].weight, tset->cell[i].sname,
		    tset->cell[i].lit, tset->cell[i].ff_val, partial_error, acc_error,
		    percent_diff);
    }
    if (nfit[CELL] > 0) {
	if (errorfunction == 0)	// SIWE
	    fprintf(fp,
		    "[%s] Total Inverse Weighted Error [CELL] = %-10.4f (|Acc. %% diff| = %10.4f)\n",
		    errorname, error[CELL], acc_percent_diff);
	else if (errorfunction == 1 || errorfunction == 4) {	// MPE and MSE
	    fprintf(fp, "[%s] Total Mean (Pecent and Square) Error [CELL] (%4.2f/%d) = ", errorname,
		    error[CELL], nfit[CELL]);
	    error[CELL] = error[CELL] / nfit[CELL];
	    fprintf(fp, "%-10.4f (|Acc. %% diff| = %10.4f)\n", error[CELL], acc_percent_diff);
	} else if (errorfunction == 2) {	// RMSE
	    fprintf(fp, "[%s] RMSE [CELL] = sqrt(%4.2f/%d) = ", errorname, error[CELL], nfit[CELL]);
	    error[CELL] = sqrt(error[CELL] / nfit[CELL]);
	    fprintf(fp, "%-10.4f (|Acc. %% diff| = %10.4f)\n", error[CELL], acc_percent_diff);
	} else if (errorfunction == 3) {	// NRMSE = RMSE/(max-min)
	    if (Emax != Emin)
		Emax_minus_Emin = Emax - Emin;
	    else
		Emax_minus_Emin = 1.0;
	    fprintf(fp, "[%s] Normalized RMSE [CELL] = sqrt(%4.2f/%d)/%4.2f = ", errorname,
		    error[CELL], nfit[CELL], Emax_minus_Emin);
	    error[CELL] = sqrt(error[CELL] / nfit[CELL]) / Emax_minus_Emin;
	    fprintf(fp, "%-10.4f (|Acc. %% diff| = %10.4f)\n", error[CELL], acc_percent_diff);
	}
    }
    // GEOMETRY ERROR CALCULATION
    for (i = 0; i < nfit[GEOMETRY]; i++) {
	if (i == 0) {
	    if (normalizederrors)
		fprintf(fp,
			"GEOMETRY %lf\nWeight\t[nWeight]\tStructure\tQM\t\tForce Field\tW. Error\tAcc. W. Error\t%% difference\n",
			secweight[GEOMETRY]);
	    else
		fprintf(fp,
			"GEOMETRY\nWeight\tStructure\tQM\t\tForce Field\tW. Error\tAcc. W. Error\t%% difference\n");
	}

	diff = tset->geom[i].lit - tset->geom[i].ff_val;
        percent_diff = 100 * fabs(diff / tset->geom[i].lit);
	acc_percent_diff += percent_diff;

	if (errorfunction == 0)	// SIWE
	    partial_error = pow(diff / tset->geom[i].weight, 2);
	else if (errorfunction == 1)	// MPE (default)
	    partial_error = tset->geom[i].nweight * percent_diff;
	else if (errorfunction == 2 || errorfunction == 3)	// RMSE or NRMSE
	    partial_error = tset->geom[i].nweight * pow(diff, 2);
	else if (errorfunction == 4)	// MSE
	    partial_error = pow(tset->geom[i].weight * diff, 2);

	error[GEOMETRY] += partial_error;
	acc_error += partial_error;

	if (errorfunction == 3) {
	    if (i == 0) {
		Emax = tset->geom[i].ff_val;
		Emin = tset->geom[i].ff_val;
	    } else if (tset->geom[i].ff_val > Emax)
		Emax = tset->geom[i].ff_val;
	    else if (tset->geom[i].ff_val < Emin)
		Emin = tset->geom[i].ff_val;
	}

	if (normalizederrors)
	    fprintf(fp, "%-4.3f\t[%4.3f]\t%-10s\t%-10.4f\t%-10.4f\t%-10.4f\t%4.4f\t\t%4.4f\n",
		    tset->geom[i].weight, tset->geom[i].nweight, tset->geom[i].sname,
		    tset->geom[i].lit, tset->geom[i].ff_val, partial_error, acc_error,
		    percent_diff);
	else
	    fprintf(fp, "%-4.3f\t%-10s\t%-10.4f\t%-10.4f\t%-10.4f\t%4.4f\t\t%4.4f\n",
		    tset->geom[i].weight, tset->geom[i].sname,
		    tset->geom[i].lit, tset->geom[i].ff_val, partial_error, acc_error,
		    percent_diff);
    }
    if (nfit[GEOMETRY] > 0) {
	if (errorfunction == 0)	// SIWE
	    fprintf(fp,
		    "[%s] Total Squared Inverse Error [GEOMETRY] = %-10.4f (|Acc. %% diff| = %10.4f)\n",
		    errorname, error[GEOMETRY], acc_percent_diff);
	else if (errorfunction == 1 || errorfunction == 4) {	// MPE
	    fprintf(fp, "[%s] Total Mean (Percent or Squared) Error [GEOMETRY] = %4.2f/%d = ",
		    errorname, error[GEOMETRY], nfit[GEOMETRY]);
	    error[GEOMETRY] = error[GEOMETRY] / nfit[GEOMETRY];
	    fprintf(fp, "%-10.4f (|Acc. %% diff| = %10.4f)\n", error[GEOMETRY], acc_percent_diff);
	} else if (errorfunction == 2) {	// RMSE
	    fprintf(fp, "[%s] Total RMSE [GEOMETRY] = sqrt(%4.2f/%d) = ", errorname,
		    error[GEOMETRY], nfit[GEOMETRY]);
	    error[GEOMETRY] = sqrt(error[GEOMETRY] / nfit[GEOMETRY]);
	    fprintf(fp, "%-10.4f (|Acc. %% diff| = %10.4f)\n", error[GEOMETRY], acc_percent_diff);
	} else if (errorfunction == 3) {	// NRMSE = RMSE/(max-min)
	    if (Emax != Emin)
		Emax_minus_Emin = Emax - Emin;
	    else
		Emax_minus_Emin = 1.0;
	    fprintf(fp,
		    "[%s] Total Normalized RMSE [GEOMETRY] = RMSE/(Emax-Emin) = sqrt(%4.2f/%d)/%4.2f = ",
		    errorname, error[GEOMETRY], nfit[GEOMETRY], Emax_minus_Emin);
	    error[GEOMETRY] = sqrt(error[GEOMETRY] / nfit[GEOMETRY]) / Emax_minus_Emin;
	    fprintf(fp, "%-10.4f (|Acc. %% diff| = %10.4f)\n", error[GEOMETRY], acc_percent_diff);
	}
    }
    // STRUCTURE ERROR CALCULATION
    for (i = 0; i < nfit[STRUCTURE]; i++) {
	if (i == 0) {
	    press = getPress();
	    if (normalizederrors)
		fprintf(fp,
			"STRUCTURE %lf\nWeight\t[nWeight]\tStructure\tQM\tForce Field\tW. Error\tAcc. W. Error\t%% difference\n",
			secweight[STRUCTURE]);
	    else
		fprintf(fp,
			"STRUCTURE\nWeight\tStructure\tQM\tForce Field\tW. Error\tAcc. W. Error\t%% difference\n");
	}

	diff = tset->struc[i].lit - tset->struc[i].ff_val;
	percent_diff = 100 * fabs(diff / tset->struc[i].lit);
	acc_percent_diff += percent_diff;

	if (errorfunction == 0)	// SIWE
	    partial_error = pow(diff / tset->struc[i].weight, 2);
	else if (errorfunction == 1)	// MPE (default)
	    partial_error = tset->struc[i].nweight * percent_diff;
	else if (errorfunction == 2 || errorfunction == 3)	// RMSE or NRMSE
	    partial_error = tset->struc[i].nweight * pow(diff, 2);
	else if (errorfunction == 4)	// MSE
	    partial_error = pow(tset->struc[i].weight * diff, 2);

	error[STRUCTURE] += partial_error;
	acc_error += partial_error;

	if (errorfunction == 3) {
	    if (i == 0) {
		Emax = tset->struc[i].ff_val;
		Emin = tset->struc[i].ff_val;
	    }
	    if (tset->struc[i].ff_val > Emax)
		Emax = tset->struc[i].ff_val;
	    else if (tset->struc[i].ff_val < Emin)
		Emin = tset->struc[i].ff_val;
	}

	if (normalizederrors)
	    fprintf(fp, "%-4.3f\t[%4.3f]\t%-10s\t%-10.4f\t%-10.4f\t%-10.4f\t[%-10.4f]\t[%-10.4f]\n",
		    tset->struc[i].weight, tset->struc[i].nweight, tset->struc[i].sname,
		    tset->struc[i].lit, tset->struc[i].ff_val, partial_error, acc_error + press,
		    percent_diff);
	else
	    fprintf(fp, "%-4.3f\t%-10s\t%-10.4f\t%-10.4f\t%-10.4f\t[%-10.4f]\t[%-10.4f]\n",
		    tset->struc[i].weight,
		    tset->struc[i].sname, tset->struc[i].lit, tset->struc[i].ff_val, partial_error,
		    acc_error + press, percent_diff);
    }
    if (nfit[STRUCTURE] > 0) {
	if (errorfunction == 0)	// SIWE 
	    fprintf(fp,
		    "[%s] Total Squared Inverse Error [STRUCTURE] = %-10.4f (|Acc. %% diff| = %10.4f)\n",
		    errorname, error[STRUCTURE], acc_percent_diff);
	if (errorfunction == 1 || errorfunction == 4) {	// MPE
	    fprintf(fp, "[%s] Total Mean (Percent and Square) Error [STRUCTURE] = %4.2f/%d = ",
		    errorname, error[STRUCTURE], nfit[STRUCTURE]);
	    error[STRUCTURE] = error[STRUCTURE] / nfit[STRUCTURE];
	    fprintf(fp, "%-10.4f (|Acc. %% diff| = %10.4f)\n", error[STRUCTURE], acc_percent_diff);
	}
	if (errorfunction == 2) {	// RMSE
	    fprintf(fp, "[%s] Total RMSE [STRUCTURE] = sqrt(%4.2f/%d) = ", errorname,
		    error[STRUCTURE], nfit[STRUCTURE]);
	    error[STRUCTURE] = sqrt(error[STRUCTURE] / nfit[STRUCTURE]);
	    fprintf(fp, "%-10.4f (|Acc. %% diff| = %10.4f)\n", error[STRUCTURE], acc_percent_diff);
	}
	if (errorfunction == 3) {	// NRMSE = RMSE/(max-min)
	    if (Emax != Emin)
		Emax_minus_Emin = Emax - Emin;
	    else
		Emax_minus_Emin = 1.0;
	    fprintf(fp,
		    "[%s] Total Normalized RMSE [STRUCTURE] = RMSE/(Emax-Emin) = sqrt(%4.2f/%d)/%4.2f = ",
		    errorname, error[STRUCTURE], nfit[STRUCTURE], Emax_minus_Emin);
	    error[STRUCTURE] = sqrt(error[STRUCTURE] / nfit[STRUCTURE]) / Emax_minus_Emin;
	    fprintf(fp, "%-10.4f (|Acc. %% diff| = %10.4f)\n", error[STRUCTURE], acc_percent_diff);
	}
    }
    // FORCE ERROR CALCULATION
    for (i = 0; i < nfit[FORCE]; i++) {
	if (i == 0) {
	    if (normalizederrors)
		fprintf(fp,
			"FORCE %lf\nWeight\t[nWeight]\tStructure\tForce Field\tW. Error\tAcc. W. Error\t%% difference\n",
			secweight[FORCE]);
	    else
		fprintf(fp,
			"FORCE\nWeight\tStructure\tForce Field\tW. Error\tAcc. W. Error\t%% difference\n");
	}

	percent_diff = fabs(tset->force[i].ff_val);	// not truly a percentage
	acc_percent_diff += percent_diff;

	if (errorfunction == 0)	// SIWE
	    partial_error = pow((tset->force[i].ff_val) / tset->force[i].weight, 2);
	else if (errorfunction == 1)	// MPE (default)
	    partial_error = tset->force[i].nweight * percent_diff;
	else if (errorfunction == 2 || errorfunction == 3)	// RMSE or NRMSE
	    partial_error = tset->force[i].nweight * pow(tset->force[i].ff_val, 2);
	else if (errorfunction == 4)	// MSE
	    partial_error = pow(tset->force[i].weight * diff, 2);

	error[FORCE] += partial_error;
	acc_error += partial_error;

	if (errorfunction == 3) {
	    if (i == 0) {
		Emax = tset->force[i].ff_val;
		Emin = tset->force[i].ff_val;
	    } else if (tset->force[i].ff_val > Emax)
		Emax = tset->force[i].ff_val;
	    else if (tset->force[i].ff_val < Emin)
		Emin = tset->force[i].ff_val;
	}

	if (normalizederrors)
	    fprintf(fp, "%-4.3f\t[%4.3f]\t%-10s\t%-10.4f\t%-10.4f\t[%-10.4f]\t[%-10.4f]\n",
		    tset->force[i].weight, tset->force[i].nweight, tset->force[i].sname,
		    tset->force[i].ff_val, partial_error, acc_error + press, percent_diff);
	else
	    fprintf(fp, "%-4.3f\t%-10s\t%-10.4f\t%-10.4f\t[%-10.4f]\t[%-10.4f]\n",
		    tset->force[i].weight,
		    tset->force[i].sname, tset->force[i].ff_val, partial_error, acc_error + press,
		    percent_diff);
    }
    if (nfit[FORCE] > 0) {
	if (errorfunction == 0)	// SIWE 
	    fprintf(fp,
		    "[%s] Total Inverse Squared Error [GEOMETRY] = %-10.4f (|Acc. %% diff| = %10.4f)\n",
		    errorname, error[GEOMETRY], acc_percent_diff);
	else if (errorfunction == 1 || errorfunction == 4) {	// MPE or MSE
	    fprintf(fp, "[%s] Total Mean (Percent or Squared) Error [FORCE] = %4.2f/%d = ",
		    errorname, error[FORCE], nfit[FORCE]);
	    error[FORCE] = error[FORCE] / nfit[FORCE];
	    fprintf(fp, "%-10.4f (|Acc. %% diff| = %10.4f)\n", error[FORCE], acc_percent_diff);
	} else if (errorfunction == 2) {	// RMSE
	    fprintf(fp, "[%s] Total RMSE [FORCE] = sqrt(%4.2f/%d) = ", errorname, error[FORCE],
		    nfit[FORCE]);
	    error[FORCE] = sqrt(error[FORCE] / nfit[FORCE]);
	    fprintf(fp, "%-10.4f (|Acc. %% diff| = %10.4f)\n", error[FORCE], acc_percent_diff);
	} else if (errorfunction == 3) {	// NRMSE = RMSE/(max-min)
	    if (Emax != Emin)
		Emax_minus_Emin = Emax - Emin;
	    else
		Emax_minus_Emin = 1.0;
	    fprintf(fp,
		    "[%s] Total Normalized RMSE [FORCE] = RMSE/(Emax-Emin) = sqrt(%4.2f/%d)/%4.2f = ",
		    errorname, error[FORCE], nfit[FORCE], Emax_minus_Emin);
	    error[FORCE] = sqrt(error[FORCE] / nfit[FORCE]) / Emax_minus_Emin;
	    fprintf(fp, "%-10.4f (|Acc. %% diff| = %10.4f)\n", error[FORCE], acc_percent_diff);
	}
    }
    // CELL STRESS ERROR CALCULATION
    for (i = 0; i < nfit[STRESS]; i++) {
	if (i == 0) {
	    if (normalizederrors)
		fprintf(fp,
			"CELL STRESS %lf\nWeight\t[nWeight]\tStructure\tQM\tForce Field\tW. Error\tAcc. W. Error\t%% difference\n",
			secweight[STRESS]);
	    else
		fprintf(fp,
			"CELL STRESS\nWeight\tStructure\tQM\tForce Field\tW. Error\tAcc. W. Error\t%% difference\n");
	}

	percent_diff = fabs(tset->stress[i].ff_val);	// not truly a percent diff
	acc_percent_diff += percent_diff;

	if (errorfunction == 0)	// SIWE
	    partial_error = pow((tset->stress[i].ff_val) / tset->stress[i].weight, 2);
	else if (errorfunction == 1)	// MPE (default)
	    partial_error = tset->stress[i].nweight * percent_diff;
	else if (errorfunction == 2 || errorfunction == 3)	// RMSE or NRMSE
	    partial_error = tset->stress[i].nweight * pow(tset->stress[i].ff_val, 2);
	else if (errorfunction == 4)	// MSE
	    partial_error = pow(tset->stress[i].weight * diff, 2);

	error[STRESS] += partial_error;
	acc_error += partial_error;

	if (errorfunction == 3) {
	    if (i == 0) {
		Emax = tset->stress[i].ff_val;
		Emin = tset->stress[i].ff_val;
	    } else if (tset->stress[i].ff_val > Emax)
		Emax = tset->stress[i].ff_val;
	    else if (tset->stress[i].ff_val < Emin)
		Emin = tset->stress[i].ff_val;
	}

	if (normalizederrors)
	    fprintf(fp, "%-4.3f[%4.3f]\t%-10s\t%-10.4f\t%-10.4f\t[%-10.4f]\t[%-10.4f]\n",
		    tset->stress[i].weight, tset->stress[i].nweight,
		    tset->stress[i].sname, tset->stress[i].ff_val, partial_error, acc_error + press,
		    percent_diff);
	else
	    fprintf(fp, "%-4.3f\t%-10s\t%-10.4f\t%-10.4f\t[%-10.4f]\t[%-10.4f]\n",
		    tset->stress[i].weight,
		    tset->stress[i].sname, tset->stress[i].ff_val, partial_error, acc_error + press,
		    percent_diff);
    }
    if (nfit[STRESS] > 0) {
	if (errorfunction == 0)	// SIWE 
	    fprintf(fp,
		    "[%s] Total Inverse Squared Error [STRESS] = %-10.4f (|Acc. %% diff| = %10.4f)\n",
		    errorname, error[STRESS], acc_percent_diff);
	else if (errorfunction == 1 || errorfunction == 4) {	// MPE
	    fprintf(fp, "[%s] Total Mean (Percent or Squared) Error [STRESS] = %4.2f/%d = ",
		    errorname, error[STRESS], nfit[STRESS]);
	    error[STRESS] = error[STRESS] / nfit[STRESS];
	    fprintf(fp, "%-10.4f (|Acc. %% diff| = %10.4f)\n", error[STRESS], acc_percent_diff);
	} else if (errorfunction == 2) {	// RMSE
	    fprintf(fp, "[%s] Total RMSE [STRESS] = sqrt(%4.2f//%d) = ", errorname, error[STRESS],
		    nfit[STRESS]);
	    error[STRESS] = sqrt(error[STRESS] / nfit[STRESS]);
	    fprintf(fp, "%-10.4f (|Acc. %% diff| = %10.4f)\n", error[STRESS], acc_percent_diff);
	} else if (errorfunction == 3) {	// NRMSE = RMSE/(max-min)
	    if (Emax != Emin)
		Emax_minus_Emin = Emax - Emin;
	    else
		Emax_minus_Emin = 1.0;
	    fprintf(fp,
		    "[%s] Total Normalized RMSE [STRESS] = RMSE/(Emax-Emin) = sqrt(%4.2f/%d)/%4.2f = ",
		    errorname, error[STRESS], nfit[STRESS], Emax_minus_Emin);
	    error[STRESS] = sqrt(error[STRESS] / nfit[STRESS]) / Emax_minus_Emin;
	    fprintf(fp, "%-10.4f (|Acc. %% diff| = %10.4f)\n", error[STRESS], acc_percent_diff);
	}
    }

    // CELL PRESSURE ERROR CALCULATION
    for (i = 0; i < nfit[PRESSURE]; i++) {
        if (i == 0) {
            if (normalizederrors)
                fprintf(fp,
                        "CELL PRESSURE %lf\nWeight\t[nWeight]\tStructure\tQM\tForce Field\tW. Error\tAcc. W. Error\t%% difference\n",
                        secweight[PRESSURE]);
            else
                fprintf(fp,
                        "CELL PRESSURE\nWeight\tStructure\tQM\tForce Field\tW. Error\tAcc. W. Error\t%% difference\n");
        }

        percent_diff = fabs(tset->press[i].ff_val);    // not truly a percent diff
        acc_percent_diff += percent_diff;

        if (errorfunction == 0) // SIWE
            partial_error = pow((tset->press[i].ff_val) / tset->press[i].weight, 2);
        else if (errorfunction == 1)    // MPE (default)
            partial_error = tset->press[i].nweight * percent_diff;
        else if (errorfunction == 2 || errorfunction == 3)      // RMSE or NRMSE
            partial_error = tset->press[i].nweight * pow(tset->press[i].ff_val, 2);
        else if (errorfunction == 4)    // MSE
            partial_error = pow(tset->press[i].weight * diff, 2);

        error[PRESSURE] += partial_error;
        acc_error += partial_error;

        if (errorfunction == 3) {
            if (i == 0) {
                Emax = tset->press[i].ff_val;
                Emin = tset->press[i].ff_val;
            } else if (tset->press[i].ff_val > Emax)
                Emax = tset->press[i].ff_val;
            else if (tset->press[i].ff_val < Emin)
                Emin = tset->press[i].ff_val;
        }
        if (normalizederrors)
            fprintf(fp, "%-4.3f[%4.3f]\t%-10s\t%-10.4f\t%-10.4f\t[%-10.4f]\t[%-10.4f]\n",
                    tset->press[i].weight, tset->press[i].nweight,
                    tset->press[i].sname, tset->press[i].ff_val, partial_error, acc_error + press,
                    percent_diff);
        else
            fprintf(fp, "%-4.3f\t%-10s\t%-10.4f\t%-10.4f\t[%-10.4f]\t[%-10.4f]\n",
                    tset->press[i].weight,
                    tset->press[i].sname, tset->press[i].ff_val, partial_error, acc_error + press,
                    percent_diff);
    }
    if (nfit[PRESSURE] > 0) {
        if (errorfunction == 0) // SIWE 
            fprintf(fp,
                    "[%s] Total Inverse Squared Error [PRESSURE] = %-10.4f (|Acc. %% diff| = %10.4f)\n",
                    errorname, error[PRESSURE], acc_percent_diff);
        else if (errorfunction == 1 || errorfunction == 4) {    // MPE
            fprintf(fp, "[%s] Total Mean (Percent or Squared) Error [PRESSURE] = %4.2f/%d = ",
                    errorname, error[PRESSURE], nfit[PRESSURE]);
            error[PRESSURE] = error[PRESSURE] / nfit[PRESSURE];
            fprintf(fp, "%-10.4f (|Acc. %% diff| = %10.4f)\n", error[PRESSURE], acc_percent_diff);
        } else if (errorfunction == 2) {        // RMSE
            fprintf(fp, "[%s] Total RMSE [PRESSURE] = sqrt(%4.2f//%d) = ", errorname, error[PRESSURE],
                    nfit[PRESSURE]);
            error[PRESSURE] = sqrt(error[PRESSURE] / nfit[PRESSURE]);
            fprintf(fp, "%-10.4f (|Acc. %% diff| = %10.4f)\n", error[PRESSURE], acc_percent_diff);
        } else if (errorfunction == 3) {        // NRMSE = RMSE/(max-min)
            if (Emax != Emin)
                Emax_minus_Emin = Emax - Emin;
            else
                Emax_minus_Emin = 1.0;
            fprintf(fp,
                    "[%s] Total Normalized RMSE [PRESSURE] = RMSE/(Emax-Emin) = sqrt(%4.2f/%d)/%4.2f = ",
                    errorname, error[PRESSURE], nfit[PRESSURE], Emax_minus_Emin);
            error[PRESSURE] = sqrt(error[PRESSURE] / nfit[PRESSURE]) / Emax_minus_Emin;
            fprintf(fp, "%-10.4f (|Acc. %% diff| = %10.4f)\n", error[PRESSURE], acc_percent_diff);
        }
    }


    // ENERGY ERROR CALCULATION
    for (i = 0; i < nfit[ENERGY]; i++) {
	eval = 0.0;
	if (i == 0) {
	    if (normalizederrors)
		fprintf(fp, "ENERGY %lf\nWeight\t[nWeight]\tStructures                   ",
			secweight[ENERGY]);
	    else
		fprintf(fp, "ENERGY\nWeight\tStructures                   ");
	    fprintf(fp, "QM\t\tForce Field\tW. Error\tAcc. W. Error\t%% difference\n");
	}

	if (normalizederrors)
	    fprintf(fp, "%-4.3f\t[%4.3f] ", tset->eng[i].weight, tset->eng[i].nweight);
	else
	    fprintf(fp, "%-4.3f\t", tset->eng[i].weight);

	for (j = 0; j < tset->eng[i].n; j++) {
	    id = get_id_from_name(tset->eng[i].sname[j]);
	    if (id == -1)
		continue;
	    if (strcmp(tset->eng[i].op[j], "+") == 0)
		eval += energies[id] / tset->eng[i].factor[j];
	    else
		eval -= energies[id] / tset->eng[i].factor[j];
	    fprintf(fp, "%-1s%-10s/%-2d", tset->eng[i].op[j], tset->eng[i].sname[j],
		    (int) tset->eng[i].factor[j]);
	}

	diff = eval - tset->eng[i].eng;
	if (tset->eng[i].eng != 0)
	    percent_diff = 100 * fabs(diff / tset->eng[i].eng);
	else
	    percent_diff = 0.0;
	acc_percent_diff += percent_diff;

	if (errorfunction == 0)	// SIWE
	    partial_error = pow(diff / tset->eng[i].weight, 2);
	else if (errorfunction == 1)	// MPE (default)
	    partial_error = tset->eng[i].nweight * percent_diff;
	else if (errorfunction == 2 || errorfunction == 3)	// RMSE or NRMSE
	    partial_error = tset->eng[i].nweight * pow(diff, 2);
	else if (errorfunction == 4)	// MSE
	    partial_error = pow(tset->eng[i].weight * diff, 2);

	error[ENERGY] += partial_error;
	acc_error += partial_error;

	if (errorfunction == 3) {
	    if (i == 0) {
		Emax = eval;
		Emin = eval;
	    } else if (eval > Emax)
		Emax = eval;
	    else if (eval < Emin)
		Emin = eval;
	}

	fprintf(fp, " %-10.4f\t%-10.4f\t%-10.4f\t%-10.4f\t%-10.4f\n",
		tset->eng[i].eng, eval, partial_error, acc_error + press, percent_diff);
    }
    if (nfit[ENERGY] > 0) {
	if (errorfunction == 0)	// SIWE 
	    fprintf(fp,
		    "[%s] Total Inverse Squared Error [ENERGY] (%d entries) = %-10.4f (|Acc. %% diff| = %10.4f)\n",
		    errorname, nfit[ENERGY], error[ENERGY], acc_percent_diff);
	if (errorfunction == 1 || errorfunction == 4) {	// MPE
	    fprintf(fp, "[%s] Total Mean (Percent or Squared) Error [ENERGY] = %4.2f/%d = ",
		    errorname, error[ENERGY], nfit[ENERGY]);
	    error[ENERGY] = error[ENERGY] / nfit[ENERGY];
	    fprintf(fp, "%-10.4f (|Acc. %% diff| = %10.4f)\n", error[ENERGY], acc_percent_diff);
	}
	if (errorfunction == 2) {	// RMSE
	    fprintf(fp, "[%s] Total RMSE [ENERGY] = sqrt(%4.2f/%d) = ", errorname, error[ENERGY],
		    nfit[ENERGY]);
	    error[ENERGY] = sqrt(error[ENERGY] / nfit[ENERGY]);
	    fprintf(fp, "%-10.4f (|Acc. %% diff| = %10.4f)\n", error[ENERGY], acc_percent_diff);
	}
	if (errorfunction == 3) {	// NRMSE = RMSE/(max-min)
	    if (Emax != Emin)
		Emax_minus_Emin = Emax - Emin;
	    else
		Emax_minus_Emin = 1.0;
	    fprintf(fp,
		    "[%s] Total Normalized RMSE [ENERGY] = RMSE/(Emax-Emin) = sqrt(%4.2f/%d)/%4.2f = ",
		    errorname, error[ENERGY], nfit[ENERGY], Emax_minus_Emin);
	    error[ENERGY] = sqrt(error[ENERGY] / nfit[ENERGY]) / Emax_minus_Emin;
	    fprintf(fp, "%-10.4f (|Acc. %% diff| = %10.4f)\n", error[ENERGY], acc_percent_diff);
	}
    }

    fprintf(fp, "Total iteration error: %4.4f\n", acc_error);
    fflush(fp);
    fclose(fp);

    total_error = 0.0;
    for (i = 0; i < 11; i++)
	total_error += error[i];

    return total_error;
}

/****************************************************************************
*  Function     : get_id_from_name
*  Description  : Returns the index entry corresponding to the molecular 
                  structure name passed as an argument
*  Parameters   : Molecular structure name
*  Effects      : Returns integer id
****************************************************************************/

int get_id_from_name(char *sname)
{
    int i, id;

    id = -1;
    for (i = 0; i < num; i++) {
	if (strcmp(fname[i], sname) == 0) {
	    id = i;
	    break;
	}
    }
    return id;
}

/****************************************************************************
*  Function     : readRDFs
*  Description  : Used for CG force field trained via radial distribution 
                  functions (RDFs). Reads and parses an RDF for every pairwise
                  interaction defined in the user-defined training set.  
                  The file names of the RDFs are taken from the parsed
                  training set values stored in the tset -> sname data structure.
                  One entry per RDF data point, each with an associated sname. 
                  Names are hardcoded in vmdRDF.pl, must have same names in target file.
 * struc has one entry per RDF data point, each of which
 * has an associated sname.  So what we do here is increment
 * i each time we read a line out of a file.  Once we get to the end
 * of a file, we go back to the top and read the next file.
 * We have to ensure that the FF RDFs exactly match the target RDFs
*  Parameters   : None
*  Effects      :  
****************************************************************************/

void readRDFs(void)
{
    int i = 0;
    char str[MAX_LINE];
    char **strtmp;
    int nwords, k;

    // allocate strtmp
    strtmp = (char **) scalloc(MAX_TOKENS, sizeof(char *), "tmp3");
    for (k = 0; k < MAX_TOKENS; k++)
	strtmp[k] = (char *) scalloc(MAX_TOKEN_LEN, sizeof(char), "tmp3i");

    while (i < nfit[STRUCTURE]) {
	FILE *fname = fopen(tset->struc[i].sname, "r");
	if (fname == NULL)
	    printf("WARNING: fopen failed for %s\n", tset->struc[i].sname);
	while (!feof(fname)) {
	    fgets(str, MAX_LINE, fname);
	    nwords = tokenize_string(str, &strtmp);
	    tset->struc[i++].ff_val = atof(strtmp[1]);
	}
	fclose(fname);
    }

    for (k = 0; k < MAX_TOKENS; k++)
	free(strtmp[k]);
    free(strtmp);
}

/****************************************************************************
*  Function     : getPress
*  Description  : Reads and parses the LAMMPS calculated system pressure from
                  the log.lammps file and writes it to p.txt
*  Parameters   : None
*  Effects      : Writes file p.txt containing pressue entry value
****************************************************************************/

double getPress(void)
{
    double ans = 0;
    char str[100];
    int ctr = 0;

    system("./getPressure.pl log.lammps p.txt");
    FILE *fname = fopen("p.txt", "r");
    if (fname == NULL)
	printf("WARNING: fopen failed for p.txt\n");
    while (!feof(fname)) {
	++ctr;
	fgets(str, MAX_LINE, fname);
	ans += fabs(atof(str));
    }
    fclose(fname);
    ans /= ctr;
    return ans;
}
