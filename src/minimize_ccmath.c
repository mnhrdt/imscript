#ifndef _MINIMIZE_CCMATH
#define _MINIMIZE_CCMATH

#include "optmiz.c"

// API
typedef double (objective_function)(double *x, int n, void *data);

// API
int minimize_objective_function(double *result, double *first, double *step,
		objective_function *f, int n, void *data);


#include <assert.h>
#include <stdbool.h>
#include <stdio.h>



#include "smapa.h"

SMART_PARAMETER(CCMATH_MIN_MAXITER,500)
SMART_PARAMETER(CCMATH_MIN_DE,1e-5)
SMART_PARAMETER(CCMATH_MIN_DELTA,1e-7)

static int global_variable_for_ccmath_optimization_n;
static void *global_variable_for_ccmath_optimization_data;
static objective_function *global_variable_for_ccmath_optimization_f;

static double func(double *x)
{
	int n = global_variable_for_ccmath_optimization_n;
	void *data = global_variable_for_ccmath_optimization_data;
	objective_function *f = global_variable_for_ccmath_optimization_f;
	return f(x, n, data);
}


// external interface to CCMATH code
// (observation: "step" is ignored)
int minimize_objective_function(double *result, double *first, double *step,
		objective_function *f, int n, void *data)
{
	global_variable_for_ccmath_optimization_n = n;
	global_variable_for_ccmath_optimization_data = data;
	global_variable_for_ccmath_optimization_f = f;

	int maxiter = CCMATH_MIN_MAXITER();
	double test = CCMATH_MIN_DELTA();
	double de = CCMATH_MIN_DE();
	for (int i = 0; i < n; i++)
		result[i] = first[i];
	return optmiz(result, n, func, de, test, maxiter);
}
#endif//_MINIMIZE_CCMATH
