// Dissociation Solver
// In MOD_SOL_aCRD
// Version 1 alpha demo
// This version should be regarded as a demonstration of function's interface.
// Reverse for further development 
// J. Fan 2019

#include "stdafx.h"
#include <stdlib.h>
#include <math.h>
#include <string>
#include "DissociationSolver.h"

using namespace std;

// Solver Version
// This solver is only returning dissociation flag
DISSOCIATIONSOLVER_API double nE11SolverVersion = 1.0;

// The Definition of Regular Parameters
#define	PI 3.141592653589793238
#define eps 1e-16

int getDissociationFlag(walker WALKER) {
	double p_dissociate = 0.0;
	double RANDOM_SEED;
	RANDOM_SEED = rand() / (double)(RAND_MAX); // Random SEED [0,1)
	p_dissociate = 1 - 1 / pow(WALKER.component, 0.5);
	if (RANDOM_SEED - p_dissociate<eps) {
		return 1;
	}
	else {
		return 0;
	}
}
