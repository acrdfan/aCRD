// Transient Diffusion Green's Function Solver
// In MOD_SOL_aCRD
// Version 1 alpha demo
// J. Fan 2018-2019

#include "stdafx.h"
#include <stdlib.h>
#include <math.h>
#include <string>
#include "TransientDiffSolver.h"

using namespace std;

// Solver Version
TRANSIENTDIFFSOLVER_API double nE02SolverVersion = 1.0;

// The Definition of Regular Parameters
#define	PI 3.141592653589793238
#define eps 1e-16

location getTransDiffRelativeLocation1D(double Dv, double timestamp, double domainRadius, double walkerRadius, location centerLocation)
{
	double RANDOM_SEED;
	double center_cover, random_location;
	RANDOM_SEED = rand() / (double)(RAND_MAX); // Random SEED [0,1)
	center_cover = 2 * (domainRadius - walkerRadius);
	random_location = center_cover*RANDOM_SEED;
	centerLocation.x = centerLocation.x - domainRadius + walkerRadius + random_location;
	return centerLocation;
}

location getTransDiffRelativeLocation3D(double Dv, double timestamp, double domainRadius, double walkerRadius, location centerLocation)
{
	double RANDOM_SEED = 0.0;
	double RANDOM_SEED_THETA = 0.0;
	double RANDOM_SEED_PHI = 0.0;
	double center_shift;
	RANDOM_SEED = rand() / (double)(RAND_MAX); // Random SEED [0,1)
	center_shift = (domainRadius - walkerRadius)*RANDOM_SEED;
	RANDOM_SEED_THETA = rand() / (double)(RAND_MAX); // Random SEED [0,1), 2*PI*R for radians
	RANDOM_SEED_PHI = rand() / (double)(RAND_MAX); // Random SEED [0,1), 2*PI*R for radians
	centerLocation.x = centerLocation.x + center_shift*sin(2 * PI*RANDOM_SEED_THETA)*cos(2 * PI*RANDOM_SEED_PHI);
	centerLocation.y = centerLocation.y * center_shift*sin(2 * PI*RANDOM_SEED_THETA)*sin(2 * PI*RANDOM_SEED_PHI);
	centerLocation.z = centerLocation.z * center_shift*cos(2 * PI*RANDOM_SEED_THETA);
	return centerLocation;
}
