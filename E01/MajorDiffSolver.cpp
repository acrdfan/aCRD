// Major Diffusion Green's Function Solver
// In MOD_SOL_aCRD
// Version 1 alpha demo
// J. Fan 2018-2019

#include "stdafx.h"
#include <stdlib.h>
#include <math.h>
#include <string>
#include "MajorDiffSolver.h"

using namespace std;

// Solver Version
// This is the Function for t_p, S(t_p) = (RANDOM_SEED)
// Also j(r,t_p). In this Version, use other Random Seed(s) as the Hypothesis of Even Distribution of Diffusion
// Need Hyper-Rectangle Method in Future Version: TODO
MAJORDIFFSOLVER_API double nE01SolverVersion = 1.0;

// The Definition of Regular Parameters
#define	PI 3.141592653589793238
#define eps 1e-16

double getMajorDiffTimeStamp1D(double Dv, double domainRadius, double walkerRadius)
{
	double FP_TIME = 0.0;
	double TIME1 = 0.0, TIME2 = 0.0;
	double reduce_l = 0.0; // L = domainRaidus
	double reduce_t = 0.0; // t = L * L / Dv
	double RANDOM_SEED;
	double S_attempt;
	double tempt;
	double check;
	double limit1;
	double limit2;
	int flagt = 1;
	S_attempt = 0.0;
	limit1 = 0.243; // > 1.0
	limit2 = 2.000; // 5.28E-8
	reduce_l = domainRadius - walkerRadius; // Adjust unit
	reduce_t = reduce_l*reduce_l / Dv; // Adjust unit
	RANDOM_SEED = rand() / (double)(RAND_MAX); // Random SEED [0,1)
	flagt = 1;
	while (flagt) {
		tempt = (limit1 + limit2) / 2.0;
		S_attempt = 2 * PI*PI*exp(-PI*PI*tempt); // Piecewise Smooth Function Approach, with Long t_p Assumption
		check = abs(S_attempt - RANDOM_SEED) / RANDOM_SEED;
		if (check - 0.0000001 < eps) {
			flagt = 0;
		}
		else
		{
			if (S_attempt - RANDOM_SEED < eps) {
				limit2 = tempt;
			}
			else
			{
				limit1 = tempt;
			}
		}
	}
	FP_TIME = tempt*reduce_t;
	return FP_TIME;
}

double getMajorDiffTimeStamp3D(double Dv, double domainRadius, double walkerRadius)
{
	double FP_TIME = 0.0;
	double TIME1 = 0.0, TIME2 = 0.0;
	double reduce_l = 0.0; // L = domainRaidus
	double reduce_t = 0.0; // t = L * L / Dv
	double RANDOM_SEED;
	double S_attempt;
	double tempt;
	double check;
	double limit1;
	double limit2;
	int flagt = 1;
	S_attempt = 0.0;
	limit1 = 0.243; // > 1.0
	limit2 = 2.000; // 5.28E-8
	reduce_l = domainRadius - walkerRadius; // Adjust unit
	reduce_t = reduce_l*reduce_l / Dv; // Adjust unit
	RANDOM_SEED = 5.30e-8 + rand() / (double)(RAND_MAX); // Random SEED (0,1)
	flagt = 1;
	while (flagt) {
		tempt = (limit1 + limit2) / 2.0;
		S_attempt = 2 * PI*PI*exp(-PI*PI*tempt); // Piecewise Smooth Function Approach, with Long t_p Assumption
		check = abs(S_attempt - RANDOM_SEED) / RANDOM_SEED;
		if (check - 0.0000001 < eps) {
			flagt = 0;
		}
		else
		{
			if (S_attempt - RANDOM_SEED < eps) {
				limit2 = tempt;
			}
			else
			{
				limit1 = tempt;
			}
		}
	}
	FP_TIME = tempt*reduce_t;
	return FP_TIME;
}

location getMajorDiffRelativeLocation1D(double Dv, double timestamp, double domainRadius, double walkerRadius, location centerLocation)
{
	double RANDOM_SEED;
	RANDOM_SEED = rand() / (double)(RAND_MAX); // Random SEED [0,1)
	if (RANDOM_SEED - 0.5 < eps) {
		centerLocation.x = centerLocation.x - domainRadius + walkerRadius;
	}
	else {
		centerLocation.x = centerLocation.x + domainRadius - walkerRadius;
	}
	centerLocation.x = centerLocation.x * 2.0;
	centerLocation.y = 0.0;
	centerLocation.z = 0.0;
	return centerLocation;
}

location getMajorDiffRelativeLocation3D(double Dv, double timestamp, double domainRadius, double walkerRadius, location centerLocation)
{
	double RANDOM_SEED_THETA = 0.0;
	double RANDOM_SEED_PHI = 0.0;
	RANDOM_SEED_THETA = rand() / (double)(RAND_MAX); // Random SEED [0,1), 2*PI*R for radians
	RANDOM_SEED_PHI = rand() / (double)(RAND_MAX); // Random SEED [0,1), 2*PI*R for radians
	centerLocation.x = centerLocation.x + domainRadius*sin(2 * PI*RANDOM_SEED_THETA)*cos(2 * PI*RANDOM_SEED_PHI) - walkerRadius;
	centerLocation.y = centerLocation.y + domainRadius*sin(2 * PI*RANDOM_SEED_THETA)*sin(2 * PI*RANDOM_SEED_PHI) - walkerRadius;
	centerLocation.z = centerLocation.z + domainRadius*cos(2 * PI*RANDOM_SEED_THETA) - walkerRadius;
	return centerLocation;
}
