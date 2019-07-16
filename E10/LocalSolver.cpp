// Local Reaction Possibility Solver
// In MOD_SOL_aCRD
// Version 1 alpha demo
// This solver is using an "adding-up" estimation to local events, you may regard this version as a demonstration to functions' interface
// J. Fan 2018-2019

#include "stdafx.h"
#include <stdlib.h>
#include <math.h>
#include <string>
#include "LocalSolver.h"

using namespace std;

// Solver Version
// This solver is using t_built+FIX_EVENTS_TIME*(1+RANDOM_SEED) as t_event
// This solver only concerns coalescence, with 1 final walker in random center
// This solver doesn't provide SPDa time update interface, only member update in main() 
LOCALSOLVER_API double nE10SolverVersion = 1.0;

// The Definition of Regular Parameters
#define	PI 3.141592653589793238
#define eps 1e-16
// This is only for first initialization estimation
// #define INITIAL_EVENTS_TIME 0.001

double getReactionTimeStamp1D(reaction_current_walker walkers, double initial_time) {
	double RX_TIME = 0.0;
	double RANDOM_SEED;
	// SPDa will trigger within INITIAL_EVENTS_TIME + INITIAL_EVENTS_TIME * Random[0,1)
	RANDOM_SEED = rand() / (double)(RAND_MAX);
	RX_TIME = initial_time + initial_time*RANDOM_SEED;
	return RX_TIME;
}

double getReactionTimeStamp3D(reaction_current_walker walkers, double initial_time) {
	double RX_TIME = 0.0;
	double RANDOM_SEED;
	// SPDa will trigger within INITIAL_EVENTS_TIME + INITIAL_EVENTS_TIME * Random[0,1)
	RANDOM_SEED = rand() / (double)(RAND_MAX);
	RX_TIME = initial_time + initial_time*RANDOM_SEED;
	return RX_TIME;
}

// Annihilation only for this version
reaction_final_walker getReactionFinalStatues1D(reaction_current_walker walkers, int init_member, domaincenter centerA, domaincenter centerB) {
	reaction_final_walker final_walkers;
	final_walkers.member[1].type = -1;
	final_walkers.member[2].type = -1;
	final_walkers.member[3].type = -1;
	final_walkers.member[4].type = -1;
	final_walkers.member[5].type = -1;
	return final_walkers;
}

// Coalescence to 1 or Annihilation in this version
// Need walker type, component, center
reaction_final_walker getReactionFinalStatues3D(reaction_current_walker walkers, int init_member, domaincenter centerA, domaincenter centerB) {
	reaction_final_walker final_walkers;
	int i_counter;
	int v_counter;
	long i_all_comp;
	long v_all_comp;
	int i;
	double RANDOM_SEED = 0.0;
	i_counter = 0;
	v_counter = 0;
	i_all_comp = 0;
	v_all_comp = 0;
	for (i = 0;i < init_member;i++) {
		if (walkers.member[i + 1].type == 0) {
			v_counter = v_counter + 1;
			v_all_comp = v_all_comp + walkers.member[i + 1].component;
		}
		if (walkers.member[i + 1].type == 1) {
			i_counter = i_counter + 1;
			i_all_comp = i_all_comp + walkers.member[i + 1].component;
		}
	}
	if (i_all_comp == v_all_comp) {
		// Annihilation
		final_walkers.member[1].type = -1;
		final_walkers.member[2].type = -1;
		final_walkers.member[3].type = -1;
		final_walkers.member[4].type = -1;
		final_walkers.member[5].type = -1;
	}
	else {
		if (i_all_comp < v_all_comp) {
			// Coalescence to vacancy cluster
			final_walkers.member[1].type = 0;
			final_walkers.member[1].component = v_all_comp - i_all_comp;
			// Select one of SPDa center as final center, try best not to break other HPD
			RANDOM_SEED = rand() / (double)(RAND_MAX); // Random SEED [0,1)
			if (RANDOM_SEED - 0.5 < eps) {
				final_walkers.member[1].center.x = centerA.x;
				final_walkers.member[1].center.y = centerA.y;
				final_walkers.member[1].center.z = centerA.z;
			}
			else {
				final_walkers.member[1].center.x = centerB.x;
				final_walkers.member[1].center.y = centerB.y;
				final_walkers.member[1].center.z = centerB.z;
			}
			final_walkers.member[2].type = -1;
			final_walkers.member[3].type = -1;
			final_walkers.member[4].type = -1;
			final_walkers.member[5].type = -1;
		}
		else {
			// Coalescence to interstitial cluster
			final_walkers.member[1].type = 1;
			final_walkers.member[1].component = i_all_comp - v_all_comp;
			// Select one of SPDa center as final center, try best not to break other HPD
			RANDOM_SEED = rand() / (double)(RAND_MAX); // Random SEED [0,1)
			if (RANDOM_SEED - 0.5 < eps) {
				final_walkers.member[1].center.x = centerA.x;
				final_walkers.member[1].center.y = centerA.y;
				final_walkers.member[1].center.z = centerA.z;
			}
			else {
				final_walkers.member[1].center.x = centerB.x;
				final_walkers.member[1].center.y = centerB.y;
				final_walkers.member[1].center.z = centerB.z;
			}
			final_walkers.member[2].type = -1;
			final_walkers.member[3].type = -1;
			final_walkers.member[4].type = -1;
			final_walkers.member[5].type = -1;
		}
	}
	return final_walkers;
}