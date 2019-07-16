#pragma once
// Version 1 alpha demo
// This file related to aCRD/Main.cpp
// This file defines several constants, parameters, structs and virtual functions
// This file is including several head files under MS Windows OS 
#ifndef MAIN_H
#define MAIN_H

#include <stdio.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string>
#include <fstream>
#include <map>
#include <vector>
#include <queue>
#include <iomanip>
#include <Windows.h>
using namespace std;

#define MAJOR_VERSION_NUMBER 1
#define MINOR_VERSION_NUMBER 0
#define REVISIONB_NUMBER 4
#define VERSION_ID "Alpha"

// The Definition of Regular Parameters
#define	PI 3.141592653589793238
#define eps 1e-16 
// Max range of one domain, unit in nm
#define EXTRA_LARGE 4999.99

// The Definition of Main Structs
struct location {
	double x;
	double y;
	double z;
};

// This structure is designed for SPD_a
struct domaincenter {
	double x;
	double y;
	double z;
	double r;
};

// The walkers follow domain time (protected) or global time (released)
struct walker {
	long ID;
	int type; // 0 V_Based; 1 I_Based; -1 for Not Using
	// int subtype; // For different defect structures, not in use for current version.
	long component; // Cluster number; 1 for single walker
	location center;
	double radius; // Use modifier for vacancies
	double energy; // In keV
	double diffusivity; // In nm^2/s
	int statue; // 0 Released; 1 Soft Protected; 2 Hard Protected 
	long regionID;
	double timestamp; // WORLD
};

// For current version, only pre-assigned final status: 1 - 5 walkers
struct dissociation_final_walker {
	walker member[5];
};

// For current version, only pre-assigned final status: 1 - 50 walkers
struct reaction_current_walker {
	walker member[50];
};

// For current version, only pre-assigned final status: 1 - 50 walkers
// TODO
struct SPDb_current_walker {
	walker member[50];
};

// For current version, only pre-assigned final status: 1 - 5 walkers
struct reaction_final_walker {
	walker member[5];
};

// For 3 Domain Design, Not for Precusor Design
// Type 1 PD, for First-Passage Event, Sealed Boundary, Single walker, Event-Driven
// Include E01, E11, and E02 from outside
struct hard_domain {
	long ID;
	double clock; // Inner timestamp
	long walker_ID;
	int walker_type;
	double radius;
	location center;
	int events_type; // Type 01 or Type 11, -1 for no event (disable)
	double events_time; // Event timestamp, relative to the initial time of domain
	int event_card; // To the index of event card in EVENT_LIST
	location e01_final_location;
	dissociation_final_walker e11_final_statue; // Walker ID will be assigned only when they are released
};

// Type 2 PD, for Local Event, None-Releasing Boundary, Multi Walker, Event-Driven
// Include E10 only, can be modified by new members
// Two overlap sphere or line centered at two key walker respectively, recording two center and two radius, outside with THRES_2
struct soft_domain_a {
	long ID;
	double clock; // Built WORLD
	int current_member;
	int final_member;
	reaction_current_walker member_list;
	domaincenter center1; // 1D: .x for location, .r for radius; 3D: .x .y .z for location, .r for radius
	domaincenter center2;
	double center_distance; // Center walkers surface distance, always less than THRES_2
	int events_type; // Type 10, -1 for no event (Disable) 
	double events_time; // Event timestamp, relative to the initial time of domain
	int event_card; // To the index of event card in EVENT_LIST
	reaction_final_walker e10_final_statue; // Walker ID will be assigned only when they are released
};

// Type 3 PD, for Local Event, Open Boundary, Multi Walker, Hybrid-Driven
// Include E10 and E20+
// TODO
struct soft_domain_b {
	long ID;
	double clock;
	int current_member;
	SPDb_current_walker member_list;
	double radius;
	location center;
	int events_type; // Type 10 or Type 20+, -1 for no event (Disable) 
	double timestep; // Hybrid factor, can be global, unit in s
};

struct event_card {
	long domainID;
	int domainType; // 0 for HPD, 1 for SPDa, 2 for SPDb
	int eventType; // Only for New Walker, 0 for Domain Event, 20 for New Walker
	double timestamp; // This is world clock stamp, not the relative stamp for solvers
};

// Event Queue
event_card EVENT_LIST[200000];
long EVENT_INDEX = 0; // For initialization
long EVENT_POINTER = 1; // For event execulation
long EVENTS_LENGTH = 0; // For sorting
// queue<long> EVENT_LIST;

long S0_INDEX = 1;
long S0_EVENT[5000]; // Max S0 at the same time, for walker ID
int S0_FULL_FLAG;

// Preassigned for test version
// Conside auto assign and memory manage in beta version: TODO
#define	INPUTS_SIZE 20000
#define DOMAIN_SIZE 100000

walker ORI_WALKERS[INPUTS_SIZE];
hard_domain HARD_PD[DOMAIN_SIZE];
soft_domain_a SOFT_PD_A[DOMAIN_SIZE];
// struct soft_domain_b SOFT_PD_B[DOMAIN_SIZE];

// Other Parameters
long i, j, k;
int inputNumber = 0;
int space_dimension = 0;

// E20
int irradiationCheck = 0;
double irradiationPulse = 0.0;
int irradiationINumber = 0;
int irradiationVNumber = 0;

double new_walker_distance = 0.0;

double X_LOWWER_LIMIT = 0.0;
double X_UPPER_LIMIT = 0.0;
double Y_LOWWER_LIMIT = 0.0;
double Y_UPPER_LIMIT = 0.0;
double Z_LOWWER_LIMIT = 0.0;
double Z_UPPER_LIMIT = 0.0;

double WORLD_CLOCK = 0.0; // Point to the last exist S0 clock
double END_CLOCK = 1000.00; // Unit in s
double EXTRA_CLOCK = 0.0; // Add to END_CLOCK
long TIME_STEP = 0;

long WALKER_INDEX = 0; // Point to the last one, storage from 1 to N
long HPDOMAIN_INDEX = 0; // Point to the last one, storage from 1 to N
long SPDOMAINA_INDEX = 0; // Point to the last one, stotage from 1 to N
long SPDOMAINB_INDEX = 0;

long WALKER_COUNTER = 0;
long HPD_COUNTER = 0;
long SPDA_COUNTER = 0;
long SPDB_COUNTER = 0;
long DUMP_COUNTER = 0;

double RANDOM_SEED = 0.0;
double TEMP_TIME = 0.0;
int CONTROL_FLAG;
int STAGE_FLAG;
int INPUT_FLAG;
int LOCATE_FLAG;
long NEARID;
double NEARDISTANCE = 0.0;
string filenameInput;
string OUTPUTNAME;
string STAMPNAME;
int DISPLAYFLAG_SYS = 1;
int DISPLAYFLAG_DEB = 0;
int DUMPFLAG_INN = 1;
int DUMPLAMMPSFLAG_INN = 1;
int DUMPCOUNTER_FLAG = 1;
int DUMPCOUNTER_INN = 0;
string ININOTES;

// Pure virtual function member list
// DO NOT include the functions from other module

void releaseModule();
// Toolbox
double getDistance(walker WALKER_1, walker WALKER_2);

double getDistance1D(walker WALKER_1, walker WALKER_2);

double getSoftDomainASurfaceDistance(walker WALKER_NEW, soft_domain_a SPDA);

double getSoftDomainASurfaceDistance1D(walker WALKER_NEW, soft_domain_a SPDA);

double getSoftDomainACentertoCenterDistance(walker WALKER_NEW, soft_domain_a SPDA);

double getSoftDomainACentertoCenterDistance1D(walker WALKER_NEW, soft_domain_a SPDA);

double getHardDomainCenterDistance(hard_domain DOMAIN_1, hard_domain DOMAIN_2);

double getEstimateRadius(walker WALKER);

long getNeighborDomainID(location DOMAIN_CENTER);

double getHardDomainSurfaceDistance1D(walker WALKER, hard_domain HPD);

double getHardDomainSurfaceDistance3D(walker WALKER, hard_domain HPD);
// Main
void checkWalkerOverlapping();

void initialSoftDomainA3D(long ID, walker WALKERS_1, walker WALKERS_2);

void initialUpdateSoftDomainA3D(long ID, walker WALKERS_NEW);

void checkInitialSoftDomains3D();

void initialSoftDomainA1D(long ID, walker WALKERS_1, walker WALKERS_2);

void initialUpdateSoftDomainA1D(long ID, walker WALKERS_NEW);

void checkInitialSoftDomains1D();

void initialHardDomain3D(long ID, walker WALKERS, double DOMAIN_R);

void setInitialHardDomains3D();

void initialHardDomain1D(long ID, walker WALKERS, double DOMAIN_R);

void setInitialHardDomains1D();

double getHardDomainModifier(hard_domain DOMAIN_1, hard_domain DOMAIN_2);

void checkInitialDomain1D();

void checkInitialDomain3D();

int recheckInitialDomain();

void buildHardDomain(long ID, walker WALKERS);

void buildSoftDomainA(long ID, walker WALKERS_1, walker WALKERS_2);

void updateSoftDomainA(long ID, walker NEW_WALKER);

// Virtual dll list
HMODULE dllLibDiffusivity;
HMODULE dllSolverMajorDiff;
HMODULE dllSolverTransientDiff;
HMODULE dllSolverLocal;
HMODULE dllSolverDissociation;

#endif