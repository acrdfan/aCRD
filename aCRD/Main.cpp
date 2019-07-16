// Version 1 alpha demo
// Please note, all code related to this prototype is under GNUv3 license
// Please refer to related documents for more details
// J. Fan 2018-2019
#include "Main.h"
#include "Parameters.h"
#include "Toolbox.h"
using namespace std;

// The debugging namespace is NOT included in release version
namespace currentEvent {
	long domainID;
	int domain_type;
	double timestamp;
	int event_type;
	location afterward_location;
	int flag_s0_in_SPDa;
	double walker_distance;
	double min_distance;
	double neighbour_allow_r;
	double O1_allow_r;
	double O1_neighour_allow_r;
	long O1_neighbout_ID;
	double O2_allow_r;
	double O3_allow_r;
	long near_domain_ID;
	int near_domain_type; // 0 HPD; 1 SPDa
	double update_spda_time_temp;
}

// The Walker Overlapping Check is specific for Initialization, if Overlap, delete one of the walker
// This procedure will only be run for once.
// This function can be applied in both 1D and 3D
void checkWalkerOverlapping() {
	double walkerDistance;
	double walkerOccupation;
	for (i = 0;i < inputNumber; i++) {
		for (j = i + 1;j < inputNumber; j++) {
			walkerDistance = getDistance(ORI_WALKERS[i + 1], ORI_WALKERS[j + 1]);
			walkerOccupation = ORI_WALKERS[i + 1].radius + ORI_WALKERS[j + 1].radius;
			if (walkerDistance - walkerOccupation < eps) {
				ORI_WALKERS[i + 1].type = -1;
				cout << "Find " << i + 1 << "th Walker Overlap with " << j + 1 << "th Walker, " << i + 1 << "th Walker has been noted as Not Using (-1)." << endl;
				WALKER_COUNTER = WALKER_COUNTER - 1;
				break;
			}
		}
	}
}

// For new SPDa
// Do not call solver
// Solver should be called only in main()
void initialSoftDomainA3D(long ID, walker WALKERS_1, walker WALKERS_2) {
	double surfaceDistance;
	cout << "Initialize SPDa for " << WALKERS_1.ID << " and " << WALKERS_2.ID << " walkers." << endl;
	cout << "Center 1 at (" << WALKERS_1.center.x << ", " << WALKERS_1.center.y << ", " << WALKERS_1.center.z << "), with raidus " << WALKERS_1.radius << endl;
	cout << "Center 2 at (" << WALKERS_2.center.x << ", " << WALKERS_2.center.y << ", " << WALKERS_2.center.z << "), with raidus " << WALKERS_2.radius << endl;
	surfaceDistance = getDistance(WALKERS_1, WALKERS_2) - WALKERS_1.radius - WALKERS_2.radius;
	SOFT_PD_A[ID].ID = ID;
	SOFT_PD_A[ID].center_distance = surfaceDistance;
	SOFT_PD_A[ID].clock = 0.0;
	SOFT_PD_A[ID].events_type = 10;
	SOFT_PD_A[ID].current_member = 2;
	SOFT_PD_A[ID].member_list.member[1] = WALKERS_1;
	SOFT_PD_A[ID].member_list.member[2] = WALKERS_2;
	SOFT_PD_A[ID].center1.x = WALKERS_1.center.x;
	SOFT_PD_A[ID].center1.y = WALKERS_1.center.y;
	SOFT_PD_A[ID].center1.z = WALKERS_1.center.z;
	SOFT_PD_A[ID].center1.r = WALKERS_1.radius;
	SOFT_PD_A[ID].center2.x = WALKERS_2.center.x;
	SOFT_PD_A[ID].center2.y = WALKERS_2.center.y;
	SOFT_PD_A[ID].center2.z = WALKERS_2.center.z;
	SOFT_PD_A[ID].center2.r = WALKERS_2.radius;
	ORI_WALKERS[WALKERS_1.ID].statue = 1;
	ORI_WALKERS[WALKERS_2.ID].statue = 1;
	ORI_WALKERS[WALKERS_1.ID].regionID = ID;
	ORI_WALKERS[WALKERS_2.ID].regionID = ID;
	SPDA_COUNTER = SPDA_COUNTER + 1; // Optional
}

// Do not call solver
// Solver should be called only in main()
void initialUpdateSoftDomainA3D(long ID, walker WALKERS_NEW) {
	double l1 = 0.0;
	double l2 = 0.0;
	double surfaceDistance;
	double domainEdgeDistance;
	double newEdgeDistance;
	// Update Center
	domainEdgeDistance = SOFT_PD_A[ID].center_distance + 2.0*SOFT_PD_A[ID].center1.r + 2.0*SOFT_PD_A[ID].center2.r;
	l1 = sqrt((WALKERS_NEW.center.x - SOFT_PD_A[ID].center1.x)*(WALKERS_NEW.center.x - SOFT_PD_A[ID].center1.x) + (WALKERS_NEW.center.y - SOFT_PD_A[ID].center1.y)*(WALKERS_NEW.center.y - SOFT_PD_A[ID].center1.y) + (WALKERS_NEW.center.z - SOFT_PD_A[ID].center1.z)*(WALKERS_NEW.center.z - SOFT_PD_A[ID].center1.z));
	l1 = l1 - SOFT_PD_A[ID].center1.r - WALKERS_NEW.radius;
	l2 = sqrt((WALKERS_NEW.center.x - SOFT_PD_A[ID].center2.x)*(WALKERS_NEW.center.x - SOFT_PD_A[ID].center2.x) + (WALKERS_NEW.center.y - SOFT_PD_A[ID].center2.y)*(WALKERS_NEW.center.y - SOFT_PD_A[ID].center2.y) + (WALKERS_NEW.center.z - SOFT_PD_A[ID].center2.z)*(WALKERS_NEW.center.z - SOFT_PD_A[ID].center2.z));
	l2 = l2 - SOFT_PD_A[ID].center2.r - WALKERS_NEW.radius;
	if (l1 - l2 < eps) {
		surfaceDistance = l1;
		newEdgeDistance = surfaceDistance + 2.0*WALKERS_NEW.radius + 2.0*SOFT_PD_A[ID].center1.r;
		if (newEdgeDistance - domainEdgeDistance > eps) {
			// Change center 2 to New
			cout << "SPD_a change its domain center." << endl;
			SOFT_PD_A[ID].center2.x = WALKERS_NEW.center.x;
			SOFT_PD_A[ID].center2.y = WALKERS_NEW.center.y;
			SOFT_PD_A[ID].center2.z = WALKERS_NEW.center.z;
			SOFT_PD_A[ID].center2.r = WALKERS_NEW.radius;
			SOFT_PD_A[ID].center_distance = surfaceDistance;
		}
	}
	else {
		surfaceDistance = l2;
		newEdgeDistance = surfaceDistance + 2.0*WALKERS_NEW.radius + 2.0*SOFT_PD_A[ID].center2.r;
		if (newEdgeDistance - domainEdgeDistance > eps) {
			// Change center 1 to New
			cout << "SPD_a change its domain center." << endl;
			SOFT_PD_A[ID].center1.x = WALKERS_NEW.center.x;
			SOFT_PD_A[ID].center1.y = WALKERS_NEW.center.y;
			SOFT_PD_A[ID].center1.z = WALKERS_NEW.center.z;
			SOFT_PD_A[ID].center1.r = WALKERS_NEW.radius;
			SOFT_PD_A[ID].center_distance = surfaceDistance;
		}
	}
	// Update Domain
	SOFT_PD_A[ID].current_member = SOFT_PD_A[ID].current_member + 1;
	SOFT_PD_A[ID].member_list.member[SOFT_PD_A[ID].current_member] = WALKERS_NEW;
	ORI_WALKERS[WALKERS_NEW.ID].statue = 1;
	ORI_WALKERS[WALKERS_NEW.ID].regionID = ID;
}


// This function is designed for SPD_a check during initialization, if walker surface within THRES_2, initial SPD_a
// SPD_a two centerc based on the d_edge = (d_surface + 2 * r1 + 2 * r2), find the largest
// Any walker distance to one of the center less than r_center + THRES_2 + r_walker will trigger E10_update
void checkInitialSoftDomains3D() {
	double surfaceDistance;
	for (i = 0;i < inputNumber;i++) {
		if (ORI_WALKERS[i + 1].type == -1 || ORI_WALKERS[i + 1].statue == 1) { continue; }
		for (j = i + 1;j < inputNumber; j++) {
			if (ORI_WALKERS[j + 1].type == -1 || ORI_WALKERS[j + 1].statue == 1) { continue; }
			// Begin to check location
			surfaceDistance = getDistance(ORI_WALKERS[i + 1], ORI_WALKERS[j + 1]) - ORI_WALKERS[i + 1].radius - ORI_WALKERS[j + 1].radius;
			//cout << "Distance between " << i + 1 << " and " << j + 1 << " is " << surfaceDistance << endl;
			if (surfaceDistance - THRES_2 < eps) {
				cout << "Find one type 1 SPD during initialization." << endl;
				// Create New Type 1 SPD
				SPDOMAINA_INDEX = SPDOMAINA_INDEX + 1;
				initialSoftDomainA3D(SPDOMAINA_INDEX, ORI_WALKERS[i + 1], ORI_WALKERS[j + 1]);
				// Searching the rest, directly compare to the SPD 
				for (k = j + 1;k < inputNumber; k++) {
					if (ORI_WALKERS[k + 1].type == -1 || ORI_WALKERS[k + 1].statue == 1) { continue; }
					surfaceDistance = getSoftDomainASurfaceDistance(ORI_WALKERS[k + 1], SOFT_PD_A[SPDOMAINA_INDEX]);
					if (surfaceDistance - THRES_2 < eps) {
						cout << "Find update to current type 1 SPD." << endl;
						// Update Current Type 1 SPD
						initialUpdateSoftDomainA3D(SPDOMAINA_INDEX, ORI_WALKERS[k + 1]);
					}
				}
			}
		}
	}
}

// For new SPDa
// Do not call solver
// Solver should be called only in main()
void initialSoftDomainA1D(long ID, walker WALKERS_1, walker WALKERS_2) {
	double surfaceDistance;
	cout << "Initialize SPDa for " << WALKERS_1.ID << " and " << WALKERS_2.ID << " walkers." << endl;
	cout << "Center at " << WALKERS_1.center.x << " and " << WALKERS_2.center.x << endl;
	cout << "Range of " << WALKERS_1.center.x - WALKERS_1.radius - THRES_2 << " to " << WALKERS_1.center.x + WALKERS_1.radius + THRES_2 << " and " << WALKERS_2.center.x - WALKERS_2.radius - THRES_2 << " to " << WALKERS_2.center.x + WALKERS_2.radius + THRES_2 << endl;
	surfaceDistance = getDistance1D(WALKERS_1, WALKERS_2) - WALKERS_1.radius - WALKERS_2.radius;
	SOFT_PD_A[ID].ID = ID;
	SOFT_PD_A[ID].center_distance = surfaceDistance;
	SOFT_PD_A[ID].clock = 0.0;
	SOFT_PD_A[ID].events_type = 10;
	SOFT_PD_A[ID].current_member = 2;
	SOFT_PD_A[ID].member_list.member[1] = WALKERS_1;
	SOFT_PD_A[ID].member_list.member[2] = WALKERS_2;
	SOFT_PD_A[ID].center1.x = WALKERS_1.center.x;
	SOFT_PD_A[ID].center1.y = WALKERS_1.center.y;
	SOFT_PD_A[ID].center1.z = WALKERS_1.center.z;
	SOFT_PD_A[ID].center1.r = WALKERS_1.radius;
	SOFT_PD_A[ID].center2.x = WALKERS_2.center.x;
	SOFT_PD_A[ID].center2.y = WALKERS_2.center.y;
	SOFT_PD_A[ID].center2.z = WALKERS_2.center.z;
	SOFT_PD_A[ID].center2.r = WALKERS_2.radius;
	ORI_WALKERS[WALKERS_1.ID].statue = 1;
	ORI_WALKERS[WALKERS_2.ID].statue = 1;
	ORI_WALKERS[WALKERS_1.ID].regionID = ID;
	ORI_WALKERS[WALKERS_2.ID].regionID = ID;
	SPDA_COUNTER = SPDA_COUNTER + 1; // Optional
}

// For SPDa update
// Do not call solver
// Solver should be called only in main()
void initialUpdateSoftDomainA1D(long ID, walker WALKERS_NEW) {
	double l1 = 0.0;
	double l2 = 0.0;
	double surfaceDistance;
	double domainEdgeDistance;
	double newEdgeDistance;
	// Update Center
	domainEdgeDistance = SOFT_PD_A[ID].center_distance + 2.0*SOFT_PD_A[ID].center1.r + 2.0*SOFT_PD_A[ID].center2.r;
	l1 = sqrt((WALKERS_NEW.center.x - SOFT_PD_A[ID].center1.x)*(WALKERS_NEW.center.x - SOFT_PD_A[ID].center1.x));
	l1 = l1 - SOFT_PD_A[ID].center1.r - WALKERS_NEW.radius;
	l2 = sqrt((WALKERS_NEW.center.x - SOFT_PD_A[ID].center2.x)*(WALKERS_NEW.center.x - SOFT_PD_A[ID].center2.x));
	l2 = l2 - SOFT_PD_A[ID].center2.r - WALKERS_NEW.radius;
	if (l1 - l2 < eps) {
		surfaceDistance = l1;
		newEdgeDistance = surfaceDistance + 2.0*WALKERS_NEW.radius + 2.0*SOFT_PD_A[ID].center1.r;
		if (newEdgeDistance - domainEdgeDistance > eps) {
			// Change center 2 to New
			cout << "SPD_a change its domain center." << endl;
			SOFT_PD_A[ID].center2.x = WALKERS_NEW.center.x;
			SOFT_PD_A[ID].center2.r = WALKERS_NEW.radius;
			SOFT_PD_A[ID].center_distance = surfaceDistance;
		}
	}
	else {
		surfaceDistance = l2;
		newEdgeDistance = surfaceDistance + 2.0*WALKERS_NEW.radius + 2.0*SOFT_PD_A[ID].center2.r;
		if (newEdgeDistance - domainEdgeDistance > eps) {
			// Change center 1 to New
			cout << "SPD_a change its domain center." << endl;
			SOFT_PD_A[ID].center1.x = WALKERS_NEW.center.x;
			SOFT_PD_A[ID].center1.r = WALKERS_NEW.radius;
			SOFT_PD_A[ID].center_distance = surfaceDistance;
		}
	}
	// Update Domain
	SOFT_PD_A[ID].current_member = SOFT_PD_A[ID].current_member + 1;
	SOFT_PD_A[ID].member_list.member[SOFT_PD_A[ID].current_member] = WALKERS_NEW;
	ORI_WALKERS[WALKERS_NEW.ID].statue = 1;
	ORI_WALKERS[WALKERS_NEW.ID].regionID = ID;
}

// This function is designed for SPD_a check during initialization, if walker surface within THRES_2, initial SPD_a
// SPD_a two centers based on the d_edge = (d_surface + 2 * r1 + 2 * r2), find the largest
// Any walker distance to one of the center less than r_center + THRES_2 + r_walker will trigger E10_update
void checkInitialSoftDomains1D() {
	double surfaceDistance;
	for (i = 0;i < inputNumber;i++) {
		if (ORI_WALKERS[i + 1].type != -1 && ORI_WALKERS[i + 1].statue != 1) {
			for (j = i + 1;j < inputNumber; j++) {
				if (ORI_WALKERS[j + 1].type != -1 && ORI_WALKERS[j + 1].statue != 1)
				{ // Begin to check location
					surfaceDistance = getDistance1D(ORI_WALKERS[i + 1], ORI_WALKERS[j + 1]) - ORI_WALKERS[i + 1].radius - ORI_WALKERS[j + 1].radius;
					if (surfaceDistance - THRES_2 < eps) {
						cout << "Find one type 1 SPD during initialization." << endl;
						// Create New Type 1 SPD
						SPDOMAINA_INDEX = SPDOMAINA_INDEX + 1;
						initialSoftDomainA1D(SPDOMAINA_INDEX, ORI_WALKERS[i + 1], ORI_WALKERS[j + 1]);
						// Searching the rest, directly compare to the SPD 
						for (k = j + 1;k < inputNumber; k++) {
							if (ORI_WALKERS[k + 1].type != -1 && ORI_WALKERS[k + 1].statue != 1)
							{
								surfaceDistance = getSoftDomainASurfaceDistance1D(ORI_WALKERS[k + 1], SOFT_PD_A[SPDOMAINA_INDEX]);
								if (surfaceDistance - THRES_2 < eps) {
									cout << "Find update to current type 1 SPD." << endl;
									// Update Current Type 1 SPD
									initialUpdateSoftDomainA1D(SPDOMAINA_INDEX, ORI_WALKERS[k + 1]);
								}
							}
						}
					}
				}
			}
		}
	}
}

// For new HPD
// Do not call solver
// Solver should be called only in main()
void initialHardDomain3D(long ID, walker WALKER, double DOMAIN_R) {
	cout << "Initialize HPD for walker " << WALKER.ID << "." << endl;
	cout << "Center at (" << WALKER.center.x << ", " << WALKER.center.y << ", " << WALKER.center.z << "), with a radius of " << DOMAIN_R << endl;
	HARD_PD[ID].ID = ID;
	HARD_PD[ID].center.x = WALKER.center.x;
	HARD_PD[ID].center.y = WALKER.center.y;
	HARD_PD[ID].center.z = WALKER.center.z;
	HARD_PD[ID].clock = 0.0;
	HARD_PD[ID].events_type = 1;
	HARD_PD[ID].radius = DOMAIN_R;
	HARD_PD[ID].walker_ID = WALKER.ID;
	HARD_PD[ID].walker_type = WALKER.type;
	ORI_WALKERS[WALKER.ID].statue = 2;
	ORI_WALKERS[WALKER.ID].regionID = ID;
	HPD_COUNTER = HPD_COUNTER + 1;
}

// This function is designed for HPD setup during initialization
// Preassigned with sqrt ( FACTOR * DV * T) in current version
// For future version, apply space index with NNL / Hilbert curve to skip overlapping check and to enable dynamic modifier : TODO
void setInitialHardDomains3D() {
	double domainRadius;
	for (i = 0;i < inputNumber;i++) {
		if (ORI_WALKERS[i + 1].type != -1 && ORI_WALKERS[i + 1].statue == 0) {
			// TODO with k-d tree index with NNL / Hc
			domainRadius = sqrt(INIDOMAIN_MODIFIER*ORI_WALKERS[i + 1].diffusivity*INITIAL_EVENTS_TIME) + ORI_WALKERS[i + 1].radius;
			HPDOMAIN_INDEX = HPDOMAIN_INDEX + 1;
			initialHardDomain3D(HPDOMAIN_INDEX, ORI_WALKERS[i + 1], domainRadius);
		}
	}
}

// For new HPD
// Do not call solver
// Solver should be called only in main()
void initialHardDomain1D(long ID, walker WALKER, double DOMAIN_RS) {
	cout << "Initialize HPD for walker " << WALKER.ID << "." << endl;
	cout << "Range from " << WALKER.center.x - DOMAIN_RS << " to " << WALKER.center.x + DOMAIN_RS << endl;
	HARD_PD[ID].ID = ID;
	HARD_PD[ID].center.x = WALKER.center.x;
	HARD_PD[ID].center.y = 0.0;
	HARD_PD[ID].center.z = 0.0;
	HARD_PD[ID].clock = 0.0;
	HARD_PD[ID].events_type = 1;
	HARD_PD[ID].radius = DOMAIN_RS;
	HARD_PD[ID].walker_ID = WALKER.ID;
	HARD_PD[ID].walker_type = WALKER.type;
	ORI_WALKERS[WALKER.ID].statue = 2;
	ORI_WALKERS[WALKER.ID].regionID = ID;
	HPD_COUNTER = HPD_COUNTER + 1;
}

// This function is designed for HPD setup during initialization
// Preassigned with sqrt ( FACTOR * DV * T) in current version
// For future version, apply space index with NNL / Hilbert curve to skip overlapping check and to enable dynamic modifier : TODO
void setInitialHardDomains1D() {
	double domainRadius;
	for (i = 0;i < inputNumber;i++) {
		if (ORI_WALKERS[i + 1].type != -1 && ORI_WALKERS[i + 1].statue == 0) {
			// TODO with k-d tree index with NNL / Hc
			domainRadius = sqrt(INIDOMAIN_MODIFIER*ORI_WALKERS[i + 1].diffusivity*INITIAL_EVENTS_TIME) + ORI_WALKERS[i + 1].radius;
			HPDOMAIN_INDEX = HPDOMAIN_INDEX + 1;
			initialHardDomain1D(HPDOMAIN_INDEX,ORI_WALKERS[i+1],domainRadius);
		}
	}
}

// The Overlapping Checking is specific for Initialization, if overlap, set the radius of both domain to the fraction of original
// If overlap, return adjust fraction for each domain to multiply, if not, return 1.0
// This function can be applied in both 1D and 3D
// Return the domain radius fix factor 
double getHardDomainModifier(hard_domain DOMAIN_1, hard_domain DOMAIN_2) {
	double l, ld;
	double fraction = 1.0;
	l = getHardDomainCenterDistance(DOMAIN_1, DOMAIN_2);
	ld = DOMAIN_1.radius + DOMAIN_2.radius;
	if (ld - l > eps) {
		fraction = l / ld;
		cout << "HPD overlap detected." << DOMAIN_1.ID << " and " << DOMAIN_2.ID << endl;
	}
	return fraction;
}

// First pass, fully check all walkers with statue == 2
// Avoid HPD overlap and SPD fully cover
// Find the most closed HPD, check overlap
// If overlap with HPD, use getHardDomainModifier to modify both
// Find the most closed SPDa, check fully cover
// If cover with SPDa, do edge check, modifier HPD region to the location of close SPDa center walker
// 1D version
void checkInitialDomain1D() {
	double centerDistance, centerDistance1, centerDistance2;
	double minDistance = 999999999.0;
	double fix;
	int minIndex;
	for (i = 0;i < HPDOMAIN_INDEX;i++) {
		// For all the rest HPD, modify their domains, after this HPD_i will not be overlapped with any HPD 
		for (j = i + 1;j < HPDOMAIN_INDEX;j++) {
			fix = getHardDomainModifier(HARD_PD[i + 1], HARD_PD[j + 1]);
			HARD_PD[i + 1].radius = fix*HARD_PD[i + 1].radius;
			HARD_PD[j + 1].radius = fix*HARD_PD[j + 1].radius;
			if (HARD_PD[i + 1].radius - ORI_WALKERS[HARD_PD[i + 1].walker_ID].radius < eps) {
				cout << "Fatal Error 900: Cannot maintain a valid HPD. Exit." << endl;
				cout << "Walker ID " << HARD_PD[i + 1].walker_ID << ", HPD ID " << HARD_PD[i + 1].ID << endl;
				system("pause");
				exit(0);
			}
			if (HARD_PD[j + 1].radius - ORI_WALKERS[HARD_PD[j + 1].walker_ID].radius < eps) {
				cout << "Fatal Error 900: Cannot maintain a valid HPD. Exit." << endl;
				cout << "Walker ID " << HARD_PD[j + 1].walker_ID << ", HPD ID " << HARD_PD[j + 1].ID << endl;
				system("pause");
				exit(0);
			}
		}
		// Find the most close SPD center, search all SPDa
		for (j = 0;j < SPDOMAINA_INDEX;j++) {
			centerDistance1 = sqrt((HARD_PD[i + 1].center.x - SOFT_PD_A[j + 1].center1.x)*(HARD_PD[i + 1].center.x - SOFT_PD_A[j + 1].center1.x));
			centerDistance2 = sqrt((HARD_PD[i + 1].center.x - SOFT_PD_A[j + 1].center2.x)*(HARD_PD[i + 1].center.x - SOFT_PD_A[j + 1].center2.x));
			if (centerDistance1 - centerDistance2 < eps) { centerDistance = centerDistance1; }
			else { centerDistance = centerDistance2; }
			if (centerDistance - minDistance < eps) { 
				minDistance = centerDistance;
				minIndex = j + 1;
			}
		}
		// Adjust HPD
		if (minDistance - HARD_PD[i + 1].radius < eps) {
			centerDistance1 = sqrt((HARD_PD[i + 1].center.x - SOFT_PD_A[minIndex].center1.x)*(HARD_PD[i + 1].center.x - SOFT_PD_A[minIndex].center1.x));
			centerDistance2 = sqrt((HARD_PD[i + 1].center.x - SOFT_PD_A[minIndex].center2.x)*(HARD_PD[i + 1].center.x - SOFT_PD_A[minIndex].center2.x));
			if (centerDistance1 - centerDistance2 < eps) { HARD_PD[i + 1].radius = centerDistance1;}
			else { HARD_PD[i + 1].radius = centerDistance2; }
			if (HARD_PD[i + 1].radius - ORI_WALKERS[HARD_PD[i + 1].walker_ID].radius < eps) {
				cout << "Fatal Error 900: Cannot maintain a valid HPD. Exit." << endl;
				cout << "Walker ID " << HARD_PD[i + 1].walker_ID << ", HPD ID " << HARD_PD[i + 1].ID << endl;
				system("pause");
				exit(0);
			}
		}
	}
}

// First pass, fully check all walkers with statue == 2
// Avoid HPD overlap and SPD fully cover
// Find the most closed HPD, check overlap
// If overlap with HPD, use getHardDomainModifier to modify both
// Find the most closed SPDa, check fully cover
// If cover with SPDa, do edge check, modifier HPD region to the location of close SPDa center walker
// 3D version
void checkInitialDomain3D() {
	double centerDistance, centerDistance1, centerDistance2;
	double minDistance = 999999999.0;
	double fix;
	int minIndex;
	for (i = 0;i < HPDOMAIN_INDEX;i++) {
		// For all the rest HPD, modify their domains, after this HPD_i will not be overlapped with any HPD
		for (j = i + 1;j < HPDOMAIN_INDEX;j++) {
			fix = getHardDomainModifier(HARD_PD[i + 1], HARD_PD[j + 1]);
			HARD_PD[i + 1].radius = fix*HARD_PD[i + 1].radius;
			HARD_PD[j + 1].radius = fix*HARD_PD[j + 1].radius;
			if (HARD_PD[i + 1].radius - ORI_WALKERS[HARD_PD[i + 1].walker_ID].radius < eps) {
				cout << "Fatal Error 900: Cannot maintain a valid HPD. Exit." << endl;
				cout << "Walker ID " << HARD_PD[i + 1].walker_ID << ", HPD ID " << HARD_PD[i + 1].ID << endl;
				system("pause");
				exit(0);
			}
			if (HARD_PD[j + 1].radius - ORI_WALKERS[HARD_PD[j + 1].walker_ID].radius < eps) {
				cout << "Fatal Error 900: Cannot maintain a valid HPD. Exit." << endl;
				cout << "Walker ID " << HARD_PD[j + 1].walker_ID << ", HPD ID " << HARD_PD[j + 1].ID << endl;
				system("pause");
				exit(0);
			}
		}
		// Find the most close SPD center, search all SPDa
		for (j = 0;j < SPDOMAINA_INDEX;j++) {
			centerDistance1 = sqrt((HARD_PD[i + 1].center.x - SOFT_PD_A[j + 1].center1.x)*(HARD_PD[i + 1].center.x - SOFT_PD_A[j + 1].center1.x)+ (HARD_PD[i + 1].center.y - SOFT_PD_A[j + 1].center1.y)*(HARD_PD[i + 1].center.y - SOFT_PD_A[j + 1].center1.y)+ (HARD_PD[i + 1].center.z - SOFT_PD_A[j + 1].center1.z)*(HARD_PD[i + 1].center.z - SOFT_PD_A[j + 1].center1.z));
			centerDistance2 = sqrt((HARD_PD[i + 1].center.x - SOFT_PD_A[j + 1].center2.x)*(HARD_PD[i + 1].center.x - SOFT_PD_A[j + 1].center2.x)+ (HARD_PD[i + 1].center.y - SOFT_PD_A[j + 1].center2.y)*(HARD_PD[i + 1].center.y - SOFT_PD_A[j + 1].center2.y)+ (HARD_PD[i + 1].center.z - SOFT_PD_A[j + 1].center2.z)*(HARD_PD[i + 1].center.z - SOFT_PD_A[j + 1].center2.z));
			if (centerDistance1 - centerDistance2 < eps) { centerDistance = centerDistance1; }
			else { centerDistance = centerDistance2; }
			if (centerDistance - minDistance < eps) {
				minDistance = centerDistance;
				minIndex = j + 1;
			}
		}
		// Adjust HPD
		if (minDistance - HARD_PD[i + 1].radius < eps) {
			centerDistance1 = sqrt((HARD_PD[i + 1].center.x - SOFT_PD_A[minIndex].center1.x)*(HARD_PD[i + 1].center.x - SOFT_PD_A[minIndex].center1.x)+ (HARD_PD[i + 1].center.y - SOFT_PD_A[minIndex].center1.y)*(HARD_PD[i + 1].center.y - SOFT_PD_A[minIndex].center1.y)+ (HARD_PD[i + 1].center.z - SOFT_PD_A[minIndex].center1.z)*(HARD_PD[i + 1].center.z - SOFT_PD_A[minIndex].center1.z));
			centerDistance2 = sqrt((HARD_PD[i + 1].center.x - SOFT_PD_A[minIndex].center2.x)*(HARD_PD[i + 1].center.x - SOFT_PD_A[minIndex].center2.x)+ (HARD_PD[i + 1].center.y - SOFT_PD_A[minIndex].center2.y)*(HARD_PD[i + 1].center.y - SOFT_PD_A[minIndex].center2.y)+ (HARD_PD[i + 1].center.z - SOFT_PD_A[minIndex].center2.z)*(HARD_PD[i + 1].center.z - SOFT_PD_A[minIndex].center2.z));
			if (centerDistance1 - centerDistance2 < eps) { HARD_PD[i + 1].radius = centerDistance1; }
			else { HARD_PD[i + 1].radius = centerDistance2; }
			if (HARD_PD[i + 1].radius - ORI_WALKERS[HARD_PD[i + 1].walker_ID].radius < eps) {
				cout << "Fatal Error 900: Cannot maintain a valid HPD. Exit." << endl;
				cout << "Walker ID " << HARD_PD[i + 1].walker_ID << ", HPD ID " << HARD_PD[i + 1].ID << endl;
				system("pause");
				exit(0);
			}
		}
	}
}

// Second pass on all HPD and SPDa: TODO
// This function is not necessary, may functional apply to 3D only
// Last function in initialization
// Return reCheck flag, if 1, run checkInitialDomain() again
int recheckInitialDomain() {
	int checkFlag;
	checkFlag = 0;
	return checkFlag;
}

// For new HPD during events
// Do not call solver
// Solver should be called only in main()
// No long in need, replaced by ASMS
void buildHardDomain(long ID, walker WALKERS) {
	cout << "Refer to ASMS" << endl;
}

// For new SPDa during events
// Do non call solver
// Solver should be called only in main()
// No long in need, replaced by ASMS
void buildSoftDomainA(long ID, walker WALKERS_1, walker WALKERS_2) {
	cout << "Refer to ASMS" << endl;
}

// As a solver-related function, this should be a part of main(), remove this definition
// No long in need, replaced by ASMS
void updateSoftDomainA(long ID, walker NEW_WALKER) {
	cout << "Refer to ASMS" << endl;
}

int main()
{
	// Welcome to PROJECT aCRD
	cout << "//////////////////////////////////////////////////////////////////" << endl;
	cout << "The following Caluculation is Based on Version " << MAJOR_VERSION_NUMBER << "." << MINOR_VERSION_NUMBER << "." << REVISIONB_NUMBER << " [" << VERSION_ID << "]" << endl;
	cout << "PPOJECT advanced Cluster Reaction Dynamics" << endl;
	cout << "PROJECT aCRD" << endl;
	cout << "//////////////////////////////////////////////////////////////////" << endl;

	// Initilization
	cout << "System Initilization" << endl;

	// Module Loading and Integrity Check
	HMODULE dllLibDiffusivity = LoadLibrary("DiffusivityLib.dll");
	HMODULE dllSolverMajorDiff = LoadLibrary("MajorDiffSolver.dll");
	HMODULE dllSolverTransientDiff = LoadLibrary("TransientDiffSolver.dll");
	HMODULE dllSolverLocal = LoadLibrary("LocalSolver.dll");
	HMODULE dllSolverDissociation = LoadLibrary("DissociationSolver.dll");
	if (dllLibDiffusivity == NULL) {
		cerr << "Cannot find DiffusivityLib.dll. Please check framework integrity." << endl;
		releaseModule();
		return -1;
	}
	if (dllSolverMajorDiff == NULL) {
		cout << "Cannot find MajorDiffSolver.dll. Please check framework integrity." << endl;
		releaseModule();
		return -1;
	}
	if (dllSolverTransientDiff == NULL) {
		cout << "Cannot find TransientDiffSolver.dll. Please check framework integrity." << endl;
		releaseModule();
		return -1;
	}
	if (dllSolverLocal == NULL) {
		cout << "Cannot find LocalSolver.dll. Please check framework integrity." << endl;
		releaseModule();
		return -1;
	}
	if (dllSolverDissociation == NULL) {
		cout << "Cannot find DissociationSolver.dll. Please check framework integrity." << endl;
		releaseModule();
		return -1;
	}
	cout << "Finish Loading Framework." << endl;
	cout << endl;

	// Function Initilization
	typedef double(*nModuleVersion);

	// Loading Library
	cout << "Begin to Load Library Module." << endl;
	// Loading Standard Diffusivity Library
	// Rely on Diffusivity Lib (DiffusivityLib.dll)
	// Contain 1 Function: getMajorDiffTimeStamp and getDvValue
	nModuleVersion getDiffusivityLibVer = (nModuleVersion)GetProcAddress(dllLibDiffusivity, TEXT("nDvLibVersion"));
	cout << "Framework Load Version " << *getDiffusivityLibVer << " Diffusivity Library." << endl;
	typedef double(*functionLibDv)(double, int, long);
	// Function double getDvValue(double E, int type, long size)
	// Diffusivity: nm^2/s
	// Energy unit in keV, Type [int] (0 for V, 1 for I), Size [Long] (Cluster)
	functionLibDv getDvValue = (functionLibDv)GetProcAddress(dllLibDiffusivity, TEXT("getDvValue"));
	if (getDvValue == NULL) {
		cout << "Cannot link to the function getDvValue. Please check the version of DiffusivityLib.dll." << endl;
		releaseModule();
		return -1;
	}
	cout << "Loading Diffusivity Library... Complete." << endl;
	// Full Lib for all parameter (radius, crystal, etc): TODO
	cout << "Finish Loading Library Module." << endl;
	cout << endl;

	// Loading Solver
	cout << "Begin to Load Solver Module." << endl;
	// Loading Event 01 Solver from dllSolverMajorDiff
	// Rely on Major Diffusion Solver (MajorDiffSolver.dll)
	// Contain 4 Functions: getMajorDiffTimeStamp1D, getMajorDiffTimeStamp3D, getMajorDiffRelativeLocation1D and getMajorDiffRelativeLocation3D
	nModuleVersion getE01SolverVer = (nModuleVersion)GetProcAddress(dllSolverMajorDiff, TEXT("nE01SolverVersion"));
	cout << "Framework Load Version " << *getE01SolverVer << " Major Diffusion Solver." << endl;
	typedef double(*functionE01_1)(double, double, double);
	typedef struct location(*functionE01_2)(double, double, double, double, location);
	// Function double getMajorDiffTimeStamp1D(double Dv, double domainRadius, double walkerRadius)
	// TimeStamp: s
	// Dv unit in nm^2/s, Radius unit in nm
	functionE01_1 getMajorDiffTimeStamp1D = (functionE01_1)GetProcAddress(dllSolverMajorDiff, TEXT("getMajorDiffTimeStamp1D"));
	// Function double getMajorDiffTimeStamp3D(double Dv, double domainRadius, double walkerRadius)
	// TimeStamp: s
	// Dv unit in nm^2/s, Radius unit in nm
	functionE01_1 getMajorDiffTimeStamp3D = (functionE01_1)GetProcAddress(dllSolverMajorDiff, TEXT("getMajorDiffTimeStamp3D"));
	// Function location getMajorDiffRelativeLocation1D(double Dv, double timestamp, double domainRadius, double walkerRadius, location centerLocation)
	// Location unit in nm
	// Dv unit in nm^2/s, time unit in s, Radius unit in nm, location unit in nm
	functionE01_2 getMajorDiffRelativeLocation1D = (functionE01_2)GetProcAddress(dllSolverMajorDiff, TEXT("getMajorDiffRelativeLocation1D"));
	// Function location getMajorDiffRelativeLocation3D(double Dv, double timestamp, double domainRadius, double walkerRadius, location centerLocation)
	// Location unit in nm
	// Dv unit in nm^2/s, times unit is s, Radius unit in nm, location unit in nm
	functionE01_2 getMajorDiffRelativeLocation3D = (functionE01_2)GetProcAddress(dllSolverMajorDiff, TEXT("getMajorDiffRelativeLocation3D"));
	if (getMajorDiffTimeStamp1D == NULL) {
		cout << "Cannot link to the function getMajorDiffTimeStamp1D. Please check the version of MajorDiffSolver.dll." << endl;
		releaseModule();
		return -1;
	}
	if (getMajorDiffTimeStamp3D == NULL) {
		cout << "Cannot link to the function getMajorDiffTimeStamp3D. Please check the version of MajorDiffSolver.dll." << endl;
		releaseModule();
		return -1;
	}
	if (getMajorDiffRelativeLocation1D == NULL) {
		cout << "Cannot link to the function getMajorDiffRelativeLocation1D. Please check the version of MajorDiffSolver.dll." << endl;
		releaseModule();
		return -1;
	}
	if (getMajorDiffRelativeLocation3D == NULL) {
		cout << "Cannot link to the function getMajorDiffRelativeLocation3D. Please check the version of MajorDiffSolver.dll." << endl;
		releaseModule();
		return -1;
	}
	cout << "Loading Major Diffusion (E01) Solver Module... Complete." << endl;

	// Loading Event 02 Solver from dllSolverTransientDiff
	// Rely on Transient Diffusion Solver (TransientDiffSolver.dll)
	// Contain 2 Functions: getTransDiffRelativeLocation1D and getTransDiffRelativeLocation3D
	nModuleVersion getE02SolverVer = (nModuleVersion)GetProcAddress(dllSolverTransientDiff, TEXT("nE02SolverVersion"));
	cout << "Framework Load Version " << *getE02SolverVer << " Transient Diffusion Solver." << endl;
	typedef struct location(*functionE02_1)(double, double, double, double, location);
	// Function location getTransDiffRelativeLocation1D(double Dv, double timestamp, double domainRadius, double walkerRadius, location centerLocation)
	// Location unit in nm
	// Time unit in s
	// Dv unit in nm^2/s, Radius unit in nm, location unit in nm
	functionE02_1 getTransDiffRelativeLocation1D = (functionE02_1)GetProcAddress(dllSolverTransientDiff, TEXT("getTransDiffRelativeLocation1D"));
	// Function location getTransDiffRelativeLocation3D(double Dv, double timestamp, double domainRadius, double walkerRadius, location centerLocation)
	// Location unit in nm
	// Time unit in s
	// Dv unit in nm^2/s, Radius unit in nm, location unit in nm
	functionE02_1 getTransDiffRelativeLocation3D = (functionE02_1)GetProcAddress(dllSolverTransientDiff, TEXT("getTransDiffRelativeLocation3D"));
	if (getTransDiffRelativeLocation1D == NULL) {
		cout << "Cannot link to the function getTransDiffRelativeLocation1D. Please check the version of TransientDiffSolver.dll." << endl;
		releaseModule();
		return -1;
	}
	if (getTransDiffRelativeLocation3D == NULL) {
		cout << "Cannot link to the function getTransDiffRelativeLocation3D. Please check the version of TransientDiffSolver.dll." << endl;
		releaseModule();
		return -1;
	}
	cout << "Loading Transient Diffusion (E02) Solver Module... Complete." << endl;

	// Loading Event 10 Solver from dllSolverLocal
	// Rely on Local Reaction Solver (LocalSolver.dll)
	// Annihilation only version in 0.5 while coalescence version in 1.0
	// Contain 4 Functions: getReactionTimeStamp1D, getReactionTimeStamp3D, getReactionFinalStatues1D and getReactionFinalStatues3D
	nModuleVersion getE10SolverVer = (nModuleVersion)GetProcAddress(dllSolverLocal, TEXT("nE10SolverVersion"));
	cout << "Framework Load Version " << *getE10SolverVer << " Local Reaction Solver." << endl;
	typedef double(*functionE10_1)(reaction_current_walker, double);
	typedef struct reaction_final_walker(*functionE10_2)(reaction_current_walker, int, domaincenter, domaincenter);
	// Function double getReactionTimeStamp1D(reaction_current_walker walkers, double initial_time)
	// Input with current walker list, storaged in SPDa
	// Return with SPDa break time: s
	functionE10_1 getReactionTimeStamp1D = (functionE10_1)GetProcAddress(dllSolverLocal, TEXT("getReactionTimeStamp1D"));
	// Function double getReactionTimeStamp3D(reaction_current_walker walkers, double initial_time)
	// Input with current walker list and overall number, storaged in SPDa
	// Return with SPDa break time: s
	functionE10_1 getReactionTimeStamp3D = (functionE10_1)GetProcAddress(dllSolverLocal, TEXT("getReactionTimeStamp3D"));
	// Function reaction_final_walker getReactionFinalStatues1D(reaction_current_walker walkers, int init_member, domaincenter centerA, domaincenter centerB)
	// Input with current walker list, overall number and center info, storaged in SPDa
	// Return with final walker list, storage to SPDa
	functionE10_2 getReactionFinalStatues1D = (functionE10_2)GetProcAddress(dllSolverLocal, TEXT("getReactionFinalStatues1D"));
	// Function reaction_final_walker getReactionFinalStatues3D(reaction_current_walker walkers, int init_member, domaincenter centerA, domaincenter centerB)
	// Input with current walker list, overall number and center info, storaged in SPDa
	// Return with final walker list, storage to SPDa
	functionE10_2 getReactionFinalStatues3D = (functionE10_2)GetProcAddress(dllSolverLocal, TEXT("getReactionFinalStatues3D"));
	if (getReactionTimeStamp1D == NULL) {
		cout << "Cannot link to the function getReactionTimeStamp1D. Please check the version of LocalSolver.dll." << endl;
		releaseModule();
		return -1;
	}
	if (getReactionTimeStamp3D == NULL) {
		cout << "Cannot link to the function getReactionTimeStamp3D. Please check the version of LocalSolver.dll." << endl;
		releaseModule();
		return -1;
	}
	if (getReactionFinalStatues1D == NULL) {
		cout << "Cannot link to the function getReactionFinalStatues1D. Please check the version of LocalSolver.dll." << endl;
		releaseModule();
		return -1;
	}
	if (getReactionFinalStatues3D == NULL) {
		cout << "Cannot link to the function getReactionFinalStatues3D. Please check the version of LocalSolver.dll." << endl;
		releaseModule();
		return -1;
	}
	cout << "Loading Local Reaction (E10) Solver Module... Complete." << endl;

	// Loading Event 11 Solver from dllSolverDissociation
	// Rely on Dissociation Solver (DissociationSolver.dll)
	// Contain 1 function: getDissociationFlag
	// Future contain 2 more functions: getDissociationFinalStatues1D and getDissociationFinalStatues3D
	// May contain other 2 functions: getDissociationTimestamp1D and getDissociationTimestamp3D
	nModuleVersion getE11SolverVer = (nModuleVersion)GetProcAddress(dllSolverDissociation, TEXT("nE11SolverVersion"));
	cout << "Framework Load Version " << *getE11SolverVer << " Dissociation Solver." << endl;
	typedef int(*functionE11_1)(walker);
	// Function int getDissociationFlag(walker WALKER)
	// Input with walker
	// Return with dissociation flag, 1: E11, 0: Nothing
	functionE11_1 getDissociationFlag = (functionE11_1)GetProcAddress(dllSolverDissociation, TEXT("getDissociationFlag"));
	if (getDissociationFlag == NULL) {
		cout << "Cannot link to the function getDissociationFlag. Please check the version of DissociationSolver.dll." << endl;
		releaseModule();
		return -1;
	}
	cout << "Loading Dissociation (E11) Solver Module... Complete." << endl;

	// E20+ & Other Events
	// Currently no solver in need, all embedded with ASMS logic

	cout << "Finish Loading Solver Module." << endl;
	cout << "//////////////////////////////////////////////////////////////////" << endl;

	// Loading User Data
Test_Start_1:
	cout << "Begin to load User Data." << endl;
	// System Parameter
	filenameInput = "System.ini";
	fstream fileInput1(filenameInput, ios::in);
	if (!fileInput1)
	{
		cerr << "Cannot Open System.ini" << endl;
		system("pause");
		return -1;
	}
	fileInput1 >> DISPLAYFLAG_SYS >> ININOTES;
	fileInput1 >> DISPLAYFLAG_DEB >> ININOTES;
	fileInput1 >> END_CLOCK >> ININOTES;
	fileInput1 >> DUMPFLAG_INN >> ININOTES;
	fileInput1 >> DUMPLAMMPSFLAG_INN >> ININOTES;
	fileInput1 >> DUMPCOUNTER_FLAG >> ININOTES;
	fileInput1.close();
	// Documents Input Initial
	// First Line to be the totol number of input
	// Second Line to be the dimesion number
	// Third line to be the boundary (Third to Fifth line for 3D)
	// Each line (3D): Type, Cluster Number, X, Y, Z
	// ofstream fileOutput("Output.txt")
	cout << "Please Input the Filename. Located in the Working Direction. Location Unit in nm." << endl;
	cin >> filenameInput;
	fstream fileInput(filenameInput, ios::in);
	if (!fileInput)
	{
		cerr << "Cannot Open the File." << endl;
		system("pause");
		return -1;
	}
	fileInput >> inputNumber;
	fileInput >> space_dimension;
	switch (space_dimension)
	{
	case 1:
		fileInput >> X_LOWWER_LIMIT >> X_UPPER_LIMIT;
		if (X_LOWWER_LIMIT - X_UPPER_LIMIT > eps) {
			cerr << "Error in boundary setting. Please check your input file." << endl;
			system("pause");
			return(-1);
		}
		for (i = 0;i < inputNumber; ++i) {
			fileInput >> ORI_WALKERS[i + 1].type >> ORI_WALKERS[i + 1].component >> ORI_WALKERS[i + 1].center.x;
			if (ORI_WALKERS[i + 1].center.x - X_LOWWER_LIMIT<eps || ORI_WALKERS[i + 1].center.x - X_UPPER_LIMIT>eps) {
				cerr << "Error in data. Please check your input file." << endl;
				system("pause");
				return(-1);
			}
			ORI_WALKERS[i + 1].radius = getEstimateRadius(ORI_WALKERS[i + 1]);
			ORI_WALKERS[i + 1].ID = i + 1;
			ORI_WALKERS[i + 1].statue = 0;
			ORI_WALKERS[i + 1].energy = 0.0; // Not in use for current version
			ORI_WALKERS[i + 1].diffusivity = getDvValue(ORI_WALKERS[i + 1].energy, ORI_WALKERS[i + 1].type, ORI_WALKERS[i + 1].component);
			ORI_WALKERS[i + 1].timestamp = 0.0;
		}
		break;
	case 3:
		fileInput >> X_LOWWER_LIMIT >> X_UPPER_LIMIT;
		if (X_LOWWER_LIMIT - X_UPPER_LIMIT > eps) {
			cerr << "Error in boundary setting. Please check your input file." << endl;
			system("pause");
			return(-1);
		}
		fileInput >> Y_LOWWER_LIMIT >> Y_UPPER_LIMIT;
		if (Y_LOWWER_LIMIT - Y_UPPER_LIMIT > eps) {
			cerr << "Error in boundary setting. Please check your input file." << endl;
			system("pause");
			return(-1);
		}
		fileInput >> Z_LOWWER_LIMIT >> Z_UPPER_LIMIT;
		if (Z_LOWWER_LIMIT - Z_UPPER_LIMIT > eps) {
			cerr << "Error in boundary setting. Please check your input file." << endl;
			system("pause");
			return(-1);
		}
		for (i = 0;i < inputNumber; ++i) {
			fileInput >> ORI_WALKERS[i + 1].type >> ORI_WALKERS[i + 1].component >> ORI_WALKERS[i + 1].center.x >> ORI_WALKERS[i + 1].center.y >> ORI_WALKERS[i + 1].center.z;
			if (ORI_WALKERS[i + 1].center.x - X_LOWWER_LIMIT<eps || ORI_WALKERS[i + 1].center.x - X_UPPER_LIMIT>eps) {
				cerr << "Error in data. Please check your input file." << endl;
				system("pause");
				return(-1);
			}
			if (ORI_WALKERS[i + 1].center.y - Y_LOWWER_LIMIT<eps || ORI_WALKERS[i + 1].center.y - Y_UPPER_LIMIT>eps) {
				cerr << "Error in data. Please check your input file." << endl;
				system("pause");
				return(-1);
			}
			if (ORI_WALKERS[i + 1].center.z - Z_LOWWER_LIMIT<eps || ORI_WALKERS[i + 1].center.z - Z_UPPER_LIMIT>eps) {
				cerr << "Error in data. Please check your input file." << endl;
				system("pause");
				return(-1);
			}
			ORI_WALKERS[i + 1].radius = getEstimateRadius(ORI_WALKERS[i + 1]);
			ORI_WALKERS[i + 1].ID = i + 1;
			ORI_WALKERS[i + 1].statue = 0;
			ORI_WALKERS[i + 1].energy = 0.0; // Not in use for current version
			ORI_WALKERS[i + 1].diffusivity = getDvValue(ORI_WALKERS[i + 1].energy, ORI_WALKERS[i + 1].type, ORI_WALKERS[i + 1].component);
			ORI_WALKERS[i + 1].timestamp = 0.0;
		}
		break;
	default:
		cerr << "Error in dimension setting. Please check your input file." << endl;
		system("pause");
		return(-1);
		break;
	}
	fileInput.close();
	WALKER_COUNTER = inputNumber;
	HPDOMAIN_INDEX = 0;
	SPDOMAINA_INDEX = 0;
	SPDOMAINB_INDEX = 0;
	cout << "Total " << WALKER_COUNTER << " walkers are read from current initial data." << endl;
	// Finish User Data. ORI_WALKERS[], WALKER_INDEX, X/Y/Z_LIMITS, space_dimension are useful data.
	double dump_min;
	double dump_max;
	dump_min = X_LOWWER_LIMIT;
	dump_max = X_UPPER_LIMIT;
	cout << "Finish Loading User Data." << endl;
	cout << "//////////////////////////////////////////////////////////////////" << endl;

	WALKER_INDEX = inputNumber; // Reserve for New Walkers, Storage to ORI_WALKERS, inputNumber is not supposed to be changed.
	// Test reading
	//for (i = 0;i < inputNumber; ++i) {
	//	cout << ORI_WALKERS[i + 1].ID << " " << ORI_WALKERS[i + 1].type << " " << ORI_WALKERS[i + 1].component << " " << ORI_WALKERS[i + 1].center.x << " " << ORI_WALKERS[i + 1].center.y << " " << ORI_WALKERS[i + 1].center.z << " with R = " << ORI_WALKERS[i + 1].radius << endl;
	//	cout << "Dv " << ORI_WALKERS[i + 1].diffusivity << endl;
	//}
	// WALKER_INDEX++; // WALKER_INDEX is now point to the new blank walker in array, this is for debug.

	// Initialization
	// Now all walkers are released
	// Build Initial Domains: delete walker overlap, allow HPD (E01+E11) and SPD_a (E10). Able to fit with restart or sync.
	checkWalkerOverlapping();

	// Build SPD_a, then build HPD for the rest, run overlap check for HPD
	// In current Version, SPD_a will be set to fixed break time and annihilation only
	// The DOMAIN_INDEX is also updated, point to the new blank domain in array
	if (space_dimension == 1) {
		cout << "Begin to initial 1D space." << endl;
		checkInitialSoftDomains1D();
		setInitialHardDomains1D();
		int Flag_Initial_HPD = 1;
		while (Flag_Initial_HPD) {
			if (DISPLAYFLAG_DEB == 1) {
				cout << "Begin a new HPD check round." << endl;
			}
			checkInitialDomain1D();
			Flag_Initial_HPD = recheckInitialDomain();
		}
	}

	if (space_dimension == 3) {
		cout << "Begin to initial 3D space." << endl;
		checkInitialSoftDomains3D();
		setInitialHardDomains3D();
		int Flag_Initial_HPD = 1;
		while (Flag_Initial_HPD) {
			if (DISPLAYFLAG_DEB == 1) {
				cout << "Begin a new HPD check round." << endl;
			}
			checkInitialDomain3D();
			// Recheck NOT in use in current version
			Flag_Initial_HPD = recheckInitialDomain();
		}
	}

	cout << endl;
	cout << "Current PROJECT aCRD main stage contains " << WALKER_COUNTER << " walkers in " << HPD_COUNTER << " HPD, " << SPDA_COUNTER << " Type 1 SPD and " << SPDB_COUNTER << " Type 2 SPD." << endl;
	cout << "Finish Space Initialization." << endl;
	cout << "//////////////////////////////////////////////////////////////////" << endl;
	// Introduce Solvers and Build Event Queue
Test_Start_2:
	cout << "Begin to apply solvers and set event queue." << endl;
	WORLD_CLOCK = 0.0;
	// HPD timestamp initialization and event set
	// E01 or E11 (Not in this version)
	if (space_dimension == 1) {
		for (i = 0;i < HPDOMAIN_INDEX;i++) {
			EVENT_INDEX = EVENT_INDEX + 1;
			HARD_PD[i + 1].events_time = getMajorDiffTimeStamp1D(ORI_WALKERS[HARD_PD[i + 1].walker_ID].diffusivity, HARD_PD[i + 1].radius, ORI_WALKERS[HARD_PD[i + 1].walker_ID].radius);
			EVENT_LIST[EVENT_INDEX].domainID = HARD_PD[i + 1].ID;
			EVENT_LIST[EVENT_INDEX].domainType = 0;
			EVENT_LIST[EVENT_INDEX].timestamp = HARD_PD[i + 1].events_time + WORLD_CLOCK;
			EVENT_LIST[EVENT_INDEX].eventType = 1;
			HARD_PD[i + 1].event_card = EVENT_INDEX;
		}
	}
	if (space_dimension == 3) {
		for (i = 0;i < HPDOMAIN_INDEX;i++) {
			EVENT_INDEX = EVENT_INDEX + 1;
			HARD_PD[i + 1].events_time = getMajorDiffTimeStamp3D(ORI_WALKERS[HARD_PD[i + 1].walker_ID].diffusivity, HARD_PD[i + 1].radius, ORI_WALKERS[HARD_PD[i + 1].walker_ID].radius);
			EVENT_LIST[EVENT_INDEX].domainID = HARD_PD[i + 1].ID;
			EVENT_LIST[EVENT_INDEX].domainType = 0;
			EVENT_LIST[EVENT_INDEX].timestamp = HARD_PD[i + 1].events_time + WORLD_CLOCK;
			EVENT_LIST[EVENT_INDEX].eventType = 1;
			HARD_PD[i + 1].event_card = EVENT_INDEX;
		}
	}
	
	// SPDa timestamp initialization and event set
	// E10
	if (space_dimension == 1) {
		for (i = 0;i < SPDOMAINA_INDEX;i++) {
			EVENT_INDEX = EVENT_INDEX + 1;
			// SPDa will trigger within INITIAL_EVENTS_TIME + INITIAL_EVENTS_TIME * Random[0,1)
			// RANDOM_SEED = rand() / (double)(RAND_MAX);
			// SOFT_PD_A[i + 1].events_time = INITIAL_EVENTS_TIME + INITIAL_EVENTS_TIME*RANDOM_SEED; // Replaced with solver interface once LocalSolver is provided
			SOFT_PD_A[i + 1].events_time = getReactionTimeStamp1D(SOFT_PD_A[i + 1].member_list, INITIAL_EVENTS_TIME); // With solver interface
			EVENT_LIST[EVENT_INDEX].domainID = SOFT_PD_A[i + 1].ID;
			EVENT_LIST[EVENT_INDEX].domainType = 1;
			EVENT_LIST[EVENT_INDEX].timestamp = SOFT_PD_A[i + 1].events_time + WORLD_CLOCK;
			EVENT_LIST[EVENT_INDEX].eventType = 10;
			SOFT_PD_A[i + 1].event_card = EVENT_INDEX;
		}
	}
	if (space_dimension == 3) {
		for (i = 0;i < SPDOMAINA_INDEX;i++) {
			EVENT_INDEX = EVENT_INDEX + 1;
			// This version will not call local solver, SPDa will trigger within INITIAL_EVENTS_TIME + INITIAL_EVENTS_TIME * Random[0,1)
			// RANDOM_SEED = rand() / (double)(RAND_MAX);
			// SOFT_PD_A[i + 1].events_time = INITIAL_EVENTS_TIME + INITIAL_EVENTS_TIME*RANDOM_SEED; // Replaced with solver interface once LocalSolver is provided
			SOFT_PD_A[i + 1].events_time = getReactionTimeStamp3D(SOFT_PD_A[i + 1].member_list, INITIAL_EVENTS_TIME); // With solver interface
			EVENT_LIST[EVENT_INDEX].domainID = SOFT_PD_A[i + 1].ID;
			EVENT_LIST[EVENT_INDEX].domainType = 1;
			EVENT_LIST[EVENT_INDEX].timestamp = SOFT_PD_A[i + 1].events_time + WORLD_CLOCK;
			EVENT_LIST[EVENT_INDEX].eventType = 10;
			SOFT_PD_A[i + 1].event_card = EVENT_INDEX;
		}
	}
	// Debug solver from solver module

	// SPDb is hybrid-driven, use world clock control or fix step events
	// Not in this version

	// E20 New Walker Event
	// Only Event Information would be applied with type=20, No need for domain
	// Use extra file to control the input infomation
	// New walker will be introduced in event queue followed by ASMS
	// TODO
	cout << "Begin to Load Extra Walker (irradiation) information." << endl;
	filenameInput = "Irradiation.ini";
	fstream fileInput2(filenameInput, ios::in);
	if (!fileInput2)
	{
		cerr << "Cannot Open Irradiation.ini" << endl;
		system("pause");
		return -1;
	}
	fileInput2 >> irradiationCheck >> ININOTES;
	if (irradiationCheck == 0) {
		cout << "No Extra Walker (irradiation) introduced in current calculation." << endl;
	}
	else {
		cout << "Loading Extra Walker (irradiation) information." << endl;
		fileInput2 >> irradiationPulse >> ININOTES;
		fileInput2 >> irradiationINumber >> ININOTES;
		fileInput2 >> irradiationVNumber >> ININOTES;
		cout << "Every " << irradiationPulse << "s, introducing " << irradiationINumber << " single I and " << irradiationVNumber << " single V." << endl;
		//system("pause");
		// Initial first event card
		// No need to decide where to insert, just get a card with type 20
		EVENT_INDEX = EVENT_INDEX + 1;
		EVENT_LIST[EVENT_INDEX].domainID = 0;
		EVENT_LIST[EVENT_INDEX].domainType = 3;
		EVENT_LIST[EVENT_INDEX].timestamp = WORLD_CLOCK + irradiationPulse;
		EVENT_LIST[EVENT_INDEX].eventType = 20;
	}
	fileInput2.close();

	// Sort event list
	// Use malloc with pointer instead of assigned array: Not in this version
	// Use type queue instead of type array: Not in this version
	for (i = 1; i < EVENT_INDEX;i++) {
		for (j = 1;j < EVENT_INDEX + 1 - i;j++) {
			if (EVENT_LIST[j].timestamp>EVENT_LIST[j + 1].timestamp) {
				swap(EVENT_LIST[j], EVENT_LIST[j + 1]);
			}
		}
	}
	// Update domain's card info
	for (i = 1;i < EVENT_INDEX;i++) {
		if (EVENT_LIST[i].domainType == 0) {
			HARD_PD[EVENT_LIST[i].domainID].event_card = i;
		}
		if (EVENT_LIST[i].domainType == 1) {
			SOFT_PD_A[EVENT_LIST[i].domainID].event_card = i;
		}
	}
	// Debug output
	if (DISPLAYFLAG_DEB == 1) {
		for (i = 0;i < EVENT_INDEX;i++) {
			if (EVENT_LIST[i + 1].domainType == 0) {
				cout << "ID " << EVENT_LIST[i + 1].domainID << " HPD will have an event at " << EVENT_LIST[i + 1].timestamp << endl;
			}
			if (EVENT_LIST[i + 1].domainType == 1) {
				cout << "ID " << EVENT_LIST[i + 1].domainID << " Type 1 SPD will have an event at " << EVENT_LIST[i + 1].timestamp << endl;
			}
			if (EVENT_LIST[i + 1].domainType == 3) {
				cout << "New walker joins at " << EVENT_LIST[i + 1].timestamp << endl;
			}
		}
	}
	cout << endl;
	cout << "Finish Initialization." << endl;
	cout << "//////////////////////////////////////////////////////////////////" << endl;
	//system("pause");
	// Event Loop
Test_Start_3:
	cout << "Entering main event stage." << endl;
	if (*getE10SolverVer - 0.9 < eps) {
		cout << "Current Version aCRD in (a) + (a) - > (0) dynamics." << endl;
	}
	else {
		cout << "Current Version aCRD in (a) + (b) - > (a + b) dynamics." << endl;
	}
	
	// 1D prototype concerns a + a - > 0 annihilation
	// S2	-	E01	-	S0	-	Rebuild	-	S1/S2
	// Rebuild - E02 check
	// S0	-	E02	-	S0	-	Rebuild	-	S1/S2
	// S1	-	E10	-	Delete

	// 3D prototype with (a) + (b) - > (a + b) coalescence
	// S2	-	E01	-	S0	-	ASMS1(E02/E10 Update)
	// S2	-	E02	-	S0	-	ASMS1(E02/E10 Update)	
	// S1	-	E10	-	S0	-	ASMS1(E02/E10 Update)	
	// S0	-	ASMS2	-	S1/S2

	// Require full solver
	// Make sure identify all walkers which change their statue. No S0 is allowed after event creator.
	STAGE_FLAG = 1; // Manual
	CONTROL_FLAG = 0;
	EVENT_POINTER = 1;
	// 1 Dimension begins here
	// No further upgrade plan
	if (space_dimension == 1) {
		while (STAGE_FLAG)
		{
			cout << "Running 1 dimension events." << endl;
			if (WORLD_CLOCK - END_CLOCK < eps) {
				CONTROL_FLAG = 1; // Check with clock
			}
			if (EVENT_POINTER > EVENT_INDEX) {
				CONTROL_FLAG = 0; // All events finished
			}
			while (CONTROL_FLAG)
			{
				// Begin to execute event
				// Read event card
				// Extra judge with E20 Event: TODO
				if (DISPLAYFLAG_SYS == 1) {
					cout << "Begin a new event." << endl;
				}
				currentEvent::domain_type = EVENT_LIST[EVENT_POINTER].domainType;
				currentEvent::domainID = EVENT_LIST[EVENT_POINTER].domainID;
				currentEvent::timestamp = EVENT_LIST[EVENT_POINTER].timestamp;
				WORLD_CLOCK = currentEvent::timestamp;
				// Execute E20 first, if not E20, do the rest

				// Locate domain
				if (currentEvent::domain_type == 0) {
					// HPD event
					currentEvent::event_type = HARD_PD[currentEvent::domainID].events_type;
					// HPD concern E01 and E11, only E01 for this version
					// May be cancelled by other event
					if (currentEvent::event_type == -1) {
						EVENT_POINTER = EVENT_POINTER + 1;
					}
					if (currentEvent::event_type == 1) {
						// Break HPD
						// Disable the domain
						HARD_PD[currentEvent::domainID].events_type = -1;
						// Call solver to get location, set walker to S0, send walker
						// getMajorDiffRelativeLocation1D(double Dv, double timestamp, double domainRadius, double walkerRadius, location centerLocation);
						currentEvent::afterward_location = getMajorDiffRelativeLocation1D(ORI_WALKERS[HARD_PD[currentEvent::domainID].walker_ID].diffusivity, HARD_PD[currentEvent::domainID].events_time, HARD_PD[currentEvent::domainID].radius, ORI_WALKERS[HARD_PD[currentEvent::domainID].walker_ID].radius, HARD_PD[currentEvent::domainID].center);
						ORI_WALKERS[HARD_PD[currentEvent::domainID].walker_ID].statue = 0; // S0 release
						ORI_WALKERS[HARD_PD[currentEvent::domainID].walker_ID].center = currentEvent::afterward_location; // Move walker
						// Move EVENT_POINTER to next event
						EVENT_POINTER = EVENT_POINTER + 1;
						// Enter ASMS 1st Stage
						S0_INDEX = 1;
						S0_EVENT[1] = HARD_PD[currentEvent::domainID].walker_ID; // Log the IDs of all walkers with S0
						S0_FULL_FLAG = 1;
						while (S0_FULL_FLAG)
						{
							// Search all exist S0 walkers
							for (i = 0;i < S0_INDEX;i++) {
								S0_FULL_FLAG = 0; // Default no need for next round, if any E02 is triggered, then set back to 1
												  // walker: ORI_WALKERS[S0_EVENT[i+1]]
								for (j = 0;j < HPDOMAIN_INDEX;j++) {
									if (HARD_PD[j + 1].events_type != -1) {
										// Get distance from current S0 walker to each HPD surface, check whether need E02 or not (TH1)
										if (getHardDomainSurfaceDistance1D(ORI_WALKERS[S0_EVENT[i + 1]], HARD_PD[j + 1]) - THRES_1 < eps) {
											// E02 is trigged for this HPD
											ORI_WALKERS[HARD_PD[j + 1].walker_ID].center = getTransDiffRelativeLocation1D(ORI_WALKERS[HARD_PD[j + 1].walker_ID].diffusivity, WORLD_CLOCK - ORI_WALKERS[HARD_PD[j + 1].walker_ID].timestamp, HARD_PD[j + 1].radius, ORI_WALKERS[HARD_PD[j + 1].walker_ID].radius, HARD_PD[j + 1].center);
											// Break, walker to S0, index + 1
											HARD_PD[j + 1].events_type = -1;
											ORI_WALKERS[HARD_PD[j + 1].walker_ID].statue = 0;
											S0_INDEX = S0_INDEX + 1;
											S0_EVENT[S0_INDEX] = HARD_PD[j + 1].walker_ID;
											S0_FULL_FLAG = 1;
										}
									}
									// Then check next HPD
								}
								// Finish check for the triggered S0, no HPD is still close, to next round: check new full S0 list
							}
						}
						cout << "Overall " << S0_INDEX << " Walkers Released in ASMS 1st Stage" << endl;
						for (i = 0;i < S0_INDEX;i++) {
							if (ORI_WALKERS[S0_EVENT[i + 1]].statue == 0 && ORI_WALKERS[S0_EVENT[i + 1]].type != -1) {
								if (ORI_WALKERS[S0_EVENT[i + 1]].center.x + ORI_WALKERS[S0_EVENT[i + 1]].radius - X_UPPER_LIMIT>eps || ORI_WALKERS[S0_EVENT[i + 1]].center.x - ORI_WALKERS[S0_EVENT[i + 1]].radius - X_LOWWER_LIMIT < eps) {
									ORI_WALKERS[S0_EVENT[i + 1]].type = -1; // Inable this particle
									cout << "1 walker exits the boundary" << endl;
								}
							}
						}
						// Enter ASMS 2nd Stage
						// Fast pass: N = 1
						if (S0_INDEX == 1) {
							// Avoid the walkers out of boundary
							if (ORI_WALKERS[S0_EVENT[1]].type != -1) {
								// Whether this S0 walker located in SPDa
								currentEvent::flag_s0_in_SPDa = 0;
								for (i = 0;i < SPDOMAINA_INDEX;i++) {
									if (SOFT_PD_A[i + 1].events_type != -1) {
										currentEvent::walker_distance = getSoftDomainASurfaceDistance1D(ORI_WALKERS[S0_EVENT[1]], SOFT_PD_A[i + 1]);
										if (currentEvent::walker_distance - THRES_2 < eps) {
											currentEvent::flag_s0_in_SPDa = 1;
											// Solver, new event
											// Update SPDa
											// Location adjust is not in this version: TODO
											SOFT_PD_A[i + 1].current_member = SOFT_PD_A[i + 1].current_member + 1;
											SOFT_PD_A[i + 1].member_list.member[SOFT_PD_A[i + 1].current_member] = ORI_WALKERS[S0_EVENT[1]];
											currentEvent::update_spda_time_temp = SOFT_PD_A[i + 1].events_time;
											SOFT_PD_A[i + 1].events_time = SOFT_PD_A[i + 1].events_time; // Replaced with solver update interface once LocalSolver is provided: TODO
											currentEvent::update_spda_time_temp = SOFT_PD_A[i + 1].events_time - currentEvent::update_spda_time_temp;
											ORI_WALKERS[S0_EVENT[1]].statue = 1;
											ORI_WALKERS[S0_EVENT[1]].timestamp = WORLD_CLOCK;
											// No new event, just timestamp update
											EVENT_LIST[SOFT_PD_A[i + 1].event_card].timestamp = EVENT_LIST[SOFT_PD_A[i + 1].event_card].timestamp + currentEvent::update_spda_time_temp;
											cout << "1/1 Walker enters SPDa" << endl;
										}
									}
								}
								// S0 not in SPDa
								if (currentEvent::flag_s0_in_SPDa != 1) {
									currentEvent::min_distance = 999999999.99;
									currentEvent::near_domain_ID = 0;
									currentEvent::near_domain_type = 0; // Default on HPD, if SPD more close, modify to 1 during searching
																		// Check all enabled HPD
									for (i = 0;i < HPDOMAIN_INDEX;i++) {
										if (HARD_PD[i + 1].events_type != -1) {
											currentEvent::walker_distance = getHardDomainSurfaceDistance1D(ORI_WALKERS[S0_EVENT[1]], HARD_PD[i + 1]);
											if (currentEvent::walker_distance - currentEvent::min_distance < eps) {
												currentEvent::min_distance = currentEvent::walker_distance;
												currentEvent::near_domain_ID = i + 1;
											}
										}
									}
									// Check all enabled SPDa
									for (i = 0;i < SPDOMAINA_INDEX;i++) {
										if (SOFT_PD_A[i + 1].events_type != -1) {
											currentEvent::walker_distance = getSoftDomainASurfaceDistance1D(ORI_WALKERS[S0_EVENT[1]], SOFT_PD_A[i + 1]);
											if (currentEvent::walker_distance - currentEvent::min_distance < eps) {
												currentEvent::min_distance = currentEvent::walker_distance;
												currentEvent::near_domain_ID = i + 1;
												currentEvent::near_domain_type = 1;
											}
										}
									}
									// Domain, Solver, new event
									// near_domain_type is not a required way to set this part, reserve for other design
									cout << "Create a new HPD during this event." << endl;
									HPDOMAIN_INDEX = HPDOMAIN_INDEX + 1;
									HARD_PD[HPDOMAIN_INDEX].ID = HPDOMAIN_INDEX;
									HARD_PD[HPDOMAIN_INDEX].center.x = ORI_WALKERS[S0_EVENT[1]].center.x;
									HARD_PD[HPDOMAIN_INDEX].center.y = 0.0;
									HARD_PD[HPDOMAIN_INDEX].center.z = 0.0;
									HARD_PD[HPDOMAIN_INDEX].clock = WORLD_CLOCK;
									HARD_PD[HPDOMAIN_INDEX].events_type = 1;
									HARD_PD[HPDOMAIN_INDEX].radius = currentEvent::min_distance;
									HARD_PD[HPDOMAIN_INDEX].walker_ID = ORI_WALKERS[S0_EVENT[1]].ID;
									HARD_PD[HPDOMAIN_INDEX].walker_type = ORI_WALKERS[S0_EVENT[1]].type;
									ORI_WALKERS[S0_EVENT[1]].statue = 2;
									ORI_WALKERS[S0_EVENT[1]].regionID = HPDOMAIN_INDEX;
									// Solver
									HARD_PD[HPDOMAIN_INDEX].events_time = getMajorDiffTimeStamp1D(ORI_WALKERS[HARD_PD[HPDOMAIN_INDEX].walker_ID].diffusivity, HARD_PD[HPDOMAIN_INDEX].radius, ORI_WALKERS[HARD_PD[HPDOMAIN_INDEX].walker_ID].radius);
									// Event
									EVENT_INDEX = EVENT_INDEX + 1;
									EVENT_LIST[EVENT_INDEX].domainID = HARD_PD[HPDOMAIN_INDEX].ID;
									EVENT_LIST[EVENT_INDEX].domainType = 0;
									EVENT_LIST[EVENT_INDEX].timestamp = WORLD_CLOCK + HARD_PD[HPDOMAIN_INDEX].events_time;
									HARD_PD[HPDOMAIN_INDEX].event_card = EVENT_INDEX;
									cout << "1/1 Walker has been protected in HPD" << endl;
								}
							}
						}
						else {
							// General Pass: N > 1, different rebuid method
							// For all S0 walkers in original S0Index
							for (i = 0;i < S0_INDEX;i++) {
								// Avoid the walkers out of boundary
								if (ORI_WALKERS[S0_EVENT[i + 1]].statue == 0 && ORI_WALKERS[S0_EVENT[i + 1]].type != -1) {
									currentEvent::flag_s0_in_SPDa = 0;
									for (j = 0;j < SPDOMAINA_INDEX;j++) {
										if (SOFT_PD_A[j + 1].events_type != -1) {
											currentEvent::walker_distance = getSoftDomainASurfaceDistance1D(ORI_WALKERS[S0_EVENT[1]], SOFT_PD_A[j + 1]);
											if (currentEvent::walker_distance - THRES_2 < eps) {
												currentEvent::flag_s0_in_SPDa = 1;
												// Solver, new event
												// Update SPDa
												// Location adjust is not in this version: TODO
												SOFT_PD_A[j + 1].current_member = SOFT_PD_A[j + 1].current_member + 1;
												SOFT_PD_A[j + 1].member_list.member[SOFT_PD_A[j + 1].current_member] = ORI_WALKERS[S0_EVENT[i + 1]];
												currentEvent::update_spda_time_temp = SOFT_PD_A[j + 1].events_time;
												SOFT_PD_A[j + 1].events_time = SOFT_PD_A[j + 1].events_time; // Replaced with solver update interface once LocalSolver is provided: TODO
												currentEvent::update_spda_time_temp = SOFT_PD_A[j + 1].events_time - currentEvent::update_spda_time_temp;
												ORI_WALKERS[S0_EVENT[i + 1]].statue = 1;
												ORI_WALKERS[S0_EVENT[i + 1]].timestamp = WORLD_CLOCK;
												// No new event, just timestamp update
												EVENT_LIST[SOFT_PD_A[j + 1].event_card].timestamp = EVENT_LIST[SOFT_PD_A[j + 1].event_card].timestamp + currentEvent::update_spda_time_temp;
												cout << "1 Walker enters SPDa" << endl;
											}
										}
									}
								}
							}
							// S0 not in SPDa
							for (i = 0;i < S0_INDEX;i++) {
								// Avoid the walkers out of boundary
								if (ORI_WALKERS[S0_EVENT[i + 1]].statue == 0 && ORI_WALKERS[S0_EVENT[i + 1]].type != -1) {
									for (j = i + 1;j < S0_INDEX;j++) {
										// Avoid the walkers out of boundary
										if (ORI_WALKERS[S0_EVENT[j + 1]].statue == 0 && ORI_WALKERS[S0_EVENT[i + 1]].statue == 0 && ORI_WALKERS[S0_EVENT[j + 1]].type != -1) {
											// surface
											currentEvent::walker_distance = getDistance1D(ORI_WALKERS[S0_EVENT[i + 1]], ORI_WALKERS[S0_EVENT[j + 1]]) - ORI_WALKERS[S0_EVENT[i + 1]].radius - ORI_WALKERS[S0_EVENT[j + 1]].radius;
											if (currentEvent::walker_distance - THRES_2 < eps) {
												// New SPDa: Domain, Solver, new event
												cout << "Create a new Type 1 SPD during this event." << endl;
												SPDOMAINA_INDEX = SPDOMAINA_INDEX + 1;
												SOFT_PD_A[SPDOMAINA_INDEX].ID = SPDOMAINA_INDEX;
												SOFT_PD_A[SPDOMAINA_INDEX].center_distance = currentEvent::walker_distance;
												SOFT_PD_A[SPDOMAINA_INDEX].clock = WORLD_CLOCK;
												SOFT_PD_A[SPDOMAINA_INDEX].events_type = 10;
												SOFT_PD_A[SPDOMAINA_INDEX].current_member = 2;
												SOFT_PD_A[SPDOMAINA_INDEX].member_list.member[1] = ORI_WALKERS[S0_EVENT[i + 1]];
												SOFT_PD_A[SPDOMAINA_INDEX].member_list.member[2] = ORI_WALKERS[S0_EVENT[j + 1]];
												SOFT_PD_A[SPDOMAINA_INDEX].center1.x = ORI_WALKERS[S0_EVENT[i + 1]].center.x;
												SOFT_PD_A[SPDOMAINA_INDEX].center1.y = ORI_WALKERS[S0_EVENT[i + 1]].center.y;
												SOFT_PD_A[SPDOMAINA_INDEX].center1.z = ORI_WALKERS[S0_EVENT[i + 1]].center.z;
												SOFT_PD_A[SPDOMAINA_INDEX].center1.r = ORI_WALKERS[S0_EVENT[i + 1]].radius;
												SOFT_PD_A[SPDOMAINA_INDEX].center2.x = ORI_WALKERS[S0_EVENT[j + 1]].center.x;
												SOFT_PD_A[SPDOMAINA_INDEX].center2.y = ORI_WALKERS[S0_EVENT[j + 1]].center.y;
												SOFT_PD_A[SPDOMAINA_INDEX].center2.z = ORI_WALKERS[S0_EVENT[j + 1]].center.z;
												SOFT_PD_A[SPDOMAINA_INDEX].center2.r = ORI_WALKERS[S0_EVENT[j + 1]].radius;
												ORI_WALKERS[S0_EVENT[i + 1]].statue = 1;
												ORI_WALKERS[S0_EVENT[j + 1]].statue = 1;
												ORI_WALKERS[S0_EVENT[i + 1]].regionID = SPDOMAINA_INDEX;
												ORI_WALKERS[S0_EVENT[j + 1]].regionID = SPDOMAINA_INDEX;
												// Solver
												// RANDOM_SEED = rand() / (double)(RAND_MAX);
												// SOFT_PD_A[SPDOMAINA_INDEX].events_time = INITIAL_EVENTS_TIME + INITIAL_EVENTS_TIME*RANDOM_SEED; // Replaced with solver interface once LocalSolver is provided
												SOFT_PD_A[SPDOMAINA_INDEX].events_time = getReactionTimeStamp1D(SOFT_PD_A[SPDOMAINA_INDEX].member_list, INITIAL_EVENTS_TIME);
												// Event
												EVENT_INDEX = EVENT_INDEX + 1;
												EVENT_LIST[EVENT_INDEX].domainID = SOFT_PD_A[SPDOMAINA_INDEX].ID;
												EVENT_LIST[EVENT_INDEX].domainType = 1;
												EVENT_LIST[EVENT_INDEX].timestamp = WORLD_CLOCK + SOFT_PD_A[SPDOMAINA_INDEX].events_time;
												SOFT_PD_A[SPDOMAINA_INDEX].event_card = EVENT_INDEX;
												cout << "2 walkers built a new SPDa" << endl;
												cout << ORI_WALKERS[S0_EVENT[i + 1]].ID << " and " << ORI_WALKERS[S0_EVENT[j + 1]].ID << endl;
											}
										}
									}
								}
							}
							for (i = 0;i < S0_INDEX;i++) {
								// Avoid the walkers out of boundary
								if (ORI_WALKERS[S0_EVENT[i + 1]].statue == 0 && ORI_WALKERS[S0_EVENT[i + 1]].type != -1) {
									// Three routes
									// Close S0 walkers: O1
									// Not close but maximum and uncover case (dynamic) is not included in this version
									currentEvent::O1_allow_r = 0.0;
									currentEvent::O1_neighbout_ID = 0;
									currentEvent::O2_allow_r = 0.0;
									currentEvent::O3_allow_r = 0.0;
									currentEvent::min_distance = EXTRA_LARGE;
									for (j = i + 1;j < S0_INDEX;j++) {
										// Avoid the walkers out of boundary
										if (ORI_WALKERS[S0_EVENT[j + 1]].statue == 0 && ORI_WALKERS[S0_EVENT[j + 1]].type != -1 && ORI_WALKERS[S0_EVENT[j + 1]].type != -1) {
											currentEvent::walker_distance = getDistance1D(ORI_WALKERS[S0_EVENT[i + 1]], ORI_WALKERS[S0_EVENT[j + 1]]) - ORI_WALKERS[S0_EVENT[i + 1]].radius - ORI_WALKERS[S0_EVENT[j + 1]].radius;
											if (currentEvent::walker_distance - currentEvent::min_distance < eps) {
												currentEvent::min_distance = currentEvent::walker_distance;
												currentEvent::O1_neighbout_ID = S0_EVENT[j + 1];
												currentEvent::O1_allow_r = getDistance1D(ORI_WALKERS[S0_EVENT[i + 1]], ORI_WALKERS[S0_EVENT[j + 1]])*(ORI_WALKERS[S0_EVENT[i + 1]].diffusivity) / (ORI_WALKERS[S0_EVENT[i + 1]].diffusivity + ORI_WALKERS[S0_EVENT[j + 1]].diffusivity) - ORI_WALKERS[S0_EVENT[i + 1]].radius;
												//currentEvent::O1_neighour_allow_r = getDistance1D(ORI_WALKERS[S0_EVENT[i + 1]], ORI_WALKERS[S0_EVENT[j + 1]])*(ORI_WALKERS[S0_EVENT[j + 1]].diffusivity) / (ORI_WALKERS[S0_EVENT[i + 1]].diffusivity + ORI_WALKERS[S0_EVENT[j + 1]].diffusivity) - ORI_WALKERS[S0_EVENT[j + 1]].radius;
											}
										}
									}
									// HPD surface: O2
									currentEvent::min_distance = EXTRA_LARGE;
									for (j = 0;j < HPDOMAIN_INDEX;j++) {
										if (HARD_PD[j + 1].events_type != -1) {
											currentEvent::walker_distance = getHardDomainSurfaceDistance1D(ORI_WALKERS[S0_EVENT[i + 1]], HARD_PD[j + 1]);
											if (currentEvent::walker_distance - currentEvent::min_distance < eps) {
												currentEvent::min_distance = currentEvent::walker_distance;
												currentEvent::O2_allow_r = currentEvent::walker_distance; // Surface to surface = center to surface - r
											}
										}
									}
									// SPDa close center: 03
									currentEvent::min_distance = EXTRA_LARGE;
									for (j = 0;j < SPDOMAINA_INDEX;j++) {
										if (SOFT_PD_A[j + 1].events_type != -1) {
											currentEvent::walker_distance = getSoftDomainACentertoCenterDistance1D(ORI_WALKERS[S0_EVENT[i + 1]], SOFT_PD_A[j + 1]) - ORI_WALKERS[S0_EVENT[i + 1]].radius;
											if (currentEvent::walker_distance - currentEvent::min_distance < eps) {
												currentEvent::min_distance = currentEvent::walker_distance;
												currentEvent::O3_allow_r = currentEvent::walker_distance; // Walker surface to SPDa center A
											}
										}
									}
									// Use the minimum none 0 one to get 3 different routes of domains, etc
									if (currentEvent::O1_allow_r - currentEvent::O2_allow_r < eps && currentEvent::O1_allow_r > eps) {
										// O1 < O2 - > O1/O3
										if (currentEvent::O1_allow_r - currentEvent::O3_allow_r < eps) {
											// O1 < O3 - > O1
											// Domain
											cout << "Create a new HPD during this event." << endl;
											HPDOMAIN_INDEX = HPDOMAIN_INDEX + 1;
											HARD_PD[HPDOMAIN_INDEX].ID = HPDOMAIN_INDEX;
											HARD_PD[HPDOMAIN_INDEX].center.x = ORI_WALKERS[S0_EVENT[i + 1]].center.x;
											HARD_PD[HPDOMAIN_INDEX].center.y = 0.0;
											HARD_PD[HPDOMAIN_INDEX].center.z = 0.0;
											HARD_PD[HPDOMAIN_INDEX].clock = WORLD_CLOCK;
											HARD_PD[HPDOMAIN_INDEX].events_type = 1;
											HARD_PD[HPDOMAIN_INDEX].radius = currentEvent::O1_allow_r + ORI_WALKERS[S0_EVENT[i + 1]].radius;
											HARD_PD[HPDOMAIN_INDEX].walker_ID = ORI_WALKERS[S0_EVENT[i + 1]].ID;
											HARD_PD[HPDOMAIN_INDEX].walker_type = ORI_WALKERS[S0_EVENT[i + 1]].type;
											ORI_WALKERS[S0_EVENT[i + 1]].statue = 2;
											ORI_WALKERS[S0_EVENT[i + 1]].regionID = HPDOMAIN_INDEX;
											// Solver
											HARD_PD[HPDOMAIN_INDEX].events_time = getMajorDiffTimeStamp1D(ORI_WALKERS[HARD_PD[HPDOMAIN_INDEX].walker_ID].diffusivity, HARD_PD[HPDOMAIN_INDEX].radius, ORI_WALKERS[HARD_PD[HPDOMAIN_INDEX].walker_ID].radius);
											// Event
											EVENT_INDEX = EVENT_INDEX + 1;
											EVENT_LIST[EVENT_INDEX].domainID = HARD_PD[HPDOMAIN_INDEX].ID;
											EVENT_LIST[EVENT_INDEX].domainType = 0;
											EVENT_LIST[EVENT_INDEX].timestamp = WORLD_CLOCK + HARD_PD[HPDOMAIN_INDEX].events_time;
											HARD_PD[HPDOMAIN_INDEX].event_card = EVENT_INDEX;
											cout << "1 Walker has been protected in HPD" << endl;
										}
										else {
											// O3 < O1 - > O3 or O1 = 0
											// Domain
											cout << "Create a new HPD during this event." << endl;
											HPDOMAIN_INDEX = HPDOMAIN_INDEX + 1;
											HARD_PD[HPDOMAIN_INDEX].ID = HPDOMAIN_INDEX;
											HARD_PD[HPDOMAIN_INDEX].center.x = ORI_WALKERS[S0_EVENT[i + 1]].center.x;
											HARD_PD[HPDOMAIN_INDEX].center.y = 0.0;
											HARD_PD[HPDOMAIN_INDEX].center.z = 0.0;
											HARD_PD[HPDOMAIN_INDEX].clock = WORLD_CLOCK;
											HARD_PD[HPDOMAIN_INDEX].events_type = 1;
											HARD_PD[HPDOMAIN_INDEX].radius = currentEvent::O3_allow_r + ORI_WALKERS[S0_EVENT[i + 1]].radius;
											HARD_PD[HPDOMAIN_INDEX].walker_ID = ORI_WALKERS[S0_EVENT[i + 1]].ID;
											HARD_PD[HPDOMAIN_INDEX].walker_type = ORI_WALKERS[S0_EVENT[i + 1]].type;
											ORI_WALKERS[S0_EVENT[i + 1]].statue = 2;
											ORI_WALKERS[S0_EVENT[i + 1]].regionID = HPDOMAIN_INDEX;
											// Solver
											HARD_PD[HPDOMAIN_INDEX].events_time = getMajorDiffTimeStamp1D(ORI_WALKERS[HARD_PD[HPDOMAIN_INDEX].walker_ID].diffusivity, HARD_PD[HPDOMAIN_INDEX].radius, ORI_WALKERS[HARD_PD[HPDOMAIN_INDEX].walker_ID].radius);
											// Event
											EVENT_INDEX = EVENT_INDEX + 1;
											EVENT_LIST[EVENT_INDEX].domainID = HARD_PD[HPDOMAIN_INDEX].ID;
											EVENT_LIST[EVENT_INDEX].domainType = 0;
											EVENT_LIST[EVENT_INDEX].timestamp = WORLD_CLOCK + HARD_PD[HPDOMAIN_INDEX].events_time;
											HARD_PD[HPDOMAIN_INDEX].event_card = EVENT_INDEX;
											cout << "1 Walker has been protected in HPD" << endl;
										}
									}
									else {
										// O2 < O1 - > O2/O3
										if (currentEvent::O2_allow_r - currentEvent::O3_allow_r < eps) {
											// O2 < O3 - > O2
											// Domain
											cout << "Create a new HPD during this event." << endl;
											HPDOMAIN_INDEX = HPDOMAIN_INDEX + 1;
											HARD_PD[HPDOMAIN_INDEX].ID = HPDOMAIN_INDEX;
											HARD_PD[HPDOMAIN_INDEX].center.x = ORI_WALKERS[S0_EVENT[i + 1]].center.x;
											HARD_PD[HPDOMAIN_INDEX].center.y = 0.0;
											HARD_PD[HPDOMAIN_INDEX].center.z = 0.0;
											HARD_PD[HPDOMAIN_INDEX].clock = WORLD_CLOCK;
											HARD_PD[HPDOMAIN_INDEX].events_type = 1;
											HARD_PD[HPDOMAIN_INDEX].radius = currentEvent::O2_allow_r + ORI_WALKERS[S0_EVENT[i + 1]].radius;
											HARD_PD[HPDOMAIN_INDEX].walker_ID = ORI_WALKERS[S0_EVENT[i + 1]].ID;
											HARD_PD[HPDOMAIN_INDEX].walker_type = ORI_WALKERS[S0_EVENT[i + 1]].type;
											ORI_WALKERS[S0_EVENT[i + 1]].statue = 2;
											ORI_WALKERS[S0_EVENT[i + 1]].regionID = HPDOMAIN_INDEX;
											// Solver
											HARD_PD[HPDOMAIN_INDEX].events_time = getMajorDiffTimeStamp1D(ORI_WALKERS[HARD_PD[HPDOMAIN_INDEX].walker_ID].diffusivity, HARD_PD[HPDOMAIN_INDEX].radius, ORI_WALKERS[HARD_PD[HPDOMAIN_INDEX].walker_ID].radius);
											// Event
											EVENT_INDEX = EVENT_INDEX + 1;
											EVENT_LIST[EVENT_INDEX].domainID = HARD_PD[HPDOMAIN_INDEX].ID;
											EVENT_LIST[EVENT_INDEX].domainType = 0;
											EVENT_LIST[EVENT_INDEX].timestamp = WORLD_CLOCK + HARD_PD[HPDOMAIN_INDEX].events_time;
											HARD_PD[HPDOMAIN_INDEX].event_card = EVENT_INDEX;
											cout << "1 Walker has been protected in HPD" << endl;
										}
										else {
											// O3 < O2 - > O3
											// Domain
											cout << "Create a new HPD during this event." << endl;
											HPDOMAIN_INDEX = HPDOMAIN_INDEX + 1;
											HARD_PD[HPDOMAIN_INDEX].ID = HPDOMAIN_INDEX;
											HARD_PD[HPDOMAIN_INDEX].center.x = ORI_WALKERS[S0_EVENT[i + 1]].center.x;
											HARD_PD[HPDOMAIN_INDEX].center.y = 0.0;
											HARD_PD[HPDOMAIN_INDEX].center.z = 0.0;
											HARD_PD[HPDOMAIN_INDEX].clock = WORLD_CLOCK;
											HARD_PD[HPDOMAIN_INDEX].events_type = 1;
											HARD_PD[HPDOMAIN_INDEX].radius = currentEvent::O3_allow_r + ORI_WALKERS[S0_EVENT[i + 1]].radius;
											HARD_PD[HPDOMAIN_INDEX].walker_ID = ORI_WALKERS[S0_EVENT[i + 1]].ID;
											HARD_PD[HPDOMAIN_INDEX].walker_type = ORI_WALKERS[S0_EVENT[i + 1]].type;
											ORI_WALKERS[S0_EVENT[i + 1]].statue = 2;
											ORI_WALKERS[S0_EVENT[i + 1]].regionID = HPDOMAIN_INDEX;
											// Solver
											HARD_PD[HPDOMAIN_INDEX].events_time = getMajorDiffTimeStamp1D(ORI_WALKERS[HARD_PD[HPDOMAIN_INDEX].walker_ID].diffusivity, HARD_PD[HPDOMAIN_INDEX].radius, ORI_WALKERS[HARD_PD[HPDOMAIN_INDEX].walker_ID].radius);
											// Event
											EVENT_INDEX = EVENT_INDEX + 1;
											EVENT_LIST[EVENT_INDEX].domainID = HARD_PD[HPDOMAIN_INDEX].ID;
											EVENT_LIST[EVENT_INDEX].domainType = 0;
											EVENT_LIST[EVENT_INDEX].timestamp = WORLD_CLOCK + HARD_PD[HPDOMAIN_INDEX].events_time;
											HARD_PD[HPDOMAIN_INDEX].event_card = EVENT_INDEX;
											cout << "1 Walker has been protected in HPD" << endl;
										}
									}
								}
							}
						} // End of N > 1 else and all N if...else
					}
					if (currentEvent::event_type == 11) {
						// Move EVENT_POINTER to next event
						EVENT_POINTER = EVENT_POINTER + 1;
						cout << "FATAL ERROR 700: Unable to process E11 for HPD in this version. Exit." << endl;
						system("pause");
						exit(0);
					}
				}
				if (currentEvent::domain_type == 1) {
					// SPDa event
					currentEvent::event_type = SOFT_PD_A[currentEvent::domainID].events_type;
					// SPD concern E10 only, this is not in input of update
					// Only annihilation in this version, 
					// Move EVENT_POINTER to next event
					EVENT_POINTER = EVENT_POINTER + 1;
					// Break SPDa
					// Disable the domain
					SOFT_PD_A[currentEvent::domainID].events_type = -1;
					// Call solver to get location, send walker, set walker to S0 (No need)
					// Update clock
					// Disable all walker (type to -1)
					for (i = 0;i < SOFT_PD_A[currentEvent::domainID].current_member;i++) {
						ORI_WALKERS[SOFT_PD_A[currentEvent::domainID].member_list.member[i + 1].ID].statue = 0;
						ORI_WALKERS[SOFT_PD_A[currentEvent::domainID].member_list.member[i + 1].ID].type = -1;
					}
					WORLD_CLOCK = currentEvent::timestamp;

					// Enter ASMS 1st Stage
					// Check location and neighbour (No need)
					// Protect walker again, build new domain or update SPDa (No need)
					// Set new events and run sort (no nned)
				}
				if (currentEvent::domain_type == 2) {
					// Move EVENT_POINTER to next event
					EVENT_POINTER = EVENT_POINTER + 1;
					cout << "FATAL ERROR 700: Unable to process Type 2 SPD in this version. Exit." << endl;
					system("pause");
					exit(0);
				}
				// Run event list sort
				// Use malloc with pointer instead of assigned array: TODO
				// Use type queue instead of type array: TODO
				for (i = 1; i < EVENT_INDEX;i++) {
					for (j = 1;j < EVENT_INDEX + 1 - i;j++) {
						if (EVENT_LIST[j].timestamp>EVENT_LIST[j + 1].timestamp) {
							swap(EVENT_LIST[j], EVENT_LIST[j + 1]);
						}
					}
				}
				for (i = 1;i < EVENT_INDEX;i++) {
					if (EVENT_LIST[i].domainType == 0) {
						HARD_PD[EVENT_LIST[i].domainID].event_card = i;
					}
					if (EVENT_LIST[i].domainType == 1) {
						SOFT_PD_A[EVENT_LIST[i].domainID].event_card = i;
					}
				}
				// One event has been processed with its executor and creator
				// Check the clock
				if (WORLD_CLOCK - END_CLOCK > eps) {
					CONTROL_FLAG = 0; // Check with clock
				}
				// Inner Dump Debug
				if (DUMPFLAG_INN) {
					ofstream file2("InnerDump.txt", ios::app);
					streambuf *f2 = cout.rdbuf(file2.rdbuf());
					DUMP_COUNTER = 0;
					cout << WORLD_CLOCK << endl;
					for (i = 1;i < WALKER_INDEX + 1;i++) {
						if (ORI_WALKERS[i].type != -1) {
							DUMP_COUNTER = DUMP_COUNTER + 1;
							if (space_dimension == 3) {
								cout << ORI_WALKERS[i].type << " " << ORI_WALKERS[i].component << " " << ORI_WALKERS[i].center.x << " " << ORI_WALKERS[i].center.y << " " << ORI_WALKERS[i].center.z << endl;
							}
							if (space_dimension == 1) {
								//cout << ORI_WALKERS[i].type << " " << ORI_WALKERS[i].component << " " << ORI_WALKERS[i].center.x << endl;
								cout << ORI_WALKERS[i].ID << " " << ORI_WALKERS[i].statue << " " << ORI_WALKERS[i].center.x << " 0.0 0.0" << endl;
								if (ORI_WALKERS[i].center.x - dump_min < eps) {
									dump_min = ORI_WALKERS[i].center.x;
								}
								if (ORI_WALKERS[i].center.x - dump_max > eps) {
									dump_max = ORI_WALKERS[i].center.x;
								}
							}
						}
					}
					//cout << "Remain " << DUMP_COUNTER << " walkers." << endl;
					//cout << "  " << DUMP_COUNTER << endl;
					cout << endl;
					cout.rdbuf(f2);
				}
				// OS2
				ofstream file3("timelog.txt", ios::app);
				streambuf *f3 = cout.rdbuf(file3.rdbuf());
				cout << WORLD_CLOCK << " " << DUMP_COUNTER << " " << dump_min << " " << dump_max << endl;
				cout.rdbuf(f3);
				// End of inner dump
			}
			STAGE_FLAG = 0;
			cout << "Current calculation has finished at " << WORLD_CLOCK << "s, aCRD has reached the end of clock, extend the calculation? 1:0." << endl;
			cin >> INPUT_FLAG;
			if (INPUT_FLAG) {
				cout << "Please input the extension time in the unit of second. " << endl;
				cin >> EXTRA_CLOCK;
				END_CLOCK = END_CLOCK + EXTRA_CLOCK;
				STAGE_FLAG = 1;
			}
		}
	}
	// 3 Dimension begins here
	if (space_dimension == 3) {
		while (STAGE_FLAG) {
			cout << "Running 3 dimension events." << endl;
			if (WORLD_CLOCK - END_CLOCK < eps) {
				CONTROL_FLAG = 1; // Check with clock
			}
			if (EVENT_POINTER > EVENT_INDEX) {
				CONTROL_FLAG = 0; // All events finished
			}
			while (CONTROL_FLAG) {
				// Begin to execute event
				// Read event card
				// Extra judge with E20 Event
				if (DISPLAYFLAG_SYS == 1) {
					cout << "Begin a new event." << endl;
				}
				currentEvent::domain_type = EVENT_LIST[EVENT_POINTER].domainType;
				currentEvent::domainID = EVENT_LIST[EVENT_POINTER].domainID;
				currentEvent::timestamp = EVENT_LIST[EVENT_POINTER].timestamp;
				currentEvent::event_type = EVENT_LIST[EVENT_POINTER].eventType;
				WORLD_CLOCK = currentEvent::timestamp;
				// Execute E20 first, if not E20, do the rest
				if (currentEvent::event_type == 20) {
					if (DISPLAYFLAG_SYS == 1) {
						cout << "Begin E20 event." << endl;
					}
					EVENT_POINTER = EVENT_POINTER + 1;
					//system("pause");
					// Create new walkers at a point out of HPD, can be in SPD
					// New walkers all in S0
					// Update WALKER_COUNTER and S0 log
					S0_INDEX = 0;
					for (i = 0;i < irradiationINumber;i++) {
						if (DISPLAYFLAG_SYS) {
							cout << "Create 1 new walker in this process." << endl;
						}
						WALKER_COUNTER = WALKER_COUNTER + 1;
						ORI_WALKERS[WALKER_COUNTER].ID = WALKER_COUNTER;
						ORI_WALKERS[WALKER_COUNTER].timestamp = WORLD_CLOCK;
						ORI_WALKERS[WALKER_COUNTER].type = 1;
						ORI_WALKERS[WALKER_COUNTER].component = 1;
						ORI_WALKERS[WALKER_COUNTER].statue = 0;
						ORI_WALKERS[WALKER_COUNTER].energy = 0.0; // Not in use for current version
						ORI_WALKERS[WALKER_COUNTER].radius = getEstimateRadius(ORI_WALKERS[WALKER_COUNTER]);
						ORI_WALKERS[WALKER_COUNTER].diffusivity = getDvValue(ORI_WALKERS[WALKER_COUNTER].energy, ORI_WALKERS[WALKER_COUNTER].type, ORI_WALKERS[WALKER_COUNTER].component);
						// Out of HPD, can be in SPD
						// Search all walkers statue 0 or 1, type != -1, prevent overlap
						// Bypass ID=WALKER_COUNTER, this is yourself
						LOCATE_FLAG = 1;
						while (LOCATE_FLAG)
						{
							// Default regard as no overlap, may change during check
							LOCATE_FLAG = 0;
							// Get random location
							RANDOM_SEED = rand() / (double)(RAND_MAX); // Random SEED [0,1)
							ORI_WALKERS[WALKER_COUNTER].center.x = X_LOWWER_LIMIT + (X_UPPER_LIMIT - X_LOWWER_LIMIT)*RANDOM_SEED;
							RANDOM_SEED = rand() / (double)(RAND_MAX); // Random SEED [0,1)
							ORI_WALKERS[WALKER_COUNTER].center.y = Y_LOWWER_LIMIT + (Y_UPPER_LIMIT - Y_LOWWER_LIMIT)*RANDOM_SEED;
							RANDOM_SEED = rand() / (double)(RAND_MAX); // Random SEED [0,1)
							ORI_WALKERS[WALKER_COUNTER].center.z = Z_LOWWER_LIMIT + (Z_UPPER_LIMIT - Z_LOWWER_LIMIT)*RANDOM_SEED;
							if (DISPLAYFLAG_DEB) {
								cout << "Randomly get 1 interstitial at " << ORI_WALKERS[WALKER_COUNTER].center.x << " " << ORI_WALKERS[WALKER_COUNTER].center.y << " " << ORI_WALKERS[WALKER_COUNTER].center.z << endl;
							}
							// Begin check, if overlap, LOCATE_FLAG=1
							for (j = 0;j < HPDOMAIN_INDEX;j++) {
								// check HPD, no overlap
								if (HARD_PD[j + 1].events_type != -1) {
									new_walker_distance = getHardDomainSurfaceDistance3D(ORI_WALKERS[WALKER_COUNTER], HARD_PD[j + 1]);
									if (new_walker_distance < eps) {
										LOCATE_FLAG = 1;
										if (DISPLAYFLAG_DEB) {
											cout << "Surface distance " << new_walker_distance << ", in side of a HPD, resample." << endl;
										}
									}
								}
							}
							for (j = 1;j < WALKER_COUNTER;j++) {
								// check walkers with S0 and S1, no overlap
								if (ORI_WALKERS[j].type != -1) {
									if (ORI_WALKERS[j].statue == 1 || ORI_WALKERS[j].statue == 0) {
										new_walker_distance = getDistance(ORI_WALKERS[WALKER_COUNTER], ORI_WALKERS[j]) - ORI_WALKERS[WALKER_COUNTER].radius - ORI_WALKERS[j].radius;
										if (new_walker_distance < eps) {
											LOCATE_FLAG = 1;
											if (DISPLAYFLAG_DEB) {
												cout << "Surface distance " << new_walker_distance << ", walker overlap, resample." << endl;
											}
										}
									}
								}
							}
							// End of check
						}// Get good location
						 // Update S0 log						
						S0_INDEX = S0_INDEX + 1;
						S0_EVENT[S0_INDEX] = ORI_WALKERS[WALKER_COUNTER].ID; // Log the IDs of all walkers with S0
					} // End of I
					for (i = 0;i < irradiationVNumber;i++) {
						if (DISPLAYFLAG_SYS) {
							cout << "Create 1 new walker in this process." << endl;
						}
						WALKER_COUNTER = WALKER_COUNTER + 1;
						ORI_WALKERS[WALKER_COUNTER].ID = WALKER_COUNTER;
						ORI_WALKERS[WALKER_COUNTER].timestamp = WORLD_CLOCK;
						ORI_WALKERS[WALKER_COUNTER].type = 0;
						ORI_WALKERS[WALKER_COUNTER].component = 1;
						ORI_WALKERS[WALKER_COUNTER].statue = 0;
						ORI_WALKERS[WALKER_COUNTER].energy = 0.0; // Not in use for current version
						ORI_WALKERS[WALKER_COUNTER].radius = getEstimateRadius(ORI_WALKERS[WALKER_COUNTER]);
						ORI_WALKERS[WALKER_COUNTER].diffusivity = getDvValue(ORI_WALKERS[WALKER_COUNTER].energy, ORI_WALKERS[WALKER_COUNTER].type, ORI_WALKERS[WALKER_COUNTER].component);
						// Out of HPD, can be in SPD
						// Search all walkers statue 0 or 1, type != -1, prevent overlap
						// Bypass ID=WALKER_COUNTER, this is yourself
						LOCATE_FLAG = 1;
						while (LOCATE_FLAG)
						{
							// Default regard as no overlap, may change during check
							LOCATE_FLAG = 0;
							// Get random location
							RANDOM_SEED = rand() / (double)(RAND_MAX); // Random SEED [0,1)
							ORI_WALKERS[WALKER_COUNTER].center.x = X_LOWWER_LIMIT + (X_UPPER_LIMIT - X_LOWWER_LIMIT)*RANDOM_SEED;
							RANDOM_SEED = rand() / (double)(RAND_MAX); // Random SEED [0,1)
							ORI_WALKERS[WALKER_COUNTER].center.y = Y_LOWWER_LIMIT + (Y_UPPER_LIMIT - Y_LOWWER_LIMIT)*RANDOM_SEED;
							RANDOM_SEED = rand() / (double)(RAND_MAX); // Random SEED [0,1)
							ORI_WALKERS[WALKER_COUNTER].center.z = Z_LOWWER_LIMIT + (Z_UPPER_LIMIT - Z_LOWWER_LIMIT)*RANDOM_SEED;
							if (DISPLAYFLAG_DEB) {
								cout << "Randomly get 1 vacancy at " << ORI_WALKERS[WALKER_COUNTER].center.x << " " << ORI_WALKERS[WALKER_COUNTER].center.y << " " << ORI_WALKERS[WALKER_COUNTER].center.z << endl;
							}
							// Begin check, if overlap, LOCATE_FLAG=1
							for (j = 0;j < HPDOMAIN_INDEX;j++) {
								// check HPD, no overlap
								if (HARD_PD[j + 1].events_type != -1) {
									new_walker_distance = getHardDomainSurfaceDistance3D(ORI_WALKERS[WALKER_COUNTER], HARD_PD[j + 1]);
									if (new_walker_distance < eps) {
										LOCATE_FLAG = 1;
										if (DISPLAYFLAG_DEB) {
											cout << "Surface distance " << new_walker_distance << ", in side of a HPD, resample." << endl;
										}
									}
								}
							}
							for (j = 1;j < WALKER_COUNTER;j++) {
								// check walkers with S0 and S1, no overlap
								if (ORI_WALKERS[j].type != -1) {
									if (ORI_WALKERS[j].statue == 1 || ORI_WALKERS[j].statue == 0) {
										new_walker_distance = getDistance(ORI_WALKERS[WALKER_COUNTER], ORI_WALKERS[j]) - ORI_WALKERS[WALKER_COUNTER].radius - ORI_WALKERS[j].radius;
										if (new_walker_distance < eps) {
											LOCATE_FLAG = 1;
											if (DISPLAYFLAG_DEB) {
												cout << "Surface distance " << new_walker_distance << ", walker overlap, resample." << endl;
											}
										}
									}
								}
							}
							// End of check
						}// Get good location
						S0_INDEX = S0_INDEX + 1;
						S0_EVENT[S0_INDEX] = ORI_WALKERS[WALKER_COUNTER].ID; // Log the IDs of all walkers with S0
					} // End of V
					  // End of create walkers
					  // Enter ASMS 1st Stage
					S0_FULL_FLAG = 1;
					while (S0_FULL_FLAG) {
						// Search all exist S0 walkers
						for (i = 0;i < S0_INDEX;i++) {
							S0_FULL_FLAG = 0; // Default no need for next round, if any E02 is triggered, then set back to 1
											  // walker: ORI_WALKERS[S0_EVENT[i+1]]
							for (j = 0;j < HPDOMAIN_INDEX;j++) {
								if (HARD_PD[j + 1].events_type != -1) {
									// Get distance from current S0 walker to each HPD surface, check whether need E02 or not (TH1)
									if (getHardDomainSurfaceDistance3D(ORI_WALKERS[S0_EVENT[i + 1]], HARD_PD[j + 1]) - THRES_1 < eps) {
										// E02 is trigged for this HPD
										ORI_WALKERS[HARD_PD[j + 1].walker_ID].center = getTransDiffRelativeLocation3D(ORI_WALKERS[HARD_PD[j + 1].walker_ID].diffusivity, WORLD_CLOCK - ORI_WALKERS[HARD_PD[j + 1].walker_ID].timestamp, HARD_PD[j + 1].radius, ORI_WALKERS[HARD_PD[j + 1].walker_ID].radius, HARD_PD[j + 1].center);
										// Break, walker to S0, index + 1
										HARD_PD[j + 1].events_type = -1;
										ORI_WALKERS[HARD_PD[j + 1].walker_ID].statue = 0;
										S0_INDEX = S0_INDEX + 1;
										S0_EVENT[S0_INDEX] = HARD_PD[j + 1].walker_ID;
										S0_FULL_FLAG = 1;
									}
								}
								// Then check next HPD
							}
							// Finish check for the triggered S0, no HPD is still close, to next round: check new full S0 list
						}
					}
					cout << "Overall " << S0_INDEX << " Walkers Released in ASMS 1st Stage" << endl;
					//system("pause");
					for (i = 0;i < S0_INDEX;i++) {
						// X boundary
						if (ORI_WALKERS[S0_EVENT[i + 1]].statue == 0 && ORI_WALKERS[S0_EVENT[i + 1]].type != -1) {
							if (ORI_WALKERS[S0_EVENT[i + 1]].center.x + ORI_WALKERS[S0_EVENT[i + 1]].radius - X_UPPER_LIMIT>eps || ORI_WALKERS[S0_EVENT[i + 1]].center.x - ORI_WALKERS[S0_EVENT[i + 1]].radius - X_LOWWER_LIMIT < eps) {
								ORI_WALKERS[S0_EVENT[i + 1]].type = -1; // Inable this particle
								cout << "1 walker exits the X-axis boundary" << endl;
							}
						}
						// Y boundary
						if (ORI_WALKERS[S0_EVENT[i + 1]].statue == 0 && ORI_WALKERS[S0_EVENT[i + 1]].type != -1) {
							if (ORI_WALKERS[S0_EVENT[i + 1]].center.y + ORI_WALKERS[S0_EVENT[i + 1]].radius - Y_UPPER_LIMIT>eps || ORI_WALKERS[S0_EVENT[i + 1]].center.y - ORI_WALKERS[S0_EVENT[i + 1]].radius - Y_LOWWER_LIMIT < eps) {
								ORI_WALKERS[S0_EVENT[i + 1]].type = -1; // Inable this particle
								cout << "1 walker exits the Y-axis boundary" << endl;
							}
						}
						// Z boundary
						if (ORI_WALKERS[S0_EVENT[i + 1]].statue == 0 && ORI_WALKERS[S0_EVENT[i + 1]].type != -1) {
							if (ORI_WALKERS[S0_EVENT[i + 1]].center.z + ORI_WALKERS[S0_EVENT[i + 1]].radius - Z_UPPER_LIMIT>eps || ORI_WALKERS[S0_EVENT[i + 1]].center.z - ORI_WALKERS[S0_EVENT[i + 1]].radius - Z_LOWWER_LIMIT < eps) {
								ORI_WALKERS[S0_EVENT[i + 1]].type = -1; // Inable this particle
								cout << "1 walker exits the Z-axis boundary" << endl;
							}
						}
					}
					// Enter ASMS 2nd Stage
					// Fast pass: N = 1
					if (S0_INDEX == 1) {
						// Avoid the walkers out of boundary
						if (ORI_WALKERS[S0_EVENT[1]].type != -1) {
							// Whether this S0 walker located in SPDa
							currentEvent::flag_s0_in_SPDa = 0;
							for (i = 0;i < SPDOMAINA_INDEX;i++) {
								if (SOFT_PD_A[i + 1].events_type != -1) {
									currentEvent::walker_distance = getSoftDomainASurfaceDistance(ORI_WALKERS[S0_EVENT[1]], SOFT_PD_A[i + 1]);
									if (currentEvent::walker_distance - THRES_2 < eps) {
										currentEvent::flag_s0_in_SPDa = 1;
										// Solver, new event
										// Update SPDa
										// Location adjust is not in this version
										SOFT_PD_A[i + 1].current_member = SOFT_PD_A[i + 1].current_member + 1;
										SOFT_PD_A[i + 1].member_list.member[SOFT_PD_A[i + 1].current_member] = ORI_WALKERS[S0_EVENT[1]];
										currentEvent::update_spda_time_temp = SOFT_PD_A[i + 1].events_time;
										SOFT_PD_A[i + 1].events_time = SOFT_PD_A[i + 1].events_time; // Replaced with solver update interface once LocalSolver is provided: Not in this version
										currentEvent::update_spda_time_temp = SOFT_PD_A[i + 1].events_time - currentEvent::update_spda_time_temp;
										ORI_WALKERS[S0_EVENT[1]].statue = 1;
										ORI_WALKERS[S0_EVENT[1]].timestamp = WORLD_CLOCK;
										// No new event, just timestamp update
										EVENT_LIST[SOFT_PD_A[i + 1].event_card].timestamp = EVENT_LIST[SOFT_PD_A[i + 1].event_card].timestamp + currentEvent::update_spda_time_temp;
										cout << "1/1 Walker enters SPDa" << endl;
									}
								}
							}
							// S0 not in SPDa
							if (currentEvent::flag_s0_in_SPDa != 1) {
								currentEvent::min_distance = EXTRA_LARGE;
								currentEvent::near_domain_ID = 0;
								currentEvent::near_domain_type = 0; // Default on HPD, if SPD more close, modify to 1 during searching
																	// Check all enabled HPD
								for (i = 0;i < HPDOMAIN_INDEX;i++) {
									if (HARD_PD[i + 1].events_type != -1) {
										currentEvent::walker_distance = getHardDomainSurfaceDistance3D(ORI_WALKERS[S0_EVENT[1]], HARD_PD[i + 1]);
										if (currentEvent::walker_distance - currentEvent::min_distance < eps) {
											currentEvent::min_distance = currentEvent::walker_distance;
											currentEvent::near_domain_ID = i + 1;
										}
									}
								}
								// Check all enabled SPDa
								for (i = 0;i < SPDOMAINA_INDEX;i++) {
									if (SOFT_PD_A[i + 1].events_type != -1) {
										currentEvent::walker_distance = getSoftDomainASurfaceDistance(ORI_WALKERS[S0_EVENT[1]], SOFT_PD_A[i + 1]);
										if (currentEvent::walker_distance - currentEvent::min_distance < eps) {
											currentEvent::min_distance = currentEvent::walker_distance;
											currentEvent::near_domain_ID = i + 1;
											currentEvent::near_domain_type = 1;
										}
									}
								}
								// Domain, Solver, new event
								// near_domain_type is not a required way to set this part, reserve for other design
								if (DISPLAYFLAG_DEB == 1) {
									cout << currentEvent::min_distance << endl;
								}
								cout << "Create a new HPD during this event." << endl;
								HPDOMAIN_INDEX = HPDOMAIN_INDEX + 1;
								HARD_PD[HPDOMAIN_INDEX].ID = HPDOMAIN_INDEX;
								HARD_PD[HPDOMAIN_INDEX].center.x = ORI_WALKERS[S0_EVENT[1]].center.x;
								HARD_PD[HPDOMAIN_INDEX].center.y = ORI_WALKERS[S0_EVENT[1]].center.y;
								HARD_PD[HPDOMAIN_INDEX].center.z = ORI_WALKERS[S0_EVENT[1]].center.z;
								HARD_PD[HPDOMAIN_INDEX].clock = WORLD_CLOCK;
								HARD_PD[HPDOMAIN_INDEX].events_type = 1;
								HARD_PD[HPDOMAIN_INDEX].radius = currentEvent::min_distance + ORI_WALKERS[S0_EVENT[1]].radius;
								HARD_PD[HPDOMAIN_INDEX].walker_ID = ORI_WALKERS[S0_EVENT[1]].ID;
								HARD_PD[HPDOMAIN_INDEX].walker_type = ORI_WALKERS[S0_EVENT[1]].type;
								ORI_WALKERS[S0_EVENT[1]].statue = 2;
								ORI_WALKERS[S0_EVENT[1]].regionID = HPDOMAIN_INDEX;
								// Solver
								HARD_PD[HPDOMAIN_INDEX].events_time = getMajorDiffTimeStamp3D(ORI_WALKERS[HARD_PD[HPDOMAIN_INDEX].walker_ID].diffusivity, HARD_PD[HPDOMAIN_INDEX].radius, ORI_WALKERS[HARD_PD[HPDOMAIN_INDEX].walker_ID].radius);
								// Event
								EVENT_INDEX = EVENT_INDEX + 1;
								EVENT_LIST[EVENT_INDEX].domainID = HARD_PD[HPDOMAIN_INDEX].ID;
								EVENT_LIST[EVENT_INDEX].domainType = 0;
								EVENT_LIST[EVENT_INDEX].timestamp = WORLD_CLOCK + HARD_PD[HPDOMAIN_INDEX].events_time;
								HARD_PD[HPDOMAIN_INDEX].event_card = EVENT_INDEX;
								cout << "1/1 Walker has been protected in HPD" << endl;
							}
						}
					}
					else {
						// General Pass: N > 1, different rebuid method
						// For all S0 walkers in original S0Index
						for (i = 0;i < S0_INDEX;i++) {
							// Avoid the walkers out of boundary
							if (ORI_WALKERS[S0_EVENT[i + 1]].statue == 0 && ORI_WALKERS[S0_EVENT[i + 1]].type != -1) {
								currentEvent::flag_s0_in_SPDa = 0;
								for (j = 0;j < SPDOMAINA_INDEX;j++) {
									if (SOFT_PD_A[j + 1].events_type != -1) {
										currentEvent::walker_distance = getSoftDomainASurfaceDistance(ORI_WALKERS[S0_EVENT[1]], SOFT_PD_A[j + 1]);
										if (currentEvent::walker_distance - THRES_2 < eps) {
											currentEvent::flag_s0_in_SPDa = 1;
											// Solver, new event
											// Update SPDa
											// Location adjust is not in this version
											SOFT_PD_A[j + 1].current_member = SOFT_PD_A[j + 1].current_member + 1;
											SOFT_PD_A[j + 1].member_list.member[SOFT_PD_A[j + 1].current_member] = ORI_WALKERS[S0_EVENT[i + 1]];
											currentEvent::update_spda_time_temp = SOFT_PD_A[j + 1].events_time;
											SOFT_PD_A[j + 1].events_time = SOFT_PD_A[j + 1].events_time; // Replaced with solver update interface once LocalSolver is provided: TODO
											currentEvent::update_spda_time_temp = SOFT_PD_A[j + 1].events_time - currentEvent::update_spda_time_temp;
											ORI_WALKERS[S0_EVENT[i + 1]].statue = 1;
											ORI_WALKERS[S0_EVENT[i + 1]].timestamp = WORLD_CLOCK;
											// No new event, just timestamp update
											EVENT_LIST[SOFT_PD_A[j + 1].event_card].timestamp = EVENT_LIST[SOFT_PD_A[j + 1].event_card].timestamp + currentEvent::update_spda_time_temp;
											cout << "1 Walker enters SPDa" << endl;
										}
									}
								}
							}
						}
						// S0 not in SPDa
						for (i = 0;i < S0_INDEX;i++) {
							// Avoid the walkers out of boundary
							if (ORI_WALKERS[S0_EVENT[i + 1]].statue == 0 && ORI_WALKERS[S0_EVENT[i + 1]].type != -1) {
								for (j = i + 1;j < S0_INDEX;j++) {
									// Avoid the walkers out of boundary
									if (ORI_WALKERS[S0_EVENT[j + 1]].statue == 0 && ORI_WALKERS[S0_EVENT[i + 1]].statue == 0 && ORI_WALKERS[S0_EVENT[j + 1]].type != -1) {
										// surface
										currentEvent::walker_distance = getDistance(ORI_WALKERS[S0_EVENT[i + 1]], ORI_WALKERS[S0_EVENT[j + 1]]) - ORI_WALKERS[S0_EVENT[i + 1]].radius - ORI_WALKERS[S0_EVENT[j + 1]].radius;
										if (currentEvent::walker_distance - THRES_2 < eps) {
											// New SPDa: Domain, Solver, new event
											cout << "Create a new Type 1 SPD during this event." << endl;
											SPDOMAINA_INDEX = SPDOMAINA_INDEX + 1;
											SOFT_PD_A[SPDOMAINA_INDEX].ID = SPDOMAINA_INDEX;
											SOFT_PD_A[SPDOMAINA_INDEX].center_distance = currentEvent::walker_distance;
											SOFT_PD_A[SPDOMAINA_INDEX].clock = WORLD_CLOCK;
											SOFT_PD_A[SPDOMAINA_INDEX].events_type = 10;
											SOFT_PD_A[SPDOMAINA_INDEX].current_member = 2;
											SOFT_PD_A[SPDOMAINA_INDEX].member_list.member[1] = ORI_WALKERS[S0_EVENT[i + 1]];
											SOFT_PD_A[SPDOMAINA_INDEX].member_list.member[2] = ORI_WALKERS[S0_EVENT[j + 1]];
											SOFT_PD_A[SPDOMAINA_INDEX].center1.x = ORI_WALKERS[S0_EVENT[i + 1]].center.x;
											SOFT_PD_A[SPDOMAINA_INDEX].center1.y = ORI_WALKERS[S0_EVENT[i + 1]].center.y;
											SOFT_PD_A[SPDOMAINA_INDEX].center1.z = ORI_WALKERS[S0_EVENT[i + 1]].center.z;
											SOFT_PD_A[SPDOMAINA_INDEX].center1.r = ORI_WALKERS[S0_EVENT[i + 1]].radius;
											SOFT_PD_A[SPDOMAINA_INDEX].center2.x = ORI_WALKERS[S0_EVENT[j + 1]].center.x;
											SOFT_PD_A[SPDOMAINA_INDEX].center2.y = ORI_WALKERS[S0_EVENT[j + 1]].center.y;
											SOFT_PD_A[SPDOMAINA_INDEX].center2.z = ORI_WALKERS[S0_EVENT[j + 1]].center.z;
											SOFT_PD_A[SPDOMAINA_INDEX].center2.r = ORI_WALKERS[S0_EVENT[j + 1]].radius;
											ORI_WALKERS[S0_EVENT[i + 1]].statue = 1;
											ORI_WALKERS[S0_EVENT[j + 1]].statue = 1;
											ORI_WALKERS[S0_EVENT[i + 1]].regionID = SPDOMAINA_INDEX;
											ORI_WALKERS[S0_EVENT[j + 1]].regionID = SPDOMAINA_INDEX;
											// Solver
											// RANDOM_SEED = rand() / (double)(RAND_MAX);
											// SOFT_PD_A[SPDOMAINA_INDEX].events_time = INITIAL_EVENTS_TIME + INITIAL_EVENTS_TIME*RANDOM_SEED; // Replaced with solver interface once LocalSolver is provided
											SOFT_PD_A[SPDOMAINA_INDEX].events_time = getReactionTimeStamp3D(SOFT_PD_A[SPDOMAINA_INDEX].member_list, INITIAL_EVENTS_TIME);
											// Event
											EVENT_INDEX = EVENT_INDEX + 1;
											EVENT_LIST[EVENT_INDEX].domainID = SOFT_PD_A[SPDOMAINA_INDEX].ID;
											EVENT_LIST[EVENT_INDEX].domainType = 1;
											EVENT_LIST[EVENT_INDEX].timestamp = WORLD_CLOCK + SOFT_PD_A[SPDOMAINA_INDEX].events_time;
											SOFT_PD_A[SPDOMAINA_INDEX].event_card = EVENT_INDEX;
											cout << "2 walkers built a new SPDa" << endl;
											cout << ORI_WALKERS[S0_EVENT[i + 1]].ID << " and " << ORI_WALKERS[S0_EVENT[j + 1]].ID << endl;
										}
									}
								}
							}
						}
						for (i = 0;i < S0_INDEX;i++) {
							// Avoid the walkers out of boundary
							if (ORI_WALKERS[S0_EVENT[i + 1]].statue == 0 && ORI_WALKERS[S0_EVENT[i + 1]].type != -1) {
								// Three routes
								// Close S0 walkers: O1
								// Not close but maximum and uncover case (dynamic) is not included in this version
								currentEvent::O1_allow_r = 0.0;
								currentEvent::O1_neighbout_ID = 0;
								currentEvent::O2_allow_r = 0.0;
								currentEvent::O3_allow_r = 0.0;
								currentEvent::min_distance = EXTRA_LARGE;
								for (j = i + 1;j < S0_INDEX;j++) {
									// Avoid the walkers out of boundary
									if (ORI_WALKERS[S0_EVENT[j + 1]].statue == 0 && ORI_WALKERS[S0_EVENT[j + 1]].type != -1 && ORI_WALKERS[S0_EVENT[j + 1]].type != -1) {
										currentEvent::walker_distance = getDistance(ORI_WALKERS[S0_EVENT[i + 1]], ORI_WALKERS[S0_EVENT[j + 1]]) - ORI_WALKERS[S0_EVENT[i + 1]].radius - ORI_WALKERS[S0_EVENT[j + 1]].radius;
										if (currentEvent::walker_distance - currentEvent::min_distance < eps) {
											currentEvent::min_distance = currentEvent::walker_distance;
											currentEvent::O1_neighbout_ID = S0_EVENT[j + 1];
											currentEvent::O1_allow_r = getDistance(ORI_WALKERS[S0_EVENT[i + 1]], ORI_WALKERS[S0_EVENT[j + 1]])*(ORI_WALKERS[S0_EVENT[i + 1]].diffusivity) / (ORI_WALKERS[S0_EVENT[i + 1]].diffusivity + ORI_WALKERS[S0_EVENT[j + 1]].diffusivity) - ORI_WALKERS[S0_EVENT[i + 1]].radius;
											//currentEvent::O1_neighour_allow_r = getDistance1D(ORI_WALKERS[S0_EVENT[i + 1]], ORI_WALKERS[S0_EVENT[j + 1]])*(ORI_WALKERS[S0_EVENT[j + 1]].diffusivity) / (ORI_WALKERS[S0_EVENT[i + 1]].diffusivity + ORI_WALKERS[S0_EVENT[j + 1]].diffusivity) - ORI_WALKERS[S0_EVENT[j + 1]].radius;
										}
									}
								}
								// HPD surface: O2
								currentEvent::min_distance = EXTRA_LARGE;
								for (j = 0;j < HPDOMAIN_INDEX;j++) {
									if (HARD_PD[j + 1].events_type != -1) {
										currentEvent::walker_distance = getHardDomainSurfaceDistance3D(ORI_WALKERS[S0_EVENT[i + 1]], HARD_PD[j + 1]);
										if (currentEvent::walker_distance - currentEvent::min_distance < eps) {
											currentEvent::min_distance = currentEvent::walker_distance;
											currentEvent::O2_allow_r = currentEvent::walker_distance; // Surface to surface = center to surface - r
										}
									}
								}
								// SPDa close center: 03
								currentEvent::min_distance = EXTRA_LARGE;
								for (j = 0;j < SPDOMAINA_INDEX;j++) {
									if (SOFT_PD_A[j + 1].events_type != -1) {
										currentEvent::walker_distance = getSoftDomainACentertoCenterDistance(ORI_WALKERS[S0_EVENT[i + 1]], SOFT_PD_A[j + 1]) - ORI_WALKERS[S0_EVENT[i + 1]].radius;
										if (currentEvent::walker_distance - currentEvent::min_distance < eps) {
											currentEvent::min_distance = currentEvent::walker_distance;
											currentEvent::O3_allow_r = currentEvent::walker_distance; // Walker surface to SPDa center A
										}
									}
								}
								// Use the minimum none 0 one to get 3 different routes of domains, etc
								if (currentEvent::O1_allow_r - currentEvent::O2_allow_r < eps && currentEvent::O1_allow_r > eps) {
									// O1 < O2 - > O1/O3
									if (currentEvent::O1_allow_r - currentEvent::O3_allow_r < eps) {
										// O1 < O3 - > O1
										// Domain
										if (DISPLAYFLAG_DEB == 1) {
											cout << currentEvent::O1_allow_r << endl;
										}
										cout << "Create a new HPD during this event." << endl;
										HPDOMAIN_INDEX = HPDOMAIN_INDEX + 1;
										HARD_PD[HPDOMAIN_INDEX].ID = HPDOMAIN_INDEX;
										HARD_PD[HPDOMAIN_INDEX].center.x = ORI_WALKERS[S0_EVENT[i + 1]].center.x;
										HARD_PD[HPDOMAIN_INDEX].center.y = ORI_WALKERS[S0_EVENT[i + 1]].center.y;
										HARD_PD[HPDOMAIN_INDEX].center.z = ORI_WALKERS[S0_EVENT[i + 1]].center.z;
										HARD_PD[HPDOMAIN_INDEX].clock = WORLD_CLOCK;
										HARD_PD[HPDOMAIN_INDEX].events_type = 1;
										HARD_PD[HPDOMAIN_INDEX].radius = currentEvent::O1_allow_r + ORI_WALKERS[S0_EVENT[i + 1]].radius;
										HARD_PD[HPDOMAIN_INDEX].walker_ID = ORI_WALKERS[S0_EVENT[i + 1]].ID;
										HARD_PD[HPDOMAIN_INDEX].walker_type = ORI_WALKERS[S0_EVENT[i + 1]].type;
										ORI_WALKERS[S0_EVENT[i + 1]].statue = 2;
										ORI_WALKERS[S0_EVENT[i + 1]].regionID = HPDOMAIN_INDEX;
										// Solver
										HARD_PD[HPDOMAIN_INDEX].events_time = getMajorDiffTimeStamp3D(ORI_WALKERS[HARD_PD[HPDOMAIN_INDEX].walker_ID].diffusivity, HARD_PD[HPDOMAIN_INDEX].radius, ORI_WALKERS[HARD_PD[HPDOMAIN_INDEX].walker_ID].radius);
										// Event
										EVENT_INDEX = EVENT_INDEX + 1;
										EVENT_LIST[EVENT_INDEX].domainID = HARD_PD[HPDOMAIN_INDEX].ID;
										EVENT_LIST[EVENT_INDEX].domainType = 0;
										EVENT_LIST[EVENT_INDEX].timestamp = WORLD_CLOCK + HARD_PD[HPDOMAIN_INDEX].events_time;
										HARD_PD[HPDOMAIN_INDEX].event_card = EVENT_INDEX;
										cout << "1 Walker has been protected in HPD" << endl;
									}
									else {
										// O3 < O1 - > O3 or O1 = 0
										// Domain
										if (DISPLAYFLAG_DEB == 1) {
											cout << currentEvent::O3_allow_r << endl;
										}
										cout << "Create a new HPD during this event." << endl;
										HPDOMAIN_INDEX = HPDOMAIN_INDEX + 1;
										HARD_PD[HPDOMAIN_INDEX].ID = HPDOMAIN_INDEX;
										HARD_PD[HPDOMAIN_INDEX].center.x = ORI_WALKERS[S0_EVENT[i + 1]].center.x;
										HARD_PD[HPDOMAIN_INDEX].center.y = ORI_WALKERS[S0_EVENT[i + 1]].center.y;
										HARD_PD[HPDOMAIN_INDEX].center.z = ORI_WALKERS[S0_EVENT[i + 1]].center.z;
										HARD_PD[HPDOMAIN_INDEX].clock = WORLD_CLOCK;
										HARD_PD[HPDOMAIN_INDEX].events_type = 1;
										HARD_PD[HPDOMAIN_INDEX].radius = currentEvent::O3_allow_r + ORI_WALKERS[S0_EVENT[i + 1]].radius;
										HARD_PD[HPDOMAIN_INDEX].walker_ID = ORI_WALKERS[S0_EVENT[i + 1]].ID;
										HARD_PD[HPDOMAIN_INDEX].walker_type = ORI_WALKERS[S0_EVENT[i + 1]].type;
										ORI_WALKERS[S0_EVENT[i + 1]].statue = 2;
										ORI_WALKERS[S0_EVENT[i + 1]].regionID = HPDOMAIN_INDEX;
										// Solver
										HARD_PD[HPDOMAIN_INDEX].events_time = getMajorDiffTimeStamp3D(ORI_WALKERS[HARD_PD[HPDOMAIN_INDEX].walker_ID].diffusivity, HARD_PD[HPDOMAIN_INDEX].radius, ORI_WALKERS[HARD_PD[HPDOMAIN_INDEX].walker_ID].radius);
										// Event
										EVENT_INDEX = EVENT_INDEX + 1;
										EVENT_LIST[EVENT_INDEX].domainID = HARD_PD[HPDOMAIN_INDEX].ID;
										EVENT_LIST[EVENT_INDEX].domainType = 0;
										EVENT_LIST[EVENT_INDEX].timestamp = WORLD_CLOCK + HARD_PD[HPDOMAIN_INDEX].events_time;
										HARD_PD[HPDOMAIN_INDEX].event_card = EVENT_INDEX;
										cout << "1 Walker has been protected in HPD" << endl;
									}
								}
								else {
									// O2 < O1 - > O2/O3
									if (currentEvent::O2_allow_r - currentEvent::O3_allow_r < eps) {
										// O2 < O3 - > O2
										// Domain
										if (DISPLAYFLAG_DEB == 1) {
											cout << currentEvent::O2_allow_r << endl;
										}
										cout << "Create a new HPD during this event." << endl;
										HPDOMAIN_INDEX = HPDOMAIN_INDEX + 1;
										HARD_PD[HPDOMAIN_INDEX].ID = HPDOMAIN_INDEX;
										HARD_PD[HPDOMAIN_INDEX].center.x = ORI_WALKERS[S0_EVENT[i + 1]].center.x;
										HARD_PD[HPDOMAIN_INDEX].center.y = ORI_WALKERS[S0_EVENT[i + 1]].center.y;
										HARD_PD[HPDOMAIN_INDEX].center.z = ORI_WALKERS[S0_EVENT[i + 1]].center.z;
										HARD_PD[HPDOMAIN_INDEX].clock = WORLD_CLOCK;
										HARD_PD[HPDOMAIN_INDEX].events_type = 1;
										HARD_PD[HPDOMAIN_INDEX].radius = currentEvent::O2_allow_r + ORI_WALKERS[S0_EVENT[i + 1]].radius;
										HARD_PD[HPDOMAIN_INDEX].walker_ID = ORI_WALKERS[S0_EVENT[i + 1]].ID;
										HARD_PD[HPDOMAIN_INDEX].walker_type = ORI_WALKERS[S0_EVENT[i + 1]].type;
										ORI_WALKERS[S0_EVENT[i + 1]].statue = 2;
										ORI_WALKERS[S0_EVENT[i + 1]].regionID = HPDOMAIN_INDEX;
										// Solver
										HARD_PD[HPDOMAIN_INDEX].events_time = getMajorDiffTimeStamp3D(ORI_WALKERS[HARD_PD[HPDOMAIN_INDEX].walker_ID].diffusivity, HARD_PD[HPDOMAIN_INDEX].radius, ORI_WALKERS[HARD_PD[HPDOMAIN_INDEX].walker_ID].radius);
										// Event
										EVENT_INDEX = EVENT_INDEX + 1;
										EVENT_LIST[EVENT_INDEX].domainID = HARD_PD[HPDOMAIN_INDEX].ID;
										EVENT_LIST[EVENT_INDEX].domainType = 0;
										EVENT_LIST[EVENT_INDEX].timestamp = WORLD_CLOCK + HARD_PD[HPDOMAIN_INDEX].events_time;
										HARD_PD[HPDOMAIN_INDEX].event_card = EVENT_INDEX;
										cout << "1 Walker has been protected in HPD" << endl;
									}
									else {
										// O3 < O2 - > O3
										// Domain
										if (DISPLAYFLAG_DEB == 1) {
											cout << currentEvent::O3_allow_r << endl;
										}
										cout << "Create a new HPD during this event." << endl;
										HPDOMAIN_INDEX = HPDOMAIN_INDEX + 1;
										HARD_PD[HPDOMAIN_INDEX].ID = HPDOMAIN_INDEX;
										HARD_PD[HPDOMAIN_INDEX].center.x = ORI_WALKERS[S0_EVENT[i + 1]].center.x;
										HARD_PD[HPDOMAIN_INDEX].center.y = ORI_WALKERS[S0_EVENT[i + 1]].center.y;
										HARD_PD[HPDOMAIN_INDEX].center.z = ORI_WALKERS[S0_EVENT[i + 1]].center.z;
										HARD_PD[HPDOMAIN_INDEX].clock = WORLD_CLOCK;
										HARD_PD[HPDOMAIN_INDEX].events_type = 1;
										HARD_PD[HPDOMAIN_INDEX].radius = currentEvent::O3_allow_r + ORI_WALKERS[S0_EVENT[i + 1]].radius;
										HARD_PD[HPDOMAIN_INDEX].walker_ID = ORI_WALKERS[S0_EVENT[i + 1]].ID;
										HARD_PD[HPDOMAIN_INDEX].walker_type = ORI_WALKERS[S0_EVENT[i + 1]].type;
										ORI_WALKERS[S0_EVENT[i + 1]].statue = 2;
										ORI_WALKERS[S0_EVENT[i + 1]].regionID = HPDOMAIN_INDEX;
										// Solver
										HARD_PD[HPDOMAIN_INDEX].events_time = getMajorDiffTimeStamp3D(ORI_WALKERS[HARD_PD[HPDOMAIN_INDEX].walker_ID].diffusivity, HARD_PD[HPDOMAIN_INDEX].radius, ORI_WALKERS[HARD_PD[HPDOMAIN_INDEX].walker_ID].radius);
										// Event
										EVENT_INDEX = EVENT_INDEX + 1;
										EVENT_LIST[EVENT_INDEX].domainID = HARD_PD[HPDOMAIN_INDEX].ID;
										EVENT_LIST[EVENT_INDEX].domainType = 0;
										EVENT_LIST[EVENT_INDEX].timestamp = WORLD_CLOCK + HARD_PD[HPDOMAIN_INDEX].events_time;
										HARD_PD[HPDOMAIN_INDEX].event_card = EVENT_INDEX;
										cout << "1 Walker has been protected in HPD" << endl;
									}
								}
							}
						}
					}// End of N > 1 else and all N if...else
					 // Create an extra new E20 card
					 // No need to decide where to insert, just get a card with type 20
					EVENT_INDEX = EVENT_INDEX + 1;
					EVENT_LIST[EVENT_INDEX].domainID = 0;
					EVENT_LIST[EVENT_INDEX].domainType = 3;
					EVENT_LIST[EVENT_INDEX].timestamp = WORLD_CLOCK + irradiationPulse;
					EVENT_LIST[EVENT_INDEX].eventType = 20;
				}
				else {
					// Locate domain
					if (currentEvent::domain_type == 0) {
						// HPD event
						currentEvent::event_type = HARD_PD[currentEvent::domainID].events_type;
						// HPD concern E01 and E11, only E01 for this version
						// May be cancelled by other event
						if (currentEvent::event_type == -1) {
							EVENT_POINTER = EVENT_POINTER + 1;
						}
						if (currentEvent::event_type == 1) {
							// Break HPD
							// Disable the domain
							HARD_PD[currentEvent::domainID].events_type = -1;
							// Call solver to get location, set walker to S0, send walker
							// getMajorDiffRelativeLocation3D(double Dv, double timestamp, double domainRadius, double walkerRadius, location centerLocation);
							currentEvent::afterward_location = getMajorDiffRelativeLocation3D(ORI_WALKERS[HARD_PD[currentEvent::domainID].walker_ID].diffusivity, HARD_PD[currentEvent::domainID].events_time, HARD_PD[currentEvent::domainID].radius, ORI_WALKERS[HARD_PD[currentEvent::domainID].walker_ID].radius, HARD_PD[currentEvent::domainID].center);
							// Move EVENT_POINTER to next event
							EVENT_POINTER = EVENT_POINTER + 1;
							// Temp E11 starts from here
							if (getDissociationFlag(ORI_WALKERS[HARD_PD[currentEvent::domainID].walker_ID])) {
								// Create one new walker
								WALKER_COUNTER = WALKER_COUNTER + 1;
								ORI_WALKERS[WALKER_COUNTER].ID = WALKER_COUNTER;
								ORI_WALKERS[WALKER_COUNTER].timestamp = WORLD_CLOCK;
								ORI_WALKERS[WALKER_COUNTER].type = ORI_WALKERS[HARD_PD[currentEvent::domainID].walker_ID].type;
								ORI_WALKERS[WALKER_COUNTER].component = 1;
								ORI_WALKERS[WALKER_COUNTER].statue = 0;
								ORI_WALKERS[WALKER_COUNTER].energy = 0.0; // Not in use for current version
								ORI_WALKERS[WALKER_COUNTER].radius = getEstimateRadius(ORI_WALKERS[WALKER_COUNTER]);
								ORI_WALKERS[WALKER_COUNTER].diffusivity = getDvValue(ORI_WALKERS[WALKER_COUNTER].energy, ORI_WALKERS[WALKER_COUNTER].type, ORI_WALKERS[WALKER_COUNTER].component);
								ORI_WALKERS[WALKER_COUNTER].center = ORI_WALKERS[HARD_PD[currentEvent::domainID].walker_ID].center; // Replace the last center
								S0_INDEX = 1;
								S0_EVENT[S0_INDEX] = ORI_WALKERS[WALKER_COUNTER].ID; // Log the IDs of all walkers with S0
								// Update old walker
								ORI_WALKERS[HARD_PD[currentEvent::domainID].walker_ID].statue = 0; // S0 release
								ORI_WALKERS[HARD_PD[currentEvent::domainID].walker_ID].center = currentEvent::afterward_location; // Move walker
								ORI_WALKERS[HARD_PD[currentEvent::domainID].walker_ID].component = ORI_WALKERS[HARD_PD[currentEvent::domainID].walker_ID].component - 1;
								ORI_WALKERS[HARD_PD[currentEvent::domainID].walker_ID].energy = ORI_WALKERS[HARD_PD[currentEvent::domainID].walker_ID].energy; // Not in use for current version
								ORI_WALKERS[HARD_PD[currentEvent::domainID].walker_ID].radius = getEstimateRadius(ORI_WALKERS[HARD_PD[currentEvent::domainID].walker_ID]);
								ORI_WALKERS[HARD_PD[currentEvent::domainID].walker_ID].diffusivity = getDvValue(ORI_WALKERS[HARD_PD[currentEvent::domainID].walker_ID].energy, ORI_WALKERS[HARD_PD[currentEvent::domainID].walker_ID].type, ORI_WALKERS[HARD_PD[currentEvent::domainID].walker_ID].component);
								S0_INDEX = S0_INDEX + 1;
								S0_EVENT[S0_INDEX] = HARD_PD[currentEvent::domainID].walker_ID; // Log the IDs of all walkers with S0
							}
							else {
								// Normal
								ORI_WALKERS[HARD_PD[currentEvent::domainID].walker_ID].statue = 0; // S0 release
								ORI_WALKERS[HARD_PD[currentEvent::domainID].walker_ID].center = currentEvent::afterward_location; // Move walker
								S0_INDEX = 1;
								S0_EVENT[1] = HARD_PD[currentEvent::domainID].walker_ID; // Log the IDs of all walkers with S0
							}
							// Temp E11 ends here
							// Enter ASMS 1st Stage
							S0_FULL_FLAG = 1;
							while (S0_FULL_FLAG) {
								// Search all exist S0 walkers
								for (i = 0;i < S0_INDEX;i++) {
									S0_FULL_FLAG = 0; // Default no need for next round, if any E02 is triggered, then set back to 1
													  // walker: ORI_WALKERS[S0_EVENT[i+1]]
									for (j = 0;j < HPDOMAIN_INDEX;j++) {
										if (HARD_PD[j + 1].events_type != -1) {
											// Get distance from current S0 walker to each HPD surface, check whether need E02 or not (TH1)
											if (getHardDomainSurfaceDistance3D(ORI_WALKERS[S0_EVENT[i + 1]], HARD_PD[j + 1]) - THRES_1 < eps) {
												// E02 is trigged for this HPD
												ORI_WALKERS[HARD_PD[j + 1].walker_ID].center = getTransDiffRelativeLocation3D(ORI_WALKERS[HARD_PD[j + 1].walker_ID].diffusivity, WORLD_CLOCK - ORI_WALKERS[HARD_PD[j + 1].walker_ID].timestamp, HARD_PD[j + 1].radius, ORI_WALKERS[HARD_PD[j + 1].walker_ID].radius, HARD_PD[j + 1].center);
												// Break, walker to S0, index + 1
												HARD_PD[j + 1].events_type = -1;
												ORI_WALKERS[HARD_PD[j + 1].walker_ID].statue = 0;
												S0_INDEX = S0_INDEX + 1;
												S0_EVENT[S0_INDEX] = HARD_PD[j + 1].walker_ID;
												S0_FULL_FLAG = 1;
											}
										}
										// Then check next HPD
									}
									// Finish check for the triggered S0, no HPD is still close, to next round: check new full S0 list
								}
							}
							cout << "Overall " << S0_INDEX << " Walkers Released in ASMS 1st Stage" << endl;
							for (i = 0;i < S0_INDEX;i++) {
								// X boundary
								if (ORI_WALKERS[S0_EVENT[i + 1]].statue == 0 && ORI_WALKERS[S0_EVENT[i + 1]].type != -1) {
									if (ORI_WALKERS[S0_EVENT[i + 1]].center.x + ORI_WALKERS[S0_EVENT[i + 1]].radius - X_UPPER_LIMIT>eps || ORI_WALKERS[S0_EVENT[i + 1]].center.x - ORI_WALKERS[S0_EVENT[i + 1]].radius - X_LOWWER_LIMIT < eps) {
										ORI_WALKERS[S0_EVENT[i + 1]].type = -1; // Inable this particle
										cout << "1 walker exits the X-axis boundary" << endl;
									}
								}
								// Y boundary
								if (ORI_WALKERS[S0_EVENT[i + 1]].statue == 0 && ORI_WALKERS[S0_EVENT[i + 1]].type != -1) {
									if (ORI_WALKERS[S0_EVENT[i + 1]].center.y + ORI_WALKERS[S0_EVENT[i + 1]].radius - Y_UPPER_LIMIT>eps || ORI_WALKERS[S0_EVENT[i + 1]].center.y - ORI_WALKERS[S0_EVENT[i + 1]].radius - Y_LOWWER_LIMIT < eps) {
										ORI_WALKERS[S0_EVENT[i + 1]].type = -1; // Inable this particle
										cout << "1 walker exits the Y-axis boundary" << endl;
									}
								}
								// Z boundary
								if (ORI_WALKERS[S0_EVENT[i + 1]].statue == 0 && ORI_WALKERS[S0_EVENT[i + 1]].type != -1) {
									if (ORI_WALKERS[S0_EVENT[i + 1]].center.z + ORI_WALKERS[S0_EVENT[i + 1]].radius - Z_UPPER_LIMIT>eps || ORI_WALKERS[S0_EVENT[i + 1]].center.z - ORI_WALKERS[S0_EVENT[i + 1]].radius - Z_LOWWER_LIMIT < eps) {
										ORI_WALKERS[S0_EVENT[i + 1]].type = -1; // Inable this particle
										cout << "1 walker exits the Z-axis boundary" << endl;
									}
								}
							}
							// Enter ASMS 2nd Stage
							// Fast pass: N = 1
							if (S0_INDEX == 1) {
								// Avoid the walkers out of boundary
								if (ORI_WALKERS[S0_EVENT[1]].type != -1) {
									// Whether this S0 walker located in SPDa
									currentEvent::flag_s0_in_SPDa = 0;
									for (i = 0;i < SPDOMAINA_INDEX;i++) {
										if (SOFT_PD_A[i + 1].events_type != -1) {
											currentEvent::walker_distance = getSoftDomainASurfaceDistance(ORI_WALKERS[S0_EVENT[1]], SOFT_PD_A[i + 1]);
											if (currentEvent::walker_distance - THRES_2 < eps) {
												currentEvent::flag_s0_in_SPDa = 1;
												// Solver, new event
												// Update SPDa
												// Location adjust is not in this version
												SOFT_PD_A[i + 1].current_member = SOFT_PD_A[i + 1].current_member + 1;
												SOFT_PD_A[i + 1].member_list.member[SOFT_PD_A[i + 1].current_member] = ORI_WALKERS[S0_EVENT[1]];
												currentEvent::update_spda_time_temp = SOFT_PD_A[i + 1].events_time;
												SOFT_PD_A[i + 1].events_time = SOFT_PD_A[i + 1].events_time; // Replaced with solver update interface once LocalSolver is provided: Not in this version
												currentEvent::update_spda_time_temp = SOFT_PD_A[i + 1].events_time - currentEvent::update_spda_time_temp;
												ORI_WALKERS[S0_EVENT[1]].statue = 1;
												ORI_WALKERS[S0_EVENT[1]].timestamp = WORLD_CLOCK;
												// No new event, just timestamp update
												EVENT_LIST[SOFT_PD_A[i + 1].event_card].timestamp = EVENT_LIST[SOFT_PD_A[i + 1].event_card].timestamp + currentEvent::update_spda_time_temp;
												cout << "1/1 Walker enters SPDa" << endl;
											}
										}
									}
									// S0 not in SPDa
									if (currentEvent::flag_s0_in_SPDa != 1) {
										currentEvent::min_distance = EXTRA_LARGE;
										currentEvent::near_domain_ID = 0;
										currentEvent::near_domain_type = 0; // Default on HPD, if SPD more close, modify to 1 during searching
																			// Check all enabled HPD
										for (i = 0;i < HPDOMAIN_INDEX;i++) {
											if (HARD_PD[i + 1].events_type != -1) {
												currentEvent::walker_distance = getHardDomainSurfaceDistance3D(ORI_WALKERS[S0_EVENT[1]], HARD_PD[i + 1]);
												if (currentEvent::walker_distance - currentEvent::min_distance < eps) {
													currentEvent::min_distance = currentEvent::walker_distance;
													currentEvent::near_domain_ID = i + 1;
												}
											}
										}
										// Check all enabled SPDa
										for (i = 0;i < SPDOMAINA_INDEX;i++) {
											if (SOFT_PD_A[i + 1].events_type != -1) {
												currentEvent::walker_distance = getSoftDomainASurfaceDistance(ORI_WALKERS[S0_EVENT[1]], SOFT_PD_A[i + 1]);
												if (currentEvent::walker_distance - currentEvent::min_distance < eps) {
													currentEvent::min_distance = currentEvent::walker_distance;
													currentEvent::near_domain_ID = i + 1;
													currentEvent::near_domain_type = 1;
												}
											}
										}
										// Domain, Solver, new event
										// near_domain_type is not a required way to set this part, reserve for other design
										if (DISPLAYFLAG_DEB == 1) {
											cout << currentEvent::min_distance << endl;
										}
										cout << "Create a new HPD during this event." << endl;
										HPDOMAIN_INDEX = HPDOMAIN_INDEX + 1;
										HARD_PD[HPDOMAIN_INDEX].ID = HPDOMAIN_INDEX;
										HARD_PD[HPDOMAIN_INDEX].center.x = ORI_WALKERS[S0_EVENT[1]].center.x;
										HARD_PD[HPDOMAIN_INDEX].center.y = ORI_WALKERS[S0_EVENT[1]].center.y;
										HARD_PD[HPDOMAIN_INDEX].center.z = ORI_WALKERS[S0_EVENT[1]].center.z;
										HARD_PD[HPDOMAIN_INDEX].clock = WORLD_CLOCK;
										HARD_PD[HPDOMAIN_INDEX].events_type = 1;
										HARD_PD[HPDOMAIN_INDEX].radius = currentEvent::min_distance + ORI_WALKERS[S0_EVENT[1]].radius;
										HARD_PD[HPDOMAIN_INDEX].walker_ID = ORI_WALKERS[S0_EVENT[1]].ID;
										HARD_PD[HPDOMAIN_INDEX].walker_type = ORI_WALKERS[S0_EVENT[1]].type;
										ORI_WALKERS[S0_EVENT[1]].statue = 2;
										ORI_WALKERS[S0_EVENT[1]].regionID = HPDOMAIN_INDEX;
										// Solver
										HARD_PD[HPDOMAIN_INDEX].events_time = getMajorDiffTimeStamp3D(ORI_WALKERS[HARD_PD[HPDOMAIN_INDEX].walker_ID].diffusivity, HARD_PD[HPDOMAIN_INDEX].radius, ORI_WALKERS[HARD_PD[HPDOMAIN_INDEX].walker_ID].radius);
										// Event
										EVENT_INDEX = EVENT_INDEX + 1;
										EVENT_LIST[EVENT_INDEX].domainID = HARD_PD[HPDOMAIN_INDEX].ID;
										EVENT_LIST[EVENT_INDEX].domainType = 0;
										EVENT_LIST[EVENT_INDEX].timestamp = WORLD_CLOCK + HARD_PD[HPDOMAIN_INDEX].events_time;
										HARD_PD[HPDOMAIN_INDEX].event_card = EVENT_INDEX;
										cout << "1/1 Walker has been protected in HPD" << endl;
									}
								}
							}
							else {
								// General Pass: N > 1, different rebuid method
								// For all S0 walkers in original S0Index
								for (i = 0;i < S0_INDEX;i++) {
									// Avoid the walkers out of boundary
									if (ORI_WALKERS[S0_EVENT[i + 1]].statue == 0 && ORI_WALKERS[S0_EVENT[i + 1]].type != -1) {
										currentEvent::flag_s0_in_SPDa = 0;
										for (j = 0;j < SPDOMAINA_INDEX;j++) {
											if (SOFT_PD_A[j + 1].events_type != -1) {
												currentEvent::walker_distance = getSoftDomainASurfaceDistance(ORI_WALKERS[S0_EVENT[1]], SOFT_PD_A[j + 1]);
												if (currentEvent::walker_distance - THRES_2 < eps) {
													currentEvent::flag_s0_in_SPDa = 1;
													// Solver, new event
													// Update SPDa
													// Location adjust is not in this version
													SOFT_PD_A[j + 1].current_member = SOFT_PD_A[j + 1].current_member + 1;
													SOFT_PD_A[j + 1].member_list.member[SOFT_PD_A[j + 1].current_member] = ORI_WALKERS[S0_EVENT[i + 1]];
													currentEvent::update_spda_time_temp = SOFT_PD_A[j + 1].events_time;
													SOFT_PD_A[j + 1].events_time = SOFT_PD_A[j + 1].events_time; // Replaced with solver update interface once LocalSolver is provided: TODO
													currentEvent::update_spda_time_temp = SOFT_PD_A[j + 1].events_time - currentEvent::update_spda_time_temp;
													ORI_WALKERS[S0_EVENT[i + 1]].statue = 1;
													ORI_WALKERS[S0_EVENT[i + 1]].timestamp = WORLD_CLOCK;
													// No new event, just timestamp update
													EVENT_LIST[SOFT_PD_A[j + 1].event_card].timestamp = EVENT_LIST[SOFT_PD_A[j + 1].event_card].timestamp + currentEvent::update_spda_time_temp;
													cout << "1 Walker enters SPDa" << endl;
												}
											}
										}
									}
								}
								// S0 not in SPDa
								for (i = 0;i < S0_INDEX;i++) {
									// Avoid the walkers out of boundary
									if (ORI_WALKERS[S0_EVENT[i + 1]].statue == 0 && ORI_WALKERS[S0_EVENT[i + 1]].type != -1) {
										for (j = i + 1;j < S0_INDEX;j++) {
											// Avoid the walkers out of boundary
											if (ORI_WALKERS[S0_EVENT[j + 1]].statue == 0 && ORI_WALKERS[S0_EVENT[i + 1]].statue == 0 && ORI_WALKERS[S0_EVENT[j + 1]].type != -1) {
												// surface
												currentEvent::walker_distance = getDistance(ORI_WALKERS[S0_EVENT[i + 1]], ORI_WALKERS[S0_EVENT[j + 1]]) - ORI_WALKERS[S0_EVENT[i + 1]].radius - ORI_WALKERS[S0_EVENT[j + 1]].radius;
												if (currentEvent::walker_distance - THRES_2 < eps) {
													// New SPDa: Domain, Solver, new event
													cout << "Create a new Type 1 SPD during this event." << endl;
													SPDOMAINA_INDEX = SPDOMAINA_INDEX + 1;
													SOFT_PD_A[SPDOMAINA_INDEX].ID = SPDOMAINA_INDEX;
													SOFT_PD_A[SPDOMAINA_INDEX].center_distance = currentEvent::walker_distance;
													SOFT_PD_A[SPDOMAINA_INDEX].clock = WORLD_CLOCK;
													SOFT_PD_A[SPDOMAINA_INDEX].events_type = 10;
													SOFT_PD_A[SPDOMAINA_INDEX].current_member = 2;
													SOFT_PD_A[SPDOMAINA_INDEX].member_list.member[1] = ORI_WALKERS[S0_EVENT[i + 1]];
													SOFT_PD_A[SPDOMAINA_INDEX].member_list.member[2] = ORI_WALKERS[S0_EVENT[j + 1]];
													SOFT_PD_A[SPDOMAINA_INDEX].center1.x = ORI_WALKERS[S0_EVENT[i + 1]].center.x;
													SOFT_PD_A[SPDOMAINA_INDEX].center1.y = ORI_WALKERS[S0_EVENT[i + 1]].center.y;
													SOFT_PD_A[SPDOMAINA_INDEX].center1.z = ORI_WALKERS[S0_EVENT[i + 1]].center.z;
													SOFT_PD_A[SPDOMAINA_INDEX].center1.r = ORI_WALKERS[S0_EVENT[i + 1]].radius;
													SOFT_PD_A[SPDOMAINA_INDEX].center2.x = ORI_WALKERS[S0_EVENT[j + 1]].center.x;
													SOFT_PD_A[SPDOMAINA_INDEX].center2.y = ORI_WALKERS[S0_EVENT[j + 1]].center.y;
													SOFT_PD_A[SPDOMAINA_INDEX].center2.z = ORI_WALKERS[S0_EVENT[j + 1]].center.z;
													SOFT_PD_A[SPDOMAINA_INDEX].center2.r = ORI_WALKERS[S0_EVENT[j + 1]].radius;
													ORI_WALKERS[S0_EVENT[i + 1]].statue = 1;
													ORI_WALKERS[S0_EVENT[j + 1]].statue = 1;
													ORI_WALKERS[S0_EVENT[i + 1]].regionID = SPDOMAINA_INDEX;
													ORI_WALKERS[S0_EVENT[j + 1]].regionID = SPDOMAINA_INDEX;
													// Solver
													// RANDOM_SEED = rand() / (double)(RAND_MAX);
													// SOFT_PD_A[SPDOMAINA_INDEX].events_time = INITIAL_EVENTS_TIME + INITIAL_EVENTS_TIME*RANDOM_SEED; // Replaced with solver interface once LocalSolver is provided
													SOFT_PD_A[SPDOMAINA_INDEX].events_time = getReactionTimeStamp3D(SOFT_PD_A[SPDOMAINA_INDEX].member_list, INITIAL_EVENTS_TIME);
													// Event
													EVENT_INDEX = EVENT_INDEX + 1;
													EVENT_LIST[EVENT_INDEX].domainID = SOFT_PD_A[SPDOMAINA_INDEX].ID;
													EVENT_LIST[EVENT_INDEX].domainType = 1;
													EVENT_LIST[EVENT_INDEX].timestamp = WORLD_CLOCK + SOFT_PD_A[SPDOMAINA_INDEX].events_time;
													SOFT_PD_A[SPDOMAINA_INDEX].event_card = EVENT_INDEX;
													cout << "2 walkers built a new SPDa" << endl;
													cout << ORI_WALKERS[S0_EVENT[i + 1]].ID << " and " << ORI_WALKERS[S0_EVENT[j + 1]].ID << endl;
												}
											}
										}
									}
								}
								for (i = 0;i < S0_INDEX;i++) {
									// Avoid the walkers out of boundary
									if (ORI_WALKERS[S0_EVENT[i + 1]].statue == 0 && ORI_WALKERS[S0_EVENT[i + 1]].type != -1) {
										// Three routes
										// Close S0 walkers: O1
										// Not close but maximum and uncover case (dynamic) is not included in this version
										currentEvent::O1_allow_r = 0.0;
										currentEvent::O1_neighbout_ID = 0;
										currentEvent::O2_allow_r = 0.0;
										currentEvent::O3_allow_r = 0.0;
										currentEvent::min_distance = EXTRA_LARGE;
										for (j = i + 1;j < S0_INDEX;j++) {
											// Avoid the walkers out of boundary
											if (ORI_WALKERS[S0_EVENT[j + 1]].statue == 0 && ORI_WALKERS[S0_EVENT[j + 1]].type != -1 && ORI_WALKERS[S0_EVENT[j + 1]].type != -1) {
												currentEvent::walker_distance = getDistance(ORI_WALKERS[S0_EVENT[i + 1]], ORI_WALKERS[S0_EVENT[j + 1]]) - ORI_WALKERS[S0_EVENT[i + 1]].radius - ORI_WALKERS[S0_EVENT[j + 1]].radius;
												if (currentEvent::walker_distance - currentEvent::min_distance < eps) {
													currentEvent::min_distance = currentEvent::walker_distance;
													currentEvent::O1_neighbout_ID = S0_EVENT[j + 1];
													currentEvent::O1_allow_r = getDistance(ORI_WALKERS[S0_EVENT[i + 1]], ORI_WALKERS[S0_EVENT[j + 1]])*(ORI_WALKERS[S0_EVENT[i + 1]].diffusivity) / (ORI_WALKERS[S0_EVENT[i + 1]].diffusivity + ORI_WALKERS[S0_EVENT[j + 1]].diffusivity) - ORI_WALKERS[S0_EVENT[i + 1]].radius;
													//currentEvent::O1_neighour_allow_r = getDistance1D(ORI_WALKERS[S0_EVENT[i + 1]], ORI_WALKERS[S0_EVENT[j + 1]])*(ORI_WALKERS[S0_EVENT[j + 1]].diffusivity) / (ORI_WALKERS[S0_EVENT[i + 1]].diffusivity + ORI_WALKERS[S0_EVENT[j + 1]].diffusivity) - ORI_WALKERS[S0_EVENT[j + 1]].radius;
												}
											}
										}
										// HPD surface: O2
										currentEvent::min_distance = EXTRA_LARGE;
										for (j = 0;j < HPDOMAIN_INDEX;j++) {
											if (HARD_PD[j + 1].events_type != -1) {
												currentEvent::walker_distance = getHardDomainSurfaceDistance3D(ORI_WALKERS[S0_EVENT[i + 1]], HARD_PD[j + 1]);
												if (currentEvent::walker_distance - currentEvent::min_distance < eps) {
													currentEvent::min_distance = currentEvent::walker_distance;
													currentEvent::O2_allow_r = currentEvent::walker_distance; // Surface to surface = center to surface - r
												}
											}
										}
										// SPDa close center: 03
										currentEvent::min_distance = EXTRA_LARGE;
										for (j = 0;j < SPDOMAINA_INDEX;j++) {
											if (SOFT_PD_A[j + 1].events_type != -1) {
												currentEvent::walker_distance = getSoftDomainACentertoCenterDistance(ORI_WALKERS[S0_EVENT[i + 1]], SOFT_PD_A[j + 1]) - ORI_WALKERS[S0_EVENT[i + 1]].radius;
												if (currentEvent::walker_distance - currentEvent::min_distance < eps) {
													currentEvent::min_distance = currentEvent::walker_distance;
													currentEvent::O3_allow_r = currentEvent::walker_distance; // Walker surface to SPDa center A
												}
											}
										}
										// Use the minimum none 0 one to get 3 different routes of domains, etc
										if (currentEvent::O1_allow_r - currentEvent::O2_allow_r < eps && currentEvent::O1_allow_r > eps) {
											// O1 < O2 - > O1/O3
											if (currentEvent::O1_allow_r - currentEvent::O3_allow_r < eps) {
												// O1 < O3 - > O1
												// Domain
												if (DISPLAYFLAG_DEB == 1) {
													cout << currentEvent::O1_allow_r << endl;
												}
												cout << "Create a new HPD during this event." << endl;
												HPDOMAIN_INDEX = HPDOMAIN_INDEX + 1;
												HARD_PD[HPDOMAIN_INDEX].ID = HPDOMAIN_INDEX;
												HARD_PD[HPDOMAIN_INDEX].center.x = ORI_WALKERS[S0_EVENT[i + 1]].center.x;
												HARD_PD[HPDOMAIN_INDEX].center.y = ORI_WALKERS[S0_EVENT[i + 1]].center.y;
												HARD_PD[HPDOMAIN_INDEX].center.z = ORI_WALKERS[S0_EVENT[i + 1]].center.z;
												HARD_PD[HPDOMAIN_INDEX].clock = WORLD_CLOCK;
												HARD_PD[HPDOMAIN_INDEX].events_type = 1;
												HARD_PD[HPDOMAIN_INDEX].radius = currentEvent::O1_allow_r + ORI_WALKERS[S0_EVENT[i + 1]].radius;
												HARD_PD[HPDOMAIN_INDEX].walker_ID = ORI_WALKERS[S0_EVENT[i + 1]].ID;
												HARD_PD[HPDOMAIN_INDEX].walker_type = ORI_WALKERS[S0_EVENT[i + 1]].type;
												ORI_WALKERS[S0_EVENT[i + 1]].statue = 2;
												ORI_WALKERS[S0_EVENT[i + 1]].regionID = HPDOMAIN_INDEX;
												// Solver
												HARD_PD[HPDOMAIN_INDEX].events_time = getMajorDiffTimeStamp3D(ORI_WALKERS[HARD_PD[HPDOMAIN_INDEX].walker_ID].diffusivity, HARD_PD[HPDOMAIN_INDEX].radius, ORI_WALKERS[HARD_PD[HPDOMAIN_INDEX].walker_ID].radius);
												// Event
												EVENT_INDEX = EVENT_INDEX + 1;
												EVENT_LIST[EVENT_INDEX].domainID = HARD_PD[HPDOMAIN_INDEX].ID;
												EVENT_LIST[EVENT_INDEX].domainType = 0;
												EVENT_LIST[EVENT_INDEX].timestamp = WORLD_CLOCK + HARD_PD[HPDOMAIN_INDEX].events_time;
												HARD_PD[HPDOMAIN_INDEX].event_card = EVENT_INDEX;
												cout << "1 Walker has been protected in HPD" << endl;
											}
											else {
												// O3 < O1 - > O3 or O1 = 0
												// Domain
												if (DISPLAYFLAG_DEB == 1) {
													cout << currentEvent::O3_allow_r << endl;
												}
												cout << "Create a new HPD during this event." << endl;
												HPDOMAIN_INDEX = HPDOMAIN_INDEX + 1;
												HARD_PD[HPDOMAIN_INDEX].ID = HPDOMAIN_INDEX;
												HARD_PD[HPDOMAIN_INDEX].center.x = ORI_WALKERS[S0_EVENT[i + 1]].center.x;
												HARD_PD[HPDOMAIN_INDEX].center.y = ORI_WALKERS[S0_EVENT[i + 1]].center.y;
												HARD_PD[HPDOMAIN_INDEX].center.z = ORI_WALKERS[S0_EVENT[i + 1]].center.z;
												HARD_PD[HPDOMAIN_INDEX].clock = WORLD_CLOCK;
												HARD_PD[HPDOMAIN_INDEX].events_type = 1;
												HARD_PD[HPDOMAIN_INDEX].radius = currentEvent::O3_allow_r + ORI_WALKERS[S0_EVENT[i + 1]].radius;
												HARD_PD[HPDOMAIN_INDEX].walker_ID = ORI_WALKERS[S0_EVENT[i + 1]].ID;
												HARD_PD[HPDOMAIN_INDEX].walker_type = ORI_WALKERS[S0_EVENT[i + 1]].type;
												ORI_WALKERS[S0_EVENT[i + 1]].statue = 2;
												ORI_WALKERS[S0_EVENT[i + 1]].regionID = HPDOMAIN_INDEX;
												// Solver
												HARD_PD[HPDOMAIN_INDEX].events_time = getMajorDiffTimeStamp3D(ORI_WALKERS[HARD_PD[HPDOMAIN_INDEX].walker_ID].diffusivity, HARD_PD[HPDOMAIN_INDEX].radius, ORI_WALKERS[HARD_PD[HPDOMAIN_INDEX].walker_ID].radius);
												// Event
												EVENT_INDEX = EVENT_INDEX + 1;
												EVENT_LIST[EVENT_INDEX].domainID = HARD_PD[HPDOMAIN_INDEX].ID;
												EVENT_LIST[EVENT_INDEX].domainType = 0;
												EVENT_LIST[EVENT_INDEX].timestamp = WORLD_CLOCK + HARD_PD[HPDOMAIN_INDEX].events_time;
												HARD_PD[HPDOMAIN_INDEX].event_card = EVENT_INDEX;
												cout << "1 Walker has been protected in HPD" << endl;
											}
										}
										else {
											// O2 < O1 - > O2/O3
											if (currentEvent::O2_allow_r - currentEvent::O3_allow_r < eps) {
												// O2 < O3 - > O2
												// Domain
												if (DISPLAYFLAG_DEB == 1) {
													cout << currentEvent::O2_allow_r << endl;
												}
												cout << "Create a new HPD during this event." << endl;
												HPDOMAIN_INDEX = HPDOMAIN_INDEX + 1;
												HARD_PD[HPDOMAIN_INDEX].ID = HPDOMAIN_INDEX;
												HARD_PD[HPDOMAIN_INDEX].center.x = ORI_WALKERS[S0_EVENT[i + 1]].center.x;
												HARD_PD[HPDOMAIN_INDEX].center.y = ORI_WALKERS[S0_EVENT[i + 1]].center.y;
												HARD_PD[HPDOMAIN_INDEX].center.z = ORI_WALKERS[S0_EVENT[i + 1]].center.z;
												HARD_PD[HPDOMAIN_INDEX].clock = WORLD_CLOCK;
												HARD_PD[HPDOMAIN_INDEX].events_type = 1;
												HARD_PD[HPDOMAIN_INDEX].radius = currentEvent::O2_allow_r + ORI_WALKERS[S0_EVENT[i + 1]].radius;
												HARD_PD[HPDOMAIN_INDEX].walker_ID = ORI_WALKERS[S0_EVENT[i + 1]].ID;
												HARD_PD[HPDOMAIN_INDEX].walker_type = ORI_WALKERS[S0_EVENT[i + 1]].type;
												ORI_WALKERS[S0_EVENT[i + 1]].statue = 2;
												ORI_WALKERS[S0_EVENT[i + 1]].regionID = HPDOMAIN_INDEX;
												// Solver
												HARD_PD[HPDOMAIN_INDEX].events_time = getMajorDiffTimeStamp3D(ORI_WALKERS[HARD_PD[HPDOMAIN_INDEX].walker_ID].diffusivity, HARD_PD[HPDOMAIN_INDEX].radius, ORI_WALKERS[HARD_PD[HPDOMAIN_INDEX].walker_ID].radius);
												// Event
												EVENT_INDEX = EVENT_INDEX + 1;
												EVENT_LIST[EVENT_INDEX].domainID = HARD_PD[HPDOMAIN_INDEX].ID;
												EVENT_LIST[EVENT_INDEX].domainType = 0;
												EVENT_LIST[EVENT_INDEX].timestamp = WORLD_CLOCK + HARD_PD[HPDOMAIN_INDEX].events_time;
												HARD_PD[HPDOMAIN_INDEX].event_card = EVENT_INDEX;
												cout << "1 Walker has been protected in HPD" << endl;
											}
											else {
												// O3 < O2 - > O3
												// Domain
												if (DISPLAYFLAG_DEB == 1) {
													cout << currentEvent::O3_allow_r << endl;
												}
												cout << "Create a new HPD during this event." << endl;
												HPDOMAIN_INDEX = HPDOMAIN_INDEX + 1;
												HARD_PD[HPDOMAIN_INDEX].ID = HPDOMAIN_INDEX;
												HARD_PD[HPDOMAIN_INDEX].center.x = ORI_WALKERS[S0_EVENT[i + 1]].center.x;
												HARD_PD[HPDOMAIN_INDEX].center.y = ORI_WALKERS[S0_EVENT[i + 1]].center.y;
												HARD_PD[HPDOMAIN_INDEX].center.z = ORI_WALKERS[S0_EVENT[i + 1]].center.z;
												HARD_PD[HPDOMAIN_INDEX].clock = WORLD_CLOCK;
												HARD_PD[HPDOMAIN_INDEX].events_type = 1;
												HARD_PD[HPDOMAIN_INDEX].radius = currentEvent::O3_allow_r + ORI_WALKERS[S0_EVENT[i + 1]].radius;
												HARD_PD[HPDOMAIN_INDEX].walker_ID = ORI_WALKERS[S0_EVENT[i + 1]].ID;
												HARD_PD[HPDOMAIN_INDEX].walker_type = ORI_WALKERS[S0_EVENT[i + 1]].type;
												ORI_WALKERS[S0_EVENT[i + 1]].statue = 2;
												ORI_WALKERS[S0_EVENT[i + 1]].regionID = HPDOMAIN_INDEX;
												// Solver
												HARD_PD[HPDOMAIN_INDEX].events_time = getMajorDiffTimeStamp3D(ORI_WALKERS[HARD_PD[HPDOMAIN_INDEX].walker_ID].diffusivity, HARD_PD[HPDOMAIN_INDEX].radius, ORI_WALKERS[HARD_PD[HPDOMAIN_INDEX].walker_ID].radius);
												// Event
												EVENT_INDEX = EVENT_INDEX + 1;
												EVENT_LIST[EVENT_INDEX].domainID = HARD_PD[HPDOMAIN_INDEX].ID;
												EVENT_LIST[EVENT_INDEX].domainType = 0;
												EVENT_LIST[EVENT_INDEX].timestamp = WORLD_CLOCK + HARD_PD[HPDOMAIN_INDEX].events_time;
												HARD_PD[HPDOMAIN_INDEX].event_card = EVENT_INDEX;
												cout << "1 Walker has been protected in HPD" << endl;
											}
										}
									}
								}
							} // End of N > 1 else and all N if...else
						}
						if (currentEvent::event_type == 11) {
							// Move EVENT_POINTER to next event
							EVENT_POINTER = EVENT_POINTER + 1;
							cout << "FATAL ERROR 700: Unable to process E11 for HPD in this version. Exit." << endl;
							system("pause");
							exit(0);
						}
					}
					if (currentEvent::domain_type == 1) {
						// SPDa event
						currentEvent::event_type = SOFT_PD_A[currentEvent::domainID].events_type;
						// SPD concern E10 only, this is not in input of update
						// Coalescence enabled
						// Move EVENT_POINTER to next event
						EVENT_POINTER = EVENT_POINTER + 1;
						// Break SPDa
						// Disable the domain
						SOFT_PD_A[currentEvent::domainID].events_type = -1;
						// Call solver to get location, send walker, set walker to S0					
						// Disable all current walker (type to -1)
						for (i = 0;i < SOFT_PD_A[currentEvent::domainID].current_member;i++) {
							ORI_WALKERS[SOFT_PD_A[currentEvent::domainID].member_list.member[i + 1].ID].statue = 0;
							ORI_WALKERS[SOFT_PD_A[currentEvent::domainID].member_list.member[i + 1].ID].type = -1;
						}
						// Call Solver to find the find statue										
						// getReactionFinalStatues3D(reaction_current_walker walkers, ini init_member)
						SOFT_PD_A[currentEvent::domainID].e10_final_statue = getReactionFinalStatues3D(SOFT_PD_A[currentEvent::domainID].member_list, SOFT_PD_A[currentEvent::domainID].current_member, SOFT_PD_A[currentEvent::domainID].center1, SOFT_PD_A[currentEvent::domainID].center2);
						// Find how many final walkers from the solver, max with 5 in current version					
						SOFT_PD_A[currentEvent::domainID].final_member = 0;
						for (i = 1;i < 6;i++) {
							if (SOFT_PD_A[currentEvent::domainID].e10_final_statue.member[i].type != -1) {
								SOFT_PD_A[currentEvent::domainID].final_member = SOFT_PD_A[currentEvent::domainID].final_member + 1;
							}
						}
						cout << "Encount with " << SOFT_PD_A[currentEvent::domainID].final_member << " new walker(s) in this event." << endl;
						// Create new walkers, all in S0
						// Update WALKER_COUNTER and S0 log
						S0_INDEX = 0;
						for (i = 0;i < SOFT_PD_A[currentEvent::domainID].final_member;i++) {
							// Update walker counter, point to new walker
							// type, component, center_location, radius(func), ID, statue(S0), energy(reserve), diffusivity(func), timestamp(WORLD_CLOCK)
							if (DISPLAYFLAG_SYS) {
								cout << "Create 1 new walker in this process." << endl;
							}
							WALKER_COUNTER = WALKER_COUNTER + 1;
							ORI_WALKERS[WALKER_COUNTER].ID = WALKER_COUNTER;
							ORI_WALKERS[WALKER_COUNTER].timestamp = WORLD_CLOCK;
							ORI_WALKERS[WALKER_COUNTER].type = SOFT_PD_A[currentEvent::domainID].e10_final_statue.member[i + 1].type;
							ORI_WALKERS[WALKER_COUNTER].component = SOFT_PD_A[currentEvent::domainID].e10_final_statue.member[i + 1].component;
							ORI_WALKERS[WALKER_COUNTER].statue = 0;
							ORI_WALKERS[WALKER_COUNTER].energy = 0.0; // Not in use for current version
							ORI_WALKERS[WALKER_COUNTER].center = SOFT_PD_A[currentEvent::domainID].e10_final_statue.member[i + 1].center;
							ORI_WALKERS[WALKER_COUNTER].radius = getEstimateRadius(ORI_WALKERS[WALKER_COUNTER]);
							ORI_WALKERS[WALKER_COUNTER].diffusivity = getDvValue(ORI_WALKERS[WALKER_COUNTER].energy, ORI_WALKERS[WALKER_COUNTER].type, ORI_WALKERS[WALKER_COUNTER].component);
							// Debug
							cout << "1 type " << ORI_WALKERS[WALKER_COUNTER].type << " walker with " << ORI_WALKERS[WALKER_COUNTER].component << " members created." << endl;
							// Update S0 log
							S0_INDEX = S0_INDEX + 1;
							S0_EVENT[S0_INDEX] = ORI_WALKERS[WALKER_COUNTER].ID; // Log the IDs of all walkers with S0
						}
						if (S0_INDEX == 0) {
							// No new walker, bypass ASMS
							cout << "Annihilation in SPDa." << endl;
						}
						else {
							// Begin ASMS, full check around
							// Don't touch EVENT_POINTER anymore
							// Enter ASMS 1st Stage
							S0_FULL_FLAG = 1;
							while (S0_FULL_FLAG) {
								// Search all exist S0 walkers
								for (i = 0;i < S0_INDEX;i++) {
									S0_FULL_FLAG = 0; // Default no need for next round, if any E02 is triggered, then set back to 1
													  // walker: ORI_WALKERS[S0_EVENT[i+1]]
									for (j = 0;j < HPDOMAIN_INDEX;j++) {
										if (HARD_PD[j + 1].events_type != -1) {
											// Get distance from current S0 walker to each HPD surface, check whether need E02 or not (TH1)
											if (getHardDomainSurfaceDistance3D(ORI_WALKERS[S0_EVENT[i + 1]], HARD_PD[j + 1]) - THRES_1 < eps) {
												// E02 is trigged for this HPD
												ORI_WALKERS[HARD_PD[j + 1].walker_ID].center = getTransDiffRelativeLocation3D(ORI_WALKERS[HARD_PD[j + 1].walker_ID].diffusivity, WORLD_CLOCK - ORI_WALKERS[HARD_PD[j + 1].walker_ID].timestamp, HARD_PD[j + 1].radius, ORI_WALKERS[HARD_PD[j + 1].walker_ID].radius, HARD_PD[j + 1].center);
												// Break, walker to S0, index + 1
												HARD_PD[j + 1].events_type = -1;
												ORI_WALKERS[HARD_PD[j + 1].walker_ID].statue = 0;
												S0_INDEX = S0_INDEX + 1;
												S0_EVENT[S0_INDEX] = HARD_PD[j + 1].walker_ID;
												S0_FULL_FLAG = 1;
											}
										}
										// Then check next HPD
									}
									// Finish check for the triggered S0, no HPD is still close, to next round: check new full S0 list
								}
							}
							cout << "Overall " << S0_INDEX << " Walkers Released in ASMS 1st Stage" << endl;
							for (i = 0;i < S0_INDEX;i++) {
								// X boundary
								if (ORI_WALKERS[S0_EVENT[i + 1]].statue == 0 && ORI_WALKERS[S0_EVENT[i + 1]].type != -1) {
									if (ORI_WALKERS[S0_EVENT[i + 1]].center.x + ORI_WALKERS[S0_EVENT[i + 1]].radius - X_UPPER_LIMIT>eps || ORI_WALKERS[S0_EVENT[i + 1]].center.x - ORI_WALKERS[S0_EVENT[i + 1]].radius - X_LOWWER_LIMIT < eps) {
										ORI_WALKERS[S0_EVENT[i + 1]].type = -1; // Inable this particle
										cout << "1 walker exits the X-axis boundary" << endl;
									}
								}
								// Y boundary
								if (ORI_WALKERS[S0_EVENT[i + 1]].statue == 0 && ORI_WALKERS[S0_EVENT[i + 1]].type != -1) {
									if (ORI_WALKERS[S0_EVENT[i + 1]].center.y + ORI_WALKERS[S0_EVENT[i + 1]].radius - Y_UPPER_LIMIT>eps || ORI_WALKERS[S0_EVENT[i + 1]].center.y - ORI_WALKERS[S0_EVENT[i + 1]].radius - Y_LOWWER_LIMIT < eps) {
										ORI_WALKERS[S0_EVENT[i + 1]].type = -1; // Inable this particle
										cout << "1 walker exits the Y-axis boundary" << endl;
									}
								}
								// Z boundary
								if (ORI_WALKERS[S0_EVENT[i + 1]].statue == 0 && ORI_WALKERS[S0_EVENT[i + 1]].type != -1) {
									if (ORI_WALKERS[S0_EVENT[i + 1]].center.z + ORI_WALKERS[S0_EVENT[i + 1]].radius - Z_UPPER_LIMIT>eps || ORI_WALKERS[S0_EVENT[i + 1]].center.z - ORI_WALKERS[S0_EVENT[i + 1]].radius - Z_LOWWER_LIMIT < eps) {
										ORI_WALKERS[S0_EVENT[i + 1]].type = -1; // Inable this particle
										cout << "1 walker exits the Z-axis boundary" << endl;
									}
								}
							}
							// Enter ASMS 2nd Stage
							// Protect walker again, build new domain or update SPDa
							// Fast pass: N = 1
							if (S0_INDEX == 1) {
								// Avoid the walkers out of boundary
								if (ORI_WALKERS[S0_EVENT[1]].type != -1) {
									// Whether this S0 walker located in SPDa
									currentEvent::flag_s0_in_SPDa = 0;
									for (i = 0;i < SPDOMAINA_INDEX;i++) {
										// Avoid the SPDa has been disabled at the beginning of this event
										if (SOFT_PD_A[i + 1].events_type != -1) {
											currentEvent::walker_distance = getSoftDomainASurfaceDistance(ORI_WALKERS[S0_EVENT[1]], SOFT_PD_A[i + 1]);
											if (currentEvent::walker_distance - THRES_2 < eps) {
												currentEvent::flag_s0_in_SPDa = 1;
												// Solver, new event
												// Update SPDa
												// Location adjust is not in this version
												SOFT_PD_A[i + 1].current_member = SOFT_PD_A[i + 1].current_member + 1;
												SOFT_PD_A[i + 1].member_list.member[SOFT_PD_A[i + 1].current_member] = ORI_WALKERS[S0_EVENT[1]];
												currentEvent::update_spda_time_temp = SOFT_PD_A[i + 1].events_time;
												SOFT_PD_A[i + 1].events_time = SOFT_PD_A[i + 1].events_time; // Replaced with solver update interface once LocalSolver is provided: Not in this version
												currentEvent::update_spda_time_temp = SOFT_PD_A[i + 1].events_time - currentEvent::update_spda_time_temp;
												ORI_WALKERS[S0_EVENT[1]].statue = 1;
												ORI_WALKERS[S0_EVENT[1]].timestamp = WORLD_CLOCK;
												// No new event, just timestamp update
												EVENT_LIST[SOFT_PD_A[i + 1].event_card].timestamp = EVENT_LIST[SOFT_PD_A[i + 1].event_card].timestamp + currentEvent::update_spda_time_temp;
												cout << "1/1 Walker enters SPDa" << endl;
											}
										}
									}
									// S0 not in SPDa
									if (currentEvent::flag_s0_in_SPDa != 1) {
										currentEvent::min_distance = EXTRA_LARGE;
										currentEvent::near_domain_ID = 0;
										currentEvent::near_domain_type = 0; // Default on HPD, if SPD more close, modify to 1 during searching
																			// Check all enabled HPD
										for (i = 0;i < HPDOMAIN_INDEX;i++) {
											if (HARD_PD[i + 1].events_type != -1) {
												currentEvent::walker_distance = getHardDomainSurfaceDistance3D(ORI_WALKERS[S0_EVENT[1]], HARD_PD[i + 1]);
												if (currentEvent::walker_distance - currentEvent::min_distance < eps) {
													currentEvent::min_distance = currentEvent::walker_distance;
													currentEvent::near_domain_ID = i + 1;
												}
											}
										}
										// Check all enabled SPDa
										for (i = 0;i < SPDOMAINA_INDEX;i++) {
											if (SOFT_PD_A[i + 1].events_type != -1) {
												currentEvent::walker_distance = getSoftDomainASurfaceDistance(ORI_WALKERS[S0_EVENT[1]], SOFT_PD_A[i + 1]);
												if (currentEvent::walker_distance - currentEvent::min_distance < eps) {
													currentEvent::min_distance = currentEvent::walker_distance;
													currentEvent::near_domain_ID = i + 1;
													currentEvent::near_domain_type = 1;
												}
											}
										}
										// Domain, Solver, new event
										// near_domain_type is not a required way to set this part, reserve for other design
										cout << "Create a new HPD during this event." << endl;
										HPDOMAIN_INDEX = HPDOMAIN_INDEX + 1;
										HARD_PD[HPDOMAIN_INDEX].ID = HPDOMAIN_INDEX;
										HARD_PD[HPDOMAIN_INDEX].center.x = ORI_WALKERS[S0_EVENT[1]].center.x;
										HARD_PD[HPDOMAIN_INDEX].center.y = ORI_WALKERS[S0_EVENT[1]].center.y;
										HARD_PD[HPDOMAIN_INDEX].center.z = ORI_WALKERS[S0_EVENT[1]].center.z;
										HARD_PD[HPDOMAIN_INDEX].clock = WORLD_CLOCK;
										HARD_PD[HPDOMAIN_INDEX].events_type = 1;
										HARD_PD[HPDOMAIN_INDEX].radius = currentEvent::min_distance + ORI_WALKERS[S0_EVENT[1]].radius;
										HARD_PD[HPDOMAIN_INDEX].walker_ID = ORI_WALKERS[S0_EVENT[1]].ID;
										HARD_PD[HPDOMAIN_INDEX].walker_type = ORI_WALKERS[S0_EVENT[1]].type;
										ORI_WALKERS[S0_EVENT[1]].statue = 2;
										ORI_WALKERS[S0_EVENT[1]].regionID = HPDOMAIN_INDEX;
										// Solver
										HARD_PD[HPDOMAIN_INDEX].events_time = getMajorDiffTimeStamp3D(ORI_WALKERS[HARD_PD[HPDOMAIN_INDEX].walker_ID].diffusivity, HARD_PD[HPDOMAIN_INDEX].radius, ORI_WALKERS[HARD_PD[HPDOMAIN_INDEX].walker_ID].radius);
										// Event
										EVENT_INDEX = EVENT_INDEX + 1;
										EVENT_LIST[EVENT_INDEX].domainID = HARD_PD[HPDOMAIN_INDEX].ID;
										EVENT_LIST[EVENT_INDEX].domainType = 0;
										EVENT_LIST[EVENT_INDEX].timestamp = WORLD_CLOCK + HARD_PD[HPDOMAIN_INDEX].events_time;
										HARD_PD[HPDOMAIN_INDEX].event_card = EVENT_INDEX;
										cout << "1/1 Walker has been protected in HPD" << endl;
									}
								}
							}
							else {
								// General Pass: N > 1, different rebuid method
								// For all S0 walkers in original S0Index
								for (i = 0;i < S0_INDEX;i++) {
									// Avoid the walkers out of boundary
									if (ORI_WALKERS[S0_EVENT[i + 1]].statue == 0 && ORI_WALKERS[S0_EVENT[i + 1]].type != -1) {
										currentEvent::flag_s0_in_SPDa = 0;
										for (j = 0;j < SPDOMAINA_INDEX;j++) {
											if (SOFT_PD_A[j + 1].events_type != -1) {
												currentEvent::walker_distance = getSoftDomainASurfaceDistance(ORI_WALKERS[S0_EVENT[1]], SOFT_PD_A[j + 1]);
												if (currentEvent::walker_distance - THRES_2 < eps) {
													currentEvent::flag_s0_in_SPDa = 1;
													// Solver, new event
													// Update SPDa
													// Location adjust is not in this version
													SOFT_PD_A[j + 1].current_member = SOFT_PD_A[j + 1].current_member + 1;
													SOFT_PD_A[j + 1].member_list.member[SOFT_PD_A[j + 1].current_member] = ORI_WALKERS[S0_EVENT[i + 1]];
													currentEvent::update_spda_time_temp = SOFT_PD_A[j + 1].events_time;
													SOFT_PD_A[j + 1].events_time = SOFT_PD_A[j + 1].events_time; // Replaced with solver update interface once LocalSolver is provided: TODO
													currentEvent::update_spda_time_temp = SOFT_PD_A[j + 1].events_time - currentEvent::update_spda_time_temp;
													ORI_WALKERS[S0_EVENT[i + 1]].statue = 1;
													ORI_WALKERS[S0_EVENT[i + 1]].timestamp = WORLD_CLOCK;
													// No new event, just timestamp update
													EVENT_LIST[SOFT_PD_A[j + 1].event_card].timestamp = EVENT_LIST[SOFT_PD_A[j + 1].event_card].timestamp + currentEvent::update_spda_time_temp;
													cout << "1 Walker enters SPDa" << endl;
												}
											}
										}
									}
								}
								// S0 not in SPDa
								for (i = 0;i < S0_INDEX;i++) {
									// Avoid the walkers out of boundary
									if (ORI_WALKERS[S0_EVENT[i + 1]].statue == 0 && ORI_WALKERS[S0_EVENT[i + 1]].type != -1) {
										for (j = i + 1;j < S0_INDEX;j++) {
											// Avoid the walkers out of boundary
											if (ORI_WALKERS[S0_EVENT[j + 1]].statue == 0 && ORI_WALKERS[S0_EVENT[i + 1]].statue == 0 && ORI_WALKERS[S0_EVENT[j + 1]].type != -1) {
												// surface
												currentEvent::walker_distance = getDistance(ORI_WALKERS[S0_EVENT[i + 1]], ORI_WALKERS[S0_EVENT[j + 1]]) - ORI_WALKERS[S0_EVENT[i + 1]].radius - ORI_WALKERS[S0_EVENT[j + 1]].radius;
												if (currentEvent::walker_distance - THRES_2 < eps) {
													// New SPDa: Domain, Solver, new event
													cout << "Create a new Type 1 SPD during this event." << endl;
													SPDOMAINA_INDEX = SPDOMAINA_INDEX + 1;
													SOFT_PD_A[SPDOMAINA_INDEX].ID = SPDOMAINA_INDEX;
													SOFT_PD_A[SPDOMAINA_INDEX].center_distance = currentEvent::walker_distance;
													SOFT_PD_A[SPDOMAINA_INDEX].clock = WORLD_CLOCK;
													SOFT_PD_A[SPDOMAINA_INDEX].events_type = 10;
													SOFT_PD_A[SPDOMAINA_INDEX].current_member = 2;
													SOFT_PD_A[SPDOMAINA_INDEX].member_list.member[1] = ORI_WALKERS[S0_EVENT[i + 1]];
													SOFT_PD_A[SPDOMAINA_INDEX].member_list.member[2] = ORI_WALKERS[S0_EVENT[j + 1]];
													SOFT_PD_A[SPDOMAINA_INDEX].center1.x = ORI_WALKERS[S0_EVENT[i + 1]].center.x;
													SOFT_PD_A[SPDOMAINA_INDEX].center1.y = ORI_WALKERS[S0_EVENT[i + 1]].center.y;
													SOFT_PD_A[SPDOMAINA_INDEX].center1.z = ORI_WALKERS[S0_EVENT[i + 1]].center.z;
													SOFT_PD_A[SPDOMAINA_INDEX].center1.r = ORI_WALKERS[S0_EVENT[i + 1]].radius;
													SOFT_PD_A[SPDOMAINA_INDEX].center2.x = ORI_WALKERS[S0_EVENT[j + 1]].center.x;
													SOFT_PD_A[SPDOMAINA_INDEX].center2.y = ORI_WALKERS[S0_EVENT[j + 1]].center.y;
													SOFT_PD_A[SPDOMAINA_INDEX].center2.z = ORI_WALKERS[S0_EVENT[j + 1]].center.z;
													SOFT_PD_A[SPDOMAINA_INDEX].center2.r = ORI_WALKERS[S0_EVENT[j + 1]].radius;
													ORI_WALKERS[S0_EVENT[i + 1]].statue = 1;
													ORI_WALKERS[S0_EVENT[j + 1]].statue = 1;
													ORI_WALKERS[S0_EVENT[i + 1]].regionID = SPDOMAINA_INDEX;
													ORI_WALKERS[S0_EVENT[j + 1]].regionID = SPDOMAINA_INDEX;
													// Solver
													// RANDOM_SEED = rand() / (double)(RAND_MAX);
													// SOFT_PD_A[SPDOMAINA_INDEX].events_time = INITIAL_EVENTS_TIME + INITIAL_EVENTS_TIME*RANDOM_SEED; // Replaced with solver interface once LocalSolver is provided
													SOFT_PD_A[SPDOMAINA_INDEX].events_time = getReactionTimeStamp3D(SOFT_PD_A[SPDOMAINA_INDEX].member_list, INITIAL_EVENTS_TIME);
													// Event
													EVENT_INDEX = EVENT_INDEX + 1;
													EVENT_LIST[EVENT_INDEX].domainID = SOFT_PD_A[SPDOMAINA_INDEX].ID;
													EVENT_LIST[EVENT_INDEX].domainType = 1;
													EVENT_LIST[EVENT_INDEX].timestamp = WORLD_CLOCK + SOFT_PD_A[SPDOMAINA_INDEX].events_time;
													SOFT_PD_A[SPDOMAINA_INDEX].event_card = EVENT_INDEX;
													cout << "2 walkers built a new SPDa" << endl;
													cout << ORI_WALKERS[S0_EVENT[i + 1]].ID << " and " << ORI_WALKERS[S0_EVENT[j + 1]].ID << endl;
												}
											}
										}
									}
								}
								for (i = 0;i < S0_INDEX;i++) {
									// Avoid the walkers out of boundary
									if (ORI_WALKERS[S0_EVENT[i + 1]].statue == 0 && ORI_WALKERS[S0_EVENT[i + 1]].type != -1) {
										// Three routes
										// Close S0 walkers: O1
										// Not close but maximum and uncover case (dynamic) is not included in this version
										currentEvent::O1_allow_r = 0.0;
										currentEvent::O1_neighbout_ID = 0;
										currentEvent::O2_allow_r = 0.0;
										currentEvent::O3_allow_r = 0.0;
										currentEvent::min_distance = EXTRA_LARGE;
										for (j = i + 1;j < S0_INDEX;j++) {
											// Avoid the walkers out of boundary
											if (ORI_WALKERS[S0_EVENT[j + 1]].statue == 0 && ORI_WALKERS[S0_EVENT[j + 1]].type != -1 && ORI_WALKERS[S0_EVENT[j + 1]].type != -1) {
												currentEvent::walker_distance = getDistance(ORI_WALKERS[S0_EVENT[i + 1]], ORI_WALKERS[S0_EVENT[j + 1]]) - ORI_WALKERS[S0_EVENT[i + 1]].radius - ORI_WALKERS[S0_EVENT[j + 1]].radius;
												if (currentEvent::walker_distance - currentEvent::min_distance < eps) {
													currentEvent::min_distance = currentEvent::walker_distance;
													currentEvent::O1_neighbout_ID = S0_EVENT[j + 1];
													currentEvent::O1_allow_r = getDistance(ORI_WALKERS[S0_EVENT[i + 1]], ORI_WALKERS[S0_EVENT[j + 1]])*(ORI_WALKERS[S0_EVENT[i + 1]].diffusivity) / (ORI_WALKERS[S0_EVENT[i + 1]].diffusivity + ORI_WALKERS[S0_EVENT[j + 1]].diffusivity) - ORI_WALKERS[S0_EVENT[i + 1]].radius;
													//currentEvent::O1_neighour_allow_r = getDistance1D(ORI_WALKERS[S0_EVENT[i + 1]], ORI_WALKERS[S0_EVENT[j + 1]])*(ORI_WALKERS[S0_EVENT[j + 1]].diffusivity) / (ORI_WALKERS[S0_EVENT[i + 1]].diffusivity + ORI_WALKERS[S0_EVENT[j + 1]].diffusivity) - ORI_WALKERS[S0_EVENT[j + 1]].radius;
												}
											}
										}
										// HPD surface: O2
										currentEvent::min_distance = EXTRA_LARGE;
										for (j = 0;j < HPDOMAIN_INDEX;j++) {
											if (HARD_PD[j + 1].events_type != -1) {
												currentEvent::walker_distance = getHardDomainSurfaceDistance3D(ORI_WALKERS[S0_EVENT[i + 1]], HARD_PD[j + 1]);
												if (currentEvent::walker_distance - currentEvent::min_distance < eps) {
													currentEvent::min_distance = currentEvent::walker_distance;
													currentEvent::O2_allow_r = currentEvent::walker_distance; // Surface to surface = center to surface - r
												}
											}
										}
										// SPDa close center: 03
										currentEvent::min_distance = EXTRA_LARGE;
										for (j = 0;j < SPDOMAINA_INDEX;j++) {
											if (SOFT_PD_A[j + 1].events_type != -1) {
												currentEvent::walker_distance = getSoftDomainACentertoCenterDistance(ORI_WALKERS[S0_EVENT[i + 1]], SOFT_PD_A[j + 1]) - ORI_WALKERS[S0_EVENT[i + 1]].radius;
												if (currentEvent::walker_distance - currentEvent::min_distance < eps) {
													currentEvent::min_distance = currentEvent::walker_distance;
													currentEvent::O3_allow_r = currentEvent::walker_distance; // Walker surface to SPDa center A
												}
											}
										}
										// Use the minimum none 0 one to get 3 different routes of domains, etc
										if (currentEvent::O1_allow_r - currentEvent::O2_allow_r < eps && currentEvent::O1_allow_r > eps) {
											// O1 < O2 - > O1/O3
											if (currentEvent::O1_allow_r - currentEvent::O3_allow_r < eps) {
												// O1 < O3 - > O1
												// Domain
												cout << "Create a new HPD during this event." << endl;
												HPDOMAIN_INDEX = HPDOMAIN_INDEX + 1;
												HARD_PD[HPDOMAIN_INDEX].ID = HPDOMAIN_INDEX;
												HARD_PD[HPDOMAIN_INDEX].center.x = ORI_WALKERS[S0_EVENT[i + 1]].center.x;
												HARD_PD[HPDOMAIN_INDEX].center.y = ORI_WALKERS[S0_EVENT[i + 1]].center.y;
												HARD_PD[HPDOMAIN_INDEX].center.z = ORI_WALKERS[S0_EVENT[i + 1]].center.z;
												HARD_PD[HPDOMAIN_INDEX].clock = WORLD_CLOCK;
												HARD_PD[HPDOMAIN_INDEX].events_type = 1;
												HARD_PD[HPDOMAIN_INDEX].radius = currentEvent::O1_allow_r + ORI_WALKERS[S0_EVENT[i + 1]].radius;
												HARD_PD[HPDOMAIN_INDEX].walker_ID = ORI_WALKERS[S0_EVENT[i + 1]].ID;
												HARD_PD[HPDOMAIN_INDEX].walker_type = ORI_WALKERS[S0_EVENT[i + 1]].type;
												ORI_WALKERS[S0_EVENT[i + 1]].statue = 2;
												ORI_WALKERS[S0_EVENT[i + 1]].regionID = HPDOMAIN_INDEX;
												// Solver
												HARD_PD[HPDOMAIN_INDEX].events_time = getMajorDiffTimeStamp3D(ORI_WALKERS[HARD_PD[HPDOMAIN_INDEX].walker_ID].diffusivity, HARD_PD[HPDOMAIN_INDEX].radius, ORI_WALKERS[HARD_PD[HPDOMAIN_INDEX].walker_ID].radius);
												// Event
												EVENT_INDEX = EVENT_INDEX + 1;
												EVENT_LIST[EVENT_INDEX].domainID = HARD_PD[HPDOMAIN_INDEX].ID;
												EVENT_LIST[EVENT_INDEX].domainType = 0;
												EVENT_LIST[EVENT_INDEX].timestamp = WORLD_CLOCK + HARD_PD[HPDOMAIN_INDEX].events_time;
												HARD_PD[HPDOMAIN_INDEX].event_card = EVENT_INDEX;
												cout << "1 Walker has been protected in HPD" << endl;
											}
											else {
												// O3 < O1 - > O3 or O1 = 0
												// Domain
												cout << "Create a new HPD during this event." << endl;
												HPDOMAIN_INDEX = HPDOMAIN_INDEX + 1;
												HARD_PD[HPDOMAIN_INDEX].ID = HPDOMAIN_INDEX;
												HARD_PD[HPDOMAIN_INDEX].center.x = ORI_WALKERS[S0_EVENT[i + 1]].center.x;
												HARD_PD[HPDOMAIN_INDEX].center.y = ORI_WALKERS[S0_EVENT[i + 1]].center.y;
												HARD_PD[HPDOMAIN_INDEX].center.z = ORI_WALKERS[S0_EVENT[i + 1]].center.z;
												HARD_PD[HPDOMAIN_INDEX].clock = WORLD_CLOCK;
												HARD_PD[HPDOMAIN_INDEX].events_type = 1;
												HARD_PD[HPDOMAIN_INDEX].radius = currentEvent::O3_allow_r + ORI_WALKERS[S0_EVENT[i + 1]].radius;
												HARD_PD[HPDOMAIN_INDEX].walker_ID = ORI_WALKERS[S0_EVENT[i + 1]].ID;
												HARD_PD[HPDOMAIN_INDEX].walker_type = ORI_WALKERS[S0_EVENT[i + 1]].type;
												ORI_WALKERS[S0_EVENT[i + 1]].statue = 2;
												ORI_WALKERS[S0_EVENT[i + 1]].regionID = HPDOMAIN_INDEX;
												// Solver
												HARD_PD[HPDOMAIN_INDEX].events_time = getMajorDiffTimeStamp3D(ORI_WALKERS[HARD_PD[HPDOMAIN_INDEX].walker_ID].diffusivity, HARD_PD[HPDOMAIN_INDEX].radius, ORI_WALKERS[HARD_PD[HPDOMAIN_INDEX].walker_ID].radius);
												// Event
												EVENT_INDEX = EVENT_INDEX + 1;
												EVENT_LIST[EVENT_INDEX].domainID = HARD_PD[HPDOMAIN_INDEX].ID;
												EVENT_LIST[EVENT_INDEX].domainType = 0;
												EVENT_LIST[EVENT_INDEX].timestamp = WORLD_CLOCK + HARD_PD[HPDOMAIN_INDEX].events_time;
												HARD_PD[HPDOMAIN_INDEX].event_card = EVENT_INDEX;
												cout << "1 Walker has been protected in HPD" << endl;
											}
										}
										else {
											// O2 < O1 - > O2/O3
											if (currentEvent::O2_allow_r - currentEvent::O3_allow_r < eps) {
												// O2 < O3 - > O2
												// Domain
												cout << "Create a new HPD during this event." << endl;
												HPDOMAIN_INDEX = HPDOMAIN_INDEX + 1;
												HARD_PD[HPDOMAIN_INDEX].ID = HPDOMAIN_INDEX;
												HARD_PD[HPDOMAIN_INDEX].center.x = ORI_WALKERS[S0_EVENT[i + 1]].center.x;
												HARD_PD[HPDOMAIN_INDEX].center.y = ORI_WALKERS[S0_EVENT[i + 1]].center.y;
												HARD_PD[HPDOMAIN_INDEX].center.z = ORI_WALKERS[S0_EVENT[i + 1]].center.z;
												HARD_PD[HPDOMAIN_INDEX].clock = WORLD_CLOCK;
												HARD_PD[HPDOMAIN_INDEX].events_type = 1;
												HARD_PD[HPDOMAIN_INDEX].radius = currentEvent::O2_allow_r + ORI_WALKERS[S0_EVENT[i + 1]].radius;
												HARD_PD[HPDOMAIN_INDEX].walker_ID = ORI_WALKERS[S0_EVENT[i + 1]].ID;
												HARD_PD[HPDOMAIN_INDEX].walker_type = ORI_WALKERS[S0_EVENT[i + 1]].type;
												ORI_WALKERS[S0_EVENT[i + 1]].statue = 2;
												ORI_WALKERS[S0_EVENT[i + 1]].regionID = HPDOMAIN_INDEX;
												// Solver
												HARD_PD[HPDOMAIN_INDEX].events_time = getMajorDiffTimeStamp3D(ORI_WALKERS[HARD_PD[HPDOMAIN_INDEX].walker_ID].diffusivity, HARD_PD[HPDOMAIN_INDEX].radius, ORI_WALKERS[HARD_PD[HPDOMAIN_INDEX].walker_ID].radius);
												// Event
												EVENT_INDEX = EVENT_INDEX + 1;
												EVENT_LIST[EVENT_INDEX].domainID = HARD_PD[HPDOMAIN_INDEX].ID;
												EVENT_LIST[EVENT_INDEX].domainType = 0;
												EVENT_LIST[EVENT_INDEX].timestamp = WORLD_CLOCK + HARD_PD[HPDOMAIN_INDEX].events_time;
												HARD_PD[HPDOMAIN_INDEX].event_card = EVENT_INDEX;
												cout << "1 Walker has been protected in HPD" << endl;
											}
											else {
												// O3 < O2 - > O3
												// Domain
												cout << "Create a new HPD during this event." << endl;
												HPDOMAIN_INDEX = HPDOMAIN_INDEX + 1;
												HARD_PD[HPDOMAIN_INDEX].ID = HPDOMAIN_INDEX;
												HARD_PD[HPDOMAIN_INDEX].center.x = ORI_WALKERS[S0_EVENT[i + 1]].center.x;
												HARD_PD[HPDOMAIN_INDEX].center.y = ORI_WALKERS[S0_EVENT[i + 1]].center.y;
												HARD_PD[HPDOMAIN_INDEX].center.z = ORI_WALKERS[S0_EVENT[i + 1]].center.z;
												HARD_PD[HPDOMAIN_INDEX].clock = WORLD_CLOCK;
												HARD_PD[HPDOMAIN_INDEX].events_type = 1;
												HARD_PD[HPDOMAIN_INDEX].radius = currentEvent::O3_allow_r + ORI_WALKERS[S0_EVENT[i + 1]].radius;
												HARD_PD[HPDOMAIN_INDEX].walker_ID = ORI_WALKERS[S0_EVENT[i + 1]].ID;
												HARD_PD[HPDOMAIN_INDEX].walker_type = ORI_WALKERS[S0_EVENT[i + 1]].type;
												ORI_WALKERS[S0_EVENT[i + 1]].statue = 2;
												ORI_WALKERS[S0_EVENT[i + 1]].regionID = HPDOMAIN_INDEX;
												// Solver
												HARD_PD[HPDOMAIN_INDEX].events_time = getMajorDiffTimeStamp3D(ORI_WALKERS[HARD_PD[HPDOMAIN_INDEX].walker_ID].diffusivity, HARD_PD[HPDOMAIN_INDEX].radius, ORI_WALKERS[HARD_PD[HPDOMAIN_INDEX].walker_ID].radius);
												// Event
												EVENT_INDEX = EVENT_INDEX + 1;
												EVENT_LIST[EVENT_INDEX].domainID = HARD_PD[HPDOMAIN_INDEX].ID;
												EVENT_LIST[EVENT_INDEX].domainType = 0;
												EVENT_LIST[EVENT_INDEX].timestamp = WORLD_CLOCK + HARD_PD[HPDOMAIN_INDEX].events_time;
												HARD_PD[HPDOMAIN_INDEX].event_card = EVENT_INDEX;
												cout << "1 Walker has been protected in HPD" << endl;
											}
										}
									}
								}
							}// End of N > 1 else and all N if...else
						}
					}
					if (currentEvent::domain_type == 2) {
						// Move EVENT_POINTER to next event
						EVENT_POINTER = EVENT_POINTER + 1;
						cout << "FATAL ERROR 700: Unable to process Type 2 SPD in this version. Exit." << endl;
						system("pause");
						exit(0);
					}
				}
				// End of event
				// Run event list sort
				// Use malloc with pointer instead of assigned array: Not in this version
				// Use type queue instead of type array: Not in this version
				for (i = 1; i < EVENT_INDEX;i++) {
					for (j = 1;j < EVENT_INDEX + 1 - i;j++) {
						if (EVENT_LIST[j].timestamp>EVENT_LIST[j + 1].timestamp) {
							swap(EVENT_LIST[j], EVENT_LIST[j + 1]);
						}
					}
				}
				for (i = 1;i < EVENT_INDEX;i++) {
					if (EVENT_LIST[i].domainType == 0) {
						HARD_PD[EVENT_LIST[i].domainID].event_card = i;
					}
					if (EVENT_LIST[i].domainType == 1) {
						SOFT_PD_A[EVENT_LIST[i].domainID].event_card = i;
					}
				}
				// One event has been processed with its executor and creator
				// Check the clock
				if (WORLD_CLOCK - END_CLOCK > eps) {
					CONTROL_FLAG = 0; // Check with clock
				}
				// Inner Dump Debug
				if (DUMPFLAG_INN) {
					ofstream file2("InnerDump.txt", ios::app);
					streambuf *f2 = cout.rdbuf(file2.rdbuf());
					DUMP_COUNTER = 0;
					cout << WORLD_CLOCK << endl;
					for (i = 1;i < WALKER_INDEX + 1;i++) {
						if (ORI_WALKERS[i].type != -1) {
							DUMP_COUNTER = DUMP_COUNTER + 1;
							if (space_dimension == 3) {
								//cout << ORI_WALKERS[i].ID << " " << ORI_WALKERS[i].type << " " << ORI_WALKERS[i].component << " " << ORI_WALKERS[i].center.x << " " << ORI_WALKERS[i].center.y << " " << ORI_WALKERS[i].center.z << endl;
								cout << " " << ORI_WALKERS[i].type << " " << ORI_WALKERS[i].component << " " << ORI_WALKERS[i].center.x << " " << ORI_WALKERS[i].center.y << " " << ORI_WALKERS[i].center.z << endl;
							}
						}
					}
					//cout << "Remain " << DUMP_COUNTER << " walkers." << endl;
					//cout << "  " << DUMP_COUNTER << endl;
					cout << endl;
					cout.rdbuf(f2);
				}
				// OS2
				ofstream file3("timelog.txt", ios::app);
				streambuf *f3 = cout.rdbuf(file3.rdbuf());
				cout << WORLD_CLOCK << " " << DUMP_COUNTER << endl;
				cout.rdbuf(f3);
				// OS3
				ofstream file4("Original_voids.txt", ios::app);
				streambuf *f4 = cout.rdbuf(file4.rdbuf());
				long dump_voids;
				dump_voids = 0;
				for (i = 1;i < WALKER_INDEX + 1;i++) {
					if (ORI_WALKERS[i].type == 0) {
						if (space_dimension == 3) {
							dump_voids = dump_voids + ORI_WALKERS[i].component;
						}
					}
				}
				cout << WORLD_CLOCK << " " << dump_voids << endl;
				cout.rdbuf(f4);
				// OS4
				ofstream file5("Total_voids.txt", ios::app);
				streambuf *f5 = cout.rdbuf(file5.rdbuf());
				dump_voids = 0;
				for (i = 1;i < WALKER_COUNTER + 1;i++) {
					if (ORI_WALKERS[i].type == 0) {
						if (space_dimension == 3) {
							dump_voids = dump_voids + ORI_WALKERS[i].component;
						}
					}
				}
				cout << WORLD_CLOCK << " " << dump_voids << endl;
				cout.rdbuf(f5);
				// OS5
				if (DUMPLAMMPSFLAG_INN) {
					DUMPCOUNTER_INN = DUMPCOUNTER_INN + 1;
					if (DUMPCOUNTER_INN == DUMPCOUNTER_FLAG) {
						DUMPCOUNTER_INN = 0;
						ofstream file6("LAMMPS_dump.txt", ios::app);
						streambuf *f6 = cout.rdbuf(file6.rdbuf());
						long lammps_type = 0;
						long lammps_number = 0;
						cout << "ITEM: TIMESTEP" << endl;
						//TIME_STEP = TIME_STEP + 1;
						cout << WORLD_CLOCK << endl;
						cout << "ITEM: NUMBER OF ATOMS" << endl;
						for (i = 1;i < WALKER_COUNTER + 1;i++) {
							if (ORI_WALKERS[i].type != -1) {
								lammps_number = lammps_number + 1;
							}
						}
						cout << lammps_number << endl;
						cout << "ITEM: BOX BOUNDS pp pp pp" << endl;
						cout << X_LOWWER_LIMIT << " " << X_UPPER_LIMIT << endl;
						cout << Y_LOWWER_LIMIT << " " << Y_UPPER_LIMIT << endl;
						cout << Z_LOWWER_LIMIT << " " << Z_UPPER_LIMIT << endl;
						cout << "ITEM: ATOMS id type x y z c_ke c_pe " << endl;
						for (i = 1;i < WALKER_COUNTER + 1;i++) {
							if (ORI_WALKERS[i].type != -1) {
								if (ORI_WALKERS[i].type == 0) {
									lammps_type = ORI_WALKERS[i].component;
								}
								if (ORI_WALKERS[i].type == 1) {
									lammps_type = ORI_WALKERS[i].component + 100;
								}
								cout << ORI_WALKERS[i].ID << " " << lammps_type << " " << ORI_WALKERS[i].center.x << " " << ORI_WALKERS[i].center.y << " " << ORI_WALKERS[i].center.z << " 0.0 " << ORI_WALKERS[i].energy << endl;
							}
						}
						cout.rdbuf(f6);
					}	
				}
				// OS6
				// End of inner dump
			}
			STAGE_FLAG = 0;
			cout << "Current calculation has finished at " << WORLD_CLOCK << "s, aCRD has reached the end of clock, extend the calculation? 1:0." << endl;
			cin >> INPUT_FLAG;
			if (INPUT_FLAG) {
				cout << "Please input the extension time in the unit of second. " << endl;
				cin >> EXTRA_CLOCK;
				END_CLOCK = END_CLOCK + EXTRA_CLOCK;
				STAGE_FLAG = 1;
			}
		}
	}

	cout << "Finish Event Loop." << endl;
	cout << "//////////////////////////////////////////////////////////////////" << endl;
	// Final Sync
	cout << "Begin to do final synchronization." << endl;

	// Run E02 to all HPDs with 
	// getSynchronization();
	// HPDs run E02
	// SPD run E10_END, in current version, do not handle SPDa as Trigger time is very short after SPDa built
	for (i = 0;i < WALKER_COUNTER;i++) {
		if (ORI_WALKERS[i + 1].statue == 2) {
			// HPD
			if (space_dimension == 1) {
				ORI_WALKERS[i + 1].statue = 0;
				ORI_WALKERS[i + 1].center = getTransDiffRelativeLocation1D(ORI_WALKERS[i + 1].diffusivity, WORLD_CLOCK - ORI_WALKERS[i + 1].timestamp, HARD_PD[ORI_WALKERS[i + 1].regionID].radius, ORI_WALKERS[i + 1].radius, ORI_WALKERS[i + 1].center);
			}
			if (space_dimension == 3) {
				ORI_WALKERS[i + 1].statue = 0;
				ORI_WALKERS[i + 1].center = getTransDiffRelativeLocation3D(ORI_WALKERS[i + 1].diffusivity, WORLD_CLOCK - ORI_WALKERS[i + 1].timestamp, HARD_PD[ORI_WALKERS[i + 1].regionID].radius, ORI_WALKERS[i + 1].radius, ORI_WALKERS[i + 1].center);
			}
			
		}
		if (ORI_WALKERS[i + 1].statue == 1) {
			// SPDa
			if (space_dimension == 1) {
				for (j = 0;j < SOFT_PD_A[ORI_WALKERS[i + 1].regionID].current_member;j++) {
					ORI_WALKERS[SOFT_PD_A[ORI_WALKERS[i + 1].regionID].member_list.member[j + 1].ID].statue = 0;
					ORI_WALKERS[SOFT_PD_A[ORI_WALKERS[i + 1].regionID].member_list.member[j + 1].ID].type = -1;
				}
			}
			if (space_dimension == 3) {
				// Don't do anything
			}
		}
	}
	cout << "Finish synchronization." << endl;
	cout << "//////////////////////////////////////////////////////////////////" << endl;

	// Debug Goto
	cout << "Test Again? 1 to Restart at File Loading, 2 to Restart at Introduce Solvers and Event Queue, 3 to Restart at Event Loop, 0 to continue." << endl;
	int Test_Flag;
	cin >> Test_Flag;
	if (Test_Flag == 1) goto Test_Start_1;
	if (Test_Flag == 2) goto Test_Start_2;
	if (Test_Flag == 3) goto Test_Start_3;
	//if (Test_Flag == 4) goto Test_Start_4;


//Test_Start_4:
	// Dump Data
	ofstream file("Dump.txt", ios::app);
	streambuf *f = cout.rdbuf(file.rdbuf());
	cout << WORLD_CLOCK << endl;
	for (i = 1;i < WALKER_INDEX;i++) {
		if (ORI_WALKERS[i].type != -1) {
			if (space_dimension == 3) {
				cout << ORI_WALKERS[i].type << " " << ORI_WALKERS[i].component << " " << ORI_WALKERS[i].center.x << " " << ORI_WALKERS[i].center.y << " " << ORI_WALKERS[i].center.z << endl;
			}
			if (space_dimension == 1) {
				cout << ORI_WALKERS[i].type << " " << ORI_WALKERS[i].component << " " << ORI_WALKERS[i].center.x << endl;
			}
		}
	}
	cout << endl;
	cout.rdbuf(f);
	cout << "Code End." << endl;
	// Exit aCRD
	cout << "Exit aCRD Framework." << endl;
	releaseModule();
	system("pause");
}