#pragma once
// Version 1 alpha demo
// This file related to aCRD/Main.cpp
// This file defines several basic functions

#ifndef TOOLBOX_H
#define TOOLBOX_H

// System function
void releaseModule(){
	cout << "Release Module from System." << endl;
	// Module Releasing
	if (dllSolverDissociation != NULL)
	{
		FreeLibrary(dllSolverDissociation);
	}
	if (dllSolverLocal != NULL)
	{
		FreeLibrary(dllSolverLocal);
	}
	if (dllSolverTransientDiff != NULL)
	{
		FreeLibrary(dllSolverTransientDiff);
	}
	if (dllSolverMajorDiff != NULL)
	{
		FreeLibrary(dllSolverMajorDiff);
	}
	if (dllLibDiffusivity != NULL)
	{
		FreeLibrary(dllLibDiffusivity);
	}
	cout << "Have released all Modules." << endl;
} // Done

// Center to Center, 3D
double getDistance(walker WALKER_1, walker WALKER_2) {
	double l = 0.0;
	l = sqrt((WALKER_1.center.x - WALKER_2.center.x)*(WALKER_1.center.x - WALKER_2.center.x) + (WALKER_1.center.y - WALKER_2.center.y)*(WALKER_1.center.y - WALKER_2.center.y) + (WALKER_1.center.z - WALKER_2.center.z)*(WALKER_1.center.z - WALKER_2.center.z));
	return l;
} // Done

// Center to Center
double getDistance1D(walker WALKER_1, walker WALKER_2) {
	double l = 0.0;
	l = sqrt((WALKER_1.center.x - WALKER_2.center.x)*(WALKER_1.center.x - WALKER_2.center.x));
	return l;
}

// Return the close walker surface distance f: walker surface and SPDA center walker surface, 3D
double getSoftDomainASurfaceDistance(walker WALKER_NEW, soft_domain_a SPDA) {
	double l1 = 0.0;
	double l2 = 0.0;
	double l = 0.0;
	l1 = sqrt((WALKER_NEW.center.x - SPDA.center1.x)*(WALKER_NEW.center.x - SPDA.center1.x) + (WALKER_NEW.center.y - SPDA.center1.y)*(WALKER_NEW.center.y - SPDA.center1.y) + (WALKER_NEW.center.z - SPDA.center1.z)*(WALKER_NEW.center.z - SPDA.center1.z));
	l1 = l1 - SPDA.center1.r - WALKER_NEW.radius;
	l2 = sqrt((WALKER_NEW.center.x - SPDA.center2.x)*(WALKER_NEW.center.x - SPDA.center2.x) + (WALKER_NEW.center.y - SPDA.center2.y)*(WALKER_NEW.center.y - SPDA.center2.y) + (WALKER_NEW.center.z - SPDA.center2.z)*(WALKER_NEW.center.z - SPDA.center2.z));
	l2 = l2 - SPDA.center2.r - WALKER_NEW.radius;
	if (l1 - l2 < eps) {
		l = l1;
	}
	else {
		l = l2;
	}
	return l;
}

// Return the close walker surface distance f: walker surface and SPDA center walker surface, 1D
double getSoftDomainASurfaceDistance1D(walker WALKER_NEW, soft_domain_a SPDA) {
	double l1 = 0.0;
	double l2 = 0.0;
	double l = 0.0;
	l1 = sqrt((WALKER_NEW.center.x - SPDA.center1.x)*(WALKER_NEW.center.x - SPDA.center1.x));
	l1 = l1 - SPDA.center1.r - WALKER_NEW.radius;
	l2 = sqrt((WALKER_NEW.center.x - SPDA.center2.x)*(WALKER_NEW.center.x - SPDA.center2.x));
	l2 = l2 - SPDA.center2.r - WALKER_NEW.radius;
	if (l1 - l2 < eps) {
		l = l1;
	}
	else {
		l = l2;
	}
	return l;
}

// Return the close walker center to SPD center walker center, 3D
double getSoftDomainACentertoCenterDistance(walker WALKER_NEW, soft_domain_a SPDA) {
	double l1 = 0.0;
	double l2 = 0.0;
	double l = 0.0;
	l1 = sqrt((WALKER_NEW.center.x - SPDA.center1.x)*(WALKER_NEW.center.x - SPDA.center1.x) + (WALKER_NEW.center.y - SPDA.center1.y)*(WALKER_NEW.center.y - SPDA.center1.y) + (WALKER_NEW.center.z - SPDA.center1.z)*(WALKER_NEW.center.z - SPDA.center1.z));
	l2 = sqrt((WALKER_NEW.center.x - SPDA.center2.x)*(WALKER_NEW.center.x - SPDA.center2.x) + (WALKER_NEW.center.y - SPDA.center2.y)*(WALKER_NEW.center.y - SPDA.center2.y) + (WALKER_NEW.center.z - SPDA.center2.z)*(WALKER_NEW.center.z - SPDA.center2.z));
	if (l1 - l2 < eps) {
		l = l1;
	}
	else {
		l = l2;
	}
	return l;
}

// Return the close walker center to SPD center walker center
double getSoftDomainACentertoCenterDistance1D(walker WALKER_NEW, soft_domain_a SPDA) {
	double l1 = 0.0;
	double l2 = 0.0;
	double l = 0.0;
	l1 = sqrt((WALKER_NEW.center.x - SPDA.center1.x)*(WALKER_NEW.center.x - SPDA.center1.x));
	l2 = sqrt((WALKER_NEW.center.x - SPDA.center2.x)*(WALKER_NEW.center.x - SPDA.center2.x));
	if (l1 - l2 < eps) {
		l = l1;
	}
	else {
		l = l2;
	}
	return l;
}

// The current version use uniform approx.
double getEstimateRadius(walker WALKER) {
	double radius = 0.0;
	if (WALKER.type == 0) {
		if (WALKER.component == 1) {
			radius = FE_BCC_A;
		}
		if (WALKER.component == 2) {
			radius = FE_BCC_A * 2.0;
		}
		if (WALKER.component == 3) {
			radius = FE_BCC_A * 2.5;
		}
		if (WALKER.component > 3) {
			radius = VACANCY_MODIFIER * pow((WALKER.component*pow(FE_BCC_A, 3.0)), 1.0 / 3.0);
		}
	}
	if (WALKER.type == 1) {
		if (WALKER.component == 1) {
			radius = FE_RADIUS;
		}
		if (WALKER.component == 2) {
			radius = FE_RADIUS * 2.0;
		}
		if (WALKER.component == 3) {
			radius = FE_RADIUS * 2.5;
		}
		if (WALKER.component > 3) {
			radius = INTERST_MODIFIER * pow((WALKER.component*pow(FE_RADIUS, 3.0)), 1.0 / 3.0);
		}
	}
	return radius;
} // Done

double getHardDomainCenterDistance(hard_domain DOMAIN_1, hard_domain DOMAIN_2) {
	double distance = 0.0;
	distance = sqrt((DOMAIN_1.center.x - DOMAIN_2.center.x)*(DOMAIN_1.center.x - DOMAIN_2.center.x) + (DOMAIN_1.center.y - DOMAIN_2.center.y)*(DOMAIN_1.center.y - DOMAIN_2.center.y) + (DOMAIN_1.center.z - DOMAIN_2.center.z)*(DOMAIN_1.center.z - DOMAIN_2.center.z));
	return distance;
}

// Part of NNS function, remove once applies: Not in this version
long getNeighborDomainID(location DOMAIN_CENTER) {
	return 0;
}

// Surface to Surface
double getHardDomainSurfaceDistance1D(walker WALKER, hard_domain HPD) {
	double distance = 0.0;
	distance = sqrt((WALKER.center.x - HPD.center.x)*(WALKER.center.x - HPD.center.x)) - HPD.radius - WALKER.radius;
	return distance;
}

// Surface to Surface
double getHardDomainSurfaceDistance3D(walker WALKER, hard_domain HPD) {
	double distance = 0.0;
	distance = sqrt((WALKER.center.x - HPD.center.x)*(WALKER.center.x - HPD.center.x) + (WALKER.center.y - HPD.center.y)*(WALKER.center.y - HPD.center.y) + (WALKER.center.z - HPD.center.z)*(WALKER.center.z - HPD.center.z)) - HPD.radius - WALKER.radius;
	return distance;
}

#endif // !TOOLBOX_H
