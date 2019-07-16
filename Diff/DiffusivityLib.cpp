// Diffusivity Library
// In MOD_LIB_aCRD
// Version 0 alpha demo
// This version should only be regarded as a demonstration of function's interface, all the parameters may be way off. 
// J. Fan 2018

#include "stdafx.h"
#include "DiffusivityLib.h"

using namespace std;

// Lib Version
DIFFUSIVITYLIB_API double nDvLibVersion = 0.5;

// Input Signiture: Energy (Unit in keV), Type (0 for V, 1 for I), Size (Cluster)
// Diffusivity: nm^2/s
double getDvValue(double E, int type, long size)
{
	double WalkerDiffusivity = 0;
	switch (type) {
	case 0:
		// Estimated Single 1000K, Phys. Rev. B 80, 144111
		// TODO
		switch (size) {
		case 1:
			WalkerDiffusivity = 17.22;
			break;
		case 2:
			WalkerDiffusivity = 8.610;
			break;
		case 3:
			WalkerDiffusivity = 5.740;
			break;
		case 4:
			WalkerDiffusivity = 4.305;
			break;
		case 5:
			WalkerDiffusivity = 3.444;
			break;
		case 6:
			WalkerDiffusivity = 2.870;
			break;
		case 7:
			WalkerDiffusivity = 2.460;
			break;
		case 8:
			WalkerDiffusivity = 0.359;
			break;
		default:
			WalkerDiffusivity = 0.003;
		}
		break;
	case 1:
		// Estimated Single 1000K, A Summary Report. Tech. Reports, Tohoku Univ. Vol. 47, No.2 (1982) p. 215
		// Self-diffusion in alpha-iron,Acta Metall. mater. Vol. 38, No. 2(1990), p. 283-292
		// TODO
		switch (size) {
		case 1:
			WalkerDiffusivity = 6.170;
			break;
		case 2:
			WalkerDiffusivity = 1.480;
			break;
		case 3:
			WalkerDiffusivity = 0.678;
			break;
		default:
			WalkerDiffusivity = 0.050;
		}
		break;
	default:
		WalkerDiffusivity = 0.000;
	}
	return WalkerDiffusivity;
}


