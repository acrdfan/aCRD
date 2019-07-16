// Major Diffusion Green's Function Solver
// In MOD_SOL_aCRD
// Version 1 alpha demo
// J. Fan 2018-2019

#ifdef MAJORDIFFSOLVER_EXPORTS
#define MAJORDIFFSOLVER_API __declspec(dllexport)
#else
#define MAJORDIFFSOLVER_API __declspec(dllimport)
#endif

struct location {
	double x;
	double y;
	double z;
};

#ifdef __cplusplus
extern "C" {
#endif
	MAJORDIFFSOLVER_API double getMajorDiffTimeStamp1D(double Dv, double domainRadius, double walkerRadius);
	MAJORDIFFSOLVER_API double getMajorDiffTimeStamp3D(double Dv, double domainRadius, double walkerRadius);
	MAJORDIFFSOLVER_API location getMajorDiffRelativeLocation1D(double Dv, double timestamp, double domainRadius, double walkerRadius, location centerLocation);
	MAJORDIFFSOLVER_API location getMajorDiffRelativeLocation3D(double Dv, double timestamp, double domainRadius, double walkerRadius, location centerLocation);
	extern MAJORDIFFSOLVER_API double nE01SolverVersion;
#ifdef __cplusplus
	}
#endif