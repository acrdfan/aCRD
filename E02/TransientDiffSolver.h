// Transient Diffusion Green's Function Solver
// In MOD_SOL_aCRD
// Version 1 alpha demo
// J. Fan 2018-2019

#ifdef TRANSIENTDIFFSOLVER_EXPORTS
#define TRANSIENTDIFFSOLVER_API __declspec(dllexport)
#else
#define TRANSIENTDIFFSOLVER_API __declspec(dllimport)
#endif

struct location {
	double x;
	double y;
	double z;
};

#ifdef __cplusplus
extern "C" {
#endif
	TRANSIENTDIFFSOLVER_API location getTransDiffRelativeLocation1D(double Dv, double timestamp, double domainRadius, double walkerRadius, location centerLocation);
	TRANSIENTDIFFSOLVER_API location getTransDiffRelativeLocation3D(double Dv, double timestamp, double domainRadius, double walkerRadius, location centerLocation);
	extern TRANSIENTDIFFSOLVER_API double nE02SolverVersion;
#ifdef __cplusplus
}
#endif