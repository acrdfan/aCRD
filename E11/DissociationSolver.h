// Dissociation Solver
// In MOD_SOL_aCRD
// Version 1 alpha demo
// J. Fan 2019

#ifdef DISSOCIATIONSOLVER_EXPORTS
#define DISSOCIATIONSOLVER_API __declspec(dllexport)
#else
#define DISSOCIATIONSOLVER_API __declspec(dllimport)
#endif

// The Definition of Main Structs
struct location {
	double x;
	double y;
	double z;
};

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

#ifdef __cplusplus
extern "C" {
#endif
	DISSOCIATIONSOLVER_API int getDissociationFlag(walker WALKER);
	extern DISSOCIATIONSOLVER_API double nE11SolverVersion;
#ifdef __cplusplus
}
#endif