// Local Reaction Possibility Solver
// In MOD_SOL_aCRD
// Version 1 alpha demo
// J. Fan 2018-2019

#ifdef LOCALSOLVER_EXPORTS
#define LOCALSOLVER_API __declspec(dllexport)
#else
#define LOCALSOLVER_API __declspec(dllimport)
#endif

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

// For current version, only pre-assigned final status: 1 - 50 walkers
struct reaction_current_walker {
	walker member[50];
};

// For current version, only pre-assigned final status: 1 - 5 walkers
struct reaction_final_walker {
	walker member[6];
};

#ifdef __cplusplus
extern "C" {
#endif
	LOCALSOLVER_API double getReactionTimeStamp1D(reaction_current_walker walkers, double initial_time);
	LOCALSOLVER_API double getReactionTimeStamp3D(reaction_current_walker walkers, double initial_time);
	LOCALSOLVER_API reaction_final_walker getReactionFinalStatues1D(reaction_current_walker walkers, int init_member, domaincenter centerA, domaincenter centerB);
	LOCALSOLVER_API reaction_final_walker getReactionFinalStatues3D(reaction_current_walker walkers, int init_member, domaincenter centerA, domaincenter centerB);
	extern LOCALSOLVER_API double nE10SolverVersion;
#ifdef __cplusplus
}
#endif