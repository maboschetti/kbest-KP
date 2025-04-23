//-----------------------------------------------------------------------------
//  File: CoinObj.h                                                       
//                                                                           
//  Project: Coin-OR (CBC) solver interface (header)                   
//                                                                           
//  Author: M.A. Boschetti                                                   
//                                                                           
//  Last update: 06.11.2012                                                  
//-----------------------------------------------------------------------------
//
//extern "C" {
//#include "gurobi_c.h"
//}

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <cassert>
#include <iomanip>

// For Branch and bound
#include "OsiSolverInterface.hpp"
#include "CbcModel.hpp"
#include "CoinModel.hpp"

// For all different
#include "CbcBranchCut.hpp"
#include "CbcBranchActual.hpp"
#include "CbcBranchAllDifferent.hpp"
#include "CbcCutGenerator.hpp"
#include "OsiClpSolverInterface.hpp"
#include "CglStored.hpp"

// For Branch and bound
//#include "OsiSolverInterface.hpp"
//#include "CbcModel.hpp"
//#include "CbcBranchUser.hpp"
//#include "CbcCompareUser.hpp"
#include "CbcCutGenerator.hpp"
#include "OsiClpSolverInterface.hpp"

// Cuts
#include "CglGomory/CglGomory.hpp"
#include "CglProbing/CglProbing.hpp"
#include "CglKnapsackCover/CglKnapsackCover.hpp"
#include "CglOddHole/CglOddHole.hpp"
#include "CglClique/CglClique.hpp"
#include "CglFlowCover/CglFlowCover.hpp"
#include "CglMixedIntegerRounding/CglMixedIntegerRounding.hpp"
#include "CglMixedIntegerRounding2/CglMixedIntegerRounding2.hpp"
#include "CglAllDifferent/CglAllDifferent.hpp"
#include "CglLiftAndProject/CglLiftAndProject.hpp"
#include "CglTwomir/CglTwomir.hpp"

// Node selection
#include "CbcCompareDepth.hpp"
#include "CbcCompareObjective.hpp"
#include "CbcCompareDefault.hpp"
#include "CbcCompareEstimate.hpp"

// Branching
#include "CbcBranchDefaultDecision.hpp"
#include "CbcBranchDynamic.hpp"


// Heuristics
#include "CbcHeuristic.hpp"
#include "CbcHeuristicLocal.hpp"
#include "CbcHeuristicRINS.hpp"
#include "CbcHeuristicFPump.hpp"
#include "CbcHeuristicGreedy.hpp"
#include "CbcHeuristicPivotAndFix.hpp"
#include "CbcHeuristicDiveCoefficient.hpp"
#include "CbcHeuristicDiveFractional.hpp"
#include "CbcHeuristicDiveVectorLength.hpp"

#include  "CoinTime.hpp"


struct Coincutinfo {
	int FlagLP;   // FlagLP=0 if we are solving the original model; FlagLP=1 if we are solving the knapsack subproblem; 
	OsiClpSolverInterface env;  
	CoinModel lp;
	int      numcols;
	int      numtoadd;
	int      maxcuts;
	int      num;
	int      nz;
	double   *x;
	int      *beg;
	int      *ind; 
	double   *val;
	double   *rhs;

};

typedef struct Coincutinfo CoinCUTINFO, *CoinCUTINFOptr;

// Function implementing a custom cutting plane procedure called by callback
int __stdcall myNewcutcallback (OsiClpSolverInterface env, void *cbdata, int wherefrom, void *cbhandle);

class CoinObj
{
private:

	int status;
	int bolfree;  // Is memory allocated inside the object?
	double TimeLimit;

public:

	OsiClpSolverInterface env;  // Solver Interface  
	CoinModel lp;               // Model
	CoinCUTINFO cutinfo;        // Data structure used to add cutting planes
	
	// Input
	int minmax;      // Mininum or Maximum (1=minimize; -1=maximize) 
	int ncols;       // Number of columns/variables
	int nrows;       // Number of rows/constraints
	int nz;          // Number of non-zero coefficients into the constraint matrix
	int nbin;        // Number of binary variables
	int nint;        // Number of integer variables
	int objsen;      // Minimization/Maximization Flag (+1=Min; -1=Max)
	int probtype;    // Problem Type (-1: err; CPXPROB_LP=0; CPXPROB_MILP=1; etc.)
	double *obj;     // Objective Function coefficients
	double *rhs;     // Right Hand Side
	char *sense;     // Constraint Senses
	int *matbeg;     // Pointers to the beginning of each columns of the constraint matrix
	int *matcnt;     // Cardinality of each columns of the constraint matrix
	int *matind;     // Row index associated to each coefficient of the constraint matrix
	double *matval;  // Value of each coefficient of the constraint matrix
	int *rmatbeg;    // Pointers to the beginning of each rows of the constraint matrix
	int *rmatcnt;    // Cardinality of each rows of the constraint matrix
	int *rmatind;    // Column index associated to each coefficient of the constraint matrix
	double *rmatval; // Value of each coefficient of the constraint matrix
	double *lb;      // Lower Bound value of each variable
	double *ub;      // Upper Bound value of each variable
	int *indices;
	int *priority;
	int *direction;
	char *xctype;
	double *rngval;

	int nsos;       // Number of SOS (Special Ordered Set) 
	int nsosnz;     // Total number of membres in all SOS added
	char *typesos;  // Array containing SOS type (CPX_TYPE_SOS1='1', CPX_TYPE_SOS2='2')
	int *sosbeg;    // Pointers to the beginning of each SOS
	int *sosind;    // Index associated to the SOS member
	double *soswt;  // Weight associated to the SOS member

	// Output
	double *xt;     // Variable Values
	double objval;  // Objective Function Value

	// Constructor and Destructor
	CoinObj(void);
	~CoinObj(void);

	//  Open: open an environment and open a new problem 
	int CoinObj::Open(void);

	// Reset: close the existing problem and create a new one 
	int CoinObj::Reset(void);

	// MallocCols: Allocate the "column" data structure
	int CoinObj::MallocCols(void);

	// MallocRows: Allocate the "row" data structure
	int CoinObj::MallocRows(void);

	// FreeCols: Deallocate the "column" data structure
	int CoinObj::FreeCols(void);

	// FreeRows: Deallocate the "row" data structure
	int CoinObj::FreeRows(void);

	// ReadMPS
	int CoinObj::ReadMPS(const char * fname);

	// CopyLP
	int CoinObj::CopyLPbyCols(void);
	int CoinObj::CopyLPbyRows(void);

	// Set a Lower Bound for the MIP
	int CoinObj::SetLB(double LB);

	// SetMIP
	int CoinObj::SetMIP(double TLim);

	// Add Fixed Cost Miner Constraints
	int CoinObj::AddFCM(int nrange);

	// Sentinel Constraint
	int CoinObj::Sentinel(long n, long m);

	// SolveMIP
	int CoinObj::SolveMIP(double *Zopt, double *Zlb, double *Gap, long *Nodes, long *Cuts);

    // Setup Custom Cutting Plane
	int CoinObj::CuttingPlane(long n, long m, double* u, double* p, double* pmax, 
		                       long *nY, long **Y, long *nUOR, long **UOR, long *nUAND, long **UAND, long *FDepA, long *FDepP, 
							   int Pred, int Cover, int Lifting, int KnapSol, int MaxCuts, long *Cuts);
	int CoinObj::CuttingPlaneHeu(long n, long m, double* u, double* p, double* pmax, 
  	                          long *nY, long **Y, long *nUOR, long **UOR, long *nUAND, long **UAND, long *FDepA, long *FDepP,
	                          int Pred, int Cover, int Lifting, int KnapSol, int MaxCuts, long *Cuts, long fact);

};