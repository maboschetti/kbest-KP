//-----------------------------------------------------------------------------
// File: subgradiente.h
//
// Dual Ascent based on a Lagrangian Column Generation for 
// Bin Packing Problem and its variants
//
// Authors: Marco A. Boschetti and S. Novellani
//
// Last update: 20.04.2025
//-----------------------------------------------------------------------------
//
// Solving Bin Packing Problem and some of its variants
//
// 1. Basic Bin Packing Problem
//    Type = 0
//    Note: No constraints
//
// 2. Cardinality constrained bin packing problem
//    Type = 1
//    Note : at most k items can be assigned to one bin
//
// 3. Class constrained bin packing problem
//    Type = 2
//    Note : each bin must contain at most k different classes (at each item is assigned a class)
//
// 4. Bin packing problem with conflicts
//    Type = 3
//    Note : some pairs of items may be in conflict and cannot be assigned to the same bin (we generate a conflict matrix)
//
// 5. Two-dimensional vector packing problem
//    Type = 4
//    Note : two capacity constraints must be satisfied in a feasible packing solution for each bin
//
//-----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "..\kbest-KP-Library\kBest-KP-Library.h"

// Gold Parameters
//#define  Prec       0.01
//#define  Inf        1000000.
//#define  FixPen     10000.
#define  Prec        0.0001
//#define  Infi       1000000000  (redefined in kBest-KP-Library.h... temporary fix)
#define  Inf         1e+10
#define  MaxVar      10000
#define  MaxCon      500
#define  MaxNZ       1000000
#define  Iprinx      1
#define  MaxColCard  100
#define  MaxColGen   100000

struct InstanceData
{
	long Type = 0;  // Type of Bin Packing Problem
	long k;         // Parameter k to setup constrants 
	long W;         // Capacity of the bin
	long n;         // Number of items
	long *w;        // Weights of items
	long W2;        // Second capacity of the bin (for 2D vector packing)
	long *w2;       // Second Weights of items (for 2D vector packing)
	long CL;        // Maximum number of classes for each bin
	long *cl;       // Classes of items
	long **CMatrix; // Conflit matrix
};

struct Parameters
{
	long Model = 1;               // Model: 0 = Set Partitioning; 1 = Set Covering
	float TimeLimit = 300.0;      // Time limit for every algorithm
	long MaxNumColumns = 1000000; // Maximum number of columns
	long MaxNumColMIP = 100000;   // Maximum number of columns for MIP
	double alpha = 10.0;          // Starting alpha for the subgradient step rule
	long MaxItMinGap = 200;       // Maximum number of iterations without improvement for the subgradient
	double MinGap = 0.001;        // MipGap for the subgradient (the minimum imporvement detected)
	double PercK = 0.1;           // k = PercK * n (where n is the number of items)
	double PercDelta = 0.1;       // Delta = PercDelta * LB (that is the maximum reduced cost)
};

class Subgradient
{
private:

	// Prototypes
	int Malloc(void);
	int Free(void);
	long GetList(FILE *fin, long *nlist, int *list);
	double random(void);
	int Generate_SR3(double *mip);

public:

	// Global Variable Declarations

	// Input: 
	struct InstanceData Inst;
	struct Parameters Param;

	// Knapsack Problem Object
	KPSolver KPSol;

	// Working variables
	double FixPen;
	double MaxM;
	double zub;
	int n, n1, n2; // n: Number of Columns (original+artificial); n1: Number of original columns; n2: Number of heuristic solution columns
	int m, mc;     // M: numner of Rows (Customers+Depots+Parks+etc);  mc: Number of cardinality constraints  
	int nmax;      // nmax: number of allocated columns
	int offset;    // Rows between 0 and offset-1 are not subject to the set-covering constraint.
	long maxnz;
	double *obj;
	double *matval;
	int *matind;
	int *matbeg;
	int *sense;
	double *rhs;
	double *cardrhs;
	double *lb;
	double *ub;
	double *rmatval;
	int *rmatind;
	int *rmatbeg;
	long *val;
	short *ind;
	long *indl;
	long *indl1;
	long *jmin;
	double *Qmin;
	double *CPmin;
	double *objpen;
	int *cardmat;

	// New to manage valid iniequalities (SR3)
	int Bol_SR3 = 0;  // If Bol_SR3 = 1, then SR3 are added
	long nSR3;
	int *rsense;
	double *rrhs;
	int rmax;  // rmax: number of allocated rows used for valid inequalities
	long rmaxnz;

	//OutPut
	FILE *ferr;
	int solstat;
	double objval;
	long *cflag;
	long *rflag;
	double *x;
	double *pen;
	double *pent;
	double *sub;
	double *subo;
	double *q;
	double *u;
	double *ubest;
	double *cred;
	double objdual;

	// Output DFF 
	double ZLBDFF;

	// Output Dual Ascent 
	double ZLBDA;
	double ZUBDA;
	double alphaDA;
	long IterDA;
	long nCoreDA;
	double TimeDA;

	// Output Dual Ascent + Lagrangian Heuristic 
	double ZLBDAH;
	double ZUBDAH;
	double alphaDAH;
	long IterDAH;
	long nCoreDAH;

	// Output Dual Ascent + Lagrangian Heuristic (BolLagrHeu==1)
	double ZLBDAH1;
	double ZUBDAH1;
	double alphaDAH1;
	long IterDAH1;
	long nCoreDAH1;

	// Output Dual Ascent + Lagrangian Heuristic (BolLagrHeu==2)
	double ZLBDAH2;
	double ZUBDAH2;
	double alphaDAH2;
	long IterDAH2;
	long nCoreDAH2;

	// Output Dual Ascent + Lagrangian Heuristic + MIP 
	double ZLBlagr;
	double ZUBlagr;
	double alpha;
	long Iter;
	long nCore;

	// Output LP-relaxation 
	double ZLBlp;
	double ZUBlp;
	double Gap;
	long Nodes;
	long Cuts;

	// Output Greedy
	double ZGreedy;

	Subgradient(void);
	~Subgradient(void);

	// Prototipi delle funzioni impiegate
	//int GenInstance(long type, long nitems, long capbin, long kcostr);
	int AllocateDataStructure(void);
	int DeAllocateDataStructure(void);

	int Dual_Ascent_SetPar(long BolLagrHeu, long *BolOpt);
	int Dual_Ascent_SetCov(int FlagHeu);
	int Solve_Model(long flaglp, double *dual, long *niter);
	int Solve_Incremental_Core(long flaglp, long flaggreedy, double *soldx, long *iter, long *BolOpt);

	long DFF(long i, long alpha, long w, long W);
	long Reductions(void);
	long LBDFF(void);
	long LBDFFC(void);
	long LBDFFC1(void);

	int Close(void);

};