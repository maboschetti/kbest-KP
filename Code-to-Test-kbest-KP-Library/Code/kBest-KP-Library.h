//-----------------------------------------------------------------------------
// Here we must insert our "license boilerplate"
// (we need to investigate what is the most suitable licence for us
//-----------------------------------------------------------------------------
// File: kBest-KP-Library.h
// Version 1.0
// Last update: 01.12.2022
// Authors: Boschetti M.A., Novellani, S.
//-----------------------------------------------------------------------------
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "Utilities.h"

// Infinity
#define  Infi  100000
#define  Infr  100000.0

// Minimum and maximum size increments for dynamic memory allocation
#define  PercMinDim  100
#define  MinDim      5
#define  PercMaxDim  10
#define  MaxDim      10

struct Instance
{
	long n;    // Number of items
	long W;    // Knapsack capacity
	long *w;   // Item weights (i.e., w[i] is the weight of item i)
	long *p;   // Item profits (i.e., p[i] is the profit of item i)
	float *pf; // Item profits float (i.e., p[i] is the profit of item i)
	long *b;   // Item bounds (i.e., b[i] id the maximum number of copies of item i)
	long k;    // Number of best solutions to find (i.e., we want find the k-best solutions) 
};

struct Solution
{
	long z;  // Optimal solution value
	long *x; // Optimal solution: x[i]=0 item i is not in solution; x[i]=1 item i is in solution
};

struct DSolution
{
	double z;  // Optimal solution value
	double *x; // Optimal solution: x[i]=0 item i is not in solution; x[i]=1 item i is in solution
};

struct kBSolutions
{
	long nk;   // Number of solutions found, 0 <= nk <= k.
	long *z;   // k-best solution values
	float *zf; // k-best solution values
	long **x;  // k-best solutions: x[k][i]=0 item i is not in the k-best solution; x[k][i]=1 item i is in k-best solution
	long **id; // k-best solutions: id[k][i] is the index of the item in the k-best solution in position i
	long *nx;  // k-best solutions: nx[k] is the number of items in the k-best solution
};

struct Constraints
{
	long Type;     // Type of contraint: 0=No Constraint; 1=Cardinality Constraints; etc.
	long MaxCard;  // Maximum cardinality of each cutting pattern
	// Etc. (Other parameters)
	// We may also consider the possibility to have a flag for each kind of contraints. 
	// This option allows combinations of different constraints
};

struct KPIs
{
	float Avg_k;       // Average k for each active state (j,w)
	float Avg_Dim_k;   // Average dimension of the vector allocated to save the k-best state (j,w)
	float Perc_Dim_k; // Percentage of memory allocated to save the k-best state (j,w)
	long Num_States_Used;    // Number of states (j,w) used (with k>0)
	float Perc_States_Used;  // Percentage of states (j,w) used (with k>0)
	long Num_States_Generated;    // Number of states (j,w,k) generated 
	float Perc_States_Generated;  // Percentage of states (j,w,k) generated
};

class KPSolver
{
	// Private data structures for computing the LP-relaxation solution
	double *ratio;
	long *ind;

	// Private Variables fo computing the optimal solution
	long **Zkp1;   // Zkp(j,w) is the cost of state w of stage j
	long **Pred1;  // Pred(j,w) is the precedessor of state w of stage j
	long **Kard;   // Kard(j,w) indicates the number of generated states (it is k(j,w) reported in the paper) 
	long **Kard1;  // Kard1(j,w) indicates the number of generated states (it is k1(j,w) reported in the paper) 
	long **Kard2;  // Kard2(j,w) indicates the number of generated states (it is k2(j,w) reported in the paper) 
	long wmin;  // Minimum item weight   

	// Private Variables for computing the k-best solutions
	long ***Zkp;   // Zkp(j,w,k) is the cost of state (w,k) of stage j
	long ***Pred;  // Pred(j,w,k) is the precedessor of state (w,k) of stage j

	// Private Variables for computing the k-best solutions
	long *ZkpL;     // Zkp(j,w,k) is the cost of state (w,k) of stage j
	long *PredL;    // Pred(j,w,k) is the precedessor of state (w,k) of stage j
	long **ZkpLD;   // Zkp(j,w,k) is the cost of state (w,k) of stage j
	float **ZkpLDF; // Zkp(j,w,k) is the cost of state (w,k) of stage j (the use of float is usefull for column generation)
	long **PredLD;  // Pred(j,w,k) is the precedessor of state (w,k) of stage j
	long **PredALD; // Preda(j,w,k) is the number of copies of item j in solution at state (w,k) of stage j
	long *Dim;      // Dimension of the array allocated for each state w of each stage j
	long *KardL;    // Kard(j,w) indicates the number of generated states (it is k(j,w) reported in the paper) 
	long *KardL1;   // Kard1(j,w) indicates the number of generated states (it is k1(j,w) reported in the paper) 
	long *KardL2;   // Kard2(j,w) indicates the number of generated states (it is k2(j,w) reported in the paper) 
	long **KardB;   // KardB(j,w,a) indicates the number of generated states (it is k_a(j,w) reported in the paper) 
	long *Bound;    // Bound(j,w) indicates the maximum number of copies of item j at state (j,w) 
	long nwk, nk, nw;  // nwk = Inst.W*Inst.k, nk=Inst.k, and nw=Inst.W
	long nnw, nn;      // nnw = Inst.n*Inst.w, nk=Inst.k, and nw=Inst.W (k,j,w)
	long dmax;         // Maximum size for state k
	long ddmin, ddmax; // Minimum and maximum size increment for state k

	// Private functions
	int Generate_State(long j, long w);
	int Generate_State_L(long j, long w);
	void Setup_KPBW(void);
	void Setup_KPBW(long j, long w);
	long KPBW(long j, long w);
	void Setup_KPBW_kBest(void);
	void Setup_KPBW_kBest(long j, long w);
	long KPBW_kBest(long j, long w, long k);
	size_t MapFW(long j, long w);
	size_t MapFW(long j, long w, long k);
	size_t MapBW(long j, long w);
	size_t MapBW(long j, long w, long k);
	void Setup_KPBW_kBest_L(void);
	long KPBW_kBest_L(long j, long w, long k);
	void Setup_KPBW_kBest_LD(void);
	long KPBW_kBest_LD(long j, long w, long k);
	int BuildSolution_kBest_LD(long k, long k1);
	void Setup_KPBW_kBest_LDC(void);
	int BuildSolution_kBest_LDC(long k, long k1);
	int FeasSol_kBest_LDC(long k);
	float KPBW_kBest_LDCF(long j, long w, long k, int *status);

	int Setup_KPBW_Bound_kBest_LDCF(void);
	int BuildSolution_Bound_kBest_LDCF(long k, long k1);
	int FeasSol_Bound_kBest_LDCF(long k);
	float KPBW_Bound_kBest_LDCF(long j, long w, long k, int *status);

public:

	// Public Variables
	struct Instance Inst;
	struct Solution OptSol;
	struct DSolution LPSol;
	struct kBSolutions kBSol;
	struct Constraints Cons;
	struct KPIs KPI;

	int PrintOn = 0;

	// Constructor and Destructor
	KPSolver(void);
	~KPSolver(void);

	// Member Functions
	int Allocate_Sol(void);
	int Allocate_LPKP(void);
	int Allocate_KPFW(void);
	int Allocate_KPBW(void);
	int Allocate_KPFW_kBest(void);
	int Allocate_KPFW_kBest_L(void);
	int Allocate_KPBW_kBest(void);
	int Allocate_KPBW_kBest_L(void);
	int Allocate_KPBW_kBest_LD(void);
	int Allocate_KPBW_kBest_LDC(void);
	int Allocate_KPBW_kBest_LDCF(void);
	int Allocate_KPBW_Bound_kBest_LDCF(void);
	int MAllocate_KPBW_kBest(void);

	int DeAllocate_Sol(void);
	int DeAllocate_LPKP(void);
	int DeAllocate_KPFW(void);
	int DeAllocate_KPBW(void);
	int DeAllocate_KPFW_kBest(void);
	int DeAllocate_KPFW_kBest_L(void);
	int DeAllocate_KPBW_kBest(void);
	int DeAllocate_KPBW_kBest_L(void);
	int DeAllocate_KPBW_kBest_LD(void);
	int DeAllocate_KPBW_kBest_LDC(void);
	int DeAllocate_KPBW_kBest_LDCF(void);
	int DeAllocate_KPBW_Bound_kBest_LDCF(void);
	int DeMAllocate_KPBW_kBest(void);

	int Solve_LPKP(void);
	int Solve_KPFW(void);
	int Solve_KPBW(void);
	int Solve_KPFW_kBest(void);
	int Solve_KPFW_kBest_L(void);
	int Solve_KPBW_kBest(void);
	int Solve_KPBW_kBest_L(void);
	int Solve_KPBW_kBest_LD(void);
	int Solve_KPBW_kBest_LDC(void);
	int Solve_KPBW_kBest_LDCF(void);
	int Solve_KPBW_Bound_kBest_LDCF(void);

	int KPI_KPBW_kBest_L(void);
	int KPI_KPBW_kBest_LDCF(void);
	int KPI_KPBW_Bound_kBest_LDCF(void);
};
