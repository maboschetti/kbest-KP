//-----------------------------------------------------------------------------
//  File: CplexObj.h                                                       
//                                                                           
//  Project: Cplex solver interface (header)                   
//                                                                           
//  Author: M.A. Boschetti and S. Novellani                                                   
//                                                                           
//  Last update: 21.04.2025                                                  
//-----------------------------------------------------------------------------
//
#include "cplex.h"

struct tupla {
	long j1, j2, j3;
};

struct cutinfo {
	int FlagLP;
	CPXENVptr env;
	CPXLPptr lp;
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

typedef struct cutinfo CUTINFO, *CUTINFOptr;

// Function implementing a custom cutting plane procedure called by callback
static int CPXPUBLIC mycutcallback (CPXCENVptr env, void *cbdata, int wherefrom,
                void *cbhandle, int *useraction_p);
static int CPXPUBLIC myNewcutcallback (CPXCENVptr env, void *cbdata, int wherefrom,
                void *cbhandle, int *useraction_p);

class CplexObj
{
private:

	int status;
	int bolfree;  // Is memory allocated inside the object?

public:

	CPXENVptr env;   // Pointer to Cplex Environment 
	CPXLPptr lp;     // Pointer to LP instance
	CUTINFO cutinfo; // Data structure used to add cutting planes
	
	// Input
	int bolmip;      // Boolean (0 = LP; 1 = MIP) 
	int minmax;      // Mininum or Maximum (1=minimize; -1=maximize) 
	int ncols;       // Number of columns/variables
	int nrows;       // Number of rows/constraints
	int nz;          // Number of non-zero coefficients into the constraint matrix
	int nbin;        // Number of binary variables
	int nint;        // Number of integer variables
	int objsen;      // Minimization/Maximization Flag (+1=Min; -1=Max)
	int probtype;    // Problem Type (-1: err; CPXPROB_LP=0; CPXPROB_MILP=1; etc.)
	int lpsolver;    // LP Solver (0: automatic; 1: Primal; 2: Dual; 3: Net; 4: Barrier)
	double *obj;     // Objective Function coefficients
	double *rhs;     // Right Hand Side
	char *sense;     // Constraint Senses
	int *matbeg;     // Pointers to the beginning of each columns of the constraint matrix
	int *matcnt;     // Cardinality of each columns of the constraint matrix
	int *matind;     // Row index associated to each coefficient of the constraint matrix
	double *matval;  // Value of each coefficient of the constraint matrix
	int rncols;      // Number of new columns/variables in the new rows
	int rnrows;      // Number of new rows/constraints
	int rnz;         // Number of non-zero coefficients into the new rows
	char *rsense;     // Constraint sense for the new rows
	double *rrhs;    // Right Hand Side for the new rows
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
	CplexObj(void);
	~CplexObj(void);

	// Open: open an environment and open a new problem 
	int Open(void);

	// Close: close environment and problem 
	int Close(void);

	// Reset: close the existing problem and create a new one 
	int Reset(void);

	// MallocCols: Allocate the "column" data structure
	int MallocCols(void);

	// MallocRows: Allocate the "row" data structure
	int MallocRows(void);

	// FreeCols: Deallocate the "column" data structure
	int FreeCols(void);

	// FreeRows: Deallocate the "row" data structure
	int FreeRows(void);

	// ReadMPS
	int ReadMPS(const char * fname);

	// CopyLP
	int CopyLPbyCols(void);
	int CopyLPbyRows(void);

    // Add New Rows
	int AddNewRows(void);

	// Set a Lower Bound for the MIP
	int SetLB(double LB);

	// SetMIP
	int SetMIP(double TLim);

	// SolveMIP
	int SolveMIP(double *Zopt, double *Zlb, double *Gap, long *Nodes, long *Cuts);

	// SolveLP
	int SolveLP(double *Zopt, double *dual, long *niter);

    // Setup Custom Cutting Plane
	int CuttingPlane(long n, long m, double* u, double* p, double* pmax, 
		             long *nY, long **Y, long *nUOR, long **UOR, long *nUAND, long **UAND, long *FDepA, long *FDepP, 
					 int Pred, int Cover, int Lifting, int KnapSol, int MaxCuts, long *Cuts);

};