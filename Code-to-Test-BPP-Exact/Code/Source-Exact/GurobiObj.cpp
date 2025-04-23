//-----------------------------------------------------------------------------
//  File: GurobiObj.cpp                                                       
//                                                                           
//  Project: Gurobi solver interface (reduced version)                 
//                                                                           
//  Authors: M.A. Boschetti and S. Novellani                                                  
//                                                                           
//  Last update: 20.04.2025                                                  
//-----------------------------------------------------------------------------
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "GurobiObj.h"

//-----------------------------------------------------------------------------
//  Constructor                                    
//-----------------------------------------------------------------------------
//
GurobiObj::GurobiObj(void)
{
	env = NULL;
	lp = NULL;

	// Setup variables
	bolfree = 0;    // Is memory allocated inside the object?
	minmax = 1;     // Default: minimize
	ncols = 0;      // Number of columns/variables
	nrows = 0;      // Number of rows/constraints
	nz = 0;         // Number of non-zero coefficients into the constraint matrix
	objsen = +1;    // Minimization/Maximization Flag (+1=Min; -1=Max)
	obj = NULL;     // Objective Function coefficients
	rhs = NULL;     // Right Hand Side
	sense = NULL;   // Constraint Senses
	matbeg = NULL;  // Pointers to the beginning of each columns of the constraint matrix
	matcnt = NULL;  // Cardinality of each columns of the constraint matrix
	matind = NULL;  // Row index associated to each coefficient of the constraint matrix
	matval = NULL;  // Value of each coefficient of the constraint matrix
	lb = NULL;      // Lower Bound value of each variable
	ub = NULL;      // Upper Bound value of each variable
	indices = NULL;
	priority = NULL;
	direction = NULL;
	xctype = NULL;
	rngval = NULL;

	rmatbeg = NULL;  // Pointers to the beginning of each columns of the constraint matrix
	rmatcnt = NULL;  // Cardinality of each columns of the constraint matrix
	rmatind = NULL;  // Row index associated to each coefficient of the constraint matrix
	rmatval = NULL;  // Value of each coefficient of the constraint matrix

	nsos = 0;        // Number of SOS (Special Ordered Set) 
	nsosnz = 0;      // Total number of membres in all SOS added
	typesos = NULL;  // Array containing SOS type (CPX_TYPE_SOS1='1', CPX_TYPE_SOS2='2')
	sosbeg = NULL;   // Pointers to the beginning of each SOS
	sosind = NULL;   // Index associated to the SOS member
	soswt = NULL;    // Weight associated to the SOS member

	// Output
	xt = NULL;      // Variable Values

	// Data structure used to add cutting planes
	cutinfo.FlagLP = 0;
	cutinfo.env = env;
	cutinfo.lp  = lp;
	cutinfo.x   = NULL;
	cutinfo.beg = NULL;
	cutinfo.ind = NULL; 
	cutinfo.val = NULL;
	cutinfo.rhs = NULL;

}

//-----------------------------------------------------------------------------
//  Destructor                                    
//-----------------------------------------------------------------------------
//
GurobiObj::~GurobiObj(void)
{
	// Close the Gurobi Enviroment
	if (lp != NULL)
	{
		GRBfreemodel(lp);
		lp = NULL;
	}
	if (env != NULL)
	{
		GRBfreeenv(env);
		env = NULL;
	}

	// Free data structures
	if (bolfree)
	{
		free(obj);
		free(rhs);
		free(sense);
		free(matbeg);
		free(matcnt);
		free(matind);
		free(matval);
		free(lb);
		free(ub);
		free(xt);
		free(indices);
		free(priority);
		free(direction);
		free(xctype);

		free(typesos);
		free(sosbeg);
		free(sosind);
		free(soswt);

		free(cutinfo.x);
		free(cutinfo.beg);
		free(cutinfo.ind); 
		free(cutinfo.val);
		free(cutinfo.rhs);
	}
}

//-----------------------------------------------------------------------------
//  Open: open an environment and open a new problem 
//-----------------------------------------------------------------------------
//
int GurobiObj::Open(void)
{
	int status;

	status = GRBloadenv(&env, NULL);
	if (status || env == NULL) 
	{
		printf("\nError... GRBloadenv...\n");
		env = NULL;
		return status;
	}

	status = GRBsetintparam(env, GRB_INT_PAR_OUTPUTFLAG, 1);
	if (status>0) return status;

	return 0;
}

//-----------------------------------------------------------------------------
//  close: close environment and problem 
//-----------------------------------------------------------------------------
//
int GurobiObj::Close(void)
{
	if (lp != NULL)
	{
		GRBfreemodel(lp);
		lp = NULL;
	}

	if (env != NULL)
	{
		GRBfreeenv(env);
		env = NULL;
	}

	// Free data structures
	if (bolfree)
	{
		FreeCols();
	}

	return 0;
}

//-----------------------------------------------------------------------------
//  Reset: close the existing problem and create a new one 
//-----------------------------------------------------------------------------
//
int GurobiObj::Reset(void)
{
	if (lp != NULL)
	{
		GRBfreemodel(lp);
		lp = NULL;
	}

	FreeCols();

	return 0;
}

//-----------------------------------------------------------------------------
//  MallocCols: Allocate the "column" data structure
//-----------------------------------------------------------------------------
//
int GurobiObj::MallocCols(void)
{
	bolfree = 1;
	obj = (double*)malloc(ncols * sizeof(double));
	rhs = (double*)malloc(nrows * sizeof(double));
	sense = (char*)malloc(nrows * sizeof(char));
	matbeg = (int*)malloc((ncols + 4) * sizeof(int));
	matcnt = (int*)malloc((ncols + 4) * sizeof(int));
	matind = (int*)malloc(nz * sizeof(int));
	matval = (double*)malloc(nz * sizeof(double));
	lb = (double*)malloc(ncols * sizeof(double));
	ub = (double*)malloc(ncols * sizeof(double));
	xt = (double*)malloc(ncols * sizeof(double));
	indices = (int*)malloc((ncols + 4) * sizeof(int));
	priority = (int*)malloc((ncols + 4) * sizeof(int));
	direction = (int*)malloc((ncols + 4) * sizeof(int));
	xctype = (char*)malloc((ncols + 4) * sizeof(char));

	if (nsos > 0)
	{
		typesos = (char*)malloc(nsos * sizeof(char));
		sosbeg = (int*)malloc((nsos + 4) * sizeof(int));
		sosind = (int*)malloc(ncols * sizeof(int));
		soswt = (double*)malloc(ncols * sizeof(double));
	}

	return 0;
}

//-----------------------------------------------------------------------------
//  MallocRows: Allocate the "row" data structure
//-----------------------------------------------------------------------------
//
int GurobiObj::MallocRows(void)
{
	bolfree = 1;
	obj = (double*)malloc(ncols * sizeof(double));
	rhs = (double*)malloc(nrows * sizeof(double));
	sense = (char*)malloc(nrows * sizeof(char));
	rmatbeg = (int*)malloc((nrows + 4) * sizeof(int));
	rmatcnt = (int*)malloc((nrows + 4) * sizeof(int));
	rmatind = (int*)malloc(nz * sizeof(int));
	rmatval = (double*)malloc(nz * sizeof(double));
	lb = (double*)malloc(ncols * sizeof(double));
	ub = (double*)malloc(ncols * sizeof(double));
	xt = (double*)malloc(ncols * sizeof(double));
	indices = (int*)malloc((ncols + 4) * sizeof(int));
	priority = (int*)malloc((ncols + 4) * sizeof(int));
	direction = (int*)malloc((ncols + 4) * sizeof(int));
	xctype = (char*)malloc((ncols + 4) * sizeof(char));

	if (nsos > 0)
	{
		typesos = (char*)malloc(nsos * sizeof(char));
		sosbeg = (int*)malloc((nsos + 4) * sizeof(int));
		sosind = (int*)malloc(ncols * sizeof(int));
		soswt = (double*)malloc(ncols * sizeof(double));
	}

	return 0;
}

//-----------------------------------------------------------------------------
//  FreeCols: Deallocate the "column" data structure
//-----------------------------------------------------------------------------
//
int GurobiObj::FreeCols(void)
{
	if (obj != NULL) free(obj);
	if (rhs != NULL) free(rhs);
	if (sense != NULL) free(sense);
	if (matbeg != NULL) free(matbeg);
	if (matcnt != NULL) free(matcnt);
	if (matind != NULL) free(matind);
	if (matval != NULL) free(matval);
	if (lb != NULL) free(lb);
	if (ub != NULL) free(ub);
	if (xt != NULL) free(xt);
	if (indices != NULL) free(indices);
	if (priority != NULL) free(priority);
	if (direction != NULL) free(direction);
	if (xctype != NULL) free(xctype);

	if (typesos != NULL) free(typesos);
	if (sosbeg != NULL) free(sosbeg);
	if (sosind != NULL) free(sosind);
	if (soswt != NULL) free(soswt);

	if (cutinfo.x != NULL) free(cutinfo.x);
	if (cutinfo.beg != NULL) free(cutinfo.beg);
	if (cutinfo.ind != NULL) free(cutinfo.ind);
	if (cutinfo.val != NULL) free(cutinfo.val);
	if (cutinfo.rhs != NULL) free(cutinfo.rhs);

	obj = NULL;
	rhs = NULL;
	sense = NULL;
	matbeg = NULL;
	matcnt = NULL;
	matind = NULL;
	matval = NULL;
	lb = NULL;
	ub = NULL;
	indices = NULL;
	priority = NULL;
	direction = NULL;
	xctype = NULL;
	rngval = NULL;

	typesos = NULL;
	sosbeg = NULL;
	sosind = NULL;
	soswt = NULL;

	xt = NULL;

	cutinfo.x = NULL;
	cutinfo.beg = NULL;
	cutinfo.ind = NULL;
	cutinfo.val = NULL;
	cutinfo.rhs = NULL;

	return 0;
}

//-----------------------------------------------------------------------------
//  FreeRows: Deallocate the "row" data structure
//-----------------------------------------------------------------------------
//
int GurobiObj::FreeRows(void)
{
	if (obj != NULL) free(obj);
	if (rhs != NULL) free(rhs);
	if (sense != NULL) free(sense);
	if (rmatbeg != NULL) free(rmatbeg);
	if (rmatcnt != NULL) free(rmatcnt);
	if (rmatind != NULL) free(rmatind);
	if (rmatval != NULL) free(rmatval);
	if (lb != NULL) free(lb);
	if (ub != NULL) free(ub);
	if (xt != NULL) free(xt);
	if (indices != NULL) free(indices);
	if (priority != NULL) free(priority);
	if (direction != NULL) free(direction);
	if (xctype != NULL) free(xctype);

	if (typesos != NULL) free(typesos);
	if (sosbeg != NULL) free(sosbeg);
	if (sosind != NULL) free(sosind);
	if (soswt != NULL) free(soswt);

	if (cutinfo.x != NULL) free(cutinfo.x);
	if (cutinfo.beg != NULL) free(cutinfo.beg);
	if (cutinfo.ind != NULL) free(cutinfo.ind);
	if (cutinfo.val != NULL) free(cutinfo.val);
	if (cutinfo.rhs != NULL) free(cutinfo.rhs);

	obj = NULL;
	rhs = NULL;
	sense = NULL;
	rmatbeg = NULL;
	rmatcnt = NULL;
	rmatind = NULL;
	rmatval = NULL;
	lb = NULL;
	ub = NULL;
	indices = NULL;
	priority = NULL;
	direction = NULL;
	xctype = NULL;
	rngval = NULL;

	typesos = NULL;
	sosbeg = NULL;
	sosind = NULL;
	soswt = NULL;

	xt = NULL;

	cutinfo.x = NULL;
	cutinfo.beg = NULL;
	cutinfo.ind = NULL;
	cutinfo.val = NULL;
	cutinfo.rhs = NULL;

	return 0;
}

//-----------------------------------------------------------------------------
//  CopyLPbyCols
//-----------------------------------------------------------------------------
//
int GurobiObj::CopyLPbyCols(void)
{
	status = GRBloadmodel(env, &lp, "lp", ncols, nrows, minmax, 0.0, obj, sense, rhs, matbeg, matcnt, matind, matval, lb, ub, xctype, NULL, NULL);

	return status;
}

//-----------------------------------------------------------------------------
//  CopyLPbyRows
//-----------------------------------------------------------------------------
//
int GurobiObj::CopyLPbyRows(void)
{
	int status;

	status = GRBnewmodel(env, &lp, "lp", ncols, obj, lb, ub, xctype, NULL);
	if (status > 0) return status;

	status = GRBaddconstrs(lp, nrows, nz, rmatbeg, rmatind, rmatval, sense, rhs, NULL);
	if (status > 0) return status;

	status = GRBsetintattr(lp, "ModelSense", minmax);
	if (status > 0) return status;

	status = GRBupdatemodel(lp);
	if (status > 0) return status;

	return status;
}

//-----------------------------------------------------------------------------
//  ReadMPS
//-----------------------------------------------------------------------------
//
int GurobiObj::ReadMPS(const char * fname)
{
	int status;
	int nz1;
	long i;

	// Load the LP instance
	status = GRBreadmodel(env, fname, &lp);
	if (status) return status;

	// Load in memory the instance

	// Sizes, Types, etc.
	status = GRBgetintattr(lp, "NumVars", &ncols);
	if (status) return status;

	status = GRBgetintattr(lp, "NumConstrs", &nrows);
	if (status) return status;

	status = GRBgetintattr(lp, "NumBinVars", &nbin);
	if (status) return status;

	status = GRBgetintattr(lp, "NumIntVars", &nint);
	if (status) return status;

	status = GRBgetintattr(lp, "NumNZs", &nz);
	if (status) return status;

	status = GRBgetintattr(lp, "ModelSense", &minmax);
	if (status) return status;

	status = GRBgetintattr(lp, "IsMIP", &probtype);
	if (status) return status;

	// Allocate memory
	MallocCols();

	// Obj, Cols, etc.
	status = GRBgetdblattrarray(lp, "Obj", 0, ncols, obj);
	if (status) return status;

	status = GRBgetvars(lp, &nz1, matbeg, matind, matval, 0, ncols);
	if (status) return status;
	if (nz1 != nz)
	{
		status = 999;
		return status;
	}
	matbeg[ncols] = nz;
	for (i = 0; i < ncols; i++)
		matcnt[i] = matbeg[i + 1] - matbeg[i];

	status = GRBgetdblattrarray(lp, "LB", 0, ncols, lb);
	if (status) return status;

	status = GRBgetdblattrarray(lp, "UB", 0, ncols, ub);
	if (status) return status;

	status = GRBgetcharattrarray(lp, "VType", 0, ncols, xctype);
	if (status) return status;

	status = GRBgetdblattrarray(lp, "RHS", 0, nrows, rhs);
	if (status) return status;

	status = GRBgetcharattrarray(lp, "Sense", 0, nrows, sense);
	if (status) return status;

	return status;
}

//-----------------------------------------------------------------------------
//  Set a Lower Bound for the MIP
//-----------------------------------------------------------------------------
//
int GurobiObj::SetLB(double LB)
{
	// Set Parameter
	status = GRBsetdblparam(env, GRB_DBL_PAR_CUTOFF, LB);

	return status;
}

//-----------------------------------------------------------------------------
//  SetMIP
//-----------------------------------------------------------------------------
//
int GurobiObj::SetMIP(double TLim)
{
	//int Param;
	int Value;
	double DValue;

	// Set Parameter
	status = GRBsetdblparam(env, GRB_DBL_PAR_TIMELIMIT, TLim);
	if (status) return status;

	DValue = 10e-9;
	status = GRBsetdblparam(env, GRB_DBL_PAR_MIPGAP, DValue);
	if (status) return status;

	DValue = 10e-10;
	status = GRBsetdblparam(env, GRB_DBL_PAR_MIPGAPABS, DValue);
	if (status) return status;

	Value = 0;
	status = GRBsetintparam(env, GRB_INT_PAR_COVERCUTS, Value);
	if (status) return status;

	return status;
}

//-----------------------------------------------------------------------------
//  Add New Rows
//-----------------------------------------------------------------------------
//
int GurobiObj::AddNewRows(void)
{
	int status = 1;

	//// Load new Rows
	////type = xctype - nrange;
	////status = GRBaddvars(lp, nrange, 0, NULL, NULL, NULL, NULL, NULL, NULL, type, NULL);
	////if (status) return status;

	//status = GRBaddconstrs(lp, nrange, k, beg, ind, val, sense, rhs, NULL);
	//if (status) return status;

	//cnt = ncols + nrange;
	//status = GRBsetdblattrarray(lp, "Obj", 0, cnt, objt);

	//status = GRBupdatemodel(lp);
	//if (status) return status;

	return status;
}

//-----------------------------------------------------------------------------
//  SolveMIP
//-----------------------------------------------------------------------------
//
int GurobiObj::SolveMIP(double *Zopt, double *Zlb, double *Gap, long *Nodes, long *Cuts)
{
	int ncuts;
	int solstat;
	double RGap;
	double aux;

	status = GRBoptimize(lp);
	if (status) return status;

	status = GRBgetintattr(lp, "Status", &solstat);
	if (status) return status;

	if ((solstat != GRB_OPTIMAL) && (status != GRB_NODE_LIMIT) && (solstat != GRB_TIME_LIMIT))
	{
		printf("Non Ammissibile... %d\n", status);
		objval = +0.0;
		RGap = 1.0;
		goto TERMINATE;
	}

	// Read the solution
	status = GRBgetdblattr(lp, "ObjVal", &objval);
	if (status) return status;

	status = GRBgetdblattr(lp, "ObjBound", Zlb);
	if (status) return status;
	RGap = minmax * (objval - (*Zlb)) / objval;

	status = GRBgetdblattrarray(lp, "X", 0, ncols, xt);
	if (status) return status;

TERMINATE:

	// Ouput
	*Zopt = objval;
	*Gap = RGap * 100.;

	status = GRBgetdblattr(lp, "NodeCount", &aux);
	if (status) return status;
	*Nodes = (long)aux;

	//status = CPXgetnumcuts(env,lp,CPX_CUT_USER,&ncuts);
	//status = CPXgetnumcuts(env,lp,CPX_CUT_COVER,&ncuts);
	ncuts = 0;
	*Cuts = ncuts;

	return status;
}

//-----------------------------------------------------------------------------
//  SolveLP
//-----------------------------------------------------------------------------
//
int GurobiObj::SolveLP(double *Zopt)
{
	int solstat;
	double RGap;

	status = GRBoptimize(lp);
	if (status) return status;

	status = GRBgetintattr(lp, "Status", &solstat);
	if (status) return status;

	if ((solstat != GRB_OPTIMAL) && (status != GRB_NODE_LIMIT) && (solstat != GRB_TIME_LIMIT))
	{
		printf("Non Ammissibile... %d\n", status);
		objval = +0.0;
		RGap = 1.0;
		goto TERMINATE;
	}

	// Read the solution
	status = GRBgetdblattr(lp, "ObjVal", &objval);
	if (status) return status;

	status = GRBgetdblattrarray(lp, "X", 0, ncols, xt);
	if (status) return status;

TERMINATE:

	// Ouput
	*Zopt = objval;

	return status;
}

//-----------------------------------------------------------------------------
//  Setup Custom Cutting Plane
//-----------------------------------------------------------------------------
//
int GurobiObj::CuttingPlane(long n, long m, double* u, double* p, double* pmax,
	long *nY, long **Y, long *nUOR, long **UOR, long *nUAND, long **UAND, long *FDepA, long *FDepP,
	int Pred, int Cover, int Lifting, int KnapSol, int MaxCuts, long *Cuts)
{
	int aux;
	long cur_numcols;

	// Set parameters

	//// Assure linear mappings between the presolved and original models 
	//status = CPXsetintparam (env, CPX_PARAM_PRELINEAR, 0);
	//if ( status )  goto TERMINATE;
 //  
	//// Turn on traditional search for use with control callbacks 
	//status = CPXsetintparam (env, CPX_PARAM_MIPSEARCH, CPX_MIPSEARCH_TRADITIONAL);
	//if ( status )  goto TERMINATE;

	//// Let MIP callbacks work on the original model
	//status = CPXsetintparam (env, CPX_PARAM_MIPCBREDLP, CPX_OFF);
	//if ( status )  goto TERMINATE;

	//// Let MIP do not reduce LP
	//status = CPXsetintparam (env, CPX_PARAM_REDUCE, CPX_PREREDUCE_NOPRIMALORDUAL);
	//if ( status )  goto TERMINATE;

	//// Let MIP do not presolve MIP
	////status = CPXsetintparam (env, CPX_PARAM_PREIND, CPX_OFF);
	////if ( status )  goto TERMINATE;

	*Cuts = 0;
	status = GRBgetintattr(lp, "NumVars", &aux);
	if (status) return status;
	cur_numcols = aux;
	cutinfo.FlagLP = 0;
	cutinfo.lp = lp;
	cutinfo.numcols = cur_numcols;
	cutinfo.numtoadd = 0;
	cutinfo.maxcuts = MaxCuts;

	cutinfo.x = (double *)malloc(cur_numcols * sizeof(double));
	if (cutinfo.x == NULL) {
		fprintf(stderr, "No memory for solution values.\n");
		goto TERMINATE;
	}
	cutinfo.ind = (int *)malloc(cur_numcols * sizeof(int));
	if (cutinfo.ind == NULL) {
		fprintf(stderr, "No memory for solution values.\n");
		goto TERMINATE;
	}
	cutinfo.val = (double *)malloc(cur_numcols * sizeof(double));
	if (cutinfo.val == NULL) {
		fprintf(stderr, "No memory for solution values.\n");
		goto TERMINATE;
	}
	cutinfo.beg = (int *)malloc(2 * sizeof(int));
	if (cutinfo.beg == NULL) {
		fprintf(stderr, "No memory for solution values.\n");
		goto TERMINATE;
	}
	cutinfo.rhs = (double *)malloc(2 * sizeof(double));
	if (cutinfo.rhs == NULL) {
		fprintf(stderr, "No memory for solution values.\n");
		goto TERMINATE;
	}

	// Set up to use MIP callback 
	status = GRBsetcallbackfunc(lp, myNewcutcallback, &cutinfo);
	if (status)  goto TERMINATE;

TERMINATE:

	return 0;
}

//-----------------------------------------------------------------------------
// Function implementing a custom cutting plane procedure called by callback
//-----------------------------------------------------------------------------
//
int __stdcall myNewcutcallback(GRBmodel *lp,
	void       *cbdata,
	int        wherefrom,
	void       *cbhandle)
{
	long i;
	int status = 0;
	//int BolAdd;
	//char sense;
	//int purgeable;
	//int param;
	//int ivalue;

	GRBCUTINFOptr cutinfo = (GRBCUTINFOptr)cbhandle;

	//*useraction_p = CPX_CALLBACK_DEFAULT; 

	// Are we solving a subproblem?
	if (cutinfo->FlagLP)
		return (status);

	// Do we have reached the maximum number of cuts?  
	if (cutinfo->numtoadd > cutinfo->maxcuts)
		return (status);

	if (wherefrom == GRB_CB_MIPSOL_SOL)
	{
		status = GRBcbget(cbdata, wherefrom, GRB_CB_MIPSOL_SOL, cutinfo->x);
		if (status)
		{
			fprintf(stderr, "Failed to get node solution.\n");
			goto TERMINATE;
		}
	}
	else
		return 0;

	printf("\nSolution:\n");
	for (i = 0; i < cutinfo->numcols; i++)
		if (cutinfo->x[i] > 0.00001)
			printf("x(%d)=%lf\n", i, cutinfo->x[i]);

	// Use a cut violation tolerance of 0.01
	//if (BolAdd==1) 
	//{
	//	cutinfo->numtoadd++;
	//	cutinfo->Cuts[0]++;
	//	if ((cutinfo->numtoadd%1000)==0) printf("\n  User Cuts added: %d\n\n",cutinfo->numtoadd);

	//	//if ((cutinfo->numtoadd%5000)==0) 
	//	//	printf("\n  User Cuts added: %d\n\n",cutinfo->numtoadd);

	//	// status = CPXcutcallbackadd(env,cbdata,wherefrom,k1,rhs[0],'G',ind,val,0);
	//	status = CPXcutcallbackadd(env,cbdata,wherefrom,cutinfo->nz,cutinfo->rhs[0],sense,cutinfo->ind,cutinfo->val,purgeable);
	//	if (status)
	//	{
	//		fprintf (stderr, "Failed to add cut.\n");
	//		goto TERMINATE;
	//	}

	//	// Tell CPLEX that cuts have been created
	//	*useraction_p = CPX_CALLBACK_SET; 

	//}

TERMINATE:

	return (status);

}




