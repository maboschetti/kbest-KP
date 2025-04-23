//-----------------------------------------------------------------------------
//  File: CplexObj.cpp                                                       
//                                                                           
//  Project: Cplex solver interface (reduced version)         
//                                                                           
//  Author: M.A. Boschetti and S. Novellani                                                   
//                                                                           
//  Last update: 20.04.2025                                                  
//-----------------------------------------------------------------------------
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "CplexObj.h"

//-----------------------------------------------------------------------------
//  Constructor                                    
//-----------------------------------------------------------------------------
//
CplexObj::CplexObj(void)
{
	env = NULL;
	lp = NULL;

	bolfree = 0;    // Is memory allocated inside the object?
	bolmip = 0;     // Deafault: LP
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
	xt = NULL;       // Variable Values

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
CplexObj::~CplexObj(void)
{
	//long i,j;

	if (lp != NULL)
	{
		status = CPXfreeprob(env, &lp);
		if (status)
		{
			printf("\nError... CPXfreeprob...\n");
			return;
		}
	}

	if (env != NULL)
	{
		status = CPXcloseCPLEX(&env);
		if (status)
		{
			printf("\nError... CPXcloseCPLEX...\n");
			return;
		}
	}

	if (bolfree)
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
	}
}

//-----------------------------------------------------------------------------
//  Open: open an environment and open a new problem 
//-----------------------------------------------------------------------------
//
int CplexObj::Open(void)
{
	int param;
	int ivalue;

	env = CPXopenCPLEX(&status);
	if (status)
	{
		printf("\nError... CPXopenCPLEX...\n");
		env = NULL;
		return status;
	}

	lp = CPXcreateprob(env,&status,"lp");
	if (status)
	{
		printf("\nError... CPXcreateprob...\n");
		lp = NULL;
		return status;
	}

	param=CPX_PARAM_DATACHECK;
	ivalue=1;
	//ivalue=0;
	status = CPXsetintparam (env, param, ivalue);
	if (status>0) return status;

	param=CPX_PARAM_SCRIND;
	ivalue=1;
	// ivalue=0;
	status = CPXsetintparam (env, param, ivalue);
	if (status>0) return status;

	param=CPX_PARAM_MIPDISPLAY;
	ivalue=2;
	// ivalue=0;
	status = CPXsetintparam (env, param, ivalue);
	if (status>0) return status;

	param=CPX_PARAM_SIMDISPLAY;
	//ivalue=2;
	ivalue=0;
	status = CPXsetintparam (env, param, ivalue);
	if (status>0) return status;

	return 0;
}

//-----------------------------------------------------------------------------
//  Close: close the existing problem 
//-----------------------------------------------------------------------------
//
int CplexObj::Close(void)
{
	status = CPXfreeprob(env, &lp);
	if (status)
	{
		printf("\nError... CPXfreeprob...\n");
		return 1;
	}

	lp = NULL;

	if (env != NULL)
	{
		status = CPXcloseCPLEX(&env);
		if (status)
		{
			printf("\nError... CPXcloseCPLEX...\n");
			return 1;
		}
	}
	env = NULL;

	if (bolfree)
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
	}

	return 0;
}

//-----------------------------------------------------------------------------
//  Reset: close the existing problem and create a new one 
//-----------------------------------------------------------------------------
//
int CplexObj::Reset(void)
{
	status = CPXfreeprob(env, &lp);
	if (status)
	{
		printf("\nError... CPXfreeprob...\n");
		return 1;
	}

	lp = CPXcreateprob(env,&status,"lp");
	if (status)
	{
		printf("\nError... CPXcreateprob...\n");
		lp = NULL;
		return 2;
	}

	FreeCols();

	return 0;
}

//-----------------------------------------------------------------------------
//  MallocCols: Allocate the "column" data structure
//-----------------------------------------------------------------------------
//
int CplexObj::MallocCols(void)
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
int CplexObj::MallocRows(void)
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
int CplexObj::FreeCols(void)
{
	//long i;

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
int CplexObj::FreeRows(void)
{
	//long i;

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
int CplexObj::CopyLPbyCols(void)
{
	long col;

	// Load the LP instance
	status = CPXcopylp(env, lp, ncols, nrows, objsen, obj, rhs, sense, matbeg, matcnt, matind, matval, lb, ub, rngval);

	// Set problem: minimization or maximization
	CPXchgobjsen(env, lp, objsen);

	// Is a MIP?
	bolmip = 0;
	if (xctype != NULL)
	{
		for (col = 0; col < ncols; col++)
			if (xctype[col] != 'C')
			{
				bolmip = 1;
				break;
			}
	}
	if (bolmip)
	{
		status = CPXchgprobtype(env, lp, CPXPROB_MILP);
		if (status) return status;

		status = CPXchgctype(env, lp, ncols, indices, xctype);
		if (status) return status;
	}

	return status;
}

//-----------------------------------------------------------------------------
//  CopyLPbyRows
//-----------------------------------------------------------------------------
//
int CplexObj::CopyLPbyRows(void)
{
	long col;

	// Load the LP instance
	status = CPXnewcols(env, lp, ncols, obj, lb, ub, NULL, NULL);
	if (status) return status;

	status = CPXaddrows(env, lp, 0, nrows, nz, rhs, sense, rmatbeg, rmatind, rmatval, NULL, NULL);
	if (status) return status;

	CPXchgobjsen(env, lp, objsen);

	// Is a MIP?
	bolmip = 0;
	if (xctype != NULL)
	{
		for (col = 0; col < ncols; col++)
			if (xctype[col] != 'C')
			{
				bolmip = 1;
				break;
			}
	}
	if (bolmip)
	{
		status = CPXchgprobtype(env, lp, CPXPROB_MILP);
		if (status) return status;

		status = CPXchgctype(env, lp, ncols, indices, xctype);
		if (status) return status;
	}

	return status;
}

//-----------------------------------------------------------------------------
//  ReadMPS
//-----------------------------------------------------------------------------
//
int CplexObj::ReadMPS(const char * fname)
{
	int status;
	int nz1, surplus;
	long i;

	// Load the LP instance
	status = CPXreadcopyprob(env, lp, fname, NULL);
	if (status) return status;

	// Load in memory the instance

	// Sizes, Types, etc.
	ncols = CPXgetnumcols(env, lp);
	nrows = CPXgetnumrows(env, lp);
	nbin = CPXgetnumbin(env, lp);
	nint = CPXgetnumint(env, lp);
	nz = CPXgetnumnz(env, lp);
	objsen = CPXgetobjsen(env, lp);  // CPX_MIN=1; CPX_MAX=-1
	probtype = CPXgetprobtype(env, lp);  // CPXPROB_LP=0; CPXPROB_MILP=1

	// Allocate memory
	MallocCols();

	// Obj, Cols, etc.
	status = CPXgetobj(env, lp, obj, 0, ncols - 1);
	if (status) return status;

	status = CPXgetcols(env, lp, &nz1, matbeg, matind, matval, nz, &surplus, 0, ncols - 1);
	if (status) return status;
	if ((nz1 != nz) || (surplus < 0))
	{
		status = 999;
		return status;
	}
	matbeg[ncols] = nz;
	for (i = 0; i < ncols; i++)
	{
		indices[i] = i;
		matcnt[i] = matbeg[i + 1] - matbeg[i];
	}

	status = CPXgetlb(env, lp, lb, 0, ncols - 1);
	if (status) return status;

	status = CPXgetub(env, lp, ub, 0, ncols - 1);
	if (status) return status;

	status = CPXgetctype(env, lp, xctype, 0, ncols - 1);
	if (status) return status;

	status = CPXgetrhs(env, lp, rhs, 0, nrows - 1);
	if (status) return status;

	status = CPXgetsense(env, lp, sense, 0, nrows - 1);
	if (status) return status;

	return status;
}

//-----------------------------------------------------------------------------
//  Set a Lower Bound for the MIP
//-----------------------------------------------------------------------------
//
int CplexObj::SetLB(double LB)
{
	int Param;
	double DValue;

	// Set Parameter
	Param = CPX_PARAM_CUTLO;
	DValue = LB - 1.;
	status = CPXsetdblparam(env, Param, DValue);

	return 0;
}

//-----------------------------------------------------------------------------
//  SetMIP
//-----------------------------------------------------------------------------
//
int CplexObj::SetMIP(double TLim)
{
	int Param;
	int Value;
	double DValue;

	// Set Parameter
	Param = CPX_PARAM_TILIM;
	DValue = TLim;
	status = CPXsetdblparam(env, Param, DValue);

	Param = CPX_PARAM_EPAGAP;
	DValue = 10e-9;
	status = CPXsetdblparam(env, Param, DValue);

	Param = CPX_PARAM_EPGAP;
	DValue = 10e-9;
	status = CPXsetdblparam(env, Param, DValue);

	Param = CPX_PARAM_COVERS;
	Value = 0;
	status = CPXsetintparam(env, Param, Value);

	return status;
}

//-----------------------------------------------------------------------------
//  Add New Rows
//-----------------------------------------------------------------------------
//
int CplexObj::AddNewRows(void)
{
	// Load the new rows
	status = CPXaddrows(env, lp, 0, rnrows, rnz, rrhs, rsense, rmatbeg, rmatind, rmatval, NULL, NULL);

	return status;
}

//-----------------------------------------------------------------------------
//  SolveMIP
//-----------------------------------------------------------------------------
//
int CplexObj::SolveMIP(double *Zopt, double *Zlb, double *Gap, long *Nodes, long *Cuts)
{
	int ncuts;
	double RGap;
	int Param;
	int Value;
	double DValue;

	// Set Parameter: number of threads
	Param = CPX_PARAM_THREADS;
	Value = 1; // 0 automatic; n: n. threads
	status = CPXsetintparam(env, Param, Value);

	// Set Parameter: Sets the upper cutoff tolerance
	Param = CPX_PARAM_CUTUP;
	DValue = (*Zopt) + 0.5;
	status = CPXsetdblparam(env, Param, DValue);

	status = CPXmipopt(env, lp);
	if (status)
	{
		printf("Cplex fail... %d\n", status);
		objval = +10e+15;
		RGap = 1.0;
		goto TERMINATE;
	}

	status = CPXgetstat(env, lp);
	//if ((status!=CPXMIP_OPTIMAL)&&(status!=CPXMIP_OPTIMAL_TOL)&&(status!=CPXMIP_TIME_LIM_FEAS)) // Cplex 11.2
	if ((status != CPXMIP_OPTIMAL) && (status != CPXMIP_OPTIMAL_TOL) && (status != CPXMIP_FEASIBLE) && (status != CPXMIP_TIME_LIM_FEAS)) // Cplex 12.2
	{
		printf("Not feasible... %d\n", status);
		objval = +10e+15;
		RGap = 1.0;
		goto TERMINATE;
	}

	// Read the solution
	status = CPXgetobjval(env, lp, &objval);
	status = CPXgetx(env, lp, xt, 0, ncols - 1);
	status = CPXgetmiprelgap(env, lp, &RGap);
	status = CPXgetbestobjval(env, lp, Zlb);

TERMINATE:

	// Ouput
	*Zopt = objval;
	*Gap = RGap * 100.;
	*Nodes = CPXgetnodecnt(env, lp);
	status = CPXgetnumcuts(env, lp, CPX_CUT_USER, &ncuts);
	//status = CPXgetnumcuts(env,lp,CPX_CUT_COVER,&ncuts);
	*Cuts = ncuts;

	return status;
}

//-----------------------------------------------------------------------------
//  SolveLP
//-----------------------------------------------------------------------------
//
int CplexObj::SolveLP(double *Zopt, double *dual, long *niter)
{
	//int ncuts;
	double RGap;
	int Param;
	int Value;
	//double DValue;

	// Set Parameter: number of threads
	Param = CPX_PARAM_THREADS;
	Value = 1;  // 0 automatic; n: n. threads
	status = CPXsetintparam(env, Param, Value);

	// Set Parameter: Screen
	Param = CPX_PARAM_SCRIND;
	Value = 0;  // 0 Turn off display of messages to screen
	status = CPXsetintparam(env, Param, Value);

	// Set Parameter: Screen
	Param = CPX_PARAM_TUNINGDISPLAY;
	Value = 0;  // 0 Turn off display
	status = CPXsetintparam(env, Param, Value);

	// Set Parameter: Screen Barrier
	Param = CPX_PARAM_BARDISPLAY;
	Value = 0;  // 0 No progress information
	status = CPXsetintparam(env, Param, Value);

	// Set Parameter: Screen Barrier
	Param = CPX_PARAM_BARALG;
	Value = 0;  // 0=default, 1=infeasibility-estimate start, 2=Infeasibility-constant start, 3=Standard Barrier
	status = CPXsetintparam(env, Param, Value);

	//// Set Parameter: optimality tolerance
	//Param = CPX_PARAM_EPOPT;
	//DValue = 0.00000001;  // 10e-1 - 10e-9
	//status = CPXsetdblparam(env, Param, DValue);

	//// Set Parameter: feasibility tolerance
	//Param = CPX_PARAM_EPRHS;
	//DValue = 0.00000001;  // 10e-1 - 10e-9
	//status = CPXsetdblparam(env, Param, DValue);

	//// Set Parameter: Emphasis_Numerical
	//Param = CPXPARAM_Emphasis_Numerical;
	//Value = 1;  // 0=OFF, 1=ON
	//status = CPXsetintparam(env, Param, Value);

	if (dual != NULL)
	{
		// Set the dual variables using the values contained in "rdual"
		status = CPXcopystart(env, lp, NULL, NULL, NULL, NULL, NULL, dual);
	}

	if (lpsolver == 0)
		status = CPXlpopt(env, lp);
	else if (lpsolver == 1)
		status = CPXprimopt(env, lp);
	else if (lpsolver == 2)
		status = CPXdualopt(env, lp);
	else if (lpsolver == 3)
		status = CPXhybnetopt(env, lp, CPX_ALG_DUAL);
	else if (lpsolver == 4)
		status = CPXbaropt(env, lp);
	else if (lpsolver == 5)
		status = CPXsiftopt(env, lp);

	status = CPXgetstat(env, lp);
	//if ((status!=CPXMIP_OPTIMAL)&&(status!=CPXMIP_OPTIMAL_TOL)&&(status!=CPXMIP_TIME_LIM_FEAS)) // Cplex 11.2
	if ((status != CPX_STAT_OPTIMAL)) // Cplex 12.2
	{
		printf("Non Ammissibile... %d\n", status);
		objval = +0.0;
		RGap = 1.0;
		goto TERMINATE;
	}

	// Read the solution
	status = CPXgetobjval(env, lp, &objval);
	status = CPXgetx(env, lp, xt, 0, ncols - 1);
	if (dual != NULL)
		status = CPXgetpi(env, lp, dual, 0, nrows - 1);

	// Iteration
	if (lpsolver == 4)
	{
		*niter = CPXgetbaritcnt(env, lp);
		printf("\n N.Iter = %d \n", (*niter));
	}

TERMINATE:

	// Ouput
	*Zopt = objval;

	return status;
}

//-----------------------------------------------------------------------------
//  Setup Custom Cutting Plane
//-----------------------------------------------------------------------------
//
int CplexObj::CuttingPlane(long n, long m, double* u, double* p, double* pmax,
	long *nY, long **Y, long *nUOR, long **UOR, long *nUAND, long **UAND, long *FDepA, long *FDepP,
	int Pred, int Cover, int Lifting, int KnapSol, int MaxCuts, long *Cuts)
{
	long cur_numcols;

	// Set parameters

	// Assure linear mappings between the presolved and original models 
	status = CPXsetintparam(env, CPX_PARAM_PRELINEAR, 0);
	if (status)  goto TERMINATE;

	// Turn on traditional search for use with control callbacks 
	status = CPXsetintparam(env, CPX_PARAM_MIPSEARCH, CPX_MIPSEARCH_TRADITIONAL);
	if (status)  goto TERMINATE;

	// Let MIP callbacks work on the original model
	status = CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF);
	if (status)  goto TERMINATE;

	// Let MIP do not reduce LP
	status = CPXsetintparam(env, CPX_PARAM_REDUCE, CPX_PREREDUCE_NOPRIMALORDUAL);
	if (status)  goto TERMINATE;

	// Let MIP do not presolve MIP
	//status = CPXsetintparam (env, CPX_PARAM_PREIND, CPX_OFF);
	//if ( status )  goto TERMINATE;

	*Cuts = 0;
	cur_numcols = CPXgetnumcols(env, lp);
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

TERMINATE:

	return 0;
}

//-----------------------------------------------------------------------------
// Function implementing a custom cutting plane procedure called by callback
//-----------------------------------------------------------------------------
//
static int CPXPUBLIC mycutcallback(CPXCENVptr env,
	void       *cbdata,
	int        wherefrom,
	void       *cbhandle,
	int        *useraction_p)
{
	int status = 0;

	CUTINFOptr cutinfo = (CUTINFOptr)cbhandle;

	int      numcols = cutinfo->numcols;
	int      numtoadd = cutinfo->numtoadd;
	int      numcuts = cutinfo->num;
	double   *x = cutinfo->x;
	int      *beg = cutinfo->beg;
	int      *ind = cutinfo->ind;
	double   *val = cutinfo->val;
	double   *rhs = cutinfo->rhs;
	//double   cutvio;
	//int      i, j, k, cutnz;
	//long     i1, j1, k1;

	*useraction_p = CPX_CALLBACK_DEFAULT;

	status = CPXgetcallbacknodex(env, cbdata, wherefrom, x, 0, numcols - 1);
	if (status)
	{
		fprintf(stderr, "Failed to get node solution.\n");
		goto TERMINATE;
	}

	// Separate
	// ...

 //ADDCUT:
	// Use a cut violation tolerance of 0.01
	//if (cutvio<-0.001) 
	//{
	   // status = CPXcutcallbackadd(env,cbdata,wherefrom,k1,rhs[0],'G',ind,val,0);
	   // if (status)
	   // {
		  //  fprintf (stderr, "Failed to add cut.\n");
		  //  goto TERMINATE;
	   // }
	//}

	// Tell CPLEX that cuts have been created
	*useraction_p = CPX_CALLBACK_SET;

TERMINATE:

	return (status);

}

//-----------------------------------------------------------------------------
// Function implementing a custom cutting plane procedure called by callback
//-----------------------------------------------------------------------------
//
static int CPXPUBLIC myNewcutcallback(CPXCENVptr env,
	void       *cbdata,
	int        wherefrom,
	void       *cbhandle,
	int        *useraction_p)
{
	//long i;
	int status = 0;
	//int BolAdd;
	//char sense;
	//int purgeable;
	//int param;
	//int ivalue;

	CUTINFOptr cutinfo = (CUTINFOptr)cbhandle;

	*useraction_p = CPX_CALLBACK_DEFAULT;

	// Are we solving a subproblem?
	if (cutinfo->FlagLP)
		return (status);

	// Do we have reached the maximum number of cuts?  
	if (cutinfo->numtoadd > cutinfo->maxcuts)
		return (status);

	status = CPXgetcallbacknodex(env, cbdata, wherefrom, cutinfo->x, 0, cutinfo->numcols - 1);
	if (status)
	{
		fprintf(stderr, "Failed to get node solution.\n");
		goto TERMINATE;
	}

	// Separate
	// ...

//AddCut:

	// Use a cut violation tolerance of 0.01
	//if (BolAdd==1) 
	//{
	//	cutinfo->numtoadd++;
	//	cutinfo->Cuts[0]++;
	//	if ((cutinfo->numtoadd%1000)==0) printf("\n  User Cuts added: %d\n\n",cutinfo->numtoadd);

	//	//if ((cutinfo->numtoadd%5000)==0) 
	//	//	printf("\n  User Cuts added: %d\n\n",cutinfo->numtoadd);

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


