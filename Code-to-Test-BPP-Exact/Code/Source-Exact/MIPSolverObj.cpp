//-----------------------------------------------------------------------------
//  File: MIPSolverObj.cpp                                                       
//                                                                           
//  Project: MIP Solver Wrapper (reduced version)                  
//                                                                           
//  Author: M.A. Boschetti and S. Novellani                                                   
//                                                                           
//  Last update: 21.04.2025                                                  
//-----------------------------------------------------------------------------
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "MIPSolverObj.h"

//-----------------------------------------------------------------------------
//  Constructor                                    
//-----------------------------------------------------------------------------
//
MIPSolverObj::MIPSolverObj(void)
{
	Solver=1;

	minmax = 1;     // Default: minimize
	ncols = 0;      // Number of columns/variables
	nrows = 0;      // Number of rows/constraints
	nz = 0;         // Number of non-zero coefficients into the constraint matrix
	probtype = 0;   // Problem Type (-1: err; CPXPROB_LP=0; CPXPROB_MILP=1; etc.)
	lpsolver = 0;   // LP Solver (0: automatic; 1: Primal; 2: Dual; 3: Net; 4: Barrier)
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

	rmatbeg = NULL;  // Pointers to the beginning of each rows of the constraint matrix
	rmatcnt = NULL;  // Cardinality of each rows of the constraint matrix
	rmatind = NULL;  // Column index associated to each coefficient of the constraint matrix
	rmatval = NULL;  // Value of each coefficient of the constraint matrix

	nsos = 0;        // Number of SOS (Special Ordered Set) 
	nsosnz = 0;      // Total number of membres in all SOS added
	typesos = NULL;  // Array containing SOS type (CPX_TYPE_SOS1='1', CPX_TYPE_SOS2='2')
	sosbeg = NULL;   // Pointers to the beginning of each SOS
	sosind = NULL;   // Index associated to the SOS member
	soswt = NULL;    // Weight associated to the SOS member

	xt = NULL;       // Variable Values

}

//-----------------------------------------------------------------------------
//  Destructor                                    
//-----------------------------------------------------------------------------
//
MIPSolverObj::~MIPSolverObj(void)
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
	if (indices != NULL) free(indices);
	if (priority != NULL) free(priority);
	if (direction != NULL) free(direction);
	if (xctype != NULL) free(xctype);

	if (rmatbeg != NULL) free(rmatbeg);
	if (rmatcnt != NULL) free(rmatcnt);
	if (rmatind != NULL) free(rmatind);
	if (rmatval != NULL) free(rmatval);

	if (typesos != NULL) free(typesos);
	if (sosbeg != NULL) free(sosbeg);
	if (sosind != NULL) free(sosind);
	if (soswt != NULL) free(soswt);

	if (xt != NULL) free(xt);
}

//-----------------------------------------------------------------------------
//  Open: open an environment and open a new problem 
//-----------------------------------------------------------------------------
//
int MIPSolverObj::Open(void)
{
	int status;

	if (Solver == 0)
		status = MIPCplex.Open();
	else if (Solver == 1)
		status = MIPGurobi.Open();
	//else if (Solver==2)
	//	status = MIPCoin.Open();
	else
		status = 9999;

	return status;
}

//-----------------------------------------------------------------------------
//  Open: open an environment and open a new problem 
//-----------------------------------------------------------------------------
//
int MIPSolverObj::Open(long s)
{
	int status;

	if ((s == 0) || (s == 1) || (s == 2))
	{
		Solver = s;
		status = Open();
	}
	else
		status = 8888;

	return status;
}

//-----------------------------------------------------------------------------
//  Close: close environment and problem 
//-----------------------------------------------------------------------------
//
int MIPSolverObj::Close(void)
{
	int status;

	if (Solver == 0)
		status = MIPCplex.Close();
	else if (Solver == 1)
		status = MIPGurobi.Close();
	//else if (Solver==2)
	//	status = MIPCoin.Open();
	else
		status = 9999;

	return status;
}

//-----------------------------------------------------------------------------
//  Close: close environment and problem 
//-----------------------------------------------------------------------------
//
int MIPSolverObj::Close(long s)
{
	int status;

	if ((s == 0) || (s == 1) || (s == 2))
	{
		Solver = s;
		status = Close();
	}
	else
		status = 8888;

	return status;
}

//-----------------------------------------------------------------------------
//  Reset: close the existing problem and create a new one 
//-----------------------------------------------------------------------------
//
int MIPSolverObj::Reset(void)
{
	int status;

	if (Solver == 0)
		status = MIPCplex.Reset();
	else if (Solver == 1)
		status = MIPGurobi.Reset();
	//else if (Solver==2)
	//	status = MIPCoin.Reset();
	else
		status = 9999;

	return status;
}

//-----------------------------------------------------------------------------
//  Reset: close the existing problem and create a new one 
//-----------------------------------------------------------------------------
//
int MIPSolverObj::Reset(long s)
{
	int status;

	if ((s == 0) || (s == 1) || (s == 2))
	{
		Solver = s;
		status = Reset();
	}
	else
		status = 8888;

	return status;
}

//-----------------------------------------------------------------------------
//  MallocCols: Allocate the "column" data structure
//-----------------------------------------------------------------------------
//
int MIPSolverObj::MallocCols(void)
{
	int status = 0;

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

	return status;
}

//-----------------------------------------------------------------------------
//  MallocRows: Allocate the "row" data structure
//-----------------------------------------------------------------------------
//
int MIPSolverObj::MallocRows(void)
{
	int status = 0;

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

	return status;
}

//-----------------------------------------------------------------------------
//  Malloc_Valid_Inequalities: Allocate the "row" data structure to add 
//  valid inequalities
//-----------------------------------------------------------------------------
//
int MIPSolverObj::Malloc_Valid_Inequalities(long newrows)
{
	int status = 0;
	long maxnz;

	maxnz = newrows * 1000;

	rrhs = (double*)malloc((newrows + 4) * sizeof(double));
	rsense = (char*)malloc((newrows + 4) * sizeof(char));
	rmatbeg = (int*)malloc((newrows + 4) * sizeof(int));
	rmatcnt = (int*)malloc((newrows + 4) * sizeof(int));
	rmatind = (int*)malloc(maxnz * sizeof(int));
	rmatval = (double*)malloc(maxnz * sizeof(double));

	return status;
}

//-----------------------------------------------------------------------------
//  Free_Valid_Inequalities: Deallocate the "row" data structure to add 
//  valid inequalities
//-----------------------------------------------------------------------------
//
int MIPSolverObj::Free_Valid_Inequalities(void)
{
	int status = 0;

	if (rrhs != NULL)
	{
		free(rrhs);
		rrhs = NULL;
	}

	if (rsense != NULL)
	{
		free(rsense);
		rsense = NULL;
	}

	if (rmatbeg != NULL)
	{
		free(rmatbeg);
		rmatbeg = NULL;
	}

	if (rmatcnt != NULL)
	{
		free(rmatcnt);
		rmatcnt = NULL;
	}

	if (rmatind != NULL)
	{
		free(rmatind);
		rmatind = NULL;
	}

	if (rmatval != NULL)
	{
		free(rmatval);
		rmatval = NULL;
	}

	return status;
}

//-----------------------------------------------------------------------------
//  FreeCols: Deallocate the "column" data structure
//-----------------------------------------------------------------------------
//
int MIPSolverObj::FreeCols(void)
{
	int status = 0;

	if (obj != NULL)
	{
		free(obj);
		obj = NULL;
	}
	if (rhs != NULL)
	{
		free(rhs);
		rhs = NULL;
	}
	if (sense != NULL)
	{
		free(sense);
		sense = NULL;
	}
	if (matbeg != NULL)
	{
		free(matbeg);
		matbeg = NULL;
	}
	if (matcnt != NULL)
	{
		free(matcnt);
		matcnt = NULL;
	}
	if (matind != NULL)
	{
		free(matind);
		matind = NULL;
	}
	if (matval != NULL)
	{
		free(matval);
		matval = NULL;
	}
	if (lb != NULL)
	{
		free(lb);
		lb = NULL;
	}
	if (ub != NULL)
	{
		free(ub);
		ub = NULL;
	}
	if (xt != NULL)
	{
		free(xt);
		xt = NULL;
	}
	if (indices != NULL)
	{
		free(indices);
		indices = NULL;
	}
	if (priority != NULL)
	{
		free(priority);
		priority = NULL;
	}
	if (direction != NULL)
	{
		free(direction);
		direction = NULL;
	}
	if (xctype != NULL)
	{
		free(xctype);
		xctype = NULL;
	}

	if (typesos != NULL)
	{
		free(typesos);
		typesos = NULL;
	}
	if (sosbeg != NULL)
	{
		free(sosbeg);
		sosbeg = NULL;
	}
	if (sosind != NULL)
	{
		free(sosind);
		sosind = NULL;
	}
	if (soswt != NULL)
	{
		free(soswt);
		soswt = NULL;
	}

	return status;
}

//-----------------------------------------------------------------------------
//  FreeRows: Deallocate the "row" data structure
//-----------------------------------------------------------------------------
//
int MIPSolverObj::FreeRows(void)
{
	int status = 0;

	free(obj);
	free(rhs);
	free(sense);
	free(rmatbeg);
	free(rmatcnt);
	free(rmatind);
	free(rmatval);
	free(lb);
	free(ub);
	free(xt);
	free(indices);
	free(priority);
	free(direction);
	free(xctype);

	if (typesos != NULL) free(typesos);
	if (sosbeg != NULL) free(sosbeg);
	if (sosind != NULL) free(sosind);
	if (soswt != NULL) free(soswt);

	return status;
}

//-----------------------------------------------------------------------------
//  CopyLPbyCols
//-----------------------------------------------------------------------------
//
int MIPSolverObj::CopyLPbyCols(void)
{
	int status = 0;

	if (Solver == 0)
	{
		MIPCplex.ncols = ncols;
		MIPCplex.nrows = nrows;
		MIPCplex.nz = nz;
		MIPCplex.objsen = objsen;
		MIPCplex.obj = obj;
		MIPCplex.rhs = rhs;
		MIPCplex.sense = sense;
		MIPCplex.matbeg = matbeg;
		MIPCplex.matcnt = matcnt;
		MIPCplex.matind = matind;
		MIPCplex.matval = matval;
		MIPCplex.lb = lb;
		MIPCplex.ub = ub;
		MIPCplex.xctype = xctype;
		MIPCplex.indices = indices;
		MIPCplex.xt = xt;
		MIPCplex.rngval = rngval;
		MIPCplex.minmax = minmax;

		status = MIPCplex.CopyLPbyCols();
	}
	else if (Solver == 1)
	{
		MIPGurobi.ncols = ncols;
		MIPGurobi.nrows = nrows;
		MIPGurobi.nz = nz;
		MIPGurobi.objsen = objsen;
		MIPGurobi.obj = obj;
		MIPGurobi.rhs = rhs;
		MIPGurobi.sense = sense;
		MIPGurobi.matbeg = matbeg;
		MIPGurobi.matcnt = matcnt;
		MIPGurobi.matind = matind;
		MIPGurobi.matval = matval;
		MIPGurobi.lb = lb;
		MIPGurobi.ub = ub;
		MIPGurobi.xctype = xctype;
		MIPGurobi.indices = indices;
		MIPGurobi.xt = xt;
		MIPGurobi.rngval = rngval;
		MIPGurobi.minmax = minmax;

		status = MIPGurobi.CopyLPbyCols();
	}
	else if (Solver == 2)
	{
		//MIPCoin.ncols = ncols;
		//MIPCoin.nrows = nrows;
		//MIPCoin.nz = nz;
		//MIPCoin.objsen = objsen;
		//MIPCoin.obj = obj;
		//MIPCoin.rhs = rhs;
		//MIPCoin.sense = sense;
		//MIPCoin.matbeg = matbeg;
		//MIPCoin.matcnt = matcnt;
		//MIPCoin.matind = matind;
		//MIPCoin.matval = matval;
		//MIPCoin.lb = lb;
		//MIPCoin.ub = ub;
		//MIPCoin.xctype = xctype;
		//MIPCoin.indices = indices;
		//MIPCoin.xt = xt;
		//MIPCoin.rngval = rngval;
		//MIPCoin.minmax = minmax;

		//status = MIPCoin.CopyLPbyCols();
	}
	else
		status = 9999;

	return status;
}

//-----------------------------------------------------------------------------
//  CopyLPbyRows
//-----------------------------------------------------------------------------
//
int MIPSolverObj::CopyLPbyRows(void)
{
	int status = 0;

	if (Solver == 0)
	{
		MIPCplex.ncols = ncols;
		MIPCplex.nrows = nrows;
		MIPCplex.nz = nz;
		MIPCplex.objsen = objsen;
		MIPCplex.obj = obj;
		MIPCplex.rhs = rhs;
		MIPCplex.sense = sense;
		MIPCplex.rmatbeg = rmatbeg;
		MIPCplex.rmatcnt = rmatcnt;
		MIPCplex.rmatind = rmatind;
		MIPCplex.rmatval = rmatval;
		MIPCplex.lb = lb;
		MIPCplex.ub = ub;
		MIPCplex.xctype = xctype;
		MIPCplex.indices = indices;
		MIPCplex.xt = xt;
		MIPCplex.rngval = rngval;
		MIPCplex.minmax = minmax;

		status = MIPCplex.CopyLPbyRows();
	}
	else if (Solver == 1)
	{
		MIPGurobi.ncols = ncols;
		MIPGurobi.nrows = nrows;
		MIPGurobi.nz = nz;
		MIPGurobi.objsen = objsen;
		MIPGurobi.obj = obj;
		MIPGurobi.rhs = rhs;
		MIPGurobi.sense = sense;
		MIPGurobi.rmatbeg = rmatbeg;
		MIPGurobi.rmatcnt = rmatcnt;
		MIPGurobi.rmatind = rmatind;
		MIPGurobi.rmatval = rmatval;
		MIPGurobi.lb = lb;
		MIPGurobi.ub = ub;
		MIPGurobi.xctype = xctype;
		MIPGurobi.indices = indices;
		MIPGurobi.xt = xt;
		MIPGurobi.rngval = rngval;
		MIPGurobi.minmax = minmax;

		status = MIPGurobi.CopyLPbyRows();
	}
	else if (Solver == 2)
	{
		//MIPCoin.ncols = ncols;
		//MIPCoin.nrows = nrows;
		//MIPCoin.nz = nz;
		//MIPCoin.objsen = objsen;
		//MIPCoin.obj = obj;
		//MIPCoin.rhs = rhs;
		//MIPCoin.sense = sense;
		//MIPCoin.rmatbeg = rmatbeg;
		//MIPCoin.rmatcnt = rmatcnt;
		//MIPCoin.rmatind = rmatind;
		//MIPCoin.rmatval = rmatval;
		//MIPCoin.lb = lb;
		//MIPCoin.ub = ub;
		//MIPCoin.xctype = xctype;
		//MIPCoin.indices = indices;
		//MIPCoin.xt = xt;
		//MIPCoin.rngval = rngval;
		//MIPCoin.minmax = minmax;

		//status = MIPCoin.CopyLPbyRows();
	}
	else
		status = 9999;

	return status;
}

//-----------------------------------------------------------------------------
//  ReadMPS
//-----------------------------------------------------------------------------
//
int MIPSolverObj::ReadMPS(const char * fname)
{
	int status = 0;

	if (Solver == 0)
		status = MIPCplex.ReadMPS(fname);
	else if (Solver == 1)
		status = MIPGurobi.ReadMPS(fname);
	//else if (Solver==2)
	//	status = MIPCoin.ReadMPS(fname);
	else
		status = 9999;

	return status;
}

//-----------------------------------------------------------------------------
//  Set a Lower Bound for the MIP
//-----------------------------------------------------------------------------
//
int MIPSolverObj::SetLB(double LB)
{
	int status;

	if (Solver == 0)
		status = MIPCplex.SetLB(LB);
	else if (Solver == 1)
		status = MIPGurobi.SetLB(LB);
	//else if (Solver==2)
	//	status = MIPCoin.SetLB(LB);
	else
		status = 9999;

	return status;
}

//-----------------------------------------------------------------------------
//  SetMIP
//-----------------------------------------------------------------------------
//
int MIPSolverObj::SetMIP(double TLim)
{
	int status;

	if (Solver == 0)
		status = MIPCplex.SetMIP(TLim);
	else if (Solver == 1)
		status = MIPGurobi.SetMIP(TLim);
	//else if (Solver==2)
	//	status = MIPCoin.SetMIP(TLim);
	else
		status = 9999;

	return status;
}

//-----------------------------------------------------------------------------
//  Add New Rows
//-----------------------------------------------------------------------------
//
int MIPSolverObj::AddNewRows(void)
{
	int status;

	if (Solver == 0)
	{
		MIPCplex.rncols = rncols;
		MIPCplex.rnrows = rnrows;
		MIPCplex.rnz = rnz;
		MIPCplex.rrhs = rrhs;
		MIPCplex.rsense = rsense;
		MIPCplex.rmatbeg = rmatbeg;
		MIPCplex.rmatcnt = rmatcnt;
		MIPCplex.rmatind = rmatind;
		MIPCplex.rmatval = rmatval;
		status = MIPCplex.AddNewRows();
	}
	else if (Solver == 1)
	{
		MIPGurobi.rncols = rncols;
		MIPGurobi.rnrows = rnrows;
		MIPGurobi.rnz = rnz;
		MIPGurobi.rrhs = rrhs;
		MIPGurobi.rsense = rsense;
		MIPGurobi.rmatbeg = rmatbeg;
		MIPGurobi.rmatcnt = rmatcnt;
		MIPGurobi.rmatind = rmatind;
		MIPGurobi.rmatval = rmatval;
		status = MIPGurobi.AddNewRows();
	}
	//else if (Solver==1)
	//	status = MIPCoin.AddNewRows();
	else
		status = 9999;

	return status;
}

//-----------------------------------------------------------------------------
//  SolveMIP
//-----------------------------------------------------------------------------
//
int MIPSolverObj::SolveMIP(double *Zopt, double *Zlb, double *Gap, long *Nodes, long *Cuts)
{
	int status;

	if (Solver == 0)
		status = MIPCplex.SolveMIP(Zopt, Zlb, Gap, Nodes, Cuts);
	else if (Solver == 1)
		status = MIPGurobi.SolveMIP(Zopt, Zlb, Gap, Nodes, Cuts);
	//else if (Solver==2)
	//	status = MIPCoin.SolveMIP(Zopt,Zlb,Gap,Nodes,Cuts);
	else
		status = 9999;

	return status;
}

//-----------------------------------------------------------------------------
//  SolveLP
//-----------------------------------------------------------------------------
//
int MIPSolverObj::SolveLP(double *Zopt, double *dual, long *niter)
{
	int status;

	*niter = 0;

	if (Solver == 0)
	{
		MIPCplex.lpsolver = lpsolver;
		status = MIPCplex.SolveLP(Zopt, dual, niter);
	}
	else if (Solver == 1)
		status = MIPGurobi.SolveLP(Zopt);
	//else if (Solver==2)
	//	status = MIPCoin.SolveMIP(Zopt,Zlb,Gap,Nodes,Cuts);
	else
		status = 9999;

	return status;
}

//-----------------------------------------------------------------------------
//  Setup Custom Cutting Plane
//-----------------------------------------------------------------------------
//
int MIPSolverObj::CuttingPlane(long n, long m, double* u, double* p, double* pmax,
	long *nY, long **Y, long *nUOR, long **UOR, long *nUAND, long **UAND, long *FDepA, long *FDepP,
	int Pred, int Cover, int Lifting, int KnapSol, int MaxCuts, long *Cuts)
{
	int status;

	if (Solver == 0)
		status = MIPCplex.CuttingPlane(n, m, u, p, pmax, nY, Y, nUOR, UOR, nUAND, UAND, FDepA, FDepP, Pred, Cover, Lifting, KnapSol, MaxCuts, Cuts);
	else if (Solver == 1)
		status = MIPGurobi.CuttingPlane(n, m, u, p, pmax, nY, Y, nUOR, UOR, nUAND, UAND, FDepA, FDepP, Pred, Cover, Lifting, KnapSol, MaxCuts, Cuts);
	//else if (Solver==2)
	//	status = MIPCoin.CuttingPlane(n,m,u,p,pmax,nY,Y,nUOR,UOR,nUAND,UAND,FDepA,FDepP,Pred,Cover,Lifting,KnapSol,MaxCuts,Cuts);
	else
		status = 9999;

	return status;
}

