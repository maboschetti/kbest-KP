//-----------------------------------------------------------------------------
//  File: CoinObj.cpp                                                       
//                                                                           
//  Project: Coin-OR (CBC) solver interface                   
//                                                                           
//  Author: M.A. Boschetti                                                   
//                                                                           
//  Last update: 06.11.2012                                                  
//-----------------------------------------------------------------------------
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "CoinObj.h"
//#include "Knapsack.h"

//-----------------------------------------------------------------------------
//  Constructor                                    
//-----------------------------------------------------------------------------
//
CoinObj::CoinObj(void)
{
	int param;
	int ivalue;
	double rvalue;

	// Open Coin-OR Environment (Solver Interface)
	//env = NULL;
	//status = GRBloadenv(&env, NULL);
	//if (status || env == NULL) 
	//{
	//	printf("\nError... GRBloadenv...\n");
	//	env = NULL;
	//	return;
	//}

	// Coin-OR creates the model when an instance is loaded
	//lp = NULL;

	// Set on/off the output on the screen
	//status = GRBsetintparam(env, GRB_INT_PAR_OUTPUTFLAG, 1);
	//if (status>0) return;

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
CoinObj::~CoinObj(void)
{
	long i,j;

	// Close the Coin-OR Enviroment (Solver Interface)... Not required by Coin-OR!
	//if (lp != NULL)
	//{
	//	GRBfreemodel(lp);
	//	lp = NULL;
	//}
	//if (env != NULL)
	//{
	//	GRBfreeenv(env);
	//	env = NULL;
	//}

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
//  Open: open an environment (Solver Interface) and open a new problem (Model)
//-----------------------------------------------------------------------------
//
int CoinObj::Open(void)
{
	int status;

	// Not required by Coin-OR!
	//status = GRBloadenv(&env, NULL);
	//if (status || env == NULL) 
	//{
	//	printf("\nError... GRBloadenv...\n");
	//	env = NULL;
	//	return status;
	//}

	//if (lp != NULL)
	//{
	//	GRBfreemodel(lp);
	//	lp = NULL;
	//}

	//status = GRBnewmodel(env, &lp, "lp", 0,	NULL, NULL, NULL, NULL, NULL);
	//if (status)
	//{
	//	printf("\nError... GRBnewmodel...\n");
	//	lp = NULL;
	//	return status;
	//}

	//FreeCols();

//00	status = GRBsetintparam(env, GRB_INT_PAR_OUTPUTFLAG, 1);
//00	if (status>0) return status;

	return 0;
}


//-----------------------------------------------------------------------------
//  Reset: close the existing problem and create a new one 
//-----------------------------------------------------------------------------
//
int CoinObj::Reset(void)
{
	// Not required by Coin-OR!
	//if (lp != NULL)
	//{
	//	GRBfreemodel(lp);
	//	lp = NULL;
	//}

	//status = GRBnewmodel(env, &lp, "lp", 0,	NULL, NULL, NULL, NULL, NULL);
	//if (status)
	//{
	//	printf("\nError... GRBnewmodel...\n");
	//	lp = NULL;
	//	return status;
	//}

	FreeCols();

	return 0;
}


//-----------------------------------------------------------------------------
//  MallocCols: Allocate the "column" data structure
//-----------------------------------------------------------------------------
//
int CoinObj::MallocCols(void)
{
	bolfree = 1;
	obj = (double*)malloc(ncols*sizeof(double));
	rhs = (double*)malloc(nrows*sizeof(double));
	sense = (char*)malloc(nrows*sizeof(char));
	matbeg = (int*)malloc((ncols+4)*sizeof(int));
	matcnt = (int*)malloc((ncols+4)*sizeof(int));
	matind = (int*)malloc(nz*sizeof(int));          
	matval = (double*)malloc(nz*sizeof(double));       
	lb = (double*)malloc(ncols*sizeof(double));
	ub = (double*)malloc(ncols*sizeof(double));
	xt = (double*)malloc(ncols*sizeof(double));
	indices = (int*)malloc((ncols+4)*sizeof(int));
	priority = (int*)malloc((ncols+4)*sizeof(int));
	direction = (int*)malloc((ncols+4)*sizeof(int));
	xctype = (char*)malloc((ncols+4)*sizeof(char));

	if (nsos>0)
	{
		typesos = (char*)malloc(nsos*sizeof(char));
		sosbeg = (int*)malloc((nsos+4)*sizeof(int));
		sosind = (int*)malloc(ncols*sizeof(int));
		soswt = (double*)malloc(ncols*sizeof(double));
	}

	//rmatbeg = (int*)malloc((nrows+4)*sizeof(int));
	//rmatcnt = (int*)malloc((nrows+4)*sizeof(int));
	//rmatind = (int*)malloc(nz*sizeof(int));          
	//rmatval = (double*)malloc(nz*sizeof(double));       

	return 0;
}


//-----------------------------------------------------------------------------
//  MallocRows: Allocate the "row" data structure
//-----------------------------------------------------------------------------
//
int CoinObj::MallocRows(void)
{
	bolfree = 1;
	obj = (double*)malloc(ncols*sizeof(double));
	rhs = (double*)malloc(nrows*sizeof(double));
	sense = (char*)malloc(nrows*sizeof(char));
	rmatbeg = (int*)malloc((nrows+4)*sizeof(int));
	rmatcnt = (int*)malloc((nrows+4)*sizeof(int));
	rmatind = (int*)malloc(nz*sizeof(int));          
	rmatval = (double*)malloc(nz*sizeof(double));       
	lb = (double*)malloc(ncols*sizeof(double));
	ub = (double*)malloc(ncols*sizeof(double));
	xt = (double*)malloc(ncols*sizeof(double));
	indices = (int*)malloc((ncols+4)*sizeof(int));
	priority = (int*)malloc((ncols+4)*sizeof(int));
	direction = (int*)malloc((ncols+4)*sizeof(int));
	xctype = (char*)malloc((ncols+4)*sizeof(char));

	if (nsos>0)
	{
		typesos = (char*)malloc(nsos*sizeof(char));
		sosbeg = (int*)malloc((nsos+4)*sizeof(int));
		sosind = (int*)malloc(ncols*sizeof(int));
		soswt = (double*)malloc(ncols*sizeof(double));
	}

	return 0;
}


//-----------------------------------------------------------------------------
//  FreeCols: Deallocate the "column" data structure
//-----------------------------------------------------------------------------
//
int CoinObj::FreeCols(void)
{
	long i;

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

	return 0;
}


//-----------------------------------------------------------------------------
//  FreeRows: Deallocate the "row" data structure
//-----------------------------------------------------------------------------
//
int CoinObj::FreeRows(void)
{
	long i;

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

	free(typesos);
	free(sosbeg);
	free(sosind);
	free(soswt);

	free(cutinfo.x);
	free(cutinfo.beg);
	free(cutinfo.ind); 
	free(cutinfo.val);
	free(cutinfo.rhs);

	return 0;
}


//-----------------------------------------------------------------------------
//  CopyLPbyCols
//-----------------------------------------------------------------------------
//
int CoinObj::CopyLPbyCols(void)
{
	double dvalue;

	status = 0;

	//status = GRBaddvars(lp,ncols,nz,matbeg,matind,matval,obj,lb,ub,NULL,NULL);
	//if (status>0)
	//	return status;

	//status = GRBaddconstrs(lp,nrows,0,NULL,NULL,NULL,sense,rhs,NULL);
	//if (status>0)
	//	return status;

	//// Set problem: minimization or maximization
	//status = GRBsetintattr(lp,"ModelSense",minmax);
	//if (status>0)
	//	return status;

	//status = GRBupdatemodel(lp);

	//status = GRBloadmodel(env,&lp,"lp",minmax,0.0,obj,sense,rhs,matbeg,matcnt,matind,matval,lb,ub,xctype,NULL,NULL);

	lp.loadBlock(ncols,nrows,matbeg,matind,matval,lb,ub,obj,sense,rhs,NULL);

	dvalue = (double)minmax;
	lp.setOptimizationDirection(dvalue);

	return status;
}


//-----------------------------------------------------------------------------
//  CopyLPbyRows
//-----------------------------------------------------------------------------
//
int CoinObj::CopyLPbyRows(void)
{
	int i,j,status;
	int start,end;
	double dvalue;

	status = 0;

	//status = GRBnewmodel(env,&lp,"lp",ncols,obj,lb,ub,xctype,NULL);
	//if (status>0) return status;

	// Add columns
	for (j=0; j<ncols; j++)
	{
		lp.setColumnBounds(j,lb[j],ub[j]);
		lp.setObjective(j,obj[j]);
		if (xctype[j]!='C')
			lp.setInteger(j);
	}

	//status = GRBaddconstrs(lp,nrows,nz,rmatbeg,rmatind,rmatval,sense,rhs,NULL);
	//if (status>0) return status;

	// Add rows
	for (i=0; i<nrows; i++)
	{
		start = rmatbeg[i];
		end = rmatbeg[i+1];
		if (sense[i]=='E')
			lp.addRow(end-start,rmatind+start,rmatval+start,rhs[i],rhs[i]);
		else if (sense[i]=='L')
			lp.addRow(end-start,rmatind+start,rmatval+start,-COIN_DBL_MAX,rhs[i]);
		else if (sense[i]=='G')
			lp.addRow(end-start,rmatind+start,rmatval+start,rhs[i],+COIN_DBL_MAX);
		else
			return 1;
	}

	//status = GRBsetintattr(lp,"ModelSense",minmax);
	//if (status>0) return status;

	dvalue = (double)minmax;
	lp.setOptimizationDirection(dvalue);

	//dvalue = lp.optimizationDirection();

	//status = GRBupdatemodel(lp);
	//if (status>0) return status;

	return status;
}


//-----------------------------------------------------------------------------
//  ReadMPS
//-----------------------------------------------------------------------------
//
int CoinObj::ReadMPS(const char * fname)
{
	int status;
	int nz1,surplus;
	//double zopt,zlb,gap;
	//long node,cuts;
	long i,j,k;

	// Load the LP instance
	//status = GRBreadmodel(env,fname,&lp);
	//if (status) return status;
	status = env.readMps(fname,"");
	if (status) return status;

	// Load in memory the instance

	// Sizes, Types, etc.
	//status = GRBgetintattr(lp,"NumVars",&ncols);
	//if (status) return status;
	ncols = env.getNumCols();

	//status = GRBgetintattr(lp,"NumConstrs",&nrows);
	//if (status) return status;
	nrows = env.getNumRows();

	//status = GRBgetintattr(lp,"NumBinVars",&nbin);
	//if (status) return status;
	nbin = -1;

	//status = GRBgetintattr(lp,"NumIntVars",&nint);
	//if (status) return status;
	nint = -1;

	//status = GRBgetintattr(lp,"NumNZs",&nz);
	//if (status) return status;
	nz = env.getNumElements();

	//status = GRBgetintattr(lp,"ModelSense",&minmax);
	//if (status) return status;
	minmax = (int)env.getObjSense();

	//status = GRBgetintattr(lp,"IsMIP",&probtype);
	//if (status) return status;
	probtype = -1;

	// Allocate memory
	MallocCols();

	// Auxiliary
	const CoinPackedMatrix *matrix = env.getMatrixByCol();
	const double *element = matrix->getElements();
	const int *row = matrix->getIndices();
	const CoinBigIndex *colstart = matrix->getVectorStarts();
	const int *colLength = matrix->getVectorLengths();
	const double *rowLower = env.getRowLower();
	const double *rowUpper = env.getRowUpper();
	const double *colLower = env.getColLower();
	const double *colUpper = env.getColUpper();
	const double *objcoeff = env.getObjCoefficients();

	// Obj, Cols, etc.
	for (j=0; j<ncols; j++)
	{
		obj[j] = objcoeff[j];

		matbeg[j] = colstart[j];
		for (k=colstart[j]; k<colstart[j+1]; k++)
		{
			matind[k] = row[k];
			matval[k] = element[k];
		}

		lb[j] = colLower[j];
		ub[j] = colUpper[j];
		if (env.isContinuous(j))
			xctype[j] = 'C';
		else
			xctype[j] = 'I';
	}
	matbeg[ncols] = nz;

	//status = GRBgetdblattrarray(lp,"Obj",0,ncols,obj);
	//if (status) return status;

	//status = GRBgetvars(lp,&nz1,matbeg,matind,matval,0,ncols);
	//if (status) return status;
	//if (nz1!=nz) 
	//{
	//	status=999;
	//	return status;
	//}
	//matbeg[ncols] = nz;
	for (i=0; i<ncols; i++)
		matcnt[i] = matbeg[i+1]-matbeg[i];

	//status = GRBgetconstrs(lp,&nz1,rmatbeg,rmatind,rmatval,0,nrows);
	//if (status) return status;
	//if (nz1!=nz) 
	//{
	//	status=999;
	//	return status;
	//}
	//rmatbeg[nrows] = nz;
	//for (i=0; i<nrows; i++)
	//	rmatcnt[i] = rmatbeg[i+1]-rmatbeg[i];

	//status = GRBgetdblattrarray(lp,"LB",0,ncols,lb);
	//if (status) return status;

	//status = GRBgetdblattrarray(lp,"UB",0,ncols,ub);
	//if (status) return status;

	//status = GRBgetcharattrarray(lp,"VType",0,ncols,xctype);
	//if (status) return status;

	for (i=0; i<nrows; i++)
	{
		if (rowLower[i]-0.001<-COIN_DBL_MAX)
		{
			rhs[i] = rowUpper[i];
			sense[i] = 'L';
		}
		else if (rowUpper[i]+0.001>+COIN_DBL_MAX)
		{
			rhs[i] = rowLower[i];
			sense[i] = 'G';
		}
		else  // We do not consider rngval constraints!
		{
			rhs[i] = rowUpper[i];
			sense[i] = 'E';
		}
	}

	//status = GRBgetdblattrarray(lp,"RHS",0,nrows,rhs);
	//if (status) return status;

	//status = GRBgetcharattrarray(lp,"Sense",0,nrows,sense);
	//if (status) return status;

	//Reset();

	//for (long i=0; i<nrows; i++)
	//{
	//	if (rmatbeg[i]>nz)
	//		nz=nz;
	//	for (long k=rmatbeg[i]; k<rmatbeg[i+1]; k++)
	//		if (rmatind[k]>=ncols)
	//			nz=nz;
	//}

	//CopyLPbyCols();
	////CopyLPbyRows();

	//SolveMIP(&zopt,&zlb,&gap,&node,&cuts);

	return status;
}


//-----------------------------------------------------------------------------
//  Set a Lower Bound for the MIP
//-----------------------------------------------------------------------------
//
int CoinObj::SetLB(double LB)
{
	// Set Parameter
	//status = GRBsetdblparam(env, GRB_DBL_PAR_CUTOFF,LB);

	return status;
}


//-----------------------------------------------------------------------------
//  SetMIP
//-----------------------------------------------------------------------------
//
int CoinObj::SetMIP(double TLim)
{
	int Param;
	int Value;
	double DValue;

	TimeLimit = TLim;

	// Change the problem from LP to MIP
	// status = CPXchgprobtype(env,lp,CPXPROB_MILP);

	// Charge the variable type (the default is "continuous")
	// status = CPXchgctype(env,lp,ncols,indices,xctype);

	//// Set Parameter
	//status = GRBsetdblparam(env,GRB_DBL_PAR_TIMELIMIT,TLim);
	//if (status) return status;
	//
	//DValue=10e-9;
	//status = GRBsetdblparam(env,GRB_DBL_PAR_MIPGAP,DValue);
	//if (status) return status;
	//
	//DValue=10e-10;
	//status = GRBsetdblparam(env,GRB_DBL_PAR_MIPGAPABS,DValue);
	//if (status) return status;
	//
	//// Load SOS 
	////status = CPXaddsos(env,lp,nsos,nsosnz,typesos,sosbeg,sosind,soswt,NULL);

	//// Set Parameter
	//Value=0;
	//status = GRBsetintparam(env, GRB_INT_PAR_COVERCUTS,Value);
	//if (status) return status;

	return status;
}


//-----------------------------------------------------------------------------
//  Sentinel Constraint
//-----------------------------------------------------------------------------
//
int CoinObj::Sentinel(long n, long m)
{
//	int i,j,k,k1;
//	int *beg=NULL;
//	int *ind=NULL;
//	char *sense=NULL;
//	double *rhs=NULL;
//	double *val=NULL;
//	char *type;
//	double rr;
//
//	// Build the Integer Sentinel Costraint
//	ind = (int *) malloc (ncols * sizeof (int));		if (ind==NULL) goto ESCI; 
//	val = (double *) malloc (ncols * sizeof (double));	if (val==NULL) goto ESCI;
//	beg = (int *) malloc ((m+4) * sizeof (int));		if (beg==NULL) goto ESCI;
//	sense = (char *) malloc (m * sizeof (char));		if (sense==NULL) goto ESCI;
//	rhs = (double *) malloc (m * sizeof (double));		if (rhs==NULL) goto ESCI;
//
//	//srand(0);
//
//	k=0;
//	for (i=0; i<m; i++)
//	{
//		rhs[i]=0.0;
//		sense[i]='E';
//		beg[i]=k;
//		for (j=0; j<n; j++)
//		{
//			k1=i*n+j;
//			//rr = (double)rand()/RAND_MAX;
//			//rr = 1;
//			//if ((xctype[k1]=='B')&&(rr<0.5))
//			if (xctype[k1]=='B')
//			{
//				ind[k]=indices[k1];
//				val[k]=1.0;
//				k++;
//			}
//		}
//		ind[k]=ncols;
//		val[k]=-1.0;
//		k++;
//		indices[ncols]=ncols;
//		priority[ncols]=1000;
////		direction[ncols] = CPX_BRANCH_DOWN; 
//		direction[ncols] = 0; 
//		xctype[ncols]='I';
//		ncols++;
//	}
//	beg[m]=k;
//
//	// Load the Integer Sentinel Constraints
//	type = xctype-m;
//	status = GRBaddvars(lp,m,0,NULL,NULL,NULL,NULL,NULL,NULL,type,NULL);
//	if (status) return status;
//
//	status = GRBaddconstrs(lp,m,k,beg,ind,val,sense,rhs,NULL);
//	if (status) return status;
//
//	status = GRBupdatemodel(lp);
//	if (status) return status;
//
//	goto TERMINATE;
//
//ESCI:
//	fprintf (stderr, "No memory for solution values.\n");
//
//TERMINATE:
//	free(beg);
//	free(ind);
//	free(val);
//	free(sense);
//	free(rhs);
//
	return status;
}


//-----------------------------------------------------------------------------
//  Add Fixed Cost Miner Constraints
//-----------------------------------------------------------------------------
//
int CoinObj::AddFCM(int nrange)
{
//	int i,j,k,k1;
//	int cnt;
//	double mincost;
//	double maxcost;
//	double delta;
//	double level,level0;
//	int *flag;
//	int *beg=NULL;
//	int *ind=NULL;
//	int *indt=NULL;
//	char *sense=NULL;
//	char *ctype=NULL;
//	char *type;
//	double *rhs=NULL;
//	double *val=NULL;
//	double *objt=NULL;
//	double rr;
//	int *indp,*pri,*dir;
//
//	// parameter
//	nrange = 10;
//
//	// Check how to split the variable set
//	//mincost = +CPX_INFBOUND;
//	//maxcost = -CPX_INFBOUND;
//	//for (j=0; j<ncols; j++)
//	//{
//	//	if (obj[j]<mincost)	mincost = obj[j];
//	//	if (obj[j]>maxcost)	maxcost = obj[j];
//	//}
//
//	//// Linear approach
//	//delta = (maxcost - mincost) / nrange;
//
//	//// Logaritmic approach
//	////...
//
//	// Build FCM Costraints
//	flag = (int *) malloc (ncols * sizeof (int));			if (flag==NULL) goto ESCI; 
//	ind = (int *) malloc ((ncols+nrange) * sizeof (int));			if (ind==NULL) goto ESCI; 
//	val = (double *) malloc ((ncols+nrange) * sizeof (double));		if (val==NULL) goto ESCI;
//	beg = (int *) malloc ((nrange+4) * sizeof (int));		if (beg==NULL) goto ESCI;
//	sense = (char *) malloc (nrange * sizeof (char));		if (sense==NULL) goto ESCI;
//	rhs = (double *) malloc (nrange * sizeof (double));		if (rhs==NULL) goto ESCI;
//	objt = (double *) malloc ((ncols+nrange) * sizeof (double));	if (objt==NULL) goto ESCI;
//	ctype = (char *) malloc ((ncols+nrange) * sizeof (char));		if (ctype==NULL) goto ESCI;
//	indt = (int *) malloc ((ncols+nrange) * sizeof (int));	if (indt==NULL) goto ESCI; 
//	indp = (int *) malloc ((ncols+nrange) * sizeof (int));	if (ind==NULL) goto ESCI; 
//	pri = (int *) malloc ((ncols+nrange) * sizeof (int));	if (ind==NULL) goto ESCI; 
//	dir = (int *) malloc ((ncols+nrange) * sizeof (int));	if (ind==NULL) goto ESCI; 
//
//	//srand(0);
//	for (j=0; j<ncols; j++)
//	{
//		flag[j] = 0;
//		objt[j] = obj[j];
//		ctype[j] = xctype[j];
//		indt[j] = j;
//
//		indp[j] = j;
//		pri[j] = 1;
//		dir[j] = 0;
//		//dir[j] = CPX_BRANCH_GLOBAL;
//		//dir[j] = CPX_BRANCH_UP; 
//
//	}
//
//	//level0 = mincost;
//	//level = mincost + delta;
//	k = 0;
//	for (i=0; i<nrange; i++)
//	{
//		mincost = +GRB_INFINITY;
//		maxcost = -GRB_INFINITY;
//		for (j=0; j<ncols; j++)
//		{
//			if (flag[j]>0)
//				continue;
//
//			if (obj[j]<mincost)	mincost = obj[j];
//			if (obj[j]>maxcost)	maxcost = obj[j];
//		}
//
//		// Linear approach
//		delta = (maxcost - mincost) / (nrange-i);
//
//		// Logaritmic approach
//		//...
//
//		level0 = mincost;
//		level = mincost + delta + 0.001;
//
//		rhs[i] = 0.0;
//		sense[i] = 'E';
//		beg[i] = k;
//		for (j=0; j<ncols; j++)
//		{
//			if (flag[j]>0)
//				continue;
//
//			if (((xctype[j]=='I')||(xctype[j]=='B'))&&(obj[j]<level))
//			{
//				flag[j] = 1;
//				objt[j] = obj[j] - level0;
//				ind[k] = j;
//				val[k] = 1.0;
//				k++;
//
//				indp[j] = j;
//				pri[ncols+i] = i+2;
//				//dir[ncols+i]=CPX_BRANCH_GLOBAL;
//				//dir[ncols+i] = CPX_BRANCH_DOWN; 
//				dir[ncols+i] = 0; 
//
//			}
//		}
//		ind[k] = ncols+i;
//		val[k] = -1.0;
//		k++;
//		objt[ncols+i] = level0;
//		indt[ncols+i] = ncols+i;
//		ctype[ncols+i] = 'I';
//		indp[ncols+i] = ncols+i;
//		pri[ncols+i] = 2*nrange;
//		//dir[ncols+i]=CPX_BRANCH_GLOBAL;
//		//dir[ncols+i] = CPX_BRANCH_DOWN; 
//		dir[ncols+i] = 0; 
//
//		level0 = level;
//		level += delta;
//
//	}
//	beg[nrange]=k;
//
//	// Load the FCM constraint
//	type = xctype-nrange;
//	status = GRBaddvars(lp,nrange,0,NULL,NULL,NULL,NULL,NULL,NULL,type,NULL);
//	if (status) return status;
//
//	status = GRBaddconstrs(lp,nrange,k,beg,ind,val,sense,rhs,NULL);
//	if (status) return status;
//
//	// Change the type for the new variable
//	//cnt = ncols + nrange;
//	//status = CPXchgctype(env,lp,cnt,indt,ctype);
//	//status = GRBsetcharattrarray(lp,"VType",ncols,nrange,ctype);
//
//	cnt = ncols + nrange;
//	status = GRBsetdblattrarray(lp,"Obj",0,cnt,objt);
//
//	status = GRBupdatemodel(lp);
//	if (status) return status;
//
//	// Load Priority Order 
//	//status = CPXcopyorder(env,lp,nrange,indp,pri,dir);
//
//	goto TERMINATE;
//
//ESCI:
//	fprintf (stderr, "No memory for solution values.\n");
//
//TERMINATE:
//	free(beg);
//	free(ind);
//	free(val);
//	free(sense);
//	free(rhs);
//	free(objt);
//	free(flag);
//	free(indt);
//	free(ctype);
//	free(indp);
//	free(pri);
//	free(dir);

	return status;
}


//-----------------------------------------------------------------------------
//  SolveMIP
//-----------------------------------------------------------------------------
//
int CoinObj::SolveMIP(double *Zopt, double *Zlb, double *Gap, long *Nodes, long *Cuts)
{
	int i,j;
	int ncuts;
	int solstat;
	double RGap;
	double aux;

	//status = CPXchgprobtype(env,lp,CPXPROB_MILP);
	//status = CPXchgctype(env,lp,ncols,indices,xctype);

	//status = CPXwriteprob(env,lp,"myprob.lp","LP");
	//status = CPXmipopt(env,lp);

	// Print MPS
	//lp.writeMps("Coin.mps",0,0,2,false);

	//aux=lp.optimizationDirection();

	// Load and Solve the model
	env.loadFromCoinModel(lp);
	CbcModel model(env);

	//***
	// Setup

	//model.solver()->setHintParam(OsiDoReducePrint,true,OsiHintTry);
	model.solver()->setHintParam(OsiDoReducePrint,true,OsiHintDo);

	// Set up some cut generators and defaults
	// Probing first as gets tight bounds on continuous

	CglProbing generator1;
	generator1.setUsingObjective(true);
	//generator1.setMaxPass(3);
	generator1.setMaxPass(30);
	generator1.setMaxProbe(100);
	generator1.setMaxLook(50);
	generator1.setRowCuts(3);
	// generator1.snapshot(*model.solver());
	// generator1.createCliques(*model.solver(),2,1000,true);
	// generator1.setMode(0);

	CglGomory generator2;
	// try larger limit
	//generator2.setLimit(300);
	generator2.setLimit(100);

	CglKnapsackCover generator3;

	CglOddHole generator4;
	generator4.setMinimumViolation(0.005);
	generator4.setMinimumViolationPer(0.00002);
	// try larger limit
	generator4.setMaximumEntries(200);

	CglClique generator5;
	generator5.setStarCliqueReport(false);
	generator5.setRowCliqueReport(false);

	CglFlowCover flowGen;
	CglMixedIntegerRounding mixedGen;
	//CglMixedIntegerRounding2 mixedGen2;
    //CglTwomir twomir;

	//CglAllDifferent generator6;
	//generator6.setMaxLook(50);
	//CglRedSplit generator6;
	//CglLiftAndProject generator6;
	//generator6.setBeta(-1);

	// Add in generators
	model.addCutGenerator(&generator1,-1,"Probing");
	model.addCutGenerator(&generator2,-1,"Gomory");
	//model.addCutGenerator(&generator3,-1,"Knapsack");
	model.addCutGenerator(&generator4,-1,"OddHole");
	model.addCutGenerator(&generator5,-1,"Clique");
	//model.addCutGenerator(&flowGen,-1,"FlowCover");
	//model.addCutGenerator(&mixedGen,-1,"MixedIntegerRounding");
	//model.addCutGenerator(&mixedGen2,-1,"MixedIntegerRounding2");
	//model.addCutGenerator(&twomir,-1,"TwoMir");

	//model.addCutGenerator(&generator6,-1,"LiftandProject");
	
	OsiClpSolverInterface * osiclp = dynamic_cast< OsiClpSolverInterface*> (model.solver());
	// go faster stripes
	if (osiclp->getNumRows()<300&&osiclp->getNumCols()<500) 
	{
		osiclp->setupForRepeatedUse(2,0);
		printf("trying slightly less reliable but faster version (? Gomory cuts okay?)\n");
		printf("may not be safe if doing cuts in tree which need accuracy (level 2 anyway)\n");
	}

	// Allow rounding heuristic
	//CbcRounding heuristic1(model);
	//model.addHeuristic(&heuristic1);

	// And local search when new solution found
	//CbcHeuristicLocal heuristic2(model);
	//model.addHeuristic(&heuristic2);

	//CbcHeuristicRINS heuristic3(model);
	//model.addHeuristic(&heuristic3);

	//CbcHeuristicFPump heuristic4(model);
	//model.addHeuristic(&heuristic4);

	//CbcHeuristicGreedyCover heuristic5(model);
	//CbcHeuristicGreedyEquality heuristic5(model);
	//CbcHeuristicGreedySOS heuristic5(model);
	//model.addHeuristic(&heuristic5);

	//CbcHeuristicPivotAndFix heuristic6(model);
	//model.addHeuristic(&heuristic6);

	//CbcHeuristicDiveCoefficient heuristic7(model);
	//CbcHeuristicDiveFractional heuristic7(model);
	//CbcHeuristicDiveVectorLength heuristic7(model);
	//model.addHeuristic(&heuristic7);

	// Redundant definition of default branching (as Default == User)
	//CbcBranchUserDecision branch;
	//CbcBranchDefaultDecision branch;
	//CbcBranchDynamicDecision branch;
	//model.setBranchingMethod(&branch);

	// Definition of node choice
	//CbcCompareUser compare;
	//model.setNodeComparison(compare);
	//CbcCompareDepth compare;
	//CbcCompareObjective compare;
	CbcCompareDefault compare;
	//CbcCompareEstimate compare;
	model.setNodeComparison(compare);

	// Do initial solve to continuous
	model.initialSolve();

	// Could tune more
	model.setMinimumDrop(CoinMin(1.0,fabs(model.getMinimizationObjValue())*1.0e-3+1.0e-4));

	if (model.getNumCols()<500)
		model.setMaximumCutPassesAtRoot(-100); // always do 100 if possible
	else if (model.getNumCols()<5000)
		model.setMaximumCutPassesAtRoot(100); // use minimum drop
	else
		model.setMaximumCutPassesAtRoot(20);
		//model.setMaximumCutPasses(5);

	// Switch off strong branching if wanted
	// model.setNumberStrong(0);
	// Do more strong branching if small
	//if (model.getNumCols()<5000)
	//	model.setNumberStrong(10);
	model.setNumberStrong(10);

	//model.solver()->setIntParam(OsiMaxNumIterationHotStart,100);
	model.solver()->setIntParam(OsiMaxNumIterationHotStart,10);

    model.setDblParam(CbcModel::CbcMaximumSeconds,TimeLimit);

	model.setPrintFrequency(50);
  
	// Switch off most output
	if (model.getNumCols()<3000) 
	{
		model.messageHandler()->setLogLevel(1);
		//model.solver()->messageHandler()->setLogLevel(0);
	} 
	else 
	{
		model.messageHandler()->setLogLevel(2);
		model.solver()->messageHandler()->setLogLevel(1);
	}

	double time1 = CoinCpuTime();

	// Do complete search

	//***

	// Solve
	model.initialSolve();
	model.branchAndBound();

	//status = GRBoptimize(lp);
	//if (status) return status;

	solstat = model.status();

	//status = GRBgetintattr(lp,"Status",&solstat);
	//if (status) return status;

	//if ((solstat!=GRB_OPTIMAL)&&(status!=GRB_NODE_LIMIT)&&(solstat!=GRB_TIME_LIMIT)) 
	if (solstat!=0)  // Da migliorare!!! 
	{
		printf("Non Ammissibile...\n",status);
		objval=+0.0;
		RGap=1.0;
		goto TERMINATE;
	}

	objval = model.getObjValue();

	// Read the solution
	//status = GRBgetdblattr(lp,"ObjVal",&objval);
	//if (status) return status;

	//status = GRBgetdblattr(lp,"ObjBound",Zlb);
	//if (status) return status;
	// Per ora
	*Zlb = objval;
	RGap = minmax*(objval-(*Zlb))/objval;

	const double *solution = model.solver()->getColSolution();
	for (j=0; j<ncols; j++)
	{
		xt[j] = solution[j];
	}

	//status = GRBgetdblattrarray(lp,"X",0,ncols,xt);
	//if (status) return status;

TERMINATE:

	// Ouput
	*Zopt = objval;
	*Gap = RGap * 100.;

	//status = GRBgetdblattr(lp,"NodeCount",&aux);
	//if (status) return status;
	//*Nodes=(long)aux;
	*Nodes = model.getNodeCount();

	//status = CPXgetnumcuts(env,lp,CPX_CUT_USER,&ncuts);
	//status = CPXgetnumcuts(env,lp,CPX_CUT_COVER,&ncuts);
	ncuts = 0;
	*Cuts = ncuts;

	return status;
}


//-----------------------------------------------------------------------------
//  Setup Custom Cutting Plane
//-----------------------------------------------------------------------------
//
int CoinObj::CuttingPlane(long n, long m, double* u, double* p, double* pmax, 
	                       long *nY, long **Y, long *nUOR, long **UOR, long *nUAND, long **UAND, long *FDepA, long *FDepP,
	                       int Pred, int Cover, int Lifting, int KnapSol, int MaxCuts, long *Cuts)
{
	int aux;
	long cur_numcols;

//	// Set parameters
//
//	//// Assure linear mappings between the presolved and original models 
//	//status = CPXsetintparam (env, CPX_PARAM_PRELINEAR, 0);
//	//if ( status )  goto TERMINATE;
// //  
//	//// Turn on traditional search for use with control callbacks 
//	//status = CPXsetintparam (env, CPX_PARAM_MIPSEARCH, CPX_MIPSEARCH_TRADITIONAL);
//	//if ( status )  goto TERMINATE;
//
//	//// Let MIP callbacks work on the original model
//	//status = CPXsetintparam (env, CPX_PARAM_MIPCBREDLP, CPX_OFF);
//	//if ( status )  goto TERMINATE;
//
//	//// Let MIP do not reduce LP
//	//status = CPXsetintparam (env, CPX_PARAM_REDUCE, CPX_PREREDUCE_NOPRIMALORDUAL);
//	//if ( status )  goto TERMINATE;
//
//	//// Let MIP do not presolve MIP
//	////status = CPXsetintparam (env, CPX_PARAM_PREIND, CPX_OFF);
//	////if ( status )  goto TERMINATE;
//
//	*Cuts=0;
//	status = GRBgetintattr(lp,"NumVars",&aux);
//	if (status) return status;
//	cur_numcols = aux;
//	cutinfo.FlagLP = 0;
//	cutinfo.lp = lp;
//	cutinfo.numcols = cur_numcols;
//	cutinfo.numtoadd = 0;
//	cutinfo.maxcuts = MaxCuts;
//
//	cutinfo.x = (double *) malloc (cur_numcols * sizeof (double));
//	if ( cutinfo.x == NULL ) {
//		fprintf (stderr, "No memory for solution values.\n");
//		goto TERMINATE;
//	}
//	cutinfo.ind = (int *) malloc (cur_numcols * sizeof (int));
//	if ( cutinfo.ind == NULL ) {
//		fprintf (stderr, "No memory for solution values.\n");
//		goto TERMINATE;
//	}
//	cutinfo.val = (double *) malloc (cur_numcols * sizeof (double));
//	if ( cutinfo.val == NULL ) {
//		fprintf (stderr, "No memory for solution values.\n");
//		goto TERMINATE;
//	}
//	cutinfo.beg = (int *) malloc (2 * sizeof (int));
//	if ( cutinfo.beg == NULL ) {
//		fprintf (stderr, "No memory for solution values.\n");
//		goto TERMINATE;
//	}
//	cutinfo.rhs = (double *) malloc (2 * sizeof (double));
//	if ( cutinfo.rhs == NULL ) {
//		fprintf (stderr, "No memory for solution values.\n");
//		goto TERMINATE;
//	}
//
//	// Set up to use MIP callback 
//	status = GRBsetcallbackfunc(lp,myNewcutcallback,&cutinfo);
//	if ( status )  goto TERMINATE;
//
//TERMINATE:
//
	return 0;
}


//-----------------------------------------------------------------------------
// Function implementing a custom cutting plane procedure called by callback
//-----------------------------------------------------------------------------
//
int __stdcall myNewcutcallback (OsiClpSolverInterface env,
               void       *cbdata,
               int        wherefrom,
               void       *cbhandle)
{
	long i;
	int status = 0;
	int BolAdd;
	char sense;
	int purgeable;
	int param;
	int ivalue;

//	GRBCUTINFOptr cutinfo = (GRBCUTINFOptr) cbhandle;
//
//	//*useraction_p = CPX_CALLBACK_DEFAULT; 
//
//	// Are we solving a subproblem?
//	if (cutinfo->FlagLP)
//		return (status);
//
//	// Do we have reached the maximum number of cuts?  
//	if (cutinfo->numtoadd>cutinfo->maxcuts)
//		return (status);
//
//	if (wherefrom == GRB_CB_MIPSOL_SOL) 
//	{
//		status = GRBcbget(cbdata,wherefrom,GRB_CB_MIPSOL_SOL,cutinfo->x);
//		if (status) 
//		{
//			fprintf(stderr, "Failed to get node solution.\n");
//			goto TERMINATE;
//		}
//	}
//	else 
//		return 0;
//
//	printf("\nSolution:\n");
//	for (i=0; i<cutinfo->numcols; i++)
//		if (cutinfo->x[i]>0.00001)
//			printf("x(%d)=%lf\n",i,cutinfo->x[i]);
//
//
//	// Use a cut violation tolerance of 0.01
//	//if (BolAdd==1) 
//	//{
//	//	cutinfo->numtoadd++;
//	//	cutinfo->Cuts[0]++;
//	//	if ((cutinfo->numtoadd%1000)==0) printf("\n  User Cuts added: %d\n\n",cutinfo->numtoadd);
//
//	//	//if ((cutinfo->numtoadd%5000)==0) 
//	//	//	printf("\n  User Cuts added: %d\n\n",cutinfo->numtoadd);
//
//	//	// status = CPXcutcallbackadd(env,cbdata,wherefrom,k1,rhs[0],'G',ind,val,0);
//	//	status = CPXcutcallbackadd(env,cbdata,wherefrom,cutinfo->nz,cutinfo->rhs[0],sense,cutinfo->ind,cutinfo->val,purgeable);
//	//	if (status)
//	//	{
//	//		fprintf (stderr, "Failed to add cut.\n");
//	//		goto TERMINATE;
//	//	}
//
//	//	// Tell CPLEX that cuts have been created
//	//	*useraction_p = CPX_CALLBACK_SET; 
//
//	//}
//
//TERMINATE:
//
	return (status);

} 




