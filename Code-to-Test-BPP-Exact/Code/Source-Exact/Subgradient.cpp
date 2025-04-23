//-----------------------------------------------------------------------------
// File: Subgradient.cpp
//
// Dual Ascent based on a Lagrangian Column Generation
// for Bin Packing Problems (BPPs)
//
// Authors: Marco A. Boschetti and S. Novellani
//
// Last update: 20.04.2025
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
// Include
//-----------------------------------------------------------------------------
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include "Utilities-2.h" 
#include "Subgradient.h"
#include "MIPSolverObj.h"

//-----------------------------------------------------------------------------
//  Constructor                                    
//-----------------------------------------------------------------------------
//
Subgradient::Subgradient(void)
{
	maxnz = 0;
	solstat = 0;

	obj = NULL;
	matval = NULL;
	matind = NULL;
	matbeg = NULL;
	sense = NULL;
	rhs = NULL;
	cardrhs = NULL;
	lb = NULL;
	ub = NULL;
	rmatval = NULL;
	rmatind = NULL;
	rmatbeg = NULL;
	val = NULL;
	ind = NULL;
	indl = NULL;
	indl1 = NULL;
	jmin = NULL;
	Qmin = NULL;
	CPmin = NULL;
	objpen = NULL;
	cardmat = NULL;

	rsense = NULL;
	rrhs = NULL;

	cflag = NULL;
	rflag = NULL;
	x = NULL;
	pen = NULL;
	pent = NULL;
	sub = NULL;
	subo = NULL;
	q = NULL;
	u = NULL;
	ubest = NULL;

	ferr = fopen("solheu.err", "w");
}

//-----------------------------------------------------------------------------
//  Destructor                                    
//-----------------------------------------------------------------------------
//
Subgradient::~Subgradient(void)
{
	if (obj != NULL) delete[] obj;
	if (matval != NULL) delete[] matval;
	if (matind != NULL) delete[] matind;
	if (matbeg != NULL) delete[] matbeg;
	if (sense != NULL) delete[] sense;
	if (rhs != NULL) delete[] rhs;
	if (cardrhs != NULL) delete[] cardrhs;
	if (lb != NULL) delete[] lb;
	if (ub != NULL) delete[] ub;
	if (rmatval != NULL) delete[] rmatval;
	if (rmatind != NULL) delete[] rmatind;
	if (rmatbeg != NULL) delete[] rmatbeg;
	if (val != NULL) delete[] val;
	if (ind != NULL) delete[] ind;
	if (indl != NULL) delete[] indl;
	if (indl1 != NULL) delete[] indl1;
	if (jmin != NULL) delete[] jmin;
	if (Qmin != NULL) delete[] Qmin;
	if (CPmin != NULL) delete[] CPmin;
	if (objpen != NULL) delete[] objpen;
	if (cardmat != NULL) delete[] cardmat;

	if (cflag != NULL) { delete[] cflag; cflag = NULL; }
	if (rflag != NULL) delete[] rflag;
	if (x != NULL) delete[] x;
	if (pen != NULL) delete[] pen;
	if (pent != NULL) delete[] pent;
	if (sub != NULL) delete[] sub;
	if (subo != NULL) delete[] subo;
	if (q != NULL) delete[] q;
	if (u != NULL) delete[] u;
	if (ubest != NULL) delete[] ubest;

	// Close lof file "logerr"
	fclose(ferr);

}

//-----------------------------------------------------------------------------
//  Malloc: allocate the data structures                                    
//-----------------------------------------------------------------------------
//
int Subgradient::Malloc(void)
{
	if (maxnz <= 0)
		maxnz = nmax * MaxColCard;

	obj = new double[nmax];
	matval = new double[maxnz];
	matind = new int[maxnz];
	matbeg = new int[nmax + 1];
	sense = new int[m];
	rhs = new double[m];
	if (mc > 0)
		cardrhs = new double[mc];
	lb = new double[nmax];
	ub = new double[nmax];
	val = new long[m];
	ind = new short[nmax];
	indl = new long[nmax];
	indl1 = new long[m];
	jmin = new long[m];
	Qmin = new double[m];
	CPmin = new double[m];
	objpen = new double[nmax];
	if (mc > 0)
		cardmat = new int[nmax*mc];

	// An option not used for the BPP
	if (Bol_SR3)
	{
		rmaxnz = rmax * 10000; // Estimate of non-zero coefficients
		rmatval = new double[rmaxnz];
		rmatind = new int[rmaxnz];
		rmatbeg = new int[rmax + 1];
		rsense = new int[rmax];
		rrhs = new double[rmax];
	}

	solstat = 0;
	objval = +Inf;
	cflag = new long[nmax];
	rflag = new long[m + mc];
	x = new double[nmax];
	pen = new double[m + mc];
	pent = new double[m + mc];
	sub = new double[m + mc];
	subo = new double[m + mc];
	q = new double[m + mc];
	u = new double[m + mc];
	ubest = new double[m + mc];

	return 0;
}

//-----------------------------------------------------------------------------
//  Free: free the memory allocated to data structures                                    
//-----------------------------------------------------------------------------
//
int Subgradient::Free(void)
{
	maxnz = 0;
	solstat = 0;

	if (obj != NULL) { delete[] obj; obj = NULL; }
	if (matval != NULL) { delete[] matval; matval = NULL; }
	if (matind != NULL) { delete[] matind; matind = NULL; }
	if (matbeg != NULL) { delete[] matbeg; matbeg = NULL; }
	if (sense != NULL) { delete[] sense; sense = NULL; }
	if (rhs != NULL) { delete[] rhs; rhs = NULL; }
	if (cardrhs != NULL) { delete[] cardrhs; cardrhs = NULL; }
	if (lb != NULL) { delete[] lb; lb = NULL; }
	if (ub != NULL) { delete[] ub; ub = NULL; }
	if (rmatval != NULL) { delete[] rmatval; rmatval = NULL; }
	if (rmatind != NULL) { delete[] rmatind; rmatind = NULL; }
	if (rmatbeg != NULL) { delete[] rmatbeg; rmatbeg = NULL; }
	if (val != NULL) { delete[] val; val = NULL; }
	if (ind != NULL) { delete[] ind; ind = NULL; }
	if (indl != NULL) { delete[] indl; indl = NULL; }
	if (indl1 != NULL) { delete[] indl1; indl1 = NULL; }
	if (jmin != NULL) { delete[] jmin; jmin = NULL; }
	if (Qmin != NULL) { delete[] Qmin; Qmin = NULL; }
	if (CPmin != NULL) { delete[] CPmin; CPmin = NULL; }
	if (objpen != NULL) { delete[] objpen; objpen = NULL; }
	if (cardmat != NULL) { delete[] cardmat; cardmat = NULL; }

	if (rsense != NULL) { delete[] rsense; rsense = NULL; }
	if (rrhs != NULL) { delete[] rrhs; rrhs = NULL; }

	if (cflag != NULL) { delete[] cflag; cflag = NULL; }
	if (rflag != NULL) { delete[] rflag; rflag = NULL; }
	if (x != NULL) { delete[] x; x = NULL; }
	if (pen != NULL) { delete[] pen; pen = NULL; }
	if (pent != NULL) { delete[] pent; pent = NULL; }
	if (sub != NULL) { delete[] sub; sub = NULL; }
	if (subo != NULL) { delete[] subo; subo = NULL; }
	if (q != NULL) { delete[] q; q = NULL; }
	if (u != NULL) { delete[] u; u = NULL; }
	if (ubest != NULL) { delete[] ubest; ubest = NULL; }

	return 0;
}

//-----------------------------------------------------------------------------
//  Close: free the memory allocated to data structures                                    
//-----------------------------------------------------------------------------
//
int Subgradient::Close(void)
{
	long i;

	//Free();

	// DeAllocate memory
	if (Inst.w)
	{
		delete[] Inst.w;
		Inst.w = NULL;
	}

	if (Inst.cl)
	{
		delete[] Inst.cl;
		Inst.cl = NULL;
	}

	if (Inst.CMatrix)
	{
		for (i = 0; i < Inst.n; i++)
			delete[] Inst.CMatrix[i];
		delete[] Inst.CMatrix;
		Inst.CMatrix = NULL;
	}

	if (Inst.w2)
	{
		delete[] Inst.w2;
		Inst.w2 = NULL;
	}

	return 0;
}


//-----------------------------------------------------------------------------
//  Function to allocate the Data Structure to save the model 
//  (also during column generation)
//-----------------------------------------------------------------------------
//
int Subgradient::AllocateDataStructure(void)
{
	long i, j, k;
	long kk;
	int status;

	// Setup the KP model
	KPSol.Inst.n = Inst.n;
	KPSol.Inst.W = Inst.W;
	KPSol.Inst.k = Param.MaxNumColumns;  // To allocate the maximum number of column generated by "k-best KP"
	KPSol.Inst.w = new long[KPSol.Inst.n];
	KPSol.Inst.p = new long[KPSol.Inst.n];      // It's not really necessary!
	KPSol.Inst.pf = new float[KPSol.Inst.n];    // To setup using dual variables when we need to generate columns 
	KPSol.Inst.b = new long[KPSol.Inst.n];      // It's not really necessary!
	KPSol.Cons.Type = Inst.Type;

	for (i = 0; i < KPSol.Inst.n; i++)
	{
		KPSol.Inst.w[i] = Inst.w[i];
		KPSol.Inst.p[i] = 0;       // It's not really necessary!
		KPSol.Inst.pf[i] = 0.0;
		KPSol.Inst.b[i] = 1;       // It's not really necessary!
	}

	// Additional constraints: initialize the required data structure
	if (Inst.Type == 1)
	{
		// Cardinality constrained bin packing problem
		KPSol.Cons.MaxCard = Inst.k;
	}
	else if (Inst.Type == 2)
	{
		// Class constrained bin packing problem
		KPSol.Cons.CL = Inst.CL;  // Maximum number of classes for each bin
		KPSol.Cons.cl = new long[KPSol.Inst.n];  // Classes of items
		for (i = 0; i < KPSol.Inst.n; i++)
		{
			KPSol.Cons.cl[i] = Inst.cl[i];
		}
	}
	else if (Inst.Type == 3)
	{
		// Bin packing problem with conflicts
		KPSol.Cons.CMatrix = new long*[KPSol.Inst.n];  // Conflit matrix
		for (i = 0; i < KPSol.Inst.n; i++)
		{
			KPSol.Cons.CMatrix[i] = new long[KPSol.Inst.n];
			for (j = 0; j < KPSol.Inst.n; j++)
				KPSol.Cons.CMatrix[i][j] = Inst.CMatrix[i][j];
		}
	}
	else if (Inst.Type == 4)
	{
		// Two-dimensional vector packing problem
		KPSol.Cons.W2 = Inst.W2;  // Second capacity of the bin
		KPSol.Cons.w2 = new long[KPSol.Inst.n];  // Second Weights of items
		for (i = 0; i < KPSol.Inst.n; i++)
		{
			KPSol.Cons.w2[i] = Inst.w2[i];
		}
	}

	KPSol.Allocate_KPBW_kBest_LDCF();

	// Setup model for BPP

	// Setup the numbers of columns and rows 
	zub = +Inf;
	n = Inst.n;  // We add the trivial solution (i.e., one item for each bin)
	m = Inst.n;  // We have a row for each item

	// Allocate the data structure
	n1 = n;  // There are no artificial variables
	n2 = 0;  // There are no heuristic solutions
	mc = 0;  // There are no cardinality constraints
	offset = 0; // All rows must be covered at least once
	nmax = Param.MaxNumColumns;
	rmax = 1000; // Maximum number of valid inequalities (not included yet... possible further development) 
	status = Malloc();
	if (status)
		return status;

	// Setup RHS and sense
	for (i = 0; i < m; i++)
	{
		rhs[i] = 1.0;
		if (Param.Model == 0)
			sense[i] = 0;  // Set partitioning constraint
		else
			sense[i] = -1;  // Set covering constraint
	}

	// Setup objective function and constraint matrix columns 
	kk = 0;
	for (j = 0; j < n; j++)
	{
		// objective function
		obj[j] = 1.0;  // Minimize bins

		// Setup the non-zero coefficient of column/cluster j
		matbeg[j] = kk;
		for (k = 0; k < 1; k++)
		{
			matind[kk] = j;
			matval[kk] = 1.0;
			kk++;
		}

		lb[j] = 0.0;
		ub[j] = 1.0;
	}

	matbeg[n] = kk;

	return 0;
}

//-----------------------------------------------------------------------------
//  Function to deallocate the data structures
//-----------------------------------------------------------------------------
//
int Subgradient::DeAllocateDataStructure(void)
{
	long i;
	int status;

	// Setup the KP model
	if (KPSol.Inst.w != NULL)
	{
		delete[] KPSol.Inst.w;
		KPSol.Inst.w = NULL;
	}
	if (KPSol.Inst.p != NULL)
	{
		delete[] KPSol.Inst.p;
		KPSol.Inst.p = NULL;
	}
	if (KPSol.Inst.pf != NULL)
	{
		delete[] KPSol.Inst.pf;
		KPSol.Inst.pf = NULL;
	}
	if (KPSol.Inst.b != NULL)
	{
		delete[] KPSol.Inst.b;
		KPSol.Inst.b = NULL;
	}

	// Additional constraints: deallocate memory for the required data structure
	if (Inst.Type == 2)
	{
		// Class constrained bin packing problem
		delete[] KPSol.Cons.cl;
		KPSol.Cons.cl = NULL;
	}
	else if (Inst.Type == 3)
	{
		// Bin packing problem with conflicts
		for (i = 0; i < KPSol.Inst.n; i++)
			delete[] KPSol.Cons.CMatrix[i];
		delete[] KPSol.Cons.CMatrix;
		KPSol.Cons.CMatrix = NULL;
	}
	else if (Inst.Type == 4)
	{
		// Two-dimensional vector packing problem
		delete[] KPSol.Cons.w2;
		KPSol.Cons.w2 = NULL;
	}

	// DeAllocate memory for the data structure for column generation
	KPSol.Inst.k = Param.MaxNumColumns;  // To deallocate the maximum number of column generated by "k-best KP"
	KPSol.DeAllocate_KPBW_kBest_LDCF();

	// DeAllocate the data structure for the model
	status = Free();
	if (status)
		return status;

	return 0;
}

//-----------------------------------------------------------------------------
//  Compare two values
//-----------------------------------------------------------------------------
//
int compare(const void * a, const void * b)
{
	return (*(int*)a - *(int*)b);
}

//-----------------------------------------------------------------------------
// Subgradient for the Set Partitioning with the Parametric Lagrangian 
// Relaxation
//-----------------------------------------------------------------------------
//
int Subgradient::Dual_Ascent_SetPar(long BolLagrHeu, long *BolOpt)
{
	int status = 0;
	int bol;
	long i, j, k, ii, jj, kk;
	long ncolred;
	long maxit, it;
	long maxiti, iti;
	long maxitf, itf;
	long NCore;
	long *Core;
	long *icov;
	long pstep;
	double objvalt;
	double den, alfa, beta, delta;
	double epsilon;
	double rho;
	double MinGap;
	double cpen;
	double fbest;
	double t1, t2, dt;
	double objf;
	double Qj;
	double gamma1, gamma2;
	double *val;
	double *dcov;
	double *preddir;
	double zgbest;
	double HGap;

	long maxKPk;
	long fbit = 0;

	int FlagCorollary = 0;
	double Mindual;

	// Get the starting time
	t1 = ((double)clock()) / CLOCKS_PER_SEC;

	// Parameters
	(*BolOpt) = 0;
	pstep = 100;  // 100
	fbest = -Inf;
	objf = -Inf;
	zgbest = ZGreedy;
	ncolred = 0;

	// Setup parameters
	maxit = 100000;  // Maximum number of iterations
	maxiti = 5;      // Maximum number of iterations without improvement before to modify the step
	maxitf = Param.MaxItMinGap;    // Maximum number of iterations without improvement larger than MinGap before to stop
	alfa = Param.alpha;  // Starting step
	//beta = 0.005;    // 0.005
	epsilon = 0.05;  // Penalized cost gap to eliminate variable from core 
	gamma1 = 0.90;   // Reduction factor for the step
	gamma2 = 1.005;  // Increment factor for the step
	delta = 0.01;    // Factor to approximate the gap using the lower bound (when updating the step)
	MinGap = Param.MinGap;  // The value of MinGap
	if (BolLagrHeu > 1)
		HGap = 0.01;  // Factor to decide when apply heuristic algorithms
	else
		HGap = 0.00001;  // Factor to decide when apply heuristic algorithms  
	rho = 0.5;  // Coefficient used "to correct" the direction of the subgradient
	maxKPk = (long)(Param.PercK * Inst.n);  // Maximum number of columns generated at each iteration
	if (maxKPk < 1)
		maxKPk = 1;

	// Allocate the array containing the "core"
	Core = new long[nmax];

	// Allocate a data structure used in the greedy heuristic
	icov = new long[m];
	dcov = new double[m];

	// Allocate the array used to sort columns
	val = new double[nmax];
	cred = new double[nmax];

	// Allocate the array used to save the direction computed in the previous iteration
	preddir = new double[m];

	// Set flag (cflag[j]=1 if j is generated) and index (indl[j]=jj if jj is th j-th column) for each column
	for (j = 0; j < n; j++)
	{
		cflag[j] = 0;
		indl[j] = j;
		val[j] = obj[j] / double(matbeg[j + 1] - matbeg[j] + 1);
	}

	// Setup weight q, penalty pen and flag associated to each row 
	srand(0);
	for (i = 0; i < m; i++)
	{
		ubest[i] = 0.0;
		q[i] = 1.0;  
		pen[i] = 0.0;
		rflag[i] = 0;
		preddir[i] = 0.0;
	}

	// Setup the first core of columns for starting the Lagrangian Column Generation
	NCore = 0;
	for (jj = 0; jj < n; jj++)
	{
		bol = 0;
		j = indl[jj];

		Core[NCore] = j;  
		NCore++;
		cflag[j] = 1;
		for (k = matbeg[j]; k < matbeg[j + 1]; k++)
		{
			i = matind[k];
			rflag[i]++;
		}
	}

	// The subgradient main loop
	itf = 0;
	for (it = 0; it < maxit; it++)
	{
		for (j = 0; j < n; j++)
		{
			x[j] = 0.0;
		}

		// For each row select the column which minimize u[i]=min{c'[j]/Q[j] : j \in Core}
		objvalt = 0.0;
		for (i = 0; i < m - mc; i++)
		{
			jmin[i] = -1;
			CPmin[i] = +Inf;
			Qmin[i] = 0.0;
		}

		// Evaluate the penalized cost of each variable in the core
		for (jj = 0; jj < NCore; jj++)
		{
			j = Core[jj];
			if (cflag[j] <= 0)
				continue;

			Qj = 0.0;
			cpen = 0.0;
			for (k = matbeg[j]; k < matbeg[j + 1]; k++)
			{
				ii = matind[k];
				cpen = cpen - matval[k] * pen[ii];
				if (ii < m - mc)
					Qj = Qj + matval[k] * q[ii];
			}
			cpen = cpen + obj[j];
			cpen = cpen / Qj;

			for (k = matbeg[j]; k < matbeg[j + 1]; k++)
			{
				ii = matind[k];
				if (ii < m - mc)
				{
					if (cpen < CPmin[ii])
					{
						jmin[ii] = j;
						CPmin[ii] = cpen;
						Qmin[ii] = Qj;
					}
				}
			}
		}

		// Compute dual variables (not necessarily feasible) and the value of the solution of the Lagrangian problem 
		for (i = 0; i < m - mc; i++)
		{
			if (jmin[i] != -1)
				u[i] = q[i] * CPmin[i] + pen[i];
			else
				u[i] = 0.0;
			objvalt += u[i] * rhs[i];
		}
		// "mc==1" iff we add the cardinality constraint (for BPP is useless but necessary for other problems)
		if (mc == 1)
		{
			u[m - 1] = +pen[m - 1]; // Side constraint
			objvalt += u[m - 1] * rhs[m - 1];
		}

		// Compute the subgradient 
		for (i = 0; i < m; i++)
			sub[i] = rhs[i];

		for (i = 0; i < m - mc; i++)
		{
			// If i is empty, then skip
			if (jmin[i] == -1)
				continue;

			j = jmin[i];
			for (k = matbeg[j]; k < matbeg[j + 1]; k++)
			{
				ii = matind[k];
				sub[ii] = sub[ii] - matval[k] * (q[i] / Qmin[i]);
			}
		}

		// Print the log
		if (it%pstep == 0)
		{
			if (BolLagrHeu == 0)
				printf(" %d) alfa = %8.5lf   objval = %12.5lf   fbest = %12.5lf    NCore = %d/%d \n", it, alfa, objvalt, fbest, NCore, n);
			else
				printf(" %d) alfa = %8.5lf   objval = %12.5lf   fbest = %12.5lf   zgbest = %12.5lf    NCore = %d/%d \n", it, alfa, objvalt, fbest, zgbest, NCore, n);
		}

		// Heuristic algorithms
		if ((BolLagrHeu > 0) && (it > 10) && ((fbest - objvalt) / fbest < HGap))
		{
			/*
			... add Heuristic Algorithms to obtain a "good" feasible solution
			*/
		}

		// Is the lower bound improved?
		if (fbest < objvalt)
		{
			bol = 1;

			// Setup KP to generate columns
			for (ii = 0; ii < Inst.n; ii++)
				KPSol.Inst.pf[ii] = (float)u[ii];
			KPSol.Inst.k = maxKPk; // 10;
			KPSol.Solve_KPBW_kBest_LDCF_CG(0, (float)1.0001);
			if (KPSol.kBSol.nk > 0)
			{
				// Check nmax
				if (n >= nmax - 1)
					return 1;

				// Add columns
				jj = matbeg[n];
				for (kk = 0; kk < KPSol.kBSol.nk; kk++)
				{
					if (KPSol.kBSol.zf[kk] < 1.0 + Prec)
						break;
					bol = 0;
					Core[NCore] = n;
					NCore++;
					cflag[n] = 1;
					indl[n] = n;
					obj[n] = 1.0;
					lb[n] = 0.0;
					ub[n] = 1.0;
					for (ii = KPSol.kBSol.nx[kk] - 1; ii >= 0; ii--)
					{
						matind[jj] = KPSol.kBSol.x[kk][ii];
						matval[jj] = 1.0;
						jj++;
					}
					n++;
					matbeg[n] = jj;
				}
			}

			// Check if the dual solution is feasible (bol==0 -> not feasible)
			if (bol == 0)
				continue;
			else
			{
				fbest = objvalt;
				alfa *= gamma2;
				iti = 0;

				// Save the best dual solution
				objdual = fbest;
				for (i = 0; i < m; i++)
					ubest[i] = u[i];

				// Delete columns with a too higher reduced cost (delete means flag "cflag[j]=0"... do not consider)
				for (jj = 0; jj < NCore; jj++)
				{
					j = Core[jj];
					if (cflag[j] <= 0)
						continue;

					bol = 1;
					cpen = obj[j];
					for (k = matbeg[j]; k < matbeg[j + 1]; k++)
					{
						ii = matind[k];
						cpen = cpen - matval[k] * u[ii];
						if (rflag[ii] <= 1)
						{
							bol = 0;
							break;
						}
					}
					objpen[j] = cpen;

					if ((bol) && (cpen > epsilon*fbest))
					{
						cflag[j] = 0;
						for (k = matbeg[j]; k < matbeg[j + 1]; k++)
						{
							ii = matind[k];
							rflag[ii]--;
						}
					}
				}

				kk = 0;
				for (jj = 0; jj < NCore; jj++)
				{
					j = Core[jj];
					if (cflag[j])
					{
						Core[kk] = Core[jj];
						kk++;
					}
				}
				NCore = kk;

				// Apply Corollary if enabled 
				if (FlagCorollary)
				{
					// Save the best dual solution
					Mindual = ubest[0];
					for (i = 1; i < m; i++)
						if (ubest[i] < Mindual)
							ubest[i] = Mindual;
					Mindual = Mindual - 1;
					for (i = 0; i < m; i++)
					{
						q[i] = ubest[i] - Mindual;
						pen[i] = Mindual;
					}
				}
			}
		}
		else
		{
			if (iti > maxiti)
			{
				alfa *= gamma1;
				iti = 0;
			}
			else
				iti++;
		}

		// Stop criteria n.1 (i.e., (z_UB - z_LB) < MinGap)
		if ((BolLagrHeu > 0) && ((((zgbest - fbest) / fbest) < MinGap) || ((zgbest - fbest) < 0.999)))
			goto TERMINATE;

		// Stop criteria n.2 (Lagrangian conditions)
		if (((fbest - objf) / fbest) < MinGap)
		{
			if (itf > maxitf)
				goto TERMINATE;
			else
				itf++;
		}
		else
		{
			itf = 0;
			objf = fbest;
		}

		// Update penalties
		den = 0.0;
		for (i = 0; i < m; i++)
		{
			den += sub[i] * sub[i];
		}

		if (den < Prec)
		{
			(*BolOpt) = 1;
			break;
		}

		beta = alfa * (delta*objvalt) / den;
		if (beta < -Prec)
			beta = -beta;
		if (beta < Prec)
			beta = Prec;

		for (i = 0; i < m; i++)
		{
			// If i is empty, then skip
			if (i < m - mc)
			{
				if (jmin[i] == -1)
					continue;
			}

			pen[i] = pen[i] + rho * beta * sub[i] + (1.0 - rho)*preddir[i];  
			if (sense[i] > 0)
				if (pen[i] > Prec)
					pen[i] = 0.0;
				else if (sense[i] < 0)
					if (pen[i] < -Prec)
						pen[i] = 0.0;

			preddir[i] = beta * sub[i];  
		}
	}

TERMINATE:

	t2 = ((double)clock()) / CLOCKS_PER_SEC;
	dt = t2 - t1;

	printf(" NcolRed = %d \n", ncolred);
	printf(" Best objval = %lf \n", fbest);
	printf(" Time = %lf sec.\n", dt);

	// Setup output
	ZLBlagr = fbest;
	ZUBlagr = zgbest;
	alpha = alfa;
	Iter = it;
	nCore = 0;
	for (j = 0; j < n; j++)
		if (cflag[j] > 0)
			nCore++;

	// Deallocate the array containing the "core"
	delete[] Core;
	delete[] icov;
	delete[] dcov;
	delete[] val;
	delete[] cred;
	delete[] preddir;

	return (status);
}

//-----------------------------------------------------------------------------
// Solve the LP-relaxation or MIP for the basic BPP
//-----------------------------------------------------------------------------
//
int Subgradient::Solve_Model(long flaglp, double *soldx, long *niter)
{
	int status;
	long i, j;
	long ii, jj, kk;

	MIPSolverObj mip; // MIP Instance

	ZUBlp = 100000; // +Inf;
	ZLBlp = 0.0;
	status = 0;
	*niter = 0;

	// Reset the solver
	status = mip.Open(0);
	if (status)
		return status;

	// Allocate data structures to enter the instance by rows
	mip.ncols = Inst.n * Inst.n + Inst.n;
	mip.nrows = Inst.n + Inst.n;
	mip.nz = 2 * Inst.n * Inst.n + Inst.n;
	status = mip.MallocRows();
	if (status)
		return status;

	// Load data structure
	mip.minmax = 1;   // Min
	mip.objsen = 1;   // Min
	mip.probtype = 1;  // 0=LP; 1=MIP

	//Setup counters
	ii = 0;  // Row counter
	jj = 0;  // Column counter
	kk = 0;  // NonZero counter

	// Variables x
	for (i = 0; i < Inst.n; i++)
		for (j = 0; j < Inst.n; j++)
		{
			mip.obj[jj] = 0.0;
			mip.lb[jj] = 0;
			mip.ub[jj] = 1;
			mip.indices[jj] = jj;
			if (flaglp > 0)
				mip.xctype[jj] = 'B';
			else
				mip.xctype[jj] = 'C';
			jj++;
		}

	// Variables y
	for (j = 0; j < Inst.n; j++)
	{
		mip.obj[jj] = 1.0;
		mip.lb[jj] = 0;
		mip.ub[jj] = 1;
		mip.indices[jj] = jj;
		if (flaglp > 0)
			mip.xctype[jj] = 'B';
		else
			mip.xctype[jj] = 'C';
		jj++;
	}

	// Constraints Capacity
	for (i = 0; i < Inst.n; i++)
	{
		mip.rhs[ii] = 0.0;
		mip.sense[ii] = 'L';
		mip.rmatbeg[ii] = kk;

		for (j = 0; j < Inst.n; j++)
		{
			mip.rmatind[kk] = i * Inst.n + j;
			mip.rmatval[kk] = Inst.w[j];
			kk++;
		}
		mip.rmatind[kk] = Inst.n * Inst.n + i;
		mip.rmatval[kk] = -Inst.W;
		kk++;

		ii++;
	}

	// Constraints Assignment
	for (j = 0; j < Inst.n; j++)
	{
		mip.rhs[ii] = 1;
		mip.sense[ii] = 'E';
		mip.rmatbeg[ii] = kk;

		for (i = 0; i < Inst.n; i++)
		{
			mip.rmatind[kk] = i * Inst.n + j;
			mip.rmatval[kk] = 1.0;
			kk++;
		}

		ii++;
	}
	mip.rmatbeg[ii] = kk;
	mip.nz = kk;

	// Load Problem
	status = mip.CopyLPbyRows();
	if (status)
		return status;

	// Optimize
	if (flaglp > 0)
	{
		mip.SetMIP(Param.TimeLimit);  // Set MIP
		status = mip.SolveMIP(&ZUBlp, &ZLBlp, &Gap, &Nodes, &Cuts);
		if (status)
			return status;
		if ((soldx != NULL) && (mip.ncols < MaxColGen))
		{
			for (j = 0; j < mip.ncols; j++)
				soldx[j] = mip.xt[j];
		}
		else
			status = 1;
	}
	else
	{
		mip.lpsolver = -flaglp;
		status = mip.SolveLP(&ZLBlp, soldx, niter);
		if (status)
			return status;
	}

	if (ZUBlp > 100000)
		ZUBlp = 100000; // +Inf;

	return status;
}

//-----------------------------------------------------------------------------
// Solve the LP-relaxation
//-----------------------------------------------------------------------------
//
int Subgradient::Solve_Incremental_Core(long flaglp, long flaggreedy, double *soldx, long *iter, long *BolOpt)
{
	int status;
	long i, j, k;
	long ii, jj, kk, nn;
	long nz;
	long NoOpt;
	double cpen;
	double Delta;
	double RDelta;
	double LBound;
	double Gamma1 = 0.1; 
	double Gamma2 = 10.;
	long BolCRed = 0;
	long MaxIter = 2;
	long MaxCol = Param.MaxNumColMIP;

	long t1, t2;
	double dt;
	double MaxTime;

	long kini = 0;
	int bol;

	MIPSolverObj mip; // MIP Instance

	//Setup executon time
	MaxTime = Param.TimeLimit - TimeDA;
	if (MaxTime < 0.1)
		MaxTime = 0.0;

	// Start clock
	t1 = clock();

	ZUBlp = +Inf;
	ZLBlp = 0.0;
	status = 0;
	(*iter) = 1;
	(*BolOpt) = 0;
	nCore = 0;

	// Compute Gap
	Gamma1 = Param.PercDelta;
	Delta = 0.0;
	double minu = +1000000.;
	double maxu = -1000000.;
	for (i = 0; i < m; i++)
	{
		if (minu > ubest[i]) minu = ubest[i];
		if (maxu < ubest[i]) maxu = ubest[i];
		Delta += rhs[i] * ubest[i];
	}
	LBound = Delta;
	if (flaglp > 0)
	{
		Delta = Gamma1 * Delta;
	}
	else
	{
		Delta = -0.00001;
	}

loop:

	if (flaggreedy)
	{
		if (ZUBlp > ZGreedy)
			ZUBlp = ZGreedy;
		if (ZUBlp > ZUBlagr)
			ZUBlp = ZUBlagr;
	}

	// Delta must not exceed the upper bound
	if (Delta > ZUBlp + Prec)
		Delta = ZUBlp + Prec;

	// If ZLBlp-ZLBlp<1, then stop... ZUBlp is the value of the optimal solution!!!
	if (ZUBlp < 1000)
	{
		if ((ZUBlp - LBound) < 1 - Prec)
		{
			(*BolOpt) = 1;
			return status;
		}
	}
	else
	{
		if ((fabs(ZUBlp - LBound) / ZUBlp) < Prec)
		{
			(*BolOpt) = 1;
			return status;
		}
	}

	// Compute core
	bol = 1;
	for (j = 0; j < n; j++)
	{
		indl[j] = j;
		cflag[j] = -1;
		cpen = 0.0;
		for (k = matbeg[j]; k < matbeg[j + 1]; k++)
		{
			ii = matind[k];
			cpen = cpen - matval[k] * ubest[ii];
		}
		cpen = cpen + obj[j];
		objpen[j] = cpen;
	}

	QSortD(objpen, indl, 0, n - 1, -1);

	nn = 0;
	nz = 0;
	for (jj = 0; jj < n; jj++)
	{
		j = indl[jj];
		if (nn >= MaxCol)
			break;
		//if (objpen[j] > Delta + 1000.)
		if ((objpen[j] > Delta) && (flaglp > 1))  // Temp 12.2023
		{
			RDelta = Delta;
			bol = 0;
			break;
		}
		RDelta = objpen[j];
		cflag[j] = nn;
		nn++;
		nz += (matbeg[j + 1] - matbeg[j]);
	}
	RDelta = Delta;

	// Generate new columns
	// Check if any column must be added to the core 
	// Here we must include "column generation"
	// Setup KP to generate columns
	if ((nn < MaxCol) && (flaglp > 0))
	{
		for (ii = 0; ii < Inst.n; ii++)
			KPSol.Inst.pf[ii] = (float)ubest[ii];
		KPSol.Inst.k = MaxCol - n;
		KPSol.Solve_KPBW_kBest_LDCF_CG(kini, (float)(1.0 + Prec - Delta));
		if (KPSol.kBSol.nk > 0)
		{
			// Check if all columns are generated
			if (KPSol.kBSol.nk == KPSol.Inst.k)
				bol = 0;

			// Add columns
			jj = matbeg[n];
			for (kk = 0; kk < KPSol.kBSol.nk; kk++)
			{
				if (1.0 + Prec - KPSol.kBSol.zf[kk] > Delta)
				{
					bol = 0;
					break;
				}
				RDelta = 1.0 - KPSol.kBSol.zf[kk];
				//Core[NCore] = n;
				//NCore++;
				cflag[n] = nn;
				nn++;
				obj[n] = 1.0;
				lb[n] = 0.0;
				ub[n] = 1.0;
				for (ii = KPSol.kBSol.nx[kk] - 1; ii >= 0; ii--)
				{
					matind[jj] = KPSol.kBSol.x[kk][ii];
					matval[jj] = 1.0;
					jj++;
				}
				n++;
				matbeg[n] = jj;
				nz += (matbeg[n] - matbeg[n - 1]);
			}
		}
	}
	//RDelta = double(long(LBound + RDelta + 0.9999)) + 0.0001;

	printf("  Colonne=%d \n", nn);
	printf("  LBound=%lf \n", LBound);
	printf("  Delta=%lf \n", Delta);
	printf("  RDelta=%lf \n", RDelta);

	// Do we generate all columns?
	if (bol == 1)
		RDelta = +(double)Infi;

	// Reset the solver
	mip.Open(0);

	// Allocate data structures to enter the instance by rows
	nCore = nn;
	mip.ncols = nn;
	mip.nrows = m;
	mip.nz = nz;
	mip.MallocCols();

	if (Bol_SR3)
	{
		mip.Malloc_Valid_Inequalities(1000);
	}

	// Load data structure
	mip.minmax = 1;    // Min
	mip.objsen = 1;    // Min
	mip.probtype = 0;  // 0=LP; 1=MIP

	// Variables 
	kk = 0;
	j = 0;
	for (jj = 0; jj < n; jj++)
	{
		if (cflag[jj] > -1)
		{
			cflag[jj] = j;
			if (BolCRed)
				mip.obj[j] = objpen[jj];
			else
				mip.obj[j] = obj[jj];
			mip.lb[j] = 0;
			mip.ub[j] = 1;
			mip.indices[j] = j;
			if (flaglp > 0)
				mip.xctype[j] = 'B';
			else
				mip.xctype[j] = 'C';

			mip.matbeg[j] = kk;
			mip.matcnt[j] = matbeg[jj + 1] - matbeg[jj];
			for (k = matbeg[jj]; k < matbeg[jj + 1]; k++)
			{
				mip.matind[kk] = matind[k];
				mip.matval[kk] = matval[k];
				kk++;
			}
			j++;
		}
	}
	mip.matbeg[j] = kk;

	// Constraints
	for (i = 0; i < m; i++)
	{
		mip.rhs[i] = rhs[i];
		if (sense[i] == 0)
			mip.sense[i] = 'E';
		else if (sense[i] == 1)
			mip.sense[i] = 'L';
		else if (sense[i] == -1)
			mip.sense[i] = 'G';
		else
			return 1215;
		if (flaglp)
			mip.sense[i] = 'G'; // Apply the Set Covering everytime... Temp!!! !!! !!!
		else
			mip.sense[i] = 'E'; // Apply the Set Partitining everytime... Temp!!! !!! !!!
	}

	// Load problem
	status = mip.CopyLPbyCols();

	// Optimize
	if (flaglp > 0)
	{
		NoOpt = 1;

		if (MaxTime > Prec)
		{
			mip.SetMIP(MaxTime);  // Set MIP
			mip.SolveMIP(&ZUBlp, &ZLBlp, &Gap, &Nodes, &Cuts);
		}
		else
		{
			ZUBlp = +Inf;
			ZLBlp = -Inf;
			Nodes = 0;
			Cuts = 0;
		}
		printf("\n\n MIP:\n");
		printf("   ZUBlp = %lf \n", ZUBlp);
		printf("   LBound = %lf \n", LBound);
		printf("   ZUBlp-LBound = %lf \n", ZUBlp - LBound);
		printf("   RDelta = %lf \n\n", RDelta);

		if ((fabs(ZUBlp - ZLBlp) / ZUBlp) < Prec)
			if (BolCRed)
			{
				if ((ZUBlp / ZUBlp) < ((RDelta / ZUBlp) + Prec))  // Puo' essere anche molto negativo se RDelta>>ZUBlp
				{
					NoOpt = 0;
					(*BolOpt) = 1;
				}
			}
			else
			{
				if (ZUBlp < (Prec + (double)(long(LBound + RDelta + 1.0 - Prec))))
				{
					NoOpt = 0;
					(*BolOpt) = 1;
				}
			}

		if ((NoOpt) && ((*iter) < MaxIter))
		{
			t2 = clock();
			dt = (double)(t2 - t1) / CLOCKS_PER_SEC;
			if (dt < MaxTime)
			{
				Delta = Gamma2 * Delta;
				MaxCol = long(Gamma2*MaxCol);
				(*iter)++;
				goto loop;
			}
		}

		if (BolCRed)
		{
			ZUBlp = ZUBlp + LBound;
			ZLBlp = ZLBlp + LBound;
			Gap = 100 * fabs(ZUBlp - ZLBlp) / ZUBlp;
		}

		if (soldx != NULL)
		{
			for (j = 0; j < n; j++)
			{
				if (cflag[j] > -1)
					soldx[j] = mip.xt[cflag[j]];
				else
					soldx[j] = 0.0;
			}
		}
	}
	else
	{
		long NIter;
		mip.lpsolver = -flaglp;
		mip.SolveLP(&ZLBlp, ubest, &NIter);
		printf("\n\n LP:\n");
		printf("   ZLBlp = %lf \n\n", ZLBlp);

		if (Bol_SR3)
		{
		ripeti:
			bol = Generate_SR3(mip.xt);
			if (bol)
			{
				mip.rncols = 0;          // Number of new columns/variables in the new rows
				mip.rnrows = nSR3;       // Number of new rows/constraints
				mip.rnz = rmatbeg[nSR3]; // Number of non-zero coefficients into the new rows
				for (i = 0; i < nSR3; i++)
				{
					mip.rrhs[i] = rrhs[i];
					mip.rsense[i] = rsense[i];
					mip.rmatbeg[i] = rmatbeg[i];    // Pointers to the beginning of each rows of the constraint matrix
					mip.rmatcnt[i] = rmatbeg[i + 1] - rmatbeg[i];    // Cardinality of each rows of the constraint matrix
					for (j = rmatbeg[i]; j < rmatbeg[i + 1]; j++)
					{
						mip.rmatind[j] = rmatind[j];    // Column index associated to each coefficient of the constraint matrix
						mip.rmatval[j] = rmatval[j]; // Value of each coefficient of the constraint matrix
					}
				}
				mip.rmatbeg[nSR3] = rmatbeg[nSR3];

				mip.AddNewRows();
				//mip.SolveLP(&ZLBlp, ubest, &NIter);
				mip.SolveLP(&ZLBlp, NULL, &NIter);
				printf("\n\n LP:\n");
				printf("   ZLBlp = %lf \n\n", ZLBlp);

				goto ripeti;
			}
		}
	}

	mip.FreeCols();

	if (Bol_SR3)
	{
		mip.Free_Valid_Inequalities();
	}

	mip.Close();

	return status;
}
 
//-----------------------------------------------------------------------------
// Generate SR3 
//-----------------------------------------------------------------------------
//
int Subgradient::Generate_SR3(double *xt)
{
	long i1, i2, i3;
	long j, k;
	long kk;
	long konta;
	double xsum;
	int bol = 0;

	nSR3 = 0;
	kk = 0;
	rmatbeg[0] = 0;
	for (i1 = 0; i1 < Inst.n - 2; i1++)
		for (i2 = i1 + 1; i2 < Inst.n - 1; i2++)
			for (i3 = i2 + 1; i3 < Inst.n; i3++)
			{
				xsum = 0;
				for (j = 0; j < n; j++)
				{
					if ((xt[j] > Prec) && (xt[j] < 1. - Prec))
					{
						konta = 0;
						for (k = matbeg[j]; k < matbeg[j + 1]; k++)
							if ((matind[k] == i1) || (matind[k] == i2) || (matind[k] == i3))
								konta++;
						if (konta >= 2)
						{
							xsum += xt[j];
						}
					}
				}
				if (xsum > 1.0 + Prec)
				{
					bol = 1;
					printf(" SR3: (%d,%d,%d) LHS=%lf \n", i1, i2, i3, xsum);
					kk = rmatbeg[nSR3];
					for (j = 0; j < n; j++)
					{
						konta = 0;
						for (k = matbeg[j]; k < matbeg[j + 1]; k++)
							if ((matind[k] == i1) || (matind[k] == i2) || (matind[k] == i3))
								konta++;
						if (konta >= 2)
						{
							rmatind[kk] = j;
							rmatval[kk] = 1.0;
							kk++;
						}
					}
					rsense[nSR3] = 'L';
					rrhs[nSR3] = 1.0;
					nSR3++;
					rmatbeg[nSR3] = kk;
				}
			}

	return bol;

}

//-----------------------------------------------------------------------------
// Subgradient for the Set Covering with cardinality constraint using a 
// the Parametric Lagrangian Relaxation
//-----------------------------------------------------------------------------
//
int Subgradient::Dual_Ascent_SetCov(int FlagHeu)
{
	int status = 0;
	int bol;
	long i, j, k, ii, jj, i1, kk;
	long maxit, it;
	long maxiti, iti;
	long maxitf, itf;
	long maxitcg, itcg;
	double den, alfa, beta, delta, gamma1, gamma2, gamma3;
	double MinGap;
	double cpen;
	double *cpmin;
	double fbest;
	double t1, t2, dt;
	double objf;
	double Qj;

	double zheu, bestzheu;
	long *indheu;
	double *objlgr;
	double *objheu;
	double *subheu;
	double *xheu, *bestxheu;

	double *preddir;
	double rho;

	long maxKPk;

	// Get the starting time
	t1 = ((double)clock()) / CLOCKS_PER_SEC;

	// Allocate auxiliary data structures
	xheu = new double[nmax];
	bestxheu = new double[nmax];
	objlgr = new double[nmax];
	objheu = new double[nmax];
	indheu = new long[nmax];
	cpmin = new double[m + mc];
	subheu = new double[m + mc];

	// Allocate the array used to save the direction computed in the previous iteration
	preddir = new double[m + mc];
	rho = 0.5;

	// Best heuristic solution
	zheu = 10000.0;
	bestzheu = 10000.0;

	// Parameters
	maxit = 100000; // Maximum number of iterations
	maxiti = 10;    // Maximum number of iteration without improvement before to reduce the step
	maxitf = Param.MaxItMinGap;   // Maximum number of iteration without improvement larger than MinGap before to stop
	maxitcg = 1;    // Maximum number of iteration for the column generation
	alfa = Param.alpha;    // Starting step
	delta = 0.01;   // 0.01
	gamma1 = 0.50;  // Reduction factor for the step
	gamma2 = 1.20;  // Increment factor for the step
	gamma3 = 1.00;  // Factor for generating new columns (fbest<gamma3*objval)
	MinGap = Param.MinGap;  // The value of MinGap
	maxKPk = (long)(Param.PercK * Inst.n);  // Maximum number of columns generated at each iteration
	if (maxKPk < 1)
		maxKPk = 1;

	// Setup initial values
	fbest = -Inf;
	objf = -Inf;

	// Set flag (cflag[j]=1 if j is generated) and index (indl[j]=jj if jj is th j-th column) for each column
	for (j = 0; j < n; j++)
	{
		cflag[j] = 0;
		indl[j] = j;
		objpen[j] = obj[j];
		objlgr[j] = obj[j];
		bestxheu[j] = 0.0;
	}

	// Sort columns for non decreasing cost
	QSortD(obj, indl, 0, n - 1, -1);

	// Setup weight q, penalty pen and flag associated to each row 
	for (i = 0; i < m + mc; i++)  // Note: m+mc
	{
		q[i] = 1.0;
		pen[i] = 1.0;
		rflag[i] = 0;
		preddir[i] = 0.0;
	}

	// Setup the first core of columns for starting the Lagrangian Column Generation
	for (jj = 0; jj < n; jj++)
	{
		bol = 0;
		j = indl[jj];
		for (k = matbeg[j]; k < matbeg[j + 1]; k++)
		{
			i = matind[k];
			if (rflag[i] == 0)
			{
				bol = 1;
				break;
			}
		}
		if (bol)
		{
			cflag[j] = 1;
			for (k = matbeg[j]; k < matbeg[j + 1]; k++)
			{
				i = matind[k];
				rflag[i]++;
			}
		}
		else
			cflag[j] = 0;
	}

	// The subgradient main loop
	iti = 0;
	itf = 0;
	itcg = maxitcg;
	for (it = 0; it < maxit; it++)
	{
		for (j = 0; j < n; j++)
		{
			x[j] = 0.0;
		}

		// For each row select the columns as follows 
		// - if c'[j]>=0 for every j \in Core: select column j' such that c'[j']/Q[j']=min{c'[j]/Q[j] : j \in Core}
		// - if at least c'[j]<0: select every column of negative reduced cost 
		objval = 0.0;
		for (i = 0; i < m; i++)
		{
			jmin[i] = -1;
			cpmin[i] = +Inf;
			Qmin[i] = 0.0;
		}
		for (j = 0; j < n; j++)
		{
			if (cflag[j] == 0)  // Is column j into the core?
				continue;
			if (matbeg[j + 1] - matbeg[j] < 1)  // Is column j empty?
				continue;

			// Evaluate dual variable u[i]
			Qj = 0.0;
			cpen = 0.0;
			for (k = matbeg[j]; k < matbeg[j + 1]; k++)
			{
				ii = matind[k];
				Qj = Qj + matval[k] * pen[ii];
			}

			objlgr[j] = 0.0;
			for (k = matbeg[j]; k < matbeg[j + 1]; k++)
			{
				i = matind[k];

				cpen = obj[j];
				if ((Qj > -0.0001) && (Qj < 0.0001))
					cpen = cpen / (double)(matbeg[j + 1] - matbeg[j]);
				else
					cpen = (pen[i] * cpen) / Qj;

				if (cpen < cpmin[i])
				{
					jmin[i] = j;
					cpmin[i] = cpen;
					Qmin[i] = Qj;
				}
			}
		}

		for (i = 0; i < m; i++)
		{
			u[i] = cpmin[i];
			objval += u[i] * rhs[i];
		}

		// Compute the subgradient vector
		for (i = 0; i < m; i++)
			sub[i] = rhs[i];
		//for (i1=0;i1<mc;i1++)
		//	sub[m+i1]=cardrhs[i1];

		// Define the current solution
		for (i = 0; i < m; i++)
		{
			j = jmin[i];

			if ((Qmin[i] > -0.001) && (Qmin[i] < 0.001))
			{
				x[j] += 1.0 / (double)(matbeg[j + 1] - matbeg[j]);
			}
			else
			{
				x[j] += pen[i] / Qmin[i];
			}
		}

		// Evaluate the subgradient
		for (j = 0; j < n; j++)
		{
			if (cflag[j] == 0)  // Is column j into the core?
				continue;

			if (x[j] < Prec)  // Is column j into the Lagrangian solution?
				continue;

			for (k = matbeg[j]; k < matbeg[j + 1]; k++)
			{
				ii = matind[k];
				sub[ii] = sub[ii] - matval[k] * x[j];
			}

		}

		// Compute a heuristic solution (a very simple heuristic... not used for BPP)
		//zheu = Inf;
		zheu = 10000.0;
		if (FlagHeu > 0)
		{
			/*
			... add Heuristic Algorithms to obtain a "good" feasible solution
			*/
		}

		// Print log
		if ((it % 100 == 0) && (Iprinx))
		{
			printf(" %d) alfa = %8.5lf   objval = %12.2lf   fbest = %12.2lf   ncol = %7d   fbestheu = %12.2lf \n", it, alfa, objval, fbest, n, bestzheu);
		}

		// Is the lower bound improved?
		if (fbest < objval)
		{
			itcg++;
			if (itcg < maxitcg)
				continue;
			else
			{
				itcg = 0;
				bol = 1;

				// Setup KP to generate columns
				for (ii = 0; ii < Inst.n; ii++)
					KPSol.Inst.pf[ii] = (float)u[ii];
				KPSol.Inst.k = maxKPk; // 10;
				KPSol.Solve_KPBW_kBest_LDCF_CG(0, (float)1.0001);
				if (KPSol.kBSol.nk > 0)
				{
					// Check nmax
					if (n >= nmax - 1)
						return 1;

					// Add columns
					jj = matbeg[n];
					for (kk = 0; kk < KPSol.kBSol.nk; kk++)
					{
						if (KPSol.kBSol.zf[kk] < 1.0001)
							break;
						bol = 0;
						objpen[n] = 1.0;
						objlgr[n] = 1.0;
						bestxheu[n] = 0.0;
						cflag[n] = 1;
						indl[n] = n;
						obj[n] = 1.0;
						lb[n] = 0.0;
						ub[n] = 1.0;
						for (ii = KPSol.kBSol.nx[kk] - 1; ii >= 0; ii--)
						{
							matind[jj] = KPSol.kBSol.x[kk][ii];
							matval[jj] = 1.0;
							jj++;
						}
						n++;
						matbeg[n] = jj;
					}
				}

				// Check if the dual solution is feasible (bol==0 -> not feasible)
				if (bol == 0)
					continue;
				else
				{
					alfa *= gamma2; // 1.02
					iti = 0;

					if (fbest < objval)
					{
						fbest = objval;
						for (i = 0; i < m; i++)
							ubest[i] = u[i];
					}

					// Eliminate columns with high penalized cost
					if (it % 1 == 0)
					{
						for (j = 0; j < n; j++)
						{
							if (cflag[j] == 0)
								continue;

							cpen = 0.0;
							for (k = matbeg[j]; k < matbeg[j + 1]; k++)
							{
								ii = matind[k];
								cpen = cpen - matval[k] * u[ii];
							}
							cpen = cpen + obj[j];
							if (j < n1 + n2) // For artificial variables the carinality constraint penalty is not added!
								for (i1 = 0; i1 < mc; i1++)
								{
									k = j * mc + i1;
									if (cardmat[k] > 0)
										cpen = cpen + pen[m + i1]; // Also the carinality constraint penalty is added!
								}
						}
					}
				}
			}
		}
		else
		{
			if (iti > maxiti)
			{
				alfa *= gamma1; 
				iti = 0;
			}
			else
				iti++;
		}

		// Stop criteria
		if (((fbest - objf) / fbest) < MinGap)
		{
			if (itf > maxitf)
				goto TERMINATE;
			else
				itf++;
		}
		else
		{
			itf = 0;
			objf = fbest;
		}

		// Update Penalties
		den = 0.0;
		bol = 1;
		for (i = 0; i < m + mc; i++)
		{
			if (((pen[i] < -Prec) || (pen[i] > +Prec)) && ((sub[i] < -Prec) || (sub[i] > +Prec)))
				bol = 0;
			den += sub[i] * sub[i];
		}

		// Check optimality and feasibility
		if (den < Prec)
		{
			continue;
		}
		beta = alfa * (delta*objval) / den;

		if (beta < -Prec)
			beta = -beta;
		if (beta < Prec)
			beta = Prec;

		for (i = 0; i < m; i++)
		{
			pen[i] = pen[i] + rho * beta * sub[i] + (1.0 - rho)*preddir[i];  
			preddir[i] = beta * sub[i];  

			if (pen[i] < Prec)
				pen[i] = 0.0;
		}
	}

TERMINATE:

	// Free auxiliary data structures
	delete[] bestxheu;
	delete[] indheu;
	delete[] objlgr;
	delete[] objheu;
	delete[] subheu;
	delete[] cpmin;
	delete[] xheu;

	t2 = ((double)clock()) / CLOCKS_PER_SEC;
	dt = t2 - t1;

	if (Iprinx)
	{
		printf(" Best objval = %lf \n", fbest);
		printf(" Time = %lf sec.\n", dt);
	}

	// Setup output
	ZLBlagr = fbest;
	ZUBlagr = bestzheu;
	alpha = alfa;
	Iter = it;
	nCore = 0;
	for (j = 0; j < n; j++)
		if (cflag[j] > 0)
			nCore++;

	return (status);
}

//-----------------------------------------------------------------------------
//  Utility for reading a list of integer values of unknown lenght                                    
//-----------------------------------------------------------------------------
//
long Subgradient::GetList(FILE *fin, long *nlist, int *list)
{
	static char line[255];
	char cc;
	long k, h;
	int bol;

	h = 0;
loop0:

	k = 0;
	bol = 0;
loop1:
	cc = fgetc(fin);
	if ((cc >= '0') && (cc <= '9'))
	{
		bol = 1;
		line[k] = cc;
		k++;
		goto loop1;
	}

	if (bol < 1)
	{
		if (cc != '\n')
			goto loop0;
	}
	else
	{
		line[k] = '\0';
		list[h] = atoi(line);
		h++;
		if (cc != '\n')
			goto loop0;
	}

	*nlist = h;

	return k;
}

//-----------------------------------------------------------------------------
//  Random                                    
//-----------------------------------------------------------------------------
//
double Subgradient::random(void)
{
	return (double)rand()/RAND_MAX;
}

//-----------------------------------------------------------------------------
//  Compute a lower bound for the BPP by Dual Feasible Functions                                    
//-----------------------------------------------------------------------------
//
long Subgradient::LBDFF(void)
{
	long i;
	long alpha;
	long LBBest;
	long LB1;
	double LB;
	double aux;
	long auxl;

	LBBest = 0;

	// DFF 1
	for (alpha = 1; alpha <= Inst.W; alpha++)
	{
		LB = 0;
		for (i = 0; i < Inst.n; i++)
		{
			aux = (double)(alpha + 1)*(double)((double)Inst.w[i] / Inst.W);
			auxl = (long)(aux + 0.0001);

			if (fabs((double)aux - (double)auxl) < 0.0001)
				LB += Inst.w[i];
			else
				LB += (long)aux*(Inst.W / alpha);
		}

		LB1 = long(((double)LB / (double)Inst.W) + 0.999);
		if (LB1 > LBBest)
			LBBest = LB1;
	}

	// DFF 1-b
	for (alpha = 1; alpha <= Inst.W; alpha++)
	{
		LB = 0;
		for (i = 0; i < Inst.n; i++)
		{
			aux = (double)(alpha + 1)*(double)((double)Inst.w[i] / Inst.W);
			auxl = (long)(aux + 0.0001);

			if (fabs((double)aux - (double)auxl) < 0.0001)
				LB += Inst.w[i] * alpha;
			else
				LB += (long)aux*(Inst.W);
		}

		LB1 = long(((double)LB / (double)((double)(alpha + 1)*Inst.W) + 0.999));
		if (LB1 > LBBest)
			LBBest = LB1;
	}

	// DFF 2
	for (alpha = 0; alpha <= Inst.W / 2; alpha++)
	{
		LB = 0;
		for (i = 0; i < Inst.n; i++)
		{
			if (Inst.w[i] > Inst.W - alpha)
				LB += Inst.W;
			else if (Inst.w[i] >= alpha)
				LB += Inst.w[i];
		}

		LB1 = long(((double)LB / (double)Inst.W) + 0.999);
		if (LB1 > LBBest)
			LBBest = LB1;
	}

	// DFF 3
	for (alpha = 1; alpha <= Inst.W / 2; alpha++)
	{
		LB = 0;
		for (i = 0; i < Inst.n; i++)
		{
			if (Inst.w[i] > Inst.W / 2)
				LB += 2 * (Inst.W / alpha - (Inst.W - Inst.w[i]) / alpha);
			else if (Inst.w[i] == Inst.W / 2)
				LB += Inst.W / alpha;
			else
				LB += 2 * (Inst.w[i] / alpha);
		}

		LB1 = long(((double)LB / (double)(2 * ((double)Inst.W / alpha))) + 0.999);
		if (LB1 > LBBest)
			LBBest = LB1;
	}

	printf("LB-DFF = %ld \n", LBBest);

	ZLBDFF = LBBest;

	return LBBest;
}

//-----------------------------------------------------------------------------
//  Compute a lower bound for the BPP by Dual Feasible Functions                                    
//-----------------------------------------------------------------------------
//
long Subgradient::LBDFFC(void)
{
	long i, j, k;
	long ai, aj;
	long maxai, maxaj;
	long LBBest;
	long LB;
	long auxw, auxW;

	LBBest = 0;

	for (i = 0; i <= 2; i++)
	{
		if (i == 0)
			maxai = Inst.W;
		else
			maxai = Inst.W / 2;
		for (ai = 1; ai <= maxai; ai++)
		{
			for (j = 0; j <= 2; j++)
			{
				if (j == 0)
					maxaj = Inst.W;
				else
					maxaj = Inst.W / 2;
				for (aj = 1; aj <= maxaj; aj++)
				{
					LB = 0;
					for (k = 0; k < Inst.n; k++)
					{
						auxw = DFF(j, aj, Inst.w[k], Inst.W);
						auxW = DFF(j, aj, Inst.W, Inst.W);
						LB += DFF(i, ai, auxw, auxW);
					}

					auxw = DFF(j, aj, Inst.W, Inst.W);
					auxW = DFF(j, aj, Inst.W, Inst.W);
					LB = long(((double)LB / (double)DFF(i, ai, auxw, auxW)) + 0.999);
					if (LB > LBBest)
						LBBest = LB;
				}
			}
		}
	}

	printf("LB-DFFC = %ld \n", LBBest);

	return LBBest;
}

//-----------------------------------------------------------------------------
//  Compute a lower bound for the BPP by Dual Feasible Functions                                    
//-----------------------------------------------------------------------------
//
long Subgradient::LBDFFC1(void)
{
	long i, j, k;
	long ai, aj;
	long maxai, maxaj;
	long LBBest;
	long LB;
	long auxw, auxW;

	LBBest = 0;

	for (i = 0; i <= 2; i++)
	{
		if (i == 0)
			maxai = Inst.W;
		else
			maxai = Inst.W / 2;
		for (ai = 1; ai <= maxai; ai++)
		{
			for (j = 0; j <= 2; j++)
			{
				if (j == 0)
					maxaj = DFF(i, ai, Inst.W, Inst.W);
				else
					maxaj = DFF(i, ai, Inst.W, Inst.W) / 2;
				for (aj = 1; aj <= maxaj; aj++)
				{
					LB = 0;
					for (k = 0; k < Inst.n; k++)
					{
						auxw = DFF(i, ai, Inst.w[k], Inst.W);
						auxW = DFF(i, ai, Inst.W, Inst.W);
						LB += DFF(j, aj, auxw, auxW);
					}

					auxw = DFF(i, ai, Inst.W, Inst.W);
					auxW = DFF(i, ai, Inst.W, Inst.W);
					LB = long(((double)LB / (double)DFF(j, aj, auxw, auxW)) + 0.999);
					if (LB > LBBest)
						LBBest = LB;
				}
			}
		}
	}

	printf("LB-DFFC = %ld \n", LBBest);

	return LBBest;
}

//-----------------------------------------------------------------------------
//  Compute a lower bound for the BPP by Dual Feasible Functions                                    
//-----------------------------------------------------------------------------
//
long Subgradient::DFF(long i, long alpha, long w, long W)
{
	double aux;
	long auxl;

	if (i == 0)
	{
		// DFF 1
		aux = (double)(alpha + 1)*(double)((double)w / W);
		auxl = (long)(aux + 0.0001);

		if (fabs((double)aux - (double)auxl) < 0.0001)
			return w;
		else
			return (long)aux*(W / alpha);
	}
	else if (i == 1)
	{
		// DFF 2
		if (w > W - alpha)
			return W;
		else if (w >= alpha)
			return w;
		else
			return 0;
	}
	else if (i == 2)
	{
		if (w > W / 2)
			return 2 * (W / alpha - (W - w) / alpha);
		else if (w == W / 2)
			return W / alpha;
		else
			return 2 * (w / alpha);
	}

	return 0;
}

//-----------------------------------------------------------------------------
//  Reductions      
//-----------------------------------------------------------------------------
//
long Subgradient::Reductions(void)
{
	long i;
	long bol = 0;

	// Knapsack Problem Object
	KPSolver KPSolT;

	KPSolT.Inst.n = Inst.n;
	KPSolT.Inst.W = Inst.W;
	KPSolT.Inst.w = new long[KPSolT.Inst.n];
	KPSolT.Inst.p = new long[KPSolT.Inst.n];

	for (i = 0; i < KPSolT.Inst.n; i++)
	{
		KPSolT.Inst.w[i] = Inst.w[i];
		KPSolT.Inst.p[i] = Inst.w[i];
	}

	KPSolT.Allocate_KPFW();

doagain:

	for (i = 0; i < Inst.n; i++)
	{
		if ((bol == 0) && (Inst.w[i] <= Inst.W / 2))
			continue;

		KPSolT.Inst.W = Inst.W - Inst.w[i];
		KPSolT.Inst.w[i] = 0;
		KPSolT.Inst.p[i] = 0;
		KPSolT.Solve_KPFW();

		if (KPSolT.OptSol.z < KPSolT.Inst.W)
			Inst.w[i] = Inst.w[i] + (KPSolT.Inst.W - KPSolT.OptSol.z);
		KPSolT.Inst.w[i] = Inst.w[i];
		KPSolT.Inst.p[i] = Inst.w[i];
	}
	if (bol == 0)
	{
		bol = 1;
		goto doagain;
	}

	KPSolT.DeAllocate_KPFW();

	return 0;
}
