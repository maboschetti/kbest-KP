//-----------------------------------------------------------------------------
// Here we must insert our "license boilerplate"
// (... to be defined)
//-----------------------------------------------------------------------------
// File: kBest-KP-Library.cpp
// Version 1.1
// Last update: 20.04.2025
// Authors: Boschetti M.A., Novellani, S.
//-----------------------------------------------------------------------------
//
#include <stdlib.h>
#include "Utilities.h"
#include "kBest-KP-Library.h"

/**
 *	KPSolver: Constructor
 */
KPSolver::KPSolver(void)
{
	Inst.n = 0;
	Inst.W = 0;
	Inst.k = 0;

	Inst.w = NULL;
	Inst.p = NULL;
	Inst.pf = NULL;
	Inst.b = NULL;

	LPSol.x = NULL;
	OptSol.x = NULL;
	kBSol.x = NULL;
	kBSol.z = NULL;
	kBSol.nx = NULL;

	Zkp1 = NULL;
	Pred1 = NULL;
	Kard = NULL;
	Kard1 = NULL;
	Kard2 = NULL;

	Zkp = NULL;
	Pred = NULL;

	ZkpL = NULL;
	PredL = NULL;
	KardL = NULL;
	KardL1 = NULL;
	KardL2 = NULL;
	KardB = NULL;

	ZkpLD = NULL;
	PredLD = NULL;
	PredALD = NULL;
	Dim = NULL;
	Bound = NULL;

	wmin = 1;
}

/**
 *  KPSolver: Destructor
 */
KPSolver::~KPSolver(void)
{
	long j, w, k;

	if (Inst.w != NULL)
		delete[] Inst.w;

	if (Inst.p != NULL)
		delete[] Inst.p;

	if (Inst.pf != NULL)
		delete[] Inst.pf;

	if (Inst.b != NULL)
		delete[] Inst.b;

	if (LPSol.x != NULL)
		delete[] LPSol.x;

	if (OptSol.x != NULL)
		delete[] OptSol.x;

	if (Zkp1 != NULL)
	{
		for (j = 0; j <= Inst.n; j++)
			if (Zkp1[j] != NULL)
				delete[] Zkp1[j];
		delete[] Zkp1;
	}

	if (Pred1 != NULL)
	{
		for (j = 0; j <= Inst.n; j++)
			if (Pred1[j] != NULL)
				delete[] Pred1[j];
		delete[] Pred1;
	}

	if (Kard != NULL)
	{
		for (j = 0; j <= Inst.n; j++)
			if (Kard[j] != NULL)
				delete[] Kard[j];
		delete[] Kard;
	}

	if (Kard1 != NULL)
	{
		for (j = 0; j <= Inst.n; j++)
			if (Kard1[j] != NULL)
				delete[] Kard1[j];
		delete[] Kard1;
	}

	if (Kard2 != NULL)
	{
		for (j = 0; j <= Inst.n; j++)
			if (Kard2[j] != NULL)
				delete[] Kard2[j];
		delete[] Kard2;
	}

	if (Zkp != NULL)
	{
		for (j = 0; j <= Inst.n; j++)
		{
			if (Zkp[j] != NULL)
			{
				for (w = 0; w <= Inst.W; w++)
					if (Zkp[j][w] != NULL)
						delete[] Zkp[j][w];
				delete[] Zkp[j];
			}
		}
		delete[] Zkp;
	}

	if (Pred != NULL)
	{
		for (j = 0; j <= Inst.n; j++)
		{
			if (Pred[j] != NULL)
			{
				for (w = 0; w <= Inst.W; w++)
					if (Pred[j][w] != NULL)
						delete[] Pred[j][w];
				delete[] Pred[j];
			}
		}
		delete[] Pred;
	}

	if (ZkpL != NULL)
		delete[] ZkpL;

	if (PredL != NULL)
		delete[] PredL;

	if (KardL != NULL)
		delete[] KardL;

	if (KardL1 != NULL)
		delete[] KardL1;

	if (KardL2 != NULL)
		delete[] KardL2;

	if (ZkpLD != NULL)
	{
		for (j = 0; j < (Inst.n + 1)*(Inst.W + 1); j++)
			if (ZkpLD[j] != NULL)
				delete[] ZkpLD[j];
		delete[] ZkpLD;
	}

	if (PredLD != NULL)
	{
		for (j = 0; j < (Inst.n + 1)*(Inst.W + 1); j++)
			if (PredLD[j] != NULL)
				delete[] PredLD[j];
		delete[] PredLD;
	}

	if (PredALD != NULL)
	{
		for (j = 0; j < (Inst.n + 1)*(Inst.W + 1); j++)
			if (PredALD[j] != NULL)
				delete[] PredALD[j];
		delete[] PredALD;
	}

	if (KardB != NULL)
	{
		for (j = 0; j < (Inst.n + 1)*(Inst.W + 1); j++)
			if (KardB[j] != NULL)
				delete[] KardB[j];
		delete[] KardB;
	}

	if (Dim != NULL)
		delete[] Dim;

	if (Bound != NULL)
		delete[] Bound;

	if (kBSol.x != NULL)
	{
		for (k = 0; k < Inst.k; k++)
			if (kBSol.x[k] != NULL)
				delete[] kBSol.x[k];
		delete[] kBSol.x;
	}

	if (kBSol.z != NULL)
		delete[] kBSol.z;

	if (kBSol.nx != NULL)
		delete[] kBSol.nx;

}

//-----------------------------------------------------------------------------
// Functions to allocate and deallocate data structures to save solutions
//-----------------------------------------------------------------------------

/**
 *  KPSolver: Allocate_Sol
 */
int KPSolver::Allocate_Sol(void)
{
	long k;

	LPSol.x = new double[Inst.n];
	if (LPSol.x == NULL)
		return 1;

	OptSol.x = new long[Inst.n];
	if (OptSol.x == NULL)
		return 1;

	kBSol.nk = 0;
	kBSol.z = new long[Inst.k];
	if (kBSol.z == NULL)
		return 1;

	kBSol.x = new long*[Inst.k];
	if (kBSol.x == NULL)
		return 1;

	for (k = 0; k < Inst.k; k++)
	{
		kBSol.x[k] = new long[Inst.n];
		if (kBSol.x[k] == NULL)
			return 1;
	}

	return 0;
};

/**
 *  KPSolver: DeAllocate_Sol
 */
int KPSolver::DeAllocate_Sol(void)
{
	long k;

	if (LPSol.x != NULL)
	{
		delete[] LPSol.x;
		LPSol.x = NULL;
	}

	if (OptSol.x != NULL)
	{
		delete[] OptSol.x;
		OptSol.x = NULL;
	}

	if (kBSol.z != NULL)
	{
		delete[] kBSol.z;
		kBSol.z = NULL;
	}

	if (kBSol.x != NULL)
	{
		for (k = 0; k < Inst.k; k++)
			if (kBSol.x[k] != NULL)
				delete[] kBSol.x[k];
		delete[] kBSol.x;
		kBSol.x = NULL;
	}

	return 0;
};

//-----------------------------------------------------------------------------
// Functions to solve the LP-relaxation
//-----------------------------------------------------------------------------

/**
 *  KPSolver: Allocate_LPKP
 */
int KPSolver::Allocate_LPKP(void)
{
	ratio = new double[Inst.n];
	if (ratio == NULL)
		return 1;

	ind = new long[Inst.n];
	if (ind == NULL)
		return 1;

	LPSol.x = new double[Inst.n];
	if (LPSol.x == NULL)
		return 1;

	return 0;
}

/**
 *  KPSolver: DeAllocate_LPKP
 */
int KPSolver::DeAllocate_LPKP(void)
{
	if (ratio != NULL)
	{
		delete[] ratio;
		ratio = NULL;
	}

	if (ind != NULL)
	{
		delete[] ind;
		ind = NULL;
	}

	if (LPSol.x != NULL)
	{
		delete[] LPSol.x;
		LPSol.x = NULL;
	}

	return 0;
}

/**
 *  KPSolver: Solve_LPKP
 */
int KPSolver::Solve_LPKP(void)
{
	long i, ii;
	long W1;
	double zub;

	// Compute the ratio values and setup an empty solution
	for (i = 0; i < Inst.n; i++)
	{
		ind[i] = i;
		ratio[i] = (double)Inst.p[i] / (double)Inst.w[i];
		LPSol.x[i] = 0;
	}

	// Sort the items according to the ratio values
	SortValues(Inst.n, ratio, ind, +1);

	// Build the LP-relaxation optimal solution
	zub = 0.0;
	W1 = Inst.W;
	for (ii = 0; ii < Inst.n; ii++)
	{
		i = ind[ii];
		if (Inst.w[i] <= W1)
		{
			LPSol.x[i] = 1.0;
			zub += Inst.p[i];
			W1 -= Inst.w[i];
		}
		else
		{
			LPSol.x[i] = (double)W1 / (double)Inst.w[i];
			zub += Inst.p[i] * (double)W1 / (double)Inst.w[i];
			break;
		}
	}
	LPSol.z = zub;

	if (PrintOn)
		printf("LP-relaxation: %f \n\n", LPSol.z);

	return 0;
}

//-----------------------------------------------------------------------------
// Functions to generate optimal solution with the forward recursion
//-----------------------------------------------------------------------------

/**
 *  KPSolver: Allocate_KPFW
 */
int KPSolver::Allocate_KPFW(void)
{
	long i;

	Zkp1 = new long*[Inst.n + 1];
	if (Zkp1 == NULL) 
		return 1;

	Pred1 = new long*[Inst.n + 1];
	if (Pred1 == NULL) 
		return 1;

	for (i = 0; i <= Inst.n; i++)
	{
		Zkp1[i] = new long[Inst.W + 1];
		if (Zkp1[i] == NULL) 
			return 1;

		Pred1[i] = new long[Inst.W + 1];
		if (Pred1[i] == NULL) 
			return 1;
	}

	OptSol.x = new long[Inst.n];
	if (OptSol.x == NULL) 
		return 1;

	return 0;
}

/**
 *  KPSolver: DeAllocate_KPFW
 */
int KPSolver::DeAllocate_KPFW(void)
{
	long i;

	if (Zkp1 != NULL)
	{
		for (i = 0; i <= Inst.n; i++)
		{
			if (Zkp1[i] != NULL)
				delete[] Zkp1[i];
		}
		delete[] Zkp1;
		Zkp1 = NULL;
	}

	if (Pred1 != NULL)
	{
		for (i = 0; i <= Inst.n; i++)
		{
			if (Pred1[i] != NULL)
				delete[] Pred1[i];
		}
		delete[] Pred1;
		Pred1 = NULL;
	}

	if (OptSol.x != NULL)
	{
		delete[] OptSol.x;
		OptSol.x = NULL;
	}

	return 0;
}

/**
 *  KPSolver: Solve_KPFW
 */
int KPSolver::Solve_KPFW(void)
{
	long i;
	long w;
	long aux;

	// Basic stage (i.e., i=0)
	for (w = 0; w <= Inst.W; w++)
	{
		Zkp1[0][w] = 0;
		Pred1[0][w] = 0;
	}

	// Compute the stages from i=1 to n
	for (i = 0; i < Inst.n; i++)
	{
		// Compute stage i+1
		for (w = 0; w <= Inst.W; w++)
		{
			// Compute state w of stage i+1
			if (w < Inst.w[i])
			{
				Zkp1[i + 1][w] = Zkp1[i][w];
				Pred1[i + 1][w] = 0;
			}
			else
			{
				aux = Zkp1[i][w - Inst.w[i]] + Inst.p[i];
				if (aux > Zkp1[i][w])
				{
					Zkp1[i + 1][w] = aux;
					Pred1[i + 1][w] = 1;
				}
				else
				{
					Zkp1[i + 1][w] = Zkp1[i][w];
					Pred1[i + 1][w] = 0;
				}
			}
		}
	}

	// Build the optimal solution
	OptSol.z = Zkp1[Inst.n][Inst.W];
	w = Inst.W;
	for (i = Inst.n; i > 0; i--)
	{
		if (Pred1[i][w] > 0)
		{
			// The item i is in the solution
			OptSol.x[i - 1] = 1;
			w = w - Inst.w[i - 1];
		}
		else
		{
			// The item i is not in the solution
			OptSol.x[i - 1] = 0;
		}
	}

	// Print the Zkp matrix
	if (PrintOn > 3)
	{
		for (i = 0; i <= Inst.n; i++)
		{
			printf("Stage (%3d): ", i);
			for (w = 0; w <= Inst.W; w++)
			{
				printf("%5d ", Zkp1[i][w]);
			}
			printf("\n");
		}
	}

	// Print the optimal solution
	if (PrintOn > 0)
	{
		printf("ZOpt = %5d \n", OptSol.z);
		if (PrintOn > 2)
		{
			for (i = 0; i < Inst.n; i++)
				printf("x(%d)=%d ", i, OptSol.x[i]);
			printf("\n");
		}
	}

	return 0;
}

//-----------------------------------------------------------------------------
// Functions to generate optimal solution with the backward recursion
//-----------------------------------------------------------------------------

/**
 *  KPSolver: Setup_KPBW
 */
void KPSolver::Setup_KPBW(void)
{
	long j, w;

	// Define the minimum item width (used to stop the search in advance)
	wmin = Inst.W;
	for (j = 0; j < Inst.n; j++)
		if (Inst.w[j] < wmin)
			wmin = Inst.w[j];

	// Setup each state as unexplored
	for (j = 0; j <= Inst.n; j++)
		for (w = 0; w <= Inst.W; w++)
			Kard[j][w] = 0;
}

/**
 *  KPSolver: Setup_KPBW (with recursive function)
 */
void KPSolver::Setup_KPBW(long j, long w)
{
	long w2;

	// Setup the state as unexplored
	Kard[j][w] = 0;

	// Stop criteria
	if (j == 0)
		return;

	// Recursion
	j = j - 1;
	Setup_KPBW(j,w);
	w2 = w - Inst.w[j];
	if (w2 >= 0)
		Setup_KPBW(j, w2);
}

/**
 *  KPSolver: KPBW
 */
long KPSolver::KPBW(long j, long w)
{
	long z1, z2;
	long w2;

	// Stop criteria for the recursion
	//if ((j == 0) || (w == 0)) // The basic version
	if ((j == 0) || (w < wmin)) 
	{
		Zkp1[j][w] = 0;
		Pred1[j][w] = -1;
		Kard[j][w] = 1;
		return 0;
	}

	// Evaluate the alternative (include or not the item j in the optimal solution)
	if (Kard[j - 1][w] == 0)
		z1 = KPBW(j - 1, w);
	else
		z1 = Zkp1[j - 1][w];
	w2 = w - Inst.w[j - 1];
	if (w2 >= 0)
		if (Kard[j - 1][w2] == 0)
			z2 = KPBW(j - 1, w2) + Inst.p[j - 1];
		else
			z2 = Zkp1[j - 1][w2] + Inst.p[j - 1];
	else
	{
		Zkp1[j][w] = z1;
		Pred1[j][w] = 0;
		Kard[j][w] = 1;
		return z1;
	}
	if (z1 >= z2)
	{
		Zkp1[j][w] = z1;
		Pred1[j][w] = 0;
		Kard[j][w] = 1;
		return z1;
	}
	Zkp1[j][w] = z2;
	Pred1[j][w] = 1;
	Kard[j][w] = 1;
	return z2;
}

/**
 *  KPSolver: Allocate_KPBW
 */
int KPSolver::Allocate_KPBW(void)
{
	long i;

	Zkp1 = new long*[Inst.n + 1];
	if (Zkp1 == NULL)
		return 1;

	Pred1 = new long*[Inst.n + 1];
	if (Pred1 == NULL)
		return 1;

	Kard = new long*[Inst.n + 1];
	if (Kard == NULL)
		return 1;

	for (i = 0; i <= Inst.n; i++)
	{
		Zkp1[i] = new long[Inst.W + 1];
		if (Zkp1[i] == NULL)
			return 1;

		Pred1[i] = new long[Inst.W + 1];
		if (Pred1[i] == NULL)
			return 1;

		Kard[i] = new long[Inst.W + 1];
		if (Kard[i] == NULL)
			return 1;
	}

	OptSol.x = new long[Inst.n];
	if (OptSol.x == NULL)
		return 1;

	return 0;
}

/**
 *  KPSolver: DeAllocate_KPBW
 */
int KPSolver::DeAllocate_KPBW(void)
{
	long i;

	if (Zkp1 != NULL)
	{
		for (i = 0; i <= Inst.n; i++)
		{
			if (Zkp1[i] != NULL)
				delete[] Zkp1[i];
		}
		delete[] Zkp1;
		Zkp1 = NULL;
	}

	if (Pred1 != NULL)
	{
		for (i = 0; i <= Inst.n; i++)
		{
			if (Pred1[i] != NULL)
				delete[] Pred1[i];
		}
		delete[] Pred1;
		Pred1 = NULL;
	}

	if (Kard != NULL)
	{
		for (i = 0; i <= Inst.n; i++)
		{
			if (Kard[i] != NULL)
				delete[] Kard[i];
		}
		delete[] Kard;
		Kard = NULL;
	}

	if (OptSol.x != NULL)
	{
		delete[] OptSol.x;
		OptSol.x = NULL;
	}

	return 0;
}

/**
 *  KPSolver: Solve_KPBW
 */
int KPSolver::Solve_KPBW(void)
{
	long i;
	long w;
	//float t1, t2, dt;

	// Setup data structures
	Setup_KPBW();
	//t1 = GetTime();
	//SetupKPBW(Inst.n, Inst.W);
	//t2 = GetTime();
	//dt = t2 - t1;
	//printf("Setup Computing time: %7.4f \n", dt);

	// Build the optimal solution
	OptSol.z = KPBW(Inst.n, Inst.W);
	w = Inst.W;
	for (i = Inst.n; i > 0; i--)
	{
		if (Pred1[i][w] > 0)
		{
			// The item i is in the solution
			OptSol.x[i - 1] = 1;
			w = w - Inst.w[i - 1];
		}
		else
		{
			// The item i is not in the solution
			OptSol.x[i - 1] = 0;
		}
	}

	// Print the optimal solution
	if (PrintOn > 0)
	{
		printf("ZOpt = %5d \n", OptSol.z);
		if (PrintOn > 2)
		{
			for (i = 0; i < Inst.n; i++)
				printf("x(%d)=%d ", i, OptSol.x[i]);
			printf("\n");
		}
	}

	return 0;
}

//-----------------------------------------------------------------------------
// Functions to generate k-best solutions with forward recursion
//-----------------------------------------------------------------------------

/**
 *  KPSolver: Generate_State
 */
int KPSolver::Generate_State(long j, long w)
{
	long kk, k1, k2, kl;
	long kmax1, kmax2;
	long w2;

	kk = 0;
	k1 = 0;
	k2 = 0;
	kmax1 = Kard[j - 1][w];
	w2 = w - Inst.w[j - 1];
	if (w2 >= 0)
		kmax2 = Kard[j - 1][w2];
	else
		kmax2 = 0;

	while (((k1 < kmax1) || (k2 < kmax2)) && (kk < Inst.k))
	{
		kl = 0;
		Zkp[j][w][kk] = -1000;
		if (k1 < kmax1)
		{
			if (Zkp[j][w][kk] < Zkp[j - 1][w][k1])
			{
				kl = -k1-1;
				Zkp[j][w][kk] = Zkp[j - 1][w][k1];
			}
		}
		if (k2 < kmax2)
		{
			w2 = w - Inst.w[j - 1];
			if (Zkp[j][w][kk] < Zkp[j - 1][w2][k2] + Inst.p[j - 1])
			{
				kl = k2+1;
				Zkp[j][w][kk] = Zkp[j - 1][w2][k2] + Inst.p[j - 1];
			}
		}
		Pred[j][w][kk] = kl;
		if (kl < 0)
			k1++;
		else
			k2++;
		kk++;
	}
	Kard[j][w] = kk;

	return 0;
}

/**
 *  KPSolver: Allocate_KPFW_kBest
 */
int KPSolver::Allocate_KPFW_kBest(void)
{
	long j, w, k;

	Zkp = new long**[Inst.n + 1];
	if (Zkp == NULL)
		return 1;

	Pred = new long**[Inst.n + 1];
	if (Pred == NULL)
		return 1;

	Kard = new long*[Inst.n + 1];
	if (Kard == NULL)
		return 1;

	for (j = 0; j <= Inst.n; j++)
	{
		Zkp[j] = new long*[Inst.W + 1];
		if (Zkp[j] == NULL)
			return 1;

		Pred[j] = new long*[Inst.W + 1];
		if (Pred[j] == NULL)
			return 1;

		Kard[j] = new long[Inst.W + 1];
		if (Kard[j] == NULL)
			return 1;

		for (w = 0; w <= Inst.W; w++)
		{
			Zkp[j][w] = new long[Inst.k];
			if (Zkp[j][w] == NULL)
				return 1;

			Pred[j][w] = new long[Inst.k];
			if (Pred[j][w] == NULL)
				return 1;
		}
	}

	kBSol.nk = 0;
	kBSol.z = new long[Inst.k];
	if (kBSol.z == NULL)
		return 1;

	kBSol.x = new long*[Inst.k];
	if (kBSol.x == NULL)
		return 1;

	for (k = 0; k < Inst.k; k++)
	{
		kBSol.x[k] = new long[Inst.n];
		if (kBSol.x[k] == NULL)
			return 1;
	}

	return 0;
}

/**
 *  KPSolver: DeAllocate_KPFW_kBest
 */
int KPSolver::DeAllocate_KPFW_kBest(void)
{
	long j, w, k;

	if (Zkp != NULL)
	{
		for (j = 0; j <= Inst.n; j++)
		{
			if (Zkp[j] != NULL)
			{
				for (w = 0; w <= Inst.W; w++)
				{
					if (Zkp[j][w] != NULL)
						delete[] Zkp[j][w];
				}
				delete[] Zkp[j];
			}
		}
		delete[] Zkp;
		Zkp = NULL;
	}

	if (Pred != NULL)
	{
		for (j = 0; j <= Inst.n; j++)
		{
			if (Pred[j] != NULL)
			{
				for (w = 0; w <= Inst.W; w++)
				{
					if (Pred[j][w] != NULL)
						delete[] Pred[j][w];
				}
				delete[] Pred[j];
			}
		}
		delete[] Pred;
		Pred = NULL;
	}

	if (Kard != NULL)
	{
		for (j = 0; j <= Inst.n; j++)
		{
			if (Kard[j] != NULL)
				delete[]  Kard[j];
		}
		delete[]  Kard;
		Kard = NULL;
	}

	if (kBSol.z != NULL)
	{
		delete[] kBSol.z;
		kBSol.z = NULL;
	}

	if (kBSol.x != NULL)
	{
		for (k = 0; k < Inst.k; k++)
			if (kBSol.x[k] != NULL)
				delete[] kBSol.x[k];
		delete[] kBSol.x;
		kBSol.x = NULL;
	}

	return 0;
}

/**
 *  KPSolver: Solve_KPFW_kBest
 */
int KPSolver::Solve_KPFW_kBest(void)
{
	long j;
	long w;
	long k, k1, ks;
	//float t1, t2, dt;

	// Initialize "Kard"
	for (j = 0; j <= Inst.n; j++)
		for (w = 0; w <= Inst.W; w++)
			Kard[j][w] = 0;

	// Basic stage (i.e., i=0)
	for (w = 0; w <= Inst.W; w++)
	{
		Zkp[0][w][0] = 0;
		Pred[0][w][0] = -1;
		Kard[0][w] = 1;
	}

	// Compute the stages from i=1 to n
	//t1 = GetTime();
	for (j = 1; j <= Inst.n; j++)
	{
		// Compute stage i+1
		for (w = 0; w <= Inst.W; w++)
		{
			Generate_State(j, w);
		}
	}
	//t2 = GetTime();
	//dt = t2 - t1;
	//printf("Generate_State Computing time: %7.4f \n", dt);

	// Build the k-Best solutions
	//t1 = t2;
	ks = Kard[Inst.n][Inst.W];
	kBSol.nk = ks;
	for (k = 0; k < ks; k++)
	{
		kBSol.z[k] = Zkp[Inst.n][Inst.W][k];
		k1 = k;
		w = Inst.W;
		for (j = Inst.n; j > 0; j--)
		{
			if (Pred[j][w][k1] > 0)
			{
				// The item i is in the solution
				kBSol.x[k][j - 1] = 1;
				k1 = Pred[j][w][k1]-1;
				w = w - Inst.w[j - 1];
			}
			else
			{
				// The item i is not in the solution
				kBSol.x[k][j - 1] = 0;
				k1 = -Pred[j][w][k1]-1;
			}
		}
	}
	//t2 = GetTime();
	//dt = t2 - t1;
	//printf("Build Computing time: %7.4f \n", dt);

	// Print the Zkp matrix
	if (PrintOn > 3)
	{
		for (k = 0; k < ks; k++)
			for (j = 0; j <= Inst.n; j++)
			{
				printf("Stage (%3d): ", j);
				for (w = 0; w <= Inst.W; w++)
				{
					printf("%5d ", Zkp[j][w][k]);
				}
				printf("\n");
			}
	}

	// Print the optimal solution
	if (PrintOn > 0)
	{
		printf("Zkp(%d) = %d \n", 0, kBSol.z[0]);
		printf("Zkp(%d) = %d \n", ks - 1, kBSol.z[ks - 1]);
	}
	if (PrintOn > 1)
	{
		for (k = 0; k < ks; k++)
		{
			printf("Zkp(%d) = %5d \n", k, kBSol.z[k]);
			if (PrintOn > 2)
			{
				for (j = 0; j < Inst.n; j++)
					printf("x(%d,%d)=%d ", k, j, kBSol.x[k][j]);
				printf("\n");
			}
		}
	}

	return 0;
}

//-----------------------------------------------------------------------------
// Functions to generate k-best solutions with backward recursion
//-----------------------------------------------------------------------------

/**
 *  KPSolver: Setup_KPBW_kBest
 */
void KPSolver::Setup_KPBW_kBest(void)
{
	long j, w;

	// Define the minimum item width (esed to stop the search in advance)
	wmin = Inst.W;
	for (j = 0; j < Inst.n; j++)
		if (Inst.w[j] < wmin)
			wmin = Inst.w[j];

	// Setup each state as unexplored
	for (j = 0; j <= Inst.n; j++)
		for (w = 0; w <= Inst.W; w++)
		{
			Pred[j][w][0] = 0;
			Kard[j][w] = 0;
			Kard1[j][w] = 0;  
			Kard2[j][w] = 0;  
		}
}

/**
 *  KPSolver: Setup_KPBW_kBest
 */
void KPSolver::Setup_KPBW_kBest(long j, long w)
{
	long w2;

	// Setup the state as unexplored
	Kard[j][w] = 0;

	// Stop criteria
	if (j == 0)
		return;

	// Recursion
	j = j - 1;
	Setup_KPBW_kBest(j, w);
	w2 = w - Inst.w[j];
	if (w2 >= 0)
		Setup_KPBW_kBest(j, w2);
}

/**
 *  KPSolver: KPBW_kBest
 */
long KPSolver::KPBW_kBest(long j, long w, long k)
{
	long z1, z2;
	long w2;

	// The current state has been already generated
	if (k < Kard[j][w])
		return Zkp[j][w][k];

	// Stop criteria for the recursion
	//if ((j == 0) || (w == 0)) // Basic version
	if ((j == 0) || (w < wmin)) 
	{
		if (k == 0)
		{
			Zkp[j][w][k] = 0;
			Pred[j][w][k] = 0;
			Kard[j][w] = 1;
			return 0;
		}
		else
			return -Infi;
	}

	// Evaluate alternatives (include or not the item j in the optimal solution)
	z1 = KPBW_kBest(j - 1, w, Kard1[j][w]);
	w2 = w - Inst.w[j - 1];
	if (w2 >= 0)
	{
		//z2 = KPBW_kBest(j - 1, w2, Kard2[j][w]) + Inst.p[j - 1];
		z2 = KPBW_kBest(j - 1, w2, Kard2[j][w]);
		if (z2 > -Infi)
			z2 = z2 + Inst.p[j - 1];
	}
	else
		z2 = -Infi;
	if (z1 >= z2)
	{
		if (z1 == -Infi)
			return -Infi;

		Zkp[j][w][k] = z1;
		Pred[j][w][k] = -(Kard1[j][w] + 1);
		Kard1[j][w]++;
	}
	else
	{
		Zkp[j][w][k] = z2;
		Pred[j][w][k] = +(Kard2[j][w] + 1);
		Kard2[j][w]++;
	}
	Kard[j][w] = k + 1;
	return Zkp[j][w][k];
}

/**
 *  KPSolver: Allocate_KPBW_kBest
 */
int KPSolver::Allocate_KPBW_kBest(void)
{
	long j, w, k;

	Zkp = new long**[Inst.n + 1];
	if (Zkp == NULL)
		return 1;

	Pred = new long**[Inst.n + 1];
	if (Pred == NULL)
		return 1;

	Kard = new long*[Inst.n + 1];
	if (Kard == NULL)
		return 1;

	Kard1 = new long*[Inst.n + 1];
	if (Kard1 == NULL)
		return 1;

	Kard2 = new long*[Inst.n + 1];
	if (Kard2 == NULL)
		return 1;

	for (j = 0; j <= Inst.n; j++)
	{
		Zkp[j] = new long*[Inst.W + 1];
		if (Zkp[j] == NULL)
			return 1;

		Pred[j] = new long*[Inst.W + 1];
		if (Pred[j] == NULL)
			return 1;

		Kard[j] = new long[Inst.W + 1];
		if (Kard[j] == NULL)
			return 1;

		Kard1[j] = new long[Inst.W + 1];
		if (Kard1[j] == NULL)
			return 1;

		Kard2[j] = new long[Inst.W + 1];
		if (Kard2[j] == NULL)
			return 1;

		for (w = 0; w <= Inst.W; w++)
		{
			Zkp[j][w] = new long[Inst.k];
			if (Zkp[j][w] == NULL)
				return 1;

			Pred[j][w] = new long[Inst.k];
			if (Pred[j][w] == NULL)
				return 1;
		}
	}

	kBSol.nk = 0;
	kBSol.z = new long[Inst.k];
	if (kBSol.z == NULL)
		return 1;

	kBSol.x = new long*[Inst.k];
	if (kBSol.x == NULL)
		return 1;

	for (k = 0; k < Inst.k; k++)
	{
		kBSol.x[k] = new long[Inst.n];
		if (kBSol.x[k] == NULL)
			return 1;
	}

	return 0;
}

/**
 *  KPSolver: DeAllocate_KPBW_kBest
 */
int KPSolver::DeAllocate_KPBW_kBest(void)
{
	long j, w, k;

	if (Zkp != NULL)
	{
		for (j = 0; j <= Inst.n; j++)
		{
			if (Zkp[j] != NULL)
			{
				for (w = 0; w <= Inst.W; w++)
					if (Zkp[j][w] != NULL)
						delete[] Zkp[j][w];
				delete[] Zkp[j];
			}
		}
		delete[] Zkp;
		Zkp = NULL;
	}

	if (Pred != NULL)
	{
		for (j = 0; j <= Inst.n; j++)
		{
			if (Pred[j] != NULL)
			{
				for (w = 0; w <= Inst.W; w++)
					if (Pred[j][w] != NULL)
						delete[] Pred[j][w];
				delete[] Pred[j];
			}
		}
		delete[] Pred;
		Pred = NULL;
	}

	if (Kard != NULL)
	{
		for (j = 0; j <= Inst.n; j++)
		{
			if (Kard[j] != NULL)
				delete[] Kard[j];
		}
		delete[] Kard;
		Kard = NULL;
	}

	if (Kard1 != NULL)
	{
		for (j = 0; j <= Inst.n; j++)
		{
			if (Kard1[j] != NULL)
				delete[] Kard1[j];
		}
		delete[] Kard1;
		Kard1 = NULL;
	}

	if (Kard2 != NULL)
	{
		for (j = 0; j <= Inst.n; j++)
		{
			if (Kard2[j] != NULL)
				delete[] Kard2[j];
		}
		delete[] Kard2;
		Kard2 = NULL;
	}

	if (kBSol.z != NULL)
	{
		delete[] kBSol.z;
		kBSol.z = NULL;
	}

	if (kBSol.x != NULL)
	{
		for (k = 0; k < Inst.k; k++)
			if (kBSol.x[k] != NULL)
				delete[] kBSol.x[k];
		delete[] kBSol.x;
		kBSol.x = NULL;
	}

	return 0;
}

/**
 *  KPSolver: MAllocate_KPBW_kBest
 */
int KPSolver::MAllocate_KPBW_kBest(void)
{
	long j, w, k;

	Zkp = (long***)malloc((Inst.n + 1) * sizeof(long**));
	if (Zkp == NULL)
		return 1;

	Pred = (long***)malloc((Inst.n + 1) * sizeof(long**));
	if (Pred == NULL)
		return 1;

	Kard = (long**)malloc((Inst.n + 1) * sizeof(long*));
	if (Kard == NULL)
		return 1;

	Kard1 = (long**)malloc((Inst.n + 1) * sizeof(long*));
	if (Kard1 == NULL)
		return 1;

	Kard2 = (long**)malloc((Inst.n + 1) * sizeof(long*));
	if (Kard2 == NULL)
		return 1;

	for (j = 0; j <= Inst.n; j++)
	{
		Zkp[j] = (long**)malloc((Inst.W + 1) * sizeof(long*));
		if (Zkp[j] == NULL)
			return 1;

		Pred[j] = (long**)malloc((Inst.W + 1) * sizeof(long*));
		if (Pred[j] == NULL)
			return 1;

		Kard[j] = (long*)malloc((Inst.W + 1) * sizeof(long));
		if (Kard[j] == NULL)
			return 1;

		Kard1[j] = (long*)malloc((Inst.W + 1) * sizeof(long));
		if (Kard1[j] == NULL)
			return 1;

		Kard2[j] = (long*)malloc((Inst.W + 1) * sizeof(long));
		if (Kard2[j] == NULL)
			return 1;

		for (w = 0; w <= Inst.W; w++)
		{
			Zkp[j][w] = (long*)malloc((Inst.k) * sizeof(long));
			if (Zkp[j][w] == NULL)
				return 1;

			Pred[j][w] = (long*)malloc((Inst.k) * sizeof(long));
			if (Pred[j][w] == NULL)
				return 1;
		}
	}

	kBSol.nk = 0;
	kBSol.z = new long[Inst.k];
	if (kBSol.z == NULL)
		return 1;

	kBSol.x = new long*[Inst.k];
	if (kBSol.x == NULL)
		return 1;

	for (k = 0; k < Inst.k; k++)
	{
		kBSol.x[k] = new long[Inst.n];
		if (kBSol.x[k] == NULL)
			return 1;
	}

	return 0;
}

/**
 *  KPSolver: DeMAllocate_KPBW_kBest
 */
int KPSolver::DeMAllocate_KPBW_kBest(void)
{
	long j, w, k;

	if (Zkp != NULL)
	{
		for (j = 0; j <= Inst.n; j++)
			if (Zkp[j] != NULL)
			{
				for (w = 0; w <= Inst.W; w++)
					if (Zkp[j][w] != NULL)
						free(Zkp[j][w]);
				free(Zkp[j]);
			}
		free(Zkp);
		Zkp = NULL;
	}

	if (Pred != NULL)
	{
		for (j = 0; j <= Inst.n; j++)
			if (Pred[j] != NULL)
			{
				for (w = 0; w <= Inst.W; w++)
					if (Pred[j][w] != NULL)
						free(Pred[j][w]);
				free(Pred[j]);
			}
		free(Pred);
		Pred = NULL;
	}

	if (Kard != NULL)
	{
		for (j = 0; j <= Inst.n; j++)
			if (Kard[j] != NULL)
				free(Kard[j]);
		free(Kard);
		Kard = NULL;
	}

	if (Kard1 != NULL)
	{
		for (j = 0; j <= Inst.n; j++)
			if (Kard1[j] != NULL)
				free(Kard1[j]);
		free(Kard1);
		Kard1 = NULL;
	}

	if (Kard2 != NULL)
	{
		for (j = 0; j <= Inst.n; j++)
			if (Kard2[j] != NULL)
				free(Kard2[j]);
		free(Kard2);
		Kard2 = NULL;
	}

	if (kBSol.z != NULL)
	{
		delete[] kBSol.z;
		kBSol.z = NULL;
	}

	if (kBSol.x != NULL)
	{
		for (k = 0; k < Inst.k; k++)
			if (kBSol.x[k] != NULL)
				delete[] kBSol.x[k];
		delete[] kBSol.x;
		kBSol.x = NULL;
	}

	return 0;
}

/**
 *  KPSolver: Solve_KPBW_kBest
 */
int KPSolver::Solve_KPBW_kBest(void)
{
	long i, j;
	long k, kk;
	long w;
	//float t1, t2, dt;

	// Setup data structures
	Setup_KPBW_kBest();
	//t1 = GetTime();
	//Setup_KPBW_kBest(Inst.n, Inst.W);
	//t2 = GetTime();
	//dt = t2 - t1;
	//printf("Setup Computing time: %7.4f \n", dt);

	// Build the k-best solutions
	//t1 = GetTime();
	k = 0;
	while (k < Inst.k)
	{
		kBSol.z[k] = KPBW_kBest(Inst.n, Inst.W, k);
		if (kBSol.z[k] <= -Infi)
			break;
		w = Inst.W;
		kk = k;
		for (i = Inst.n; i > 0; i--)
		{
			if (Pred[i][w][kk] > 0)
			{
				// The item i is in the solution
				kBSol.x[k][i - 1] = 1;
				kk = Pred[i][w][kk] - 1;
				w = w - Inst.w[i - 1];
			}
			else if (Pred[i][w][kk] < 0)
			{
				// The item i is not in the solution
				kBSol.x[k][i - 1] = 0;
				kk = -Pred[i][w][kk] - 1;
			}
			else
			{
				kBSol.x[k][i - 1] = 0;
				break;
			}
		}
		for (j = i; j > 0; j--)
			kBSol.x[k][j - 1] = 0;

		k++;
	}
	kBSol.nk = k;
	//t2 = GetTime();
	//dt = t2 - t1;
	//printf("Generate Computing time: %7.4f \n", dt);

	// Print the optimal solution
	if (PrintOn > 0)
	{
		printf("Zkp(%d) = %d \n", 0, kBSol.z[0]);
		printf("Zkp(%d) = %d \n", k - 1, kBSol.z[k - 1]);
	}
	if (PrintOn > 1)
	{
		for (kk = 0; kk < k; kk++)
		{
			printf("Zkp(%d) = %5d \n", kk, kBSol.z[kk]);
			if (PrintOn > 2)
			{
				for (i = 0; i < Inst.n; i++)
					printf("x(%d,%d)=%d ", kk, i, kBSol.x[kk][i]);
				printf("\n");
			}
		}
	}

	return 0;
}

//-----------------------------------------------------------------------------
// Functions to map bi- and tri-dimensional arrays on one dimensional arrays 
//-----------------------------------------------------------------------------

/**
 *  KPSolver: MapFW
 */
size_t KPSolver::MapFW(long j, long w)
{
	return j * nw + w;
}

/**
 *  KPSolver: MapFW
 */
size_t KPSolver::MapFW(long j, long w, long k)
{
	return j * nwk + w * nk + k;  // First choice
	//return k * nnw + j * nw + w;
	//return k * nnw + w * nn + j;
}

//-----------------------------------------------------------------------------
// Functions to generate k-best solutions with forward recursion 
// based on one dimensional arrays 
//-----------------------------------------------------------------------------

/**
 *  KPSolver: Generate_State_L
 */
int KPSolver::Generate_State_L(long j, long w)
{
	long kk, k1, k2, kl;
	long kmax1, kmax2;
	long w2;
	long z1, z2;
	size_t index;

	kk = 0;
	k1 = 0;
	k2 = 0;
	kmax1 = KardL[MapFW(j - 1, w)];
	w2 = w - Inst.w[j - 1];
	if (w2 >= 0)
		kmax2 = KardL[MapFW(j - 1, w2)];
	else
		kmax2 = 0;

	while (((k1 < kmax1) || (k2 < kmax2)) && (kk < Inst.k))
	{
		index = MapFW(j, w, kk);
		kl = 0;
		ZkpL[index] = -1000;
		if (k1 < kmax1)
		{
			z1 = ZkpL[MapFW(j - 1, w, k1)];
			if (ZkpL[index] < z1)
			{
				kl = -k1 - 1;
				ZkpL[index] = z1;
			}
		}
		if (k2 < kmax2)
		{
			z2 = ZkpL[MapFW(j - 1, w2, k2)] + Inst.p[j - 1];
			w2 = w - Inst.w[j - 1];
			if (ZkpL[index] < z2)
			{
				kl = k2 + 1;
				ZkpL[index] = z2;
			}
		}
		PredL[index] = kl;
		if (kl < 0)
			k1++;
		else
			k2++;
		kk++;
	}
	KardL[MapFW(j, w)] = kk;

	return 0;
}

/**
 *  KPSolver: Allocate_KPFW_kBest_L
 */
int KPSolver::Allocate_KPFW_kBest_L(void)
{
	long k;

	// Setup parameters
	nnw = (Inst.n + 1)*(Inst.W + 1);
	nwk = (Inst.W + 1)*Inst.k;
	nk = Inst.k;
	nw = Inst.W + 1;
	nn = Inst.n + 1;

	ZkpL = new long[(Inst.n + 1)*(Inst.W + 1)*Inst.k];
	if (ZkpL == NULL)
		return 1;

	PredL = new long[(Inst.n + 1)*(Inst.W + 1)*Inst.k];
	if (PredL == NULL)
		return 1;

	KardL = new long[(Inst.n + 1)*(Inst.W + 1)];
	if (KardL == NULL)
		return 1;

	kBSol.nk = 0;
	kBSol.z = new long[Inst.k];
	if (kBSol.z == NULL)
		return 1;

	kBSol.x = new long*[Inst.k];
	if (kBSol.x == NULL)
		return 1;

	for (k = 0; k < Inst.k; k++)
	{
		kBSol.x[k] = new long[Inst.n];
		if (kBSol.x[k] == NULL)
			return 1;
	}

	return 0;
}

/**
 *  KPSolver: DeAllocate_KPFW_kBest_L
 */
int KPSolver::DeAllocate_KPFW_kBest_L(void)
{
	long k;

	if (ZkpL != NULL)
	{
		delete[] ZkpL;
		ZkpL = NULL;
	}

	if (PredL != NULL)
	{
		delete[] PredL;
		PredL = NULL;
	}

	if (KardL != NULL)
	{
		delete[] KardL;
		KardL = NULL;
	}

	if (kBSol.z != NULL)
	{
		delete[] kBSol.z;
		kBSol.z = NULL;
	}

	if (kBSol.x != NULL)
	{
		for (k = 0; k < Inst.k; k++)
			if (kBSol.x[k] != NULL)
				delete[] kBSol.x[k];
		delete[] kBSol.x;
		kBSol.x = NULL;
	}

	return 0;
}

/**
 *  KPSolver: Solve_KPFW_kBest_L
 */
int KPSolver::Solve_KPFW_kBest_L(void)
{
	long i, j;
	long w;
	long k, k1, ks;
	size_t index;
	//float t1, t2, dt;

	// Initialize "KardL"
	//for (j = 0; j <= Inst.n; j++)
	//	for (w = 0; w <= Inst.W; w++)
	//		KardL[MapFW(j, w)] = 0;
	//
	// The equivalent initialization
	for (i = 0; i < (Inst.n + 1)*(Inst.W + 1); i++)
		KardL[i] = 0;

	// Basic stage (i.e., i=0)
	for (w = 0; w <= Inst.W; w++)
	{
		index = MapFW(0, w, 0);
		ZkpL[index] = 0;
		PredL[index] = -1;
		KardL[MapFW(0, w)] = 1;
	}

	// Compute the stages from i=1 to n
	//t1 = GetTime();
	for (j = 1; j <= Inst.n; j++)
	{
		// Compute stage i+1
		for (w = 0; w <= Inst.W; w++)
		{
			Generate_State_L(j, w);
		}
	}
	//t2 = GetTime();
	//dt = t2 - t1;
	//printf("Generate_State Computing time: %7.4f \n", dt);

	// Build the k-Best solutions
	//t1 = t2;
	ks = KardL[MapFW(Inst.n, Inst.W)];
	kBSol.nk = ks;
	for (k = 0; k < ks; k++)
	{
		kBSol.z[k] = ZkpL[MapFW(Inst.n, Inst.W, k)];
		k1 = k;
		w = Inst.W;
		for (j = Inst.n; j > 0; j--)
		{
			index = MapFW(j, w, k1);
			if (PredL[index] > 0)
			{
				// The item i is in the solution
				kBSol.x[k][j - 1] = 1;
				k1 = PredL[index] - 1;
				w = w - Inst.w[j - 1];
			}
			else
			{
				// The item i is not in the solution
				kBSol.x[k][j - 1] = 0;
				k1 = -PredL[index] - 1;
			}
		}
	}
	//t2 = GetTime();
	//dt = t2 - t1;
	//printf("Build Computing time: %7.4f \n", dt);

	// Print the Zkp matrix
	if (PrintOn > 3)
	{
		for (k = 0; k < ks; k++)
			for (j = 0; j <= Inst.n; j++)
			{
				printf("Stage (%3d): ", j);
				for (w = 0; w <= Inst.W; w++)
				{
					printf("%5d ", ZkpL[MapFW(j, w, k)]);
				}
				printf("\n");
			}
	}

	// Print the optimal solution
	if (PrintOn > 0)
	{
		printf("Zkp(%d) = %d \n", 0, kBSol.z[0]);
		printf("Zkp(%d) = %d \n", ks - 1, kBSol.z[ks - 1]);
	}
	if (PrintOn > 1)
	{
		for (k = 0; k < ks; k++)
		{
			printf("Zkp(%d) = %5d \n", k, kBSol.z[k]);
			if (PrintOn > 2)
			{
				for (j = 0; j < Inst.n; j++)
					printf("x(%d,%d)=%d ", k, j, kBSol.x[k][j]);
				printf("\n");
			}
		}
	}

	return 0;
}

//-----------------------------------------------------------------------------
// Functions to map bi- and tri-dimensional arrays on one dimensional arrays 
//-----------------------------------------------------------------------------

/**
 *  KPSolver: MapBW
 */
size_t KPSolver::MapBW(long j, long w)
{
	return j * nw + w;
}

/**
 *  KPSolver: MapBW
 */
size_t KPSolver::MapBW(long j, long w, long k)
{
	//return j * nwk + w * nk + k;
	//return k * nnw + j * nw + w;
	return k * nnw + w * nn + j;  // Best!
}

//-----------------------------------------------------------------------------
// Functions to generate k-best solutions with backward recursion 
// based on one dimensional arrays 
//-----------------------------------------------------------------------------

/**
 *  KPSolver: Setup_KPBW_kBest_L
 */
void KPSolver::Setup_KPBW_kBest_L(void)
{
	//long w;
	long j;
	size_t index;

	// Define the minimum item width (esed to stop the search in advance)
	wmin = Inst.W;
	for (j = 0; j < Inst.n; j++)
		if (Inst.w[j] < wmin)
			wmin = Inst.w[j];

	//for (j = 0; j <= Inst.n; j++)
	//	for (w = 0; w <= Inst.W; w++)
	//	{
	//		PredL[MapBW(j,w,0)] = 0;
	//		index = MapBW(j, w);
	//		KardL[index] = 0;
	//		KardL1[index] = 0;
	//		KardL2[index] = 0;
	//	}
	//
	// The equivalent setup 
	for (index = 0; index < (Inst.n + 1)*(Inst.W + 1); index++)
	{
		KardL[index] = 0;
		KardL1[index] = 0;
		KardL2[index] = 0;
	}
}

/**
 *  KPSolver: KPBW_kBest_L
 */
long KPSolver::KPBW_kBest_L(long j, long w, long k)
{
	long z1, z2;
	long w2;
	size_t index2, index3;

	// Setup index
	index2 = MapBW(j, w);
	index3 = MapBW(j, w, k);

	// The current state has been already generated
	if (k < KardL[index2])
		return ZkpL[index3];

	// Stop criteria for the recursion
	//if ((j == 0) || (w == 0)) // Basic version
	if ((j == 0) || (w < wmin)) 
	{
		if (k == 0)
		{
			ZkpL[index3] = 0;
			PredL[index3] = 0;
			KardL[index2] = 1;
			return 0;
		}
		else
			return -Infi;
	}

	// Evaluate the alternative (include or not the item j in the optimal solution)
	z1 = KPBW_kBest_L(j - 1, w, KardL1[index2]);
	w2 = w - Inst.w[j - 1];
	if (w2 >= 0)
	{
		//z2 = KPBW_kBest_L(j - 1, w2, KardL2[MapBW(j, w)]) + Inst.p[j - 1];
		z2 = KPBW_kBest_L(j - 1, w2, KardL2[index2]);
		if (z2 > -Infi)
			z2 = z2 + Inst.p[j - 1];
	}
	else
		z2 = -Infi;
	if (z1 >= z2)
	{
		if (z1 == -Infi)
			return -Infi;

		ZkpL[index3] = z1;
		PredL[index3] = -(KardL1[index2] + 1);
		KardL1[index2]++;
	}
	else
	{
		ZkpL[index3] = z2;
		PredL[index3] = +(KardL2[index2] + 1);
		KardL2[index2]++;
	}
	KardL[index2] = k + 1;

	return ZkpL[index3];
}

/**
 *  KPSolver: Allocate_KPBW_kBest_L
 */
int KPSolver::Allocate_KPBW_kBest_L(void)
{
	long k;

	// Setup parameters
	nnw = (Inst.n + 1)*(Inst.W + 1); 
	nwk = (Inst.W + 1)*Inst.k;
	nk = Inst.k;
	nw = Inst.W + 1;
	nn = Inst.n + 1;

	ZkpL = new long[(Inst.n + 1) * (Inst.W + 1) * Inst.k];
	if (ZkpL == NULL)
		return 1;

	PredL = new long[(Inst.n + 1) * (Inst.W + 1) * Inst.k];
	if (PredL == NULL)
		return 1;

	KardL = new long[(Inst.n + 1) * (Inst.W + 1)];
	if (KardL == NULL)
		return 1;

	KardL1 = new long[(Inst.n + 1) * (Inst.W + 1)];
	if (KardL1 == NULL)
		return 1;

	KardL2 = new long[(Inst.n + 1) * (Inst.W + 1)];
	if (KardL2 == NULL)
		return 1;

	kBSol.nk = 0;
	kBSol.z = new long[Inst.k];
	if (kBSol.z == NULL)
		return 1;

	kBSol.x = new long*[Inst.k];
	if (kBSol.x == NULL)
		return 1;

	for (k = 0; k < Inst.k; k++)
	{
		kBSol.x[k] = new long[Inst.n];
		if (kBSol.x[k] == NULL)
			return 1;
	}

	return 0;
}

/**
 *  KPSolver: DeAllocate_KPBW_kBest_L
 */
int KPSolver::DeAllocate_KPBW_kBest_L(void)
{
	long k;

	if (ZkpL != NULL)
	{
		delete[] ZkpL;
		ZkpL = NULL;
	}

	if (PredL != NULL)
	{
		delete[] PredL;
		PredL = NULL;
	}

	if (KardL != NULL)
	{
		delete[] KardL;
		KardL = NULL;
	}

	if (KardL1 != NULL)
	{
		delete[] KardL1;
		KardL1 = NULL;
	}

	if (KardL2 != NULL)
	{
		delete[] KardL2;
		KardL2 = NULL;
	}

	if (kBSol.z != NULL)
	{
		delete[] kBSol.z;
		kBSol.z = NULL;
	}

	if (kBSol.x != NULL)
	{
		for (k = 0; k < Inst.k; k++)
			if (kBSol.x[k] != NULL)
				delete[] kBSol.x[k];
		delete[] kBSol.x;
		kBSol.x = NULL;
	}

	return 0;
}

/**
 *  KPSolver: Solve_KPBW_kBest_L
 */
int KPSolver::Solve_KPBW_kBest_L(void)
{
	long i, j;
	long k, kk;
	long w;
	size_t index;
	//float t1, t2, dt;

	// Setup data structures
	//t1 = GetTime();
	Setup_KPBW_kBest_L();
	//Setup_KPBW_kBest(Inst.n, Inst.W);
	//t2 = GetTime();
	//dt = t2 - t1;
	//printf("Setup Computing time: %7.4f \n", dt);

	// Build the k-best solutions
	k = 0;
	while (k < Inst.k)
	{
		kBSol.z[k] = KPBW_kBest_L(Inst.n, Inst.W, k);
		if (kBSol.z[k] <= -Infi)
			break;
		w = Inst.W;
		kk = k;
		for (i = Inst.n; i > 0; i--)
		{
			index = MapBW(i, w, kk);
			if (PredL[index] > 0)
			{
				// The item i is in the solution
				kBSol.x[k][i - 1] = 1;
				kk = PredL[index] - 1;
				w = w - Inst.w[i - 1];
			}
			else if (PredL[index] < 0)
			{
				// The item i is not in the solution
				kBSol.x[k][i - 1] = 0;
				kk = -PredL[index] - 1;
			}
			else
			{
				kBSol.x[k][i - 1] = 0;
				break;
			}
		}
		for (j = i; j > 0; j--)
			kBSol.x[k][j - 1] = 0;

		k++;
	}
	kBSol.nk = k;

	// Print the optimal solution
	if (PrintOn > 0)
	{
		printf("Zkp(%d) = %d \n", 0, kBSol.z[0]);
		printf("Zkp(%d) = %d \n", k - 1, kBSol.z[k - 1]);
	}
	if (PrintOn > 1)
	{
		for (kk = 0; kk < k; kk++)
		{
			printf("Zkp(%d) = %5d \n", kk, kBSol.z[kk]);
			if (PrintOn > 2)
			{
				for (i = 0; i < Inst.n; i++)
					printf("x(%d,%d)=%d ", kk, i, kBSol.x[kk][i]);
				printf("\n");
			}
		}
	}

	return 0;
}

/**
 *  KPSolver: KPI_KPBW_kBest_L
 */
int KPSolver::KPI_KPBW_kBest_L(void)
{
	long j, w;
	long Konta;
	long Konta1;
	long Konta2;
	long Konta3;
	long Konta4;
	long Konta5;
	long Konta6;
	long Konta7;
	float AvgK;
	float Perc, Perc1;

	Konta = 0;
	Konta1 = 0;
	Konta2 = 0;
	Konta3 = 0;
	Konta4 = 0;
	Konta5 = 0;
	Konta6 = 0;
	Konta7 = 0;
	for (j = 0; j <= Inst.n; j++)
		for (w = 0; w <= Inst.W; w++)
		{
			if (KardL[MapBW(j, w)] > 0)
				Konta++;
			if (KardL[MapBW(j, w)] > 10)
				Konta3++;
			if (KardL[MapBW(j, w)] > 50)
				Konta4++;
			if (KardL[MapBW(j, w)] > 100)
				Konta5++;
			if ((j > Inst.n / 2) && (KardL[MapBW(j, w)] > 0))
				Konta6++;
			if ((j > Inst.n / 4) && (KardL[MapBW(j, w)] > 0))
				Konta7++;
			Konta1 += KardL[MapBW(j, w)];
		}

	Konta2 = Konta1 / Konta;
	AvgK = (float)Konta1 / Konta;

	// To avoid long integer overflow we split each computation
	// If (Inst.n + 1), (Inst.W + 1), and (Inst.k) are too large we can have an overflow

	//Perc = (float)Konta / ((float)(Inst.n + 1)*(Inst.W + 1)) * 100.;
	Perc = (float)Konta / ((float)(Inst.n + 1));
	Perc = Perc / ((float)(Inst.W + 1));
	Perc = Perc * (float)100.;

	//Perc1 = (float)Konta1 / ((float)(Inst.n + 1)*(Inst.W + 1)*(Inst.k)) * 100.;
	Perc1 = (float)Konta1 / ((float)(Inst.n + 1));
	Perc1 = Perc1 / ((float)(Inst.W + 1));
	Perc1 = Perc1 / ((float)(Inst.k));
	Perc1 = Perc1 * (float)100.;

	KPI.Avg_k = AvgK;
	KPI.Num_States_Used = Konta;
	KPI.Perc_States_Used = Perc;
	KPI.Num_States_Generated = Konta1;
	KPI.Perc_States_Generated = Perc1;
	KPI.Avg_Dim_k = (float)Inst.k;  // It is not dynamic
	KPI.Perc_Dim_k = 100.0;  // It is not dynamic

	if (PrintOn == 1)
	{
		//printf("Number of states (j,w) with k>10: %d \n", Konta3);
		//printf("Number of states (j,w) with k>50: %d \n", Konta4);
		//printf("Number of states (j,w) with k>100: %d \n", Konta5);
		//printf("Number of states (j,w) with k>0 and j>n/2: %d \n", Konta6);
		//printf("Number of states (j,w) with k>0 and j>n/4: %d \n", Konta7);
		printf("States (j,w) used (active): %d // %d (%f%%) \n", Konta, (Inst.n + 1)*(Inst.W + 1), Perc);
		printf("States (j,w,k) generated: %d // %d  (%f%%) \n", Konta1, (Inst.n + 1)*(Inst.W + 1)*(Inst.k), Perc1);
		//printf("Average k for each active state (j,w): %d \n", Konta2);
		printf("Average k for each active state (j,w): %f \n", AvgK);
	}

	return 0;
}

//-----------------------------------------------------------------------------
// Functions to generate k-best solutions with backward recursion 
// based on one dimensional arrays with dynamic reallocation of states k
//-----------------------------------------------------------------------------

/**
 *  KPSolver: Setup_KPBW_kBest_LD
 */
void KPSolver::Setup_KPBW_kBest_LD(void)
{
	//long w;
	long j;
	size_t index;

	// Define the minimum item width (esed to stop the search in advance)
	wmin = Inst.W;
	for (j = 0; j < Inst.n; j++)
		if (Inst.w[j] < wmin)
			wmin = Inst.w[j];

	// Setup sizes for allocation
	dmax = Inst.k;  // Maximum size
	ddmin = Inst.n / PercMinDim;  //Minimum size increment
	if (ddmin < MinDim)
		ddmin = MinDim;
	ddmax = Inst.n / PercMaxDim;  //Maxmum size increment
	if (ddmax < MaxDim)
		ddmax = MaxDim;

	// Setup each state as unexplored
	//for (j = 0; j <= Inst.n; j++)
	//	for (w = 0; w <= Inst.W; w++)
	//	{
	//		PredL[MapBW(j,w,0)] = 0;
	//		index = MapBW(j, w);
	//		KardL[index] = 0;
	//		KardL1[index] = 0;
	//		KardL2[index] = 0;
	//	}
	//
	// The equivalent setup 
	for (index = 0; index < (Inst.n + 1)*(Inst.W + 1); index++)
	{
		KardL[index] = 0;
		KardL1[index] = 0;
		KardL2[index] = 0;
	}
}

/**
 *  KPSolver: KPBW_kBest_LD
 */
long KPSolver::KPBW_kBest_LD(long j, long w, long k)
{
	long z1, z2;
	long w2;
	size_t index2;

	// Setup index
	index2 = MapBW(j, w);

	// Check if array for states k are already allocated
	if (Dim[index2] == 0)
	{
		Dim[index2] = ddmin;
		ZkpLD[index2] = (long*)malloc(ddmin * sizeof(long));
		PredLD[index2] = (long*)malloc(ddmin * sizeof(long));
	}
	if (Dim[index2] <= k)
	{
		Dim[index2] = (k / ddmax + 1) *ddmax;
		if (Dim[index2] > dmax)
			Dim[index2] = dmax;
		ZkpLD[index2] = (long*)realloc(ZkpLD[index2], Dim[index2] * sizeof(long));
		PredLD[index2] = (long*)realloc(PredLD[index2], Dim[index2] * sizeof(long));
	}

	// The current state has been already generated
	if (k < KardL[index2])
		return ZkpLD[index2][k];

	// Stop criteria for the recursion
	//if ((j == 0) || (w == 0)) // Basic version
	if ((j == 0) || (w < wmin)) 
	{
		if (k == 0)
		{
			ZkpLD[index2][k] = 0;
			PredLD[index2][k] = 0;
			KardL[index2] = 1;
			return 0;
		}
		else
			return -Infi;
	}

	// Evaluate the alternative (include or not the item j in the optimal solution)
	z1 = KPBW_kBest_LD(j - 1, w, KardL1[index2]);
	w2 = w - Inst.w[j - 1];
	if (w2 >= 0)
	{
		//z2 = KPBW_kBest_L(j - 1, w2, KardL2[MapBW(j, w)]) + Inst.p[j - 1];
		z2 = KPBW_kBest_LD(j - 1, w2, KardL2[index2]);
		if (z2 > -Infi)
			z2 = z2 + Inst.p[j - 1];
	}
	else
		z2 = -Infi;
	if (z1 >= z2)
	{
		if (z1 == -Infi)
			return -Infi;

		ZkpLD[index2][k] = z1;
		PredLD[index2][k] = -(KardL1[index2] + 1);
		KardL1[index2]++;
	}
	else
	{
		ZkpLD[index2][k] = z2;
		PredLD[index2][k] = +(KardL2[index2] + 1);
		KardL2[index2]++;
	}
	KardL[index2] = k + 1;

	return ZkpLD[index2][k];
}

/**
 *  KPSolver: Allocate_KPBW_kBest_LD
 */
int KPSolver::Allocate_KPBW_kBest_LD(void)
{
	long k;
	size_t index;

	// Setup parameters
	nnw = (Inst.n + 1)*(Inst.W + 1); 
	nwk = (Inst.W + 1)*Inst.k;
	nk = Inst.k;
	nw = Inst.W + 1;
	nn = Inst.n + 1;

	ZkpLD = (long**)malloc((Inst.n + 1) * (Inst.W + 1) * sizeof(new long*));
	if (ZkpLD == NULL)
		return 1;

	PredLD = (long**)malloc((Inst.n + 1) * (Inst.W + 1) * sizeof(new long*));
	if (PredLD == NULL)
		return 1;

	Dim = (long*)malloc((Inst.n + 1) * (Inst.W + 1) * sizeof(new long));
	if (Dim == NULL)
		return 1;

	KardL = (long*)malloc((Inst.n + 1) * (Inst.W + 1) * sizeof(new long));
	if (KardL == NULL)
		return 1;

	KardL1 = (long*)malloc((Inst.n + 1) * (Inst.W + 1) * sizeof(new long));
	if (KardL1 == NULL)
		return 1;

	KardL2 = (long*)malloc((Inst.n + 1) * (Inst.W + 1) * sizeof(new long));
	if (KardL2 == NULL)
		return 1;

	for (index = 0; index < (Inst.n + 1)*(Inst.W + 1); index++)
	{
		ZkpLD[index] = NULL;
		PredLD[index] = NULL;
		Dim[index] = 0;
	}

	kBSol.nk = 0;
	kBSol.z = new long[Inst.k];
	if (kBSol.z == NULL)
		return 1;

	kBSol.x = new long*[Inst.k];
	if (kBSol.x == NULL)
		return 1;

	for (k = 0; k < Inst.k; k++)
	{
		kBSol.x[k] = new long[Inst.n];
		if (kBSol.x[k] == NULL)
			return 1;
	}

	return 0;
}

/**
 *  KPSolver: DeAllocate_KPBW_kBest_LD
 */
int KPSolver::DeAllocate_KPBW_kBest_LD(void)
{
	long k;
	size_t index;

	if (ZkpLD != NULL)
	{
		for (index = 0; index < (Inst.n + 1)*(Inst.W + 1); index++)
			if (ZkpLD[index] != NULL)
				free(ZkpLD[index]);
		free(ZkpLD);
		ZkpLD = NULL;
	}

	if (PredLD != NULL)
	{
		for (index = 0; index < (Inst.n + 1)*(Inst.W + 1); index++)
			if (PredLD[index] != NULL)
				free(PredLD[index]);
		free(PredLD);
		PredLD = NULL;
	}

	if (Dim != NULL)
	{
		free(Dim);
		Dim = NULL;
	}
	
	if (KardL != NULL)
	{
		free(KardL);
		KardL = NULL;
	}

	if (KardL1 != NULL)
	{
		free(KardL1);
		KardL1 = NULL;
	}

	if (KardL2 != NULL)
	{
		free(KardL2);
		KardL2 = NULL;
	}

	if (kBSol.z != NULL)
	{
		delete[] kBSol.z;
		kBSol.z = NULL;
	}

	if (kBSol.x != NULL)
	{
		for (k = 0; k < Inst.k; k++)
			if (kBSol.x[k] != NULL)
				delete[] kBSol.x[k];
		delete[] kBSol.x;
		kBSol.x = NULL;
	}

	return 0;
}

/**
 *  KPSolver: Solve_KPBW_kBest_LD
 */
int KPSolver::Solve_KPBW_kBest_LD(void)
{
	long i, j;
	long k, kk;
	long w;
	size_t index;
	//float t1, t2, dt;

	// Setup data structures
	//t1 = GetTime();
	Setup_KPBW_kBest_LD();
	//Setup_KPBW_kBest(Inst.n, Inst.W);
	//t2 = GetTime();
	//dt = t2 - t1;
	//printf("Setup Computing time: %7.4f \n", dt);

	// Build the k-best solutions
	k = 0;
	while (k < Inst.k)
	{
		kBSol.z[k] = KPBW_kBest_LD(Inst.n, Inst.W, k);
		if (kBSol.z[k] <= -Infi)
			break;
		w = Inst.W;
		kk = k;
		for (i = Inst.n; i > 0; i--)
		{
			index = MapBW(i, w);
			if (PredLD[index][kk] > 0)
			{
				// The item i is in the solution
				kBSol.x[k][i - 1] = 1;
				kk = PredLD[index][kk] - 1;
				w = w - Inst.w[i - 1];
			}
			else if (PredLD[index][kk] < 0)
			{
				// The item i is not in the solution
				kBSol.x[k][i - 1] = 0;
				kk = -PredLD[index][kk] - 1;
			}
			else
			{
				kBSol.x[k][i - 1] = 0;
				break;
			}
		}
		for (j = i; j > 0; j--)
			kBSol.x[k][j - 1] = 0;

		k++;
	}
	kBSol.nk = k;

	// Print the optimal solution
	if (PrintOn > 0)
	{
		printf("Zkp(%d) = %d \n", 0, kBSol.z[0]);
		printf("Zkp(%d) = %d \n", k - 1, kBSol.z[k - 1]);
	}
	if (PrintOn > 1)
	{
		for (kk = 0; kk < k; kk++)
		{
			printf("Zkp(%d) = %5d \n", kk, kBSol.z[kk]);
			if (PrintOn > 2)
			{
				for (i = 0; i < Inst.n; i++)
					printf("x(%d,%d)=%d ", kk, i, kBSol.x[kk][i]);
				printf("\n");
			}
		}
	}

	return 0;
}

//-----------------------------------------------------------------------------
// Functions to generate k-best solutions with backward recursion 
// based on one dimensional arrays with dynamic reallocation of states k
// This version manages further constraints
//-----------------------------------------------------------------------------

/**
 *  KPSolver: Setup_KPBW_kBest_LDC
 */
void KPSolver::Setup_KPBW_kBest_LDC(void)
{
	long j;
	size_t index;

	// Define the minimum item width (esed to stop the search in advance)
	wmin = Inst.W;
	for (j = 0; j < Inst.n; j++)
		if (Inst.w[j] < wmin)
			wmin = Inst.w[j];

	// Setup sizes for allocation
	// Maximum size
	if (Cons.Type == 0)
		dmax = Inst.k;
	else
		dmax = Infi;  // If we have further constraints is not limited

	//Minimum size increment
	ddmin = Inst.n / PercMinDim;
	if (ddmin < MinDim)
		ddmin = MinDim;

	//Maxmum size increment
	ddmax = Inst.n / PercMaxDim;
	if (ddmax < MaxDim)
		ddmax = MaxDim;

	for (index = 0; index < (Inst.n + 1)*(Inst.W + 1); index++)
	{
		KardL[index] = 0;
		KardL1[index] = 0;
		KardL2[index] = 0;
	}
}

/**
 *  KPSolver: FeasSol_kBest_LDC
 */
int KPSolver::FeasSol_kBest_LDC(long k)
{
	long i, j;
	long i1, j1;
	long W;
	long bol;

	// Here put the code to check the feasibility
	if (Cons.Type == 0)  // Included to save time for the unconstrained KP (and avoid 4 if statements)
	{
		return 1;
	}
	else if (Cons.Type == 1)
	{
		// Type 1: Cardinality constrained bin packing problem
		if (kBSol.nx[k] > Cons.MaxCard)
			return 0;
	}
	else if (Cons.Type == 2)
	{
		// Type 2: Class constrained bin packing problem
		W = 1;  // we count the last one in advance
		for (i = 0; i < kBSol.nx[k] - 1; i++)
		{
			i1 = kBSol.x[k][i];
			bol = 1;
			for (j = i + 1; j < kBSol.nx[k]; j++)
			{
				j1 = kBSol.x[k][j];
				if (Cons.cl[i] == Cons.cl[j])
				{
					bol = 0;
					break;
				}
			}
			if (bol)
				W++;
		}
		if (W > Cons.CL)
			return 0;
	}
	else if (Cons.Type == 3)
	{
		// Type 3: Bin packing problem with conflicts
		for (i = 0; i < kBSol.nx[k] - 1; i++)
		{
			i1 = kBSol.x[k][i];
			for (j = i + 1; j < kBSol.nx[k]; j++)
			{
				j1 = kBSol.x[k][j];
				if (Cons.CMatrix[i][j] > 0)
					return 0;
			}
		}
	}
	else if (Cons.Type == 4)
	{
		// Type 4: Two-dimensional vector packing problem
		W = 0;
		for (i = 0; i < kBSol.nx[k]; i++)
			W += Cons.w2[kBSol.x[k][i]];
		if (W > Cons.W2)
			return 0;
	}

	return 1;
}

/**
 *  KPSolver: BuildSolution_kBest_LD
 */
int KPSolver::BuildSolution_kBest_LD(long k, long k1)
{
	long i, j, kk;
	long w;
	size_t index;

	w = Inst.W;
	kk = k1;
	for (i = Inst.n; i > 0; i--)
	{
		index = MapBW(i, w);
		if (PredLD[index][kk] > 0)
		{
			// The item i is in the solution
			kBSol.x[k][i - 1] = 1;
			kk = PredLD[index][kk] - 1;
			w = w - Inst.w[i - 1];
		}
		else if (PredLD[index][kk] < 0)
		{
			// The item i is not in the solution
			kBSol.x[k][i - 1] = 0;
			kk = -PredLD[index][kk] - 1;
		}
		else
		{
			kBSol.x[k][i - 1] = 0;
			break;
		}
	}
	for (j = i; j > 0; j--)
		kBSol.x[k][j - 1] = 0;

	return 0;
}

/**
 *  KPSolver: BuildSolution_kBest_LDC
 */
int KPSolver::BuildSolution_kBest_LDC(long k, long k1)
{
	long i, kk;
	long w;
	long konta;
	size_t index;

	w = Inst.W;
	kk = k1;
	konta = 0;
	for (i = Inst.n; i > 0; i--)
	{
		index = MapBW(i, w);
		if (PredLD[index][kk] > 0)
		{
			// The item i is in the solution
			kBSol.x[k][konta] = i - 1;
			kk = PredLD[index][kk] - 1;
			w = w - Inst.w[i - 1];
			konta++;
		}
		else if (PredLD[index][kk] < 0)
		{
			// The item i is not in the solution
			kk = -PredLD[index][kk] - 1;
		}
		else
		{
			break;
		}
	}

	kBSol.nx[k] = konta;

	return 0;
}

/**
 *  KPSolver: Allocate_KPBW_kBest_LDC
 */
int KPSolver::Allocate_KPBW_kBest_LDC(void)
{
	long k;
	size_t index;

	// Setup Parameters
	nnw = (Inst.n + 1)*(Inst.W + 1); 
	nwk = (Inst.W + 1)*Inst.k;
	nk = Inst.k;
	nw = Inst.W + 1;
	nn = Inst.n + 1;

	ZkpLD = (long**)malloc((Inst.n + 1) * (Inst.W + 1) * sizeof(long*));
	if (ZkpLD == NULL)
		return 1;

	PredLD = (long**)malloc((Inst.n + 1) * (Inst.W + 1) * sizeof(long*));
	if (PredLD == NULL)
		return 1;

	Dim = (long*)malloc((Inst.n + 1) * (Inst.W + 1) * sizeof(long));
	if (Dim == NULL)
		return 1;

	KardL = (long*)malloc((Inst.n + 1) * (Inst.W + 1) * sizeof(long));
	if (KardL == NULL)
		return 1;

	KardL1 = (long*)malloc((Inst.n + 1) * (Inst.W + 1) * sizeof(long));
	if (KardL1 == NULL)
		return 1;

	KardL2 = (long*)malloc((Inst.n + 1) * (Inst.W + 1) * sizeof(long));
	if (KardL2 == NULL)
		return 1;

	for (index = 0; index < (Inst.n + 1)*(Inst.W + 1); index++)
	{
		ZkpLD[index] = NULL;
		PredLD[index] = NULL;
		Dim[index] = 0;
	}

	kBSol.nk = 0;
	kBSol.z = new long[Inst.k];
	if (kBSol.z == NULL)
		return 1;

	kBSol.nx = new long[Inst.k];
	if (kBSol.nx == NULL)
		return 1;

	kBSol.x = new long*[Inst.k];
	if (kBSol.x == NULL)
		return 1;

	for (k = 0; k < Inst.k; k++)
	{
		kBSol.x[k] = new long[Inst.n];
		if (kBSol.x[k] == NULL)
			return 1;
	}

	return 0;
}

/**
 *  KPSolver: DeAllocate_KPBW_kBest_LDC
 */
int KPSolver::DeAllocate_KPBW_kBest_LDC(void)
{
	long k;
	size_t index;

	if (ZkpLD != NULL)
	{
		for (index = 0; index < (Inst.n + 1)*(Inst.W + 1); index++)
			if (ZkpLD[index] != NULL)
				free(ZkpLD[index]);
		free(ZkpLD);
		ZkpLD = NULL;
	}

	if (PredLD != NULL)
	{
		for (index = 0; index < (Inst.n + 1)*(Inst.W + 1); index++)
			if (PredLD[index] != NULL)
				free(PredLD[index]);
		free(PredLD);
		PredLD = NULL;
	}

	if (Dim != NULL)
	{
		free(Dim);
		Dim = NULL;
	}

	if (KardL)
	{
		free(KardL);
		KardL = NULL;
	}

	if (KardL1 != NULL)
	{
		free(KardL1);
		KardL1 = NULL;
	}

	if (KardL2 != NULL)
	{
		free(KardL2);
		KardL2 = NULL;
	}

	if (kBSol.z != NULL)
	{
		delete[] kBSol.z;
		kBSol.z = NULL;
	}

	if (kBSol.nx != NULL)
	{
		delete[] kBSol.nx;
		kBSol.nx = NULL;
	}

	if (kBSol.x != NULL)
	{
		for (k = 0; k < Inst.k; k++)
			if (kBSol.x[k] != NULL)
				delete[] kBSol.x[k];
		delete[] kBSol.x;
		kBSol.x = NULL;
	}

	return 0;
}

/**
 *  KPSolver: Solve_KPBW_kBest_LDC
 */
int KPSolver::Solve_KPBW_kBest_LDC(void)
{
	long i;
	long k, kk;
	//float t1, t2, dt;

	// Setup data structures
	//t1 = GetTime();
	Setup_KPBW_kBest_LDC();
	//t2 = GetTime();
	//dt = t2 - t1;
	//printf("Setup Computing time: %7.4f \n", dt);

	// Build the k-best solutions
	k = 0;
	kk = 0;
	while (k < Inst.k)
	{
		kBSol.z[k] = KPBW_kBest_LD(Inst.n, Inst.W, kk);
		if (kBSol.z[k] <= -Infi)
			break;
		BuildSolution_kBest_LDC(k, kk);
		if (FeasSol_kBest_LDC(k))
			k++;
		kk++;
	}
	kBSol.nk = k;

	// Print the optimal solution
	if (PrintOn > 0)
	{
		printf("Zkp(%d) = %d \n", 0, kBSol.z[0]);
		printf("Zkp(%d) = %d \n", k - 1, kBSol.z[k - 1]);
	}
	if (PrintOn > 1)
	{
		printf("Numero soluzioni generate = %5d \n", kk);
		for (kk = 0; kk < k; kk++)
		{
			printf("Zkp(%d) = %5d \n", kk, kBSol.z[kk]);
			if (PrintOn > 2)
			{
				for (i = 0; i < kBSol.nx[kk]; i++)
					printf("x(%d,%d)=%d ", kk, i, kBSol.x[kk][i]);
				printf("\n");
			}
		}
	}

	return 0;
}

//-----------------------------------------------------------------------------
// Functions to generate k-best solutions with backward recursion 
// based on one dimensional arrays with dynamic reallocation of states k
// This version manages further constraints and use float for Zkp
//-----------------------------------------------------------------------------

/**
 *  KPSolver: KPBW_kBest_LDCF
 */
float KPSolver::KPBW_kBest_LDCF(long j, long w, long k, int *status)
{
	float z1, z2;
	long w2;
	size_t index2;
	int status1, status2;

	// Setup status
	*status = 0;

	// Setup index
	index2 = MapBW(j, w);

	// Check if array for states k are already allocated
	if (Dim[index2] == 0)
	{
		Dim[index2] = ddmin;
		ZkpLDF[index2] = (float*)malloc(ddmin * sizeof(float));
		if (ZkpLDF[index2] == NULL)
		{
			*status = 1;
			return -Infr;
		}
		PredLD[index2] = (long*)malloc(ddmin * sizeof(long));
		if (PredLD[index2] == NULL)
		{
			*status = 1;
			return -Infr;
		}
	}
	if (Dim[index2] <= k)
	{
		Dim[index2] = (k / ddmax + 1) *ddmax;
		if (Dim[index2] > dmax)
			Dim[index2] = dmax;
		ZkpLDF[index2] = (float*)realloc(ZkpLDF[index2], Dim[index2] * sizeof(float));
		if (ZkpLDF[index2] == NULL)
		{
			*status = 1;
			return -Infr;
		}
		PredLD[index2] = (long*)realloc(PredLD[index2], Dim[index2] * sizeof(long));
		if (PredLD[index2] == NULL)
		{
			*status = 1;
			return -Infr;
		}
	}

	// The current state has been already generated
	if (k < KardL[index2])
		return ZkpLDF[index2][k];

	// Stop criteria for the recursion
	//if ((j == 0) || (w == 0)) // Basic version
	if ((j == 0) || (w < wmin))
	{
		if (k == 0)
		{
			ZkpLDF[index2][k] = 0.0;
			PredLD[index2][k] = 0;
			KardL[index2] = 1;
			return 0.0;
		}
		else
		{
			*status = -1;
			return -Infr;
		}
	}

	// Evaluate the alternative (include or not the item j in the optimal solution)
	z1 = KPBW_kBest_LDCF(j - 1, w, KardL1[index2], status);
	if (*status > 0)
		return -Infr;
	status1 = *status;

	w2 = w - Inst.w[j - 1];
	if (w2 >= 0)
	{
		//z2 = KPBW_kBest_L(j - 1, w2, KardL2[MapBW(j, w)]) + Inst.p[j - 1];
		z2 = KPBW_kBest_LDCF(j - 1, w2, KardL2[index2], status);
		if (*status > 0)
			return -Infr;
		if (*status == 0)
			z2 = z2 + Inst.pf[j - 1];
		status2 = *status;
	}
	else
	{
		z2 = -Infr;
		status2 = -1;
	}

	if (z1 >= z2)
	{
		if (status1 < 0)
		{
			*status = -1;
			return -Infr;
		}
		ZkpLDF[index2][k] = z1;
		PredLD[index2][k] = -(KardL1[index2] + 1);
		KardL1[index2]++;
	}
	else
	{
		if (status2 < 0)
		{
			*status = -1;
			return -Infr;
		}
		ZkpLDF[index2][k] = z2;
		PredLD[index2][k] = +(KardL2[index2] + 1);
		KardL2[index2]++;
	}
	KardL[index2] = k + 1;
	*status = 0;

	return ZkpLDF[index2][k];
}

/**
 *  KPSolver: Allocate_KPBW_kBest_LDCF
 */
int KPSolver::Allocate_KPBW_kBest_LDCF(void)
{
	long k;
	size_t index;

	// Setup parameters
	nnw = (Inst.n + 1)*(Inst.W + 1); 
	nwk = (Inst.W + 1)*Inst.k;
	nk = Inst.k;
	nw = Inst.W + 1;
	nn = Inst.n + 1;

	ZkpLDF = (float**)malloc((Inst.n + 1) * (Inst.W + 1) * sizeof(float*));
	if (ZkpLDF == NULL)
		return 1;

	PredLD = (long**)malloc((Inst.n + 1) * (Inst.W + 1) * sizeof(long*));
	if (PredLD == NULL)
		return 1;

	Dim = (long*)malloc((Inst.n + 1) * (Inst.W + 1) * sizeof(long));
	if (Dim == NULL)
		return 1;

	KardL = (long*)malloc((Inst.n + 1) * (Inst.W + 1) * sizeof(long));
	if (KardL == NULL)
		return 1;

	KardL1 = (long*)malloc((Inst.n + 1) * (Inst.W + 1) * sizeof(long));
	if (KardL1 == NULL)
		return 1;

	KardL2 = (long*)malloc((Inst.n + 1) * (Inst.W + 1) * sizeof(long));
	if (KardL2 == NULL)
		return 1;

	for (index = 0; index < (Inst.n + 1)*(Inst.W + 1); index++)
	{
		ZkpLDF[index] = NULL;
		PredLD[index] = NULL;
		Dim[index] = 0;
	}

	kBSol.nk = 0;
	kBSol.zf = (float*)malloc((Inst.k) * sizeof(float));
	if (kBSol.zf == NULL)
		return 1;

	kBSol.nx = (long*)malloc((Inst.k) * sizeof(long));
	if (kBSol.nx == NULL)
		return 1;

	kBSol.x = (long**)malloc((Inst.k) * sizeof(long*));
	if (kBSol.x == NULL)
		return 1;
	for (k = 0; k < Inst.k; k++)
	{
		kBSol.x[k] = (long*)malloc((Inst.n) * sizeof(long));
		if (kBSol.x[k] == NULL)
			return 1;
	}

	return 0;
}

/**
 *  KPSolver: DeAllocate_KPBW_kBest_LDCF
 */
int KPSolver::DeAllocate_KPBW_kBest_LDCF(void)
{
	long k;
	size_t index;

	if (ZkpLDF != NULL)
	{
		for (index = 0; index < (Inst.n + 1)*(Inst.W + 1); index++)
		{
			if (ZkpLDF[index] != NULL)
				free(ZkpLDF[index]);
		}
		free(ZkpLDF);
		ZkpLDF = NULL;
	}

	if (PredLD != NULL)
	{
		for (index = 0; index < (Inst.n + 1)*(Inst.W + 1); index++)
		{
			if (PredLD[index] != NULL)
				free(PredLD[index]);
		}
		free(PredLD);
		PredLD = NULL;
	}

	if (Dim != NULL)
	{
		free(Dim);
		Dim = NULL;
	}

	if (KardL != NULL)
	{
		free(KardL);
		KardL = NULL;
	}

	if (KardL1 != NULL)
	{
		free(KardL1);
		KardL1 = NULL;
	}

	if (KardL2 != NULL)
	{
		free(KardL2);
		KardL2 = NULL;
	}

	if (kBSol.zf != NULL)
	{
		free(kBSol.zf);
		kBSol.zf = NULL;
	}

	if (kBSol.nx != NULL)
	{
		free(kBSol.nx);
		kBSol.nx = NULL;
	}

	if (kBSol.x != NULL)
	{
		for (k = 0; k < Inst.k; k++)
			free(kBSol.x[k]);
		free(kBSol.x);
		kBSol.x = NULL;
	}

	return 0;
}

/**
 *  KPSolver: Solve_KPBW_kBest_LDCF
 */
int KPSolver::Solve_KPBW_kBest_LDCF(void)
{
	long i;
	long k, kk;
	int status;
	//float t1, t2, dt;

	// Setup data structures
	//t1 = GetTime();
	Setup_KPBW_kBest_LDC();
	//t2 = GetTime();
	//dt = t2 - t1;
	//printf("Setup Computing time: %7.4f \n", dt);

	// Build the k-best solutions
	k = 0;
	kk = 0;
	while (k < Inst.k)
	{
		kBSol.zf[k] = KPBW_kBest_LDCF(Inst.n, Inst.W, kk, &status);
		if (status > 0)
			return status;
		if (status < 0)
			break;
		BuildSolution_kBest_LDC(k, kk);
		if (FeasSol_kBest_LDC(k))
			k++;
		kk++;
	}
	kBSol.nk = k;

	// Print the optimal solution
	if (PrintOn > 0)
	{
		printf("Zkp(%d) = %7.2f \n", 0, kBSol.zf[0]);
		printf("Zkp(%d) = %7.2f \n", k - 1, kBSol.zf[k - 1]);
	}
	if (PrintOn > 1)
	{
		printf("Numero soluzioni generate = %5d \n", kk);
		for (kk = 0; kk < k; kk++)
		{
			printf("Zkp(%d) = %7.2f \n", kk, kBSol.zf[kk]);
			if (PrintOn > 2)
			{
				for (i = 0; i < kBSol.nx[kk]; i++)
					printf("x(%d,%d)=%d ", kk, i, kBSol.x[kk][i]);
				printf("\n");
			}
		}
	}

	return 0;
}

/**
 *  KPSolver: Solve_KPBW_kBest_LDCF_CG
 *  (specifically designed for Column Generation)
 */
int KPSolver::Solve_KPBW_kBest_LDCF_CG(long kini, float delta)
{
	long i;
	long k, kk;
	int status;
	//float t1, t2, dt;

	// Setup data structures
	//t1 = GetTime();
	if (kini == 0)  // Otherwise we continue the previous generation
		Setup_KPBW_kBest_LDC();
	//t2 = GetTime();
	//dt = t2 - t1;
	//printf("Setup Computing time: %7.4f \n", dt);

	// Build the k-best solutions
	k = 0;
	kk = kini;  // Restart column generation from kini (this feature allows the core extension
	while (k < Inst.k)
	{
		kBSol.zf[k] = KPBW_kBest_LDCF(Inst.n, Inst.W, kk, &status);
		if (status > 0)
			return status;
		if (status < 0)
			break;
		BuildSolution_kBest_LDC(k, kk);
		if (FeasSol_kBest_LDC(k))
			k++;
		
		// Check delta
		if (kBSol.zf[kk] < delta)
			break;

		kk++;
	}
	kBSol.nk = k;

	// Print the optimal solution
	if (PrintOn > 0)
	{
		printf("Zkp(%d) = %7.2f \n", 0, kBSol.zf[0]);
		printf("Zkp(%d) = %7.2f \n", k - 1, kBSol.zf[k - 1]);
	}
	if (PrintOn > 1)
	{
		printf("Numero soluzioni generate = %5d \n", kk);
		for (kk = 0; kk < k; kk++)
		{
			printf("Zkp(%d) = %7.2f \n", kk, kBSol.zf[kk]);
			if (PrintOn > 2)
			{
				for (i = 0; i < kBSol.nx[kk]; i++)
					printf("x(%d,%d)=%d ", kk, i, kBSol.x[kk][i]);
				printf("\n");
			}
		}
	}

	return 0;
}

/**
 *  KPSolver: KPI_KPBW_kBest_LDCF
 */
int KPSolver::KPI_KPBW_kBest_LDCF(void)
{
	long j, w;
	size_t index2;
	long Konta;
	long Konta1;
	long Konta2;
	long Konta3;
	long Konta4;
	long Konta5;
	long Konta6;
	long Konta7;
	float AvgK;
	float TotDim;
	float AvgDim;
	float Perc, Perc1, Perc2;

	Konta = 0;
	Konta1 = 0;
	Konta2 = 0;
	Konta3 = 0;
	Konta4 = 0;
	Konta5 = 0;
	Konta6 = 0;
	Konta7 = 0;
	TotDim = 0.0;
	for (j = 0; j <= Inst.n; j++)
		for (w = 0; w <= Inst.W; w++)
		{
			// Setup index
			index2 = MapBW(j, w);

			if (KardL[index2] > 0)
				Konta++;
			if (KardL[index2] > 10)
				Konta3++;
			if (KardL[index2] > 50)
				Konta4++;
			if (KardL[index2] > 100)
				Konta5++;
			if ((j > Inst.n / 2) && (KardL[index2] > 0))
				Konta6++;
			if ((j > Inst.n / 4) && (KardL[index2] > 0))
				Konta7++;
			Konta1 += KardL[index2];

			TotDim += (float)Dim[index2];
		}

	Konta2 = Konta1 / Konta;
	AvgK = (float) Konta1 / Konta;

	// To avoid long integer overflow we split each computation
	// If (Inst.n + 1), (Inst.W + 1), and (Inst.k) are too large we can have an overflow

	//Perc = (float)Konta / ((Inst.n + 1)*(Inst.W + 1)) * 100.;
	Perc = (float)Konta / (Inst.n + 1);
	Perc = Perc / (Inst.W + 1);
	Perc = Perc * (float)100.;

	//Perc1 = (float)Konta1 / ((Inst.n + 1)*(Inst.W + 1)*(Inst.k)) * 100.;
	Perc1 = (float)Konta1 / (Inst.n + 1);
	Perc1 = Perc1 / (Inst.W + 1);
	Perc1 = Perc1 / (Inst.k);
	Perc1 = Perc1 * (float)100.;

	//AvgDim = (float)TotDim / ((Inst.n + 1)*(Inst.W + 1));
	AvgDim = (float)TotDim / (Inst.n + 1);
	AvgDim = AvgDim / (Inst.W + 1);

	//Perc2 = (float)TotDim / ((Inst.n + 1)*(Inst.W + 1)*(Inst.k)) * 100.;
	Perc2 = (float)TotDim / (Inst.n + 1);
	Perc2 = Perc2 / (Inst.W + 1);
	Perc2 = Perc2 / (Inst.k);
	Perc2 = Perc2 * (float)100.;

	// Setup KPI
	KPI.Avg_k = AvgK;
	KPI.Num_States_Used = Konta;
	KPI.Perc_States_Used = Perc;
	KPI.Num_States_Generated = Konta1;
	KPI.Perc_States_Generated = Perc1;
	KPI.Avg_Dim_k = AvgDim;
	KPI.Perc_Dim_k = Perc2;

	// If (Inst.n + 1), (Inst.W + 1), and (Inst.k) are too large we can have an overflow
	if (PrintOn == 1)
	{
		//printf("Numbero of states (j,w) with k>10: %d \n", Konta3);
		//printf("Numbero of states (j,w) with k>50: %d \n", Konta4);
		//printf("Numbero of states (j,w) with k>100: %d \n", Konta5);
		//printf("Numbero of states (j,w) with k>0 and j>n/2: %d \n", Konta6);
		//printf("Numbero of states (j,w) with k>0 and j>n/4: %d \n", Konta7);
		printf("States (j,w) used (active): %d // %d (%f%%) \n", Konta, (Inst.n + 1)*(Inst.W + 1), Perc);
		printf("States (j,w,k) generated: %d // %d  (%f%%) \n", Konta1, (Inst.n + 1)*(Inst.W + 1)*(Inst.k), Perc1);
		//printf("Average k for each active state (j,w): %d \n", Konta2);
		printf("Average k for each active state (j,w): %f \n", AvgK);
		printf("Average Dim(j,w): %f \n", AvgDim);
		printf("Percentage of memory allocated to save the k-best: %f%% \n", Perc2);
	}

	return 0;
}

//-----------------------------------------------------------------------------
// Knapsack Problem Bounded and Unbounded version
// Functions to generate k-best solutions with backward recursion 
// based on one dimensional arrays with dynamic reallocation of states k
// This version manage further constraints and use float for Zkp
//-----------------------------------------------------------------------------

/**
 *  KPSolver: Setup_KPBW_Bound_kBest_LDCF
 */
int KPSolver::Setup_KPBW_Bound_kBest_LDCF(void)
{
	long j;
	size_t index;

	// Define the minimum item width (esed to stop the search in advance)
	wmin = Inst.W;
	for (j = 0; j < Inst.n; j++)
		if (Inst.w[j] < wmin)
			wmin = Inst.w[j];
	if (wmin <= 0)
		return 1;

	// Setup sizes for allocation

	// Maximum size
	if (Cons.Type == 0)
		dmax = Inst.k;
	else
		dmax = Infi;  // If we have further constraints is not limited

	//Minimum size increment
	ddmin = Inst.n / 100;
	if (ddmin < 5)
		ddmin = 5;

	//Maxmum size increment
	ddmax = Inst.n / 10;
	if (ddmax < 10)
		ddmax = 10;

	for (index = 0; index < (Inst.n + 1)*(Inst.W + 1); index++)
	{
		KardL[index] = 0;
	}

	return 0;
}

/**
 *  KPSolver: KPBW_Bound_kBest_LDCF
 */
float KPSolver::KPBW_Bound_kBest_LDCF(long j, long w, long k, int *status)
{
	float z1, z2;
	long w2;
	long q, q1;
	size_t index2;
	int status1, status2;

	// Setup status
	*status = 0;

	// Setup index
	index2 = MapBW(j, w);

	// Check if array for states k are already allocated
	if (Dim[index2] == 0)
	{
		Dim[index2] = ddmin;
		ZkpLDF[index2] = (float*)malloc(ddmin * sizeof(float));
		if (ZkpLDF[index2] == NULL)
		{
			*status = 1;
			return -Infr;
		}
		PredLD[index2] = (long*)malloc(ddmin * sizeof(long));
		if (PredLD[index2] == NULL)
		{
			*status = 1;
			return -Infr;
		}
		PredALD[index2] = (long*)malloc(ddmin * sizeof(long));
		if (PredLD[index2] == NULL)
		{
			*status = 1;
			return -Infr;
		}
		Bound[index2] = 0;
		if (j > 0)
		{
			Bound[index2] = w / Inst.w[j - 1];
			if (Bound[index2] > Inst.b[j - 1])
				Bound[index2] = Inst.b[j - 1];
		}
		if (Bound[index2] >= 0)
		{
			KardB[index2] = (long*)malloc((Bound[index2] + 1) * sizeof(long));
			if (KardB[index2] == NULL)
			{
				*status = 1;
				return -Infr;
			}
			for (q = 0; q <= Bound[index2]; q++)
				KardB[index2][q] = 0;
		}
	}
	if (Dim[index2] <= k)
	{
		Dim[index2] = (k / ddmax + 1) *ddmax;
		if (Dim[index2] > dmax)
			Dim[index2] = dmax;
		ZkpLDF[index2] = (float*)realloc(ZkpLDF[index2], Dim[index2] * sizeof(float));
		if (ZkpLDF[index2] == NULL)
		{
			*status = 1;
			return -Infr;
		}
		PredLD[index2] = (long*)realloc(PredLD[index2], Dim[index2] * sizeof(long));
		if (PredLD[index2] == NULL)
		{
			*status = 1;
			return -Infr;
		}
		PredALD[index2] = (long*)realloc(PredALD[index2], Dim[index2] * sizeof(long));
		if (PredLD[index2] == NULL)
		{
			*status = 1;
			return -Infr;
		}
	}

	// The current state has been already generated
	if (k < KardL[index2])
		return ZkpLDF[index2][k];

	// Stop criteria for the recursion
	//if ((j == 0) || (w == 0)) // Basic version
	if ((j == 0) || (w < wmin))
	{
		if (k == 0)
		{
			ZkpLDF[index2][k] = 0.0;
			PredLD[index2][k] = 0;
			KardL[index2] = 1;
			return 0.0;
		}
		else
		{
			*status = -1;
			return -Infr;
		}
	}

	// Evaluate the alternative (include or not the item j in the optimal solution)
	z1 = -Infr;
	status1 = -1;
	q1 = -1;

	for (q = 0; q <= Bound[index2]; q++)
	{
		w2 = w - q * Inst.w[j - 1];
		if (w2 >= 0)
		{
			//z2 = KPBW_kBest_L(j - 1, w2, KardL2[MapBW(j, w)]) + Inst.p[j - 1];
			z2 = KPBW_Bound_kBest_LDCF(j - 1, w2, KardB[index2][q], status);
			if (*status > 0)
				return -Infr;
			if (*status == 0)
			{
				z2 = z2 + q * Inst.pf[j - 1];
				if ((z2 > z1) || (status1 == -1))
				{
					z1 = z2;
					status1 = 0;
					q1 = q;
				}
			}
			status2 = *status;
		}
	}

	if (status1 < 0)
	{
		*status = -1;
		return -Infr;
	}
	ZkpLDF[index2][k] = z1;
	PredLD[index2][k] = KardB[index2][q1] + 1;
	PredALD[index2][k] = q1;
	KardB[index2][q1]++;

	KardL[index2] = k + 1;
	*status = 0;

	return ZkpLDF[index2][k];
}

/**
 *  KPSolver: FeasSol_Bound_kBest_LDCF
 */
int KPSolver::FeasSol_Bound_kBest_LDCF(long k)
{
	long i;
	long Konta;

	// Here put the code to check the feasibility
	// Type 1: cardinality constraint (example)
	if (Cons.Type == 1)
	{
		Konta = 0;
		for (i = 0; i < kBSol.nx[k]; i++)
			Konta += kBSol.x[k][i];
		if (Konta > Cons.MaxCard)
			return 0;
	}

	return 1;
}

/**
 *  KPSolver: BuildSolution_Bound_kBest_LDCF
 */
int KPSolver::BuildSolution_Bound_kBest_LDCF(long k, long k1)
{
	long i, kk;
	long w;
	long q;
	long konta;
	size_t index;

	w = Inst.W;
	kk = k1;
	konta = 0;
	for (i = Inst.n; i > 0; i--)
	{
		index = MapBW(i, w);
		if (PredLD[index][kk] > 0)
		{
			// The item i is in the solution
			q = PredALD[index][kk];
			if (q > 0)
			{
				w = w - q * Inst.w[i - 1];
				kBSol.id[k][konta] = i - 1;
				kBSol.x[k][konta] = q;
				konta++;
			}
			kk = PredLD[index][kk] - 1;
		}
		else
		{
			break;
		}
	}

	kBSol.nx[k] = konta;

	return 0;
}

/**
 *  KPSolver: Allocate_KPBW_Bound_kBest_LDCF
 */
int KPSolver::Allocate_KPBW_Bound_kBest_LDCF(void)
{
	long k;
	size_t index;

	// Setup parameters
	nnw = (Inst.n + 1)*(Inst.W + 1); 
	nwk = (Inst.W + 1)*Inst.k;
	nk = Inst.k;
	nw = Inst.W + 1;
	nn = Inst.n + 1;

	ZkpLDF = (float**)malloc((Inst.n + 1) * (Inst.W + 1) * sizeof(float*));
	if (ZkpLDF == NULL)
		return 1;

	PredLD = (long**)malloc((Inst.n + 1) * (Inst.W + 1) * sizeof(long*));
	if (PredLD == NULL)
		return 1;

	PredALD = (long**)malloc((Inst.n + 1) * (Inst.W + 1) * sizeof(long*));
	if (PredALD == NULL)
		return 1;

	Dim = (long*)malloc((Inst.n + 1) * (Inst.W + 1) * sizeof(long));
	if (Dim == NULL)
		return 1;

	KardL = (long*)malloc((Inst.n + 1) * (Inst.W + 1) * sizeof(long));
	if (KardL == NULL)
		return 1;

	KardB = (long**)malloc((Inst.n + 1) * (Inst.W + 1) * sizeof(long*));
	if (KardB == NULL)
		return 1;

	Bound = (long*)malloc((Inst.n + 1) * (Inst.W + 1) * sizeof(long));
	if (Bound == NULL)
		return 1;

	for (index = 0; index < (Inst.n + 1) * (Inst.W + 1); index++)
	{
		ZkpLDF[index] = NULL;
		PredLD[index] = NULL;
		PredALD[index] = NULL;
		KardB[index] = NULL;
		Dim[index] = 0;
	}

	kBSol.nk = 0;
	kBSol.zf = (float*)malloc((Inst.k) * sizeof(float));
	if (kBSol.zf == NULL)
		return 1;

	kBSol.nx = (long*)malloc((Inst.k) * sizeof(long));
	if (kBSol.nx == NULL)
		return 1;

	kBSol.x = (long**)malloc((Inst.k) * sizeof(long*));
	if (kBSol.x == NULL)
		return 1;
	kBSol.id = (long**)malloc((Inst.k) * sizeof(long*));
	if (kBSol.id == NULL)
		return 1;
	for (k = 0; k < Inst.k; k++)
	{
		kBSol.x[k] = (long*)malloc((Inst.n) * sizeof(long));
		if (kBSol.x[k] == NULL)
			return 1;
		kBSol.id[k] = (long*)malloc((Inst.n) * sizeof(long));
		if (kBSol.id[k] == NULL)
			return 1;
	}

	return 0;
}

/**
 *  KPSolver: DeAllocate_KPBW_Bound_kBest_LDCF
 */
int KPSolver::DeAllocate_KPBW_Bound_kBest_LDCF(void)
{
	long k;
	size_t index;

	if (ZkpLDF != NULL)
	{
		for (index = 0; index < (Inst.n + 1)*(Inst.W + 1); index++)
		{
			if (ZkpLDF[index] != NULL)
				free(ZkpLDF[index]);
		}
		free(ZkpLDF);
		ZkpLDF = NULL;
	}

	if (PredLD != NULL)
	{
		for (index = 0; index < (Inst.n + 1)*(Inst.W + 1); index++)
		{
			if (PredLD[index] != NULL)
				free(PredLD[index]);
		}
		free(PredLD);
		PredLD = NULL;
	}

	if (PredALD != NULL)
	{
		for (index = 0; index < (Inst.n + 1)*(Inst.W + 1); index++)
		{
			if (PredALD[index] != NULL)
				free(PredALD[index]);
		}
		free(PredALD);
		PredALD = NULL;
	}

	if (KardB != NULL)
	{
		for (index = 0; index < (Inst.n + 1)*(Inst.W + 1); index++)
		{
			if (KardB[index] != NULL)
				free(KardB[index]);
		}
		free(KardB);
		KardB = NULL;
	}

	if (Dim != NULL)
	{
		free(Dim);
		Dim = NULL;
	}

	if (KardL != NULL)
	{
		free(KardL);
		KardL = NULL;
	}

	if (Bound != NULL)
	{
		free(Bound);
		Bound = NULL;
	}

	if (kBSol.zf != NULL)
	{
		free(kBSol.zf);
		kBSol.zf = NULL;
	}

	if (kBSol.nx != NULL)
	{
		free(kBSol.nx);
		kBSol.nx = NULL;
	}

	if (kBSol.x != NULL)
	{
		for (k = 0; k < Inst.k; k++)
			free(kBSol.x[k]);
		free(kBSol.x);
		kBSol.x = NULL;
	}

	if (kBSol.id != NULL)
	{
		for (k = 0; k < Inst.k; k++)
			free(kBSol.id[k]);
		free(kBSol.id);
		kBSol.id = NULL;
	}

	return 0;
}

/**
 *  KPSolver: Solve_KPBW_Bound_kBest_LDCF
 */
int KPSolver::Solve_KPBW_Bound_kBest_LDCF(void)
{
	long i;
	long k, kk;
	int status;
	//float t1, t2, dt;

	// Setup data structures
	//t1 = GetTime();
	status = Setup_KPBW_Bound_kBest_LDCF();
	if (status > 0)
		return status;
	//t2 = GetTime();
	//dt = t2 - t1;
	//printf("Setup Computing time: %7.4f \n", dt);

	// Build the k-best solutions
	k = 0;
	kk = 0;
	while (k < Inst.k)
	{
		kBSol.zf[k] = KPBW_Bound_kBest_LDCF(Inst.n, Inst.W, kk, &status);
		if (status > 0)
			return status;
		if (status < 0)
			break;
		BuildSolution_Bound_kBest_LDCF(k, kk);
		if (FeasSol_Bound_kBest_LDCF(k))
			k++;
		kk++;
	}
	kBSol.nk = k;

	// Print the optimal solution
	if (PrintOn > 0)
	{
		printf("Zkp(%d) = %7.2f \n", 0, kBSol.zf[0]);
		printf("Zkp(%d) = %7.2f \n", k - 1, kBSol.zf[k - 1]);
	}
	if (PrintOn > 1)
	{
		printf("Numero soluzioni generate = %5d \n", kk);
		for (kk = 0; kk < k; kk++)
		{
			printf("Zkp(%d) = %7.2f \n", kk, kBSol.zf[kk]);
			if (PrintOn > 2)
			{
				for (i = 0; i < kBSol.nx[kk]; i++)
					printf("id(%d,%d)=%d ", kk, i, kBSol.id[kk][i]);
				printf("\n");
				for (i = 0; i < kBSol.nx[kk]; i++)
					printf(" x(%d,%d)=%d ", kk, i, kBSol.x[kk][i]);
				printf("\n");
			}
		}
	}

	return 0;
}

/**
 *  KPSolver: KPI_KPBW_Bound_kBest_LDCF
 */
int KPSolver::KPI_KPBW_Bound_kBest_LDCF(void)
{
	long j, w;
	size_t index2;
	long Konta;
	long Konta1;
	long Konta2;
	long Konta3;
	long Konta4;
	long Konta5;
	long Konta6;
	long Konta7;
	float AvgK;
	float TotDim;
	float AvgDim;
	float Perc, Perc1, Perc2;

	Konta = 0;
	Konta1 = 0;
	Konta2 = 0;
	Konta3 = 0;
	Konta4 = 0;
	Konta5 = 0;
	Konta6 = 0;
	Konta7 = 0;
	TotDim = 0.0;
	for (j = 0; j <= Inst.n; j++)
		for (w = 0; w <= Inst.W; w++)
		{
			// Setup index
			index2 = MapBW(j, w);

			if (KardL[index2] > 0)
				Konta++;
			if (KardL[index2] > 10)
				Konta3++;
			if (KardL[index2] > 50)
				Konta4++;
			if (KardL[index2] > 100)
				Konta5++;
			if ((j > Inst.n / 2) && (KardL[index2] > 0))
				Konta6++;
			if ((j > Inst.n / 4) && (KardL[index2] > 0))
				Konta7++;
			Konta1 += KardL[index2];

			TotDim += (float)Dim[index2];
		}

	Konta2 = Konta1 / Konta;
	AvgK = (float)Konta1 / Konta;

	// To avoid long integer overflow we split each computation
	// If (Inst.n + 1), (Inst.W + 1), and (Inst.k) are too large we can have an overflow

	//Perc = (float)Konta / ((Inst.n + 1)*(Inst.W + 1)) * 100.;
	Perc = (float)Konta / (Inst.n + 1);
	Perc = Perc / (Inst.W + 1);
	Perc = Perc * (float)100.;

	//Perc1 = (float)Konta1 / ((Inst.n + 1)*(Inst.W + 1)*(Inst.k)) * 100.;
	Perc1 = (float)Konta1 / (Inst.n + 1);
	Perc1 = Perc1 / (Inst.W + 1);
	Perc1 = Perc1 / (Inst.k);
	Perc1 = Perc1 * (float)100.;

	//AvgDim = (float)TotDim / ((Inst.n + 1)*(Inst.W + 1));
	AvgDim = (float)TotDim / (Inst.n + 1);
	AvgDim = AvgDim / (Inst.W + 1);

	//Perc2 = (float)TotDim / ((Inst.n + 1)*(Inst.W + 1)*(Inst.k)) * 100.;
	Perc2 = (float)TotDim / (Inst.n + 1);
	Perc2 = Perc2 / (Inst.W + 1);
	Perc2 = Perc2 / (Inst.k);
	Perc2 = Perc2 * (float)100.;

	// Setup KPI
	KPI.Avg_k = AvgK;
	KPI.Num_States_Used = Konta;
	KPI.Perc_States_Used = Perc;
	KPI.Num_States_Generated = Konta1;
	KPI.Perc_States_Generated = Perc1;
	KPI.Avg_Dim_k = AvgDim;
	KPI.Perc_Dim_k = Perc2;

	if (PrintOn == 1)
	{
		//printf("Numbero of states (j,w) with k>10: %d \n", Konta3);
		//printf("Numbero of states (j,w) with k>50: %d \n", Konta4);
		//printf("Numbero of states (j,w) with k>100: %d \n", Konta5);
		//printf("Numbero of states (j,w) with k>0 and j>n/2: %d \n", Konta6);
		//printf("Numbero of states (j,w) with k>0 and j>n/4: %d \n", Konta7);
		printf("States (j,w) used (active): %d // %d (%f%%) \n", Konta, (Inst.n + 1)*(Inst.W + 1), Perc);
		printf("States (j,w,k) generated: %d // %d  (%f%%) \n", Konta1, (Inst.n + 1)*(Inst.W + 1)*(Inst.k), Perc1);
		//printf("Average k for each active state (j,w): %d \n", Konta2);
		printf("Average k for each active state (j,w): %f \n", AvgK);
		printf("Average Dim(j,w): %f \n", AvgDim);
		printf("Percentage of memory allocated to save the k-best: %f%% \n", Perc2);
	}

	return 0;
}



