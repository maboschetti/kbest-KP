//-----------------------------------------------------------------------------
// Here we must insert our "license boilerplate"
// (we need to investigate what is the most suitable licence for us
//-----------------------------------------------------------------------------
// File: Main-KP-Test.cpp
// Version 1.0
// Last update: 01.12.2022
// Authors: Boschetti M.A., Novellani, S.
//-----------------------------------------------------------------------------
//
#include <time.h>
#include "kBest-KP-Library.h"

int ReadInstance(struct Instance *Inst)
{
	long i;
	FILE *fin;

	fin = fopen("..\\..\\Instances\\input.txt", "r");
	if (fin == NULL)
	{
		printf("File not found!\n");
		exit(1);
	}

	fscanf(fin, "%d", &(Inst->n));
	fscanf(fin, "%d", &(Inst->W));

	Inst->w = new long[Inst->n];
	Inst->p = new long[Inst->n];

	for (i = 0; i < Inst->n; i++)
		fscanf(fin, "%d", &(Inst->w[i]));

	for (i = 0; i < Inst->n; i++)
		fscanf(fin, "%d", &(Inst->p[i]));

	Inst->k = 4;

	fclose(fin);

	return 0;
}

int GenerateInstance(struct Instance *Inst)
{
	long i;
	FILE *fin;
	
	Inst->n = 1000;
	Inst->W = 1000;
	Inst->k = 1000;

	printf("Data: n=%d; W=%d; k=%d \n\n", Inst->n, Inst->W, Inst->k);

	Inst->w = new long[Inst->n];
	Inst->p = new long[Inst->n];
	Inst->pf = new float[Inst->n];
	Inst->b = new long[Inst->n];

	for (i = 0; i < Inst->n; i++)
	{
		Inst->w[i] = long(((double)(rand()) / RAND_MAX) * Inst->W) + 1;
		Inst->w[i] = 100 * long(Inst->w[i] / 100) + 1;
		//Inst->w[i] = 10 * long(Inst->w[i] / 10) + 1;
		if (Inst->w[i] > Inst->W)
			Inst->w[i] = Inst->W;		
		if (Inst->w[i] < 1)
			Inst->w[i] = 1;
	}

	for (i = 0; i < Inst->n; i++)
	{
		Inst->p[i] = long(((double)(rand()) / RAND_MAX) * Inst->W) + 1;
		Inst->pf[i] = (float)Inst->p[i];
		Inst->b[i] = 1;
	}

	return 0;
}


int main(void)
{
	float t1, t2, dt;
	float tt1, tt2, dtt;
	int status;

	KPSolver KPSol;

	//ReadInstance(&(KPSol.Inst));
	GenerateInstance(&(KPSol.Inst));

	// Compute the LP-Relaxation
	KPSol.Allocate_LPKP();
	KPSol.Solve_LPKP();
	KPSol.DeAllocate_LPKP();

	goto salto1;

	// Solve the 1-best KP by the Forward Recursion 
	tt1 = GetTime();
	status = KPSol.Allocate_KPFW();
	if (status > 0)
	{
		printf("Error...\n");
		exit(1);
	}
	t1 = GetTime();
	status = KPSol.Solve_KPFW();
	if (status > 0)
	{
		printf("Error...\n");
		exit(1);
	}
	t2 = GetTime();
	status = KPSol.DeAllocate_KPFW();
	if (status > 0)
	{
		printf("Error...\n");
		exit(1);
	}
	tt2 = GetTime();
	dt = t2 - t1;
	printf("Solve_KPFW: Generation Computing Time: %7.4f \n", dt);
	dtt = tt2 - tt1;
	printf("Solve_KPFW: Total Computing Time: %7.4f \n\n", dtt);

	// Solve the 1-best KP by the Backward Recursion 
	tt1 = GetTime();
	status = KPSol.Allocate_KPBW();
	if (status > 0)
	{
		printf("Error...\n");
		exit(1);
	}
	t1 = GetTime();
	status = KPSol.Solve_KPBW();
	if (status > 0)
	{
		printf("Error...\n");
		exit(1);
	}
	t2 = GetTime();
	status = KPSol.DeAllocate_KPBW();
	if (status > 0)
	{
		printf("Error...\n");
		exit(1);
	}
	tt2 = GetTime();
	dt = t2 - t1;
	printf("Solve_KPBW: Generation Computing Time: %7.4f \n", dt);
	dtt = tt2 - tt1;
	printf("Solve_KPBW: Total Computing Time: %7.4f \n\n", dtt);

	// Solve the k-best KP by the Forward Recursion 
	tt1 = GetTime();
	status = KPSol.Allocate_KPFW_kBest();
	if (status > 0)
	{
		printf("Error...\n");
		exit(1);
	}
	t1 = GetTime();
	status = KPSol.Solve_KPFW_kBest();
	if (status > 0)
	{
		printf("Error...\n");
		exit(1);
	}
	t2 = GetTime();
	status = KPSol.DeAllocate_KPFW_kBest();
	if (status > 0)
	{
		printf("Error...\n");
		exit(1);
	}
	tt2 = GetTime();
	dt = t2 - t1;
	printf("Solve_KPFW_kBest: Generation Computing Time: %7.4f \n", dt);
	dtt = tt2 - tt1;
	printf("Solve_KPFW_kBest: Total Computing Time: %7.4f \n\n", dtt);

salto1:

	// Solve the k-best KP by the Forward Recursion 
	tt1 = GetTime();
	status = KPSol.Allocate_KPFW_kBest_L();
	if (status > 0)
	{
		printf("Error...\n");
		exit(1);
	}
	t1 = GetTime();
	status = KPSol.Solve_KPFW_kBest_L();
	if (status > 0)
	{
		printf("Error...\n");
		exit(1);
	}
	t2 = GetTime();
	status = KPSol.DeAllocate_KPFW_kBest_L();
	if (status > 0)
	{
		printf("Error...\n");
		exit(1);
	}
	tt2 = GetTime();
	dt = t2 - t1;
	printf("Solve_KPFW_kBest_L: Generation Computing Time: %7.4f \n", dt);
	dtt = tt2 - tt1;
	printf("Solve_KPFW_kBest_L: Total Computing Time: %7.4f \n\n", dtt);

	goto salto2;

	// Solve the k-best KP by the Backward Recursion 
	tt1 = GetTime();
	status = KPSol.Allocate_KPBW_kBest();
	if (status > 0)
	{
		printf("Error...\n");
		exit(1);
	}
	t1 = GetTime();
	status = KPSol.Solve_KPBW_kBest();
	if (status > 0)
	{
		printf("Error...\n");
		exit(1);
	}
	t2 = GetTime();
	status = KPSol.DeAllocate_KPBW_kBest();
	if (status > 0)
	{
		printf("Error...\n");
		exit(1);
	}
	tt2 = GetTime();
	dt = t2 - t1;
	printf("Solve_KPBW_kBest: Generation Computing Time: %7.4f \n", dt);
	dtt = tt2 - tt1;
	printf("Solve_KPBW_kBest: Total Computing Time: %7.4f \n\n", dtt);

	// Solve the k-best KP by the Backward Recursion with Malloc and Free 
	tt1 = GetTime();
	status = KPSol.MAllocate_KPBW_kBest();
	if (status > 0)
	{
		printf("Error...\n");
		exit(1);
	}
	t1 = GetTime();
	status = KPSol.Solve_KPBW_kBest();
	if (status > 0)
	{
		printf("Error...\n");
		exit(1);
	}
	t2 = GetTime();
	status = KPSol.DeMAllocate_KPBW_kBest();
	if (status > 0)
	{
		printf("Error...\n");
		exit(1);
	}
	tt2 = GetTime();
	dt = t2 - t1;
	printf("Solve_KPBW_kBest (M): Generation Computing Time: %7.4f \n", dt);
	dtt = tt2 - tt1;
	printf("Solve_KPBW_kBest (M): Total Computing Time: %7.4f \n\n", dtt);

salto2:

	// Solve the k-best KP by the Backward Recursion 
	tt1 = GetTime();
	status =KPSol.Allocate_KPBW_kBest_L();
	if (status > 0)
	{
		printf("Error...\n");
		exit(1);
	}
	t1 = GetTime();
	status = KPSol.Solve_KPBW_kBest_L();
	if (status > 0)
	{
		printf("Error...\n");
		exit(1);
	}
	t2 = GetTime();
	KPSol.KPI_KPBW_kBest_L();
	status = KPSol.DeAllocate_KPBW_kBest_L();
	if (status > 0)
	{
		printf("Error...\n");
		exit(1);
	}
	tt2 = GetTime();
	dt = t2 - t1;
	printf("Solve_KPBW_kBest_L: Generation Computing Time: %7.4f \n", dt);
	dtt = tt2 - tt1;
	printf("Solve_KPBW_kBest_L: Total Computing Time: %7.4f \n\n", dtt);

	// Solve the k-best KP by the Backward Recursion with dynamic reallocation of states k 
	tt1 = GetTime();
	status = KPSol.Allocate_KPBW_kBest_LD();
	if (status > 0)
	{
		printf("Error...\n");
		exit(1);
	}
	t1 = GetTime();
	status = KPSol.Solve_KPBW_kBest_LD();
	if (status > 0)
	{
		printf("Error...\n");
		exit(1);
	}
	t2 = GetTime();
	KPSol.KPI_KPBW_kBest_L();
	status = KPSol.DeAllocate_KPBW_kBest_LD();
	if (status > 0)
	{
		printf("Error...\n");
		exit(1);
	}
	tt2 = GetTime();
	dt = t2 - t1;
	printf("Solve_KPBW_kBest_LD (M): Generation Computing Time: %7.4f \n", dt);
	dtt = tt2 - tt1;
	printf("Solve_KPBW_kBest_LD (M): Total Computing Time: %7.4f \n\n", dtt);

	// Solve the k-best KP by the Backward Recursion with dynamic reallocation of states k
	// This version manage further constraints
	KPSol.Cons.Type = 0;  // Set to 1 if you want to check maximum cardinality constraint
	KPSol.Cons.MaxCard = KPSol.Inst.n / 5;
	tt1 = GetTime();
	status = KPSol.Allocate_KPBW_kBest_LDC();
	if (status > 0)
	{
		printf("Error...\n");
		exit(1);
	}
	t1 = GetTime();
	status = KPSol.Solve_KPBW_kBest_LDC();
	if (status > 0)
	{
		printf("Error...\n");
		exit(1);
	}
	t2 = GetTime();
	KPSol.KPI_KPBW_kBest_L();
	status = KPSol.DeAllocate_KPBW_kBest_LDC();
	if (status > 0)
	{
		printf("Error...\n");
		exit(1);
	}
	tt2 = GetTime();
	dt = t2 - t1;
	printf("Solve_KPBW_kBest_LDC (M): Generation Computing Time: %7.4f \n", dt);
	dtt = tt2 - tt1;
	printf("Solve_KPBW_kBest_LDC (M): Total Computing Time: %7.4f \n\n", dtt);

	// Solve the k-best KP by the Backward Recursion with dynamic reallocation of states k
	// This version manage further constraints and use float for Zkp
	KPSol.Cons.Type = 0;  // Set to 1 if you want to check maximum cardinality constraint
	KPSol.Cons.MaxCard = KPSol.Inst.n / 5;
	tt1 = GetTime();
	status = KPSol.Allocate_KPBW_kBest_LDCF();
	if (status > 0)
	{
		printf("Error...\n");
		exit(1);
	}
	t1 = GetTime();
	status = KPSol.Solve_KPBW_kBest_LDCF();
	if (status > 0)
	{
		printf("Error...\n");
		exit(1);
	}
	t2 = GetTime();
	KPSol.KPI_KPBW_kBest_LDCF();
	status = KPSol.DeAllocate_KPBW_kBest_LDCF();
	if (status > 0)
	{
		printf("Error...\n");
		exit(1);
	}
	tt2 = GetTime();
	dt = t2 - t1;
	printf("Solve_KPBW_kBest_LDCF (M): Generation Computing Time: %7.4f \n", dt);
	dtt = tt2 - tt1;
	printf("Solve_KPBW_kBest_LDCF (M): Total Computing Time: %7.4f \n\n", dtt);

	// Solve the k-best KP Bounded by the Backward Recursion with dynamic reallocation of states k
	// This version manage further constraints and use float for Zkp
	KPSol.Cons.Type = 0;  // Set to 1 if you want to check maximum cardinality constraint
	KPSol.Cons.MaxCard = KPSol.Inst.n / 5;
	tt1 = GetTime();
	status = KPSol.Allocate_KPBW_Bound_kBest_LDCF();
	if (status > 0)
	{
		printf("Error...\n");
		exit(1);
	}
	t1 = GetTime();
	status = KPSol.Solve_KPBW_Bound_kBest_LDCF();
	if (status > 0)
	{
		printf("Error...\n");
		exit(1);
	}
	t2 = GetTime();
	KPSol.KPI_KPBW_Bound_kBest_LDCF();
	status = KPSol.DeAllocate_KPBW_Bound_kBest_LDCF();
	if (status > 0)
	{
		printf("Error...\n");
		exit(1);
	}
	tt2 = GetTime();
	dt = t2 - t1;
	printf("Solve_KPBW_Bound_kBest_LDCF (M): Generation Computing Time: %7.4f \n", dt);
	dtt = tt2 - tt1;
	printf("Solve_KPBW_Bound_kBest_LDCF (M): Total Computing Time: %7.4f \n\n", dtt);

	return 0;
}
 
  