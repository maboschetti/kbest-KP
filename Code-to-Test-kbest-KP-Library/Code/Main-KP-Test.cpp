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

#define READ_INSTANCE //if defined reads instance from files, otherwise a random instance is generated
//#define TESTS //different folder position
#define CSV_OUTPUT //if output in csv file otherwise text file

int ReadInstance(struct Instance* Inst, char* namefile, char* num_k)
{
	long i;
	FILE* fin;

	char new_namefile[100];

#ifndef TESTS
	sprintf(new_namefile, "..\\Instances/");
#else
	sprintf(new_namefile, "..\\..\\Instances/");
#endif
	strcat(new_namefile, namefile);
	strcat(new_namefile, ".txt");
	//

	//fin = fopen("..\\..\\Instances\\input.txt", "r");
	fin = fopen(new_namefile, "r");
	if (fin == NULL)
	{
		printf("File not found!\n");
		exit(1);
	}

	fscanf(fin, "%d", &(Inst->n));
	fscanf(fin, "%d", &(Inst->W));
	//printf("n=%d\n", (Inst->n));
	//printf("W=%d\n", (Inst->W));

	Inst->w = new long[Inst->n];
	Inst->p = new long[Inst->n];

	for (i = 0; i < Inst->n; i++) {
		fscanf(fin, "%d", &(Inst->w[i]));
		fscanf(fin, "%d", &(Inst->p[i]));
		//printf("%d ", (Inst->w[i]));
		//printf("%d\n", (Inst->p[i]));
	}

	//for (i = 0; i < Inst->n; i++)
		//fscanf(fin, "%d", &(Inst->p[i]));

	Inst->k = atoi(num_k);

	fclose(fin);

	Inst->pf = new float[Inst->n];
	Inst->b = new long[Inst->n];

	for (i = 0; i < Inst->n; i++)
	{
		Inst->pf[i] = (float)Inst->p[i];
		Inst->b[i] = 1;
	}


	return 0;
}


int GenerateInstance(struct Instance *Inst)
{
	long i;
	FILE *fin;
	
	Inst->n = 1000;
	Inst->W = 1000;
	Inst->k = 1000;

	Inst->w = new long[Inst->n];
	Inst->p = new long[Inst->n];
	Inst->pf = new float[Inst->n];
	Inst->b = new long[Inst->n];

	for (i = 0; i < Inst->n; i++)
	{
		Inst->w[i] = long(((double)(rand()) / RAND_MAX) * Inst->W) + 1;
		//Inst->w[i] = 100 * long(Inst->w[i] / 100);
		Inst->w[i] = 10 * long(Inst->w[i] / 10) + 1;
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


int main(int argc, char** argv)
{

		float t1, t2, dt;
		float tt1, tt2, dtt;
		int status;

		KPSolver KPSol;

#ifdef READ_INSTANCE

		//argv[1] name input file
		//argv[2] number of k-solutions
		//argv[3] random number for ordering instances
		//argv[4] instance ID
		//argv[5] algorithm type
		//argv[6] run ID (from 0 to 9)

		ReadInstance(&(KPSol.Inst), argv[1], argv[2]); 

		//get information on the instance
		int size_n = trunc(log10(KPSol.Inst.n)) + 1;
		int size_W = trunc(log10(KPSol.Inst.W)) + 1;
		//printf("%d %d %c\n", size_n, size_W, argv[1][size_n + size_W + 3]);
		char buffer[20] = "";
		int ii = size_n + size_W + 3;
		while (argv[1][ii] != '_') {
			buffer[ii - (size_n + size_W + 3)] = argv[1][ii];
			//printf("%c\n", argv[1][ii]);
			ii++;
		}
		int Multiplier = atoi(buffer);
		//printf("M%d\n", Multiplier);
		char buffer2[20] = "";
		int jj = ii + 2;
		while (argv[1][jj] != '_') {
			buffer2[jj - (ii + 2)] = argv[1][jj];
			//printf("%c\n", argv[1][jj]);
			jj++;
		}
		int Divisor = atoi(buffer2);
		//printf("D%d\n", Divisor);


		//create an output file
		char new_namefile[100];
#ifndef TESTS
//		sprintf(new_namefile, "..\\Results/");
		sprintf(new_namefile, "Results/");
#else
		sprintf(new_namefile, "Results/");
#endif
		strcat(new_namefile, "AlgoType_");
		strcat(new_namefile, argv[5]);
		strcat(new_namefile, "_K");
		strcat(new_namefile, argv[2]);
		strcat(new_namefile, "_run");
		strcat(new_namefile, argv[6]);
#ifndef CSV_OUTPUT
		strcat(new_namefile, ".txt");
#else
		strcat(new_namefile, ".csv");
#endif

		FILE* fout;
		//fout = fopen("resultsKPFW.txt", "a");
		fout = fopen(new_namefile, "a");
		if (fout == NULL)
		{
			printf("Error opening the output file!\n");
			exit(1);
		}
#ifndef CSV_OUTPUT
		fprintf(fout, "%s\t%d\t%d\t%d\t%d\t%d\t\%d\t%d\t%d\t%d", argv[1], atoi(argv[5]), atoi(argv[4]), atoi(argv[3]), atoi(argv[6]), KPSol.Inst.n, KPSol.Inst.W, KPSol.Inst.k, Multiplier, Divisor);
#else
		fprintf(fout, "%s;%d;%d;%d;%d;%d;%d;%d;%d;%d", argv[1], atoi(argv[5]), atoi(argv[4]), atoi(argv[3]), atoi(argv[6]), KPSol.Inst.n, KPSol.Inst.W, KPSol.Inst.k, Multiplier, Divisor);
#endif
		//
#else
		GenerateInstance(&(KPSol.Inst));
#endif

	//// Allocate the data structure to save the solutions
	//KPSol.Allocate_Sol();

		
	switch (atoi(argv[5]))
	{
	case 0:
		// Compute the LP-Relaxation
		KPSol.Allocate_LPKP();
		KPSol.Solve_LPKP();
		KPSol.DeAllocate_LPKP();

		break;

	case 1:
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
#ifndef CSV_OUTPUT
		fprintf(fout, "\t%d", KPSol.OptSol.z);
#else
		fprintf(fout, ";%d", KPSol.OptSol.z);
#endif
		status = KPSol.DeAllocate_KPFW();
		if (status > 0)
		{
			printf("Error...\n");
			exit(1);
		}
		tt2 = GetTime();
		dt = t2 - t1;
		//printf("Solve_KPFW: Generation Computing Time: %7.6f \n", dt);
#ifndef CSV_OUTPUT
		fprintf(fout, "\t%7.6f", dt);
#else
		fprintf(fout, ";%7.6f", dt);
#endif
		dtt = tt2 - tt1;
		//printf("Solve_KPFW: Total Computing Time: %7.6f \n\n", dtt);
#ifndef CSV_OUTPUT
		fprintf(fout, "\t%7.6f", dtt);
#else
		fprintf(fout, ";%7.6f", dtt);
#endif

		break;

	case 2:
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
#ifndef CSV_OUTPUT
		fprintf(fout, "\t%d", KPSol.OptSol.z);
#else
		fprintf(fout, ";%d", KPSol.OptSol.z);
#endif
		status = KPSol.DeAllocate_KPBW();
		if (status > 0)
		{
			printf("Error...\n");
			exit(1);
		}
		tt2 = GetTime();
		dt = t2 - t1;
		//printf("Solve_KPBW: Generation Computing Time: %7.6f \n", dt);
#ifndef CSV_OUTPUT
		fprintf(fout, "\t%7.6f", dt);
#else
		fprintf(fout, ";%7.6f", dt);
#endif
		dtt = tt2 - tt1;
		//printf("Solve_KPBW: Total Computing Time: %7.6f \n\n", dtt);
#ifndef CSV_OUTPUT
		fprintf(fout, "\t%7.6f", dtt);
#else
		fprintf(fout, ";%7.6f", dtt);
#endif
		
		break;

	case 3:
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
#ifndef CSV_OUTPUT
		fprintf(fout, "\t%d", KPSol.kBSol.z[0]);
		fprintf(fout, "\t%d", KPSol.kBSol.z[KPSol.kBSol.nk - 1]);
#else
		fprintf(fout, ";%d", KPSol.kBSol.z[0]);
		fprintf(fout, ";%d", KPSol.kBSol.z[KPSol.kBSol.nk - 1]);
#endif
		status = KPSol.DeAllocate_KPFW_kBest();
		if (status > 0)
		{
			printf("Error...\n");
			exit(1);
		}
		tt2 = GetTime();
		dt = t2 - t1;
		//printf("Solve_KPFW_kBest: Generation Computing Time: %7.6f \n", dt);
#ifndef CSV_OUTPUT
		fprintf(fout, "\t%7.6f", dt);
#else
		fprintf(fout, ";%7.6f", dt);
#endif
		dtt = tt2 - tt1;
		//printf("Solve_KPFW_kBest: Total Computing Time: %7.6f \n\n", dtt);
#ifndef CSV_OUTPUT
		fprintf(fout, "\t%7.6f", dtt);
#else
		fprintf(fout, ";%7.6f", dtt);
#endif
		if (KPSol.kBSol.z[0] < KPSol.kBSol.z[KPSol.kBSol.nk - 1]) {
			//fprintf(fout, "\tError!");
			printf("Error...\n");
			exit(1);
		}

		break; 

	case 4:
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
#ifndef CSV_OUTPUT
		fprintf(fout, "\t%d", KPSol.kBSol.z[0]);
		fprintf(fout, "\t%d", KPSol.kBSol.z[KPSol.kBSol.nk - 1]);
#else
		fprintf(fout, ";%d", KPSol.kBSol.z[0]);
		fprintf(fout, ";%d", KPSol.kBSol.z[KPSol.kBSol.nk - 1]);
#endif
		status = KPSol.DeAllocate_KPFW_kBest_L();
		if (status > 0)
		{
			printf("Error...\n");
			exit(1);
		}
		tt2 = GetTime();
		dt = t2 - t1;
		//printf("Solve_KPFW_kBest_L: Generation Computing Time: %7.6f \n", dt);
#ifndef CSV_OUTPUT
		fprintf(fout, "\t%7.6f", dt);
#else
		fprintf(fout, ";%7.6f", dt);
#endif
		dtt = tt2 - tt1;
		//printf("Solve_KPFW_kBest_L: Total Computing Time: %7.6f \n\n", dtt);
#ifndef CSV_OUTPUT
		fprintf(fout, "\t%7.6f", dtt);
#else
		fprintf(fout, ";%7.6f", dtt);
#endif
		/*if (KPSol.kBSol.z[0] < KPSol.kBSol.z[KPSol.kBSol.nk - 1]) {
			//fprintf(fout, "\tError!");
			printf("Error...\n");
			exit(1);
		}*/

		break;

	case 5:
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
#ifndef CSV_OUTPUT
		fprintf(fout, "\t%d", KPSol.kBSol.z[0]);
		fprintf(fout, "\t%d", KPSol.kBSol.z[KPSol.kBSol.nk - 1]);
#else
		fprintf(fout, ";%d", KPSol.kBSol.z[0]);
		fprintf(fout, ";%d", KPSol.kBSol.z[KPSol.kBSol.nk - 1]);
#endif
		status = KPSol.DeAllocate_KPBW_kBest();
		if (status > 0)
		{
			printf("Error...\n");
			exit(1);
		}
		tt2 = GetTime();
		dt = t2 - t1;
		//printf("Solve_KPBW_kBest: Generation Computing Time: %7.6f \n", dt);
#ifndef CSV_OUTPUT
		fprintf(fout, "\t%7.6f", dt);
#else
		fprintf(fout, ";%7.6f", dt);
#endif
		dtt = tt2 - tt1;
		//printf("Solve_KPBW_kBest: Total Computing Time: %7.6f \n\n", dtt);
#ifndef CSV_OUTPUT
		fprintf(fout, "\t%7.6f", dtt);
#else
		fprintf(fout, ";%7.6f", dtt);
#endif
		/*if (KPSol.kBSol.z[0] < KPSol.kBSol.z[KPSol.kBSol.nk - 1]) {
			//fprintf(fout, "\tError!");
			printf("Error...\n");
			exit(1);
		}*/

		break;

	case 6:
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
#ifndef CSV_OUTPUT
		fprintf(fout, "\t%d", KPSol.kBSol.z[0]);
		fprintf(fout, "\t%d", KPSol.kBSol.z[KPSol.kBSol.nk - 1]);
#else
		fprintf(fout, ";%d", KPSol.kBSol.z[0]);
		fprintf(fout, ";%d", KPSol.kBSol.z[KPSol.kBSol.nk - 1]);
#endif
		status = KPSol.DeMAllocate_KPBW_kBest();
		if (status > 0)
		{
			printf("Error...\n");
			exit(1);
		}
		tt2 = GetTime();
		dt = t2 - t1;
		//printf("Solve_KPBW_kBest (M): Generation Computing Time: %7.6f \n", dt);
#ifndef CSV_OUTPUT
		fprintf(fout, "\t%7.6f", dt);
#else
		fprintf(fout, ";%7.6f", dt);
#endif
		dtt = tt2 - tt1;
		//printf("Solve_KPBW_kBest (M): Total Computing Time: %7.6f \n\n", dtt);
#ifndef CSV_OUTPUT
		fprintf(fout, "\t%7.6f", dt);
#else
		fprintf(fout, ";%7.6f", dt);
#endif
		/*if (KPSol.kBSol.z[0] < KPSol.kBSol.z[KPSol.kBSol.nk - 1]) {
			//fprintf(fout, "\tError!");
			printf("Error...\n");
			exit(1);
		}*/

		break;

	case 7:
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
#ifndef CSV_OUTPUT
		fprintf(fout, "\t%d", KPSol.kBSol.z[0]);
		fprintf(fout, "\t%d", KPSol.kBSol.z[KPSol.kBSol.nk - 1]);
#else
#endif
		KPSol.KPI_KPBW_kBest_L();
#ifndef CSV_OUTPUT
		fprintf(fout, "\t%7.4f", KPSol.KPI.Avg_k);
		fprintf(fout, "\t%d", KPSol.KPI.Num_States_Used);
		fprintf(fout, "\t%7.4f", KPSol.KPI.Perc_States_Used);
		fprintf(fout, "\t%d", KPSol.KPI.Num_States_Generated);
		fprintf(fout, "\t%7.4f", KPSol.KPI.Perc_States_Generated);
		fprintf(fout, "\t%7.4f", KPSol.KPI.Avg_Dim_k);
		fprintf(fout, "\t%7.4f", KPSol.KPI.Perc_Dim_k);
#else
		fprintf(fout, ";%7.4f", KPSol.KPI.Avg_k);
		fprintf(fout, ";%d", KPSol.KPI.Num_States_Used);
		fprintf(fout, ";%7.4f", KPSol.KPI.Perc_States_Used);
		fprintf(fout, ";%d", KPSol.KPI.Num_States_Generated);
		fprintf(fout, ";%7.4f", KPSol.KPI.Perc_States_Generated);
		fprintf(fout, ";%7.4f", KPSol.KPI.Avg_Dim_k);
		fprintf(fout, ";%7.4f", KPSol.KPI.Perc_Dim_k);
#endif
		status = KPSol.DeAllocate_KPBW_kBest_L();
		if (status > 0)
		{
			printf("Error...\n");
			exit(1);
		}
		tt2 = GetTime();
		dt = t2 - t1;
		//printf("Solve_KPBW_kBest_L: Generation Computing Time: %7.6f \n", dt);
#ifndef CSV_OUTPUT
		fprintf(fout, "\t%7.6f", dt);
#else
		fprintf(fout, ";%7.6f", dt);
#endif
		dtt = tt2 - tt1;
		//printf("Solve_KPBW_kBest_L: Total Computing Time: %7.6f \n\n", dtt);
#ifndef CSV_OUTPUT
		fprintf(fout, "\t%7.6f", dtt);
#else
		fprintf(fout, ";%7.6f", dtt);
#endif
		/*if (KPSol.kBSol.z[0] < KPSol.kBSol.z[KPSol.kBSol.nk - 1]) {
			//fprintf(fout, "\tError!");
			printf("Error...\n");
			exit(1);
		}*/

		break;

	case 8:
		// Solve the k-best KP by the Backward Recursion with dynamic reallocation of states k 
		//KPSol.Inst.k = KPSol.Inst.k * 10;
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
#ifndef CSV_OUTPUT
		fprintf(fout, "\t%d", KPSol.kBSol.z[0]);
		fprintf(fout, "\t%d", KPSol.kBSol.z[KPSol.kBSol.nk - 1]);
#else
		fprintf(fout, ";%d", KPSol.kBSol.z[0]);
		fprintf(fout, ";%d", KPSol.kBSol.z[KPSol.kBSol.nk - 1]);
#endif
		KPSol.KPI_KPBW_kBest_L();
#ifndef CSV_OUTPUT
		fprintf(fout, "\t%7.4f", KPSol.KPI.Avg_k);
		fprintf(fout, "\t%d", KPSol.KPI.Num_States_Used);
		fprintf(fout, "\t%7.4f", KPSol.KPI.Perc_States_Used);
		fprintf(fout, "\t%d", KPSol.KPI.Num_States_Generated);
		fprintf(fout, "\t%7.4f", KPSol.KPI.Perc_States_Generated);
		fprintf(fout, "\t%7.4f", KPSol.KPI.Avg_Dim_k);
		fprintf(fout, "\t%7.4f", KPSol.KPI.Perc_Dim_k);
#else
		fprintf(fout, ";%7.4f", KPSol.KPI.Avg_k);
		fprintf(fout, ";%d", KPSol.KPI.Num_States_Used);
		fprintf(fout, ";%7.4f", KPSol.KPI.Perc_States_Used);
		fprintf(fout, ";%d", KPSol.KPI.Num_States_Generated);
		fprintf(fout, ";%7.4f", KPSol.KPI.Perc_States_Generated);
		fprintf(fout, ";%7.4f", KPSol.KPI.Avg_Dim_k);
		fprintf(fout, ";%7.4f", KPSol.KPI.Perc_Dim_k);
#endif
		status = KPSol.DeAllocate_KPBW_kBest_LD();
		if (status > 0)
		{
			printf("Error...\n");
			exit(1);
		}
		tt2 = GetTime();
		dt = t2 - t1;
		//printf("Solve_KPBW_kBest_LD (M): Generation Computing Time: %7.6f \n", dt);
#ifndef CSV_OUTPUT
		fprintf(fout, "\t%7.6f", dt);
#else
		fprintf(fout, ";%7.6f", dt);
#endif
		dtt = tt2 - tt1;
		//printf("Solve_KPBW_kBest_LD (M): Total Computing Time: %7.6f \n\n", dtt);
#ifndef CSV_OUTPUT
		fprintf(fout, "\t%7.6f", dtt);
#else
		fprintf(fout, ";%7.6f", dtt);
#endif

		/*if (KPSol.kBSol.z[0] < KPSol.kBSol.z[KPSol.kBSol.nk - 1]) {
			//fprintf(fout, "\tError!");
			printf("Error...\n");
			exit(1);
		}*/

		break;
	
	case 9:
		// Solve the k-best KP by the Backward Recursion with dynamic reallocation of states k
		// This version manage further constraints
		KPSol.Cons.Type = 0;  // Set to 1 is we want check constraints
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
#ifndef CSV_OUTPUT
		fprintf(fout, "\t%d", KPSol.kBSol.z[0]);
		fprintf(fout, "\t%d", KPSol.kBSol.z[KPSol.kBSol.nk - 1]);
#else
		fprintf(fout, ";%d", KPSol.kBSol.z[0]);
		fprintf(fout, ";%d", KPSol.kBSol.z[KPSol.kBSol.nk - 1]);
#endif
		KPSol.KPI_KPBW_kBest_L();
#ifndef CSV_OUTPUT
		fprintf(fout, "\t%7.4f", KPSol.KPI.Avg_k);
		fprintf(fout, "\t%d", KPSol.KPI.Num_States_Used);
		fprintf(fout, "\t%7.4f", KPSol.KPI.Perc_States_Used);
		fprintf(fout, "\t%d", KPSol.KPI.Num_States_Generated);
		fprintf(fout, "\t%7.4f", KPSol.KPI.Perc_States_Generated);
		fprintf(fout, "\t%7.4f", KPSol.KPI.Avg_Dim_k);
		fprintf(fout, "\t%7.4f", KPSol.KPI.Perc_Dim_k);
#else
		fprintf(fout, ";%d", KPSol.kBSol.z[0]);
		fprintf(fout, ";%d", KPSol.kBSol.z[KPSol.kBSol.nk - 1]);
#endif
		status = KPSol.DeAllocate_KPBW_kBest_LDC();
		if (status > 0)
		{
			printf("Error...\n");
			exit(1);
		}
		tt2 = GetTime();
		dt = t2 - t1;
		//printf("Solve_KPBW_kBest_LDC (M): Generation Computing Time: %7.6f \n", dt);
#ifndef CSV_OUTPUT
		fprintf(fout, "\t%7.6f", dt);
#else
		fprintf(fout, ";%7.6f", dt);
#endif
		dtt = tt2 - tt1;
		//printf("Solve_KPBW_kBest_LDC (M): Total Computing Time: %7.6f \n\n", dtt);
#ifndef CSV_OUTPUT
		fprintf(fout, "\t%7.6f", dtt);
#else
		fprintf(fout, ";%7.6f", dtt);
#endif
		/*if (KPSol.kBSol.z[0] < KPSol.kBSol.z[KPSol.kBSol.nk - 1]) {
			//fprintf(fout, "\tError!");
			printf("Error...\n");
			exit(1);
		}*/

		break;

	case 10:
		// Solve the k-best KP by the Backward Recursion with dynamic reallocation of states k
		// This version manage further constraints and use float for Zkp
		KPSol.Cons.Type = 0;  // Set to 1 is we want check constraints
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
#ifndef CSV_OUTPUT
		fprintf(fout, "\t%7.4f", KPSol.kBSol.zf[0]);
		fprintf(fout, "\t%7.4f", KPSol.kBSol.zf[KPSol.kBSol.nk - 1]);
#else
		fprintf(fout, ";%7.4f", KPSol.kBSol.zf[0]);
		fprintf(fout, ";%7.4f", KPSol.kBSol.zf[KPSol.kBSol.nk - 1]);
#endif
		KPSol.KPI_KPBW_kBest_LDCF();
#ifndef CSV_OUTPUT
		fprintf(fout, "\t%7.4f", KPSol.KPI.Avg_k);
		fprintf(fout, "\t%d", KPSol.KPI.Num_States_Used);
		fprintf(fout, "\t%7.4f", KPSol.KPI.Perc_States_Used);
		fprintf(fout, "\t%d", KPSol.KPI.Num_States_Generated);
		fprintf(fout, "\t%7.4f", KPSol.KPI.Perc_States_Generated);
		fprintf(fout, "\t%7.4f", KPSol.KPI.Avg_Dim_k);
		fprintf(fout, "\t%7.4f", KPSol.KPI.Perc_Dim_k);
#else
		fprintf(fout, ";%7.4f", KPSol.KPI.Avg_k);
		fprintf(fout, ";%d", KPSol.KPI.Num_States_Used);
		fprintf(fout, ";%7.4f", KPSol.KPI.Perc_States_Used);
		fprintf(fout, ";%d", KPSol.KPI.Num_States_Generated);
		fprintf(fout, ";%7.4f", KPSol.KPI.Perc_States_Generated);
		fprintf(fout, ";%7.4f", KPSol.KPI.Avg_Dim_k);
		fprintf(fout, ";%7.4f", KPSol.KPI.Perc_Dim_k);
#endif
		status = KPSol.DeAllocate_KPBW_kBest_LDCF();
		if (status > 0)
		{
			printf("Error...\n");
			exit(1);
		}
		tt2 = GetTime();
		dt = t2 - t1;
		//printf("Solve_KPBW_kBest_LDCF (M): Generation Computing Time: %7.6f \n", dt);
#ifndef CSV_OUTPUT
		fprintf(fout, "\t%7.6f", dt);
#else
		fprintf(fout, ";%7.6f", dt);
#endif
		dtt = tt2 - tt1;
		//printf("Solve_KPBW_kBest_LDCF (M): Total Computing Time: %7.6f \n\n", dtt);
#ifndef CSV_OUTPUT
		fprintf(fout, "\t%7.6f", dtt);
#else
		fprintf(fout, ";%7.6f", dtt);
#endif
		/*if (KPSol.kBSol.z[0]  < KPSol.kBSol.z[KPSol.kBSol.nk - 1]) {
			//fprintf(fout, "\tError!");
			printf("Error...\n");
			exit(1);
		}*/

		break;

	case 11:
		// Solve the k-best KP Bounded by the Backward Recursion with dynamic reallocation of states k
		// This version manage further constraints and use float for Zkp
		KPSol.Cons.Type = 0;  // Set to 1 is we want check constraints
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
#ifndef CSV_OUTPUT
		fprintf(fout, "\t%7.4f", KPSol.kBSol.zf[0]);
		fprintf(fout, "\t%7.4f", KPSol.kBSol.zf[KPSol.kBSol.nk - 1]);
#else
		fprintf(fout, ";%7.4f", KPSol.kBSol.zf[0]);
		fprintf(fout, ";%7.4f", KPSol.kBSol.zf[KPSol.kBSol.nk - 1]);
#endif
		KPSol.KPI_KPBW_Bound_kBest_LDCF();
#ifndef CSV_OUTPUT
		fprintf(fout, "\t%7.4f", KPSol.KPI.Avg_k);
		fprintf(fout, "\t%d", KPSol.KPI.Num_States_Used);
		fprintf(fout, "\t%7.4f", KPSol.KPI.Perc_States_Used);
		fprintf(fout, "\t%d", KPSol.KPI.Num_States_Generated);
		fprintf(fout, "\t%7.4f", KPSol.KPI.Perc_States_Generated);
		fprintf(fout, "\t%7.4f", KPSol.KPI.Avg_Dim_k);
		fprintf(fout, "\t%7.4f", KPSol.KPI.Perc_Dim_k);
#else
		fprintf(fout, ";%7.4f", KPSol.KPI.Avg_k);
		fprintf(fout, ";%d", KPSol.KPI.Num_States_Used);
		fprintf(fout, ";%7.4f", KPSol.KPI.Perc_States_Used);
		fprintf(fout, ";%d", KPSol.KPI.Num_States_Generated);
		fprintf(fout, ";%7.4f", KPSol.KPI.Perc_States_Generated);
		fprintf(fout, ";%7.4f", KPSol.KPI.Avg_Dim_k);
		fprintf(fout, ";%7.4f", KPSol.KPI.Perc_Dim_k);
#endif
		status = KPSol.DeAllocate_KPBW_Bound_kBest_LDCF();
		if (status > 0)
		{
			printf("Error...\n");
			exit(1);
		}
		tt2 = GetTime();
		dt = t2 - t1;
		//printf("Solve_KPBW_Bound_kBest_LDCF (M): Generation Computing Time: %7.4f \n", dt);
#ifndef CSV_OUTPUT
		fprintf(fout, "\t%7.4f", dt);
#else
		fprintf(fout, ";%7.4f", dt);
#endif
		dtt = tt2 - tt1;
		//printf("Solve_KPBW_Bound_kBest_LDCF (M): Total Computing Time: %7.4f \n\n", dtt);
#ifndef CSV_OUTPUT
		fprintf(fout, "\t%7.4f", dtt);
#else
		fprintf(fout, ";%7.4f", dtt);
#endif
		/*if (KPSol.kBSol.z[0] < KPSol.kBSol.z[KPSol.kBSol.nk - 1]) {
			//fprintf(fout, "\tError!");
			printf("Error...\n");
			exit(1);
		}*/

		break;

	default:
		break;
	}

	fprintf(fout, "\n");
	//fprintf(fout, "%s", "\n");
	fclose(fout);



	return 0;
}
 
  