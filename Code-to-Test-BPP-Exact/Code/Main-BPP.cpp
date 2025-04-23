//-----------------------------------------------------------------------------
// Here we must insert our "license boilerplate"
// (... to be defined)
//-----------------------------------------------------------------------------
// File: Main-BPP.cpp
// Version 1.0
// Last update: 20.04.2025
// Authors: Boschetti M.A., Novellani, S.
//-----------------------------------------------------------------------------
//
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "Source-Exact/Subgradient.h"

// Prototypes
long EliminatePath(char *fname);
long ReadRow(FILE *fin, char *row);
void Delay(long milliseconds);
int GenInstance(long irun, long type, long nitems, long capbin, long kcons, Subgradient *SubLR);
int ReadInstance(char *filename, Subgradient *SubLR);

/**
 *	Main
 */
int main(int argc, char *argv[])
{
	int status;        // Status used to handle exceptions
	Subgradient SubLR; // Solver

	FILE *fpro;  // Input file
	FILE *fout;  // Output file
	char fname[255];
	char fname1[255];
	char row[255];
	long i, j, npro;
	long t1, t2;
	long iter;
	long BolOpt;
	double dt, dtlr, dtot;
	double dtlp;
	double *xt;
	double Mega = 1000000.;
	double Giga = 1000000000.;
	double Zlp, Zmip, Zopt;
	double ZGreedy, tGreedy;
	double GapGreedy, GapLagrHeu, GapLR0, GapLR, GapOpt;
	long BolBPP = 0;  // If BolBPP=1, solve the classical model for the BPP
	long BolMIP = 1;  
	long BolLR = 1;
	long BolLP = 1;  // If BolLP=1, solve the LP-relaxation by the Incremental Core approach
	long BolEsatto = 1;
	long BolGreedy = 1;
	long BolLagrHeu = 1;
	long NRip = 1;     // Number of times you run the same instances
	long Tdelay = 50;  // Time in milliseconds between two run
	int ProType=0;
	char filename[255];
	char filename1[255];

	double Best_dt;
	double Best_Zlb;
	double Best_Zub;
	double Best_Gap;
	double Best_Alfa;
	long Best_Nodes;
	long Best_Iter;
	long Best_NCore;

	FILE *ftab2;

	double TimeDA;
	double TimeMip;
	double TimeTot;

	long nrun;
	long type;
	long nitems;
	long capbin;
	long kcons;
	long irun;

	// Setup parameters
	SubLR.FixPen = 100000.0;
	SubLR.MaxM = 20000.0;

	fpro = fopen("InputTests.pro", "r");
	//fpro = fopen("InputTests2.pro", "r");
	//fpro = fopen("InputTestsF.pro", "r");
	//fpro = fopen("InputTestsF2.pro", "r");
	if (fpro == NULL)
		exit(1);
	fout = fopen("Exact-Tab.out", "w");
	ftab2 = fopen("Exact-Tab.csv", "w");

	fprintf(fout, "\n                               ");
	fprintf(fout, "      Size        ");
	fprintf(fout, "       DFF      ");
	if (BolBPP)
	{
		fprintf(fout, "                       Cplex                ");
	}
	fprintf(fout, "                            LRStandard                            ");
	if (BolLP)
	{
		fprintf(fout, "       LP-relaxation             ");
	}
	fprintf(fout, "                         Exact                                  ");
	fprintf(fout, "\n");

	fprintf(fout, "Istances                 ");
	fprintf(fout, "    n         W         k   ");
	fprintf(fout, "  ZLBDFF  ");
	if (BolBPP)
	{
		fprintf(fout, "     Zlp        Zopt     Gap       Nodes      Time  ");
	}
	fprintf(fout, "   ZLBlagr    ZUBlagr      Alfa       Iter      Time     nCore   ");
	if (BolLP)
	{
		fprintf(fout, "    Zlp        Time     nCore   ");
	}
	fprintf(fout, "    Zlp         Zopt       Gap      Nodes      Time     nCore   Iter   T.Time ");
	fprintf(fout, "  Opt  ");
	fprintf(fout, "\n");
	fflush(fout);

	fprintf(ftab2, "Istances ; n ; W ; k ; ZLBDFF; z_LBlagr ; z_UBlagr ; Alfa; Iter; Time ; nCore; z_LP; z_opt; Gap ; Nodes; Time ; nCore ; Iter; T.Time; Opt \n");
	fflush(ftab2);

	// Read Parameters

	// Elimininate header
	ReadRow(fpro, row);

	// Read the model to use: Set Partitioning or Set Covering
	fscanf(fpro, "%d", &SubLR.Param.Model);
	// Eliminate description
	ReadRow(fpro, row);

	// Read the number of instrances
	fscanf(fpro, "%f", &SubLR.Param.TimeLimit);
	// Eliminate description
	ReadRow(fpro, row);

	// Read the maximum number of columns
	fscanf(fpro, "%d", &SubLR.Param.MaxNumColumns);
	// Eliminate description
	ReadRow(fpro, row);

	// Read the initial maximum number of columns
	fscanf(fpro, "%d", &SubLR.Param.MaxNumColMIP);
	// Eliminate description
	ReadRow(fpro, row);

	// Read the initial step for the subgradient
	fscanf(fpro, "%lf", &SubLR.Param.alpha);
	// Eliminate description
	ReadRow(fpro, row);

	// Read the maximum number of iterations without improvement for the subgradient
	fscanf(fpro, "%ld", &SubLR.Param.MaxItMinGap);
	// Eliminate description
	ReadRow(fpro, row);

	// Read the MipGap for the subgradient
	fscanf(fpro, "%lf", &SubLR.Param.MinGap);
	// Eliminate description
	ReadRow(fpro, row);

	// Read PerK: k = PercK * n (where n is the number of items)
	fscanf(fpro, "%lf", &SubLR.Param.PercK);
	// Eliminate description
	ReadRow(fpro, row);

	// Read PercDelta: Delta = PercDelta * LB (that is the maximum reduced cost to include a column in the MIP)
	fscanf(fpro, "%lf", &SubLR.Param.PercDelta);
	// Eliminate description
	ReadRow(fpro, row);

	// Read Problem Type: 0=generate (using "Pro_Type  Num_items  Bin_Size  Constraint"; 1=read from "filename"
	fscanf(fpro, "%d", &ProType);
	// Eliminate description
	ReadRow(fpro, row);

	// Read Instances
	// Elimininate header
	ReadRow(fpro, row);

	// Read the number of instrances
	fscanf(fpro, "%d", &npro);

	// Elimininate header
	ReadRow(fpro, row);
	ReadRow(fpro, row);

	for (i = 0; i < npro; i++)
	{
		if (ProType == 0)
		{
			// Read Data
			fscanf(fpro, "%d %d %d %d %d", &nrun, &type, &nitems, &capbin, &kcons);
		}
		else
		{
			// Read Data
			fscanf(fpro, "%d", &nrun);
			fscanf(fpro, "%s", filename);
		}

		for (irun = 0; irun < nrun; irun++)
		{
			if (ProType == 0)
			{
				// Generate the instance
				GenInstance(irun, type, nitems, capbin, kcons, &SubLR);

				// Define name
				strcpy(fname, "BPP-");
				itoa(i, fname1, 10);
				strcat(fname, fname1);
				strcat(fname, "-");
				itoa(irun, fname1, 10);
				strcat(fname, fname1);
			}
			else
			{
				// Generate filename
				strcpy(filename1, filename);
				strcat(filename1, "_");
				itoa(irun, fname1, 10);
				strcat(filename1, fname1);
				strcat(filename1, ".txt");

				// Generate the instance
				status = ReadInstance(filename1, &SubLR);
				if (status)
				{
					printf("ERROR... reading file...\n");
					break;
				}
				strcpy(fname, filename);
				j = EliminatePath(fname);
				strcpy(fname, fname + j);
				strcat(fname, "_");
				itoa(irun, fname1, 10);
				strcat(fname, fname1);
			}

			// Allocate the data structures
			status = SubLR.AllocateDataStructure();
			if (status)
				exit(status);

			// Print name instance
			fprintf(fout, "%-20s", fname);

			// Instance sizes
			fprintf(fout, "%10d", SubLR.Inst.n);
			fprintf(fout, "%10d", SubLR.Inst.W);
			fprintf(fout, "%10d", SubLR.Inst.k);
			fflush(fout);

			// Setup of the variables needed to generate the tables
			Zlp = 0.0;
			Zmip = +Inf;
			GapLR = 100.;
			GapOpt = 100.;

			// Greedy
			tGreedy = 0.0;
			ZGreedy = +Inf;
			SubLR.ZGreedy = +Inf;
			// Include a heuristic algorithm 
			// ...

			// Reductions
			//SubLR.Reductions();

			// Lower Bound based on Dual Feasible Functions
			SubLR.LBDFF();
			//SubLR.LBDFFC1();
			fprintf(fout, "%11.2lf", SubLR.ZLBDFF);

			// Solve the classical MIP model for the BPP 
			SubLR.ZLBlp = 0.0;
			SubLR.ZUBlp = 0.0;
			SubLR.Gap = 0.0;
			SubLR.Nodes = 0;
			xt = new double[SubLR.Param.MaxNumColumns];

			// Solve the the classical MIP model if (BolBPP > 0)
			if (BolBPP)
			{
				t1 = clock();
				//SubLR.Solve_Model(1, xt, &iter); 
				t2 = clock();
				dt = (double)(t2 - t1) / CLOCKS_PER_SEC;

				fprintf(fout, "%11.2lf", SubLR.ZLBlp);
				fprintf(fout, "%11.2lf", SubLR.ZUBlp);
				fprintf(fout, "%10.4lf", SubLR.Gap);
				fprintf(fout, "%10d", SubLR.Nodes);
				fprintf(fout, "%10.2lf", dt);
				fflush(fout);
			}

			// Solve the Dual Ascent based on a Lagrangian Column Generation without Lagrangian Heuristic
			BolLagrHeu = 0;
			SubLR.ZLBlagr = 0.0;
			SubLR.ZUBlagr = 10000.0; // +Inf;
			SubLR.ZGreedy = 10000.0; // +Inf;
			SubLR.alpha = 0.0;
			SubLR.Iter = 0;
			SubLR.nCore = 0;
			BolOpt = 0;
			Best_dt = +Inf;
			for (j = 0; j < NRip; j++)
			{
				t1 = clock();
				if (BolLR)
				{
					if (SubLR.Param.Model == 0)
						SubLR.Dual_Ascent_SetPar(BolLagrHeu, &BolOpt);
					else
						SubLR.Dual_Ascent_SetCov(BolLagrHeu);
				}
				t2 = clock();
				dt = (double)(t2 - t1) / CLOCKS_PER_SEC;
				if (dt < Best_dt)
				{
					Best_dt = dt;
					Best_Zlb = SubLR.ZLBlagr;
					Best_Zub = SubLR.ZUBlagr;
					Best_Alfa = SubLR.alpha;
					Best_Iter = SubLR.Iter;
					Best_NCore = SubLR.nCore;
				}
				Delay(Tdelay);
			}
			dt = Best_dt;
			TimeDA = Best_dt;
			SubLR.TimeDA = Best_dt;
			SubLR.ZLBlagr = Best_Zlb;
			SubLR.ZUBlagr = Best_Zub;
			SubLR.alpha = Best_Alfa;
			SubLR.Iter = Best_Iter;
			SubLR.nCore = Best_NCore;
			SubLR.ZLBDA = Best_Zlb;
			SubLR.ZUBDA = Best_Zub;
			SubLR.alphaDA = Best_Alfa;
			SubLR.IterDA = Best_Iter;
			SubLR.nCoreDA = Best_NCore;

			dtlr = dt;
			//GapLR = 100. * abs(Zlp - SubLR.ZLBlagr)/Zlp;
			GapLR0 = 100. * (Zlp - SubLR.ZLBlagr) / Zlp;

			//if (SubLR.ZUBlagr > ZGreedy)
			//	SubLR.ZUBlagr = ZGreedy;

			fprintf(fout, "%11.2lf", SubLR.ZLBlagr);
			if (BolOpt)
				fprintf(fout, "*");
			else
				fprintf(fout, " ");
			fprintf(fout, "%11.2lf ", SubLR.ZUBlagr);
			fprintf(fout, "%10.4lf", SubLR.alpha);
			fprintf(fout, "%10d", SubLR.Iter);
			fprintf(fout, "%10.2lf", dtlr);
			fprintf(fout, "%10d", SubLR.nCore);
			fflush(fout);

			// Solve LP by Incremental Core to obtain optimal dual solution
			dtlp = 0;
			if (BolLP)
			{
				t1 = clock();
				SubLR.Solve_Incremental_Core(0, BolGreedy, NULL, &iter, &BolOpt);
				t2 = clock();
				dtlp = (double)(t2 - t1) / CLOCKS_PER_SEC;
				dtot = dtlr + dtlp;
				TimeTot = dtot;

				if (SubLR.ZLBlp < -Giga)
					fprintf(fout, "%11.2lfG", SubLR.ZLBlp / Giga);
				else if (SubLR.ZLBlp < -Mega)
					fprintf(fout, "%11.2lfM", SubLR.ZLBlp / Mega);
				else
					fprintf(fout, "%11.2lf ", SubLR.ZLBlp);

				fprintf(fout, "%10.2lf", dtlp);
				fprintf(fout, "%10d", SubLR.nCore);
				fflush(fout);
			}

			// Solve by Incremental Core using MIP
			xt = new double[SubLR.nmax];
			Best_dt = +Inf;
			for (j = 0; j < NRip; j++)
			{
				t1 = clock();
				if (BolEsatto)
					SubLR.Solve_Incremental_Core(1, BolGreedy, xt, &iter, &BolOpt);
				t2 = clock();
				dt = (double)(t2 - t1) / CLOCKS_PER_SEC;
				if (dt < Best_dt)
				{
					Best_dt = dt;
					Best_Zlb = SubLR.ZLBlp;
					Best_Zub = SubLR.ZUBlp;
					Best_Gap = SubLR.Gap;
					Best_Nodes = SubLR.Nodes;
					Best_NCore = SubLR.nCore;
					Best_Iter = iter;
				}
				Delay(Tdelay);
			}
			dt = Best_dt;
			TimeMip = Best_dt;
			SubLR.ZLBlp = Best_Zlb;
			SubLR.ZUBlp = Best_Zub;
			SubLR.Gap = Best_Gap;
			SubLR.Nodes = Best_Nodes;
			SubLR.nCore = Best_NCore;
			iter = Best_Iter;

			dtot = dtlr + dtlp + dt;
			TimeTot = dtot;
			Zopt = Zmip;
			if (SubLR.ZUBlp < Zmip)
				Zopt = SubLR.ZUBlp;
			//GapOpt = 100. * abs(SubLR.ZUBlp - Zmip)/Zmip;
			GapOpt = 100. * (SubLR.ZUBlp - Zmip) / Zmip;
			GapGreedy = 100. * (SubLR.ZGreedy - Zopt) / Zopt;
			GapLagrHeu = 100. * (SubLR.ZUBlagr - Zopt) / Zopt;
			if (BolGreedy)
				dtot += tGreedy;

			if (SubLR.ZLBlp < -Giga)
				fprintf(fout, "%11.2lfG", SubLR.ZLBlp / Giga);
			else if (SubLR.ZLBlp < -Mega)
				fprintf(fout, "%11.2lfM", SubLR.ZLBlp / Mega);
			else
				fprintf(fout, "%11.2lf ", SubLR.ZLBlp);
			if (SubLR.ZUBlp > Giga)
				fprintf(fout, "%11.2lfG", SubLR.ZUBlp / Giga);
			else if (SubLR.ZUBlp > Mega)
				fprintf(fout, "%11.2lfM", SubLR.ZUBlp / Mega);
			else
				fprintf(fout, "%11.2lf ", SubLR.ZUBlp);
			fprintf(fout, "%10.4lf", SubLR.Gap);
			fprintf(fout, "%10d", SubLR.Nodes);
			fprintf(fout, "%10.2lf", dt);
			fprintf(fout, "%10d", SubLR.nCore);
			fprintf(fout, "%6d", iter);
			fprintf(fout, "%10.2lf", dtot);
			fprintf(fout, "%5d", BolOpt);
			fflush(fout);
			delete[] xt;

			fprintf(fout, "\n");

			// Print Tables
			fprintf(ftab2, "%-20s ; ", fname);
			fprintf(ftab2, "%10d ; ", SubLR.Inst.n);
			fprintf(ftab2, "%10d ; ", SubLR.Inst.W);
			fprintf(ftab2, "%10d ; ", SubLR.Inst.k);
			fprintf(ftab2, "%10.2lf ; ", SubLR.ZLBDFF);
			fprintf(ftab2, "%11.2lf ; ", SubLR.ZUBlagr);
			fprintf(ftab2, "%10.4lf ; ", SubLR.alpha);
			fprintf(ftab2, "%10d ; ", SubLR.Iter);
			fprintf(ftab2, "%10.2lf ; ", dtlr);
			fprintf(ftab2, "%10d ; ", SubLR.nCore);
			fprintf(ftab2, "%11.2lf ; ", SubLR.ZLBlp);
			fprintf(ftab2, "%11.2lf ; ", SubLR.ZUBlp);
			fprintf(ftab2, "%10.4lf ; ", SubLR.Gap);
			fprintf(ftab2, "%10d ; ", SubLR.Nodes);
			fprintf(ftab2, "%10.2lf ; ", dt);
			fprintf(ftab2, "%10d ; ", SubLR.nCore);
			fprintf(ftab2, "%6d ; ", iter);
			fprintf(ftab2, "%10.2lf ; ", dtot);
			fprintf(ftab2, "%5d ; \n", BolOpt);
			fflush(ftab2);

			// DeAllocate the data structures
			status = SubLR.DeAllocateDataStructure();
			if (status)
				exit(status);

			SubLR.Close();

		}
	}

	fclose(fpro);
	fclose(fout);
	fclose(ftab2);

	return 0;
}

/**
 *  EliminatePath: remove the path in the string "fname" containing the filename
 */
long EliminatePath(char *fname)
{
	long i, j;

	j = 0;
	for (i=0; fname[i]!='\0'; i++)
		if ((fname[i] == 92) || (fname[i] == 47) || (fname[i] == 46))
		{
			j = i + 1;
		}

	return j;
}

/**
 *  ReadRow: read a row from a text file
 */
long ReadRow(FILE *fin, char *row)
{
	long i = 0;
	while (true)
	{
		fscanf(fin, "%c", &(row[i]));
		if (row[i] == '\n')
		{
			row[i] = '\0';
			break;
		}
		i++;
	}

	return i;
}

/**
 *  Delay: wait for the given amount of milliseconds 
 */
void Delay(long milliseconds)
{
    long pause;
    clock_t now,then;

    pause = milliseconds*(CLOCKS_PER_SEC/1000);
    now = then = clock();
    while( (now-then) < pause )
        now = clock();
}

/**
 *  GenInstance: generate a BPP instance
 */
int GenInstance(long irun, long type, long nitems, long capbin, long kcons, Subgradient *SubLR)
{
	long i, j;
	double coef;

	// Assign general data
	SubLR->Inst.Type = type;
	SubLR->Inst.n = nitems;
	SubLR->Inst.W = capbin;
	SubLR->Inst.k = kcons;  

	// Use of "kcons" according the instance type
	// Type 0 - Basic BPP: No constraints and No used
	// Type 1 - Cardinality constrained BPP: in the same bin cannot be allocated more than "kcons" items
	// Type 2 - Class constrained BPP: in the same bin cannot be allocated more than "kcons" classes 
	// Type 3 - BPP with conflicts: "kcons" represents the percentage of conflicts [0,...,100]
	// Type 4 - Two-dimensional vector packing problem: "kcons" is the capacityfor the second weight

	// Allocate memory
	SubLR->Inst.w = new long[nitems];
	SubLR->Inst.cl = NULL;
	SubLR->Inst.CMatrix = NULL;
	SubLR->Inst.w2 = NULL;
	if (type == 2)
	{
		SubLR->Inst.cl = new long[nitems];
		SubLR->Inst.CL = kcons;
	}
	else if (type == 3)
	{
		SubLR->Inst.CMatrix = new long*[nitems];
		for (i = 0; i < nitems; i++)
			SubLR->Inst.CMatrix[i] = new long[nitems];
	}
	else if (type == 4)
	{
		SubLR->Inst.w2 = new long[nitems];
		SubLR->Inst.W2 = kcons;
	}

	// For each item generate weight et al.
	long seed = SubLR->Inst.n + SubLR->Inst.W + SubLR->Inst.k + 17 * (irun + 1);
	srand(seed);
	for (i = 0; i < nitems; i++)
	{
		// Generate weight
		SubLR->Inst.w[i] = long(((((double)(rand()) / RAND_MAX)*(0.5)) + 0.2) * capbin) + 1;
		//SubLR->Inst.w[i] = long(((((double)(rand()) / RAND_MAX)*(0.4)) + 0.1) * capbin) + 1;
		//SubLR->Inst.w[i] = long(((((double)(rand()) / RAND_MAX)*(1.0)) + 0.0) * capbin) + 1;
		//SubLR->Inst.w[i] = 100 * long(SubLR->Inst.w[i] / 100) + 1;
		if (SubLR->Inst.w[i] > capbin)
			SubLR->Inst.w[i] = capbin;
		if (SubLR->Inst.w[i] < 1)
			SubLR->Inst.w[i] = 1;

		if (type == 2)
		{
			// Generate classes
			coef = ((double)(rand()) / RAND_MAX);
			SubLR->Inst.cl[i] = long(coef * 2 * kcons) + 1; // The number of classes is the double of the maximum number of classes for each bin
			//if (SubLR->Inst.cl[i] > kcons)
			//	SubLR->Inst.cl[i] = kcons;
			if (SubLR->Inst.cl[i] < 1)
				SubLR->Inst.cl[i] = 1;
		}
		else if (type == 3)
		{
			// Generate conflit matrix (kcons is the pencentage of conflicts)
			SubLR->Inst.CMatrix[i][i] = 0;
			for (j = i + 1; j < nitems; j++)
			{
				SubLR->Inst.CMatrix[i][j] = long(((double)(rand()) / RAND_MAX * ((double)(kcons) / 100.)) + 0.5);
				if (SubLR->Inst.CMatrix[i][j] > 1)
					SubLR->Inst.CMatrix[i][j] = 1;
				if (SubLR->Inst.CMatrix[i][j] < 0)
					SubLR->Inst.CMatrix[i][j] = 0;
				SubLR->Inst.CMatrix[j][i] = SubLR->Inst.CMatrix[i][j];
			}
		}
		else if (type == 4)
		{
			// Generate weight
			SubLR->Inst.w2[i] = long(((double)(rand()) / RAND_MAX) * kcons/2.) + 1;
			//SubLR.Inst.w2[i] = 100 * long(SubLR.Inst.w[i] / 100) + 1;
			if (SubLR->Inst.w2[i] > kcons)
				SubLR->Inst.w2[i] = kcons;
			if (SubLR->Inst.w2[i] < 1)
				SubLR->Inst.w2[i] = 1;
		}
	}

	return 0;
}

/**
 *  ReadInstance: read an instance from a text file in a given format
 */
int ReadInstance(char *filename, Subgradient *SubLR)
{
	long i;
	FILE *fin;

	// Open text file containing the instance
	fin = fopen(filename, "r");
	if (fin == NULL)
		return 1;

	// Assign general data
	SubLR->Inst.Type = 0;  // Basic BPP

	// Read Number of Items
	fscanf(fin, "%ld", &(SubLR->Inst.n));

	// Read Bin Capacity
	fscanf(fin, "%ld", &(SubLR->Inst.W));

	SubLR->Inst.k = 0;

	// Allocate memory
	SubLR->Inst.w = new long[SubLR->Inst.n];
	SubLR->Inst.cl = NULL;
	SubLR->Inst.CMatrix = NULL;
	SubLR->Inst.w2 = NULL;

	// For each item generate weight et al.
	for (i = 0; i < SubLR->Inst.n; i++)
	{
		// Read weight
		fscanf(fin, "%ld", &(SubLR->Inst.w[i]));
	}

	fclose(fin);

	return 0;
}

