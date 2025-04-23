//-----------------------------------------------------------------------------
// Here we must insert our "license boilerplate"
// (... to be defined)
//-----------------------------------------------------------------------------
// File: Utilities.cpp
// Version 1.0
// Last update: 18.04.2025
// Authors: Boschetti M.A., Novellani, S.
//-----------------------------------------------------------------------------
//
#include "Utilities.h"

float GetTime(void)
{
	clock_t t;

	t = clock();

	return ((float)t) / CLOCKS_PER_SEC;
}

/*-----------------------------------------------------------------------------
//  Generate a delay (used in the sync-mechanism)
//-----------------------------------------------------------------------------
*/
void delay(int milli_seconds)
{
	// Storing start time
	clock_t start_time = clock();

	// looping till required time is not achieved
	while (clock() < start_time + milli_seconds);
}

/**
 *  Sort Value(one dimensional)
 *
 *	@param  nval	number of values
 *	@param  val		values
 *	@param  ind		index of the i-th ordered value
 *	@param	sense	-1 increasing; +1 decreasing
 */
int SortValues(long nval, double *val, long *ind, int sense)
{
	long i;
	Sorter Data;

	Data.NLevel = 1;
	Data.NItem = nval;
	Data.Init();

	for (i = 0; i < nval; i++)
	{
		Data.Item[i].Pos = i;
		Data.Item[i].Val[0] = val[i];
	}

	Data.Sort(sense);

	for (i = 0; i < nval; i++)
	{
		ind[i] = Data.Item[i].Pos;
	}

	return 0;
}

/**
 *	Sorter: Constructor
 */
Sorter::Sorter(void)
{
	NLevel = 0;
	NItem = 0;
	Item = NULL_PTR;
}

/**
 *  Sorter: Destructor
 */
Sorter::~Sorter(void)
{
	NLevel = 0;
	NItem = 0;
	if (Item != NULL_PTR)
		delete[] Item;
}

/**
 *  Sorter: Initialize
 */
int Sorter::Init(void)
{
	Item = new QSData_STR[NItem];

	return 0;
}

/**
 *	Sort item:
 *
 *	@param sense -1 increasing; +1 decreasing.
 */
int Sorter::Sort(int sense)
{
	QSort(Item, 0, NItem - 1, NLevel, sense);

	return 0;
}

/**
 *  Quick sort for Double Variables
 *
 *  @param sense -1 increasing; +1 decreasing
 */
void Sorter::QSort(QSData_STR *a, long first, long last, long QSNVal, int sense)
{
	long mid;

	if (last >= first + 1)
	{
		mid = Split(a, first, last, QSNVal, sense);
		QSort(a, first, mid - 1, QSNVal, sense);
		QSort(a, mid + 1, last, QSNVal, sense);
	}
}

/**
 * 	Quick sort for Double Varibles
 *	Fuction for perfoming the splitting
 */
int Sorter::Split(QSData_STR *a, long first, long last, long QSNVal, int sense)
{
	long left, right;
	QSData_STR temp;

	left = first + 1;
	right = last;

	for (;;)
	{
		while ((left <= last) && (QSCmp(&a[left], &a[first], QSNVal, sense) < 0))
			left++;
		while ((QSCmp(&a[right], &a[first], QSNVal, sense) > 0))
			right--;

		if (left < right)
		{
			temp = a[left];
			a[left] = a[right];
			a[right] = temp;
			left++;
			right--;
		}
		else
		{
			temp = a[first];
			a[first] = a[right];
			a[right] = temp;
			return right;
		}
	}
}

/**
 *  Quick sort for Double Varibles
 *  Fuction for comparing two items
 */
int Sorter::QSCmp(QSData_STR *pi, QSData_STR *pj, long QSNVal, int sense)
{
	long i;

	for (i = 0; i < QSNVal; i++)
	{
		if (pi->Val[i] > pj->Val[i])
			return -sense;
		if (pi->Val[i] < pj->Val[i])
			return sense;
	}

	return 0;
}


