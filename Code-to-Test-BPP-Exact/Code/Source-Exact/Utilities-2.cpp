//-----------------------------------------------------------------------------
//  File: Utilities-2.cpp
//
//  Some sorting functions
//
//  Authors: M.A.Boschetti and S. Novellani         
//
//  Last update: 18.04.2025 
//-----------------------------------------------------------------------------
//

#include "Utilities-2.h"

//--------------------------------------------------------------------
//  Quick sort for "long" Varibles 
//  (sense: -1 increasing; +1 decreasing)
//--------------------------------------------------------------------

void QSort(long *a, long *ind, int first, int last, int sense)
{
	int mid;

	if (last >= first + 1)
	{
		mid = Split(a, ind, first, last, sense);
		QSort(a, ind, first, mid - 1, sense);
		QSort(a, ind, mid + 1, last, sense);
	}
}

int Split(long *a, long *ind, int first, int last, int sense)
{
	int left, right;
	long temp;

	left = first + 1;
	right = last;

	for (;;)
	{
		while ((left <= last) && (QSCmp(&a[ind[left]], &a[ind[first]], sense) < 0))
			left++;
		while (QSCmp(&a[ind[right]], &a[ind[first]], sense) > 0)
			right--;

		if (left < right)
		{
			temp = ind[left];
			ind[left] = ind[right];
			ind[right] = temp;
			left++;
			right--;
		}
		else
		{
			temp = ind[first];
			ind[first] = ind[right];
			ind[right] = temp;
			return right;
		}
	}
}

int QSCmp(long *pi, long *pj, int sense)
{
	if (*pi>*pj)
		return -sense;
	if (*pi<*pj)
		return sense;

	return 0;
}

//--------------------------------------------------------------------
//  Quick sort for "double" Varibles 
//  (sense: -1 increasing; +1 decreasing)
//--------------------------------------------------------------------

void QSortD(double *a, long *ind, int first, int last, int sense)
{
	int mid;

	if (last >= first + 1)
	{
		mid = SplitD(a, ind, first, last, sense);
		QSortD(a, ind, first, mid - 1, sense);
		QSortD(a, ind, mid + 1, last, sense);
	}
}

int SplitD(double *a, long *ind, int first, int last, int sense)
{
	int left, right;
	long temp;

	left = first + 1;
	right = last;

	for (;;)
	{
		while ((left <= last) && (QSCmpD(&a[ind[left]], &a[ind[first]], sense) < 0))
			left++;
		while (QSCmpD(&a[ind[right]], &a[ind[first]], sense) > 0)
			right--;

		if (left < right)
		{
			temp = ind[left];
			ind[left] = ind[right];
			ind[right] = temp;
			left++;
			right--;
		}
		else
		{
			temp = ind[first];
			ind[first] = ind[right];
			ind[right] = temp;
			return right;
		}
	}
}

int QSCmpD(double *pi, double *pj, int sense)
{
	if (*pi > *pj)
		return -sense;
	if (*pi < *pj)
		return sense;

	return 0;
}

