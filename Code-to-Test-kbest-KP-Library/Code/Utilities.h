//-----------------------------------------------------------------------------
// Here we must insert our "license boilerplate"
// (we need to investigate what is the most suitable licence for us
//-----------------------------------------------------------------------------
// File: Utilities.h
// Version 1.0
// Last update: 01.12.2022
// Authors: Boschetti M.A., Novellani, S.
//-----------------------------------------------------------------------------
//
#ifndef __UTILITIES_DEFS_H
#define __UTILITIES_DEFS_H

#include <time.h>

#define NULL_PTR 0 

float GetTime(void);
void delay(int milli_seconds);

int SortValues(long nval, double *val, long *ind, int sense);

/**
 *	@brief @b Class for sorting items.
 */
class Sorter
{

	/// @b Struct: Data structure used in class Sorter (for sorting items)
	struct QSData_STR
	{
		long Pos; 		/*!< Original position of the value that should be sorted */
		double Val[8]; 	/*!< Value that should be lexicographically sorted */
		long Ind[4]; 	/*!< index of the original data that should be retrieved */
	};

private:

	void QSort(QSData_STR *a, long first, long last, long QSNVal, int sense);
	int Split(QSData_STR *a, long first, long last, long QSNVal, int sense);
	int QSCmp(QSData_STR *pi, QSData_STR *pj, long QSNVal, int sense);

public:

	long NLevel;   		/*!< Number of levels that should be considered in the sorting */
	long NItem;    		/*!< number of items to be sorted */
	QSData_STR *Item;  	/*!< Item to be sorted */

	// Constructor and Destructorr
	Sorter(void);
	~Sorter(void);

	// Initialize
	int Init(void);

	// Sort item
	int Sort(int sense);

};

#endif /* __UTILITIES_DEFS_H */

