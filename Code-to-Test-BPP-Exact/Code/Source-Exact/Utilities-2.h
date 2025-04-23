//-----------------------------------------------------------------------------
//  File: Utilities-2.h
//
//  Some sorting functions
//
//  Authors: M.A.Boschetti and S. Novellani         
//
//  Last update: 18.04.2025 
//-----------------------------------------------------------------------------
//
void QSort(long *a, long *ind, int first, int last, int sense);
int Split(long *a, long *ind, int first, int last, int sense);
int QSCmp(long *pi, long *pj, int sense);
void QSortD(double *a, long *ind, int first, int last, int sense);
int SplitD(double *a, long *ind, int first, int last, int sense);
int QSCmpD(double *pi, double *pj, int sense);
