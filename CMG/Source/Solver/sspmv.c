/* ========================================================================== */
/* === CMG/Source/Solver/spmv.c         ===================================== */
/* ========================================================================== */

/* 	
 *   Sparse matrix vector multiplication
 *   Matrix is assumed to be symmetric and only upper triangular part is stored
 *
 */

/* CMG, Copyright (c) 2008-2010  Ioannis Koutis, Gary Miller                */
/* 
/* The CMG solver is distributed under the terms of the GNU General Public  */
/* Lincense Version 3.0 of the Free Software Foundation.                    */
/* CMG is also available under other licenses; contact authors for details. */


/* 
 * spmv
 * 
 * Input:  --precision-- vectors: x,y
 *      :  size of vectors :n 
 *     
 * Output: --precision-- vector: z = x-y
 *
 *
 */

#include "cmg.h"


void sspmv (mSize n, precision *a, mIndex *ia, mIndex *ja, precision *x, precision *y)
{
    mIndex i,j,k;
    precision sum;
    
    /* Initialize y vector */
    for (i=0;i<n;i++){
        y[i]=0;}
    
    
   for (i=0;i<n;i++) {
        sum = a[ia[i]]*x[i];
        for (j=ia[i]+1; j<ia[i+1]; j++) {
            k = ja[j];
            sum += a[j]*x[k]; 
            y[k] += a[j]*x[i];}
        y[i] +=sum;
    }

    return; 
   
}