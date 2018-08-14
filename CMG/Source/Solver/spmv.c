/* ========================================================================== */
/* === CMG/Source/Solver/spmv.c         ===================================== */
/* ========================================================================== */

/* 	
 *   Sparse matrix vector multiplication
 *
 */

/* CMG, Copyright (c) 2008-2010  Ioannis Koutis, Gary Miller                */
/* 
/* The CMG solver is distributd under the terms of the GNU General Public   */
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


void spmv (mSize n, precision *a, mIndex *ia, mIndex *ja, precision *x, precision *y)
{
    mIndex i,j;
    
    /* Initialize y vector */
    for (i=0;i<n;i++){
        y[i]=0;}
    
   for (i=0;i<n;i++) {
        for (j=ia[i]; j<ia[i+1]; j++) 
            y[i] = y[i]+ a[j]*x[ja[j]]; 
    }

    return; 
   
}