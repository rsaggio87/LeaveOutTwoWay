/* ========================================================================== */
/* === CMG/Source/dspmv.c               ===================================== */
/* ========================================================================== */

/* 	
 *   Sparse matrix vector multiplication
 *
 */


/* 
 * spmv
 * 
 * Input:  double vectors: x,y
 *      :  size of vectors :n 
 *     
 * Output: double vector: z = x-y
 *
 *
 */

#include "cmg.h"
#include "malloc.h"

void dspmv (mSize n, double *a, mIndex *ia, mIndex *ja, double *x, double *y)
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