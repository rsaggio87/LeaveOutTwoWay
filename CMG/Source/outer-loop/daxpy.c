/* ========================================================================== */
/* === CMG/Source/daxpy.c               ===================================== */
/* ========================================================================== */



/* 
 * daxpy
 * 
 * Input:  double vectors: x,y
 *      :  double a
 *      :  size of vectors :n 
 *     
 * Output: double vector: z = a*x+y
 *
 *
 */

#include "cmg.h"
#include "malloc.h"

void daxpy (double a, double *x, double *y, double *z, mSize n)
{
    mIndex i;
    
    for (i=0;i<n;i++) 
        z[i] = a*x[i]+y[i];
    
}