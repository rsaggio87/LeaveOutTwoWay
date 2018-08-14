/* ========================================================================== */
/* === CMG/Source/Solver/Double/fvtof.c ===================================== */
/* ========================================================================== */

/* 	
 *   Cast double vector to float
 *
 */


/* 
 * vinw
 * 
 * Input:  double vector v
 *      :  float vector w
 *      :  size of vectors :n 
 *     
 * Output: w = (float) v
 *
 *
 */

#include "cmg.h"
#include "malloc.h"

void dvtof (double *v, float *w, mSize n)
{
    mIndex i;
    
    for (i=0;i<n;i++) 
        w[i] = (float) v[i];
    
    
}