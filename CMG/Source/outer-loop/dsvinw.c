/* ========================================================================== */
/* === CMG/Source/Solver/Double/dsvinw.c ===================================== */
/* ========================================================================== */

/* 	
 *   Inner product of vectors
 *
 */


/* 
 * vinw
 * 
 * Input:  double vector v
 *      :  float vector w
 *      :  size of vectors :n 
 *     
 * Output: double scalar: z = v^T w
 *
 *
 */

#include "cmg.h"
#include "malloc.h"

double dsvinw (double *v, float *w, mSize n)
{
    mIndex i;
    double z;
    
    z = 0;
    for (i=0;i<n;i++) 
        z = z+v[i]*( (double) w[i]);
    
    return z;
    
}