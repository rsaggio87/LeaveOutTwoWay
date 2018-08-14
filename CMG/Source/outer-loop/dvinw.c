/* ========================================================================== */
/* === CMG/Source/Solver/Double/dvinw.c ===================================== */
/* ========================================================================== */

/* 	
 *   Inner product of vectors
 *
 */


/* 
 * dvinw
 * 
 * Input:  double vectors: v,w
 *      :  size of vectors :n 
 *     
 * Output: double scalar: z = v^T w
 *
 *
 */

#include "cmg.h"
#include "malloc.h"

double dvinw (double *v, double *w, mSize n)
{
    mIndex i;
    double z;
    
    z = 0;
    for (i=0;i<n;i++) 
        z = z+v[i]*w[i];
    
    return z;
    
}