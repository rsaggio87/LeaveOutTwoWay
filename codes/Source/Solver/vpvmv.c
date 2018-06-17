/* ========================================================================== */
/* === CMG/Source/Solver/vpvmv.c      =====================================   */
/* ========================================================================== */

/* 	
 *   Pointwise addition of three vectors
 *
 */


/* CMG, Copyright (c) 2008-2010  Ioannis Koutis, Gary Miller                */
/* 
/* The CMG solver is distributed under the terms of the GNU General Public  */
/* Lincense Version 3.0 of the Free Software Foundation.                    */
/* CMG is also available under other licenses; contact authors for details. */


/* 
 * vpvmv
 * 
 * Input:  precision vectors: x,y,z
 *      :  size of vectors :n 
 *     
 * Output: precision vector: w = x+y-z
 *
 *
 */

#include "cmg.h"

void vpvmv (precision *x, precision *y, precision *z, precision *w, mSize n)
{
    mIndex i;
    
    for (i=0;i<n;i++) 
        w[i] = x[i]+y[i]-z[i];
    
}