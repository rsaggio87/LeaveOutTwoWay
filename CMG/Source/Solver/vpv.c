/* ========================================================================== */
/* === CMG/Source/Solver/vpv.c          ===================================== */
/* ========================================================================== */

/* 	
 *   Pointwise subtraction of vectors
 *
 */

/* CMG, Copyright (c) 2008-2010  Ioannis Koutis, Gary Miller                */
/* 
/* The CMG solver is distributed under the terms of the GNU General Public  */
/* Lincense Version 3.0 of the Free Software Foundation.                    */
/* CMG is also available under other licenses; contact authors for details. */


/* 
 * vpv
 * 
 * Input:  --precision-- vectors: x,y
 *      :  size of vectors :n 
 *     
 * Output: --precision-- vector: z = x+y
 *
 *
 */

#include "cmg.h"


void vpv (precision *x, precision *y, precision *z, mSize n)
{
    mIndex i;
    
    for (i=0;i<n;i++) 
        z[i] = x[i]+y[i];
    
}