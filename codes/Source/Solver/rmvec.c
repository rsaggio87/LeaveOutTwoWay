/* ========================================================================== */
/* === CMG/Source/Solver/rmvec.c        ===================================== */
/* ========================================================================== */


/* 	
 *  projection and projection-transpose operators
 */


/* CMG, Copyright (c) 2008-2010  Ioannis Koutis, Gary Miller                */
/* 
/* The CMG solver is distributed under the terms of the GNU General Public  */
/* Lincense Version 3.0 of the Free Software Foundation.                    */
/* CMG is also available under other licenses; contact authors for details. */


/* 
 * rmvecmul
 * 
 * Input: cluster-membership vector: *CI
 *      : --precision-- vector: *x
 *      : number of clusters: m
 *      : dimension of x: n
 *
 * Output: --precision-- vector: *y
 *
 * Function: If R is the mxn 0-1 representing CI, y<-R*x
 *
 */

 /* constraint: max(CI) = m */


#include "cmg.h"

void rmvecmul (mIndex *CI, precision *x, mSize  n, precision *y, mSize m, boolean dc)
{
    mIndex i;
    
    /* initialize y */
    for (i=0;i<m;i++) 
        y[i]= (precision) 0;
    
    if (!dc) {
        for (i=0;i<n;i++)
        y[CI[i]] = y[CI[i]]+x[i];
    }
    else { 
        for (i=0;i<n;i++)
            if (CI[i]<(m-1))
                y[CI[i]] = y[CI[i]]+x[i];
    }
        
}



/* 
 * trmvecmul
 * 
 * Input: cluster-membership vector: *CI
 *      : --precision-- vector: *x
 *      : number of clusters: m
 *      : dimension of x: n
 *
 * Output: --precision-- vector: *y
 *
 * Function: If R is the mxn 0-1 representing CI, y<-R'*x
 */
 
void trmvecmul (mIndex *CI, precision *x, mSize  m, precision *y, mSize n, boolean dc)
{
    mIndex i;
    
	
    if (!dc) {
        for (i=0;i<n;i++)
            y[i] = x[CI[i]];
    }
    else {
        for (i=0;i<n;i++) {
            if (CI[i]<(m-1)) 
			{y[i] = x[CI[i]];}
            else 
			{y[i] = (precision) 0.0;}
        }
    }
}
