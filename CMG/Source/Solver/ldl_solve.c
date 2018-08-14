/* ========================================================================== */
/* === CMG/Source/Solver/ldl_solve.c    ===================================== */
/* ========================================================================== */

/* 	
 *   Solve given Cholesky factorization 
 *
 */

/* CMG, Copyright (c) 2008-2010  Ioannis Koutis, Gary Miller                */
/* 
/* The CMG solver is distributed under the terms of the GNU General Public  */
/* Lincense Version 3.0 of the Free Software Foundation.                    */
/* CMG is also available under other licenses; contact authors for details. */

/* 
 * ldl_solve
 * 
 * Input:  cholesky factorization: ldl_p
 *      :  --precision-- vector: b
 *
 * Output:  --precision-- vector: x 
 *
 *
 */

#include "cmg.h"
#include "stdlib.h"

void ldl_solve (ldl_p *chol, precision *b, precision *x)
{
    mSize n;
    mIndex j;
    double *y;
    d_matrix ld, ldT;
    int i;
 

    n = chol->ld.n;
    ld = chol->ld;
    ldT = chol->ldT;
    y = malloc(n*sizeof(double));
    


    /* initialize y */
    for (i=0;i<n;i++) 
        y[i] = (double) b[chol->p[i]];
    
   
    for (i=1;i<n;i++) {
        for (j=ld.ia[i]; j<(ld.ia[i+1]-1); j++) 
           y[i] = y[i]-(ld.a[j])*y[ld.ja[j]];     
     }
    
    for (i=0;i<n;i++) {
       y[i]=y[i]/(ld.a[ld.ia[i+1]-1]);
    }

    for (i=(n-2);i>=0;i--) {
       for (j=(ldT.ia[i+1]-1); j>ldT.ia[i]; j--) {
          y[i] = y[i]-y[ldT.ja[j]]*(ldT.a[j]);}
    }
    

    for (i=0;i<n;i++)
        x[i]= (precision) y[chol->invp[i]];
    
    
    free(y);
 
    
       
    /* sspmv(n, ld.a, ld.ia, ld.ja, b, x); */
      
    
}