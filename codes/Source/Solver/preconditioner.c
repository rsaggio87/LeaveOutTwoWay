/* ========================================================================== */
/* === CMG/Source/Solver/preconditioner.c  ================================== */
/* ========================================================================== */

/* 	
 *   Solve given the hierarchy of graphs
 *
 */

/* CMG, Copyright (c) 2008-2010  Ioannis Koutis, Gary Miller                */
/* 
/* The CMG solver is distributed under the terms of the GNU General Public  */
/* Lincense Version 3.0 of the Free Software Foundation.                    */
/* CMG is also available under other licenses; contact authors for details. */



/* 
 * preconditioner
 * 
 * Input:  pointer to array of hierarchy levels: H
 *      :  --precision-- vector: b
 *      :  integer: l
 *
 * Output:  --precision-- vector: x 
 *
 *
 */

#include "cmg.h"



void preconditioner(s_hlevel *H, precision *b, int level, int iter, precision *x)
{
    mSize n,m,nc;
    int i,j;
    precision *y, *Bb, *r, *z;
    precision *b_small;
    precision s;
    /* char uplo,transa; */
  
    
    n = H[level].A.n;
    
    if ((H[level].islast ==1) && !(H[level].iterative)) {
        ldl_solve(&(H[level].chol) , b, x); 
        /* ssspmv(n, H[level].A.a, H[level].A.ia, H[level].A.ja, b, x); */
        if ((level == 0) && (H[level].laplacian==1))
          x[n]= (precision) 0.0;
        return;
    }
    
	

 
    if (H[level].iterative == 1){
        Bb = H[level].lws1;
        y = H[level].lws2;; 
        vvmul(H[level].invD,b,Bb,n);                                                 /*  Bb=H{level}.invD.*b */
        sspmv(n, H[level].A.a, H[level].A.ia, H[level].A.ja, Bb, y);
                                                                                     /* r = b-A*Bb */
                                                                                     /* x = Bb-r */
        vpvmv(Bb,y,b,x,n);                                                           /* x = Bb+y-b */
		vvmul(H[level].invD,b,x,n);
        return; 
    }

    /* main cycle */
    nc = H[level].nc;
	m = H[level+1].A.n;
    
    /* initialize the solution vector x */
    for (i=0;i<n;i++)
        x[i] = (precision) 0.0;
    
    Bb = H[level].lws1;                           /* fixed working space */
    vvmul(H[level].invD,b,Bb,n);                  /* Bb=H{level}.invD.*b */
 
    y  = H[level].lws2;                       /* y is used for temp storage */
    r  = H[level].lws2;
    for (j=1;j<=iter; j++){
         
        
        /* Jacobi pre-smooth */
        if (j==1){
                for (i=0;i<n;i++)
                x[i] = Bb[i];}
        else{
                   
        sspmv(n, H[level].A.a, H[level].A.ia, H[level].A.ja, x, y);
        vmv(b,y,r,n);                                                                   
        vvmul(H[level].invD,r,y,n);               /*  y = invD.*(b-A x)   */
        vpv(x,y,x,n);                             /*  x = x+invD.*(b-Ax) */        
        }
     
 
        b_small = H[level].sws1;  
        z = H[level].sws2;
        
        
        sspmv(n, H[level].A.a, H[level].A.ia, H[level].A.ja, x, y);
        vmv(b,y,r,n);                           /* r= b - (H{level}.A*x) */
        
        rmvecmul(H[level].cI , r, n, b_small, nc, H[level].dc);                       /*  b_small = Rt*r;  */
        preconditioner(H,b_small,level+1, H[level].repeat,z);                         /*  z = preconditioner(r) */
        trmvecmul(H[level].cI, z, nc, y, n, H[level].dc);                             /*  y = R*z */
        
        vpv(y,x,x,n);  
        
    
        /* Jacobi post-smooth */
        sspmv(n, H[level].A.a, H[level].A.ia, H[level].A.ja, x, y);
        vmv(b,y,r,n);                                                                   
        vvmul(H[level].invD,r,y,n);               /*  y = invD.*(b-A x)   */
        vpv(x,y,x,n);                             /*  x = x+invD.*(b-Ax) */        
    }
    
    if ((level == 0) && (H[level].laplacian==1)){
        x[n]= (precision)	 0.0; 
    }
    
    return;
 

}