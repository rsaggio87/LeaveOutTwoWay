/* ========================================================================== */
/* === CMG/Source/Solver/Souble/dpcg.c  ===================================== */
/* ========================================================================== */

/* 	
 *   Implementation of preconditioned conjugate gradients
 *
 */


/* 
 * dpcg
 * 
 * Input:  pointer to array of hierarchy levels: H
 *      :  float vector: b
 *      :  integer: l
 *
 * Output:  float vector: x 
 *
 *
 */

#include "cmg.h"
#include "malloc.h"


void dpcg(mSize n, d_matrix *A, double *b, s_hlevel *H, unsigned int iter, double *xe)
{
    mSize m;
    unsigned int i;
    double *re, *ro, *pe, *po,  *xo;
    float *ze, *zo;
    double *y;
    float *yf;
    double r_in_zo, r_in_ze;
    double a,c;
  
    y = malloc(n*sizeof(double));
    re = malloc(n*sizeof(double));
    ro = malloc(n*sizeof(double));
    pe = malloc(n*sizeof(double));
    po = malloc(n*sizeof(double));
    xe = malloc(n*sizeof(double));
    xo = malloc(n*sizeof(double));
    ze = malloc(n*sizeof(float));
    zo = malloc(n*sizeof(float));
    
    /* initial vector is always taken to be zero */
    
    /* initialize */
    for (i=0; i<n; i++)  {
        xe[i]=0;
        re[i]=b[i];   /* residual = b - A x */
    }
    dvtof(re,yf,n);
    preconditioner(H,yf,0,ze);
    for (i=0; i<n; i++)  
         pe[i]=ze[i];
    
    r_in_ze = dsvinw(re,ze,n);
       
    
    for (i=1; i<iter; i++){
        
    if (i%2 == 1) {
        dspmv(n, A->a, A->ia, A->ja, pe, y); /* y holds A*p_e */
        a = r_in_ze/dvinw(y,pe,n); 
        daxpy(a,pe,xe,xo,n); 
        daxpy(-a,y,re,ro,n);
        
        /* here we need to add tol condition */
        dvtof(ro,yf,n);        
        preconditioner(H,yf,0,zo);        
        r_in_zo = dsvinw(re,ze,n);        
        c = r_in_zo/r_in_ze;
        daxpy(c,pe,ro,po,n);
    } 
    if (i%2 == 0) {
        dspmv(n, A->a, A->ia, A->ja, po, y); /* y holds A*p_o */
        a = r_in_zo/dvinw(y,po,n); 
        daxpy(a,po,xo,xe,n); 
        daxpy(-a,y,ro,re,n);
        
        /* here we need to add tol condition */
        dvtof(re,yf,n);        
        preconditioner(H,yf,0,ze);        
        r_in_ze = dsvinw(ro,zo,n);        
        c = r_in_ze/r_in_zo;
        daxpy(c,po,re,pe,n);
    } 
    }
    
    if (iter %2 == 1) 
        for (i=0;i<n; i++)
            xe[i] = xo[i];
            
        

  
    
    free(y);
    free(re);
    free(ro);
    free(pe);
    free(po);
    free(ze);
    free(zo);
    free(xo);


    
}