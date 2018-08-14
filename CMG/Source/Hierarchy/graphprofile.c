/* ========================================================================== */
/* === CMG/Source/Hierarchy/graphprofile.c      ============================= */
/* ========================================================================== */



/* CMG, Copyright (c) 2008-2010  Ioannis Koutis, Gary Miller                */
/* 
/* The CMG solver is distributd under the terms of the GNU General Public   */
/* Lincense Version 3.0 of the Free Software Foundation.                    */
/* CMG is also available under other licenses; contact authors for details. */


#include "mex.h"


#define A_IN prhs[0]
#define P_IN prhs[1]
#define B_OUT plhs[0]
#define D_OUT plhs[1]
#define M_OUT plhs[2]




void mexFunction (int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{

    mwIndex *B_irow, *B_jcol, *A_irow, *A_jcol;
    mwIndex i,j,cmax_i;
    mwSize n;
    uint32_T *B;
    double *A, *P, *D, *M, cmax, candmax, colsum;
      

    if ((nrhs !=1) && (nrhs !=2))
        mexErrMsgTxt("Wrong number of input arguments");
    if ((nlhs > 3) || (nlhs < 1))
        mexErrMsgTxt("Wrong number of output arguments");
    
    if (!mxIsSparse(A_IN))
        mexErrMsgTxt("Input must be a sparse matrix");

    n = mxGetN(A_IN);
    if (( n!=mxGetM(A_IN)) || ( n!=mxGetN(A_IN)))
        mexErrMsgTxt("Matrix must be symmetric"); 
    
    
    if (nrhs == 2) {
    if ((!mxIsSparse(P_IN)) || (mxGetM(P_IN)!=n) || (mxGetN(P_IN) !=n))
        mexErrMsgTxt("Matrix dimensions must be identical"); }

    A_irow = mxGetIr(A_IN);
    A_jcol = mxGetJc(A_IN);
    A = mxGetPr(A_IN);

    
    /* Create output */
    B_OUT = mxCreateNumericMatrix( 1, n, mxUINT32_CLASS, mxREAL);
    D_OUT = mxCreateDoubleMatrix(n,1,mxREAL);
    M_OUT = mxCreateDoubleMatrix(n,1,mxREAL);
    if ((D_OUT == NULL) || (B_OUT == NULL) || (M_OUT == NULL))
         mexErrMsgTxt("Unsufficient memory");
    B = (uint32_T *) mxGetPr(B_OUT);
    D = mxGetPr(D_OUT);
    M = mxGetPr(M_OUT);
    
   
    if (nrhs == 1){
    for (j=0 ; j<n ; j++){
		cmax_i = n+1;  /* initialiaze */
        cmax = 0.0;
        colsum = 0.0;
        for (i=A_jcol[j]; i<A_jcol[j+1]; i++) {
            colsum = colsum+A[i];
            if (A[i]>cmax){
             cmax_i = A_irow[i];
             cmax = A[i];}}
        D[j] = colsum;
        B[j] = (uint32_T) cmax_i; 
        M[j] = cmax;
    }}
    
    if (nrhs == 2){
    P = mxGetPr(P_IN);
    for (j=0 ; j<n ; j++){
        cmax_i = n+1;  /* initialiaze */
        cmax = 0.0;
        colsum = 0.0;
        for (i=A_jcol[j]; i<A_jcol[j+1]; i++) {
            colsum = colsum+A[i];
            candmax = P[i]*A[i];
            if (candmax>cmax){
             cmax_i = A_irow[i];
             cmax = candmax;}}
        D[j] = colsum;
        B[j] = (uint32_T) cmax_i;
        M[j] = cmax;
    }}    
}   
              
