/* ========================================================================== */
/* === CMG/Source/Hierarchy/laplacian2.c        ============================= */
/* ========================================================================== */



/* CMG, Copyright (c) 2008-2010  Ioannis Koutis, Gary Miller                */
/* 
/* The CMG solver is distributd under the terms of the GNU General Public   */
/* Lincense Version 3.0 of the Free Software Foundation.                    */
/* CMG is also available under other licenses; contact authors for details. */



#include "mex.h"


#define A_IN prhs[0]
#define D_IN prhs[1]
#define B_OUT plhs[0]



void mexFunction (int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{

    mwIndex *B_irow, *B_jcol, *A_jcol, *A_irow;
    mwIndex i,j, ind, flag;
    mwSize n;
    double *B, *A, *D;
      

    if ((nrhs != 2))
        mexErrMsgTxt("Wrong number of input arguments");
    if (nlhs > 1)
        mexErrMsgTxt("Wrong number of output arguments");
    
    if (!mxIsSparse(A_IN))
        mexErrMsgTxt("Input must be a sparse matrix");
    
    if (!(mxIsDouble(D_IN)) || mxIsSparse(D_IN)  ) {
      mexErrMsgTxt("Second argument must be a non-sparse vector");
    }
    
  
           
    n = (mwSize) mxGetN(A_IN);

    if (n!=mxGetM(A_IN))
        mexErrMsgTxt("Input must be symmetric");
    
    if (((n!=mxGetN(D_IN)) || (mxGetM(D_IN)!=1)) && ((n!=mxGetM(D_IN)) || (mxGetN(D_IN)!=1)))
        mexErrMsgTxt("Inconsistent argument dimensions");
    
    /* Get Input */
    A_irow = mxGetIr(A_IN);
    A_jcol = mxGetJc(A_IN);
    A = mxGetPr(A_IN);
    D = mxGetPr(D_IN);
    

    /* Create output */
    B_OUT = mxCreateSparse(n,n,A_jcol[n]+n,mxREAL);
    B_irow = mxGetIr(B_OUT);
    B_jcol = mxGetJc(B_OUT);
    B = mxGetPr(B_OUT);
    
    for (j=0;j<=n;j++)
     B_jcol[j] = A_jcol[j]+j;

 
    ind = 0;  
    for (j=0;j<n;j++) {
       i=A_jcol[j];
       while((A_irow[i]<=j) && (i<A_jcol[j+1])){
           B[ind] = -A[i];
           B_irow[ind] = A_irow[i];
           ind = ind+1;
           i = i+1; }
       if (A_irow[i]==j)
           mexErrMsgTxt("Adjacency contains a diagonal element");
       B[ind] = D[j];
       B_irow[ind] = j;       
       ind = ind+1;
       while(i<A_jcol[j+1]){
           B[ind] = -A[i];
           B_irow[ind] = A_irow[i];
           ind = ind+1;
           i = i+1;}}

 
}   
              
