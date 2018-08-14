/* ========================================================================== */
/* === CMG/Source/Hierarchy/adjacency.c   ==================================  */
/* ========================================================================== */


/* CMG, Copyright (c) 2008-2010  Ioannis Koutis, Gary Miller                */
/* 
/* The CMG solver is distributd under the terms of the GNU General Public   */
/* Lincense Version 3.0 of the Free Software Foundation.                    */
/* CMG is also available under other licenses; contact authors for details. */


#include "mex.h"


#define A_IN prhs[0]
#define B_OUT plhs[0]



void mexFunction (int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{

    mwIndex *B_irow, *B_jcol, *A_jcol, *A_irow;
    mwIndex i,j, ind, flag;
    mwSize m,n;
    double *B, *A;
      

    if ((nrhs != 1))
        mexErrMsgTxt("Wrong number of input arguments");
    if (nlhs > 1)
        mexErrMsgTxt("Wrong number of output arguments");
    
    if (!mxIsSparse(A_IN))
        mexErrMsgTxt("Input must be a sparse matrix");
  
           
    m = (mwSize) mxGetM(A_IN);
    n = (mwSize) mxGetN(A_IN);

    if (m!=n)
        mexErrMsgTxt("Input must be symmetric");
    
    /* Get Input */
    A_irow = mxGetIr(A_IN);
    A_jcol = mxGetJc(A_IN);
    A = mxGetPr(A_IN);
    

    /* Create output */
    B_OUT = mxCreateSparse(n,n,A_jcol[n]-n,mxREAL);
    B_irow = mxGetIr(B_OUT);
    B_jcol = mxGetJc(B_OUT);
    B = mxGetPr(B_OUT);
    
    for (j=0;j<=n;j++)
     B_jcol[j] = A_jcol[j]-j;

 
    ind = 0;  
    for (j=0;j<n;j++) {
       flag = 0;
       for (i=A_jcol[j]; i<A_jcol[j+1]; i++) { 
           if ( j!= A_irow[i] ){
               B[ind] = -A[i];
               B_irow[ind] = A_irow[i];
               ind = ind+1;}
           else 
               flag = 1;}
       if (flag == 0)
           mexErrMsgTxt("Laplacian missing a diagonal element");  }
 
}   
              
