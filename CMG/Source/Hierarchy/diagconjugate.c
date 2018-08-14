
/* ========================================================================== */
/* === CMG/Source/Solver/Hierarchy/diagconjugate.c ========================   */
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

    mwIndex *B_irow, *B_jcol;
    mwIndex i,j;
    mwSize m,n;
    double *B, *D,  mult;
     

    /* Input validation */
    if ((nrhs != 2))
        mexErrMsgTxt("Wrong number of input arguments");
    if (nlhs > 1)
        mexErrMsgTxt("Wrong number of output arguments");
    
    if (!mxIsSparse(A_IN))
        mexErrMsgTxt("Input must be a sparse matrix");

    if ( !mxIsDouble(D_IN) || mxIsSparse(D_IN) ) 
      mexErrMsgTxt("Second arguments must be a full vector");
   
    m = (mwSize) mxGetM(A_IN);
    n = (mwSize) mxGetN(A_IN);

    if (m != n)
        mexErrMsgTxt("Input matrix must be square");
    
    if ( (mxGetM(D_IN) !=1) && (mxGetN(D_IN)!=1)) 
        mexErrMsgTxt("Second arguments must be a vector");

    if (( mxGetM(D_IN) !=n) && (mxGetN(D_IN)!=n)) 
        mexErrMsgTxt("Vector dimensions must agree with matrix dimensions");

    
 
    B_OUT = mxDuplicateArray(A_IN);
    B_irow = mxGetIr(B_OUT);
    B_jcol = mxGetJc(B_OUT);
    B = mxGetPr(B_OUT);
    D = mxGetPr(D_IN);
    
    
    for (j=0;j<n;j++)
       for (i=B_jcol[j]; i<B_jcol[j+1]; i++) {
            if (j< B_irow[i])
                mult = D[j]*D[B_irow[i]];
            else
                mult = D[B_irow[i]]*D[j];
            B[i] = B[i]*mult;
       }                 
    
 
}   
              
