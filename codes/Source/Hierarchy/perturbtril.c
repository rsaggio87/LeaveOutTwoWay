/* ========================================================================== */
/* === CMG/Source/Hierarchy/perturbtril.c       ============================= */
/* ========================================================================== */



/* CMG, Copyright (c) 2008-2010  Ioannis Koutis, Gary Miller                */
/* 
/* The CMG solver is distributd under the terms of the GNU General Public   */
/* Lincense Version 3.0 of the Free Software Foundation.                    */
/* CMG is also available under other licenses; contact authors for details. */


#include "mex.h"


#define A_IN prhs[0]
#define lb_IN prhs[1]
#define ub_IN prhs[2]
#define B_OUT plhs[0]



void mexFunction (int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{

    mwIndex *B_irow, *B_jcol;
    double *B, lb, ub;
    double  rf; 
    mwIndex m,n,i,j;
    mxArray *random_GET[1], *random_SEND[2];
    double *random_source;

    
       

    if ((nrhs != 3))
        mexErrMsgTxt("Wrong number of input arguments");
    if (nlhs > 1)
        mexErrMsgTxt("Wrong number of output arguments");
    
    if (!mxIsSparse(A_IN))
        mexErrMsgTxt("Input must be a sparse matrix");

    if ( (nrhs>=2) && (!(mxIsDouble(lb_IN)) || mxIsSparse(lb_IN) )) 
      mexErrMsgTxt("Second and third arguments must be scalars");
    
    
    if ( (nrhs==3) && (!(mxIsDouble(ub_IN)) || mxIsSparse(ub_IN) )) 
      mexErrMsgTxt("Second and third arguments must be scalars");
    
    
    if ( ((nrhs >=2) && (mxGetN(lb_IN)!=1)) || ((nrhs ==3) && (mxGetN(ub_IN)!=1)) ) 
      mexErrMsgTxt("Second and third arguments must be scalars");
    
     
    
    lb = mxGetPr(lb_IN)[0];
    ub = mxGetPr(ub_IN)[0];
    if (ub <= lb)
        mexErrMsgTxt("Third argument must be larger than second");
    
        
    m = (mwSize) mxGetM(A_IN);
    n = (mwSize) mxGetN(A_IN);

    if (m!=n)
        mexErrMsgTxt("Input must be symmetric");
    

    /* Create output by duplicating input */
    B_OUT = mxDuplicateArray(A_IN);
    B_irow = mxGetIr(B_OUT);
    B_jcol = mxGetJc(B_OUT);
    B = mxGetPr(B_OUT);

    
    
 
    /* Get random sourse */ 
    random_SEND[0] = mxCreateDoubleScalar(B_jcol[n]);
    random_SEND[1] = mxCreateDoubleScalar(1);
    mexCallMATLAB(1,random_GET,2,random_SEND,"rand");
    random_source = mxGetPr(random_GET[0]); 
    
    
 
    for (j=0;j<n;j++)
       for (i=B_jcol[j]; i<B_jcol[j+1]; i++) 
           if ( j< B_irow[i] )
               B[i] = (lb +(ub-lb)*random_source[i])*B[i];
           else 
               B[i] = 0.0;
    
    /* Deallocate memory */
    mxDestroyArray(random_GET[0]);
    mxDestroyArray(random_SEND[0]);
    mxDestroyArray(random_SEND[1]); 
}   
              
