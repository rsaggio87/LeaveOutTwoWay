/* ========================================================================== */
/* === CMG/Source/Hierarchy/update_groups.c     ============================= */
/* ========================================================================== */


/* CMG, Copyright (c) 2008-2010  Ioannis Koutis, Gary Miller                */
/* 
/* The CMG solver is distributed under the terms of the GNU General Public   */
/* Lincense Version 3.0 of the Free Software Foundation.                    */
/* CMG is also available under other licenses; contact authors for details. */


#include "mex.h"

#define A_IN prhs[0]
#define C_IN prhs[1]
#define C_OUT plhs[0]
#define D_OUT plhs[1]


void mexFunction (int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{

    mwIndex *A_irow, *A_jcol;
    mwIndex iA,j;
    mwSize n,m;
    double *D, *A, colsum, *Bcolsum;
    uint32_T *Cin, *Cout;
      

    if (nrhs !=2)
        mexErrMsgTxt("Wrong number of input arguments");
    if ((nlhs != 2))
        mexErrMsgTxt("Wrong number of output arguments");
    
    if (!mxIsSparse(A_IN))
        mexErrMsgTxt("First input must be a sparse matrix");
    
    n = mxGetN(C_IN);
    m = mxGetM(C_IN);
    if ((m!=1) && (n!=1))
        mexErrMsgTxt("Second input must be a vector");
    if (mxGetClassID(C_IN) != mxUINT32_CLASS)
        mexErrMsgTxt("Input vector must be of type uint32");
    Cin = (uint32_T *) mxGetPr(C_IN);
  
   
    A_irow = mxGetIr(A_IN);
    A_jcol = mxGetJc(A_IN);
    A =  mxGetPr(A_IN);
  
    /* Create output */
    C_OUT = mxCreateNumericMatrix(m, n, mxUINT32_CLASS, mxREAL);
    D_OUT = mxCreateDoubleMatrix(n,1,mxREAL);
    if ((D_OUT == NULL) || (C_OUT == NULL)) 
         mexErrMsgTxt("Unsufficient memory");
    Cout = (uint32_T *) mxGetPr(C_OUT);
    D = mxGetPr(D_OUT);
    

    Bcolsum = mxMalloc(n*sizeof(double));
 
    
    n = mxGetN(A_IN);
    for (j=0 ; j<n ; j++){
        colsum = 0.0;
        Bcolsum[j] = 0.0;
        for (iA=A_jcol[j]; iA<A_jcol[j+1]; iA++) {
            colsum = colsum+A[iA];
            if (Cin[j] == (uint32_T) A_irow[iA]){
                Bcolsum[j] = Bcolsum[j]+A[iA];
                Bcolsum[A_irow[iA]] = Bcolsum[A_irow[iA]] + A[iA];}
            }
        D[j] = colsum;
    }
    
    for (j=0 ; j<n; j++){
        if ((Bcolsum[j]/D[j])<0.125) 
            Cout[j] = (uint32_T) j;
        else
            Cout[j] =  Cin[j];}
    

    mxFree(Bcolsum);
}    
        
              
