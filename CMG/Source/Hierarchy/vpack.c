/* ========================================================================== */
/* === CMG/Source/Hierarchy/vpack.c             ============================= */
/* ========================================================================== */



/* CMG, Copyright (c) 2008-2010  Ioannis Koutis, Gary Miller                */
/* 
/* The CMG solver is distributd under the terms of the GNU General Public   */
/* Lincense Version 3.0 of the Free Software Foundation.                    */
/* CMG is also available under other licenses; contact authors for details. */



#include "mex.h"

#define MAX(a,b)  (((a) > (b)) ? (a) : (b))
#define MIN(a,b)  (((a) < (b)) ? (a) : (b))

#define v_IN prhs[0]
#define w_OUT plhs[0]



void mexFunction (int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{

    mwIndex i, first_not_used, last_used;
    mwSize m,n;
    bool *v , vc;
    uint32_T *w ,wc;

      

    if ( nrhs != 1 )
        mexErrMsgTxt("Wrong number of input arguments");
    if (nlhs > 1)
        mexErrMsgTxt("Wrong number of output arguments");
    
    if (!mxIsClass(v_IN,"logical"))
        mexErrMsgTxt("Input must be a boolean vector");
              
    m = (mwSize) mxGetM(v_IN);
    n = (mwSize) mxGetN(v_IN);
    
    if ( (m!=1) && (n!=1) )
        mexErrMsgTxt("Input must be a boolean vector");
 
    v = (bool *) mxGetPr(v_IN);
    w_OUT = mxCreateNumericMatrix(m,n,mxUINT32_CLASS,mxREAL);

    if (w_OUT == NULL)
        mexErrMsgTxt("Out of memory");
    w = (uint32_T *) mxGetPr(w_OUT);
 
    first_not_used=0;
    last_used = 1;
    for (i=0; i<MAX(m,n); i++){
        if ( v[i] ) {
            w[i] = (uint32_T) last_used;
            last_used = last_used+1;
            first_not_used = first_not_used+1;
            while ( first_not_used<=(i+1)&& v[first_not_used])
                first_not_used = first_not_used+1;
        }
    }    
     
        
}   
              
