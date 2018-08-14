
/* CMG, Copyright (c) 2008-2010  Ioannis Koutis, Gary Miller                */
/* 
/* The CMG solver is distributed under the terms of the GNU General Public   */
/* Lincense Version 3.0 of the Free Software Foundation.                    */
/* CMG is also available under other licenses; contact authors for details. */

#include "mex.h"
#include "cmg.h"
#include "stdlib.h"

#define MAX(a,b)  (((a) > (b)) ? (a) : (b))
#define MIN(a,b)  (((a) < (b)) ? (a) : (b))

#define H_IN prhs[0]
#define b_IN prhs[1]
#define x_OUT plhs[0]


void mexFunction (int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{

    mwSize n;
    float *b , *x;
    ldl_p chol;
    mwSize i;
    mwIndex *jld;
    s_hlevel *H;
    int nlevels;
    mxArray *C;


    /* Input validation */
    if (nrhs !=2)
        mexErrMsgTxt("Wrong number of input arguments");
    
    if (!mxIsCell(H_IN))
        mexErrMsgTxt("First argument must be a Hierarchy cell");

    if (!mxIsSingle(b_IN))   
      mexErrMsgTxt("Second argument must be a non-sparse single precision vector");

    if ( (mxGetN(b_IN)!=1)) 
        mexErrMsgTxt("Second argument must be a column vector");

    n = mxGetM(b_IN);
    nlevels = (int) MAX(mxGetM(H_IN),mxGetN(H_IN));
    H = malloc(nlevels*sizeof(s_hlevel));

	if (nlevels>1) {
		C = mxGetCell(H_IN,nlevels-2);
		H[nlevels-2].islast = (boolean) *(mxGetPr(mxGetField(C,0,"islast")));
        if (H[nlevels-2].islast)
			nlevels = nlevels-1; 
	}

  
    for (i=0; i<nlevels; i++) {
        
        C = mxGetCell(H_IN,i);
       
        H[i].islast = (boolean) *(mxGetPr(mxGetField(C,0,"islast")));
        H[i].iterative = (boolean) *(mxGetPr(mxGetField(C,0,"iterative")));
        
                
        if (i<(nlevels-1)) {
        H[i].cI = (mIndex *) mxGetPr(mxGetField(C,0,"cI"));
        H[i].nc =  (mSize) *(mxGetPr(mxGetField(C,0,"nc")));
        H[i].invD = (float *) mxGetPr(mxGetField(C,0,"invD"));
        H[i].dc = (boolean) *(mxGetPr(mxGetField(C,0,"dc")));
        H[i].repeat = (int) *(mxGetPr(mxGetField(C,0,"repeat")));
        H[i].lws1 = (float *) mxGetPr(mxGetField(C,0,"lws1"));
        H[i].lws2 = (float *) mxGetPr(mxGetField(C,0,"lws2"));
        H[i].sws1 = (float *) mxGetPr(mxGetField(C,0,"sws1"));
        H[i].sws2 = (float *) mxGetPr(mxGetField(C,0,"sws2"));
        H[i].sws3 = (float *) mxGetPr(mxGetField(C,0,"sws3"));
        }

        if (i==0)
        H[i].laplacian = (boolean) *(mxGetPr(mxGetField(C,0,"laplacian")));
        
        H[i].A.a = (float *) mxGetPr(mxGetField(mxGetField(C,0,"A"),0,"a"));
        H[i].A.ia = (mIndex *) mxGetPr(mxGetField(mxGetField(C,0,"A"),0,"ia"));
        H[i].A.ja = (mIndex *) mxGetPr(mxGetField(mxGetField(C,0,"A"),0,"ja"));  
        H[i].A.n = (mSize) *(mxGetPr(mxGetField(mxGetField(C,0,"A"),0,"n")));    
        H[i].A.issym = (boolean) *(mxGetPr(mxGetField(mxGetField(C,0,"A"),0,"issym")));    

                                                                    
        if (i==(nlevels-1)){
            
            if (!H[i].iterative){
                H[i].chol.ld.a = (double *) mxGetPr(mxGetField(mxGetField(mxGetField(C,0,"chol"),0,"ld"),0,"a"));
                H[i].chol.ld.ia = (mIndex *) mxGetPr(mxGetField(mxGetField(mxGetField(C,0,"chol"),0,"ld"),0,"ia"));        
                H[i].chol.ld.ja = (mIndex *) mxGetPr(mxGetField(mxGetField(mxGetField(C,0,"chol"),0,"ld"),0,"ja"));       
                H[i].chol.ld.n = (mSize ) *(mxGetPr(mxGetField(mxGetField(mxGetField(C,0,"chol"),0,"ld"),0,"n")));   
                
                H[i].chol.ldT.a = (double *) mxGetPr(mxGetField(mxGetField(mxGetField(C,0,"chol"),0,"ldT"),0,"a"));
                H[i].chol.ldT.ia = (mIndex *) mxGetPr(mxGetField(mxGetField(mxGetField(C,0,"chol"),0,"ldT"),0,"ia"));        
                H[i].chol.ldT.ja = (mIndex *) mxGetPr(mxGetField(mxGetField(mxGetField(C,0,"chol"),0,"ldT"),0,"ja"));       
                H[i].chol.ldT.n = (mSize ) *(mxGetPr(mxGetField(mxGetField(mxGetField(C,0,"chol"),0,"ldT"),0,"n"))); 
                
                H[i].chol.p   = (mIndex *) mxGetPr(mxGetField(mxGetField(C,0,"chol"),0,"p"));
                H[i].chol.invp   = (mIndex *) mxGetPr(mxGetField(mxGetField(C,0,"chol"),0,"invp"));
            }
            else {
                H[i].invD = (float *) mxGetPr(mxGetField(C,0,"invD"));
                H[i].lws1 = (float *) mxGetPr(mxGetField(C,0,"lws1"));
                H[i].lws2 = (float *) mxGetPr(mxGetField(C,0,"lws2"));
                H[i].sws1 = (float *) mxGetPr(mxGetField(C,0,"sws1"));
                H[i].sws2 = (float *) mxGetPr(mxGetField(C,0,"sws2"));
                H[i].sws3 = (float *) mxGetPr(mxGetField(C,0,"sws3"));
            }
        
        }


 
    }

    
    b = (float *) mxGetPr(b_IN);
    x_OUT = mxCreateNumericArray(1,&n,mxSINGLE_CLASS,mxREAL);
    x = (float *) mxGetPr(x_OUT);
    
    preconditioner( H, b, 0, 1,x);
    
    free(H);
    
    
    return;

}   
              
