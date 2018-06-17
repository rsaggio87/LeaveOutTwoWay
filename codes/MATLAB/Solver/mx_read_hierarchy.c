/* CMG, Copyright (c) 2008-2010  Ioannis Koutis, Gary Miller                */
/* 
/* The CMG solver is distributd under the terms of the GNU General Public   */
/* Lincense Version 3.0 of the Free Software Foundation.                    */
/* CMG is also available under other licenses; contact authors for details. */

#include "mex.h"
#include "cmg.h"
#include "malloc.h"

#define MAX(a,b)  (((a) > (b)) ? (a) : (b))
#define MIN(a,b)  (((a) < (b)) ? (a) : (b))

#define NAME_IN prhs[0]
#define C_OUT plhs[0]


void mexFunction (int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{

    mwSize n,nelem;
    float *fw;
    double *bin;
    ldl_p chol;
    mwSize i;
    mwIndex *jld;
    s_hlevel *H;
    int nlevels;
    mxArray *C,*S, *S2, *UIA, *FIA;
    char *str;
    int s, nlevels;
    mIndex *iw;
    const char *shlevel_fields[] = {"laplacian", "islast","iterative","dc","nc","repeat","invD","A","cI"};
    const char *matrix_fields[] = {"a", "ia","ja","n","issym"};
    mwSize *dims;



    /* Input validation */
    if (nrhs !=1)
        mexErrMsgTxt("Wrong number of input arguments");
    
    if (!mxIsChar(NAME_IN))
        mexErrMsgTxt("Argument must be a string");

    str = mxArrayToString(prhs[i]); 
    H = allocate_hierarchy(str,&nlevels);
    s = read_hierarchy(str, & H);
    if (s == 0)
        mexErrMsgTxt("Reading hierarchy failed");
    
    C_OUT = mxCreateCellArray(1, &nlevels);

  
    for (i=0; i<(nlevels-1); i++) {
        
        *dims = 1;
        S = mxCreateStructArray(1,dims,9,shlevel_fields);
        S2 = mxCreateStructArray(1,dims,5,shlevel_fields);
        FIA = mxNumericArray(1,dims,mxClassIDFromClassName('uint32'),mxREAL)
        iw = (mIndex *) mxGetPr(UIA);
        
        
        if (i==0)
        *iw = H[i].laplacian;
        mxSetField(S,1,"laplacian",UIA);
        *iw = H[i].islast;
        mxSetField(S,1,"islast",UIA);
        *iw = H[i].iterative;        
        mxSetField(S,1,"iterative",UIA);        
        *iw = H[i].dc;        
        mxSetField(S,1,"iterative",dc);
        
        
        *iw = H[i].repeat;        
        mxSetField(S,1,"repeat",UIA);
        *iw = H[i].nc;        
        mxSetField(S,1,"nc",UIA);
        
        *dims = H[i].A.n;
        FIA = mxNumericArray(1,dims,mxClassIDFromClassName('single'),mxREAL)
        *fw = H[i].nc;
        
        
        
        C = mxGetCell(C_OUT,i);
       
        H[i].islast = (boolean) *(mxGetPr(mxGetField(C,0,"islast")));
        H[i].iterative = (boolean) *(mxGetPr(mxGetField(C,0,"iterative")));
     
     
                
        if (i<(nlevels-1)) {
        H[i].cI = (mIndex *) mxGetPr(mxGetField(C,0,"cI"));
        H[i].nc =  (mSize) *(mxGetPr(mxGetField(C,0,"nc")));
        H[i].invD = (float *) mxGetPr(mxGetField(C,0,"invD"));
        H[i].dc = (boolean) *(mxGetPr(mxGetField(C,0,"dc")));
        H[i].repeat = (int) *(mxGetPr(mxGetField(C,0,"repeat")));
        H[i].ws1 = (float *) mxGetPr(mxGetField(C,0,"ws1"));
        H[i].ws2 = (float *) mxGetPr(mxGetField(C,0,"ws2"));
        H[i].ws3 = (float *) mxGetPr(mxGetField(C,0,"ws2"));
        }

        if (i==0)
        H[i].laplacian = (boolean) *(mxGetPr(mxGetField(C,0,"laplacian")));
        
        H[i].A.a = (float *) mxGetPr(mxGetField(mxGetField(C,0,"A"),0,"a"));
        H[i].A.ia = (mIndex *) mxGetPr(mxGetField(mxGetField(C,0,"A"),0,"ia"));
        H[i].A.ja = (mIndex *) mxGetPr(mxGetField(mxGetField(C,0,"A"),0,"ja"));  
        H[i].A.n = (mSize) *(mxGetPr(mxGetField(mxGetField(C,0,"A"),0,"n")));    
        H[i].A.issym = (boolean) *(mxGetPr(mxGetField(mxGetField(C,0,"A"),0,"issym")));    

                                                                    
        if (i==(nlevels-1)){
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


 
    }


    b = (float *) mxGetPr(b_IN);
    x_OUT = mxCreateNumericArray(1,&n,mxSINGLE_CLASS,mxREAL);
    x = (float *) mxGetPr(x_OUT);
    
    preconditioner( H, b, 0, x);
    
    free(H);
    
    
    return;

}   
              
