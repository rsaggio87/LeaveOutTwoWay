/* ========================================================================== */
/* === CMG/Source/Hierarchy/splitforest.c       ============================= */
/* ========================================================================== */


/* CMG, Copyright (c) 2008-2010  Ioannis Koutis, Gary Miller                */
/* 
/* The CMG solver is distributd under the terms of the GNU General Public   */
/* Lincense Version 3.0 of the Free Software Foundation.                    */
/* CMG is also available under other licenses; contact authors for details. */



#include "mex.h"
#include <stdlib.h>

#define A_IN prhs[0]
#define C_OUT plhs[0]


void mexFunction (int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{

    uint32_T *C, *A;
    mwIndex i,j, k, ik, middlek, *walkbuffer,jwalk, jwalkb, jwalka, new_front;
    mwSize n, m,*indegree;
    int *ancestors, ancestors_in_path, removed_ancestors, ncut, *newancestorbuff;
    double *Icut,*Jcut, *Icutp, *Jcutp;
    bool *visited, new_ancestor_flag, walkterminated, startwalk, cut_mode;
    
    if (nrhs != 1) 
        mexErrMsgTxt("Wrong number of input arguments");
    if ((nlhs > 1))
        mexErrMsgTxt("Wrong number of output arguments");
    

    n = mxGetN(A_IN);
    m = mxGetM(A_IN);
    if ((m!=1) && (n!=1))
        mexErrMsgTxt("Input must be a vector");
    if (mxGetClassID(A_IN) != mxUINT32_CLASS)
        mexErrMsgTxt("Input vector must be of type uint32");
   
    /* Create output */
    C_OUT = mxCreateNumericMatrix(m, n, mxUINT32_CLASS, mxREAL);
    if (C_OUT == NULL) 
         mexErrMsgTxt("Unsufficient memory");
    C = (uint32_T *) mxGetPr(C_OUT);   
    A = (uint32_T *) mxGetPr(A_IN);
    if (m>n)
        n = m;
    for (j=0; j<n; j++)
        C[j]= A[j];

    
    
    ancestors = mxMalloc(n*sizeof(int));
    indegree = mxMalloc((n+2)*sizeof(mwSize)); /* +2 for handling nodes with no neighbors */
    visited = mxMalloc(n*sizeof(bool));
    walkbuffer = mxMalloc(20*sizeof(mwIndex));
    newancestorbuff = mxMalloc(20*sizeof(int));
    if ((indegree == NULL) || (walkbuffer == NULL) || (visited == NULL) || (newancestorbuff==NULL))
        mexErrMsgTxt("Unsufficient memory");
  


    
    /* compute indegrees*/
    for (j=0; j< n; j++){
        indegree[j]  =(mwSize) 0; 
        ancestors[j] = (int) 0; 
        visited[j]=false;}

    
    for (j=0; j<n; j++)
            indegree[C[j]]=indegree[C[j]]+1;


    
    ncut=0;
    /* partition into clusters of small diameter */
    for (j=0; j<n; j++) {
        jwalk=j;
        startwalk=true;
        while (startwalk && (indegree[jwalk] == 0) && !visited[jwalk] && C[j]!=(n+1) ){
            startwalk=false;
            ancestors_in_path=1;
            k=0;
            walkbuffer[k]=jwalk;
            newancestorbuff[k]=0;
            while (k<=5 || visited[jwalk]){
                jwalk = C[jwalk];
                walkterminated=(jwalk == walkbuffer[k]) || ((k>0) && (jwalk == walkbuffer[k-1]));
                if (walkterminated)
                    break; /* while */
                k = k+1;
                walkbuffer[k] = jwalk;
                if (visited[jwalk]){
                    newancestorbuff[k] = ancestors_in_path;}
                else {
                    ancestors_in_path = ancestors_in_path+1;
                    newancestorbuff[k] = ancestors_in_path;}} 
            if (k>5){ /* large diamater - make cut */
                ncut=ncut+1;
                middlek = k/2;
                C[walkbuffer[middlek]]=walkbuffer[middlek];  /* cut middle edge */
                indegree[walkbuffer[middlek+1]]=indegree[walkbuffer[middlek+1]]-1; /* update indegree */
                for (ik=(middlek+1); ik<=k; ik++){
                    ancestors[walkbuffer[ik]] = ancestors[walkbuffer[ik]]-ancestors[walkbuffer[middlek]];}
                for (ik=0; ik<=middlek; ik++) {   /* update ancestors and visited flag */
                    visited[walkbuffer[ik]]=true;
                    ancestors[walkbuffer[ik]] = ancestors[walkbuffer[ik]]+newancestorbuff[ik];}
                jwalk = walkbuffer[middlek+1]; /* set first vertex in new walk */
                startwalk=true;} /* end cut procedure */
            if (!startwalk)
                for (ik=0; ik<=k; ik++){
                    ancestors[walkbuffer[ik]] = ancestors[walkbuffer[ik]]+newancestorbuff[ik];                    
                    visited[walkbuffer[ik]]=true;}}
    } 
    
    
        
    /* partition into clusters of high conductance */
    for (j=0; j<n; j++) {
        jwalk=j;
        startwalk=true;
        while (startwalk && (indegree[jwalk] == 0) && C[j]!=(n+1)){
            startwalk=false;            
            jwalkb=jwalk;
            cut_mode = false;
            while (true){
                jwalka = C[jwalk];
                walkterminated= ((jwalka == jwalk) || (jwalka == jwalkb));
                if (walkterminated)
                    break; /* while */

                if (!cut_mode && (ancestors[jwalk] >2) && ((ancestors[jwalka]-ancestors[jwalk])>2)) { /*low conductance - make cut*/
                    C[jwalk] = jwalk; /* cut edge */
                    ncut = ncut+1;
                    indegree[jwalka] = indegree[jwalka]-1;
                    removed_ancestors = ancestors[jwalk];
                    new_front = jwalka;
                    cut_mode = true;} /* end making cut */
                jwalkb = jwalk;
                jwalk = jwalka;
                if (cut_mode)
                    ancestors[jwalk] = ancestors[jwalk] -removed_ancestors;}
            if (cut_mode){
                startwalk=true;
                jwalk=new_front;}
        } 
    } 
    
       
        mxFree(indegree);
        mxFree(visited);
        mxFree(ancestors);
        mxFree(walkbuffer);
        mxFree(newancestorbuff);
        
}

    