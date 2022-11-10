/*==========================================================
 * dtwm_costmatrix.c - example in MATLAB External Interfaces
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2007-2012 The MathWorks, Inc.
 *
 *========================================================*/

#include "mex.h"
#include "math.h"
#include "matrix.h"
#include "limits.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>		// mkt: for memcpy()

#define EM(ii,jj)((jj)*m+(ii))

/*
 * Auxiliary function: return the arg min, ignoring NANs, -1 if all NANs
 * TODO: remove isnan and explain, check, time
*/
static inline
int argmin(const double *list, int n)
{
    int ii=-1;
    double vv=INFINITY;
    for(int i=0; i<n; i++) {
        /* The following is a faster equivalent to
         *    if(!isnan(list[i]) && list[i]<vv)
         * because   (NAN < x) is false
         */
        if(list[i]<vv) {
            ii=i;
            vv=list[i];
        }
    }
    return ii;
}    

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    double *cm_in, *cm_out, *lm, *clist, *sc;
    int *pn, *di, *dj, *sm_in, *sm_out;
    int m,n,ee;
    bool *wm;
    
    double nan;
    nan = mxGetNaN();
    
    m = mxGetM(prhs[0]); /*rows*/
    n = mxGetN(prhs[0]); /*columns*/
    
    /*mexPrintf("\nrows: %i, columns: %i\n", m,n);*/

//mkt: mxGetPr deprecated (and generated errors under Xcode/Clang) 
	cm_in = (double *)mxGetData(prhs[0]); 	// cm_in = mxGetPr(prhs[0]);
	pn = (int *)mxGetData(prhs[1]);			// pn = mxGetPr(prhs[1]);
	di = (int *)mxGetData(prhs[2]);			// di = mxGetPr(prhs[2]);
	dj = (int *)mxGetData(prhs[3]);			// dj = mxGetPr(prhs[3]);
	sc = (double *)mxGetData(prhs[4]);		// sc = mxGetPr(prhs[4]);    
	lm = (double *)mxGetData(prhs[5]);		// lm = mxGetPr(prhs[5]);
	clist = (double *)mxGetData(prhs[6]);	// double = mxGetPr(prhs[6]);
	sm_in = (int *)mxGetData(prhs[7]);		// sm_in = mxGetPr(prhs[7]);
	wm = (bool *)mxGetData(prhs[8]);		// wm = mxGetPr(prhs[8]);
                
    int nsteps = mxGetNumberOfElements(prhs[1]);
    int npats=pn[nsteps-1]+1;
    
    ee = mxGetNumberOfElements(prhs[0]);
        
    plhs[0] = mxCreateDoubleMatrix(m,n,mxREAL);
    plhs[1] = mxCreateNumericMatrix(m,n,mxINT32_CLASS,mxREAL);
    
    cm_out = (double *)mxGetData(plhs[0]);	// cm_out = mxGetPr(plhs[0]);
    sm_out = (int *)mxGetData(plhs[1]);		// sm_out = mxGetPr(plhs[1]);
    
    memcpy(cm_out,cm_in,sizeof(double)*ee);
    memcpy(sm_out,sm_in,sizeof(int)*ee);
 
    int ii, jj, p;
    double cc;
    
    /*j: columns, i: rows
     */
    
    for (int j=0; j<n; j++) {
        for (int i=0; i<m; i++) {
           
            /* out of window? */
            if(!wm[EM(i,j)]){
                continue;
            } 
            
            for(int z=0; z<npats; z++) {
                clist[z]=nan;             
            }

            for(int s=0; s<nsteps; s++) {
                p = pn[s];
                ii = i-di[s];
                jj = j-dj[s];
                
                if(ii>=0 && jj>=0) {
                    cc = sc[s];
                    if(cc==-1.0) {
                        clist[p]=cm_out[EM(ii,jj)];
                    } else {
                        clist[p]+=cc*lm[EM(ii,jj)];
                    }
                }
            }

            int minc=argmin(clist,npats);
            if(minc>-1) {
                cm_out[EM(i,j)]=clist[minc];
                sm_out[EM(i,j)]=minc+1;
            }

        }
    }
    
    /*mxDestroyArray(*clist);*/
}
