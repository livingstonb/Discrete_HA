/*==========================================================
 * draw_random.c
 *
 * Draws random numbers from a probability distribution
 *
 * The calling syntax is:
 *
 *		outMatrix = draw_random(cumdistVec, nOut,mOut)
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2007-2012 The MathWorks, Inc.
 *
 *========================================================*/

#include "mex.h"

// recursively find index for which target < Vec[i] and target >= Vec[i-1]
int search( double *Vec, double target, int start, int end)
{
    int mid;
            
    if ( start == end )
        return start;
    else
        mid = (end + start) / 2;
    
    if ( target > Vec[mid] )
        return search ( Vec, target, mid, end);
    else if ( mid == 0 )
        return mid;
    else if ( target < Vec[mid-1] )
        return search ( Vec, target, start, mid);
    else
        return mid;
}
 
// main computational routine
void draw_random( double *inVec, double *cumdistVec, double *outVec)
{
    int n, ndist;
    
    n = (int) sizeof(inVec) / sizeof(double);
    ndist = (int) sizeof(cumdistVec) / sizeof(double);
    
    for (i=0, i<n; i++) {
        /* find output index of element of cumdistVec */
        outVec[i] = search(cumdistVec, invVec[i], 0, ndist-1) + 1;
    }
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double *inVec;                  /* input 1xn input matrix */
    double *cumdistVec;             /* 1xndist input matrix */
    size_t ncols;                   /* size of matrix */
    double *outVec   ;              /* output matrix */

    /* check for proper number of arguments */
    if(nrhs!=2) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Two inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","One output required.");
    }
    /* make sure the first input argument is scalar */
    if( !mxIsDouble(prhs[0]) || 
         mxIsComplex(prhs[0]) ) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar","Input vector 1 must be of type double.");
    }
    
    /* make sure the second input argument is type double */
    if( !mxIsDouble(prhs[1]) || 
         mxIsComplex(prhs[1])) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble","Input vector 2 must be type double.");
    }
    
    /* check that number of columns in input arguments are 1 */
    if ( mxGetM(prhs[0])!=1 ) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Input 1 must be row vector.");
    }
    
    if ( mxGetM(phrs[1]!=1 ) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Input 2 must be row vector.");
    }
    
    /* create pointers to the real data in the input matrices  */
    inVec = mxGetPr(prhs[0]);
    cumdistVec = mxGetPr(prhs[1]);
    
    // get sizes
    ncols_inVec = mxGetN(prhs[0]);
    ncols_cumdistVec = mxgetN(phrs[0]);

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(1, (mwSize) ncols_inVec, mxREAL);

    /* get a pointer to the real data in the output matrix */
    outVec = mxGetPr(plhs[0]);

    /* call the computational routine */
    get_random(inVec, cumdistVec, outVec);
}
