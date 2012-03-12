#include "mex.h"

void oak_assim(int ntime,int n, int r,double* Sf,double* Sa);

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{

/*local variables */

  int ntime;
  double* Ef;
  double* Ea;
  int n,r;

/*convert variables */

  /* Check for proper number of arguments. */
  if (nrhs != 2) {
    mexErrMsgTxt("2 inputs required.");
  } else if(nlhs > 2) {
    mexErrMsgTxt("Too many output arguments.");
  }

  if (mxGetNumberOfElements(prhs[0]) != 1) {
    mexErrMsgTxt("First parameter should be a scalar.");    
  }
  ntime = mxGetScalar(prhs[0]);

  n = mxGetM(prhs[1]);
  r = mxGetN(prhs[1]);

  Ef = mxGetPr(prhs[1]);

  plhs[0] = mxCreateDoubleMatrix(n,r, mxREAL);

  Ea = mxGetPr(plhs[0]);

  mexPrintf("Hello, world! %d %d %d\n",ntime,n,r);

/* function call */

  oak_assim(ntime,n,r,Ef,Ea);

/*free local variables */



}
