#include "mex.h"

void oak_done();

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{


/*convert variables */

  if (nrhs != 0) {
    mexErrMsgTxt("0 inputs required.");
  } else if(nlhs > 0) {
    mexErrMsgTxt("Too many output arguments.");
  }

/* function call */

oak_done();


}
