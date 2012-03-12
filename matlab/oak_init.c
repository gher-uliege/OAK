#include "mex.h"

void oak_init(char* fname);

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{

/*local variables */

mwSize strlen;
char *fname;
int retval;

/*convert variables */

strlen = mxGetN(prhs[0])+1;
fname = malloc(strlen * sizeof(char));
retval = mxGetString(prhs[0], fname, strlen);

  if (retval != 0) {
    mexErrMsgTxt("Convertion to string failed");
  }

/* function call */

oak_init(fname);

/*free local variables */

free(fname);

}
