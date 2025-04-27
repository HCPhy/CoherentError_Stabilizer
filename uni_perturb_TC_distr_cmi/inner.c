#include "mex.h"
#include <stdint.h> // For int8_t, int32_t

void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{
    /* Check for proper number of arguments */
    if (nrhs != 2) {
        mexErrMsgIdAndTxt("MATLAB:inner:nargin",
            "Two input vectors required.");
    }
    if (nlhs > 1) {
        mexErrMsgIdAndTxt("MATLAB:inner:nargout",
            "Too many output arguments.");
    }

    /* Check that both inputs are int8 vectors */
    if (!mxIsInt8(prhs[0]) ||
        mxGetNumberOfDimensions(prhs[0]) > 2) {
        mexErrMsgIdAndTxt("MATLAB:inner:inputNotInt8Vector",
            "Input v must be an int8 vector.");
    }
    if (!mxIsInt8(prhs[1]) ||
        mxGetNumberOfDimensions(prhs[1]) > 2) {
        mexErrMsgIdAndTxt("MATLAB:inner:inputNotInt8Vector",
            "Input w must be an int8 vector.");
    }

    /* Get the number of elements in each input */
    size_t n_v = mxGetNumberOfElements(prhs[0]);
    size_t n_w = mxGetNumberOfElements(prhs[1]);

    /* Check that the vectors are of the same length and even */
    if (n_v != n_w) {
        mexErrMsgIdAndTxt("MATLAB:inner:lengthMismatch",
            "Vectors must be of the same length.");
    }
    if (n_v % 2 != 0) {
        mexErrMsgIdAndTxt("MATLAB:inner:oddLength",
            "Vectors must have even length.");
    }

    /* Get pointers to the data in the input vectors */
    int8_t *v = (int8_t *)mxGetData(prhs[0]);
    int8_t *w = (int8_t *)mxGetData(prhs[1]);

    /* Compute the symplectic inner product */
    int32_t t = 0;
    for (size_t i = 0; i < n_v; i += 2) {
        t += (int32_t)v[i] * (int32_t)w[i+1] + (int32_t)w[i] * (int32_t)v[i+1];
    }

    /* Compute t modulo 2 */
    int32_t t_mod = t % 2;
    if (t_mod < 0)
        t_mod += 2; // Ensure the result is non-negative

    /* Create the output scalar */
    plhs[0] = mxCreateDoubleScalar((double)t_mod);
}