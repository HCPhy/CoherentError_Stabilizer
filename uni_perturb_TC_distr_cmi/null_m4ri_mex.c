#include "mex.h"
#include <m4ri/m4ri.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs != 1) {
        mexErrMsgTxt("One input required.");
    }
    if (!mxIsLogical(prhs[0])) {
        mexErrMsgTxt("Input must be a logical (binary) matrix.");
    }

    size_t m = mxGetM(prhs[0]);
    size_t n = mxGetN(prhs[0]);
    bool *H_data = mxGetLogicals(prhs[0]);

    // Initialize A
    mzd_t *A = mzd_init(m, n);
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n; j++) {
            if (H_data[i + m*j]) {
                mzd_write_bit(A, i, j, 1);
            }
        }
    }

    // To find the null space (right kernel) of A:
    // Compute the left kernel of A^T
    mzd_t *AT = mzd_init(n, m);
    mzd_transpose(AT, A);

    // Compute left kernel of AT
    mzd_t *kernel_left = mzd_kernel_left(AT, 0);

    // kernel_left->nrows = dimension of the left kernel of AT
    // Each row of kernel_left is a basis vector of the left kernel of AT
    // Transpose logic:
    // If K is in the left kernel of AT, it means K * AT = 0
    // Interpreted as columns, A * K^T = 0, so K^T is in the right kernel of A.
    //
    // The returned kernel_left matrix is of dimension (k, n) for some k = dimension of kernel.
    // We actually want vectors of length n forming the null space of A.
    // The rows of kernel_left correspond to solutions K, so their transpose forms solution vectors in the null space of A.

    // Convert kernel_left rows into MATLAB array
    // kernel_left: k rows, n cols. Each row corresponds to a vector in GF(2).
    plhs[0] = mxCreateLogicalMatrix(n, kernel_left->nrows);
    bool *N_data = mxGetLogicals(plhs[0]);

    for (size_t r = 0; r < kernel_left->nrows; r++) {
        for (size_t c = 0; c < kernel_left->ncols; c++) {
            N_data[c + n*r] = mzd_read_bit(kernel_left, r, c) ? true : false;
        }
    }

    mzd_free(A);
    mzd_free(AT);
    mzd_free(kernel_left);
}