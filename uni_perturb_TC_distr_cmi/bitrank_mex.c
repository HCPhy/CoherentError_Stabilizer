#include "mex.h"
#include "matrix.h"
#include <stdint.h>
#include <string.h>

/**
 * Function: bitrank_packed
 * ------------------------
 * Computes the rank of a binary matrix using Gaussian elimination over GF(2) with bit-packing.
 * This version does not use AVX2 intrinsics and relies on standard C operations.
 *
 * Parameters:
 *   mat_packed    - Pointer to the bit-packed binary matrix.
 *   n_row         - Number of rows in the matrix.
 *   n_col         - Number of columns in the matrix.
 *   words_per_row - Number of 64-bit words per row.
 *   rk            - Pointer to store the computed rank.
 */
void bitrank_packed(uint64_t *mat_packed, mwSize n_row, mwSize n_col, mwSize words_per_row, double *rk) {
    *rk = 0.0;
    mwSize col = 0;
    mwSize bits_per_word = 64;

    for (mwSize i = 0; i < n_row; i++) {
        // Declare and initialize pivot_row before the while loop
        mwSize pivot_row = i;

        // Find the next pivot column
        while (col < n_col) {
            mwSize word_idx = col / bits_per_word;
            mwSize bit_idx = col % bits_per_word;
            uint64_t mask = ((uint64_t)1) << bit_idx;

            // Search for a pivot in the current column
            int found = 0;
            pivot_row = i; // Initialize pivot_row for this column
            for (mwSize j = i; j < n_row; j++) {
                if (mat_packed[j * words_per_row + word_idx] & mask) {
                    pivot_row = j;
                    found = 1;
                    break;
                }
            }

            if (found) {
                break;
            }
            col++;
        }

        if (col >= n_col) {
            break; // No more pivot columns
        }

        // Swap rows if necessary
        if (pivot_row != i) {
            for (mwSize w = 0; w < words_per_row; w++) {
                uint64_t temp = mat_packed[i * words_per_row + w];
                mat_packed[i * words_per_row + w] = mat_packed[pivot_row * words_per_row + w];
                mat_packed[pivot_row * words_per_row + w] = temp;
            }
        }

        (*rk)++; // Increment rank

        mwSize word_idx = col / bits_per_word;
        mwSize bit_idx = col % bits_per_word;
        uint64_t mask = ((uint64_t)1) << bit_idx;

        // Eliminate all other rows with a 1 in the pivot column
        for (mwSize j = 0; j < n_row; j++) {
            if (j != i && (mat_packed[j * words_per_row + word_idx] & mask)) {
                // Perform XOR operation across the entire row
                for (mwSize w = 0; w < words_per_row; w++) {
                    mat_packed[j * words_per_row + w] ^= mat_packed[i * words_per_row + w];
                }
            }
        }

        col++; // Move to the next column
    }
}

/**
 * MEX Gateway Function
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Input Validation
    if (nrhs != 1) {
        mexErrMsgIdAndTxt("bitrank_mex:nargin", "One input required.");
    }
    if (!mxIsLogical(prhs[0])) {
        mexErrMsgIdAndTxt("bitrank_mex:inputNotLogical", "Input must be a logical matrix.");
    }

    // Retrieve Matrix Dimensions
    mwSize n_row = mxGetM(prhs[0]);
    mwSize n_col = mxGetN(prhs[0]);

    // Calculate the number of 64-bit words per row
    mwSize bits_per_word = 64;
    mwSize words_per_row = (n_col + bits_per_word - 1) / bits_per_word;

    // Access Input Matrix Data
    const unsigned char *input_mat = (const unsigned char*)mxGetLogicals(prhs[0]);

    // Allocate Bit-Packed Matrix
    uint64_t *mat_packed = (uint64_t*)mxCalloc(n_row * words_per_row, sizeof(uint64_t));
    if (mat_packed == NULL) {
        mexErrMsgIdAndTxt("bitrank_mex:memory", "Unable to allocate memory for packed matrix.");
    }

    // Pack the input matrix into mat_packed
    for (mwSize row = 0; row < n_row; row++) {
        for (mwSize col_idx = 0; col_idx < n_col; col_idx++) {
            if (input_mat[row + col_idx * n_row]) { // MATLAB is column-major
                mwSize word_idx = col_idx / bits_per_word;
                mwSize bit_idx = col_idx % bits_per_word;
                mat_packed[row * words_per_row + word_idx] |= ((uint64_t)1) << bit_idx;
            }
        }
    }

    // Create Output for Rank
    plhs[0] = mxCreateDoubleScalar(0);
    double *rk = mxGetPr(plhs[0]);

    // Compute the Rank
    bitrank_packed(mat_packed, n_row, n_col, words_per_row, rk);

    // Free Allocated Memory
    mxFree(mat_packed);
}