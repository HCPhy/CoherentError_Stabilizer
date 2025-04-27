#include "mex.h"
#include <stdint.h>
#include <string.h>
#include <stdlib.h>

/*
 * null_gf2_bitmask_mex(H)
 *
 * Compute a basis for the null space of the binary matrix H (over GF(2))
 * using a bitmask representation for fast row operations.
 */

/* Get the number of 64-bit words needed to store n bits */
static inline mwSize word_count(mwSize n) {
    return (n + 63U) >> 6; /* ceil(n/64) */
}

/* Set the bit at position col in the row */
static inline void set_bit(uint64_t *row, mwSize col) {
    row[col >> 6] |= (uint64_t)1 << (col & 63);
}

/* Check if a bit at 'col' is set in 'row' */
static inline int get_bit(const uint64_t *row, mwSize col) {
    return (row[col >> 6] >> (col & 63)) & 1U;
}

/* XOR two rows: dest ^= src (word-wise) */
static inline void xor_row(uint64_t *dest, const uint64_t *src, mwSize wc) {
    for (mwSize i = 0; i < wc; i++) {
        dest[i] ^= src[i];
    }
}

/* Swap two rows of wc words each */
static inline void swap_rows(uint64_t *A, uint64_t *B, mwSize wc) {
    for (mwSize i = 0; i < wc; i++) {
        uint64_t temp = A[i];
        A[i] = B[i];
        B[i] = temp;
    }
}

/* 
 * GF(2) RREF on a bitmask-based matrix
 * Input: 
 *   M: matrix data, m rows, n columns
 * Output:
 *   M will be transformed into RREF
 *   pivots: array of pivot column indices
 * Returns number of pivots
 */
static int gf2_rref_bitmask(uint64_t *M, mwSize m, mwSize n, int *pivots) {
    mwSize wc = word_count(n);
    mwSize pivot_row = 0;
    int pivot_count = 0;

    for (mwSize col = 0; col < n && pivot_row < m; col++) {
        /* Find a pivot in this column at or below pivot_row */
        mwSize pivot_r = m; /* sentinel */
        for (mwSize r = pivot_row; r < m; r++) {
            if (get_bit(M + r*wc, col)) {
                pivot_r = r;
                break;
            }
        }

        if (pivot_r == m) {
            /* no pivot in this column */
            continue;
        }

        /* If pivot_r != pivot_row, swap */
        if (pivot_r != pivot_row) {
            swap_rows(M + pivot_r*wc, M + pivot_row*wc, wc);
        }

        /* Record pivot column */
        pivots[pivot_count++] = (int)col;

        /* Eliminate below and above pivot_row */
        for (mwSize r = 0; r < m; r++) {
            if (r != pivot_row && get_bit(M + r*wc, col)) {
                xor_row(M + r*wc, M + pivot_row*wc, wc);
            }
        }

        pivot_row++;
    }

    return pivot_count;
}

/*
 * null_gf2_bitmask_mex:
 *
 * Input:
 *   H: m-by-n matrix over GF(2) (logical or double mod 2)
 *
 * Output:
 *   N: n-by-(num_free_vars) logical matrix representing a basis for the nullspace over GF(2)
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs != 1) {
        mexErrMsgTxt("One input required: H.");
    }
    if (nlhs > 1) {
        mexErrMsgTxt("One output required.");
    }

    const mxArray *H_in = prhs[0];
    if (!mxIsDouble(H_in) && !mxIsLogical(H_in)) {
        mexErrMsgTxt("H must be double or logical.");
    }

    mwSize m = mxGetM(H_in);
    mwSize n = mxGetN(H_in);
    mwSize wc = word_count(n);

    /* Allocate matrix M as bitmasks */
    uint64_t *M = (uint64_t *)mxCalloc(m*wc, sizeof(uint64_t));
    if (M == NULL) {
        mexErrMsgTxt("Out of memory.");
    }

    /* Convert H to bitmask form */
    if (mxIsDouble(H_in)) {
        double *H_d = (double *)mxGetData(H_in);
        for (mwSize j = 0; j < n; j++) {
            for (mwSize i = 0; i < m; i++) {
                if (((int)H_d[i + j*m]) & 1) {
                    set_bit(M + i*wc, j);
                }
            }
        }
    } else {
        bool *H_b = (bool *)mxGetData(H_in);
        for (mwSize j = 0; j < n; j++) {
            for (mwSize i = 0; i < m; i++) {
                if (H_b[i + j*m]) {
                    set_bit(M + i*wc, j);
                }
            }
        }
    }

    /* Array to hold pivot positions */
    int *pivots = (int *)mxCalloc((m < n ? m : n), sizeof(int));
    int pivot_count = gf2_rref_bitmask(M, m, n, pivots);

    /* Determine which columns are pivots */
    bool *isPivot = (bool *)mxCalloc(n, sizeof(bool));
    for (int i = 0; i < pivot_count; i++) {
        isPivot[pivots[i]] = true;
    }

    /* Count free variables */
    int free_count = 0;
    for (mwSize col = 0; col < n; col++) {
        if (!isPivot[col]) {
            free_count++;
        }
    }

    /* Create output N: n-by-free_count (logical) */
    plhs[0] = mxCreateLogicalMatrix(n, free_count);
    bool *N = (bool *)mxGetData(plhs[0]);

    /* Extract the free variables */
    int *free_vars = (int *)mxCalloc(free_count, sizeof(int));
    {
        int idx = 0;
        for (mwSize col = 0; col < n; col++) {
            if (!isPivot[col]) {
                free_vars[idx++] = (int)col;
            }
        }
    }

    /* Back-substitution:  
       For each free variable f:
         - Initialize v with v[f] = 1 and others = 0
         - Solve for pivot columns using the RREF
         
       The RREF form after gf2_rref_bitmask is such that:
       Each pivot row corresponds to one pivot. Pivot i in pivots[i].
       The pivot row is i (0-based) since we processed them in order.
    */

    /* We'll represent v as a bitmask as well for speed */
    uint64_t *v = (uint64_t *)mxCalloc(wc, sizeof(uint64_t));
    for (int fv_i = 0; fv_i < free_count; fv_i++) {
        /* Reset v to zero */
        for (mwSize w = 0; w < wc; w++) {
            v[w] = 0ULL;
        }
        /* Set the free variable bit */
        set_bit(v, free_vars[fv_i]);

        /* Back-substitute: 
           We know pivot rows: i-th pivot in row i.
           Solve from bottom pivot up:
        */
        for (int pi = pivot_count - 1; pi >= 0; pi--) {
            int c = pivots[pi]; /* pivot column */
            /* Compute the sum of known variables in that pivot row */
            /* M[row=pi, :] represents equation:
               x[c] = sum(M[pi, all other set bits]*x[that col])
             */
            uint64_t sumvec[wc];
            memcpy(sumvec, M + pi*wc, wc*sizeof(uint64_t));

            /* Clear the pivot bit itself to get the other variables */
            sumvec[c >> 6] &= ~((uint64_t)1 << (c & 63));

            /* The value of x[c] is the XOR of all x[cols] where M[pi,cols] = 1 */
            /* XOR sum of v & sumvec (AND first to select relevant bits) */
            /* We just XOR bit by bit: 
               effectively x[c] = parity of intersection of set bits in v and sumvec.
             */
            uint64_t acc = 0ULL;
            for (mwSize w = 0; w < wc; w++) {
                acc ^= (v[w] & sumvec[w]);
            }

            /* Compute parity of acc */
            /* Brian Kernighan’s trick or builtin __builtin_popcountll could be used.
               We'll do a simple popcount:
            */
            int parity = 0;
            while (acc) {
                acc &= (acc - 1); 
                parity = !parity;
            }

            /* Set v[c] = parity */
            if (parity) {
                v[c >> 6] ^= ((uint64_t)1 << (c & 63));
            } else {
                /* ensure it's cleared */
                v[c >> 6] &= ~((uint64_t)1 << (c & 63));
            }
        }

        /* Store v into N(:, fv_i) */
        /* v is a bitmask of length n, N is n-by-free_count column major */
        for (mwSize col = 0; col < n; col++) {
            N[col + fv_i*n] = get_bit(v, col);
        }
    }

    mxFree(M);
    mxFree(pivots);
    mxFree(isPivot);
    mxFree(free_vars);
    mxFree(v);
}