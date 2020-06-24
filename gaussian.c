/*--------------------------- gaussian.c ----------------------------
 *
 * This file contains forward and backward substitution functions of
 * Guassian elimination solving
 *
 *      A x = B
 *
 * No specific form of A and B is assumed. Operations on A and B are
 * performed simultaneously.
 -------------------------------------------------------------------*/
#include "galois.h"

// perform forward substitution on a matrix to transform it to a upper triangular structure
//static long long forward_substitute(int nrow, int ncolA, int ncolB, GF_ELEMENT A[][ncolA], GF_ELEMENT B[][ncolB])
long long forward_substitute(int nrow, int ncolA, int ncolB, GF_ELEMENT **A, GF_ELEMENT **B)
{
    //printf("entering foward_substitute()...\n");
    long operations = 0;
    int i, j, k, m, n, p;
    int pivot;
    GF_ELEMENT quotient;

    // transform A into upper triangular structure by row operation
    int boundary = nrow >= ncolA ? ncolA : nrow;

    int has_a_dimension;
    for (i=0; i<boundary; i++) {
        has_a_dimension = 1;            // whether this column is all-zero

        if (A[i][i] == 0) {
            has_a_dimension = 0;
            /* Look for nonzero element below diagonal */
            for (pivot=i+1; pivot<nrow; pivot++) {
                if (A[pivot][i] != 0) {
                    has_a_dimension = 1;
                    break;
                }
            }
            // if this column is an zero column, skip this column
            if (!has_a_dimension)
                continue;
            else {
                // swap row
                GF_ELEMENT tmp2;
                for (m=0; m<ncolA; m++) {
                    tmp2 = A[i][m];
                    A[i][m] = A[pivot][m];
                    A[pivot][m] = tmp2;
                }
                // swap B accordingly
                // rows of B are in block memory, so only exchanging pointers
                GF_ELEMENT *temp_p;
                temp_p = B[i];
                B[i] = B[pivot];
                B[pivot] = temp_p;
            }
        }
        // Eliminate nonzero elements beow diagonal
        for (j=i+1; j<nrow; j++) {
            if (A[j][i] == 0)
                continue;   // skip zeros
            quotient = galois_divide(A[j][i], A[i][i]);
            operations += 1;
            // eliminate the items under row i at col i
            galois_multiply_add_region(&(A[j][i]), &(A[i][i]), quotient, ncolA-i);
            operations += (ncolA-i);
            // simultaneously do the same thing on right matrix B
            galois_multiply_add_region(B[j], B[i], quotient, ncolB);
            operations += ncolB;
        }
    }
    return operations;
}

// perform back-substitution on full-rank upper trianguler matrix A
long long back_substitute(int nrow, int ncolA, int ncolB, GF_ELEMENT *A[], GF_ELEMENT *B[])
{
    //printf("entering back_substitute()...\n");
    long operations = 0;

    // Transform the upper triangular matrix A into diagonal.
    int i, j, k, l;
    for (i=ncolA-1; i>=0; i--) {
        // eliminate all items above A[i][i]
        for (j=0; j<i; j++) {
            if (A[j][i] == 0)
                continue;       // skip zeros
            GF_ELEMENT quotient = galois_divide(A[j][i], A[i][i]);
            operations += 1;
            A[j][i] = 0;
            // doing accordingly to B
            galois_multiply_add_region(B[j], B[i], quotient, ncolB);
            operations += ncolB;
        }
        // diagonalize diagonal element
        if (A[i][i] != 1) {
            galois_multiply_region(B[i], galois_divide(1, A[i][i]), ncolB);
            operations += ncolB;
            A[i][i] = 1;
        }

    }
    return operations;
}
