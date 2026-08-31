// Small dense linear algebra for the Levenberg-Marquardt fitter.
//
// Transcribed from the SMAP / GPUmleFit_LM sources (Ries lab, EMBL; original
// template by Keith Lidke) so that the numerical behaviour is preserved:
//   MatInvLib.h        -> mat_inv_n
//   CPUsplineLib.cpp   -> cholesky, lu_evaluate
#pragma once

#include <cmath>
#include <cstring>

namespace smappy {

// Cholesky-like decomposition of a symmetric matrix A (n x n, row major).
// Returns 1 if A is not positive definite, in which case L/U are unusable.
inline int cholesky(const float* A, int n, float* L, float* U) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < (i + 1); ++j) {
            float s = 0.0f;
            for (int k = 0; k < j; ++k) s += U[i * n + k] * U[j * n + k];
            if (i == j) {
                if (A[i * n + i] - s < 0.0f) return 1;
                U[i * n + j] = std::sqrt(A[i * n + i] - s);
                L[j * n + i] = U[i * n + j];
            } else {
                U[i * n + j] = (1.0f / U[j * n + j] * (A[i * n + j] - s));
                L[j * n + i] = U[i * n + j];
            }
        }
    }
    return 0;
}

// Solve A x = b given the decomposition from cholesky().
inline void lu_evaluate(const float* L, const float* U, const float* b, int n,
                        float* x) {
    float y[8] = {0};
    for (int i = 0; i < n; ++i) {
        y[i] = b[i];
        for (int j = 0; j < i; ++j) y[i] -= L[j * n + i] * y[j];
        y[i] /= L[i * n + i];
    }
    for (int i = n - 1; i >= 0; --i) {
        x[i] = y[i];
        for (int j = i + 1; j < n; ++j) x[i] -= U[j * n + i] * x[j];
        x[i] /= U[i * n + i];
    }
}

// n x n matrix inversion; only the diagonal of the inverse is returned.
// M is overwritten. Used for the Cramer-Rao lower bounds.
inline void mat_inv_n(float* M, float* Minv, float* diag, int n) {
    float tmp1 = 0.0f;
    float yy[8];

    for (int jj = 0; jj < n; ++jj) {
        for (int ii = 0; ii <= jj; ++ii) {
            if (ii > 0) {
                for (int kk = 0; kk <= ii - 1; ++kk)
                    tmp1 += M[ii + kk * n] * M[kk + jj * n];
                M[ii + jj * n] -= tmp1;
                tmp1 = 0.0f;
            }
        }
        for (int ii = jj + 1; ii < n; ++ii) {
            if (jj > 0) {
                for (int kk = 0; kk <= jj - 1; ++kk)
                    tmp1 += M[ii + kk * n] * M[kk + jj * n];
                M[ii + jj * n] = (1.0f / M[jj + jj * n]) * (M[ii + jj * n] - tmp1);
                tmp1 = 0.0f;
            } else {
                M[ii + jj * n] = (1.0f / M[jj + jj * n]) * M[ii + jj * n];
            }
        }
    }

    tmp1 = 0.0f;
    for (int num = 0; num < n; ++num) {
        yy[0] = (num == 0) ? 1.0f : 0.0f;
        for (int ii = 1; ii < n; ++ii) {
            float b = (ii == num) ? 1.0f : 0.0f;
            for (int jj = 0; jj <= ii - 1; ++jj) tmp1 += M[ii + jj * n] * yy[jj];
            yy[ii] = b - tmp1;
            tmp1 = 0.0f;
        }
        Minv[n - 1 + num * n] = yy[n - 1] / M[(n - 1) + (n - 1) * n];
        for (int ii = n - 2; ii >= 0; --ii) {
            for (int jj = ii + 1; jj < n; ++jj)
                tmp1 += M[ii + jj * n] * Minv[jj + num * n];
            Minv[ii + num * n] = (1.0f / M[ii + ii * n]) * (yy[ii] - tmp1);
            tmp1 = 0.0f;
        }
    }

    if (diag)
        for (int ii = 0; ii < n; ++ii) diag[ii] = Minv[ii * n + ii];
}

}  // namespace smappy
