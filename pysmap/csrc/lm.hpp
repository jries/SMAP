// Levenberg-Marquardt maximum-likelihood fit of one ROI, for Poisson noise.
//
// Transcribed from the SMAP / GPUmleFit_LM sources (Ries lab, EMBL), which
// implement the algorithm of Smith et al., Nat. Methods 7, 373 (2010) with the
// LM damping of Li et al., Nat. Methods 15, 367 (2018).  The ten near-identical
// kernels of the original (five PSF models x two noise models) collapse into
// this one loop: the PSF model is a template parameter and only the Poisson
// noise model remains, EM excess noise being handled outside the fitter.
#pragma once

#include <cmath>
#include <cstring>

#include "linalg.hpp"

namespace smapfit {

constexpr float TOLERANCE = 1e-6f;
constexpr float INIT_LAMBDA = 0.1f;
constexpr float SCALE_UP = 10.0f;
constexpr float SCALE_DOWN = 0.1f;
constexpr float ACCEPTANCE = 1.5f;

// Accumulate the Poisson log-likelihood error, its gradient and the
// (Gauss-Newton) Hessian over all pixels of the ROI.
template <class Model>
inline void accumulate(const Model& model, const float* data, int sz,
                       const float* theta, float* err, float* jacobian,
                       float* hessian) {
    constexpr int NV = Model::NV;
    float dudt[NV], mu;

    *err = 0.0f;
    std::memset(jacobian, 0, NV * sizeof(float));
    std::memset(hessian, 0, NV * NV * sizeof(float));

    model.prepare(theta, sz);
    for (int ix = 0; ix < sz; ++ix)
        for (int iy = 0; iy < sz; ++iy) {
            model.value(ix, iy, theta, dudt, &mu);
            float d = data[iy * sz + ix];

            if (d > 0.0f)
                *err += 2.0f * ((mu - d) - d * std::log(mu / d));
            else {
                *err += 2.0f * mu;
                d = 0.0f;
            }

            const float t1 = 1.0f - d / mu;
            for (int l = 0; l < NV; ++l) jacobian[l] += t1 * dudt[l];

            const float t2 = d / (mu * mu);
            for (int l = 0; l < NV; ++l)
                for (int m = l; m < NV; ++m) {
                    hessian[l * NV + m] += t2 * dudt[l] * dudt[m];
                    hessian[m * NV + l] = hessian[l * NV + m];
                }
        }
}

// Fisher information -> Cramer-Rao lower bounds, plus the log-likelihood.
template <class Model>
inline void crlb_and_logl(const Model& model, const float* data, int sz,
                          const float* theta, float* crlb, float* logl) {
    constexpr int NV = Model::NV;
    float M[NV * NV] = {0}, Minv[NV * NV] = {0}, dudt[NV], mu;
    float div = 0.0f;

    model.prepare(theta, sz);
    for (int ix = 0; ix < sz; ++ix)
        for (int iy = 0; iy < sz; ++iy) {
            model.value(ix, iy, theta, dudt, &mu);
            const float d = data[iy * sz + ix];

            for (int k = 0; k < NV; ++k)
                for (int l = k; l < NV; ++l) {
                    M[k * NV + l] += dudt[l] * dudt[k] / mu;
                    M[l * NV + k] = M[k * NV + l];
                }

            if (mu > 0.0f) {
                if (d > 0.0f)
                    div += d * std::log(mu) - mu - d * std::log(d) + d;
                else
                    div += -mu;
            }
        }

    mat_inv_n(M, Minv, crlb, NV);
    *logl = div;
}

// Fit one ROI.  `theta` and `crlb` must have room for Model::NV values.
template <class Model>
inline void lm_fit(const Model& model, const float* data, int sz, int iterations,
                   float* theta, float* crlb, float* logl, int* used_iterations) {
    constexpr int NV = Model::NV;

    float old_theta[NV], maxjump[NV];
    float new_update[NV], old_update[NV];
    float jacobian[NV] = {0}, hessian[NV * NV] = {0};
    float L[NV * NV] = {0}, U[NV * NV] = {0};

    for (int i = 0; i < NV; ++i) {
        new_update[i] = 1e13f;
        old_update[i] = 1e13f;
        maxjump[i] = 1.0f;
    }

    model.init(data, sz, theta, maxjump);
    for (int i = 0; i < NV; ++i) old_theta[i] = theta[i];

    float new_lambda = INIT_LAMBDA, old_lambda = INIT_LAMBDA, mu = 1.0f;
    float new_err = 0.0f, old_err = 1e13f;
    int err_flag = 0;

    accumulate(model, data, sz, theta, &new_err, jacobian, hessian);

    int k = 0;
    for (; k < iterations; ++k) {
        if (std::fabs((new_err - old_err) / new_err) < TOLERANCE) break;  // converged

        if (new_err > ACCEPTANCE * old_err) {
            // the step made things much worse: go back and damp harder
            for (int i = 0; i < NV; ++i) {
                theta[i] = old_theta[i];
                new_update[i] = old_update[i];
            }
            new_lambda = old_lambda;
            new_err = old_err;
            mu = std::max((1 + new_lambda * SCALE_UP) / (1 + new_lambda), 1.3f);
            new_lambda *= SCALE_UP;
        } else if (new_err < old_err && err_flag == 0) {
            new_lambda *= SCALE_DOWN;
            mu = 1 + new_lambda;
        }

        for (int i = 0; i < NV; ++i) hessian[i * NV + i] *= mu;

        std::memset(L, 0, NV * NV * sizeof(float));
        std::memset(U, 0, NV * NV * sizeof(float));
        err_flag = cholesky(hessian, NV, L, U);

        if (err_flag != 0) {
            mu = std::max((1 + new_lambda * SCALE_UP) / (1 + new_lambda), 1.3f);
            new_lambda *= SCALE_UP;
            continue;
        }

        for (int i = 0; i < NV; ++i) {
            old_theta[i] = theta[i];
            old_update[i] = new_update[i];
        }
        old_lambda = new_lambda;
        old_err = new_err;

        lu_evaluate(L, U, jacobian, NV, new_update);

        for (int i = 0; i < NV; ++i) {
            // a reversal of direction means we overshot: halve the step limit
            if (new_update[i] / old_update[i] < -0.5f) maxjump[i] *= 0.5f;
            new_update[i] = new_update[i] / (1 + std::fabs(new_update[i] / maxjump[i]));
            theta[i] -= new_update[i];
        }
        model.clamp(theta, sz);

        accumulate(model, data, sz, theta, &new_err, jacobian, hessian);
    }

    *used_iterations = k;
    crlb_and_logl(model, data, sz, theta, crlb, logl);
}

}  // namespace smapfit
