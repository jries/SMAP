/*!
 * \file GPUgaussLib.h
//author Yiming Li
//date 20170301
 * \brief Prototypes for all the Cuda helper functions.
 */

// This code provides a set of functions that can be called from inside
// NVIDIA CUDA Kernels.
#ifndef GPUGAUSSLIB_H
#define GPUGAUSSLIB_H

__device__ char d_assert[1024];

__device__ float kernel_IntGauss1D(const int ii, const float x, const float sigma);

__device__ float kernel_alpha(const float z, const float Ax, const float Bx, const float d);

__device__ float kernel_dalphadz(const float z, const float Ax, const float Bx, const float d);

__device__ float kernel_d2alphadz2(const float z, const float Ax, const float Bx, const float d);

__device__ void kernel_DerivativeIntGauss1D(const int ii, const float x, const float sigma, const float N,
        const float PSFy, float *dudt, float *d2udt2);

__device__ void kernel_DerivativeIntGauss1DSigma(const int ii, const float x, 
        const float Sx, const float N, const float PSFy, float *dudt, float *d2udt2);

__device__ void kernel_DerivativeIntGauss2DSigma(const int ii, const int jj, const float x, const float y,
        const float S, const float N, const float PSFx, const float PSFy, float *dudt, float *d2udt2);

__device__ void kernel_CenterofMass2D(const int sz, const float *data, float *x, float *y);

__device__ void kernel_GaussFMaxMin2D(const int sz, const float sigma, const float * data, float *MaxN, float *MinBG);

__device__ void kernel_DerivativeGauss2D(int ii, int jj, float PSFSigma,float* theta, float *dudt, float *model);

__device__ void kernel_DerivativeGauss2D_sigma(int ii, int jj, float* theta, float *dudt, float *model);

__device__ inline void kernel_DerivativeIntGauss2Dz(const int ii, const int jj, const float *theta,
	const float PSFSigma_x, const float PSFSigma_y, const float Ax, const float Ay, 
	const float Bx, const float By, const float gamma, const float d, float *pPSFx, float *pPSFy, float *dudt, float *d2udt2,float *model);

__device__ void kernel_DerivativeGauss2D_sigmaxy(int ii, int jj,float* theta, float *dudt, float *model);

#endif