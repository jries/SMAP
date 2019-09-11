/*!
 * \file CPUgaussLib.h
//author Yiming Li
//date 20170301
 */


#ifndef CPUGAUSSLIB_H
#define CPUGAUSSLIB_H

 float erf(float x);
float erfc(float x);

 float kernel_IntGauss1D(const int ii, const float x, const float sigma);

 float kernel_alpha(const float z, const float Ax, const float Bx, const float d);

 float kernel_dalphadz(const float z, const float Ax, const float Bx, const float d);

 float kernel_d2alphadz2(const float z, const float Ax, const float Bx, const float d);

 void kernel_DerivativeIntGauss1D(const int ii, const float x, const float sigma, const float N,
        const float PSFy, float *dudt, float *d2udt2);

 void kernel_DerivativeIntGauss1DSigma(const int ii, const float x, 
        const float Sx, const float N, const float PSFy, float *dudt, float *d2udt2);

 void kernel_DerivativeIntGauss2DSigma(const int ii, const int jj, const float x, const float y,
        const float S, const float N, const float PSFx, const float PSFy, float *dudt, float *d2udt2);

 void kernel_CenterofMass2D(const int sz, const float *data, float *x, float *y);

 void kernel_GaussFMaxMin2D(const int sz, const float sigma, const float * data, float *MaxN, float *MinBG);

 void kernel_DerivativeGauss2D(int ii, int jj, float PSFSigma,float* theta, float *dudt, float *model);

 void kernel_DerivativeGauss2D_sigma(int ii, int jj, float* theta, float *dudt, float *model);

 void kernel_DerivativeIntGauss2Dz(const int ii, const int jj, const float *theta,
	const float PSFSigma_x, const float PSFSigma_y, const float Ax, const float Ay, 
	const float Bx, const float By, const float gamma, const float d, float *pPSFx, float *pPSFy, float *dudt, float *d2udt2,float *model);

 void kernel_DerivativeGauss2D_sigmaxy(int ii, int jj,float* theta, float *dudt, float *model);

#endif