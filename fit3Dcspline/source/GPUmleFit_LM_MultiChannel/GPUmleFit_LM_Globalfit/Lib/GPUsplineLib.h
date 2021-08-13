#ifndef GPUSPLINELIB_H
#define GPUSPLINELIB_H

__device__ void kernel_computeDelta3D(float x_delta, float y_delta, float z_delta, float *delta_f, float *delta_dxf, float *delta_dyf, float *delta_dzf);


__device__ float kernal_fAt3D(int zc, int yc, int xc, int xsize, int ysize, int zsize, float *delta_f, float *coeff);

__device__ int kernel_cholesky(float *A,int n, float *L, float*U);

__device__ void kernel_luEvaluate(float *L,float *U, float *b, int n, float *x);

__device__ void kernel_DerivativeSpline(int xc, int yc, int zc, int xsize, int ysize, int zsize, float *delta_f, float *delta_dxf, float *delta_dyf, float *delta_dzf,const float *coeff,float *theta, float*dudt,float*model);

__device__ void kernel_DerivativeSpline_multiChannel(int xc, int yc, int zc, int xsize, int ysize, int zsize, float *delta_f, float *delta_dxf, float *delta_dyf, float *delta_dzf,const float *coeff,float *theta, float*dudt,float*model);



__device__ void kernel_evalBSpline(float xi, int deg, float *bi);

__device__ void kernel_computeDelta3D_bSpline(float x_delta, float y_delta, float z_delta, float *delta_f, float *delta_dxf, float *delta_dyf, float *delta_dzf);


__device__ void kernel_Derivative_bSpline(float xi, float yi, float zi, int xsize, int nDatax, int nDatay, int nDataz, float *delta_f, float *delta_dfx, float *delta_dfy, float *delta_dfz,const float *cMat,float *theta, float*dudt,float*model);

#endif