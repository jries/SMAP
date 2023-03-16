#ifndef CPUSPLINELIB_H
#define CPUSPLINELIB_H

 void kernel_computeDelta3D(float x_delta, float y_delta, float z_delta, float *delta_f, float *delta_dxf, float *delta_dyf, float *delta_dzf);

 int kernel_cholesky(float *A,int n, float *L, float*U);

 void kernel_luEvaluate(float *L,float *U, float *b, int n, float *x);


 void kernel_DerivativeSpline(int xc, int yc, int zc, int xsize, int ysize, int zsize, float *delta_f, float *delta_dxf, float *delta_dyf, float *delta_dzf,const float *coeff,float *theta, float*dudt,float*model);


#endif