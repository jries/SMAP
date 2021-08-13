#include "GPUsplineLib.h"
//#include <cuda_runtime.h>
#include "definitions.h"
//#define pi 3.141592f



//**************************************************************************************************************************************************
// This function for calculation of the common term for Cspline is adpopted from
//"Analyzing Single Molecule Localization Microscopy Data Using Cubic Splines", Hazen Babcok, Xiaowei Zhuang,Scientific Report, 1, 552 , 2017.
__device__ inline void kernel_computeDelta3D(float x_delta, float y_delta, float z_delta, float *delta_f, float *delta_dxf, float *delta_dyf, float *delta_dzf) {
    
	int i,j,k;
	float cx,cy,cz;

	cz = 1.0;
	for(i=0;i<4;i++){
		cy = 1.0;
		for(j=0;j<4;j++){
			cx = 1.0;
			for(k=0;k<4;k++){
				delta_f[i*16+j*4+k] = cz * cy * cx;
				if(k<3){
					delta_dxf[i*16+j*4+k+1] = ((float)k+1) * cz * cy * cx;
				}
				
				if(j<3){
					delta_dyf[i*16+(j+1)*4+k] = ((float)j+1) * cz * cy * cx;
				}
				
				if(i<3){
					delta_dzf[(i+1)*16+j*4+k] = ((float)i+1) * cz * cy * cx;
				}
				
				cx = cx * x_delta;
			}
			cy = cy * y_delta;
		}
		cz= cz * z_delta;
	}
}


//***********************************************************************************************************
__device__ inline int kernel_cholesky(float *A,int n, float *L, float*U) {
	int info = 0;
	for (int i = 0; i < n; i++)
		for (int j = 0; j < (i+1); j++) {
			float s = 0;
			for (int k = 0; k < j; k++)
				s += U[i * n + k] * U[j * n + k];

			if (i==j){
				if (A[i*n+i]-s>=0){
					U[i * n + j] = sqrt(A[i * n + i] - s);
					L[j*n+i]=U[i * n + j];
				}
				else{
					info =1;
					return info;
				}
			}
			else{
				U[i * n + j] = (1.0 / U[j * n + j] * (A[i * n + j] - s));
				L[j*n+i]=U[i * n + j];
			}

		}
	return info;
}
//******************************************************************************************************
__device__ inline void kernel_luEvaluate(float *L,float *U, float *b, const int n, float *x) {
	//Ax = b -> LUx = b. Then y is defined to be Ux
	//for sigmaxy, we have 6 parameters
	float y[5*Max_No_Channel] = {0};
	int i = 0;
	int j = 0;
	// Forward solve Ly = b
	for (i = 0; i < n; i++)
	{
		y[i] = b[i];
		for (j = 0; j < i; j++)
		{
			y[i] -= L[j*n+i] * y[j];
		}
		y[i] /= L[i*n+i];
	}
	// Backward solve Ux = y
	for (i = n - 1; i >= 0; i--)
	{
		x[i] = y[i];
		for (j = i + 1; j < n; j++)
		{
			x[i] -= U[j*n+i] * x[j];
		}
		x[i] /= U[i*n+i];
	}

}


//**************************************************************************************************************

__device__ inline void kernel_DerivativeSpline(int xc, int yc, int zc, int xsize, int ysize, int zsize, float *delta_f, float *delta_dxf, float *delta_dyf, float *delta_dzf,const float *coeff,float *theta, float*dudt,float*model) {
	int i;
	float temp =0;
	//float dudt_temp[NV_PSP] = {0};//,temp;
	memset(dudt,0,NV_PSP*sizeof(float));
	//for (i=0;i<NV_PSP;i++) dudt[i]=0;
	
	xc = max(xc,0);
	xc = min(xc,xsize-1);

	yc = max(yc,0);
	yc = min(yc,ysize-1);

	zc = max(zc,0);
	zc = min(zc,zsize-1);
	
	

	for (i=0;i<64;i++){		
		temp+=delta_f[i]*coeff[i*(xsize*ysize*zsize)+zc*(xsize*ysize)+yc*xsize+xc];
		dudt[0]+=delta_dxf[i]*coeff[i*(xsize*ysize*zsize)+zc*(xsize*ysize)+yc*xsize+xc];
		dudt[1]+=delta_dyf[i]*coeff[i*(xsize*ysize*zsize)+zc*(xsize*ysize)+yc*xsize+xc];
		dudt[4]+=delta_dzf[i]*coeff[i*(xsize*ysize*zsize)+zc*(xsize*ysize)+yc*xsize+xc];
				//temp = tex1Dfetch(tex_test, i*(xsize*ysize*zsize)+zc*(xsize*ysize)+yc*xsize+xc);
				//pd+=delta_f[i]*temp;
	}
	dudt[0]*=-1.0f*theta[2];
	dudt[1]*=-1.0f*theta[2];
	dudt[4]*=theta[2];
	dudt[2]=temp;
	dudt[3]=1.0f;
	*model = theta[3]+theta[2]*temp;
	
	//return pd;
}


//**************************************************************************************************************

__device__ inline void kernel_DerivativeSpline_multiChannel(int xc, int yc, int zc, int xsize, int ysize, int zsize, float *delta_f, float *delta_dxf, float *delta_dyf, float *delta_dzf,const float *coeff,float *theta, float*dudt,float*model) {
	int i;
	float temp =0;
	//float dudt_temp[NV_PSP] = {0};//,temp;
	memset(dudt,0,5*sizeof(float));
	//for (i=0;i<NV_PSP;i++) dudt[i]=0;
	
	xc = max(xc,0);
	xc = min(xc,xsize-1);

	yc = max(yc,0);
	yc = min(yc,ysize-1);

	zc = max(zc,0);
	zc = min(zc,zsize-1);
	
	

	for (i=0;i<64;i++){		
		temp+=delta_f[i]*coeff[i*(xsize*ysize*zsize)+zc*(xsize*ysize)+yc*xsize+xc];
		dudt[0]+=delta_dxf[i]*coeff[i*(xsize*ysize*zsize)+zc*(xsize*ysize)+yc*xsize+xc];
		dudt[1]+=delta_dyf[i]*coeff[i*(xsize*ysize*zsize)+zc*(xsize*ysize)+yc*xsize+xc];
		dudt[2]+=delta_dzf[i]*coeff[i*(xsize*ysize*zsize)+zc*(xsize*ysize)+yc*xsize+xc];
				//temp = tex1Dfetch(tex_test, i*(xsize*ysize*zsize)+zc*(xsize*ysize)+yc*xsize+xc);
				//pd+=delta_f[i]*temp;
	}
	dudt[0]*=-1.0f*theta[3];
	dudt[1]*=-1.0f*theta[3];
	dudt[2]*=theta[3];
	dudt[3]=temp;
	dudt[4]=1.0f;
	*model = theta[4]+theta[3]*temp;
	
	//return pd;
}


//**************************************************************************************************************
__device__ inline void kernel_evalBSpline(float xi, int deg, float *bi) {
	float x2,x3;



	//switch(deg){
	//case 2:
	//	if (xi>=1.0f/2.0f&&xi<3.0f/2.0f){
	//		*bi = 9.0f/8.0f-3.0f/2.0f*xi+1.0f/2.0f*x2;
	//	}
	//	if (xi>=-1.0f/2.0f&&xi<1.0f/2.0f){
	//		*bi = 3.0f/4.0f-x2;
	//	}
	//	if (xi>=-3.0f/2.0f&&xi<-1.0f/2.0f){
	//		*bi = 9.0f/8.0f+3.0f/2.0f*xi+1.0f/2.0f*x2;
	//	}
	//case 3:
	//	if (xi>=1.0f&&xi<2.0f){
	//		*bi = 4.0f/3.0f-2.0f*xi+x2-1.0f/6.0f*x3;
	//	}
	//	if (xi>=0.0f&&xi<1.0f){
	//		*bi = 2.0f/3.0f-x2+1.0f/2.0f*x3;
	//	}
	//	if (xi>=-1.0f&&xi<0.0f){
	//		*bi = 2.0f/3.0f-x2-1.0f/2.0f*x3;
	//	}
	//	if (xi>=-2.0f&&xi<-1.0f){
	//		*bi = 4.0f/3.0f+2.0f*xi+x2+1.0f/6.0f*x3;
	//	}
	//}
	*bi=0;
	switch(deg){
	case 2:
		x2 = xi*xi;
		if (xi>=0.5f&&xi<1.5f){
			*bi = 1.125f-1.5f*xi+0.5f*x2;
		}
		else if (xi>=-0.5f&&xi<0.5f){
			*bi = 0.75f-x2;
		}
		else if (xi>=-1.5f&&xi<-0.5f){
			*bi = 1.125f+1.5*xi+0.5f*x2;
		}

		break;
	case 3:
		x2 = xi*xi;
		x3 = x2*xi;
		if (xi>=1.0f&&xi<2.0f){
			*bi = 4.0f/3.0f-2.0f*xi+x2-1.0f/6.0f*x3;
		}
		else if (xi>=0.0f&&xi<1.0f){
			*bi = 2.0f/3.0f-x2+0.5f*x3;
		}
		else if (xi>=-1.0f&&xi<0.0f){
			*bi = 2.0f/3.0f-x2-0.5f*x3;
		}
		else if (xi>=-2.0f&&xi<-1.0f){
			*bi = 4.0f/3.0f+2.0f*xi+x2+1.0f/6.0f*x3;
		}

	}
}





//**************************************************************************************************************
__device__ inline void kernel_computeDelta3D_bSpline(float x_delta, float y_delta, float z_delta, float *delta_f, float *delta_dxf, float *delta_dyf, float *delta_dzf) {
	float dx_delta, dy_delta,dz_delta;
	int nIndex = 0;
	int q,j,i;
	float Bz1,Bz2,dBz1,dBz2,Bx1,Bx2,dBx1,dBx2,By1,By2,dBy1,dBy2;

	dx_delta = x_delta-1.0f;
	dy_delta = y_delta-1.0f;
	dz_delta = z_delta-1.0f;

	for (q=1;q<=2;q++){
		kernel_evalBSpline(z_delta+q-1.0f, 3, &Bz1);
		kernel_evalBSpline(z_delta-q, 3, &Bz2);
		kernel_evalBSpline(dz_delta+q-1.0f/2.0f, 2, &dBz1);
		kernel_evalBSpline(dz_delta-q+1.0f/2.0f, 2, &dBz2);

		for (j=1;j<=2;j++){
			kernel_evalBSpline(x_delta+j-1.0f, 3, &Bx1);
			kernel_evalBSpline(x_delta-j, 3, &Bx2);
			kernel_evalBSpline(dx_delta+j-1.0f/2.0f, 2, &dBx1);
			kernel_evalBSpline(dx_delta-j+1.0f/2.0f, 2, &dBx2);
			for (i=1;i<=2;i++){
				kernel_evalBSpline(y_delta+i-1.0f, 3, &By1);
				kernel_evalBSpline(y_delta-i, 3, &By2);
				kernel_evalBSpline(dy_delta+i-1.0f/2.0f, 2, &dBy1);
				kernel_evalBSpline(dy_delta-i+1.0f/2.0f, 2, &dBy2);

				delta_f[8*nIndex] = By1*Bx1*Bz1;
				delta_f[8*nIndex+1]= By2*Bx1*Bz1;
				delta_f[8*nIndex+2] = By1*Bx2*Bz1;
				delta_f[8*nIndex+3] = By2*Bx2*Bz1;
				delta_f[8*nIndex+4] = By1*Bx1*Bz2;
				delta_f[8*nIndex+5] = By2*Bx1*Bz2;
				delta_f[8*nIndex+6] = By1*Bx2*Bz2;
				delta_f[8*nIndex+7] = By2*Bx2*Bz2;

				delta_dxf[8*nIndex] = dBy1*Bx1*Bz1;
				delta_dxf[8*nIndex+1]  = dBy2*Bx1*Bz1;
				delta_dxf[8*nIndex+2]  = dBy1*Bx2*Bz1;
				delta_dxf[8*nIndex+3]  = dBy2*Bx2*Bz1;
				delta_dxf[8*nIndex+4] = dBy1*Bx1*Bz2;
				delta_dxf[8*nIndex+5]  = dBy2*Bx1*Bz2;
				delta_dxf[8*nIndex+6] = dBy1*Bx2*Bz2;
				delta_dxf[8*nIndex+7]  = dBy2*Bx2*Bz2;

				delta_dyf[8*nIndex] = By1*dBx1*Bz1;
				delta_dyf[8*nIndex+1] = By2*dBx1*Bz1;
				delta_dyf[8*nIndex+2] = By1*dBx2*Bz1;
				delta_dyf[8*nIndex+3] = By2*dBx2*Bz1;
				delta_dyf[8*nIndex+4] = By1*dBx1*Bz2;
				delta_dyf[8*nIndex+5] = By2*dBx1*Bz2;
				delta_dyf[8*nIndex+6] = By1*dBx2*Bz2;
				delta_dyf[8*nIndex+7] = By2*dBx2*Bz2;

				delta_dzf[8*nIndex] = By1*Bx1*dBz1;
				delta_dzf[8*nIndex+1] = By2*Bx1*dBz1;
				delta_dzf[8*nIndex+2] = By1*Bx2*dBz1;
				delta_dzf[8*nIndex+3] = By2*Bx2*dBz1;
				delta_dzf[8*nIndex+4] = By1*Bx1*dBz2;
				delta_dzf[8*nIndex+5] = By2*Bx1*dBz2;
				delta_dzf[8*nIndex+6] = By1*Bx2*dBz2;
				delta_dzf[8*nIndex+7] = By2*Bx2*dBz2;

				nIndex = nIndex+1;
			}
		}
	}

}






//**************************************************************************************************************

__device__ inline void kernel_Derivative_bSpline(float xi, float yi, float zi, int nDatax, int nDatay, int nDataz, float *delta_f, float *delta_dfx, float *delta_dfy, float *delta_dfz,const float *cMat,float *theta, float*dudt,float*model) {
	int xc,yc,zc,kx,ky,kz,q,j,i;
	int nIndex = 0.0f;
	float f=0.0f, dfx=0.0f, dfy = 0.0f, dfz = 0.0f;

	xc = floor(xi);
	yc = floor(yi);
	zc = floor(zi);

	xc = max(xc,1);
	xc = min(xc,nDatax-1);

	yc = max(yc,1);
	yc = min(yc,nDatay-1);

	zc = max(zc,1);
	zc = min(zc,nDataz-1);

	if (xi==float(nDatax))
		xc = nDatax;

	if (yi==float(nDatay))
		yc = nDatay;

	if (zi==float(nDataz))
		zc = nDataz;

	kx = xc+1;
	ky = yc+1;
	kz = zc+1;

	for (q=1;q<=2;q++){
		for (j=1;j<=2;j++){
			for (i=1;i<=2;i++){
				f += cMat[(kz-q)*(nDatax+4)*(nDatay+4)+(kx-j)*(nDatax+4)+(ky-i)]*delta_f[8*nIndex]+
					 cMat[(kz-q)*(nDatax+4)*(nDatay+4)+(kx-j)*(nDatax+4)+(ky+i-1)]*delta_f[8*nIndex+1]+
					 cMat[(kz-q)*(nDatax+4)*(nDatay+4)+(kx+j-1)*(nDatax+4)+(ky-i)]*delta_f[8*nIndex+2]+
					 cMat[(kz-q)*(nDatax+4)*(nDatay+4)+(kx+j-1)*(nDatax+4)+(ky+i-1)]*delta_f[8*nIndex+3]+
					 cMat[(kz+q-1)*(nDatax+4)*(nDatay+4)+(kx-j)*(nDatax+4)+(ky-i)]*delta_f[8*nIndex+4]+
					 cMat[(kz+q-1)*(nDatax+4)*(nDatay+4)+(kx-j)*(nDatax+4)+(ky+i-1)]*delta_f[8*nIndex+5]+
					 cMat[(kz+q-1)*(nDatax+4)*(nDatay+4)+(kx+j-1)*(nDatax+4)+(ky-i)]*delta_f[8*nIndex+6]+
					 cMat[(kz+q-1)*(nDatax+4)*(nDatay+4)+(kx+j-1)*(nDatax+4)+(ky+i-1)]*delta_f[8*nIndex+7];

				dfx += (-1.0*cMat[(kz-q)*(nDatax+4)*(nDatay+4)+(kx-j)*(nDatax+4)+(ky-i)]+cMat[(kz-q)*(nDatax+4)*(nDatay+4)+(kx-j)*(nDatax+4)+(ky-i+1)])*delta_dfx[8*nIndex]+
					 (-1.0*cMat[(kz-q)*(nDatax+4)*(nDatay+4)+(kx-j)*(nDatax+4)+(ky+i-1)]+cMat[(kz-q)*(nDatax+4)*(nDatay+4)+(kx-j)*(nDatax+4)+(ky+i)])*delta_dfx[8*nIndex+1]+
					 (-1.0*cMat[(kz-q)*(nDatax+4)*(nDatay+4)+(kx+j-1)*(nDatax+4)+(ky-i)]+cMat[(kz-q)*(nDatax+4)*(nDatay+4)+(kx+j-1)*(nDatax+4)+(ky-i+1)])*delta_dfx[8*nIndex+2]+
					 (-1.0*cMat[(kz-q)*(nDatax+4)*(nDatay+4)+(kx+j-1)*(nDatax+4)+(ky+i-1)]+cMat[(kz-q)*(nDatax+4)*(nDatay+4)+(kx+j-1)*(nDatax+4)+(ky+i)])*delta_dfx[8*nIndex+3]+
					 (-1.0*cMat[(kz+q-1)*(nDatax+4)*(nDatay+4)+(kx-j)*(nDatax+4)+(ky-i)]+ cMat[(kz+q-1)*(nDatax+4)*(nDatay+4)+(kx-j)*(nDatax+4)+(ky-i+1)])*delta_dfx[8*nIndex+4]+
					 (-1.0*cMat[(kz+q-1)*(nDatax+4)*(nDatay+4)+(kx-j)*(nDatax+4)+(ky+i-1)]+cMat[(kz+q-1)*(nDatax+4)*(nDatay+4)+(kx-j)*(nDatax+4)+(ky+i)])*delta_dfx[8*nIndex+5]+
					 (-1.0*cMat[(kz+q-1)*(nDatax+4)*(nDatay+4)+(kx+j-1)*(nDatax+4)+(ky-i)]+cMat[(kz+q-1)*(nDatax+4)*(nDatay+4)+(kx+j-1)*(nDatax+4)+(ky-i+1)])*delta_dfx[8*nIndex+6]+
					 (-1.0*cMat[(kz+q-1)*(nDatax+4)*(nDatay+4)+(kx+j-1)*(nDatax+4)+(ky+i-1)]+cMat[(kz+q-1)*(nDatax+4)*(nDatay+4)+(kx+j-1)*(nDatax+4)+(ky+i)])*delta_dfx[8*nIndex+7];

				dfy += (-1.0*cMat[(kz-q)*(nDatax+4)*(nDatay+4)+(kx-j)*(nDatax+4)+(ky-i)]+cMat[(kz-q)*(nDatax+4)*(nDatay+4)+(kx-j+1)*(nDatax+4)+(ky-i)])*delta_dfy[8*nIndex]+
					 (-1.0*cMat[(kz-q)*(nDatax+4)*(nDatay+4)+(kx-j)*(nDatax+4)+(ky+i-1)]+cMat[(kz-q)*(nDatax+4)*(nDatay+4)+(kx-j+1)*(nDatax+4)+(ky+i-1)])*delta_dfy[8*nIndex+1]+
					 (-1.0*cMat[(kz-q)*(nDatax+4)*(nDatay+4)+(kx+j-1)*(nDatax+4)+(ky-i)]+cMat[(kz-q)*(nDatax+4)*(nDatay+4)+(kx+j)*(nDatax+4)+(ky-i)])*delta_dfy[8*nIndex+2]+
					 (-1.0*cMat[(kz-q)*(nDatax+4)*(nDatay+4)+(kx+j-1)*(nDatax+4)+(ky+i-1)]+cMat[(kz-q)*(nDatax+4)*(nDatay+4)+(kx+j)*(nDatax+4)+(ky+i-1)])*delta_dfy[8*nIndex+3]+
					 (-1.0*cMat[(kz+q-1)*(nDatax+4)*(nDatay+4)+(kx-j)*(nDatax+4)+(ky-i)]+ cMat[(kz+q-1)*(nDatax+4)*(nDatay+4)+(kx-j+1)*(nDatax+4)+(ky-i)])*delta_dfy[8*nIndex+4]+
					 (-1.0*cMat[(kz+q-1)*(nDatax+4)*(nDatay+4)+(kx-j)*(nDatax+4)+(ky+i-1)]+cMat[(kz+q-1)*(nDatax+4)*(nDatay+4)+(kx-j+1)*(nDatax+4)+(ky+i-1)])*delta_dfy[8*nIndex+5]+
					 (-1.0*cMat[(kz+q-1)*(nDatax+4)*(nDatay+4)+(kx+j-1)*(nDatax+4)+(ky-i)]+cMat[(kz+q-1)*(nDatax+4)*(nDatay+4)+(kx+j)*(nDatax+4)+(ky-i)])*delta_dfy[8*nIndex+6]+
					 (-1.0*cMat[(kz+q-1)*(nDatax+4)*(nDatay+4)+(kx+j-1)*(nDatax+4)+(ky+i-1)]+cMat[(kz+q-1)*(nDatax+4)*(nDatay+4)+(kx+j)*(nDatax+4)+(ky+i-1)])*delta_dfy[8*nIndex+7];
				     
				     
				dfz += (-1.0*cMat[(kz-q)*(nDatax+4)*(nDatay+4)+(kx-j)*(nDatax+4)+(ky-i)]+cMat[(kz-q+1)*(nDatax+4)*(nDatay+4)+(kx-j)*(nDatax+4)+(ky-i)])*delta_dfz[8*nIndex]+
					 (-1.0*cMat[(kz-q)*(nDatax+4)*(nDatay+4)+(kx-j)*(nDatax+4)+(ky+i-1)]+cMat[(kz-q+1)*(nDatax+4)*(nDatay+4)+(kx-j)*(nDatax+4)+(ky+i-1)])*delta_dfz[8*nIndex+1]+
					 (-1.0*cMat[(kz-q)*(nDatax+4)*(nDatay+4)+(kx+j-1)*(nDatax+4)+(ky-i)]+cMat[(kz-q+1)*(nDatax+4)*(nDatay+4)+(kx+j-1)*(nDatax+4)+(ky-i)])*delta_dfz[8*nIndex+2]+
					 (-1.0*cMat[(kz-q)*(nDatax+4)*(nDatay+4)+(kx+j-1)*(nDatax+4)+(ky+i-1)]+cMat[(kz-q+1)*(nDatax+4)*(nDatay+4)+(kx+j-1)*(nDatax+4)+(ky+i-1)])*delta_dfz[8*nIndex+3]+
					 (-1.0*cMat[(kz+q-1)*(nDatax+4)*(nDatay+4)+(kx-j)*(nDatax+4)+(ky-i)]+ cMat[(kz+q)*(nDatax+4)*(nDatay+4)+(kx-j)*(nDatax+4)+(ky-i)])*delta_dfz[8*nIndex+4]+
					 (-1.0*cMat[(kz+q-1)*(nDatax+4)*(nDatay+4)+(kx-j)*(nDatax+4)+(ky+i-1)]+cMat[(kz+q)*(nDatax+4)*(nDatay+4)+(kx-j)*(nDatax+4)+(ky+i-1)])*delta_dfz[8*nIndex+5]+
					 (-1.0*cMat[(kz+q-1)*(nDatax+4)*(nDatay+4)+(kx+j-1)*(nDatax+4)+(ky-i)]+cMat[(kz+q)*(nDatax+4)*(nDatay+4)+(kx+j-1)*(nDatax+4)+(ky-i)])*delta_dfz[8*nIndex+6]+
					 (-1.0*cMat[(kz+q-1)*(nDatax+4)*(nDatay+4)+(kx+j-1)*(nDatax+4)+(ky+i-1)]+cMat[(kz+q)*(nDatax+4)*(nDatay+4)+(kx+j-1)*(nDatax+4)+(ky+i-1)])*delta_dfz[8*nIndex+7];

				nIndex = nIndex+1;
			}
		}
	}


	//dudt[0]=-1.0f*theta[2]*dfy;
	//dudt[1]=-1.0f*theta[2]*dfx;
	//flipped in bSpline
	dudt[0]=-1.0f*theta[2]*dfx;
	dudt[1]=-1.0f*theta[2]*dfy;
	dudt[4]=theta[2]*dfz;
	dudt[2]=f;
	dudt[3]=1.0f;
	*model = theta[3]+theta[2]*f;
	
	//return pd;
}
