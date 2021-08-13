#include "CPUsplineLib.h"
//#include <cuda_runtime.h>
#include "definitions.h"
#include <math.h>
#include <string.h>
//#define pi 3.141592f
#include "mex.h"


 //**************************************************************************************************************************************************
// This function for calculation of the common term for Cspline is adpopted from
//"Analyzing Single Molecule Localization Microscopy Data Using Cubic Splines", Hazen Babcok, Xiaowei Zhuang,Scientific Report, 1, 552 , 2017.
 void kernel_computeDelta3D(float x_delta, float y_delta, float z_delta, float *delta_f, float *delta_dxf, float *delta_dyf, float *delta_dzf) {
    
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
 int kernel_cholesky(float *A,int n, float *L, float*U) {
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
 void kernel_luEvaluate(float *L,float *U, float *b, const int n, float *x) {
	//Ax = b -> LUx = b. Then y is defined to be Ux

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

	//mexPrintf("after first for \n");
	//return;
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

 void kernel_DerivativeSpline(int xc, int yc, int zc, int xsize, int ysize, int zsize, float *delta_f, float *delta_dxf, float *delta_dyf, float *delta_dzf,const float *coeff,float *theta, float*dudt,float*model) {
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
	}
	dudt[0]*=-1*theta[2];
	dudt[1]*=-1*theta[2];
	dudt[4]*=theta[2];
	dudt[2]=temp;
	dudt[3]=1;
	*model = theta[3]+theta[2]*temp;
}

 //**************************************************************************************************************

 void kernel_DerivativeSpline_multiChannel(int xc, int yc, int zc, int xsize, int ysize, int zsize, float *delta_f, float *delta_dxf, float *delta_dyf, float *delta_dzf,const float *coeff,float *theta, float*dudt,float*model) {
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