//Copyright (c) 2017 Ries Lab, European Molecular Biology Laboratory, Heidelberg.
//author: Yiming Li
//email: yiming.li@embl.de
//date: 2018.10.08


/*! 
* [out] = simSplinePSF_c(Npixels,coeff,I,bg,cor)
//************varargin**********************
*  \varargin:
*  \Npxiels: size of output ROI
*  \coeff: spline coeff
*  \I: intensity
*  \bg: background
*  \cor: coordinates

//************Output**********************
*  \out: image stack

*/
#ifdef _WIN32
    #include <windows.h>
#endif

#pragma comment(lib, "kernel32.lib")

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mex.h>
//#include "definitions.h"
//#include "CPUmleFit_LM.h"

void kernel_computeDelta3D(float x_delta, float y_delta, float z_delta, float *delta_f) {
    
	int i,j,k;
	float cx,cy,cz;

	cz = 1.0;
	for(i=0;i<4;i++){
		cy = 1.0;
		for(j=0;j<4;j++){
			cx = 1.0;
			for(k=0;k<4;k++){
				delta_f[i*16+j*4+k] = cz * cy * cx;
				cx = cx * x_delta;
			}
			cy = cy * y_delta;
		}
		cz= cz * z_delta;
	}
}

 float fAt3Dj(int xc, int yc, int zc, int xsize, int ysize, int zsize, float*delta_f, float *coeff) {
	int i;
	float temp =0;
	//float dudt_temp[NV_PSP] = {0};//,temp;
	
	//for (i=0;i<NV_PSP;i++) dudt[i]=0;
	
	xc = fmax(xc,0);
	xc = fmin(xc,xsize-1);

	yc = fmax(yc,0);
	yc = fmin(yc,ysize-1);

	zc = fmax(zc,0);
	zc = fmin(zc,zsize-1);
	
	

	for (i=0;i<64;i++){		
		temp+=delta_f[i]*coeff[i*(xsize*ysize*zsize)+zc*(xsize*ysize)+yc*xsize+xc];
	}
	return temp;
}

//*******************************************************************************************
void mexFunction(int nlhs, mxArray *plhs[],	int	nrhs, const	mxArray	*prhs[]) {
	/*!
	*  \brief Entry point in the code for Matlab.  Equivalent to main().
	*  \param nlhs number of left hand mxArrays to return
	*  \param plhs array of pointers to the output mxArrays
	*  \param nrhs number of input mxArrays
	*  \param prhs array of pointers to the input mxArrays.
	*/

	//int blockx=0;
	//int threadx=0;

	//float PSFSigma=1, Ax=0, Ay=0, Bx=0, By=0, gamma=0.5, d=0.5, PSFSigma_y=1,initZ=-1;
	//int iterations=50;
	//float * startParameters = 0;
	//float *data=0, *d_data=0;
	//float *Parameters=0,*CRLBs=0,*LogLikelihood=0;
	//size_t Ndim, Nfitraw, fittype, sz, cameratype=0;//default cameratype 0 EMCCD
	//mwSize noParameters;
	////sCMOS
	//float *varim=0;

	////Spline

	//int silent = 0;


	// simPSF
	int Npixels,Nfitraw, xstart, ystart, zstart;
	float *coeff,*cor, *data;
	float I, bg, xcenter, ycenter, zcenter,xc, yc, zc,off,temp;
	const mwSize *datasize_spline;
	const mwSize *datasize=0;
	int spline_xsize, spline_ysize, spline_zsize,ii,jj,kk;
	size_t Ndim;
	mwSize dataDim[3]={0};
	float delta_f[64]={0};



	Npixels = (int)mxGetScalar(prhs[0]);
	coeff =(float *) mxGetData(prhs[1]);
	I = (float)mxGetScalar(prhs[2]);
	bg = (float)mxGetScalar(prhs[3]);
	cor = (float *) mxGetData(prhs[4]);

	datasize_spline=mxGetDimensions(prhs[1]);
	spline_xsize = datasize_spline[0];
	spline_ysize = datasize_spline[1];
	spline_zsize = datasize_spline[2];

	datasize=mxGetDimensions(prhs[4]);
	Ndim=mxGetNumberOfDimensions(prhs[0]);
    Nfitraw=datasize[0];
	dataDim[0]= Npixels;
	dataDim[1]= Npixels;
	dataDim[2]= Nfitraw;
	plhs[0] = mxCreateNumericArray(3, dataDim, mxSINGLE_CLASS, mxREAL);
	off = floor(((spline_xsize+1)-Npixels)/2.0);
	data = (float*) mxGetData(plhs[0]);

	for (int kk = 0; kk<Nfitraw; kk++){
		xcenter = cor[kk];
		ycenter = cor[Nfitraw+kk];
		zcenter = cor[2*Nfitraw+kk];

		xc = -1.0*((xcenter-float(Npixels)/2)+0.5);
		yc = -1.0*((ycenter-float(Npixels)/2)+0.5);

		xstart = floor(xc);
		xc = xc-xstart;

		ystart = floor(yc);
		yc = yc-ystart;

		zstart = floor(zcenter);
		zc = zcenter -zstart;

		kernel_computeDelta3D(xc, yc, zc, delta_f);
		for (ii=0;ii<Npixels;ii++) for(jj=0;jj<Npixels;jj++) {
			temp = fAt3Dj(ii+xstart+off,jj+ystart+off,zstart,spline_xsize,spline_ysize,spline_zsize,delta_f,coeff);
			data[Npixels*Npixels*kk+Npixels*jj+ii] = temp*I + bg;		
		}
	}

	return;
}
