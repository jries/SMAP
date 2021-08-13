//Copyright (c) 2021 Ries Lab, European Molecular Biology Laboratory, Heidelberg &Li Lab, Southern University of Science and Technology, Shenzhen.
//author: Yiming Li
//email: liym2019@sustech.edu.cn
//date: 2021.05.02

#ifdef _WIN32
    #include <windows.h>
#endif

#pragma comment(lib, "kernel32.lib")

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mex.h>
#include "definitions.h"
#include "CPUmleFit_LM.h"

//*******************************************************************************************
void mexFunction(int nlhs, mxArray *plhs[],	int	nrhs, const	mxArray	*prhs[]) {
	/*!
	*  \brief Entry point in the code for Matlab.  Equivalent to main().
	*  \param nlhs number of left hand mxArrays to return
	*  \param plhs array of pointers to the output mxArrays
	*  \param nrhs number of input mxArrays
	*  \param prhs array of pointers to the input mxArrays.
	*/
    float * test;
	int blockx=0;
	int threadx=0;
	float PSFSigma = 1.0;
	const mwSize *datasize=0;
	//const int *datasize=0;
	int iterations=50,i;
	float * startParameters = 0;
	float *data=0, *d_data=0, *dTAll=0, *initZall=0;
	float *Parameters=0,*CRLBs=0,*LogLikelihood=0;
	//size_t Ndim, fittype, cameratype=0,noParameters,Nfitraw,noChannels,sz;;//default cameratype 0 EMCCD
	int Ndim, fittype, cameratype=0,noParameters,Nfitraw,noChannels,sz;;//default cameratype 0 EMCCD
	//int noParameters,Nfitraw,noChannels,sz;
	int sumShared=0;
	int *shared, *d_shared;
	//sCMOS
	float *varim=0;

	//Spline
	int spline_xsize, spline_ysize, spline_zsize;
	float *coeff=0,*d_coeff=0;
	const mwSize *datasize_spline;
	//const int *datasize_spline;
	int silent = 0;



	//fun with timing
	//mwSize freq;

	//mwSize start,stop;
#ifdef _WIN32
	LARGE_INTEGER freq;
	LARGE_INTEGER start,stop;
	double togpu=0,fit=0,fromgpu=0,cleanup=0;
	QueryPerformanceFrequency( &freq );
	QueryPerformanceCounter(&start);
#endif


	//input checks
	//if (nrhs<2)
	//	mexErrMsgIdAndTxt("CPUmleFit_LM:WrongNoArgs","Inputs must include: CPUmleFit_LM(data,fittype,iterations,startparameter,var,silent) !\n");
	//else
	//{
	if (mxGetClassID(prhs[0])!=mxSINGLE_CLASS)
		mexErrMsgIdAndTxt("CPUmleFit_LM_MultiChannel:WrongDataType","Data must be comprised of single floats!\n");

	datasize=mxGetDimensions(prhs[0]);
	Ndim=mxGetNumberOfDimensions(prhs[0]);

	Nfitraw=datasize[2];

	if (datasize[0] > IMSZBIG)
		mexErrMsgIdAndTxt("CPUmleFit_LM_MultiChannel:SubregionSize","X,Y dimension of data must be smaller than 51.\n");
	if (datasize[1]!=datasize[0])
		mexErrMsgIdAndTxt("CPUmleFit_LM_MultiChannel:SubregionShape","Fit Box must be square");
	sz=datasize[0];
	noChannels = datasize[3];

	
	data=(float *) mxGetData(prhs[0]);


	fittype = (int) mxGetScalar(prhs[1]);
	shared = (int *) mxGetData(prhs[2]);

	iterations=(int) mxGetScalar(prhs[3]);

	switch (fittype){
		case 1: 
			PSFSigma = (float)mxGetScalar(prhs[4]);
			break;
		case 2:
			coeff =(float *) mxGetData(prhs[4]);
			datasize_spline=mxGetDimensions(prhs[4]);
			spline_xsize = datasize_spline[0];
			spline_ysize = datasize_spline[1];
			spline_zsize = datasize_spline[2];
			break;
		}

	if (mxGetClassID(prhs[4])!=mxSINGLE_CLASS)
		mexErrMsgIdAndTxt("CPUmleFit_LM_MultiChannel:WrongDataType","coeff must be comprised of single floats!\n");

	/*coeff =(float *) mxGetData(prhs[3]);
	datasize_spline=mxGetDimensions(prhs[3]);
	spline_xsize = datasize_spline[0];
	spline_ysize = datasize_spline[1];
	spline_zsize = datasize_spline[2];*/

	dTAll = (float *) mxGetData(prhs[5]);

	if (nrhs > 6) {
		if (mxGetNumberOfElements(prhs[6]) != mxGetNumberOfElements(prhs[0]))
			cameratype = 0;//EMCCD
		else
		{
			if (mxGetClassID(prhs[6]) != mxSINGLE_CLASS)
				mexErrMsgIdAndTxt("CPUmleFit_LM_MultiChannel:WrongDataType", "var must be comprised of single floats!\n");
			varim = (float*)mxGetData(prhs[6]);
			cameratype = 1;//sCMOS
		}
	}


	if (nrhs > 7)
		silent = (int)mxGetScalar(prhs[7]);


	if (nrhs > 8) {
		initZall = (float*)mxGetData(prhs[8]);
	}
	else {
		initZall = (float*)malloc(Nfitraw * sizeof(float));
		if (initZall) {
			for (int i = 0; i < Nfitraw; i++) {
				initZall[i - 1] = spline_zsize / 2.0f;
			}
		}
		//memset(initZall, spline_zsize / 2.0f, Nfitraw * sizeof(float));

	}

	if (silent == 0)
	{
		mexPrintf("subregion size: %d fittype: %d\n", sz, fittype);
		//mexPrintf("Number of fits: %d, Number padded fits: %d\n", Nfitraw, Nfitspad);
		mexPrintf("Number of fits: %d no padding\n", Nfitraw);
		//mexPrintf("Memory required: %dMB Memory available: %dMB\n",requiredMemory/(1024*1024),availableMemory);
	}

	sumShared=0;
	for (i=0;i<5;i++){
		sumShared+=shared[i];
	}
	noParameters = (5*noChannels-sumShared*(noChannels-1));

	/*mexPrintf("noChannels %f noParameters %f sumShared %f Nfitraw%f test%f \n", noChannels,noParameters,sumShared, Nfitraw, 1);*/


	plhs[0]=mxCreateNumericMatrix(Nfitraw, noParameters+1, mxSINGLE_CLASS, mxREAL);
	plhs[1]=mxCreateNumericMatrix(Nfitraw, noParameters, mxSINGLE_CLASS, mxREAL);

	plhs[2]=mxCreateNumericMatrix(Nfitraw, 1, mxSINGLE_CLASS, mxREAL);

	//****************************test*****************************
	//plhs[3]=mxCreateNumericMatrix(100, 1, mxSINGLE_CLASS, mxREAL);
	//test = (float*) mxGetData(plhs[3]);
	//test[0]=1;
	//test[1]=noChannels;
	//test[2]=noParameters;
	//test[3]=sumShared;
	//test[4]=Nfitraw;
	//test[5]=shared[0];
	//test[6]=shared[1];
	//test[7]=shared[2];
	//test[8]=shared[3];
	//test[9]=shared[4];
	//****************************test*****************************

	
	Parameters = (float*) mxGetData(plhs[0]);
	CRLBs = (float*) mxGetData(plhs[1]);
	LogLikelihood = (float*) mxGetData(plhs[2]);

	//check allocations
	if (plhs[0] == NULL || plhs[1] == NULL || plhs[2] == NULL)
		mexErrMsgIdAndTxt("CPUmleFit_LM_MultiChannel:NotEnoughMemory","Could not allocate memory for results.\n");
#ifdef _WIN32
	QueryPerformanceCounter(&stop);
	togpu = (double) (stop.QuadPart-start.QuadPart)/freq.QuadPart;
	QueryPerformanceCounter(&start);
#endif


	//setup kernel
	//	blockx = (int) ceil( (float)Nfitraw/(float)BlockSize);
	//	threadx= BlockSize;
	//
	//#ifdef _DEBUG
	//	mexPrintf("threadx: %d,blockx: %d\n", threadx, blockx);
	//#endif
	/*mexPrintf("before iterations %f\n",(double)iterations);*/


switch (fittype){
case 1:
	if (cameratype == 0) {
		for (int ii = 0; ii < Nfitraw; ii++)
			kernel_MLEFit_sigma_EMCCD_multi(ii, data, PSFSigma, dTAll, (int)sz, iterations, noParameters, noChannels, Parameters, CRLBs, LogLikelihood, (const int)Nfitraw, shared);
	}
	else if (cameratype == 1) {
		for (int ii = 0; ii < Nfitraw; ii++)
			kernel_MLEFit_sigma_sCMOS_multi(ii, data, PSFSigma, dTAll, (int)sz, iterations, noParameters, noChannels, Parameters, CRLBs, LogLikelihood, (const int)Nfitraw, shared, varim);
	}

	//extern void kernel_MLEFit_sigma_EMCCD_multi_wrapper(dim3 dimGrid, dim3 dimBlock, const float *d_data, const float PSFSigma, const int sz, const int iterations, 
	//const int NV, const int noChannels, float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits, const int*d_shared);


	break;
case 2:
	//if (initZ <0)
	//	initZ = spline_zsize/2.0f;
	if (cameratype == 0) {
		for (int ii = 0; ii < Nfitraw; ii++)
			kernel_splineMLEFit_z_EMCCD_multi(ii, data, coeff, dTAll, spline_xsize, spline_ysize, spline_zsize, (int)sz, iterations, noParameters, noChannels, Parameters, CRLBs, LogLikelihood, initZall, (const int)Nfitraw, shared);
	}
	else if (cameratype == 1) {
		for (int ii = 0; ii < Nfitraw; ii++)
			kernel_splineMLEFit_z_sCMOS_multi(ii, data, coeff, dTAll, spline_xsize, spline_ysize, spline_zsize, (int)sz, iterations, noParameters, noChannels, Parameters, CRLBs, LogLikelihood, initZall, (const int)Nfitraw, shared, varim);

	}
}



	//if (initZ <0)
	//initZ = spline_zsize/2.0f;

	//for (int ii = 0; ii<Nfitraw; ii++)
	//	kernel_splineMLEFit_z_EMCCD_multi(ii,data,coeff,dTAll,spline_xsize,spline_ysize,spline_zsize,(int) sz,iterations,noParameters,noChannels,Parameters,CRLBs,LogLikelihood,initZ,(const int) Nfitraw,shared);
			
    #ifdef _WIN32
	QueryPerformanceCounter(&stop);    
	fit = (double)(stop.QuadPart-start.QuadPart)/freq.QuadPart;
	QueryPerformanceCounter(&start);

	QueryPerformanceCounter(&stop);
	fromgpu = (double)(stop.QuadPart-start.QuadPart)/freq.QuadPart;
	QueryPerformanceCounter(&start);

	//cleanup
	QueryPerformanceCounter(&stop);
     

	if (silent==0)
	{//mexPrintf("Memory copies to GPU %f seconds\n", togpu);
		mexPrintf("Actual fitting time %f Actual fits per second %f\n", fit,
			Nfitraw/fit);
		//mexPrintf("Memory copies from GPU %f seconds\n", fromgpu);
		//mexPrintf("Clean up from GPU %f seconds\n", (double)(stop.QuadPart-start.QuadPart)/freq.QuadPart);
	}
    #endif 
	return;
}
