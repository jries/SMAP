//Copyright (c) 2021 Ries Lab, European Molecular Biology Laboratory, Heidelberg &Li Lab, Southern University of Science and Technology, Shenzhen.
//author: Yiming Li
//email: liym2019@sustech.edu.cn
//date: 2021.05.02


#include <windows.h>
#include <tchar.h>
#include <io.h>
#pragma comment(lib, "kernel32.lib")

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mex.h>
#include <cuda_runtime.h>
#include <definitions.h>

#ifndef max
//! not defined in the C standard used by visual studio
#define max(a,b) (((a) > (b)) ? (a) : (b))
#endif
#ifndef min
//! not defined in the C standard used by visual studio
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif

extern void kernel_splineMLEFit_z_EMCCD_multi_wrapper(dim3 dimGrid, dim3 dimBlock, const float* d_data, const float* d_coeff, const float* dTAll, const int spline_xsize, const int spline_ysize, const int spline_zsize, const int sz, const int iterations,
	const int NV, const int noChannels, float* d_Parameters, float* d_CRLBs, float* d_LogLikelihood, float* d_initZall, const int Nfits, const int* d_shared);

extern void kernel_MLEFit_sigma_EMCCD_multi_wrapper(dim3 dimGrid, dim3 dimBlock, const float* d_data, const float PSFSigma, const float* dTAll, const int sz, const int iterations,
	const int NV, const int noChannels, float* d_Parameters, float* d_CRLBs, float* d_LogLikelihood, const int Nfits, const int* d_shared);

extern void kernel_splineMLEFit_z_sCMOS_multi_wrapper(dim3 dimGrid, dim3 dimBlock, const float* d_data, const float* d_coeff, const float* dTAll, const int spline_xsize, const int spline_ysize, const int spline_zsize, const int sz, const int iterations,
	const int NV, const int noChannels, float* d_Parameters, float* d_CRLBs, float* d_LogLikelihood, float* d_initZall, const int Nfits, const int* d_shared, const float* d_varim);

extern void kernel_MLEFit_sigma_sCMOS_multi_wrapper(dim3 dimGrid, dim3 dimBlock, const float* d_data, const float PSFSigma, const float* dTAll, const int sz, const int iterations,
	const int NV, const int noChannels, float* d_Parameters, float* d_CRLBs, float* d_LogLikelihood, const int Nfits, const int* d_shared, const float* d_varim);


void CUDAERROR(const char* instr, int lineNumber);
void cudasafe(cudaError_t err, char* str, int lineNumber);

//*******************************************************************************************
void CUDAERROR(const char* instr, int lineNumber) {
	cudaError_t errornum;
	const char* str;
	if (errornum = cudaGetLastError()) {
		//reset all cuda devices
		int deviceCount = 0;
		int ii = 0;
		cudasafe(cudaGetDeviceCount(&deviceCount), "cudaGetDeviceCount", __LINE__); //query number of GPUs
		for (ii = 0; ii < deviceCount; ii++) {
			cudaSetDevice(ii);
			cudaDeviceReset();
		}
		str = cudaGetErrorString(errornum);
		//mexPrintf("gpuGaussNDmle(line %i): %s in %s\n",lineNumber, str, instr);
		cudaDeviceReset();
		mexErrMsgIdAndTxt("GPUmleFit_LM_MultiChannel:cudaFail", "GPUmleFit_LM_MultiChannel(line %i): %s in %s\n", lineNumber, str, instr);
		exit(1); // might not stop matlab
	}
}

//*******************************************************************************************
void cudasafe(cudaError_t err, char* str, int lineNumber)
{
	if (err != cudaSuccess)
	{
		//reset all cuda devices
		int deviceCount = 0;
		int ii = 0;
		cudasafe(cudaGetDeviceCount(&deviceCount), "cudaGetDeviceCount", __LINE__); //query number of GPUs
		for (ii = 0; ii < deviceCount; ii++) {
			cudaSetDevice(ii);
			cudaDeviceReset();
		}
		mexErrMsgIdAndTxt("GPUmleFit_LM_MultiChannel:cudaFail", "%s failed with error code %i at line %d\n", str, err, lineNumber);
		exit(1); // might not stop matlab
	}
}

//*******************************************************************************************
void cudaavailable(int silent) {
	int driverVersion = 0, runtimeVersion = 0, deviceCount = 0;
	if (cudaSuccess == cudaDriverGetVersion(&driverVersion)) {
		if (silent == 0)
			mexPrintf("CUDA driver version: %d\n", driverVersion);

	}
	else {
		mexErrMsgIdAndTxt("GPUmleFit_LM_MultiChannel:nodriver", "Could not query CUDA driver version\n");
	}
	if (cudaSuccess == cudaRuntimeGetVersion(&runtimeVersion)) {
		if (silent == 0)
			mexPrintf("CUDA driver version: %d\n", runtimeVersion);

	}
	else {
		mexErrMsgIdAndTxt("GPUmleFit_LM_MultiChannel:noruntime", "Could not query CUDA runtime version\n");
	}
	if (cudaSuccess == cudaGetDeviceCount(&deviceCount)) {
		if (silent == 0)
			mexPrintf("CUDA devices detected: %d\n", deviceCount);

	}
	else {
		mexErrMsgIdAndTxt("GPUmleFit_LM_MultiChannel:nodevices", "Could not query CUDA device count\n", runtimeVersion);
	}
	if (deviceCount < 1) {
		mexErrMsgIdAndTxt("GPUmleFit_LM_MultiChannel:NoDevice", "No CUDA capable devices were detected");
	}
}


//*******************************************************************************************
void mexFunction(int nlhs, mxArray* plhs[], int	nrhs, const	mxArray* prhs[]) {
	/*!
	*  \brief Entry point in the code for Matlab.  Equivalent to main().
	*  \param nlhs number of left hand mxArrays to return
	*  \param plhs array of pointers to the output mxArrays
	*  \param nrhs number of input mxArrays
	*  \param prhs array of pointers to the input mxArrays.
	*/
	//#define _DEBUG 1
	int blockx = 0;
	int threadx = 0;
	const mwSize* datasize = 0;
	float PSFSigma = 1.0;

	int iterations = 50, i;
	float* startParameters = 0;
	float* data = 0, * d_data = 0, * dTAll = 0, * d_dTAll = 0, * initZall = 0, * d_initZall = 0;
	float* d_Parameters = 0, * d_CRLBs = 0, * d_LogLikelihood = 0;
	int Ndim, fittype, cameratype = 0, noParameters, Nfitraw, noChannels, sz;;//default cameratype 0 EMCCD
	//int noParameters,Nfitraw,noChannels,sz;
	int sumShared = 0;
	int* shared, * d_shared;
	//sCMOS
	float* varim = 0, * d_varim = 0;

	//Spline
	int spline_xsize, spline_ysize, spline_zsize;
	float* coeff = 0, * d_coeff = 0;
	const mwSize* datasize_spline;
	int silent = 0;

	//fun with timing
	LARGE_INTEGER freq;
	double togpu = 0, fit = 0, fromgpu = 0, cleanup = 0;
	LARGE_INTEGER start, stop;
	QueryPerformanceFrequency(&freq);
	QueryPerformanceCounter(&start);


	//input checks
	//if (nrhs<2)
	//	mexErrMsgIdAndTxt("GPUmleFit_LM:WrongNoArgs","Inputs must include: GPUmleFit_LM(data,fittype,iterations,startparameter,var,silent) !\n");
	//else
	//{
	if (mxGetClassID(prhs[0]) != mxSINGLE_CLASS)
		mexErrMsgIdAndTxt("GPUmleFit_LM_MultiChannel:WrongDataType", "Data must be comprised of single floats!\n");

	datasize = mxGetDimensions(prhs[0]);
	Ndim = mxGetNumberOfDimensions(prhs[0]);

	//if (Ndim==2)Nfitraw=1;else 
	Nfitraw = datasize[2];

	if (datasize[0] > IMSZBIG)
		mexErrMsgIdAndTxt("GPUmleFit_LM_MultiChannel:SubregionSize", "X,Y dimension of data must be smaller than 51.\n");
	if (datasize[1] != datasize[0])
		mexErrMsgIdAndTxt("GPUmleFit_LM_MultiChannel:SubregionShape", "Fit Box must be square");
	sz = datasize[0];
	noChannels = datasize[3];

	data = (float*)mxGetData(prhs[0]);


	fittype = (int)mxGetScalar(prhs[1]);
	shared = (int*)mxGetData(prhs[2]);

	iterations = (int)mxGetScalar(prhs[3]);

	switch (fittype) {
	case 1:
		PSFSigma = (float)mxGetScalar(prhs[4]);
		break;
	case 2:
		coeff = (float*)mxGetData(prhs[4]);
		datasize_spline = mxGetDimensions(prhs[4]);
		spline_xsize = datasize_spline[0];
		spline_ysize = datasize_spline[1];
		spline_zsize = datasize_spline[2];
		break;
	}

	if (mxGetClassID(prhs[4]) != mxSINGLE_CLASS)
		mexErrMsgIdAndTxt("GPUmleFit_LM_MultiChannel:WrongDataType", "coeff must be comprised of single floats!\n");



	dTAll = (float*)mxGetData(prhs[5]);

	if (nrhs > 6) {
		if (mxGetNumberOfElements(prhs[6]) != mxGetNumberOfElements(prhs[0]))
			cameratype = 0;//EMCCD
		else
		{
			if (mxGetClassID(prhs[6]) != mxSINGLE_CLASS)
				mexErrMsgIdAndTxt("GPUmleFit_LM_MultiChannel:WrongDataType", "var must be comprised of single floats!\n");
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


	//query GPUs
	int deviceCount = 0;
	int driverVersion = 0;
	int runtimeVersion = 0;
	cudaDeviceProp deviceProp;

	cudaavailable(silent);
	//mexPrintf("after cudaavailable\n");


	cudasafe(cudaGetDeviceProperties(&deviceProp, 0), "Could not get properties for device 0.", __LINE__);

	//cudasafe(cudaSetDevice(0),"Could not select GPU 0.",__LINE__);
	//cudasafe(cudaDeviceReset(),"Error on cudaDeviceReset().",__LINE__); //release context so future cudaSetDevice calls work

	if (silent == 0)
		mexPrintf("Using GPU %s\n", deviceProp.name);

	int BlockSize = BSZ; //for no shared scenario, we don't care!

	//check if we are doing something silly and allocating more memory than the card has
	const size_t availableMemory = deviceProp.totalGlobalMem / (1024 * 1024);
	size_t requiredMemory = 0;
	requiredMemory += 3 * sz * sz * Nfitraw * sizeof(float) * noChannels;
	requiredMemory += spline_xsize * spline_ysize * spline_zsize * 64 * sizeof(float) * noChannels;

	if (requiredMemory > 0.95 * deviceProp.totalGlobalMem)
		mexErrMsgIdAndTxt("GPUmleFit_LM_MultiChannel:NotEnoughMemory", "Trying to allocation %dMB. GPU only has %dMB. Please break your fitting into multiple smaller runs.\n",
			requiredMemory / (1024 * 1024), availableMemory);
	if (silent == 0)
	{
		mexPrintf("subregion size: %d fittype: %d\n", sz, fittype);
		//mexPrintf("Number of fits: %d, Number padded fits: %d\n", Nfitraw, Nfitspad);
		mexPrintf("Number of fits: %d no padding\n", Nfitraw);
		mexPrintf("Memory required: %dMB Memory available: %dMB\n", requiredMemory / (1024 * 1024), availableMemory);
	}

	//create device variable for data and copy to device
	cudasafe(cudaMalloc((void**)&d_data, sz * sz * Nfitraw * noChannels * sizeof(float)), "Failed cudaMalloc on d_data.", __LINE__);
	cudasafe(cudaMemset(d_data, 0, sz * sz * Nfitraw * noChannels * sizeof(float)), "Failed cudaMemset on d_data.", __LINE__);
	cudasafe(cudaMemcpy(d_data, data, sz * sz * Nfitraw * noChannels * sizeof(float), cudaMemcpyHostToDevice), "Failed cudaMemcpy on d_data.", __LINE__);

	//sCMOS
	if (cameratype == 1) {
		cudasafe(cudaMalloc((void**)&d_varim, sz * sz * Nfitraw * sizeof(float)), "Failed cudaMalloc on d_varim.", __LINE__);
		cudasafe(cudaMemset(d_varim, 0, sz * sz * Nfitraw * sizeof(float)), "Failed cudaMemset on d_varim.", __LINE__);
		cudasafe(cudaMemcpy(d_varim, varim, sz * sz * Nfitraw * sizeof(float), cudaMemcpyHostToDevice), "Failed cudaMemcpy on d_varim.", __LINE__);
	}

	//create device variable for coefficient
	if (fittype == 2) {
		cudasafe(cudaMalloc((void**)&d_coeff, spline_xsize * spline_ysize * spline_zsize * 64 * noChannels * sizeof(float)), "Failed cudaMalloc on d_coeff.", __LINE__);
		cudasafe(cudaMemset(d_coeff, 0, spline_xsize * spline_ysize * spline_zsize * 64 * noChannels * sizeof(float)), "Failed cudaMemset on d_coeff.", __LINE__);
		cudasafe(cudaMemcpy(d_coeff, coeff, spline_xsize * spline_ysize * spline_zsize * 64 * noChannels * sizeof(float), cudaMemcpyHostToDevice), "Failed cudaMemcpy on d_coeff.", __LINE__);
	}

	//create device variable for dTAll
	cudasafe(cudaMalloc((void**)&d_dTAll, NV_PS * Nfitraw * noChannels * 2 * sizeof(float)), "Failed cudaMalloc on d_dTAll.", __LINE__);
	cudasafe(cudaMemset(d_dTAll, 0, NV_PS * Nfitraw * noChannels * 2 * sizeof(float)), "Failed cudaMemset on d_dTAll.", __LINE__);
	cudasafe(cudaMemcpy(d_dTAll, dTAll, NV_PS * Nfitraw * noChannels * 2 * sizeof(float), cudaMemcpyHostToDevice), "Failed cudaMemcpy on d_dTAll.", __LINE__);

	//create device variable for initZall
	cudasafe(cudaMalloc((void**)&d_initZall, Nfitraw * sizeof(float)), "Failed cudaMalloc on d_initZall.", __LINE__);
	cudasafe(cudaMemset(d_initZall, 0, Nfitraw * sizeof(float)), "Failed cudaMemset on d_initZall.", __LINE__);
	cudasafe(cudaMemcpy(d_initZall, initZall, Nfitraw * sizeof(float), cudaMemcpyHostToDevice), "Failed cudaMemcpy on d_initZall.", __LINE__);


	//create device variable for shared
	cudasafe(cudaMalloc((void**)&d_shared, 5 * Nfitraw * sizeof(int)), "Failed cudaMalloc on d_dTAll.", __LINE__);
	cudasafe(cudaMemset(d_shared, 0, NV_PS * Nfitraw * sizeof(int)), "Failed cudaMemset on d_dTAll.", __LINE__);
	cudasafe(cudaMemcpy(d_shared, shared, NV_PS * Nfitraw * sizeof(int), cudaMemcpyHostToDevice), "Failed cudaMemcpy on d_dTAll.", __LINE__);

	//calculate # of fitting parameters
	sumShared = 0;
	for (i = 0; i < 5; i++) {
		sumShared += shared[i];
	}
	noParameters = (5 * noChannels - sumShared * (noChannels - 1));


	//create output for parameters and CRLBs
	switch (fittype) {
	case 1:
		cudasafe(cudaMalloc((void**)&d_Parameters, (noParameters + 1) * Nfitraw * sizeof(float)), "Failed cudaMalloc on d_Parameters", __LINE__);
		cudasafe(cudaMemset(d_Parameters, 0, (noParameters + 1) * Nfitraw * sizeof(float)), "Failed cudaMemset on d_Parameters.", __LINE__);
		cudasafe(cudaMalloc((void**)&d_CRLBs, noParameters * Nfitraw * sizeof(float)), "Failed cudaMalloc on d_CRLBs.", __LINE__);
		cudasafe(cudaMemset(d_CRLBs, 0, noParameters * Nfitraw * sizeof(float)), "Failed cudaMemset on d_CRLBs.", __LINE__);
		plhs[0] = mxCreateNumericMatrix(Nfitraw, (noParameters + 1), mxSINGLE_CLASS, mxREAL);
		plhs[1] = mxCreateNumericMatrix(Nfitraw, noParameters, mxSINGLE_CLASS, mxREAL);
		break;
	case 2:
		cudasafe(cudaMalloc((void**)&d_Parameters, (noParameters + 1) * Nfitraw * sizeof(float)), "Failed cudaMalloc on d_Parameters", __LINE__);
		cudasafe(cudaMemset(d_Parameters, 0, (noParameters + 1) * Nfitraw * sizeof(float)), "Failed cudaMemset on d_Parameters.", __LINE__);
		cudasafe(cudaMalloc((void**)&d_CRLBs, noParameters * Nfitraw * sizeof(float)), "Failed cudaMalloc on d_CRLBs.", __LINE__);
		cudasafe(cudaMemset(d_CRLBs, 0, noParameters * Nfitraw * sizeof(float)), "Failed cudaMemset on d_CRLBs.", __LINE__);
		plhs[0] = mxCreateNumericMatrix(Nfitraw, (noParameters + 1), mxSINGLE_CLASS, mxREAL);
		plhs[1] = mxCreateNumericMatrix(Nfitraw, noParameters, mxSINGLE_CLASS, mxREAL);
		break;
	}



	cudasafe(cudaMalloc((void**)&d_LogLikelihood, Nfitraw * sizeof(float)), "Failed cudaMalloc on d_LogLikelihood.", __LINE__);
	plhs[2] = mxCreateNumericMatrix(Nfitraw, 1, mxSINGLE_CLASS, mxREAL);



	//check allocations
	if (plhs[0] == NULL || plhs[1] == NULL || plhs[2] == NULL)
		mexErrMsgIdAndTxt("GPUmleFit_LM_MultiChannel:NotEnoughMemory", "Could not allocate memory for results.\n");

	QueryPerformanceCounter(&stop);
	togpu = (double)(stop.QuadPart - start.QuadPart) / freq.QuadPart;
	QueryPerformanceCounter(&start);

	//setup kernel
	blockx = (int)ceil((float)Nfitraw / (float)BlockSize);
	threadx = BlockSize;

	dim3 dimBlock(threadx);
	dim3 dimGrid(blockx);

	if (silent == 0)
		mexPrintf("threadx: %d,blockx: %d\n", threadx, blockx);

	switch (fittype) {
	case 1:
		if (cameratype == 0) {
			kernel_MLEFit_sigma_EMCCD_multi_wrapper(dimGrid, dimBlock, d_data, PSFSigma, d_dTAll, (int)sz, iterations, noParameters, noChannels, d_Parameters, d_CRLBs, d_LogLikelihood, (const int)Nfitraw, d_shared);
		}
		else if (cameratype == 1) {
			kernel_MLEFit_sigma_sCMOS_multi_wrapper(dimGrid, dimBlock, d_data, PSFSigma, d_dTAll, (int)sz, iterations, noParameters, noChannels, d_Parameters, d_CRLBs, d_LogLikelihood, (const int)Nfitraw, d_shared, d_varim);
		}

		//extern void kernel_MLEFit_sigma_EMCCD_multi_wrapper(dim3 dimGrid, dim3 dimBlock, const float *d_data, const float PSFSigma, const int sz, const int iterations, 
		//const int NV, const int noChannels, float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits, const int*d_shared);


		break;
	case 2:
		/*if (initZ <0)
			initZ = spline_zsize/2.0f;*/
		if (cameratype == 0) {
			kernel_splineMLEFit_z_EMCCD_multi_wrapper(dimGrid, dimBlock, d_data, d_coeff, d_dTAll, spline_xsize, spline_ysize, spline_zsize, sz, iterations,
				noParameters, noChannels, d_Parameters, d_CRLBs, d_LogLikelihood, d_initZall, Nfitraw, d_shared);
		}
		else if (cameratype == 1) {
			kernel_splineMLEFit_z_sCMOS_multi_wrapper(dimGrid, dimBlock, d_data, d_coeff, d_dTAll, spline_xsize, spline_ysize, spline_zsize, sz, iterations,
				noParameters, noChannels, d_Parameters, d_CRLBs, d_LogLikelihood, d_initZall, Nfitraw, d_shared, d_varim);
		}

	}

	CUDAERROR("kernel", __LINE__);
	cudasafe(cudaDeviceSynchronize(), "sync failed", __LINE__);
	QueryPerformanceCounter(&stop);
	fit = (double)(stop.QuadPart - start.QuadPart) / freq.QuadPart;
	QueryPerformanceCounter(&start);



	cudasafe(cudaMemcpy((void*)mxGetData(plhs[0]), d_Parameters, (noParameters + 1) * Nfitraw * sizeof(float), cudaMemcpyDeviceToHost),
		"cudaMemcpy failed for d_Parameters.", __LINE__);
	cudasafe(cudaMemcpy((void*)mxGetData(plhs[1]), d_CRLBs, noParameters * Nfitraw * sizeof(float), cudaMemcpyDeviceToHost),
		"cudaMemcpy failed for d_CRLBs.", __LINE__);
	cudasafe(cudaMemcpy((void*)mxGetData(plhs[2]), d_LogLikelihood, Nfitraw * sizeof(float), cudaMemcpyDeviceToHost),
		"cudaMemcpy failed for d_LogLikelihood.", __LINE__);

	QueryPerformanceCounter(&stop);
	fromgpu = (double)(stop.QuadPart - start.QuadPart) / freq.QuadPart;
	QueryPerformanceCounter(&start);

	//cleanup
	cudasafe(cudaFree(d_data), "cudaFree failed on d_data.", __LINE__);
	cudasafe(cudaFree(d_Parameters), "cudaFree failed on d_Parameters.", __LINE__);
	cudasafe(cudaFree(d_CRLBs), "cudaFree failed on d_CRLBs.", __LINE__);
	cudasafe(cudaFree(d_LogLikelihood), "cudaFree failed on d_LogLikelihood.", __LINE__);
	cudasafe(cudaFree(d_coeff), "cudaFree failed on d_coeff.", __LINE__);
	cudasafe(cudaFree(d_shared), "cudaFree failed on d_coeff.", __LINE__);
	cudasafe(cudaFree(d_dTAll), "cudaFree failed on d_coeff.", __LINE__);
	cudasafe(cudaFree(d_varim), "cudaFree failed on d_varim.", __LINE__);
	//cudasafe(cudaFree(d_gainim),"cudaFree failed on d_varim.",__LINE__);

	//cudasafe(cudaDeviceReset(),"Error on cudaDeviceReset().",__LINE__); //release context so future cudaSetDevice calls work
	QueryPerformanceCounter(&stop);
	if (silent == 0)
	{
		mexPrintf("Memory copies to GPU %f seconds\n", togpu);
		mexPrintf("Actual fitting time %f Actual fits per second %f\n", fit,
			Nfitraw / fit);
		mexPrintf("Memory copies from GPU %f seconds\n", fromgpu);
		mexPrintf("Clean up from GPU %f seconds\n", (double)(stop.QuadPart - start.QuadPart) / freq.QuadPart);
	}
	return;
}
