
#include "definitions.h"
#include "GPUmleFit_LM_EMCCD.h"
#include "GPUgaussLib.h"
#include "GPUmleFit_LM_sCMOS.h"
//#include "GPUgaussLib.cuh"
//#include "GPUgaussMLE_LM.h"
//#include "GPUgaussMLE_sCMOS.h"
//#include "GPUsplineMLE.h"

//EMCCD wrapper

extern void kernel_splineMLEFit_z_EMCCD_multi_wrapper(dim3 dimGrid, dim3 dimBlock, const float *d_data,const float *d_coeff, const float *d_dTAll, const int spline_xsize, const int spline_ysize, const int spline_zsize, const int sz, const int iterations, 
        const int NV, const int noChannels, float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,float * d_initZall, const int Nfits, const int *d_shared)
{
	kernel_splineMLEFit_z_EMCCD_multi<<<dimGrid, dimBlock>>>(d_data,d_coeff, d_dTAll, spline_xsize, spline_ysize, spline_zsize,  sz, iterations, NV, noChannels, d_Parameters, d_CRLBs, d_LogLikelihood, d_initZall, Nfits, d_shared);
}


extern void kernel_MLEFit_sigma_EMCCD_multi_wrapper(dim3 dimGrid, dim3 dimBlock, const float *d_data, const float PSFSigma, const float *d_dTAll,const int sz, const int iterations, 
	const int NV, const int noChannels, float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits, const int*d_shared)
{
	kernel_MLEFit_sigma_EMCCD_multi<<<dimGrid, dimBlock>>>(d_data,PSFSigma,  d_dTAll,sz, iterations, NV, noChannels, d_Parameters, d_CRLBs, d_LogLikelihood, Nfits, d_shared);
}

//sCMOS wrapper
//*******************************************************************************************
extern void kernel_splineMLEFit_z_sCMOS_multi_wrapper(dim3 dimGrid, dim3 dimBlock, const float* d_data, const float* d_coeff, const float* d_dTAll, const int spline_xsize, const int spline_ysize, const int spline_zsize, const int sz, const int iterations,
	const int NV, const int noChannels, float* d_Parameters, float* d_CRLBs, float* d_LogLikelihood, float* d_initZall, const int Nfits, const int* d_shared, const float* d_varim)
{
	kernel_splineMLEFit_z_sCMOS_multi << <dimGrid, dimBlock >> > (d_data, d_coeff, d_dTAll, spline_xsize, spline_ysize, spline_zsize, sz, iterations, NV, noChannels, d_Parameters, d_CRLBs, d_LogLikelihood, d_initZall, Nfits, d_shared, d_varim);
}


extern void kernel_MLEFit_sigma_sCMOS_multi_wrapper(dim3 dimGrid, dim3 dimBlock, const float* d_data, const float PSFSigma, const float* d_dTAll, const int sz, const int iterations,
	const int NV, const int noChannels, float* d_Parameters, float* d_CRLBs, float* d_LogLikelihood, const int Nfits, const int* d_shared, const float* d_varim)
{
	kernel_MLEFit_sigma_sCMOS_multi << <dimGrid, dimBlock >> > (d_data, PSFSigma, d_dTAll, sz, iterations, NV, noChannels, d_Parameters, d_CRLBs, d_LogLikelihood, Nfits, d_shared, d_varim);
}