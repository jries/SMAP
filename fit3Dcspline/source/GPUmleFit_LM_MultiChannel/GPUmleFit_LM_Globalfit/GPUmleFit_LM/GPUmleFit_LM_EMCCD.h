/*!
 * \file GPUmleFit_LM_EMCCD.h
 //author Yiming Li
//date 20170301
 * \brief Prototypes for the actual Cuda kernels.  
 */
#include "definitions.h"
#ifndef GPUMLEFIT_LM_EMCCD_H
#define GPUMLEFIT_LM_EMCCD_H

__global__ void kernel_splineMLEFit_z_EMCCD_multi(const float *d_data,const float *d_coeff, const float *d_dTAll, const int spline_xsize, const int spline_ysize, const int spline_zsize, const int sz, const int iterations, 
        const int NV, const int noChannels, float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,float *initZ, const int Nfits, const int *d_shared);

__global__ void kernel_MLEFit_sigma_EMCCD_multi(const float *d_data, const float PSFSigma, const float *d_dTAll,const int sz, const int iterations, 
	const int NV, const int noChannels, float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits, const int*d_shared);

#endif