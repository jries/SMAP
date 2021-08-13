#include "definitions.h"
#ifndef CPUMLEFIT_LM_EMCCD_H
#define CPUMLEFIT_LM_EMCCD_H


void kernel_splineMLEFit_z_EMCCD_multi(const int subregion, const float *d_data,const float *d_coeff, const float *d_dTAll, const int spline_xsize, const int spline_ysize, const int spline_zsize, const int sz, const int iterations, 
        const int NV, const int noChannels, float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,float *initZ, const int Nfits, const int *d_shared);

void kernel_MLEFit_sigma_EMCCD_multi(const int subregion,const float *d_data, const float PSFSigma, const float *d_dTAll,const int sz, const int iterations, const int NV, const int noChannels, float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits, const int*d_shared);

void kernel_splineMLEFit_z_sCMOS_multi(const int subregion, const float* d_data, const float* d_coeff, const float* d_dTAll, const int spline_xsize, const int spline_ysize, const int spline_zsize, const int sz, const int iterations,
	const int NV, const int noChannels, float* d_Parameters, float* d_CRLBs, float* d_LogLikelihood, float* initZ, const int Nfits, const int* d_shared, const float* d_varim);

void kernel_MLEFit_sigma_sCMOS_multi(const int subregion, const float* d_data, const float PSFSigma, const float* d_dTAll, const int sz, const int iterations, const int NV, const int noChannels, float* d_Parameters, float* d_CRLBs, float* d_LogLikelihood, const int Nfits, const int* d_shared, const float* d_varim);



#endif