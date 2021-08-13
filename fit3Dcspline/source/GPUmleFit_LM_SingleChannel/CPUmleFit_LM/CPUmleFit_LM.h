#include "definitions.h"
#ifndef CPUMLEFIT_LM_EMCCD_H
#define CPUMLEFIT_LM_EMCCD_H
 void kernel_MLEFit_LM_EMCCD(const int subregion,const float *d_data,const float PSFSigma, const int sz, const int iterations, 
	float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits);

 void kernel_MLEFit_LM_Sigma_EMCCD(const int subregion,const float *d_data,const float PSFSigma, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits);

 void kernel_MLEFit_LM_z_EMCCD(const int subregion,const float *d_data, const float PSFSigma_x, const float Ax, const float Ay, const float Bx, 
	const float By, const float gamma, const float d, const float PSFSigma_y, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits);

 void kernel_MLEFit_LM_sigmaxy_EMCCD(const int subregion,const float *d_data, const float PSFSigma, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits);

 void kernel_splineMLEFit_z_EMCCD(const int subregion,const float *d_data,const float *d_coeff, const int spline_xsize, const int spline_ysize, const int spline_zsize, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,float initZ,const int Nfits);


  void kernel_MLEFit_LM_sCMOS(const int subregion,const float *d_data,const float PSFSigma, const int sz, const int iterations, 
	float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits,const float *d_varim);

 void kernel_MLEFit_LM_Sigma_sCMOS(const int subregion,const float *d_data,const float PSFSigma, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits,const float *d_varim);

 void kernel_MLEFit_LM_z_sCMOS(const int subregion,const float *d_data, const float PSFSigma_x, const float Ax, const float Ay, const float Bx, 
	const float By, const float gamma, const float d, const float PSFSigma_y, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits,const float *d_varim);

 void kernel_MLEFit_LM_sigmaxy_sCMOS(const int subregion,const float *d_data, const float PSFSigma, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits,const float *d_varim);

 void kernel_splineMLEFit_z_sCMOS(const int subregion,const float *d_data,const float *d_coeff, const int spline_xsize, const int spline_ysize, const int spline_zsize, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,float initZ,const int Nfits,const float *d_varim);

#endif