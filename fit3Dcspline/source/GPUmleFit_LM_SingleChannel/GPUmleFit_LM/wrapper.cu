//file GPUmleFit_LM_sCMOS.cu
//author Yiming Li
//date 20170301

 //brief Wrap the Cuda kernel calls as standard external C functions.  This allows the kernels to be
 // called without doing anything special in the C code and simplifies building the code.
 //

//Terms of Use 
//
//This file is part of GPUmleFit_LM. 
//
//GPUmleFit_LM Fitter is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. 
//
//GPUmleFit_LM Fitter is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. 
//
//You should have received a copy of the GNU General Public License along with GPUmleFit_LM Fitter. If not, see <http://www.gnu.org/licenses/>. 
//
//Additional permission under GNU GPL version 3 section 7 

#include "definitions.h"
#include "GPUmleFit_LM_EMCCD.h"
#include "GPUgaussLib.h"
#include "GPUmleFit_LM_sCMOS.h"

//EMCCD wrapper

//*******************************************************************************************
extern void kernel_MLEFit_EMCCD_wrapper(dim3 dimGrid, dim3 dimBlock, const float *d_data, const float PSFSigma, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits) 
{
/*!
 *  \brief Basic maximum likelihood estimator fit based kernel
 *  \param dimGrid number of blocks 
 *  \param dimBlock number of threads per block
 *  \param d_data an array of subregions to be processed copied into video memory
 *  \param PSFSigma the sigma value to use for the point spread function
 *  \param sz nxn size of the subregion
 *  \param iterations maximum allowed iterations before aborting fitting
 *  \param d_Parameters pointer to result array of fitted parameters, x, y coords, etc.
 *  \param d_CRLBs pointer to result array of Cramer-Rao lower bound estimates 
 *  \param d_LogLikelihood pointer to result array of loglikelihood estimates of fitting
 *  \param Nfits number of subregions to fit
 */

	kernel_MLEFit_LM_EMCCD<<<dimGrid, dimBlock>>>(d_data, PSFSigma, sz, iterations, d_Parameters, d_CRLBs, d_LogLikelihood, Nfits);
	
}

//*******************************************************************************************
extern void kernel_MLEFit_sigma_EMCCD_wrapper(dim3 dimGrid, dim3 dimBlock, const float *d_data, const float PSFSigma, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits) 
{
/*!
 *  \brief Basic maximum likelihood estimator fit based kernel
 *  \param dimGrid number of blocks 
 *  \param dimBlock number of threads per block
 *  \param d_data an array of subregions to be processed copied into video memory
 *  \param PSFSigma the sigma value to use for the point spread function
 *  \param sz nxn size of the subregion
 *  \param iterations maximum allowed iterations before aborting fitting
 *  \param d_Parameters pointer to result array of fitted parameters, x, y coords, etc.
 *  \param d_CRLBs pointer to result array of Cramer-Rao lower bound estimates 
 *  \param d_LogLikelihood pointer to result array of loglikelihood estimates of fitting
 *  \param Nfits number of subregions to fit
 */
	kernel_MLEFit_LM_Sigma_EMCCD<<<dimGrid, dimBlock>>>(d_data, PSFSigma, sz, iterations, d_Parameters, d_CRLBs, d_LogLikelihood, Nfits);
	
}

//*******************************************************************************************
extern void kernel_MLEFit_z_EMCCD_wrapper(dim3 dimGrid, dim3 dimBlock, const float *d_data, const float PSFSigma_x, const float Ax, const float Ay, const float Bx, 
		const float By, const float gamma, const float d, const float PSFSigma_y, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits) 
{
/*!
 *  \brief Basic maximum likelihood estimator fit based kernel
 *  \param dimGrid number of blocks 
 *  \param dimBlock number of threads per block
 *  \param d_data an array of subregions to be processed copied into video memory
 *  \param PSFSigma_x the sigma value to use for the point spread function on the x axis
 *  \param Ax ???
 *  \param Ay ???
 *  \param Bx ???
 *  \param By ???
 *  \param gamma ???
 *  \param d ???
 *  \param PSFSigma_y the sigma value to use for the point spread function on the y axis
 *  \param sz nxn size of the subregion
 *  \param iterations maximum allowed iterations before aborting fitting
 *  \param d_Parameters pointer to result array of fitted parameters, x, y coords, etc.
 *  \param d_CRLBs pointer to result array of Cramer-Rao lower bound estimates 
 *  \param d_LogLikelihood pointer to result array of loglikelihood estimates of fitting
 *  \param Nfits number of subregions to fit
 */
	kernel_MLEFit_LM_z_EMCCD<<<dimGrid, dimBlock>>>(d_data, PSFSigma_x, Ax, Ay, Bx, By, gamma, d, PSFSigma_y, sz, iterations, d_Parameters, d_CRLBs, d_LogLikelihood, Nfits);
}

//*******************************************************************************************
extern void kernel_MLEFit_sigmaxy_EMCCD_wrapper(dim3 dimGrid, dim3 dimBlock, const float *d_data, const float PSFSigma, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits) 
{
	/*!
 *  \brief Basic maximum likelihood estimator fit based kernel
 *  \param dimGrid number of blocks 
 *  \param dimBlock number of threads per block
 *  \param d_data an array of subregions to be processed copied into video memory
 *  \param PSFSigma the sigma value to use for the point spread function
 *  \param sz nxn size of the subregion
 *  \param iterations maximum allowed iterations before aborting fitting
 *  \param d_Parameters pointer to result array of fitted parameters, x, y coords, etc.
 *  \param d_CRLBs pointer to result array of Cramer-Rao lower bound estimates 
 *  \param d_LogLikelihood pointer to result array of loglikelihood estimates of fitting
 *  \param Nfits number of subregions to fit
 */
	kernel_MLEFit_LM_sigmaxy_EMCCD<<<dimGrid, dimBlock>>>(d_data, PSFSigma, sz, iterations, d_Parameters, d_CRLBs, d_LogLikelihood, Nfits);
}


extern void kernel_splineMLEFit_z_EMCCD_wrapper(dim3 dimGrid, dim3 dimBlock, const float *d_data,const float *d_coeff, const int spline_xsize, const int spline_ysize, const int spline_zsize, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,float initZ, const int Nfits) 
/*!
 *  \brief Basic maximum likelihood estimator fit based kernel
 *  \param dimGrid number of blocks 
 *  \param dimBlock number of threads per block
 *  \param d_data an array of subregions to be processed copied into video memory
 *  \param d_coeff spline coefficient
 *  \param sz nxn size of the subregion
 *  \param iterations maximum allowed iterations before aborting fitting
 *  \param d_Parameters pointer to result array of fitted parameters, x, y coords, etc.
 *  \param d_CRLBs pointer to result array of Cramer-Rao lower bound estimates 
 *  \param d_LogLikelihood pointer to result array of loglikelihood estimates of fitting
 *  \param Nfits number of subregions to fit
 */
{
	kernel_splineMLEFit_z_EMCCD<<<dimGrid, dimBlock>>>(d_data, d_coeff, spline_xsize, spline_ysize, spline_zsize, sz, iterations, d_Parameters, d_CRLBs, d_LogLikelihood,initZ, Nfits);
}


//sCMOS wrapper
//*******************************************************************************************
extern void kernel_MLEFit_sCMOS_wrapper(dim3 dimGrid, dim3 dimBlock, const float *d_data, const float PSFSigma, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits, const float *d_varim) 
{
/*!
 *  \brief Basic maximum likelihood estimator fit based kernel
 *  \param dimGrid number of blocks 
 *  \param dimBlock number of threads per block
 *  \param d_data an array of subregions to be processed copied into video memory
 *  \param PSFSigma the sigma value to use for the point spread function
 *  \param sz nxn size of the subregion
 *  \param iterations maximum allowed iterations before aborting fitting
 *  \param d_Parameters pointer to result array of fitted parameters, x, y coords, etc.
 *  \param d_CRLBs pointer to result array of Cramer-Rao lower bound estimates 
 *  \param d_LogLikelihood pointer to result array of loglikelihood estimates of fitting
 *  \param Nfits number of subregions to fit
 *  \param d_varim variance map for sCMOS
 */

	kernel_MLEFit_LM_sCMOS<<<dimGrid, dimBlock>>>(d_data, PSFSigma, sz, iterations, d_Parameters, d_CRLBs, d_LogLikelihood, Nfits,d_varim);
}

//*******************************************************************************************
extern void kernel_MLEFit_sigma_sCMOS_wrapper(dim3 dimGrid, dim3 dimBlock, const float *d_data, const float PSFSigma, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits, const float *d_varim) 
{
/*!
 *  \brief Basic maximum likelihood estimator fit based kernel
 *  \param dimGrid number of blocks 
 *  \param dimBlock number of threads per block
 *  \param d_data an array of subregions to be processed copied into video memory
 *  \param PSFSigma the sigma value to use for the point spread function
 *  \param sz nxn size of the subregion
 *  \param iterations maximum allowed iterations before aborting fitting
 *  \param d_Parameters pointer to result array of fitted parameters, x, y coords, etc.
 *  \param d_CRLBs pointer to result array of Cramer-Rao lower bound estimates 
 *  \param d_LogLikelihood pointer to result array of loglikelihood estimates of fitting
 *  \param Nfits number of subregions to fit
 *  \param d_varim variance map for sCMOS
 */
	kernel_MLEFit_LM_Sigma_sCMOS<<<dimGrid, dimBlock>>>(d_data, PSFSigma, sz, iterations, d_Parameters, d_CRLBs, d_LogLikelihood, Nfits,d_varim);
}

//*******************************************************************************************
extern void kernel_MLEFit_z_sCMOS_wrapper(dim3 dimGrid, dim3 dimBlock, const float *d_data, const float PSFSigma_x, const float Ax, const float Ay, const float Bx, 
		const float By, const float gamma, const float d, const float PSFSigma_y, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits, const float *d_varim) 
{
/*!
 *  \brief Basic maximum likelihood estimator fit based kernel
 *  \param dimGrid number of blocks 
 *  \param dimBlock number of threads per block
 *  \param d_data an array of subregions to be processed copied into video memory
 *  \param PSFSigma_x the sigma value to use for the point spread function on the x axis
 *  \param Ax ???
 *  \param Ay ???
 *  \param Bx ???
 *  \param By ???
 *  \param gamma ???
 *  \param d ???
 *  \param PSFSigma_y the sigma value to use for the point spread function on the y axis
 *  \param sz nxn size of the subregion
 *  \param iterations maximum allowed iterations before aborting fitting
 *  \param d_Parameters pointer to result array of fitted parameters, x, y coords, etc.
 *  \param d_CRLBs pointer to result array of Cramer-Rao lower bound estimates 
 *  \param d_LogLikelihood pointer to result array of loglikelihood estimates of fitting
 *  \param Nfits number of subregions to fit
 *  \param d_varim variance map for sCMOS
 */
	kernel_MLEFit_LM_z_sCMOS<<<dimGrid, dimBlock>>>(d_data, PSFSigma_x, Ax, Ay, Bx, By, gamma, d, PSFSigma_y, sz, iterations, d_Parameters, d_CRLBs, d_LogLikelihood, Nfits,d_varim);
}

//*******************************************************************************************
extern void kernel_MLEFit_sigmaxy_sCMOS_wrapper(dim3 dimGrid, dim3 dimBlock, const float *d_data, const float PSFSigma, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits, const float *d_varim) 
{
	/*!
 *  \brief Basic maximum likelihood estimator fit based kernel
 *  \param dimGrid number of blocks 
 *  \param dimBlock number of threads per block
 *  \param d_data an array of subregions to be processed copied into video memory
 *  \param PSFSigma the sigma value to use for the point spread function
 *  \param sz nxn size of the subregion
 *  \param iterations maximum allowed iterations before aborting fitting
 *  \param d_Parameters pointer to result array of fitted parameters, x, y coords, etc.
 *  \param d_CRLBs pointer to result array of Cramer-Rao lower bound estimates 
 *  \param d_LogLikelihood pointer to result array of loglikelihood estimates of fitting
 *  \param Nfits number of subregions to fit
 *  \param d_varim variance map for sCMOS
 */
	kernel_MLEFit_LM_sigmaxy_sCMOS<<<dimGrid, dimBlock>>>(d_data, PSFSigma, sz, iterations, d_Parameters, d_CRLBs, d_LogLikelihood, Nfits,d_varim);
}


extern void kernel_splineMLEFit_z_sCMOS_wrapper(dim3 dimGrid, dim3 dimBlock, const float *d_data,const float *d_coeff, const int spline_xsize, const int spline_ysize, const int spline_zsize, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,float initZ, const int Nfits, const float *d_varim) 
{
	/*!
 *  \brief Basic maximum likelihood estimator fit based kernel
 *  \param dimGrid number of blocks 
 *  \param dimBlock number of threads per block
 *  \param d_data an array of subregions to be processed copied into video memory
 *  \param d_coeff spline coefficient
 *  \param sz nxn size of the subregion
 *  \param iterations maximum allowed iterations before aborting fitting
 *  \param d_Parameters pointer to result array of fitted parameters, x, y coords, etc.
 *  \param d_CRLBs pointer to result array of Cramer-Rao lower bound estimates 
 *  \param d_LogLikelihood pointer to result array of loglikelihood estimates of fitting
 *  \param Nfits number of subregions to fit
 *  \param d_varim variance map for sCMOS
 */
	kernel_splineMLEFit_z_sCMOS<<<dimGrid, dimBlock>>>(d_data, d_coeff, spline_xsize, spline_ysize, spline_zsize, sz, iterations, d_Parameters, d_CRLBs, d_LogLikelihood,initZ, Nfits, d_varim);
}
