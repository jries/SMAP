//Copyright (c) 2017 Ries Lab, European Molecular Biology Laboratory, Heidelberg.
//author: Yiming Li
//email: yiming.li@embl.de
//date: 2017.07.24

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



/*!
 * \file CPUmleFit_LM.cpp
 * \brief This contains the definitions for all the fitting mode.  
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mex.h"
#include "definitions.h"
#include "MatInvLib.h"
#include <math.h>
#include "CPUsplineLib.h"
#include "CPUgaussLib.h"
#include "CPUmleFit_LM.h"

//*******************************************************************************************
//theta is: {x,y,N,bg}
 void kernel_MLEFit_LM_EMCCD(const int subregion, const float *d_data,const float PSFSigma, const int sz, const int iterations, 
	float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits){
		/*! 
	 * \brief basic MLE fitting kernel.  No additional parameters are computed.
	 * \param d_data array of subregions to fit copied to GPU
	 * \param PSFSigma sigma of the point spread function
	 * \param sz nxn size of the subregion to fit
	 * \param iterations number of iterations for solution to converge
	 * \param d_Parameters array of fitting parameters to return for each subregion
	 * \param d_CRLBs array of Cramer-Rao lower bound estimates to return for each subregion
	 * \param d_LogLikelihood array of loglikelihood estimates to return for each subregion
	 * \param Nfits number of subregions to fit
	 */

		const int NV=NV_P;
		float M[NV*NV],Diag[NV], Minv[NV*NV];
		//int tx = threadIdx.x;
		//int bx = blockIdx.x;
		//int BlockSize = blockDim.x;
		int ii, jj, kk, ll, l, m, i;


		float model, data;
		float Div;

		float newTheta[NV],oldTheta[NV];
		float newLambda = INIT_LAMBDA, oldLambda = INIT_LAMBDA, mu;
	    float newUpdate[NV] = {1e13, 1e13, 1e13, 1e13},oldUpdate[NV] = {1e13, 1e13, 1e13, 1e13};
		float maxJump[NV]={1.0,1.0,100,20};
		float newDudt[NV] ={0};

		float newErr = 1e12, oldErr = 1e13;

		float jacobian[NV]={0};
		float hessian[NV*NV]={0};
		float t1,t2;

		float Nmax;
		int errFlag=0;
		float L[NV*NV] = {0}, U[NV*NV] = {0};


		//Prevent read/write past end of array
		if (subregion>=Nfits) return;

		for (ii=0;ii<NV*NV;ii++)M[ii]=0;
		for (ii=0;ii<NV*NV;ii++)Minv[ii]=0;

		//copy in data
		const float *s_data = d_data+(sz*sz*subregion);

		//initial values
		kernel_CenterofMass2D(sz, s_data, &newTheta[0], &newTheta[1]);
		kernel_GaussFMaxMin2D(sz, PSFSigma, s_data, &Nmax, &newTheta[3]);
		newTheta[2]=max(0.0, (Nmax-newTheta[3])*2*PI*PSFSigma*PSFSigma);
		newTheta[3] = max(newTheta[3],0.01);

		maxJump[2]=max(newTheta[2],maxJump[2]);

		maxJump[3]=max(newTheta[3],maxJump[3]);

		for (ii=0;ii<NV;ii++)oldTheta[ii]=newTheta[ii];

		//updateFitValues
		newErr = 0;
		memset(jacobian,0,NV*sizeof(float));
		memset(hessian,0,NV*NV*sizeof(float));
		for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
			kernel_DerivativeGauss2D(ii,jj,PSFSigma,newTheta,newDudt,&model);
			data=s_data[sz*jj+ii];

			if (data>0)
				newErr = newErr + 2*((model-data)-data*log(model/data));
			else
			{
				newErr = newErr + 2*model;
				data = 0;
			}

			t1 = 1-data/model;
			for (l=0;l<NV;l++){
				jacobian[l]+=t1*newDudt[l];
			}

			t2 = data/pow(model,2);
			for (l=0;l<NV;l++) for(m=l;m<NV;m++) {
				hessian[l*NV+m] +=t2*newDudt[l]*newDudt[m];
				hessian[m*NV+l] = hessian[l*NV+m];
			}
		}

		for (kk=0;kk<iterations;kk++) {//main iterative loop

			if(fabs((newErr-oldErr)/newErr)<TOLERANCE){
				//CONVERGED;
				break;
			}
			else{
				if(newErr>ACCEPTANCE*oldErr){
					//copy Fitdata

					for (i=0;i<NV;i++){
						newTheta[i]=oldTheta[i];
						newUpdate[i]=oldUpdate[i];
					}
					newLambda = oldLambda;
					newErr = oldErr;
					mu = max( (1 + newLambda*SCALE_UP)/(1 + newLambda),1.3f);         
					newLambda = SCALE_UP*newLambda;

				}
				else if(newErr<oldErr&&errFlag==0){
					newLambda = SCALE_DOWN*newLambda;
				    mu = 1+newLambda;
				}


				for (i=0;i<NV;i++){
					hessian[i*NV+i]=hessian[i*NV+i]*mu;
				}
				memset(L,0,NV*sizeof(float));
				memset(U,0,NV*sizeof(float));
				errFlag = kernel_cholesky(hessian,NV,L,U);
				if (errFlag ==0){
					for (i=0;i<NV;i++){
						oldTheta[i]=newTheta[i];
						oldUpdate[i] = newUpdate[i];
					}
					oldLambda = newLambda;
					oldErr=newErr;

					kernel_luEvaluate(L,U,jacobian,NV,newUpdate);	
					
					//updateFitParameters
					for (ll=0;ll<NV;ll++){
						if (newUpdate[ll]/oldUpdate[ll]< -0.5f){
							maxJump[ll] = maxJump[ll]*0.5;
						}
					    newUpdate[ll] = newUpdate[ll]/(1+fabs(newUpdate[ll]/maxJump[ll]));
						newTheta[ll] = newTheta[ll]-newUpdate[ll];
					}
					//restrict range
					newTheta[0] = max(newTheta[0],(float(sz)-1)/2-sz/4.0);
					newTheta[0] = min(newTheta[0],(float(sz)-1)/2+sz/4.0);
					newTheta[1] = max(newTheta[1],(float(sz)-1)/2-sz/4.0);
					newTheta[1] = min(newTheta[1],(float(sz)-1)/2+sz/4.0);
					newTheta[2] = max(newTheta[2],1.0);
					newTheta[3] = max(newTheta[3],0.01);
					


					newErr = 0;
					memset(jacobian,0,NV*sizeof(float));
					memset(hessian,0,NV*NV*sizeof(float));
					for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
						//calculating derivatives
						kernel_DerivativeGauss2D(ii,jj,PSFSigma,newTheta,newDudt,&model);
						data=s_data[sz*jj+ii];			

						if (data>0)
							newErr = newErr + 2*((model-data)-data*log(model/data));
						else
						{
							newErr = newErr + 2*model;
							data = 0;
						}

						t1 = 1-data/model;
						for (l=0;l<NV;l++){
							jacobian[l]+=t1*newDudt[l];
						}

						t2 = data/pow(model,2);
						for (l=0;l<NV;l++) for(m=l;m<NV;m++) {
							hessian[l*NV+m] +=t2*newDudt[l]*newDudt[m];
							hessian[m*NV+l] = hessian[l*NV+m];
						}
					}
				}
				else
				{
					mu = max( (1 + newLambda*SCALE_UP)/(1 + newLambda),1.3f);         
					newLambda = SCALE_UP*newLambda;
				}
			}
		}
		//output iteration
		d_Parameters[Nfits*NV+subregion]=kk;
		// Calculating the CRLB and LogLikelihood
		Div=0.0f;
		for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
			kernel_DerivativeGauss2D(ii,jj,PSFSigma,newTheta,newDudt,&model);
			data=s_data[sz*jj+ii];

			//Building the Fisher Information Matrix
			for (kk=0;kk<NV;kk++)for (ll=kk;ll<NV;ll++){
				M[kk*NV+ll]+= newDudt[ll]*newDudt[kk]/model;
				M[ll*NV+kk]=M[kk*NV+ll];
			}

			//LogLikelyhood
			if (model>0)
				if (data>0)Div+=data*log(model)-model-data*log(data)+data;
				else
					Div+=-model;
		}

		// Matrix inverse (CRLB=F^-1) and output assigments
		kernel_MatInvN(M, Minv, Diag, NV);
		//write to global arrays
		for (kk=0;kk<NV;kk++) d_Parameters[Nfits*kk+subregion]=newTheta[kk];
		for (kk=0;kk<NV;kk++) d_CRLBs[Nfits*kk+subregion]=Diag[kk];
		d_LogLikelihood[subregion] = Div;

		return;
}

//*********************************************************************************************************************************************

 void kernel_MLEFit_LM_Sigma_EMCCD(const int subregion, const float *d_data,const float PSFSigma, const int sz, const int iterations, 
	float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits){
		/*! 
	 * \brief basic MLE fitting kernel.  No additional parameters are computed.
	 * \param d_data array of subregions to fit copied to GPU
	 * \param PSFSigma sigma of the point spread function
	 * \param sz nxn size of the subregion to fit
	 * \param iterations number of iterations for solution to converge
	 * \param d_Parameters array of fitting parameters to return for each subregion
	 * \param d_CRLBs array of Cramer-Rao lower bound estimates to return for each subregion
	 * \param d_LogLikelihood array of loglikelihood estimates to return for each subregion
	 * \param Nfits number of subregions to fit
	 */

		const int NV=NV_PS;
		float M[NV*NV],Diag[NV], Minv[NV*NV];
		/*int tx = threadIdx.x;
		int bx = blockIdx.x;
		int BlockSize = blockDim.x;*/
		int ii, jj, kk, ll, l, m, i;


		float model, data;
		float Div;

		float newTheta[NV],oldTheta[NV];
		float newLambda = INIT_LAMBDA, oldLambda = INIT_LAMBDA, mu;
		float newUpdate[NV] = {1e13, 1e13, 1e13, 1e13, 1e13},oldUpdate[NV] = {1e13, 1e13, 1e13, 1e13, 1e13};
		float maxJump[NV]={1.0,1.0,100,20,0.5};
		float newDudt[NV] ={0};

		float newErr = 1e12, oldErr = 1e13;

		float jacobian[NV]={0};
		float hessian[NV*NV]={0};
		float t1,t2;

		float Nmax;
		int errFlag=0;
		float L[NV*NV] = {0}, U[NV*NV] = {0};


		//Prevent read/write past end of array
		if (subregion>=Nfits) return;

		for (ii=0;ii<NV*NV;ii++)M[ii]=0;
		for (ii=0;ii<NV*NV;ii++)Minv[ii]=0;

		//copy in data
		const float *s_data = d_data+(sz*sz*subregion);

		//initial values
		kernel_CenterofMass2D(sz, s_data, &newTheta[0], &newTheta[1]);
		kernel_GaussFMaxMin2D(sz, PSFSigma, s_data, &Nmax, &newTheta[3]);
		newTheta[2]=max(0.0, (Nmax-newTheta[3])*2*PI*PSFSigma*PSFSigma);
		newTheta[3] = max(newTheta[3],0.01);
		newTheta[4]=PSFSigma;

		maxJump[2]=max(newTheta[2],maxJump[2]);

		maxJump[3]=max(newTheta[3],maxJump[3]);;


		for (ii=0;ii<NV;ii++)oldTheta[ii]=newTheta[ii];

		//updateFitValues
		newErr = 0;
		memset(jacobian,0,NV*sizeof(float));
		memset(hessian,0,NV*NV*sizeof(float));
		for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
			kernel_DerivativeGauss2D_sigma(ii,jj,newTheta,newDudt,&model);
			data=s_data[sz*jj+ii];

			if (data>0)
				newErr = newErr + 2*((model-data)-data*log(model/data));
			else
			{
				newErr = newErr + 2*model;
				data = 0;
			}

			t1 = 1-data/model;
			for (l=0;l<NV;l++){
				jacobian[l]+=t1*newDudt[l];
			}

			t2 = data/pow(model,2);
			for (l=0;l<NV;l++) for(m=l;m<NV;m++) {
				hessian[l*NV+m] +=t2*newDudt[l]*newDudt[m];
				hessian[m*NV+l] = hessian[l*NV+m];
			}
		}

		for (kk=0;kk<iterations;kk++) {//main iterative loop

			if(fabs((newErr-oldErr)/newErr)<TOLERANCE){
				//CONVERGED;
				break;
			}
			else{
				if(newErr>ACCEPTANCE*oldErr){
					//copy Fitdata

					for (i=0;i<NV;i++){
						newTheta[i]=oldTheta[i];
						newUpdate[i]=oldUpdate[i];
					}
					newLambda = oldLambda;
					newErr = oldErr;
					mu = max( (1 + newLambda*SCALE_UP)/(1 + newLambda),1.3f);         
					newLambda = SCALE_UP*newLambda;
				}
				else if(newErr<oldErr&&errFlag==0){
					newLambda = SCALE_DOWN*newLambda;
				    mu = 1+newLambda;
				}


				for (i=0;i<NV;i++){
					hessian[i*NV+i]=hessian[i*NV+i]*mu;
				}
				memset(L,0,NV*sizeof(float));
				memset(U,0,NV*sizeof(float));
				errFlag = kernel_cholesky(hessian,NV,L,U);
				if (errFlag ==0){
					for (i=0;i<NV;i++){
						oldTheta[i]=newTheta[i];
						oldUpdate[i] = newUpdate[i];
					}
					oldLambda = newLambda;
					oldErr=newErr;

					kernel_luEvaluate(L,U,jacobian,NV,newUpdate);	
					
					//updateFitParameters
					for (ll=0;ll<NV;ll++){
						if (newUpdate[ll]/oldUpdate[ll]< -0.5f){
							maxJump[ll] = maxJump[ll]*0.5;
						}
					    newUpdate[ll] = newUpdate[ll]/(1+fabs(newUpdate[ll]/maxJump[ll]));
						newTheta[ll] = newTheta[ll]-newUpdate[ll];
					}
					//restrict range
					newTheta[0] = max(newTheta[0],(float(sz)-1)/2-sz/4.0);
					newTheta[0] = min(newTheta[0],(float(sz)-1)/2+sz/4.0);
					newTheta[1] = max(newTheta[1],(float(sz)-1)/2-sz/4.0);
					newTheta[1] = min(newTheta[1],(float(sz)-1)/2+sz/4.0);
					newTheta[2] = max(newTheta[2],1.0);
					newTheta[3] = max(newTheta[3],0.01);
					newTheta[4] = max(newTheta[4],0.0);
					newTheta[4] = min(newTheta[4],sz/2.0f);


					newErr = 0;
					memset(jacobian,0,NV*sizeof(float));
					memset(hessian,0,NV*NV*sizeof(float));
					for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
						//calculating derivatives
						kernel_DerivativeGauss2D_sigma(ii,jj,newTheta,newDudt,&model);
						data=s_data[sz*jj+ii];			

						if (data>0)
							newErr = newErr + 2*((model-data)-data*log(model/data));
						else
						{
							newErr = newErr + 2*model;
							data = 0;
						}

						t1 = 1-data/model;
						for (l=0;l<NV;l++){
							jacobian[l]+=t1*newDudt[l];
						}

						t2 = data/pow(model,2);
						for (l=0;l<NV;l++) for(m=l;m<NV;m++) {
							hessian[l*NV+m] +=t2*newDudt[l]*newDudt[m];
							hessian[m*NV+l] = hessian[l*NV+m];
						}
					}
				}
				else
				{
					mu = max( (1 + newLambda*SCALE_UP)/(1 + newLambda),1.3f);         
					newLambda = SCALE_UP*newLambda;
				}
			}
		}
		//output iteration
		d_Parameters[Nfits*NV+subregion]=kk;
		// Calculating the CRLB and LogLikelihood
		Div=0.0f;
		for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
			
			kernel_DerivativeGauss2D_sigma(ii,jj,newTheta,newDudt,&model);
			data=s_data[sz*jj+ii];

			//Building the Fisher Information Matrix
			for (kk=0;kk<NV;kk++)for (ll=kk;ll<NV;ll++){
				M[kk*NV+ll]+= newDudt[ll]*newDudt[kk]/model;
				M[ll*NV+kk]=M[kk*NV+ll];
			}

			//LogLikelyhood
			if (model>0)
				if (data>0)Div+=data*log(model)-model-data*log(data)+data;
				else
					Div+=-model;
		}

		// Matrix inverse (CRLB=F^-1) and output assigments
		kernel_MatInvN(M, Minv, Diag, NV);
		//write to global arrays
		for (kk=0;kk<NV;kk++) d_Parameters[Nfits*kk+subregion]=newTheta[kk];
		for (kk=0;kk<NV;kk++) d_CRLBs[Nfits*kk+subregion]=Diag[kk];
		d_LogLikelihood[subregion] = Div;

		return;
}

//*********************************************************************************************************************************************

 void kernel_MLEFit_LM_z_EMCCD(const int subregion,const float *d_data, const float PSFSigma_x, const float Ax, const float Ay, const float Bx, 
	const float By, const float gamma, const float d, const float PSFSigma_y, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits){
			/*! 
	 * \brief basic MLE fitting kernel.  No additional parameters are computed.
	 * \param d_data array of subregions to fit copied to GPU
	 * \param PSFSigma_x sigma of the point spread function on the x axis
	 * \param Ax ???
	 * \param Ay ???
	 * \param Bx ???
	 * \param By ???
	 * \param gamma ???
	 * \param d ???
	 * \param PSFSigma_y sigma of the point spread function on the y axis
	 * \param sz nxn size of the subregion to fit
	 * \param iterations number of iterations for solution to converge
	 * \param d_Parameters array of fitting parameters to return for each subregion
	 * \param d_CRLBs array of Cramer-Rao lower bound estimates to return for each subregion
	 * \param d_LogLikelihood array of loglikelihood estimates to return for each subregion
	 * \param Nfits number of subregions to fit
	 */

		const int NV=NV_PZ;
		float M[NV*NV],Diag[NV], Minv[NV*NV];
		//int tx = threadIdx.x;
		//int bx = blockIdx.x;
		//int BlockSize = blockDim.x;
		int ii, jj, kk, ll, l, m, i;


		float model, data;
		float Div;
		float PSFy, PSFx;

		float newTheta[NV],oldTheta[NV];
		float newLambda = INIT_LAMBDA, oldLambda = INIT_LAMBDA, mu;
		float newUpdate[NV] = {1e13, 1e13, 1e13, 1e13, 1e13},oldUpdate[NV] = {1e13, 1e13, 1e13, 1e13, 1e13};
		float maxJump[NV]={1.0,1.0,100,20,2};
		float newDudt[NV] ={0};

		float newErr = 1e12, oldErr = 1e13;

		float jacobian[NV]={0};
		float hessian[NV*NV]={0};
		float t1,t2;

		float Nmax;
//		float temp;
		int errFlag=0;
		float L[NV*NV] = {0}, U[NV*NV] = {0};


		//Prevent read/write past end of array
		if (subregion>=Nfits) return;

		for (ii=0;ii<NV*NV;ii++)M[ii]=0;
		for (ii=0;ii<NV*NV;ii++)Minv[ii]=0;

		//copy in data
		const float *s_data = d_data+(sz*sz*subregion);

		//initial values
		kernel_CenterofMass2D(sz, s_data, &newTheta[0], &newTheta[1]);
		kernel_GaussFMaxMin2D(sz, PSFSigma_x, s_data, &Nmax, &newTheta[3]);
		newTheta[2]=max(0.0, (Nmax-newTheta[3])*2*PI*PSFSigma_x*PSFSigma_y*sqrt(2.0f));
		newTheta[3] = max(newTheta[3],0.01);
		newTheta[4]=0;

		maxJump[2]=max(newTheta[2],maxJump[2]);

		maxJump[3]=max(newTheta[3],maxJump[3]);

		for (ii=0;ii<NV;ii++)oldTheta[ii]=newTheta[ii];

		//updateFitValues
		newErr = 0;
		memset(jacobian,0,NV*sizeof(float));
		memset(hessian,0,NV*NV*sizeof(float));
		for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
			 kernel_DerivativeIntGauss2Dz(ii, jj, newTheta, PSFSigma_x,PSFSigma_y, Ax,Ay,Bx,By, gamma, d, &PSFx, &PSFy, newDudt, NULL,&model);
			data=s_data[sz*jj+ii];

			if (data>0)
				newErr = newErr + 2*((model-data)-data*log(model/data));
			else
			{
				newErr = newErr + 2*model;
				data = 0;
			}

			t1 = 1-data/model;
			for (l=0;l<NV;l++){
				jacobian[l]+=t1*newDudt[l];
			}

			t2 = data/pow(model,2);
			for (l=0;l<NV;l++) for(m=l;m<NV;m++) {
				hessian[l*NV+m] +=t2*newDudt[l]*newDudt[m];
				hessian[m*NV+l] = hessian[l*NV+m];
			}
		}

		for (kk=0;kk<iterations;kk++) {//main iterative loop

			if(fabs((newErr-oldErr)/newErr)<TOLERANCE){
				//CONVERGED;
				break;
			}
			else{
				if(newErr>ACCEPTANCE*oldErr){
					//copy Fitdata

					for (i=0;i<NV;i++){
						newTheta[i]=oldTheta[i];
						newUpdate[i]=oldUpdate[i];
					}
					newLambda = oldLambda;
					newErr = oldErr;
					mu = max( (1 + newLambda*SCALE_UP)/(1 + newLambda),1.3f);         
					newLambda = SCALE_UP*newLambda;
				}
				else if(newErr<oldErr&&errFlag==0){
					newLambda = SCALE_DOWN*newLambda;
				    mu = 1+newLambda;
				}


				for (i=0;i<NV;i++){
					hessian[i*NV+i]=hessian[i*NV+i]*mu;
				}
				memset(L,0,NV*sizeof(float));
				memset(U,0,NV*sizeof(float));
				errFlag = kernel_cholesky(hessian,NV,L,U);
				if (errFlag ==0){
					for (i=0;i<NV;i++){
						oldTheta[i]=newTheta[i];
						oldUpdate[i] = newUpdate[i];
					}
					oldLambda = newLambda;
					oldErr=newErr;

					kernel_luEvaluate(L,U,jacobian,NV,newUpdate);	
					
					//updateFitParameters
					for (ll=0;ll<NV;ll++){
						if (newUpdate[ll]/oldUpdate[ll]< -0.5f){
							maxJump[ll] = maxJump[ll]*0.5;
						}
					    newUpdate[ll] = newUpdate[ll]/(1+fabs(newUpdate[ll]/maxJump[ll]));
						newTheta[ll] = newTheta[ll]-newUpdate[ll];
					}

					//restrict range
					newTheta[0] = max(newTheta[0],(float(sz)-1)/2-sz/4.0);
					newTheta[0] = min(newTheta[0],(float(sz)-1)/2+sz/4.0);
					newTheta[1] = max(newTheta[1],(float(sz)-1)/2-sz/4.0);
					newTheta[1] = min(newTheta[1],(float(sz)-1)/2+sz/4.0);
					newTheta[2] = max(newTheta[2],1.0);
					newTheta[3] = max(newTheta[3],0.01);


					newErr = 0;
					memset(jacobian,0,NV*sizeof(float));
					memset(hessian,0,NV*NV*sizeof(float));
					for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
						//calculating derivatives
						kernel_DerivativeIntGauss2Dz(ii, jj, newTheta, PSFSigma_x,PSFSigma_y, Ax,Ay,Bx,By, gamma, d, &PSFx, &PSFy, newDudt, NULL,&model);
						data=s_data[sz*jj+ii];			

						if (data>0)
							newErr = newErr + 2*((model-data)-data*log(model/data));
						else
						{
							newErr = newErr + 2*model;
							data = 0;
						}

						t1 = 1-data/model;
						for (l=0;l<NV;l++){
							jacobian[l]+=t1*newDudt[l];
						}

						t2 = data/pow(model,2);
						for (l=0;l<NV;l++) for(m=l;m<NV;m++) {
							hessian[l*NV+m] +=t2*newDudt[l]*newDudt[m];
							hessian[m*NV+l] = hessian[l*NV+m];
						}
					}
				}
				else
				{
					mu = max( (1 + newLambda*SCALE_UP)/(1 + newLambda),1.3f);         
					newLambda = SCALE_UP*newLambda;
				}
			}
		}
		//output iteration
		d_Parameters[Nfits*NV+subregion]=kk;
		// Calculating the CRLB and LogLikelihood
		Div=0.0f;
		for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
		    kernel_DerivativeIntGauss2Dz(ii, jj, newTheta, PSFSigma_x,PSFSigma_y, Ax,Ay,Bx,By, gamma, d, &PSFx, &PSFy, newDudt, NULL,&model);
			data=s_data[sz*jj+ii];

			//Building the Fisher Information Matrix
			for (kk=0;kk<NV;kk++)for (ll=kk;ll<NV;ll++){
				M[kk*NV+ll]+= newDudt[ll]*newDudt[kk]/model;
				M[ll*NV+kk]=M[kk*NV+ll];
			}

			//LogLikelyhood
			if (model>0)
				if (data>0)Div+=data*log(model)-model-data*log(data)+data;
				else
					Div+=-model;
		}

		// Matrix inverse (CRLB=F^-1) and output assigments
		kernel_MatInvN(M, Minv, Diag, NV);
		//write to global arrays
		for (kk=0;kk<NV;kk++) d_Parameters[Nfits*kk+subregion]=newTheta[kk];
		for (kk=0;kk<NV;kk++) d_CRLBs[Nfits*kk+subregion]=Diag[kk];
		d_LogLikelihood[subregion] = Div;

		return;
}

//*********************************************************************************************************************************************

 void kernel_MLEFit_LM_sigmaxy_EMCCD(const int subregion,const float *d_data, const float PSFSigma, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits){
			/*! 
	 * \brief basic MLE fitting kernel.  No additional parameters are computed.
	 * \param d_data array of subregions to fit copied to GPU
	 * \param PSFSigma sigma of the point spread function
	 * \param sz nxn size of the subregion to fit
	 * \param iterations number of iterations for solution to converge
	 * \param d_Parameters array of fitting parameters to return for each subregion
	 * \param d_CRLBs array of Cramer-Rao lower bound estimates to return for each subregion
	 * \param d_LogLikelihood array of loglikelihood estimates to return for each subregion
	 * \param Nfits number of subregions to fit
	 */

		const int NV=NV_PS2;
		float M[NV*NV],Diag[NV], Minv[NV*NV];
		//int tx = threadIdx.x;
		//int bx = blockIdx.x;
		//int BlockSize = blockDim.x;
		int ii, jj, kk, ll, l, m, i;


		float model, data;
		float Div;

		float newTheta[NV],oldTheta[NV];
		float newLambda = INIT_LAMBDA, oldLambda = INIT_LAMBDA, mu;
		float newUpdate[NV] = {1e13, 1e13, 1e13, 1e13, 1e13, 1e13},oldUpdate[NV] = {1e13, 1e13, 1e13, 1e13, 1e13, 1e13};
		float maxJump[NV]={1.0,1.0,100,20,0.5,0.5};
		float newDudt[NV] ={0};

		float newErr = 1e12, oldErr = 1e13;

		float jacobian[NV]={0};
		float hessian[NV*NV]={0};
		float t1,t2;

		float Nmax;
		//float temp;
		int errFlag=0;
		float L[NV*NV] = {0}, U[NV*NV] = {0};


		//Prevent read/write past end of array
		if (subregion>=Nfits) return;

		for (ii=0;ii<NV*NV;ii++)M[ii]=0;
		for (ii=0;ii<NV*NV;ii++)Minv[ii]=0;

		//copy in data
		const float *s_data = d_data+(sz*sz*subregion);

		//initial values
		kernel_CenterofMass2D(sz, s_data, &newTheta[0], &newTheta[1]);
		kernel_GaussFMaxMin2D(sz, PSFSigma, s_data, &Nmax, &newTheta[3]);
		newTheta[2]=max(0.0, (Nmax-newTheta[3])*2*PI*PSFSigma*PSFSigma);
		newTheta[3] = max(newTheta[3],0.01);
		newTheta[4]=PSFSigma;
		newTheta[5]=PSFSigma;

		maxJump[2]=max(newTheta[2],maxJump[2]);

		maxJump[3]=max(newTheta[3],maxJump[3]);

		for (ii=0;ii<NV;ii++)oldTheta[ii]=newTheta[ii];

		//updateFitValues
		newErr = 0;
		memset(jacobian,0,NV*sizeof(float));
		memset(hessian,0,NV*NV*sizeof(float));
		for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
			kernel_DerivativeGauss2D_sigmaxy( ii,  jj, newTheta, newDudt, &model);
			data=s_data[sz*jj+ii];

			if (data>0)
				newErr = newErr + 2*((model-data)-data*log(model/data));
			else
			{
				newErr = newErr + 2*model;
				data = 0;
			}

			t1 = 1-data/model;
			for (l=0;l<NV;l++){
				jacobian[l]+=t1*newDudt[l];
			}

			t2 = data/pow(model,2);
			for (l=0;l<NV;l++) for(m=l;m<NV;m++) {
				hessian[l*NV+m] +=t2*newDudt[l]*newDudt[m];
				hessian[m*NV+l] = hessian[l*NV+m];
			}
		}

		for (kk=0;kk<iterations;kk++) {//main iterative loop

			if(fabs((newErr-oldErr)/newErr)<TOLERANCE){
				//CONVERGED;
				break;
			}
			else{
				if(newErr>ACCEPTANCE*oldErr){
					//copy Fitdata

					for (i=0;i<NV;i++){
						newTheta[i]=oldTheta[i];
						newUpdate[i]=oldUpdate[i];
					}
					newLambda = oldLambda;
					newErr = oldErr;
					mu = max( (1 + newLambda*SCALE_UP)/(1 + newLambda),1.3f);         
					newLambda = SCALE_UP*newLambda;
				}
				else if(newErr<oldErr&&errFlag==0){
					newLambda = SCALE_DOWN*newLambda;
				    mu = 1+newLambda;
				}


				for (i=0;i<NV;i++){
					hessian[i*NV+i]=hessian[i*NV+i]*mu;
				}
				memset(L,0,NV*sizeof(float));
				memset(U,0,NV*sizeof(float));
				errFlag = kernel_cholesky(hessian,NV,L,U);
				if (errFlag ==0){
					for (i=0;i<NV;i++){
						oldTheta[i]=newTheta[i];
						oldUpdate[i] = newUpdate[i];
					}
					oldLambda = newLambda;
					oldErr=newErr;

					kernel_luEvaluate(L,U,jacobian,NV,newUpdate);	
					
					//updateFitParameters
					for (ll=0;ll<NV;ll++){
						if (newUpdate[ll]/oldUpdate[ll]< -0.5f){
							maxJump[ll] = maxJump[ll]*0.5;
						}
					    newUpdate[ll] = newUpdate[ll]/(1+fabs(newUpdate[ll]/maxJump[ll]));
						newTheta[ll] = newTheta[ll]-newUpdate[ll];
					}
					//restrict range
					newTheta[0] = max(newTheta[0],(float(sz)-1)/2-sz/4.0);
					newTheta[0] = min(newTheta[0],(float(sz)-1)/2+sz/4.0);
					newTheta[1] = max(newTheta[1],(float(sz)-1)/2-sz/4.0);
					newTheta[1] = min(newTheta[1],(float(sz)-1)/2+sz/4.0);
					newTheta[2] = max(newTheta[2],1.0);
					newTheta[3] = max(newTheta[3],0.01);
					newTheta[4] = max(newTheta[4],PSFSigma/10.0f);
					newTheta[5] = max(newTheta[5],PSFSigma/10.0f);

					newErr = 0;
					memset(jacobian,0,NV*sizeof(float));
					memset(hessian,0,NV*NV*sizeof(float));
					for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
						//calculating derivatives
						kernel_DerivativeGauss2D_sigmaxy( ii,  jj, newTheta, newDudt, &model);
						data=s_data[sz*jj+ii];			

						if (data>0)
							newErr = newErr + 2*((model-data)-data*log(model/data));
						else
						{
							newErr = newErr + 2*model;
							data = 0;
						}

						t1 = 1-data/model;
						for (l=0;l<NV;l++){
							jacobian[l]+=t1*newDudt[l];
						}

						t2 = data/pow(model,2);
						for (l=0;l<NV;l++) for(m=l;m<NV;m++) {
							hessian[l*NV+m] +=t2*newDudt[l]*newDudt[m];
							hessian[m*NV+l] = hessian[l*NV+m];
						}
					}
				}
				else
				{
					mu = max( (1 + newLambda*SCALE_UP)/(1 + newLambda),1.3f);         
					newLambda = SCALE_UP*newLambda;
				}
			}
		}
		d_Parameters[Nfits*NV+subregion]=kk;
		// Calculating the CRLB and LogLikelihood
		Div=0.0f;
		for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
			//need to check why don't use newTheta[4] instead of PSFSigma!!!
			kernel_DerivativeGauss2D_sigmaxy( ii,  jj, newTheta, newDudt, &model);
			data=s_data[sz*jj+ii];

			//Building the Fisher Information Matrix
			for (kk=0;kk<NV;kk++)for (ll=kk;ll<NV;ll++){
				M[kk*NV+ll]+= newDudt[ll]*newDudt[kk]/model;
				M[ll*NV+kk]=M[kk*NV+ll];
			}

			//LogLikelyhood
			if (model>0)
				if (data>0)Div+=data*log(model)-model-data*log(data)+data;
				else
					Div+=-model;
		}

		// Matrix inverse (CRLB=F^-1) and output assigments
		kernel_MatInvN(M, Minv, Diag, NV);
		//write to global arrays
		for (kk=0;kk<NV;kk++) d_Parameters[Nfits*kk+subregion]=newTheta[kk];
		for (kk=0;kk<NV;kk++) d_CRLBs[Nfits*kk+subregion]=Diag[kk];
		d_LogLikelihood[subregion] = Div;

		return;
}

//******************************************************************************************************

void kernel_splineMLEFit_z_EMCCD(const int subregion,const float *d_data,const float *d_coeff, const int spline_xsize, const int spline_ysize, const int spline_zsize, const int sz, const int iterations, 
	float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,float initZ, const int Nfits){
		/*! 
	 * \brief basic MLE fitting kernel.  No additional parameters are computed.
	 * \param d_data array of subregions to fit copied to GPU
	 * \param d_coeff array of spline coefficients of the PSF model
	 * \param spline_xsize,spline_ysize,spline_zsize, x,y,z size of spline coefficients
	 * \param sz nxn size of the subregion to fit
	 * \param iterations number of iterations for solution to converge
	 * \param d_Parameters array of fitting parameters to return for each subregion
	 * \param d_CRLBs array of Cramer-Rao lower bound estimates to return for each subregion
	 * \param d_LogLikelihood array of loglikelihood estimates to return for each subregion
	 * \param initZ intial z position used for fitting
	 * \param Nfits number of subregions to fit
	 */
	
	
   const int NV=NV_PSP;
    float M[NV*NV],Diag[NV], Minv[NV*NV];
    //int tx = threadIdx.x;
    //int bx = blockIdx.x;
    //int BlockSize = blockDim.x;
    int ii, jj, kk, ll, l, m, i;
	int xstart, ystart, zstart;

	const float *s_coeff;
	s_coeff = d_coeff;

    float model, data;
    float Div;
    //float dudt[NV_PS];
    float newTheta[NV],oldTheta[NV];
	float newLambda = INIT_LAMBDA, oldLambda = INIT_LAMBDA, mu;
	float newUpdate[NV] = {1e13, 1e13, 1e13, 1e13, 1e13},oldUpdate[NV] = {1e13, 1e13, 1e13, 1e13, 1e13};
	float maxJump[NV]={1.0,1.0,100,20,2};
	float newDudt[NV] ={0};

	float newErr = 1e12, oldErr = 1e13;

	float off;
	float jacobian[NV]={0};
	float hessian[NV*NV]={0};
	float t1,t2;

	float Nmax;
	float xc,yc,zc;
	float delta_f[64]={0}, delta_dxf[64]={0}, delta_dyf[64]={0}, delta_dzf[64]={0};
	int errFlag=0;
	float L[NV*NV] = {0}, U[NV*NV] = {0};

    
    //Prevent read/write past end of array
    if (subregion>=Nfits) return;
    
    for (ii=0;ii<NV*NV;ii++)M[ii]=0;
    for (ii=0;ii<NV*NV;ii++)Minv[ii]=0;

    //copy in data
      const float *s_data = d_data+(sz*sz*subregion);
    
    //initial values
    kernel_CenterofMass2D(sz, s_data, &newTheta[0], &newTheta[1]);
    kernel_GaussFMaxMin2D(sz, 1.5, s_data, &Nmax, &newTheta[3]);
	//central pixel of spline model
	newTheta[3] = max(newTheta[3],0.01);
	newTheta[2]= (Nmax-newTheta[3])/d_coeff[(int)(spline_zsize/2)*(spline_xsize*spline_ysize)+(int)(spline_ysize/2)*spline_xsize+(int)(spline_xsize/2)]*4;

    newTheta[4]=initZ;

	maxJump[2]=max(newTheta[2],maxJump[2]);

	maxJump[3]=max(newTheta[3],maxJump[3]);

	maxJump[4]= max(spline_zsize/3.0f,maxJump[4]);


	for (ii=0;ii<NV;ii++)oldTheta[ii]=newTheta[ii];

	//updateFitValues
	xc = -1.0*((newTheta[0]-float(sz)/2)+0.5);
	yc = -1.0*((newTheta[1]-float(sz)/2)+0.5);

	off = floor((float(spline_xsize)+1.0-float(sz))/2);

	xstart = floor(xc);
	xc = xc-xstart;

	ystart = floor(yc);
	yc = yc-ystart;

	//zstart = floor(newTheta[4]);
	zstart = floor(newTheta[4]);
	zc = newTheta[4] -zstart;

	newErr = 0;
	memset(jacobian,0,NV*sizeof(float));
	memset(hessian,0,NV*NV*sizeof(float));
	kernel_computeDelta3D(xc, yc, zc, delta_f, delta_dxf, delta_dyf, delta_dzf);

	for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
		kernel_DerivativeSpline(ii+xstart+off,jj+ystart+off,zstart,spline_xsize,spline_ysize,spline_zsize,delta_f,delta_dxf,delta_dyf,delta_dzf,s_coeff,newTheta,newDudt,&model);
		data=s_data[sz*jj+ii];

		if (data>0)
			newErr = newErr + 2*((model-data)-data*log(model/data));
		else
		{
			newErr = newErr + 2*model;
			data = 0;
		}

		t1 = 1-data/model;
		for (l=0;l<NV;l++){
			jacobian[l]+=t1*newDudt[l];
		}

		t2 = data/pow(model,2);
		for (l=0;l<NV;l++) for(m=l;m<NV;m++) {
			hessian[l*NV+m] +=t2*newDudt[l]*newDudt[m];
			hessian[m*NV+l] = hessian[l*NV+m];
		}
	}

	for (kk=0;kk<iterations;kk++) {//main iterative loop

			if(fabs((newErr-oldErr)/newErr)<TOLERANCE){
				//CONVERGED;
				break;
			}
			else{
				if(newErr>ACCEPTANCE*oldErr){
					//copy Fitdata
					
					for (i=0;i<NV;i++){
						newTheta[i]=oldTheta[i];
						newUpdate[i]=oldUpdate[i];
					}
					newLambda = oldLambda;
					newErr = oldErr;
					mu = max( (1 + newLambda*SCALE_UP)/(1 + newLambda),1.3f);         
					newLambda = SCALE_UP*newLambda;
				}
				else if(newErr<oldErr&&errFlag==0){
					newLambda = SCALE_DOWN*newLambda;
				    mu = 1+newLambda;
				}
				

				for (i=0;i<NV;i++){
					hessian[i*NV+i]=hessian[i*NV+i]*mu;
				}
				memset(L,0,NV*sizeof(float));
				memset(U,0,NV*sizeof(float));
				errFlag = kernel_cholesky(hessian,NV,L,U);
				if (errFlag ==0){
					for (i=0;i<NV;i++){
						oldTheta[i]=newTheta[i];
						oldUpdate[i] = newUpdate[i];
					}
					oldLambda = newLambda;
					oldErr=newErr;

					kernel_luEvaluate(L,U,jacobian,NV,newUpdate);	
					
					//updateFitParameters
					for (ll=0;ll<NV;ll++){
						if (newUpdate[ll]/oldUpdate[ll]< -0.5f){
							maxJump[ll] = maxJump[ll]*0.5f;
						}
					    newUpdate[ll] = newUpdate[ll]/(1+fabs(newUpdate[ll]/maxJump[ll]));
						newTheta[ll] = newTheta[ll]-newUpdate[ll];
					}
					//restrict range
					newTheta[0] = max(newTheta[0],(float(sz)-1)/2.0f-sz/4.0f);
					newTheta[0] = min(newTheta[0],(float(sz)-1)/2.0f+sz/4.0f);
					newTheta[1] = max(newTheta[1],(float(sz)-1)/2.0f-sz/4.0f);
					newTheta[1] = min(newTheta[1],(float(sz)-1)/2.0f+sz/4.0f);
					newTheta[2] = max(newTheta[2],1.0f);
					newTheta[3] = max(newTheta[3],0.01f);
					newTheta[4] = max(newTheta[4],0.0f);
					newTheta[4] = min(newTheta[4],float(spline_zsize));

					//updateFitValues
					xc = -1.0*((newTheta[0]-float(sz)/2)+0.5f);
					yc = -1.0*((newTheta[1]-float(sz)/2)+0.5f);

					xstart = floor(xc);
					xc = xc-xstart;

					ystart = floor(yc);
					yc = yc-ystart;

					//zstart = floor(newTheta[4]);
					zstart = floor(newTheta[4]);
					zc = newTheta[4] -zstart;


					newErr = 0;
					memset(jacobian,0,NV*sizeof(float));
					memset(hessian,0,NV*NV*sizeof(float));
					kernel_computeDelta3D(xc, yc, zc, delta_f, delta_dxf, delta_dyf, delta_dzf);
					for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
						kernel_DerivativeSpline(ii+xstart+off,jj+ystart+off,zstart,spline_xsize,spline_ysize,spline_zsize,delta_f,delta_dxf,delta_dyf,delta_dzf,s_coeff,newTheta,newDudt,&model);
						data=s_data[sz*jj+ii];

						if (data>0)
							newErr = newErr + 2*((model-data)-data*log(model/data));
						else
						{
							newErr = newErr + 2*model;
							data = 0;
						}

						t1 = 1-data/model;
						for (l=0;l<NV;l++){
							jacobian[l]+=t1*newDudt[l];
						}

						t2 = data/pow(model,2);
						for (l=0;l<NV;l++) for(m=l;m<NV;m++) {
							hessian[l*NV+m] +=t2*newDudt[l]*newDudt[m];
							hessian[m*NV+l] = hessian[l*NV+m];
						}
					}
				}
				else
				{
					mu = max( (1 + newLambda*SCALE_UP)/(1 + newLambda),1.3f);         
					newLambda = SCALE_UP*newLambda;
				}
			}	
	}
	//output iteration time
	d_Parameters[Nfits*NV+subregion]=kk;
    
    // Calculating the CRLB and LogLikelihood
	Div=0.0;

	xc = -1.0*((newTheta[0]-float(sz)/2)+0.5);
	yc = -1.0*((newTheta[1]-float(sz)/2)+0.5);

	xstart = floor(xc);
	xc = xc-xstart;

	ystart = floor(yc);
	yc = yc-ystart;

	//zstart = floor(newTheta[4]);
	zstart = floor(newTheta[4]);
	zc = newTheta[4] -zstart;

	kernel_computeDelta3D(xc, yc, zc, delta_f, delta_dxf, delta_dyf, delta_dzf);

    for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
		kernel_DerivativeSpline(ii+xstart+off,jj+ystart+off,zstart,spline_xsize,spline_ysize,spline_zsize,delta_f,delta_dxf,delta_dyf,delta_dzf,s_coeff,newTheta,newDudt,&model);
		data=s_data[sz*jj+ii];
        
        //Building the Fisher Information Matrix
        for (kk=0;kk<NV;kk++)for (ll=kk;ll<NV;ll++){
            M[kk*NV+ll]+= newDudt[ll]*newDudt[kk]/model;
            M[ll*NV+kk]=M[kk*NV+ll];
        }
        
        //LogLikelyhood
        if (model>0)
            if (data>0)Div+=data*log(model)-model-data*log(data)+data;
            else
                Div+=-model;
    }
    
    // Matrix inverse (CRLB=F^-1) and output assigments
    kernel_MatInvN(M, Minv, Diag, NV);
  
    
    //write to global arrays
    for (kk=0;kk<NV;kk++) d_Parameters[Nfits*kk+subregion]=newTheta[kk];
   for (kk=0;kk<NV;kk++) d_CRLBs[Nfits*kk+subregion]=Diag[kk];
   d_LogLikelihood[subregion] = Div;
    return;
}


//*********************************************************************************************************************************************
 void kernel_MLEFit_LM_sCMOS(const int subregion, const float *d_data,const float PSFSigma, const int sz, const int iterations, 
	float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits,const float *d_varim){
		/*! 
	 * \brief basic MLE fitting kernel.  No additional parameters are computed.
	 * \param d_data array of subregions to fit copied to GPU
	 * \param PSFSigma sigma of the point spread function
	 * \param sz nxn size of the subregion to fit
	 * \param iterations number of iterations for solution to converge
	 * \param d_Parameters array of fitting parameters to return for each subregion
	 * \param d_CRLBs array of Cramer-Rao lower bound estimates to return for each subregion
	 * \param d_LogLikelihood array of loglikelihood estimates to return for each subregion
	 * \param Nfits number of subregions to fit
	 * \d_varim variance map of scmos
	 */

		const int NV=NV_P;
		float M[NV*NV],Diag[NV], Minv[NV*NV];
		//int tx = threadIdx.x;
		//int bx = blockIdx.x;
		//int BlockSize = blockDim.x;
		int ii, jj, kk, ll, l, m, i;


		float model, data;
		float Div;

		float newTheta[NV],oldTheta[NV];
		float newLambda = INIT_LAMBDA, oldLambda = INIT_LAMBDA, mu;
	    float newUpdate[NV] = {1e13, 1e13, 1e13, 1e13},oldUpdate[NV] = {1e13, 1e13, 1e13, 1e13};
		float maxJump[NV]={1.0,1.0,100,20};
		float newDudt[NV] ={0};

		float newErr = 1e12, oldErr = 1e13;

		float jacobian[NV]={0};
		float hessian[NV*NV]={0};
		float t1,t2;

		float Nmax;
		int errFlag=0;
		float L[NV*NV] = {0}, U[NV*NV] = {0};


		//Prevent read/write past end of array
		if (subregion>=Nfits) return;

		for (ii=0;ii<NV*NV;ii++)M[ii]=0;
		for (ii=0;ii<NV*NV;ii++)Minv[ii]=0;

		//copy in data
		const float *s_data = d_data+(sz*sz*subregion);
		const float *s_varim = d_varim+(sz*sz*subregion);

		//initial values
		kernel_CenterofMass2D(sz, s_data, &newTheta[0], &newTheta[1]);
		kernel_GaussFMaxMin2D(sz, PSFSigma, s_data, &Nmax, &newTheta[3]);
		newTheta[2]=max(0.0, (Nmax-newTheta[3])*2*PI*PSFSigma*PSFSigma);
		newTheta[3] = max(newTheta[3],0.01);

		maxJump[2]=max(newTheta[2],maxJump[2]);

		maxJump[3]=max(newTheta[3],maxJump[3]);

		for (ii=0;ii<NV;ii++)oldTheta[ii]=newTheta[ii];

		//updateFitValues
		newErr = 0;
		memset(jacobian,0,NV*sizeof(float));
		memset(hessian,0,NV*NV*sizeof(float));
		for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
			kernel_DerivativeGauss2D(ii,jj,PSFSigma,newTheta,newDudt,&model);
			model +=s_varim[sz*jj+ii];
			data=s_data[sz*jj+ii]+s_varim[sz*jj+ii];

			if (data>0)
				newErr = newErr + 2*((model-data)-data*log(model/data));
			else
			{
				newErr = newErr + 2*model;
				data = 0;
			}

			t1 = 1-data/model;
			for (l=0;l<NV;l++){
				jacobian[l]+=t1*newDudt[l];
			}

			t2 = data/pow(model,2);
			for (l=0;l<NV;l++) for(m=l;m<NV;m++) {
				hessian[l*NV+m] +=t2*newDudt[l]*newDudt[m];
				hessian[m*NV+l] = hessian[l*NV+m];
			}
		}

		for (kk=0;kk<iterations;kk++) {//main iterative loop

			if(fabs((newErr-oldErr)/newErr)<TOLERANCE){
				//CONVERGED;
				break;
			}
			else{
				if(newErr>ACCEPTANCE*oldErr){
					//copy Fitdata

					for (i=0;i<NV;i++){
						newTheta[i]=oldTheta[i];
						newUpdate[i]=oldUpdate[i];
					}
					newLambda = oldLambda;
					newErr = oldErr;
					mu = max( (1 + newLambda*SCALE_UP)/(1 + newLambda),1.3f);         
					newLambda = SCALE_UP*newLambda;

				}
				else if(newErr<oldErr&&errFlag==0){
					newLambda = SCALE_DOWN*newLambda;
				    mu = 1+newLambda;
				}


				for (i=0;i<NV;i++){
					hessian[i*NV+i]=hessian[i*NV+i]*newLambda;
				}
				memset(L,0,NV*sizeof(float));
				memset(U,0,NV*sizeof(float));
				errFlag = kernel_cholesky(hessian,NV,L,U);
				if (errFlag ==0){
					for (i=0;i<NV;i++){
						oldTheta[i]=newTheta[i];
						oldUpdate[i] = newUpdate[i];
					}
					oldLambda = newLambda;
					oldErr=newErr;

					kernel_luEvaluate(L,U,jacobian,NV,newUpdate);	
					
					//updateFitParameters
					for (ll=0;ll<NV;ll++){
						if (newUpdate[ll]/oldUpdate[ll]< -0.5f){
							maxJump[ll] = maxJump[ll]*0.5;
						}
					    newUpdate[ll] = newUpdate[ll]/(1+fabs(newUpdate[ll]/maxJump[ll]));
						newTheta[ll] = newTheta[ll]-newUpdate[ll];
					}
					//restrict range
					newTheta[0] = max(newTheta[0],(float(sz)-1)/2-sz/4.0);
					newTheta[0] = min(newTheta[0],(float(sz)-1)/2+sz/4.0);
					newTheta[1] = max(newTheta[1],(float(sz)-1)/2-sz/4.0);
					newTheta[1] = min(newTheta[1],(float(sz)-1)/2+sz/4.0);
					newTheta[2] = max(newTheta[2],1.0);
					newTheta[3] = max(newTheta[3],0.01);
					


					newErr = 0;
					memset(jacobian,0,NV*sizeof(float));
					memset(hessian,0,NV*NV*sizeof(float));
					for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
						//calculating derivatives
						kernel_DerivativeGauss2D(ii,jj,PSFSigma,newTheta,newDudt,&model);
						model +=s_varim[sz*jj+ii];
						data=s_data[sz*jj+ii]+s_varim[sz*jj+ii];		

						if (data>0)
							newErr = newErr + 2*((model-data)-data*log(model/data));
						else
						{
							newErr = newErr + 2*model;
							data = 0;
						}

						t1 = 1-data/model;
						for (l=0;l<NV;l++){
							jacobian[l]+=t1*newDudt[l];
						}

						t2 = data/pow(model,2);
						for (l=0;l<NV;l++) for(m=l;m<NV;m++) {
							hessian[l*NV+m] +=t2*newDudt[l]*newDudt[m];
							hessian[m*NV+l] = hessian[l*NV+m];
						}
					}
				}
				else
				{
					mu = max( (1 + newLambda*SCALE_UP)/(1 + newLambda),1.3f);         
					newLambda = SCALE_UP*newLambda;
				}
			}
		}
		//output iteration
		d_Parameters[Nfits*NV+subregion]=kk;
		// Calculating the CRLB and LogLikelihood
		Div=0.0f;
		for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
			kernel_DerivativeGauss2D(ii,jj,PSFSigma,newTheta,newDudt,&model);
			model +=s_varim[sz*jj+ii];
			data=s_data[sz*jj+ii]+s_varim[sz*jj+ii];

			//Building the Fisher Information Matrix
			for (kk=0;kk<NV;kk++)for (ll=kk;ll<NV;ll++){
				M[kk*NV+ll]+= newDudt[ll]*newDudt[kk]/model;
				M[ll*NV+kk]=M[kk*NV+ll];
			}

			//LogLikelyhood
			if (model>0)
				if (data>0)Div+=data*log(model)-model-data*log(data)+data;
				else
					Div+=-model;
		}

		// Matrix inverse (CRLB=F^-1) and output assigments
		kernel_MatInvN(M, Minv, Diag, NV);
		//write to global arrays
		for (kk=0;kk<NV;kk++) d_Parameters[Nfits*kk+subregion]=newTheta[kk];
		for (kk=0;kk<NV;kk++) d_CRLBs[Nfits*kk+subregion]=Diag[kk];
		d_LogLikelihood[subregion] = Div;

		return;
}

//*********************************************************************************************************************************************

 void kernel_MLEFit_LM_Sigma_sCMOS(const int subregion, const float *d_data,const float PSFSigma, const int sz, const int iterations, 
	float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits,const float *d_varim){
		/*! 
	 * \brief basic MLE fitting kernel.  No additional parameters are computed.
	 * \param d_data array of subregions to fit copied to GPU
	 * \param PSFSigma sigma of the point spread function
	 * \param sz nxn size of the subregion to fit
	 * \param iterations number of iterations for solution to converge
	 * \param d_Parameters array of fitting parameters to return for each subregion
	 * \param d_CRLBs array of Cramer-Rao lower bound estimates to return for each subregion
	 * \param d_LogLikelihood array of loglikelihood estimates to return for each subregion
	 * \param Nfits number of subregions to fit
	 * \d_varim variance map of scmos
	 */
		const int NV=NV_PS;
		float M[NV*NV],Diag[NV], Minv[NV*NV];
		/*int tx = threadIdx.x;
		int bx = blockIdx.x;
		int BlockSize = blockDim.x;*/
		int ii, jj, kk, ll, l, m, i;


		float model, data;
		float Div;

		float newTheta[NV],oldTheta[NV];
		float newLambda = INIT_LAMBDA, oldLambda = INIT_LAMBDA, mu;
		float newUpdate[NV] = {1e13, 1e13, 1e13, 1e13, 1e13},oldUpdate[NV] = {1e13, 1e13, 1e13, 1e13, 1e13};
		float maxJump[NV]={1.0,1.0,100,20,0.5};
		float newDudt[NV] ={0};

		float newErr = 1e12, oldErr = 1e13;

		float jacobian[NV]={0};
		float hessian[NV*NV]={0};
		float t1,t2;

		float Nmax;
		int errFlag=0;
		float L[NV*NV] = {0}, U[NV*NV] = {0};


		//Prevent read/write past end of array
		if (subregion>=Nfits) return;

		for (ii=0;ii<NV*NV;ii++)M[ii]=0;
		for (ii=0;ii<NV*NV;ii++)Minv[ii]=0;

		//copy in data
		const float *s_data = d_data+(sz*sz*subregion);
		const float *s_varim = d_varim+(sz*sz*subregion);

		//initial values
		kernel_CenterofMass2D(sz, s_data, &newTheta[0], &newTheta[1]);
		kernel_GaussFMaxMin2D(sz, PSFSigma, s_data, &Nmax, &newTheta[3]);
		newTheta[2]=max(0.0, (Nmax-newTheta[3])*2*PI*PSFSigma*PSFSigma);
		newTheta[3] = max(newTheta[3],0.01);
		newTheta[4]=PSFSigma;

		maxJump[2]=max(newTheta[2],maxJump[2]);

		maxJump[3]=max(newTheta[3],maxJump[3]);


		for (ii=0;ii<NV;ii++)oldTheta[ii]=newTheta[ii];

		//updateFitValues
		newErr = 0;
		memset(jacobian,0,NV*sizeof(float));
		memset(hessian,0,NV*NV*sizeof(float));
		for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
			kernel_DerivativeGauss2D_sigma(ii,jj,newTheta,newDudt,&model);
			model +=s_varim[sz*jj+ii];
			data=s_data[sz*jj+ii]+s_varim[sz*jj+ii];

			if (data>0)
				newErr = newErr + 2*((model-data)-data*log(model/data));
			else
			{
				newErr = newErr + 2*model;
				data = 0;
			}

			t1 = 1-data/model;
			for (l=0;l<NV;l++){
				jacobian[l]+=t1*newDudt[l];
			}

			t2 = data/pow(model,2);
			for (l=0;l<NV;l++) for(m=l;m<NV;m++) {
				hessian[l*NV+m] +=t2*newDudt[l]*newDudt[m];
				hessian[m*NV+l] = hessian[l*NV+m];
			}
		}

		for (kk=0;kk<iterations;kk++) {//main iterative loop

			if(fabs((newErr-oldErr)/newErr)<TOLERANCE){
				//CONVERGED;
				break;
			}
			else{
				if(newErr>ACCEPTANCE*oldErr){
					//copy Fitdata

					for (i=0;i<NV;i++){
						newTheta[i]=oldTheta[i];
						newUpdate[i]=oldUpdate[i];
					}
					newLambda = oldLambda;
					newErr = oldErr;
					mu = max( (1 + newLambda*SCALE_UP)/(1 + newLambda),1.3f);         
					newLambda = SCALE_UP*newLambda;

				}
				else if(newErr<oldErr&&errFlag==0){
					newLambda = SCALE_DOWN*newLambda;
				    mu = 1+newLambda;
				}


				for (i=0;i<NV;i++){
					hessian[i*NV+i]=hessian[i*NV+i]*mu;
				}
				memset(L,0,NV*sizeof(float));
				memset(U,0,NV*sizeof(float));
				errFlag = kernel_cholesky(hessian,NV,L,U);
				if (errFlag ==0){
					for (i=0;i<NV;i++){
						oldTheta[i]=newTheta[i];
						oldUpdate[i] = newUpdate[i];
					}
					oldLambda = newLambda;
					oldErr=newErr;

					kernel_luEvaluate(L,U,jacobian,NV,newUpdate);	
					
					//updateFitParameters
					for (ll=0;ll<NV;ll++){
						if (newUpdate[ll]/oldUpdate[ll]< -0.5f){
							maxJump[ll] = maxJump[ll]*0.5;
						}
					    newUpdate[ll] = newUpdate[ll]/(1+fabs(newUpdate[ll]/maxJump[ll]));
						newTheta[ll] = newTheta[ll]-newUpdate[ll];
					}
					//restrict range
					newTheta[0] = max(newTheta[0],(float(sz)-1)/2-sz/4.0);
					newTheta[0] = min(newTheta[0],(float(sz)-1)/2+sz/4.0);
					newTheta[1] = max(newTheta[1],(float(sz)-1)/2-sz/4.0);
					newTheta[1] = min(newTheta[1],(float(sz)-1)/2+sz/4.0);
					newTheta[2] = max(newTheta[2],1.0);
					newTheta[3] = max(newTheta[3],0.01);
					newTheta[4] = max(newTheta[4],0.0);
					newTheta[4] = min(newTheta[4],sz/2.0f);


					newErr = 0;
					memset(jacobian,0,NV*sizeof(float));
					memset(hessian,0,NV*NV*sizeof(float));
					for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
						//calculating derivatives
						kernel_DerivativeGauss2D_sigma(ii,jj,newTheta,newDudt,&model);
						model +=s_varim[sz*jj+ii];
						data=s_data[sz*jj+ii]+s_varim[sz*jj+ii];		

						if (data>0)
							newErr = newErr + 2*((model-data)-data*log(model/data));
						else
						{
							newErr = newErr + 2*model;
							data = 0;
						}

						t1 = 1-data/model;
						for (l=0;l<NV;l++){
							jacobian[l]+=t1*newDudt[l];
						}

						t2 = data/pow(model,2);
						for (l=0;l<NV;l++) for(m=l;m<NV;m++) {
							hessian[l*NV+m] +=t2*newDudt[l]*newDudt[m];
							hessian[m*NV+l] = hessian[l*NV+m];
						}
					}
				}
				else
				{
					mu = max( (1 + newLambda*SCALE_UP)/(1 + newLambda),1.3f);         
					newLambda = SCALE_UP*newLambda;
				}
			}
		}
		d_Parameters[Nfits*NV+subregion]=kk;
		// Calculating the CRLB and LogLikelihood
		Div=0.0f;
		for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
			kernel_DerivativeGauss2D_sigma(ii,jj,newTheta,newDudt,&model);
			model +=s_varim[sz*jj+ii];
			data=s_data[sz*jj+ii]+s_varim[sz*jj+ii];

			//Building the Fisher Information Matrix
			for (kk=0;kk<NV;kk++)for (ll=kk;ll<NV;ll++){
				M[kk*NV+ll]+= newDudt[ll]*newDudt[kk]/model;
				M[ll*NV+kk]=M[kk*NV+ll];
			}

			//LogLikelyhood
			if (model>0)
				if (data>0)Div+=data*log(model)-model-data*log(data)+data;
				else
					Div+=-model;
		}

		// Matrix inverse (CRLB=F^-1) and output assigments
		kernel_MatInvN(M, Minv, Diag, NV);
		//write to global arrays
		for (kk=0;kk<NV;kk++) d_Parameters[Nfits*kk+subregion]=newTheta[kk];
		for (kk=0;kk<NV;kk++) d_CRLBs[Nfits*kk+subregion]=Diag[kk];
		d_LogLikelihood[subregion] = Div;

		return;
}

//*********************************************************************************************************************************************

 void kernel_MLEFit_LM_z_sCMOS(const int subregion,const float *d_data, const float PSFSigma_x, const float Ax, const float Ay, const float Bx, 
	const float By, const float gamma, const float d, const float PSFSigma_y, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits,const float *d_varim){
			/*! 
	 * \brief basic MLE fitting kernel.  No additional parameters are computed.
	 * \param d_data array of subregions to fit copied to GPU
	 * \param PSFSigma_x sigma of the point spread function on the x axis
	 * \param Ax ???
	 * \param Ay ???
	 * \param Bx ???
	 * \param By ???
	 * \param gamma ???
	 * \param d ???
	 * \param PSFSigma_y sigma of the point spread function on the y axis
	 * \param sz nxn size of the subregion to fit
	 * \param iterations number of iterations for solution to converge
	 * \param d_Parameters array of fitting parameters to return for each subregion
	 * \param d_CRLBs array of Cramer-Rao lower bound estimates to return for each subregion
	 * \param d_LogLikelihood array of loglikelihood estimates to return for each subregion
	 * \param Nfits number of subregions to fit
	 * \d_varim variance map of scmos
	 */

		const int NV=NV_PZ;
		float M[NV*NV],Diag[NV], Minv[NV*NV];
		//int tx = threadIdx.x;
		//int bx = blockIdx.x;
		//int BlockSize = blockDim.x;
		int ii, jj, kk, ll, l, m, i;


		float model, data;
		float Div;
		float PSFy, PSFx;

		float newTheta[NV],oldTheta[NV];
		float newLambda = INIT_LAMBDA, oldLambda = INIT_LAMBDA, mu;
		float newUpdate[NV] = {1e13, 1e13, 1e13, 1e13, 1e13},oldUpdate[NV] = {1e13, 1e13, 1e13, 1e13, 1e13};
		float maxJump[NV]={1.0,1.0,100,20,2};
		float newDudt[NV] ={0};

		float newErr = 1e12, oldErr = 1e13;

		float jacobian[NV]={0};
		float hessian[NV*NV]={0};
		float t1,t2;

		float Nmax;
		int errFlag=0;
		float L[NV*NV] = {0}, U[NV*NV] = {0};


		//Prevent read/write past end of array
		if (subregion>=Nfits) return;

		for (ii=0;ii<NV*NV;ii++)M[ii]=0;
		for (ii=0;ii<NV*NV;ii++)Minv[ii]=0;

		//copy in data
		const float *s_data = d_data+(sz*sz*subregion);
		const float *s_varim = d_varim+(sz*sz*subregion);
		//initial values
		kernel_CenterofMass2D(sz, s_data, &newTheta[0], &newTheta[1]);
		kernel_GaussFMaxMin2D(sz, PSFSigma_x, s_data, &Nmax, &newTheta[3]);
		newTheta[2]=max(0.0, (Nmax-newTheta[3])*2*PI*PSFSigma_x*PSFSigma_y*sqrt(2.0f));
		newTheta[3] = max(newTheta[3],0.01);
		newTheta[4]=0;

		maxJump[2]=max(newTheta[2],maxJump[2]);

		maxJump[3]=max(newTheta[3],maxJump[3]);

		for (ii=0;ii<NV;ii++)oldTheta[ii]=newTheta[ii];

		//updateFitValues
		newErr = 0;
		memset(jacobian,0,NV*sizeof(float));
		memset(hessian,0,NV*NV*sizeof(float));
		for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
			 kernel_DerivativeIntGauss2Dz(ii, jj, newTheta, PSFSigma_x,PSFSigma_y, Ax,Ay,Bx,By, gamma, d, &PSFx, &PSFy, newDudt, NULL,&model);
			model +=s_varim[sz*jj+ii];
			data=s_data[sz*jj+ii]+s_varim[sz*jj+ii];

			if (data>0)
				newErr = newErr + 2*((model-data)-data*log(model/data));
			else
			{
				newErr = newErr + 2*model;
				data = 0;
			}

			t1 = 1-data/model;
			for (l=0;l<NV;l++){
				jacobian[l]+=t1*newDudt[l];
			}

			t2 = data/pow(model,2);
			for (l=0;l<NV;l++) for(m=l;m<NV;m++) {
				hessian[l*NV+m] +=t2*newDudt[l]*newDudt[m];
				hessian[m*NV+l] = hessian[l*NV+m];
			}
		}

		for (kk=0;kk<iterations;kk++) {//main iterative loop

			if(fabs((newErr-oldErr)/newErr)<TOLERANCE){
				//CONVERGED;
				break;
			}
			else{
				if(newErr>ACCEPTANCE*oldErr){
					//copy Fitdata

					for (i=0;i<NV;i++){
						newTheta[i]=oldTheta[i];
						newUpdate[i]=oldUpdate[i];
					}
					newLambda = oldLambda;
					newErr = oldErr;
					mu = max( (1 + newLambda*SCALE_UP)/(1 + newLambda),1.3f);         
					newLambda = SCALE_UP*newLambda;
				}
				else if(newErr<oldErr&&errFlag==0){
					newLambda = SCALE_DOWN*newLambda;
				    mu = 1+newLambda;
				}

				for (i=0;i<NV;i++){
					hessian[i*NV+i]=hessian[i*NV+i]*mu;
				}
				memset(L,0,NV*sizeof(float));
				memset(U,0,NV*sizeof(float));
				errFlag = kernel_cholesky(hessian,NV,L,U);
				if (errFlag ==0){
					for (i=0;i<NV;i++){
						oldTheta[i]=newTheta[i];
						oldUpdate[i] = newUpdate[i];
					}
					oldLambda = newLambda;
					oldErr=newErr;

					kernel_luEvaluate(L,U,jacobian,NV,newUpdate);	
					
					//updateFitParameters
					for (ll=0;ll<NV;ll++){
						if (newUpdate[ll]/oldUpdate[ll]< -0.5f){
							maxJump[ll] = maxJump[ll]*0.5;
						}
					    newUpdate[ll] = newUpdate[ll]/(1+fabs(newUpdate[ll]/maxJump[ll]));
						newTheta[ll] = newTheta[ll]-newUpdate[ll];
					}
					//restrict range
					newTheta[0] = max(newTheta[0],(float(sz)-1)/2-sz/4.0);
					newTheta[0] = min(newTheta[0],(float(sz)-1)/2+sz/4.0);
					newTheta[1] = max(newTheta[1],(float(sz)-1)/2-sz/4.0);
					newTheta[1] = min(newTheta[1],(float(sz)-1)/2+sz/4.0);
					newTheta[2] = max(newTheta[2],1.0);
					newTheta[3] = max(newTheta[3],0.01);
					//newTheta[4] = max(newTheta[4],0.0);
					//newTheta[4] = min(newTheta[4],sz/2.0f);


					newErr = 0;
					memset(jacobian,0,NV*sizeof(float));
					memset(hessian,0,NV*NV*sizeof(float));
					for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
						//calculating derivatives
						kernel_DerivativeIntGauss2Dz(ii, jj, newTheta, PSFSigma_x,PSFSigma_y, Ax,Ay,Bx,By, gamma, d, &PSFx, &PSFy, newDudt, NULL,&model);
						model +=s_varim[sz*jj+ii];
						data=s_data[sz*jj+ii]+s_varim[sz*jj+ii];		

						if (data>0)
							newErr = newErr + 2*((model-data)-data*log(model/data));
						else
						{
							newErr = newErr + 2*model;
							data = 0;
						}

						t1 = 1-data/model;
						for (l=0;l<NV;l++){
							jacobian[l]+=t1*newDudt[l];
						}

						t2 = data/pow(model,2);
						for (l=0;l<NV;l++) for(m=l;m<NV;m++) {
							hessian[l*NV+m] +=t2*newDudt[l]*newDudt[m];
							hessian[m*NV+l] = hessian[l*NV+m];
						}
					}
				}
				else
				{
					mu = max( (1 + newLambda*SCALE_UP)/(1 + newLambda),1.3f);         
					newLambda = SCALE_UP*newLambda;
				}
			}
		}
		d_Parameters[Nfits*NV+subregion]=kk;
		// Calculating the CRLB and LogLikelihood
		Div=0.0f;
		for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
		    kernel_DerivativeIntGauss2Dz(ii, jj, newTheta, PSFSigma_x,PSFSigma_y, Ax,Ay,Bx,By, gamma, d, &PSFx, &PSFy, newDudt, NULL,&model);
			model +=s_varim[sz*jj+ii];
			data=s_data[sz*jj+ii]+s_varim[sz*jj+ii];

			//Building the Fisher Information Matrix
			for (kk=0;kk<NV;kk++)for (ll=kk;ll<NV;ll++){
				M[kk*NV+ll]+= newDudt[ll]*newDudt[kk]/model;
				M[ll*NV+kk]=M[kk*NV+ll];
			}

			//LogLikelyhood
			if (model>0)
				if (data>0)Div+=data*log(model)-model-data*log(data)+data;
				else
					Div+=-model;
		}

		// Matrix inverse (CRLB=F^-1) and output assigments
		kernel_MatInvN(M, Minv, Diag, NV);
		//write to global arrays
		for (kk=0;kk<NV;kk++) d_Parameters[Nfits*kk+subregion]=newTheta[kk];
		for (kk=0;kk<NV;kk++) d_CRLBs[Nfits*kk+subregion]=Diag[kk];
		d_LogLikelihood[subregion] = Div;

		return;
}

//*********************************************************************************************************************************************

 void kernel_MLEFit_LM_sigmaxy_sCMOS(const int subregion,const float *d_data, const float PSFSigma, const int sz, const int iterations, 
        float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits,const float *d_varim){
			/*! 
	 * \brief basic MLE fitting kernel.  No additional parameters are computed.
	 * \param d_data array of subregions to fit copied to GPU
	 * \param PSFSigma sigma of the point spread function
	 * \param sz nxn size of the subregion to fit
	 * \param iterations number of iterations for solution to converge
	 * \param d_Parameters array of fitting parameters to return for each subregion
	 * \param d_CRLBs array of Cramer-Rao lower bound estimates to return for each subregion
	 * \param d_LogLikelihood array of loglikelihood estimates to return for each subregion
	 * \param Nfits number of subregions to fit
	  * \d_varim variance map of scmos
	 */

		const int NV=NV_PS2;
		float M[NV*NV],Diag[NV], Minv[NV*NV];
		//int tx = threadIdx.x;
		//int bx = blockIdx.x;
		//int BlockSize = blockDim.x;
		int ii, jj, kk, ll, l, m, i;


		float model, data;
		float Div;

		float newTheta[NV],oldTheta[NV];
		float newLambda = INIT_LAMBDA, oldLambda = INIT_LAMBDA, mu;
		float newUpdate[NV] = {1e13, 1e13, 1e13, 1e13, 1e13, 1e13},oldUpdate[NV] = {1e13, 1e13, 1e13, 1e13, 1e13, 1e13};
		float maxJump[NV]={1.0,1.0,100,20,0.5,0.5};
		float newDudt[NV] ={0};

		float newErr = 1e12, oldErr = 1e13;

		float jacobian[NV]={0};
		float hessian[NV*NV]={0};
		float t1,t2;

		float Nmax;
		int errFlag=0;
		float L[NV*NV] = {0}, U[NV*NV] = {0};


		//Prevent read/write past end of array
		if (subregion>=Nfits) return;

		for (ii=0;ii<NV*NV;ii++)M[ii]=0;
		for (ii=0;ii<NV*NV;ii++)Minv[ii]=0;

		//copy in data
		const float *s_data = d_data+(sz*sz*subregion);
		const float *s_varim = d_varim+(sz*sz*subregion);

		//initial values
		kernel_CenterofMass2D(sz, s_data, &newTheta[0], &newTheta[1]);
		kernel_GaussFMaxMin2D(sz, PSFSigma, s_data, &Nmax, &newTheta[3]);
		newTheta[2]=max(0.0, (Nmax-newTheta[3])*2*PI*PSFSigma*PSFSigma);
		newTheta[3] = max(newTheta[3],0.01);
		newTheta[4]=PSFSigma;
		newTheta[5]=PSFSigma;

		maxJump[2]=max(newTheta[2],maxJump[2]);

		maxJump[3]=max(newTheta[3],maxJump[3]);

		for (ii=0;ii<NV;ii++)oldTheta[ii]=newTheta[ii];

		//updateFitValues
		newErr = 0;
		memset(jacobian,0,NV*sizeof(float));
		memset(hessian,0,NV*NV*sizeof(float));
		for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
			kernel_DerivativeGauss2D_sigmaxy( ii,  jj, newTheta, newDudt, &model);
			model +=s_varim[sz*jj+ii];
			data=s_data[sz*jj+ii]+s_varim[sz*jj+ii];

			if (data>0)
				newErr = newErr + 2*((model-data)-data*log(model/data));
			else
			{
				newErr = newErr + 2*model;
				data = 0;
			}

			t1 = 1-data/model;
			for (l=0;l<NV;l++){
				jacobian[l]+=t1*newDudt[l];
			}

			t2 = data/pow(model,2);
			for (l=0;l<NV;l++) for(m=l;m<NV;m++) {
				hessian[l*NV+m] +=t2*newDudt[l]*newDudt[m];
				hessian[m*NV+l] = hessian[l*NV+m];
			}
		}

		for (kk=0;kk<iterations;kk++) {//main iterative loop

			if(fabs((newErr-oldErr)/newErr)<TOLERANCE){
				//CONVERGED;
				break;
			}
			else{
				if(newErr>ACCEPTANCE*oldErr){
					//copy Fitdata

					for (i=0;i<NV;i++){
						newTheta[i]=oldTheta[i];
						newUpdate[i]=oldUpdate[i];
					}
					newLambda = oldLambda;
					newErr = oldErr;
					mu = max( (1 + newLambda*SCALE_UP)/(1 + newLambda),1.3f);         
					newLambda = SCALE_UP*newLambda;
				}
				else if(newErr<oldErr&&errFlag==0){
					newLambda = SCALE_DOWN*newLambda;
				    mu = 1+newLambda;
				}


				for (i=0;i<NV;i++){
					hessian[i*NV+i]=hessian[i*NV+i]*mu;
				}
				memset(L,0,NV*sizeof(float));
				memset(U,0,NV*sizeof(float));
				errFlag = kernel_cholesky(hessian,NV,L,U);
				if (errFlag ==0){
					for (i=0;i<NV;i++){
						oldTheta[i]=newTheta[i];
						oldUpdate[i] = newUpdate[i];
					}
					oldLambda = newLambda;
					oldErr=newErr;

					kernel_luEvaluate(L,U,jacobian,NV,newUpdate);	
					
					//updateFitParameters
					for (ll=0;ll<NV;ll++){
						if (newUpdate[ll]/oldUpdate[ll]< -0.5f){
							maxJump[ll] = maxJump[ll]*0.5;
						}
					    newUpdate[ll] = newUpdate[ll]/(1+fabs(newUpdate[ll]/maxJump[ll]));
						newTheta[ll] = newTheta[ll]-newUpdate[ll];
					}
					//restrict range
					newTheta[0] = max(newTheta[0],(float(sz)-1)/2-sz/4.0);
					newTheta[0] = min(newTheta[0],(float(sz)-1)/2+sz/4.0);
					newTheta[1] = max(newTheta[1],(float(sz)-1)/2-sz/4.0);
					newTheta[1] = min(newTheta[1],(float(sz)-1)/2+sz/4.0);
					newTheta[2] = max(newTheta[2],1.0);
					newTheta[3] = max(newTheta[3],0.01);
					newTheta[4] = max(newTheta[4],PSFSigma/10.0f);
					newTheta[5] = max(newTheta[5],PSFSigma/10.0f);

					newErr = 0;
					memset(jacobian,0,NV*sizeof(float));
					memset(hessian,0,NV*NV*sizeof(float));
					for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
						//calculating derivatives
						kernel_DerivativeGauss2D_sigmaxy( ii,  jj, newTheta, newDudt, &model);
						model +=s_varim[sz*jj+ii];
						data=s_data[sz*jj+ii]+s_varim[sz*jj+ii];			

						if (data>0)
							newErr = newErr + 2*((model-data)-data*log(model/data));
						else
						{
							newErr = newErr + 2*model;
							data = 0;
						}

						t1 = 1-data/model;
						for (l=0;l<NV;l++){
							jacobian[l]+=t1*newDudt[l];
						}

						t2 = data/pow(model,2);
						for (l=0;l<NV;l++) for(m=l;m<NV;m++) {
							hessian[l*NV+m] +=t2*newDudt[l]*newDudt[m];
							hessian[m*NV+l] = hessian[l*NV+m];
						}
					}
				}
				else
				{
					mu = max( (1 + newLambda*SCALE_UP)/(1 + newLambda),1.3f);         
					newLambda = SCALE_UP*newLambda;
				}
			}
		}
		//output iteration
		d_Parameters[Nfits*NV+subregion]=kk;
		// Calculating the CRLB and LogLikelihood
		Div=0.0f;
		for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
			kernel_DerivativeGauss2D_sigmaxy( ii,  jj, newTheta, newDudt, &model);
			model +=s_varim[sz*jj+ii];
			data=s_data[sz*jj+ii]+s_varim[sz*jj+ii];

			//Building the Fisher Information Matrix
			for (kk=0;kk<NV;kk++)for (ll=kk;ll<NV;ll++){
				M[kk*NV+ll]+= newDudt[ll]*newDudt[kk]/model;
				M[ll*NV+kk]=M[kk*NV+ll];
			}

			//LogLikelyhood
			if (model>0)
				if (data>0)Div+=data*log(model)-model-data*log(data)+data;
				else
					Div+=-model;
		}

		// Matrix inverse (CRLB=F^-1) and output assigments
		kernel_MatInvN(M, Minv, Diag, NV);
		//write to global arrays
		for (kk=0;kk<NV;kk++) d_Parameters[Nfits*kk+subregion]=newTheta[kk];
		for (kk=0;kk<NV;kk++) d_CRLBs[Nfits*kk+subregion]=Diag[kk];
		d_LogLikelihood[subregion] = Div;

		return;
}

//******************************************************************************************************

void kernel_splineMLEFit_z_sCMOS(const int subregion,const float *d_data,const float *d_coeff, const int spline_xsize, const int spline_ysize, const int spline_zsize, const int sz, const int iterations, 
	float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,float initZ, const int Nfits,const float *d_varim){
		/*! 
	 * \brief basic MLE fitting kernel.  No additional parameters are computed.
	 * \param d_data array of subregions to fit copied to GPU
	 * \param d_coeff array of spline coefficients of the PSF model
	 * \param spline_xsize,spline_ysize,spline_zsize, x,y,z size of spline coefficients
	 * \param sz nxn size of the subregion to fit
	 * \param iterations number of iterations for solution to converge
	 * \param d_Parameters array of fitting parameters to return for each subregion
	 * \param d_CRLBs array of Cramer-Rao lower bound estimates to return for each subregion
	 * \param d_LogLikelihood array of loglikelihood estimates to return for each subregion
	 * \param initZ intial z position used for fitting
	 * \param Nfits number of subregions to fit
	 * \param d_varim variance map of sCMOS
	 */
	
   const int NV=NV_PSP;
    float M[NV*NV],Diag[NV], Minv[NV*NV];
    //int tx = threadIdx.x;
    //int bx = blockIdx.x;
    //int BlockSize = blockDim.x;
    int ii, jj, kk, ll, l, m, i;
	int xstart, ystart, zstart;

	const float *s_coeff;
	s_coeff = d_coeff;

    float model, data;
    float Div;
    float newTheta[NV],oldTheta[NV];
	float newLambda = INIT_LAMBDA, oldLambda = INIT_LAMBDA, mu;
	float newUpdate[NV] = {1e13, 1e13, 1e13, 1e13, 1e13},oldUpdate[NV] = {1e13, 1e13, 1e13, 1e13, 1e13};
	float maxJump[NV]={1.0,1.0,100,20,2};
	float newDudt[NV] ={0};

	float newErr = 1e12, oldErr = 1e13;

	float off;
	float jacobian[NV]={0};
	float hessian[NV*NV]={0};
	float t1,t2;

	float Nmax;
	float xc,yc,zc;
	float delta_f[64]={0}, delta_dxf[64]={0}, delta_dyf[64]={0}, delta_dzf[64]={0};
	int errFlag=0;
	float L[NV*NV] = {0}, U[NV*NV] = {0};

    
    //Prevent read/write past end of array
    if (subregion>=Nfits) return;
    
    for (ii=0;ii<NV*NV;ii++)M[ii]=0;
    for (ii=0;ii<NV*NV;ii++)Minv[ii]=0;

    //copy in data
      const float *s_data = d_data+(sz*sz*subregion);
	  const float *s_varim = d_varim+(sz*sz*subregion);
    
    //initial values
    kernel_CenterofMass2D(sz, s_data, &newTheta[0], &newTheta[1]);
    kernel_GaussFMaxMin2D(sz, 1.5, s_data, &Nmax, &newTheta[3]);

	//central pixel of spline model
	newTheta[3] = max(newTheta[3],0.01);
	newTheta[2]= (Nmax-newTheta[3])/d_coeff[(int)(spline_zsize/2)*(spline_xsize*spline_ysize)+(int)(spline_ysize/2)*spline_xsize+(int)(spline_xsize/2)]*4;

    newTheta[4]=initZ;

	maxJump[2]=max(newTheta[2],maxJump[2]);

	maxJump[3]=max(newTheta[3],maxJump[3]);

	maxJump[4]= max(spline_zsize/3.0f,maxJump[4]);


	for (ii=0;ii<NV;ii++)oldTheta[ii]=newTheta[ii];

	//updateFitValues
	xc = -1.0*((newTheta[0]-float(sz)/2)+0.5);
	yc = -1.0*((newTheta[1]-float(sz)/2)+0.5);

	off = floor((float(spline_xsize)+1.0-float(sz))/2);

	xstart = floor(xc);
	xc = xc-xstart;

	ystart = floor(yc);
	yc = yc-ystart;

	//zstart = floor(newTheta[4]);
	zstart = floor(newTheta[4]);
	zc = newTheta[4] -zstart;

	newErr = 0;
	memset(jacobian,0,NV*sizeof(float));
	memset(hessian,0,NV*NV*sizeof(float));
	kernel_computeDelta3D(xc, yc, zc, delta_f, delta_dxf, delta_dyf, delta_dzf);

	for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
		kernel_DerivativeSpline(ii+xstart+off,jj+ystart+off,zstart,spline_xsize,spline_ysize,spline_zsize,delta_f,delta_dxf,delta_dyf,delta_dzf,s_coeff,newTheta,newDudt,&model);
		model +=s_varim[sz*jj+ii];
		data=s_data[sz*jj+ii]+s_varim[sz*jj+ii];	

		if (data>0)
			newErr = newErr + 2*((model-data)-data*log(model/data));
		else
		{
			newErr = newErr + 2*model;
			data = 0;
		}

		t1 = 1-data/model;
		for (l=0;l<NV;l++){
			jacobian[l]+=t1*newDudt[l];
		}

		t2 = data/pow(model,2);
		for (l=0;l<NV;l++) for(m=l;m<NV;m++) {
			hessian[l*NV+m] +=t2*newDudt[l]*newDudt[m];
			hessian[m*NV+l] = hessian[l*NV+m];
		}
	}

	for (kk=0;kk<iterations;kk++) {//main iterative loop

			if(fabs((newErr-oldErr)/newErr)<TOLERANCE){
				//CONVERGED;
				break;
			}
			else{
				if(newErr>ACCEPTANCE*oldErr){
					//copy Fitdata
					for (i=0;i<NV;i++){
						newTheta[i]=oldTheta[i];
						newUpdate[i]=oldUpdate[i];
					}
					newLambda = oldLambda;
					newErr = oldErr;
					mu = max( (1 + newLambda*SCALE_UP)/(1 + newLambda),1.3f);         
					newLambda = SCALE_UP*newLambda;
				}
				else if(newErr<oldErr&&errFlag==0){
					newLambda = SCALE_DOWN*newLambda;
				    mu = 1+newLambda;
				}
				

				for (i=0;i<NV;i++){
					hessian[i*NV+i]=hessian[i*NV+i]*mu;
				}
				memset(L,0,NV*sizeof(float));
				memset(U,0,NV*sizeof(float));
				errFlag = kernel_cholesky(hessian,NV,L,U);
				if (errFlag ==0){
					for (i=0;i<NV;i++){
						oldTheta[i]=newTheta[i];
						oldUpdate[i] = newUpdate[i];
					}
					oldLambda = newLambda;
					oldErr=newErr;

					kernel_luEvaluate(L,U,jacobian,NV,newUpdate);	
					
					//updateFitParameters
					for (ll=0;ll<NV;ll++){
						if (newUpdate[ll]/oldUpdate[ll]< -0.5f){
							maxJump[ll] = maxJump[ll]*0.5;
						}
					    newUpdate[ll] = newUpdate[ll]/(1+fabs(newUpdate[ll]/maxJump[ll]));
						newTheta[ll] = newTheta[ll]-newUpdate[ll];
					}
					//restrict range
					newTheta[0] = max(newTheta[0],(float(sz)-1)/2-sz/4.0);
					newTheta[0] = min(newTheta[0],(float(sz)-1)/2+sz/4.0);
					newTheta[1] = max(newTheta[1],(float(sz)-1)/2-sz/4.0);
					newTheta[1] = min(newTheta[1],(float(sz)-1)/2+sz/4.0);
					newTheta[2] = max(newTheta[2],1.0);
					newTheta[3] = max(newTheta[3],0.01);
					newTheta[4] = max(newTheta[4],0.0);
					newTheta[4] = min(newTheta[4],float(spline_zsize));

					//updateFitValues
					xc = -1.0*((newTheta[0]-float(sz)/2)+0.5);
					yc = -1.0*((newTheta[1]-float(sz)/2)+0.5);

					xstart = floor(xc);
					xc = xc-xstart;

					ystart = floor(yc);
					yc = yc-ystart;

					zstart = floor(newTheta[4]);
					zc = newTheta[4] -zstart;


					newErr = 0;
					memset(jacobian,0,NV*sizeof(float));
					memset(hessian,0,NV*NV*sizeof(float));
					kernel_computeDelta3D(xc, yc, zc, delta_f, delta_dxf, delta_dyf, delta_dzf);
					for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
						kernel_DerivativeSpline(ii+xstart+off,jj+ystart+off,zstart,spline_xsize,spline_ysize,spline_zsize,delta_f,delta_dxf,delta_dyf,delta_dzf,s_coeff,newTheta,newDudt,&model);
						model +=s_varim[sz*jj+ii];
						data=s_data[sz*jj+ii]+s_varim[sz*jj+ii];	

						if (data>0)
							newErr = newErr + 2*((model-data)-data*log(model/data));
						else
						{
							newErr = newErr + 2*model;
							data = 0;
						}

						t1 = 1-data/model;
						for (l=0;l<NV;l++){
							jacobian[l]+=t1*newDudt[l];
						}

						t2 = data/pow(model,2);
						for (l=0;l<NV;l++) for(m=l;m<NV;m++) {
							hessian[l*NV+m] +=t2*newDudt[l]*newDudt[m];
							hessian[m*NV+l] = hessian[l*NV+m];
						}
					}
				}
				else
				{
					mu = max( (1 + newLambda*SCALE_UP)/(1 + newLambda),1.3f);         
					newLambda = SCALE_UP*newLambda;
				}

			}


		
	}
	//output iteration time
	d_Parameters[Nfits*NV+subregion]=kk;
    
    // Calculating the CRLB and LogLikelihood
	Div=0.0;

	xc = -1.0*((newTheta[0]-float(sz)/2)+0.5);
	yc = -1.0*((newTheta[1]-float(sz)/2)+0.5);

	//off = (float(spline_xsize)+1.0-2*float(sz))/2;

	xstart = floor(xc);
	xc = xc-xstart;

	ystart = floor(yc);
	yc = yc-ystart;

	//zstart = floor(newTheta[4]);
	zstart = floor(newTheta[4]);
	zc = newTheta[4] -zstart;

	kernel_computeDelta3D(xc, yc, zc, delta_f, delta_dxf, delta_dyf, delta_dzf);

    for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
		kernel_DerivativeSpline(ii+xstart+off,jj+ystart+off,zstart,spline_xsize,spline_ysize,spline_zsize,delta_f,delta_dxf,delta_dyf,delta_dzf,s_coeff,newTheta,newDudt,&model);
		model +=s_varim[sz*jj+ii];
		data=s_data[sz*jj+ii]+s_varim[sz*jj+ii];	
        
        //Building the Fisher Information Matrix
        for (kk=0;kk<NV;kk++)for (ll=kk;ll<NV;ll++){
            M[kk*NV+ll]+= newDudt[ll]*newDudt[kk]/model;
            M[ll*NV+kk]=M[kk*NV+ll];
        }
        
        //LogLikelyhood
        if (model>0)
            if (data>0)Div+=data*log(model)-model-data*log(data)+data;
            else
                Div+=-model;
    }
    
    // Matrix inverse (CRLB=F^-1) and output assigments
    kernel_MatInvN(M, Minv, Diag, NV);
  
    
    //write to global arrays
    for (kk=0;kk<NV;kk++) d_Parameters[Nfits*kk+subregion]=newTheta[kk];
   for (kk=0;kk<NV;kk++) d_CRLBs[Nfits*kk+subregion]=Diag[kk];
   d_LogLikelihood[subregion] = Div;
	//d_LogLikelihood[BlockSize*bx+tx] = 1;
    
    
    return;
}


