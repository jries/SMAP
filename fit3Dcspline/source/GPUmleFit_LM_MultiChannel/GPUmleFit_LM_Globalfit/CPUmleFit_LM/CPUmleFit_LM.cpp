
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

void kernel_MLEFit_sigma_EMCCD_multi(const int subregion, const float *d_data, const float PSFSigma, const float *d_dTAll,const int sz, const int iterations, const int NV, const int noChannels, float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,const int Nfits, const int*d_shared){
			//float test = 0;
			float M[25*Max_No_Channel*Max_No_Channel]={0},Diag[5*Max_No_Channel]={0},Minv[25*Max_No_Channel*Max_No_Channel]={0};
			//int tx = threadIdx.x;
			//int bx = blockIdx.x;
			//int BlockSize = blockDim.x;
			int ii,jj,kk,ll,l,m,i,j,n;
			//int xstart[Max_No_Channel]={0},ystart[Max_No_Channel]={0},zstart[Max_No_Channel]={0};


			float model[Max_No_Channel]={0},data[Max_No_Channel];
			float Div=0;

			float newTheta[5*Max_No_Channel]={0},oldTheta[5*Max_No_Channel]={0},newThetaAll[5*Max_No_Channel]={0};
			float newLambda = INIT_LAMBDA, oldLambda = INIT_LAMBDA, mu;
			float newUpdate[5*Max_No_Channel] ,oldUpdate[5*Max_No_Channel];
			float maxJump_Init[5*Max_No_Channel]={0},maxJump[5*Max_No_Channel]={0};
			float newDudt[5*Max_No_Channel]={0},newDudtAll[5*Max_No_Channel*Max_No_Channel];

			float newErr = 1e12, oldErr = 1e13;

			int off;
			float jacobian[5*Max_No_Channel]={0};
			float hessian[25*Max_No_Channel*Max_No_Channel]={0};
			float t1[Max_No_Channel]={0}, t2[Max_No_Channel]={0};

			float Nmax;
			/*float xc[Max_No_Channel]={0},yc[Max_No_Channel]={0},zc[Max_No_Channel]={0};
			float delta_f[64*Max_No_Channel]={0},delta_dxf[64*Max_No_Channel]={0},delta_dyf[64*Max_No_Channel]={0},delta_dzf[64*Max_No_Channel]={0};*/
			int errFlag = 0;
			float L[25*Max_No_Channel*Max_No_Channel]={0},U[25*Max_No_Channel*Max_No_Channel]={0};

			//Prevent read/write past end of array
			//if ((bx*BlockSize+tx)>=Nfits) return;
			if ((subregion)>=Nfits) return;

			for (i=0;i<NV;i++){
				newUpdate[i]=1e13;
				oldUpdate[i]=1e13;
			}
			//memset(newUpdate,1e13,NV*NV*sizeof(float));
			//memset(oldUpdate,1e13,NV*NV*sizeof(float));
			
			for (i=0;i<noChannels;i++){
				maxJump_Init[0+i*5] = 1.0f;
				maxJump_Init[1+i*5] = 1.0f;
				maxJump_Init[2+i*5] = 100.0f;
				maxJump_Init[3+i*5] = 20.0f;
				maxJump_Init[4+i*5] = 0.5f;
			}



			//copy in data
			const float *s_data[Max_No_Channel];
			const float *dT;
			const int *cshared;
			for (i=0;i<noChannels;i++){
				s_data[i]= d_data+(sz*sz*(Nfits*i+subregion));
			}
			dT = d_dTAll+ (subregion)*5*noChannels*2;
			cshared = d_shared+ (subregion)*5;

			//initial values
			for (i=0;i<noChannels;i++){
				kernel_CenterofMass2D(sz, s_data[i], &newTheta[0+5*i], &newTheta[1+5*i]);
				kernel_GaussFMaxMin2D(sz, 1.5, s_data[i], &Nmax, &newTheta[4+5*i]);
				newTheta[4+5*i] = PSFSigma;

				newTheta[2+5*i]= max(0.0, (Nmax-newTheta[3+5*i])*2*PI*PSFSigma*PSFSigma);
				newTheta[3+5*i] = max(newTheta[3+5*i],0.01);

				maxJump_Init[2+5*i] = max(maxJump_Init[2+5*i],newTheta[2+5*i]);
				maxJump_Init[3+5*i] = max(maxJump_Init[3+5*i],newTheta[3+5*i]);
			} 

			// map paramters
			for (i=0;i<5;i++){
				if (cshared[i]==1){
					for (j=1;j<noChannels;j++){
						newTheta[i+j*5]=newTheta[i] * dT[noChannels * 5 + i + j * 5] +dT[i+j*5];
					}
				}
			}

			for (ii=0;ii<5*noChannels;ii++)oldTheta[ii]=newTheta[ii];

			//combine all channel
			n = 0;
			for (i=0;i<5;i++){
				if (cshared[i]==1){
					maxJump[n]=maxJump_Init[i];
					newThetaAll[n]=newTheta[i];
				}
				else
				{
					for (j=0;j<noChannels;j++){
						maxJump[n+j]=maxJump_Init[i+j*5];
						newThetaAll[n+j]=newTheta[i+j*5];
					}
					n = n+j-1;
				}
				n=n+1;
			}

			newErr=0;

			memset(jacobian,0,NV*sizeof(float));
			memset(hessian,0,NV*NV*sizeof(float));
			memset(newDudtAll,0,NV*noChannels*sizeof(float));


			for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
				for (i=0;i<noChannels;i++){
					data[i] = *(s_data[i]+sz*jj+ii);
					kernel_DerivativeGauss2D_sigma(ii,jj,&newTheta[5*i],&newDudt[5*i],&model[i]);
					//s_data[i]= d_data+(sz*sz*bx*BlockSize+sz*sz*tx+sz*sz*Nfits*i);
				}

				n = 0;
				for (i=0;i<5;i++){
					if (cshared[i]==1){
						for (j=0;j<noChannels;j++){
							//newDudtAll[n+j*NV]=newDudt[i+j*5];
							newDudtAll[n + j * NV] = newDudt[i + j * 5] * dT[noChannels * 5 + i + j * 5];
						}	
					}
					else
					{
						for (j=0;j<noChannels;j++){
							newDudtAll[n+j+j*NV]=newDudt[i+j*5];
						}
						n = n+j-1;
					}
					n++;
				}

				for (i=0;i<noChannels;i++){ 
					if (data[i]>0)
						newErr = newErr + 2*((model[i]-data[i])-data[i]*log(model[i]/data[i]));
					else
					{
						newErr = newErr + 2*model[i];
						data[i] = 0;
					}
					t1[i] =1-data[i]/model[i];
					t2[i] = data[i]/pow(model[i],2);
				}

				for (l=0;l<NV;l++){
					for (i=0;i<noChannels;i++){ 
						jacobian[l]+=t1[i]*newDudtAll[l+i*NV];
					}
				}

				for (l=0;l<NV;l++) for(m=l;m<NV;m++) {
					for (i=0;i<noChannels;i++){ 
						hessian[l*NV+m] +=t2[i]*newDudtAll[l+i*NV]*newDudtAll[m+i*NV];
					}
					hessian[m*NV+l] = hessian[l*NV+m];
				}
			}
			
			for (kk=0;kk<iterations;kk++) {//main iterative loop

				if(fabs((newErr-oldErr)/newErr)<TOLERANCE){
					//newStatus = CONVERGED;
					break;
				}
				else{
					if(newErr>ACCEPTANCE*oldErr){
						//copy Fitdata
						for (ii=0;ii<5*noChannels;ii++)newTheta[ii]=oldTheta[ii];
						for (i=0;i<NV;i++){
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
					memset(L,0,NV*NV*sizeof(float));
					memset(U,0,NV*NV*sizeof(float));
					errFlag = kernel_cholesky(hessian,NV,L,U);
					if (errFlag ==0){
						for (ii=0;ii<5*noChannels;ii++)oldTheta[ii]=newTheta[ii];
						for (i=0;i<NV;i++){
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
							//test = oldUpdate[1];
							newUpdate[ll] = newUpdate[ll]/(1+fabs(newUpdate[ll]/maxJump[ll]));
							newThetaAll[ll] = newThetaAll[ll]-newUpdate[ll];
						}
						
						n=0;
						for (i=0;i<5;i++){
							if (cshared[i]==1){
								switch (i){
								case 0:
								case 1:
									/*newThetaAll[n]=max(newThetaAll[n],(float(sz)-1)/2-sz/4.0);
									newThetaAll[n] = min(newThetaAll[n],(float(sz)-1)/2+sz/4.0);*/
	/*								test = 1;*/
									break;
								case 2:
									/*newThetaAll[n]=max(newThetaAll[n],0.0);
									newThetaAll[n]=min(newThetaAll[n],float(spline_zsize));*/
									newThetaAll[n]=max(newThetaAll[n],1.0);
									break;
								case 3:
									/*newThetaAll[n]=max(newThetaAll[n],1.0);*/
									newThetaAll[n]=max(newThetaAll[n],0.01);
									break;
								case 4:
									//newThetaAll[n]=max(newThetaAll[n],0.01);
									newThetaAll[n]=max(newThetaAll[n],0.5);
									newThetaAll[n]=min(newThetaAll[n],sz/2.0f);
									/*newTheta[4] = max(newTheta[4],0.5);
									newTheta[4] = min(newTheta[4],sz/2.0f);*/
									
									break;
								}
								for (j=0;j<noChannels;j++){
									//newTheta[i+5*j]=newThetaAll[n]+dT[i+j*5];
									newTheta[i + 5 * j] = newThetaAll[n] * dT[noChannels * 5 + i + j * 5] + dT[i + j * 5];
								}
							}

							else
							{
								for (j=0;j<noChannels;j++){
									switch (i){
									case 0:
									case 1:

										/*newThetaAll[n+j]=max(newThetaAll[n+j],(float(sz)-1)/2-sz/4.0);
										newThetaAll[n+j] = min(newThetaAll[n+j],(float(sz)-1)/2+sz/4.0);*/


										break;
									case 2:

										/*newThetaAll[n+j]=max(newThetaAll[n+j],0.0);
										newThetaAll[n+j]=min(newThetaAll[n+j],float(spline_zsize));*/
										newThetaAll[n+j]=max(newThetaAll[n+j],1.0);

										break;
									case 3:

										/*newThetaAll[n+j]=max(newThetaAll[n+j],1.0);*/
										newThetaAll[n+j]=max(newThetaAll[n+j],0.01);
										break;
									case 4:
										//newThetaAll[n+j]=max(newThetaAll[n+j],0.01);

										newThetaAll[n+j]=max(newThetaAll[n+j],0.5);
										newThetaAll[n+j]=min(newThetaAll[n+j],sz/2.0f);

										break;
									}
									newTheta[i+5*j]=newThetaAll[n+j]+dT[i+j*5];
								}
								n = n+j-1;
							}
							n = n+1;
						}


						//prepare for next iteration
						newErr=0;

						memset(jacobian,0,NV*sizeof(float));
						memset(hessian,0,NV*NV*sizeof(float));
						memset(newDudtAll,0,NV*noChannels*sizeof(float));

						for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
							for (i=0;i<noChannels;i++){
								data[i] = *(s_data[i]+sz*jj+ii);
								kernel_DerivativeGauss2D_sigma(ii,jj,&newTheta[5*i],&newDudt[5*i],&model[i]);
							}

							n = 0;
							for (i=0;i<5;i++){
								if (cshared[i]==1){
									for (j=0;j<noChannels;j++){
										newDudtAll[n+j*NV]=newDudt[i+j*5];
									}	
								}
								else
								{
									for (j=0;j<noChannels;j++){
										newDudtAll[n+j+j*NV]=newDudt[i+j*5];
									}
									n = n+j-1;
								}
								n++;
							}

							for (i=0;i<noChannels;i++){ 
								if (data[i]>0)
									newErr = newErr + 2*((model[i]-data[i])-data[i]*log(model[i]/data[i]));
								else
								{
									newErr = newErr + 2*model[i];
									data[i] = 0;
								}
								t1[i] =1-data[i]/model[i];
								t2[i] = data[i]/pow(model[i],2);
							}

							for (l=0;l<NV;l++){
								for (i=0;i<noChannels;i++){ 
									jacobian[l]+=t1[i]*newDudtAll[l+i*NV];
								}
							}

							for (l=0;l<NV;l++) for(m=l;m<NV;m++) {
								for (i=0;i<noChannels;i++){ 
									hessian[l*NV+m] +=t2[i]*newDudtAll[l+i*NV]*newDudtAll[m+i*NV];
								}
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
			/*d_Parameters[Nfits*(NV+1)+BlockSize*bx+tx]=maxJump[0];
			d_Parameters[Nfits*(NV+2)+BlockSize*bx+tx]=maxJump[1];
			d_Parameters[Nfits*(NV+3)+BlockSize*bx+tx]=maxJump[2];
			d_Parameters[Nfits*(NV+4)+BlockSize*bx+tx]=maxJump[3];
			d_Parameters[Nfits*(NV+5)+BlockSize*bx+tx]=maxJump[4];
			d_Parameters[Nfits*(NV+6)+BlockSize*bx+tx]=maxJump[5];
			d_Parameters[Nfits*(NV+7)+BlockSize*bx+tx]=maxJump[6];*/
			//d_Parameters[Nfits*(NV+4)+BlockSize*bx+tx]=maxJump[3];


			// Calculating the CRLB and LogLikelihood
			Div=0.0;
			
			for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
				for (i=0;i<noChannels;i++){
					data[i] = *(s_data[i]+sz*jj+ii);
					kernel_DerivativeGauss2D_sigma(ii,jj,&newTheta[5*i],&newDudt[5*i],&model[i]);
				}

				n = 0;
				for (i=0;i<5;i++){
					if (cshared[i]==1){
						for (j=0;j<noChannels;j++){
							//newDudtAll[n+j*NV]=newDudt[i+j*5];
							newDudtAll[n + j * NV] = newDudt[i + j * 5] * dT[noChannels * 5 + i + j * 5];
						}	
					}
					else
					{
						for (j=0;j<noChannels;j++){
							newDudtAll[n+j+j*NV]=newDudt[i+j*5];
						}
						n = n+j-1;
					}
					n++;
				}

				//Building the Fisher Information Matrix
				for (kk=0;kk<NV;kk++)for (ll=kk;ll<NV;ll++){
					for (j=0;j<noChannels;j++){
						M[kk*NV+ll]+= newDudtAll[ll+j*NV]*newDudtAll[kk+j*NV]/model[j];
					}
					M[ll*NV+kk]=M[kk*NV+ll];
				}
        
				//LogLikelyhood
				for(i=0;i<noChannels;i++){
				if (model[i]>0)
					if (data[i]>0)Div+=data[i]*log(model[i])-model[i]-data[i]*log(data[i])+data[i];
					else
						Div+=-model[i];
				}
			}
			// Matrix inverse (CRLB=F^-1) and output assigments
			kernel_MatInvN(M, Minv, Diag, NV);

			/*n = 0;
			for (i=0;i<5;i++){
				if (cshared[i]!=1)
				{
					for (j=0;j<noChannels;j++){
						newThetaAll[n+j]=newThetaAll[n+j]+dT[i+j*5];
					}
					n = n+j-1;
				}
				n=n+1;
			}*/

			//write to global arrays
			for (kk=0;kk<NV;kk++) d_Parameters[Nfits*kk+subregion]=newThetaAll[kk];
			for (kk=0;kk<NV;kk++) d_CRLBs[Nfits*kk+subregion]=Diag[kk];
			//d_LogLikelihood[BlockSize*bx+tx] = newUpdate[0];
			d_LogLikelihood[subregion] = Div;
			//d_LogLikelihood[BlockSize*bx+tx] = 1;


			return;


}
//**********************************************************************************************************************************
//multichannel fit
void kernel_splineMLEFit_z_EMCCD_multi(const int subregion, const float *d_data,const float *d_coeff, const float *d_dTAll, const int spline_xsize, const int spline_ysize, const int spline_zsize, const int sz, const int iterations, const int NV, const int noChannels, float *d_Parameters, float *d_CRLBs, float *d_LogLikelihood,float *initZ, const int Nfits, const int *d_shared){
	float M[25*Max_No_Channel*Max_No_Channel]={0},Diag[5*Max_No_Channel]={0},Minv[25*Max_No_Channel*Max_No_Channel]={0};
	//int tx = threadIdx.x;
	//int bx = blockIdx.x;
	//int BlockSize = blockDim.x;
	int ii,jj,kk,ll,l,m,i,j,n;
	int xstart[Max_No_Channel]={0},ystart[Max_No_Channel]={0},zstart[Max_No_Channel]={0};


	float model[Max_No_Channel]={0},data[Max_No_Channel];
	float Div=0;

	float newTheta[5*Max_No_Channel]={0},oldTheta[5*Max_No_Channel]={0},newThetaAll[5*Max_No_Channel]={0};
	float newLambda = INIT_LAMBDA, oldLambda = INIT_LAMBDA, mu;
	float newUpdate[5*Max_No_Channel] ,oldUpdate[5*Max_No_Channel];
	float maxJump_Init[5*Max_No_Channel]={0},maxJump[5*Max_No_Channel]={0};
	float newDudt[5*Max_No_Channel]={0},newDudtAll[5*Max_No_Channel*Max_No_Channel];

	float newErr = 1e12, oldErr = 1e13;

	int off;
	float jacobian[5*Max_No_Channel]={0};
	float hessian[25*Max_No_Channel*Max_No_Channel]={0};
	float t1[Max_No_Channel]={0}, t2[Max_No_Channel]={0};

	float Nmax;
	float xc[Max_No_Channel]={0},yc[Max_No_Channel]={0},zc[Max_No_Channel]={0};
	float delta_f[64*Max_No_Channel]={0},delta_dxf[64*Max_No_Channel]={0},delta_dyf[64*Max_No_Channel]={0},delta_dzf[64*Max_No_Channel]={0};
	int errFlag = 0;
	float L[25*Max_No_Channel*Max_No_Channel]={0},U[25*Max_No_Channel*Max_No_Channel]={0};

	//Prevent read/write past end of array
	if ((subregion)>=Nfits) return;

	for (i=0;i<NV;i++){
		newUpdate[i]=1e13;
		oldUpdate[i]=1e13;
	}
	//memset(newUpdate,1e13,NV*NV*sizeof(float));
	//memset(oldUpdate,1e13,NV*NV*sizeof(float));

	for (i=0;i<noChannels;i++){
		maxJump_Init[0+i*5] = 1.0f;
		maxJump_Init[1+i*5] = 1.0f;
		maxJump_Init[2+i*5] = spline_zsize/5.0f;
		maxJump_Init[3+i*5] = 100.0f;
		maxJump_Init[4+i*5] = 20.0f;
	}



	//copy in data
	const float *s_data[Max_No_Channel];
	const float *dT;
	const int *cshared;
	for (i=0;i<noChannels;i++){
		s_data[i]= d_data+(sz*sz*(Nfits*i+subregion));
	}
	dT = d_dTAll+ (subregion)*5*noChannels*2;
	cshared = d_shared+ (subregion)*5;

	//initial values
	for (i=0;i<noChannels;i++){
		kernel_CenterofMass2D(sz, s_data[i], &newTheta[0+5*i], &newTheta[1+5*i]);
		kernel_GaussFMaxMin2D(sz, 1.5, s_data[i], &Nmax, &newTheta[4+5*i]);
		newTheta[2+5*i] = initZ[subregion];

		newTheta[3+5*i]= (Nmax-newTheta[3+5*i])/d_coeff[(int)(spline_zsize/2)*(spline_xsize*spline_ysize)+(int)(spline_ysize/2)*spline_xsize+(int)(spline_xsize/2)+i*spline_xsize*spline_ysize*spline_zsize*64]*4;
		newTheta[4+5*i] = max(newTheta[4+5*i],0.01);

		maxJump_Init[3+5*i] = max(maxJump_Init[3+5*i],newTheta[3+5*i]);
		maxJump_Init[4+5*i] = max(maxJump_Init[4+5*i],newTheta[4+5*i]);
	} 

	// map paramters
	for (i=0;i<5;i++){
		if (cshared[i]==1){
			for (j=1;j<noChannels;j++){
				//newTheta[i+j*5]=newTheta[i]+dT[i+j*5];
				newTheta[i + j * 5] = newTheta[i] * dT[noChannels * 5 + i + j * 5] + dT[i + j * 5];

			}
		}
	}

	for (ii=0;ii<5*noChannels;ii++)oldTheta[ii]=newTheta[ii];

	//combine all channel
	n = 0;
	for (i=0;i<5;i++){
		if (cshared[i]==1){
			maxJump[n]=maxJump_Init[i];
			newThetaAll[n]=newTheta[i];
		}
		else
		{
			for (j=0;j<noChannels;j++){
				maxJump[n+j]=maxJump_Init[i+j*5];
				newThetaAll[n+j]=newTheta[i+j*5];
			}
			n = n+j-1;
		}
		n=n+1;
	}



	//prepare values for iteration
	off = (int)((float(spline_xsize)+1.0-float(sz))/2);

	for (i=0;i<noChannels;i++){ 
		xc[i] = -1.0*((newTheta[0+5*i]-float(sz)/2)+0.5);
		yc[i] = -1.0*((newTheta[1+5*i]-float(sz)/2)+0.5);

		xstart[i] = floor(xc[i]);
		xc[i] = xc[i]-xstart[i];

		ystart[i] = floor(yc[i]);
		yc[i] = yc[i]-ystart[i];

		//zstart = floor(newTheta[4]);
		zstart[i] = floor(newTheta[2+5*i]);
		zc[i] = newTheta[2+5*i] -zstart[i];
	}

	newErr=0;

	memset(jacobian,0,NV*sizeof(float));
	memset(hessian,0,NV*NV*sizeof(float));
	memset(newDudtAll,0,NV*noChannels*sizeof(float));

	for (i=0;i<noChannels;i++){ 
		kernel_computeDelta3D(xc[i], yc[i], zc[i], &delta_f[64*i], &delta_dxf[64*i], &delta_dyf[64*i], &delta_dzf[64*i]);
	}

	for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
		for (i=0;i<noChannels;i++){
			data[i] = *(s_data[i]+sz*jj+ii);
			kernel_DerivativeSpline_multiChannel(ii+xstart[i]+off,jj+ystart[i]+off,zstart[i],spline_xsize,spline_ysize,spline_zsize,&delta_f[64*i],&delta_dxf[64*i],&delta_dyf[64*i],&delta_dzf[64*i],&d_coeff[spline_xsize*spline_ysize*spline_zsize*64*i],&newTheta[5*i],&newDudt[5*i],&model[i]);
			//s_data[i]= d_data+(sz*sz*bx*BlockSize+sz*sz*tx+sz*sz*Nfits*i);
		}

		n = 0;
		for (i=0;i<5;i++){
			if (cshared[i]==1){
				for (j=0;j<noChannels;j++){
					//newDudtAll[n+j*NV]=newDudt[i+j*5];
					newDudtAll[n + j * NV] = newDudt[i + j * 5] * dT[noChannels * 5 + i + j * 5];
				}	
			}
			else
			{
				for (j=0;j<noChannels;j++){
					newDudtAll[n+j+j*NV]=newDudt[i+j*5];
				}
				n = n+j-1;
			}
			n++;
		}

		for (i=0;i<noChannels;i++){ 
			if (data[i]>0)
				newErr = newErr + 2*((model[i]-data[i])-data[i]*log(model[i]/data[i]));
			else
			{
				newErr = newErr + 2*model[i];
				data[i] = 0;
			}
			t1[i] =1-data[i]/model[i];
			t2[i] = data[i]/pow(model[i],2);
		}

		for (l=0;l<NV;l++){
			for (i=0;i<noChannels;i++){ 
				jacobian[l]+=t1[i]*newDudtAll[l+i*NV];
			}
		}

		for (l=0;l<NV;l++) for(m=l;m<NV;m++) {
			for (i=0;i<noChannels;i++){ 
				hessian[l*NV+m] +=t2[i]*newDudtAll[l+i*NV]*newDudtAll[m+i*NV];
			}
			hessian[m*NV+l] = hessian[l*NV+m];
		}
	}

	for (kk=0;kk<iterations;kk++) {//main iterative loop

		if(fabs((newErr-oldErr)/newErr)<TOLERANCE){
			//newStatus = CONVERGED;
			break;
		}
		else{
			if(newErr>ACCEPTANCE*oldErr){
				//copy Fitdata
				for (ii=0;ii<5*noChannels;ii++)newTheta[ii]=oldTheta[ii];
				for (i=0;i<NV;i++){
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
			memset(L,0,NV*NV*sizeof(float));
			memset(U,0,NV*NV*sizeof(float));
			
			errFlag = kernel_cholesky(hessian,NV,L,U);
			if (errFlag ==0){
				for (ii=0;ii<5*noChannels;ii++)oldTheta[ii]=newTheta[ii];
				for (i=0;i<NV;i++){
					oldUpdate[i] = newUpdate[i];
				}
				oldLambda = newLambda;
				oldErr=newErr;
				
				
				//mexPrintf("before kernel_luEvaluate\n");
				//mexPrintf("before kernel_luEvaluate NV %f\n",(double)NV);
				//for (ii=0;ii<100;ii++){
				//	mexPrintf("before L is %f, i is  %f\n",(double)L[ii],(double)ii);
				//	mexPrintf("before U is %f, i is  %f\n",(double)U[ii],(double)ii);
				//}
				
				kernel_luEvaluate(L,U,jacobian,NV,newUpdate);
				//return;

				//updateFitParameters
				for (ll=0;ll<NV;ll++){
					if (newUpdate[ll]/oldUpdate[ll]< -0.5f){
						maxJump[ll] = maxJump[ll]*0.5;

					}
					//test = oldUpdate[1];
					newUpdate[ll] = newUpdate[ll]/(1+fabs(newUpdate[ll]/maxJump[ll]));
					newThetaAll[ll] = newThetaAll[ll]-newUpdate[ll];
				}
				
				n=0;
				for (i=0;i<5;i++){
					if (cshared[i]==1){
						switch (i){
						case 0:
						case 1:
							newThetaAll[n]=max(newThetaAll[n],(float(sz)-1)/2-sz/4.0);
							newThetaAll[n] = min(newThetaAll[n],(float(sz)-1)/2+sz/4.0);
							/*								test = 1;*/
							break;
						case 2:
							newThetaAll[n]=max(newThetaAll[n],0.0);
							newThetaAll[n]=min(newThetaAll[n],float(spline_zsize));
							break;
						case 3:
							newThetaAll[n]=max(newThetaAll[n],1.0);
							break;
						case 4:
							newThetaAll[n]=max(newThetaAll[n],0.01);
							break;
						}
						for (j=0;j<noChannels;j++){
							//newTheta[i+5*j]=newThetaAll[n]+dT[i+j*5];
							newTheta[i + 5 * j] = newThetaAll[n] * dT[noChannels * 5 + i + j * 5] + dT[i + j * 5];
						}
					}

					else
					{
						for (j=0;j<noChannels;j++){
							switch (i){
							case 0:
							case 1:
								/*newThetaAll[n]=max(newThetaAll[n],(float(sz)-1)/2-sz/4.0);
								newThetaAll[n] = min(newThetaAll[n],(float(sz)-1)/2+sz/4.0);*/

								newThetaAll[n+j]=max(newThetaAll[n+j],(float(sz)-1)/2-sz/4.0);
								newThetaAll[n+j] = min(newThetaAll[n+j],(float(sz)-1)/2+sz/4.0);


								break;
							case 2:
								/*newThetaAll[n]=max(newThetaAll[n],0.0);
								newThetaAll[n]=min(newThetaAll[n],float(spline_zsize));*/

								newThetaAll[n+j]=max(newThetaAll[n+j],0.0);
								newThetaAll[n+j]=min(newThetaAll[n+j],float(spline_zsize));

								break;
							case 3:
								/*newThetaAll[n]=max(newThetaAll[n],1.0);*/

								newThetaAll[n+j]=max(newThetaAll[n+j],1.0);
								break;
							case 4:
								/*newThetaAll[n]=max(newThetaAll[n],0.01);*/

								newThetaAll[n+j]=max(newThetaAll[n+j],0.01);
								break;
							}
							newTheta[i+5*j]=newThetaAll[n+j]+dT[i+j*5];
						}
						n = n+j-1;
					}
					n = n+1;
				}


				//prepare for next iteration
				for (i=0;i<noChannels;i++){ 
					xc[i] = -1.0*((newTheta[0+5*i]-float(sz)/2)+0.5);
					yc[i] = -1.0*((newTheta[1+5*i]-float(sz)/2)+0.5);

					xstart[i] = floor(xc[i]);
					xc[i] = xc[i]-xstart[i];

					ystart[i] = floor(yc[i]);
					yc[i] = yc[i]-ystart[i];

					//zstart = floor(newTheta[4]);
					zstart[i] = floor(newTheta[2+5*i]);
					zc[i] = newTheta[2+5*i] -zstart[i];
				}

				newErr=0;

				memset(jacobian,0,NV*sizeof(float));
				memset(hessian,0,NV*NV*sizeof(float));
				memset(newDudtAll,0,NV*noChannels*sizeof(float));

				for (i=0;i<noChannels;i++){ 
					kernel_computeDelta3D(xc[i], yc[i], zc[i], &delta_f[64*i], &delta_dxf[64*i], &delta_dyf[64*i], &delta_dzf[64*i]);
				}


				for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
					for (i=0;i<noChannels;i++){
						data[i] = *(s_data[i]+sz*jj+ii);
						kernel_DerivativeSpline_multiChannel(ii+xstart[i]+off,jj+ystart[i]+off,zstart[i],spline_xsize,spline_ysize,spline_zsize,&delta_f[64*i],&delta_dxf[64*i],&delta_dyf[64*i],&delta_dzf[64*i],&d_coeff[spline_xsize*spline_ysize*spline_zsize*64*i],&newTheta[5*i],&newDudt[5*i],&model[i]);
						//s_data[i]= d_data+(sz*sz*bx*BlockSize+sz*sz*tx+sz*sz*Nfits*i);
					}

					n = 0;
					for (i=0;i<5;i++){
						if (cshared[i]==1){
							for (j=0;j<noChannels;j++){
								//newDudtAll[n+j*NV]=newDudt[i+j*5];
								newDudtAll[n + j * NV] = newDudt[i + j * 5] * dT[noChannels * 5 + i + j * 5];
							}	
						}
						else
						{
							for (j=0;j<noChannels;j++){
								newDudtAll[n+j+j*NV]=newDudt[i+j*5];
							}
							n = n+j-1;
						}
						n++;
					}

					for (i=0;i<noChannels;i++){ 
						if (data[i]>0)
							newErr = newErr + 2*((model[i]-data[i])-data[i]*log(model[i]/data[i]));
						else
						{
							newErr = newErr + 2*model[i];
							data[i] = 0;
						}
						t1[i] =1-data[i]/model[i];
						t2[i] = data[i]/pow(model[i],2);
					}

					for (l=0;l<NV;l++){
						for (i=0;i<noChannels;i++){ 
							jacobian[l]+=t1[i]*newDudtAll[l+i*NV];
						}
					}

					for (l=0;l<NV;l++) for(m=l;m<NV;m++) {
						for (i=0;i<noChannels;i++){ 
							hessian[l*NV+m] +=t2[i]*newDudtAll[l+i*NV]*newDudtAll[m+i*NV];
						}
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
	/*d_Parameters[Nfits*(NV+1)+BlockSize*bx+tx]=maxJump[0];
	d_Parameters[Nfits*(NV+2)+BlockSize*bx+tx]=maxJump[1];
	d_Parameters[Nfits*(NV+3)+BlockSize*bx+tx]=maxJump[2];
	d_Parameters[Nfits*(NV+4)+BlockSize*bx+tx]=maxJump[3];
	d_Parameters[Nfits*(NV+5)+BlockSize*bx+tx]=maxJump[4];
	d_Parameters[Nfits*(NV+6)+BlockSize*bx+tx]=maxJump[5];
	d_Parameters[Nfits*(NV+7)+BlockSize*bx+tx]=maxJump[6];*/
	//d_Parameters[Nfits*(NV+4)+BlockSize*bx+tx]=maxJump[3];


	// Calculating the CRLB and LogLikelihood
	Div=0.0;
	for (i=0;i<noChannels;i++){ 
		xc[i] = -1.0*((newTheta[0+5*i]-float(sz)/2)+0.5);
		yc[i] = -1.0*((newTheta[1+5*i]-float(sz)/2)+0.5);

		xstart[i] = floor(xc[i]);
		xc[i] = xc[i]-xstart[i];

		ystart[i] = floor(yc[i]);
		yc[i] = yc[i]-ystart[i];

		//zstart = floor(newTheta[4]);
		zstart[i] = floor(newTheta[2+5*i]);
		zc[i] = newTheta[2+5*i] -zstart[i];
	}


	for (i=0;i<noChannels;i++){ 
		kernel_computeDelta3D(xc[i], yc[i], zc[i], &delta_f[64*i], &delta_dxf[64*i], &delta_dyf[64*i], &delta_dzf[64*i]);
	}

	for (ii=0;ii<sz;ii++) for(jj=0;jj<sz;jj++) {
		for (i=0;i<noChannels;i++){
			data[i] = *(s_data[i]+sz*jj+ii);
			kernel_DerivativeSpline_multiChannel(ii+xstart[i]+off,jj+ystart[i]+off,zstart[i],spline_xsize,spline_ysize,spline_zsize,&delta_f[64*i],&delta_dxf[64*i],&delta_dyf[64*i],&delta_dzf[64*i],&d_coeff[spline_xsize*spline_ysize*spline_zsize*64*i],&newTheta[5*i],&newDudt[5*i],&model[i]);
			//s_data[i]= d_data+(sz*sz*bx*BlockSize+sz*sz*tx+sz*sz*Nfits*i);
		}

		n = 0;
		for (i=0;i<5;i++){
			if (cshared[i]==1){
				for (j=0;j<noChannels;j++){
					//newDudtAll[n+j*NV]=newDudt[i+j*5];
					newDudtAll[n + j * NV] = newDudt[i + j * 5] * dT[noChannels * 5 + i + j * 5];

				}	
			}
			else
			{
				for (j=0;j<noChannels;j++){
					newDudtAll[n+j+j*NV]=newDudt[i+j*5];
				}
				n = n+j-1;
			}
			n++;
		}

		//Building the Fisher Information Matrix
		for (kk=0;kk<NV;kk++)for (ll=kk;ll<NV;ll++){
			for (j=0;j<noChannels;j++){
				M[kk*NV+ll]+= newDudtAll[ll+j*NV]*newDudtAll[kk+j*NV]/model[j];
			}
			M[ll*NV+kk]=M[kk*NV+ll];
		}

		//LogLikelyhood
		for(i=0;i<noChannels;i++){
			if (model[i]>0)
				if (data[i]>0)Div+=data[i]*log(model[i])-model[i]-data[i]*log(data[i])+data[i];
				else
					Div+=-model[i];
		}
	}
	// Matrix inverse (CRLB=F^-1) and output assigments
	kernel_MatInvN(M, Minv, Diag, NV);

	//write to global arrays
	for (kk=0;kk<NV;kk++) d_Parameters[Nfits*kk+subregion]=newThetaAll[kk];
	for (kk=0;kk<NV;kk++) d_CRLBs[Nfits*kk+subregion]=Diag[kk];
	//d_LogLikelihood[BlockSize*bx+tx] = newUpdate[0];
	d_LogLikelihood[subregion] = Div;
	//d_LogLikelihood[BlockSize*bx+tx] = 1;


	return;
}



//sCMOS
//*******************************************************************************************
void kernel_MLEFit_sigma_sCMOS_multi(const int subregion, const float* d_data, const float PSFSigma, const float* d_dTAll, const int sz, const int iterations, const int NV, const int noChannels, float* d_Parameters, float* d_CRLBs, float* d_LogLikelihood, const int Nfits, const int* d_shared, const float* d_varim) {
	//float test = 0;
	float M[25 * Max_No_Channel * Max_No_Channel] = { 0 }, Diag[5 * Max_No_Channel] = { 0 }, Minv[25 * Max_No_Channel * Max_No_Channel] = { 0 };
	//int tx = threadIdx.x;
	//int bx = blockIdx.x;
	//int BlockSize = blockDim.x;
	int ii, jj, kk, ll, l, m, i, j, n;
	//int xstart[Max_No_Channel]={0},ystart[Max_No_Channel]={0},zstart[Max_No_Channel]={0};


	float model[Max_No_Channel] = { 0 }, data[Max_No_Channel];
	float Div = 0;

	float newTheta[5 * Max_No_Channel] = { 0 }, oldTheta[5 * Max_No_Channel] = { 0 }, newThetaAll[5 * Max_No_Channel] = { 0 };
	float newLambda = INIT_LAMBDA, oldLambda = INIT_LAMBDA, mu;
	float newUpdate[5 * Max_No_Channel], oldUpdate[5 * Max_No_Channel];
	float maxJump_Init[5 * Max_No_Channel] = { 0 }, maxJump[5 * Max_No_Channel] = { 0 };
	float newDudt[5 * Max_No_Channel] = { 0 }, newDudtAll[5 * Max_No_Channel * Max_No_Channel];

	float newErr = 1e12, oldErr = 1e13;

	int off;
	float jacobian[5 * Max_No_Channel] = { 0 };
	float hessian[25 * Max_No_Channel * Max_No_Channel] = { 0 };
	float t1[Max_No_Channel] = { 0 }, t2[Max_No_Channel] = { 0 };

	float Nmax;
	/*float xc[Max_No_Channel]={0},yc[Max_No_Channel]={0},zc[Max_No_Channel]={0};
	float delta_f[64*Max_No_Channel]={0},delta_dxf[64*Max_No_Channel]={0},delta_dyf[64*Max_No_Channel]={0},delta_dzf[64*Max_No_Channel]={0};*/
	int errFlag = 0;
	float L[25 * Max_No_Channel * Max_No_Channel] = { 0 }, U[25 * Max_No_Channel * Max_No_Channel] = { 0 };

	//Prevent read/write past end of array
	//if ((bx*BlockSize+tx)>=Nfits) return;
	if ((subregion) >= Nfits) return;

	for (i = 0; i < NV; i++) {
		newUpdate[i] = 1e13;
		oldUpdate[i] = 1e13;
	}
	//memset(newUpdate,1e13,NV*NV*sizeof(float));
	//memset(oldUpdate,1e13,NV*NV*sizeof(float));

	for (i = 0; i < noChannels; i++) {
		maxJump_Init[0 + i * 5] = 1.0f;
		maxJump_Init[1 + i * 5] = 1.0f;
		maxJump_Init[2 + i * 5] = 100.0f;
		maxJump_Init[3 + i * 5] = 20.0f;
		maxJump_Init[4 + i * 5] = 0.5f;
	}



	//copy in data
	const float* s_data[Max_No_Channel];
	const float* s_varim[Max_No_Channel];
	const float* dT;
	const int* cshared;
	for (i = 0; i < noChannels; i++) {
		s_data[i] = d_data + (sz * sz * (Nfits * i + subregion));
	}

	for (i = 0; i < noChannels; i++) {
		s_varim[i] = d_varim + (sz * sz * (Nfits * i + subregion));
	}

	dT = d_dTAll + (subregion) * 5 * noChannels * 2;
	cshared = d_shared + (subregion) * 5;

	//initial values
	for (i = 0; i < noChannels; i++) {
		kernel_CenterofMass2D(sz, s_data[i], &newTheta[0 + 5 * i], &newTheta[1 + 5 * i]);
		kernel_GaussFMaxMin2D(sz, 1.5, s_data[i], &Nmax, &newTheta[4 + 5 * i]);
		newTheta[4 + 5 * i] = PSFSigma;

		newTheta[2 + 5 * i] = max(0.0, (Nmax - newTheta[3 + 5 * i]) * 2 * PI * PSFSigma * PSFSigma);
		newTheta[3 + 5 * i] = max(newTheta[3 + 5 * i], 0.01);

		maxJump_Init[2 + 5 * i] = max(maxJump_Init[2 + 5 * i], newTheta[2 + 5 * i]);
		maxJump_Init[3 + 5 * i] = max(maxJump_Init[3 + 5 * i], newTheta[3 + 5 * i]);
	}

	// map paramters
	for (i = 0; i < 5; i++) {
		if (cshared[i] == 1) {
			for (j = 1; j < noChannels; j++) {
				newTheta[i + j * 5] = newTheta[i] * dT[noChannels * 5 + i + j * 5] + dT[i + j * 5];
			}
		}
	}

	for (ii = 0; ii < 5 * noChannels; ii++)oldTheta[ii] = newTheta[ii];

	//combine all channel
	n = 0;
	for (i = 0; i < 5; i++) {
		if (cshared[i] == 1) {
			maxJump[n] = maxJump_Init[i];
			newThetaAll[n] = newTheta[i];
		}
		else
		{
			for (j = 0; j < noChannels; j++) {
				maxJump[n + j] = maxJump_Init[i + j * 5];
				newThetaAll[n + j] = newTheta[i + j * 5];
			}
			n = n + j - 1;
		}
		n = n + 1;
	}

	newErr = 0;

	memset(jacobian, 0, NV * sizeof(float));
	memset(hessian, 0, NV * NV * sizeof(float));
	memset(newDudtAll, 0, NV * noChannels * sizeof(float));


	for (ii = 0; ii < sz; ii++) for (jj = 0; jj < sz; jj++) {
		for (i = 0; i < noChannels; i++) {
			kernel_DerivativeGauss2D_sigma(ii, jj, &newTheta[5 * i], &newDudt[5 * i], &model[i]);
			data[i] = *(s_data[i] + sz * jj + ii) + *(s_varim[i] + sz * jj + ii);
			model[i] += *(s_varim[i] + sz * jj + ii);
			//s_data[i]= d_data+(sz*sz*bx*BlockSize+sz*sz*tx+sz*sz*Nfits*i);
		}

		n = 0;
		for (i = 0; i < 5; i++) {
			if (cshared[i] == 1) {
				for (j = 0; j < noChannels; j++) {
					//newDudtAll[n+j*NV]=newDudt[i+j*5];
					newDudtAll[n + j * NV] = newDudt[i + j * 5] * dT[noChannels * 5 + i + j * 5];
				}
			}
			else
			{
				for (j = 0; j < noChannels; j++) {
					newDudtAll[n + j + j * NV] = newDudt[i + j * 5];
				}
				n = n + j - 1;
			}
			n++;
		}

		for (i = 0; i < noChannels; i++) {
			if (data[i] > 0)
				newErr = newErr + 2 * ((model[i] - data[i]) - data[i] * log(model[i] / data[i]));
			else
			{
				newErr = newErr + 2 * model[i];
				data[i] = 0;
			}
			t1[i] = 1 - data[i] / model[i];
			t2[i] = data[i] / pow(model[i], 2);
		}

		for (l = 0; l < NV; l++) {
			for (i = 0; i < noChannels; i++) {
				jacobian[l] += t1[i] * newDudtAll[l + i * NV];
			}
		}

		for (l = 0; l < NV; l++) for (m = l; m < NV; m++) {
			for (i = 0; i < noChannels; i++) {
				hessian[l * NV + m] += t2[i] * newDudtAll[l + i * NV] * newDudtAll[m + i * NV];
			}
			hessian[m * NV + l] = hessian[l * NV + m];
		}
	}

	for (kk = 0; kk < iterations; kk++) {//main iterative loop

		if (fabs((newErr - oldErr) / newErr) < TOLERANCE) {
			//newStatus = CONVERGED;
			break;
		}
		else {
			if (newErr > ACCEPTANCE* oldErr) {
				//copy Fitdata
				for (ii = 0; ii < 5 * noChannels; ii++)newTheta[ii] = oldTheta[ii];
				for (i = 0; i < NV; i++) {
					newUpdate[i] = oldUpdate[i];
				}
				newLambda = oldLambda;
				newErr = oldErr;
				mu = max((1 + newLambda * SCALE_UP) / (1 + newLambda), 1.3f);
				newLambda = SCALE_UP * newLambda;
			}
			else if (newErr < oldErr && errFlag == 0) {
				newLambda = SCALE_DOWN * newLambda;
				mu = 1 + newLambda;
			}


			for (i = 0; i < NV; i++) {
				hessian[i * NV + i] = hessian[i * NV + i] * mu;
			}
			memset(L, 0, NV * NV * sizeof(float));
			memset(U, 0, NV * NV * sizeof(float));
			errFlag = kernel_cholesky(hessian, NV, L, U);
			if (errFlag == 0) {
				for (ii = 0; ii < 5 * noChannels; ii++)oldTheta[ii] = newTheta[ii];
				for (i = 0; i < NV; i++) {
					oldUpdate[i] = newUpdate[i];
				}
				oldLambda = newLambda;
				oldErr = newErr;

				kernel_luEvaluate(L, U, jacobian, NV, newUpdate);

				//updateFitParameters
				for (ll = 0; ll < NV; ll++) {
					if (newUpdate[ll] / oldUpdate[ll] < -0.5f) {
						maxJump[ll] = maxJump[ll] * 0.5;

					}
					//test = oldUpdate[1];
					newUpdate[ll] = newUpdate[ll] / (1 + fabs(newUpdate[ll] / maxJump[ll]));
					newThetaAll[ll] = newThetaAll[ll] - newUpdate[ll];
				}

				n = 0;
				for (i = 0; i < 5; i++) {
					if (cshared[i] == 1) {
						switch (i) {
						case 0:
						case 1:
							/*newThetaAll[n]=max(newThetaAll[n],(float(sz)-1)/2-sz/4.0);
							newThetaAll[n] = min(newThetaAll[n],(float(sz)-1)/2+sz/4.0);*/
							/*								test = 1;*/
							break;
						case 2:
							/*newThetaAll[n]=max(newThetaAll[n],0.0);
							newThetaAll[n]=min(newThetaAll[n],float(spline_zsize));*/
							newThetaAll[n] = max(newThetaAll[n], 1.0);
							break;
						case 3:
							/*newThetaAll[n]=max(newThetaAll[n],1.0);*/
							newThetaAll[n] = max(newThetaAll[n], 0.01);
							break;
						case 4:
							//newThetaAll[n]=max(newThetaAll[n],0.01);
							newThetaAll[n] = max(newThetaAll[n], 0.5);
							newThetaAll[n] = min(newThetaAll[n], sz / 2.0f);
							/*newTheta[4] = max(newTheta[4],0.5);
							newTheta[4] = min(newTheta[4],sz/2.0f);*/

							break;
						}
						for (j = 0; j < noChannels; j++) {
							//newTheta[i+5*j]=newThetaAll[n]+dT[i+j*5];
							newTheta[i + 5 * j] = newThetaAll[n] * dT[noChannels * 5 + i + j * 5] + dT[i + j * 5];
						}
					}

					else
					{
						for (j = 0; j < noChannels; j++) {
							switch (i) {
							case 0:
							case 1:

								/*newThetaAll[n+j]=max(newThetaAll[n+j],(float(sz)-1)/2-sz/4.0);
								newThetaAll[n+j] = min(newThetaAll[n+j],(float(sz)-1)/2+sz/4.0);*/


								break;
							case 2:

								/*newThetaAll[n+j]=max(newThetaAll[n+j],0.0);
								newThetaAll[n+j]=min(newThetaAll[n+j],float(spline_zsize));*/
								newThetaAll[n + j] = max(newThetaAll[n + j], 1.0);

								break;
							case 3:

								/*newThetaAll[n+j]=max(newThetaAll[n+j],1.0);*/
								newThetaAll[n + j] = max(newThetaAll[n + j], 0.01);
								break;
							case 4:
								//newThetaAll[n+j]=max(newThetaAll[n+j],0.01);

								newThetaAll[n + j] = max(newThetaAll[n + j], 0.5);
								newThetaAll[n + j] = min(newThetaAll[n + j], sz / 2.0f);

								break;
							}
							newTheta[i + 5 * j] = newThetaAll[n + j] + dT[i + j * 5];
						}
						n = n + j - 1;
					}
					n = n + 1;
				}


				//prepare for next iteration
				newErr = 0;

				memset(jacobian, 0, NV * sizeof(float));
				memset(hessian, 0, NV * NV * sizeof(float));
				memset(newDudtAll, 0, NV * noChannels * sizeof(float));

				for (ii = 0; ii < sz; ii++) for (jj = 0; jj < sz; jj++) {
					for (i = 0; i < noChannels; i++) {
						kernel_DerivativeGauss2D_sigma(ii, jj, &newTheta[5 * i], &newDudt[5 * i], &model[i]);
						data[i] = *(s_data[i] + sz * jj + ii) + *(s_varim[i] + sz * jj + ii);
						model[i] += *(s_varim[i] + sz * jj + ii);
					}

					n = 0;
					for (i = 0; i < 5; i++) {
						if (cshared[i] == 1) {
							for (j = 0; j < noChannels; j++) {
								newDudtAll[n + j * NV] = newDudt[i + j * 5];
							}
						}
						else
						{
							for (j = 0; j < noChannels; j++) {
								newDudtAll[n + j + j * NV] = newDudt[i + j * 5];
							}
							n = n + j - 1;
						}
						n++;
					}

					for (i = 0; i < noChannels; i++) {
						if (data[i] > 0)
							newErr = newErr + 2 * ((model[i] - data[i]) - data[i] * log(model[i] / data[i]));
						else
						{
							newErr = newErr + 2 * model[i];
							data[i] = 0;
						}
						t1[i] = 1 - data[i] / model[i];
						t2[i] = data[i] / pow(model[i], 2);
					}

					for (l = 0; l < NV; l++) {
						for (i = 0; i < noChannels; i++) {
							jacobian[l] += t1[i] * newDudtAll[l + i * NV];
						}
					}

					for (l = 0; l < NV; l++) for (m = l; m < NV; m++) {
						for (i = 0; i < noChannels; i++) {
							hessian[l * NV + m] += t2[i] * newDudtAll[l + i * NV] * newDudtAll[m + i * NV];
						}
						hessian[m * NV + l] = hessian[l * NV + m];
					}
				}
			}
			else
			{
				mu = max((1 + newLambda * SCALE_UP) / (1 + newLambda), 1.3f);
				newLambda = SCALE_UP * newLambda;
			}
		}
	}

	//output iteration time
	d_Parameters[Nfits * NV + subregion] = kk;
	/*d_Parameters[Nfits*(NV+1)+BlockSize*bx+tx]=maxJump[0];
	d_Parameters[Nfits*(NV+2)+BlockSize*bx+tx]=maxJump[1];
	d_Parameters[Nfits*(NV+3)+BlockSize*bx+tx]=maxJump[2];
	d_Parameters[Nfits*(NV+4)+BlockSize*bx+tx]=maxJump[3];
	d_Parameters[Nfits*(NV+5)+BlockSize*bx+tx]=maxJump[4];
	d_Parameters[Nfits*(NV+6)+BlockSize*bx+tx]=maxJump[5];
	d_Parameters[Nfits*(NV+7)+BlockSize*bx+tx]=maxJump[6];*/
	//d_Parameters[Nfits*(NV+4)+BlockSize*bx+tx]=maxJump[3];


	// Calculating the CRLB and LogLikelihood
	Div = 0.0;

	for (ii = 0; ii < sz; ii++) for (jj = 0; jj < sz; jj++) {
		for (i = 0; i < noChannels; i++) {
			kernel_DerivativeGauss2D_sigma(ii, jj, &newTheta[5 * i], &newDudt[5 * i], &model[i]);
			data[i] = *(s_data[i] + sz * jj + ii) + *(s_varim[i] + sz * jj + ii);
			model[i] += *(s_varim[i] + sz * jj + ii);
		}

		n = 0;
		for (i = 0; i < 5; i++) {
			if (cshared[i] == 1) {
				for (j = 0; j < noChannels; j++) {
					//newDudtAll[n+j*NV]=newDudt[i+j*5];
					newDudtAll[n + j * NV] = newDudt[i + j * 5] * dT[noChannels * 5 + i + j * 5];
				}
			}
			else
			{
				for (j = 0; j < noChannels; j++) {
					newDudtAll[n + j + j * NV] = newDudt[i + j * 5];
				}
				n = n + j - 1;
			}
			n++;
		}

		//Building the Fisher Information Matrix
		for (kk = 0; kk < NV; kk++)for (ll = kk; ll < NV; ll++) {
			for (j = 0; j < noChannels; j++) {
				M[kk * NV + ll] += newDudtAll[ll + j * NV] * newDudtAll[kk + j * NV] / model[j];
			}
			M[ll * NV + kk] = M[kk * NV + ll];
		}

		//LogLikelyhood
		for (i = 0; i < noChannels; i++) {
			if (model[i] > 0)
				if (data[i] > 0)Div += data[i] * log(model[i]) - model[i] - data[i] * log(data[i]) + data[i];
				else
					Div += -model[i];
		}
	}
	// Matrix inverse (CRLB=F^-1) and output assigments
	kernel_MatInvN(M, Minv, Diag, NV);

	/*n = 0;
	for (i=0;i<5;i++){
		if (cshared[i]!=1)
		{
			for (j=0;j<noChannels;j++){
				newThetaAll[n+j]=newThetaAll[n+j]+dT[i+j*5];
			}
			n = n+j-1;
		}
		n=n+1;
	}*/

	//write to global arrays
	for (kk = 0; kk < NV; kk++) d_Parameters[Nfits * kk + subregion] = newThetaAll[kk];
	for (kk = 0; kk < NV; kk++) d_CRLBs[Nfits * kk + subregion] = Diag[kk];
	//d_LogLikelihood[BlockSize*bx+tx] = newUpdate[0];
	d_LogLikelihood[subregion] = Div;
	//d_LogLikelihood[BlockSize*bx+tx] = 1;


	return;


}
//**********************************************************************************************************************************
//multichannel fit
void kernel_splineMLEFit_z_sCMOS_multi(const int subregion, const float* d_data, const float* d_coeff, const float* d_dTAll, const int spline_xsize, const int spline_ysize, const int spline_zsize, const int sz, const int iterations, const int NV, const int noChannels, float* d_Parameters, float* d_CRLBs, float* d_LogLikelihood, float *initZ, const int Nfits, const int* d_shared, const float* d_varim) {
	float M[25 * Max_No_Channel * Max_No_Channel] = { 0 }, Diag[5 * Max_No_Channel] = { 0 }, Minv[25 * Max_No_Channel * Max_No_Channel] = { 0 };
	//int tx = threadIdx.x;
	//int bx = blockIdx.x;
	//int BlockSize = blockDim.x;
	int ii, jj, kk, ll, l, m, i, j, n;
	int xstart[Max_No_Channel] = { 0 }, ystart[Max_No_Channel] = { 0 }, zstart[Max_No_Channel] = { 0 };


	float model[Max_No_Channel] = { 0 }, data[Max_No_Channel];
	float Div = 0;

	float newTheta[5 * Max_No_Channel] = { 0 }, oldTheta[5 * Max_No_Channel] = { 0 }, newThetaAll[5 * Max_No_Channel] = { 0 };
	float newLambda = INIT_LAMBDA, oldLambda = INIT_LAMBDA, mu;
	float newUpdate[5 * Max_No_Channel], oldUpdate[5 * Max_No_Channel];
	float maxJump_Init[5 * Max_No_Channel] = { 0 }, maxJump[5 * Max_No_Channel] = { 0 };
	float newDudt[5 * Max_No_Channel] = { 0 }, newDudtAll[5 * Max_No_Channel * Max_No_Channel];

	float newErr = 1e12, oldErr = 1e13;

	int off;
	float jacobian[5 * Max_No_Channel] = { 0 };
	float hessian[25 * Max_No_Channel * Max_No_Channel] = { 0 };
	float t1[Max_No_Channel] = { 0 }, t2[Max_No_Channel] = { 0 };

	float Nmax;
	float xc[Max_No_Channel] = { 0 }, yc[Max_No_Channel] = { 0 }, zc[Max_No_Channel] = { 0 };
	float delta_f[64 * Max_No_Channel] = { 0 }, delta_dxf[64 * Max_No_Channel] = { 0 }, delta_dyf[64 * Max_No_Channel] = { 0 }, delta_dzf[64 * Max_No_Channel] = { 0 };
	int errFlag = 0;
	float L[25 * Max_No_Channel * Max_No_Channel] = { 0 }, U[25 * Max_No_Channel * Max_No_Channel] = { 0 };

	//Prevent read/write past end of array
	if ((subregion) >= Nfits) return;

	for (i = 0; i < NV; i++) {
		newUpdate[i] = 1e13;
		oldUpdate[i] = 1e13;
	}
	//memset(newUpdate,1e13,NV*NV*sizeof(float));
	//memset(oldUpdate,1e13,NV*NV*sizeof(float));

	for (i = 0; i < noChannels; i++) {
		maxJump_Init[0 + i * 5] = 1.0f;
		maxJump_Init[1 + i * 5] = 1.0f;
		maxJump_Init[2 + i * 5] = spline_zsize / 5.0f;
		maxJump_Init[3 + i * 5] = 100.0f;
		maxJump_Init[4 + i * 5] = 20.0f;
	}



	//copy in data
	const float* s_data[Max_No_Channel];
	const float* s_varim[Max_No_Channel];
	const float* dT;
	const int* cshared;
	for (i = 0; i < noChannels; i++) {
		s_data[i] = d_data + (sz * sz * (Nfits * i + subregion));
	}

	for (i = 0; i < noChannels; i++) {
		s_varim[i] = d_varim + (sz * sz * (Nfits * i + subregion));
	}


	dT = d_dTAll + (subregion) * 5 * noChannels * 2;
	cshared = d_shared + (subregion) * 5;

	//initial values
	for (i = 0; i < noChannels; i++) {
		kernel_CenterofMass2D(sz, s_data[i], &newTheta[0 + 5 * i], &newTheta[1 + 5 * i]);
		kernel_GaussFMaxMin2D(sz, 1.5, s_data[i], &Nmax, &newTheta[4 + 5 * i]);
		newTheta[2 + 5 * i] = initZ[subregion];

		newTheta[3 + 5 * i] = (Nmax - newTheta[3 + 5 * i]) / d_coeff[(int)(spline_zsize / 2) * (spline_xsize * spline_ysize) + (int)(spline_ysize / 2) * spline_xsize + (int)(spline_xsize / 2) + i * spline_xsize * spline_ysize * spline_zsize * 64] * 4;
		newTheta[4 + 5 * i] = max(newTheta[4 + 5 * i], 0.01);

		maxJump_Init[3 + 5 * i] = max(maxJump_Init[3 + 5 * i], newTheta[3 + 5 * i]);
		maxJump_Init[4 + 5 * i] = max(maxJump_Init[4 + 5 * i], newTheta[4 + 5 * i]);
	}

	// map paramters
	for (i = 0; i < 5; i++) {
		if (cshared[i] == 1) {
			for (j = 1; j < noChannels; j++) {
				//newTheta[i+j*5]=newTheta[i]+dT[i+j*5];
				newTheta[i + j * 5] = newTheta[i] * dT[noChannels * 5 + i + j * 5] + dT[i + j * 5];

			}
		}
	}

	for (ii = 0; ii < 5 * noChannels; ii++)oldTheta[ii] = newTheta[ii];

	//combine all channel
	n = 0;
	for (i = 0; i < 5; i++) {
		if (cshared[i] == 1) {
			maxJump[n] = maxJump_Init[i];
			newThetaAll[n] = newTheta[i];
		}
		else
		{
			for (j = 0; j < noChannels; j++) {
				maxJump[n + j] = maxJump_Init[i + j * 5];
				newThetaAll[n + j] = newTheta[i + j * 5];
			}
			n = n + j - 1;
		}
		n = n + 1;
	}



	//prepare values for iteration
	off = (int)((float(spline_xsize) + 1.0 - float(sz)) / 2);

	for (i = 0; i < noChannels; i++) {
		xc[i] = -1.0 * ((newTheta[0 + 5 * i] - float(sz) / 2) + 0.5);
		yc[i] = -1.0 * ((newTheta[1 + 5 * i] - float(sz) / 2) + 0.5);

		xstart[i] = floor(xc[i]);
		xc[i] = xc[i] - xstart[i];

		ystart[i] = floor(yc[i]);
		yc[i] = yc[i] - ystart[i];

		//zstart = floor(newTheta[4]);
		zstart[i] = floor(newTheta[2 + 5 * i]);
		zc[i] = newTheta[2 + 5 * i] - zstart[i];
	}

	newErr = 0;

	memset(jacobian, 0, NV * sizeof(float));
	memset(hessian, 0, NV * NV * sizeof(float));
	memset(newDudtAll, 0, NV * noChannels * sizeof(float));

	for (i = 0; i < noChannels; i++) {
		kernel_computeDelta3D(xc[i], yc[i], zc[i], &delta_f[64 * i], &delta_dxf[64 * i], &delta_dyf[64 * i], &delta_dzf[64 * i]);
	}

	for (ii = 0; ii < sz; ii++) for (jj = 0; jj < sz; jj++) {
		for (i = 0; i < noChannels; i++) {
			kernel_DerivativeSpline_multiChannel(ii + xstart[i] + off, jj + ystart[i] + off, zstart[i], spline_xsize, spline_ysize, spline_zsize, &delta_f[64 * i], &delta_dxf[64 * i], &delta_dyf[64 * i], &delta_dzf[64 * i], &d_coeff[spline_xsize * spline_ysize * spline_zsize * 64 * i], &newTheta[5 * i], &newDudt[5 * i], &model[i]);
			data[i] = *(s_data[i] + sz * jj + ii) + *(s_varim[i] + sz * jj + ii);
			model[i] += *(s_varim[i] + sz * jj + ii);
			//s_data[i]= d_data+(sz*sz*bx*BlockSize+sz*sz*tx+sz*sz*Nfits*i);
		}

		n = 0;
		for (i = 0; i < 5; i++) {
			if (cshared[i] == 1) {
				for (j = 0; j < noChannels; j++) {
					//newDudtAll[n+j*NV]=newDudt[i+j*5];
					newDudtAll[n + j * NV] = newDudt[i + j * 5] * dT[noChannels * 5 + i + j * 5];
				}
			}
			else
			{
				for (j = 0; j < noChannels; j++) {
					newDudtAll[n + j + j * NV] = newDudt[i + j * 5];
				}
				n = n + j - 1;
			}
			n++;
		}

		for (i = 0; i < noChannels; i++) {
			if (data[i] > 0)
				newErr = newErr + 2 * ((model[i] - data[i]) - data[i] * log(model[i] / data[i]));
			else
			{
				newErr = newErr + 2 * model[i];
				data[i] = 0;
			}
			t1[i] = 1 - data[i] / model[i];
			t2[i] = data[i] / pow(model[i], 2);
		}

		for (l = 0; l < NV; l++) {
			for (i = 0; i < noChannels; i++) {
				jacobian[l] += t1[i] * newDudtAll[l + i * NV];
			}
		}

		for (l = 0; l < NV; l++) for (m = l; m < NV; m++) {
			for (i = 0; i < noChannels; i++) {
				hessian[l * NV + m] += t2[i] * newDudtAll[l + i * NV] * newDudtAll[m + i * NV];
			}
			hessian[m * NV + l] = hessian[l * NV + m];
		}
	}

	for (kk = 0; kk < iterations; kk++) {//main iterative loop

		if (fabs((newErr - oldErr) / newErr) < TOLERANCE) {
			//newStatus = CONVERGED;
			break;
		}
		else {
			if (newErr > ACCEPTANCE* oldErr) {
				//copy Fitdata
				for (ii = 0; ii < 5 * noChannels; ii++)newTheta[ii] = oldTheta[ii];
				for (i = 0; i < NV; i++) {
					newUpdate[i] = oldUpdate[i];
				}
				newLambda = oldLambda;
				newErr = oldErr;
				mu = max((1 + newLambda * SCALE_UP) / (1 + newLambda), 1.3f);
				newLambda = SCALE_UP * newLambda;
			}
			else if (newErr < oldErr && errFlag == 0) {
				newLambda = SCALE_DOWN * newLambda;
				mu = 1 + newLambda;
			}


			for (i = 0; i < NV; i++) {
				hessian[i * NV + i] = hessian[i * NV + i] * mu;
			}
			memset(L, 0, NV * NV * sizeof(float));
			memset(U, 0, NV * NV * sizeof(float));

			errFlag = kernel_cholesky(hessian, NV, L, U);
			if (errFlag == 0) {
				for (ii = 0; ii < 5 * noChannels; ii++)oldTheta[ii] = newTheta[ii];
				for (i = 0; i < NV; i++) {
					oldUpdate[i] = newUpdate[i];
				}
				oldLambda = newLambda;
				oldErr = newErr;


				//mexPrintf("before kernel_luEvaluate\n");
				//mexPrintf("before kernel_luEvaluate NV %f\n",(double)NV);
				//for (ii=0;ii<100;ii++){
				//	mexPrintf("before L is %f, i is  %f\n",(double)L[ii],(double)ii);
				//	mexPrintf("before U is %f, i is  %f\n",(double)U[ii],(double)ii);
				//}

				kernel_luEvaluate(L, U, jacobian, NV, newUpdate);
				//return;

				//updateFitParameters
				for (ll = 0; ll < NV; ll++) {
					if (newUpdate[ll] / oldUpdate[ll] < -0.5f) {
						maxJump[ll] = maxJump[ll] * 0.5;

					}
					//test = oldUpdate[1];
					newUpdate[ll] = newUpdate[ll] / (1 + fabs(newUpdate[ll] / maxJump[ll]));
					newThetaAll[ll] = newThetaAll[ll] - newUpdate[ll];
				}

				n = 0;
				for (i = 0; i < 5; i++) {
					if (cshared[i] == 1) {
						switch (i) {
						case 0:
						case 1:
							newThetaAll[n] = max(newThetaAll[n], (float(sz) - 1) / 2 - sz / 4.0);
							newThetaAll[n] = min(newThetaAll[n], (float(sz) - 1) / 2 + sz / 4.0);
							/*								test = 1;*/
							break;
						case 2:
							newThetaAll[n] = max(newThetaAll[n], 0.0);
							newThetaAll[n] = min(newThetaAll[n], float(spline_zsize));
							break;
						case 3:
							newThetaAll[n] = max(newThetaAll[n], 1.0);
							break;
						case 4:
							newThetaAll[n] = max(newThetaAll[n], 0.01);
							break;
						}
						for (j = 0; j < noChannels; j++) {
							//newTheta[i+5*j]=newThetaAll[n]+dT[i+j*5];
							newTheta[i + 5 * j] = newThetaAll[n] * dT[noChannels * 5 + i + j * 5] + dT[i + j * 5];
						}
					}

					else
					{
						for (j = 0; j < noChannels; j++) {
							switch (i) {
							case 0:
							case 1:
								/*newThetaAll[n]=max(newThetaAll[n],(float(sz)-1)/2-sz/4.0);
								newThetaAll[n] = min(newThetaAll[n],(float(sz)-1)/2+sz/4.0);*/

								newThetaAll[n + j] = max(newThetaAll[n + j], (float(sz) - 1) / 2 - sz / 4.0);
								newThetaAll[n + j] = min(newThetaAll[n + j], (float(sz) - 1) / 2 + sz / 4.0);


								break;
							case 2:
								/*newThetaAll[n]=max(newThetaAll[n],0.0);
								newThetaAll[n]=min(newThetaAll[n],float(spline_zsize));*/

								newThetaAll[n + j] = max(newThetaAll[n + j], 0.0);
								newThetaAll[n + j] = min(newThetaAll[n + j], float(spline_zsize));

								break;
							case 3:
								/*newThetaAll[n]=max(newThetaAll[n],1.0);*/

								newThetaAll[n + j] = max(newThetaAll[n + j], 1.0);
								break;
							case 4:
								/*newThetaAll[n]=max(newThetaAll[n],0.01);*/

								newThetaAll[n + j] = max(newThetaAll[n + j], 0.01);
								break;
							}
							newTheta[i + 5 * j] = newThetaAll[n + j] + dT[i + j * 5];
						}
						n = n + j - 1;
					}
					n = n + 1;
				}


				//prepare for next iteration
				for (i = 0; i < noChannels; i++) {
					xc[i] = -1.0 * ((newTheta[0 + 5 * i] - float(sz) / 2) + 0.5);
					yc[i] = -1.0 * ((newTheta[1 + 5 * i] - float(sz) / 2) + 0.5);

					xstart[i] = floor(xc[i]);
					xc[i] = xc[i] - xstart[i];

					ystart[i] = floor(yc[i]);
					yc[i] = yc[i] - ystart[i];

					//zstart = floor(newTheta[4]);
					zstart[i] = floor(newTheta[2 + 5 * i]);
					zc[i] = newTheta[2 + 5 * i] - zstart[i];
				}

				newErr = 0;

				memset(jacobian, 0, NV * sizeof(float));
				memset(hessian, 0, NV * NV * sizeof(float));
				memset(newDudtAll, 0, NV * noChannels * sizeof(float));

				for (i = 0; i < noChannels; i++) {
					kernel_computeDelta3D(xc[i], yc[i], zc[i], &delta_f[64 * i], &delta_dxf[64 * i], &delta_dyf[64 * i], &delta_dzf[64 * i]);
				}


				for (ii = 0; ii < sz; ii++) for (jj = 0; jj < sz; jj++) {
					for (i = 0; i < noChannels; i++) {
						kernel_DerivativeSpline_multiChannel(ii + xstart[i] + off, jj + ystart[i] + off, zstart[i], spline_xsize, spline_ysize, spline_zsize, &delta_f[64 * i], &delta_dxf[64 * i], &delta_dyf[64 * i], &delta_dzf[64 * i], &d_coeff[spline_xsize * spline_ysize * spline_zsize * 64 * i], &newTheta[5 * i], &newDudt[5 * i], &model[i]);
						data[i] = *(s_data[i] + sz * jj + ii) + *(s_varim[i] + sz * jj + ii);
						model[i] += *(s_varim[i] + sz * jj + ii);
						//s_data[i]= d_data+(sz*sz*bx*BlockSize+sz*sz*tx+sz*sz*Nfits*i);
					}

					n = 0;
					for (i = 0; i < 5; i++) {
						if (cshared[i] == 1) {
							for (j = 0; j < noChannels; j++) {
								//newDudtAll[n+j*NV]=newDudt[i+j*5];
								newDudtAll[n + j * NV] = newDudt[i + j * 5] * dT[noChannels * 5 + i + j * 5];
							}
						}
						else
						{
							for (j = 0; j < noChannels; j++) {
								newDudtAll[n + j + j * NV] = newDudt[i + j * 5];
							}
							n = n + j - 1;
						}
						n++;
					}

					for (i = 0; i < noChannels; i++) {
						if (data[i] > 0)
							newErr = newErr + 2 * ((model[i] - data[i]) - data[i] * log(model[i] / data[i]));
						else
						{
							newErr = newErr + 2 * model[i];
							data[i] = 0;
						}
						t1[i] = 1 - data[i] / model[i];
						t2[i] = data[i] / pow(model[i], 2);
					}

					for (l = 0; l < NV; l++) {
						for (i = 0; i < noChannels; i++) {
							jacobian[l] += t1[i] * newDudtAll[l + i * NV];
						}
					}

					for (l = 0; l < NV; l++) for (m = l; m < NV; m++) {
						for (i = 0; i < noChannels; i++) {
							hessian[l * NV + m] += t2[i] * newDudtAll[l + i * NV] * newDudtAll[m + i * NV];
						}
						hessian[m * NV + l] = hessian[l * NV + m];
					}
				}
			}
			else
			{
				mu = max((1 + newLambda * SCALE_UP) / (1 + newLambda), 1.3f);
				newLambda = SCALE_UP * newLambda;
			}
		}
	}

	//output iteration time
	d_Parameters[Nfits * NV + subregion] = kk;
	/*d_Parameters[Nfits*(NV+1)+BlockSize*bx+tx]=maxJump[0];
	d_Parameters[Nfits*(NV+2)+BlockSize*bx+tx]=maxJump[1];
	d_Parameters[Nfits*(NV+3)+BlockSize*bx+tx]=maxJump[2];
	d_Parameters[Nfits*(NV+4)+BlockSize*bx+tx]=maxJump[3];
	d_Parameters[Nfits*(NV+5)+BlockSize*bx+tx]=maxJump[4];
	d_Parameters[Nfits*(NV+6)+BlockSize*bx+tx]=maxJump[5];
	d_Parameters[Nfits*(NV+7)+BlockSize*bx+tx]=maxJump[6];*/
	//d_Parameters[Nfits*(NV+4)+BlockSize*bx+tx]=maxJump[3];


	// Calculating the CRLB and LogLikelihood
	Div = 0.0;
	for (i = 0; i < noChannels; i++) {
		xc[i] = -1.0 * ((newTheta[0 + 5 * i] - float(sz) / 2) + 0.5);
		yc[i] = -1.0 * ((newTheta[1 + 5 * i] - float(sz) / 2) + 0.5);

		xstart[i] = floor(xc[i]);
		xc[i] = xc[i] - xstart[i];

		ystart[i] = floor(yc[i]);
		yc[i] = yc[i] - ystart[i];

		//zstart = floor(newTheta[4]);
		zstart[i] = floor(newTheta[2 + 5 * i]);
		zc[i] = newTheta[2 + 5 * i] - zstart[i];
	}


	for (i = 0; i < noChannels; i++) {
		kernel_computeDelta3D(xc[i], yc[i], zc[i], &delta_f[64 * i], &delta_dxf[64 * i], &delta_dyf[64 * i], &delta_dzf[64 * i]);
	}

	for (ii = 0; ii < sz; ii++) for (jj = 0; jj < sz; jj++) {
		for (i = 0; i < noChannels; i++) {
			kernel_DerivativeSpline_multiChannel(ii + xstart[i] + off, jj + ystart[i] + off, zstart[i], spline_xsize, spline_ysize, spline_zsize, &delta_f[64 * i], &delta_dxf[64 * i], &delta_dyf[64 * i], &delta_dzf[64 * i], &d_coeff[spline_xsize * spline_ysize * spline_zsize * 64 * i], &newTheta[5 * i], &newDudt[5 * i], &model[i]);
			data[i] = *(s_data[i] + sz * jj + ii) + *(s_varim[i] + sz * jj + ii);
			model[i] += *(s_varim[i] + sz * jj + ii);
			//s_data[i]= d_data+(sz*sz*bx*BlockSize+sz*sz*tx+sz*sz*Nfits*i);
		}

		n = 0;
		for (i = 0; i < 5; i++) {
			if (cshared[i] == 1) {
				for (j = 0; j < noChannels; j++) {
					//newDudtAll[n+j*NV]=newDudt[i+j*5];
					newDudtAll[n + j * NV] = newDudt[i + j * 5] * dT[noChannels * 5 + i + j * 5];

				}
			}
			else
			{
				for (j = 0; j < noChannels; j++) {
					newDudtAll[n + j + j * NV] = newDudt[i + j * 5];
				}
				n = n + j - 1;
			}
			n++;
		}

		//Building the Fisher Information Matrix
		for (kk = 0; kk < NV; kk++)for (ll = kk; ll < NV; ll++) {
			for (j = 0; j < noChannels; j++) {
				M[kk * NV + ll] += newDudtAll[ll + j * NV] * newDudtAll[kk + j * NV] / model[j];
			}
			M[ll * NV + kk] = M[kk * NV + ll];
		}

		//LogLikelyhood
		for (i = 0; i < noChannels; i++) {
			if (model[i] > 0)
				if (data[i] > 0)Div += data[i] * log(model[i]) - model[i] - data[i] * log(data[i]) + data[i];
				else
					Div += -model[i];
		}
	}
	// Matrix inverse (CRLB=F^-1) and output assigments
	kernel_MatInvN(M, Minv, Diag, NV);

	//write to global arrays
	for (kk = 0; kk < NV; kk++) d_Parameters[Nfits * kk + subregion] = newThetaAll[kk];
	for (kk = 0; kk < NV; kk++) d_CRLBs[Nfits * kk + subregion] = Diag[kk];
	//d_LogLikelihood[BlockSize*bx+tx] = newUpdate[0];
	d_LogLikelihood[subregion] = Div;
	//d_LogLikelihood[BlockSize*bx+tx] = 1;


	return;
}







