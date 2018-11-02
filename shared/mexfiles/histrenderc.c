#include "mex.h"
#include <stdio.h>
//#include <math.h>




/*void correlate(double *n1, double *G, mwSize lenG, mwSize lenn)*/
double gaussrender(float *srim,float *xpix, float *ypix, mwSize *srec, float *sigma, float *Gtemplate, float Gsigma, float roiks,  float *N, int uselut, float *c, float *lut, float *rangec, mwSize Gx,mwSize numlocs, mwSize sl)
{
mwSize Gsizegauss,indc,xt,yt,col,srimindlin;
float dx,dy,intcorrection,gaussnorm;
long k,dn,xr,yr,xax,yax,xp,yp;
double numberOfLocs;

Gsizegauss=(Gx-1)/2; /**/
for(k=0;k<numlocs;k++)
    {

    xr=xpix[k]+0.5;
    yr=ypix[k]+0.5;


    gaussnorm=N[k];
    if(uselut==1)
        
    indc=(c[k]-rangec[0])/(rangec[1]-rangec[0])*(sl-1);
    
    
    if(xr>=0&&xr<srec[0]&&yr>=0&&yr<srec[1])
        {
        numberOfLocs++;
        if(uselut==1)
            {
            if(indc>=0&indc<sl)
                {
                for(col=0;col<3;col++)
                    {
                    srimindlin=col*srec[1]*srec[0]+yr*srec[0]+xr;
                    srim[srimindlin]+=gaussnorm*lut[indc+col*sl];
                    }
                }
            }
         else

            {
            srim[xr+yr*srec[0]]+=gaussnorm;
            }

        }
    }
return numberOfLocs;
}           

                    


/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
   float *srim;
  float *Gtemplate,*xpix,*ypix,*sigma,*N,*c;
  float *lut,*rangec;
  float Gsigma,roiks;
  int uselut;
  mwSize *srec,srec3[]={10,10,3};
 double numberOfLocs;
    mwSize Gx,Gy,numlocs,sl,sz;
  /* xpix, ypix, srec, sigma, 4.Gtemplate, 5. Gsigma, roiks, 7. N, uselut, 9. c, 10. lut, 11. rangec */

  /*  create a pointer to the input matrix y */
  xpix = (float*) mxGetData(prhs[0]);
 ypix = (float*)mxGetData(prhs[1]);
 srec = (mwSize*)mxGetData(prhs[2]);
 sigma = (float*)mxGetData(prhs[3]);
 Gtemplate = (float*)mxGetData(prhs[4]);
 N = (float*)mxGetData(prhs[7]);
 c = (float*)mxGetData(prhs[9]);
 lut = (float*)mxGetData(prhs[10]);
 rangec = (float*)mxGetData(prhs[11]);
  

  
  /*  get the dimensions of the matrix input y */
  Gx = mxGetM(prhs[4]);

  numlocs=mxGetM(prhs[0]);
  sl=mxGetM(prhs[10]);
  
  Gsigma=mxGetScalar(prhs[5]);
  roiks=mxGetScalar(prhs[6]);
  uselut=mxGetScalar(prhs[8]);

  /*  set the output pointer to the output matrix */
 if(uselut==0)
     {
      srec3[2]=1;
      }

     srec3[0]=srec[0];
     srec3[1]=srec[1];
     plhs[0] = mxCreateNumericArray(3,srec3,mxSINGLE_CLASS,mxREAL);
     /*printf("output size %i,%i,%i\n",srec3[0],srec3[1],sl);*/

  
  /*  create a C pointer to a copy of the output matrix */
  srim = mxGetData(plhs[0]);
  
  /*  call the C subroutine */
  numberOfLocs=gaussrender(srim,xpix, ypix, srec, sigma, Gtemplate, Gsigma, roiks,  N, uselut, c, lut, rangec,Gx,numlocs,sl);
  plhs[1]=mxCreateDoubleScalar(numberOfLocs);
}
