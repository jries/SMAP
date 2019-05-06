#include "mex.h"
#include <stdio.h>
#include <math.h>
#include "matrix.h"

double erf(double x)
{
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;
    double t,y;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x);

    // A&S formula 7.1.26
     t = 1.0/(1.0 + p*x);
     y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return sign*y;
}


/*void correlate(double *n1, double *G, mwSize lenG, mwSize lenn)*/
double gaussrender(float *srim,float *xpix, float *ypix,float *zpix, mwSize *srec, float *sigmax, float *sigmay, float roiks,  float *N,mwSize numlocs)
{
mwSize srimindlin;
float dx,dy,dz,intcorrectionx,intcorrectiony,gaussnorm,inttemp,oldfac,sx2,sy2;
long k,dnx,dny,dnz,xr,yr,zr,xax,yax,zax,xp,yp,zp;
double numberOfLocs;

numberOfLocs=0;


for(k=0;k<numlocs;k++)
    {
    dnx=roiks*sigmax[k]+1;
    dny=roiks*sigmax[k]+1;
    dnz=roiks*sigmay[k]+1;
    
    xr=xpix[k]+0.5;
    yr=ypix[k]+0.5;
    zr=zpix[k]+0.5;
    if(xr>=-dnx&&xr<srec[0]+dnx&&yr>=-dny&&yr<srec[1]+dny)
        {
        if(xr>=0&&xr<srec[0]&&yr>=0&&yr<srec[1]&&zr>=0&&zr<srec[2])
            {    
            numberOfLocs++;
        }
        sx2=sigmax[k]*sigmax[k]*2;
        sy2=sigmay[k]*sigmay[k]*2;
                
        dx=xpix[k]-xr;
        dy=ypix[k]-yr;
        dz=zpix[k]-zr;
        intcorrectionx=(erf((dnx+0.5)/sigmax[k]/1.4142135624));
        intcorrectiony=(erf((dny+0.5)/sigmay[k]/1.4142135624));
        //intcorrection=1;
        gaussnorm=N[k]/(2*3.1415926536*sigmax[k]*sigmax[k]*sigmay[k]*intcorrectionx*intcorrectionx*intcorrectiony);
        for(xax=-dnx;xax<=dnx;xax++)
            {
            xp=xr+xax;
//             xt=(xax-dx)*Gsigma/sigmax[k]+Gsizegauss+0.5; /* careful: indexing moves from 1 to zero. How to take this into account??? minus 1. seems to be right. wrong in matlab?*/
            for(yax=-dny;yax<=dny;yax++)
                {
                yp=yr+yax;
//                 yt=(yax-dy)*Gsigma/sigmax[k]+Gsizegauss+0.5;
                for (zax=-dnz;zax<=dnz;zax++)
                    {
                    zp=zr+zax;
//                     zt=(zax-dz)*Gsigma/sigmay[k]+Gsizegauss+0.5;
                    if(xp>=0 && xp<srec[0] && yp>=0 && yp<srec[1]&& zp>=0 && zp<srec[2])
                        {
//                             inttemp=(Gtemplate[xt+yt*Gx+zt*Gx*Gx]*gaussnorm);
                            inttemp=gaussnorm*exp(-(xax-dx)*(xax-dx)/sx2-(yax-dy)*(yax-dy)/sx2-(zax-dz)*(zax-dz)/sy2);
//                             srim[xp][yp][zp]=inttemp;
                            srimindlin=zp*srec[1]*srec[0]+yp*srec[0]+xp;
//                             srimindlin=1;
                            srim[srimindlin]=srim[srimindlin]+inttemp;

                        }
                    }
                }
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
  float *Gtemplate,*xpix,*ypix,*sigmax,*sigmay,*N,*zpix;
  float Gsigma,roiks;
  double numberOfLocs;
  unsigned long *srec;

    mwSize Gx,Gy,numlocs,sl,sz;
/* xpix, ypix,zpix, srec, sigmax,sigmay, 6.Gtemplate, 7. Gsigma, 8.roiks, 9. N */


  /*  create a pointer to the input matrix y */
  xpix = (float*) mxGetData(prhs[0]);
 ypix = (float*)mxGetData(prhs[1]);
  zpix = (float*)mxGetData(prhs[2]);
 srec = (unsigned long*)mxGetData(prhs[3]);
 sigmax = (float*)mxGetData(prhs[4]);
 sigmay = (float*)mxGetData(prhs[5]);
 N = (float*)mxGetData(prhs[7]);
   
  /*  get the dimensions of the matrix input y */
  numlocs=mxGetM(prhs[0]);
// printf("output size %i,%i,%i\n",srec[2],srec[1],srec[0]);
  

  roiks=mxGetScalar(prhs[6]);
  plhs[0] = mxCreateNumericArray(3,srec,mxSINGLE_CLASS,mxREAL);

  /*  create a C pointer to a copy of the output matrix */
  srim = mxGetData(plhs[0]);
  
  /*  call the C subroutine */
  numberOfLocs=gaussrender(srim,xpix, ypix,zpix, srec, sigmax,sigmay, roiks,  N,numlocs);
//   numberOfLocs=1;
  plhs[1]=mxCreateDoubleScalar(numberOfLocs);
  
}
