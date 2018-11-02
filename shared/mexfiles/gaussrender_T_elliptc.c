#include "mex.h"
#include <stdio.h>
#include <math.h>

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
double gaussrender(float *srim,float *xpix, float *ypix, mwSize *srec, float *sigmax, float *sigmay, float *Gtemplate, float Gsigma, float roiks,  float *N, int uselut, float *c, float *lut, float *rangec, mwSize Gx,mwSize numlocs, mwSize sl, float transparency)
{
mwSize Gsizegauss,indc,xt,yt,col,srimindlin;
float dx,dy,intcorrectionx,intcorrectiony,gaussnorm,inttemp,oldfac,oldfac2;
long k,dnx,dny,xr,yr,xax,yax,xp,yp;
double numberOfLocs;

numberOfLocs=0;
Gsizegauss=(Gx-1)/2; /**/

for(k=0;k<numlocs;k++)
    {
    dnx=roiks*sigmax[k]+1;
    dny=roiks*sigmay[k]+1;
    
    xr=xpix[k]+0.5;
    yr=ypix[k]+0.5;
    if(xr>=-dnx&&xr<srec[0]+dnx&&yr>=-dny&&yr<srec[1]+dny)
        {
        if(xr>=0&&xr<srec[0]&&yr>=0&&yr<srec[1])
            {    
            numberOfLocs++;
        }
        dx=xpix[k]-xr;
        dy=ypix[k]-yr;
        intcorrectionx=(erf((dnx+0.5)/sigmax[k]/1.4142135624));
        intcorrectiony=(erf((dny+0.5)/sigmay[k]/1.4142135624));
        //intcorrection=1;
        gaussnorm=N[k]/(2*3.1415926536*sigmax[k]*sigmay[k]*intcorrectionx*intcorrectiony);
        if(uselut==1)
            {
            indc=(c[k]-rangec[0])/(rangec[1]-rangec[0])*(sl-1);
        }
        for(xax=-dnx;xax<=dnx;xax++)
            {
            xt=(xax-dx)*Gsigma/sigmax[k]+Gsizegauss+0.5; /* careful: indexing moves from 1 to zero. How to take this into account??? minus 1. seems to be right. wrong in matlab?*/
            for(yax=-dny;yax<=dny;yax++)
                {
                yt=(yax-dy)*Gsigma/sigmay[k]+Gsizegauss+0.5;
                xp=xr+xax; 
                yp=yr+yax;
                if(xp>=0&&xp<srec[0]&&yp>=0&&yp<srec[1]&&xt>=0&&xt<Gx&&yt>=0&&yt<Gx)
                    {
                    if(uselut==1)
                        {
                        if(indc>=0&indc<sl)
                            {
                            inttemp=(Gtemplate[xt+yt*Gx]*gaussnorm);
                            oldfac=transparency*inttemp;
                            oldfac2=oldfac*5;
                            if(oldfac2>1) oldfac2=1;
                            
                            for(col=0;col<3;col++)
                                {
                                srimindlin=col*srec[1]*srec[0]+yp*srec[0]+xp;
                                
                                /* srim[srimindlin]+=Gtemplate[xt+yt*Gx]*gaussnorm*lut[indc+col*sl];*/

                                
                                srim[srimindlin]=oldfac*lut[indc+col*sl]+(1-oldfac2)*srim[srimindlin];
                                }
                            }
                        }
                     else

                        {
                        inttemp=(Gtemplate[xt+yt*Gx]*gaussnorm);
                        //oldfac=(1-transparency*inttemp);
                        //if(oldfac<0) oldfac=0;
                        //srim[xp+yp*srec[0]]=transparency*inttemp*inttemp+oldfac*srim[xp+yp*srec[0]];
                        
                        oldfac=transparency*inttemp;
                        oldfac2=oldfac;
                        if(oldfac2>1) oldfac2=1;
                        srimindlin=yp*srec[0]+xp;
                        srim[srimindlin]=oldfac+(1-oldfac2)*srim[srimindlin];
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
  float *Gtemplate,*xpix,*ypix,*sigmax,*sigmay,*N,*c;
  float *lut,*rangec;
  float Gsigma,roiks,transparency;
  int uselut;
  double numberOfLocs;
  mwSize *srec,srec3[]={10,10,3};

    mwSize Gx,Gy,numlocs,sl,sz;
  /* xpix, ypix, srec, sigma, 4.Gtemplate, 5. Gsigma, roiks, 7. N, uselut, 9. c, 10. lut, 11. rangec */

  /*  create a pointer to the input matrix y */
  xpix = (float*) mxGetData(prhs[0]);
 ypix = (float*)mxGetData(prhs[1]);
 srec = (mwSize*)mxGetData(prhs[2]);
 sigmax = (float*)mxGetData(prhs[3]);
 sigmay = (float*)mxGetData(prhs[4]);
 Gtemplate = (float*)mxGetData(prhs[5]);
 N = (float*)mxGetData(prhs[8]);
 c = (float*)mxGetData(prhs[10]);
 lut = (float*)mxGetData(prhs[11]);
 rangec = (float*)mxGetData(prhs[12]);
  
   
  /*  get the dimensions of the matrix input y */
  Gx = mxGetM(prhs[5]);

  numlocs=mxGetM(prhs[0]);
  sl=mxGetM(prhs[11]);
  
  Gsigma=mxGetScalar(prhs[6]);
  roiks=mxGetScalar(prhs[7]);
  uselut=mxGetScalar(prhs[9]);
 transparency = mxGetScalar(prhs[13]);

  /*  set the output pointer to the output matrix */
 if(uselut==0)
     {
      srec3[2]=1;
      }

     srec3[0]=srec[0];
     srec3[1]=srec[1];
         //  printf("output size %i,%i,%i\n",srec3[0],srec3[1],srec3[2]);
     plhs[0] = mxCreateNumericArray(3,srec3,mxSINGLE_CLASS,mxREAL);


  
  /*  create a C pointer to a copy of the output matrix */
  srim = mxGetData(plhs[0]);
  
  /*  call the C subroutine */
  numberOfLocs=gaussrender(srim,xpix, ypix, srec, sigmax,sigmay, Gtemplate, Gsigma, roiks,  N, uselut, c, lut, rangec,Gx,numlocs,sl,transparency);
  plhs[1]=mxCreateDoubleScalar(numberOfLocs);
  
}
