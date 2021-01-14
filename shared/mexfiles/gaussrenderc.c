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
double gaussrender(float *srim,float *xpix, float *ypix,  unsigned int *srec, float *sigma, float *Gtemplate, float Gsigma, float roiks,  float *N, int uselut, float *c, float *lut, float *rangec,  unsigned int Gx, unsigned int numlocs,  unsigned int sl)
{
 unsigned int Gsizegauss,indc,xt,yt,col,srimindlin;
float dx,dy,intcorrection,gaussnorm;
long k,dn,xr,yr,xax,yax,xp,yp;
double numberOfLocs;

Gsizegauss=(Gx-1)/2; /**/
numberOfLocs=0;
for(k=0;k<numlocs;k++)
    {
    dn=roiks*sigma[k]+1;
    xr=xpix[k]+0.5;
    yr=ypix[k]+0.5;
    if(xr>=-dn&&xr<srec[0]+dn&&yr>=-dn&&yr<srec[1]+dn)
        {
    if(xr>=0&&xr<srec[0]&&yr>=0&&yr<srec[1])
        {    
    numberOfLocs++;
    }
    
    dx=xpix[k]-xr;
    dy=ypix[k]-yr;
    intcorrection=(erf((dn+0.5)/sigma[k]/1.4142135624));
    //intcorrection=1;
    gaussnorm=N[k]/(2*3.1415926536*sigma[k]*sigma[k]*intcorrection*intcorrection);
    if(uselut==1)
        {
        indc=(c[k]-rangec[0])/(rangec[1]-rangec[0])*(sl-1);
    }

        
//     }
    for(xax=-dn;xax<=dn;xax++)
        {
        xt=(xax-dx)*Gsigma/sigma[k]+Gsizegauss+0.5; /* careful: indexing moves from 1 to zero. How to take this into account??? minus 1. seems to be right. wrong in matlab?*/
        for(yax=-dn;yax<=dn;yax++)
            {
            yt=(yax-dy)*Gsigma/sigma[k]+Gsizegauss+0.5;
            xp=xr+xax; 
            yp=yr+yax;
            if(xp>=0&&xp<srec[0]&&yp>=0&&yp<srec[1]&&xt>=0&&xt<Gx&&yt>=0&&yt<Gx)
                {
                if(uselut==1)
                    {
                    if(indc>=0&indc<sl)
                        {
                        
                        for(col=0;col<3;col++)
                            {
                            
                            srimindlin=col*srec[1]*srec[0]+yp*srec[0]+xp;
                            srim[srimindlin]+=Gtemplate[xt+yt*Gx]*gaussnorm*lut[indc+col*sl];
                            }
                        }
                    }
                 else
                        
                    {
                    srim[xp+yp*srec[0]]+=Gtemplate[xt+yt*Gx]*gaussnorm;
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
  float *Gtemplate,*xpix,*ypix,*sigma,*N,*c;
  float *lut,*rangec;
  float Gsigma,roiks;
  int uselut;
  double numberOfLocs;
  unsigned int *srec;
  mwSize srec3[]={10,10,3};

     unsigned int Gx,Gy,numlocs,sl,sz;
  /* xpix, ypix, srec, sigma, 4.Gtemplate, 5. Gsigma, roiks, 7. N, uselut, 9. c, 10. lut, 11. rangec */

  /*  create a pointer to the input matrix y */
  xpix = (float*) mxGetData(prhs[0]);
 ypix = (float*)mxGetData(prhs[1]);
 srec = ( unsigned int*)mxGetData(prhs[2]);
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

     srec3[0]=(mwSize)srec[0];
     srec3[1]=(mwSize)srec[1];
     plhs[0] = mxCreateNumericArray(3,srec3,mxSINGLE_CLASS,mxREAL);
     /*printf("output size %i,%i,%i\n",srec3[0],srec3[1],sl);*/
//      plhs[1]=mxCreateDoubleScalar();
  
  /*  create a C pointer to a copy of the output matrix */
  srim = mxGetData(plhs[0]);
//  numberOfLocs=(double*)mxGetData(plhs[1]);
  
  /*  call the C subroutine */
  numberOfLocs=gaussrender(srim,xpix, ypix, srec, sigma, Gtemplate, Gsigma, roiks,  N, uselut, c, lut, rangec,Gx,numlocs,sl);
   plhs[1]=mxCreateDoubleScalar(numberOfLocs);
  
}
