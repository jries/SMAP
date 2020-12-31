#include "mex.h"
#include <stdio.h>
#include <math.h>

/*
 * xtimesy.c - example found in API guide
 *
 * multiplies an input scalar times an input matrix and outputs a
 * matrix 
 *
 * This is a MEX-file for MATLAB.
 * Copyright 1984-2006 The MathWorks, Inc.
 */

/* $Revision: 1.10.6.2 $ */

/*void correlate(double *n1, double *G, mwSize lenG, mwSize lenn)*/
  //  cs(d list,d x,d y,l frames,d dX,l dT,mws maxactive,mws lenx);
void   findn(double *nb,double *x,double *y,double *z,double *x2,double *y2,double *z2, double dx, double dz, long lenx, long lenx2)
{
double dx2,dtest,dtest2,r2,rz2, dtestz,dtestz2,dz2,ry2,rx2;
long indx2,this,testind2;

dtest=dx*2.5;
dtestz=dz*2.5;
dx2=dx*dx;
dz2=dz*dz;
dtest2=dtest*dtest;
dtestz2=dtestz*dtestz;
indx2=0;
for(this=0;this<lenx;this++) 
{
    if (this%100000==0)
        {
        mexPrintf(" %i of %i k\n",this/1000,lenx/1000);
        mexEvalString("drawnow");
        }
    while((x2[indx2]<x[this]-dtest)&&(indx2<lenx2-1))
    {
        indx2++;
    }
    testind2=indx2;
    while((x2[testind2]<x[this]+dtest)&&(testind2<lenx2-1))
    {
        ry2=(y2[testind2]-y[this])*(y2[testind2]-y[this]);
        if (ry2<dtest2)
            {
            rz2=(z2[testind2]-z[this])*(z2[testind2]-z[this]);
            if (rz2<dtestz2)
                {
                rx2=(x2[testind2]-x[this])*(x2[testind2]-x[this]); 
                nb[this]+=exp(-(rx2+ry2)/(2*dx2)-rz2/(2*dz2));
                }
            }
        testind2++;
        
    }
    
}
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
double *x,*y,*z,*x2,*y2,*z2,*nb,dx,dz;
long  lenx,lenx2;


  
  /*  check for proper number of arguments */
  /* NOTE: You do not need an else statement when using mexErrMsgTxt
     within an if statement, because it will never get to the else
     statement if mexErrMsgTxt is executed. (mexErrMsgTxt breaks you out of
     the MEX-file) */
  if(nrhs!=8) 
    mexErrMsgTxt("8 inputs required.");
  if(nlhs!=1) 
    mexErrMsgTxt("1 output required.");
  

  /*  get the scalar input x */
 
  /*  create a pointer to the input matrix y */
 x = mxGetPr(prhs[0]);
 y = mxGetPr(prhs[1]);
 z = mxGetPr(prhs[2]);
  
  x2 = mxGetPr(prhs[3]);
 y2 = mxGetPr(prhs[4]);
 z2 = mxGetPr(prhs[5]);
  
  /*  get the dimensions of the matrix input y */
  lenx = mxGetM(prhs[0]);
 lenx2 = mxGetM(prhs[3]);
 // printf("length of positons: %i\n",lenx);
  
dx=mxGetScalar(prhs[6]);
dz=mxGetScalar(prhs[7]);

  /*  set the output pointer to the output matrix */
//  printf("inputs associated\n");
 
  plhs[0] = mxCreateNumericMatrix(lenx,1,mxDOUBLE_CLASS,mxREAL);


//     printf("outputmatrix created\n");
  
  
  /*  create a C pointer to a copy of the output matrix */
  nb = mxGetData(plhs[0]);

//    printf("outputmatrix created 2, call routine\n");
  /*  call the C subroutine */
  findn(nb,x,y,z,x2,y2,z2,dx,dz,lenx,lenx2);

}
