#include "mex.h"
#include <stdio.h>


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
void   findn(double *nb,double *x,double *y,double *x2,double *y2,double dx,long lenx,long lenx2)
{
double dx2,rx2,ry2;
long indx2,this,testind2;

dx2=dx*dx;
indx2=0;
for(this=0;this<lenx;this++) 
{
    if (this%10000==0)
        {
        mexPrintf(" %i of %i k\n",this/1000,lenx/1000);
        mexEvalString("drawnow");
        }    
    while((x2[indx2]<x[this]-dx)&&(indx2<lenx2-1))
    {
        indx2++;
    }
    testind2=indx2;
    
//     while((x2[testind2]>x[this]-dx)&&(testind2>0))
//     {
//         testind2--;
//     }
    while((x2[testind2]<x[this]+dx)&&(testind2<lenx2-1))
    {
        
        ry2=(y2[testind2]-y[this])*(y2[testind2]-y[this]);
        if (ry2<dx2)
            {
                rx2=(x2[testind2]-x[this])*(x2[testind2]-x[this]); 
                if ((rx2+ry2)<dx2)
                nb[this]++;
                
            }
        testind2++;        
    }    
}
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
double *x,*y,*x2,*y2,*nb,dx;
long  lenx,lenx2;

  
  /*  check for proper number of arguments */
  /* NOTE: You do not need an else statement when using mexErrMsgTxt
     within an if statement, because it will never get to the else
     statement if mexErrMsgTxt is executed. (mexErrMsgTxt breaks you out of
     the MEX-file) */
  if(nrhs!=5) 
    mexErrMsgTxt("5 inputs required.");
  if(nlhs!=1) 
    mexErrMsgTxt("1 output required.");
  

  /*  get the scalar input x */
 
  /*  create a pointer to the input matrix y */
 x = mxGetPr(prhs[0]);
 y = mxGetPr(prhs[1]);
 x2 = mxGetPr(prhs[2]);
 y2 = mxGetPr(prhs[3]);
  
  /*  get the dimensions of the matrix input y */
  lenx = mxGetM(prhs[0]);
lenx2 = mxGetM(prhs[2]);
 // printf("length of positons: %i\n",lenx);
  
dx=mxGetScalar(prhs[4]);


  /*  set the output pointer to the output matrix */
//  printf("inputs associated\n");
 
  plhs[0] = mxCreateNumericMatrix(lenx,1,mxDOUBLE_CLASS,mxREAL);


//     printf("outputmatrix created\n");
  
  
  /*  create a C pointer to a copy of the output matrix */
  nb = mxGetData(plhs[0]);

//    printf("outputmatrix created 2, call routine\n");
  /*  call the C subroutine */
  findn(nb,x,y,x2,y2,dx,lenx,lenx2);

}
