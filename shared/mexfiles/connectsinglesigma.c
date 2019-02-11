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
void  cs(double *list, double *x, double *y, double *frames,double dX,long dT,double *dxsigma,mwSize lenx)
{
    long thisentry,stopnow,particlenumber,numpart;
    long particlefound,testentry,numdark; 
    double xh,yh,frh,frtest,lph,dxsigmac; //coordinates of main particle
    
    
    thisentry=0;stopnow=0;
    particlenumber=0;
    numpart=0;
  
    
    while(thisentry<=lenx-1)
    {       
        
        
    while((list[thisentry]>0)&& (thisentry<=lenx-1)) //find next particle which is not connected: list(of this)=0
    {
        thisentry+=1;
    }  
    
    //store coordinates
    particlenumber++;
    xh=x[thisentry];
    yh=y[thisentry];
    lph=dxsigma[thisentry];
    frh=frames[thisentry];
    list[thisentry]=particlenumber;
    
    numdark=0;
    
    testentry=thisentry;
    while ( (numdark<=dT) && (testentry <=lenx-1) && (frames[testentry]>=frh) )//search for dT frames in the future
    {
        particlefound=0;
        
        
        // find index of next frame
        while((frames[testentry]==frh) && (testentry<=lenx-1)) 
        {
            testentry+=1;
        } 
    
        
//         frtest=frames[testentry];
        frtest=frh+1;
        
        //find first entry with fitting x in frame
        while((x[testentry]<xh-dX) && (frames[testentry]==frtest) && (testentry<=lenx-1)) 
        {
            testentry+=1;
        } 

        // compare y for all possible x
        while((x[testentry]<xh+dX) && (frames[testentry]==frtest) && (testentry<=lenx-1))
        {
            dxsigmac=dxsigma[thisentry]*dxsigma[thisentry]+lph*lph;
            if (dxsigmac>dX*dX)
            {
                dxsigmac=dX*dX;
                }
            if ( (y[testentry]-yh)*(y[testentry]-yh)+(x[testentry]-xh)*(x[testentry]-xh)< dxsigmac && (list[testentry]==0))
//             if ((y[testentry]>yh-dX) && (y[testentry]<yh+dX) && (list[testentry]==0)) 
            {
                particlefound=1;
                list[testentry]=particlenumber;
                
                numpart+=1;
                xh=(x[testentry]+xh)/2;
                yh=(y[testentry]+yh)/2;
                frh=frames[testentry];
                numdark=0;
                break;
            }
            else
            {
                testentry+=1;
            }
        }
        if (particlefound==0) 
            {
//             frh=frames[testentry];  //this already went to the next frame...
            frh=frtest; //current frame
            numdark++; //no particle found: add up numdark
            }

    }
    }
    
}
    


/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  double *list,*x,*y,dX,*frames,*appearencelist, *dxsigma; 
  long int  dT;
  mwSize lenx; 


  
  /*  check for proper number of arguments */
  /* NOTE: You do not need an else statement when using mexErrMsgTxt
     within an if statement, because it will never get to the else
     statement if mexErrMsgTxt is executed. (mexErrMsgTxt breaks you out of
     the MEX-file) */
  if(nrhs!=6) 
    mexErrMsgTxt("6 inputs required.");
  if(nlhs!=1) 
    mexErrMsgTxt("One output required.");
  
  /* check to make sure the first input argument is a scalar */
//   if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
//       mxGetN(prhs[0])*mxGetM(prhs[0])!=1 ) {
//     mexErrMsgTxt("Input x must be a scalar.");
//   }
  
  /*  get the scalar input x */
 
  /*  create a pointer to the input matrix y */
 x = mxGetPr(prhs[0]);
 y = mxGetPr(prhs[1]);
 frames = mxGetPr(prhs[2]);
  
  /*  get the dimensions of the matrix input y */
  lenx = mxGetM(prhs[0]);

 /* printf("length of positons: %i\n",lenx);*/
  
dX=mxGetScalar(prhs[3]);
dT=mxGetScalar(prhs[4]);
dxsigma=mxGetPr(prhs[5]);
  /*  set the output pointer to the output matrix */
//  printf("inputs associated\n");

  plhs[0] = mxCreateNumericMatrix(lenx,1,mxDOUBLE_CLASS,mxREAL);

//     printf("outputmatrix created\n");
  
  
  /*  create a C pointer to a copy of the output matrix */
  list = mxGetData(plhs[0]);
//    printf("outputmatrix created 2, call routine\n");
  /*  call the C subroutine */
  cs(list,x,y,frames,dX,dT,dxsigma,lenx);

}
