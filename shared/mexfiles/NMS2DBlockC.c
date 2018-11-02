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
void  NMS(double *maxima,float *img,short sx,short sy,short n, mwSize maxmax)
{
short i1,i2,j1,j2,mi,mj,stop;
long maxind;

short *where;


maxind=0;

for(i1=n;i1<sy-2*n;i1+=n) //carful if sx and sy are exchanged 
{
    for(j1=n;j1<sx-2*n;j1+=n) 
    {
        mi=i1; mj=j1;
        for(i2=i1;i2<i1+n;i2+=1)
        {
            for(j2=j1;j2<j1+n;j2+=1)
            {
                if(*(img+i2*sx+j2) > *(img+mi*sx+mj)) 
                {
                    mi=i2;
                    mj=j2;
                }
            }
        }
        
        
//     for i1=n+1:n+1:simg(1)-n-n2;
// %     i1
//     for j1=n+1:n+1:simg(2)-n-n2;
        
//         for i2=i1:i1+n
//             for j2=j1:j1+n
//                 if imin(i2,j2)>imin(mi,mj)
//                     mi=i2;mj=j2;
//                 end
//             end
//         end
        stop=0;
        for(i2=mi-n;i2<mi+n;i2+=1) //+1???
        {
            for(j2=mj-n;j2<mj+n;j2+=1)
            {
               
                if(*(img+i2*sx+j2) > *(img+mi*sx+mj)) 
                {
                    stop=1;
                    break;
                }
            }
            if (stop==1) break;
        }
        if ((stop==0)&(maxind<maxmax))
        {
            
            *(maxima+maxind)=mj+1;
//             maxind++;
            *(maxima+maxind+maxmax)=mi+1;
            *(maxima+maxind+2*maxmax)=*(img+mi*sx+mj);
            maxind++;
        }
    }
//     *(maxima)=5;
}
            
        

//         for i2=mi-n2:mi+n2 %future: exclude block already tested. not done here!
//             for j2=mj-n2:mj+n2 %or go backwards, then  already tested ones last %think about investigating only mj-1 and mj+1: neares neighbours: image already filtered
//                 if imin(i2,j2)>imin(mi,mj) 
//                     stop=1;
//                     break
//                 end
//             end
//             if stop
//                 break
//             end
//         end

//         if ~stop
//             maxima(maxind,:)=[mi mj imin(mi,mj)];
//             
//             maxind=maxind+1;
//         end

//     end
// end
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  float *img; //filtered image
  short sx,sy, n; //size of image 
  double *maxima;
mwSize maxmax;
mwSize dims[2];

  
  /*  check for proper number of arguments */
  /* NOTE: You do not need an else statement when using mexErrMsgTxt
     within an if statement, because it will never get to the else
     statement if mexErrMsgTxt is executed. (mexErrMsgTxt breaks you out of
     the MEX-file) */
  if(nrhs!=3) 
    mexErrMsgTxt("Three inputs required.");
  if(nlhs!=1) 
    mexErrMsgTxt("One output required.");
  
  /* check to make sure the first input argument is a scalar */
//   if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
//       mxGetN(prhs[0])*mxGetM(prhs[0])!=1 ) {
//     mexErrMsgTxt("Input x must be a scalar.");
//   }
  
  /*  get the scalar input x */
 
  /*  create a pointer to the input matrix y */
 img = mxGetPr(prhs[0]);
  
  /*  get the dimensions of the matrix input y */
  sx = mxGetM(prhs[0]);
  sy = mxGetN(prhs[0]);
  
  
  n=mxGetScalar(prhs[1]);
  maxmax=mxGetScalar(prhs[2]);
  /*  set the output pointer to the output matrix */
dims[0]=maxmax;
dims[1]=2;

  plhs[0] = mxCreateNumericMatrix(maxmax,3,mxDOUBLE_CLASS,mxREAL);

  
  /*  create a C pointer to a copy of the output matrix */
  maxima = mxGetData(plhs[0]);
  
  /*  call the C subroutine */
  NMS(maxima,img,sx,sy,n,maxmax);
  
}
