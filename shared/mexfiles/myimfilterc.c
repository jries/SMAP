#include "mex.h"
#include <stdio.h>


/*void correlate(double *n1, double *G, mwSize lenG, mwSize lenn)*/
void imfilter(float *in, float *out,mwSize sx,mwSize sy,short dn)
{
short x1,x2,y1,y2;
float norm;
norm=(2*dn+1)*(2*dn+1);
       
for(x1=dn;x1<sy-dn;x1++)
{
    for(y1=dn;y1<sx-dn;y1++)
    {
        for(x2=x1-dn;x2<=x1+dn;x2++)
        {
            for(y2=y1-dn;y2<=y1+dn;y2++)
            {
                *(out+x1*sx+y1)+=*(in+x2*sx+y2);
            }
        }
        *(out+x1*sx+y1)/=norm; 
    }
}
}           

                    


/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  float *in,*out; 
  short dn; //already dn here


mwSize sx,sy;
  
  /*  check for proper number of arguments */
  /* NOTE: You do not need an else statement when using mexErrMsgTxt
     within an if statement, because it will never get to the else
     statement if mexErrMsgTxt is executed. (mexErrMsgTxt breaks you out of
     the MEX-file) */
  if(nrhs!=2) 
    mexErrMsgTxt("Two inputs required.");
  if(nlhs!=1) 
    mexErrMsgTxt("One output required.");
  
  /* check to make sure the first input argument is a scalar */
//   if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
//       mxGetN(prhs[0])*mxGetM(prhs[0])!=1 ) {
//     mexErrMsgTxt("Input x must be a scalar.");
//   }
  
  /*  get the scalar input x */
 
  /*  create a pointer to the input matrix y */
 in = mxGetPr(prhs[0]);
  
  /*  get the dimensions of the matrix input y */
  sx = mxGetM(prhs[0]);
  sy = mxGetN(prhs[0]);
  
  
  dn=mxGetScalar(prhs[1]);

  /*  set the output pointer to the output matrix */

  plhs[0] = mxCreateNumericMatrix(sx,sy,mxSINGLE_CLASS,mxREAL);

  
  /*  create a C pointer to a copy of the output matrix */
  out = mxGetData(plhs[0]);
  
  /*  call the C subroutine */
  imfilter(in,out,sx,sy,dn);
  
}
