#include "mex.h"
#include <stdio.h>
#include "math.h"

void sumcombine(double *vout,double *vin,double *list,mwSize numlocs, mwSize numcombined)
{
    double v, thisparticle,oldparticle,v0=0;
    mwSize part=0,k;
    
    oldparticle=list[0];
    for(k=0;k<numlocs;k++)
        {
        v=vin[k];
        thisparticle=list[k];
        if (oldparticle==thisparticle)
            {
            v0=v0+v;
            }
        else
            {
            vout[part]=v0;
            part++;
            oldparticle=thisparticle;
            
            v0=v;
            }       
        }
    vout[part]=v0;
}


/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  double *vin,*list,*vout;
  mwSize numlocs, numcombined;    
  if(nrhs!=2) 
    mexErrMsgTxt("2 inputs required.");
  if(nlhs!=1) 
    mexErrMsgTxt("1 output required.");
  
  /*  create a pointer to the input matrix y */
 vin = mxGetPr(prhs[0]);
 list = mxGetPr(prhs[1]);

  numlocs = mxGetM(prhs[0]);
numcombined=list[numlocs-1]-list[0]+1;
  plhs[0] = mxCreateNumericMatrix(numcombined,1,mxDOUBLE_CLASS,mxREAL);
  vout = mxGetData(plhs[0]);
  sumcombine(vout,vin,list,numlocs,numcombined);

}
