#include "mex.h"
#include <stdio.h>
#include "math.h"

void countlocs(double *vout,double *list,mwSize numlocs)
{
    double ng=0,groupc;
    mwSize k,l,inds=0;
    
    groupc=list[0];
    for(k=0;k<numlocs;k++)
        {
        if (list[k]!=groupc)
            {
            groupc=list[k];
            for (l=inds;l<k;l++)
                {
                vout[l]=ng;
                }
            ng=0;
            inds=k;
            }
        ng++;
          
        }
    vout[numlocs-1]=ng;
}



//  groupc=listsort(1);
//             ng=0;
//             inds=1;
//             for k=1:length(numbergroup)
//                 if listsort(k)~=groupc
//                     groupc=listsort(k);
//                     numbergroup(inds:k-1)=ng;
//                     ng=0;
//                     inds=k;
//                 end
//                 ng=ng+1;
//             end


/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  double *list,*vout;
  mwSize numlocs;    
  if(nrhs!=1) 
    mexErrMsgTxt("1 inputs required.");
  if(nlhs!=1) 
    mexErrMsgTxt("1 output required.");
  
  /*  create a pointer to the input matrix y */
 list = mxGetPr(prhs[0]);

  numlocs = mxGetM(prhs[0]);

  plhs[0] = mxCreateNumericMatrix(numlocs,1,mxDOUBLE_CLASS,mxREAL);
  vout = mxGetData(plhs[0]);
  countlocs(vout,list,numlocs);

}
