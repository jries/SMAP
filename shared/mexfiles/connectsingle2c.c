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
void  cs(double *list, double *x, double *y, double *frames,double dX,long dT,mwSize maxactive,mwSize lenx)
{
    long thisentry,stopnow,particlenumber,numpart;
    long particlefound,testentry,numdark; 
    double xh,yh,frh,frtest; //coordinates of main particle
    
    
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
            if ((y[testentry]>yh-dX) && (y[testentry]<yh+dX) && (list[testentry]==0)) 
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
    
//mxArray *xac,*yac; 
// double xc, yc;
//         
// long numpartcurrent,entry,fr,particlefound,ap,np,error,api;
// 
// double *xac,*yac;
// long *lastframe,*particlenumber,*appearenceparticle ;

// xac = mxCalloc(lenx , sizeof (double));
// printf("1");
// yac = mxCalloc(lenx , sizeof (double));
// printf("2");
// lastframe = mxCalloc(lenx , sizeof (long));
// printf("3");
// particlenumber =  mxCalloc(lenx , sizeof (long));
// printf("4");
// appearenceparticle =  mxCalloc(lenx , sizeof (long));
// printf("5 ");


// printf("lastf %i ", lastframe[2]);

// printf("start routine\n");

// error=0;
// numpartcurrent=0;
// if (xac == NULL|yac==NULL|lastframe==NULL|particlenumber==NULL|appearenceparticle==NULL) {
// printf("not sufficient memory\n");
// } else {
//     /* Allocation succeeded.  Do something. */
// printf("b: sufficient memory allocated\n");
// 
// for(entry=0;entry<lenx;entry++) //look at all entry locs in list. entry counts position
// {
// //     printf("\n %i,",entry);
//     xc=x[entry]; //look now at 1 loc at position xc yc in frame fr
//     yc=y[entry];
//     fr=frames[entry];
//     particlefound=0;  //now there is no corresponding particle later
// //     printf("before api, frame: %i", fr);
//     for(api=0;api<maxactive;api++) //look from here to maxactive entries later. api counts entries
//     {
// //         printf(" :%i! ",lastframe[api]);
//         if(lastframe[api]>0&&lastframe[api]<fr) //lastframe initialized to 0
//         {
//             if(xac[api]-dX<xc && xac[api]+dX>xc && yac[api]-dX<yc && yac[api]+dX>yc)
//                 
//             {
// //                 printf("pf ");
//                 particlefound=1;
//                 ap=api;
//             }
//             else
//             {
//                 if(lastframe[api]<fr-dT)
//                 {
// //                     printf("df ");
//                     lastframe[api]=0;
//                 }
//             }
//         }
//     }
// //     printf("\n pos 1");
//     if(particlefound==1)
//     {
// //         printf("ap ");
//         list[entry]=particlenumber[ap];
//         xac[ap]=(xac[ap]+xc)/2;yac[ap]=(yc+yac[ap])/2;lastframe[ap]=fr;
//        appearenceparticle[ap]++;
//        appearencelist[entry]=appearenceparticle[ap];
//     }
//     else
//     {
// //         printf("np ");
//         numpartcurrent++;
//         list[entry]=numpartcurrent;
//         //find next entry whichis 0
//         for(np=0;np<maxactive;np++)
//         {
//             if(lastframe[np]==0) break;
//         }
//         if(np==maxactive-1)
//         {
//             break;
//         printf("stopped 2");
//         }
//         else
//         {
// //             printf("np: %i",np);
// //             printf("fr: %i",fr);
//             
//             xac[np]=xc;yac[np]=yc;lastframe[np]=fr;
//             particlenumber[np]=numpartcurrent;
//                    appearenceparticle[np]=1;
//        appearencelist[entry]=appearenceparticle[np];
// //             printf("lf(np): %i",lastframe[np]);
//         }
//   
// }
// }
// 
//  mxFree(xac);
//  mxFree(yac);
//  mxFree(lastframe);
//  mxFree(particlenumber);
//  mxFree(appearenceparticle);
// 
// }
// }

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  double *list,*x,*y,dX,*frames,*appearencelist; 
  long int  dT;
  mwSize lenx, maxactive; 


  
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
maxactive=mxGetScalar(prhs[5]);
  /*  set the output pointer to the output matrix */
//  printf("inputs associated\n");

  plhs[0] = mxCreateNumericMatrix(lenx,1,mxDOUBLE_CLASS,mxREAL);

//     printf("outputmatrix created\n");
  
  
  /*  create a C pointer to a copy of the output matrix */
  list = mxGetData(plhs[0]);
//    printf("outputmatrix created 2, call routine\n");
  /*  call the C subroutine */
  cs(list,x,y,frames,dX,dT,maxactive,lenx);

}
