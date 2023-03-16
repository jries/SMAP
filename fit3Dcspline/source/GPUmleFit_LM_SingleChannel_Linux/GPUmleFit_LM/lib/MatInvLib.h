/*!
 * \file MatInvLib.h
 * \author Keith Lidke
 * \date January 10, 2010
 * \brief This provides the matrix inversion helper functions.
 */

#ifndef MATINVLIB_H
#define MATINVLIB_H
//__device__ void kernel_MatInv4(float *, float *, float *,float *);
//__device__ void kernel_MatInv5(float *, float *, float *,float *);
//__device__ void kernel_MatInvN(float *, float *, float *,float *);

//*******************************************************************************************
__device__ inline void kernel_MatInvN(float * M, float * Minv, float * DiagMinv, int sz) {
	/*! 
	 * \brief nxn partial matrix inversion
	 * \param M matrix to inverted
	 * \param Minv inverted matrix result
	 * \param DiagMinv just the inverted diagonal
	 * \param sz size of the matrix
	 */
    int ii, jj, kk, num, b;
    float tmp1=0;
    float yy[25];
    
    for (jj = 0; jj < sz; jj++) {
        //calculate upper matrix
        for (ii=0;ii<=jj;ii++)
            //deal with ii-1 in the sum, set sum(kk=0->ii-1) when ii=0 to zero
            if (ii>0) {
            for (kk=0;kk<=ii-1;kk++) tmp1+=M[ii+kk*sz]*M[kk+jj*sz];
            M[ii+jj*sz]-=tmp1;
            tmp1=0;
            }
  
        for (ii=jj+1;ii<sz;ii++)
            if (jj>0) {
            for (kk=0;kk<=jj-1;kk++) tmp1+=M[ii+kk*sz]*M[kk+jj*sz];
            M[ii+jj*sz]=(1/M[jj+jj*sz])*(M[ii+jj*sz]-tmp1);
            tmp1=0;
            }
            else { M[ii+jj*sz]=(1/M[jj+jj*sz])*M[ii+jj*sz]; }
    }
 
    tmp1=0;
    
    for (num=0;num<sz;num++) {
        // calculate yy
        if (num==0) yy[0]=1;
        else yy[0]=0;
        
        for (ii=1;ii<sz;ii++) {
            if (ii==num) b=1;
            else b=0;
            for (jj=0;jj<=ii-1;jj++) tmp1+=M[ii+jj*sz]*yy[jj];
            yy[ii]=b-tmp1;
            tmp1=0;    
        }
        
        // calculate Minv
        Minv[sz-1+num*sz]=yy[sz-1]/M[(sz-1)+(sz-1)*sz];
        
        for (ii=sz-2;ii>=0;ii--) {
            for (jj=ii+1;jj<sz;jj++) tmp1+=M[ii+jj*sz]*Minv[jj+num*sz];
            Minv[ii+num*sz]=(1/M[ii+ii*sz])*(yy[ii]-tmp1);
            tmp1=0;
        }
    }
    
    if (DiagMinv) for (ii=0;ii<sz;ii++) DiagMinv[ii]=Minv[ii*sz+ii];
    
    return;
    
}

#endif