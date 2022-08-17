/*
Fclosepairs_st.c

Copyright (C) 2022 Thomas Shaw, Frank Fazekas and Sarah Veatch
This file is part of SPACETIME RESOLUTION.
SPACETIME RESOLUTION is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
SPACETIME RESOLUTION is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
%along with SPACETIME RESOLUTION.  If not, see <https://www.gnu.org/licenses/>.
ADAPTED FROM:

closepair.c

$Revision : 1.34 $     $Date : 2018 / 12 / 18 02 : 43 : 11 $
Copyright(C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001 - 2018
Licence : GNU Public Licence >= 2
*/

#include "mex.h"
#include "matrix.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#define OK 0
#define ERR_OVERFLOW 1

double sqrt();

extern bool utIsInterruptPending();

void Fclosepairs_ts(n, x, y, t, rmax, taumin, taumax, noutmax,
        nout, dxout, dyout, dtout, err)
     /* inputs */
     int n, noutmax;
     double *x, *y, rmax, *t, taumin, taumax;
     /* outputs */
     int *nout, *err;
     double *dxout, *dyout, *dtout;
{
  int k, i, j, jfirst; /* note points are ordered by t, not x */
  double xi, yi, ti, r2max, tfirst, dx, dy, dx2, d2, tau;
 
  r2max = rmax * rmax;
 
  *nout = 0;
  *err = 0;
  if(n == 0) 
    return;

  jfirst = 0;
  
  i = 0;
  k = 0; /* k indexes the output arrays */
  
  for(; i < n; i++) {
      /* check interrupt (CTRL-C) status every 1024 pts and abort if pending */
      if ( (!(i & 0x3ff)) && utIsInterruptPending() ) break;
      
      xi = x[i];
      yi = y[i];
      ti = t[i];

      /* adjust starting position jfirst */
      tfirst = ti + taumin;
      while((t[jfirst] < tfirst) && (jfirst + 1 < n)) /* +1 is because we are going to add 1 if true */
        ++jfirst;

      /* process from j=jfirst until dt > taumax */
      for(j=jfirst; j < n; j++) {
        if (j==i)
          continue; /* skip the i,i pair */

        tau = t[j] - ti;
        if (tau >= taumax)
            break;

        dx = x[j] - xi;
        /* check spatial constraint */
        if (dx < rmax) {
          dx2 = dx * dx;
          dy = y[j] - yi;
          d2 = dx2 + dy * dy;
          if(d2 < r2max) {
            /* add this (i, j) pair to output */
            dxout[k] = dx;
            dyout[k] = dy;
            dtout[k] = tau;

            ++k;
            if (k == noutmax) {
                *nout = k;
                *err = 1;
                return;
            }
          }
        }  
      }
  }
  *nout = k;
  return;
}

void Fclosepairs_st(n, x, y, t, rmax, taumin, taumax, noutmax,
        nout, dxout, dyout, dtout, err)
     /* inputs */
     int n, noutmax;
     double *x, *y, rmax, *t, taumin, taumax;
     /* outputs */
     int *nout, *err;
     double *dxout, *dyout, *dtout;
{
  int k, i, j, jfirst; /* note points are ordered by x, not t */
  double xi, yi, ti, r2max, xfirst, dx, dy, dx2, d2, tau;
 
  r2max = rmax * rmax;
 
  *nout = 0;
  *err = 0;
  if(n == 0) 
    return;

  jfirst = 0;
  
  i = 0;
  k = 0; /* k indexes the output arrays */
  
  for(; i < n; i++) {
      /* check interrupt (CTRL-C) status every 1024 pts and abort if pending */
      if ( (!(i & 0x3ff)) && utIsInterruptPending() ) break;
      
      xi = x[i];
      yi = y[i];
      ti = t[i];

      /* adjust starting position jfirst */
      xfirst = xi - rmax;
      while((x[jfirst] < xfirst) && (jfirst + 1 < n))
        ++jfirst;

      /* process from j=jfirst until dx > rmax */
      for(j=jfirst; j < n; j++) {
        if (j==i)
          continue; /* skip the i,i pair */

        dx = x[j] - xi;
        if (dx >= rmax)
          break;

        tau = t[j] - ti;
        if (tau >= taumax || tau < taumin)
            continue;

        dx2 = dx * dx;
        /* check spatial constraint */
        dy = y[j] - yi;
        d2 = dx2 + dy * dy;
        if(d2 < r2max) {
          /* add this (i, j) pair to output */
          dxout[k] = dx;
          dyout[k] = dy;
          dtout[k] = tau;

          ++k;
          if (k == noutmax) {
              *nout = k;
              *err = 1;
              return;
          }
        }
      }
  }
  *nout = k;
  return;
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    int n;
    double *x, *y, *t;
    double rmax, taumin, taumax;
    double *dxout, *dyout, *dtout;
    int *tmp, sortbyt, noutmax, *err, nout;
    int i;

    /* RHS args are:
     *  0   x // coords
     *  1   y
     *  2   t
     *  3   rmax
     *  4   taumin
     *  5   taumax
     *  6   noutmax
     *  7   sortbyt
     */

    /* LHS args are:
     * 0    dxout
     * 1    dyout
     * 2    dtout
     * 3    err
     */

    /*
    // TODO: clearly this needs some more safety checks...
    // For now we just assume stuff is in the right formats.
    // */
    n = mxGetNumberOfElements(prhs[0]);

    x = (double *) mxGetPr(prhs[0]);
    y = (double *) mxGetPr(prhs[1]);
    t = (double *) mxGetPr(prhs[2]);

    rmax = mxGetScalar(prhs[3]);

    taumin = mxGetScalar(prhs[4]);
    taumax = mxGetScalar(prhs[5]);

    tmp = (int *) mxGetPr(prhs[6]);
    noutmax = *tmp;
    tmp = (int *) mxGetPr(prhs[7]);
    sortbyt = *tmp;

    plhs[0] = mxCreateNumericMatrix(noutmax, 1, mxDOUBLE_CLASS, mxREAL);
    plhs[1] = mxCreateNumericMatrix(noutmax, 1, mxDOUBLE_CLASS, mxREAL);
    plhs[2] = mxCreateNumericMatrix(noutmax, 1, mxDOUBLE_CLASS, mxREAL);
    plhs[3] = mxCreateNumericMatrix(1,1, mxINT64_CLASS, mxREAL);

    dxout = (double *) mxGetPr(plhs[0]);
    dyout = (double *) mxGetPr(plhs[1]);
    dtout = (double *) mxGetPr(plhs[2]);
    err = (int *) mxGetPr(plhs[3]);

    if (sortbyt) {
        Fclosepairs_ts(n, x, y, t, rmax, taumin, taumax, noutmax,
                &nout, dxout, dyout, dtout, err);
    } else {
        Fclosepairs_st(n, x, y, t, rmax, taumin, taumax, noutmax,
                &nout, dxout, dyout, dtout, err);
    }

   dxout = mxRealloc(dxout, nout * sizeof(double));
   dyout = mxRealloc(dyout, nout * sizeof(double));
   dtout = mxRealloc(dtout, nout * sizeof(double));

   mxSetPr(plhs[0], dxout);
   mxSetPr(plhs[1], dyout);
   mxSetPr(plhs[2], dtout);

   for (i=0; i<3; i++)
       mxSetM(plhs[i], nout);

   return;
}
