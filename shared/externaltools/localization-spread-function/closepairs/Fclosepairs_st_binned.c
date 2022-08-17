/*
Fclosepairs_ts_binned.c

Copyright (C) 2021 Frank Fazekas, Thomas Shaw, and Sarah Veatch
This file is part of MEAN SHIFT DRIFT CORRECTION.
MEAN SHIFT DRIFT CORRECTION is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
MEAN SHIFT DRIFT CORRECTION is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
%along with MEAN SHIFT DRIFT CORRECTION.  If not, see <https://www.gnu.org/licenses/>.
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

void Fclosepairs_ts_binned(n, x, y, t, rmax, nrout, taumin, taumax, ntout,
        counts)
     /* inputs */
     int n;
     double *x, *y, rmax, *t, taumin, taumax;
     /* outputs */
     int nrout, ntout;
     uint64_t *counts;
{
  int k, i, j, jfirst; /* note points are ordered by t, not x */
  double xi, yi, ti, r2max, tfirst, dx, dy, dx2, d2, dt, dr, d, tau;
 
  r2max = rmax * rmax;
  dt = (taumax - taumin)/((double) ntout);
  dr = rmax / ((double) nrout);
 
  if(n == 0) 
    return;

  jfirst = 0;
  
  i = 0;
  
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

      /* process from j=jfirst until dx > rmax */
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
            d = sqrt(d2);

            k = floor((tau - taumin)/dt);
            counts[k*nrout + ((int) (d/dr))] += 1;
          }
        }  
      }
  }
}

void Fclosepairs_st_binned(n, x, y, t, rmax, nrout, taumin, taumax, ntout,
        counts)
     /* inputs */
     int n;
     double *x, *y, rmax, *t, taumin, taumax;
     /* outputs */
     int nrout, ntout;
     uint64_t *counts;
{
  int k, i, j, jfirst; /* note points are ordered by x, not t */
  double xi, yi, ti, r2max, xfirst, dx, dy, dx2, d2, dt, dr, d, tau;
 
  r2max = rmax * rmax;
  dt = (taumax - taumin)/((double) ntout);
  dr = rmax / ((double) nrout);
 
  if(n == 0) 
    return;

  jfirst = 0;
  
  i = 0;
  
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

        k = (int) ((tau - taumin)/dt);

        dx2 = dx * dx;
        /* check spatial constraint */
        dy = y[j] - yi;
        d2 = dx2 + dy * dy;
        if(d2 < r2max) {
          /* add this (i, j) pair to output */
          d = sqrt(d2);
          counts[k*nrout + ((int) (d/dr))] += 1;
        }
      }
  }
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    int n;
    double *x, *y, *t;
    double rmax, taumin, taumax;
    uint64_t *counts;
    int nrout, ntout, *tmp, sortbyt;
    int i;

    /* RHS args are:
     *  0   x // coords
     *  1   y
     *  2   t
     *  3   rmax
     *  4   nrout
     *  5   taumin
     *  6   taumax
     *  7   ntout
     *  8   sortbyt
     */

    /* LHS args are:
     * 0    (r,t)-hist
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

    tmp = (int *) mxGetPr(prhs[4]);
    nrout = *tmp;
    
    taumin = mxGetScalar(prhs[5]);
    taumax = mxGetScalar(prhs[6]);

    tmp = (int *) mxGetPr(prhs[7]);
    ntout = *tmp;
    tmp = (int *) mxGetPr(prhs[8]);
    sortbyt = *tmp;

    plhs[0] = mxCreateNumericMatrix(nrout, ntout, mxUINT64_CLASS, mxREAL);

    counts = (uint64_t *) mxGetPr(plhs[0]);

    if (sortbyt) {
        Fclosepairs_ts_binned(n, x, y, t, rmax, nrout, taumin, taumax, ntout,
                counts);
    } else {
        Fclosepairs_st_binned(n, x, y, t, rmax, nrout, taumin, taumax, ntout,
                counts);
    }

    return;
}
