/*
Fcrosspairs_ts_binned.c

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
#include <stdio.h>

#define OK 0
#define ERR_OVERFLOW 1

double sqrt();
extern bool utIsInterruptPending();

/* this one is for inputs ordered by x (hence _st_ for space-time) */
void Fcrosspairs_st_binned(n1, x1, y1, t1, n2, x2, y2, t2, rmax, nrout, taumin, taumax, ntout,
        counts)
     /* inputs */
     int n1, n2;
     double *x1, *y1, *x2, *y2, rmax, *t1, *t2, taumin, taumax;
     /* outputs */
     int nrout, ntout;
     uint64_t *counts;
{
  int k, i, j, jfirst; /* note points are ordered by x, not t */
  double x1i, y1i, t1i, r2max, xfirst, dx, dy, dx2, d2, dt, dr, d, tau;
 
  r2max = rmax * rmax;
  dt = (taumax - taumin)/((double) ntout);
  dr = rmax / ((double) nrout);
 
  if(n1 == 0 || n2 == 0) 
    return;

  jfirst = 0;
  
  i = 0;
  
  for(; i < n1; i++) {
      /* check interrupt (CTRL-C) status every 1024 pts and abort if pending */
      if ( (!(i & 0x3ff)) && utIsInterruptPending() ) break;
      
      x1i = x1[i];
      y1i = y1[i];
      t1i = t1[i];

      /* adjust starting position jfirst */
      xfirst = x1i - rmax;
      while((x2[jfirst] < xfirst) && (jfirst + 1 < n2))
        ++jfirst;

      /* process from j=jfirst until dx > rmax */
      for(j=jfirst; j < n2; j++) {
        dx = x2[j] - x1i;
        if (dx >= rmax)
          break;

        tau = t2[j] - t1i;

        k = floor((tau - taumin)/dt);
        if (k >= ntout || k < 0)
            continue;

        dx2 = dx * dx;
        /* check spatial constraint */
        dy = y2[j] - y1i;
        d2 = dx2 + dy * dy;
        if(d2 < r2max) {
          /* add this (i, j) pair to output */
          d = sqrt(d2);
          counts[k*nrout + ((int) (d/dr))] += 1;
        }
      }
  }
}

/* this one is ordered by t (hence _ts_ for time-space) */
void Fcrosspairs_ts_binned(n1, x1, y1, t1, n2, x2, y2, t2, rmax, nrout, taumin, taumax, ntout,
        counts)
     /* inputs */
     int n1, n2;
     double *x1, *y1, *x2, *y2, rmax, *t1, *t2, taumin, taumax;
     /* outputs */
     int nrout, ntout;
     uint64_t *counts;
{
  int k, i, j, jfirst; /* note points are ordered by t */
  double x1i, y1i, t1i, r2max, tfirst, dx, dy, dx2, d2, dt, dr, d, tau;
 
  r2max = rmax * rmax;
  dt = (taumax - taumin)/((double) ntout);
  dr = rmax / ((double) nrout);
 
  if(n1 == 0 || n2 == 0) 
    return;

  jfirst = 0;
  
  i = 0;
  
  for(; i < n1; i++) {
      /* check interrupt (CTRL-C) status every 1024 pts and abort if pending */
      if ( (!(i & 0x3ff)) && utIsInterruptPending() ) break;
      
      x1i = x1[i];
      y1i = y1[i];
      t1i = t1[i];

      /* adjust starting position jfirst */
      tfirst = t1i + taumin;
      while((t2[jfirst] < tfirst) && (jfirst + 1 < n2))
        ++jfirst;

      /* process from j=jfirst until dx > rmax */
      for(j=jfirst; j < n2; j++) {
        tau = t2[j] - t1i;
        if (tau >= taumax)
            break;

        dx = x2[j] - x1i;

        /* check spatial constraint */
        if (dx < rmax) {
          dx2 = dx * dx;
          dy = y2[j] - y1i;
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

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    int n1, n2;
    double *x1, *y1, *t1, *x2, *y2, *t2;
    double rmax, taumin, taumax;
    uint64_t *counts, *counts_allr;
    int nrout, ntout, *tmp, sortbyt;
    int i;

    /* RHS args are:
     *  0   x1 // coords for first pp
     *  1   y1
     *  2   t1
     *  3   x2 // coords for second pp
     *  4   y2
     *  5   t2
     *  6   rmax
     *  7   nrout
     *  8   taumin
     *  9   taumax
     * 10   ntout
     * 11   sortbyt = 0 for space, 1 for time
     */

    /* LHS args are:
     * 0,1    (r,t)-hist, (t)-hist at all r.
     */

    /*
    // TODO: clearly this needs some more safety checks...
    // For now we just assume stuff is in the right formats.
    // */
    n1 = mxGetNumberOfElements(prhs[0]);
    n2 = mxGetNumberOfElements(prhs[3]);

    x1 = (double *) mxGetPr(prhs[0]);
    y1 = (double *) mxGetPr(prhs[1]);
    t1 = (double *) mxGetPr(prhs[2]);
    x2 = (double *) mxGetPr(prhs[3]);
    y2 = (double *) mxGetPr(prhs[4]);
    t2 = (double *) mxGetPr(prhs[5]);

    rmax = mxGetScalar(prhs[6]);

    tmp = (int *) mxGetPr(prhs[7]);
    nrout = *tmp;
    
    taumin = mxGetScalar(prhs[8]);
    taumax = mxGetScalar(prhs[9]);

    tmp = (int *) mxGetPr(prhs[10]);
    ntout = *tmp;
    tmp = (int *) mxGetPr(prhs[11]);
    sortbyt = *tmp;

    plhs[0] = mxCreateNumericMatrix(nrout, ntout, mxUINT64_CLASS, mxREAL);

    counts = (uint64_t *) mxGetPr(plhs[0]);

    if (sortbyt) {
        Fcrosspairs_ts_binned(n1, x1, y1, t1, n2, x2, y2, t2,
                rmax, nrout, taumin, taumax, ntout, counts);
    } else {
        Fcrosspairs_st_binned(n1, x1, y1, t1, n2, x2, y2, t2,
                rmax, nrout, taumin, taumax, ntout, counts);
    }

    return;
}
