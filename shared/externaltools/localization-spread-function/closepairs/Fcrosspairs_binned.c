/*
Fcrosspairs_binned.c

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

void Fcrosspairs_st_binned(n1, x1, y1, n2, x2, y2, rmax, nrout,
        counts)
     /* inputs */
     int n1, n2;
     double *x1, *y1, *x2, *y2, rmax;
     /* outputs */
     int nrout;
     uint64_t *counts;
{
  int k, i, j, jfirst; /* note points are ordered by x, not t */
  double x1i, y1i, r2max, xfirst, dx, dy, dx2, d2, dr, d;
 
  r2max = rmax * rmax;
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

      /* adjust starting position jfirst */
      xfirst = x1i - rmax;
      while((x2[jfirst] < xfirst) && (jfirst + 1 < n2))
        ++jfirst;

      /* process from j=jfirst until dx > rmax */
      for(j=jfirst; j < n2; j++) {
        dx = x2[j] - x1i;
        if (dx >= rmax)
          break;

        dx2 = dx * dx;
        /* check spatial constraint */
        dy = y2[j] - y1i;
        d2 = dx2 + dy * dy;
        if(d2 < r2max) {
          /* add this (i, j) pair to output */
          d = sqrt(d2);
          counts[((int) (d/dr))] += 1;
        }
      }
  }
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    int n1, n2;
    double *x1, *y1, *x2, *y2;
    double rmax;
    uint64_t *counts, *counts_allr;
    int nrout, *tmp;
    int i;

    /* RHS args are:
     *  0   x1 // coords for first pp
     *  1   y1
     *  2   x2 // coords for second pp
     *  3   y2
     *  4   rmax
     *  5   nrout
     */

    /* LHS args are:
     * 0    (r)-hist
     */

    /*
    // TODO: clearly this needs some more safety checks...
    // For now we just assume stuff is in the right formats.
    // */
    n1 = mxGetNumberOfElements(prhs[0]);
    n2 = mxGetNumberOfElements(prhs[2]);

    x1 = (double *) mxGetPr(prhs[0]);
    y1 = (double *) mxGetPr(prhs[1]);
    x2 = (double *) mxGetPr(prhs[2]);
    y2 = (double *) mxGetPr(prhs[3]);

    rmax = mxGetScalar(prhs[4]);

    tmp = (int *) mxGetPr(prhs[5]);
    nrout = *tmp;
    
    plhs[0] = mxCreateNumericMatrix(nrout, 1, mxUINT64_CLASS, mxREAL);

    counts = (uint64_t *) mxGetPr(plhs[0]);

    Fcrosspairs_st_binned(n1, x1, y1, n2, x2, y2,
                rmax, nrout, counts);

    return;
}
