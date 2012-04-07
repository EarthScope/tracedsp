/*********************************************************************
 * envelope.c
 *
 * Calculate envelope and Hilber transforms of time series data.
 *
 * The Hilbert transform and supporting convolution code used below
 * were originally written by Dave Hale and extracted from the CWP
 * library source code.  All legal notices and most of the notes
 * remain intact.  The source has been modified to use doubles, do
 * error checking, simplified for the needs of this code base and to
 * fit the style.
 *
 * Modified: 2011.136
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "envelope.h"

/* Pi, with comical precision */
#define PI 3.1415926535897932384626433832795

static void conv (int lx, int ifx, double *x,
		  int ly, int ify, double *y,
		  int lz, int ifz, double *z);

/********************************************************************
 * envelope:
 *
 * Calculate envelope of a time series as square root of x(n)^2 + y(n)^2
 * where x(n) is the original sample and y(n) is the Hilbert transform.
 *
 * Arguments:
 *   data       : array of input and output data samples as doubles
 *   npts       : number of samples
 *
 * Returns samples in output on success and -1 on error
 ********************************************************************/
int
envelope (double *data, int npts)
{
  double *hil;
  int idx;
  
  if ( ! data || npts < 1 )
    {    
      fprintf (stderr, "envelope(): No data buffer specified\n");
      return -1;
    }
  
  if ( (hil = (double *) malloc (npts * sizeof(double))) == NULL )
    {
      fprintf (stderr, "envelope(): Cannot allocate buffer for Hilbert transform\n");
      return -1;
    }
  
  /* Calculate Hilbert transform */
  if ( hilbert (data, hil, npts) < 0 )
    {    
      fprintf (stderr, "envelope(): Error calculating Hilbert transform\n");
      free (hil);
      return -1;
    }
  
  /* Calculate envelope */
  for (idx=0; idx < npts; idx++)
    {
      data[idx] = sqrt(data[idx]*data[idx] + hil[idx]*hil[idx]);
    }
  
  free (hil);
  
  return npts;
}  /* End of envelope()*/


/* Originally from hilbert.c */
/* Copyright (c) Colorado School of Mines, 2010.*/
/* All rights reserved.                       */

/*********************** self documentation **********************/
/*****************************************************************************
HILBERT - Compute Hilbert transform y of x

hilbert		compute the Hilbert transform

******************************************************************************
Function Prototype:
void hilbert (int n, float x[], float y[]);

******************************************************************************
Input:
n		length of x and y
x		array[n] to be Hilbert transformed

Output:
y		array[n] containing Hilbert transform of x

******************************************************************************
Notes:
The Hilbert transform is computed by convolving x with a
windowed (approximate) version of the ideal Hilbert transformer.

******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/

#define LHHALF 30	/* half-length of Hilbert transform filter*/
#define LH 2*LHHALF+1	/* filter length must be odd */

/*****************************************************************************
Compute Hilbert transform y of x
******************************************************************************
Input:
n		length of x and y
x		array[n] to be Hilbert transformed

Output:
y		array[n] containing Hilbert transform of x
******************************************************************************
Notes:
The Hilbert transform is computed by convolving x with a
windowed (approximate) version of the ideal Hilbert transformer.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/
int
hilbert (double *x, double *y, int npts)
{
  static int madeh=0;
  static double h[LH];
  int i;
  double taper;
  
  if ( ! x || ! y || npts < 1 )
    {
      fprintf (stderr, "hilbert(): No data buffer specified\n");
      return -1;
    }
  
  /* If not made, make Hilbert transform filter; use Hamming window */
  if ( ! madeh )
    {
      h[LHHALF] = 0.0;
      for (i=1; i<=LHHALF; i++)
	{
	  taper = 0.54 + 0.46 * cos(PI*(double)i/(double)(LHHALF));
	  h[LHHALF+i] = taper * (-(double)(i%2)*2.0/(PI*(double)(i)));
	  h[LHHALF-i] = -h[LHHALF+i];
	}
      madeh = 1;
    }
  
  /* convolve Hilbert transform with input array */
  conv (LH,-LHHALF,h,npts,0,x,npts,0,y);
  
  return npts;
}  /* End of hilbert() */


/* Originally from convolution.c */
/* Copyright (c) Colorado School of Mines, 2010.*/
/* All rights reserved.                       */

/*********************** self documentation **********************/
/*****************************************************************************
CONVOLUTION - Compute z = x convolved with y

conv	compute the convolution of two input vector arrays

******************************************************************************
Input:
lx		length of x array
ifx		sample index of first x
x		array[lx] to be convolved with y
ly		length of y array
ify		sample index of first y
y		array[ly] with which x is to be convolved
lz		length of z array
ifz		sample index of first z

Output:
z		array[lz] containing x convolved with y

******************************************************************************
Function Prototype:
void conv (int lx, int ifx, float *x, int ly, int ify, float *y,
	int lz, int ifz, float *z);

******************************************************************************
Notes:
The operation z = x convolved with y is defined to be
           ifx+lx-1
    z[i] =   sum    x[j]*y[i-j]  ;  i = ifz,...,ifz+lz-1
            j=ifx
The x samples are contained in x[0], x[1], ..., x[lx-1]; likewise for
the y and z samples.  The sample indices of the first x, y, and z values
determine the location of the origin for each array.  For example, if
z is to be a weighted average of the nearest 5 samples of y, one might
use 
	...
	x[0] = x[1] = x[2] = x[3] = x[4] = 1.0/5.0;
	conv(5,-2,x,lx,0,y,ly,0,z);
	...
In this example, the filter x is symmetric, with index of first sample = -2.

This function is optimized for architectures that can simultaneously perform
a multiply, add, and one load from memory; e.g., the IBM RISC System/6000.
Because, for each value of i, it accumulates the convolution sum z[i] in a
scalar, this function is not likely to be optimal for vector architectures.

******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 11/23/91
*****************************************************************************/
/**************** end self doc ********************************/

/*****************************************************************************
Compute z = x convolved with y; i.e.,

           ifx+lx-1
    z[i] =   sum    x[j]*y[i-j]  ;  i = ifz,...,ifz+lz-1
            j=ifx
******************************************************************************
Input:
lx		length of x array
ifx		sample index of first x
x		array[lx] to be convolved with y
ly		length of y array
ify		sample index of first y
y		array[ly] with which x is to be convolved
lz		length of z array
ifz		sample index of first z

Output:
z		array[lz] containing x convolved with y
******************************************************************************
Notes:
The x samples are contained in x[0], x[1], ..., x[lx-1]; likewise for
the y and z samples.  The sample indices of the first x, y, and z values
determine the location of the origin for each array.  For example, if
z is to be a weighted average of the nearest 5 samples of y, one might
use 
	...
	x[0] = x[1] = x[2] = x[3] = x[4] = 1.0/5.0;
	conv(5,-2,x,lx,0,y,ly,0,z);
	...
In this example, the filter x is symmetric, with index of first sample = -2.

******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 11/23/91
*****************************************************************************/
void conv (int lx, int ifx, double *x,
	   int ly, int ify, double *y,
	   int lz, int ifz, double *z)
{
  int ilx=ifx+lx-1;
  int ily=ify+ly-1;
  int ilz=ifz+lz-1;
  int i,j,ilow,ihigh,jlow,jhigh;
  double sa,sb,xa,xb,ya,yb,*t;
  
  /* if x is longer than y, swap x and y */
  if ( lx > ly )
    {
      i = ifx;  ifx = ify;  ify = i;
      i = ilx;  ilx = ily;  ily = i;
      i = lx;  lx = ly;  ly = i;
      t = x;  x = y;  y = t;
    }
  
  /* adjust pointers for indices of first samples */
  x -= ifx;
  y -= ify;
  z -= ifz;
  
  /* OFF LEFT:  i < ify+ifx */
  
  /* zero output for all i */
  ilow = ifz;
  ihigh = ify+ifx-1;
  if ( ihigh > ilz ) ihigh = ilz;
  for ( i=ilow; i<=ihigh; ++i)
    z[i] = 0.0;
  
  /* ROLLING ON:  ify+ifx <= i < ify+ilx */
	
  /* if necessary, do one i so that number of j in overlap is odd */
  if ( i < (ify+ilx) && i <= ilz )
    {
      jlow = ifx;
      jhigh = i-ify;
      if ( (jhigh-jlow)%2 )
	{
	  sa = 0.0;
	  for (j=jlow; j<=jhigh; ++j)
	    sa += x[j]*y[i-j];
	  z[i++] = sa;
	}
    }
	
  /* loop over pairs of i and j */
  ilow = i;
  ihigh = ilx+ify-1;
  if ( ihigh > ilz) ihigh = ilz;
  jlow = ifx;
  jhigh = ilow-ify;
  for (i=ilow; i < ihigh; i+=2,jhigh+=2)
    {
      sa = sb = 0.0;
      xb = x[jhigh+1];
      yb = 0.0;
      for (j=jhigh; j >= jlow; j-=2)
	{
	  sa += xb*yb;
	  ya = y[i-j];
	  sb += xb*ya;
	  xa = x[j];
	  sa += xa*ya;
	  yb = y[i+1-j];
	  sb += xa*yb;
	  xb = x[j-1];
	}
      z[i] = sa;
      z[i+1] = sb;
    }
  
  /* if number of i is odd */
  if ( i == ihigh )
    {
      jlow = ifx;
      jhigh = i-ify;
      sa = 0.0;
      for (j=jlow; j <= jhigh; ++j)
	sa += x[j]*y[i-j];
      z[i++] = sa;
    }
  
  /* MIDDLE:  ify+ilx <= i <= ily+ifx */
  
  /* determine limits for i and j */
  ilow = i;
  ihigh = ily+ifx;
  if ( ihigh > ilz ) ihigh = ilz;
  jlow = ifx;
  jhigh = ilx;
  
  /* if number of j is even, do j in pairs with no leftover */
  if ( (jhigh-jlow)%2 )
    {
      for (i=ilow; i < ihigh; i+=2)
	{
	  sa = sb = 0.0;
	  yb = y[i+1-jlow];
	  xa = x[jlow];
	  for (j=jlow; j < jhigh; j+=2)
	    {
	      sb += xa*yb;
	      ya = y[i-j];
	      sa += xa*ya;
	      xb = x[j+1];
	      sb += xb*ya;
	      yb = y[i-1-j];
	      sa += xb*yb;
	      xa = x[j+2];
	    }
	  z[i] = sa;
	  z[i+1] = sb;
	}
      
      /* else, number of j is odd, so do j in pairs with leftover */
    }
  else
    {
      for (i=ilow; i < ihigh; i+=2)
	{
	  sa = sb = 0.0;
	  yb = y[i+1-jlow];
	  xa = x[jlow];
	  for (j=jlow; j<jhigh; j+=2) {
	    sb += xa*yb;
	    ya = y[i-j];
	    sa += xa*ya;
	    xb = x[j+1];
	    sb += xb*ya;
	    yb = y[i-1-j];
	    sa += xb*yb;
	    xa = x[j+2];
	  }
	  z[i] = sa+x[jhigh]*y[i-jhigh];
	  z[i+1] = sb+x[jhigh]*y[i+1-jhigh];
	}
    }
  
  /* if number of i is odd */
  if ( i == ihigh )
    {
      sa = 0.0;
      for (j=jlow; j<=jhigh; ++j)
	sa += x[j]*y[i-j];
      z[i++] = sa;
    }
  
  /* ROLLING OFF:  ily+ifx < i <= ily+ilx */
  
  /* if necessary, do one i so that number of j in overlap is even */
  if (i <= (ily+ilx) && i <= ilz)
    {
      jlow = i-ily;
      jhigh = ilx;
      if ( !((jhigh-jlow)%2) )
	{
	  sa = 0.0;
	  for (j=jlow; j <= jhigh; ++j)
	    sa += x[j]*y[i-j];
	  z[i++] = sa;
	}
    }
  
  /* number of j is now even, so loop over both i and j in pairs */
  ilow = i;
  ihigh = ily+ilx;
  if ( ihigh > ilz ) ihigh = ilz;
  jlow = ilow-ily;
  jhigh = ilx-2; /* Dave's new patch */
  for (i=ilow; i < ihigh; i+=2,jlow+=2)
    {
      sa = sb = 0.0;
      xa = x[jlow];
      yb = 0.0;
      for (j=jlow; j<jhigh; j+=2) {
	sb += xa*yb;
	ya = y[i-j];
	sa += xa*ya;
	xb = x[j+1];
	sb += xb*ya;
	yb = y[i-1-j];
	sa += xb*yb;
	xa = x[j+2];
      }
      sb += xa*yb;
      ya = y[i-j];
      sa += xa*ya;
      xb = x[j+1];
      sb += xb*ya;
      yb = y[i-1-j];
      sa += xb*yb;
      z[i] = sa;
      z[i+1] = sb;
    }
  
  /* if number of i is odd */
  if ( i == ihigh )
    {
      jlow = i-ily;
      jhigh = ilx;
      sa = 0.0;
      for (j=jlow; j <= jhigh; ++j)
	sa += x[j]*y[i-j];
      z[i++] = sa;
    }
  
  /* OFF RIGHT:  ily+ilx < i */
  
  /* zero output for all i */
  ilow = i;
  ihigh = ilz;
  for (i=ilow; i <= ihigh; ++i)
    z[i] = 0.0;
  
} /* End of conv() */
