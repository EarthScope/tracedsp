/**********************************************************************
 * iirfilter.c:
 *
 * Routines to filter time-series data using a Butterworth IIR filter
 * derived from low and/or high filter orders and cutoff frequencies.
 *
 * The iirfilter() routine is the entry point and only routine exposed
 * in the iirfilter.h include.
 *
 * Originally taken from the PQLII source code from PASSCAL with the
 * following attribution: written by jcf feb 1988; and heavily
 * modified after that.
 *
 * Chad Trabant - IRIS Data Management Center
 * modified: 2011.151
 **********************************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

/* Pi, with comical precision */
#define PI 3.1415926535897932384626433

typedef struct {
  double real;
  double imag;
} iircomplex;

iircomplex add_c (iircomplex u, iircomplex v);
iircomplex mul_c (iircomplex u, iircomplex v);
iircomplex cmul_c (double a, iircomplex u);
iircomplex sub_c (iircomplex u, iircomplex v);
iircomplex div_c (iircomplex u, iircomplex v);
iircomplex conj_c (iircomplex u);

void applyfilter (double a1, double a2, double b1, double b2,
		  int npts, int reverse, double *input, double *output);
void lowpass (double fc, double dt, int n, iircomplex *p, double *b);
void highpass (double fc, double dt, int n, iircomplex *p, double *b);


/**********************************************************************
 * iirfilter: filter a timeseries with a derived IIR definition
 *
 * Filter a timeseries using a Butterworth IIR digital filter determined
 * from high and low filter orders and cutoff frequencies.  Either the 
 * high or low order filter is optional or both may be used for a bandpass
 * filter.
 *
 * To avoid end effects the mean is removed from the series prior to
 * filtering, after the filter operation is complete the mean value is
 * added back to each sample.
 *
 * The *output data array will be allocated if necessary and it is up
 * to the calling routine to free that memory.  If *output is non-zero
 * it will be re-allocated for reuse.
 *
 * input      : Pointer to input data array
 * inputtype  : Input sample type ('i'=INT32, 'f'=FLOAT32, 'd'=FLOAT64)
 * samplecnt  : Number of samples in input data array
 * reverse    : Apply 2-pass filter, one forward, one reverse
 * output     : Pointer to pointer to output data array
 * outputtype : Output sample type ('i'=INT32, 'f'=FLOAT32, 'd'=FLOAT64)
 * highorder  : Order of high-pass filter, must be even, 0 to disable high-pass filter
 * highcutoff : Cutoff frequency of high-pass filter
 * loworder   : Order of low-pass filter, must be even, 0 to disable low-pass filter
 * lowcutoff  : Cutoff frequency of low-pass filter
 * samprate   : Sampling rate in Hz
 * verbose    : Verbosity level
 *
 * Returns 0 on success and 1 on error.
 **********************************************************************/
int
iirfilter ( void *input, char inputtype, int samplecnt, int reverse,
	    void **output, char outputtype,
	    int highorder, double highcutoff,
	    int loworder, double lowcutoff,
	    double samprate, int verbose )
{
  int         *iptr;
  float       *fptr;
  double      *dptr;
  int          idx;
  int          insampsize;        /* Input sample size in bytes */
  int          outsampsize;       /* Output sample size in bytes */
  iircomplex   pl[12], ph[12];    /* Filter poles */
  double       b0l, b0h;          /* Filter gains */
  double      *doutput;           /* Output samples as doubles */
  double       a1, a2, b1, b2;
  double       mean;
  
  /* Sanity check of input parameters */
  if ( ! input || ! output )
    return 1;
  
  /* Determine sample sizes */
  switch ( inputtype )
    {
    case 'i': insampsize = 4; break;
    case 'f': insampsize = 4; break;
    case 'd': insampsize = 8; break;
    default:
      fprintf (stderr, "[iirfilter] Unknown input data type: %d\n", inputtype);
      return 1;
    }
  switch ( outputtype )
    {
    case 'i': outsampsize = 4; break;
    case 'f': outsampsize = 4; break;
    case 'd': outsampsize = 8; break;
    default:
      fprintf (stderr, "[iirfilter] Unknown outputput data type: %d\n", outputtype);
      return 1;
    }
  
  /* Verify high-pass filter parameters */
  if ( (highorder%2) != 0 || highorder < 0 || highorder > 13 )
    {
      fprintf (stderr, "[iirfilter] Order of high-pass filter cannot be %d\n", highorder);
      fprintf (stderr, "[iirfilter]   Order must be even and less than 13\n");
      return 1;
    }
  else if ( highorder )
    {
      if ( highcutoff == 0.0 )
	{
	  fprintf (stderr, "[iirfilter] High-pass cutoff frequency cannot be 0\n");
	  return 1;
	}
      
      if ( highcutoff > samprate/2.0 ) 
	{
	  fprintf (stderr, "[iirfilter] High-pass cutoff frequency cannot be higher than Nyquist (%f)\n", samprate/2.0);
	  return 1;
	}
    }
  
  /* Verify low-pass filter parameters */
  if ( (loworder%2) != 0 || loworder < 0 || loworder > 13 )
    {
      fprintf (stderr, "[iirfilter] Order of low-pass filter cannot be %d\n", loworder);
      fprintf (stderr, "[iirfilter]   Order must be even and less than 13\n");
      return 1;
    }
  else if ( loworder )
    {
      if ( lowcutoff == 0.0 )
	{
	  fprintf (stderr, "[iirfilter] Low-pass cutoff frequency cannot be 0\n");
	  return 1;
	}
      
      if ( lowcutoff > samprate/2.0 ) 
	{
	  fprintf (stderr, "[iirfilter] Low-pass cutoff frequency cannot be higher than Nyquist (%f)\n", samprate/2.0);
	  return 1;
	}
      
      if ( highorder && lowcutoff <= highcutoff ) 
	{
	  fprintf (stderr, "[iirfilter] Low-pass cutoff frequency cannot be <= high-pass cutoff\n");
	  return 1;
	}
    }
  
  /* Verify sampling rate */
  if ( samprate <= 0.0 )
    {
      fprintf (stderr, "[iirfilter] Sampling rate must be greater than 0.0 (%f)\n", samprate);
      return 1;
    }
  
  /* Get high-pass filter poles if necessary */
  if ( highorder )
    {
      highpass (highcutoff, 1.0/samprate, highorder, ph, &b0h);
      
      if ( verbose > 1 )
	fprintf (stderr, "[iirfilter] Gain of high-pass filter: %f\n", b0h);
    }
  
  /* Get low-pass filter poles if necessary */
  if ( loworder )
    {
      lowpass (lowcutoff, 1.0/samprate, loworder, pl, &b0l);
      
      if ( verbose > 1 )
	fprintf (stderr, "[iirfilter] Gain of low-pass filter: %f\n", b0l);
    }
  
  /* Allocate internal filter output buffer */
  if ( (doutput = (double *) malloc (sizeof(double) * samplecnt)) == NULL )
    {
      fprintf (stderr, "[iirfilter] Error allocating filter buffer\n");
      return 1;
    }
  
  /* Demean the input data while copying to internal filter buffer (doutput) */
  switch ( inputtype )
    {
    case 'i':
      iptr = (int *) input;
      
      mean = (double) (*iptr);
      for (idx=1; idx < samplecnt; idx++)
	mean = mean + ( *(iptr+idx) - mean) / (idx + 1);
      
      if ( verbose > 1 )
	fprintf (stderr, "[iirfilter] Removing mean value of %f\n", mean);
      
      for (idx=0; idx < samplecnt; idx++)
        *(doutput+idx) = *(iptr+idx) - mean;
      
      break;
      
    case 'f':
      fptr = (float *) input;
      
      mean = (double) (*fptr);
      for (idx=1; idx < samplecnt; idx++)
	mean = mean + ( *(fptr+idx) - mean) / (idx + 1);
      
      if ( verbose > 1 )
	fprintf (stderr, "[iirfilter] Removing mean value of %f\n", mean);
      
      for (idx=0; idx < samplecnt; idx++)
        *(doutput+idx) = *(fptr+idx) - mean;
      
      break;
      
    case 'd':
      dptr = (double *) input;
      
      mean = *dptr;
      for (idx=1; idx < samplecnt; idx++)
	mean = mean + ( *(dptr+idx) - mean) / (idx + 1);
      
      if ( verbose > 1 )
	fprintf (stderr, "[iirfilter] Removing mean value of %f\n", mean);
      
      for (idx=0; idx < samplecnt; idx++)
        *(doutput+idx) = *(dptr+idx) - mean;
      
      break;
    }
  
  /* Filter the data: filtering is implemented as a cascade of second order filters */

  /*  Apply high-pass filter using poles ph
   *    Numerator polynomial is z**2 - 2*z + 1 */
  if ( highorder )
    {
      for ( idx=0 ; idx < highorder ; idx += 2 )
	{
	  a1 = -2*ph[idx].real;
	  a2 = ph[idx].real*ph[idx].real + ph[idx].imag*ph[idx].imag;
	  b1 = -2;
	  b2 = 1;
	  
	  applyfilter (a1, a2, b1, b2, samplecnt, 0, doutput, doutput);
	}
      
      /* Compensate for filter gain */
      for ( idx=0 ; idx < samplecnt ; idx++ )
	doutput[idx] = b0h * doutput[idx];
      
      /* Apply filter in reverse order to remove phase distortion */
      if ( reverse )
	{
	  for ( idx=0 ; idx < highorder ; idx += 2 )
	    {
	      a1 = -2*ph[idx].real;
	      a2 = ph[idx].real*ph[idx].real + ph[idx].imag*ph[idx].imag;
	      b1 = -2;
	      b2 = 1;
	      
	      applyfilter (a1, a2, b1, b2, samplecnt, 1, doutput, doutput);
	    }
	  
	  /* Compensate for filter gain */
	  for ( idx=0 ; idx < samplecnt ; idx++ )
	    doutput[idx] = b0h * doutput[idx];
	}
    }
  
  /* Apply low-pass filter using poles pl
   *   Numerator polynomial is z**2 + 2*z + 1 */
  if ( loworder )
    {
      for ( idx=0 ; idx < loworder ; idx += 2 )
	{
	  a1 = -2*pl[idx].real;
	  a2 = pl[idx].real*pl[idx].real + pl[idx].imag*pl[idx].imag;
	  b1 = 2;
	  b2 = 1;
	  
	  applyfilter (a1, a2, b1, b2, samplecnt, 0, doutput, doutput);
	}
      
      /* Compensate for filter gain */
      for ( idx=0 ; idx < samplecnt ; idx++ )
	doutput[idx] = b0l * doutput[idx];
      
      /* Apply filter in reverse order to remove phase distortion */
      if ( reverse )
	{
	  for ( idx=0 ; idx < loworder ; idx += 2 )
	    {
	      a1 = -2*pl[idx].real;
	      a2 = pl[idx].real*pl[idx].real + pl[idx].imag*pl[idx].imag;
	      b1 = 2;
	      b2 = 1;
	      
	      applyfilter (a1, a2, b1, b2, samplecnt, 1, doutput, doutput);
	    }
	  
	  /* Compensate for filter gain */
	  for ( idx=0 ; idx < samplecnt ; idx++ )
	    doutput[idx] = b0l * doutput[idx];
	}
    }
  
  /* (Re)Allocate output buffer */
  if ( ! *output )
    {
      if ( verbose > 1 )
	fprintf (stderr, "[iirfilter] Allocating output data buffer\n");
      
      if ( (*output = (void *) malloc (outsampsize * samplecnt)) == NULL )
	{
	  fprintf (stderr, "[iirfilter] Error allocating output buffer\n");
	  return 1;
	}
    }
  else
    {
      if ( verbose > 1 )
	fprintf (stderr, "[iirfilter] Reallocating output data buffer\n");
      
      if ( (*output = (void *) realloc (*output, outsampsize * samplecnt)) == NULL )
	{
	  fprintf (stderr, "[iirfilter] Error allocating output buffer\n");
	  return 1;
	}
    }

  /* Copy samples into output buffer rounding integer data types, no test for overflow */
  /* Add mean back into sample values while copying */
  switch ( outputtype )
    {
    case 'i':
      for ( idx=0; idx < samplecnt; idx++ ) 
	{
	  *((int *)*output+idx) = (int) (*(doutput+idx) + mean + 0.5);
	}
      break;
    case 'f':
      for ( idx=0; idx < samplecnt; idx++ )
	{
	  *((float *)*output+idx) = (float) (*(doutput+idx) + mean);
	}
      break;
    case 'd':
      for ( idx=0; idx < samplecnt; idx++ ) 
	{
	  *((double *)*output+idx) = (double) (*(doutput+idx) + mean);
	}
      break;
    }
  
  if ( doutput )
    free (doutput);
  
  return 0;
}  /* End of iirfilter() */


/**********************************************************************
 * add_c: Add two iircomplex numbers: u + v
 **********************************************************************/
iircomplex
add_c (iircomplex u, iircomplex v)
{
  iircomplex w;
  
  w.real = u.real + v.real;
  w.imag = u.imag + v.imag;
  
  return (w);
}


/**********************************************************************
 * mul_c: Multiply two iircomplex numbers: u * v
 **********************************************************************/
iircomplex
mul_c (iircomplex u, iircomplex v)
{
  iircomplex w;
  
  w.real = u.real*v.real - u.imag*v.imag;
  w.imag = u.real*v.imag + u.imag*v.real;
  
  return (w);
}


/**********************************************************************
 * cmul_c: Multiply a real number and an iircomplex number
 **********************************************************************/
iircomplex
cmul_c (double a, iircomplex u)
{
  iircomplex w;
  
  w.real = a * u.real;
  w.imag = a * u.imag;
  
  return (w);
}


/**********************************************************************
 * sub_c: Subtract two iircomplex numbers: u - v
 **********************************************************************/
iircomplex
sub_c (iircomplex u, iircomplex v)
{
  iircomplex w;
  
  w.real = u.real - v.real;
  w.imag = u.imag - v.imag;
  
  return (w);
}


/**********************************************************************
 * div_c: Divide two iircomplex numbers: u / v
 **********************************************************************/
iircomplex
div_c (iircomplex u, iircomplex v)
{
  iircomplex w;

  w.real = w.imag = 0;
  
  /* Check for divide by 0 */
  if ( v.real != 0 && v.imag != 0 )
    {
      w.real = ((u.real * v.real) + (u.imag * v.imag)) /
	((v.real * v.real) + (v.imag * v.imag));
      w.imag = ((u.imag * v.real) - (u.real * v.imag)) /
	((v.real * v.real) + (v.imag * v.imag));
    }
  else
    {
      fprintf (stderr, "ERROR: iircomplex division by 0 in div_c\n");
    }
  
  return (w);
}


/**********************************************************************
 * conj_c: Calculate the iircomplex conjugate
 **********************************************************************/
iircomplex
conj_c (iircomplex u)
{
  iircomplex w;
  
  w.real = u.real;
  w.imag = -u.imag;
  
  return (w);
}


/**********************************************************************
 * applyfilter:
 *
 * Apply a second order recursive filter to a data array.
 * Denomonator polynomial is z**2 + a1*z + a2
 * Numerator polynomial is z**2 + b1*z + b2
 *
 * input   : Array of input data
 * output  : Array of output data (may be the same as input)
 * npts    : Number of points in input/output arrays
 * reverse : Apply filter in reverse order if true
 *
 **********************************************************************/
void
applyfilter (double a1, double a2, double b1, double b2,
	     int npts, int reverse, double *input, double *output)
{
  double d1 = 0.0;
  double d2 = 0.0;
  double out;
  int i;
  
  if ( ! reverse )
    {
      for ( i=0 ; i<npts ; i++ )
	{
	  out = input[i] + d1;
	  d1 = b1*input[i] - a1*out + d2;
	  d2 = b2*input[i] - a2*out;
	  output[i] = out;
	}
    }
  else
    {
      for ( i=npts-1 ; i>=0 ; i-- )
	{
	  out = input[i] + d1;
	  d1 = b1*input[i] - a1*out + d2;
	  d2 = b2*input[i] - a2*out;
	  output[i] = out;
	}
    }
}


/**********************************************************************
 * lowpass:
 *
 * Compute low-pass filter poles for Butterworth filter
 *   fc = desired cutoff frequency
 *   dt = sample period in seconds
 *   n = number of poles (MUST BE EVEN)
 *   p = pole locations (RETURNED)
 *   b = gain factor for filter (RETURNED)
 *
 * This routine calculates a continuous Butterworth low-pass IIR with
 * specified cutoff frequency.  A discrete filter is then calculated
 * utilizing the bilinear transform.  Methods used here follow those
 * in Digital Filters and Signal Processing by Leland B. Jackson.
 *
 **********************************************************************/
void
lowpass (double fc, double dt, int n, iircomplex *p, double *b)
{  
  iircomplex one, x, y;
  double wcp, wc, b0;
  int i, i1;
  
  wcp = 2 * fc * PI;
  wc = (2.0/dt) * tan(wcp*dt/2.0);
  one.real = 1.0;
  one.imag = 0.0;
  
  /* Calculate position of poles for continuous filter */
  for ( i=0 ; i < n ; i += 2 )
    {
      i1 = i + 1;
      p[i].real = -wc * cos(i1*PI/(2*n));
      p[i].imag = wc * sin(i1*PI/(2*n));
      p[i+1] = conj_c(p[i]);
    }
  
  /* Calculate position of poles for discrete filter using
   * the bilinear transformation */
  for ( i=0 ; i < n ; i += 2 )
    {
      p[i] = cmul_c(dt/2,p[i]);
      x = add_c(one,p[i]);
      y = sub_c(one,p[i]);
      p[i] = div_c(x,y);
      p[i+1] = conj_c(p[i]);
    }
  
  /* Calculate filter gain */
  b0 = 1.0;
  for ( i=0 ; i < n ; i +=2 )
    {
      x = sub_c(one,p[i]);
      y = sub_c(one,p[i+1]);
      x = mul_c(x,y);
      b0 = b0 * 4.0 / x.real;
    }
  *b = 1.0 / b0;
}


/**********************************************************************
 * highpass:
 *
 * Compute high-pass filter poles for Butterworth filter
 *   fc = desired cutoff frequency
 *   dt = sample period in seconds
 *   n = number of poles (MUST BE EVEN)
 *   p = pole locations (RETURNED)
 *   b = gain factor for filter (RETURNED)
 *
 * This routine calculates a low-pass IIR with lowpass().  Then this
 * filter is converted to a high pass filter.  Methods used here
 * follow those in Digital Filters and Signal Processing by Leland
 * B. Jackson.
 *
 **********************************************************************/
void
highpass (double fc, double dt, int n, iircomplex *p, double *b)
{  
  iircomplex one, x, y;
  double wcp, wc,  alpha, b0;
  int i;
  
  wcp = 2 * fc * PI;
  wc = (2.0/dt) * tan(wcp*dt/2.0);
  alpha = cos(wc*dt);
  one.real = 1.0;
  one.imag = 0.0;
  
  /* Compute poles for low-pass filter  */
  lowpass (fc, dt, n, p, &b0) ;
  
  /* Now find poles for high-pass filter */
  for ( i=0 ; i < n ; i += 2 )
    {
      x = cmul_c (alpha,one) ;
      x = sub_c (x,p[i]) ;
      y = cmul_c (alpha,p[i]) ;
      y = sub_c(one,y) ;
      p[i] = div_c(x,y) ;
      p[i+1] = conj_c(p[i]) ;
    }
  
  /* Calculate gain for high-pass filter */
  b0 = 1.0;
  for ( i=0 ; i < n ; i += 2 )
    {
      x = add_c(one,p[i]);
      y = add_c(one,p[i+1]);
      x = mul_c(x,y);
      b0 = b0 * 4.0 / x.real;
    }
  *b = 1.0 / b0;
}
