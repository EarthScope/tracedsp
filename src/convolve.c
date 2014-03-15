/*********************************************************************
 * convolve.c
 *
 * Time-series convolution routines, the primary entry points are the
 * convolve_sac() and convolve_resp() routines which (de)convolve a
 * time-series with the frequency response described either in a SAC
 * formatted poles and zeros file or a SEED RESP file.
 *
 * The (de)convolution scheme used herein is more or less standard
 * with the following parameters:
 *
 * number of points in the FFT (nfft): next power of 2 of # samples
 * number of frequencies (nfreqs):  nfft/2 + 1
 * frequencies from 0 (DC) to Nyquist at intervals of samprate/nfft.
 *
 * Modified: 2014.073
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <float.h>

#include <evresp.h>

#include "convolve.h"
#include "whiten.h"
#include "getpzfr.h"

#define PI 3.1415926535897932384626433832795

/*********************************************************************
 * convolve:
 *
 * Convole a time-series with the frequency response described by a
 * complex frequency response.
 *
 * This is a modified version of SAC 2000's transfer() routine.
 *
 * Arguments:
 *   data       : array of input data samples
 *   npts       : number of samples
 *   delta      : sampling period in seconds
 *   nfreqs     : number of unique points in FFT (FFT mag. is symmetric)
 *   nfft       : number of all points used in FFT (must be 2^x)
 *   creal      : real components of FR to convolve (nfreqs)
 *   cimag      : imaginary components of FR to convolve (nfreqs)
 *   dreal      : real components of FR to deconvolve (nfreqs)
 *   dimag      : imaginary components of FR to deconvolve (nfreqs)
 *   taperfreq  : spectrum taper filter definition
 *                  f0,f1 = frequency range for high-pass taper
 *                  f2,f3 = frequency range for low-pass taper
 *   prewhiten  : order of predictive filter to prewhiten the data,
 *                a value of 0 (not the pointer itself) indicates none
 *
 * Returns 0 on success and -1 on error.
 *********************************************************************/
int
convolve (double data[], int npts, double delta, int nfreqs, int nfft,
	  double *creal, double *cimag, double *dreal, double *dimag,
	  double taperfreq[], int *prewhiten,
	  int verbose)
{
  int i, j;
  int wlwarning = 0;
  double whitecoef[21];
  double delfreq, denr, fac, freq;
  
  double *sreal, *simag;
  
  delfreq = 1.0 / (delta * nfft);
  
  if ( ! creal && ! cimag && ! dreal && ! dimag )
    {
      fprintf (stderr, "convolve(): No frequency responses supplied\n");
      return -1;
    }
  
  if ( verbose > 2 )
    {
      fprintf (stderr, "Convolution parameters:\n");
      fprintf (stderr, "  npts: %d\n", npts);
      fprintf (stderr, "  nfft: %d\n", nfft);
      fprintf (stderr, "  nfreqs: %d\n", nfreqs);
      fprintf (stderr, "  delfreq: %g\n", delfreq);
      if ( taperfreq )
	fprintf (stderr, "  Spectral tapering: %g/%g => %g/%g\n",
		 taperfreq[0], taperfreq[1], taperfreq[2], taperfreq[3]);
      else
	fprintf (stderr, "  Spectral tapering: none\n");
    }
  
  /* Allocate complex number arrays */
  if ( (sreal = (double *) malloc (2*nfft*sizeof(double))) == NULL )
    {
      fprintf (stderr, "convolve(): Error allocating memory\n");
      return -1;
    }
  if ( (simag = (double *) malloc (2*nfft*sizeof(double))) == NULL )
    {
      fprintf (stderr, "convolve(): Error allocating memory\n");
      return -1;
    }
  
  /* Prewhiten the data if requested */
  if ( *prewhiten > 0 )
    {
      if ( verbose > 1 )
	fprintf (stderr, "Prewhitening the data\n");
      
      if ( whiten (data, npts, prewhiten, whitecoef) )
	{
	  fprintf (stderr, "convolve(): Error prewhitening the data\n");
	  return -1;
	}
    }
  
  /* For deconvolution need the inverse response function:
   * Compute 1 / Function, applying waterlevel of FLT_MIN */
  if ( dreal && dimag )
    {
      for ( i = 0; i < nfreqs; i++ )
	{
	  denr = (dreal[i]*dreal[i] + dimag[i]*dimag[i]);
	  
	  if ( denr <= FLT_MIN )
	    {
	      dreal[i] = 0.0;
	      dimag[i] = 0.0;
	      wlwarning++;
	    }
	  else
	    {
	      denr = 1.0 / denr;
	      dreal[i] = dreal[i] * denr;
	      dimag[i] = -dimag[i] * denr;
	    }
	}
      
      if ( verbose && wlwarning > 0 )
	fprintf (stderr, "WARNING: Water level reached for %d coeffients during response inversion\n",
		 wlwarning);
    }
  
  /* Determine if performing convolution, deconvolution or both (aka transfer).
     Calculate composite if needed, target response will be in creal & cimag */
  if ( creal && cimag && dreal && dimag )
    {
      /* Multiply response functions to get composite */
      double temp;
      for ( i = 0; i < nfreqs; i++ )
	{
	  temp = creal[i]*dreal[i] - cimag[i]*dimag[i];
	  cimag[i] = creal[i]*dimag[i] + cimag[i]*dreal[i];
	  creal[i] = temp;
	}
    }
  else if ( dreal && dimag )
    {
      /* If only deconvolution, redirect the creal & cimag pointers */
      creal = dreal;
      cimag = dimag;
    }
  
  /* Print tapered composite response */
  if ( verbose > 3 )
    {
      double amp;
      double phase;
      double real;
      double imag;
      
      fprintf (stderr, "Composite convolution operator (FAP):\n");
      for ( i = 0; i < nfreqs; i++ )
	{
	  freq = i * delfreq;
	  
	  if ( taperfreq )
	    fac = (spectraltaper (freq, taperfreq[1], taperfreq[0]) *
		   spectraltaper (freq, taperfreq[2], taperfreq[3]));
	  else
	    fac = 1.0;
	  
	  real = creal[i] * fac;
	  imag = cimag[i] * fac;
	  
	  amp = sqrt (real*real + imag*imag);
	  phase = atan2 (imag, real + 1.0e-200) * 180/PI;
	  
	  fprintf (stderr, "%.6E  %.6E  %.6E\n", freq, amp, phase);
	}
    }
  
  /* Scale and optionally taper the frequency response function */
  for ( i = 0; i < nfreqs; i++ )
    {
      freq = i * delfreq;
      
      if ( taperfreq )
	fac = delfreq * (spectraltaper (freq, taperfreq[1], taperfreq[0]) *
			 spectraltaper (freq, taperfreq[2], taperfreq[3]));
      else
	fac = delfreq;
      
      creal[i] *= fac;
      cimag[i] *= fac;
    }
  
  /* Fill a complex vector with data samples, zero-padding as needed */
  for ( i = 0; i < npts; i++ )
    {
      sreal[i] = data[i] * delta;
      simag[i] = 0.0;
    }
  for ( i = npts; i < nfft; i++ )
    {
      sreal[i] = 0.0;
      simag[i] = 0.0;
    }
  
  /* Transform the data */
  if ( verbose > 1 )
    fprintf (stderr, "Calculating FFT\n");
  fft ( sreal, simag, nfft, 1 );
  
  /* Multiply transformed data by frequency response operator */
  for ( i = 0; i < nfreqs; ++i )
    {
      double tempR = creal[i]*sreal[i] - cimag[i]*simag[i];
      double tempI = creal[i]*simag[i] + cimag[i]*sreal[i];
      sreal[i] = tempR;
      simag[i] = tempI;
      /* Input data are real so F(N-j) = (F(j))*   */
      if ( i > 0 && i < (nfreqs-1) )
	{
	  j = nfft - i;
	  sreal[j] = tempR;
	  simag[j] = -tempI;
	}
    }
  
  /* Perform the inverse transform */
  if ( verbose > 1 )
    fprintf (stderr, "Calculating inverse FFT\n");
  fft ( sreal, simag, nfft, -1 );
  
  /* Copy the transformed data back into the original data array */
  for ( i = 0; i < npts; ++i )
    data[i] = sreal[i];
  
  /* Undo the effects of prewhitening if necessary */
  if ( *prewhiten > 0 )
    {
      if ( verbose > 1 )
	fprintf (stderr, "Dewhitening the data\n");
      
      if ( dewhiten (data, npts, *prewhiten, whitecoef) )
	{
	  fprintf (stderr, "convolve(): Error dewhitening data\n");
	  return -1;
	}
    }
  
  if ( sreal )
    free (sreal);
  if ( simag )
    free (simag);
  
  return 0;
}  /* End of convolve() */


/*********************************************************************
 * calcfr_sac:
 *
 * Calculate the frequency response described by a set of poles and
 * zeros.  The poles, zeros and constant must be in the same format
 * used by SAC 2000.
 *
 * Arguments:
 *   nfreqs     : number of frequencies
 *   delfreq    : frequency step
 *   sacpzfilename : name of file containing poles and zeros in SAC format
 *   xreal      : pointer to array for FR real values, (re)allocated
 *   ximag      : pointer to array for FR imaginary values, (re)allocated
 *   verbose    : controls level of diagnostic output
 *
 * The frequency response arrays xreal and ximag will be (re)allocated
 * and it is up to the caller to free these arrays after use.
 *
 * Returns 0 on success and -1 on error.
 *********************************************************************/
int
calcfr_sac (int nfreqs, double delfreq, char *sacpzfilename,
	    double **xreal, double **ximag, int verbose)
{
  /* Sanity checks */
  if ( ! sacpzfilename || ! xreal || ! ximag )
    return -1;
  
  /* Allocate complex number arrays */
  if ( (*xreal = (double *) realloc (*xreal, 2*nfreqs*sizeof(double))) == NULL )
    {
      fprintf (stderr, "calcfr_sac(): Error allocating memory\n");
      return -1;
    }
  if ( (*ximag = (double *) realloc (*ximag, 2*nfreqs*sizeof(double))) == NULL )
    {
      fprintf (stderr, "calcfr_sac(): Error allocating memory\n");
      return -1;
    }
  
  /* Read poles & zeros and calculate frequency response */
  if ( verbose > 1 )
    fprintf (stderr, "Determinig frequency response from SAC Poles & Zeros\n");
  
  if ( getpzfr (nfreqs, delfreq, *xreal, *ximag, sacpzfilename) )
    {
      fprintf (stderr, "calcfr_sac(): Error determining frequency response for: %s\n",
	       sacpzfilename);
      return -1;
    }
  
  return 0;
}  /* End of calcfr_sac() */


/*********************************************************************
 * calcfr_resp:
 *
 * Calculate the frequency response described by a SEED RESP file.
 * The complex response is evaluated using evalresp.
 *
 * Arguments:
 *   nfreqs     : number of frequencies
 *   delfreq    : frequency step
 *   net        : SEED network code to match, can be '*'
 *   sta        : SEED station code to match, can be '*'
 *   loc        : SEED location id to match, can be '*'
 *   chan       : SEED channel code to match, can be '*'
 *   startstage : starting stage, can be -1 for first
 *   stopstage  : stoping stage, can be -1 for last
 *   units      : response units, can be: DIS, VEL, ACC or DEF
 *   resptime   : time of response to match as epoch time
 *   usedelay   : flag to control usage of estimated delay in RESP
 *   respfilename : name of RESP file
 *   totalsensflag : use total sensitivity intead of gain product
 *   xreal      : pointer to array for FR real values, (re)allocated
 *   ximag      : pointer to array for FR imaginary values, (re)allocated
 *   verbose    : controls level of diagnostic output
 *
 * The frequency response arrays xreal and ximag will be (re)allocated
 * and it is up to the caller to free these arrays after use.
 *
 * Returns 0 on success and -1 on error.
 *********************************************************************/
int
calcfr_resp (int nfreqs, double delfreq, char *net, char *sta,
	     char *loc, char *chan, int startstage, int stopstage,
	     char *units, time_t resptime, int usedelay,
	     char *respfilename, int totalsensflag, 
	     double **xreal, double **ximag, int verbose)
{
  int i;
  struct tm *rtime;
  char datetime[20];
  struct response *resp = NULL;
  double *freqs = NULL;
  
  /* Sanity checks */
  if ( ! respfilename || ! xreal || ! ximag )
    return -1;
  
  /* Allocate complex number arrays */
  if ( (*xreal = (double *) realloc (*xreal, 2*nfreqs*sizeof(double))) == NULL )
    {
      fprintf (stderr, "calcfr_resp(): Error allocating memory\n");
      return -1;
    }
  if ( (*ximag = (double *) realloc (*ximag, 2*nfreqs*sizeof(double))) == NULL )
    {
      fprintf (stderr, "calcfr_resp(): Error allocating memory\n");
      return -1;
    }
  
  /* Construct an appropriate time string: YYYY,JJJ,HH:MM:SS */
  rtime = gmtime(&resptime);
  snprintf (datetime, sizeof(datetime), "%04d,%03d,%02d:%02d:%02d",
	    rtime->tm_year+1900, rtime->tm_yday+1,
	    rtime->tm_hour, rtime->tm_min, rtime->tm_sec);
  
  /* Allocate frequencies array */
  if ( (freqs = (double *) malloc(nfreqs * sizeof(double))) == NULL )
    {
      fprintf (stderr, "calcfr_resp(): Cannot allocate memory\n");
      return -1;
    }
  
  /* Calculate frequencies */
  for (i=0 ; i < nfreqs ; i++)
    {
      freqs[i] = i * delfreq;
    }
  
  /* Control usage of any delay specified in the SEED response information */
  if ( usedelay )
    use_delay (TRUE);
  else
    use_delay (FALSE);
  
  if ( verbose > 1 )
    fprintf (stderr, "Determinig frequency response from SEED RESP\n");
  
  /* evresp (evalresp) arguments in order:
   * char *stalst = List of stations (space or comma separated), glob matching
   * char *chalst = List of channels (space or comma separated), glob matching
   * char *net_code = Network code, glob matching
   * char *locidlst = List of locations (comma separated), glob matching
   * char *date_time = Date and Time in YYYY,DDD,HH:MM:SS format
   * char *units = Output units: DEF, DIS, VEL, ACC (DEF=no conversion)
   * char *file = Input RESP file name
   * double *freqs = Array of desired frequencies for evaluation
   * int nfreqs = Number of frequencies in freqs array
   * char *rtype = Response type ("cs" = complex spectrum), unused I think.
   * char *verbose = Verbosity string (e.g. NULL or "-v")
   * int start_stage = Starting stage for evaluation, -1 for beginning
   * int stop_stage = Ending stage for evaluation, -1 for end
   * int stdio_flag = Set to true if input is from stdin
   * int listinterp_out_flag = Interpolate List blockette output, 0=no
   * int listinterp_in_flag = Interpolate List blockette input, 0=no
   * double ilstinterp_tension = Interpolation tension
   * int useTotalSensitivityFlag = Use of total sensitivity instead of gain product
   */
  
  /* Read RESP information and calculate frequency response using evalresp */
  resp = evresp_itp (sta, chan, net, loc, datetime, units, respfilename, freqs, nfreqs,
		     "cs", (verbose>1)?"-v":NULL, startstage, stopstage, 0, 0, 0, 0, totalsensflag);
  
  if ( ! resp )
    {
      fprintf (stderr, "calcfr_resp(): Error with evresp()\n");
      return -1;
    }
  else if ( resp->nfreqs != nfreqs )
    {
      fprintf (stderr, "calcfr_resp(): Number of frequencies requested (%d) not equal to returned (%d), FAP response?\n",
	       nfreqs, resp->nfreqs);
      return -1;
    }
  else if ( verbose )
    {
      fprintf (stderr, "Found response for %s.%s.%s.%s at %s (stages %d to %d)\n",
	       resp->network, resp->station, resp->locid, resp->channel,
	       datetime, startstage, stopstage);
    }
  
  /* Split evalresp's complex response into real & imaginary arrays */
  for (i=0 ; i < nfreqs ; i++)
    {
      (*xreal)[i] = resp->rvec[i].real;
      (*ximag)[i] = resp->rvec[i].imag;
    }
  
  if ( resp )
    free_response (resp);
  if ( freqs )
    free (freqs);
  
  return 0;
}  /* End of calcfr_resp() */


/*********************************************************************
 * convolve_sac:
 *
 * Convole a time-series with the frequency response described by a
 * set of poles and zeros.  The poles, zeros and constant must be in
 * the same format used by SAC 2000.
 *
 * Arguments:
 *   data       : array of input data samples
 *   npts       : number of samples
 *   delta      : sampling period in seconds
 *   taperfreq  : spectrum taper filter definition
 *                  f0,f1 = frequency range for high-pass taper
 *                  f2,f3 = frequency range for low-pass taper
 *   prewhiten  : order of predictive filter to prewhiten the data,
 *                a value of 0 (not the pointer itself) indicates none
 *   deconvflag : 0=Convolution and 1=Deconvolution
 *   sacpzfilename : name of file containing poles and zeros in SAC format
 *   verbose    : controls level of diagnostic output
 *
 * Returns 0 on success and -1 on error.
 *********************************************************************/
int
convolve_sac (double data[], int npts, double delta, double taperfreq[], 
	      int *prewhiten, int deconvflag, char *sacpzfilename,
	      int verbose)
{
  int retval;
  int nfft;
  int nfreqs;
  double delfreq;

  double *xreal = NULL;
  double *ximag = NULL;
  
  nfft = next2 (npts);
  nfreqs = nfft/2 + 1;
  
  delfreq = 1.0 / (delta * nfft);
  
  /* Calculate frequency response from SAC P&Zs */
  if ( calcfr_sac (nfreqs, delfreq, sacpzfilename, &xreal, &ximag, verbose) )
    {
      return -1;
    }
  
  /* Call main convolution routine */
  if ( deconvflag )
    retval = convolve (data, npts, delta, nfreqs, nfft, NULL, NULL,
		       xreal, ximag, taperfreq, prewhiten, verbose);
  else
    retval = convolve (data, npts, delta, nfreqs, nfft, xreal, ximag,
		       NULL, NULL, taperfreq, prewhiten, verbose);
  
  if ( xreal )
    free (xreal);
  if ( ximag )
    free (ximag);
  
  return retval;
}  /* End of convolve_sac() */


/*********************************************************************
 * convolve_resp:
 *
 * Convole a time-series with the frequency response described by a
 * SEED RESP file.  The complex response is evaluated using evalresp.
 *
 * Arguments:
 *   data       : array of input data samples
 *   npts       : number of samples
 *   delta      : sampling period in seconds
 *   net        : SEED network code to match, can be '*'
 *   sta        : SEED station code to match, can be '*'
 *   loc        : SEED location id to match, can be '*'
 *   chan       : SEED channel code to match, can be '*'
 *   startstage : starting stage, can be -1 for first
 *   stopstage  : stoping stage, can be -1 for last
 *   units      : response units, can be: DIS, VEL, ACC or DEF
 *   resptime   : time of response to match as epoch time
 *   usedelay   : flag to control usage of estimated delay in RESP
 *   taperfreq  : spectrum taper filter definition
 *                  f0,f1 = frequency range for high-pass taper
 *                  f2,f3 = frequency range for low-pass taper
 *   prewhiten  : order of predictive filter to prewhiten the data,
 *                a value of 0 (not the pointer itself) indicates none
 *   deconvflag : 0=Convolution and 1=Deconvolution
 *   respfilename : name of RESP file
 *   totalsensflag : use total sensitivity intead of gain product
 *   verbose    : controls level of diagnostic output
 *
 * Returns 0 on success and -1 on error.
 *********************************************************************/
int
convolve_resp (double data[], int npts, double delta,
	       char *net, char *sta, char *loc, char *chan,
	       int startstage, int stopstage,
	       char *units, time_t resptime, int usedelay,
	       double taperfreq[], int *prewhiten,
	       int deconvflag, char *respfilename,
	       int totalsensflag, int verbose)
{
  int retval;
  int nfft;
  int nfreqs;
  double delfreq;

  double *xreal = NULL;
  double *ximag = NULL;
  
  nfft = next2 (npts);
  nfreqs = nfft/2 + 1;
  
  delfreq = 1.0 / (delta * nfft);
  
  /* Calculate frequency response from SEED RESP */
  if ( calcfr_resp (nfreqs, delfreq, net, sta, loc, chan, startstage,
		    stopstage, units, resptime, usedelay, respfilename,
		    totalsensflag, &xreal, &ximag, verbose) )
    {
      return -1;
    }
  
  /* Call main convolution routine */
  if ( deconvflag )
    retval = convolve (data, npts, delta, nfreqs, nfft, NULL, NULL,
		       xreal, ximag, taperfreq, prewhiten, verbose);
  else
    retval = convolve (data, npts, delta, nfreqs, nfft, xreal, ximag,
		       NULL, NULL, taperfreq, prewhiten, verbose);
  
  if ( xreal )
    free (xreal);
  if ( ximag )
    free (ximag);
  
  return retval;
}  /* End of convolve_resp() */


/*********************************************************************
 * next2:
 *
 * Return the power of 2 greater than or equal to value or 0 on error.
 *********************************************************************/
int
next2 (int value)
{
  int npower = 2;
  int nextpow2 = 0;
  
  while ( nextpow2 < value )
    {
      /* 2^X = (unsigned int) 1 << X; */
      nextpow2 = (unsigned int) 1 << npower++;
      
      /* Saftey check, stop at ~ 1e30 */
      if ( npower > 100 )
	return 0;
    }
  
  return (nextpow2);
}  /* End of next2() */


/*********************************************************************
 * spectraltaper:
 *
 * Taper spectra by a cosine.
 *
 * freq = frequency in question
 * fqh  = transition frequency between unity and the taper
 * fql  = transition frequency between zero and the taper
 *
 * if fql > fqh : lo-pass
 * if fqh > fql : hi-pass
 *
 * Return tapered value.
 ********************************************************************/
double
spectraltaper (double freq, double fqh, double fql)
{
  const double twopi = 6.283185307179586;
  double dblepi;
  double taper_v = 0;
  
  dblepi = 0.5*twopi;
  
  if ( fql > fqh )
    {
      if ( freq < fqh )
	taper_v = 1.0;
      if ( freq >= fqh && freq <= fql )
	taper_v = 0.5 * (1.0 + cos (dblepi*(freq - fqh)/(fql - fqh)));
      if ( freq > fql )
	taper_v = 0.0;
    }
  else if ( fqh > fql )
    {
      if ( freq < fql )
	taper_v = 0.0;
      if ( freq >= fql && freq <= fqh )
	taper_v = 0.5 * (1.0 - cos (dblepi*(freq - fql)/(fqh - fql)));
      if ( freq > fqh )
	taper_v = 1.0;
    }
  else
    {
      fprintf (stderr, "spectraltaper(): Invalid window specified\n");
      fprintf (stderr, "    freq: %g, fqh: %g, fql: %g\n", freq, fqh, fql);
    }
  
  return (taper_v);
}  /* End of spectraltaper() */


/*********************************************************************
 * findtaper:
 *
 * Determine spectral taper parameters (frequecies) from the amplitude
 * response.  This routine was designed primarily for responses with a
 * "flat" pass band.
 *
 * Any taper cutoff or pass values equal to -1.0 will be filled in
 * using the following criterion:
 *
 * taperfreq[0] - lower cutoff frequency
 *   -> lowest frequency with amplitude greater than dB down from maximum
 *
 * taperfreq[1] - lower pass frequency
 *   -> 110% of the lower cutoff frequency
 *
 * taperfreq[2] - upper pass frequency
 *   -> 95% of the upper cutoff frequency
 *
 * taperfreq[3] - upper cutoff frequency
 *   -> highest frequency with amplitude greater than dB down from maximum
 *
 * taperfreq = input taper frequencies, -1.0 values to be replaced
 * xreal    = real portion of the response
 * ximag    = imaginary portion of the response
 * nfreqs   = number of frequencies
 * delfreq  = frequency step
 * lcdBdown = lower corner cutoff specified as dB down from maximum.
 * ucdBdown = upper corner cutoff specified as dB down from maximum.
 *              if dB down is negative a default of 3dB will be used.
 *
 * Determined taper frequencies are stored directly in the taperfreq
 * parameters.
 *
 * This routine could probably benefit from some optimization at the
 * cost of clarity and simplicity.
 *
 * Return 0 on success and -1 on error.
 ********************************************************************/
int
findtaper (double *taperfreq, double *xreal, double *ximag,
	   int nfreqs, double delfreq,
	   double lcdBdown, double ucdBdown)
{
  int i;
  double amp;
  double maxamp = 0.0;
  double lowestfreq = 0.0;
  double highestfreq = 0.0;
  double dBdelta = 0.0;
  
  /*
  double lowestamp = 0.0;
  double highestamp = 0.0;
  double maxfreq = 0.0;
  */

  /* Set default dBdown if needed */
  if ( lcdBdown < 0.0 )
    lcdBdown = 3.0;
  if ( ucdBdown < 0.0 )
    ucdBdown = 3.0;
  
  /* Search for a lower frequency taper cutoff */
  if ( taperfreq[0] == -1.0 )
    {
      maxamp = 0.0;
      lowestfreq = 0.0;
      
      /* Determine frequency of maximum amplitude while searching
       * for lowest frequency with amplitude above the specified dB down.
       * The spectra in the first bin for DC (i==0) is specifically excluded */
      for ( i = (nfreqs-1); i > 0; i-- )
	{
	  amp = sqrt (xreal[i] * xreal[i] + ximag[i] * ximag[i]);
	  
	  /* Track the maximum amplitude/frequency and calculate minimum amplitue of 1% */
	  if ( maxamp == 0.0 || amp > maxamp )
	    {
	      maxamp = amp;
	      //maxfreq = i * delfreq;
	    }
	  
	  dBdelta = 20 * log10 (amp / maxamp);
	  
	  /* Track amplitude/frequency above the minimum amplitude */
	  if ( -dBdelta <= lcdBdown )
	    {
	      //lowestamp = amp;
	      lowestfreq = i * delfreq;
	    }
	}
      
      /* Set lower taper band to lowest frequency with amplitude at least 1% of maximum */
      taperfreq[0] = lowestfreq;
    }
  
  /* Set lower taper band to taper pass to 110% of cutoff */
  if ( taperfreq[1] == -1.0 )
    {
      taperfreq[1] = lowestfreq * 1.1;
    }
  
  /* Search for a upper frequency taper cutoff */
  if ( taperfreq[3] == -1.0 )
    {
      maxamp = 0.0;
      highestfreq = 0.0;
      
      /* Determine frequency of maximum amplitude while searching
       * for highest frequency with amplitude above the specified dB down. */
      for ( i = 0; i < nfreqs-1; i++ )
	{
	  amp = sqrt (xreal[i] * xreal[i] + ximag[i] * ximag[i]);
	  
	  /* Track the maximum amplitude/frequency and calculate minimum amplitue of 1% */
	  if ( maxamp == 0.0 || amp > maxamp )
	    {
	      maxamp = amp;
	      //maxfreq = i * delfreq;
	    }
	  
	  dBdelta = 20 * log10 (amp / maxamp);
	  
	  /* Track amplitude/frequency above the minimum amplitude */
	  if ( -dBdelta <= ucdBdown )
	    {
	      //highestamp = amp;
	      highestfreq = i * delfreq;
	    }
	}
      
      /* Set upper taper band to lowest frequency with amplitude at least 1% of maximum */
      taperfreq[3] = highestfreq;
    }
  
  /* Set upper taper band to taper pass to 95% of cutoff */
  if ( taperfreq[2] == -1.0 )
    {
      taperfreq[2] = taperfreq[3] * 0.95;
    }
  
  return 0;
}  /* End of findtaper() */


/*********************************************************************
 * fft:
 * Discrete FFT for perfect powers of two.
 *
 * The sine tables are calcualted when the fuction is first used and
 * extended when needed for subsequent calls where nfft has
 * increased.
 *
 * direction: (+1)=Forward DFT, (-1)=Inverse DFT
 *
 * There is a minor optimization if this routine is called with the
 * same nfft repeatedly, the sine tables will be not recalculated.
 *
 * This routine can be called with the first two arguments set to 0 to
 * request clean-up of allocated sine table buffers.
 *
 * This is a combination of modified versions of routines from
 * Numerical Utilities (NumUtils), www.xgraph.org/numutil.html
 *
 * Returns 0 on success and -1 on error.
 *********************************************************************/
int
fft (double real[], double imag[], int nfft, int direction)
{
  static double *sintab1 = 0;
  static double *sintab2 = 0;
  static int sinallocated = 0;
  
  int i, step, j, m, mmax;
  int n2, k;
  double sinphi, phi, coef_i, coef_r, sci, scr;
  double tmpr, tmpi, pi2;
  double treal, timag;
  
  /* Perform clean-up if requested */
  if ( ! real && ! imag )
    {
      if ( sintab1 )
	free (sintab1);
      if ( sintab2 )
	free (sintab2);
      sintab1 = 0;
      sintab2 = 0;
      sinallocated = 0;
    }
  
  if ( direction != 1 && direction != -1 )
    {
      fprintf (stderr, "fft(): direction argument must be 1 (forward) or -1 (inverse)\n");
      return -1;
    }
  
  /* Sanity check fft length */
  m = 2;
  while ( m < nfft )
    {
      m = m << 1;
    }
  if ( m != nfft )
    {
      fprintf (stderr, "fft(): nfft (%d) is not a perfect power of two\n", nfft);
      return -1;
    }
  
  /* Allocate sine tables if needed */
  if ( sinallocated != nfft )
    {
      if ( (sintab1 = (double *) realloc (sintab1, nfft * sizeof(double))) == NULL )
	{
	  fprintf (stderr, "Error allocating memory\n");
	  return -1;
	}
      if ( (sintab2 = (double *) realloc (sintab2, nfft * sizeof(double))) == NULL )
	{
	  fprintf (stderr, "Error allocating memory\n");
	  return -1;
	}
      
      pi2 = 0.5 * PI;
      mmax = 1;
      k = 0;
      
      while ( mmax < nfft )
	{
	  step = mmax << 1;
	  phi = pi2 / (double) mmax;
	  sintab1[k] = sin (phi);
	  sintab2[k] = sin (phi * 2.0);
	  k++;
	  m = 0;
	  while ( m < mmax )
	    {
	      i = m;
	      while ( i < nfft )
		{
		  j = i + mmax;
		  i += step;
		}
	      m++;
	    }
	  
	  mmax = step;
	}
      
      sinallocated = nfft;
    }
  
  n2 = nfft >> 1;
  j = 0;
  i = 0;
  
  /* Perform decimation in time by bit reversing up front */
  while ( i < nfft )
    {
      if ( i < j )
	{
	  treal = real[j]; timag = imag[j];
	  real[j] = real[i]; imag[j] = imag[i];
	  real[i] = treal; imag[i] = timag;
	}
      m = n2;
      while ( (j > m-1) && (m >= 2) )
	{ j -= m;  m = m >> 1; }
      j += m;
      i++;
    }
  
  mmax = 1;
  k = 0;
  while ( mmax < nfft )
    {
      step = mmax << 1;
      sinphi = -1 * (double)direction * sintab1[k];
      sci = -1 * (double)direction * sintab2[k];
      scr = -2.0 * sinphi * sinphi;
      k++;
      m = 0;
      coef_r = 1.0;
      coef_i = 0.0;
      
      while ( m < mmax )
	{
	  i = m;
	  while ( i < nfft )
	    {
	      j = i + mmax;
	      tmpr = (coef_r * real[j]) - (coef_i * imag[j]);
	      tmpi = (coef_r * imag[j]) + (coef_i * real[j]);
	      real[j] = real[i] - tmpr;
	      imag[j] = imag[i] - tmpi;
	      real[i] += tmpr;
	      imag[i] += tmpi;
	      i += step;
	    }
	  
	  tmpr = coef_r;
	  m++;
	  coef_r = coef_r + (coef_r * scr) - (coef_i * sci);
	  coef_i = coef_i + (coef_i * scr) + (tmpr * sci);
	}
      
      mmax = step;
    }
  
  /* A forward transform will scale the values by nfft, compensate. */
  /* This particular code base does not require this compensation
     if ( direction > 0 )
     {
      double scale = 1.0 / (double) nfft;
      for ( i=0; i < nfft; i++ )
	{
	  real[i] *= scale;
	  imag[i] *= scale;
	}
    }
  */
  
  return 0;
}  /* End of fft() */
