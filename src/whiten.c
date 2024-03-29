/*********************************************************************
 * whiten.c
 *
 * Routines for whitening and dewhitening time-series data.  The
 * whitening process includes designing a prediction error filter
 * based on the autocorrelation function of the input sequence and
 * then applying this filter.  The dewhitening process then applies a
 * prediction filter using the same coefficients.
 *
 * modified: 2012.098
 *********************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "whiten.h"

/* Maximum number of linear prediction coefficients (filter order) */
#define NCMAX 12

#define MAX(A, B) ((A) > (B) ? (A) : (B)) /* Max macro definition */
#define MIN(A, B) ((A) > (B) ? (B) : (A)) /* Min macro definition */

/*********************************************************************
 * whiten:
 *
 * Whiten an input sequence in-place using a prediction error filter.
 * The linear prediction coefficients for the prediction error filter
 * are determined directly from the input sequence and return for
 * later use (presumably for use with dewhiten()).  Modified version
 * of SAC 2000's prewit() routine.
 *
 * Arguments:
 *   data  :  array of input data samples
 *   npts  :  number of data samples in the data array
 *  *order :  int pointer to predictor order
 *   coef  :  array of linear prediction coefficients (returned)
 *
 * Returns 0 on success and -1 on error.
 *********************************************************************/
int
whiten (double data[], int npts, int *order, double coef[])
{
  int idx, jdx, j2, kdx, kb, torder;
  double at, q;
  double ref[NCMAX + 1];
  double sa[NCMAX + 1];
  double pefcoef[NCMAX];

  /* Range check for predictor order */
  if (*order < 1 || *order > NCMAX)
  {
    fprintf (stderr, "Error, prewit(): Predictor order out of bounds (1-12): %d", *order);
    return -1;
  }

  /* Design the prewhitening filter */
  LPCautocorr (data, npts, *order, coef, ref);

  /* DEAD CODE: The code below that tests for filter stability is in
   * the original routine and not enabled for unknown reasons.  It is
   * left in this version for future reference as the capability might
   * be resurrected at some point. */

  /* Check for stability, and truncate filter if necessary
   * so the recursive de-whitener is safely stable
   * Use an adaptation of LPTRN code from Markel and Gray
   * found in the IEEE ASSP book of signal processing codes */
  torder = *order;

  /* do 100 i = NCMAX+1, 1, -1
   *   if ( abs( ref(i) ) .gt. 0.95 ) torder = i - 1
   * 100  continue
   *
   * If any reflection coefficients are too close to +1 or -1
   * then the system is getting dangerous, so truncate, and
   * then regenerate the filter coeffs for the truncated filter */
  if (torder != *order)
  {
    for (idx = 0; idx < torder; idx++)
    {
      sa[idx] = ref[idx];
    }
    for (jdx = 1; jdx < torder; jdx++)
    {
      j2 = (jdx + 1) / 2;
      q  = ref[jdx];
      for (kdx = 0; kdx < j2; kdx++)
      {
        kb      = jdx - kdx;
        at      = sa[kdx] + q * sa[kb];
        sa[kb]  = sa[kb] + q * sa[kdx];
        sa[kdx] = at;
      }
    }
    for (idx = 0; idx < torder; idx++)
    {
      coef[idx + 1] = sa[idx];
    }
    for (idx = torder + 1; idx < (*order + 1); idx++)
    {
      coef[idx] = 0.0;
    }
    coef[0] = 1.0;
    *order  = torder;
  }

  /*  Prewhiten the data
   *
   *    Negate coefficients of prediction error filter generated by
   *    LPCautocorr().  For historical reasons, different storage
   *    modes are used in LPCautocorr() and applypef() and
   *    applypf(). */
  for (idx = 0; idx < *order; idx++)
  {
    pefcoef[idx] = -coef[idx];
  }

  /* Apply prediction error filter */
  if (applypef (data, npts, pefcoef, *order, data))
    return -1;

  return 0;
} /* End of whiten() */

/*********************************************************************
 * dewhiten:
 *
 * Dewhiten an input sequence in-place using a prediction filter
 * calculated with whiten().  Modified version of SAC 2000's dewit()
 * routine.
 *
 * Arguments:
 *   data  :  array of input data samples
 *   npts  :  number of data samples in the data array
 *   order :  predictor order
 *   coef  :  array of linear prediction coefficients
 *
 * Returns 0 on success and -1 on error.
 *********************************************************************/
int
dewhiten (double data[], int npts, int order, double coef[])
{
  long int i;
  double precoef[NCMAX];

  /* Range check for predictor order */
  if (order < 1 || order > 12)
  {
    fprintf (stderr, "Error, dewhiten(): Predictor order out of bounds (1-12): %d", order);
    return -1;
  }

  /*  Dewhiten the data
   *
   *    Negate coefficients of prediction error filter generated by
   *    LPCautocorr().  For historical reasons, different storage
   *    modes are used in LPCautocorr() and applypef() and
   *    applypf(). */
  for (i = 0; i < order; i++)
  {
    precoef[i] = -coef[i];
  }

  /* Apply prediction filter */
  if (applypf (data, npts, precoef, order, data))
    return -1;

  return 0;
} /* End of dewhiten() */

/*********************************************************************
 * LPCautocorr:
 *
 * Calculate Linear Prediction Coefficients and related reflection
 * coefficients by applying the Levinson-Durbin algorithm to the
 * autocorrelation function of the input data.
 *
 * Input arguments:
 *   data  :  array of input data samples
 *   npts  :  number of data samples in the data array
 *   order :  predictor order
 *   lpc   :  array of linear prediction coefficients
 *   ref   :  array of reflection coefficients
 *
 * Returns the minimum mean square error
 *********************************************************************/
double
LPCautocorr (double data[], int npts, int order,
             double lpc[], double ref[])
{
  double corr[NCMAX];
  double sum;
  double r;
  double error;
  double ctemp;
  int idx, jdx;
  int aclag = order + 1;

  /* Compute autocorrelation */
  while (aclag--)
  {
    for (idx = aclag, sum = 0; idx < npts; idx++)
      sum += data[idx] * data[idx - aclag];

    corr[aclag] = sum;
  }

  if (corr[0] == 0.0)
  {
    for (idx = 0; idx < order; idx++)
      ref[idx] = 0.0;

    return 0.0;
  }

  error = corr[0];

  for (idx = 0; idx < order; idx++)
  {
    /* Sum up this iteration's reflection coefficient */
    r = -corr[idx + 1];

    for (jdx = 0; jdx < idx; jdx++)
      r -= lpc[jdx] * corr[idx - jdx];

    ref[idx] = r /= error;

    /* Update LPC coefficients and total error */
    lpc[idx] = r;

    for (jdx = 0; jdx < idx / 2; jdx++)
    {
      ctemp = lpc[jdx];
      lpc[jdx] += r * lpc[idx - 1 - jdx];
      lpc[idx - 1 - jdx] += r * ctemp;
    }

    if (idx % 2)
      lpc[jdx] += lpc[jdx] * r;

    error *= 1.0 - r * r;
  }

  return error;
} /* End of LPCautocorr() */

/*********************************************************************
 * applypef:
 *
 * Apply prediction error filter to data given linear prediction
 * coefficients.  Modified version of SAC 2000's pef() routine.
 *
 * data   :  array of input data samples
 * npts   :  number of data samples in the data array
 * coef   :  prediction error filter coefficients
 * nc     :  number of pef coefficients
 * result :  array of output data samples, can be same as input.
 *
 * Returns the minimum mean square error
 *********************************************************************/
int
applypef (double data[], int npts, double coef[], int nc, double result[])
{
  int bufptr, datptr;
  int idx, lsamp, ncmp, nextra;
  double e;

  double inbuf[2000];
  double outbuf[2000];

  /* Initializations */
  datptr = 0;
  nextra = 2000 - nc;

  if (nextra < 1)
  {
    fprintf (stderr, "Error, applypef(): Filter too large: %d\n", nc);
    return -1;
  }

  for (idx = 0; idx < 2000; idx++)
    inbuf[idx] = outbuf[idx] = 0.0;

  while (datptr < npts)
  {
    /* Shift input buffer points back */
    memmove (inbuf, &inbuf[nextra], nc * sizeof (double));

    /* Calculate start and stop points for buffer index */
    bufptr = nc;
    lsamp  = MIN (2000, npts - datptr + nc);
    ncmp   = lsamp - nc;

    /* Load new data into input buffer */
    memcpy (&inbuf[nc], &data[datptr], ncmp * sizeof (double));

    /* Filter data */
    while (bufptr <= lsamp)
    {
      e = inbuf[bufptr];

      for (idx = 0; idx < nc; idx++)
      {
        e -= coef[idx] * inbuf[bufptr - idx - 1];
      }

      outbuf[bufptr] = e;

      bufptr++;
    }

    /* Store newly filtered data and update data pointer */
    memcpy (&result[datptr], &outbuf[nc], ncmp * sizeof (double));

    datptr += ncmp;
  }

  return 0;
} /* End of applypef() */

/*********************************************************************
 * applypf:
 *
 * Apply prediction filter to data given linear prediction
 * coefficients.  Modified version of SAC 2000's predflt() routine.
 *
 * data   :  array of input data samples
 * npts   :  number of data samples in the data array
 * coef   :  prediction error filter coefficients
 * nc     :  number of pef coefficients
 * result :  array of output data samples, can be same as input.
 *
 * Returns the minimum mean square error
 *********************************************************************/
int
applypf (double data[], int npts, double coef[], int nc, double result[])
{
  int bufptr, dataptr;
  int idx, lsamp, ncmp;
  double history[NCMAX];
  double eo;

  double inbuf[2000];
  double outbuf[2000];

  if (nc > NCMAX)
  {
    fprintf (stderr, "Error, applypf(): Filter too large: %d\n", nc);
    return -1;
  }

  for (idx = 0; idx < 2000; idx++)
    inbuf[idx] = outbuf[idx] = 0.0;

  for (idx = 0; idx < NCMAX; idx++)
    history[idx] = 0.0;

  dataptr = 0;
  while (dataptr < npts)
  {
    /* Calculate start and stop points for buffer index */
    bufptr = 0;
    lsamp  = MIN (2000, (npts - dataptr));
    ncmp   = lsamp;

    /* Load new data into input buffer */
    memcpy (inbuf, &data[dataptr], ncmp * sizeof (double));

    /* Filter data */
    while (bufptr < lsamp)
    {
      eo = inbuf[bufptr];

      for (idx = (nc - 1); idx >= 0; idx--)
      {
        eo += coef[idx] * history[idx];

        if (idx > 0)
          history[idx] = history[idx - 1];
      }

      outbuf[bufptr] = eo;
      history[0]     = eo;

      bufptr++;
    }

    /* Store newly filtered data */
    memcpy (&result[dataptr], outbuf, ncmp * sizeof (double));

    /* Increment data pointer */
    dataptr += ncmp;
  }

  return 0;
} /* End of applypf() */
