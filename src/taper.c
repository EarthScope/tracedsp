/*********************************************************************
 * taper.c
 *
 * Symmetric series tapering.
 *
 * Modified: 2011.024
 *********************************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "taper.h"

/* Pi, with comical precision */
#define PI 3.1415926535897932384626433

/********************************************************************
 * taper:
 *
 * -- notes from SAC help --
 *
 * The general form for the taper is:
 *
 * DATA(J)=DATA(J)*(F0-F1*COS(OMEGA*(J-1))
 *
 * This equation would be applied to the left hand side of each
 * signal.  A symmetric one is applied to the right hand side.  The
 * following table defines the various parameters used in the
 * different tapers.  In this table N is the length of the taper on
 * each end.
 *
 * TYPE     OMEGA     F0    F1
 * HANNING  PI/N      0.50  0.50
 * HAMMING  PI/N      0.54  0.46
 * COSINE   PI/(2*N)  1.00  1.00  // Not used in SAC
 *
 * Arguments:
 *   data       : array of input data samples as doubles
 *   npts       : number of samples
 *   width      : taper window width in percent (0 - 0.5)
 *   type       : TAPER_HANNING, TAPER_HAMMING or TAPER_COSINE
 *
 * Returns samples in window length on success and -1 on error
 ********************************************************************/
int
taper (double *data, int npts, double width, int type)
{
  double omega;
  double factor;
  int widthpts;
  int idx;
  
  /* Sanity checks */
  if ( ! data )
    {
      fprintf (stderr, "taper(): No data buffer specified\n");
      return -1;
    }
  
  if ( width > 0.5 )
    {
      fprintf (stderr, "taper(): Window width (%g) cannot be larger than 0.5\n",
	       width);
      return -1;
    }
  
  /* Calculate number of samples in window width */
  widthpts = (npts * width);
  
  if ( npts <= 0 )
    return 0;
  
  /* Calculate omega */
  switch ( type )
    {
    case TAPER_HANNING:
    case TAPER_HAMMING:
      omega = PI / widthpts;
      break;
    case TAPER_COSINE:
      omega = PI / (2 * widthpts);
      break;
    default:
      fprintf (stderr, "taper(): Unrecognized taper type: %d\n", type);
      return -1;
      break;
    }
  
  /* Apply taper to both ends of trace */
  switch ( type )
    {
    case TAPER_HANNING:
      for ( idx=0 ; idx < widthpts ; idx++ )
	{
	  factor = (0.50 - 0.50 * cos(omega * idx));
	  data[idx] *= factor;
	  data[npts - 1 - idx] *= factor;
	}
      break;
    case TAPER_HAMMING:
      for ( idx=0 ; idx < widthpts ; idx++ )
	{
	  factor = (0.54 - 0.46 * cos(omega * idx));
	  data[idx] *= factor;
	  data[npts - 1 - idx] *= factor;
	}
      break;
    case TAPER_COSINE:
      for ( idx=0 ; idx < widthpts ; idx++ )
	{
	  factor = sin(omega * (idx)); /* Matches SAC "cosine" tapering */
	  data[idx] *= factor;
	  data[npts - 1 - idx] *= factor;
	}
      break;
    }
  
  return widthpts;
}  /* End of taper() */
