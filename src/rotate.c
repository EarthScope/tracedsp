/*********************************************************************
 * rotate.c
 *
 * 2-D and 3-D series rotation routines.
 *
 * Modified: 2012.306
 *********************************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "rotate.h"

/* Pi, with comical precision */
#define PI 3.1415926535897932384626433832795


/*********************************************************************
 * rotate2:
 *
 * Perform 2-D rotation.
 *
 * If the N and E input data are positive to the North and East
 * respectively the resulting R component will be positive along the
 * azimuth rotated to and the T component will be positive along an
 * axis 90 clockwise of the R component.
 *
 * Arguments:
 *   n       : The "north" component, 90 degress counter-clockwise from "east"
 *   e       : The "east" component, 90 degress clockwise from "north"
 *   npts    : number of samples in each component array
 *   azimuth : rotation angle (clockwise) in degrees
 *   r       : radial output component, may be the same as input array
 *   t       : transverse output component, may be the same as input array
 *
 *********************************************************************/
void
rotate2 ( double *n, double *e, long int npts, double azimuth,
	  double *r, double *t )
{
  double sina, cosa;
  double N, E;
  long int i;
  
  azimuth *= (PI / 180.0);
  sina = sin ( azimuth );
  cosa = cos ( azimuth );
  
  /* Apply rotation matrix */
  for ( i=0; i < npts; i++ )
    {
      N = n[i];
      E = e[i];
      
      r[i] = sina * E + cosa * N;
      t[i] = cosa * E - sina * N;
    }
  
} /* End of rotate2 */


/*********************************************************************
 * rotate3:
 *
 * Perform 3-D rotation to the LQT, ray-oriented system.
 *
 * If the N and E input data are positive to the North and East
 * respectively and Z is positive up the resulting components are:
 *
 *   L: Longitudinal, in the direction of the ray away from the source
 *        as defined by the incidence and azimuth (P wave energy).
 *   Q: 90 degrees from L in the vertical plane pointing upward (SV energy).
 *   T: 90 degrees from L and Q in the horizontal plane (SH energy).
 *
 * Arguments:
 *   z         : The "vertical" component, perpendicular to "north" and "east"
 *   n         : The "north" component, 90 degress counter-clockwise from "east"
 *   e         : The "east" component, 90 degress clockwise from "north"
 *   npts      : number of samples in each component array
 *   azimuth   : rotation angle (degrees clockwise) for "north" and "east"
 *   incidence : rotation angle (degress from vertical) for radial and "vertical"
 *   l         : ray direction output component, may be the same as input array
 *   q         : transverse output component, may be the same as input array
 *   t         : transverse output component, may be the same as input array
 *
 *********************************************************************/
void
rotate3 ( double *z, double *n, double *e, long int npts,
	  double azimuth, double incidence,
	  double *l, double *q, double *t )
{
  double sina, cosa;
  double sini, cosi;
  double Z, N, E;
  long int i;
  
  /* Convert azimuth to backazimuth needed for rotation matrix */
  azimuth = ( azimuth >= 180.0 ) ? azimuth - 180.0 : azimuth + 180.0;
  
  azimuth *= (PI / 180.0);
  incidence *= (PI / 180.0);
  
  sina = sin ( azimuth );
  cosa = cos ( azimuth );
  sini = sin ( incidence );
  cosi = cos ( incidence );
  
  /* Apply rotation matrix */
  for ( i=0; i < npts; i++ )
    {
      Z = z[i];
      N = n[i];
      E = e[i];
      
      l[i] = cosi * Z - sini * cosa * N - sini * sina * E;
      q[i] = sini * Z + cosi * cosa * N + cosi * sina * E;
      t[i] = sina * N - cosa * E;
    }
  
} /* End of rotate3 */
