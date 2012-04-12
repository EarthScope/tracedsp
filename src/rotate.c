/*********************************************************************
 * rotate.c
 *
 * 2-D and 3-D series rotation routines.  Originally these came from
 * the source code for SeismicHander by Klaus Stammler.  They have
 * been modified for use in this code base.
 *
 * Modified: 2007.181
 *********************************************************************/

/* Originally these routines came from the source code for SeismicHandler
 * with the notices below.
 *
 * original file name SHMATH.C
 *      ========
 *
 * version 34, 22-May-2006
 *
 * mathematical subroutines
 * K. Stammler, 16-MAY-1990
 */

/*
 *  SeismicHandler, seismic analysis software
 *  Copyright (C) 1992,  Klaus Stammler, Federal Institute for Geosciences
 *                                       and Natural Resources (BGR), Germany
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "rotate.h"

/* Pi, with comical precision */
#define PI 3.1415926535897932384626433832795

/* half of square root of 3 */
#define SQRT3_2 0.8660254

/* square root of 3 */
#define SQRT3 1.7320508

#define SQRT1_3 5.77350269e-1
#define SQRT1_2 7.07106781e-1
#define SQRT2_3 8.16496581e-1


/* 2-dimensional rotation
 *
 * parameter of routine
 * double    *n, *e;      input; north & east component
 * long      lth;         input; length of input & output traces
 * double    angle;       input; rotation angle
 * double    *r, *t;      output; rotated traces
 */
void
rotate2 ( double *n, double *e, long lth, double angle,
	  double *r, double *t )
{
  double ann, ane, aen, aee;  /* rotation matrix */
  long  i;                    /* counter */
  
  /* Determine rotation matrix */
  angle *= (PI / 180.0);
  ann = cos( angle );
  ane = -sin( angle );
  aen = -ane;
  aee = ann;
  
  /* Apply rotation matrix */
  for  (i=0; i < lth ;i++)
    {
      *r++ = ann * *n + ane * *e;
      *t++ = aen * *n + aee * *e;
      e++; n++;
    }
  
} /* End of rotate2 */


/* 3-dimensional rotation
 *
 * parameter of routine
 * double    *z, *n, *e;  input; vertical, north & east component
 * long      lth;         input; length of input & output traces
 * double    azim, inci;  input; rotation angles
 * int       type;        input; type of rotation
 * double    *l, *q, *t;  output; rotated traces
 */
void
rotate3 ( double *z, double *n, double *e, long lth,
	  double azim, double inci, int type,
	  double *l, double *q, double *t )
{
  double  azz, azn, aze;  /* rotation matrix */
  double  anz, ann, ane;
  double  aez, aen, aee;
  double  cosi, cosa, sini, sina;
  long  i;                   /* counter */
  
  /* Determine rotation matrix */
  azim *= (PI / 180.0);
  cosa = cos( azim );
  sina = sin( azim );
  if  (type == ROT_ZNE_TO_LQT)
    {
      inci *= (PI / 180.0);
      cosi = cos( inci );
      sini = sin( inci );
    }
  
  if  (type == ROT_ZNE_TO_LQT)
    {
      azz = cosi;         azn = -sini*cosa;     aze = -sini*sina;
      anz = sini;         ann =  cosi*cosa;     ane =  cosi*sina;
      aez = 0.0;          aen =  sina;          aee = -cosa;
    }
  else if  (type == ROT_ZNE_TO_UVW)
    {
#     ifdef XXX
      if  (sina != 0.0)   sina = 1.0 / sina;
      if  (cosa != 0.0)   cosa = 1.0 / cosa;
      azz = cosa/3.0;     anz = cosa/3.0;       aez = cosa/3.0;
      azn = 0.;           ann = sina/SQRT3;     aen = -sina/SQRT3;
      aze = -sina*2.0/3.0;ane = sina/3.0;       aee = sina/3.0;
#     endif
      azz = SQRT1_3;      anz = SQRT1_3;        aez = SQRT1_3;
      azn = 0.0;          ann = SQRT1_2;        aen = -SQRT1_2;
      aze = -SQRT2_3;     ane = 0.5*SQRT2_3;    aee = 0.5*SQRT2_3;
    }
  else if  (type == ROT_UVW_TO_ZNE)
    {
#     ifdef XXX
      azz = cosa;         anz = cosa;           aez = cosa;
      azn = 0.;           ann = sina*SQRT3_2;   aen = -sina*SQRT3_2;
      aze = -sina;        ane = sina*0.5;       aee = sina*0.5;
#     endif
      azz = SQRT1_3;      azn = SQRT1_3;        aze = SQRT1_3;
      anz = 0.0;          ann = SQRT1_2;        ane = -SQRT1_2;
      aez = -SQRT2_3;     aen = 0.5*SQRT2_3;    aee = 0.5*SQRT2_3;
    }
  else
    {
      azz = azn = aze = anz = ann = ane = aez = aen = aee = 0.0;
    }
  
  /* Apply rotation matrix */
  for (i=0; i < lth ;i++)
    {
      *l++ = azz * *z  +  azn * *n  +  aze * *e;
      *q++ = anz * *z  +  ann * *n  +  ane * *e;
      *t++ = aez * *z  +  aen * *n  +  aee * *e;
      z++; n++; e++;
    }
  
} /* End of rotate3 */
