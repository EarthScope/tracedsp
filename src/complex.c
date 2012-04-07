/*********************************************************************
 * complex.c:
 *
 * Complex number operation routines.
 *
 * Based on routines from the SAC 2000 source code.
 *
 * modified: 2012.098
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "complex.h"

/* Return the imaginary part of a complex number */
double
cmplximag (complexd c)
{
  return (c.imag);
}

/* Return the real part of a complex number */
double
cmplxreal (complexd c)
{
  return (c.real);
}

/* Add two complex numbers */
complexd
cmplxadd (complexd c1, complexd c2)
{
  c1.real += c2.real;
  c1.imag += c2.imag;
  return (c1);
}

/* Return complex conjugate */
complexd
cmplxconj (complexd c)
{
  c.imag = -c.imag;
  return (c);
}

/* Multiple two complex numbers */
complexd
cmplxmul (complexd c1, complexd c2)
{
  complexd c3;
  
  c3.real = (c1.real * c2.real) - (c1.imag * c2.imag);
  c3.imag = (c1.real * c2.imag) + (c1.imag * c2.real);
  return (c3);
}

/* Multiple a real number and a complex number */
complexd
cmplxdmul (double d, complexd c)
{
  c.real *= d;
  c.imag *= d;
  return (c);
}

/* Create complex number from double */
complexd
dbltocmplx (double d1, double d2)
{
  complexd c;
  c.real = d1;
  c.imag = d2;
  return (c);
}

/* Subtract two complex numbers */
complexd
cmplxsub (complexd c1, complexd c2)
{
  c1.real -= c2.real;
  c1.imag -= c2.imag;
  return (c1);
}

/* Absolute value of complex number */
double
cmplxabs (complexd c)
{
  return (sqrt((double)((c.real*c.real) + (c.imag*c.imag))));
} 

double
cmplxang (complexd c)
{
  double d = 1.0;

  if ( c.imag < 1.0e-14 && c.imag > -1.0e-14 && fabs(c.real) > 1.0e-7 )
    {
      /* d = acos(c.real/cmplxabs(c)); */
      if ( c.real > 0 ) 
	d = 0;
      else if ( c.real < 0 ) 
	d = PI;
    }
  else if (c.real < 1.0e-14 && c.real > -1.0e-14 && fabs(c.imag) > 1.0e-7 )
    {
      /* d = asin(c.imag/cmplxabs(c)); */
      if ( c.imag > 0 ) 
	d = PI / 2;
      else if ( c.real < 0 ) 
	d = -PI / 2;
    }
  else
    {
      d = atan(c.imag/c.real);
      if (c.real < 0.0)
	{
	  if (c.imag < 0.0) d -= PI;
	  else              d += PI;
	}
    }

  return (d);
}

/* Return square root of complex number */
complexd
cmplxsqrt (complexd c)
{
  double sqrtsave, angle;

  sqrtsave = sqrt (cmplxabs(c));
  angle = cmplxang(c);
  
  c.real = sqrtsave * cos (angle/2.0);
  c.imag = sqrtsave * sin (angle/2.0);

  if (c.real < 0.0)
    {
      c.real = -c.real;
      c.imag = -c.imag;
    }
  else if (c.real < 1.0e-14 && c.real > -1.0e-14 && c.imag < 0.0)
    c.imag = -c.imag;
  
  return (c);
}

/* Divide complex numbers */
complexd
cmplxdiv (complexd c1, complexd c2)
{
  complexd c;
  double d;
  
  if (c2.real == 0.0 && c2.imag == 0.0)
    {
      fprintf(stderr, "cmplxdiv(): complex divide by zero\n");
      c.real = c.imag = 0.0; 
    }
  else
    {
      d = c2.real*c2.real + c2.imag*c2.imag;
      
      c.real = (c1.real*c2.real + c1.imag*c2.imag) / d;
      c.imag = (c2.real*c1.imag - c1.real*c2.imag) / d;
    }
  
  return (c);
}

/* Logarithm of complex number */
complexd
cmplxlog (complexd c)
{
  complexd c1;
  
  c1.real = log (cmplxabs(c));
  c1.imag = cmplxang (c);
  
  return (c1);
}

/* Raise natural logarithm e to the power of the complex number */
complexd
cmplxexp (complexd c)
{
  double d;
  
  if (c.real == 0.0)
    d = 1.0;
  else
    d = exp(c.real);
  
  if (c.imag == 0.0)
    {
      c.real = d;
      return (c);
    }
  
  c.real = d * cos (c.imag);
  c.imag = d * sin (c.imag);

  return (c);
}

/* Raise complex number to specified power */
complexd
cmplxpow (complexd c, double d)
{
  if (c.real == 0.0 && c.imag == 0.0)
    return (c);

  c = cmplxlog (c);
  c.real = d*c.real;
  c.imag = d*c.imag;
  
  return (cmplxexp(c));
}

/* Negate complex number */
complexd
cmplxneg (complexd c)
{
  c.real = -c.real;
  c.imag = -c.imag;
  return (c);
}
