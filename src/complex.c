/*********************************************************************
 * complex.c:
 *
 * Complex number operation routines.
 *
 * Based on routines from the SAC 2000 source code.
 *
 * modified: 2006.160
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "complex.h"

double
cmplximag (complexf c)
{
  return (c.imag);
}

double
cmplxreal (complexf c)
{
  return (c.real);
}

complexf
cmplxadd (complexf c1, complexf c2)
{
  c1.real += c2.real;
  c1.imag += c2.imag;
  return (c1);
}

complexf
cmplxcj (complexf c)
{
  c.imag = -c.imag;
  return (c);
}

complexf
cmplxmul (complexf c1, complexf c2)
{
  complexf c3;
  
  c3.real = (c1.real * c2.real) - (c1.imag * c2.imag);
  c3.imag = (c1.real * c2.imag) + (c1.imag * c2.real);
  return (c3);
}

complexf
flttocmplx (double d1, double d2)
{
  complexf c;
  c.real = d1;
  c.imag = d2;
  return (c);
}

complexf
cmplxsub (complexf c1, complexf c2)
{
  c1.real -= c2.real;
  c1.imag -= c2.imag;
  return (c1);
}

double
cmplxabs (complexf c)
{
  return (sqrt((double)((c.real*c.real) + (c.imag*c.imag))));
} 

double
cmplxang (complexf c)
{
  double d = 1.0;

  if ( c.imag < 1.0e-14 && c.imag > -1.0e-14 && fabs(c.real) > 1.0e-7 )
    {
      /* d = acos(c.real/cmplxabs(c)); */
      if ( c.real > 0 ) 
	d = 0;
      else if ( c.real < 0 ) 
	d = 3.1415926536;
    }
  else if (c.real < 1.0e-14 && c.real > -1.0e-14 && fabs(c.imag) > 1.0e-7 )
    {
      /* d = asin(c.imag/cmplxabs(c)); */
      if ( c.imag > 0 ) 
	d = 3.1415926536/2;
      else if ( c.real < 0 ) 
	d = -3.1415926536/2;
    }
  else
    {
      d = atan(c.imag/c.real);
      if (c.real < 0.0)
	{
	  if (c.imag < 0.0) d -= 3.1415926536;
	  else           d += 3.1415926536;
	}
    }

  return (d);
}

complexf
cmplxsqrt (complexf c)
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

complexf
cmplxdiv (complexf c1, complexf c2)
{
  complexf c;
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

complexf
cmplxlog (complexf c)
{
  complexf c1;
  
  c1.real = log (cmplxabs(c));
  c1.imag = cmplxang (c);
  
  return (c1);
}

complexf
cmplxexp (complexf c)
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

complexf
cmplxpow (complexf c, double d)
{
  if (c.real == 0.0 && c.imag == 0.0)
    return (c);

  c = cmplxlog (c);
  c.real = d*c.real;
  c.imag = d*c.imag;
  
  return (cmplxexp(c));
}

complexf
cmplxneg (complexf c)
{
  c.real = -c.real;
  c.imag = -c.imag;
  return (c);
}
