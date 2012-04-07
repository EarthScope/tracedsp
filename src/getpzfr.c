
/*********************************************************************
 *
 * This code is designed to read a SAC-style poles and zeros file and
 * calculate their corresponding complex frequency response.
 *
 * The calcfr() routine is a cleaned up version of the equivalent in
 * SAC 2000.
 *
 * modified: 2012.098
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>

#include "complex.h"

#define	MAXPOLES 30
#define	MAXZEROS 30

void calcfr (int nfreq, double delfreq, double constant,
	     int nzero, complexd zeros[],
	     int npole, complexd poles[],
	     double xreal[], double ximag[]);


/*********************************************************************
 * getpzfr():
 *
 * Read a SAC format poles and zeros file and calculate their
 * frequency response using calcfr().
 *
 * nfreq   : number of points that will be used in the FFT
 * delfreq : time-interval between data points in the time-domain
 * xreal   : array of real part of frequency response
 * ximag   : array of imaginary part of frequncy response
 *           xreal and ximag must already be allocated with enough
 *           room for nfreq entries.
 * pzfilename : name of file containing poles and zeros in SAC format
 *
 * Return 0 on success and -1 on error.
 *********************************************************************/
int
getpzfr (int nfreq, double delfreq, double xreal[], double ximag[],
         char *pzfilename)
{
  int   idx;
  FILE *fp;
  
  char  line[1024];
  char *lp;
  
  char  reading = 0; /* 0 = unknown, 1 = zeros, 2 = poles */
  int   linecount = 0;
  
  complexd zeros[MAXZEROS];
  complexd poles[MAXPOLES];
  double constant;
  int rdconstant = 0;
  int nzeros  = 0;
  int zeroidx = 0;
  int rdzeros = 0;
  int npoles  = 0;
  int poleidx = 0;
  int rdpoles = 0;

  double reald;
  double imagd;
  
  /* Initialize transfer function */
  for ( idx = 0; idx < nfreq; idx++ )
    {
      xreal[idx] = 1.0e0;
      ximag[idx] = 0.0e0;
    }

  /* Initialize constant, poles and zeros */
  constant = 1.0;
  
  for ( idx = 0; idx < MAXZEROS; idx++ )
    {
      zeros[idx] = dbltocmplx (0.0, 0.0);
    }
  
  for ( idx = 0; idx < MAXPOLES; idx++ )
    {
      poles[idx] = dbltocmplx (0.0, 0.0);
    }
  
  /* Open P&Z file */
  if ( (fp = fopen(pzfilename, "rb")) == NULL )
    {
      fprintf (stderr, "Error opening %s: %s\n",
	       pzfilename, strerror(errno));
      return -1;
    }
  
  /* Read the P&Z file */
  while ( fgets (line, sizeof(line), fp) )
    {
      lp = line;
      linecount++;
      
      /* Skip leading white space */
      while ( isspace(*lp) ) lp++;
      
      if ( *lp == '*' || *lp == '#' ) /* Skip comments lines */
	continue;
      if ( strlen (lp) == 0 ) /* Skip empty lines */
	continue;
      else if ( ! strncasecmp ("ZEROS", lp, 5) )
	{
	  if ( rdzeros )
	    {
	      fprintf (stderr, "%s: Already read a ZEROS line, skipping\n",
		       pzfilename);
	      reading = 0;
	      continue;
	    }
	  if ( (sscanf (lp+5, " %d\n", &nzeros)) != 1 )
	    {
	      fprintf (stderr, "Error parsing ZEROS (line %d): %s\n",
		       linecount, lp);
	      return -1;
	    }
	  if ( nzeros > MAXZEROS )
	    {
	      fprintf (stderr, "Number of ZEROS (%d) cannot be more than %d\n",
		       nzeros, MAXZEROS);
	      return -1;
	    }
	  rdzeros = 1;
	  reading = 1;
	  continue;
	}
      else if ( ! strncasecmp ("POLES", lp, 5) )
	{
	  if ( rdpoles )
	    {
	      fprintf (stderr, "%s: Already read a POLES line, skipping\n",
		       pzfilename);
	      reading = 0;
	      continue;
	    }
	  if ( (sscanf (lp+5, " %d\n", &npoles)) != 1 )
	    {
	      fprintf (stderr, "Error parsing POLES (line %d): %s\n",
		       linecount, lp);
	      return -1;
	    }
	  if ( nzeros > MAXPOLES )
	    {
	      fprintf (stderr, "Number of POLES (%d) cannot be more than %d\n",
		       npoles, MAXPOLES);
	      return -1;
	    }
	  rdpoles = 1;
	  reading = 2;
	  continue;
	}
      else if ( ! strncasecmp ("CONSTANT", lp, 8) )
	{
	  if ( rdconstant )
	    {
	      fprintf (stderr, "%s: Already read a CONSTANT line, skipping\n",
		       pzfilename);
	      reading = 0;
	      continue;
	    }
	  if ( (sscanf (lp+8, " %lf\n", &constant)) != 1 )
	    {
	      fprintf (stderr, "Error parsing POLES (line %d): %s\n",
		       linecount, lp);
	      return -1;
	    }
	  rdconstant = 1;
	  reading = 0;
	  continue;
	}
      else if ( reading == 1 ) /* Reading ZEROS */
	{
	  if ( zeroidx >= MAXZEROS )
	    {
	      fprintf (stderr, "Too many zeros in %s, max: %d\n",
		       pzfilename, MAXZEROS);
	      return -1;
	    }
	  if ( (sscanf (lp, "%lf %lf", &reald, &imagd)) != 2 )
	    {
	      fprintf (stderr, "Error parsing ZERO (line %d): %s\n",
		       linecount, lp);
	      return -1;
	    }
	  
	  zeros[zeroidx++] = dbltocmplx (reald, imagd);
	  continue;
	}
      else if ( reading == 2 ) /* Reading POLES */
	{
	  if ( poleidx >= MAXPOLES )
	    {
	      fprintf (stderr, "Too many poles in %s, max: %d\n",
		       pzfilename, MAXPOLES);
	      return -1;
	    }
	  if ( (sscanf (lp, "%lf %lf", &reald, &imagd)) != 2 )
	    {
	      fprintf (stderr, "Error parsing POLE (line %d): %s\n",
		       linecount, lp);
	      return -1;
	    }
	  
	  poles[poleidx++] = dbltocmplx (reald, imagd);
	  continue;
	}
    } /* Done reading file */
  
  if ( fp )
    fclose (fp);
  
  if ( poleidx != npoles )
    {
      fprintf (stderr, "Read %d POLES but was expecting %d\n", poleidx, npoles);
      return -1;
    }
  
  fprintf (stderr, "Calculating frequency response for poles and zeros\n");
  
  calcfr (nfreq, delfreq, constant, nzeros, zeros, npoles, poles, xreal, ximag);
  
  return 0;
} /* end of function getpzfr() */


/*********************************************************************
 * calcfr():
 *
 * Calculate the frequency response function for specified poles and
 * zeros with a given constant and specified number of frequencies.
 *
 * nfreq   : number of points that will be used in the FFT
 * delfreq : time-interval between data points in the time-domain
 * nzero   : number of zeros
 * zeros   : array of zeros
 * npole   : number of poles
 * poles   : array of poles
 * xreal   : array of real part of frequency response
 * ximag   : array of imaginary part of frequncy response
 *           xreal and ximag must already be allocated with enough
 *           room for nfreq entries.
 *     
 *********************************************************************/
void
calcfr (int nfreq, double delfreq, double constant,
	int nzero, complexd zeros[],
	int npole, complexd poles[],
	double xreal[], double ximag[])
{
  double delomg, fac, omega;
  double ti, ti0, tid, tin;
  double tr, tr0, trd, trn;
  
  int idx, jdx;
  
  delomg = 2 * PI * delfreq;
  
  for ( jdx = 0; jdx < nfreq; jdx++ )
    {
      omega = delomg * jdx;
      trn = 1.0;
      tin = 0.0;
      
      if ( nzero > 0 )
	{
	  for ( idx = 0; idx < nzero; idx++ )
	    {
	      tr = -1.0 * cmplxreal (zeros[idx]);
	      ti = omega - cmplximag (zeros[idx]);
	      tr0 = trn*tr - tin*ti;
	      ti0 = trn*ti + tin*tr;
	      trn = tr0;
	      tin = ti0;
	    }
	}
      
      trd = 1.0;
      tid = 0.0;
      
      if ( npole > 0 )
	{
	  for ( idx = 0; idx < npole; idx++ )
	    {
	      tr = -1.0 * cmplxreal (poles[idx]);
	      ti = omega - cmplximag (poles[idx]);
	      tr0 = trd*tr - tid*ti;
	      ti0 = trd*ti + tid*tr;
	      trd = tr0;
	      tid = ti0;
	    }
	}
      
      fac = constant / (trd*trd + tid*tid);
      xreal[jdx] = fac * (trn*trd + tin*tid);
      ximag[jdx] = fac * (trd*tin - trn*tid);
    }
  
  return;
} /* end of function calcfr() */
