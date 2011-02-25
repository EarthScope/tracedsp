/***************************************************************************
 * tracedsp.c - time-series processor for Mini-SEED and SAC data
 *
 * Opens user specified files, parses the input data, applys
 * processing steps to the timeseries and writes the data.
 *
 * Written by Chad Trabant, IRIS Data Management Center.
 *
 * modified 2011.055
 ***************************************************************************/

// Add stats output.

// check if doubles can be used throughout convolution code
// convert conversions from int to default to doubles

// Add rotation code and options

// Add resampling process

// Add dbdown specification for -Ta

// Add trim/lop segment synchronization process

// Add processing log, a summary line for each operation


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <errno.h>
#include <ctype.h>
#include <math.h>

#include <libmseed.h>

#include "iirfilter.h"
#include "convolve.h"
#include "decimate.h"
#include "rotate.h"
#include "taper.h"
#include "sacformat.h"

#define VERSION "0.9.4dev"
#define PACKAGE "tracedsp"

/* Linkable structure to hold input file names */
struct filelink {
  char *filename;
  int   format;     /* File format: 1 = Mini-SEED, 2 = SAC (alpha or binary) */
  struct filelink *next;
};

/* Linkable structure to hold processing parameters, normally sparse */
struct proclink {
  int    type;
  char  *filename[2];
  int    filetype[2];
  int    respstart;
  int    respstop;
  double lpcutoff;
  int    lporder;
  double hpcutoff;
  int    hporder;
  int    reverseflag;
  double scalefactor;
  int    decimfactor;
  double taperwidth;
  int    tapertype;
  struct proclink *next;
};

#define PROC_STATS       1
#define PROC_LPFILTER    2
#define PROC_HPFILTER    3
#define PROC_BPFILTER    4
#define PROC_CONVOLVE    5
#define PROC_CONVRESP    6 // Subtype of CONVOLVE
#define PROC_DECONVRESP  7 // Subtype of CONVOLVE
#define PROC_CONVSAC     8 // Subtype of CONVOLVE
#define PROC_DECONVSAC   9 // Subtype of CONVOLVE
#define PROC_DIFF2       10
#define PROC_INTTRAP     11
#define PROC_RMEAN       12
#define PROC_SCALE       13
#define PROC_DECIMATE    14
#define PROC_TAPER       15

/* Default order of high/low pass filter */
#define DEFAULT_FILTER_ORDER 4

/* Maximum number of metadata fields per line */
#define MAXMETAFIELDS 17

/* Linkable structure to hold station metadata */
struct metalist {
  char *key;
  char *data;
  struct metalist *next;
};

/* Structure for metadata */
struct metanode {
  char *metafields[MAXMETAFIELDS];
  hptime_t starttime;
  hptime_t endtime;
};

static int procFilter (MSTraceSeg *seg, struct proclink *plp);
static int procConvolve (MSTraceID *id, MSTraceSeg *seg, struct proclink *plp);
static int procDiff2 (MSTraceSeg *seg, struct proclink *plp);
static int procIntTrap (MSTraceSeg *seg, struct proclink *plp);
static int procRMean (MSTraceSeg *seg, struct proclink *plp);
static int procScale (MSTraceSeg *seg, struct proclink *plp);
static int procDecimate (MSTraceSeg *seg, struct proclink *plp);
static int procTaper (MSTraceSeg *seg, struct proclink *plp);

static int64_t readMSEED (char *mseedfile, MSTraceList *mstl);
static int64_t readSAC (char *sacfile, MSTraceList *mstl);
static int parseSAC (FILE *ifp, struct SACHeader *sh, float **data, int format,
                     int verbose, char *sacfile);
static int readBinaryHeaderSAC (FILE *ifp, struct SACHeader *sh, int *format,
				int *swapflag, int verbose, char *sacfile);
static int readBinaryDataSAC (FILE *ifp, float *data, int datacnt,
			      int swapflag, int verbose, char *sacfile);
static int readAlphaHeaderSAC (FILE *ifp, struct SACHeader *sh);
static int readAlphaDataSAC (FILE *ifp, float *data, int datacnt);

static int writeMSEED (MSTraceID *id, MSTraceSeg *seg, char *outputfile);
static int writeSAC (MSTraceID *id, MSTraceSeg *seg, int format, char *outputfile);
static int writeBinarySAC (struct SACHeader *sh, float *fdata, int npts, char *outfile);
static int writeAlphaSAC (struct SACHeader *sh, float *fdata, int npts, char *outfile);
static int insertSACMetaData (struct SACHeader *sh, hptime_t sacstarttime);
static int swapSACHeader (struct SACHeader *sh);

static int differentiate2 (void *input, char inputtype, int length,
			   double rate, void *output);
static int integrateTrap (void *input, char inputtype, int length,
			  double halfstep, void *output);
static int parameterProc (int argcount, char **argvec);
static char *getOptVal (int argcount, char **argvec, int argopt, int dasharg);
static int readMetaData (char *metafile);
static struct metalist *addMetaNode (struct metalist **listroot,
				     void *key, int keylen,
				     void *data, int datalen);
static int addFile (char *filename);
static int addListFile (char *filename);
static void addProcess (int type, char *string1, char *string2, int ivalue1,
			int ivalue2, double dvalue1, double dvalue2);
static void recordHandler (char *record, int reclen, void *vofp);
static void usage (void);

static flag    verbose       = 0;
static flag    basicsum      = 0;    /* Controls printing of basic summary */
static double  timetol       = -1.0; /* Time tolerance for continuous traces */
static double  sampratetol   = -1.0; /* Sample rate tolerance for continuous traces */
static int     reclen        = -1;   /* SEED record length for input data */
static int     packreclen    = -1;   /* SEED record length for output data */
static int     packencoding  = 4;    /* SEED encoding format for output data */
static char   *encodingstr   = 0;    /* SEED encoding format string */
static int     byteorder     = -1;   /* Byte order of output data, use libmseed default */
static int     srateblkt     = 0;    /* Add blockette 100 to Mini-SEED */
static flag    dataformat    = 2;    /* 0 = No output, 1 = Mini-SEED, 2 = SAC */
static flag    sacinformat   = 0;    /* 0=auto, 1=alpha, 2=binary (host), 3=binary (LE), 4=binary (BE) */
static flag    sacoutformat  = 2;    /* 1=alpha, 2=binary (host), 3=binary (LE), 4=binary (BE) */
static char   *sacnet        = 0;    /* SAC network code override */
static char   *sacloc        = 0;    /* SAC location ID override */
static char   *metadatafile  = 0;    /* File containing metadata for output (SAC, etc.) */
static int     prewhiten     = 0;    /* Prewhitening for [de]convolution, predictor order */
static double *spectaperfreq = 0;    /* Spectrum taper frequencies */
static double  lcdBdown      = -1.0; /* Lower corner dB down cutoff */
static double  ucdBdown      = -1.0; /* Upper corner dB down cutoff */
static flag    resptotalsens = 0;    /* Controls evalresp's usage of total sensitivity in RESP */
static flag    respusedelay  = 0;    /* Controls evalresp's usage of estimated delay in RESP */
static flag    respusename   = 1;    /* Controls evalresp's matching of Net, Sta, Loc and Chan */
static char   *respunits     = 0;    /* Controls units for evalresp calculated responses */
static hptime_t starttime    = HPTERROR;
static hptime_t endtime      = HPTERROR;
static char   *outputfile    = 0;    /* Output file name */
static char   *outputdir     = 0;    /* Output base directory */
static int     outputbytes   = 0;    /* Bytes written to output file */
static char   *channel       = 0;    /* Forced output channel/component */

/* Root and tail of input file name list */
struct filelink *filelist = 0;
struct filelink *filelisttail = 0;

/* Root of processing actions list */
struct proclink *proclist = 0;

/* A list of station and coordinates */
struct metalist *metadata = 0;
static int seedinc = 0;

int
main (int argc, char **argv)
{
  struct filelink *flp, *nextflp;
  struct proclink *plp, *nextplp;
  MSTraceList *mstl = NULL;
  MSTraceID *id = NULL;
  MSTraceSeg *seg = NULL;
  char srcname[100];
  int errflag = 0;
  int64_t totalsamps = 0;
  int64_t sampsread;
  int totalfiles = 0;
  
  /* Process given parameters (command line and parameter file) */
  if ( parameterProc (argc, argv) < 0 )
    return -1;
  
  /* Setup encoding environment variable if specified, ugly kludge */
  if ( encodingstr )
    {
      int inputencoding = strtoul (encodingstr, NULL, 10);
      
      if ( inputencoding == 0 && errno == EINVAL )
        {
          ms_log (2, "Error parsing input encoding format: %s\n", encodingstr);
          return -1;
        }
      
      MS_UNPACKENCODINGFORMAT (inputencoding);
    }
  
  /* Initialize MSTraceList */
  mstl = mstl_init (NULL);
  
  /* Loop over input file list */
  flp = filelist;
  while ( flp != 0 )
    {
      /* Read Mini-SEED file */
      if ( flp->format == 1 )
	{
	  if ( verbose >= 2 )
	    fprintf (stderr, "Processing Mini-SEED: %s\n", flp->filename);
	  
	  sampsread = readMSEED (flp->filename, mstl);
	  
	  if ( sampsread < 0 )
	    {
	      fprintf (stderr, "Error reading '%s', skipping\n", flp->filename);
	      flp = flp->next;
	      continue;
	    }
	  
	  totalsamps += sampsread;
	}
      /* Read SAC file */
      else
	{
	  if ( verbose >= 2 )
	    fprintf (stderr, "Processing SAC: %s\n", flp->filename);
	  
	  sampsread = readSAC (flp->filename, mstl);
	  
	  if ( sampsread < 0 )
	    {
	      fprintf (stderr, "Error reading '%s', skipping\n", flp->filename);
	      flp = flp->next;
	      continue;
	    }
	  
	  totalsamps += sampsread;
	}
      
      totalfiles++;
      flp = flp->next;
    } /* End of looping over file list */
  
  if ( basicsum )
    printf ("Input Files: %d, Samples: %lld\n", totalfiles, (long long int) totalsamps);
  
  /* Loop through all read MSTraceIDs, apply processing and write resulting data */
  id = mstl->traces;
  while ( id )
    {
      seg = id->first;
      while ( seg )
	{
	  /* Loop through process list */
	  plp = proclist;
	  while ( plp && ! errflag )
	    {
	      if ( plp->type == PROC_LPFILTER || plp->type == PROC_HPFILTER || plp->type == PROC_BPFILTER )
		{
		  if ( procFilter (seg, plp) )
		    {
		      fprintf (stderr, "Error applying filter for %s\n", srcname);
		      errflag = -1;
		      break;
		    }
		}
	      else if ( plp->type == PROC_CONVOLVE )
		{
		  if ( procConvolve (id, seg, plp) )
		    {
		      fprintf (stderr, "Error (de)convolving response (%s - %s) from %s\n",
			       (plp->filename[0])?plp->filename[0]:"None",
			       (plp->filename[1])?plp->filename[1]:"None", srcname);
		      errflag = -1;
		      break;
		    }
		}
	      else if ( plp->type == PROC_DIFF2 )
		{
		  if ( procDiff2 (seg, plp) )
		    {
		      fprintf (stderr, "Error differentiating time-series for %s\n",
			       srcname);
		      errflag = -1;
		      break;
		    }
		}
	      else if ( plp->type == PROC_INTTRAP )
		{
		  if ( procIntTrap (seg, plp) )
		    {
		      fprintf (stderr, "Error integrating time-series for %s\n",
			       srcname);
		      errflag = -1;
		      break;
		    }
		}
	      else if ( plp->type == PROC_RMEAN )
		{
		  if ( procRMean (seg, plp) )
		    {
		      fprintf (stderr, "Error removing mean from time-series %s\n",
			       srcname);
		      errflag = -1;
		      break;
		    }
		}
	      else if ( plp->type == PROC_SCALE )
		{
		  if ( procScale (seg, plp) )
		    {
		      fprintf (stderr, "Error scaling values of time-series %s\n",
			       srcname);
		      errflag = -1;
		      break;
		    }
		}
	      else if ( plp->type == PROC_DECIMATE )
		{
		  if ( procDecimate (seg, plp) )
		    {
		      fprintf (stderr, "Error decimating time-series %s\n",
			       srcname);
		      errflag = -1;
		      break;
		    }
		}
	      else if ( plp->type == PROC_TAPER )
		{
		  if ( procTaper (seg, plp) )
		    {
		      fprintf (stderr, "Error tapering time-series %s\n",
			       srcname);
		      errflag = -1;
		      break;
		    }
		}
	      
	      plp = plp->next;
	    }
	  
	  /* If no processing errors occurred write the data out */
	  if ( ! errflag )
	    {
	      /* Write Mini-SEED data */
	      if ( dataformat == 1 )
		{
		  if ( writeMSEED (id, seg, outputfile) < 0 )
		    fprintf (stderr, "Error writing Mini-SEED\n");
		}
	      /* Write SAC data */
	      else if ( dataformat == 2 )
		{
		  if ( writeSAC (id, seg, sacoutformat, outputfile) < 0 )
		    fprintf (stderr, "Error writing SAC\n");
		}
	    }
	  
	  seg = seg->next;
	}
      
      id = id->next;
    }
  
  /* Make sure everything is cleaned up */
  mstl_free (&mstl, 1);
  
  flp = filelist;
  while ( flp )
    {
      nextflp = flp->next;
      if ( flp->filename )
	free (flp->filename);
      free (flp);
      flp = nextflp;
    }
  
  plp = proclist;
  while ( plp )
    {
      nextplp = plp->next;
      if ( plp->filename[0] )
	free (plp->filename[0]);
      if ( plp->filename[1] )
	free (plp->filename[1]);
      free (plp);
      plp = nextplp;
    }
  
  return errflag;
}  /* End of main() */


/***************************************************************************
 * procFilter:
 * 
 * Apply filter to the MSTraceSeg, results are requested to be
 * returned as floats.
 *
 * Returns 0 on success and -1 on error.
 ***************************************************************************/
static int
procFilter (MSTraceSeg *seg, struct proclink *plp)
{
  void *datasamples = 0;
  
  if ( ! seg || ! plp )
    return -1;
  
  /* Filter trace */
  if ( plp->lporder || plp->hporder )
    {
      if ( verbose )
	{
	  if ( plp->lporder && ! plp->hporder )
	    fprintf (stderr, "Low-pass filter cutoff: %f, order: %d\n",
		     plp->lpcutoff, plp->lporder);
	  else if ( ! plp->lporder && plp->hporder )
	    fprintf (stderr, "High-pass filter cutoff: %f, order: %d\n",
		     plp->hpcutoff, plp->hporder);
	  else
	    fprintf (stderr, "Band-pass filter HP cutoff: %f, order: %d => LP cutoff: %f, order: %d\n",
		     plp->lpcutoff, plp->lporder, plp->hpcutoff, plp->hporder);
	}
      
      /* Apply the filter */
      if ( iirfilter (seg->datasamples, seg->sampletype, seg->numsamples, plp->reverseflag,
		      &datasamples, 'f', plp->hporder, plp->hpcutoff,
		      plp->lporder, plp->lpcutoff, seg->samprate, verbose) )
	{
	  if ( datasamples )
	    free (datasamples);
	  return -1;
	}
      else
	{
	  /* Free the original buffer and replace it with the filtered */
	  if ( seg->datasamples )
	    free (seg->datasamples);
	  
	  seg->datasamples = datasamples;
	  seg->sampletype = 'f';
	}
    }
  
  return 0;
}  /* End of procFilter() */


/***************************************************************************
 * procConvolve:
 * 
 * Convolve or deconvolve supplied response(s) from signal.  The
 * convolution routines require floats so this routine will always
 * convert the supplied data to floats and return the results as
 * floats.
 *
 * Returns 0 on success and non-zero on error.
 ***************************************************************************/
static int
procConvolve (MSTraceID *id, MSTraceSeg *seg, struct proclink *plp)
{
  float *fdata = 0;
  int idx;
  int retval = 0;
  int nfft;
  int nfreqs;
  double delfreq;
  
  double *creal = NULL;
  double *cimag = NULL;
  double *dreal = NULL;
  double *dimag = NULL;
  double *xreal = NULL;
  double *ximag = NULL;
  
  if ( ! seg || ! plp )
    return -1;
  
  /* Convert samples to floats if needed */
  if ( seg->sampletype == 'f' )
    {
      fdata = seg->datasamples;
    }
  else
    {
      if ( (fdata = (float *) malloc (seg->numsamples * sizeof(float))) == NULL )
	{
	  fprintf (stderr, "Error allocating memory\n");
	  return -1;
	}
      
      if ( seg->sampletype == 'i' )
	{
	  int32_t *iptr = seg->datasamples;
	  for (idx = 0; idx < seg->numsamples; idx++)
	    fdata[idx] = (float) iptr[idx];
	  
	  free (seg->datasamples);
	  seg->datasamples = fdata;
	  seg->sampletype = 'f';
	}
      else if ( seg->sampletype == 'd' )
	{
	  double *dptr = seg->datasamples;
	  for (idx = 0; idx < seg->numsamples; idx++)
	    fdata[idx] = (float) dptr[idx];
	  
	  free (seg->datasamples);
	  seg->datasamples = fdata;
	  seg->sampletype = 'f';
	}
      else
	{
	  fprintf (stderr, "procConvolve(): unsupported sample type: %c\n",
		   seg->sampletype);
	  if ( fdata )
	    free (fdata);
	  return -1;
	}
    }
  
  if ( verbose >= 1 && plp->filetype[0] && plp->filetype[1] )
    {
      fprintf (stderr, "Deconvolution-convolution have been paired into a transfer operation\n");
    }
  
  /* Calculate common parameters */
  nfft = next2 (seg->numsamples);
  nfreqs = nfft/2 + 1;
  delfreq = seg->samprate / nfft;
  
  /* Loop over two potential response files */
  for (idx = 0; idx < 2; idx++)
    {
      if ( plp->filename[idx] == NULL )
	break;
      
      /* Calculate response functions from SEED RESP */
      if ( plp->filetype[idx] == PROC_CONVRESP || plp->filetype[idx] == PROC_DECONVRESP )
	{
	  if ( verbose )
	    {
	      if ( plp->filetype[idx] == PROC_CONVRESP )
		fprintf (stderr, "Convolving with SEED RESP response '%s'\n",
			 plp->filename[idx]);
	      else
		fprintf (stderr, "Deconvolving SEED RESP response '%s'\n",
			 plp->filename[idx]);
	    }
	  
	  if ( respusename )
	    {
	      retval = calcfr_resp (nfreqs, delfreq,
				    id->network, id->station, id->location, id->channel,
				    plp->respstart, plp->respstop,
				    respunits, MS_HPTIME2EPOCH (seg->starttime), respusedelay,
				    plp->filename[idx], resptotalsens,
				    &xreal, &ximag, verbose);
	    }
	  else
	    {
	      retval = calcfr_resp (nfreqs, delfreq, "*", "*", "*", "*",
				    plp->respstart, plp->respstop,
				    respunits, MS_HPTIME2EPOCH (seg->starttime), respusedelay,
				    plp->filename[idx], resptotalsens,
				    &xreal, &ximag, verbose);
	    }
	  
	  /* Assign to convolution or deconvolution array */
	  if ( plp->filetype[idx] == PROC_CONVRESP )
	    {
	      creal = xreal; xreal = NULL;
	      cimag = ximag; ximag = NULL;
	    }
	  else
	    {
	      dreal = xreal; xreal = NULL;
	      dimag = ximag; ximag = NULL;	  
	    }
	}
      /* Calculate response functions from SAC P&Zs */
      else if ( plp->filetype[idx] == PROC_CONVSAC || plp->filetype[idx] == PROC_DECONVSAC )
	{
	  if ( verbose )
	    {
	      if ( plp->filetype[idx] == PROC_CONVSAC )
		fprintf (stderr, "Convolving with SAC Poles & Zeros response '%s'\n",
			 plp->filename[idx]);
	      else
		fprintf (stderr, "Deconvolving SAC Poles & Zeros response '%s'\n",
			 plp->filename[idx]);
	    }
	  
	  retval = calcfr_sac (nfreqs, delfreq, plp->filename[idx],
			       &xreal, &ximag , verbose);
	  
	  /* Assign to convolution or deconvolution array */
	  if ( plp->filetype[idx] == PROC_CONVSAC )
	    {
	      creal = xreal; xreal = NULL;
	      cimag = ximag; ximag = NULL;
	    }
	  else
	    {
	      dreal = xreal; xreal = NULL;
	      dimag = ximag; ximag = NULL;
	    }
	}
    }
  
  /* Check if tapering parameters need to be calculated for deconvolution */
  if ( dreal && dimag && spectaperfreq &&
       (spectaperfreq[0] == -1.0 || spectaperfreq[1] == -1.0 ||
	spectaperfreq[2] == -1.0 || spectaperfreq[3] == -1.0) )
    {
      /* Determine taper parameters for deconvolution response */
      if ( findtaper (spectaperfreq, dreal, dimag, nfreqs, delfreq, lcdBdown, ucdBdown) )
	{
	  fprintf (stderr, "Error determing spectral taper parameters\n");
	  return -1;
	}
      
      if ( verbose )
	fprintf (stderr, "Final spectral tapering (Hz): %g/%g => %g/%g [cutoffs %g/%g]\n",
		 spectaperfreq[0], spectaperfreq[1], spectaperfreq[2], spectaperfreq[3],
		 lcdBdown, ucdBdown);
    }
  
  /* Perform convolution, deconvolution or both */
  retval = convolve (fdata, seg->numsamples, 1.0/seg->samprate, nfreqs, nfft,
		     creal, cimag, dreal, dimag, spectaperfreq, &prewhiten, verbose);
  
  /* Free response function arrays */
  if ( creal )
    free (creal);
  if ( cimag )
    free (cimag);
  if ( dreal )
    free (dreal);
  if ( dimag )
    free (dimag);
  
  return retval;
}  /* End of procConvolve() */


/***************************************************************************
 * procDiff2:
 * 
 * Prepare for an uncentered, 2-point differentiation.
 *
 * Returns 0 on success and non-zero on error.
 ***************************************************************************/
static int
procDiff2 (MSTraceSeg *seg, struct proclink *plp)
{
  int idx;
  int count;
  
  if ( ! seg || ! plp )
    return -1;
  
  /* Convert integer samples to floats in-place */
  if ( seg->sampletype == 'i' )
    {
      int32_t *iptr = seg->datasamples;
      float *fdata = seg->datasamples;
      
      for (idx = 0; idx < seg->numsamples; idx++)
	fdata[idx] = (float) iptr[idx];
      
      seg->sampletype = 'f';
    }
  else if ( seg->sampletype != 'f' && seg->sampletype != 'd' )
    {
      fprintf (stderr, "procDiff2(): unsupported sample type: %c\n",
	       seg->sampletype);
      return -1;
    }
  
  /* Perform differentiation */
  if ( verbose )
    fprintf (stderr, "Differentiating (%c) time-series\n", seg->sampletype);
  
  count = differentiate2 (seg->datasamples, seg->sampletype, seg->numsamples,
			  seg->samprate, seg->datasamples);
  
  if ( count < 0 )
    return -1;
  
  /* Update sample counts */
  seg->samplecnt = count;
  seg->numsamples = count;
  
  /* Shift the start time by 1/2 the sample interval */
  seg->starttime += (0.5 / seg->samprate) * HPTMODULUS;
  
  return 0;
}  /* End of procDiff2() */


/***************************************************************************
 * procIntTrap:
 * 
 * Prepare for integration using the trapezoidal (midpoint) method.
 *
 * Returns 0 on success and non-zero on error.
 ***************************************************************************/
static int
procIntTrap (MSTraceSeg *seg, struct proclink *plp)
{
  int idx;
  int count;
  
  if ( ! seg || ! plp )
    return -1;
  
  /* Convert integer samples to floats in-place */
  if ( seg->sampletype == 'i' )
    {
      int32_t *iptr = seg->datasamples;
      float *fdata = seg->datasamples;
      
      for (idx = 0; idx < seg->numsamples; idx++)
	fdata[idx] = (float) iptr[idx];
      
      seg->sampletype = 'f';
    }
  else if ( seg->sampletype != 'f' && seg->sampletype != 'd' )
    {
      fprintf (stderr, "procIntTrap(): unsupported sample type: %c\n",
	       seg->sampletype);
      return -1;
    }
  
  /* Perform integration */
  if ( verbose )
    fprintf (stderr, "Integrating (%c) time-series\n", seg->sampletype);
  
  count = integrateTrap (seg->datasamples, seg->sampletype, seg->numsamples,
			 (0.5/seg->samprate), seg->datasamples);
  
  if ( count < 0 )
    return -1;
  
  /* Update sample counts */
  seg->samplecnt = count;
  seg->numsamples = count;
  
  /* Shift the start time by 1/2 the sample interval */
  seg->starttime += (0.5 / seg->samprate) * HPTMODULUS;
  
  return 0;
}  /* End of procIntTrap() */


/***************************************************************************
 * procRMean:
 * 
 * Removes the mean from a time-series.  Integer samples will be
 * converted to floats.
 *
 * Returns 0 on success and non-zero on error.
 ***************************************************************************/
static int
procRMean (MSTraceSeg *seg, struct proclink *plp)
{
  int idx;
  double Mean, pM;
  int32_t *idata = seg->datasamples;
  float *fdata = seg->datasamples;
  double *ddata = seg->datasamples;
  
  if ( ! seg || ! plp )
    return -1;
  
  /* Convert integer samples to floats in-place */
  if ( seg->sampletype == 'i' )
    {
      for (idx = 0; idx < seg->numsamples; idx++)
	fdata[idx] = (float) idata[idx];
      
      seg->sampletype = 'f';
    }
  else if ( seg->sampletype != 'f' && seg->sampletype != 'd' )
    {
      fprintf (stderr, "procRMean(): unsupported sample type: %c\n",
	       seg->sampletype);
      return -1;
    }
  
  /* Find mean value */
  if ( seg->sampletype == 'f' )
    {
      /* Calculate running mean */
      pM = Mean = *fdata;
      for (idx = 1; idx < seg->numsamples; idx++)
	{
	  Mean = pM + (*(fdata+idx) - pM) / (idx + 1);
	  pM = Mean;
	}
      
      if ( verbose )
	fprintf (stderr, "Removing mean of %g from time-series\n", Mean);
      
      /* Remove mean */
      for (idx = 0; idx < seg->numsamples; idx++)
	{
	  *(fdata+idx) -= Mean;
	}
    }
  else if ( seg->sampletype == 'd' )
    {
      /* Calculate running mean */
      pM = Mean = *ddata;
      for (idx = 1; idx < seg->numsamples; idx++)
	{
	  Mean = pM + (*(ddata+idx) - pM) / (idx + 1);
	  pM = Mean;
	}
      
      if ( verbose )
	fprintf (stderr, "Removing mean of %g from time-series\n", Mean);
      
      /* Remove mean */
      for (idx = 0; idx < seg->numsamples; idx++)
	{
	  *(ddata+idx) -= Mean;
	}
    }
  
  return 0;
}  /* End of procRMean() */


/***************************************************************************
 * procScale:
 *
 * Scales all data samples in a time-series.  Integer samples will be
 * converted to floats.
 *
 * Returns 0 on success and non-zero on error.
 ***************************************************************************/
static int
procScale (MSTraceSeg *seg, struct proclink *plp)
{
  int idx;
  int32_t *idata = seg->datasamples;
  float *fdata = seg->datasamples;
  double *ddata = seg->datasamples;
  
  if ( ! seg || ! plp )
    return -1;
  
  /* Convert integer samples to floats in-place */
  if ( seg->sampletype == 'i' )
    {
      for (idx = 0; idx < seg->numsamples; idx++)
	fdata[idx] = (float) idata[idx];
      
      seg->sampletype = 'f';
    }
  else if ( seg->sampletype != 'f' && seg->sampletype != 'd' )
    {
      fprintf (stderr, "procScale(): unsupported sample type: %c\n",
	       seg->sampletype);
      return -1;
    }
  
  /* Scale sample values */
  if ( seg->sampletype == 'f' )
    {
      if ( verbose )
	fprintf (stderr, "Scaling time-series by %g\n", plp->scalefactor);
      
      /* Scale samples */
      for (idx = 0; idx < seg->numsamples; idx++)
	{
	  *(fdata+idx) *= plp->scalefactor;
	}
    }
  else if ( seg->sampletype == 'd' )
    {
      if ( verbose )
	fprintf (stderr, "Scaling time-series by %g\n", plp->scalefactor);
      
      /* Scale samples */
      for (idx = 0; idx < seg->numsamples; idx++)
	{
	  *(ddata+idx) *= plp->scalefactor;
	}
    }
  
  return 0;
}  /* End of procScale() */


/***************************************************************************
 * procDecimate:
 *
 * Decimates the time-series by a given factor.  Integer and float
 * samples will be converted to doubles.
 *
 * Returns 0 on success and non-zero on error.
 ***************************************************************************/
static int
procDecimate (MSTraceSeg *seg, struct proclink *plp)
{
  int idx;
  int numsamples;
  int32_t *idata = seg->datasamples;
  float *fdata = seg->datasamples;
  double *ddata = seg->datasamples;
  
  if ( ! seg || ! plp )
    return -1;
  
  /* Convert integer and float samples to doubles */
  if ( seg->sampletype == 'i' || seg->sampletype == 'f' )
    {
      /* Allocate memory for double samples */
      if ( ! (ddata = (double *) malloc (seg->numsamples * sizeof(double))) )
	{
	  fprintf (stderr, "procDecimate(): Cannot allocate memory\n");
	  return -1;
	}
      
      /* Convert samples to doubles */
      if ( seg->sampletype == 'i' )
	for (idx = 0; idx < seg->numsamples; idx++)
	  ddata[idx] = (double) idata[idx];
      else
	for (idx = 0; idx < seg->numsamples; idx++)
	  ddata[idx] = (double) fdata[idx];
      
      free (seg->datasamples);
      seg->datasamples = ddata;
      seg->sampletype = 'd';
    }
  
  if ( verbose )
    fprintf (stderr, "Decimating time-series by a factor of %d (%g -> %g sps)\n",
	     plp->decimfactor, seg->samprate, seg->samprate/plp->decimfactor);
  
  /* Perform the decimation and filtering */
  numsamples = decimate (seg->datasamples, seg->numsamples, plp->decimfactor, NULL, -1, -1);
  
  if ( numsamples >= 0 )
    {
      /* Adjust sample rate, sample count, end time and sample buffer size */
      seg->samprate /= plp->decimfactor;
      seg->samplecnt = numsamples;
      seg->numsamples = numsamples;
      seg->endtime = seg->starttime +
	(((double)(seg->numsamples-1)/seg->samprate * HPTMODULUS) + 0.5);
      
      if ( ! (seg->datasamples = realloc (seg->datasamples, numsamples*sizeof(double))) )
	{
	  fprintf (stderr, "procDecimate(): Error reallocating sample buffer\n");
	  return -1;
	}
    }
  else
    {
      fprintf (stderr, "procDecimate(): Error decimating time-series\n");
      return -1;      
    }
  
  return 0;
}  /* End of procDecimate() */


/***************************************************************************
 * procTaper:
 *
 * Tapers the time-series using a specified type and width.  Integer
 * and float samples will be converted to doubles.
 *
 * Returns 0 on success and non-zero on error.
 ***************************************************************************/
static int
procTaper (MSTraceSeg *seg, struct proclink *plp)
{
  int idx;
  int retval;
  int32_t *idata = seg->datasamples;
  float *fdata = seg->datasamples;
  double *ddata = seg->datasamples;
  
  if ( ! seg || ! plp )
    return -1;
  
  /* Convert integer and float samples to doubles */
  if ( seg->sampletype == 'i' || seg->sampletype == 'f' )
    {
      /* Allocate memory for double samples */
      if ( ! (ddata = (double *) malloc (seg->numsamples * sizeof(double))) )
	{
	  fprintf (stderr, "procTaper(): Cannot allocate memory\n");
	  return -1;
	}
      
      /* Convert samples to doubles */
      if ( seg->sampletype == 'i' )
	for (idx = 0; idx < seg->numsamples; idx++)
	  ddata[idx] = (double) idata[idx];
      else
	for (idx = 0; idx < seg->numsamples; idx++)
	  ddata[idx] = (double) fdata[idx];
      
      free (seg->datasamples);
      seg->datasamples = ddata;
      seg->sampletype = 'd';
    }
  
  if ( seg->sampletype != 'd' )
    {
      fprintf (stderr, "procTaper(): Unrecognized sample type: '%c'\n", seg->sampletype);
      return -1;
    }
  
  if ( verbose )
    {
      char *typestr = "";
      switch ( plp->tapertype )
	{
	case TAPER_HANNING: typestr = "Hanning"; break;
	case TAPER_HAMMING: typestr = "Hamming"; break;
	case TAPER_COSINE: typestr = "Cosine"; break;
	}
      
      fprintf (stderr, "Tapering time-series using width %g (%s)\n",
	       plp->taperwidth, typestr);
    }
  
  /* Perform the tapering */
  retval = taper (seg->datasamples, seg->numsamples, plp->taperwidth, plp->tapertype);
  
  if ( retval < 0 )
    {
      fprintf (stderr, "procTaper(): Error tapering time-series\n");
      return -1;      
    }
  
  if ( verbose )
    fprintf (stderr, "Taper window length: %d samples\n", retval);
  
  return 0;
}  /* End of procTaper() */


/***************************************************************************
 * readMSEED:
 * 
 * Read file containing Mini-SEED and add data to the supplied MSTraceList.
 *
 * Returns the number of samples read on success and -1 on error.
 ***************************************************************************/
static int64_t
readMSEED (char *mseedfile, MSTraceList *mstl)
{
  MSRecord *msr = 0;
  char srcname[100];
  int64_t totalsamps = 0;
  int retcode = MS_NOERROR;
  
  if ( ! mseedfile )
    {
      fprintf (stderr, "readMSEED(): No input file specified\n");
      return -1;
    }
  
  if ( ! mstl )
    {
      fprintf (stderr, "readMSEED(): No MSTraceList specified\n");
      return -1;
    }
  
  /* Loop over the input file reading records */
  while ( (retcode = ms_readmsr (&msr, mseedfile, reclen, NULL, NULL, 1, 1, verbose-3)) == MS_NOERROR )
    {
      /* Skip data records that do not contain time-series data */
      if ( msr->numsamples == 0 || (msr->sampletype != 'i' && msr->sampletype != 'f' && msr->sampletype != 'd') )
	{
	  if ( verbose >= 3 )
	    {
	      char stime[100];
	      msr_srcname (msr, srcname, 1);
	      ms_hptime2seedtimestr (msr->starttime, stime, 1);
	      fprintf (stderr, "Skipping (no time-series data) %s, %s\n", srcname, stime);
	    }
	  continue;
	}
      
      /* Check if record matches start/end time criteria */
      if ( starttime != HPTERROR && (msr->starttime < starttime) )
	{
	  if ( verbose >= 3 )
	    {
	      char stime[100];
	      msr_srcname (msr, srcname, 1);
	      ms_hptime2seedtimestr (msr->starttime, stime, 1);
	      fprintf (stderr, "Skipping (start time) %s, %s\n", srcname, stime);
	    }
	  continue;
	}
      
      if ( endtime != HPTERROR && (msr_endtime(msr) > endtime) )
	{
	  if ( verbose >= 3 )
	    {
	      char stime[100];
	      msr_srcname (msr, srcname, 1);
	      ms_hptime2seedtimestr (msr->starttime, stime, 1);
	      fprintf (stderr, "Skipping (end time) %s, %s\n", srcname, stime);
	    }
	  continue;
	}
      
      totalsamps += msr->samplecnt;
      
      if ( verbose >= 3 )
	msr_print (msr, verbose - 3);
      
      /* Add to the MSTraceList */
      mstl_addmsr (mstl, msr, 0, 1, timetol, sampratetol);
    }
  
  /* Print error if not EOF and not counting down records */
  if ( retcode != MS_ENDOFFILE )
    {
      fprintf (stderr, "Error reading %s: %s\n", mseedfile, ms_errorstr(retcode));
      totalsamps = -1;
    }
  
  /* Make sure everything is cleaned up */
  ms_readmsr (&msr, NULL, 0, NULL, NULL, 0, 0, 0);
  
  return totalsamps;
}  /* End of readMSEED() */


/***************************************************************************
 * readSAC:
 *
 * Read a SAC file and add data samples to a MSTraceGroup.  As the SAC
 * data is read in a MSRecord struct is used as a holder for the input
 * information.
 *
 * The SAC header contents is stored at the private pointer of the
 * MSTraceSeg structure (MSTraceSeg->prvtptr).
 *
 * Returns the number of samples read on success and -1 on failure
 ***************************************************************************/
static int64_t
readSAC (char *sacfile, MSTraceList *mstl)
{
  FILE *ifp = 0;
  MSRecord *msr = 0;
  MSTraceSeg *seg;
  
  struct SACHeader sh;
  float *fdata = 0;
  int datacnt;
  
  /* Open input file */
  if ( (ifp = fopen (sacfile, "rb")) == NULL )
    {
      fprintf (stderr, "Cannot open input file: %s (%s)\n",
	       sacfile, strerror(errno));
      return -1;
    }
  
  /* Parse input SAC file into a header structure and data buffer */
  if ( (datacnt = parseSAC (ifp, &sh, &fdata, sacinformat, verbose, sacfile)) < 0 )
    {
      fprintf (stderr, "Error parsing %s\n", sacfile);
      return -1;
    }
  
  if ( ! (msr = msr_init(msr)) )
    {
      fprintf (stderr, "Cannot initialize MSRecord strcture\n");
      return -1;
    }
  
  /* Populate MSRecord structure with header details */
  if ( strncmp (SUNDEF, sh.knetwk, 8) ) ms_strncpclean (msr->network, sh.knetwk, 2);
  if ( strncmp (SUNDEF, sh.kstnm, 8) ) ms_strncpclean (msr->station, sh.kstnm, 5);
  if ( strncmp (SUNDEF, sh.khole, 8) ) ms_strncpclean (msr->location, sh.khole, 2);
  if ( strncmp (SUNDEF, sh.kcmpnm, 8) ) ms_strncpclean (msr->channel, sh.kcmpnm, 3);
  
  if ( sacnet )
    ms_strncpclean (msr->network, sacnet, 2);
  
  if ( sacloc )
    ms_strncpclean (msr->location, sacloc, 2);
  
  msr->starttime = ms_time2hptime (sh.nzyear, sh.nzjday, sh.nzhour, sh.nzmin, sh.nzsec, sh.nzmsec * 1000);
  
  /* Adjust for Begin ('B' SAC variable) time offset */
  msr->starttime += (double) sh.b * HPTMODULUS;
  
  /* Calculate sample rate from interval(period) rounding to nearest 0.000001 Hz */
  msr->samprate = (double) ((int)((1 / sh.delta) * 100000 + 0.5)) / 100000;
  
  msr->samplecnt = msr->numsamples = datacnt;
  
  msr->sampletype = 'f';
  msr->datasamples = fdata;
  
  if ( verbose >= 1 )
    {
      fprintf (stderr, "[%s] %d samps @ %.6f Hz for N: '%s', S: '%s', L: '%s', C: '%s'\n",
	       sacfile, msr->numsamples, msr->samprate,
	       msr->network, msr->station,  msr->location, msr->channel);
    }
  
  /* Add new data to MSTraceList, do not merge with other segments */
  if ( ! (seg = mstl_addmsr (mstl, msr, 0, 0, timetol, sampratetol)) )
    {
      fprintf (stderr, "[%s] Error adding samples to MSTraceGroup\n", sacfile);
    }
  
  /* Store SAC header structure in private pointer of MSTraceSeg */
  if ( ! seg->prvtptr )
    {
      if ( ! (seg->prvtptr = malloc (sizeof(struct SACHeader))) )
	{
	  fprintf (stderr, "Error allocating memory for SACHheader\n");
	}
      else
	{
	  memcpy (seg->prvtptr, &sh, sizeof(struct SACHeader));
	}
    }
  else if ( verbose >= 1 )
    {
      fprintf (stderr, "%s_%s_%s_%s: SAC header contents not stored for merged segment",
	       msr->network, msr->station,  msr->location, msr->channel);
    }
  
  fclose (ifp);
  
  if ( fdata )
    free (fdata);
  
  msr->datasamples = 0;
  
  if ( msr )
    msr_free (&msr);
  
  return datacnt;
}  /* End of readSAC() */


/***************************************************************************
 * parseSAC:
 *
 * Parse a SAC file, autodetecting format dialect (ALPHA,
 * binary, big or little endian).  Results will be placed in the
 * supplied SAC header struct and data (float sample array in host
 * byte order).  The data array will be allocated by this routine and
 * must be free'd by the caller.  The data array will contain the
 * number of samples indicated in the SAC header (sh->npts).
 *
 * The format argument is interpreted as:
 * 0 : Unknown, detection needed
 * 1 : ALPHA
 * 2 : Binary, byte order detection needed
 * 3 : Binary, little endian
 * 4 : Binary, big endian
 *
 * Returns number of data samples in file or -1 on failure.
 ***************************************************************************/
static int
parseSAC (FILE *ifp, struct SACHeader *sh, float **data, int format,
	  int verbose, char *sacfile)
{
  char fourc[4];
  int swapflag = 0;
  int rv;
  
  /* Argument sanity */
  if ( ! ifp || ! sh || ! data )
    return -1;
  
  /* Read the first 4 characters */
  if ( fread (&fourc, 4, 1, ifp) < 1 )
    return -1;
  
  /* Determine if the file is ALPHA or binary SAC,
   * if the first 4 characters are spaces assume ALPHA SAC */
  if ( format == 0 )
    {
      if ( fourc[0] == ' ' && fourc[1] == ' ' && fourc[2] == ' ' && fourc[3] == ' ' )
	format = 1;
      else
	format = 2;  /* Byte order detection will occur below */
    }
  
  /* Rewind the file position pointer to the beginning */
  rewind (ifp);
  
  
  /* Read the header */
  if ( format == 1 )  /* Process SAC ALPHA header */
    {
      if ( (rv = readAlphaHeaderSAC (ifp, sh)) )
	{
	  fprintf (stderr, "[%s] Error parsing SAC ALPHA header at line %d\n",
		   sacfile, rv);
	  return -1;
	}
    }
  else if ( format >= 2 && format <= 4 ) /* Process SAC binary header */
    {
      if ( readBinaryHeaderSAC (ifp, sh, &format, &swapflag, verbose, sacfile) )
	{
	  fprintf (stderr, "[%s] Error parsing SAC header\n", sacfile);
	  return -1;
	}
    }
  else
    {
      fprintf (stderr, "[%s] Unrecognized format value: %d\n", sacfile, format);
      return -1;
    }
  
  /* Sanity check the start time */
  if ( sh->nzyear < 1900 || sh->nzyear >3000 ||
       sh->nzjday < 1 || sh->nzjday > 366 ||
       sh->nzhour < 0 || sh->nzhour > 23 ||
       sh->nzmin < 0 || sh->nzmin > 59 ||
       sh->nzsec < 0 || sh->nzsec > 60 ||
       sh->nzmsec < 0 || sh->nzmsec > 999999 )
    {
      fprintf (stderr, "[%s] Unrecognized format (not SAC?)\n", sacfile);
      return -1;
    }
  
  if ( verbose )
    {
      if ( format == 1 )
	fprintf (stderr, "[%s] Reading SAC ALPHA format\n", sacfile);
      if ( format == 3 )
	fprintf (stderr, "[%s] Reading SAC binary format (little-endian)\n", sacfile);
      if ( format == 4 )
	fprintf (stderr, "[%s] Reading SAC binary format (big-endian)\n", sacfile);
    }
  
  if ( verbose > 2 )
    fprintf (stderr, "[%s] SAC header version number: %d\n", sacfile, sh->nvhdr);
  
  if ( sh->nvhdr != 6 )
    fprintf (stderr, "[%s] WARNING SAC header version (%d) not expected value of 6\n",
	     sacfile, sh->nvhdr);
  
  if ( sh->npts <= 0 )
    {
      fprintf (stderr, "[%s] No data, number of samples: %d\n", sacfile, sh->npts);
      return -1;
    }
  
  if ( sh->iftype != ITIME )
    {
      fprintf (stderr, "[%s] Data is not time series (IFTYPE=%d), cannot convert other types\n",
	       sacfile, sh->iftype);
      return -1;
    }
  
  if ( ! sh->leven )
    {
      fprintf (stderr, "[%s] Data is not evenly spaced (LEVEN not true), cannot convert\n", sacfile);
      return -1;
    }
  
  /* Allocate space for data samples */
  if ( ! (*data = (float *) calloc (1, sizeof(float) * sh->npts)) )
    {
      fprintf (stderr, "Error allocating memory for data samples\n");
      return -1;
    }
  
  /* Read the data samples */
  if ( format == 1 )  /* Process SAC ALPHA data */
    {
      if ( (rv = readAlphaDataSAC (ifp, *data, sh->npts)) )
	{
	  fprintf (stderr, "[%s] Error parsing SAC ALPHA data at line %d\n",
		   sacfile, rv);
	  return -1;
	}
    }
  else if ( format >= 2 && format <= 4 ) /* Process SAC binary data */
    {
      if ( readBinaryDataSAC (ifp, *data, sh->npts, swapflag, verbose, sacfile) )
	{
	  fprintf (stderr, "[%s] Error reading SAC data samples\n", sacfile);
	  return -1;
	}
    }
  else
    {
      fprintf (stderr, "[%s] Unrecognized format value: %d\n", sacfile, format);
      return -1;
    }      
  
  return sh->npts;
}  /* End of parseSAC() */


/***************************************************************************
 * readBinaryHeaderSAC:
 *
 * Read a binary header from a file and parse into a SAC header
 * struct.  Also determines byte order and sets the swap flag unless
 * already dictated by the format.
 *
 * Returns 0 on sucess or -1 on failure.
 ***************************************************************************/
static int
readBinaryHeaderSAC (FILE *ifp, struct SACHeader *sh, int *format,
		  int *swapflag, int verbose, char *sacfile)
{
  int bigendianhost;
  int32_t hdrver;
  
  /* Read the binary header into memory */
  if ( fread (sh, sizeof(struct SACHeader), 1, ifp) != 1 )
    {
      fprintf (stderr, "[%s] Could not read SAC header from file\n", sacfile);
      
      if ( ferror (ifp) )
	fprintf (stderr, "[%s] Error reading from file\n", sacfile);
      
      return -1;
    }
  
  /* Determine if host is big-endian */
  bigendianhost = ms_bigendianhost();
  
  *swapflag = 0;
  
  /* Test byte order using the header version if unknown */
  /* Also set the swapflag appropriately */
  if ( *format == 2 )
    {
      memcpy (&hdrver, &sh->nvhdr, 4);
      if ( hdrver < 1 || hdrver > 10 )
	{
	  ms_gswap4 (&hdrver);
	  if ( hdrver < 1 || hdrver > 10 )
	    {
	      fprintf (stderr, "[%s] Cannot determine byte order (not SAC?)\n", sacfile);
	      return -1;
	    }
	  
	  *format = ( bigendianhost ) ? 3 : 4;
	  *swapflag = 1;
	}
      else
	{
	  *format =  ( bigendianhost ) ? 4 : 3;
	}
    }
  else if ( *format == 3 && bigendianhost ) *swapflag = 1;
  else if ( *format == 4 && ! bigendianhost ) *swapflag = 1;
  
  if ( verbose > 1 )
    {
      if ( *swapflag )
	fprintf (stderr, "[%s] Byte swapping required\n", sacfile);
      else
	fprintf (stderr, "[%s] Byte swapping NOT required\n", sacfile);
    }
  
  /* Byte swap all values in header */
  if ( *swapflag )
    swapSACHeader (sh);  
  
  return 0;
}  /* End of readBinaryHeaderSAC() */


/***************************************************************************
 * readBinaryDataSAC:
 *
 * Read binary data from a file and add to an array, the array
 * must already be allocated with datacnt floats.
 *
 * Returns 0 on sucess or -1 on failure.
 ***************************************************************************/
static int
readBinaryDataSAC (FILE *ifp, float *data, int datacnt, int swapflag,
		int verbose, char *sacfile)
{
  int samplesread = 0;
  int dataidx;
  
  /* Read in data samples */
  if ( (samplesread = fread (data, sizeof(float), datacnt, ifp)) != datacnt )
    {
      fprintf (stderr, "[%s] Only read %d of %d expected data samples\n",
	       sacfile, samplesread, datacnt);
      return -1;
    }
  
  /* Swap data samples */
  if ( swapflag )
    {
      for ( dataidx = 0; dataidx < datacnt; dataidx++ ) 
	{
	  ms_gswap4 (data + dataidx);
	}
    }
  
  return 0;
}   /* End of readBinaryDataSAC() */


/***************************************************************************
 * readAlphaHeaderSAC:
 *
 * Read a alphanumeric header from a file and parse into a SAC header
 * struct.
 *
 * Returns 0 on sucess or a positive number indicating line number of
 * parsing failure.
 ***************************************************************************/
static int
readAlphaHeaderSAC (FILE *ifp, struct SACHeader *sh)
{
  char line[1025];
  int linecnt = 1;  /* The header starts at line 1 */
  int lineidx;
  int count;
  int hvidx = 0;
  char *cp;
  
  if ( ! ifp || ! sh )
    return -1;
  
  /* The first 14 lines x 5 values are floats */
  for (lineidx=0; lineidx < 14; lineidx++)
    {
      if ( ! fgets(line, sizeof(line), ifp) )
	return linecnt;
      
      count = sscanf (line, " %f %f %f %f %f ", (float *) sh + hvidx,
		      (float *) sh + hvidx + 1, (float *) sh + hvidx + 2,
		      (float *) sh + hvidx + 3, (float *) sh + hvidx + 4);
      
      if ( count != 5 )
	return linecnt;
      
      hvidx += 5;
      linecnt++;
    }
  
  /* The next 8 lines x 5 values are integers */
  for (lineidx=0; lineidx < 8; lineidx++)
    {
      if ( ! fgets(line, sizeof(line), ifp) )
	return linecnt;
      
      count = sscanf (line, " %d %d %d %d %d ", (int32_t *) sh + hvidx,
		      (int32_t *) sh + hvidx + 1, (int32_t *) sh + hvidx + 2,
		      (int32_t *) sh + hvidx + 3, (int32_t *) sh + hvidx + 4);
      
      if ( count != 5 )
	return linecnt;
      
      hvidx += 5;
      linecnt++;
    }
  
  /* Set pointer to start of string variables */
  cp =  (char *) sh + (hvidx * 4);
  
  /* The next 8 lines each contain 24 bytes of string data */
  for (lineidx=0; lineidx < 8; lineidx++)
    {
      memset (line, 0, sizeof(line));
      if ( ! fgets(line, sizeof(line), ifp) )
	return linecnt;
      
      memcpy (cp, line, 24);
      cp += 24;
      
      linecnt++;
    }
  
  /* Make sure each of the 23 string variables are left justified */
  cp =  (char *) sh + (hvidx * 4);  
  for (count=0; count < 24; count++)
    {
      int ridx, widx, width;
      char *fcp;
      
      /* Each string variable is 8 characters with one exception */
      if ( count != 1 )
	{
	  width = 8;
	}
      else
	{
	  width = 16;
	  count++;
	}
      
      /* Pointer to field */
      fcp = cp + (count * 8);

      /* Find first character that is not a space */
      ridx = 0;
      while ( *(fcp + ridx) == ' ' )
	ridx++;
      
      /* Remove any leading spaces */
      if ( ridx > 0 )
	{
	  for (widx=0; widx < width; widx++, ridx++)
	    {
	      if ( ridx < width )
		*(fcp + widx) = *(fcp + ridx);
	      else
		*(fcp + widx) = ' ';
	    }
	}
    }
  
  return 0;
}  /* End of readAlphaHeaderSAC() */


/***************************************************************************
 * readAlphaDataSAC:
 *
 * Read a alphanumeric data from a file and add to an array, the array
 * must already be allocated with datacnt floats.
 *
 * Returns 0 on sucess or a positive number indicating line number of
 * parsing failure.
 ***************************************************************************/
static int
readAlphaDataSAC (FILE *ifp, float *data, int datacnt)
{
  char line[1025];
  int linecnt = 31; /* Data samples start on line 31 */
  int samplesread = 0;
  int count;
  int dataidx = 0;
  
  if ( ! ifp || ! data || ! datacnt)
    return -1;
  
  /* Each data line should contain 5 floats unless the last */
  for (;;)
    {
      if ( ! fgets(line, sizeof(line), ifp) )
	return linecnt;
      
      count = sscanf (line, " %f %f %f %f %f ", (float *) data + dataidx,
		      (float *) data + dataidx + 1, (float *) data + dataidx + 2,
		      (float *) data + dataidx + 3, (float *) data + dataidx + 4);
      
      samplesread += count;
      
      if ( samplesread >= datacnt )
	break;
      else if ( count != 5 )
	return linecnt;
      
      dataidx += 5;
      linecnt++;
    }
  
  return 0;
}  /* End of readAlphaDataSAC() */


/***************************************************************************
 * writeMSEED:
 * 
 * Write data buffer to output file as Mini-SEED, converting sample
 * type as required by the encoding format.
 *
 * Returns the number of samples written, -1 on error.
 ***************************************************************************/
static int
writeMSEED (MSTraceID *id, MSTraceSeg *seg, char *outputfile)
{
  MSTrace *mst = NULL;
  MSRecord *msr = NULL;
  BTime btime;
  struct blkt_100_s Blkt100;
  char outfile[1024];
  FILE *ofp = 0;
  int packedsamples = 0;
  int packedrecords = 0;
  int idx;
  
  int32_t *idata;
  float *fdata;
  double *ddata;
  
  /* Determine file open mode:
   * If bytes have already been written: append
   * If nothing has been written: overwrite */
  char *mode = ( outputbytes > 0 ) ? "ab" : "wb";

  if ( ! id || ! seg || ! outputfile )
    return -1;
  
  if ( seg->numsamples <= 0 || seg->samprate == 0.0 )
    return 0;
  
  /* Populate MSTrace used for packing */
  mst = mst_init (NULL);
  strncpy (mst->network, id->network, sizeof(mst->network)-1);
  mst->network[sizeof(mst->network)-1] = '\0';
  strncpy (mst->station, id->station, sizeof(mst->station)-1);
  mst->station[sizeof(mst->station)-1] = '\0';
  strncpy (mst->location, id->location, sizeof(mst->location)-1);
  mst->location[sizeof(mst->location)-1] = '\0';
  strncpy (mst->channel, (channel)?channel:id->channel, sizeof(mst->channel)-1);
  mst->channel[sizeof(mst->channel)-1] = '\0';
  mst->dataquality = id->dataquality;
  
  mst->starttime = seg->starttime;
  mst->endtime = seg->endtime;
  mst->samprate = seg->samprate;
  mst->samplecnt = seg->samplecnt;
  mst->datasamples = seg->datasamples;
  mst->numsamples = seg->numsamples;
  mst->sampletype = seg->sampletype;
  
  if ( verbose )
    fprintf (stderr, "Writing Mini-SEED for %.8s.%.8s.%.8s.%.8s\n",
	     mst->network, mst->station, mst->location, mst->channel);
  
  if ( ! outputfile )
    {
      if ( ms_hptime2btime (mst->starttime, &btime) )
	{
	  fprintf (stderr, "Cannot convert MSTrace hptime to BTime\n");
	  return -1;
	}
      
      /* Create output file name: Net.Sta.Loc.Chan.Qual.Year.Day.Hour.Min.Sec.mseed */
      snprintf (outfile, sizeof(outfile), "%s%s%s.%s.%s.%s.%.1s.%04d,%03d,%02d:%02d:%02d.mseed",
		(outputdir)?outputdir:"",(outputdir)?"/":"",
		mst->network, mst->station, mst->location, mst->channel,
		&mst->dataquality, btime.year, btime.day, btime.hour,
		btime.min, btime.sec);
    }
  else
    {
      strncpy (outfile, outputfile, sizeof(outfile));
    }
  
  /* Open output file */
  if ( strcmp (outfile, "-") == 0 )
    {
      ofp = stdout;
    }
  else if ( (ofp = fopen (outfile, mode)) == NULL )
    {
      fprintf (stderr, "Cannot open output file: %s (%s)\n",
	       outfile, strerror(errno));
      return -1;
    }
  
  /* Convert data buffer and set the data type depending on encoding format */
  if ( packencoding == 3 || packencoding == 10 || packencoding == 11 )
    {
      /* Use the same buffer (ints being equal to floats and smaller than doubles) */      
      idata = (int32_t *) mst->datasamples;
      
      if ( mst->sampletype == 'f' )
	{
	  if ( verbose )
	    fprintf (stderr, "Converting (truncating) floats to integer for Mini-SEED\n");
	  
	  fdata = (float *) mst->datasamples;
	  for (idx=0; idx < mst->numsamples; idx++)
	    idata[idx] = (int) fdata[idx];
	}
      else if ( mst->sampletype == 'd' )
	{
	  if ( verbose )
	    fprintf (stderr, "Converting (truncating) doubles to integer for Mini-SEED\n");
	  
	  ddata = (double *) mst->datasamples;
	  for (idx=0; idx < mst->numsamples; idx++)
	    idata[idx] = (int) ddata[idx];
	}
      else if ( mst->sampletype != 'i' )
	{
	  fprintf (stderr, "writeMSEED(): Unsupported sample type: %c\n", mst->sampletype);
	  return -1;
	}
      
      mst->sampletype = seg->sampletype = 'i';
    }
  else if ( packencoding == 4 ) /* 32-bit floats */
    {
      /* Use the same buffer (floats being equal to ints and smaller than doubles) */
      fdata = (float *) mst->datasamples;
      
      if ( mst->sampletype == 'i' )
	{
	  if ( verbose )
	    fprintf (stderr, "Converting integer to float for Mini-SEED\n");
	  
	  idata = (int32_t *) mst->datasamples;
	  for (idx=0; idx < mst->numsamples; idx++)
	    fdata[idx] = (float) idata[idx];
	}
      else if ( mst->sampletype == 'd' )
	{
	  if ( verbose )
	    fprintf (stderr, "Converting double to float for Mini-SEED\n");
	  
	  ddata = (double *) mst->datasamples;
	  for (idx=0; idx < mst->numsamples; idx++)
	    fdata[idx] = (float) ddata[idx];
	}
      else if ( mst->sampletype != 'f' )
	{
	  fprintf (stderr, "writeMSEED(): Unsupported sample type: %c\n", mst->sampletype);
	  return -1;
	}
      
      mst->sampletype = seg->sampletype = 'f';
    }
  else if ( packencoding == 5 ) /* 64-bit floats */
    {
      if ( mst->sampletype == 'i' )
	{
	  if ( verbose )
	    fprintf (stderr, "Converting integer to double for Mini-SEED\n");
	  
	  /* Allocate new array of doubles */
	  if ( (ddata = (double *) malloc (sizeof(double) * mst->numsamples)) == NULL )
	    {
	      fprintf (stderr, "Error allocating memory\n");
	      return -1;
	    }
	  
	  idata = (int32_t *) mst->datasamples;
	  for (idx=0; idx < mst->numsamples; idx++)
	    ddata[idx] = (double) idata[idx];
	  
	  /* Free earlier buffer and replace with new array of doubles */
	  if ( mst->datasamples )
	    free (mst->datasamples);
	  mst->datasamples = seg->datasamples = ddata;
	}
      else if ( mst->sampletype == 'f' )
	{
	  if ( verbose )
	    fprintf (stderr, "Converting float to double for Mini-SEED\n");
	  
	  /* Allocate new array of doubles */
	  if ( (ddata = (double *) malloc (sizeof(double) * mst->numsamples)) == NULL )
	    {
	      fprintf (stderr, "Error allocating memory\n");
	      return -1;
	    }
	  
	  fdata = (float *) mst->datasamples;
	  for (idx=0; idx < mst->numsamples; idx++)
	    ddata[idx] = (double) fdata[idx];

	  /* Free earlier buffer and replace with new array of doubles */
	  if ( mst->datasamples )
	    free (mst->datasamples);
	  mst->datasamples = seg->datasamples = ddata;
	}
      else if ( mst->sampletype != 'd' )
	{
	  fprintf (stderr, "writeMSEED(): Unsupported sample type: %c\n", mst->sampletype);
	  return -1;
	}
      
      mst->sampletype = seg->sampletype = 'd';
    }
  else
    {
      fprintf (stderr, "Unrecognized encoding format: %d\n", packencoding);
      return -1;
    }
  
  /* If a blockette 100 is requested add it */
  if ( srateblkt )
    {
      msr = msr_init (NULL);
      memset (&Blkt100, 0, sizeof(struct blkt_100_s));
      Blkt100.samprate = (float) msr->samprate;
      msr_addblockette (msr, (char *) &Blkt100, sizeof(struct blkt_100_s), 100, 0);
    }
  
  /* Pack data into records */
  packedrecords = mst_pack (mst, recordHandler, ofp, packreclen, packencoding,
			    byteorder, &packedsamples, 1, verbose-2, msr);
  
  /* Unless anerror occurred the sample buffer has been released, adjust */
  if ( packedrecords >= 0 )
    {
      seg->datasamples = NULL;
      seg->numsamples = 0;
    }
  
  if ( verbose )
    fprintf (stderr, "Packed %d samples into %d records\n",
	     packedsamples, packedrecords);
  
  /* Make sure everything is cleaned up */
  if ( msr )
    msr_free (&msr);
  
  /* Free temporary MSTrace */
  mst->datasamples = NULL;
  mst->prvtptr = NULL;
  mst_free (&mst);
  
  if ( ofp && ofp != stdout )
    fclose (ofp);
  
  return packedsamples;
}  /* End of writeMSEED) */


/***************************************************************************
 * writeSAC:
 * 
 * Write data buffer to output file as SAC.
 *
 * Returns the number of samples written or -1 on error.
 ***************************************************************************/
static int
writeSAC (MSTraceID *id, MSTraceSeg *seg, int format, char *outputfile)
{
  struct SACHeader sh = NullSACHeader;
  char outfile[1024];
  char sacchannel[11];
  
  float *fdata = 0;
  double *ddata = 0;
  int32_t *idata = 0;
  hptime_t submsec;
  int idx;
  
  if ( ! id || ! seg )
    return -1;
  
  if ( seg->numsamples <= 0 || seg->samprate == 0.0 )
    return 0; 
  
  /* Substitute channel name if specified */
  strncpy (sacchannel, (channel)?channel:id->channel, sizeof(sacchannel)-1);
  sacchannel[sizeof(sacchannel)-1] = '\0';
  
  if ( verbose )
    fprintf (stderr, "Writing SAC for %.8s.%.8s.%.8s.%.8s\n",
	     id->network, id->station, id->location, sacchannel);
  
  /* If an original input SAC header is available use it as a base to update */
  if ( seg->prvtptr )
    {
      struct SACHeader *osh = (struct SACHeader *) seg->prvtptr;
      hptime_t ostarttime;
      
      memcpy (&sh, osh, sizeof (struct SACHeader));
      
      /* Calculate original input start time */
      ostarttime = ms_time2hptime (osh->nzyear, osh->nzjday, osh->nzhour, osh->nzmin, osh->nzsec, osh->nzmsec * 1000);
      
      /* Adjust for Begin ('B' SAC variable) time offset */
      ostarttime += (double) osh->b * HPTMODULUS;
      
      /* If data start time has changed adjust Begin and End header variables */
      if ( ostarttime != seg->starttime )
	{
	  float shift = (float) (seg->starttime - ostarttime) / HPTMODULUS;
	  
	  if ( verbose >= 1 )
	    fprintf (stderr, "Updating data begin and end time variables (%g seconds)\n", shift);
	  
	  sh.b += shift;
	  sh.e = sh.b + (seg->numsamples - 1) * (1 / seg->samprate);
	}
    }
  else /* Otherwise set zero time and time of first and last sample */
    {
      BTime btime;
      
      /* Set time zero to start time */
      ms_hptime2btime (seg->starttime, &btime);
      sh.nzyear = btime.year;
      sh.nzjday = btime.day;
      sh.nzhour = btime.hour;
      sh.nzmin = btime.min;
      sh.nzsec = btime.sec;
      sh.nzmsec = btime.fract / 10;
      
      /* Determine any sub-millisecond portion of the start time in HP time */
      submsec = (seg->starttime - 
		 ms_time2hptime (sh.nzyear, sh.nzjday, sh.nzhour,
				 sh.nzmin, sh.nzsec, sh.nzmsec * 1000));
      
      /* Set begin and end offsets from reference time for first and last sample,
       * any sub-millisecond start time is stored in these offsets. */
      sh.b = ((float)submsec / HPTMODULUS);
      sh.e = sh.b + (seg->numsamples - 1) * (1 / seg->samprate);
    }
  
  /* Set time-series source parameters */
  if ( *id->network != '\0' )
    strncpy (sh.knetwk, id->network, 8);
  if ( *id->station != '\0' )
    strncpy (sh.kstnm, id->station, 8);
  if ( *id->location != '\0' )
    strncpy (sh.khole, id->location, 8);
  if ( *id->channel != '\0' )
    strncpy (sh.kcmpnm, id->channel, 8);
  
  /* Set misc. header variables */
  sh.nvhdr = 6;                 /* Header version = 6 */
  sh.leven = 1;                 /* Evenly spaced data */
  sh.iftype = ITIME;            /* Data is time-series */
  
  /* Insert metadata if present */
  if ( metadata )
    insertSACMetaData (&sh, seg->starttime);
  
  /* Set sampling interval (seconds), sample count */
  sh.delta = 1.0 / seg->samprate;
  sh.npts = seg->numsamples;
  
  /* Original input SAC variables that may no longer be valid:
   * scale : amplitude scale factor
   * idep  : type of amplitude
   */ 
  
  /* Convert data buffer to floats */
  if ( seg->sampletype == 'f' )
    {
      fdata = (float *) seg->datasamples;
      
      /* Determine minimum and maximum sample values */
      sh.depmin = *fdata;
      sh.depmax = *fdata;
      for (idx=1; idx < seg->numsamples; idx++)
	{
	  if ( fdata[idx] < sh.depmin )
	    sh.depmin = fdata[idx];
	  if ( fdata[idx] > sh.depmax )
	    sh.depmax = fdata[idx];
	}
    }
  else if ( seg->sampletype == 'i' )
    {
      idata = (int32_t *) seg->datasamples;
      
      fdata = (float *) malloc (seg->numsamples * sizeof(float));

      if ( fdata == NULL )
	{
	  fprintf (stderr, "Error allocating memory\n");
	  return -1;
	}
      
      /* Populate array of float samples and determine min and max */
      sh.depmin = *idata;
      sh.depmax = *idata;
      for (idx=0; idx < seg->numsamples; idx++)
	{
	  fdata[idx] = (float) idata[idx];
	  
	  if ( fdata[idx] < sh.depmin )
	    sh.depmin = fdata[idx];
	  if ( fdata[idx] > sh.depmax )
	    sh.depmax = fdata[idx];
	}
    }
  else if ( seg->sampletype == 'd' )
    {
      ddata = (double *) seg->datasamples;
      
      fdata = (float *) malloc (seg->numsamples * sizeof(float));
      
      if ( fdata == NULL )
	{
	  fprintf (stderr, "Error allocating memory\n");
	  return -1;
	}
      
      /* Populate array of float samples and determine min and max */
      sh.depmin = *ddata;
      sh.depmax = *ddata;
      for (idx=0; idx < seg->numsamples; idx++)
	{
	  fdata[idx] = (float) ddata[idx];
	  
	  if ( fdata[idx] < sh.depmin )
	    sh.depmin = fdata[idx];
	  if ( fdata[idx] > sh.depmax )
	    sh.depmax = fdata[idx];
	}
    }
  else
    {
      fprintf (stderr, "Error, unrecognized sample type: '%c'\n",
	       seg->sampletype);
      return -1;
    }
  
  if ( format >= 2 && format <= 4 )
    {
      if ( ! outputfile )
	{
	  /* Create output file name: Net.Sta.Loc.Chan.Qual.Year.Day.Hour.Min.Sec.SAC */
	  snprintf (outfile, sizeof(outfile), "%s%s%s.%s.%s.%s.%.1s.%04d,%03d,%02d:%02d:%02d.SAC",
		    (outputdir)?outputdir:"",(outputdir)?"/":"",
		    id->network, id->station, id->location, id->channel,
		    &id->dataquality, sh.nzyear, sh.nzjday, sh.nzhour,
		    sh.nzmin, sh.nzsec);
	}
      else
	{
	  if ( verbose && outputbytes )
	    fprintf (stderr, "Warning, output file will be overwritten: %s\n", outputfile);
	  
	  strncpy (outfile, outputfile, sizeof(outfile));
	}
      
      /* Byte swap the data header and data if needed */
      if ( (format == 3 && ms_bigendianhost()) ||
	   (format == 4 && ! ms_bigendianhost()) )
	{
	  if ( verbose )
	    fprintf (stderr, "Byte swapping SAC header and data\n");

	  swapSACHeader (&sh);
	  
	  for (idx=0; idx < seg->numsamples; idx++)
	    {
	      ms_gswap4 (fdata + idx);
	    }
	}
	   
      if ( verbose > 1 )
	fprintf (stderr, "Writing binary SAC file '%s'\n", outfile);

      /* Reset output bytes counter */
      outputbytes = 0;
      
      if ( writeBinarySAC (&sh, fdata, seg->numsamples, outfile) )
	return -1;
    }
  else if ( format == 1 )
    {
      if ( ! outputfile )
	{
	  /* Create output file name: Net.Sta.Loc.Chan.Qual.Year.Day.Hour.Min.Sec.SACA */
	  snprintf (outfile, sizeof(outfile), "%s%s%s.%s.%s.%s.%.1s.%04d,%03d,%02d:%02d:%02d.SACA",
		    (outputdir)?outputdir:"",(outputdir)?"/":"",
		    id->network, id->station, id->location, id->channel,
		    &id->dataquality, sh.nzyear, sh.nzjday, sh.nzhour,
		    sh.nzmin, sh.nzsec);
	}
      else
	{
	  if ( verbose && outputbytes )
	    fprintf (stderr, "Warning, output file will be overwritten: %s\n", outputfile);
	  
	  strncpy (outfile, outputfile, sizeof(outfile));
	}
      
      if ( verbose > 1 )
	fprintf (stderr, "Writing alphanumeric SAC file: '%s'\n", outfile);
      
      /* Reset output bytes counter */
      outputbytes = 0;
      
      if ( writeAlphaSAC (&sh, fdata, seg->numsamples, outfile) )
	return -1;
    }
  else
    {
      fprintf (stderr, "Error, unrecognized SAC format: '%d'\n", format);
    }
  
  /* Free any buffer allocated by this routine */
  if ( fdata && seg->sampletype != 'f' )
    free (fdata);
  
  if ( verbose )
    fprintf (stderr, "Wrote %d samples to %s (%d bytes)\n",
	     seg->numsamples, outfile, outputbytes);
  
  return seg->numsamples;
}  /* End of writeSAC() */


/***************************************************************************
 * writeBinarySAC:
 * Write binary SAC file.
 *
 * Returns 0 on success, and -1 on failure.
 ***************************************************************************/
static int
writeBinarySAC (struct SACHeader *sh, float *fdata, int npts, char *outfile)
{
  FILE *ofp;
  
  /* Open output file */
  if ( (ofp = fopen (outfile, "wb")) == NULL )
    {
      fprintf (stderr, "Cannot open output file: %s (%s)\n",
	       outfile, strerror(errno));
      return -1;
    }
  
  /* Write SAC header to output file */
  if ( fwrite (sh, sizeof(struct SACHeader), 1, ofp) != 1 )
    {
      fprintf (stderr, "Error writing SAC header to output file\n");
      return -1;
    }
  outputbytes += sizeof(struct SACHeader);
  
  /* Write float data to output file */
  if ( fwrite (fdata, sizeof(float), npts, ofp) != npts )
    {
      fprintf (stderr, "Error writing SAC data to output file\n");
      return -1;
    }
  outputbytes += sizeof(float) * npts;
  
  fclose (ofp);
  
  return 0;
}  /* End of writeBinarySAC() */


/***************************************************************************
 * writeAlphaSAC:
 * Write alphanumeric SAC file.
 *
 * Returns 0 on success, and -1 on failure.
 ***************************************************************************/
static int
writeAlphaSAC (struct SACHeader *sh, float *fdata, int npts, char *outfile)
{
  FILE *ofp;
  int idx, fidx;
  
  /* Declare and set up pointers to header variable type sections */
  float   *fhp = (float *) sh;
  int32_t *ihp = (int32_t *) sh + (NUMFLOATHDR);
  char    *shp = (char *) sh + (NUMFLOATHDR * 4 + NUMINTHDR * 4);
  
  /* Open output file */
  if ( (ofp = fopen (outfile, "wb")) == NULL )
    {
      fprintf (stderr, "Cannot open output file: %s (%s)\n",
	       outfile, strerror(errno));
      return -1;
    }
  
  /* Write SAC header float variables to output file, 5 variables per line */
  for (idx=0; idx < NUMFLOATHDR; idx += 5)
    {
      for (fidx=idx; fidx < (idx+5) && fidx < NUMFLOATHDR; fidx++)
	fprintf (ofp, "%#15.7g", *(fhp + fidx));
      
      fprintf (ofp, "\n");
    }
  outputbytes += 1064;
  
  /* Write SAC header integer variables to output file, 5 variables per line */
  for (idx=0; idx < NUMINTHDR; idx += 5)
    {
      for (fidx=idx; fidx < (idx+5) && fidx < NUMINTHDR; fidx++)
	fprintf (ofp, "%10d", *(ihp + fidx));
      
      fprintf (ofp, "\n");
    }
  outputbytes += 408;
  
  /* Write SAC header string variables to output file, 3 variables per line */
  for (idx=0; idx < NUMSTRHDR; idx += 3)
    {
      if ( idx == 0 )
	fprintf (ofp, "%-8.8s%-16.16s", shp, shp + 8);
      else
	fprintf (ofp, "%-8.8s%-8.8s%-8.8s", shp+(idx*8), shp+((idx+1)*8), shp+((idx+2)*8));
      
      fprintf (ofp, "\n");
    }
  outputbytes += 200;
  
  /* Write float data to output file, 5 values per line */
  for (idx=0; idx < npts; idx += 5)
    {
      for (fidx=idx; fidx < (idx+5) && fidx < npts && fidx >= 0; fidx++)
	{
	  fprintf (ofp, "%#15.7g", *(fdata + fidx));
	  outputbytes += 15;
	}
      
      fprintf (ofp, "\n");
      outputbytes++;
    }
  
  fclose (ofp);
  
  return 0;
}  /* End of writeAlphaSAC() */


/***************************************************************************
 * insertSACMetaData:
 *
 * Search the metadata list for the first matching source and insert
 * the metadata into the SAC header if found.  The source names (net,
 * sta, loc, chan) are used to find a match.  If metadata list entries
 * include a '*' they will match everything, for example if the
 * channel field is '*' all channels for the specified network,
 * station and location will match the list entry.
 *
 * The metadata list should be populated with an array of pointers to:
 *  0:  Network (knetwk)
 *  1:  Station (kstnm)
 *  2:  Location (khole)
 *  3:  Channel (kcmpnm)
 *  4:  Latitude (stla)
 *  5:  Longitude (stlo)
 *  6:  Elevation (stel) [not currently used by SAC]
 *  7:  Depth (stdp) [not currently used by SAC]
 *  8:  Component Azimuth (cmpaz), degrees clockwise from north
 *  9:  Component Incident Angle (cmpinc), degrees from vertical
 *  10: Instrument Name (kinst)
 *  11: Scale Factor (scale)
 *  12: Scale Frequency, unused
 *  13: Scale Units, unused
 *  14: Sampling rate, unused
 *  15: Start time, used for matching
 *  16: End time, used for matching
 *
 * Returns 0 on sucess, 1 when no matching metadata found and -1 on failure.
 ***************************************************************************/
static int
insertSACMetaData (struct SACHeader *sh, hptime_t sacstarttime)
{
  struct metalist *mlp = metadata;
  struct metanode *mn = NULL;
  hptime_t sacendtime;
  char *endptr;
  char sacnetwork[9];
  char sacstation[9];
  char saclocation[9];
  char sacchannel[9];
  int retval = 1;
  
  if ( ! mlp || ! sh )
    return -1;
  
  /* Determine source name parameters for comparison, as a special case if the 
   * location code is not set it will match '--' */
  if ( strncmp (sh->knetwk, SUNDEF, 8) ) ms_strncpclean (sacnetwork, sh->knetwk, 8);
  else sacnetwork[0] = '\0';
  if ( strncmp (sh->kstnm, SUNDEF, 8) ) ms_strncpclean (sacstation, sh->kstnm, 8);
  else sacstation[0] = '\0';
  if ( strncmp (sh->khole, SUNDEF, 8) ) ms_strncpclean (saclocation, sh->khole, 8);
  else { saclocation[0] = '-'; saclocation[1] = '-'; saclocation[2] = '\0'; }
  if ( strncmp (sh->kcmpnm, SUNDEF, 8) ) ms_strncpclean (sacchannel, sh->kcmpnm, 8);
  else sacchannel[0] = '\0';
  
  /* Calculate end time of SAC data */
  sacendtime = sacstarttime + (((sh->npts - 1) * sh->delta) * HPTMODULUS);
  
  while ( mlp )
    {
      mn = (struct metanode *) mlp->data;
      
      /* Sanity check that source name fields are present */
      if ( ! mn->metafields[0] || ! mn->metafields[1] || 
	   ! mn->metafields[2] || ! mn->metafields[3] )
	{
	  fprintf (stderr, "insertmetadata(): error, source name fields not all present\n");
	}
      /* Test if network, station, location and channel; also handle simple wildcards */
      else if ( ( ! strncmp (sacnetwork, mn->metafields[0], 8) || (*(mn->metafields[0]) == '*') ) &&
		( ! strncmp (sacstation, mn->metafields[1], 8) || (*(mn->metafields[1]) == '*') ) &&
		( ! strncmp (saclocation, mn->metafields[2], 8) || (*(mn->metafields[2]) == '*') ) &&
		( ! strncmp (sacchannel, mn->metafields[3], 8) || (*(mn->metafields[3]) == '*') ) )
	{
	  /* Check time window match */
	  if ( mn->starttime != HPTERROR || mn->endtime != HPTERROR )
	    {
	      /* Check for overlap with metadata window */
	      if ( mn->starttime != HPTERROR && mn->endtime != HPTERROR )
		{
		  if ( ! (sacendtime >= mn->starttime && sacstarttime <= mn->endtime) )
		    {
		      mlp = mlp->next;
		      continue;
		    }
		}
	      /* Check if data after start time */
	      else if ( mn->starttime != HPTERROR )
		{
		  if ( sacendtime < mn->starttime )
		    {
		      mlp = mlp->next;
		      continue;
		    }
		}
	      /* Check if data before end time */
	      else if ( mn->endtime != HPTERROR )
		{
		  if ( sacstarttime > mn->endtime )
		    {
		      mlp = mlp->next;
		      continue;
		    }
		}
	    }
	  
	  if ( verbose )
	    fprintf (stderr, "Inserting metadata for N: '%s', S: '%s', L: '%s', C: '%s' (%s - %s)\n",
		     sacnetwork, sacstation, saclocation, sacchannel,
		     (mn->metafields[15])?mn->metafields[15]:"NONE",
		     (mn->metafields[16])?mn->metafields[16]:"NONE");
	  
	  /* Insert metadata into SAC header */
	  if ( mn->metafields[4] ) sh->stla = (float) strtod (mn->metafields[4], &endptr);
	  if ( mn->metafields[5] ) sh->stlo = (float) strtod (mn->metafields[5], &endptr);
	  if ( mn->metafields[6] ) sh->stel = (float) strtod (mn->metafields[6], &endptr);
	  if ( mn->metafields[7] ) sh->stdp = (float) strtod (mn->metafields[7], &endptr);
	  if ( mn->metafields[8] ) sh->cmpaz = (float) strtod (mn->metafields[8], &endptr);
	  if ( mn->metafields[9] ) {
	    sh->cmpinc = (float) strtod (mn->metafields[9], &endptr);
	    if ( seedinc ) sh->cmpinc += 90;
	  }
	  if ( mn->metafields[10] ) strncpy (sh->kinst, mn->metafields[10], 8);
	  if ( mn->metafields[11] ) sh->scale = (float) strtod (mn->metafields[11], &endptr);
	  
	  retval = 0;
	  break;
	}
      
      mlp = mlp->next;
    }
  
  return retval;
}  /* End of insertSACMetaData() */


/***************************************************************************
 * swapSACHeader:
 *
 * Byte swap all multi-byte quantities (floats and ints) in SAC header
 * struct.
 *
 * Returns 0 on sucess and -1 on failure.
 ***************************************************************************/
static int
swapSACHeader (struct SACHeader *sh)
{
  int32_t *ip;
  int idx;
  
  if ( ! sh )
    return -1;
  
  for ( idx=0; idx < (NUMFLOATHDR + NUMINTHDR); idx++ )
    {
      ip = (int32_t *) sh + idx;
      ms_gswap4 (ip);
    }
  
  return 0;
}  /* End of swapSACHeader() */


/***************************************************************************
 * differentiate2:
 * 
 * Perform uncentered, two-point differentiation.  The output array
 * can be the same as the input array, otherwise the output array must
 * already be allocated.  The length of the output array will be one
 * less than the input array.
 *
 * This differentiation implies a shift of the time by 1/2 step.
 *
 * The input data can be either floats or doubles, the inputtype must
 * be set to either 'f' or 'd' indicating how the input data should be
 * treated.
 *
 * Returns the number of samples in the output array on success and -1
 * on failure.
 ***************************************************************************/
static int
differentiate2 (void *input, char inputtype, int length, double rate, void *output)
{
  int idx;
  
  float *fin, *fout;
  double *din, *dout;
  
  if ( ! input || ! output )
    {
      fprintf (stderr, "differentiate2(): input or output pointers are NULL\n");      
      return -1;
    }
  
  if ( inputtype != 'f' && inputtype != 'd' )
    {
      fprintf (stderr, "differentiate2(): unrecognized input sample type: %c\n",
	       inputtype);
      return -1;
    }
  
  if ( inputtype == 'f' )
    {
      fin = input;
      fout = output;
      
      for (idx=0; idx < (length-1); idx++)
	{
	  fout[idx] = rate * (fin[idx+1] - fin[idx]);
	}
    }
  else if ( inputtype == 'd' )
    {
      din = input;
      dout = output;
      
      for (idx=0; idx < (length-1); idx++)
	{
	  dout[idx] = rate * (din[idx+1] - din[idx]);
	}
    }

  return length - 1;
} /* End of differentiate2() */


/***************************************************************************
 * integrateTrap:
 * 
 * Perform integration using the trapezoidal (midpoint) method.  The
 * output array can be the same as the input array, otherwise the
 * output array must already be allocated.  The length of the output
 * array will be one less than the input array.
 *
 * This integration implies a shift of the time by 1/2 step.
 *
 * The input data can be either floats or doubles, the inputtype must
 * be set to either 'f' or 'd' indicating how the input data should be
 * treated.
 *
 * Returns the number of samples in the output array on success and -1
 * on failure.
 ***************************************************************************/
static int
integrateTrap (void *input, char inputtype, int length, double halfstep, void *output)
{
  int idx;
  double prtint;
  double totint = 0.0;
  
  float *fin, *fout;
  double *din, *dout;
  
  if ( ! input || ! output )
    {
      fprintf (stderr, "integrateTrap(): input or output pointers are NULL\n");      
      return -1;
    }
  
  if ( inputtype != 'f' && inputtype != 'd' )
    {
      fprintf (stderr, "integrateTrap(): unrecognized input sample type: %c\n",
	       inputtype);
      return -1;
    }
  
  if ( inputtype == 'f' )
    {
      fin = input;
      fout = output;
      
      for (idx=0; idx < (length-1); idx++)
	{
	  prtint = halfstep * (fin[idx] + fin[idx+1]);
	  totint = totint + prtint;
	  fout[idx] = totint;
	}
    }
  else if ( inputtype == 'd' )
    {
      din = input;
      dout = output;
      
      dout[0] = din[0];
      
      for (idx=0; idx < (length-1); idx++)
	{
	  prtint = halfstep * (din[idx] + din[idx+1]);
	  totint = totint + prtint;
	  dout[idx] = totint;
	}
    }
  
  return length - 1;
} /* End of integrateTrap() */


/***************************************************************************
 * parameterProc():
 * Process the command line parameters.
 *
 * Returns 0 on success, and -1 on failure
 ***************************************************************************/
static int
parameterProc (int argcount, char **argvec)
{
  char  *filename;
  char  *filtstr = NULL;
  char  *spectaperstr = NULL;
  char  *dBdownstr = NULL;
  char  *tptr;
  int    optind;
  
  /* Process all command line arguments */
  for (optind = 1; optind < argcount; optind++)
    {
      if (strcmp (argvec[optind], "-V") == 0)
	{
	  fprintf (stderr, "%s version: %s\n", PACKAGE, VERSION);
	  exit (0);
	}
      else if (strcmp (argvec[optind], "-h") == 0)
	{
	  usage();
	  exit (0);
	}
      else if (strncmp (argvec[optind], "-v", 2) == 0)
	{
	  verbose += strspn (&argvec[optind][1], "v");
	}
      else if (strcmp (argvec[optind], "-s") == 0)
	{
	  basicsum = -1;
	}
      else if (strcmp (argvec[optind], "-tt") == 0)
	{
	  timetol = strtod (getOptVal(argcount, argvec, optind++, 0), NULL);
	}
      else if (strcmp (argvec[optind], "-rt") == 0)
	{
	  sampratetol = strtod (getOptVal(argcount, argvec, optind++, 0), NULL);
	}
      else if (strcmp (argvec[optind], "-ts") == 0)
	{
	  starttime = ms_seedtimestr2hptime (getOptVal(argcount, argvec, optind++, 0));
	  if ( starttime == HPTERROR )
	    return -1;
	}
      else if (strcmp (argvec[optind], "-te") == 0)
	{
	  endtime = ms_seedtimestr2hptime (getOptVal(argcount, argvec, optind++, 0));
	  if ( endtime == HPTERROR )
	    return -1;
	}
      else if (strcmp (argvec[optind], "-MSEED") == 0)
        {
          dataformat = 1;
        }
      else if (strcmp (argvec[optind], "-Mr") == 0)
	{
	  reclen = strtol (getOptVal(argcount, argvec, optind++, 1), NULL, 10);
	}
      else if (strcmp (argvec[optind], "-Me") == 0)
	{
	  encodingstr = getOptVal(argcount, argvec, optind++, 0);
	}
      else if (strcmp (argvec[optind], "-MR") == 0)
	{
	  packreclen = strtol (getOptVal(argcount, argvec, optind++, 0), NULL, 10);
	}
      /*
      else if (strcmp (argvec[optind], "-S") == 0)
        {
          srateblkt = 1;
        }
      */
      else if (strcmp (argvec[optind], "-ME") == 0)
	{
	  packencoding = strtol (getOptVal(argcount, argvec, optind++, 0), NULL, 10);
	}
      else if (strcmp (argvec[optind], "-Sif") == 0)
	{
	  sacinformat = strtoul (getOptVal(argcount, argvec, optind++, 0), NULL, 10);
	}
      else if (strcmp (argvec[optind], "-Sf") == 0)
	{
	  sacoutformat = strtoul (getOptVal(argcount, argvec, optind++, 0), NULL, 10);
	}
      else if (strcmp (argvec[optind], "-m") == 0)
	{
	  metadatafile = getOptVal(argcount, argvec, optind++, 0);
	}
      else if (strcmp (argvec[optind], "-msi") == 0)
        {
          seedinc = 1;
        }
      else if (strcmp (argvec[optind], "-o") == 0)
        {
          outputfile = getOptVal(argcount, argvec, optind++, 1);
        }
      else if (strcmp (argvec[optind], "-od") == 0)
        {
          outputdir = getOptVal(argcount, argvec, optind++, 0);
        }
      else if (strcmp (argvec[optind], "-c") == 0)
        {
          channel = getOptVal(argcount, argvec, optind++, 0);
        }
      else if (strcmp (argvec[optind], "-CR") == 0)
        {
          filename = getOptVal(argcount, argvec, optind++, 0);
	  if ( filename )
	    addProcess (PROC_CONVOLVE, filename, NULL, PROC_CONVRESP, 0, 0.0, 0.0);
        }
      else if (strcmp (argvec[optind], "-DR") == 0)
        {
          filename = getOptVal(argcount, argvec, optind++, 0);
	  if ( filename )
	    addProcess (PROC_CONVOLVE, filename, NULL, PROC_DECONVRESP, 0, 0.0, 0.0);
        }
      else if (strcmp (argvec[optind], "-CS") == 0)
        {
          filename = getOptVal(argcount, argvec, optind++, 0);
	  if ( filename )
	    addProcess (PROC_CONVOLVE, filename, NULL, PROC_CONVSAC, 0, 0.0, 0.0);
        }
      else if (strcmp (argvec[optind], "-DS") == 0)
        {
          filename = getOptVal(argcount, argvec, optind++, 0);
	  if ( filename )
	    addProcess (PROC_CONVOLVE, filename, NULL, PROC_DECONVSAC, 0, 0.0, 0.0);
        }
      else if (strcmp (argvec[optind], "-ST") == 0)
	{
	  spectaperstr = getOptVal(argcount, argvec, optind++, 1);
	}
      else if (strcmp (argvec[optind], "-STa") == 0)
	{
	  dBdownstr = getOptVal(argcount, argvec, optind++, 1);
	}
      else if (strcmp (argvec[optind], "-W") == 0)
	{
	  if ( prewhiten < 0 )
	    prewhiten *= -1;
	  else
	    prewhiten = 6;
	}
      else if (strcmp (argvec[optind], "-Wo") == 0)
	{
	  if ( prewhiten )
	    prewhiten = strtol (getOptVal(argcount, argvec, optind++, 0), NULL, 10);
	  else
	    prewhiten = -1 * strtol (getOptVal(argcount, argvec, optind++, 0), NULL, 10);
	}
      else if (strcmp (argvec[optind], "-Rts") == 0)
	{
	  resptotalsens = 1;
	}
      else if (strcmp (argvec[optind], "-Red") == 0)
	{
	  respusedelay = 1;
	}
      else if (strcmp (argvec[optind], "-Rin") == 0)
	{
	  respusename = 0;
	}
      else if (strcmp (argvec[optind], "-Ru") == 0)
	{
	  respunits = getOptVal(argcount, argvec, optind++, 0);
	}
      else if (strcmp (argvec[optind], "-LP") == 0)
	{
	  filtstr = getOptVal(argcount, argvec, optind++, 0);
	  addProcess (PROC_LPFILTER, filtstr, NULL, 1, 0, 0.0, 0.0);
	}
      else if (strcmp (argvec[optind], "-LP1") == 0)
	{
	  filtstr = getOptVal(argcount, argvec, optind++, 0);
	  addProcess (PROC_LPFILTER, filtstr, NULL, 0, 0, 0.0, 0.0);
	}
      else if (strcmp (argvec[optind], "-HP") == 0)
	{
	  filtstr = getOptVal(argcount, argvec, optind++, 0);
	  addProcess (PROC_HPFILTER, filtstr, NULL, 1, 0, 0.0, 0.0);
	}
      else if (strcmp (argvec[optind], "-HP1") == 0)
	{
	  filtstr = getOptVal(argcount, argvec, optind++, 0);
	  addProcess (PROC_HPFILTER, filtstr, NULL, 0, 0, 0.0, 0.0);
	}
      else if (strcmp (argvec[optind], "-BP") == 0)
	{
	  char *hpfiltstr = getOptVal(argcount, argvec, optind++, 0);
	  char *lpfiltstr = strchr (hpfiltstr, ':');
	  if ( lpfiltstr )
	    {
	      *lpfiltstr++ ='\0';
	    }
	  else
	    {
	      fprintf (stderr, "Band-pass filter option requires a ':'\n");
	      exit (1);
	    }
	  addProcess (PROC_BPFILTER, lpfiltstr, hpfiltstr, 1, 0, 0.0, 0.0);
	}
      else if (strcmp (argvec[optind], "-BP1") == 0)
	{
	  char *hpfiltstr = getOptVal(argcount, argvec, optind++, 0);
	  char *lpfiltstr = strchr (hpfiltstr, ':');
	  if ( lpfiltstr )
	    {
	      *lpfiltstr++ ='\0';
	    }
	  else
	    {
	      fprintf (stderr, "Band-pass filter option requires a ':'\n");
	      exit (1);
	    }
	  addProcess (PROC_BPFILTER, lpfiltstr, hpfiltstr, 0, 0, 0.0, 0.0);
	}
      else if (strcmp (argvec[optind], "-D2") == 0)
        {
	  addProcess (PROC_DIFF2, NULL, NULL, 0, 0, 0.0, 0.0);
        }
      else if (strcmp (argvec[optind], "-IT") == 0)
        {
	  addProcess (PROC_INTTRAP, NULL, NULL, 0, 0, 0.0, 0.0);
        }
      else if (strcmp (argvec[optind], "-RM") == 0)
        {
	  addProcess (PROC_RMEAN, NULL, NULL, 0, 0, 0.0, 0.0);
        }
      else if (strcmp (argvec[optind], "-SC") == 0)
        {
	  double scalefactor = strtod (getOptVal(argcount, argvec, optind++, 1), NULL);
	  addProcess (PROC_SCALE, NULL, NULL, 0, 0, scalefactor, 0.0);
        }
      else if (strcmp (argvec[optind], "-SI") == 0)
        {
	  double scalefactor = strtod (getOptVal(argcount, argvec, optind++, 1), NULL);
	  addProcess (PROC_SCALE, NULL, NULL, 0, 0, (scalefactor)?1.0/scalefactor:0.0, 0.0);
        }
      else if (strcmp (argvec[optind], "-DEC") == 0)
        {
	  int decimfactor = strtol (getOptVal(argcount, argvec, optind++, 0), NULL, 10);
	  addProcess (PROC_DECIMATE, NULL, NULL, decimfactor, 0, 0.0, 0.0);
        }
      else if (strcmp (argvec[optind], "-TAP") == 0)
        {
	  double taperwidth = 0.0;
	  int tapertype = TAPER_HANNING;
	  
	  char *taperwidthstr = getOptVal(argcount, argvec, optind++, 0);
	  char *tapertypestr = strchr (taperwidthstr, ':');

	  if ( tapertypestr )
	    *tapertypestr++ ='\0';
	  
	  taperwidth = strtod (taperwidthstr, NULL);
	  
	  if ( tapertypestr && *tapertypestr )
	    {
	      if ( ! strcasecmp (tapertypestr, "HANNING") )
		tapertype = TAPER_HANNING;
	      else if ( ! strcasecmp (tapertypestr, "HAMMING") )
		tapertype = TAPER_HAMMING;
	      else if ( ! strcasecmp (tapertypestr, "COSINE") )
		tapertype = TAPER_COSINE;
	      else
		{
		  fprintf (stderr, "Unrecognized taper type: '%s'\n", tapertypestr);
		  exit (1);
		}
	    }
	  
	  addProcess (PROC_TAPER, NULL, NULL, tapertype, 0, taperwidth, 0.0);
        }
      else if (strncmp (argvec[optind], "-", 1) == 0 &&
	       strlen (argvec[optind]) > 1 )
	{
	  fprintf (stderr, "Unknown option: %s\n", argvec[optind]);
	  exit (1);
	}
      else
	{
	  tptr = argvec[optind];
          
          /* Check for an input file list */
          if ( tptr[0] == '@' )
            {
              if ( addListFile (tptr+1) < 0 )
                {
                  fprintf (stderr, "Error adding list file %s", tptr+1);
                  exit (1);
                }
            }
          /* Otherwise this is an input file */
          else
            {
              /* Add file to global file list */
              if ( addFile (tptr) )
                {
                  fprintf (stderr, "Error adding file to input list %s", tptr);
                  exit (1);
                }
            }
	}
    }
  
  /* Make sure input files were specified */
  if ( filelist == 0 )
    {
      fprintf (stderr, "No input files were specified\n\n");
      fprintf (stderr, "%s version %s\n\n", PACKAGE, VERSION);
      fprintf (stderr, "Try %s -h for usage\n", PACKAGE);
      exit (1);
    }
  
  /* Make sure output directory exists */
  if ( outputdir )
    {
      if ( access (outputdir, W_OK) )
	{
	  fprintf (stderr, "Cannot write to output diretory: %s (%s)\n",
		   outputdir, strerror(errno));
	  exit (1);
	}
    }
  
  /* Set prewhitening if only -Wo was specified without -W */
  if ( prewhiten < 0 )
    prewhiten *= -1;
  
  /* Parse taper envelope frequencies */
  if ( spectaperstr )
    {
      int parsed = 0;
      
      if ( ! spectaperfreq )
	{
	  if ( (spectaperfreq = (double *) malloc (4 * sizeof(double))) == NULL )
	    {
	      fprintf (stderr, "Error allocating memory\n");
	      exit(1);
	    }
	  
	  spectaperfreq[0] = -1.0;
	  spectaperfreq[1] = -1.0;
	  spectaperfreq[2] = -1.0;
	  spectaperfreq[3] = -1.0;
	}
      
      parsed = sscanf (spectaperstr, "%lf/%lf/%lf/%lf",
		       &spectaperfreq[0], &spectaperfreq[1],
		       &spectaperfreq[2], &spectaperfreq[3]);
      
      if ( parsed != 4 )
	{
	  fprintf (stderr, "Taper frequcies specified incorrectly: %s\n", spectaperstr);
	  fprintf (stderr, "Try %s -h for usage\n", PACKAGE);
	  exit(1);
	}
      
      if ( spectaperfreq[0] > spectaperfreq[1] )
	{
	  fprintf (stderr, "Taper frequcies specified incorrectly: %s\n", spectaperstr);
	  fprintf (stderr, "Cut-off frequency of lower taper bound (%g) cannot be greater than pass frequency (%g)\n",
		   spectaperfreq[0], spectaperfreq[1]);
	  exit(1);
	}
      
      if ( spectaperfreq[2] > spectaperfreq[3] )
	{
	  fprintf (stderr, "Taper frequcies specified incorrectly: %s\n", spectaperstr);
	  fprintf (stderr, "Cut-off frequency of upper taper bound (%g) cannot be less than pass frequency (%g)\n",
		   spectaperfreq[3], spectaperfreq[2]);
	  exit(1);
	}
    }
  
  /* Parse dB down auto taper cutoff */
  if ( dBdownstr )
    {
      if ( ! spectaperfreq )
	{
	  if ( (spectaperfreq = (double *) malloc (4 * sizeof(double))) == NULL )
	    {
	      fprintf (stderr, "Error allocating memory\n");
	      exit(1);
	    }
	  
	  spectaperfreq[0] = -1.0;
	  spectaperfreq[1] = -1.0;
	  spectaperfreq[2] = -1.0;
	  spectaperfreq[3] = -1.0;
	}
      
      /* dB down cutoff string is lowercorner[/uppercorner] */
      
      /* Check for lower & upper corner cutoff separator */
      tptr = strchr (dBdownstr, '/');
      
      /* Parse lowercorner/uppercorner cutoff pair */
      if ( tptr )
	{
	  *tptr++ = '\0';
	  lcdBdown = strtod (dBdownstr, NULL);
	  ucdBdown = strtod (tptr, NULL);
	}
      /* Parse lowercorner cutoff */
      else
	{
	  lcdBdown = strtod (dBdownstr, NULL);
	}
    }
  
  /* Sanity check RESP units */
  if ( respunits )
    {
      if ( strcmp (respunits, "DIS") && strcmp (respunits, "VEL") &&
	   strcmp (respunits, "ACC") && strcmp (respunits, "DEF") )
	{
	  fprintf (stderr, "Unrecognized SEED RESP response units: %s\n", respunits);
	  fprintf (stderr, "Should be one of DIS, VEL, ACC or DEF\n");
	  exit(1);
	}
    }
  else
    {
      /* Default is no unit conversion */
      respunits = "DEF";
    }
  
  /* Report the program version */
  if ( verbose )
    fprintf (stderr, "%s version: %s\n", PACKAGE, VERSION);
  
  /* Read metadata file if specified */
  if ( metadatafile )
    {
      if ( readMetaData (metadatafile) )
	{
	  fprintf (stderr, "Error reading metadata file\n");
	  return -1;
	}
    }
  
  return 0;
}  /* End of parameterProc() */


/***************************************************************************
 * getOptVal:
 * Return the value to a command line option; checking that the value is 
 * itself not an option (starting with '-') and is not past the end of
 * the argument list.
 *
 * argcount: total arguments in argvec
 * argvec: argument list
 * argopt: index of option to process, value is expected to be at argopt+1
 * dasharg: if the argument can potentially start with a dash
 *
 * Returns value on success and exits with error message on failure
 ***************************************************************************/
static char *
getOptVal (int argcount, char **argvec, int argopt, int dasharg)
{
  if ( argvec == NULL || argvec[argopt] == NULL ) {
    fprintf (stderr, "getOptVal(): NULL option requested\n");
    exit (1);
    return 0;
  }
  
  /* When the value potentially starts with a dash (-) */
  if ( (argopt+1) < argcount && dasharg )
    return argvec[argopt+1];
  
  /* Otherwise check that the value is not another option */
  if ( (argopt+1) < argcount && *argvec[argopt+1] != '-' )
    return argvec[argopt+1];
  
  fprintf (stderr, "Option %s requires a value\n", argvec[argopt]);
  exit (1);
  return 0;
}  /* End of getOptVal() */


/***************************************************************************
 * readMetaData:
 *
 * Read a file of metadata into a structured list, each line should
 * contain the following fields comma-separated in this order:
 *
 * The metadata list should be populated with an array of pointers to:
 *  0:  Network (knetwk)
 *  1:  Station (kstnm)
 *  2:  Location (khole)
 *  3:  Channel (kcmpnm)
 *  4:  Latitude (stla)
 *  5:  Longitude (stlo)
 *  6:  Elevation (stel) [not currently used by SAC]
 *  7:  Depth (stdp) [not currently used by SAC]
 *  8:  Component Azimuth (cmpaz), degrees clockwise from north
 *  9:  Component Incident Angle (cmpinc), degrees from vertical
 *  10: Instrument Name (kinst)
 *  11: Scale Factor (scale)
 *  12: Scale Frequency, unused
 *  13: Scale Units, unused
 *  14: Sampling rate, unused
 *  15: Start time, used for matching
 *  16: End time, used for matching
 *
 * Any lines not containing at least 3 commas are skipped.  If fields
 * are not specified the values are set to NULL with the execption
 * that the first 4 fields (net, sta, loc & chan) cannot be empty.
 *
 * Any lines beginning with '#' are skipped, think comments.
 *
 * Returns 0 on sucess and -1 on failure.
 ***************************************************************************/
static int
readMetaData (char *metafile)
{
  struct metanode mn;
  FILE *mfp;
  char line[1024];
  char *lineptr;
  char *fp;
  int idx, count;
  int linecount = 0;
  
  if ( ! metafile )
    return -1;
  
  if ( (mfp = fopen (metafile, "rb")) == NULL )
    {
      fprintf (stderr, "Cannot open metadata output file: %s (%s)\n",
	       metafile, strerror(errno));
      return -1;
    }

  if ( verbose )
    fprintf (stderr, "Reading station/channel metadata from %s\n", metafile);
  
  while ( fgets (line, sizeof(line), mfp) )
    {
      linecount++;

      /* Truncate at line return if any */
      if ( (fp = strchr (line, '\n')) )
	*fp = '\0';
      
      /* Count the number of commas */
      count = 0;
      fp = line;
      while ( (fp = strchr (fp, ',')) )
	{
	  count++;
	  fp++;
	}
      
      /* Must have at least 3 commas for Net, Sta, Loc, Chan ... */
      if ( count < 3 )
	{
	  if ( verbose > 1 )
	    fprintf (stderr, "Skipping metadata line: %s\n", line);
	  continue;
	}
      
      /* Check for comment line beginning with '#' */
      if ( line[0] == '#' )
	{
	  if ( verbose > 1 )
	    fprintf (stderr, "Skipping comment line: %s\n", line);
	  continue;
	}
      
      /* Create a copy of the line */
      lineptr = strdup (line);
      
      mn.metafields[0] = fp = lineptr;
      mn.starttime = HPTERROR;
      mn.endtime = HPTERROR;
      
      /* Separate line on commas and index in metafields array */
      for (idx = 1; idx < MAXMETAFIELDS; idx++)
	{
	  mn.metafields[idx] = NULL;
	  
	  if ( fp )
	    {
	      if ( (fp = strchr (fp, ',')) )
		{
		  *fp++ = '\0';
		  
		  if ( *fp != ',' && *fp != '\0' )
		    mn.metafields[idx] = fp;
		}
	    }
	}
      
      /* Trim last field if more fields exist */
      if ( (fp = strchr (fp, ',')) )
	*fp = '\0';
      
      /* Sanity check, source name fields must be populated */
      for (idx = 0; idx <= 3; idx++)
	{
	  if ( mn.metafields[idx] == NULL )
	    {
	      fprintf (stderr, "Error, field %d cannot be empty in metadata file line %d\n",
		       idx+1, linecount);
	      fprintf (stderr, "Perhaps a wildcard character (*) was the intention?\n");
	      
	      exit (1);
	    }
	}
      
      /* Parse and convert start time */
      if ( mn.metafields[15] )
	{
	  if ( (mn.starttime = ms_timestr2hptime (mn.metafields[15])) == HPTERROR )
	    {
	      fprintf (stderr, "Error parsing metadata start time: '%s'\n", mn.metafields[15]);
	      exit (1);
	    }
	}
      
      /* Parse and convert end time */
      if ( mn.metafields[16] )
	{
	  if ( (mn.endtime = ms_timestr2hptime (mn.metafields[16])) == HPTERROR )
	    {
	      fprintf (stderr, "Error parsing metadata end time: '%s'\n", mn.metafields[16]);
	      exit (1);
	    }
	}
      
      /* Add the metanode to the metadata list */
      if ( ! addMetaNode (&metadata, NULL, 0, &mn, sizeof(struct metanode)) )
	{
	  fprintf (stderr, "Error adding metadata fields to list\n");
	}
    }
  
  fclose (mfp);
  
  return 0;
}  /* End of readMetaData() */


/***************************************************************************
 * addMetaNode:
 *
 * Add node to the specified list.
 *
 * Return a pointer to the added node on success and NULL on error.
 ***************************************************************************/
static struct metalist *
addMetaNode (struct metalist **listroot, void *key, int keylen,
	     void *data, int datalen)
{
  struct metalist *lastlp, *newlp;
  
  if ( data == NULL )
    {
      fprintf (stderr, "addMetaNode(): No data specified\n");
      return NULL;
    }
  
  lastlp = *listroot;
  while ( lastlp != 0 )
    {
      if ( lastlp->next == 0 )
        break;
      
      lastlp = lastlp->next;
    }
  
  /* Create new sacmetanode */
  newlp = (struct metalist *) malloc (sizeof (struct metalist));
  memset (newlp, 0, sizeof (struct metalist));
  
  if ( key )
    {
      newlp->key = malloc (keylen);
      memcpy (newlp->key, key, keylen);
    }
  
  if ( data)
    {
      newlp->data = malloc (datalen);
      memcpy (newlp->data, data, datalen);
    }
  
  newlp->next = 0;
  
  if ( lastlp == 0 )
    *listroot = newlp;
  else
    lastlp->next = newlp;
  
  return newlp;
}  /* End of addMetaNode() */


/***************************************************************************
 * addFile:
 *
 * Add file to end of the global file list (filelist).
 *
 * Returns 0 on success and -1 on error.
 ***************************************************************************/
static int
addFile (char *filename)
{
  FILE *ifp = 0;
  struct filelink *newlp;
  char header[48];
  
  if ( ! filename )
    {
      fprintf (stderr, "addFile(): No file name specified\n");
      return - 1;
    }
  
  newlp = (struct filelink *) calloc (1, sizeof(struct filelink));
  
  if ( ! newlp )
    {
      fprintf (stderr, "addFile(): Error allocating memory\n");
      return -1;
    }
  
  newlp->filename = strdup(filename);
  
  if ( ! newlp->filename )
    {
      fprintf (stderr, "addFile(): Error duplicating string\n");
      return -1;      
    }
  
  /* Determine file format: open and read first 48 bytes */
  if ( (ifp = fopen (filename, "rb")) == NULL )
    {
      fprintf (stderr, "Cannot open input file: %s (%s)\n", filename, strerror(errno));
      free (newlp->filename);
      free (newlp);
      return -1;
    }
  
  if ( fread (header, 48, 1, ifp) != 1  )
    {
      fprintf (stderr, "Cannot read %s\n", filename);
      free (newlp->filename);
      free (newlp);
      return -1;
    }
  
  fclose (ifp);
  
  /* Check for Mini-SEED otherwise default to SAC */
  if ( MS_ISVALIDHEADER (header) )
    newlp->format = 1;
  else
    newlp->format = 2;
  
  /* Add new file to the end of the list */
  if ( filelisttail == 0 )
    {
      filelist = newlp;
      filelisttail = newlp;
    }
  else
    {
      filelisttail->next = newlp;
      filelisttail = newlp;
    }
  
  return 0;
}  /* End of addFile() */


/***************************************************************************
 * addListFile:
 *
 * Add files listed in the specified file to the global input file list.
 *
 * Returns count of files added on success and -1 on error.
 ***************************************************************************/
static int
addListFile (char *filename) 
{
  FILE *fp;
  char filelistent[1024];
  int filecount = 0;
  
  if ( verbose >= 1 )
    fprintf (stderr, "Reading list file '%s'\n", filename);
  
  if ( ! (fp = fopen(filename, "rb")) )
    {
      fprintf (stderr, "Cannot open list file %s: %s\n", filename, strerror(errno));
      return -1;
    }
  
  while ( fgets (filelistent, sizeof(filelistent), fp) )
    {
      char *cp;
      
      /* End string at first newline character */
      if ( (cp = strchr(filelistent, '\n')) )
        *cp = '\0';
      
      /* Skip empty lines */
      if ( ! strlen (filelistent) )
        continue;
      
      /* Skip comment lines */
      if ( *filelistent == '#' )
        continue;
      
      if ( verbose > 1 )
        fprintf (stderr, "Adding '%s' from list file\n", filelistent);
      
      if ( addFile (filelistent) )
        return -1;
      
      filecount++;
    }
  
  fclose (fp);
  
  return filecount;
}  /* End of addListFile() */


/***************************************************************************
 * addProcess:
 *
 * Add process to end of the global process list.
 ***************************************************************************/
static void
addProcess (int type, char *string1, char *string2, int ivalue1, int ivalue2,
	    double dvalue1, double dvalue2)
{
  struct proclink *lastlp, *newlp;
  char *filename = NULL;
  char *start = 0;
  char *stop = 0;
  char *tptr;
  double lpcutoff = 0.0;
  double hpcutoff = 0.0;
  double scalefactor = 0.0;
  int filetype = 0;
  int lporder = 0;
  int hporder = 0;
  int respstart = -1;
  int respstop = -1;
  int reverseflag = 0;
  int decimfactor = -1;
  int tapertype = -1;
  double taperwidth = 0.0;
  

  /* (De)Convolution */
  if ( type == PROC_CONVOLVE )
    {
      /* Parse stage limits from RESP file name */
      if ( ivalue1 == PROC_CONVRESP || ivalue1 == PROC_DECONVRESP )
	{
	  /* string1 == filename, ivalue1 == sub/filetype */
	  
	  if ( (start = strchr (string1, ':')) )
	    {
	      *start++ = '\0';
	      
	      if ( (stop = strchr (start, ':')) )
		{
		  *stop++ = '\0';
		}
	    }
	  
	  if ( start )
	    respstart = strtol (start, NULL, 10);
	  if ( stop )
	    respstop = strtol (stop, NULL, 10);
	  
	  filename = string1;
	  filetype = ivalue1;
	}
  
      /* SAC poles and zeros */
      if ( ivalue1 == PROC_CONVSAC || ivalue1 == PROC_DECONVSAC )
	{
	  /* string1 == filename, ivalue1 == sub/filetype */
	  
	  filename = string1;
	  filetype = ivalue1;
	}
    }
  /* Low-pass filter */
  if ( type == PROC_LPFILTER )
    {
      /* string1 == lpfilter, ivalue1 == reverseflag */
      
      /* Check for slash character indicating an order was specified */
      if ( (tptr = strchr (string1, '/')) )
	{
	  *tptr++ = '\0';
	  lpcutoff = strtod (string1, NULL);
	  lporder = strtoul (tptr, NULL, 10);
	}
      else
	{
	  lpcutoff = strtod (string1, NULL);
	  lporder = DEFAULT_FILTER_ORDER;
	}
      
      reverseflag = ivalue1;
    }
  
  /* High-pass filter */
  if ( type == PROC_HPFILTER )
    {
      /* string1 == hpfilter, ivalue1 == reverseflag */
      
      /* Check for slash character indicating an order was specified */
      if ( (tptr = strchr (string1, '/')) )
	{
	  *tptr++ = '\0';
	  hpcutoff = strtod (string1, NULL);
	  hporder = strtoul (tptr, NULL, 10);
	}
      else
	{
	  hpcutoff = strtod (string1, NULL);
	  hporder = DEFAULT_FILTER_ORDER;
	}
      
      reverseflag = ivalue1;
    }
  
  /* Band-pass filter */
  if ( type == PROC_BPFILTER )
    {
      /* string1 == lpfilter, string2 == hpfilter, ivalue1 == reverseflag */
      
      /* Check for slash character indicating a LP order was specified */
      if ( (tptr = strchr (string1, '/')) )
	{
	  *tptr++ = '\0';
	  lpcutoff = strtod (string1, NULL);
	  lporder = strtoul (tptr, NULL, 10);
	}
      else
	{
	  lpcutoff = strtod (string1, NULL);
	  lporder = DEFAULT_FILTER_ORDER;
	}
      
      /* Check for slash character indicating a HP order was specified */
      if ( (tptr = strchr (string2, '/')) )
	{
	  *tptr++ = '\0';
	  hpcutoff = strtod (string2, NULL);
	  hporder = strtoul (tptr, NULL, 10);
	}
      else
	{
	  hpcutoff = strtod (string2, NULL);
	  hporder = DEFAULT_FILTER_ORDER;
	}
      
      reverseflag = ivalue1;
    }
  
  /* Scale factor */
  if ( type == PROC_SCALE )
    {
      /* dvalue1 == scalefactor */
      
      scalefactor = dvalue1;
    }
  
  /* Decimation factor */
  if ( type == PROC_DECIMATE )
    {
      /* ivalue1 == decimfactor */
      
      decimfactor = ivalue1;
    }
  
  /* Taper */
  if ( type == PROC_TAPER )
    {
      /* ivalue1 == tapertype */
      /* dvalue1 == taperwidth */
      
      tapertype = ivalue1;
      taperwidth = dvalue1;
    }
  
  /* Find last process in list */
  lastlp = proclist;
  while ( lastlp != 0 )
    {
      if ( lastlp->next == 0 )
	break;
      
      lastlp = lastlp->next;
    }

  /* If this is a paired deconvolution-convolution option add operation to last entry */
  if ( type == PROC_CONVOLVE && lastlp && lastlp->type == PROC_CONVOLVE && lastlp->filetype[1] == 0 &&
       (lastlp->filetype[0] == PROC_DECONVRESP || lastlp->filetype[0] == PROC_DECONVSAC) &&
       (filetype == PROC_CONVRESP || filetype == PROC_CONVSAC) )
    {
      if ( filename )
	lastlp->filename[1] = strdup(filename);
      else
	lastlp->filename[1] = NULL;
      lastlp->filetype[1] = filetype;
    }
  /* Otherwise add a new processing entry */
  else
    {
      /* Allocate and populate new process entry */
      newlp = (struct proclink *) malloc (sizeof(struct proclink));
      
      if ( ! newlp )
	{
	  fprintf (stderr, "Error allocating memory\n");
	  exit(1);
	}
      
      newlp->type = type;
      if ( filename )
	newlp->filename[0] = strdup(filename);
      else
	newlp->filename[0] = NULL;
      newlp->filetype[0] = filetype;
      newlp->filename[1] = NULL;
      newlp->filetype[1] = 0;
      newlp->respstart = respstart;
      newlp->respstop = respstop;
      newlp->lpcutoff = lpcutoff;
      newlp->lporder = lporder;
      newlp->hpcutoff = hpcutoff;
      newlp->hporder = hporder;
      newlp->reverseflag = reverseflag;
      newlp->scalefactor = scalefactor;
      newlp->decimfactor = decimfactor;
      newlp->tapertype = tapertype;
      newlp->taperwidth = taperwidth;
      newlp->next = 0;
      
      if ( lastlp == 0 )
	proclist = newlp;
      else
	lastlp->next = newlp;
    }
  
}  /* End of addProcess() */


/***************************************************************************
 * recordHandler:
 * Saves passed records to the output file.
 ***************************************************************************/
static void
recordHandler (char *record, int reclen, void *vofp)
{
  FILE *ofp = (FILE *)vofp;
  
  if ( ofp )
    {
      if ( fwrite(record, reclen, 1, ofp) != 1 )
	{
	  fprintf (stderr, "Error writing Mini-SEED to output file\n");
	}
      
      outputbytes += reclen;
    }
}  /* End of recordHandler() */


/***************************************************************************
 * usage():
 *
 * Print the usage message.
 ***************************************************************************/
static void
usage (void)
{
  fprintf (stderr, "%s - time-series signal processor: %s\n\n", PACKAGE, VERSION);
  fprintf (stderr, "Usage: %s [options] file1 [file2] [file3] ...\n\n", PACKAGE);
  fprintf (stderr,
	   " ## Input/Output Options ##\n"
	   " -V            Report program version\n"
	   " -h            Show this usage message\n"
	   " -v            Be more verbose, multiple flags can be used\n"
	   " -s            Print a basic summary after reading all input files\n"
	   " -tt secs      Specify a time tolerance for continuous traces\n"
	   " -rt diff      Specify a sample rate tolerance for continuous traces\n"
           " -ts time      Limit Mini-SEED input to records that start after time\n"
           " -te time      Limit Mini-SEED input to records that end before time\n"
           "                 time format: 'YYYY[,DDD,HH,MM,SS,FFFFFF]' delimiters: [,:.]\n"
	   " -MSEED        Write output as Mini-SEED instead of default SAC\n"
           " -Mr bytes     Specify record length in bytes, default is autodetection\n"
           " -Me encoding  Specify encoding format of data samples for input data\n"
	   " -MR bytes     Specify record length for output Mini-SEED, default is 4096\n"
/*	   " -S            Include Blockette 1001 (hi-res sample rate) in output Mini-SEED\n"*/
	   " -ME encoding  Encoding format for SEED output data, default is 4 (floats)\n"
	   " -Sf format    Specify SAC output format (default is 2:binary)\n"
           "                 1=alpha, 2=binary (host byte order),\n"
           "                 3=binary (little-endian), 4=binary (big-endian)\n" 
	   " -m metafile   File containing station metadata for SAC output\n"
           " -msi          Convert component inclination/dip from SEED to SAC convention\n"
	   " -c channel    Set channel/component name for processed output\n"
	   " -o outfile    Specify output file instead of generating names\n"
	   " -od outdir    Specify output directory for generated file names\n"
	   "\n"
	   " ## Convolution Options ##\n"
	   " -CR respfile[:#:#] Specify SEED RESP file/dir for convolution\n"
	   " -DR respfile[:#:#] Specify SEED RESP file/dir for deconvolution\n"
	   " -CS pzfile    Specify poles and zeros file for convolution\n"
	   " -DS pzfile    Specify poles and zeros file for deconvolution\n"
	   " -ST freqs     Specify envelope for a spectrum taper for (de)convolution operations\n"
	   "                 Frequencies are specify as: 'f1/f2/f3/f4'\n"
	   " -STa dBdown   Automatically determine envelope for a spectrum taper\n"
	   " -W            Use prewhitening and dewhitening for convolution operations\n"
	   " -Wo order     Specify prediction filter order used for whitening, default: 6\n"
	   " -Rts          Use the total sensitivity in SEED RESP files instead of gains\n"
	   " -Red          Use the estimated delay in SEED RESP files\n"
	   " -Rin          Ignore channel name fields, match any net, sta, loc & chan\n"
	   " -Ru units     Specify units for SEED RESP based responses (default: DEF)\n"
	   "                 Can be DIS, VEL, ACC or DEF, DEF = no unit conversion\n"
	   "\n"
	   " ## Filter Options ##\n"
	   " -LP frequency[/order]\n"
	   " -HP frequency[/order]\n"
	   " -BP frequency[/order]:frequency[/order]\n"
	   "               Specify low-pass, high-pass and band-pass filters by defining\n"
	   "                 cut off frequencies and optionally a filter order (default 4).\n"
	   " -LP1, -HP1, -BP1: Used for single-pass filtering, possible phase distortion\n"
	   "\n"
	   " ## Other Processing Options ##\n"
	   " -D2           Differentiate using 2 point (uncentered) method\n"
	   " -IT           Integrate using trapezoidal (midpoint) method\n"
	   " -RM           Remove mean from time-series\n"
	   " -SC factor    Scale the data samples by a specified factor\n"
	   " -SI factor    Scale the data samples by inverse of specified factor\n"
	   " -DEC factor   Decimate time series by specified factor\n"
	   " -TAP width[:type]  Apply symmetric taper to time series\n"
	   "\n"
	   " file#         File of input Mini-SEED or SAC\n"
	   "\n");
}  /* End of usage() */
