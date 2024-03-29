/***************************************************************************
 * tracedsp.c - time series processor for miniSEED and SAC data
 *
 * Opens user specified files, parses the input data, applies
 * specified processing steps to the timeseries and writes the data.
 *
 * A number of processing steps are available, see the manual for details.
 *
 * Written by Chad Trabant, EarthScope Data Services
 ***************************************************************************/

#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include <libmseed.h>

#include "convolve.h"
#include "decimate.h"
#include "envelope.h"
#include "iirfilter.h"
#include "rotate.h"
#include "sacformat.h"
#include "taper.h"

#define VERSION "0.12"
#define PACKAGE "tracedsp"

/* Linkable structure to hold input file names */
struct filelink
{
  char *filename;
  int format; /* File format: 1 = miniSEED, 2 = SAC (alpha or binary) */
  struct filelink *next;
};

/* Linkable structure to hold processing parameters, normally sparse */
struct proclink
{
  int type;
  char *filename[2];
  int filetype[2];
  int respstart;
  int respstop;
  double lpcutoff;
  int lporder;
  double hpcutoff;
  int hporder;
  int reverseflag;
  double scalefactor;
  int decimfactor;
  double taperwidth;
  int tapertype;
  double *coefficients;
  int coefficientcount;
  char rotateENZ[3];   /* Component names for Z,N,E */
  char rotatedENZ[3];  /* Final component names */
  double rotations[2]; /* 0: Azimuth, 1: Incident angle */
  struct proclink *next;
};

/* Additional segment details, stored at MS3TraceSeg->prvtptr  */
struct segdetails
{
  struct SACHeader *sacheader;
  int procerror;
  int rotated;
};

#define PROC_STATS 1
#define PROC_LPFILTER 2
#define PROC_HPFILTER 3
#define PROC_BPFILTER 4
#define PROC_CONVOLVE 5
#define PROC_CONVRESP 6   // Subtype of CONVOLVE
#define PROC_DECONVRESP 7 // Subtype of CONVOLVE
#define PROC_CONVSAC 8    // Subtype of CONVOLVE
#define PROC_DECONVSAC 9  // Subtype of CONVOLVE
#define PROC_DIFF2 10
#define PROC_INTTRAP 11
#define PROC_RMEAN 12
#define PROC_SCALE 13
#define PROC_DECIMATE 14
#define PROC_TAPER 15
#define PROC_POLYNOMIALM 16
#define PROC_ENVELOPE 17
#define PROC_DATATRIM 18
#define PROC_ROTATE 19

/* Maximum processing log size in bytes, 10 MB */
#define MAXPROCLOG 10485760

/* Default order of high/low pass filter */
#define DEFAULT_FILTER_ORDER 4

/* Maximum number of metadata fields per line */
#define MAXMETAFIELDS 17

/* Linkable structure to hold station metadata */
struct metalist
{
  char *key;
  char *data;
  struct metalist *next;
};

/* Structure for metadata */
struct metanode
{
  char *metafields[MAXMETAFIELDS];
  nstime_t starttime;
  nstime_t endtime;
};

static int procFilter (MS3TraceSeg *seg, struct proclink *plp);
static int procConvolve (MS3TraceID *id, MS3TraceSeg *seg, struct proclink *plp);
static int procDiff2 (MS3TraceSeg *seg, struct proclink *plp);
static int procIntTrap (MS3TraceSeg *seg, struct proclink *plp);
static int procRMean (MS3TraceSeg *seg, struct proclink *plp);
static int procScale (MS3TraceSeg *seg, struct proclink *plp);
static int procDecimate (MS3TraceSeg *seg, struct proclink *plp);
static int procTaper (MS3TraceSeg *seg, struct proclink *plp);
static int procPolynomialM (MS3TraceSeg *seg, struct proclink *plp);
static int procEnvelope (MS3TraceSeg *seg, struct proclink *plp);
static int procDataTrim (MS3TraceSeg *seg, nstime_t lateststart, nstime_t earliestend);
static int procRotate (MS3TraceList *mstl, MS3TraceID *tid, MS3TraceSeg *tseg, struct proclink *plp);

static int64_t readMSEED (char *mseedfile, MS3TraceList *mstl);
static int64_t readSAC (char *sacfile, MS3TraceList *mstl);
static int parseSAC (FILE *ifp, struct SACHeader *sh, float **data, int format,
                     int verbose, char *sacfile);
static int readBinaryHeaderSAC (FILE *ifp, struct SACHeader *sh, int *format,
                                int *swapflag, int verbose, char *sacfile);
static int readBinaryDataSAC (FILE *ifp, float *data, int datacnt,
                              int swapflag, int verbose, char *sacfile);
static int readAlphaHeaderSAC (FILE *ifp, struct SACHeader *sh);
static int readAlphaDataSAC (FILE *ifp, float *data, int datacnt);

static int writeMSEED (MS3TraceID *id, MS3TraceSeg *seg, char *outputfile);
static int writeSAC (MS3TraceID *id, MS3TraceSeg *seg, int format, char *outputfile);
static int writeBinarySAC (struct SACHeader *sh, float *fdata, int npts, char *outfile);
static int writeAlphaSAC (struct SACHeader *sh, float *fdata, int npts, char *outfile);
static int insertSACMetaData (struct SACHeader *sh, nstime_t sacstarttime);
static int swapSACHeader (struct SACHeader *sh);

static int addToString (char **string, char *add, char *delim, int where, int maxlen);
static int addToProcLog (const char *format, ...);

static int calcStats (void *input, char inputtype, int length);
static int differentiate2 (void *input, char inputtype, int length,
                           double rate, void *output);
static int integrateTrap (void *input, char inputtype, int length,
                          double halfstep, void *output);
static int applyPolynomialM (void *input, char inputtype, int length,
                             double *coeff, int numcoeff, void *output);

static int convertSamples (MS3TraceSeg *seg, char type);
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
static char *procDescription (int proctype);

static void recordHandler (char *record, int reclen, void *vofp);
static void usage (void);

static flag verbose       = 0;
static flag basicsum      = 0;    /* Controls printing of basic summary */
static int packreclen     = -1;   /* miniSEED record length for output data */
static int packencoding   = 4;    /* miniSEED encoding format for output data */
static int msformat       = 2;    /* miniSEED format version */
static flag dataformat    = 2;    /* 0 = No output, 1 = miniSEED, 2 = SAC */
static flag informat      = 0;    /* 0 = auto detect, 1 = miniSEED, 2 = SAC */
static flag sacinformat   = 0;    /* 0=auto, 1=alpha, 2=binary (host), 3=binary (LE), 4=binary (BE) */
static flag sacoutformat  = 2;    /* 1=alpha, 2=binary (host), 3=binary (LE), 4=binary (BE) */
static char *sacnet       = 0;    /* SAC network code override */
static char *sacloc       = 0;    /* SAC location ID override */
static char *metadatafile = 0;    /* File containing metadata for output (SAC, etc.) */
static char *proclog      = 0;    /* Detailed processing log as a string */
static char *proclogfile  = 0;    /* Detailed processing log file */
static int prewhiten      = 0;    /* Prewhitening for [de]convolution, predictor order */
static double *freqlimit  = 0;    /* Frequency limits for deconvolution */
static double lcdBdown    = -1.0; /* Lower corner dB down cutoff */
static double ucdBdown    = -1.0; /* Upper corner dB down cutoff */
static flag resptotalsens = 0;    /* Controls evalresp's usage of total sensitivity in RESP */
static flag respusedelay  = 0;    /* Controls evalresp's usage of estimated delay in RESP */
static flag respusename   = 1;    /* Controls evalresp's matching of Net, Sta, Loc and Chan */
static char *respunits    = 0;    /* Controls units for evalresp calculated responses */
static nstime_t starttime = NSTUNSET;
static nstime_t endtime   = NSTUNSET;
static char *outputfile   = 0; /* Output file name */
static char *outputdir    = 0; /* Output base directory */
static int outputbytes    = 0; /* Bytes written to output file */
static char *channel      = 0; /* Forced output channel/component */

/* Root and tail of input file name list */
struct filelink *filelist     = 0;
struct filelink *filelisttail = 0;

/* Root of processing actions list */
struct proclink *proclist = 0;

/* A list of station and coordinates */
struct metalist *metadata = 0;
static int seedinc        = 0;

int
main (int argc, char **argv)
{
  struct filelink *flp, *nextflp;
  struct proclink *plp, *nextplp;
  MS3TraceList *mstl = NULL;
  MS3TraceID *id     = NULL;
  MS3TraceSeg *seg   = NULL;
  int errflag        = 0;
  int64_t totalsamps = 0;
  int64_t sampsread;
  int totalfiles       = 0;
  nstime_t lateststart = NSTUNSET;
  nstime_t earliestend = NSTUNSET;
  char stime[50];
  char etime[50];

  /* Process given parameters (command line and parameter file) */
  if (parameterProc (argc, argv) < 0)
    return -1;

  /* Initialize MS3TraceList */
  mstl = mstl3_init (NULL);

  /* Loop over input file list */
  flp = filelist;
  while (flp != 0)
  {
    /* Read miniSEED file */
    if (flp->format == 1)
    {
      if (verbose >= 2)
        fprintf (stderr, "Processing miniSEED: %s\n", flp->filename);

      sampsread = readMSEED (flp->filename, mstl);

      if (sampsread < 0)
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
      if (verbose >= 2)
        fprintf (stderr, "Processing SAC: %s\n", flp->filename);

      sampsread = readSAC (flp->filename, mstl);

      if (sampsread < 0)
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

  if (basicsum)
    printf ("Input Files: %d, Samples: %lld\n", totalfiles, (long long int)totalsamps);

  /* Determine latest start and earliest end times across all channels */
  id = mstl->traces.next[0];
  while (id)
  {
    if (lateststart == NSTUNSET || lateststart < id->earliest)
      lateststart = id->earliest;

    if (earliestend == NSTUNSET || earliestend > id->latest)
      earliestend = id->latest;

    id = id->next[0];
  }

  /* Loop through process list, apply to each time series segment */
  plp = proclist;
  while (plp && !errflag)
  {
    id = mstl->traces.next[0];
    while (id)
    {
      seg = id->first;
      while (seg)
      {
        ms_nstime2timestr (seg->starttime, stime, ISOMONTHDAY_DOY_Z, NANO_MICRO);
        ms_nstime2timestr (seg->endtime, etime, ISOMONTHDAY_DOY_Z, NANO_MICRO);
        addToProcLog ("Processing %s [%s - %s]: %s", id->sid, stime, etime, procDescription (plp->type));

        if (plp->type == PROC_STATS)
        {
          if (calcStats (seg->datasamples, seg->sampletype, seg->numsamples))
          {
            fprintf (stderr, "Error calculating statistics\n");
            errflag = -1;
            ((struct segdetails *)seg->prvtptr)->procerror = errflag;
            break;
          }
        }
        else if (plp->type == PROC_LPFILTER || plp->type == PROC_HPFILTER || plp->type == PROC_BPFILTER)
        {
          if (procFilter (seg, plp))
          {
            fprintf (stderr, "Error applying filter for %s\n", id->sid);
            errflag = -1;
            ((struct segdetails *)seg->prvtptr)->procerror = errflag;
            break;
          }
        }
        else if (plp->type == PROC_CONVOLVE)
        {
          if (procConvolve (id, seg, plp))
          {
            fprintf (stderr, "Error (de)convolving response (%s - %s) from %s\n",
                     (plp->filename[0]) ? plp->filename[0] : "None",
                     (plp->filename[1]) ? plp->filename[1] : "None", id->sid);
            errflag = -1;
            ((struct segdetails *)seg->prvtptr)->procerror = errflag;
            break;
          }
        }
        else if (plp->type == PROC_DIFF2)
        {
          if (procDiff2 (seg, plp))
          {
            fprintf (stderr, "Error differentiating time series for %s\n",
                     id->sid);
            errflag = -1;
            ((struct segdetails *)seg->prvtptr)->procerror = errflag;
            break;
          }
        }
        else if (plp->type == PROC_INTTRAP)
        {
          if (procIntTrap (seg, plp))
          {
            fprintf (stderr, "Error integrating time series for %s\n",
                     id->sid);
            errflag = -1;
            ((struct segdetails *)seg->prvtptr)->procerror = errflag;
            break;
          }
        }
        else if (plp->type == PROC_RMEAN)
        {
          if (procRMean (seg, plp))
          {
            fprintf (stderr, "Error removing mean from time series %s\n",
                     id->sid);
            errflag = -1;
            ((struct segdetails *)seg->prvtptr)->procerror = errflag;
            break;
          }
        }
        else if (plp->type == PROC_SCALE)
        {
          if (procScale (seg, plp))
          {
            fprintf (stderr, "Error scaling values of time series %s\n",
                     id->sid);
            errflag = -1;
            ((struct segdetails *)seg->prvtptr)->procerror = errflag;
            break;
          }
        }
        else if (plp->type == PROC_DECIMATE)
        {
          if (procDecimate (seg, plp))
          {
            fprintf (stderr, "Error decimating time series %s\n",
                     id->sid);
            errflag = -1;
            ((struct segdetails *)seg->prvtptr)->procerror = errflag;
            break;
          }
        }
        else if (plp->type == PROC_TAPER)
        {
          if (procTaper (seg, plp))
          {
            fprintf (stderr, "Error tapering time series %s\n",
                     id->sid);
            errflag = -1;
            ((struct segdetails *)seg->prvtptr)->procerror = errflag;
            break;
          }
        }
        else if (plp->type == PROC_POLYNOMIALM)
        {
          if (procPolynomialM (seg, plp))
          {
            fprintf (stderr, "Error applying Maclaurin polynomial to time series %s\n",
                     id->sid);
            errflag = -1;
            ((struct segdetails *)seg->prvtptr)->procerror = errflag;
            break;
          }
        }
        else if (plp->type == PROC_ENVELOPE)
        {
          if (procEnvelope (seg, plp))
          {
            fprintf (stderr, "Error calculating envelope of time series %s\n",
                     id->sid);
            errflag = -1;
            ((struct segdetails *)seg->prvtptr)->procerror = errflag;
            break;
          }
        }
        else if (plp->type == PROC_DATATRIM)
        {
          if (procDataTrim (seg, lateststart, earliestend))
          {
            fprintf (stderr, "Error synchronizing time series windows %s\n",
                     id->sid);
            errflag = -1;
            ((struct segdetails *)seg->prvtptr)->procerror = errflag;
            break;
          }
        }
        else if (plp->type == PROC_ROTATE)
        {
          if (!((struct segdetails *)(seg->prvtptr))->rotated)
            if (procRotate (mstl, id, seg, plp) < 0)
            {
              fprintf (stderr, "Error rotating time series %s\n",
                       id->sid);
              errflag = -1;
              ((struct segdetails *)seg->prvtptr)->procerror = errflag;
              break;
            }
        }

        seg = seg->next;
      }

      id = id->next[0];
    }

    plp = plp->next;
  }

  /* Write each segment with no processing errors */
  id = mstl->traces.next[0];
  while (id)
  {
    seg = id->first;
    while (seg)
    {
      if (!((struct segdetails *)seg->prvtptr)->procerror)
      {
        ms_nstime2timestr (seg->starttime, stime, ISOMONTHDAY_DOY_Z, NANO_MICRO);
        ms_nstime2timestr (seg->endtime, etime, ISOMONTHDAY_DOY_Z, NANO_MICRO);
        addToProcLog ("Writing %s [%s - %s]", id->sid, stime, etime);

        if (dataformat == 1)
        {
          if (writeMSEED (id, seg, outputfile) < 0)
            fprintf (stderr, "Error writing miniSEED\n");
        }
        else if (dataformat == 2)
        {
          if (writeSAC (id, seg, sacoutformat, outputfile) < 0)
            fprintf (stderr, "Error writing SAC\n");
        }
      }
      seg = seg->next;
    }
    id = id->next[0];
  }

  /* Write process log to output file */
  if (proclog && proclogfile)
  {
    FILE *plf;

    if (verbose)
      fprintf (stderr, "Writing process log file\n");

    if ((plf = fopen (proclogfile, "wb")) == NULL)
    {
      fprintf (stderr, "Error opening process log file: %s\n", strerror (errno));
    }
    else
    {
      if (fprintf (plf, "%s\n", proclog) < 0)
      {
        fprintf (stderr, "Error writing process log file\n");
      }

      fclose (plf);
    }
  }

  /* Make sure everything is cleaned up */
  mstl3_free (&mstl, 1);

  flp = filelist;
  while (flp)
  {
    nextflp = flp->next;
    if (flp->filename)
      free (flp->filename);
    free (flp);
    flp = nextflp;
  }

  plp = proclist;
  while (plp)
  {
    nextplp = plp->next;
    if (plp->filename[0])
      free (plp->filename[0]);
    if (plp->filename[1])
      free (plp->filename[1]);
    free (plp);
    plp = nextplp;
  }

  return errflag;
} /* End of main() */

/***************************************************************************
 * procFilter:
 *
 * Apply filter to the MS3TraceSeg, results are requested to be
 * returned as doubles.
 *
 * Returns 0 on success and -1 on error.
 ***************************************************************************/
static int
procFilter (MS3TraceSeg *seg, struct proclink *plp)
{
  void *datasamples = 0;

  if (!seg || !plp)
    return -1;

  /* Filter trace */
  if (plp->lporder || plp->hporder)
  {
    if (plp->lporder && !plp->hporder)
      addToProcLog ("Low-pass filter cutoff: %f, order: %d",
                    plp->lpcutoff, plp->lporder);
    else if (!plp->lporder && plp->hporder)
      addToProcLog ("High-pass filter cutoff: %f, order: %d",
                    plp->hpcutoff, plp->hporder);
    else
      addToProcLog ("Band-pass filter HP cutoff: %f, order: %d => LP cutoff: %f, order: %d",
                    plp->lpcutoff, plp->lporder, plp->hpcutoff, plp->hporder);

    /* Apply the filter */
    if (iirfilter (seg->datasamples, seg->sampletype, seg->numsamples, plp->reverseflag,
                   &datasamples, 'd', plp->hporder, plp->hpcutoff,
                   plp->lporder, plp->lpcutoff, seg->samprate, verbose))
    {
      if (datasamples)
        free (datasamples);
      return -1;
    }
    else
    {
      /* Free the original buffer and replace it with the filtered */
      if (seg->datasamples)
        free (seg->datasamples);

      seg->datasamples = datasamples;
      seg->sampletype  = 'd';
    }
  }

  return 0;
} /* End of procFilter() */

/***************************************************************************
 * procConvolve:
 *
 * Convolve or deconvolve supplied response(s) from signal.  The
 * convolution routines require doubless so this routine will always
 * convert the supplied data to doubles and return the results as
 * doubles.
 *
 * Returns 0 on success and non-zero on error.
 ***************************************************************************/
static int
procConvolve (MS3TraceID *id, MS3TraceSeg *seg, struct proclink *plp)
{
  int idx;
  int retval = 0;
  int nfft;
  int nfreqs;
  double delfreq;

  double *freqs = NULL;
  double *creal = NULL;
  double *cimag = NULL;
  double *dreal = NULL;
  double *dimag = NULL;
  double *xreal = NULL;
  double *ximag = NULL;

  char net[11] = {0};
  char sta[11] = {0};
  char loc[11] = {0};
  char chan[11] = {0};

  if (!seg || !plp)
    return -1;

  /* Decompose Source ID into separate codes */
  ms_sid2nslc (id->sid, net, sta, loc, chan);

  /* Convert samples to doubles if needed */
  if (seg->sampletype != 'd')
  {
    if (convertSamples (seg, 'd'))
      return -1;
  }

  /* Calculate common parameters */
  nfft    = next2 (seg->numsamples);
  nfreqs  = nfft / 2 + 1;
  delfreq = seg->samprate / nfft;

  /* Loop over two potential response files */
  for (idx = 0; idx < 2; idx++)
  {
    if (plp->filename[idx] == NULL)
      break;

    /* Calculate response functions from SEED RESP */
    if (plp->filetype[idx] == PROC_CONVRESP || plp->filetype[idx] == PROC_DECONVRESP)
    {
      if (plp->filetype[idx] == PROC_CONVRESP)
        addToProcLog ("Convolving with SEED RESP response '%s'", plp->filename[idx]);
      else
        addToProcLog ("Deconvolving SEED RESP response '%s'", plp->filename[idx]);

      if (respusename)
      {
        retval = calcfr_resp (nfreqs, delfreq,
                              net, sta, loc, chan,
                              plp->respstart, plp->respstop,
                              respunits, MS_NSTIME2EPOCH (seg->starttime), respusedelay,
                              plp->filename[idx], resptotalsens,
                              &freqs, &xreal, &ximag, verbose);
      }
      else
      {
        retval = calcfr_resp (nfreqs, delfreq, "*", "*", "*", "*",
                              plp->respstart, plp->respstop,
                              respunits, MS_NSTIME2EPOCH (seg->starttime), respusedelay,
                              plp->filename[idx], resptotalsens,
                              &freqs, &xreal, &ximag, verbose);
      }

      if (retval)
      {
        fprintf (stderr, "Error determining frequency response\n");
        return -1;
      }

      /* Assign to convolution or deconvolution array */
      if (plp->filetype[idx] == PROC_CONVRESP)
      {
        creal = xreal;
        xreal = NULL;
        cimag = ximag;
        ximag = NULL;
      }
      else
      {
        dreal = xreal;
        xreal = NULL;
        dimag = ximag;
        ximag = NULL;
      }
    }
    /* Calculate response functions from SAC P&Zs */
    else if (plp->filetype[idx] == PROC_CONVSAC || plp->filetype[idx] == PROC_DECONVSAC)
    {
      if (plp->filetype[idx] == PROC_CONVSAC)
        addToProcLog ("Convolving with SAC Poles & Zeros response '%s'", plp->filename[idx]);
      else
        addToProcLog ("Deconvolving SAC Poles & Zeros response '%s'", plp->filename[idx]);

      retval = calcfr_sac (nfreqs, delfreq, plp->filename[idx],
                           &freqs, &xreal, &ximag, verbose);

      /* Assign to convolution or deconvolution array */
      if (plp->filetype[idx] == PROC_CONVSAC)
      {
        creal = xreal;
        xreal = NULL;
        cimag = ximag;
        ximag = NULL;
      }
      else
      {
        dreal = xreal;
        xreal = NULL;
        dimag = ximag;
        ximag = NULL;
      }
    }
  }

  if (plp->filetype[0] && plp->filetype[1])
  {
    addToProcLog ("Deconvolution-convolution operations combined into a single transfer operation");
  }

  /* Check if frequency limit parameters need to be calculated for deconvolution */
  if (dreal && dimag && freqlimit)
  {
    if ((freqlimit[0] == -1.0 || freqlimit[1] == -1.0 ||
         freqlimit[2] == -1.0 || freqlimit[3] == -1.0))
    {
      /* Determine frequency limits for deconvolution response */
      if (findtaper (freqlimit, dreal, dimag, nfreqs, delfreq, lcdBdown, ucdBdown))
      {
        fprintf (stderr, "Error determining deconvolution frequency limit parameters\n");
        return -1;
      }
    }

    if (plp->filetype[0] && plp->filetype[1])
      addToProcLog ("Transfer frequency limits (Hz): %g/%g => %g/%g [cutoffs %g/%g]",
                    freqlimit[0], freqlimit[1], freqlimit[2], freqlimit[3],
                    lcdBdown, ucdBdown);
    else
      addToProcLog ("Deconvolution frequency limits (Hz): %g/%g => %g/%g [cutoffs %g/%g]",
                    freqlimit[0], freqlimit[1], freqlimit[2], freqlimit[3],
                    lcdBdown, ucdBdown);
  }

  /* Perform convolution, deconvolution or both */
  retval = convolve (seg->datasamples, seg->numsamples, 1.0 / seg->samprate, nfreqs, nfft,
                     freqs, creal, cimag, dreal, dimag, freqlimit, &prewhiten, verbose);

  /* Free response function arrays */
  if (freqs)
    free (freqs);
  if (creal)
    free (creal);
  if (cimag)
    free (cimag);
  if (dreal)
    free (dreal);
  if (dimag)
    free (dimag);

  return retval;
} /* End of procConvolve() */

/***************************************************************************
 * procDiff2:
 *
 * Prepare for an uncentered, 2-point differentiation.  Integer
 * samples will be converted to doubles.
 *
 * Returns 0 on success and non-zero on error.
 ***************************************************************************/
static int
procDiff2 (MS3TraceSeg *seg, struct proclink *plp)
{
  int count;

  if (!seg || !plp)
    return -1;

  /* Convert integer samples to doubles */
  if (seg->sampletype == 'i')
  {
    if (convertSamples (seg, 'd'))
      return -1;
  }
  else if (seg->sampletype != 'f' && seg->sampletype != 'd')
  {
    fprintf (stderr, "procDiff2(): unsupported sample type: %c\n",
             seg->sampletype);
    return -1;
  }

  /* Perform differentiation */
  addToProcLog ("Differentiating (%c) time series", seg->sampletype);

  count = differentiate2 (seg->datasamples, seg->sampletype, seg->numsamples,
                          seg->samprate, seg->datasamples);

  if (count < 0)
    return -1;

  /* Update sample counts */
  seg->samplecnt  = count;
  seg->numsamples = count;

  /* Shift the start time by 1/2 the sample interval */
  seg->starttime += (0.5 / seg->samprate) * NSTMODULUS;

  return 0;
} /* End of procDiff2() */

/***************************************************************************
 * procIntTrap:
 *
 * Prepare for integration using the trapezoidal (midpoint) method.
 * Integer samples will be converted to doubles.
 *
 * Returns 0 on success and non-zero on error.
 ***************************************************************************/
static int
procIntTrap (MS3TraceSeg *seg, struct proclink *plp)
{
  int count;

  if (!seg || !plp)
    return -1;

  /* Convert integer samples to doubles */
  if (seg->sampletype == 'i')
  {
    if (convertSamples (seg, 'd'))
      return -1;
  }
  else if (seg->sampletype != 'f' && seg->sampletype != 'd')
  {
    fprintf (stderr, "procIntTrap(): unsupported sample type: %c\n",
             seg->sampletype);
    return -1;
  }

  /* Perform integration */
  addToProcLog ("Integrating (%c) time series", seg->sampletype);

  count = integrateTrap (seg->datasamples, seg->sampletype, seg->numsamples,
                         (0.5 / seg->samprate), seg->datasamples);

  if (count < 0)
    return -1;

  /* Update sample counts */
  seg->samplecnt  = count;
  seg->numsamples = count;

  /* Shift the start time by 1/2 the sample interval */
  seg->starttime += (0.5 / seg->samprate) * NSTMODULUS;

  return 0;
} /* End of procIntTrap() */

/***************************************************************************
 * procRMean:
 *
 * Removes the mean from a time series.  Integer samples will be
 * converted to doubles.
 *
 * Returns 0 on success and non-zero on error.
 ***************************************************************************/
static int
procRMean (MS3TraceSeg *seg, struct proclink *plp)
{
  int idx;
  double Mean, pM;
  float *fdata;
  double *ddata;

  if (!seg || !plp)
    return -1;

  /* Convert integer samples to doubles */
  if (seg->sampletype == 'i')
  {
    if (convertSamples (seg, 'd'))
      return -1;
  }
  else if (seg->sampletype != 'f' && seg->sampletype != 'd')
  {
    fprintf (stderr, "procRMean(): unsupported sample type: %c\n",
             seg->sampletype);
    return -1;
  }

  fdata = seg->datasamples;
  ddata = seg->datasamples;

  /* Find mean value */
  if (seg->sampletype == 'f')
  {
    /* Calculate running mean */
    pM = Mean = *fdata;
    for (idx = 1; idx < seg->numsamples; idx++)
    {
      Mean = pM + (*(fdata + idx) - pM) / (idx + 1);
      pM   = Mean;
    }

    addToProcLog ("Removing mean of %g from time series", Mean);

    /* Remove mean */
    for (idx = 0; idx < seg->numsamples; idx++)
    {
      *(fdata + idx) -= Mean;
    }
  }
  else if (seg->sampletype == 'd')
  {
    /* Calculate running mean */
    pM = Mean = *ddata;
    for (idx = 1; idx < seg->numsamples; idx++)
    {
      Mean = pM + (*(ddata + idx) - pM) / (idx + 1);
      pM   = Mean;
    }

    addToProcLog ("Removing mean of %g from time series", Mean);

    /* Remove mean */
    for (idx = 0; idx < seg->numsamples; idx++)
    {
      *(ddata + idx) -= Mean;
    }
  }

  return 0;
} /* End of procRMean() */

/***************************************************************************
 * procScale:
 *
 * Scales all data samples in a time series.  Integer samples will be
 * converted to doubles.
 *
 * Returns 0 on success and non-zero on error.
 ***************************************************************************/
static int
procScale (MS3TraceSeg *seg, struct proclink *plp)
{
  int idx;
  float *fdata;
  double *ddata;

  if (!seg || !plp)
    return -1;

  /* Convert integer samples to doubles */
  if (seg->sampletype == 'i')
  {
    if (convertSamples (seg, 'd'))
      return -1;
  }
  else if (seg->sampletype != 'f' && seg->sampletype != 'd')
  {
    fprintf (stderr, "procScale(): unsupported sample type: %c\n",
             seg->sampletype);
    return -1;
  }

  fdata = seg->datasamples;
  ddata = seg->datasamples;

  /* Scale sample values */
  if (seg->sampletype == 'f')
  {
    addToProcLog ("Scaling time series by %g", plp->scalefactor);

    /* Scale samples */
    for (idx = 0; idx < seg->numsamples; idx++)
    {
      *(fdata + idx) *= plp->scalefactor;
    }
  }
  else if (seg->sampletype == 'd')
  {
    addToProcLog ("Scaling time series by %g", plp->scalefactor);

    /* Scale samples */
    for (idx = 0; idx < seg->numsamples; idx++)
    {
      *(ddata + idx) *= plp->scalefactor;
    }
  }

  return 0;
} /* End of procScale() */

/***************************************************************************
 * procDecimate:
 *
 * Decimates the time series by a given factor.  Integer and float
 * samples will be converted to doubles.
 *
 * Returns 0 on success and non-zero on error.
 ***************************************************************************/
static int
procDecimate (MS3TraceSeg *seg, struct proclink *plp)
{
  int numsamples;

  if (!seg || !plp)
    return -1;

  /* Convert samples to doubles */
  if (seg->sampletype != 'd')
  {
    if (convertSamples (seg, 'd'))
      return -1;
  }

  addToProcLog ("Decimating time series by a factor of %d (%g -> %g sps)",
                plp->decimfactor, seg->samprate, seg->samprate / plp->decimfactor);

  /* Perform the decimation and filtering */
  numsamples = decimate (seg->datasamples, seg->numsamples, plp->decimfactor, NULL, -1, -1);

  if (numsamples >= 0)
  {
    /* Adjust sample rate, sample count, end time and sample buffer size */
    seg->samprate /= plp->decimfactor;
    seg->samplecnt  = numsamples;
    seg->numsamples = numsamples;
    seg->endtime    = seg->starttime +
                   (((double)(seg->numsamples - 1) / seg->samprate * NSTMODULUS) + 0.5);

    if (!(seg->datasamples = realloc (seg->datasamples, numsamples * sizeof (double))))
    {
      fprintf (stderr, "procDecimate(): Error reallocating sample buffer\n");
      return -1;
    }
  }
  else
  {
    fprintf (stderr, "procDecimate(): Error decimating time series\n");
    return -1;
  }

  return 0;
} /* End of procDecimate() */

/***************************************************************************
 * procTaper:
 *
 * Tapers the time series using a specified type and width.  Integer
 * and float samples will be converted to doubles.
 *
 * Returns 0 on success and non-zero on error.
 ***************************************************************************/
static int
procTaper (MS3TraceSeg *seg, struct proclink *plp)
{
  int retval;
  char *typestr = "";

  if (!seg || !plp)
    return -1;

  /* Convert integer and float samples to doubles */
  if (seg->sampletype != 'd')
  {
    if (convertSamples (seg, 'd'))
      return -1;
  }

  /* Perform the tapering */
  retval = taper (seg->datasamples, seg->numsamples, plp->taperwidth, plp->tapertype);

  if (retval < 0)
  {
    fprintf (stderr, "procTaper(): Error tapering time series\n");
    return -1;
  }

  switch (plp->tapertype)
  {
  case TAPER_HANNING:
    typestr = "Hanning";
    break;
  case TAPER_HAMMING:
    typestr = "Hamming";
    break;
  case TAPER_COSINE:
    typestr = "Cosine";
    break;
  }

  addToProcLog ("Tapering time series using width %g (%d samples) type: %s",
                plp->taperwidth, retval, typestr);

  return 0;
} /* End of procTaper() */

/***************************************************************************
 * procPolynomialM:
 *
 * Apply a Maclaurin type polynomial to the data series.  Integer and
 * float samples will be converted to doubles.
 *
 * Returns 0 on success and non-zero on error.
 ***************************************************************************/
static int
procPolynomialM (MS3TraceSeg *seg, struct proclink *plp)
{
  int idx;
  int retval;
  char coeffval[20];
  char *coeffstr = NULL;

  if (!seg || !plp)
    return -1;

  if (!plp->coefficients)
    return -1;

  /* Convert samples to doubles */
  if (seg->sampletype != 'd')
  {
    if (convertSamples (seg, 'd'))
      return -1;
  }

  /* Build coefficient string to print */
  for (idx = 0; idx < plp->coefficientcount; idx++)
  {
    snprintf (coeffval, sizeof (coeffval), "%g", plp->coefficients[idx]);

    retval = addToString (&coeffstr, coeffval, ",", 0, 96);

    if (retval < 0)
    {
      if (retval == -2)
        addToString (&coeffstr, "...", " ", 0, 100);

      break;
    }
  }

  addToProcLog ("Applying Maclaurin type polynomial to data, cofficients: %s", coeffstr);

  if (coeffstr)
    free (coeffstr);

  /* Calculate envelope, replacing original data array */
  retval = applyPolynomialM (seg->datasamples, seg->sampletype, seg->numsamples,
                             plp->coefficients, plp->coefficientcount, seg->datasamples);

  if (retval < 0)
  {
    fprintf (stderr, "procPolynomialM(): Error applying polynomial to time series\n");
    return -1;
  }

  return 0;
} /* End of procPolynomialM() */

/***************************************************************************
 * procEnvelope:
 *
 * Calculates envelope of time series data.  Integer and float samples
 * will be converted to doubles.
 *
 * Returns 0 on success and non-zero on error.
 ***************************************************************************/
static int
procEnvelope (MS3TraceSeg *seg, struct proclink *plp)
{
  int retval;

  if (!seg || !plp)
    return -1;

  /* Convert samples to doubles */
  if (seg->sampletype != 'd')
  {
    if (convertSamples (seg, 'd'))
      return -1;
  }

  addToProcLog ("Calculating envelope of time series");

  /* Calculate envelope, replacing original data array */
  retval = envelope (seg->datasamples, seg->numsamples);

  if (retval < 0)
  {
    fprintf (stderr, "procEnvelope(): Error caluclating envelope of time series\n");
    return -1;
  }

  return 0;
} /* End of procEnvelope() */

/***************************************************************************
 * procDataTrim:
 *
 * Synchronize data start and end times by:
 *  1) trimming start times for each segement to the latest start time
 *  2) trimming end times for each segement to the earliest end time
 *
 * Returns 0 on success and non-zero on error.
 ***************************************************************************/
static int
procDataTrim (MS3TraceSeg *seg, nstime_t lateststart, nstime_t earliestend)
{
  int64_t trimcount;
  int samplesize;
  void *datasamples;

  if (!seg)
    return -1;

  if (lateststart == NSTUNSET || earliestend == NSTUNSET)
    return -1;

  /* Skip segments that do not have integer, float or double types */
  if (seg->sampletype != 'i' && seg->sampletype != 'f' && seg->sampletype != 'd')
    return 0;

  samplesize = ms_samplesize (seg->sampletype);

  /* Trim samples from beginning of segment if earlier than latest */
  if (seg->starttime < lateststart)
  {
    trimcount = (int64_t) ((double)MS_NSTIME2EPOCH ((lateststart - seg->starttime)) * seg->samprate + 0.5);

    if (trimcount > 0 && trimcount < seg->numsamples)
    {
      addToProcLog ("Trimming %lld samples from beginning of trace", (long long)trimcount);

      memmove (seg->datasamples, seg->datasamples + (trimcount * samplesize), (trimcount * samplesize));
      datasamples = realloc (seg->datasamples, (seg->numsamples - trimcount) * samplesize);

      if (!datasamples)
      {
        fprintf (stderr, "procDataTrim(): Error reallocating sample buffer\n");
        return -1;
      }

      seg->datasamples = datasamples;
      seg->starttime   = seg->starttime + MS_EPOCH2NSTIME ((trimcount / seg->samprate));
      seg->numsamples -= trimcount;
      seg->samplecnt -= trimcount;
    }
  }

  /* Trim samples from end of segment if later than earliest */
  if (seg->endtime > earliestend)
  {
    trimcount = (int64_t) ((double)MS_NSTIME2EPOCH ((seg->endtime - earliestend)) * seg->samprate + 0.5);

    if (trimcount > 0 && trimcount < seg->numsamples)
    {
      addToProcLog ("Trimming %lld samples from end of trace", (long long)trimcount);

      datasamples = realloc (seg->datasamples, (seg->numsamples - trimcount) * samplesize);

      if (!datasamples)
      {
        fprintf (stderr, "procDataTrim(): Error reallocating sample buffer\n");
        return -1;
      }

      seg->datasamples = datasamples;
      seg->endtime     = seg->endtime - MS_EPOCH2NSTIME ((trimcount / seg->samprate));
      seg->numsamples -= trimcount;
      seg->samplecnt -= trimcount;
    }
  }

  return 0;
} /* End of procDataTrim() */

/***************************************************************************
 * procRotate:
 *
 * Rotate seismogram sets.
 *
 * Returns the number of segments rotated, 0 on no operation and negative on error.
 ***************************************************************************/
static int
procRotate (MS3TraceList *mstl, MS3TraceID *tid, MS3TraceSeg *tseg, struct proclink *plp)
{
  MS3TraceID *id       = NULL;
  MS3TraceID *ENZid[3] = {NULL, NULL, NULL};

  MS3TraceSeg *seg       = NULL;
  MS3TraceSeg *ENZseg[3] = {NULL, NULL, NULL};

  struct SACHeader *ENZsacheader[3] = {NULL, NULL, NULL};

  nstime_t nstimetol;
  char stime[50];
  char etime[50];
  char tsname[50];
  char *cptr   = NULL;
  int snlength = 0;
  int idx;
  int retval = 0;

  if (!tid || !tseg || !plp)
    return -1;

  /* Cannot rotate data without coverage */
  if (!tseg->samprate || !tseg->numsamples)
    return 0;

  /* Identify trace ID's that match the requested set of components to rotate */
  snlength = snprintf (tsname, sizeof (tsname), "%s", tid->sid);

  for (idx = 0; idx < 3; idx++)
  {
    if (plp->rotateENZ[idx] == '\0')
      continue;

    tsname[snlength - 1] = plp->rotateENZ[idx];

    if ((id = mstl3_findID (mstl, tsname, 0, NULL)))
      ENZid[idx] = id;
  }

  /* Check that an appropriate channel set was found for requested 3-D or 2-D rotation */
  if (plp->rotations[1] && (!ENZid[0] || !ENZid[1] || !ENZid[2]))
  {
    if (verbose)
      fprintf (stderr, "Cannot find 3-D rotation channel set for %.*s, components %c, %c, %c\n",
               snlength - 1, tsname, plp->rotateENZ[0], plp->rotateENZ[1], plp->rotateENZ[2]);
    return 0;
  }
  if (plp->rotations[0] && (!ENZid[0] || !ENZid[1]))
  {
    if (verbose)
      fprintf (stderr, "Cannot find 2-D rotation channel set for %.*s, components %c, %c\n",
               snlength - 1, tsname, plp->rotateENZ[0], plp->rotateENZ[1]);
    return 0;
  }

  /* Determine time toleranace, as 1/2 sample period, in high precision time ticks */
  nstimetol = (tseg->samprate) ? (nstime_t) (NSTMODULUS / (2 * tseg->samprate)) : 0;

  /* Search for segments in each component group that match the target segment
   * in both time, sample count and sample rate */
  for (idx = 0; idx < 3; idx++)
  {
    if (ENZid[idx] == tid)
    {
      ENZseg[idx] = tseg;
    }
    else if (ENZid[idx])
    {
      for (seg = ENZid[idx]->first; seg; seg = seg->next)
      {
        if ((seg->starttime <= (tseg->starttime + nstimetol) && seg->starttime >= (tseg->starttime - nstimetol)) &&
            (seg->endtime <= (tseg->endtime + nstimetol) && seg->endtime >= (tseg->endtime - nstimetol)))
        {
          if (MS_ISRATETOLERABLE (seg->samprate, tseg->samprate))
            if (seg->numsamples == tseg->numsamples)
              ENZseg[idx] = seg;
        }
      }
    }
  }

  /* Check that an appropriate segment set was found for requested 3-D or 2-D rotation */
  if (plp->rotations[1] && (!ENZseg[0] || !ENZseg[1] || !ENZseg[2]))
  {
    if (verbose)
    {
      ms_nstime2timestr (tseg->starttime, stime, ISOMONTHDAY_DOY_Z, NANO_MICRO);
      ms_nstime2timestr (tseg->endtime, etime, ISOMONTHDAY_DOY_Z, NANO_MICRO);
      fprintf (stderr, "Cannot find 3-D rotation segment set matching %s [%s - %s]\n",
               tid->sid, stime, etime);
    }
    return 0;
  }
  if (plp->rotations[0] && (!ENZseg[0] || !ENZseg[1]))
  {
    if (verbose)
    {
      ms_nstime2timestr (tseg->starttime, stime, ISOMONTHDAY_DOY_Z, NANO_MICRO);
      ms_nstime2timestr (tseg->endtime, etime, ISOMONTHDAY_DOY_Z, NANO_MICRO);
      fprintf (stderr, "Cannot find 2-D rotation segment set matching %s [%s - %s]\n",
               tid->sid, stime, etime);
    }
    return 0;
  }

  if (plp->rotations[1])
    addToProcLog ("Rotating 3-D channel set %.*s, components %c,%c,%c (azimuth: %g, incidence: %g)",
                  snlength - 1, tsname, plp->rotateENZ[0], plp->rotateENZ[1], plp->rotateENZ[2],
                  plp->rotations[0], plp->rotations[1]);
  else
    addToProcLog ("Rotating 2-D channel set %.*s, components %c,%c (azimuth: %g)",
                  snlength - 1, tsname, plp->rotateENZ[0], plp->rotateENZ[1],
                  plp->rotations[0]);

  /* Convert samples to doubles */
  for (idx = 0; idx < 3; idx++)
  {
    if (ENZseg[idx])
    {
      if (ENZseg[idx]->sampletype != 'd')
      {
        if (convertSamples (ENZseg[idx], 'd'))
          return -1;
      }
    }
  }

  /* Perform 3-D or 2-D rotation */
  if (plp->rotations[1])
    rotate3 (ENZseg[2]->datasamples, ENZseg[1]->datasamples, ENZseg[0]->datasamples,
             tseg->numsamples,
             plp->rotations[0], plp->rotations[1],
             ENZseg[2]->datasamples, ENZseg[1]->datasamples, ENZseg[0]->datasamples);
  else
    rotate2 (ENZseg[1]->datasamples, ENZseg[0]->datasamples,
             tseg->numsamples,
             plp->rotations[0],
             ENZseg[1]->datasamples, ENZseg[0]->datasamples);

  /* Rename channel orientation codes and mark as rotated */
  /* idx: 0=E, 1=N, 2=Z */
  for (idx = 0; idx < 3; idx++)
  {
    if (!ENZid[idx])
      continue;

    if (plp->rotatedENZ[idx])
    {
      if ((cptr = strrchr (ENZid[idx]->sid, plp->rotateENZ[idx])))
        *cptr = plp->rotatedENZ[idx];
    }

    if (ENZseg[idx] && ENZseg[idx]->prvtptr)
    {
      struct segdetails *sd = (struct segdetails *)ENZseg[idx]->prvtptr;

      sd->rotated = 1;

      if (sd->sacheader)
      {
        /* Rename orientation code */
        if (plp->rotatedENZ[idx])
        {
          if ((cptr = strrchr (sd->sacheader->kcmpnm, plp->rotatedENZ[idx])))
            *cptr = plp->rotatedENZ[idx];
        }

        ENZsacheader[idx] = sd->sacheader;
      }
    }
  }

  /* Update SAC orientation values */
  if (plp->rotations[1] && ENZsacheader[0] && ENZsacheader[1] && ENZsacheader[2]) /* 3-D */
  {
    /* Rotate E/T component azimuth */
    if (ENZsacheader[0]->cmpaz != FUNDEF)
      ENZsacheader[0]->cmpaz += plp->rotations[0];

    /* Rotate N/Q component azimuth and reverse polarity */
    if (ENZsacheader[1]->cmpaz != FUNDEF)
      ENZsacheader[1]->cmpaz += plp->rotations[0] + 180;

    /* Rotate Z/L component azimuth, assume perpendicular to Q */
    if (ENZsacheader[2]->cmpaz != FUNDEF &&
        ENZsacheader[1]->cmpaz != FUNDEF)
      ENZsacheader[2]->cmpaz = ENZsacheader[1]->cmpaz - 180;

    for (idx = 0; idx < 3; idx++)
    {
      while (ENZsacheader[idx]->cmpaz > 360.0)
        ENZsacheader[idx]->cmpaz -= 360.0;

      while (ENZsacheader[idx]->cmpaz < 0.0)
        ENZsacheader[idx]->cmpaz += 360.0;
    }

    /* Update incidence angle for L and Q components */
    if (ENZsacheader[1]->cmpinc != FUNDEF)
      ENZsacheader[1]->cmpinc -= plp->rotations[1];

    if (ENZsacheader[2]->cmpinc != FUNDEF)
      ENZsacheader[2]->cmpinc += plp->rotations[1];

    for (idx = 1; idx < 3; idx++)
    {
      while (ENZsacheader[idx]->cmpinc > 180.0)
        ENZsacheader[idx]->cmpinc -= 180.0;

      while (ENZsacheader[idx]->cmpinc < 0.0)
        ENZsacheader[idx]->cmpinc += 180.0;
    }
  }
  else if (ENZsacheader[0] && ENZsacheader[1]) /* 2-D, horizontals only */
  {
    for (idx = 0; idx < 2; idx++)
    {
      if (ENZsacheader[idx]->cmpaz != FUNDEF)
        ENZsacheader[idx]->cmpaz += plp->rotations[0];

      while (ENZsacheader[idx]->cmpaz > 360.0)
        ENZsacheader[idx]->cmpaz -= 360.0;

      while (ENZsacheader[idx]->cmpaz < 0.0)
        ENZsacheader[idx]->cmpaz += 360.0;
    }
  }

  if (retval < 0)
  {
    fprintf (stderr, "procRotate(): Error rotating time series\n");
    return -1;
  }

  return 0;
} /* End of procRotate() */

/***************************************************************************
 * readMSEED:
 *
 * Read file containing miniSEED and add data to the supplied MS3TraceList.
 *
 * Returns the number of samples read on success and -1 on error.
 ***************************************************************************/
static int64_t
readMSEED (char *mseedfile, MS3TraceList *mstl)
{
  MS3Record *msr     = NULL;
  MS3TraceSeg *seg   = NULL;
  uint32_t flags     = 0;
  int64_t totalsamps = 0;
  int retcode        = MS_NOERROR;
  char stime[100];

  if (!mseedfile)
  {
    fprintf (stderr, "readMSEED(): No input file specified\n");
    return -1;
  }

  if (!mstl)
  {
    fprintf (stderr, "readMSEED(): No MS3TraceList specified\n");
    return -1;
  }

  /* Set flags to validate CRCs, parse byte ranges from file names and unpack data */
  flags |= MSF_VALIDATECRC;
  flags |= MSF_PNAMERANGE;
  flags |= MSF_UNPACKDATA;

  /* Loop over the input file reading records */
  while ((retcode = ms3_readmsr (&msr, mseedfile, flags, verbose - 3)) == MS_NOERROR)
  {
    /* Skip data records that do not contain time series data */
    if (msr->numsamples == 0 || (msr->sampletype != 'i' && msr->sampletype != 'f' && msr->sampletype != 'd'))
    {
      if (verbose >= 3)
      {
        ms_nstime2timestr (msr->starttime, stime, ISOMONTHDAY_DOY_Z, NANO_MICRO);
        fprintf (stderr, "Skipping (no time series data) %s, %s\n", msr->sid, stime);
      }
      continue;
    }

    /* Check if record matches start/end time criteria */
    if (starttime != NSTUNSET && (msr->starttime < starttime))
    {
      if (verbose >= 3)
      {
        ms_nstime2timestr (msr->starttime, stime, ISOMONTHDAY_DOY_Z, NANO_MICRO);
        fprintf (stderr, "Skipping (start time) %s, %s\n", msr->sid, stime);
      }
      continue;
    }

    if (endtime != NSTUNSET && (msr3_endtime (msr) > endtime))
    {
      if (verbose >= 3)
      {
        ms_nstime2timestr (msr->starttime, stime, ISOMONTHDAY_DOY_Z, NANO_MICRO);
        fprintf (stderr, "Skipping (end time) %s, %s\n", msr->sid, stime);
      }
      continue;
    }

    totalsamps += msr->samplecnt;

    if (verbose >= 3)
      msr3_print (msr, verbose - 3);

    /* Add new data to MS3TraceList, merge with other segments */
    if (!(seg = mstl3_addmsr (mstl, msr, 0, 1, flags, NULL)))
    {
      fprintf (stderr, "[%s] Error adding samples to MS3TraceList\n", mseedfile);
    }

    /* Add additional segment details structure */
    if (seg && !seg->prvtptr)
    {
      if (!(seg->prvtptr = malloc (sizeof (struct segdetails))))
      {
        fprintf (stderr, "Error allocating memory for additiona segment details\n");
      }

      ((struct segdetails *)seg->prvtptr)->procerror = 0;
      ((struct segdetails *)seg->prvtptr)->rotated   = 0;
      ((struct segdetails *)seg->prvtptr)->sacheader = NULL;
    }
  }

  /* Print error if not EOF and not counting down records */
  if (retcode != MS_ENDOFFILE)
  {
    fprintf (stderr, "Error reading %s: %s\n", mseedfile, ms_errorstr (retcode));
    totalsamps = -1;
  }

  /* Make sure everything is cleaned up */
  ms3_readmsr (&msr, NULL, flags, verbose - 3);

  return totalsamps;
} /* End of readMSEED() */

/***************************************************************************
 * readSAC:
 *
 * Read a SAC file and add data samples to a MS3TraceGroup.  As the SAC
 * data is read in a MS3Record struct is used as a holder for the input
 * information.
 *
 * The SAC header contents is stored at the private pointer of the
 * MS3TraceSeg structure (MS3TraceSeg->prvtptr->sacheader).
 *
 * Returns the number of samples read on success and -1 on failure
 ***************************************************************************/
static int64_t
readSAC (char *sacfile, MS3TraceList *mstl)
{
  FILE *ifp     = 0;
  MS3Record *msr = 0;
  MS3TraceSeg *seg;
  uint32_t flags = 0;

  struct SACHeader sh;
  float *fdata = 0;
  int datacnt;

  char net[11] = {0};
  char sta[11] = {0};
  char loc[11] = {0};
  char chan[11] = {0};

  /* Open input file */
  if ((ifp = fopen (sacfile, "rb")) == NULL)
  {
    fprintf (stderr, "Cannot open input file: %s (%s)\n",
             sacfile, strerror (errno));
    return -1;
  }

  /* Parse input SAC file into a header structure and data buffer */
  if ((datacnt = parseSAC (ifp, &sh, &fdata, sacinformat, verbose, sacfile)) < 0)
  {
    fprintf (stderr, "Error parsing %s\n", sacfile);
    return -1;
  }

  if (!(msr = msr3_init (msr)))
  {
    fprintf (stderr, "Cannot initialize MS3Record strcture\n");
    return -1;
  }

  /* Populate MS3Record structure with header details */
  if (strncmp (SUNDEF, sh.knetwk, 8))
    ms_strncpclean (net, sh.knetwk, 8);
  if (strncmp (SUNDEF, sh.kstnm, 8))
    ms_strncpclean (sta, sh.kstnm, 8);
  if (strncmp (SUNDEF, sh.khole, 8))
    ms_strncpclean (loc, sh.khole, 8);
  if (strncmp (SUNDEF, sh.kcmpnm, 8))
    ms_strncpclean (chan, sh.kcmpnm, 8);

  if (sacnet)
    ms_strncpclean (net, sacnet, sizeof (net));

  if (sacloc)
    ms_strncpclean (loc, sacloc, sizeof (loc));

  ms_sid2nslc (msr->sid, net, sta, loc, chan);

  msr->starttime = ms_time2nstime (sh.nzyear, sh.nzjday, sh.nzhour, sh.nzmin, sh.nzsec, sh.nzmsec * 1000000);

  /* Adjust for Begin ('B' SAC variable) time offset */
  msr->starttime += (double)sh.b * NSTMODULUS;

  /* Calculate sample rate from interval(period) rounding to nearest 0.000001 Hz */
  msr->samprate = (double)((int)((1 / sh.delta) * 100000 + 0.5)) / 100000;

  msr->samplecnt = msr->numsamples = datacnt;

  msr->sampletype  = 'f';
  msr->datasamples = fdata;

  if (verbose >= 1)
  {
    fprintf (stderr, "[%s] %" PRId64 " samps @ %.6f Hz for %s\n",
             sacfile, msr->numsamples, msr->samprate, msr->sid);
  }

  /* Add new data to MS3TraceList, do not merge with other segments */
  if (!(seg = mstl3_addmsr (mstl, msr, 0, 0, flags, NULL)))
  {
    fprintf (stderr, "[%s] Error adding samples to MS3TraceList\n", sacfile);
  }

  /* Add additional segment details structure */
  if (!seg->prvtptr)
  {
    if (!(seg->prvtptr = malloc (sizeof (struct segdetails))))
    {
      fprintf (stderr, "Error allocating memory for additiona segment details\n");
    }

    ((struct segdetails *)seg->prvtptr)->procerror = 0;
    ((struct segdetails *)seg->prvtptr)->rotated   = 0;
    ((struct segdetails *)seg->prvtptr)->sacheader = NULL;
  }

  /* Store SAC header structure in segment details */
  if (!((struct segdetails *)seg->prvtptr)->sacheader)
  {
    if (!(((struct segdetails *)seg->prvtptr)->sacheader = malloc (sizeof (struct SACHeader))))
    {
      fprintf (stderr, "Error allocating memory for SACHheader\n");
    }
    else
    {
      memcpy (((struct segdetails *)seg->prvtptr)->sacheader, &sh, sizeof (struct SACHeader));
    }
  }
  else if (verbose >= 1)
  {
    fprintf (stderr, "%s: SAC header contents not stored for merged segment", msr->sid);
  }

  fclose (ifp);

  if (fdata)
    free (fdata);

  msr->datasamples = 0;

  if (msr)
    msr3_free (&msr);

  return datacnt;
} /* End of readSAC() */

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
  if (!ifp || !sh || !data)
    return -1;

  /* Read the first 4 characters */
  if (fread (&fourc, 4, 1, ifp) < 1)
    return -1;

  /* Determine if the file is ALPHA or binary SAC,
   * if the first 4 characters are spaces assume ALPHA SAC */
  if (format == 0)
  {
    if (fourc[0] == ' ' && fourc[1] == ' ' && fourc[2] == ' ' && fourc[3] == ' ')
      format = 1;
    else
      format = 2; /* Byte order detection will occur below */
  }

  /* Rewind the file position pointer to the beginning */
  rewind (ifp);

  /* Read the header */
  if (format == 1) /* Process SAC ALPHA header */
  {
    if ((rv = readAlphaHeaderSAC (ifp, sh)))
    {
      fprintf (stderr, "[%s] Error parsing SAC ALPHA header at line %d\n",
               sacfile, rv);
      return -1;
    }
  }
  else if (format >= 2 && format <= 4) /* Process SAC binary header */
  {
    if (readBinaryHeaderSAC (ifp, sh, &format, &swapflag, verbose, sacfile))
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
  if (sh->nzyear < 1900 || sh->nzyear > 3000 ||
      sh->nzjday < 1 || sh->nzjday > 366 ||
      sh->nzhour < 0 || sh->nzhour > 23 ||
      sh->nzmin < 0 || sh->nzmin > 59 ||
      sh->nzsec < 0 || sh->nzsec > 60 ||
      sh->nzmsec < 0 || sh->nzmsec > 999999)
  {
    fprintf (stderr, "[%s] Unrecognized format (not SAC?)\n", sacfile);
    return -1;
  }

  if (verbose)
  {
    if (format == 1)
      fprintf (stderr, "[%s] Reading SAC ALPHA format\n", sacfile);
    if (format == 3)
      fprintf (stderr, "[%s] Reading SAC binary format (little-endian)\n", sacfile);
    if (format == 4)
      fprintf (stderr, "[%s] Reading SAC binary format (big-endian)\n", sacfile);
  }

  if (verbose > 2)
    fprintf (stderr, "[%s] SAC header version number: %d\n", sacfile, sh->nvhdr);

  if (sh->nvhdr != 6)
    fprintf (stderr, "[%s] WARNING SAC header version (%d) not expected value of 6\n",
             sacfile, sh->nvhdr);

  if (sh->npts <= 0)
  {
    fprintf (stderr, "[%s] No data, number of samples: %d\n", sacfile, sh->npts);
    return -1;
  }

  if (sh->iftype != ITIME)
  {
    fprintf (stderr, "[%s] Data is not time series (IFTYPE=%d), cannot convert other types\n",
             sacfile, sh->iftype);
    return -1;
  }

  if (!sh->leven)
  {
    fprintf (stderr, "[%s] Data is not evenly spaced (LEVEN not true), cannot convert\n", sacfile);
    return -1;
  }

  /* Allocate space for data samples */
  if (!(*data = (float *)calloc (1, sizeof (float) * sh->npts)))
  {
    fprintf (stderr, "Error allocating memory for data samples\n");
    return -1;
  }

  /* Read the data samples */
  if (format == 1) /* Process SAC ALPHA data */
  {
    if ((rv = readAlphaDataSAC (ifp, *data, sh->npts)))
    {
      fprintf (stderr, "[%s] Error parsing SAC ALPHA data at line %d\n",
               sacfile, rv);
      return -1;
    }
  }
  else if (format >= 2 && format <= 4) /* Process SAC binary data */
  {
    if (readBinaryDataSAC (ifp, *data, sh->npts, swapflag, verbose, sacfile))
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
} /* End of parseSAC() */

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
  if (fread (sh, sizeof (struct SACHeader), 1, ifp) != 1)
  {
    fprintf (stderr, "[%s] Could not read SAC header from file\n", sacfile);

    if (ferror (ifp))
      fprintf (stderr, "[%s] Error reading from file\n", sacfile);

    return -1;
  }

  /* Determine if host is big-endian */
  bigendianhost = ms_bigendianhost ();

  *swapflag = 0;

  /* Test byte order using the header version if unknown */
  /* Also set the swapflag appropriately */
  if (*format == 2)
  {
    memcpy (&hdrver, &sh->nvhdr, 4);
    if (hdrver < 1 || hdrver > 10)
    {
      ms_gswap4 (&hdrver);
      if (hdrver < 1 || hdrver > 10)
      {
        fprintf (stderr, "[%s] Cannot determine byte order (not SAC?)\n", sacfile);
        return -1;
      }

      *format   = (bigendianhost) ? 3 : 4;
      *swapflag = 1;
    }
    else
    {
      *format = (bigendianhost) ? 4 : 3;
    }
  }
  else if (*format == 3 && bigendianhost)
    *swapflag = 1;
  else if (*format == 4 && !bigendianhost)
    *swapflag = 1;

  if (verbose > 1)
  {
    if (*swapflag)
      fprintf (stderr, "[%s] Byte swapping required\n", sacfile);
    else
      fprintf (stderr, "[%s] Byte swapping NOT required\n", sacfile);
  }

  /* Byte swap all values in header */
  if (*swapflag)
    swapSACHeader (sh);

  return 0;
} /* End of readBinaryHeaderSAC() */

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
  if ((samplesread = fread (data, sizeof (float), datacnt, ifp)) != datacnt)
  {
    fprintf (stderr, "[%s] Only read %d of %d expected data samples\n",
             sacfile, samplesread, datacnt);
    return -1;
  }

  /* Swap data samples */
  if (swapflag)
  {
    for (dataidx = 0; dataidx < datacnt; dataidx++)
    {
      ms_gswap4 (data + dataidx);
    }
  }

  return 0;
} /* End of readBinaryDataSAC() */

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
  int linecnt = 1; /* The header starts at line 1 */
  int lineidx;
  int count;
  int hvidx = 0;
  char *cp;

  if (!ifp || !sh)
    return -1;

  /* The first 14 lines x 5 values are floats */
  for (lineidx = 0; lineidx < 14; lineidx++)
  {
    if (!fgets (line, sizeof (line), ifp))
      return linecnt;

    count = sscanf (line, " %f %f %f %f %f ", (float *)sh + hvidx,
                    (float *)sh + hvidx + 1, (float *)sh + hvidx + 2,
                    (float *)sh + hvidx + 3, (float *)sh + hvidx + 4);

    if (count != 5)
      return linecnt;

    hvidx += 5;
    linecnt++;
  }

  /* The next 8 lines x 5 values are integers */
  for (lineidx = 0; lineidx < 8; lineidx++)
  {
    if (!fgets (line, sizeof (line), ifp))
      return linecnt;

    count = sscanf (line, " %d %d %d %d %d ", (int32_t *)sh + hvidx,
                    (int32_t *)sh + hvidx + 1, (int32_t *)sh + hvidx + 2,
                    (int32_t *)sh + hvidx + 3, (int32_t *)sh + hvidx + 4);

    if (count != 5)
      return linecnt;

    hvidx += 5;
    linecnt++;
  }

  /* Set pointer to start of string variables */
  cp = (char *)sh + (hvidx * 4);

  /* The next 8 lines each contain 24 bytes of string data */
  for (lineidx = 0; lineidx < 8; lineidx++)
  {
    memset (line, 0, sizeof (line));
    if (!fgets (line, sizeof (line), ifp))
      return linecnt;

    memcpy (cp, line, 24);
    cp += 24;

    linecnt++;
  }

  /* Make sure each of the 23 string variables are left justified */
  cp = (char *)sh + (hvidx * 4);
  for (count = 0; count < 24; count++)
  {
    int ridx, widx, width;
    char *fcp;

    /* Each string variable is 8 characters with one exception */
    if (count != 1)
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
    while (*(fcp + ridx) == ' ')
      ridx++;

    /* Remove any leading spaces */
    if (ridx > 0)
    {
      for (widx = 0; widx < width; widx++, ridx++)
      {
        if (ridx < width)
          *(fcp + widx) = *(fcp + ridx);
        else
          *(fcp + widx) = ' ';
      }
    }
  }

  return 0;
} /* End of readAlphaHeaderSAC() */

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
  int linecnt     = 31; /* Data samples start on line 31 */
  int samplesread = 0;
  int count;
  int dataidx = 0;

  if (!ifp || !data || !datacnt)
    return -1;

  /* Each data line should contain 5 floats unless the last */
  for (;;)
  {
    if (!fgets (line, sizeof (line), ifp))
      return linecnt;

    count = sscanf (line, " %f %f %f %f %f ", (float *)data + dataidx,
                    (float *)data + dataidx + 1, (float *)data + dataidx + 2,
                    (float *)data + dataidx + 3, (float *)data + dataidx + 4);

    samplesread += count;

    if (samplesread >= datacnt)
      break;
    else if (count != 5)
      return linecnt;

    dataidx += 5;
    linecnt++;
  }

  return 0;
} /* End of readAlphaDataSAC() */

/***************************************************************************
 * writeMSEED:
 *
 * Write data buffer to output file as miniSEED, converting sample
 * type as required by the encoding format.
 *
 * Returns the number of samples written, -1 on error.
 ***************************************************************************/
static int
writeMSEED (MS3TraceID *id, MS3TraceSeg *seg, char *outputfile)
{
  MS3Record *msr = NULL;
  uint32_t flags = 0;
  char outfile[1024];
  FILE *ofp             = 0;
  int64_t packedsamples = 0;
  int packedrecords     = 0;

  char net[11] = {0};
  char sta[11] = {0};
  char loc[11] = {0};
  char chan[11] = {0};

  uint16_t year;
  uint16_t yday;
  uint8_t hour;
  uint8_t min;
  uint8_t sec;
  uint32_t nsec;

  /* Determine file open mode:
   * If bytes have already been written: append
   * If nothing has been written: overwrite */
  char *mode = (outputbytes > 0) ? "ab" : "wb";

  if (!id || !seg)
    return -1;

  if (seg->numsamples <= 0 || seg->samprate == 0.0)
    return 0;

  /* Convert data buffer and type depending on encoding format */
  if (packencoding == 3 || packencoding == 10 || packencoding == 11) /* 32-bit integers */
  {
    if (convertSamples (seg, 'i'))
      return -1;

    if (seg->sampletype != 'i')
    {
      fprintf (stderr, "writeMSEED(): Unsupported sample type: %c\n", seg->sampletype);
      return -1;
    }
  }
  else if (packencoding == 4) /* 32-bit floats */
  {
    if (convertSamples (seg, 'f'))
      return -1;

    if (seg->sampletype != 'f')
    {
      fprintf (stderr, "writeMSEED(): Unsupported sample type: %c\n", seg->sampletype);
      return -1;
    }
  }
  else if (packencoding == 5) /* 64-bit floats */
  {
    if (convertSamples (seg, 'd'))
      return -1;

    if (seg->sampletype != 'd')
    {
      fprintf (stderr, "writeMSEED(): Unsupported sample type: %c\n", seg->sampletype);
      return -1;
    }
  }
  else
  {
    fprintf (stderr, "Unrecognized encoding format: %d\n", packencoding);
    return -1;
  }

  /* Populate MS3Record used for packing */
  if (!(msr = msr3_init (msr)))
  {
    ms_log (2, "Could not allocate MS3Record, out of memory?\n");
    return -1;
  }

  msr->reclen   = packreclen;
  msr->encoding = packencoding;

  ms_sid2nslc (id->sid, net, sta, loc, chan);

  if (channel)
  {
    strncpy (chan, channel, sizeof (chan) - 1);
  }

  ms_nslc2sid (msr->sid, sizeof (msr->sid), 0, net, sta, loc, chan);

  msr->pubversion  = id->pubversion;
  msr->starttime   = seg->starttime;
  msr->samprate    = seg->samprate;
  msr->samplecnt   = seg->samplecnt;
  msr->datasamples = seg->datasamples;
  msr->numsamples  = seg->numsamples;
  msr->sampletype  = seg->sampletype;

  if (verbose)
    fprintf (stderr, "Writing miniSEED for %s\n", msr->sid);

  if (!outputfile)
  {
    if (ms_nstime2time (msr->starttime, &year, &yday, &hour, &min, &sec, &nsec))
    {
      fprintf (stderr, "Cannot decompose nstime into date-time components\n");
      return -1;
    }

    /* Create output file name: Net.Sta.Loc.Chan.Qual.Year,Day,Hour:Min:Sec.mseed */
    snprintf (outfile, sizeof (outfile), "%s%s%s.%s.%s.%s.%d.%04d,%03d,%02d:%02d:%02d.mseed",
              (outputdir) ? outputdir : "", (outputdir) ? "/" : "",
              net, sta, loc, chan, msr->pubversion,
              year, yday, hour, min, sec);
  }
  else
  {
    strncpy (outfile, outputfile, sizeof (outfile));
  }

  /* Open output file */
  if (strcmp (outfile, "-") == 0)
  {
    ofp = stdout;
  }
  else if ((ofp = fopen (outfile, mode)) == NULL)
  {
    fprintf (stderr, "Cannot open output file: %s (%s)\n",
             outfile, strerror (errno));
    return -1;
  }

  /* Set data flush flag */
  flags |= MSF_FLUSHDATA;

  /* Set miniSEED v2 if requested */
  if (msformat == 2)
    flags |= MSF_PACKVER2;

  /* Pack data into records */
  packedrecords = msr3_pack (msr, recordHandler, ofp, &packedsamples, flags, verbose - 2);

  if (verbose)
    fprintf (stderr, "Packed %" PRId64 " samples into %d records\n",
             packedsamples, packedrecords);

  /* Make sure everything is cleaned up */
  msr->datasamples = NULL;
  msr->numsamples  = 0;
  if (msr)
    msr3_free (&msr);

  if (ofp && ofp != stdout)
    fclose (ofp);

  return packedsamples;
} /* End of writeMSEED() */

/***************************************************************************
 * writeSAC:
 *
 * Write data buffer to output file as SAC.
 *
 * Returns the number of samples written or -1 on error.
 ***************************************************************************/
static int
writeSAC (MS3TraceID *id, MS3TraceSeg *seg, int format, char *outputfile)
{
  struct SACHeader sh = NullSACHeader;
  char outfile[1024];

  float *fdata = 0;

  char net[11] = {0};
  char sta[11] = {0};
  char loc[11] = {0};
  char chan[11] = {0};

  uint16_t year;
  uint16_t yday;
  uint8_t hour;
  uint8_t min;
  uint8_t sec;
  uint32_t nsec;

  nstime_t submsec;
  int idx;

  if (!id || !seg)
    return -1;

  if (seg->numsamples <= 0 || seg->samprate == 0.0)
    return 0;

  if (verbose)
    fprintf (stderr, "Writing SAC for %s\n", id->sid);

  /* If an original input SAC header is available use it as a base to update */
  if (seg->prvtptr && ((struct segdetails *)seg->prvtptr)->sacheader)
  {
    struct SACHeader *osh = (struct SACHeader *)((struct segdetails *)seg->prvtptr)->sacheader;
    nstime_t ostarttime;

    memcpy (&sh, osh, sizeof (struct SACHeader));

    /* Calculate original input start time */
    ostarttime = ms_time2nstime (osh->nzyear, osh->nzjday, osh->nzhour, osh->nzmin,
                                 osh->nzsec, osh->nzmsec * 1000000);

    /* Adjust for Begin ('B' SAC variable) time offset */
    ostarttime += (double)osh->b * NSTMODULUS;

    /* If data start time has changed adjust Begin and End header variables */
    if (ostarttime != seg->starttime)
    {
      float shift = (float)(seg->starttime - ostarttime) / NSTMODULUS;

      if (verbose >= 1)
        fprintf (stderr, "Updating data begin and end time variables (%g seconds)\n", shift);

      sh.b += shift;
      sh.e = sh.b + (seg->numsamples - 1) * (1 / seg->samprate);
    }
  }
  else /* Otherwise set zero time and time of first and last sample */
  {
    /* Set time zero to start time */
    ms_nstime2time (seg->starttime, &year, &yday, &hour, &min, &sec, &nsec);
    sh.nzyear = year;
    sh.nzjday = yday;
    sh.nzhour = hour;
    sh.nzmin  = min;
    sh.nzsec  = sec;
    sh.nzmsec = nsec / 1000000;

    /* Determine any sub-millisecond portion of the start time in HP time */
    submsec = (seg->starttime -
               ms_time2nstime (sh.nzyear, sh.nzjday, sh.nzhour,
                               sh.nzmin, sh.nzsec, sh.nzmsec * 1000000));

    /* Set begin and end offsets from reference time for first and last sample,
       * any sub-millisecond start time is stored in these offsets. */
    sh.b = ((float)submsec / NSTMODULUS);
    sh.e = sh.b + (seg->numsamples - 1) * (1 / seg->samprate);
  }

  /* Set time series source parameters */
  ms_sid2nslc (id->sid, net, sta, loc, chan);

  /* Channel override */
  if (channel)
    strncpy (chan, channel, sizeof (chan));

  if (*net != '\0')
    strncpy (sh.knetwk, net, sizeof (sh.knetwk));
  if (*sta != '\0')
    strncpy (sh.kstnm, sta, sizeof (sh.kstnm));
  if (*loc != '\0')
    strncpy (sh.khole, loc, sizeof (sh.khole));
  if (*chan != '\0')
    strncpy (sh.kcmpnm, chan, sizeof (sh.kcmpnm));

  /* Set misc. header variables */
  sh.nvhdr  = 6;     /* Header version = 6 */
  sh.leven  = 1;     /* Evenly spaced data */
  sh.iftype = ITIME; /* Data is time series */

  /* Insert metadata if present */
  if (metadata)
    insertSACMetaData (&sh, seg->starttime);

  /* Set sampling interval (seconds), sample count */
  sh.delta = 1.0 / seg->samprate;
  sh.npts  = seg->numsamples;

  /* Original input SAC variables that may no longer be valid:
   * scale : amplitude scale factor
   * idep  : type of amplitude
   */

  /* Convert samples to floats */
  if (seg->sampletype != 'f')
  {
    if (convertSamples (seg, 'f'))
      return -1;
  }

  if (seg->sampletype != 'f')
  {
    fprintf (stderr, "writeSAC: Unusable sample type: %c\n", seg->sampletype);
    return -1;
  }

  fdata = (float *)seg->datasamples;

  /* Determine minimum and maximum sample values */
  sh.depmin = *fdata;
  sh.depmax = *fdata;
  for (idx = 1; idx < seg->numsamples; idx++)
  {
    if (fdata[idx] < sh.depmin)
      sh.depmin = fdata[idx];
    if (fdata[idx] > sh.depmax)
      sh.depmax = fdata[idx];
  }

  if (format >= 2 && format <= 4)
  {
    if (!outputfile)
    {
      /* Create output file name: Net.Sta.Loc.Chan.Qual.Year.Day.Hour.Min.Sec.SAC */
      snprintf (outfile, sizeof (outfile), "%s%s%s.%s.%s.%s.%d.%04d,%03d,%02d:%02d:%02d.SAC",
                (outputdir) ? outputdir : "", (outputdir) ? "/" : "",
                net, sta, loc, chan,
                id->pubversion, sh.nzyear, sh.nzjday, sh.nzhour,
                sh.nzmin, sh.nzsec);
    }
    else
    {
      if (verbose && outputbytes)
        fprintf (stderr, "Warning, output file will be overwritten: %s\n", outputfile);

      strncpy (outfile, outputfile, sizeof (outfile));
    }

    /* Byte swap the data header and data if needed */
    if ((format == 3 && ms_bigendianhost ()) ||
        (format == 4 && !ms_bigendianhost ()))
    {
      if (verbose)
        fprintf (stderr, "Byte swapping SAC header and data\n");

      swapSACHeader (&sh);

      for (idx = 0; idx < seg->numsamples; idx++)
      {
        ms_gswap4 (fdata + idx);
      }
    }

    if (verbose > 1)
      fprintf (stderr, "Writing binary SAC file '%s'\n", outfile);

    /* Reset output bytes counter */
    outputbytes = 0;

    if (writeBinarySAC (&sh, fdata, seg->numsamples, outfile))
      return -1;
  }
  else if (format == 1)
  {
    if (!outputfile)
    {
      /* Create output file name: Net.Sta.Loc.Chan.Qual.Year.Day.Hour.Min.Sec.SACA */
      snprintf (outfile, sizeof (outfile), "%s%s%s.%s.%s.%s.%d.%04d,%03d,%02d:%02d:%02d.SACA",
                (outputdir) ? outputdir : "", (outputdir) ? "/" : "",
                net, sta, loc, chan,
                id->pubversion, sh.nzyear, sh.nzjday, sh.nzhour,
                sh.nzmin, sh.nzsec);
    }
    else
    {
      if (verbose && outputbytes)
        fprintf (stderr, "Warning, output file will be overwritten: %s\n", outputfile);

      strncpy (outfile, outputfile, sizeof (outfile));
    }

    if (verbose > 1)
      fprintf (stderr, "Writing alphanumeric SAC file: '%s'\n", outfile);

    /* Reset output bytes counter */
    outputbytes = 0;

    if (writeAlphaSAC (&sh, fdata, seg->numsamples, outfile))
      return -1;
  }
  else
  {
    fprintf (stderr, "Error, unrecognized SAC format: '%d'\n", format);
  }

  if (verbose)
    fprintf (stderr, "Wrote %lld samples to %s (%d bytes)\n",
             (long long int)seg->numsamples, outfile, outputbytes);

  return seg->numsamples;
} /* End of writeSAC() */

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
  if ((ofp = fopen (outfile, "wb")) == NULL)
  {
    fprintf (stderr, "Cannot open output file: %s (%s)\n",
             outfile, strerror (errno));
    return -1;
  }

  /* Write SAC header to output file */
  if (fwrite (sh, sizeof (struct SACHeader), 1, ofp) != 1)
  {
    fprintf (stderr, "Error writing SAC header to output file\n");
    return -1;
  }
  outputbytes += sizeof (struct SACHeader);

  /* Write float data to output file */
  if (fwrite (fdata, sizeof (float), npts, ofp) != npts)
  {
    fprintf (stderr, "Error writing SAC data to output file\n");
    return -1;
  }
  outputbytes += sizeof (float) * npts;

  fclose (ofp);

  return 0;
} /* End of writeBinarySAC() */

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
  float *fhp   = (float *)sh;
  int32_t *ihp = (int32_t *)sh + (NUMFLOATHDR);
  char *shp    = (char *)sh + (NUMFLOATHDR * 4 + NUMINTHDR * 4);

  /* Open output file */
  if ((ofp = fopen (outfile, "wb")) == NULL)
  {
    fprintf (stderr, "Cannot open output file: %s (%s)\n",
             outfile, strerror (errno));
    return -1;
  }

  /* Write SAC header float variables to output file, 5 variables per line */
  for (idx = 0; idx < NUMFLOATHDR; idx += 5)
  {
    for (fidx = idx; fidx < (idx + 5) && fidx < NUMFLOATHDR; fidx++)
      fprintf (ofp, "%#15.7g", *(fhp + fidx));

    fprintf (ofp, "\n");
  }
  outputbytes += 1064;

  /* Write SAC header integer variables to output file, 5 variables per line */
  for (idx = 0; idx < NUMINTHDR; idx += 5)
  {
    for (fidx = idx; fidx < (idx + 5) && fidx < NUMINTHDR; fidx++)
      fprintf (ofp, "%10d", *(ihp + fidx));

    fprintf (ofp, "\n");
  }
  outputbytes += 408;

  /* Write SAC header string variables to output file, 3 variables per line */
  for (idx = 0; idx < NUMSTRHDR; idx += 3)
  {
    if (idx == 0)
      fprintf (ofp, "%-8.8s%-16.16s", shp, shp + 8);
    else
      fprintf (ofp, "%-8.8s%-8.8s%-8.8s", shp + (idx * 8), shp + ((idx + 1) * 8), shp + ((idx + 2) * 8));

    fprintf (ofp, "\n");
  }
  outputbytes += 200;

  /* Write float data to output file, 5 values per line */
  for (idx = 0; idx < npts; idx += 5)
  {
    for (fidx = idx; fidx < (idx + 5) && fidx < npts && fidx >= 0; fidx++)
    {
      fprintf (ofp, "%#15.7g", *(fdata + fidx));
      outputbytes += 15;
    }

    fprintf (ofp, "\n");
    outputbytes++;
  }

  fclose (ofp);

  return 0;
} /* End of writeAlphaSAC() */

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
insertSACMetaData (struct SACHeader *sh, nstime_t sacstarttime)
{
  struct metalist *mlp = metadata;
  struct metanode *mn  = NULL;
  nstime_t sacendtime;
  char *endptr;
  char sacnetwork[9];
  char sacstation[9];
  char saclocation[9];
  char sacchannel[9];
  int retval = 1;

  if (!mlp || !sh)
    return -1;

  /* Determine source name parameters for comparison, as a special case if the
   * location code is not set it will match '--' */
  if (strncmp (sh->knetwk, SUNDEF, 8))
    ms_strncpclean (sacnetwork, sh->knetwk, 8);
  else
    sacnetwork[0] = '\0';
  if (strncmp (sh->kstnm, SUNDEF, 8))
    ms_strncpclean (sacstation, sh->kstnm, 8);
  else
    sacstation[0] = '\0';
  if (strncmp (sh->khole, SUNDEF, 8))
    ms_strncpclean (saclocation, sh->khole, 8);
  else
  {
    saclocation[0] = '-';
    saclocation[1] = '-';
    saclocation[2] = '\0';
  }
  if (strncmp (sh->kcmpnm, SUNDEF, 8))
    ms_strncpclean (sacchannel, sh->kcmpnm, 8);
  else
    sacchannel[0] = '\0';

  /* Calculate end time of SAC data */
  sacendtime = sacstarttime + (((sh->npts - 1) * sh->delta) * NSTMODULUS);

  while (mlp)
  {
    mn = (struct metanode *)mlp->data;

    /* Sanity check that source name fields are present */
    if (!mn->metafields[0] || !mn->metafields[1] ||
        !mn->metafields[2] || !mn->metafields[3])
    {
      fprintf (stderr, "insertmetadata(): error, source name fields not all present\n");
    }
    /* Test if network, station, location and channel; also handle simple wildcards */
    else if ((!strncmp (sacnetwork, mn->metafields[0], 8) || (*(mn->metafields[0]) == '*')) &&
             (!strncmp (sacstation, mn->metafields[1], 8) || (*(mn->metafields[1]) == '*')) &&
             (!strncmp (saclocation, mn->metafields[2], 8) || (*(mn->metafields[2]) == '*')) &&
             (!strncmp (sacchannel, mn->metafields[3], 8) || (*(mn->metafields[3]) == '*')))
    {
      /* Check time window match */
      if (mn->starttime != NSTUNSET || mn->endtime != NSTUNSET)
      {
        /* Check for overlap with metadata window */
        if (mn->starttime != NSTUNSET && mn->endtime != NSTUNSET)
        {
          if (!(sacendtime >= mn->starttime && sacstarttime <= mn->endtime))
          {
            mlp = mlp->next;
            continue;
          }
        }
        /* Check if data after start time */
        else if (mn->starttime != NSTUNSET)
        {
          if (sacendtime < mn->starttime)
          {
            mlp = mlp->next;
            continue;
          }
        }
        /* Check if data before end time */
        else if (mn->endtime != NSTUNSET)
        {
          if (sacstarttime > mn->endtime)
          {
            mlp = mlp->next;
            continue;
          }
        }
      }

      if (verbose)
        fprintf (stderr, "Inserting metadata for N: '%s', S: '%s', L: '%s', C: '%s' (%s - %s)\n",
                 sacnetwork, sacstation, saclocation, sacchannel,
                 (mn->metafields[15]) ? mn->metafields[15] : "NONE",
                 (mn->metafields[16]) ? mn->metafields[16] : "NONE");

      /* Insert metadata into SAC header */
      if (mn->metafields[4])
        sh->stla = (float)strtod (mn->metafields[4], &endptr);
      if (mn->metafields[5])
        sh->stlo = (float)strtod (mn->metafields[5], &endptr);
      if (mn->metafields[6])
        sh->stel = (float)strtod (mn->metafields[6], &endptr);
      if (mn->metafields[7])
        sh->stdp = (float)strtod (mn->metafields[7], &endptr);
      if (mn->metafields[8])
        sh->cmpaz = (float)strtod (mn->metafields[8], &endptr);
      if (mn->metafields[9])
      {
        sh->cmpinc = (float)strtod (mn->metafields[9], &endptr);
        if (seedinc)
          sh->cmpinc += 90;
      }
      if (mn->metafields[10])
        strncpy (sh->kinst, mn->metafields[10], 8);
      if (mn->metafields[11])
        sh->scale = (float)strtod (mn->metafields[11], &endptr);

      retval = 0;
      break;
    }

    mlp = mlp->next;
  }

  return retval;
} /* End of insertSACMetaData() */

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

  if (!sh)
    return -1;

  for (idx = 0; idx < (NUMFLOATHDR + NUMINTHDR); idx++)
  {
    ip = (int32_t *)sh + idx;
    ms_gswap4 (ip);
  }

  return 0;
} /* End of swapSACHeader() */

/***************************************************************************
 * addToString:
 *
 * Concatinate one string to another with a delimiter in-between
 * growing the target string as needed up to a maximum length.  The
 * new addition can be added to either the beggining or end of the
 * string using the where flag:
 *
 * where == 0 means add new addition to end of string
 * where != 0 means add new addition to beginning of string
 *
 * Return 0 on success, -1 on memory allocation error and -2 when
 * string would grow beyond maximum length.
 ***************************************************************************/
static int
addToString (char **string, char *add, char *delim, int where, int maxlen)
{
  int length;
  char *ptr;

  if (!string || !add)
    return -1;

  /* If string is empty, allocate space and copy the addition */
  if (!*string)
  {
    length = strlen (add) + 1;

    if (length > maxlen)
      return -2;

    if ((*string = (char *)malloc (length)) == NULL)
      return -1;

    strcpy (*string, add);
  }
  /* Otherwise add the addition with a delimiter */
  else
  {
    length = strlen (*string) + strlen (delim) + strlen (add) + 1;

    if (length > maxlen)
      return -2;

    if ((ptr = (char *)malloc (length)) == NULL)
      return -1;

    /* Put addition at beginning of the string */
    if (where)
    {
      snprintf (ptr, length, "%s%s%s",
                (add) ? add : "",
                (delim) ? delim : "",
                *string);
    }
    /* Put addition at end of the string */
    else
    {
      snprintf (ptr, length, "%s%s%s",
                *string,
                (delim) ? delim : "",
                (add) ? add : "");
    }

    /* Free previous string and set pointer to newly allocated space */
    free (*string);
    *string = ptr;
  }

  return 0;
} /* End of addToString() */

/***************************************************************************
 * addToProcLog:
 *
 * Format printf-style log message, print to stderr if verbose and add
 * to global process log if requested.
 *
 * Return 0 on success, -1 on errors.
 ***************************************************************************/
static int
addToProcLog (const char *format, ...)
{
  static char message[1024];
  va_list varlist;

  if (!format)
    return -1;

  va_start (varlist, format);
  vsnprintf (message, sizeof (message), format, varlist);
  va_end (varlist);

  if (verbose)
    fprintf (stderr, "%s\n", message);

  /* Add to stored process log buffer, delimit entries with newlines */
  if (proclogfile)
    if (addToString (&proclog, message, "\n", 0, MAXPROCLOG))
      return -1;

  return 0;
} /* End of addToProcLog() */

/***************************************************************************
 * calcStats():
 *
 * Calculate various statistics for the specified data and add the
 * stats to the process log.
 *
 * Running arithmetic mean and standard deviation calculations are
 * based on Knuth's "The Art of Computer Programming: Semi-Numerical
 * Algorithms" who further credits B.P. Welford, "Technometrics 4"
 * (1962), 419-420.
 *
 * The algorithm goes:
 * M[1] = x[1];  S[1] = 0;
 * for (i=2; i<=n; ++i) {
 *   M[i] = M[i-1] + (x[i] - M[i-1]) / i
 *   S[i] = S[i-1] + (x[i] - M[i-1]) * (x[i] - M[i])
 * }
 * SD = sqrt(S[n] / (n - 1));
 *
 * RMS (root-mean-square) is related to the standard deviation (SD)
 * and the arithmetic mean (M) in the following way:
 *
 * RMS^2 = M^2 + SD^2
 *  or
 * RMS = sqrt( M^2 + SD^2 )
 *
 * Return 0 on success and -1 on error.
 ***************************************************************************/
static int
calcStats (void *input, char inputtype, int length)
{
  int idx;
  int32_t *idata;
  float *fdata;
  double *ddata;

  double min;
  double max;
  double M   = 0.0;
  double S   = 0.0;
  double pM  = 0.0;
  double pS  = 0.0;
  double SD  = 0.0;
  double RMS = 0.0;

  if (!input)
    return -1;

  idata = (int32_t *)input;
  fdata = (float *)input;
  ddata = (double *)input;

  /* Check sample type and initialize */
  if (inputtype == 'i')
  {
    min = (double)*idata;
    max = (double)*idata;
    M = pM = (double)*idata;
  }
  else if (inputtype == 'f')
  {
    min = (double)*fdata;
    max = (double)*fdata;
    M = pM = (double)*fdata;
  }
  else if (inputtype == 'd')
  {
    min = (double)*ddata;
    max = (double)*ddata;
    M = pM = (double)*ddata;
  }
  else
  {
    fprintf (stderr, "calcStats: Unsupported sample type: '%c'\n", inputtype);
    return -1;
  }

  /* Determine running min, max, mean, standard deviation and RMS,
   * skip first value of each series, already initialized */
  for (idx = 1; idx < length; idx++)
  {
    if (inputtype == 'i')
    {
      if (*(idata + idx) > max)
        max = *(idata + idx);

      if (*(idata + idx) < min)
        min = *(idata + idx);

      M = pM + (*(idata + idx) - pM) / (idx + 1);
      S = pS + (*(idata + idx) - pM) * (*(idata + idx) - M);
    }
    else if (inputtype == 'f')
    {
      if (*(fdata + idx) > max)
        max = *(fdata + idx);

      if (*(fdata + idx) < min)
        min = *(fdata + idx);

      M = pM + (*(fdata + idx) - pM) / (idx + 1);
      S = pS + (*(fdata + idx) - pM) * (*(fdata + idx) - M);
    }
    else if (inputtype == 'd')
    {
      if (*(ddata + idx) > max)
        max = *(ddata + idx);

      if (*(ddata + idx) < min)
        min = *(ddata + idx);

      M = pM + (*(ddata + idx) - pM) / (idx + 1);
      S = pS + (*(ddata + idx) - pM) * (*(ddata + idx) - M);
    }

    pM = M;
    pS = S;
  }

  if (S < 0.0)
  {
    fprintf (stderr, "Accumulator overload, S: %g\n", S);
    return -1;
  }

  SD  = sqrt (S / (length - 1));
  RMS = sqrt ((M * M) + (SD * SD));

  addToProcLog ("Statistics min: %8.1f, max: %8.1f, mean: %8.1f, SD: %8.1f, RMS: %8.1f",
                min, max, M, SD, RMS);

  return 0;
} /* End of calcStats() */

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

  if (!input || !output)
  {
    fprintf (stderr, "differentiate2(): input or output pointers are NULL\n");
    return -1;
  }

  if (inputtype != 'f' && inputtype != 'd')
  {
    fprintf (stderr, "differentiate2(): unrecognized input sample type: %c\n",
             inputtype);
    return -1;
  }

  if (inputtype == 'f')
  {
    fin  = input;
    fout = output;

    for (idx = 0; idx < (length - 1); idx++)
    {
      fout[idx] = rate * (fin[idx + 1] - fin[idx]);
    }
  }
  else if (inputtype == 'd')
  {
    din  = input;
    dout = output;

    for (idx = 0; idx < (length - 1); idx++)
    {
      dout[idx] = rate * (din[idx + 1] - din[idx]);
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

  if (!input || !output)
  {
    fprintf (stderr, "integrateTrap(): input or output pointers are NULL\n");
    return -1;
  }

  if (inputtype != 'f' && inputtype != 'd')
  {
    fprintf (stderr, "integrateTrap(): unrecognized input sample type: %c\n",
             inputtype);
    return -1;
  }

  if (inputtype == 'f')
  {
    fin  = input;
    fout = output;

    for (idx = 0; idx < (length - 1); idx++)
    {
      prtint    = halfstep * (fin[idx] + fin[idx + 1]);
      totint    = totint + prtint;
      fout[idx] = totint;
    }
  }
  else if (inputtype == 'd')
  {
    din  = input;
    dout = output;

    dout[0] = din[0];

    for (idx = 0; idx < (length - 1); idx++)
    {
      prtint    = halfstep * (din[idx] + din[idx + 1]);
      totint    = totint + prtint;
      dout[idx] = totint;
    }
  }

  return length - 1;
} /* End of integrateTrap() */

/***************************************************************************
 * applyPolynomialM:
 *
 * Apply a Maclaurin-type polynomial to the data samples.
 *
 * output(x) = a0 + a1*x + a2*x^2 + ... + an*x^n
 *
 *  where a0, a1, a2 ... an are the coefficients of the series.
 *  where x = input sample, ouput(x) = output sample.
 *
 * The input data can be either floats or doubles, the inputtype must
 * be set to either 'f' or 'd' indicating how the input data should be
 * treated.
 *
 * The output array can be the same as the input array in which case
 * the input array is replaced.  The output array will be the same
 * type as the input array and should already be allocated.
 *
 * Returns the number of samples in the output array on success and -1
 * on failure.
 ***************************************************************************/
static int
applyPolynomialM (void *input, char inputtype, int length,
                  double *coeff, int numcoeff, void *output)
{
  int idx;
  int cdx;

  float *fin  = (float *)input;
  float *fout = (float *)output;
  float fval;
  double *din  = (double *)input;
  double *dout = (double *)output;
  double dval;

  if (!input || !coeff || !output)
  {
    fprintf (stderr, "applyPolynomialM(): input, coeff or output pointers are NULL\n");
    return -1;
  }

  if (numcoeff < 1) /* Nothing to do with no coefficients */
    return 0;

  if (inputtype != 'f' && inputtype != 'd')
  {
    fprintf (stderr, "applyPolynomialM(): unsupported input sample type: %c\n",
             inputtype);
    return -1;
  }

  if (inputtype == 'f')
  {
    for (idx = 0; idx < length; idx++)
    {
      fval = fin[idx];

      fout[idx] = coeff[0];
      for (cdx = 1; cdx < numcoeff; cdx++)
      {
        fout[idx] += (coeff[cdx] * pow (fval, cdx));
      }
    }
  }
  else if (inputtype == 'd')
  {
    for (idx = 0; idx < length; idx++)
    {
      dval = din[idx];

      dout[idx] = coeff[0];
      for (cdx = 1; cdx < numcoeff; cdx++)
      {
        dout[idx] += (coeff[cdx] * pow (dval, cdx));
      }
    }
  }

  return length;
} /* End of applyPolynomialM() */

/***************************************************************************
 * convertSamples:
 *
 * Convert samples for a specified MS3TraceSeg to specified type.
 *
 * Returns 0 on success, and -1 on failure
 ***************************************************************************/
static int
convertSamples (MS3TraceSeg *seg, char type)
{
  int32_t *idata;
  float *fdata;
  double *ddata;
  int64_t idx;

  if (!seg)
  {
    fprintf (stderr, "convertSamples: Error, no MS3TraceSeg specified!\n");
    return -1;
  }

  idata = (int32_t *)seg->datasamples;
  fdata = (float *)seg->datasamples;
  ddata = (double *)seg->datasamples;

  /* Convert sample type if needed */
  if (seg->sampletype != type)
  {
    if (seg->sampletype == 'a' || type == 'a')
    {
      fprintf (stderr, "Error, cannot convert ASCII samples to/from numeric type\n");
      return -1;
    }

    /* Convert to integers */
    else if (type == 'i')
    {
      if (seg->sampletype == 'f') /* Convert floats to integers with simple rounding */
      {
        addToProcLog ("Converting segment samples from floats to integers");

        for (idx = 0; idx < seg->numsamples; idx++)
        {
          /* Check for loss of sub-integer */
          if ((fdata[idx] - (int32_t)fdata[idx]) > 0.000001)
          {
            fprintf (stderr, "Warning, Loss of precision when converting floats to integers, loss: %g\n",
                     (fdata[idx] - (int32_t)fdata[idx]));
            return -1;
          }

          idata[idx] = (int32_t) (fdata[idx] + 0.5);
        }
      }
      else if (seg->sampletype == 'd') /* Convert doubles to integers with simple rounding */
      {
        addToProcLog ("Converting segment samples from doubles to integers");

        for (idx = 0; idx < seg->numsamples; idx++)
        {
          /* Check for loss of sub-integer */
          if ((ddata[idx] - (int32_t)ddata[idx]) > 0.000001)
          {
            fprintf (stderr, "Warning, Loss of precision when converting doubles to integers, loss: %g\n",
                     (ddata[idx] - (int32_t)ddata[idx]));
            return -1;
          }

          idata[idx] = (int32_t) (ddata[idx] + 0.5);
        }

        /* Reallocate buffer for reduced size needed */
        if (!(seg->datasamples = realloc (seg->datasamples, (seg->numsamples * sizeof (int32_t)))))
        {
          fprintf (stderr, "Error, cannot re-allocate buffer for sample conversion\n");
          return -1;
        }
      }

      seg->sampletype = 'i';
    }

    /* Convert to floats */
    else if (type == 'f')
    {
      if (seg->sampletype == 'i') /* Convert integers to floats */
      {
        addToProcLog ("Converting segment samples from integers to floats");

        for (idx = 0; idx < seg->numsamples; idx++)
          fdata[idx] = (float)idata[idx];
      }
      else if (seg->sampletype == 'd') /* Convert doubles to floats */
      {
        addToProcLog ("Converting segment samples from doubles to floats");

        for (idx = 0; idx < seg->numsamples; idx++)
          fdata[idx] = (float)ddata[idx];

        /* Reallocate buffer for reduced size needed */
        if (!(seg->datasamples = realloc (seg->datasamples, (seg->numsamples * sizeof (float)))))
        {
          fprintf (stderr, "Error, cannot re-allocate buffer for sample conversion\n");
          return -1;
        }
      }

      seg->sampletype = 'f';
    }
    /* Convert to doubles */
    else if (type == 'd')
    {
      if (!(ddata = (double *)malloc (seg->numsamples * sizeof (double))))
      {
        fprintf (stderr, "Error, cannot allocate buffer for sample conversion to doubles\n");
        return -1;
      }

      if (seg->sampletype == 'i') /* Convert integers to doubles */
      {
        addToProcLog ("Converting segment samples from integers to doubles");

        for (idx = 0; idx < seg->numsamples; idx++)
          ddata[idx] = (double)idata[idx];

        free (idata);
      }
      else if (seg->sampletype == 'f') /* Convert floats to doubles */
      {
        addToProcLog ("Converting segment samples from floats to doubles");

        for (idx = 0; idx < seg->numsamples; idx++)
          ddata[idx] = (double)fdata[idx];

        free (fdata);
      }

      seg->datasamples = ddata;
      seg->sampletype  = 'd';
    }
  }

  return 0;
} /* End of convertSamples() */

/***************************************************************************
 * parameterProc():
 * Process the command line parameters.
 *
 * Returns 0 on success, and -1 on failure
 ***************************************************************************/
static int
parameterProc (int argcount, char **argvec)
{
  char string1[10];
  char string2[10];
  char *filename;
  char *filtstr      = NULL;
  char *freqlimitstr = NULL;
  char *dBdownstr    = NULL;
  char *tptr;
  int optind;
  int idx;

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
      usage ();
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
    else if (strcmp (argvec[optind], "-ts") == 0)
    {
      starttime = ms_timestr2nstime (getOptVal (argcount, argvec, optind++, 0));
      if (starttime == NSTERROR)
        return -1;
    }
    else if (strcmp (argvec[optind], "-te") == 0)
    {
      endtime = ms_timestr2nstime (getOptVal (argcount, argvec, optind++, 0));
      if (endtime == NSTERROR)
        return -1;
    }
    else if (strcmp (argvec[optind], "-Im") == 0)
    {
      informat = 1;
    }
    else if (strcmp (argvec[optind], "-Is") == 0)
    {
      informat = 2;
    }
    else if (strcmp (argvec[optind], "-MSEED") == 0)
    {
      dataformat = 1;
    }
    else if (strcmp (argvec[optind], "-MR") == 0)
    {
      packreclen = strtol (getOptVal (argcount, argvec, optind++, 0), NULL, 10);
    }
    else if (strcmp (argvec[optind], "-ME") == 0)
    {
      packencoding = strtol (getOptVal (argcount, argvec, optind++, 0), NULL, 10);
    }
    else if (strcmp (argvec[optind], "-MF") == 0)
    {
      msformat = strtol (getOptVal (argcount, argvec, optind++, 0), NULL, 10);
    }
    else if (strcmp (argvec[optind], "-Sif") == 0)
    {
      sacinformat = strtoul (getOptVal (argcount, argvec, optind++, 0), NULL, 10);
    }
    else if (strcmp (argvec[optind], "-Sf") == 0)
    {
      sacoutformat = strtoul (getOptVal (argcount, argvec, optind++, 0), NULL, 10);
    }
    else if (strcmp (argvec[optind], "-m") == 0)
    {
      metadatafile = getOptVal (argcount, argvec, optind++, 0);
    }
    else if (strcmp (argvec[optind], "-msi") == 0)
    {
      seedinc = 1;
    }
    else if (strcmp (argvec[optind], "-o") == 0)
    {
      outputfile = getOptVal (argcount, argvec, optind++, 1);
    }
    else if (strcmp (argvec[optind], "-od") == 0)
    {
      outputdir = getOptVal (argcount, argvec, optind++, 0);
    }
    else if (strcmp (argvec[optind], "-plog") == 0)
    {
      proclogfile = getOptVal (argcount, argvec, optind++, 1);
    }
    else if (strcmp (argvec[optind], "-c") == 0)
    {
      channel = getOptVal (argcount, argvec, optind++, 0);
    }
    else if (strcmp (argvec[optind], "-CR") == 0)
    {
      filename = getOptVal (argcount, argvec, optind++, 0);
      if (filename)
        addProcess (PROC_CONVOLVE, filename, NULL, PROC_CONVRESP, 0, 0.0, 0.0);
    }
    else if (strcmp (argvec[optind], "-DR") == 0)
    {
      filename = getOptVal (argcount, argvec, optind++, 0);
      if (filename)
        addProcess (PROC_CONVOLVE, filename, NULL, PROC_DECONVRESP, 0, 0.0, 0.0);
    }
    else if (strcmp (argvec[optind], "-CS") == 0)
    {
      filename = getOptVal (argcount, argvec, optind++, 0);
      if (filename)
        addProcess (PROC_CONVOLVE, filename, NULL, PROC_CONVSAC, 0, 0.0, 0.0);
    }
    else if (strcmp (argvec[optind], "-DS") == 0)
    {
      filename = getOptVal (argcount, argvec, optind++, 0);
      if (filename)
        addProcess (PROC_CONVOLVE, filename, NULL, PROC_DECONVSAC, 0, 0.0, 0.0);
    }
    else if (strcmp (argvec[optind], "-FL") == 0)
    {
      freqlimitstr = getOptVal (argcount, argvec, optind++, 1);
    }
    else if (strcmp (argvec[optind], "-FLa") == 0)
    {
      dBdownstr = getOptVal (argcount, argvec, optind++, 1);
    }
    else if (strcmp (argvec[optind], "-W") == 0)
    {
      if (prewhiten < 0)
        prewhiten *= -1;
      else
        prewhiten = 6;
    }
    else if (strcmp (argvec[optind], "-Wo") == 0)
    {
      if (prewhiten)
        prewhiten = strtol (getOptVal (argcount, argvec, optind++, 0), NULL, 10);
      else
        prewhiten = -1 * strtol (getOptVal (argcount, argvec, optind++, 0), NULL, 10);
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
      respunits = getOptVal (argcount, argvec, optind++, 0);
    }
    else if (strcmp (argvec[optind], "-LP") == 0)
    {
      filtstr = getOptVal (argcount, argvec, optind++, 0);
      addProcess (PROC_LPFILTER, filtstr, NULL, 1, 0, 0.0, 0.0);
    }
    else if (strcmp (argvec[optind], "-LP1") == 0)
    {
      filtstr = getOptVal (argcount, argvec, optind++, 0);
      addProcess (PROC_LPFILTER, filtstr, NULL, 0, 0, 0.0, 0.0);
    }
    else if (strcmp (argvec[optind], "-HP") == 0)
    {
      filtstr = getOptVal (argcount, argvec, optind++, 0);
      addProcess (PROC_HPFILTER, filtstr, NULL, 1, 0, 0.0, 0.0);
    }
    else if (strcmp (argvec[optind], "-HP1") == 0)
    {
      filtstr = getOptVal (argcount, argvec, optind++, 0);
      addProcess (PROC_HPFILTER, filtstr, NULL, 0, 0, 0.0, 0.0);
    }
    else if (strcmp (argvec[optind], "-BP") == 0)
    {
      char *hpfiltstr = getOptVal (argcount, argvec, optind++, 0);
      char *lpfiltstr = strchr (hpfiltstr, ':');
      if (lpfiltstr)
      {
        *lpfiltstr++ = '\0';
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
      char *hpfiltstr = getOptVal (argcount, argvec, optind++, 0);
      char *lpfiltstr = strchr (hpfiltstr, ':');
      if (lpfiltstr)
      {
        *lpfiltstr++ = '\0';
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
      double scalefactor = strtod (getOptVal (argcount, argvec, optind++, 1), NULL);
      addProcess (PROC_SCALE, NULL, NULL, 0, 0, scalefactor, 0.0);
    }
    else if (strcmp (argvec[optind], "-SI") == 0)
    {
      double scalefactor = strtod (getOptVal (argcount, argvec, optind++, 1), NULL);
      addProcess (PROC_SCALE, NULL, NULL, 0, 0, (scalefactor) ? 1.0 / scalefactor : 0.0, 0.0);
    }
    else if (strcmp (argvec[optind], "-DEC") == 0)
    {
      int decimfactor = strtol (getOptVal (argcount, argvec, optind++, 0), NULL, 10);
      addProcess (PROC_DECIMATE, NULL, NULL, decimfactor, 0, 0.0, 0.0);
    }
    else if (strcmp (argvec[optind], "-TAP") == 0)
    {
      double taperwidth = 0.0;
      int tapertype     = TAPER_HANNING;

      char *taperwidthstr = getOptVal (argcount, argvec, optind++, 0);
      char *tapertypestr  = strchr (taperwidthstr, ':');

      if (tapertypestr)
        *tapertypestr++ = '\0';

      taperwidth = strtod (taperwidthstr, NULL);

      if (tapertypestr && *tapertypestr)
      {
        if (!strcasecmp (tapertypestr, "HANNING"))
          tapertype = TAPER_HANNING;
        else if (!strcasecmp (tapertypestr, "HAMMING"))
          tapertype = TAPER_HAMMING;
        else if (!strcasecmp (tapertypestr, "COSINE"))
          tapertype = TAPER_COSINE;
        else
        {
          fprintf (stderr, "Unrecognized taper type: '%s'\n", tapertypestr);
          exit (1);
        }
      }

      addProcess (PROC_TAPER, NULL, NULL, tapertype, 0, taperwidth, 0.0);
    }
    else if (strcmp (argvec[optind], "-POLYM") == 0)
    {
      char *coeffstr = getOptVal (argcount, argvec, optind++, 1);
      addProcess (PROC_POLYNOMIALM, coeffstr, NULL, 0, 0, 0.0, 0.0);
    }
    else if (strcmp (argvec[optind], "-ENV") == 0)
    {
      addProcess (PROC_ENVELOPE, NULL, NULL, 0, 0, 0.0, 0.0);
    }
    else if (strcmp (argvec[optind], "-DTRIM") == 0)
    {
      addProcess (PROC_DATATRIM, NULL, NULL, 0, 0, 0.0, 0.0);
    }
    else if (strcmp (argvec[optind], "-ROTATE") == 0)
    {
      /* -ROTATE "E/1,N/2,Z/3:azimuth,incidence" */

      /* Component codes:
	   * This should be a string of comma-separated, single-characters representing
	   * the E, N and Z components, in that order regardless of their actual codes.
	   * Each code may optionally be follow by a / and another code represting the
	   * code that should be used for the resulting segment.
	   *
	   * Examples:
	   *   "E,N,Z"
	   *   "1,2,Z"
	   *   "1/T,2/Q,Z/L"
	   *   "E/T,N/R"
	   *
	   * Rotation angles:
	   * This should be a string of two floating point values separated by a comma.
	   * The first value is the angle of horizontal rotation (azimuth).  The second,
	   * optional, value is the angle of vertical tilt (incidence).
	   */

      char *components = getOptVal (argcount, argvec, optind++, 0);
      char *angles     = strchr (components, ':');
      double azimuth   = 0.0;
      double incidence = 0.0;

      if (angles)
        *angles++ = '\0';

      /* Parse rotation components and replacement codes */
      memset (string1, 0, sizeof (string1));
      memset (string2, 0, sizeof (string2));
      tptr = components;
      for (idx = 0; tptr && idx < 3; idx++)
      {
        string1[idx] = tptr[0];

        if (strlen (tptr) >= 2)
        {
          if (tptr[1] == '/')
            string2[idx] = tptr[2];
        }

        if ((tptr = strchr (tptr, ',')))
          tptr++;
      }

      /* Parse the azimuth and incidence angles */
      if (angles)
      {
        tptr = strchr (angles, ',');

        if (tptr)
        {
          *tptr++   = '\0';
          incidence = strtod (tptr, NULL);
        }

        azimuth = strtod (angles, NULL);
      }
      else
      {
        fprintf (stderr, "Rotation operation requires at least an azimuth\n");
        exit (1);
      }

      if (!string1[1])
      {
        fprintf (stderr, "Rotation requires at least two horizontal components\n");
        exit (1);
      }
      if (incidence && !string1[2])
      {
        fprintf (stderr, "3-D Rotation requires three components\n");
        exit (1);
      }

      /* The components will be parsed in addProcess */
      addProcess (PROC_ROTATE, string1, string2, 0, 0, azimuth, incidence);
    }
    else if (strcmp (argvec[optind], "-STATS") == 0)
    {
      addProcess (PROC_STATS, NULL, NULL, 0, 0, 0.0, 0.0);
    }
    else if (strncmp (argvec[optind], "-", 1) == 0 &&
             strlen (argvec[optind]) > 1)
    {
      fprintf (stderr, "Unknown option: %s\n", argvec[optind]);
      exit (1);
    }
    else
    {
      tptr = argvec[optind];

      /* Check for an input file list */
      if (tptr[0] == '@')
      {
        if (addListFile (tptr + 1) < 0)
        {
          fprintf (stderr, "Error adding list file %s", tptr + 1);
          exit (1);
        }
      }
      /* Otherwise this is an input file */
      else
      {
        /* Add file to global file list */
        if (addFile (tptr))
        {
          fprintf (stderr, "Error adding file to input list %s", tptr);
          exit (1);
        }
      }
    }
  }

  /* Make sure input files were specified */
  if (filelist == 0)
  {
    fprintf (stderr, "No input files were specified\n\n");
    fprintf (stderr, "%s version %s\n\n", PACKAGE, VERSION);
    fprintf (stderr, "Try %s -h for usage\n", PACKAGE);
    exit (1);
  }

  /* Make sure output directory exists */
  if (outputdir)
  {
    if (access (outputdir, W_OK))
    {
      fprintf (stderr, "Cannot write to output diretory: %s (%s)\n",
               outputdir, strerror (errno));
      exit (1);
    }
  }

  /* Set prewhitening if only -Wo was specified without -W */
  if (prewhiten < 0)
    prewhiten *= -1;

  /* Parse frequency limits string */
  if (freqlimitstr)
  {
    int parsed = 0;

    if (!freqlimit)
    {
      if ((freqlimit = (double *)malloc (4 * sizeof (double))) == NULL)
      {
        fprintf (stderr, "Error allocating memory\n");
        exit (1);
      }

      freqlimit[0] = -1.0;
      freqlimit[1] = -1.0;
      freqlimit[2] = -1.0;
      freqlimit[3] = -1.0;
    }

    parsed = sscanf (freqlimitstr, "%lf/%lf/%lf/%lf",
                     &freqlimit[0], &freqlimit[1],
                     &freqlimit[2], &freqlimit[3]);

    if (parsed != 4)
    {
      fprintf (stderr, "Deconvolution frequency limits incorrectly: %s\n", freqlimitstr);
      fprintf (stderr, "Try %s -h for usage\n", PACKAGE);
      exit (1);
    }

    if (freqlimit[0] > freqlimit[1])
    {
      fprintf (stderr, "Deconvolution frequency limits specified incorrectly: %s\n", freqlimitstr);
      fprintf (stderr, "Cut-off frequency of lower bound (%g) cannot be greater than pass frequency (%g)\n",
               freqlimit[0], freqlimit[1]);
      exit (1);
    }

    if (freqlimit[2] > freqlimit[3])
    {
      fprintf (stderr, "Deconvolution freqency limits specified incorrectly: %s\n", freqlimitstr);
      fprintf (stderr, "Cut-off frequency of upper bound (%g) cannot be less than pass frequency (%g)\n",
               freqlimit[3], freqlimit[2]);
      exit (1);
    }
  }

  /* Parse dB down auto limit cutoff */
  if (dBdownstr)
  {
    if (!freqlimit)
    {
      if ((freqlimit = (double *)malloc (4 * sizeof (double))) == NULL)
      {
        fprintf (stderr, "Error allocating memory\n");
        exit (1);
      }

      freqlimit[0] = -1.0;
      freqlimit[1] = -1.0;
      freqlimit[2] = -1.0;
      freqlimit[3] = -1.0;
    }

    /* dB down cutoff string is lowercorner[/uppercorner] */

    /* Check for lower & upper corner cutoff separator */
    tptr = strchr (dBdownstr, '/');

    /* Parse lowercorner/uppercorner cutoff pair */
    if (tptr)
    {
      *tptr++  = '\0';
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
  if (respunits)
  {
    if (strcmp (respunits, "DIS") && strcmp (respunits, "VEL") &&
        strcmp (respunits, "ACC") && strcmp (respunits, "DEF"))
    {
      fprintf (stderr, "Unrecognized SEED RESP response units: %s\n", respunits);
      fprintf (stderr, "Should be one of DIS, VEL, ACC or DEF\n");
      exit (1);
    }
  }

  /* Default is no unit conversion */
  if (!respunits)
  {
    respunits = "DEF";
  }

  /* Report the program version */
  if (verbose)
    fprintf (stderr, "%s version: %s\n", PACKAGE, VERSION);

  /* Read metadata file if specified */
  if (metadatafile)
  {
    if (readMetaData (metadatafile))
    {
      fprintf (stderr, "Error reading metadata file\n");
      return -1;
    }
  }

  return 0;
} /* End of parameterProc() */

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
  if (argvec == NULL || argvec[argopt] == NULL)
  {
    fprintf (stderr, "getOptVal(): NULL option requested\n");
    exit (1);
    return 0;
  }

  /* When the value potentially starts with a dash (-) */
  if ((argopt + 1) < argcount && dasharg)
    return argvec[argopt + 1];

  /* Otherwise check that the value is not another option */
  if ((argopt + 1) < argcount && *argvec[argopt + 1] != '-')
    return argvec[argopt + 1];

  fprintf (stderr, "Option %s requires a value\n", argvec[argopt]);
  exit (1);
  return 0;
} /* End of getOptVal() */

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
  char yearstr[5] = {0};
  int idx, count;
  int linecount = 0;

  if (!metafile)
    return -1;

  if ((mfp = fopen (metafile, "rb")) == NULL)
  {
    fprintf (stderr, "Cannot open metadata output file: %s (%s)\n",
             metafile, strerror (errno));
    return -1;
  }

  if (verbose)
    fprintf (stderr, "Reading station/channel metadata from %s\n", metafile);

  while (fgets (line, sizeof (line), mfp))
  {
    linecount++;

    /* Truncate at line return if any */
    if ((fp = strchr (line, '\n')))
      *fp = '\0';

    /* Count the number of commas */
    count = 0;
    fp    = line;
    while ((fp = strchr (fp, ',')))
    {
      count++;
      fp++;
    }

    /* Must have at least 3 commas for Net, Sta, Loc, Chan ... */
    if (count < 3)
    {
      if (verbose > 1)
        fprintf (stderr, "Skipping metadata line: %s\n", line);
      continue;
    }

    /* Check for comment line beginning with '#' */
    if (line[0] == '#')
    {
      if (verbose > 1)
        fprintf (stderr, "Skipping comment line: %s\n", line);
      continue;
    }

    /* Create a copy of the line */
    lineptr = strdup (line);

    mn.metafields[0] = fp = lineptr;
    mn.starttime          = NSTUNSET;
    mn.endtime            = NSTUNSET;

    /* Separate line on commas and index in metafields array */
    for (idx = 1; idx < MAXMETAFIELDS; idx++)
    {
      mn.metafields[idx] = NULL;

      if (fp)
      {
        if ((fp = strchr (fp, ',')))
        {
          *fp++ = '\0';

          if (*fp != ',' && *fp != '\0')
            mn.metafields[idx] = fp;
        }
      }
    }

    /* Trim last field if more fields exist */
    if (fp && (fp = strchr (fp, ',')))
      *fp = '\0';

    /* Sanity check, source name fields must be populated */
    for (idx = 0; idx <= 3; idx++)
    {
      if (mn.metafields[idx] == NULL)
      {
        fprintf (stderr, "Error, field %d cannot be empty in metadata file line %d\n",
                 idx + 1, linecount);
        fprintf (stderr, "Perhaps a wildcard character (*) was the intention?\n");

        exit (1);
      }
    }

    /* Parse and convert start time */
    if (mn.metafields[15])
    {
      /* Check if first four characters are a year before 1900 and set NSTUNSET (open) time */
      strncpy (yearstr, mn.metafields[15], 4);
      if (strlen (yearstr) >= 4 &&
          isdigit (yearstr[0]) &&
          isdigit (yearstr[1]) &&
          isdigit (yearstr[2]) &&
          isdigit (yearstr[3]) &&
          strtoul (yearstr, NULL, 10) < 1900)
      {
        if (verbose > 1)
          fprintf (stderr, "Setting open start time for year: %s\n", yearstr);

        mn.starttime = NSTUNSET;
      }

      else if ((mn.starttime = ms_timestr2nstime (mn.metafields[15])) == NSTERROR)
      {
        fprintf (stderr, "Error parsing metadata start time: '%s'\n", mn.metafields[15]);
        exit (1);
      }
    }

    /* Parse and convert end time */
    if (mn.metafields[16])
    {
      /* Check if first four characters are a year beyond 2100 and set NSTUNSET (open) time */
      strncpy (yearstr, mn.metafields[16], 4);
      if (strlen (yearstr) >= 4 &&
          isdigit (yearstr[0]) &&
          isdigit (yearstr[1]) &&
          isdigit (yearstr[2]) &&
          isdigit (yearstr[3]) &&
          strtoul (yearstr, NULL, 10) > 2100)
      {
        if (verbose > 1)
          fprintf (stderr, "Setting open end time for year: %s\n", yearstr);

        mn.endtime = NSTUNSET;
      }

      else if ((mn.endtime = ms_timestr2nstime (mn.metafields[16])) == NSTERROR)
      {
        fprintf (stderr, "Error parsing metadata end time: '%s'\n", mn.metafields[16]);
        exit (1);
      }
    }

    /* Add the metanode to the metadata list */
    if (!addMetaNode (&metadata, NULL, 0, &mn, sizeof (struct metanode)))
    {
      fprintf (stderr, "Error adding metadata fields to list\n");
    }
  }

  fclose (mfp);

  return 0;
} /* End of readMetaData() */

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

  if (data == NULL)
  {
    fprintf (stderr, "addMetaNode(): No data specified\n");
    return NULL;
  }

  lastlp = *listroot;
  while (lastlp != 0)
  {
    if (lastlp->next == 0)
      break;

    lastlp = lastlp->next;
  }

  /* Create new sacmetanode */
  newlp = (struct metalist *)malloc (sizeof (struct metalist));
  memset (newlp, 0, sizeof (struct metalist));

  if (key)
  {
    newlp->key = malloc (keylen);
    memcpy (newlp->key, key, keylen);
  }

  if (data)
  {
    newlp->data = malloc (datalen);
    memcpy (newlp->data, data, datalen);
  }

  newlp->next = 0;

  if (lastlp == 0)
    *listroot = newlp;
  else
    lastlp->next = newlp;

  return newlp;
} /* End of addMetaNode() */

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
  char header[64] = {0};
  uint8_t formatversion;

  if (!filename)
  {
    fprintf (stderr, "addFile(): No file name specified\n");
    return -1;
  }

  newlp = (struct filelink *)calloc (1, sizeof (struct filelink));

  if (!newlp)
  {
    fprintf (stderr, "addFile(): Error allocating memory\n");
    return -1;
  }

  newlp->filename = strdup (filename);

  if (!newlp->filename)
  {
    fprintf (stderr, "addFile(): Error duplicating string\n");
    return -1;
  }

  /* Determine file format: open and read first some data to test */
  if ((ifp = fopen (filename, "rb")) == NULL)
  {
    fprintf (stderr, "Cannot open input file: %s (%s)\n", filename, strerror (errno));
    free (newlp->filename);
    free (newlp);
    return -1;
  }

  if (fread (header, sizeof (header), 1, ifp) != 1)
  {
    fprintf (stderr, "Cannot read %s\n", filename);
    free (newlp->filename);
    free (newlp);
    return -1;
  }

  fclose (ifp);

  /* Set input format or check for miniSEED, otherwise default to SAC */
  if (informat)
    newlp->format = informat;
  else if (ms3_detect (header, sizeof (header), &formatversion) > 0)
    newlp->format = 1;
  else
    newlp->format = 2;

  /* Add new file to the end of the list */
  if (filelisttail == 0)
  {
    filelist     = newlp;
    filelisttail = newlp;
  }
  else
  {
    filelisttail->next = newlp;
    filelisttail       = newlp;
  }

  return 0;
} /* End of addFile() */

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

  if (verbose >= 1)
    fprintf (stderr, "Reading list file '%s'\n", filename);

  if (!(fp = fopen (filename, "rb")))
  {
    fprintf (stderr, "Cannot open list file %s: %s\n", filename, strerror (errno));
    return -1;
  }

  while (fgets (filelistent, sizeof (filelistent), fp))
  {
    char *cp;

    /* End string at first newline character */
    if ((cp = strchr (filelistent, '\n')))
      *cp = '\0';

    /* Skip empty lines */
    if (!strlen (filelistent))
      continue;

    /* Skip comment lines */
    if (*filelistent == '#')
      continue;

    if (verbose > 1)
      fprintf (stderr, "Adding '%s' from list file\n", filelistent);

    if (addFile (filelistent))
      return -1;

    filecount++;
  }

  fclose (fp);

  return filecount;
} /* End of addListFile() */

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
  char *start    = 0;
  char *stop     = 0;
  char *tptr, *eptr;
  double lpcutoff      = 0.0;
  double hpcutoff      = 0.0;
  double scalefactor   = 0.0;
  int filetype         = 0;
  int lporder          = 0;
  int hporder          = 0;
  int respstart        = -1;
  int respstop         = -1;
  int reverseflag      = 0;
  int decimfactor      = -1;
  int tapertype        = -1;
  double taperwidth    = 0.0;
  double *coefficients = NULL;
  int coefficientcount = 0;
  char rotateENZ[3]    = {'\0', '\0', '\0'};
  char rotatedENZ[3]   = {'\0', '\0', '\0'};
  double rotations[2]  = {0.0, 0.0};

  /* (De)Convolution */
  if (type == PROC_CONVOLVE)
  {
    /* Parse stage limits from RESP file name */
    if (ivalue1 == PROC_CONVRESP || ivalue1 == PROC_DECONVRESP)
    {
      /* string1 == filename, ivalue1 == sub/filetype */

      if ((start = strchr (string1, ':')))
      {
        *start++ = '\0';

        if ((stop = strchr (start, ':')))
        {
          *stop++ = '\0';
        }
      }

      if (start)
        respstart = strtol (start, NULL, 10);
      if (stop)
        respstop = strtol (stop, NULL, 10);

      filename = string1;
      filetype = ivalue1;
    }

    /* SAC poles and zeros */
    if (ivalue1 == PROC_CONVSAC || ivalue1 == PROC_DECONVSAC)
    {
      /* string1 == filename, ivalue1 == sub/filetype */

      filename = string1;
      filetype = ivalue1;
    }
  }
  /* Low-pass filter */
  if (type == PROC_LPFILTER)
  {
    /* string1 == lpfilter, ivalue1 == reverseflag */

    /* Check for slash character indicating an order was specified */
    if ((tptr = strchr (string1, '/')))
    {
      *tptr++  = '\0';
      lpcutoff = strtod (string1, NULL);
      lporder  = strtoul (tptr, NULL, 10);
    }
    else
    {
      lpcutoff = strtod (string1, NULL);
      lporder  = DEFAULT_FILTER_ORDER;
    }

    reverseflag = ivalue1;
  }

  /* High-pass filter */
  if (type == PROC_HPFILTER)
  {
    /* string1 == hpfilter, ivalue1 == reverseflag */

    /* Check for slash character indicating an order was specified */
    if ((tptr = strchr (string1, '/')))
    {
      *tptr++  = '\0';
      hpcutoff = strtod (string1, NULL);
      hporder  = strtoul (tptr, NULL, 10);
    }
    else
    {
      hpcutoff = strtod (string1, NULL);
      hporder  = DEFAULT_FILTER_ORDER;
    }

    reverseflag = ivalue1;
  }

  /* Band-pass filter */
  if (type == PROC_BPFILTER)
  {
    /* string1 == lpfilter, string2 == hpfilter, ivalue1 == reverseflag */

    /* Check for slash character indicating a LP order was specified */
    if ((tptr = strchr (string1, '/')))
    {
      *tptr++  = '\0';
      lpcutoff = strtod (string1, NULL);
      lporder  = strtoul (tptr, NULL, 10);
    }
    else
    {
      lpcutoff = strtod (string1, NULL);
      lporder  = DEFAULT_FILTER_ORDER;
    }

    /* Check for slash character indicating a HP order was specified */
    if ((tptr = strchr (string2, '/')))
    {
      *tptr++  = '\0';
      hpcutoff = strtod (string2, NULL);
      hporder  = strtoul (tptr, NULL, 10);
    }
    else
    {
      hpcutoff = strtod (string2, NULL);
      hporder  = DEFAULT_FILTER_ORDER;
    }

    reverseflag = ivalue1;
  }

  /* Scale factor */
  if (type == PROC_SCALE)
  {
    /* dvalue1 == scalefactor */

    scalefactor = dvalue1;
  }

  /* Decimation factor */
  if (type == PROC_DECIMATE)
  {
    /* ivalue1 == decimfactor */

    decimfactor = ivalue1;
  }

  /* Taper */
  if (type == PROC_TAPER)
  {
    /* ivalue1 == tapertype */
    /* dvalue1 == taperwidth */

    tapertype  = ivalue1;
    taperwidth = dvalue1;
  }

  /* Polynomial */
  if (type == PROC_POLYNOMIALM)
  {
    /* string1 == c0,c1,c2,... */
    double cval;

    tptr = string1;
    while (tptr)
    {
      cval = strtod (tptr, &eptr);
      if (eptr != tptr)
      {
        coefficientcount++;
        if (!(coefficients = realloc (coefficients, coefficientcount * sizeof (double))))
        {
          fprintf (stderr, "Error allocating memory\n");
          exit (1);
        }

        coefficients[coefficientcount - 1] = cval;
      }

      if ((tptr = strchr (tptr, ',')))
        tptr++;
    }
  }

  /* Rotation */
  if (type == PROC_ROTATE)
  {
    /* string1 == pre-rotation components */
    /* string2 == post-rotation components */
    /* dvalue1 == azimuth */
    /* dvalue2 == incidence */

    rotateENZ[0]  = string1[0];
    rotateENZ[1]  = string1[1];
    rotateENZ[2]  = string1[2];
    rotatedENZ[0] = string2[0];
    rotatedENZ[1] = string2[1];
    rotatedENZ[2] = string2[2];
    rotations[0]  = dvalue1;
    rotations[1]  = dvalue2;
  }

  /* Find last process in list */
  lastlp = proclist;
  while (lastlp != 0)
  {
    if (lastlp->next == 0)
      break;

    lastlp = lastlp->next;
  }

  /* If this is a paired deconvolution-convolution option add operation to last entry */
  if (type == PROC_CONVOLVE && lastlp && lastlp->type == PROC_CONVOLVE && lastlp->filetype[1] == 0 &&
      (lastlp->filetype[0] == PROC_DECONVRESP || lastlp->filetype[0] == PROC_DECONVSAC) &&
      (filetype == PROC_CONVRESP || filetype == PROC_CONVSAC))
  {
    if (filename)
      lastlp->filename[1] = strdup (filename);
    else
      lastlp->filename[1] = NULL;
    lastlp->filetype[1] = filetype;
  }
  /* Otherwise add a new processing entry */
  else
  {
    /* Allocate and populate new process entry */
    newlp = (struct proclink *)malloc (sizeof (struct proclink));

    if (!newlp)
    {
      fprintf (stderr, "Error allocating memory\n");
      exit (1);
    }

    newlp->type = type;
    if (filename)
      newlp->filename[0] = strdup (filename);
    else
      newlp->filename[0] = NULL;
    newlp->filetype[0]      = filetype;
    newlp->filename[1]      = NULL;
    newlp->filetype[1]      = 0;
    newlp->respstart        = respstart;
    newlp->respstop         = respstop;
    newlp->lpcutoff         = lpcutoff;
    newlp->lporder          = lporder;
    newlp->hpcutoff         = hpcutoff;
    newlp->hporder          = hporder;
    newlp->reverseflag      = reverseflag;
    newlp->scalefactor      = scalefactor;
    newlp->decimfactor      = decimfactor;
    newlp->tapertype        = tapertype;
    newlp->taperwidth       = taperwidth;
    newlp->coefficients     = coefficients;
    newlp->coefficientcount = coefficientcount;
    newlp->rotateENZ[0]     = rotateENZ[0];
    newlp->rotateENZ[1]     = rotateENZ[1];
    newlp->rotateENZ[2]     = rotateENZ[2];
    newlp->rotatedENZ[0]    = rotatedENZ[0];
    newlp->rotatedENZ[1]    = rotatedENZ[1];
    newlp->rotatedENZ[2]    = rotatedENZ[2];
    newlp->rotations[0]     = rotations[0];
    newlp->rotations[1]     = rotations[1];

    newlp->next = 0;

    if (lastlp == 0)
      proclist = newlp;
    else
      lastlp->next = newlp;
  }

} /* End of addProcess() */

/***************************************************************************
 * procDescription:
 * Return pointer to a text description of a processing number.
 ***************************************************************************/
static char *
procDescription (int proctype)
{
  switch (proctype)
  {

  case PROC_STATS:
    return "Calculate statistics";
    break;
  case PROC_LPFILTER:
    return "Lowpass filter";
    break;
  case PROC_HPFILTER:
    return "Highpass filter";
    break;
  case PROC_BPFILTER:
    return "Bandpass filter";
    break;
  case PROC_CONVOLVE:
    return "Convolve/Deconvolve";
    break;
  case PROC_CONVRESP:
    return "Convole RESP";
    break;
  case PROC_DECONVRESP:
    return "Deconvolve RESP";
    break;
  case PROC_CONVSAC:
    return "Convolve SACPZ";
    break;
  case PROC_DECONVSAC:
    return "Deconvolve SACPZ";
    break;
  case PROC_DIFF2:
    return "Differentiate";
    break;
  case PROC_INTTRAP:
    return "Integrate";
    break;
  case PROC_RMEAN:
    return "Remove mean";
    break;
  case PROC_SCALE:
    return "Scale";
    break;
  case PROC_DECIMATE:
    return "Decimate";
    break;
  case PROC_TAPER:
    return "Taper";
    break;
  case PROC_POLYNOMIALM:
    return "Apply polynomial";
    break;
  case PROC_ENVELOPE:
    return "Envelope";
    break;
  case PROC_DATATRIM:
    return "Data trim";
    break;
  case PROC_ROTATE:
    return "Rotate";
    break;
  default:
    return "Unknown operation";
    break;
  }
} /* End of procDescription() */

/***************************************************************************
 * recordHandler:
 * Saves passed records to the output file.
 ***************************************************************************/
static void
recordHandler (char *record, int reclen, void *vofp)
{
  FILE *ofp = (FILE *)vofp;

  if (ofp)
  {
    if (fwrite (record, reclen, 1, ofp) != 1)
    {
      fprintf (stderr, "Error writing miniSEED to output file\n");
    }

    outputbytes += reclen;
  }
} /* End of recordHandler() */

/***************************************************************************
 * usage():
 *
 * Print the usage message.
 ***************************************************************************/
static void
usage (void)
{
  fprintf (stderr, "%s - time series signal processor: %s\n\n", PACKAGE, VERSION);
  fprintf (stderr, "Usage: %s [options] file1 [file2] [file3] ...\n\n", PACKAGE);
  fprintf (stderr,
           " ## Input/Output Options ##\n"
           " -V            Report program version\n"
           " -h            Show this usage message\n"
           " -v            Be more verbose, multiple flags can be used\n"
           " -s            Print a basic summary after reading all input files\n"
           " -ts time      Limit miniSEED input to records that start after time\n"
           " -te time      Limit miniSEED input to records that end before time\n"
           "                 time format: 'YYYY[,DDD,HH,MM,SS,FFFFFF]' delimiters: [,:.]\n"
           " -I[ms]        Force input to be considered (m)iniSEED or (s)AC\n"
           "                 By default the input format is autodetected\n"
           " -MSEED        Write output as miniSEED instead of default SAC\n"
           " -Mr bytes     Specify record length in bytes, default is autodetection\n"
           " -Me encoding  Specify encoding format of data samples for input data\n"
           " -MR bytes     Specify record length for output miniSEED, default is 4096\n"
           /*	   " -S            Include Blockette 1001 (hi-res sample rate) in output miniSEED\n"*/
           " -ME encoding  Encoding format for SEED output data, default is 4 (floats)\n"
           " -Sf format    Specify SAC output format (default is 2:binary)\n"
           "                 1=alpha, 2=binary (host byte order),\n"
           "                 3=binary (little-endian), 4=binary (big-endian)\n"
           " -m metafile   File containing station metadata for SAC output\n"
           " -msi          Convert component inclination/dip from SEED to SAC convention\n"
           " -c channel    Set channel/component name for processed output\n"
           " -o outfile    Specify output file instead of generating file names\n"
           " -od outdir    Specify output directory for generated file names\n"
           " -plog file    Write detailed processing log to specified file\n"
           "\n"
           " ## Convolution Options ##\n"
           " -CR respfile[:#:#] Specify SEED RESP file/dir for convolution\n"
           " -DR respfile[:#:#] Specify SEED RESP file/dir for deconvolution\n"
           " -CS pzfile    Specify poles and zeros file for convolution\n"
           " -DS pzfile    Specify poles and zeros file for deconvolution\n"
           " -FL freqs     Specify frequency limits for deconvolution operations\n"
           "                 Frequencies are specify as: 'f1/f2/f3/f4'\n"
           " -FLa dBdown   Automatically determine frequency limits for specified dBdown\n"
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
           " -RM           Remove mean from time series\n"
           " -SC factor    Scale the data samples by a specified factor\n"
           " -SI factor    Scale the data samples by inverse of specified factor\n"
           " -DEC factor   Decimate time series by specified factor\n"
           " -TAP width[:type]  Apply symmetric taper to time series\n"
           " -POLYM c1,c2,...   Apply a Maclaurin type polynomial with given coefficients\n"
           " -ENV          Calculate envelope of time series\n"
           " -DTRIM        Trim time series to latest start and earliest end\n"
           " -ROTATE E[/1],N[/2],Z[/3]:azimuth[,incidence]\n"
           "                 Rotate component sets, 2-D or 3-D\n"
           " -STATS        Add simple stats to processing log, verbose to print\n"
           "\n"
           " file#         File of input miniSEED or SAC\n"
           "\n");
} /* End of usage() */
