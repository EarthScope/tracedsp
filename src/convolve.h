/* Routines for convolution/deconvolution */

#ifndef CONVOLVE_H
#define CONVOLVE_H 1

#ifdef __cplusplus
extern "C" {
#endif

int convolve (double data[], int npts, double delta, int nfreqs, int nfft,
	      double *creal, double *cimag, double *dreal, double *dimag,
	      double taperfreq[], int *prewhiten,
	      int verbose);

int calcfr_sac (int nfreqs, double delfreq, char *sacpzfilename,
		double **xreal, double **ximag, int verbose);

int calcfr_resp (int nfreqs, double delfreq, char *net, char *sta,
		 char *loc, char *chan, int startstage, int stopstage,
		 char *units, time_t resptime, int usedelay,
		 char *respfilename, int totalsensflag, 
		 double **xreal, double **ximag, int verbose);

int convolve_sac (double data[], int npts, double delta, double taperfreq[],
		  int *prewhiten, int deconvflag, char *sacpzfilename,
		  int verbose);
  
int convolve_resp (double data[], int npts, double delta,
		   char *net, char *sta, char *loc, char *chan,
		   int startstage, int stopstage,
		   char *units, time_t resptime, int usedelay,
		   double taperfreq[], int *prewhiten,
		   int deconvflag, char *respfilename,
		   int totalsensflag, int verbose);

int next2 (int value);
double spectraltaper (double freq, double fqh, double fql);
int findtaper (double *taperfreq, double *xreal, double *ximag,
	       int nfreqs, double delfreq, double lcdBdown, double ucdBdown);
int fft (double real[], double imag[], int nfreq, int direction);


#ifdef __cplusplus
}
#endif

#endif /* CONVOLVE_H */
