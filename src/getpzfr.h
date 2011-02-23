
#ifndef GETPZFR_H
#define GETPZFR_H 1

#ifdef __cplusplus
extern "C" {
#endif

int getpzfr (int nfreq, double delfrq, double xreal[], double ximag[],
	     char *pzfilename);

#ifdef __cplusplus
}
#endif

#endif /* GETPZFR_H */
