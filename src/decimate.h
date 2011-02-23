/* Routine for decimation */

#ifndef DECIMATE_H
#define DECIMATE_H 1

#ifdef __cplusplus
extern "C" {
#endif

int decimate (double *data, int npts, int factor,
	      double *fir, int firc, int firsym);

#ifdef __cplusplus
}
#endif

#endif /* DECIMATE_H */
