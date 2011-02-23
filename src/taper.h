/* Routines for series taper */

#ifndef TAPER_H
#define TAPER_H 1

#ifdef __cplusplus
extern "C" {
#endif

extern int taper (double *data, int npts, double width, int type);

/* Taper types */
#define TAPER_HANNING  1
#define TAPER_HAMMING  2
#define TAPER_COSINE   3

#ifdef __cplusplus
}
#endif

#endif /* TAPER_H */
