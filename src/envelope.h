/* Routines for series envelope and Hilbert transform */

#ifndef ENVELOPE_H
#define ENVELOPE_H 1

#ifdef __cplusplus
extern "C" {
#endif

extern int envelope (double *data, int npts);
extern int hilbert (double *x, double *y, int npts);

#ifdef __cplusplus
}
#endif

#endif /* ENVELOPE_H */
