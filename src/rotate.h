/* Routines for series rotation */

#ifndef ROTATE_H
#define ROTATE_H 1

#ifdef __cplusplus
extern "C" {
#endif

void rotate2 ( double *n, double *e, long int npts, double azimuth,
	       double *r, double *t );
void rotate3 ( double *z, double *n, double *e, long int npts,
	       double azimiuth, double incidence,
	       double *l, double *q, double *t );

#ifdef __cplusplus
}
#endif

#endif /* ROTATE_H */
