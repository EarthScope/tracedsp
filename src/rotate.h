/* Routines for series rotation */

#ifndef ROTATE_H
#define ROTATE_H 1

#ifdef __cplusplus
extern "C" {
#endif

/* Rotation types for 3-D rotation */
#define ROT_ZNE_TO_LQT 0
#define ROT_ZNE_TO_UVW 1
#define ROT_UVW_TO_ZNE 2

void rotate2 ( double *n, double *e, long lth, double angle,
	       double *r, double *t );
void rotate3 ( double *z, double *n, double *e, long lth,
	       double azim, double inci, int type,
	       double *l, double *q, double *t );

#ifdef __cplusplus
}
#endif

#endif /* ROTATE_H */
