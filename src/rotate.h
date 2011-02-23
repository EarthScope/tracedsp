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

/* Define data types */
#define SAMPLE float
#define REAL double

void mt_rot2 ( SAMPLE *n, SAMPLE *e, long lth, REAL angle,
	       SAMPLE *r, SAMPLE *t );
void mt_rot3 ( SAMPLE *z, SAMPLE *n, SAMPLE *e, long lth,
	       REAL azim, REAL inci, int type,
	       SAMPLE *l, SAMPLE *q, SAMPLE *t );

#ifdef __cplusplus
}
#endif

#endif /* ROTATE_H */
