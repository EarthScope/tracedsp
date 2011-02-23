/* Routines to whiten and dewhiten data */

#ifndef WHITEN_H
#define WHITEN_H 1

#ifdef __cplusplus
extern "C" {
#endif


int whiten (float data[], int npts, int *order, double coef[]);
int dewhiten (float data[], int npts, int order, double coef[]);
double LPCautocorr (float data[], int npts, int order,
                    double lpc[], double ref[] );
int applypef (float data[], int npts, double coef[], int nc, float result[]);
int applypf (float data[], int npts, double coef[], int nc, float result[]);


#ifdef __cplusplus
}
#endif

#endif /* WHITEN_H */
