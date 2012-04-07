/* Routines to whiten and dewhiten data */

#ifndef WHITEN_H
#define WHITEN_H 1

#ifdef __cplusplus
extern "C" {
#endif


int whiten (double data[], int npts, int *order, double coef[]);
int dewhiten (double data[], int npts, int order, double coef[]);
double LPCautocorr (double data[], int npts, int order,
                    double lpc[], double ref[] );
int applypef (double data[], int npts, double coef[], int nc, double result[]);
int applypf (double data[], int npts, double coef[], int nc, double result[]);


#ifdef __cplusplus
}
#endif

#endif /* WHITEN_H */
