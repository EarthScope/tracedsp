/* Routines for handling complex numbers */

#ifndef COMPLEX_H
#define COMPLEX_H 1

#ifdef __cplusplus
extern "C" {
#endif

#define PI 3.1415926535897932384626433832795

typedef struct complexd_s {
  double real;
  double imag;
} complexd;

double cmplximag (complexd c);
double cmplxreal (complexd c);
complexd cmplxadd (complexd c1, complexd c2);
complexd cmplxconj (complexd c);
complexd cmplxmul (complexd c1, complexd c2);
complexd cmplxdmul (double d, complexd c);
complexd dbltocmplx (double d1, double d2);
complexd cmplxsub(complexd c1, complexd c2);
double cmplxabs (complexd c);
double cmplxang (complexd c);
complexd cmplxsqrt (complexd c);
complexd cmplxdiv (complexd c1, complexd c2);
complexd cmplxlog (complexd c);
complexd cmplxexp (complexd c);
complexd cmplxpow (complexd c, double d);
complexd cmplxneg (complexd c);

#ifdef __cplusplus
}
#endif

#endif /* COMPLEX_H */
