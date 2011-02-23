/* Routines for handling complex numbers */

#ifndef COMPLEX_H
#define COMPLEX_H 1

#ifdef __cplusplus
extern "C" {
#endif

typedef struct complexf_s {
  double real;
  double imag;
} complexf;

double cmplximag (complexf c);
double cmplxreal (complexf c);
complexf cmplxadd (complexf c1, complexf c2);
complexf cmplxcj (complexf c);
complexf cmplxmul (complexf c1, complexf c2);
complexf flttocmplx (double d1, double d2);
complexf cmplxsub(complexf c1, complexf c2);
double cmplxabs (complexf c);
double cmplxang (complexf c);
complexf cmplxsqrt (complexf c);
complexf cmplxdiv (complexf c1, complexf c2);
complexf cmplxlog (complexf c);
complexf cmplxexp (complexf c);
complexf cmplxpow (complexf c, double d);
complexf cmplxneg (complexf c);

#ifdef __cplusplus
}
#endif

#endif /* COMPLEX_H */
