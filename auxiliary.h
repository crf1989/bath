#ifndef AUXILIARY_H
#define AUXILIARY_H 1

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <complex.h>

const double PI = 3.141592653589793;
const double BOLTZMANN_CONSTANT = 8.6173324e-5;
const double u = 5.7883818066e-5; /* bohr magneton, eV/T */
const double g = 2;		  /* lande factor */

int pow2();
double square (double);
double cabs2 (double complex);
double drand ();
double urand ();
void box_muller (double complex*, double);
void tee (FILE*, char*, ...);
double fermi (double, double, double);
double complex sfg_solve (double, double);

FILE* flog = NULL;

/********************************************************/

int pow2 (int n)
{
  assert ((0 <= n) && (n < 32));
  int result = 1;
  while (n--)
    result <<= 1;
  return result;
}
  
double square (double x)
{
  return x*x;
}

double cabs2 (double complex x)
{
  return square (cabs(x));
}

double drand ()
{
  return (random() + 1.0)/(RAND_MAX+1.0);
}

double urand ()
{
  return 2 * (drand() - 0.5);
}

void box_muller (double complex* x, double sigma)
{
  double r1 = drand ();
  double r2 = drand ();
  *x = sigma* (sqrt(-2*log(r2))*cos(2*PI*r1) +
	       sqrt(-2*log(r2))*sin(2*PI*r1)*I);
}

void tee (FILE* fp, char* FORMAT, ...)
{
  va_list argp;
  va_start (argp, FORMAT);
  vprintf (FORMAT, argp);
  va_end (argp);

  va_start (argp, FORMAT);
  vfprintf (fp, FORMAT, argp);
  va_end (argp);
}

double fermi (double e, double u, double T)
{
  assert (T >= 0);
  if (T == 0)
    if (e == u)
      return 0.5;
    else 
      return (e > u ? 0 : 1);
  else
    return 1/(exp((e-u)/(BOLTZMANN_CONSTANT*T))+1);
}

double complex sfg_solve (double e, double t)
{
  double complex lambda1 = 0.5*(e/t-csqrt(square(e/t)-4));
  double complex lambda2 = 0.5*(e/t+csqrt(square(e/t)-4));
  if (cabs(lambda1) <= 1)
    return lambda1*t;
  else
    return lambda2*t;
}



#endif /* AUXILIARY_H */
