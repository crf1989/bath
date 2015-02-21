#ifndef CHAIN_H
#define CHAIN_H 1

#include <stdio.h>
#include <complex.h>
#include <assert.h>

typedef struct
{
  int size;
  int n;			/* number for self energy */
  double h;			/* hopping energy */
  double dt;			/* time interval */
  double complex* psi;
  double complex* cl;
  double complex* cr;
} chain;

chain* chain_alloc (int, int, double, double);
void chain_clear (chain*);
double get_norm (chain*);
void crank_nicolson (chain*);

/**********************************************************/

chain* chain_alloc (int size, int n, double h, double dt)
{
  chain* p = malloc (sizeof(chain));
  assert (size > 1);
  p->size = size;
  p->n = n;
  p->h = h;
  p->dt = dt;
  p->psi = calloc (size, sizeof(double complex));
  assert (p->psi != NULL);
  p->cl = calloc (n, sizeof(double complex));
  assert (p->cl != NULL);
  p->cr = calloc (n, sizeof(double complex));
  assert (p->cr != NULL);
  return p;
}

void chain_clear (chain* p)
{
  free (p->psi);
  free (p->cl);
  free (p->cr);
  free (p);
}

double get_norm (chain* p)
{
  double sum = 0;
  for (int i = 0; i < p->size; ++i)
    sum += cabs2 (p->psi[i]);
  return sum;
}

void crank_nicolson (chain* p)
{
  double t = p->h;
  double dt = p->dt;
  double complex* a = malloc (p->size * sizeof(double complex));
  double complex* b = malloc (p->size * sizeof(double complex));
  double complex* c = malloc (p->size * sizeof(double complex));
  double complex* d = malloc (p->size * sizeof(double complex));
  
  d[0] = p->psi[0] - I*((-t)*p->psi[1])*dt/2;
  d[p->size-1] = p->psi[p->size-1] - I*((-t)*p->psi[p->size-2])*dt/2;
  for (int i = 1; i < p->size-1; ++i)
    d[i] = p->psi[i] - I*((-t)*(p->psi[i-1]+p->psi[i+1]))*dt/2;
  
  a[0] = 0; c[p->size-1] = 0;
  for (int i = 1; i < p->size; ++i)
    a[i] = I*(-t)*dt/2;
  for (int i = 0; i < p->size; ++i)
    b[i] = 1;
  for (int i = 0; i < p->size-1; ++i)
    c[i] = I*(-t)*dt/2;

  c[0] /= b[0];
  for (int i = 1; i < p->size-1; ++i)
    c[i] /= (b[i]-a[i]*c[i-1]);
  d[0] /= b[0];
  for (int i = 1; i < p->size; ++i)
    d[i] = (d[i]-a[i]*d[i-1])/(b[i]-a[i]*c[i-1]);

  p->psi[p->size-1] = d[p->size-1];
  for (int i = p->size-2; i >= 0; --i)
    p->psi[i] = d[i] - c[i]*p->psi[i+1];

  free (a);
  free (b);
  free (c);
  free (d);
}

  
  
#endif /* CHAIN_H */
