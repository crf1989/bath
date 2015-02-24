#ifndef SPIN_H 
#define SPIN_H 1

#include "parameter.h"
#include "auxiliary.h"

typedef struct
{
  double s[3];
} vector;

void genvec (vector* a);
vector* spin_alloc (int n);
void spin_free (vector* p);
void rotate_spin (vector* p, vector H);

/******************************************************/

void genvec (vector* a)
{
  /* theta is in [0, PI] */
  double theta = asin (urand()) + PI/2;
  /* phi is in [0, 2PI] */
  double phi = 2*PI * drand ();

  a->s[0] = S * sin (theta) * cos (phi);
  a->s[1] = S * sin (theta) * sin (phi);
  a->s[2] = S * cos (theta);
}

vector* spin_alloc (int n)
{
  vector* p = malloc (n * sizeof(vector));
  for (int i = 0; i < n; ++i)
    genvec (&p[i]);
  return p;
}

void spin_free (vector* p)
{
  free (p);
}

void rotate_spin (vector* p, vector H) 
{
  double R[3][3];
  /* angular velocity */
  double omega = sqrt (square(H.s[0]) + square(H.s[1]) + square(H.s[2]));
  if (fabs (omega <1e-8))
    return;
  /* unit rotation axis */
  H.s[0] /= omega; H.s[1] /= omega; H.s[2] /= omega;
  double cth = cos (-TIME_INTERVAL*omega); /* cos (theta) */
  double sth = sin (-TIME_INTERVAL*omega); /* sin (theta) */
  
  R[0][0] = cth + H.s[0]*H.s[0]*(1-cth);
  R[0][1] = H.s[0]*H.s[1]*(1-cth) - H.s[2]*sth;
  R[0][2] = H.s[0]*H.s[2]*(1-cth) + H.s[1]*sth;
  
  R[1][0] = H.s[1]*H.s[0]*(1-cth) + H.s[2]*sth;
  R[1][1] = cth + H.s[1]*H.s[1]*(1-cth);
  R[1][2] = H.s[1]*H.s[2]*(1-cth) - H.s[0]*sth;
  
  R[2][0] = H.s[2]*H.s[0]*(1-cth) - H.s[1]*sth;
  R[2][1] = H.s[2]*H.s[1]*(1-cth) + H.s[0]*sth;
  R[2][2] = cth + H.s[2]*H.s[2]*(1-cth);
  
  double stx = p->s[0];
  double sty = p->s[1];
  double stz = p->s[2];
  
  /* matrix multiplication */
  for (int k = 0; k < 3; ++k)
    p->s[k] = R[k][0]*stx + R[k][1]*sty + R[k][2]*stz;
}


#endif	/* SPIN_H */
