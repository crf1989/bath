#ifndef SYSTEM_H
#define SYSTEM_H 1

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>

#include "auxiliary.h"
#include "parameter.h"
#include "chain.h"
#include "bath.h"
#include "spin.h"

chain* center1 = NULL;
chain* center2 = NULL;
bath* left1 = NULL;
bath* left2 = NULL;
bath* right1 = NULL;
bath* right2 = NULL;
vector* spin = NULL;

void init_system ();
void rotate_electron (int i, vector s);
/*******************************************************/

void init_system ()
{
  left1 = bath_alloc (NOISE_NUMBER, TIME_INTERVAL, LEFT_LEAD_TEMPERATURE,
		      LEFT_UP_LEAD_CHEMICAL_POTENTIAL, HOPPING_ENERGY);
  assert (left1 != NULL);
  left2 = bath_alloc (NOISE_NUMBER, TIME_INTERVAL, LEFT_LEAD_TEMPERATURE,
		      LEFT_DOWN_LEAD_CHEMICAL_POTENTIAL, HOPPING_ENERGY);
  assert (left2 != NULL);

  right1 = bath_alloc (NOISE_NUMBER, TIME_INTERVAL, RIGHT_LEAD_TEMPERATURE,
		       RIGHT_UP_LEAD_CHEMICAL_POTENTIAL, HOPPING_ENERGY);
  assert (right1 != NULL);
  right2 = bath_alloc (NOISE_NUMBER, TIME_INTERVAL, RIGHT_LEAD_TEMPERATURE,
		       RIGHT_DOWN_LEAED_CHEMICAL_POTENTIAL, HOPPING_ENERGY);
  assert (right2 != NULL);

  center1 = chain_alloc (CENTER_SITE_NUMBER, SELF_ENERGY_NUMBER,
			 HOPPING_ENERGY, TIME_INTERVAL);
  assert (center1 != NULL);
  center2 = chain_alloc (CENTER_SITE_NUMBER, SELF_ENERGY_NUMBER,
			 HOPPING_ENERGY, TIME_INTERVAL);
  assert (center2 != NULL);

  spin = spin_alloc (CENTER_SITE_NUMBER+1);

  generate_bath (left1);
  generate_bath (left2);
  generate_bath (right1);
  generate_bath (right2);

  if (WRITE_SPECTRAL)
    {
      assert (system ("mkdir -p check/spectral") != -1);
      fflush (stdout);
      write_spectral (left1,
		      "#This is the spectral for left up lead\n"
		      "#omega\tspectral\n", "check/spectral/left_up_spectral.dat");
      write_spectral (left2,
		      "#This is the spectral for left down lead\n"
		      "#omega\tspectral\n", "check/spectral/left_down_spectral.dat");

      write_spectral (right1,
		      "#This is the spectral for right up lead\n"
		      "#omega\tspectral\n", "check/spectral/right_up_spectral.dat");
      write_spectral (right2,
		      "#This is the spectral for right down lead\n"
		      "#omega\tspectral\n", "check/spectral/right_down_spectral.dat");

    }
  if (WRITE_FRIC)
    {
      assert (system ("mkdir -p check/fric") != -1);
      fflush (stdout);
      write_fric (left1,
		  "#This is the self energy for left up lead\n"
		  "#time\treal\timag\n", "check/fric/left_up_fric.dat");
      write_fric (left2,
		  "#This is the self energy for left down lead\n"
		  "#time\treal\timag\n", "check/fric/left_down_fric.dat");

      write_fric (right1,
		  "#This is the self energy for right up lead\n"
		  "#time\treal\timag\n", "check/fric/right_fric.dat");
      write_fric (right2,
		  "#This is the self energy for right down lead\n"
		  "#time\treal\timag\n", "check/fric/right_down_fric.dat");
    }
  if (WRITE_NOISE)
    {
      assert (system ("mkdir -p check/noise") != -1);
      write_noise (left1,
		   "#This is the noise in omega domain for left up lead\n"
		   "#omega\treal\timag\n", "check/noise/left_noise.dat");
      write_noise (left2,
		   "#This is the noise in omega domain for left down lead\n"
		   "#omega\treal\timag\n", "check/noise/left_noise.dat");

      write_noise (right1,
		   "#This is the noise in omega domain for right up lead\n"
		   "#omega\treal\timag\n", "check/noise/right_noise.dat");
      write_noise (right2,
		   "#This is the noise in omega domain for right down lead\n"
		   "#omega\treal\timag\n", "check/noise/right_down_noise.dat");
    }
  if (WRITE_ETA)
    {
      assert (system ("mkdir -p check/eta") != -1);
      write_eta (left1,
		 "#This is eta (noise in time domain) for left up lead\n"
		 "#time\treal\timag\n", "check/eta/left_up_eta.dat");
      write_eta (left2,
		 "#This is eta (noise in time domain) for left down lead\n"
		 "#time\treal\timag\n", "check/eta/left_down_eta.dat");
      write_eta (right1,
		 "#This is eta (noise in time domain) for right up lead\n"
		 "#time\treal\timag\n", "check/eta/right_eta.dat");
      write_eta (right2,
		 "#This is eta (noise in time domain) for right down lead\n"
		 "#time\treal\timag\n", "check/eta/right_down_eta.dat");
    }
  if (WRITE_CORRELATION)
    {
      assert (system ("mkdir -p check/correlation") != -1);
      write_correlation (left1,
			 "#This is eta correlation for left up lead\n"
			 "#time\treal\timag\n", "check/correlation/left_up_correlation.dat");
      write_correlation (left2,
			 "#This is eta correlation for left down lead\n"
			 "#time\treal\timag\n", "check/correlation/left_down_correlation.dat");
      write_correlation (right1,
			 "#This is eta correlation for right up lead\n"
			 "#time\treal\timag\n", "check/correlation/right_up_correlation.dat");
      write_correlation (right2,
			 "#This is eta correlation for right down lead\n"
			 "#time\treal\timag\n", "check/correlation/right_down_correlation.dat");
    }
}


void rotate_electron (int i, vector s)
{
  double complex sigma[2][2];
  
  double phi = sqrt (square(s.s[0]) + square(s.s[1]) + square(s.s[2]));
  
  /* if |s| is zero, do nothing. */
  if (fabs (phi) < 1e-8)
    return;
  
  /* \vec{s}.\vec{simga}/|s|, here simga is Pauli matrice. */
  for (int k = 0; k < 3; ++k)
	s.s[k] /= phi;
  
  double cth = cos (phi*TIME_INTERVAL);
  double sth = sin (phi*TIME_INTERVAL);
  
  sigma[0][0] = cth + I*s.s[2]*sth;
  sigma[0][1] = I*s.s[0]*sth + s.s[1]*sth;
  sigma[1][0] = I*s.s[0]*sth - s.s[1]*sth;
  sigma[1][1] = cth - I*s.s[2]*sth;

  /* matrice multiplication */
  double complex tmp1 = center1->psi[i];
  double complex tmp2 = center2->psi[i];
  
  center1->psi[i] = sigma[0][0]*tmp1 + sigma[0][1]*tmp2;
  center2->psi[i] = sigma[1][0]*tmp1 + sigma[1][1]*tmp2;
}

/* void md_step () */
/* { */
/*   double complex left1_fric = 0; */
/*   double complex left2_fric = 0; */
/*   double complex right1_fric = 0; */
/*   double complex right2_fric = 0; */
/*   int m = i%center1->n; */
/*   int n = (i+1)%center1->n; */

/*   for (int j = 0; j < center1->n; ++j) */
/*     { */
/*       left1_fric += left1->fric[left1->size-center1->n+j] * */
/* 	center1->cl[(n+j)%center1->n]; */
/*       left2_fric += left2->fric[left2->size-center2->n+j] * */
/* 	center2->cl[(n+j)%center2->n]; */
	  
/*       right1_fric += right1->fric[right1->size-center1->n+j] * */
/* 	center1->cr[(n+j)%center1->n]; */
/*       right2_fric += right2->fric[right2->size-center2->n+j] * */
/* 	center2->cr[(n+j)%center2->n]; */
/*     } */

/*   left1_fric *= dt; */
/*   left2_fric *= dt; */
/*   right1_fric *= dt; */
/*   right2_fric *= dt; */

/*   center1->psi[0] += -I*(left1->eta[i]+left1_fric)*dt; */
/*   center2->psi[0] += -I*(left2->eta[i]+left2_fric)*dt; */
/*   center1->psi[center1->size-1] += -I*(right1->eta[i]+right1_fric)*dt; */
/*   center2->psi[center2->size-1] += -I*(right2->eta[i]+right2_fric)*dt; */

/*   vector H; */
/*   vector s; */
/*   for (int j = 0; j < CENTER_SITE_NUMBER; ++j) */
/*     { */
/*       for (int k = 0; k < 3; ++k) */
/* 	s.s[k] = 0.5*JH*(spin[j].s[k]+spin[j+1].s[k]); */
/*       rotate_electron (j, s); */
/*     } */
/*   for (int j = 0; j < CENTER_SITE_NUMBER+1; ++j) */
/*     { */
/*       for (int k = 0; k < 3; ++k) */
/* 	if (j == 0) */
/* 	  H.s[k] = J1*spin[j+1].s[k]; */
/* 	else if (j == CENTER_SITE_NUMBER) */
/* 	  H.s[k] = J1*spin[j-1].s[k]; */
/* 	else */
/* 	  H.s[k] = J1*(spin[j-1].s[k]+spin[j-1].s[k]); */
	  
/*       if (j == 0) */
/* 	{ */
/* 	  H.s[0] += JH*creal(conj(center1->psi[j])*center2->psi[j]); */
/* 	  H.s[1] += JH*cimag(conj(center1->psi[j])*center2->psi[j]); */
/* 	  H.s[2] += JH*0.5*(cabs2(center1->psi[j])-cabs2(center2->psi[j])); */
/* 	} */
/*       else if (j == CENTER_SITE_NUMBER) */
/* 	{ */
/* 	  H.s[0] += JH*creal(conj(center1->psi[j-1])*center2->psi[j-1]); */
/* 	  H.s[1] += JH*cimag(conj(center1->psi[j-1])*center2->psi[j-1]); */
/* 	  H.s[2] += JH*0.5*(cabs2(center1->psi[j-1])-cabs2(center2->psi[j-1])); */
/* 	} */
/*       else */
/* 	{ */
/* 	  H.s[0] += JH*creal(conj(center1->psi[j-1])*center2->psi[j-1])+ */
/* 	    JH*creal(conj(center1->psi[j])*center2->psi[j]); */
/* 	  H.s[1] += JH*cimag(conj(center1->psi[j-1])*center2->psi[j-1])+ */
/* 	    JH*cimag(conj(center1->psi[j])*center2->psi[j]); */
/* 	  H.s[2] += JH*0.5*(cabs2(center1->psi[j-1])-cabs2(center2->psi[j-1]))+ */
/* 	    JH*0.5*(cabs2(center1->psi[j])-cabs2(center2->psi[j])); */
/* 	} */
/*       rotate_spin (&spin[j], H); */
/*     } */

/*   crank_nicolson (center1); */
/*   crank_nicolson (center2); */

/*   center1->cl[n] = center1->psi[0]; */
/*   center2->cl[n] = center2->psi[0]; */
/*   center1->cr[n] = center1->psi[center1->size-1]; */
/*   center2->cr[n] = center2->psi[center2->size-1]; */
/* } */
  
#endif	/* SYSTEM_H */