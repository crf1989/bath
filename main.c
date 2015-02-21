#include <stdio.h>
#include <stdlib.h>

#include "init.h"

int main () 
{
  init ();
  /* set_output (); */

  /* double dt = TIME_INTERVAL; */
  /* tee (flog, "Simulation starts ... "); */
  /* fflush(stdout); */

  double mean_left1_number_current = 0;
  double mean_left2_number_current = 0;
  double mean_right1_number_current = 0;
  double mean_right2_number_current = 0;
  double mean_left1_energy_current = 0;
  double mean_left2_energy_current = 0;
  double mean_right1_energy_current = 0;
  double mean_right2_energy_current = 0;

/*   for (int i = 0; i < NOISE_NUMBER; ++i) */
/*     { */
/*       double complex left1_fric = 0; */
/*       double complex left2_fric = 0; */
/*       double complex right1_fric = 0; */
/*       double complex right2_fric = 0; */
/*       int m = i%center1->n; */
/*       int n = (i+1)%center1->n; */

/*       for (int j = 0; j < center1->n; ++j) */
/*       	{ */
/*       	  left1_fric += left1->fric[left1->size-center1->n+j] * */
/*       	    center1->cl[(n+j)%center1->n]; */
/* 	  left2_fric += left2->fric[left2->size-center2->n+j] * */
/*       	    center2->cl[(n+j)%center2->n]; */
	  
/*       	  right1_fric += right1->fric[right1->size-center1->n+j] * */
/*       	    center1->cr[(n+j)%center1->n]; */
/* 	  right2_fric += right2->fric[right2->size-center2->n+j] * */
/*       	    center2->cr[(n+j)%center2->n]; */
/*       	} */

/*       left1_fric *= dt; */
/*       left2_fric *= dt; */
/*       right1_fric *= dt; */
/*       right2_fric *= dt; */

/*       center1->psi[0] += -I*(left1->eta[i]+left1_fric)*dt; */
/*       center2->psi[0] += -I*(left2->eta[i]+left2_fric)*dt; */
/*       center1->psi[center1->size-1] += -I*(right1->eta[i]+right1_fric)*dt; */
/*       center2->psi[center2->size-1] += -I*(right2->eta[i]+right2_fric)*dt; */

/* #ifdef ADD_SPIN */
/*       vector H; */
/*       vector s; */
/*       for (int j = 0; j < CENTER_SITE_NUMBER; ++j) */
/* 	{ */
/* 	  for (int k = 0; k < 3; ++k) */
/* 	    s.s[k] = 0.5*JH*(spin[j].s[k]+spin[j+1].s[k]); */
/* 	  rotate_electron (j, s); */
/* 	} */
/*       for (int j = 0; j < CENTER_SITE_NUMBER+1; ++j) */
/* 	{ */
/* 	  for (int k = 0; k < 3; ++k) */
/* 	    if (j == 0) */
/* 	      H.s[k] = J1*spin[j+1].s[k]; */
/* 	    else if (j == CENTER_SITE_NUMBER) */
/* 	      H.s[k] = J1*spin[j-1].s[k]; */
/* 	    else */
/* 	      H.s[k] = J1*(spin[j-1].s[k]+spin[j-1].s[k]); */
	  
/* 	  if (j == 0) */
/* 	    { */
/* 	      H.s[0] += JH*creal(conj(center1->psi[j])*center2->psi[j]); */
/* 	      H.s[1] += JH*cimag(conj(center1->psi[j])*center2->psi[j]); */
/* 	      H.s[2] += JH*0.5*(cabs2(center1->psi[j])-cabs2(center2->psi[j])); */
/* 	    } */
/* 	  else if (j == CENTER_SITE_NUMBER) */
/* 	    { */
/* 	      H.s[0] += JH*creal(conj(center1->psi[j-1])*center2->psi[j-1]); */
/* 	      H.s[1] += JH*cimag(conj(center1->psi[j-1])*center2->psi[j-1]); */
/* 	      H.s[2] += JH*0.5*(cabs2(center1->psi[j-1])-cabs2(center2->psi[j-1])); */
/* 	    } */
/* 	  else */
/* 	    { */
/* 	      H.s[0] += JH*creal(conj(center1->psi[j-1])*center2->psi[j-1])+ */
/* 		JH*creal(conj(center1->psi[j])*center2->psi[j]); */
/* 	      H.s[1] += JH*cimag(conj(center1->psi[j-1])*center2->psi[j-1])+ */
/* 		JH*cimag(conj(center1->psi[j])*center2->psi[j]); */
/* 	      H.s[2] += JH*0.5*(cabs2(center1->psi[j-1])-cabs2(center2->psi[j-1]))+ */
/* 		JH*0.5*(cabs2(center1->psi[j])-cabs2(center2->psi[j])); */
/* 	    } */
/* 	  rotate_spin (&spin[j], H); */
/* 	} */
/* #endif */

/*       crank_nicolson (center1); */
/*       crank_nicolson (center2); */

/*       center1->cl[n] = center1->psi[0]; */
/*       center2->cl[n] = center2->psi[0]; */
/*       center1->cr[n] = center1->psi[center1->size-1]; */
/*       center2->cr[n] = center2->psi[center2->size-1]; */

/*       if (i%stride != 0) */
/*   	continue; */
      
/*       if (WRITE_PSI) */
/*   	{ */
/*   	  for (int j = 0; j < center1->size; ++j) */
/*   	    fprintf (fpsi[j], "%.8g\t%.8g\t%.8g\t%.8g\t%.8g\n", */
/*   		     (i+1)*dt,  */
/* 		     creal (center1->psi[j]), cimag (center1->psi[j]), */
/* 		     creal(center2->psi[j]), cimag(center2->psi[j])); */
/*   	  } */
/*       if (WRITE_NORM) */
/*   	fprintf (fnorm, "%.8g\t%.8g\n", (i+1)*dt,  */
/* 		 get_norm (center1)+get_norm(center2)); */

/*       { */
/*   	double left1_number_current = */
/*   	  2*cimag (conj(center1->cl[m])*(left1->eta[i]+left1_fric)); */
/*   	double left2_number_current = */
/*   	  2*cimag (conj(center2->cl[m])*(left2->eta[i]+left2_fric)); */
/*   	double right1_number_current = */
/*   	  2*cimag (conj(center1->cr[m])*(right1->eta[i]+right1_fric)); */
/*   	double right2_number_current = */
/*   	  2*cimag (conj(center2->cr[m])*(right2->eta[i]+right2_fric)); */
      
/*   	mean_left1_number_current += left1_number_current; */
/*   	mean_left2_number_current += left2_number_current; */
/*   	mean_right1_number_current += right1_number_current; */
/*   	mean_right2_number_current += right2_number_current; */
	
/*   	fprintf (fnc, "%g\t%g\t%g\t%g\t%g\n", i*dt, */
/*   		 left1_number_current, left2_number_current, */
/*   		 right1_number_current, right2_number_current); */
/*   	fprintf (fmnc, "%g\t%g\t%g\t%g\t%g\n", i*dt, */
/*   		 mean_left1_number_current/(i/stride), */
/* 		 mean_left2_number_current/(i/stride), */
/*   		 mean_right1_number_current/(i/stride), */
/*   		 mean_right2_number_current/(i/stride)); */
/*       } */

/*       { */
/*   	double left1_energy_current = */
/*   	  -2*creal (conj(left1->eta[i]+left1_fric)* */
/*   			(center1->cl[n]-center1->cl[m])/dt); */
/* 	double left2_energy_current = */
/*   	  -2*creal (conj(left2->eta[i]+left2_fric)* */
/*   			(center2->cl[n]-center2->cl[m])/dt); */
/*   	double right1_energy_current = */
/*   	  -2*creal (conj(right1->eta[i]+right1_fric)* */
/*   			(center1->cr[n]-center1->cr[m])/dt); */
/*   	double right2_energy_current = */
/*   	  -2*creal (conj(right2->eta[i]+right2_fric)* */
/*   			(center2->cr[n]-center2->cr[m])/dt); */

/*   	mean_left1_energy_current += left1_energy_current; */
/* 	mean_left2_energy_current += left2_energy_current; */
/*   	mean_right1_energy_current += right1_energy_current; */
/*   	mean_right2_energy_current += right2_energy_current; */

/*   	fprintf (fec, "%g\t%g\t%g\t%g\t%g\n", i*dt, */
/*   		 left1_energy_current, left2_energy_current, */
/* 		 right1_energy_current, right2_energy_current); */
/*   	fprintf (fmec, "%g\t%g\t%g\t%g\t%g\n", i*dt, */
/*   		 mean_left1_energy_current/(i/stride), */
/* 		 mean_left2_energy_current/(i/stride), */
/*   		 mean_right1_energy_current/(i/stride), */
/* 		 mean_right2_energy_current/(i/stride)); */
/*       } */
/*   } */
  
/*   tee (flog, "done.\n\n"); */

/*   tee (flog, "mean particle number current for left lead spin up:\t%g\n", */
/*        mean_left1_number_current/(NOISE_NUMBER/stride)); */
/*   tee (flog, "mean particle number current for left lead spin down:\t%g\n", */
/*        mean_left2_number_current/(NOISE_NUMBER/stride)); */
/*   tee (flog, "mean particle number current for left lead:\t\t%g\n", */
/*        (mean_left1_number_current+mean_left2_number_current)/ */
/*        (NOISE_NUMBER/stride)); */
/*   tee (flog, "\n"); */

/*   tee (flog, "mean particle number current for right lead spin up:\t%g\n", */
/*        mean_right1_number_current/(NOISE_NUMBER/stride)); */
/*   tee (flog, "mean particle number current for right lead spin down:\t%g\n", */
/*        mean_right2_number_current/(NOISE_NUMBER/stride)); */
/*   tee (flog, "mean particle number current for right lead:\t\t%g\n", */
/*        (mean_right1_number_current+mean_right2_number_current)/ */
/*        (NOISE_NUMBER/stride)); */
/*   tee (flog, "\n"); */


/*   tee (flog, "mean energy current for left lead spin up:\t%g\n", */
/*        mean_left1_energy_current/(NOISE_NUMBER/stride)); */
/*   tee (flog, "mean energy current for left lead spin down:\t%g\n", */
/*        mean_left2_energy_current/(NOISE_NUMBER/stride)); */
/*   tee (flog, "mean energy current for left lead:\t\t%g\n", */
/*        (mean_left1_energy_current+mean_left2_energy_current)/ */
/*        (NOISE_NUMBER/stride)); */
/*   tee (flog, "\n"); */

/*   tee (flog, "mean energy current for right lead spin up:\t%g\n", */
/*        mean_right1_energy_current/(NOISE_NUMBER/stride)); */
/*   tee (flog, "mean energy current for right lead spin down:\t%g\n", */
/*        mean_right2_energy_current/(NOISE_NUMBER/stride)); */
/*   tee (flog, "mean energy current for right lead:\t\t%g\n", */
/*        (mean_right1_energy_current+mean_right2_energy_current)/ */
/*        (NOISE_NUMBER/stride)); */
/*   tee (flog, "\n"); */

  /* close_all (); */

  return 0;
}
