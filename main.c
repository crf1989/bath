#include <stdio.h>
#include <stdlib.h>

#include "init.h"

int main () 
{
  init ();
  set_output ();

  tee (flog, "Thermalizing system ...\t\t\t");
  fflush (stdout);
  time2 = time(0);
  for (int i = 0; i < NOISE_NUMBER; ++i)
    {
      md_step ();
      next (left1); next (left2);
      next (right1); next (right2);
    }
  tee (flog, "complete, %.8gs is used\n", difftime (time(0), time2));

  tee (flog, "Simulation starts ...\t\t\t");
  time2 = time(0);
  fflush (stdout);


  double dt = TIME_INTERVAL;

  double mean_left1_number_current = 0;
  double mean_left2_number_current = 0;
  double mean_right1_number_current = 0;
  double mean_right2_number_current = 0;
  double mean_left1_energy_current = 0;
  double mean_left2_energy_current = 0;
  double mean_right1_energy_current = 0;
  double mean_right2_energy_current = 0;

  for (int i = 0; i < TIME_STEP_NUMBER-1; ++i)
    {
      int m = left1->pos%center1->n;
      int n = (m+1)%center1->n;
      
      md_step ();

      if (i%stride != 0)
	{
	  next (left1); next (left2);
	  next (right1); next (right2);
	  continue;
	}
      if (WRITE_PSI)
  	{
  	  for (int j = 0; j < center1->size; ++j)
  	    fprintf (fpsi[j], "%.8g\t%.8g\t%.8g\t%.8g\t%.8g\n",
  		     (i+1)*dt,
		     creal (center1->psi[j]), cimag (center1->psi[j]),
		     creal(center2->psi[j]), cimag(center2->psi[j]));
  	  }
      if (WRITE_NORM)
  	fprintf (fnorm, "%.8g\t%.8g\n", (i+1)*dt,
		 get_norm (center1)+get_norm(center2));
      if (WRITE_SPIN)
	{
	  for (int j = 0; j < center1->size+1; ++j)
	    fprintf (fspin[j], "%.8g\t%.8g\t%.8g\t%.8g\t%.8g\n",
		     (i+1)*dt,
		     spin[j].s[0], 
		     spin[j].s[1], 
		     spin[j].s[2],
		     sqrt(square(spin[j].s[0])
			  +square(spin[j].s[1])
			  +square(spin[j].s[2])));
	}

      {
      	double left1_number_current =
      	  2*cimag (conj(center1->cl[m])*(eta(left1)+left1_fric));
      	double left2_number_current =
      	  2*cimag (conj(center2->cl[m])*(eta(left2)+left2_fric));
      	double right1_number_current =
      	  2*cimag (conj(center1->cr[m])*(eta(right1)+right1_fric));
      	double right2_number_current =
      	  2*cimag (conj(center2->cr[m])*(eta(right2)+right2_fric));
      
      	mean_left1_number_current += left1_number_current;
      	mean_left2_number_current += left2_number_current;
      	mean_right1_number_current += right1_number_current;
      	mean_right2_number_current += right2_number_current;
	
      	fprintf (fnc, "%g\t%g\t%g\t%g\t%g\n", i*dt,
      		 left1_number_current, left2_number_current,
      		 right1_number_current, right2_number_current);
      	fprintf (fmnc, "%g\t%g\t%g\t%g\t%g\n", i*dt,
      		 mean_left1_number_current/(i/stride),
      		 mean_left2_number_current/(i/stride),
      		 mean_right1_number_current/(i/stride),
      		 mean_right2_number_current/(i/stride));
      }

      {
      	double left1_energy_current =
      	  -2*creal (conj(eta(left1)+left1_fric)*
      		    (center1->cl[n]-center1->cl[m])/dt);
      	double left2_energy_current =
      	  -2*creal (conj(eta(left2)+left2_fric)*
      		    (center2->cl[n]-center2->cl[m])/dt);
      	double right1_energy_current =
      	  -2*creal (conj(eta(right1)+right1_fric)*
      		    (center1->cr[n]-center1->cr[m])/dt);
      	double right2_energy_current =
      	  -2*creal (conj(eta(right2)+right2_fric)*
      		    (center2->cr[n]-center2->cr[m])/dt);
	
      	mean_left1_energy_current += left1_energy_current;
      	mean_left2_energy_current += left2_energy_current;
      	mean_right1_energy_current += right1_energy_current;
      	mean_right2_energy_current += right2_energy_current;
	
      	fprintf (fec, "%g\t%g\t%g\t%g\t%g\n", i*dt,
      		 left1_energy_current, left2_energy_current,
      		 right1_energy_current, right2_energy_current);
      	fprintf (fmec, "%g\t%g\t%g\t%g\t%g\n", i*dt,
      		 mean_left1_energy_current/(i/stride),
      		 mean_left2_energy_current/(i/stride),
      		 mean_right1_energy_current/(i/stride),
      		 mean_right2_energy_current/(i/stride));
      }
      next (left1); next (left2);
      next (right1); next (right2);
    }
  
  tee (flog, "complete, %.8gs is used\n", difftime(time(0),time2));

  tee (flog, "mean particle number current:\n");
  tee (flog, "\tleft up:\t\t\t%g\n",
       mean_left1_number_current/(TIME_STEP_NUMBER/stride));
  tee (flog, "\tleft down:\t\t\t%g\n",
       mean_left2_number_current/(TIME_STEP_NUMBER/stride));
  tee (flog, "\tleft total:\t\t\t%g\n",
       (mean_left1_number_current+mean_left2_number_current)/
       (TIME_STEP_NUMBER/stride));
  tee (flog, "\n");

  tee (flog, "\tright up:\t\t\t%g\n",
       mean_right1_number_current/(TIME_STEP_NUMBER/stride));
  tee (flog, "\tright down:\t\t\t%g\n",
       mean_right2_number_current/(TIME_STEP_NUMBER/stride));
  tee (flog, "\tright total:\t\t\t%g\n",
       (mean_right1_number_current+mean_right2_number_current)/
       (TIME_STEP_NUMBER/stride));
  tee (flog, "\n");


  tee (flog, "mean energy current:\n");
  tee (flog, "\tleft up:\t\t\t%g\n",
       mean_left1_energy_current/(TIME_STEP_NUMBER/stride));
  tee (flog, "\tleft down:\t\t\t%g\n",
       mean_left2_energy_current/(TIME_STEP_NUMBER/stride));
  tee (flog, "\tleft total:\t\t\t%g\n",
       (mean_left1_energy_current+mean_left2_energy_current)/
       (TIME_STEP_NUMBER/stride));
  tee (flog, "\n");

  tee (flog, "\tright up:\t\t\t%g\n",
       mean_right1_energy_current/(TIME_STEP_NUMBER/stride));
  tee (flog, "\tright down:\t\t\t%g\n",
       mean_right2_energy_current/(TIME_STEP_NUMBER/stride));
  tee (flog, "\tright total:\t\t\t%g\n",
       (mean_right1_energy_current+mean_right2_energy_current)/
       (TIME_STEP_NUMBER/stride));
  tee (flog, "\n");

  close_all ();

  return 0;
}
