#ifndef INIT_H
#define INIT_H 1

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>

#include "system.h"

time_t time1, time2;

FILE** fpsi = NULL;
FILE** fspin = NULL;
FILE* fnorm = NULL;
FILE* fnc = NULL;		/* number current */
FILE* fec = NULL;		/* energy current */
FILE* fmnc = NULL;		/* mean number current */
FILE* fmec = NULL;		/* mean energy current */

void init ();
void set_output ();
void close_all ();

/*******************************************************/

void init ()
{
  flog = fopen ("log", "w");
  assert (flog != NULL);
  
  time1 = time(0);
  tee (flog, "Prgram starts at %s\n", ctime(&time1));

  srandom (time(0));

  NOISE_NUMBER = pow2 (NOISE_NUMBER_POW);
  SELF_ENERGY_NUMBER = pow2 (SELF_ENERGY_NUMBER_POW);
  TIME_STEP_NUMBER = pow2 (TIME_STEP_NUMBER_POW);

  tee (flog,
       "S\t\t\t\t\t=\t%.8g\n"
       "J1\t\t\t\t\t=\t%.8g\n"
       "JH\t\t\t\t\t=\t%.8g\n"
       "CENTER_SITE_NUMBER\t\t\t=\t%d\n"
       "NOISE_NUMBER\t\t\t\t=\t2^%d\t=\t%d\n"
       "SELF_ENERGY_NUMBER\t\t\t=\t2^%d\t=\t%d\n"
       "TIME_STEP_NUMBER\t\t\t=\t2^%d\t=\t%d\n"
       "HOPPING_ENERGY\t\t\t\t=\t%.8g eV\n"
       "TIME_INTERVAL\t\t\t\t=\t%.8g hbar/eV\n"
       "LEFT_LEAD_TEMPERATURE\t\t\t=\t%.8g K\n"
       "LEFT_UP_LEAD_CHEMICAL_POTENTIAL\t\t=\t%.8g eV\n"
       "LEFT_DOWN_LEAD_CHEMICAL_POTENTIAL\t=\t%.8g eV\n"
       "RIGHT_LEAD_TEMPERATURE\t\t\t=\t%.8g K\n"
       "RIGHT_UP_LEAD_CHEMICAL_POTENIAL\t\t=\t%.8g eV\n"
       "RIGHT_DOWN_LEAD_CHEMICAL_POTENIAL\t=\t%.8g eV\n\n",
       S, J1, JH,
       CENTER_SITE_NUMBER, 
       NOISE_NUMBER_POW, NOISE_NUMBER,
       SELF_ENERGY_NUMBER_POW, SELF_ENERGY_NUMBER,
       TIME_STEP_NUMBER_POW, TIME_STEP_NUMBER,
       HOPPING_ENERGY, TIME_INTERVAL,
       LEFT_LEAD_TEMPERATURE, 
       LEFT_UP_LEAD_CHEMICAL_POTENTIAL, LEFT_DOWN_LEAD_CHEMICAL_POTENTIAL,
       RIGHT_LEAD_TEMPERATURE, 
       RIGHT_UP_LEAD_CHEMICAL_POTENTIAL, RIGHT_DOWN_LEAED_CHEMICAL_POTENTIAL);

  if (WRITE_SPECTRAL)
    tee (flog, "WRITE_SPECTRAL\t\t\t\ton\n");
  else
    tee (flog, "WRITE_SPECTRAL\t\t\t\toff\n");
  if (WRITE_FRIC)
    tee (flog, "WRITE_FRIC\t\t\t\ton\n");
  else
    tee (flog, "WRITE_FRIC\t\t\t\toff\n");
  if (WRITE_NOISE)
    tee (flog, "WRITE_NOISE\t\t\t\ton\n");
  else
    tee (flog, "WRITE_NOISE\t\t\t\toff\n");
  if (WRITE_ETA)
    tee (flog, "WRITE_ETA\t\t\t\ton\n");
  else
    tee (flog, "WRITE_ETA\t\t\t\toff\n");
  if (WRITE_CORRELATION)
    tee (flog, "WRITE_CORRELATION\t\t\ton\n");
  else
    tee (flog, "WRITE_CORRELATION\t\t\toff\n");
  if (WRITE_PSI)
    tee (flog, "WRITE_PSI\t\t\t\ton\n");
  else
    tee (flog, "WRITE_PSI\t\t\t\toff\n");
  if (WRITE_NORM)
    tee (flog, "WRITE_NORM\t\t\t\ton\n");
  else
    tee (flog, "WRITE_NORM\t\t\t\toff\n");
  if (WRITE_SPIN)
    tee (flog, "WRITE_SPIN\t\t\t\ton\n");
  else
    tee (flog, "WRITE_SPIN\t\t\t\toff\n");
  tee (flog, "\n");
  
  tee (flog, "omega interval dw = %.8g eV, range [%.8g,%.8g] eV\n\n",
       2*PI/(NOISE_NUMBER*TIME_INTERVAL),
       -PI/TIME_INTERVAL, PI/TIME_INTERVAL);

  time2 = time(0);
  tee (flog, "Initializing system ...\t\t\t");
  init_system ();
  tee (flog, "complete, %.8gs is used\n", difftime(time(0), time2));
  fflush (stdout);

}

void set_output ()
{
  if (WRITE_PSI)
    {
      assert (system ("rm -rf psi") != -1);
      assert (system ("mkdir -p psi") != -1);
      fpsi = malloc (CENTER_SITE_NUMBER * sizeof(FILE*));
      assert (fpsi != NULL);
      char filename[128];
      for (int i = 0; i < CENTER_SITE_NUMBER; ++i)
	{
	  sprintf (filename, "psi/%d.dat", i);
	  fpsi[i] = fopen (filename, "w");
	  assert (fpsi[i] != NULL);
	  fprintf (fpsi[i], 
		   "#This is psi at site %d with respect with time\n"
		   "#time\tupreal\tupimag\tdownreal\tdownimag\n", i);
	}
    }
  if (WRITE_NORM)
    {
      assert (system ("mkdir -p check/norm") != -1);
      fnorm = fopen ("check/norm/norm.dat", "w");
      assert (fnorm != NULL);
      fprintf (fnorm, 
	       "#This is the norm of psi with respect with time\n"
	       "#time\treal\timag\n");
    }
  if (WRITE_SPIN)
    {
      assert (system ("rm -rf spin") != -1);
      assert (system ("mkdir -p spin") != -1);
      fspin = malloc ((CENTER_SITE_NUMBER+1) * sizeof(FILE*));
      assert (fspin != NULL);
      char filename[128];
      for (int i = 0; i < CENTER_SITE_NUMBER+1; ++i)
	{
	  sprintf (filename, "spin/%d.dat", i);
	  fspin[i] = fopen (filename, "w");
	  assert (fspin[i] != NULL);
	  fprintf (fspin[i],
		   "#This is the vector of spin with respect with time\n"
		   "#time\tx\ty\tz\tS\n");
	}
    }

  
  assert (system ("mkdir -p current") != -1);
  fnc = fopen ("current/number_current.dat", "w");
  assert (fnc != NULL);
  fec = fopen ("current/energy_current.dat", "w");
  assert (fec != NULL);
  fmnc = fopen ("current/mean_number_current.dat", "w");
  assert (fmnc != NULL);
  fmec = fopen ("current/mean_energy_current.dat", "w");
  assert (fmec != NULL);

  fprintf (fnc, 
	   "#This is the particle number current\n"
	   "#time\tlu\tru\tld\trd\n");
  fprintf (fmnc,
	   "#This is the mean particle number current\n"
	   "#time\tlu\tru\tld\trd\n");
  fprintf (fec,
	   "#This is energy current\n"
	   "#time\tlu\tru\tld\trd\n");
  fprintf (fmec, 
	   "#This is the mean energy current\n"
	   "#time\tlu\tru\tld\trd\n");
}

void close_all ()
{
  if (WRITE_PSI)
    {
      for (int i = 0; i < CENTER_SITE_NUMBER; ++i)
	fclose (fpsi[i]);
      free (fpsi);
    }
  if (WRITE_NORM)
    {
      fclose (fnorm);
    }
  if (WRITE_SPIN)
    {
      for (int i = 0; i < CENTER_SITE_NUMBER+1; ++i)
	fclose (fspin[i]);
      free (fspin);
    }
  free_bath (left1);
  free_bath (left2);
  free_bath (right1);
  free_bath (right2);
  chain_clear (center1);
  chain_clear (center2);
  spin_free (spin);

  fclose (fnc);
  fclose (fec);
  fclose (fmnc);
  fclose (fmec);

  time2 = time(0);
  tee (flog, "Prgram finishes at %s", ctime(&time2));
  tee (flog, "Totally %gs has been used\n", difftime (time2,time1));
}
  
#endif	/* INIT_H */
