#ifndef BATH_H
#define BATH_H 1

#include <stdio.h>
#include <complex.h>
#include <fftw3.h>
#include <assert.h>

#include "auxiliary.h"
#include "parameter.h"

typedef struct 
{
  int size;
  int pos;
  double dt;			/* time interval */
  double T;			/* temperature */
  double u;			/* chemical potential */
  double h;			/* hopping energy */
  double dw;			/* omega interval */
  double* spectral;
  double complex* fric; 	/* self energy */
  double complex* noise;
  double complex* eta;
} bath;

bath* bath_alloc (int size, double dt, double T,
		  double u, double h);
void free_bath (bath*);
void generate_spectral (bath*);
void generate_fric (bath*);
void generate_noise (bath*);
void generate_eta (bath*);
void generate_bath (bath*);
void regenerate_bath (bath*);
double complex eta (bath*);
void next (bath*);


void write_spectral (bath*, char* MESSAGE, char* filename);
void write_fric (bath*, char* MESSAGE, char* filename);
void write_noise (bath*, char* MESSAGE, char* filename);
void write_eta (bath*, char* MESSAGE, char* filename);
void write_correlation (bath*, char* MESSAGE, char* filename);


/*************************************************************/

bath* bath_alloc (int size, double dt, double T,
		  double u, double h)
{
  bath* p = malloc (sizeof(bath));
  p->size = size;
  p->pos = 0;
  p->dt = dt;
  p->T = T;
  p->u = u;
  p->h = h;
  p->dw = 2*PI/(size*dt);
  p->spectral = malloc (size * sizeof(double));
  assert (p->spectral != NULL);
  p->fric = malloc (size * sizeof(double complex));
  assert (p->fric != NULL);
  p->noise = malloc (size * sizeof(double complex));
  assert (p->noise != NULL);
  p->eta = malloc (size * sizeof(double complex));
  assert (p->eta != NULL);
  return p;
}

void free_bath (bath* p)
{
  free (p->spectral);
  free (p->fric);
  free (p->noise);
  free (p->eta);
  free (p);
}

void generate_spectral (bath* p)
{
  for (int i = 0; i <= p->size/2; ++i)
    {
      p->fric[i] = sfg_solve (p->dw*i, p->h);
      p->spectral[i] = fermi (p->dw*i, p->u, p->T) *
	(-2)*cimag (p->fric[i]);
    }
  for (int i = p->size/2+1; i < p->size; ++i)
    {
      p->fric[i] = sfg_solve (p->dw*(i-p->size), p->h);
      p->spectral[i] = fermi (p->dw*(i-p->size), p->u, p->T) *
  	(-2)*cimag (p->fric[i]);
    }
}

void generate_fric (bath* p)
{
  fftw_plan plan = fftw_plan_dft_1d (p->size,
				     p->fric,
				     p->fric,
				     FFTW_BACKWARD,
				     FFTW_ESTIMATE);
  fftw_execute (plan);
  double n = sqrt (p->size);
  for (int i = 0; i < p->size; ++i)
      p->fric[i] /= n;
  fftw_destroy_plan (plan);
}

void generate_noise (bath* p)
{
  for (int i = 0; i < p->size; ++i)
    box_muller (&p->noise[i],
  		sqrt (0.5*p->size*p->dt*p->spectral[i]));
}

void generate_eta (bath* p)
{
  fftw_plan fft = fftw_plan_dft_1d (p->size,
  				    p->noise,
  				    p->eta,
  				    FFTW_BACKWARD,
  				    FFTW_ESTIMATE);
  fftw_execute (fft);
  double n = sqrt (p->size);
  for (int i = 0; i < p->size; ++i)
    p->eta[i] /= n;

  fftw_destroy_plan (fft);
}

void generate_bath (bath* p)
{
  generate_spectral (p);
  generate_fric (p);
  generate_noise (p);
  generate_eta (p);
}

void regenerate_bath (bath* p)
{
  generate_noise (p);
  generate_eta (p);
}  

double complex eta (bath* p)
{
  return p->eta[p->pos];
}

void next (bath* p)
{
  ++p->pos;
  if (p->pos == p->size)
    {
      p->pos = 0;
      regenerate_bath (p);
    }
}

void write_spectral (bath* p, char* MESSAGE, char* filename)
{
  FILE* fp = fopen (filename, "w");
  assert (fp != NULL);
  
  fprintf (fp, "%s", MESSAGE);
  for (int i = 0; i < p->size; ++i)
    fprintf (fp, "%.8g\t%.8g\n",
	     i*p->dw, p->spectral[i]);

  fclose (fp);
}

void write_fric (bath* p, char* MESSAGE, char* filename)
{
  FILE* fp = fopen (filename, "w");
  assert (fp != NULL);

  fprintf (fp, "%s", MESSAGE);
  for (int i = 0; i < p->size; ++i)
    fprintf (fp, "%.8g\t%.8g\t%.8g\n", i*p->dt, 
	     creal (p->fric[i]), cimag (p->fric[i]));
  
  fclose (fp);
}  

void write_noise (bath* p, char* MESSAGE, char* filename)
{
  FILE* fp = fopen (filename, "w");
  assert (fp != NULL);

  fprintf (fp, "%s", MESSAGE);
  for (int i = 0; i < p->size; ++i)
    fprintf (fp, "%.8g\t%.8g\t%.8g\n", i*p->dw, 
	     creal (p->noise[i]), cimag (p->noise[i]));
  
  fclose (fp);
}

void write_eta (bath* p, char* MESSAGE, char* filename)
{
  FILE* fp = fopen (filename, "w");
  assert (fp != NULL);

  fprintf (fp, "%s", MESSAGE);
  for (int i = 0; i < p->size; ++i)
    fprintf (fp, "%.8g\t%.8g\t%.8g\n", i*p->dt, 
	     creal (p->eta[i]), cimag (p->eta[i]));
  
  fclose (fp);
}

void write_correlation (bath* p, char* MESSAGE, char* filename)
{
  FILE* fp = fopen (filename, "w");
  assert (fp != NULL);

  fprintf (fp, "%s", MESSAGE);
  for (int i = 0; i < p->size/2; ++i)
    {
      double complex corr = 0;
      for (int j = 0; j < p->size-i; ++j)
	corr += p->eta[j]*conj(p->eta[j+i]);
      corr /= p->size-i;
      fprintf (fp, "%.8g\t%.8g\t%.8g\n", i*p->dt, 
	       creal (corr), cimag (corr));
    }
}
	  
  
#endif /* BATH_H */
