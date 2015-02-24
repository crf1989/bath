#ifndef PARAMETER_H
#define PARAMETER_H 1

const int CENTER_SITE_NUMBER = 10;
const int NOISE_NUMBER_POW = 22;
const int SELF_ENERGY_NUMBER_POW = 10;
const int TIME_STEP_NUMBER_POW = 24;
int NOISE_NUMBER;
int SELF_ENERGY_NUMBER;
int TIME_STEP_NUMBER;

const double HOPPING_ENERGY = 1; /* eV */
const double TIME_INTERVAL = 1e-3;
const double S = 3.968626966596886;
const double J1 = 0;  	/* eV */
const double JH = 0.1;		/* eV */
const double B = 0;		/* magnetic field, T */

const double u = 5.7883818066e-5; /* bohr magneton, eV/T */
const double g = 2;		  /* lande factor */

const double LEFT_LEAD_TEMPERATURE = 50; /* Kelvin */
const double LEFT_UP_LEAD_CHEMICAL_POTENTIAL = -1.6; /* eV */
const double LEFT_DOWN_LEAD_CHEMICAL_POTENTIAL = -1.7; /* eV */

const double RIGHT_LEAD_TEMPERATURE = 50; /* kelvin */
const double RIGHT_UP_LEAD_CHEMICAL_POTENTIAL = -1.8; /* eV */
const double RIGHT_DOWN_LEAED_CHEMICAL_POTENTIAL = -1.9; /* eV */

const int stride = 100;

#define WRITE_SPECTRAL 1
#define WRITE_FRIC 1
#define WRITE_NOISE 1
#define WRITE_ETA 1
#define WRITE_CORRELATION 0

#define WRITE_PSI 1
#define WRITE_NORM 1
#define WRITE_SPIN 1

#endif /* PARAMETER_H */
