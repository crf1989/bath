#ifndef PARAMETER_H
#define PARAMETER_H 1

const int CENTER_SITE_NUMBER = 10;
const int NOISE_NUMBER_POW = 20;
const int SELF_ENERGY_NUMBER_POW = 10;
const int TIME_STEP_NUMBER_POW = 22;
int NOISE_NUMBER;
int SELF_ENERGY_NUMBER;
int TIME_STEP_NUMBER;

const double HOPPING_ENERGY = 1; /* eV */
const double TIME_INTERVAL = 1e-3;
const double S = 3.968626966596886;
const double J1 = 0;  	/* eV */
const double JH = 0;		/* eV */
const double B = 0;		/* magnetic field, T */

const double LEFT_LEAD_TEMPERATURE = 50; /* Kelvin */
const double LEFT_UP_LEAD_CHEMICAL_POTENTIAL = -1.9; /* eV */
const double LEFT_DOWN_LEAD_CHEMICAL_POTENTIAL = -1.9; /* eV */

const double RIGHT_LEAD_TEMPERATURE = 50; /* kelvin */
const double RIGHT_UP_LEAD_CHEMICAL_POTENTIAL = -2; /* eV */
const double RIGHT_DOWN_LEAED_CHEMICAL_POTENTIAL = -2; /* eV */

const int stride = 100;

#define WRITE_SPECTRAL 0
#define WRITE_FRIC 0
#define WRITE_NOISE 0
#define WRITE_ETA 0
#define WRITE_CORRELATION 0

#define WRITE_PSI 0
#define WRITE_NORM 0
#define WRITE_SPIN 0

#endif /* PARAMETER_H */
