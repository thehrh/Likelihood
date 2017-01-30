#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

/*define size of the grid*/
#define RESE 600
#define REST 1000
#define EMAX 60.0
#define TMAX 10.0
#define STEPE EMAX/RESE
#define STEPT TMAX/REST
#define THREADS_PER_BLOCK 128

#ifdef USE_SP
typedef float user_data_t;
#else
typedef double user_data_t;
#endif

/* normalized LL energy spectrum - 15.4: average energy, 3.8: betha, 4802: normalization factor*/
#define LL_energy_spectrum(E) pow(E, 3.8)*exp(-(1.0+3.8)*E/15.4)/4802.516160
/* normalized LL arrival time spectrum - approximated by log-normal-dist */
#define MY -1.0324
#define SIGMA 0.9134
#define LL_time_spectrum(t) exp( - (log(t)-MY)*(log(t)-MY)/(2*SIGMA*SIGMA) ) / (t*SIGMA*sqrt(2*M_PI))
#define ALPHA 2.22
/*gaussian (around 0): error for the energy smearing from photon counts:
N is proportional to E -> N = 1/alpha*E -> N=alpha*E -> factor of alpha for sigma2 (sigma2=E)
*/
#define GAUSS(i, sigma2) 1/sqrt(2.0*M_PI*sigma2*ALPHA)*exp(-(i*i)/(2.0*sigma2*ALPHA))

#ifndef SPECTRUM_H_
#define SPECTRUM_H_

void generateDist(user_data_t mass, user_data_t dist, user_data_t events, user_data_t *distribution, user_data_t *triggerEffs, bool useEnergyRes);
void createSpectrum(user_data_t *spectrum, user_data_t mass, user_data_t distance, user_data_t events, bool useEnergyRes, bool useTriggerEff, user_data_t noise);
void getEvent(int *eventEnergy, int *eventTime, double mass, double distance, double events, int filenumber);
double getLLH(double mass, double distance, double events, bool triggEff, bool energyRes, double noise, int *eventTime, int *eventEnergy);

#endif
