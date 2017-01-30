/* File needed for SWIG to build the python module "spectrum" */

%module spectrum

%include "carrays.i"
%array_class(double, doubleArray);
%array_class(int, intArray);

%{
#define SWIG_FILE_WITH_INIT
#include "spectrum.h"
%}

#ifdef USE_SP
typedef float user_data_t;
#else
typedef double user_data_t;
#endif

void createSpectrum(user_data_t *spectrum, user_data_t mass, user_data_t distance, user_data_t events, bool useEnergyRes, bool useTriggerEff, user_data_t noise);

void getEvent(int *eventEnergy, int *eventTime, double mass, double distance, double events, int filenumber);
double getLLH(double mass, double distance, double events, bool triggEff, bool energyRes, double noise, int *eventTime, int *eventEnergy);
