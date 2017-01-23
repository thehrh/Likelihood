/* File: spectrum.i */


%module spectrum

%include "carrays.i"
%array_class(double, doubleArray);
%array_class(int, intArray);

%{
#define SWIG_FILE_WITH_INIT
#include "spectrum.h"
%}

double LLSpectrumTotal (double t, double E, double mass, double dist);
void ProbFirstHitDist (double mass, double dist, double events, double *result);
void correlation(double mass, double dist, double events, double *newSpec);
void generateDist(double mass, double dist, double events, double *distribution, double *triggerEfficiency, bool energyRes);
void createSpectrum(double *spectrum, double mass, double distance, double events, bool energyRes, bool triggEff, double noise);
void getTriggerEfficiency(double *triggerEfficiency, bool triggEff);
void addNoise(double *spectrum, double noise);
void getEvent(int *eventEnergy, int *eventTime, double mass, double distance, double events, int filenumber);
double getLLH(double mass, double distance, double events, double *triggerEfficiency, bool energyRes, double noise, int *eventTime, int *eventEnergy);
