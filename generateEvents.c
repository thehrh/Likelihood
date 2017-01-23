/* 
*Author: Maike Jung
*Date: 15.11.2016

*Purpose: Draw random events from a certain mass spectrum

SN - Model: Lawrence-Livermore
    time spectrum is convoluted with the first hit distribution, to account for not knowing the absolute arrival times

UNITS: mass: eV
       energy: MeV
       distance: Mpc
       time: s

BINNING: defined in header-file (spectrum.h)
*/

#include "spectrum.h"
#include <time.h>

/*
generate random events for: E=0.1-60.0 0.1
                            t=0.0-10.0 0.001
*/

void createEvents(double mass, double distance, double events, bool triggEff, bool energyRes, int filenumber, double *spectrum, double max){

    /*storing time & energy in file*/
    char filename[sizeof "1.5eV_ideal/eventsGenerated_1.45eV_10.5Mpc_1000Events_real_1111.txt"];
    if (triggEff && energyRes){
        sprintf(filename, "events_%.2feV_%.1fMpc_%.0fEvents_real_%d.txt", mass, distance, events, filenumber);
    }
    else {
        sprintf(filename, "events_%.2feV_%.1fMpc_%.0fEvents_ideal_%d.txt", mass, distance, events, filenumber);
    }
    FILE *f = fopen(filename, "w");
    if (f == NULL){
        printf("Error opening file!\n");
        exit(1);
    }

    int eventsGenerated = 0;
    int randE, randT;
    double randCheck;
    while(eventsGenerated < events){
        randE = rand() % (RESE);
        randT = rand() % (REST);
        randCheck = rand()*max/RAND_MAX;
        if (spectrum[randT*(RESE-1)+randE] >= randCheck){
            fprintf(f, "%d %d\n", randE, randT);
            eventsGenerated += 1;
        }
    }
    fclose(f);
}

double findSpectrumMax(double *spectrum){
    // find the maximum value in the spectrum - needed to draw random events
    double max = 0.0;
    int i;
    for(i=0; i<((RESE-1)*REST); i++){
        if(spectrum[i]>max) max = spectrum[i];
    }
    return max;
}



int main(void){
    /*set parameters*/
    /*flag for trigger efficiency*/
    bool useTriggerEff = true;
    bool useEnergyRes = true;
    double mass = 1.0;
    double distance = 5.0;
    double events = 10.0;
    int filenumber, i;
    double noise = pow(10,-5);
    double max;

    // generate spectrum from which time/energy events are drawn
    double *spectrum= (double*) malloc((RESE-1) * REST * sizeof(double));
    createSpectrum(spectrum, mass, distance, events, useEnergyRes, useTriggerEff, noise);

    /*create a file from the spectrum that can then be plotted to look at the spectrum*/
    char filename[sizeof "spectrum_0.1eV_1Mpc_160events.txt"];
    sprintf(filename, "spectrum_%.2feV_%.3fMpc_%.0fevents.txt", mass, distance, events);
    FILE *f = fopen(filename, "w+");
    for(i=0; i<((RESE-1)*REST);i++){
        fprintf(f, "%e\n", spectrum[i]);
    }
    fclose(f);

    max = findSpectrumMax(spectrum);

    srand( (unsigned)time( NULL ) );
    /*calculate uncertainty for certain configuration*/
    for (filenumber=1; filenumber<2; filenumber++){
        printf("creating file %d \n", filenumber);
        createEvents(mass, distance, events, useTriggerEff, useEnergyRes, filenumber, spectrum, max);
    }

    free(spectrum);
    printf("DONE\n");
    return 0;
}
