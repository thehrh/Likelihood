/* 
*Author: Maike Jung
*Date: 26.01.2017

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
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>

/*
generate random events for: E=0.1-60.0 0.1
                            t=0.0-10.0 0.001
*/

void fillTriggerEff(user_data_t *triggerEffs, bool useTriggerEff){
    /*Read in trigger efficiency based on the chosen resolution (eff. vs. energy).
     Note that the trigger efficiency file for the chosen resolution needs
     to be located in the current directory.*/
    int i;
    user_data_t triggerEns[RESE+1]; // TODO: why do we need this at all?
    if(useTriggerEff){
        FILE *myFile;
        if (RESE==600){
            myFile = fopen("trigger_efficiency_100keV_steps.txt", "r");
        }
        else if (RESE==6000){
            myFile = fopen("trigger_efficiency_10keV_steps.txt", "r");
        }
        else if (RESE==60000){
            myFile = fopen("trigger_efficiency_1keV_steps.txt", "r");
        }
        else {
            printf("Invalid grid size for the energy resolution.");
        }
        for (i = 0; i < RESE+1; i++) {
            fscanf(myFile, "%lf %lf", &triggerEns[i], &triggerEffs[i]);
        }
        fclose(myFile);
    }
    else{
        /* initialize with 1s for ideal trigger efficiency */
        for(i = 0; i < RESE+1 ; triggerEffs[i++] = 1.0);
    }
}

void addNoise(user_data_t *spectrum, user_data_t noise){
    int i;
    // add constant noise floor to the spectrum
    for (i=0; i<(RESE-1)*REST; i++){
        //spectrum[i] *= 0.99;
        spectrum[i] += noise;
    }
}

void createSpectrum(user_data_t *spectrum, user_data_t mass, user_data_t distance, user_data_t events, bool useEnergyRes, bool useTriggerEff, user_data_t noise){
    user_data_t triggerEffs[RESE+1];

    /*get trigger efficiencies as function of energy*/
    fillTriggerEff(triggerEffs, useTriggerEff);

    /*create the spectrum from which the random events are drawn*/
    generateDist(mass, distance, events, spectrum, triggerEffs, useEnergyRes);

    /*sprinkle with some noise*/
    addNoise(spectrum, noise);

}

void createEvents(user_data_t mass, user_data_t distance, user_data_t events, bool triggEff, bool energyRes, int filenumber, user_data_t *spectrum, user_data_t max){

    /* file for storing time & energy */
    char filename[sizeof "DATA/10.00Mpc_700Events_1.57eV_event_1.45eV_10.5Mpc_1000Events_real_1111.txt"];

    if (triggEff && energyRes){
        sprintf(filename, "DATA/%.2fMpc_%.0fEvents_%.2feV/events_%.2feV_%.2fMpc_%.0fEvents_real_%d.txt",distance, events, mass, mass, distance, events, filenumber);
    }
    else {


        sprintf(filename, "DATA/events_%.2feV_%.1fMpc_%.0fEvents_ideal_%d_test.txt", mass, distance, events, filenumber);

    }
    FILE *f = fopen(filename, "w");
    if (f == NULL){
        printf("Error opening file!\n");
        exit(1);
    }

    int eventsGenerated = 0;
    int randE, randT;
    user_data_t randCheck;
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

user_data_t findSpectrumMax(user_data_t *spectrum){
    // find the maximum value in the spectrum - needed to draw random events
    user_data_t max = 0.0;
    int i;
    for(i=0; i<((RESE-1)*REST); i++){
        if(spectrum[i]>max) max = spectrum[i];
    }
    return max;
}


int main(void){
    bool triggEff = true;
    bool energyRes = true;
    bool plot = false;
    user_data_t mass = 1.0;
    user_data_t distance = 5.0;
    user_data_t events = 10.0;
    user_data_t max;
    int filenumber, i;
    double noise_rate = pow(10,-3); //Hz, expected total noise rate
    // noise that has to be added to each bin:
    user_data_t noise = noise_rate*STEPT;
    
    // generate spectrum from whicht time/energy events are drawn
    user_data_t spectrum[(RESE-1)*REST];
	//double *spectrum= (double*) malloc((RESE-1) * REST * sizeof(double));
    createSpectrum(spectrum, mass, distance, events, energyRes, triggEff, noise);

    if(plot){
        /*create a file from the spectrum that can then be ploted*/
        char filename[sizeof "spectrum_0.1eV_1Mpc_160events_test.txt"];
        sprintf(filename, "spectrum_%.2feV_%.3fMpc_%.0fevents_test.txt",mass, distance, events);
        FILE *f = fopen(filename, "w+");
        for(i=0; i<((RESE-1)*REST);i++){
            fprintf(f, "%e\n", spectrum[i]);
        }
        fclose(f);
    }

    max = findSpectrumMax(spectrum);

    srand( (unsigned)time( NULL ) );

    /*create files that contain pseudo experiments*/
    // create directory
    /*
    mkdir("DATA", S_IRWXU);
    char dirname[sizeof "DATA/10.00Mpc_700Events_1.57eV"];
    sprintf(dirname, "DATA/%.2fMpc_%.0fEvents_%.2feV", distance, events, mass);
    mkdir(dirname, S_IRWXU);

    for (filenumber=1; filenumber<1001; filenumber++){ 
        printf("creating file %d \n", filenumber);
        createEvents(mass, distance, events, useTriggerEff, useEnergyRes, filenumber, spectrum, max);
    }
    */

    //free(spectrum);
    printf("DONE\n");
    return 0;
}
