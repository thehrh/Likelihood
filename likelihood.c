/* 
*Author: Maike Jung
*Date: 15.11.2016

*Purpose: Calculate the likelihood for the random events generated with generateEvents.c to belong to a certain mass spectrum

SN - Model: Lawrence-Livermore
    time spectrum is convoluted with the first hit distribution, to account for not knowing the absolute arrival times

UNITS: mass: eV
       energy: MeV
       distance: Mpc
       time: s

BINNING: see spectrum.h

*/

#include "spectrum.h"

/// load event
void getEvent(int *eventEnergy, int *eventTime, double mass, double distance, double events, int filenumber){

    /*load events & store energy and time in arrays*/
    char filename[sizeof "1eV_ideal/eventsGenerated_1.34eV_10.5Mpc_1000Events_ideal_1111.txt"];
    sprintf(filename, "events_%.2feV_%.1fMpc_%.0fEvents_real_%d.txt", mass, distance, events, filenumber);

    FILE *f = fopen(filename, "r");
    int i;
    for(i = 0; i < events; eventEnergy[i++] = 1);
    for(i = 0; i < events; eventTime[i++] = 1);
    for (i = 0; i < events; i++){
        fscanf(f, "%d %d", &eventEnergy[i], &eventTime[i]);
    }
}
///function that calculates likelihood - and will then be minimized in python
double getLLH(double mass, double distance, double events, bool triggEff, bool energyRes, double noise, int *eventTime, int *eventEnergy){

    double llh = 0.0;
    int i;
    user_data_t spectrum[(RESE-1)*REST];
	//double *spectrum= (double*) malloc((RESE-1) * REST * sizeof(double));
    createSpectrum(spectrum, mass, distance, events, energyRes, triggEff, noise);

    for (i = 0; i < events; i++){
        if (spectrum[eventTime[i]*(RESE-1)+eventEnergy[i]] < pow(10,-200)){
            llh += -10000000;   
            printf("value of spectrum very small - check");
        }
        else llh += log(spectrum[eventTime[i]*(RESE-1)+eventEnergy[i]]);
    }
    llh*=-1;
    return llh;
}


// method that scans over range
void calcLLH(double mass, double distance, double events, bool triggEff, bool energyRes, int filenumber, double noise){
   
    /*load events & store energy and time in arrays*/
    int eventEnergy[(int) events];
    int eventTime[(int) events];
    getEvent(eventEnergy, eventTime, mass, distance, events, filenumber);

    // calculate the likelihood
    int i;
    double llh;
    double testMass;
    // store current value of the minimumLLH and the corresponding mass
    double minLLH = INFINITY;
    double massOfMinLLH = 0.0;
    // go over all the spectra around a certain range of the input mass & calculate the likelihood for each spectrum
    double *testSpectrum= (double*) malloc((RESE-1) * REST * sizeof(double));
    // first go over broad range - there are no negative entries in the spectrum!!!!!
    for (testMass = mass - 0.5; testMass <= mass + 0.5; testMass+=0.1){
        llh = 0.0;
        createSpectrum(testSpectrum, testMass, distance, events, energyRes, triggEff, noise);
        for (i = 0; i < events; i++){
            if (testSpectrum[eventTime[i]*(RESE-1)+eventEnergy[i]] < pow(10,-200)){
                llh += -10000000;   
            }
            else llh += log(testSpectrum[eventTime[i]*(RESE-1)+eventEnergy[i]]);
            //printf("tset %d %f \n", i, llh);
        }
        llh*=-1;
        printf("mass %f, llh %f\n",testMass, llh);
        if (llh < minLLH) {
            minLLH = llh;
            massOfMinLLH = testMass;
        }
    }
    double currentMinimum = massOfMinLLH;

    // check around the minimum mass in finer steps
    for (testMass = currentMinimum - 0.05; testMass <= currentMinimum + 0.05; testMass += 0.01){
        if(testMass >= 0.0){
            llh = 0.0;
            createSpectrum(testSpectrum, testMass, distance, events, energyRes, triggEff, noise);
            for (i = 0; i < events; i++){
                if (testSpectrum[eventTime[i]*(RESE-1)+eventEnergy[i]] < pow(10,-200)){
                    llh += -10000000;   
                }
                else llh += log(testSpectrum[eventTime[i]*(RESE-1)+eventEnergy[i]]);
                //printf("tset %d %f \n", i, llh);
            }
            llh*=-1;
            printf("mass %f, llh %f \n",testMass, llh);
            if (llh < minLLH) {
                minLLH = llh;
                massOfMinLLH = testMass;
            }
        }
    }

    printf("mass found: %f eV\n", massOfMinLLH);

    free(testSpectrum);

    /*write to file*/
    char filename2[sizeof "Results_Likelihood/test_masses_1.55eV_1Mpc_real.txt"];
    if (triggEff && energyRes){
        sprintf(filename2, "TEST_masses_%.2feV_%.1fMpc_real_1.txt", mass, distance);
    }
    else {
        sprintf(filename2, "Results_Likelihood/masses_%.2feV_%.1fMpc_ideal_test3.txt", mass, distance);
    }
    FILE *g = fopen(filename2, "a+");
    fprintf(g, "%f \n", massOfMinLLH);
    fclose(g);
}

int main(void){
	/*set parameters*/
    /*flag for trigger efficiency*/
    bool triggEff = true;
    bool energyRes = true;
    double mass = 1.0;
    double distance = 1.0;
    double events = 160;
    int filenumber;
    double noise = pow(10,-5);

    /*calculate uncertainty for certain configuration*/
    for (filenumber=1; filenumber<1; filenumber++){ 
        printf("evaluating file %d \n", filenumber);
        calcLLH(mass, distance, events, triggEff, energyRes, filenumber, noise);
    }

    printf("DONE\n");
    return 0;
}
