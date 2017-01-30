/* 
*Author: Maike Jung
*Date: 2.11.2016

*Purpose: Calculate the deltaLLH for the different mass spectra to determine how well a certain neutrino mass (with a certain number of events at a certain distance) could be detected 

SN - Model: Lawrence-Livermore
    time spectrum is convoluted with the first hit distribution, to account for not knowing the absolute arrival times

UNITS: mass: eV
       energy: MeV
       distance: Mpc
       time: s
*/

#include "spectrum.h"
void calculateLLH(double mass, double distance, double events, double stepSize, bool triggEff, bool energyRes, bool checkEntireSpectrum, double noise){

    /*read in trigger efficiency*/
    int i;
    double triggerEnergy[RESE+1];
    double triggerEfficiency[RESE+1];
    /* initialize with 1 in case triggerEfficiency is not used*/
    for(i = 0; i < RESE+1 ; triggerEfficiency[i++] = 1);
    if(triggEff){
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
            fscanf(myFile, "%lf %lf", &triggerEnergy[i], &triggerEfficiency[i]);
        }
        fclose(myFile);
    }

	double *myTrue= (double*) malloc((RESE-1) * REST * sizeof(double));
	double *myTest= (double*) malloc((RESE-1) * REST * sizeof(double));
    /*calculate my true*/
	generateDist(mass, distance, events, myTrue, triggerEfficiency, energyRes);

    /*calculate my test and calculate DeltaLLH for the different masses*/
	bool foundUpper = true;
	bool foundLower = true;
	double upperLimit, lowerLimit;
	int k = 0;
	double deltaLLHi;
	double testMass = mass;
	double deltaLLH[REST]; /*array in which results are stored - 1000 should be large enough but check*/
	
    /*calculate the likelihood over a larger range*/
    if (checkEntireSpectrum) {
        char filenames[sizeof "binnedLLH_Spectrum_pos_1.5eV_ideal.txt"];
        sprintf(filenames, "binnedLLH_Spectrum_neg_%feV_real3.txt", mass);
        FILE *f = fopen(filenames, "a+");
        if (f == NULL) {
            printf("Error opening file!\n");
            exit(1);
        }
        printf("looking at entire range \n");

        for (testMass = mass; testMass > mass-0.5; testMass -= 0.01){
            deltaLLHi = 0.0;
            double *myTest= (double*) malloc((RESE-1) * REST * sizeof(double));
		    generateDist(testMass, distance, events, myTest, triggerEfficiency, energyRes);
            /*calculate DLLH*/
            for (i = 1; i < REST*(RESE-1); i++){
                double mtrue = myTrue[i]*(10.0/REST)*(60.0/RESE)*events + noise; /*calculate area for each bin*/
                double mtest = myTest[i]*(10.0/REST)*(60.0/RESE)*events + noise;
                deltaLLHi += fabs(mtest - mtrue*log(mtest) - mtrue + mtrue*log(mtrue));
                //printf("deltaLLHi: %d %e %e %e\n", i, mtrue, mtest, deltaLLHi);
            }
            fprintf(f, "%f %f\n", testMass, deltaLLHi);
		    printf("mass %f, deltaLLH %f\n", testMass, deltaLLHi);
            free(myTest);
        }
        fclose(f);
    }

    if (checkEntireSpectrum) {
        char filenames[sizeof "binnedLLH_Spectrum_pos_1.5eV_ideal.txt"];
        sprintf(filenames, "binnedLLH_Spectrum_pos_%feV_real3.txt", mass);
        FILE *f = fopen(filenames, "a+");
        if (f == NULL) {
            printf("Error opening file!\n");
            exit(1);
        }
        printf("looking at entire range \n");
        double mtrue;
        double mtest;
        for (testMass = mass; testMass < mass+0.5; testMass += 0.01){
            double *myTest= (double*) malloc((RESE-1) * REST * sizeof(double));
            deltaLLHi = 0.0;
		    generateDist(testMass, distance, events, myTest, triggerEfficiency, energyRes);
            /*calculate DLLH*/
            for (i = 1; i < REST*(RESE-1); i++){
                if(myTrue[i]> 1 || myTest[i] >1 ) printf("h %e %e\n", myTrue[i], myTest[i]);
                mtrue = myTrue[i]*(10.0/REST)*(60.0/RESE)*events + noise; /*calculate area for each bin*/
                mtest = myTest[i]*(10.0/REST)*(60.0/RESE)*events + noise;
                if(mtrue> 1 || mtest >1 || mtrue< pow(10,-8) || mtest <pow(10,-8) ){ 
                    printf("h %e %e\n", mtrue, mtest);
                }
                deltaLLHi += fabs(mtest - mtrue*log(mtest) - mtrue + mtrue*log(mtrue) ); 
                //deltaLLHi += fabs(mtest - mtrue*log(mtest) - mtrue + mtrue*log(mtrue)); 
            }
            fprintf(f, "%f %f\n", testMass, deltaLLHi);
		    printf("mass %f, deltaLLH %f\n", testMass, deltaLLHi);
            free(myTest);
        }
        fclose(f);
    }

    foundUpper = false;
    foundLower = false;
	/*determine only 1sigma uncertainty*/
	/*search for upper and lower limit*/
	while (foundUpper){
		deltaLLHi = 0.0;
		generateDist(testMass, distance, events, myTest, triggerEfficiency, energyRes);

		/*calculate DLLH*/
		for (i = 1; i < REST*(RESE-1); i++){
            double mtrue = myTrue[i]*(10.0/REST)*(60.0/RESE)*events + noise; /*calculate area for each bin*/
            double mtest = myTest[i]*(10.0/REST)*(60.0/RESE)*events + noise;
			deltaLLHi += fabs(mtest - mtrue*log(mtest) - mtrue + mtrue*log(mtrue)); 
		
        }
		deltaLLH[k] = deltaLLHi;
		printf("mass %f, deltaLLH %f\n", testMass, deltaLLHi);
		if (deltaLLHi > 0.5){
			foundUpper = false;
            /*interpolate between the last two values and return the mass for which DLLH reaches 0.5*/
            double mass1 = testMass - stepSize;
            double mass2 = testMass;
            double value1 = deltaLLH[k-1];
            double value2 = deltaLLH[k];
            double a = (value1 - value2)/(mass1 - mass2);
            double b = value1 - a*mass1;
            upperLimit = (0.5 - b)/a;
            printf("upper limit %f \n", upperLimit);
		}
		testMass += stepSize;
        k++;
	}
	testMass = mass - stepSize;

	while (foundLower){
        /*if test mass is to close to 0, return 0*/
        if (testMass < stepSize) {
            foundLower = false;
            lowerLimit = 0.0;
            break;
        }
        deltaLLHi = 0.0;
		generateDist(testMass, distance, events, myTest, triggerEfficiency, energyRes);

		/*calculate DLLH*/
		for (i = 1; i < REST*(RESE-1); i++){
            double mtrue = myTrue[i]*(1/(REST*0.1))*(1/(RESE/60.0))*events; /*calculate area for each bin*/
            double mtest = myTest[i]*(1/(REST*0.1))*(1/(RESE/60.0))*events;
			deltaLLHi += fabs(mtest - mtrue*log(mtest) - mtrue + mtrue*log(mtrue)); 	
        }
		deltaLLH[k] = deltaLLHi;
		printf("mass %f, deltaLLH %f\n", testMass, deltaLLHi);
		if (deltaLLHi > 0.5){
			foundLower = false;
            /*interpolate between the last two values and return the mass for which DLLH reaches 0.5*/
            double mass1 = testMass + stepSize;
            double mass2 = testMass;
            double value1 = deltaLLH[k-1];
            double value2 = deltaLLH[k];
            double a = (value1 - value2)/(mass1 - mass2);
            double b = value1 - a*mass1;
            lowerLimit = (0.5 - b)/a;
            printf("lower limit %f \n", lowerLimit);
		}
		testMass -= stepSize;
        k++;
	}
	free(myTrue);
	free(myTest);
    if(!checkEntireSpectrum){
        printf("input mass: %f eV \n", mass);
        printf("mass found: %f - %f + %f eV\n", mass, mass-lowerLimit, upperLimit-mass);

	    /*storing mass and uncertainty in file*/
        char filename[sizeof "massUncertainty_10.5Mpc_1000Events_DLLH_trigg_eff_energy_res.txt"];
        if (triggEff && energyRes){
            sprintf(filename, "massUncertainty_%.1fMpc_%.0fEvents_DLLH_trigg_eff_energy_res_noise.txt", distance, events);
        }
        else if (energyRes){
            sprintf(filename, "massUncertainty_%.1fMpc_%.0fEvents_DLLH_energy_res.txt", distance, events);
        }
        else if (triggEff){
            sprintf(filename, "massUncertainty_%.1fMpc_%.0fEvents_DLLH_trigg_eff.txt", distance, events);
        } 
        else {
            sprintf(filename, "massUncertainty_%.1fMpc_%.0fEvents_DLLH.txt", distance, events);
        }
        FILE *f = fopen(filename, "a+");
        if (f == NULL)
        {
            printf("Error opening file!\n");
            exit(1);
        }
        fprintf(f, "%f %f %f\n", mass, mass-lowerLimit, upperLimit-mass);
        fclose(f);
    }
}

int main(void){
	/*set parameters*/
    // if checkEntireSpectrum = true: calulcate likelihood over a large range of time
    // else: checke when likelihood reaches 0.5 to find the 1sigma limits
    bool triggEff = true;
    bool energyRes = true;
    bool checkEntireSpectrum = true;
    double mass;
    double distance = 5.0;
    double events = 10;
    double stepSize = 0.01;
    double noise = pow(10,-5);
    
    /*calculate uncertainty for certain configuration*/
    for (mass = 1.0; mass < 1.6; mass+=1.1){
        calculateLLH(mass, distance, events, stepSize, triggEff, energyRes, checkEntireSpectrum, noise);
    }
    printf("DONE\n");
    return 0;
}
