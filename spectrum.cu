/*
*Author: Maike Jung
*Date: 15.11.2016

*Purpose: create the arrival time spectrum of the neutrinos, that can then be used to
    generate random events: generateEvents.c
    calculate the likelihood for these events:  likelihood.c
    calculate the binned-likelihood:    binned_likelihood.c

SN - Model: Lawrence-Livermore
    time spectrum is convoluted with the first hit distribution, to account for not knowing the absolute arrival times

UNITS: mass: eV
       energy: MeV
       distance: Mpc
       time: s

add noise of 10-5!
*/

#include "spectrum.h"


/* time shift due to neutrino mass */
double getDeltaT(double E, double mass, double dist){
    double tDelta = dist*51.4635*(mass/E)*(mass/E);
    return tDelta;
}

double getTimeDelay(double t, double E, double mass, double dist){
    return t - getDeltaT(E, mass, dist);
}

double LL_time_spectrum_shifted(double t, double E, double mass, double dist){
    double time = getTimeDelay(t, E, mass, dist);
    if (time <= 0){
        // unphysical?
        return 0.0;
    }
    return LL_time_spectrum(time);
}

/* arrival time probability for a certain mass/distance - normalized */
double LLSpectrumTotal (double t, double E, double mass, double dist){
    return LL_time_spectrum_shifted(t, E, mass, dist)*LL_energy_spectrum(E);
}

void cumSumT(double *arrayToSum, double *cumulative){
    /*calculate the cumulative sum of the arrival time distribution*/
    int k, l;
    double cum;
    for (k = 0; k < REST; k++){
        cum = 0.0;
        for (l = 0; l <= k; l++){
            cum += arrayToSum[l];
        }
        cumulative[k] = cum;
    }
}

void firstHitDistWeightedArrivalTimeDist(double *arrivalTimeDist, double *cumulative, double events, double *result){
    int m;
    double count = 0.0;
    for (m = 0; m < REST; m++){
        result[m] = arrivalTimeDist[m]*events*pow((1 - cumulative[m]), events-1);
        count += result[m];
    }
    for (m = 0; m < REST; m++){
        result[m] = result[m]/count;
    }
}

/* calculate the probability to get the first hit after a certain amount of time */
void ProbFirstHitDist (double mass, double dist, double events, double *result){
    /*arrival time distribution of all the hits (for a certain mass) - project the E,t spectrum
    on the t axis - t in 0.01 steps from 0 to 10 seconds*/
    double totalArrivalTimeDist[REST];
    int i;
    double sum;
    double y, e;
    for (i = 0; i < REST; i++){
        /* set the sum to zero for each time bin */
        sum = 0.0;
        /*Integrate over the energy part for every time bin. We move in 0.01 MeV
        steps up to 60 MeV. For each pair of time and energy, we compute the
        product of time and energy PDF ("LLSpectrumTotal"), continually 
        incrementing the sum*/
        for (e = 0.01; e < EMAX; e += 0.01) {
            y = LLSpectrumTotal(i*STEPT, e, mass, dist);
            sum += y * 0.01;
        }
        totalArrivalTimeDist[i] = sum*STEPT;
    }

    double cumulative[REST];
    cumSumT(totalArrivalTimeDist, cumulative);

    firstHitDistWeightedArrivalTimeDist(totalArrivalTimeDist, cumulative, events, result);
}

void convolveHitDistWithLLTimeSpec(double *hitDist, double *convolSpec){
    int i, j;
    double pNew;
    /*perform the convolution*/
    for (i = 0; i < REST*1.3; i++){
        pNew = 0.0;
        for (j = 0; j < REST; j++){
            if ((i-0.3*REST + j) < REST && (i-0.3*REST + j) > 0){
                pNew += hitDist[j] * LL_time_spectrum( (j+i-0.3*REST)*STEPT );
            }
        convolSpec[i] = pNew;
        }
    }
}

/*calculate the correlation - new spectrum between -3 and 10s*/
/*this is stored in an array so newSpec[0] corresponds to a time of -3s
and newSpec[1.3*REST-1] to 10s*/
void correlation(double mass, double dist, double events, double *newSpec){
    double hitDist[REST];
    ProbFirstHitDist(mass, dist, events, hitDist);
    convolveHitDistWithLLTimeSpec(hitDist, newSpec);
}

void applyEnergyRes(int t, double *distribution, double *energySpectrum){
    int f, g;
    for (f=1; f<RESE; f+=1){
        double pNew = 0.0;
        for (g=-RESE; g<RESE+1; g+=5){
            if (f-g >= 0 && f-g <= RESE){
                pNew += GAUSS(g*STEPE, f*STEPE)*energySpectrum[f-g];
            }
            distribution[t*(RESE-1)+f-1] = pNew*STEPE;
        }
    }
}

/* This becomes the function we want to run in parallel on the GPU */
__global__
void getEnergySpec(double *mass, double *dist, double *timeArray, double *distribution, double *triggerEffs, bool *useEnergyRes){
	double time, pUnsmeared;
	int e, f, g, arrayIndex;

    int t = blockIdx.x;//*blockDim.x + threadIdx.x;
    if (t < REST ){
        double energySpectrum[RESE];
        energySpectrum[0] = 0.0;
        for (e=1; e<RESE; e++){
            // make this explicit for now, until know how to call function correctly
            time =  t*STEPT - (*dist)*51.4635*(*mass/e*STEPE)*(*mass/e*STEPE);//getTimeDelay(t*STEPT, e*STEPE, mass, dist);
            arrayIndex = (int) (time/(STEPT) + 0.3*REST);
            if (arrayIndex <= 0){
                arrayIndex = 0;
            }
            pUnsmeared = LL_energy_spectrum(e*STEPE)*timeArray[arrayIndex]*triggerEffs[e];
            if (!useEnergyRes){
                distribution[t*(RESE-1) +e-1] = pUnsmeared;
            }
            energySpectrum[e] = pUnsmeared;
        }
        if (useEnergyRes){
            // also make this explicit for now
            // applyEnergyRes(t, distribution, energySpectrum);
            for (f=1; f<RESE; f+=1){
                double pNew = 0.0;
                for (g=-RESE; g<RESE+1; g+=5){
                    if (f-g >= 0 && f-g <= RESE){
                        pNew += GAUSS(g*STEPE, f*STEPE)*energySpectrum[f-g];
                    }
                distribution[t*(RESE-1)+f-1] = pNew*STEPE;
                }
            }
        }
    }
}

void normalize(double *distribution){
	// normalize the spectrum to 1
	int k;
	double normalize = 0;

	for (k=0; k<(RESE-1)*REST; k++){
		normalize += distribution[k]*STEPT*STEPE;
	}

	for (k=0; k<(RESE-1)*REST; k++){
        distribution[k] *= 1.0/normalize;
    	}

}

/*generate the proper distribution*/
void generateDist(double mass, double dist, double events, double *distribution, double *triggerEffs, bool useEnergyRes){
    double *timeArray, *d_timeArray;
    timeArray = (double*)malloc((1.3*REST)*sizeof(double));
    double *d_mass, *d_dist, *d_distribution, *d_triggerEffs;
    bool *d_useEnergyRes;
    int size = sizeof(double);

    //double timeArray[(int) (1.3*REST)];
    correlation(mass, dist, events, timeArray);
    cudaMalloc(&d_mass, size); cudaMalloc(&d_dist, size);
    cudaMalloc(&d_distribution, (RESE-1) * REST * size);
    cudaMalloc(&d_triggerEffs, RESE*size);
    cudaMalloc(&d_useEnergyRes, sizeof(bool));
    cudaMalloc(&d_timeArray, (1.3*REST)*sizeof(double));

    cudaMemcpy(d_mass, &mass, size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_dist, &dist, size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_timeArray, &timeArray, 1.3*REST*size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_distribution, &distribution, (RESE-1) * REST * size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_triggerEffs, &triggerEffs, RESE*size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_useEnergyRes, &useEnergyRes, sizeof(bool), cudaMemcpyHostToDevice);

    getEnergySpec<<< REST, 1 >>>(d_mass, d_dist, d_timeArray, d_distribution, d_triggerEffs, d_useEnergyRes);

    cudaMemcpy(&distribution, d_distribution, (RESE-1) * REST * size, cudaMemcpyDeviceToHost);

    cudaFree(d_mass); cudaFree(d_dist); cudaFree(d_timeArray); cudaFree(d_distribution);
    cudaFree(d_triggerEffs); cudaFree(d_useEnergyRes);

	normalize(distribution);
}
