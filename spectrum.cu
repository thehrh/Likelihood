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

// Define this to turn on error checking
#define CUDA_ERROR_CHECK

#define CudaSafeCall( err ) __cudaSafeCall( err, __FILE__, __LINE__ )
#define CudaCheckError()    __cudaCheckError( __FILE__, __LINE__ )

inline void __cudaSafeCall( cudaError err, const char *file, const int line )
{
#ifdef CUDA_ERROR_CHECK
    if ( cudaSuccess != err )
    {
        fprintf( stderr, "cudaSafeCall() failed at %s:%i : %s\n",
                 file, line, cudaGetErrorString( err ) );
        exit( -1 );
    }
#endif

    return;
}

inline void __cudaCheckError( const char *file, const int line )
{
#ifdef CUDA_ERROR_CHECK
    cudaError err = cudaGetLastError();
    if ( cudaSuccess != err )
    {
        fprintf( stderr, "cudaCheckError() failed at %s:%i : %s\n",
                 file, line, cudaGetErrorString( err ) );
        exit( -1 );
    }

    // More careful checking. However, this will affect performance.
    // Comment away if needed.
    err = cudaDeviceSynchronize();
    if( cudaSuccess != err )
    {
        fprintf( stderr, "cudaCheckError() with sync failed at %s:%i : %s\n",
                 file, line, cudaGetErrorString( err ) );
        exit( -1 );
    }
#endif

    return;
}


/* time shift due to neutrino mass */
user_data_t getDeltaT(user_data_t E, user_data_t mass, user_data_t dist){
    user_data_t tDelta = dist*51.4635*(mass/E)*(mass/E);
    return tDelta;
}

user_data_t getTimeDelay(user_data_t t, user_data_t E, user_data_t mass, user_data_t dist){
    return t - getDeltaT(E, mass, dist);
}

user_data_t LL_time_spectrum_shifted(user_data_t t, user_data_t E, user_data_t mass, user_data_t dist){
    user_data_t time = getTimeDelay(t, E, mass, dist);
    if (time <= 0){
        // unphysical?
        return 0.0;
    }
    return LL_time_spectrum(time);
}

/* arrival time probability for a certain mass/distance - normalized */
user_data_t LLSpectrumTotal (user_data_t t, user_data_t E, user_data_t mass, user_data_t dist){
    return LL_time_spectrum_shifted(t, E, mass, dist)*LL_energy_spectrum(E);
}

void cumSumT(user_data_t *arrayToSum, user_data_t *cumulative){
    /*calculate the cumulative sum of the arrival time distribution*/
    int k, l;
    user_data_t cum;
    for (k = 0; k < REST; k++){
        cum = 0.0;
        for (l = 0; l <= k; l++){
            cum += arrayToSum[l];
        }
        cumulative[k] = cum;
    }
}

void firstHitDistWeightedArrivalTimeDist(user_data_t *arrivalTimeDist, user_data_t *cumulative, user_data_t events, user_data_t *result){
    int m;
    user_data_t count = 0.0;
    for (m = 0; m < REST; m++){
        result[m] = arrivalTimeDist[m]*events*pow((1 - cumulative[m]), events-1);
        count += result[m];
    }
    for (m = 0; m < REST; m++){
        result[m] = result[m]/count;
    }
}

/* calculate the probability to get the first hit after a certain amount of time */
void ProbFirstHitDist (user_data_t mass, user_data_t dist, user_data_t events, user_data_t *result){
    /*arrival time distribution of all the hits (for a certain mass) - project the E,t spectrum
    on the t axis - t in 0.01 steps from 0 to 10 seconds*/
    user_data_t totalArrivalTimeDist[REST];
    int i;
    user_data_t sum;
    user_data_t y, e;
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

    user_data_t cumulative[REST];
    cumSumT(totalArrivalTimeDist, cumulative);

    firstHitDistWeightedArrivalTimeDist(totalArrivalTimeDist, cumulative, events, result);
}

void convolveHitDistWithLLTimeSpec(user_data_t *hitDist, user_data_t *convolSpec){
    int i, j;
    user_data_t pNew;
    /*perform the convolution*/
    for (i = 0; i < REST*1.3; i++){
        pNew = 0.0;
        for (j = 0; j < REST; j++){
            if ((i-0.3*REST + j) < REST && (i-0.3*REST + j) > 0){
                pNew += hitDist[j] * LL_time_spectrum( (j+i-0.3*REST)*STEPT );
            }
        }
        convolSpec[i] = pNew;
    }
}

/*calculate the correlation - new spectrum between -3 and 10s*/
/*this is stored in an array so newSpec[0] corresponds to a time of -3s
and newSpec[1.3*REST-1] to 10s*/
void correlation(user_data_t mass, user_data_t dist, user_data_t events, user_data_t *newSpec){
    user_data_t hitDist[REST];
    ProbFirstHitDist(mass, dist, events, hitDist);
    convolveHitDistWithLLTimeSpec(hitDist, newSpec);
}

void applyEnergyRes(int t, user_data_t *distribution, user_data_t *energySpectrum){
    int f, g;
    for (f=1; f<RESE; f+=1){
        user_data_t pNew = 0.0;
        for (g=-RESE; g<RESE+1; g+=5){
            if (f-g >= 0 && f-g <= RESE){
                pNew += GAUSS(g*STEPE, f*STEPE)*energySpectrum[f-g];
            }
            distribution[t*(RESE-1)+f-1] = pNew*STEPE;
        }
    }
}


void normalize(user_data_t *distribution){
	// normalize the spectrum to 1
	int k;
	user_data_t normalize = 0;

	for (k=0; k<(RESE-1)*REST; k++){
		normalize += distribution[k]*STEPT*STEPE;
	}

	for (k=0; k<(RESE-1)*REST; k++){
		distribution[k] *= 1.0/normalize;
    	}
}


__global__
void getEnergySpec(user_data_t mass, user_data_t dist, user_data_t *timeArray, user_data_t *triggerEffs, user_data_t *distribution, bool useEnergyRes){
    user_data_t time, pUnsmeared, pNew, p_E_LL, p_t_LL, triggerEff;
    int e, f, g, arrayIndex;
    int t = blockIdx.x*blockDim.x + threadIdx.x;
    if (t < REST){
        user_data_t energySpectrum[RESE];
        energySpectrum[0] = 0.0;
        for (e=1; e<RESE; e++){
            // make this explicit for now, until know how to call function correctly
            time =  t*STEPT - dist*51.4635*(mass/(e*STEPE))*(mass/(e*STEPE));//getTimeDelay(t*STEPT, e*STEPE, mass, dist);
            arrayIndex = (int) (time/(STEPT) + 0.3*REST);
            if (arrayIndex <= 0){
                arrayIndex = 0;
	        }
            p_E_LL = LL_energy_spectrum(e*STEPE);
            p_t_LL = timeArray[arrayIndex];
            triggerEff = triggerEffs[e];
            /*
            if (t==0){
                printf("\n LL energy spectrum entry %d: %.10f \n", e, p_E_LL);
                printf("\n trigger Eff %d: %.10f \n", e, triggerEff);
                printf("\n timeArray entry %d: %.10f \n", arrayIndex, p_t_LL);
            }
            */
	        pUnsmeared = p_E_LL*p_t_LL*triggerEff;
            if (!useEnergyRes){
                distribution[t*(RESE-1) +e-1] = pUnsmeared;
            }
            energySpectrum[e] = pUnsmeared;
	    }
        if (useEnergyRes){
            // also make this explicit for now
            // applyEnergyRes(t, distribution, energySpectrum);
            for (f=1; f<RESE; f+=1){
                /*
                if (t==0){
                    printf("\n For t=0, energy spectrum entry %d: %.10f \n", f, energySpectrum[f]);
                }
                */
                pNew = 0.0;
                for (g=-RESE; g<RESE+1; g+=5){
                    if (f-g >= 0 && f-g < RESE){
                        pNew += GAUSS(g*STEPE, f*STEPE)*energySpectrum[f-g];
		            }
                }
                /*
                if (t==0){
                    printf("\n For t=0, write %.10f to %d \n", pNew*STEPE, t*(RESE-1)+f-1);
                }
                */
                distribution[t*(RESE-1)+f-1] = pNew*STEPE;
            }
	    }
        /*
        if (t==0){
            printf("\n spectrum[120]: %.10f", distribution[120]);
        }
        */
    }	
}

void generateDist(user_data_t mass, user_data_t dist, user_data_t events, user_data_t *distribution, user_data_t *triggerEffs, bool useEnergyRes){
    user_data_t timeArray[(int) (1.3*REST)];
    user_data_t *d_triggerEffs, *d_distribution, *d_timeArray;
    int size = sizeof(user_data_t);

    correlation(mass, dist, events, timeArray);

    /*
    //create a file from the timeArray for debugging
    char filename[sizeof "timeArray_CUDA.txt"];
    sprintf(filename, "timeArray_CUDA.txt");
    FILE *f = fopen(filename, "w+");
    for(int i=0; i<(int)(1.3*REST); i++){
        fprintf(f, "%e\n", timeArray[i]);
    }
    fclose(f);
    */

    cudaMalloc((void **)&d_distribution, (RESE-1) * REST * size);
    cudaMalloc((void **)&d_triggerEffs, (RESE+1)*size);
    cudaMalloc((void **)&d_timeArray, 1.3*REST*size);

    cudaMemcpy(d_timeArray, &timeArray, 1.3*REST*size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_triggerEffs, triggerEffs, (RESE+1)*size, cudaMemcpyHostToDevice);

    getEnergySpec<<<(REST + 511) / 512, 512>>>(mass, dist, d_timeArray, d_triggerEffs, d_distribution, useEnergyRes);
    CudaCheckError();

    cudaMemcpy(distribution, d_distribution, (RESE-1) * REST * size, cudaMemcpyDeviceToHost);

    /*
    //create a file from the dist before norm
    char filename2[sizeof "spec_before_norm_CUDA.txt"];
    sprintf(filename2, "spec_before_norm_CUDA.txt");
    FILE *f2 = fopen(filename2, "w+");
    for(int i=0; i<(RESE-1)*REST; i++){
        fprintf(f2, "%e\n", distribution[i]);
    }
    fclose(f2);
    */

    cudaFree(d_timeArray);
    cudaFree(d_distribution);
    cudaFree(d_triggerEffs);

    normalize(distribution);
}
