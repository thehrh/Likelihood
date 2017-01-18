profile: generateEventsPG likelihoodPG

default: generateEvents likelihood

likelihood: likelihood.o spectrum.o
	gcc -o likelihood likelihood.o spectrum.o -lm

likelihoodPG: likelihood.o spectrumPG.o
	gcc -o likelihoodPG likelihood.o spectrumPG.o -lm -pg

likelihood.o: likelihood.c spectrum.h
	gcc -c likelihood.c

generateEvents: generateEvents.o spectrum.o
	gcc -o generateEvents generateEvents.o spectrum.o -lm

generateEventsPG: generateEvents.o spectrumPG.o
	gcc -o generateEventsPG generateEvents.o spectrumPG.o -lm -pg

generateEvents.o: generateEvents.c spectrum.h
	gcc -c generateEvents.c

spectrumPG.o: spectrum.c spectrum.h
	gcc -c spectrum.c -o spectrumPG.o -pg

spectrum.o: spectrum.c spectrum.h
	gcc -c spectrum.c

clean:
	rm generateEvents.o generateEvents generateEventsPG spectrum.o spectrumPG.o likelihood.o likelihood likelihoodPG spectrum*.txt events*.txt
