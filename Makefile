default: generateEvents likelihood

likelihood: likelihood.o spectrum.o
	gcc -o likelihood likelihood.o spectrum.o -lm

likelihood.o: likelihood.c spectrum.h
	gcc -c likelihood.c

generateEvents: generateEvents.o spectrum.o
	gcc -o generateEvents generateEvents.o spectrum.o -lm

generateEvents.o: generateEvents.c spectrum.h
	gcc -c generateEvents.c

spectrum.o: spectrum.c spectrum.h
	gcc -c spectrum.c

clean:
	rm generateEvents.o generateEvents generateEvents spectrum.o spectrum.o likelihood.o likelihood likelihood spectrum*.txt events*.txt
