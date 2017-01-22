# Likelihood
nvcc -c spectrum.cu
g++ -L/usr/local/cuda/lib64/ -lcuda -lcudart spectrum.o generateEvents.c -lm -o generateEvents (-pg)
(nvprof) ./generateEvents
gprof generateEvents gmon.out
