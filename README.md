# Likelihood
<b>Compile CUDA code and link:</b><br>
nvcc -c spectrum.cu <br>
g++ -L/usr/local/cuda/lib64/ -lcuda -lcudart spectrum.o generateEvents.c -lm -o generateEvents (-pg)<br>
<b>Run and profile:</b><br>
(nvprof) ./generateEvents<br>
gprof generateEvents gmon.out
