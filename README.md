# Likelihood
<b>To run pure C++ code: Compile CUDA code and link:</b><br>
`nvcc -c spectrum.cu` <br>
`g++ -L/usr/local/cuda/lib64/ -lcuda -lcudart spectrum.o generateEvents.c -lm -o generateEvents (-pg)`<br>
<br>
`g++ -L/usr/local/cuda/lib64/ -lcuda -lcudart spectrum.o likelihood.c -lm -o likelihood (-pg)`<br>
<b>Run and profile:</b><br>
`(nvprof) ./generateEvents` (or ./likelihood)<br>
`gprof generateEvents/likelihood gmon.out`<br>

<b>Interface C++ and CUDA code with Python:</b><br>
`nvcc -Xcompiler -fPIC -shared -o spectrum.so spectrum.cu`<br>
(then execute `calc_minimal_likelihood.py`)
