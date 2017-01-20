create a " .i" file that contains the name of your module, includes the header file & all the functions you want to use
 'carrays.i' and the function below let you use array

then run:
swig -python name.i

this creates a wraper file
then you run setup.py in which you include your sources, so the wraper file and all the .c programs you need, then you run:

python setup.py build_ext --inplace

and thats it, now you can use your module as a normal python module (include modulename)
