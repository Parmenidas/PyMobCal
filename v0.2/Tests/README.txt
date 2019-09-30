In this folder we compare results between PyMobCal and the orginal Mobcal.


SINGLE FUNCTION COMPARISON

The following folders compare the functions with the same name:

che
deriv
diffeq
dljpot
fcoord
gsang
mobil2
mobil4
ncoord
rantate
rotate

In this case the fortran code is compiled using f2py and then called by the python script.
The script sends the *same* data to both the fortran and python function and compares the
output. The .pyf files needed to compile the fortran module are included


EXTENDED OUTPUT

This folder compares the full programs with verbose output.
Mobcal is converted into a python module and driven by python.
The reason for that is to use the exact same random number generator
so that we can compare numerically the output of both codes.
