We compare the output of the full software. To compare the output for numerically accuracy, 
we need to use the same random number generator. To do so, mobcal is made into a python module
and the random number generator is replace by the python code.
Note that now we run *two*  python codes independently. One is PyMobCal. 
The other is a python code that sets the random number generator and then runs mobcal