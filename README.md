# PyMobCal
Software for calculation of collisional cross sections
This software is a "translation" of the software MobCal that can be found at https://www.indiana.edu/~nano/software/. 
That website contains a manual as well as examples. Remember to cite the original papers if you are using this software. 
The original software was written in Fortran 77. As this is a widely used software, 
I thought to convert it into a more modern language such as Python. 
At the current stage, the software should be 100% compatible with the original version in Fortran 77 (same input/output files, etc.). 
Thus, if you have used mobcal already, just copy all the files in one folder and type 'python PyMobCal' 
and wait for magic (or bugs) to happen! Currently, the code needs huge testing. 
Thus, always compare its results with the original mobcal and let me know how it goes. Some more technical notes. 
I have decided to keep the code as close to the original as possible so that people with experience in the Fortran77 version 
may understand it. This should also help with initial debug. This choice is not painless as it hamper performance a bit. 
However, at the current stage, I want to make sure that the code is stable before further optimization. 
Changes from the original version are reported in the source code with the keyword DIFF. 
When possible I have replaced built-in function with standard python libraries. 
Examples are the random number generator and the rotation routine.
Examples, reference and more can be found at https://www.indiana.edu/~nano/software/
