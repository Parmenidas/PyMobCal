# PyMobCal
Software for calculation of collisional cross sections.

This software is a "translation" of the software MobCal that can be found at https://www.indiana.edu/~nano/software/. 
That website contains a manual as well as examples. Remember to cite the original papers if you are using this software (see below). 
The original software was written in Fortran 77. As this is a widely used software, 
I thought to convert it into a more modern language such as Python or C. 
At the current stage, the software should be 100% compatible with the original version in Fortran 77 (same input/output files, etc.). 

There are currently two versions of the code. One written in Python and one in C. See the manual within each version folder for details. Currently, this is still a development version, but the code should be relatively stable at this stage. In any case, always compare its results with the original mobcal and let me know how it goes.
If you would like to share your files, it would be appreciated. They will be added to the samples and to a future test suite.

# Some more technical notes 
I have decided to keep the code as close to the original as possible so that people with experience in the Fortran77 version 
may understand it. Changes from the original version are reported in the source code with the keyword DIFF. 
When possible I have replaced built-in function with standard libraries. 
Examples are the random number generator and the rotation routine.

# References
When using this software, be sure to cite the following papers:

* M. F. Mesleh, J. M. Hunter, A. A. Shvartsburg, G. C. Schatz, and M. F. Jarrold, Structural Information from Ion Mobility Measurements: Effects of the Long Range Potential, J. Phys. Chem. 1996, 100, 16082-16086; Erratum, J. Phys. Chem. A 1997, 101, 968.

* A. A. Shvartsburg and M. F. Jarrold, An Exact Hard Spheres Scattering Model for the Mobilities of Polyatomic Ions, Chem. Phys. Lett. 1996, 261, 86-91.

# Major Versions
Note: Starting from version 3, a manual is available within each folder.
* v0.1 - A brand new MobCal!
* v0.2 - bug fixes
* v0.3 - improvement in the trajectory method. The original algorithm has been replaced by Velocity Verlet. Several optimizations to make the code faster. This version of the code should be used if you want to add/remove/update features in mobcal. Not the most optimized version, but the easiest to modify.
* v0.3c - This is the same code as version 0.3, but written in C. This version should be used once the testing of novel features has been done with the python code. It provides a performance boost.  
* v0.3pthread **(recommended)** - This is the version that users should choose. It takes adavantage of multi-core of modern CPUs and so outperforms the single core version (python and regular c).
* v0.3mpi - This is the version for clusters. Unless you have a huge molecule, it is probably overkill. It also requires a working MPI environment.

# Developing new features
The best way is to check both the manual and the developer's manual within each release. The python version of the code has been designed to be the easiest to modify and is well suited for testing and rapid development. Move to the C version only after testing is done.
