# PyMobCal
Software for calculation of collisional cross sections.

This project is a rewritten Python and C implementation based on MobCal, whose original Fortran 77 code, manual, examples, and algorithm descriptions are available at https://www.indiana.edu/~nano/software/. If you use this software, please cite the original papers (see below). Because the original Fortran 77 implementation is widely used but is not optimized for newer machines, I decided to rewrite it in more modern languages while keeping the published algorithms and workflow intact. At this stage, the code is intended to be 100% compatible with the original Fortran 77 version, preserving the original input and output file formats.

Two versions are currently available—one in Python and one in C. See the manual inside each version's folder for details. Although this is still a development release, the code should be relatively stable. Please compare your results against the original MobCal and let me know how they match up. Sharing your input/output files would also be greatly appreciated, as they will help expand the sample set and build a future test suite.


# Some more technical notes 
I have decided to keep the code as close to the original as possible so people familiar with the Fortran 77 version can follow it easily. Changes from the original implementation are marked in the source code with the keyword DIFF. Where possible, I have replaced built-in functions with standard library equivalents, including the random number generator and the rotation routine.

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
