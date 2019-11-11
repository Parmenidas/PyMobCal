REQUIRED: 
1) you need to have the GNU Scientific Library (GSL) installed. Otherwise the installation will fail.
2) You need to have installed and operational a MPI framework.

To compile:

make cmobcal

==========
IMPORTANT:
==========
*Be careful how you set INUM,ITN and the number of nodes.*

Both the total number of cycles (TM method) and number of Monte Carlo trajectories (EHS, PA method) depend  on the number of cores/nodes.
In addition to the number of cores/nodes, the number of cycles depends on the parameter ITN and the number of MC trajectories depends on INUM.

Thus, for *one* core:
Number of cycles = ITN = N_TM
MC trajectories = INUM = N_MC 

If you use *N_CORES* cores these formulas become:
Number of cycles = ITN * N_CORES = N_TM
MC trajectories = INUM * N_CORES = N_MC 

You may have two scenarios:
1) Your statistics is not enough, so you want more points, you may increase both ITN,INUM and NPROCS
2) Your statistics is enough, but it takes too long. In this case you want to keep the N_TM and N_MC contant.
   For instance, if you move from one to N_CORES, then you replace:
   ITN   ==> ITN/N_CORES
   INUM  ==> INUM/N_CORES
   In this way, N_TM, N_MC do not change, but the same work is shared among N_CORES, so in principle it is N_CORES times faster.
