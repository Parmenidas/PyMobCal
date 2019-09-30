'''
Test of function rotate
note that the fortran version of rotate has been changed adding 
open/close the file 

**IMPORTANT**: These function create their own output files that
should be checked for differences (out.fortrant.dat and out.pyhton.dat).
This is in *addition* to diff_file.dat that checks the return value(s)
of those functions
'''

import numpy as np
from scipy.spatial.transform import Rotation as R
import rotate_f

class printSwitch:
    '''
    Some switches used to turn on/off some debugging features.
    Unless you are modifying the code, keep defaults
    '''
    def __init__(self):
        self.ip=0
        self.it=0
        self.iu1=0
        self.iu2=0
        self.iu3=0
        self.iv=0
        self.im2=0
        self.im4=0
        self.igs=0
    
    def SetIt(self,value):
        self.it=value

    def SetIu2(self,value):
        self.iu2=value

    def SetIu3(self,value):
        self.iu3=value
    

class constants:
    '''
    Constants used within the code
    '''
    def __init__(self):
        self.pi=np.pi
        self.cang= 180/np.pi
        self.xe= 1.60217733e-19
        self.xk=1.380658e-23
        self.xn= 6.0221367e23
        self.xeo= 8.854187817e-12
        self.xmv=0.02241410
        self.eo=1.34e-03*self.xe
        self.ro=3.043e-10
        self.ro2=self.ro*self.ro
        self.ipr=1000
        # self.inp=40   
        # self.itn=10 
        # self.imp=25
        self.inp=1   
        self.itn=1 
        self.imp=1
        self.cmin=0.0005
        self.sw1=0.00005
        self.sw2=0.005 
        self.dtsf1=0.5
        self.dtsf2=0.1
        self.inwr=1
        self.ifail=100
        self.ifailc=0              #This should probably be brought outside. It is changed by the code, so not a constant
        self.inum=250000
        self.inor=30
        self.v=2309.9
        self.b=1.5654e-10
        self.ntheta=33.3
        self.nphi=64.250
        self.ngamma=134.30
        self.ehsm=0
        self.tmm=0
        self.immmax=0
        self.immmin=self.inor
        ## I have replace traj_lost with 1000 instead of 30000 just to test the error message in a reasonable amount of time
        self.traj_lost=1000        #DIFF: This was hard-coded into the function. I believe this is the right place
        self.vec_dim=100            #DIF: this is what is used within the FORTRAN function as dimension of arrays 

def rotate(vector,theta,phi,gamma,ConstList,switch,fout):
    '''
    Rotates the vector(s) provided
    '''
    if switch.iu2==1 or switch.iu3==1:  
        print('\n\n coordinates rotated by ROTATE\n\n theta={0: .4E} phi={1: .4E} gamma={2: .4E}\n'.format(theta*ConstList.cang,phi*ConstList.cang,gamma*ConstList.cang),file=fout)

    rz=R.from_rotvec(theta*np.array([0,0,1]))
    rx=R.from_rotvec(-phi*np.array([1,0,0]))
    rzb=R.from_rotvec(gamma*np.array([0,0,1]))
    newVector=rzb.apply(rx.apply(rz.apply(vector)))

    if switch.iu2==1:
        x=' '  
        print(9*x+"initial coordinates"+24*x+"new coordinates\n",file=fout)
        for v,rotv in zip(vector,newVector):
            print(' {: .4E} {: .4E} {: .4E}      {: .4E} {: .4E} {: .4E}'.format(v[0],v[1],v[2],rotv[0],rotv[1],rotv[2]),file=fout)    

    return newVector



##########################################################################
## Set variables. This is the values used by *BOTH* FORTRAN and PYTHON code
##########################################################################
v1=np.array([1,1,1])
v2=np.array([1,-1,1])
np.random.RandomState(seed=43)
n=np.random.rand((3))
theta=n[0]*2.0*np.pi
phi=np.arcsin((n[1]*2.0)-1.0)+(np.pi/2.0)
gamma=n[2]*2.0*np.pi

parmList=constants()
switch=printSwitch()
switch.iu2=1
vector=np.array([v1,v2])
##################################
## Set Fortran variables
##################################
rotate_f.printswitch.iu2=switch.iu2
rotate_f.angles.theta=theta
rotate_f.angles.phi=phi
rotate_f.angles.gamma=gamma
rotate_f.constants.cang=parmList.cang
rotate_f.constants.inatom=2
for c,v in enumerate(vector):
    rotate_f.coordinates.ox[c]=v[0]
    rotate_f.coordinates.oy[c]=v[1]
    rotate_f.coordinates.oz[c]=v[2]
##################################
## Set Python variables
##################################
fout_t_py=open("out.python.dat","w") # file with the result for python
fout_py=open("python_file.dat","w") # file with the result for python
fout_f=open("fortran_file.dat","w") # file with the result for fortran
fout_diff=open("diff_file.dat","w") # file with the difference between the two files above
newVector=rotate(vector,theta,phi,gamma,parmList,switch,fout_t_py)
for v in newVector:
    print(v[0],v[1],v[2],file=fout_py)

rotate_f.rotate()
for c,v in enumerate(newVector):
    print(rotate_f.coordinates.fx[c],rotate_f.coordinates.fy[c],rotate_f.coordinates.fz[c],file=fout_f)
    print(rotate_f.coordinates.fx[c]-v[0],rotate_f.coordinates.fy[c]-v[1],rotate_f.coordinates.fz[c]-v[2],file=fout_diff)

fout_py.close()
fout_t_py.close()
fout_f.close()
fout_diff.close()


