'''
Test of function deriv

Here we provide a known potential (see function harmonic)
to both functions and check deriv. We chose the harmonic 
as it is analytically solvable and it is easy to check the output
for mistakes
'''


import deriv 
import numpy as np

class molecule:
    '''
    Contains information about the molecules as a whole
    '''

    def __init__(self):
        self.TotalMass=0
        self.dipole=0
        self.totalCharge=0
        self.totalAbsoluteCharge=0
        self.mGas=0
        self.mu=0  #reduced mass
        self.mobility=0
        self.temperature=float(298)  #kelvin

    def SetMass(self,mass):
        self.TotalMass=mass

    def SetDipole(self,dipole):
        self.dipole=dipole

    def SetMobility(self,mobility):
        self.mobility=mobility

    def SetCharges(self,tcharge,acharge):
        self.totalCharge=tcharge
        self.totalAbsoluteCharge=acharge
    
    def SetMassGas(self,gasMass,scale=1):
        self.mGas=gasMass
        self.mu=((self.TotalMass*self.mGas)/(self.TotalMass+self.mGas))/(scale)


class LennardJones:
    '''
    It contains some variable that are needed for the calculation.
    Instead of being evaluated over and over again, we compute them once and reuse
    '''
    def __init__(self,eolj,rolj):
        self.eolj=eolj
        self.rolj=rolj
        self.eox4=4*eolj
        self.ro2lj=np.square(rolj)
        self.ro6lj=np.power(self.ro2lj,3)
        self.ro12lj=np.square(self.ro6lj)
        self.dro6=6*self.ro6lj
        self.dro12=12*self.ro12lj
        self.romax=np.amax(rolj)


def deriv_py(w,ion,coord,lj,chargeList,dipole):
    '''
c     Defines Hamilton's equations of motion as the time derivatives 
c     of the coordinates and momenta.
    '''
    dw=np.zeros(6)

#     From Hamilton's equations, the time derivatives of the coordinates
#     are the conjugates divided by the mass.
    dw[0]=w[1]/ion.mu
    dw[2]=w[3]/ion.mu
    dw[4]=w[5]/ion.mu

#     Hamilton's equations for the time derivatives of the momenta
#     evaluated by using the coordinate derivatives together with the
#     chain rule.
    x=w[0]
    y=w[2]
    z=w[4]
#    These are analytical derivatives.
    pot,dpotx,dpoty,dpotz,dmax=harmonic(x,y,z)
    dw[1]=-dpotx
    dw[3]=-dpoty
    dw[5]=-dpotz

    return dw,pot,dpotx,dpoty,dpotz,dmax


def harmonic(x,y,z):
    '''
    Harmonic potential. This potential "replaces" the original dljpot in the fortran code.
    We did *NOT* change any fortran functions, but we used the .pyf file to
    play this trick. In brief, in the .pyf file we "explain" that when running the
    code, we should use the function we provide instead of dljpot. Of course, the substitute function 
    must have the same input/output variables passed in the same order.
    '''
    k=1
    pot=-0.5*k*(x**2+y**2+z**2)
    dpotx=k*x
    dpoty=k*y
    dpotz=k*z
    dmax=0
    return pot,dpotx,dpoty,dpotz,dmax



# We need to set some variables. Firstly the general variable used by
# *BOTH* FORTRAN and PYTHON code.
mu=1

## Some variables for the python function.
## All zeros (except for mu) because they are needed for dljpot, 
## but we are not using it in this example
coord=np.array([[0,0,0]])
chargeList=np.array([0])
dipole=0
lj=LennardJones(np.array([0]),np.array([0]))
ion=molecule()
ion.mu=mu

# Here, the FORTRAN variables
deriv.constants.mu=mu

fout_f=open("fortran_file.dat","w") # file with the result for fortran
fout_py=open("python_file.dat","w") # file with the result for python
fout_diff=open("diff_file.dat","w") # file with the difference between the two files above

for t in np.arange(0,2*np.pi,0.1):
    #Harmonic motion:
    x=np.cos(2*t)
    vx=np.sin(2*t)
    y=np.cos(2*t)
    vy=np.sin(2*t)
    z=np.cos(2*t)
    vz=np.sin(2*t)

    w=np.array([x,vx,y,vy,z,vz])
    
    dw,pot,dpotx,dpoty,dpotz,dmax=deriv.deriv(w,harmonic)
    print(' '.join(map(str,dw)),"\t",pot,"\t",dpotx,"\t",dpoty,"\t",dpotz,"\t",dmax,file=fout_f)
    dw2,pot2,dpotx2,dpoty2,dpotz2,dmax2=deriv_py(w,ion,coord,lj,chargeList,dipole)
    print(' '.join(map(str,dw2)),"\t",pot2,"\t",dpotx2,"\t",dpoty2,"\t",dpotz2,"\t",dmax2,file=fout_py)
    print(' '.join(map(str,dw2-dw)),"\t",pot2-pot,"\t",dpotx2-dpotx,"\t",dpoty2-dpoty,"\t",dpotz2-dpotz,"\t",dmax2-dmax,file=fout_diff)


fout_f.close()
fout_py.close()
fout_diff.close()
