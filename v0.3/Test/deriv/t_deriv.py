'''
Test of function deriv. This version has been vectorized

Here we provide a known potential (see function harmonic)
to both functions and check deriv. We chose the harmonic 
as it is analytically solvable and it is easy to check the output
for mistakes

We add the variable invMu in the class molecules. As we use 
1/mu a lot, this reduces the operations bit.

'''

import deriv_f
import numpy as np

class molecule:
    '''
    /* 
    * Structure containing information about the system under investigation
    * 
    * Attributes
    * ----------
    * TotalMass : total mass of the molecule (excluding the gas particle)
    * dipole : dipole of the molecule (excluding the gas particle)
    * totalCharge : total charge of the molecule (excluding the gas particle)
    * totalAbsoluteCharge : sum of the absolute value of the charges (excluding the gas particle)
    * mGas : mass of the gas particle
    * mu : reduced mass of the system (molecule + gas particle)
    * invMu: 1/mu. divisions is an expensive operation, so we compute it once here and reuse
    * mobility : mobility coefficient of the molecule (excluding the gas particle)
    * temperature : temperature of the system in Kelvin
    * 
    */    

    '''

    def __init__(self):
        self.TotalMass=0.
        self.dipole=0.
        self.totalCharge=0.
        self.totalAbsoluteCharge=0.
        self.mGas=0.
        self.mu=0.  # reduced mass
        self.invMu=0. # as division is very expensive, we calculate this once 
        self.mobility=0.
        self.temperature=298.  #kelvin

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
        self.invMu=np.reciprocal(self.mu)

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

def deriv_py(position,momenta,ion,coord,lj,chargeList,dipole):
    '''
    /*
    * This function computes the derivative with respect time of position and momenta
    * 
    * Arguments
    * ---------
    * position : vector with the position of the gas particle. *must* be of dimesion = 3
    * moment1 : vector with the linear  momentum (that is, velocity * mass) of the gas particle. *must* be of dimesion = 3
    * ion : proeprties of the molecules under examination
    * coord : coordinates of the atoms
    * lj : structure with info about the Lennard-Jones parameters
    * chargeList : vector with the charges of each atom in the molecules
    * dipole : total dipole of the molecule
    * 
    * Returns
    * -------
    * pot : total potential energy (LJ + dipole)
    * dpotx,dpoty,dpotz : derivative of the potential along each direction x,y,z. NOTE: the force is -dpot (*minus* dpot)
    * dmax : minimum distance between the particle of the gas and the molecule (namely, closest atom of the molecule)
    * Dposition : vector with the derivative of the position of the gas particle. Dimension = 3
    * Dmomentum : vector with the derivative of momentum of the gas particle. Dimension = 3
    */
    '''
    # Update derivative of the positions with respect to time
    Dposition=momenta*ion.invMu

    # Update the potential and its derivatives
    pot,dpotx,dpoty,dpotz,dmax=harmonic(position[0],position[1],position[2])

    # Update derivative of the momenta with respect to time
    Dmomenta=np.array([-dpotx,-dpoty,-dpotz])

    return Dposition,Dmomenta,pot,dpotx,dpoty,dpotz,dmax

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
ion.invMu=1./mu

# Here, the FORTRAN variables
deriv_f.constants.mu=mu

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
    
    dw,pot,dpotx,dpoty,dpotz,dmax=deriv_f.deriv(w,harmonic)
    print(' '.join(map(str,dw)),"\t",pot,"\t",dpotx,"\t",dpoty,"\t",dpotz,"\t",dmax,file=fout_f)
    Dposition,Dmomenta,pot2,dpotx2,dpoty2,dpotz2,dmax2=deriv_py(np.array([x,y,z]),np.array([vx,vy,vz]),ion,coord,lj,chargeList,dipole)
    dw2=[Dposition[0],Dmomenta[0],Dposition[1],Dmomenta[1],Dposition[2],Dmomenta[2]]
    print(' '.join(map(str,dw2)),"\t",pot2,"\t",dpotx2,"\t",dpoty2,"\t",dpotz2,"\t",dmax2,file=fout_py)
    print(' '.join(map(str,dw2-dw)),"\t",pot2-pot,"\t",dpotx2-dpotx,"\t",dpoty2-dpoty,"\t",dpotz2-dpotz,"\t",dmax2-dmax,file=fout_diff)


fout_f.close()
fout_py.close()
fout_diff.close()
