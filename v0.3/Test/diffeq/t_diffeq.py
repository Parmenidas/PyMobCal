'''
Test of function diffeq. This version has been vectorized

Here we provide a known potential (see function harmonic)
to both functions. We chose the harmonic 
as it is analytically solvable and it is easy to check the output
for mistakes.

To make it work, we added a common/testing/ to the fortran function
to mimick some local variables

The main change at this point is the replacement of the old integration
algorithm (runge-kutta + predictor correct) with velocity verlet.

As we are comparing different algorithms, we calculate the differenze with the analytic solution.

'''

import diffeq_f 
import numpy as np

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
        self.traj_lost=30000        #DIFF: This was hard-coded into the function. I believe this is the right place
        self.vec_dim=100            #DIF: this is what is used within the FORTRAN function as dimension of arrays 

def diffeq_py(tim,dt,position,momenta,parmList,ion,coord,lj,chargeList,dipole,Dposition,Dmomenta):
    '''
    /*
    * This function integrates the equations of motions using the velocity verlet algorithm.
    * 
    * NOTE: IN/OUT shows the variable that are modified by the function
    * 
    * Arguments
    * ---------
    * time (IN/OUT) : the time of the previous step
    * dt (IN/OUT) : integration time step
    * position (IN/OUT) : vector with the position of the gas particle. *must* be of dimesion = 3
    * momenta (IN/OUT) : vector with the linear  momentum (that is, velocity * mass) of the gas particle. *must* be of dimesion = 3
    * coord : matrix with the position of each atom. Each row is a different atom an contains the coordinates x,y,z
    * ion : information about the molecules (dipole, mass, etc.)
    * lj : structure with info about the Lennard-Jones parameters
    * chargeList : vector with the charges of each atom in the molecules
    * dipole : total dipole of the molecule
    * n_atoms : number of atoms within the molecule
    * Dposition (IN/OUT) : vector with the derivative of the position of the gas particle. *must* be of dimesion = 3
    * Dmomentum (IN/OUT) : vector with the derivative of momentum of the gas particle. *must* be of dimesion = 3
    * 
    * Returns
    * -------
    * dt (IN/OUT) : integration time step
    * time (IN/OUT) : the time of the updated position, velocity and their derivatives
    * pot : total potential energy (LJ + dipole)
    * dpotx,dopty,dpotz : derivative of the potential along each direction x,y,z. NOTE: the force is -dpot (*minus* dpot)
    * dmax : minimum distance between the particle of the gas and the molecule (namely, closest atom of the molecule)
    * Dposition (IN/OUT): vector with the derivative of the position of the gas particle. Dimesion = 3
    * Dmomenta (IN/OUT): vector with the derivative of momentum of the gas particle. Dimesion = 3
    * position (IN/OUT) : vector with the position of the gas particle. Dimension = 3
    * momenta (IN/OUT) : vector with the linear  momentum (that is, velocity * mass) of the gas particle. Dimension = 3
    */    


    '''

    position=position+momenta*ion.invMu*dt+0.5*dt*dt*ion.invMu*Dmomenta

    #half step velocity
    momenta=momenta+0.5*dt*Dmomenta

    #compute Forces
    pot,dpotx,dpoty,dpotz,dmax=harmonic(position[0],position[1],position[2])

    Dmomenta=np.array([-dpotx,-dpoty,-dpotz])
 
    #finish step velocity
    momenta=momenta+0.5*dt*Dmomenta
    Dposition=momenta*ion.invMu

    tim=tim+dt
 
    return position,momenta,Dposition,Dmomenta,dt,tim,pot,dpotx,dpoty,dpotz,dmax

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


def harmonic(x,y,z):
    '''
    Harmonic potential. This potential "replaces" the original dljpot in the fortran code.
    We did *NOT* change any fortran functions, but we used the .pyf file to
    play this trick. In brief, in the .pyf file we "explain" that when running the
    code, we should use the function we provide instead of dljpot. Of course, the substitute function 
    must have the same input/output variables passed in the same order.
    '''
    k=1.
    pot=-0.5*k*(x**2+y**2+z**2)
    dpotx=k*x
    dpoty=k*y
    dpotz=k*z
    dmax=0
    return pot,dpotx,dpoty,dpotz,dmax

# We need to set some variables. Firstly the general variable used by
# *BOTH* FORTRAN and PYTHON code.
mu=1.
time=0.
dt=0.1 
w=np.array([0.,1.,2.,0.,1.,1.])  # initial state x=0 vx=1, y=2 vy=0, z=1 vz=1

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
parmList=constants()
tim=time
dw,pot,dpotx,dpoty,dpotz,dmax=deriv_py(w,ion,coord,lj,chargeList,dipole)
# Here, the FORTRAN variables
diffeq_f.constants.mu=mu
l=0
dw_f=dw
time_f=time
dt_f=dt
w_f=w

fout_f=open("fortran_file.dat","w") # file with the result for fortran
fout_py=open("python_file.dat","w") # file with the result for python
fout_f_diff=open("diff_f_file.dat","w") # file with the difference between fortran code and analytic model
fout_p_diff=open("diff_p_file.dat","w") # file with the difference between python code and analytic model

position=np.array([w[0],w[2],w[4]])
momenta=np.array([w[1],w[3],w[5]])
Dposition=np.array([dw[0],dw[2],dw[4]])
Dmomenta=np.array([dw[1],dw[3],dw[5]])

for steps in range(0,100):

    position,momenta,Dposition,Dmomenta,dt,tim,pot,dpotx,dpoty,dpotz,dmax=diffeq_py(tim,dt,position,momenta,parmList,ion,coord,lj,chargeList,dipole,Dposition,Dmomenta)
    dw=np.array([Dposition[0],Dmomenta[0],Dposition[1],Dmomenta[1],Dposition[2],Dmomenta[2]])
    print(round(tim,3),' '.join(map(str,dw)),file=fout_py)
    l,time_f,dt_f,w_f,dw_f,pot_f,dmax_f=diffeq_f.diffeq(l,time_f,dt_f,w_f,dw_f)
    print(round(time_f,3),' '.join(map(str,dw_f)),file=fout_f)
    dw_ideal=np.array([np.cos(tim),-position[0],-2*np.sin(tim),-position[1],-np.sin(tim)+np.cos(tim),-position[2]])
    print(round(tim,3),' '.join(map(str,dw-dw_ideal)),file=fout_p_diff)
    print(round(tim,3),' '.join(map(str,dw_f-dw_ideal)),file=fout_f_diff)

fout_f.close()
fout_py.close()
fout_f_diff.close()
fout_p_diff.close()
