'''
Test of function diffeq

Here we provide a known potential (see function harmonic)
to both functions. We chose the harmonic 
as it is analytically solvable and it is easy to check the output
for mistakes.
To make it work, we added a common/testing/ to the fortran function
to mimick some local variables
'''


import diffeq 
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

class IntegratorState:
    '''
    Retains information about the integration method
    '''
    def __init__(self,dt):
        self.var=2.97013888888
        self.cvar=0.990972222222
        self.acst=0.332866152768
        self.a=np.array([0.50,0.292893218814,1.70710678118,0.1666666666667])
        self.b=np.array([2.0,1.0,1.0,2.0])
        self.c=np.array([-0.5,-0.292893218814,-1.70710678118,-0.5])
        self.ampc=np.array([-0.111059153612,0.672667757774,-1.70633621697,2.33387888707,-1.8524668225])
        self.amcc=np.array([0.0189208128941,-0.121233356692,0.337771548703,-0.55921513665])

        self.q=np.zeros(6)
        self.hvar=dt*self.var
        self.hcvar=dt*self.cvar
        self.Dt=0.5*dt
        self.array_f=np.zeros([6,6])
        self.l=0
        self.k=0
    
    def Reset(self,dt):
        self.q=np.zeros(6)
        self.hvar=dt*self.var
        self.hcvar=dt*self.cvar
        self.Dt=0.5*dt

    def GetStatus(self):
        return self.l

    def SetStatus(self,value):
        self.l=value

    def ComputeRKG(self,w,ion,coord,lj,chargeList,dipole,tim):
        self.k=0
        while True:
            for j in range(0,4):
                if ((-1)**(j+1))>0:
                    tim=tim+0.5*self.Dt
                dw,pot,dpotx,dpoty,dpotz,dmax=deriv_py(w,ion,coord,lj,chargeList,dipole)
                dw=self.Dt*dw
                r=self.a[j]*np.subtract(dw,self.b[j]*self.q)           
                w=w+r
                self.q=self.q+3*r+self.c[j]*dw
            dw,pot,dpotx,dpoty,dpotz,dmax=deriv_py(w,ion,coord,lj,chargeList,dipole)

            if self.k>0:
                break
            else:
                self.k=1
        
        if (self.l-6)>=0:
            self.l=-1
            self.Dt=2.0*self.Dt
            return w,self.Dt,tim,dw,pot,dpotx,dpoty,dpotz,dmax
        else:           
            self.array_f[self.l-1]=dw
            return w,self.Dt,tim,dw,pot,dpotx,dpoty,dpotz,dmax

    def ComputeAMpc(self,w,dw,ion,coord,lj,chargeList,dipole,tim):
        savw=w
        savdw=dw
        self.array_f[5]=savdw
        for j in range(6):
            for i in range(5):
                self.array_f[-1,j]=self.array_f[-1,j]+self.ampc[i]*self.array_f[i,j]
        w=self.array_f[-1]*self.hvar+w
        tim=tim+self.Dt
        dw,pot,dpotx,dpoty,dpotz,dmax=deriv_py(w,ion,coord,lj,chargeList,dipole)
        
        self.array_f[-1]=self.acst*dw
        for j in range(6):
            for i in range(4):
                self.array_f[i,j]=self.array_f[i+1,j]
                self.array_f[-1,j]=self.array_f[i,j]*self.amcc[i]+self.array_f[-1,j]
        self.array_f[4]=savdw
        w=savw+self.hcvar*np.add(self.array_f[4],self.array_f[-1])
        dw,pot,dpotx,dpoty,dpotz,dmax=deriv_py(w,ion,coord,lj,chargeList,dipole)
        return w,self.Dt,tim,dw,pot,dpotx,dpoty,dpotz,dmax


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

def diffeq_py(tim,dt,w,parmList,ion,coord,lj,chargeList,dipole,Integrator,dw):
    '''
c     Integration subroutine - uses 5th order runge-kutta-gill to 
c     initiate and 5th order adams-moulton predictor-corrector to 
c     propagate. Parameter l is initially set to zero and then 
c     incremented to tell the subroutine when to switch between 
c     integration methods. DIFFEQ calls subroutine DERIV to define 
c     the equations of motion to be integrated.
    '''

    l=Integrator.GetStatus()

    if l==0:
        Integrator.Reset(dt)
    
    if l>=0:
        l=l+1
        Integrator.SetStatus(l)
#     This is the runge-kutta-gill part...the steps are broken up into
#     half steps to improve accuracy.
        w,dt,tim,dw,pot,dpotx,dpoty,dpotz,dmax=Integrator.ComputeRKG(w,ion,coord,lj,chargeList,dipole,tim)

#     This is the adams-moulton predictor-corrector part.
    if l<0:
        w,dt,tim,dw,pot,dpotx,dpoty,dpotz,dmax=Integrator.ComputeAMpc(w,dw,ion,coord,lj,chargeList,dipole,tim)

    return w,dt,tim,dw,pot,dpotx,dpoty,dpotz,dmax

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
time=0
dt=0.1 
w=np.array([0,1,2,0,1,1])  # initial state x=0 vx=1, y=2 vy=0, z=1 vz=1

## Some variables for the python function.
## All zeros (except for mu) because they are needed for dljpot, 
## but we are not using it in this example
coord=np.array([[0,0,0]])
chargeList=np.array([0])
dipole=0
lj=LennardJones(np.array([0]),np.array([0]))
ion=molecule()
ion.mu=mu
Integrator=IntegratorState(dt)
parmList=constants()
tim=time
dw,pot,dpotx,dpoty,dpotz,dmax=deriv_py(w,ion,coord,lj,chargeList,dipole)
# Here, the FORTRAN variables
diffeq.constants.mu=mu
l=0
dw_f=dw
time_f=time
dt_f=dt
w_f=w

fout_f=open("fortran_file.dat","w") # file with the result for fortran
fout_py=open("python_file.dat","w") # file with the result for python
fout_diff=open("diff_file.dat","w") # file with the difference between the two files above

for steps in range(0,100):
    w,dt,tim,dw,pot,dpotx,dpoty,dpotz,dmax=diffeq_py(tim,dt,w,parmList,ion,coord,lj,chargeList,dipole,Integrator,dw)
    print(round(tim,3),' '.join(map(str,dw)),file=fout_py)
    l,time_f,dt_f,w_f,dw_f,pot_f,dmax_f=diffeq.diffeq(l,time_f,dt_f,w_f,dw_f)
    print(round(time_f,3),' '.join(map(str,dw_f)),file=fout_f)
    print(round(tim,3),' '.join(map(str,dw-dw_f)),file=fout_diff)

fout_f.close()
fout_py.close()
fout_diff.close()
