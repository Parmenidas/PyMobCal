'''
Test of function gsang
note that the fortran version of gsang has been changed adding 
open/close the file and a check for ns>30000 has been replaced with 
a check for ns>1000 to speed up the test. The function diffeq has an
extra common /testing/ for local variables.

**IMPORTANT**: These function create their own output files that
should be checked for differences (out.fortrant.dat and out.pyhton.dat).
This is in *addition* to diff_file.dat
'''


import gsang_f
import numpy as np
import sys

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
                dw,pot,dpotx,dpoty,dpotz,dmax=deriv(w,ion,coord,lj,chargeList,dipole)
                dw=self.Dt*dw
                r=self.a[j]*np.subtract(dw,self.b[j]*self.q)           
                w=w+r
                self.q=self.q+3*r+self.c[j]*dw
            dw,pot,dpotx,dpoty,dpotz,dmax=deriv(w,ion,coord,lj,chargeList,dipole)

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
        dw,pot,dpotx,dpoty,dpotz,dmax=deriv(w,ion,coord,lj,chargeList,dipole)
        
        self.array_f[-1]=self.acst*dw
        for j in range(6):
            for i in range(4):
                self.array_f[i,j]=self.array_f[i+1,j]
                self.array_f[-1,j]=self.array_f[i,j]*self.amcc[i]+self.array_f[-1,j]
        self.array_f[4]=savdw
        w=savw+self.hcvar*np.add(self.array_f[4],self.array_f[-1])
        dw,pot,dpotx,dpoty,dpotz,dmax=deriv(w,ion,coord,lj,chargeList,dipole)
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

def MapCommonsConstants(f_pack,switch,constants,mu,M,dipol,mobility,cfact,inatom,icoord,iic,lj):
    '''
    Set the common blocks /constants/,/printswitch/ and /trajectory/
    from the corresponding python data structures
    '''
    f_pack.printswitch.ip=switch.ip
    f_pack.printswitch.it=switch.it
    f_pack.printswitch.iu1=switch.iu1
    f_pack.printswitch.iu2=switch.iu2
    f_pack.printswitch.iu3=switch.iu3
    f_pack.printswitch.iv=switch.iv
    f_pack.printswitch.im2=switch.im2
    f_pack.printswitch.im4=switch.im4
    f_pack.printswitch.igs=switch.igs
    f_pack.constants.mu=mu
    f_pack.constants.ro=constants.ro
    f_pack.constants.eo=constants.eo
    f_pack.constants.pi=constants.pi
    f_pack.constants.cang=constants.cang
    f_pack.constants.ro2=constants.ro2
    f_pack.constants.dipol=dipol
    f_pack.constants.emax=None             # Only used in the function potent() that however is NOT used in "vanilla" Mobcal
    f_pack.constants.m1=4.0026       
    f_pack.constants.m2=M                 # This is the mass of the ion
    f_pack.constants.xe=constants.xe
    f_pack.constants.xk=constants.xk
    f_pack.constants.xn=constants.xn
    f_pack.constants.xeo=constants.xeo
    f_pack.constants.xmv=constants.xmv
    f_pack.constants.mconst=mobility
    f_pack.constants.correct=cfact
    f_pack.constants.romax=lj.romax
    f_pack.constants.inatom=inatom
    f_pack.constants.icoord=icoord
    f_pack.constants.iic=iic
    f_pack.trajectory.sw1=constants.sw1
    f_pack.trajectory.sw2=constants.sw2
    f_pack.trajectory.dtsf1=constants.dtsf1
    f_pack.trajectory.dtsf2=constants.dtsf2
    f_pack.trajectory.cmin=constants.cmin
    f_pack.trajectory.ifail=constants.ifail
    f_pack.trajectory.ifailc=constants.ifailc
    f_pack.trajectory.inwr=constants.inwr

def MapLJConstants(f_pack,lj):
    '''
    Set the common block /ljparameters/from the corresponding python data structures
    Fortran has a fixed size for arrays = 1000
    '''
    f_pack.ljparameters.eolj=np.zeros(1000)
    f_pack.ljparameters.rolj=np.zeros(1000)
    f_pack.ljparameters.eox4=np.zeros(1000)
    f_pack.ljparameters.ro6lj=np.zeros(1000)
    f_pack.ljparameters.ro12lj=np.zeros(1000)
    f_pack.ljparameters.dro6=np.zeros(1000)
    f_pack.ljparameters.dro12=np.zeros(1000)
        
    for c,i in enumerate(lj.eolj):
        f_pack.ljparameters.eolj[c]=lj.eolj[c]
        f_pack.ljparameters.rolj[c]=lj.rolj[c]
        f_pack.ljparameters.eox4[c]=lj.eox4[c]
        f_pack.ljparameters.ro6lj[c]=lj.ro6lj[c]
        f_pack.ljparameters.ro12lj[c]=lj.ro12lj[c]
        f_pack.ljparameters.dro6[c]=lj.dro6[c]
        f_pack.ljparameters.dro12[c]=lj.dro12[c]

def gsang(coord,theta,phi,gamma,chargeList,dipole,lj,v,b,switch,fout,ion,parmList):
    vy=-v
    vx=0.0
    vz=0.0
    d1=0.0
    istep=0
    vxyz=np.abs(vy)       ##NOTE: not used?
    if switch.it==1:
        print("\n specific trajectory parameters\n\n v ={0: .4E}    b ={1: 0.4E}".format(v,b),file=fout)
#     determine time step
    
    top=(v/95.2381)-0.5
    if v>=1000:
        top=10.0
    if v>=2000:
        top=10.0-((v-2000)*7.5e-3)
    if v>=3000:
        top=2.5
    dt1=top*parmList.dtsf1*1.0e-11/v
    dt2=dt1*parmList.dtsf2
    dt=dt1
    if switch.it==1:
        print(' time steps, dt1 ={0: .4E} dt2 ={1: .4E}'.format(dt1,dt2),file=fout)
    
#     determine trajectory start position
    e0=0.5*ion.mu*v*v
    x=b
    z=0.0
#    ymin=0.0
#    ymax=0.0
    ymax=np.amax(coord[:,1])
    ymin=np.amin(coord[:,1])
    ymax=ymax*1e10
    ymin=ymin*1e10
    iymin=np.trunc(ymin)-1
    iymax=np.trunc(ymax)+1
    id2=iymax
    y=float(id2)*1.0e-10
    pot,dpotx,dpoty,dpotz,dmax=dljpot(x,y,z,coord,lj,chargeList,dipole)

    if np.abs(pot/e0)<=parmList.sw1:
        while True:
            id2=id2-1
            y=float(id2)*1.0e-10
            pot,dpotx,dpoty,dpotz,dmax=dljpot(x,y,z,coord,lj,chargeList,dipole)
            if id2<iymin:
                if switch.it==1:
                    print(' trajectory not started - potential too small',file=fout)
                ang=0.0
                erat=1.0
                return ang,erat,d1,istep
            y=float(id2)*1.0e-10
            pot,dpotx,dpoty,dpotz,dmax=dljpot(x,y,z,coord,lj,chargeList,dipole)
            if np.abs(pot/e0)>=parmList.sw1:
                break
    else:    
        while True:
            id2=id2+10
            y=float(id2)*1.0e-10
            pot,dpotx,dpoty,dpotz,dmax=dljpot(x,y,z,coord,lj,chargeList,dipole)
            if np.abs(pot/e0)<=parmList.sw1:
                break
        
        while True:
            id2=id2-1
            y=float(id2)*1.0e-10
            pot,dpotx,dpoty,dpotz,dmax=dljpot(x,y,z,coord,lj,chargeList,dipole)
            if np.abs(pot/e0)>=parmList.sw1:
                break

    y=float(id2)*1.0e-10
    etot=e0+pot
    if switch.it==1:
        print(' trajectory start position ={0: .4E}'.format(y*1e10),file=fout)
    d1=y
#     initial coordinates and momenta
    w=np.zeros(6)
    w[0]=x
    w[1]=vx*ion.mu
    w[2]=y
    w[3]=vy*ion.mu
    w[4]=z
    w[5]=vz*ion.mu
    tim=0.0
    if switch.it==1:
        x=' '
        print('\n\n trajectory ns, x,  y,  z,  kin e, dt,    tot e',file=fout)
        print(16*x+'vx, vy, vz, pot e, pot/e0\n',file=fout)    
    
    #     initialize the time derivatives of the coordinates and momenta
    dw,pot,dpotx,dpoty,dpotz,dmax=deriv(w,ion,coord,lj,chargeList,dipole)
    ns=0
    nw=0
    ang=0.
    erat=0.0
    #l=0
    #print(dw)
    Integrator=IntegratorState(dt)
    while True:
        while True:
            while True:
                while True:
                    w,dt,tim,dw,pot,dpotx,dpoty,dpotz,dmax=diffeq(tim,dt,w,parmList,ion,coord,lj,chargeList,dipole,Integrator,dw)
                    nw=nw+1
                    if nw==parmList.inwr:
                        break
                ns=ns+nw
                nw=0
            #     print out the trajectory coordinates and velocities
                if switch.it==1:
                    e=w[1]**2/(2.0*ion.mu)+w[3]**2/(2.0*ion.mu)+w[5]**2/(2.0*ion.mu)
                    x=' '
                    #print(" {:5d}".format(ns)," {: 11.4E}".format(w[0])," {: .4e}".format(w[2])," {: .4e}".format(w[4])," {: .4e}".format(e)," {: .4e}".format(dt)," {: .4e}".format(pot+e),file=fout)
                    print(" {:5d} {: .4E} {: .4E} {: .4E} {: .4E} {: .4E} {: .4E}".format(ns,w[0],w[2],w[4],e,dt,pot+e),file=fout)
                    print(7*x+"{: .4E} {: .4E} {: .4E} {: .4E} {: .4E}".format(dw[0],dw[2],dw[4],pot,np.abs(pot/e0)),file=fout)

            #     check if trajectory has become "lost" (too many steps)

                if ns>parmList.traj_lost:
                    print(' trajectory lost: b ={0: .4E} v={1: .4E}'.format(b,v),file=fout)
                    ang=parmList.pi/2.0
                    e=0.5*ion.mu*(dw[0]**2+dw[2]**2+dw[4]**2)
                    erat=(e+pot)/etot
                    istep=ns
                    return ang,erat,d1,istep

            #     check if the trajectory is finished
                if dmax>=lj.romax:
                    break
            if np.abs(pot/e0)>parmList.sw2 and dt==dt1:
                dt=dt2
                Integrator.SetStatus(0)
            if np.abs(pot/e0)<parmList.sw2 and dt==dt2:
                dt=dt1
                Integrator.SetStatus(0)
            if np.abs(pot/e0)<=parmList.sw1:
                break
        if ns>=50:
            break
    istep=ns
    #     determine scattering angle 
    if(dw[0]>0.0):        
        num=dw[2]*(-v)
        den=v*np.sqrt(dw[0]**2+dw[2]**2+dw[4]**2)
        ang=np.arccos(num/den)

    if(dw[0]<0.0):
        num=dw[2]*(-v)
        den=v*np.sqrt(dw[0]**2+dw[2]**2+dw[4]**2)
        ang=-1*np.arccos(num/den)
    
    #     check for energy conservation
    e=0.5*ion.mu*(dw[0]**2+dw[2]**2+dw[4]**2)
    erat=(e+pot)/etot
    if erat<1.01 and erat>0.99:
        return ang,erat,d1,istep  #ang
    print('\n energy not conserved: e ratio ={0: .4E} v ={1: .4E} b ={2: .4E}'.format(erat,v,b),file=fout)
    print(' gst2 ={0: .4E} theta ={1: .4E} phi ={2: .4E} gamma ={3: .4E}'.format(0.5*ion.mu*v*v/parmList.eo,theta*parmList.cang,phi*parmList.cang,gamma*parmList.cang),file=fout)
    if switch.ip==1:
        print('')
    parmList.ifailc=parmList.ifailc+1   ## make ifailc external to list
    if parmList.ifailc==parmList.ifail:
        fout.close()
        print('Energy not conserved')
        sys.exit(0)
    return ang,erat,d1,istep

def diffeq(tim,dt,w,parmList,ion,coord,lj,chargeList,dipole,Integrator,dw):
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

def deriv(w,ion,coord,lj,chargeList,dipole):
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
    pot,dpotx,dpoty,dpotz,dmax=dljpot(x,y,z,coord,lj,chargeList,dipole)
    dw[1]=-dpotx
    dw[3]=-dpoty
    dw[5]=-dpotz

    return dw,pot,dpotx,dpoty,dpotz,dmax


def dljpot(x,y,z,vector,lj,charge,dipol):
    rx=0.0
    ry=0.0
    rz=0.0
    e00=0.0
    de00x=0.0
    de00y=0.0
    de00z=0.0
    sum1=0.0
    sum2=0.0
    sum3=0.0
    sum4=0.0
    sum5=0.0
    sum6=0.0
    dmax=2*lj.romax

    for ivec,ie,ir12,ir6,dr6,dr12,q in zip(vector,lj.eox4,lj.ro12lj,lj.ro6lj,lj.dro6,lj.dro12,charge):
        xx=x-ivec[0]
        xx2=xx*xx
        yy=y-ivec[1]
        yy2=yy*yy
        zz=z-ivec[2]
        zz2=zz*zz
        rxyz2=xx2+yy2+zz2
        rxyz=np.sqrt(rxyz2)
        if rxyz<dmax:
            dmax=rxyz
        rxyz3=rxyz2*rxyz
        rxyz5=rxyz3*rxyz2
        rxyz6=rxyz5*rxyz
        rxyz8=rxyz5*rxyz3
        rxyz12=rxyz6*rxyz6
        rxyz14=rxyz12*rxyz2
    #     LJ potential 
        e00=e00+(ie*((ir12/rxyz12)-(ir6/rxyz6)))
    #     LJ derivative
        de00=ie*((dr6/rxyz8)-(dr12/rxyz14))
        de00x=de00x+(de00*xx)
        de00y=de00y+(de00*yy)
        de00z=de00z+(de00*zz)
    #     ion-induced dipole potential
        if(q==0):
            continue
        rxyz3i=q/rxyz3
        rxyz5i=-3*q/rxyz5
        rx=rx+(xx*rxyz3i)
        ry=ry+(yy*rxyz3i)
        rz=rz+(zz*rxyz3i)
    #     ion-induced dipole derivative
        sum1=sum1+(rxyz3i+(xx2*rxyz5i))
        sum2=sum2+(xx*yy*rxyz5i)
        sum3=sum3+(xx*zz*rxyz5i)
        sum4=sum4+(rxyz3i+(yy2*rxyz5i))
        sum5=sum5+(yy*zz*rxyz5i)
        sum6=sum6+(rxyz3i+(zz2*rxyz5i))

    pot=e00-(dipol*((rx*rx)+(ry*ry)+(rz*rz)))
    dpotx=de00x-(dipol*((2.0*rx*sum1)+(2.*ry*sum2)+(2.*rz*sum3)))
    dpoty=de00y-(dipol*((2.0*rx*sum2)+(2.*ry*sum4)+(2.*rz*sum5)))
    dpotz=de00z-(dipol*((2.0*rx*sum3)+(2.*ry*sum5)+(2.*rz*sum6)))

    return pot,dpotx,dpoty,dpotz,dmax

#####################################################################
## Test of function gsang
## note that the fortran version of gsang has been changed adding 
## open/close the file and a check for ns>30000 has been replaced with 
## a check for ns>1000 to speed up the test. The function diffeq has an
## extra common /testing/ for local variables.
## 
## **IMPORTANT**: These function create their own output files that
## should be checked for differences (out.fortrant.dat and out.pyhton.dat).
## This is in *addition* to diff_file.dat
#####################################################################



##########################################################################
## Set variables. This is the values used by *BOTH* FORTRAN and PYTHON code
##########################################################################
#Define atom postions, charges, etc... [equivalent to fcoord()]
#atoms positioned along x. LJ: sigma=1, epsilon=1
#Firstly, the position of the gas
x=0.
y=0.
z=0.
#Secondly, the position of the ion
fx=0.
fy=0.
fz=0.

romax=1.
pcharge=0.
eolj=103-10
rolj=10e-10
mu=700.
theta,phi,gamma=1.,2.,3. #no particular meaning, random

parmList=constants()
dipole=parmList.xe**2*0.204956e-30/(8*parmList.pi*parmList.xeo)
switch=printSwitch()
switch.it=1
lj=LennardJones(np.array([eolj]),np.array([rolj]))
chargeList=np.array([0])
###############################################################
# Set Commons (FORTRAN). Not all variables in commons are actually used. 
# We only set the ones used.
###############################################################
MapCommonsConstants(gsang_f,switch,parmList,mu,None,dipole,None,None,1,1,1,lj)
MapLJConstants(gsang_f,lj)
gsang_f.charge.pcharge[0]=chargeList[0]
gsang_f.coordinates.fx[0]=fx
gsang_f.coordinates.fy[0]=fy
gsang_f.coordinates.fz[0]=fz
gsang_f.angles.theta=theta
gsang_f.angles.phi=phi
gsang_f.angles.gamma=gamma
#######################################
## Python code
#######################################
coord=np.array([[fx,fy,fz]])
ion=molecule()
ion.mu=mu

fout_py=open("python_file.dat","w") # file with the result for python
fout_f=open("fortran_file.dat","w") # file with the result for fortran
fout_diff=open("diff_file.dat","w") # file with the difference of the two files above
fout_t_py=open("out.python.dat","w") # file with the result for python

#The pairs (v,b) are chosen to trigger different parts of the code. Missing an example for "trajectory too long" error
headers=((3.0,parmList.ro,'fail - ns>50'),(0.01,0.,'fail - ns>50, erat 1'),(3.0,0.5,'fail - potential too small'),(3000.,0.8*parmList.ro,'success'),(1000.,0.5*parmList.ro,'energy not conserved'))

for group in headers:
    v=group[0]
    b=group[1]
    note=group[2]
    #print (v,b,note)    
    erat_f,ang_f,d1_f,istep_f=gsang_f.gsang(v,b)
    print(note+'\n',erat_f,ang_f,d1_f,istep_f,file=fout_f)
    #print (v,b,note)    
    ang,erat,d1,istep=gsang(coord,theta,phi,gamma,chargeList,dipole,lj,v,b,switch,fout_t_py,ion,parmList)
    #gsang(coord,theta,phi,gamma,chargeList,dipole,lj,v,b,switch,fout_t_py,ion,parmList)
    print(note+'\n',erat,ang,d1,istep,file=fout_py)
    print(note+'\n',erat_f-erat,ang_f-ang,d1_f-d1,istep_f-istep,file=fout_diff)

fout_py.close()
fout_f.close()
fout_diff.close()
fout_t_py.close()
