import math
import numpy as np
import scipy
from scipy.spatial.transform import Rotation as R
import sys

# This software is a version of a prior software called mobcal (see details below)
# Mobcal is the standard in ion mobility calculation and its original code
# was developed in fortran. Here I have "translated" the same code in python.
# The code is as close as possible to the original version in FORTRAN, with only
# a few changes. This decision was taken in order to allow everybody who has previous experience
# with the FORTRAN code to easily understand this new version and develop it further or fix bugs.
# A few changes were made for performance,readiility and maitencance (see sections called DIFF).
# For instance, the software uses the numpy random generator instead of its own. 
# Every time that it makes sense I have used an existing function found in standard python packages.

# a  
# c
# c     PROGRAM MOBCAL
# c
# c     Program to Calculate Mobilities
# c
# c     See: M. F. Mesleh, J. M. Hunter, A. A. Shvartsburg, G. C. Schatz,
# c     and M. F. Jarrold, "Structural Information from Ion Mobility
# c     Measurements: Effects of the Long Range Potential" J. Phys. Chem.
# c     1996, 100, 16082-16086; A. A. Shvartsburg and M. F. Jarrold,
# c     "An Exact Hard Spheres Scattering Model for the Mobilities of
# c     Polyatomic Ions" Chem. Phys. Lett. 1996, 261, 86.
# c


#DIFF: Python does not have the concept of COMMON used in fortran. 
# I replaced it with classes, except when the common is one array/value  

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
            return w,dw,self.Dt,tim,dw,pot,dpotx,dpoty,dpotz,dmax
        else:           
            self.array_f[self.l-1]=dw
            return w,dw,self.Dt,tim,dw,pot,dpotx,dpoty,dpotz,dmax

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


    
class AtomData:
    '''
    Contains information about each atom within the molecule
    '''
    def __init__(self,mass=[],charge=[],coord=[]):
        self.Mass=mass
        self.Charges=charge
        self.Coords=coord
    
    def UpdatePosition(self,coord=[]):
        self.Coords=coord

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

class infile:
    '''
    Information found in the input file. This is the class used within the code itself.
    Ideally, this software can read different types of input, but all of those will
    populate this class. 
    In this way, we make the "reading file" phase independent from its refernce within the code itself.
    '''
    def __init__(self,inList):
        self.atoms=inList["atoms"]
        self.icoord=inList["icoord"]
        self.units=inList["units"]
        self.charges=inList["charges"]
        self.cfact=inList["cfact"]
        self.label=inList["label"]
        
class atom:
    '''
    It contains information about the atoms in the periodic table (not within the molecule!!)
    Not all atoms of the periodic table are supported, as some properties have not been estimated yet
    '''
    def __init__(self,name,mass,e,s,rhs,id):
        self.name=name
        self.mass=mass
        self.epsilonLJ=e
        self.sigmaLJ=s
        self.rhs=rhs
        self.id=id
    
    def Unpack(self):
        return self.mass,self.epsilonLJ,self.sigmaLJ,self.rhs

class constants:
    '''
    Constants used within the code
    '''
    def __init__(self):
        self.pi=math.pi
        self.cang= 180/math.pi
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
        self.vec_dim=100            #DIIF: this is what is used within the FORTRAN function as dimension of arrays 


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
    

class LennardJones:
    '''
    It contains some variable that are needed for the calculation.
    Instead of being evaluated over and over again, we compute them once and reuse
    '''
    def __init__(self,eolj,rolj):
        self.eolj=eolj
        self.rolj=rolj
        self.eox4=4*eolj
        self.ro2lj=np.multiply(rolj,rolj)
        self.ro6lj=np.power(self.ro2lj,3)
        self.ro12lj=np.multiply(self.ro6lj,self.ro6lj)
        self.dro6=6*self.ro6lj
        self.dro12=12*self.ro12lj
        self.romax=np.amax(rolj)

class HardSphere:
    '''
    Parameters used for the hard sphere model.
    '''
    def __init__(self,rhs):
        self.rhs=rhs
        self.rhs2=np.multiply(rhs,rhs)

def GetAtom(AtomList,ConstList,id):
    '''
    For the original file format of mobcal, the atoms is identifyed by its mass
    This function comapares the mass found in teh file with the database and set all
    the internal variable for that atom.
    If the atom is NOT in the database it is defaulted to a carbon atom and a warning is issued.
    Whether using a carbon atom as default is a good idea or not, I will let you decide
    '''
    flag=0
    for element in AtomList:
        if element.id==id:
            flag=1
            return element.Unpack()

    if flag==0:
        print("****** WARNING   WARNING   WARNING ******")
        print("One of the atoms listed is unknown, assuming carbon")
        print("****** WARNING   WARNING   WARNING ******")
        default=atom('Carbon',12.01,1.34e-3*ConstList.xe,3.043e-10,2.7e-10,12)
        return default.Unpack()

def rotate(vector,theta,phi,gamma,ConstList,switch,fout):
    '''
    Rotates the vector(s) provided
    '''
    if switch.iu2==1 or switch.iu3==1:  
        print('coordinates rotated by ROTATE\n theta=',theta*ConstList.cang,' phi=',phi*ConstList.cang,' gamma=',gamma*ConstList.cang,file=fout)

    rz=R.from_rotvec(theta*np.array([0,0,1]))
    rx=R.from_rotvec(-phi*np.array([1,0,0]))
    rzb=R.from_rotvec(gamma*np.array([0,0,1]))
    newVector=rzb.apply(rx.apply(rz.apply(vector)))

    if switch.iu2==1:  
        print("initial coordinates\t\t\t\t\t\t\t new coordinates",file=fout)
        for v,rotv in zip(vector,newVector):
            print(v[0],v[1],v[2],rotv[0],rotv[1],rotv[2],file=fout)    

    return newVector

def fcoord(fileIN,fileOUT,parmList,AtomList,switch):
    '''
    Read the input file. The format is the original mobcal format
    '''
    inList={
    "atoms": 0,
    "icoord": 0,
    "units": " ",
    "charges": " ",
    "cfact" :0,
    "label":" "
    }
    
    ion=molecule()

    #check the files are actually there
    try:
        fin=open(fileIN,"r")
    except IOError as e:
        print("Couldn't open file (%s)." % e)
        #return 0

    try:
        fout=open(fileOUT,"w")
    except IOError as e:
        print("Couldn't open file (%s)." % e)
        #return 0

    print("input file name = ",fileIN,file=fout)
    tmp=fin.readline()
    inList["label"]=tmp
    print("input file label =",tmp,file=fout,end="")
    tmp=fin.readline()
    inList["icoord"]=int(tmp)
    print("number of coordinate sets=", tmp,file=fout,end="")
    tmp=fin.readline()
    inList["atoms"]=int(tmp)
    print("number of atoms=",tmp,file=fout,end="")
    
    tmp=fin.readline()
    x=tmp.split()
    inList["units"]=x[0]
    if x[0]=='au':
        fout.write("coordinates in atomic units\n")
    elif x[0]=='ang':
        fout.write("coordinates in angstroms\n")
    else:
        fout.write("Coordinates selected not supported\n")
        #return 0

    tmp=fin.readline()
    x=tmp.split()
    inList["charges"]=x[0]
    if x[0]=='equal':
        fout.write("using a uniform charge distribution\n")
    elif x[0]=='calc':
        fout.write("using a calculated (non-uniform) charge distribution\n")
    elif x[0]=='none':
        fout.write("using no charge - only LJ interactions\n")
    else:
        fout.write("Charges model selected not supported\n")
        fin.close()
        fout.close()
        sys.exit(0)
        #return 0

    # DIFF: the following line is swicthed with the previous line  inthe original code
    tmp=fin.readline()
    inList["cfact"]=float(tmp)
    print("correction factor for coordinates = ", tmp,file=fout,end="")

    #Seeting up some array before reading coords. Try to be consist with original code as well as optimize a bit
    fx=np.zeros(inList["atoms"])
    fy=np.zeros(inList["atoms"])
    fz=np.zeros(inList["atoms"])
    imass=np.zeros(inList["atoms"],dtype=int)
    if inList["charges"]=='equal':
        pcharges=np.full(inList["atoms"],1/inList["atoms"])
    else:
        pcharges=np.zeros(inList["atoms"])
    xmass=np.zeros(inList["atoms"])
    eolj=np.zeros(inList["atoms"])
    rolj=np.zeros(inList["atoms"])
    rhs=np.zeros(inList["atoms"])

    #Reading atoms coordinates starts here:
    for i in range(inList["atoms"]):
        tmp=fin.readline()
        x=tmp.split()
        fx[i]=float(x[0])
        fy[i]=float(x[1])
        fz[i]=float(x[2])
        imass[i]=int(np.rint(float(x[3])))
        if inList["charges"]=='calc':
            pcharges[i]=float(x[4])

    if inList["units"]=='au':
        fx=fx*0.52917706
        fy=fy*0.52917706
        fz=fz*0.52917706

    tcharge=np.sum(pcharges) #total charge
    acharge=np.sum(np.abs(pcharges)) #total absolute charge
    #DIFF: Original code only prints the charges is model==calc (see next 2 lines). I prefer always: helps catching bugs
    # if inList["charges"]=='calc':
    #     print("total charge = ",tcharge,"\ntotal absolute charge =",acharge,file=fout)
    print("total charge = ",tcharge,"\ntotal absolute charge =",acharge,file=fout)
    ion.SetCharges(tcharge,acharge)

    for i in range(inList["atoms"]):
        xmass[i],eolj[i],rolj[i],rhs[i]=GetAtom(AtomList,parmList,imass[i])
    tMass=np.sum(xmass)
    print("mass of ion = ",tMass,file=fout)
    ion.SetMass(tMass)
    
    lj=LennardJones(eolj,rolj)
    hs=HardSphere(rhs)

    if switch.iu1==1:
        print('initial coordinates \t\t\t mass \t  charge \t\t\t LJ parameters',file=fout)
    
    #The center of mass is used as origin of the reference frame
    fxo=np.sum(np.multiply(fx,xmass))/tMass
    fyo=np.sum(np.multiply(fy,xmass))/tMass
    fzo=np.sum(np.multiply(fz,xmass))/tMass
    print('center of mass coordinates = ',fxo,' ',fyo,' ',fzo,file=fout)

    Oxyz=np.zeros([inList["atoms"],3])   # original coordinates vector
    for i in range(inList["atoms"]):
        Oxyz[i][0]=(fx[i]-fxo)*1.e-10*inList["cfact"]
        Oxyz[i][1]=(fy[i]-fyo)*1.e-10*inList["cfact"]
        Oxyz[i][2]=(fz[i]-fzo)*1.e-10*inList["cfact"]

    if switch.iu1==1:
        for i in range(inList["atoms"]):
            print(fx[i],fy[i],fz[i],imass[i],pcharges[i],eolj[i]/parmList.xe,rolj[i]*1e+10,file=fout)
        print("\n")
    
    if inList["icoord"]==1:
        fin.close()

#    determine structural asymmetry parameter
    theta=0
    asymp=0

    for igamma in range(0,360,2):
        for iphi in range(0,180,2):
            xyzsum=0
            yzsum=0
            gamma=float(igamma/parmList.cang)
            phi=float(iphi/parmList.cang)
            rotVector=rotate(Oxyz,theta,phi,gamma,parmList,switch,fout)
            for ivec in rotVector:
                xyzsum=xyzsum+np.linalg.norm(ivec)
                yzsum=yzsum+np.linalg.norm(ivec[1:])
            hold=((parmList.pi/4)*xyzsum)/yzsum
            if(hold>asymp):
                asymp=hold

    InputData=infile(inList)
    atoms=AtomData(xmass,pcharges,Oxyz)
    return atoms,ion,InputData,asymp,lj,hs,fin,fout

def ncoord(fin,fout,parmList,AtomList,switch,ion,inputList):
    '''
    Read coordinates after the first set (if any). The format is the original mobcal format
    '''
    
    fx=np.zeros(inputList.atoms)
    fy=np.zeros(inputList.atoms)
    fz=np.zeros(inputList.atoms)
    imass=np.zeros(inputList.atoms)
    xmass=np.zeros(inputList.atoms)

    tmp=fin.readline()   #unused line

    #Reading atoms coordinates starts here:
    for i in range(inputList.atoms):
        tmp=fin.readline()
        x=tmp.split()
        fx[i]=float(x[0])
        fy[i]=float(x[1])
        fz[i]=float(x[2])
        imass[i]=int(np.rint(float(x[3])))

    if inputList.units=='au':
        fx=fx*0.52917706
        fy=fy*0.52917706
        fz=fz*0.52917706

    for i in range(inputList.atoms):
        xmass[i],eolj,rolj,rhs=GetAtom(AtomList,parmList,imass[i])
    
    tMass=np.sum(xmass)
    if tMass!=ion.TotalMass:
        print('masses do not add up')
        fin.close()
        fout.close()
        sys.exit(0)    
    
    #The center of mass is used as origin of the reference frame
    fxo=np.sum(np.multiply(fx,xmass))/tMass
    fyo=np.sum(np.multiply(fy,xmass))/tMass
    fzo=np.sum(np.multiply(fz,xmass))/tMass

    Oxyz=np.zeros([inputList.atoms,3])   # original coordinates vector
    for i in range(inputList.atoms):
        Oxyz[i][0]=(fx[i]-fxo)*1.e-10*inputList.cfact
        Oxyz[i][1]=(fy[i]-fyo)*1.e-10*inputList.cfact
        Oxyz[i][2]=(fz[i]-fzo)*1.e-10*inputList.cfact

#    determine structural asymmetry parameter
    theta=0
    asymp=0

    for igamma in range(0,360,2):
        for iphi in range(0,180,2):
            xyzsum=0
            yzsum=0
            gamma=float(igamma/parmList.cang)
            phi=float(iphi/parmList.cang)
            rotVector=rotate(Oxyz,theta,phi,gamma,parmList,switch,fout)
            for ivec in rotVector:
                xyzsum=xyzsum+np.linalg.norm(ivec)
                yzsum=yzsum+np.linalg.norm(ivec[1:])
            hold=((parmList.pi/4)*xyzsum)/yzsum
            if(hold>asymp):
                asymp=hold

    return Oxyz,asymp,fin,fout


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
        w,dw,dt,tim,dw,pot,dpotx,dpoty,dpotz,dmax=Integrator.ComputeRKG(w,ion,coord,lj,chargeList,dipole,tim)

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

def gsang(coord,theta,phi,gamma,chargeList,dipole,lj,v,b,switch,fout,ion,parmList):
    vy=-v
    vx=0.0
    vz=0.0
    vxyz=np.abs(vy)       ##NOTE: not used?
    if switch.it==1:
        print("\nspecific trajectory parameters",file=fout)
        print('v =',v,'b =',b,file=fout)
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
        print('time steps, dt1 =',dt1,' dt2 =',dt2,file=fout)
    
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
                    print('trajectory not started - potential too small',file=fout)
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
        print('trajectory start position =',y*1e10,file=fout)
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
        print('\n\ntrajectory ns, x,  y,  z,  kin e, dt,    tot e',file=fout)
        print('vx, vy, vz, pot e, pot/e0',file=fout)    
    
    #     initialize the time derivatives of the coordinates and momenta
    dw,pot,dpotx,dpoty,dpotz,dmax=deriv(w,ion,coord,lj,chargeList,dipole)
    ns=0
    nw=0
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
                    print(ns,w[0],w[2],w[4],e,dt,pot+e,file=fout)
                    print(dw[0],dw[2],dw[4],pot,np.abs(pot/e0),file=fout)

            #     check if trajectory has become "lost" (too many steps)

                if ns>parmList.traj_lost:
                    print('trajectory lost: b =',b, 'v=',v,file=fout)
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
        return ang,erat,d1,istep
    print('\nenergy not conserved: e ratio =',erat,'v = ',v,'b = ',b,file=fout)
    print('gst2 =',0.5*ion.mu*v*v/parmList.eo,' theta =',theta*parmList.cang,' phi =', phi*parmList.cang,' gamma =',gamma*parmList.cang,file=fout)
    if switch.ip==1:
        print('')
    parmList.ifailc=parmList.ifailc+1   ## make ifailc external to list
    if parmList.ifailc==parmList.ifail:
        fout.close()
        print('Energy not conserved')
        sys.exit(0)

def rantate(vector,parmList,switch,fout):
    '''
c     Rotates the cluster/molecule to a random orientation.  
    '''
    rnt=np.random.rand()
    rnp=np.random.rand()
    rng=np.random.rand()
    # tmp_rng=[[0.32738494873046875,0.17474853992462158,0.96336913108825684],[ 0.84189116954803467,0.10158562660217285,5.4060041904449463E-002]]
    # rnt=tmp_rng[im][0]
    # rnp=tmp_rng[im][1]
    # rng=tmp_rng[im][2]

    theta=rnt*2.0*parmList.pi
    phi=np.arcsin((rnp*2.0)-1.0)+(parmList.pi/2.0)
    gamma=rng*2.0*parmList.pi
    newVector=rotate(vector,theta,phi,gamma,parmList,switch,fout)
    return newVector


def mobil2(coordList,chargeList,fout,parmList,switch,lj,dipole,ion,mconst,iic):
    if switch.im2==0:
        print('\nmobility calculation by MOBIL2 (trajectory method)\n',file=fout)
        print('global trajectory parameters\n\n sw1 =',parmList.sw1,'sw2 =',parmList.sw2,file=fout)
        print('dtsf1 =',parmList.dtsf1,'dtsf2 =',parmList.dtsf2,file=fout)
        print('inwr =',parmList.inwr,'ifail =',parmList.ifail,file=fout)

    switch.SetIt(0)
    switch.SetIu2(0)

    #     determine maximum extent and orientate along x axis
    
    if switch.im2==0:
        print('maximum extent orientated along x axis',file=fout)

    rmax=0
    count=0
    ihold=0
    for ivec in coordList:
        tmp=np.linalg.norm(ivec)
        if tmp>rmax:
            rmax=tmp
            vector=ivec
            ihold=count
        count=count+1
    
    rzy=np.sqrt(vector[2]**2+vector[1]**2)
    phi=np.arccos(vector[2]/rzy)
    phi=phi+(0.5*parmList.pi)  # I want to make the z component =0, thus shit to y-axis
    if vector[1]<0:
        phi=(2*parmList.pi)-phi
    phi=(2.*parmList.pi)-phi
    theta=0.0 
    gamma=0.0
    
    newVector=rotate(coordList,theta,phi,gamma,parmList,switch,fout)
    rotVector=newVector[ihold]
    rxy=np.sqrt(rotVector[0]**2+rotVector[1]**2)
    gamma=np.arccos(rotVector[0]/rxy)
    if rotVector[1]<0:
        gamma=(2*parmList.pi)-gamma
    gamma=(2.*parmList.pi)-gamma

    if(switch.im2==0):
        switch.SetIu3(1)
    if(switch.ip==1):
        switch.SetIu2(1)

    newVector=rotate(coordList,theta,phi,gamma,parmList,switch,fout)
    switch.SetIu3(0)
    switch.SetIu2(0)
    rotVector=newVector[ihold]
    hold=rotVector[0]/rmax

    if hold<0.9999999999 or hold>1.0000000001 or rotVector[1]>1.0e-20 or rotVector[2]>1.0e-20 or rotVector[1]<-1.0e-20 or rotVector[2]<-1.0e-20:
        print('\nProblem orientating along x axis\n',file=fout)
        for ivec in newVector:
            print(ivec[0],ivec[1],ivec[2],np.linalg.norm(ivec),file=fout)

    
    #     determine rmax, emax, and r00 along x, y, and z directions    
    if switch.ip==1:
        print("\n")
    irn=1000
    ddd=(rmax+lj.romax)/float(irn)

    y=0.0
    z=0.0
    emaxx=0.0
    for ir in range(irn):
        x=rmax+lj.romax-(float(ir)*ddd)
        pot,dpotx,dpoty,dpotz,dmax=dljpot(x,y,z,newVector,lj,chargeList,dipole)
        if pot>0:
            break
        r00x=x
        if  pot<emaxx:
            rmaxx=x
            emaxx=pot
    if switch.im2==0:
        print('along x axis emax =',emaxx/parmList.xe,'eV rmax =',rmaxx*1.0e10,'A r00 =',r00x*1.0e10,'A',file=fout)
    
    x=0.0
    z=0.0
    emaxy=0.0
    for ir in range(irn):
        y=rmax+lj.romax-(float(ir)*ddd)
        pot,dpotx,dpoty,dpotz,dmax=dljpot(x,y,z,newVector,lj,chargeList,dipole)
        if pot>0:
            break
        r00y=y
        if  pot<emaxy:
            rmaxy=y
            emaxy=pot
    if switch.im2==0:
        print('along x axis emax =',emaxy/parmList.xe,'eV rmax =',rmaxy*1.0e10,'A r00 =',r00y*1.0e10,'A',file=fout)

    x=0.0
    y=0.0
    emaxz=0.0
    for ir in range(irn):
        z=rmax+lj.romax-(float(ir)*ddd)
        pot,dpotx,dpoty,dpotz,dmax=dljpot(x,y,z,newVector,lj,chargeList,dipole)
        if pot>0:
            break
        r00z=z
        if  pot<emaxz:
            rmaxz=z
            emaxz=pot
    if switch.im2==0:
        print('along x axis emax =',emaxz/parmList.xe,'eV rmax =',rmaxz*1.0e10,'A r00 =',r00z*1.0e10,'A',file=fout)

#     set-up integration over gst
    tst=parmList.xk*ion.temperature/parmList.eo
    if switch.im2==0:
        print('\nt*=',tst,file=fout)
    tst3=tst*tst*tst
    dgst=5e-7*6*np.sqrt(tst)
    gst=dgst
    sum=0.0
    sum2=0.0
    sum1=np.sum(np.sqrt(np.arange(1,parmList.inp+1)))
    if switch.im2==0:
        print('\nset-up gst integration - integration over velocity',file=fout)
        print('\npgst \t\t wgst \t\t\t v \t\t\t ke/kt \t\t gst^5*frac of exp(gst^2/tst)  \t\t sum ',file=fout)
    
    wgst=np.zeros(parmList.vec_dim)
    pgst=np.zeros(parmList.vec_dim)
    #count=0
    for i in np.arange(1,parmList.inp+1):
        hold1=np.sqrt(float(i))
        hold2=np.sqrt(float(i-1))
        sum2=sum2+hold2
        wgst[i-1]=hold1/sum1      # CAREFUL if number. python is 0-based, but fortran is 1-based
        gstt=tst3*(sum2+(hold1/2))/sum1

        while True:
            sum=sum+(np.exp(-gst*gst/tst)*gst*gst*gst*gst*gst*dgst)
            gst=gst+dgst
            if sum>gstt:
                pgst[i-1]=gst-(dgst/2.0)
                break
            if sum==gstt:
                break

        hold1=np.sqrt((pgst[i-1]*pgst[i-1]*parmList.eo)/(0.5*ion.mu))
        hold2=0.5*ion.mu*hold1*hold1/(parmList.xk*ion.temperature)
        hold3=np.exp(-pgst[i-1]*pgst[i-1]/tst)*pgst[i-1]**5

        if switch.im2==0:            
            print(pgst[i-1],wgst[i-1],hold1,hold2,hold3,sum/tst3,file=fout)
            # print(count,"-",pgst[i-1],wgst[i-1],hold1,hold2,hold3,sum/tst3)
            # count=count+1

#     determine b2max
    b2max=np.zeros(parmList.vec_dim)
    cosx=np.zeros(500)
    dbst2=1.0
    dbst22=dbst2/10.0
    cmin=0.0005
    if switch.im2==0:
        print('\nset up b2 integration - integration over impact parameter','\n\nminimum value of (1-cosX) =',cmin,file=fout)
    gst2=np.square(pgst)
    vList=np.sqrt((gst2*parmList.eo)/(0.5*ion.mu))
    for ig in range(parmList.inp-1,-1,-1):
        ibst=np.trunc(rmaxx/parmList.ro)-6
        if ig<parmList.inp-1:
            ibst=int(np.trunc(b2max[ig+1]/dbst2))-6
        if ibst<0:
            ibst=0
        if switch.ip==1:
            print('gst2 =',gst2,'v =',vList[ig],'\nb \t\t bst2 \t\t X ang \t\t cos(X) \t\t e ratio')
        while True:        
            bst2=dbst2*float(ibst)
            b=parmList.ro*np.sqrt(bst2)
            
            ang,erat,d1,istep=gsang(newVector,theta,phi,gamma,chargeList,dipole,lj,vList[ig],b,switch,fout,ion,parmList)
            #print (ibst,b,ang,np.cos(ang))
            cosx[ibst]=1.0-np.cos(ang)
            if(switch.ip==1):
                print(b,bst2,ang,cosx[ibst],erat,file=fout)

            if ibst<5 or ( ibst>=5 and (not (cosx[ibst]<cmin and cosx[ibst-1]<cmin and cosx[ibst-2] < cmin and cosx[ibst-3]<cmin and cosx[ibst-4]<cmin))):
                ibst=ibst+1
                if ibst>=500:
                    print('ibst greater than 500')
                    fout.close()
                    sys.exit(0)
            else:
                break
        b2max[ig]=float(ibst-5)*dbst2
        while True:
            b2max[ig]=b2max[ig]+dbst22
            b=parmList.ro*np.sqrt(b2max[ig])
            
            ang,erat,d1,istep=gsang(newVector,theta,phi,gamma,chargeList,dipole,lj,vList[ig],b,switch,fout,ion,parmList)
            if (1.0-np.cos(ang)<=cmin):
                break
    if switch.im2==0:
        print('\ngst \t\t\t\t b2max/ro2 \t\t\t b/A',file=fout)
        for ig in range(parmList.inp):
            print(pgst[ig],b2max[ig],parmList.ro*np.sqrt(b2max[ig])*1.0e10,file=fout)

#     Calculate Omega(1,1)*, Omega(1,2)*, Omega(1,3)*, and Omega(2,2)*
#     by integrating Q(1)* or Q(2)* over all orientations, and initial 
#     relative velocities.
    if(switch.im2==0):
        print('\n\nnumber of complete cycles (itn) =',parmList.itn,'\nnumber of velocity points (inp) =',parmList.inp,'\nnumber of random points (imp) =',parmList.imp,'\ntotal number of points =',parmList.itn*parmList.inp*parmList.imp,file=fout)
    if switch.ip==1:
        print('start mobility calculation',file=fout)
    
    q1st=np.zeros(parmList.inp)
    q2st=np.zeros(parmList.inp)
    om11st=np.zeros(parmList.itn)
    om12st=np.zeros(parmList.itn)
    om13st=np.zeros(parmList.itn)
    om22st=np.zeros(parmList.itn)

    for ic in range(parmList.itn):
        if switch.ip==1:
            print('\ncycle number, ic =',ic,file=fout)
        # DIFF : next 2 lines unecessary in original code?
        #   gst2=pgst(ig)*pgst(ig)
        #   v=dsqrt((gst2*eo)/(0.5d0*mu))
        for ig in range(parmList.inp):
            if switch.ip==1:
                print('\n\nic =',ic,'ig=',ig,'gst2=',gst2[ig],'v=',vList[ig],file=fout)
            temp1=0.0
            temp2=0.0
            for im in range(parmList.imp):
                if switch.ip==1 and im==0:
                    print('b/A \t\t ang \t \t (1-cosX) \t\t e ratio \t\t theta \t\t phi \t\t gamma',file=fout)
                rnb=np.random.rand()
                newVector=rantate(coordList,parmList,switch,fout)
                bst2=rnb*b2max[ig]
                b=parmList.ro*np.sqrt(bst2)
                if switch.igs==1:
                    ftemp=open('hold')
                    print(iic,ic,ig,im,file=ftemp)
                    print(vList[ig],b,file=ftemp)
                    print(theta*parmList.cang,phi*parmList.cang,gamma*parmList.cang,file=ftemp)
                    ftemp.close()
                ang,erat,d1,istep=gsang(newVector,theta,phi,gamma,chargeList,dipole,lj,vList[ig],b,switch,fout,ion,parmList)
                #print("ang:",ang)
                hold1=1.0-np.cos(ang)
                hold2=np.sin(ang)*np.sin(ang)
                if switch.ip==1:
                    print(b*1.e10,ang*parmList.cang,hold1,erat,theta*parmList.cang,phi*parmList.cang,gamma*parmList.cang,file=fout)
                temp1=temp1+(hold1*b2max[ig]/float(parmList.imp))
                temp2=temp2+(1.5*hold2*b2max[ig]/float(parmList.imp))
            om11st[ic]=om11st[ic]+(temp1*wgst[ig])
            om12st[ic]=om12st[ic]+(temp1*pgst[ig]*pgst[ig]*wgst[ig]*(1.0/(3.0*tst)))
            om13st[ic]=om13st[ic]+(temp1*(pgst[ig]**4)*wgst[ig]*(1.0/(12.0*tst*tst)))
            om22st[ic]=om22st[ic]+(temp2*pgst[ig]*pgst[ig]*wgst[ig]*(1.0/(3.0*tst)))
            q1st[ig]=q1st[ig]+temp1
            q2st[ig]=q2st[ig]+temp2
            if switch.ip==1:
                print('\nv=',vList[ig],'q1st=',q1st[ig],file=fout)
        if switch.ip==1:
            print('\nOMEGA(1,1)*=',om11st[ic],file=fout)

    #     calculate running averages
    temp=1.0/(mconst/(np.sqrt(ion.temperature)*om11st*parmList.pi*parmList.ro2))
    hold1=np.sum(om11st)
    hold2=np.sum(temp)

    if switch.im2==0:
        print('\nsummary of mobility calculations \n\ncycle \t cs/A^2 \t\t avge cs/A^2 \t\t\t Ko^-1 \t \t avge Ko^-1',file=fout)
        for icc in range(1,parmList.itn+1):
            print(icc,om11st[icc-1]*parmList.pi*parmList.ro2*1.e20,hold1*parmList.pi*parmList.ro2*1e20/float(icc),temp[icc-1],hold2/float(icc),file=fout)
        print('\n\naverage values for q1st\n','gst2\t\t\twgst\t\t\tq1st',file=fout)
        for ig in range(parmList.inp):
            print(pgst[ig]*pgst[ig],wgst[ig],q1st[ig]/float(parmList.inp),file=fout)
    
    mom11st=np.sum(om11st)
    mom12st=np.sum(om12st)
    mom13st=np.sum(om13st)
    mom22st=np.sum(om22st)
    mom11st=mom11st/float(parmList.itn)
    mom12st=mom12st/float(parmList.itn)
    mom13st=mom13st/float(parmList.itn)
    mom22st=mom22st/float(parmList.itn)
    sdom11st=0.0
    for item in om11st:
        hold=mom11st-item
        sdom11st=sdom11st+(hold*hold)
    sdom11st=np.sqrt(sdom11st/float(parmList.itn))
    sterr=sdom11st/np.sqrt(float(parmList.itn))

    if switch.im2==0:
        print('\n\nmean OMEGA*(1,1) =',mom11st,'\nstandard deviation =',sdom11st,'\nstandard error of mean =',sterr,file=fout)
    cs=mom11st*parmList.pi*parmList.ro2
    sdevpc=100.0*sdom11st/mom11st
    
    #     Use omegas to obtain higher order correction factor to mobility

    ayst=mom22st/mom11st
    best=((5.0*mom12st)-(4.0*mom13st))/mom11st
    cest=mom12st/mom11st
    term=((4.0*ayst)/(15.0))+(.50*((ion.TotalMass-ion.mGas)**2.0)/(ion.mGas*ion.TotalMass))
    u2=term-(.083330*(2.40*best+1.0)*(ion.mGas/ion.TotalMass))
    w=(ion.mGas/ion.TotalMass)
    delta=((((6.0*cest)-5.0)**2.0)*w)/(60.0*(1.0+u2))
    f=1.0/(1.0-delta)
    if switch.im2==0:
        print('\n\nf value for second order correction=',f,file=fout)
        print('(integrations for second order correction are not',file=fout)
        print('accurate, check if correction becomes significant)',file=fout)
        print('\nomega*12 =',mom12st,'\tomega*13 =',mom13st,'\tomega*22 =',mom22st,file=fout)
        print('    u2 =',u2,'\t       w =',w,'\t   delta =',delta,file=fout)
    mob=(mconst*f)/(np.sqrt(ion.temperature)*cs)
    print('\n\naverage (second order) TM mobility =',mob,file=fout)
    print('inverse average (second order) TM mobility =',1.0/mob,file=fout)
    print('average TM cross section =',cs*1.e20,file=fout)

    return mob,cs,sdevpc

def che(im,coordList,lj,cop,versor,parmList,yr,zr,kp):
    '''
c     Guides hard sphere scattering trajectory. Adapted from code 
c     written by Alexandre Shvartsburg.

    '''

# c     If this is a secondary collision, the incident vector lies on the
# c     x axis in transformed coordinates. 
    if im!=1:
        yr=0
        zr=0
    
#c     xl is a large number greater than any cluster dimension.
    xl=1e+6
    #DIFF: in the original coda ki is the index of the vector that meets the criteria below.
    #Here ki is 0 or 1 for whether it meets the criteria or not.
    # the vector that meets the criteria is stored directly in hs_vector
    ki=0
    hs_vector=np.zeros(3)
    for item,hs,hs2 in zip(coordList,lj.rhs,lj.rhs2):
        if item[0]>1e-16 or im==1:

            #     yd and zd are the coordinates of the impact points for atom in
            #    with respect to its own coordinates (if such a point exists).
            #     dev is the impact parameter.   

            yd=yr-item[1]
            zd=zr-item[2]
            ras=(yd*yd)+(zd*zd)
            dev=np.sqrt(ras)
            # If the collision with in-th atom occurs, then
            if dev<=hs:
                #Find xc - the x coordinate of collision point with the in-th atom.
                xc=item[0]-np.sqrt(hs2-ras)
                # Find the smallest xc (earliest collision).
                if xc<xl:
                    xl=xc
                    ki=1
                    hs_vector=item
                    hs_hs2=hs2

    #If a collisions took place, then xv, yv, and zv are the vectors
    #going from the center of ki-th atom to the collision point.
    if ki!=0:
        kp=1
        xv=xl-hs_vector[0]
        yv=yr-hs_vector[1]
        zv=zr-hs_vector[2]
        #Transform all coordinates by a parallel move such that the 
        # collision point is at (0,0,0) in new coordinates.
        for item in coordList:
            item[0]=item[0]-xl
            item[1]=item[1]-yv
            item[2]=item[2]-zv
    # c     Transform the coordinates of all atoms such that the direction
    # c     vector of the reflected particle is collinear to axis x. (Note
    # c     that the number of such transformations is infinite, thus the
    # c     direction of the y or z axis is arbitrary.) The direction cosines
    # c     of the incoming ray are also transformed. xve1, xve2, and xve3 
    # c     are the direction vectors of reflected ray in the coordinate 
    # c     system of incoming ray. 

        #Evaluate the transformation matrix elements.
        xxv=2.0*xv*xv
        xyv=2.0*xv*yv
        xzv=2.0*xv*zv
        xyz=xyv*xzv
        rad2=hs_hs2-xxv
        ad1=(rad2*rad2)+(xyv*xyv)
        adr1=np.sqrt(ad1)
        adr2=np.sqrt(ad1*ad1+xyz*xyz+xzv*xzv*rad2*rad2)
        xve1=1.0-2.0*xv*xv/hs_hs2
        xve2=-xyv/hs_hs2
        xve3=-xzv/hs_hs2
        yve1=xyv/adr1
        yve2=rad2/adr1
        yve3=0.0
        zve1=rad2*xzv/adr2
        zve2=-xyz/adr2
        zve3=ad1/adr2
        # NOTE: 997 format(15x,1pe12.5,1x,e12.5,1x,e12.5) --Unecessary?
        # Transform the coordinates and direction cosines of incoming ray.
        for item in coordList:
            xne=xve1*item[0]+xve2*item[1]+xve3*item[2]
            yne=yve1*item[0]+yve2*item[1]+yve3*item[2]
            item[2]=zve1*item[0]+zve2*item[1]+zve3*item[2]
            item[0]=xne
            item[1]=yne
        
        xne=xve1*versor[0]+xve2*versor[1]+xve3*versor[2]
        yne=yve1*versor[0]+yve2*versor[1]+yve3*versor[2]
        versor[2]=zve1*versor[0]+zve2*versor[1]+zve3*versor[2]
        versor[0]=xne
        versor[1]=yne
        cof=np.cos(0.50*(parmList.pi-np.arccos(versor[0])))
    else:
        cof=cop

    return cof,coordList,versor,yr,zr,kp

def mobil4(coordList,chargeList,fout,parmList,switch,lj,ion,mconst):
    '''
c     Calculates the collision cross-section using the projection 
c     approximation and the exact hard spheres scattering model. Adapted
c     from code written by Alexandre Shvartsburg.

    '''

# c     inum is the number of Monte Carlo trajectories employed.
# c     inor is the maximum number of successive reflections followed 
# c     for each trajectory. cof stores the cosines of halves of angles
# c     between the incidence and reflected vectors (for successive 
# c     collisions along a single trajectory). The corresponding collision
# c     cross sections are accumulated in crof.
# c     Initialize the cross-section arrays to zero.
# c     crb is the projection (i.e. "hard-sphere cross-section").
# c     imm is the highest collision order encountered for any trajectory

    crof=np.zeros(parmList.inor)
    crb=0.0
    imm=1
    newVector=coordList

# Monte Carlo starts here

    for i in range(parmList.inum):
        newVector=rantate(newVector,parmList,switch,fout)
        y1=newVector[:,1]+lj.rhs
        y2=newVector[:,1]-lj.rhs
        z1=newVector[:,2]+lj.rhs
        z2=newVector[:,2]-lj.rhs
        ymin=np.amin(y2)
        ymax=np.amax(y1)
        zmin=np.amin(z2)
        zmax=np.amax(z1)
        
        # ydi and zdi are the lengths of the box sides.
        ydi=ymax-ymin
        zdi=zmax-zmin
        # pls is the area of the box.
        pls=ydi*zdi
        # yr and zr are the random coordinates inside the box.
        yr=ymin+ydi*np.random.rand()
        zr=zmin+zdi*np.random.rand()
        # Write the direction cosines of the initial incidence vector 
        # (collinear to the x axis).
        # DIFF: The original code added this vector as an extra coordinate.

        versor=[1,0,0]
        kp=0
        cof=np.zeros(parmList.inor+2)
        for im in range(1,parmList.inor+1):
            cof[im+1],newVector,versor,yr,zr,kp=che(im,newVector,lj,cof[im],versor,parmList,yr,zr,kp)
        # c     If the inclusion of the next order of collision had not changed
        # c     incidence/reflection angle (i.e., no collision occurred) stop 
        # c     following this trajectory.
            if cof[im+1]==cof[im]:
                if im-1 >imm:
                    imm=im-1
                    #     Set all further incidence/reflection angle cosines along this
                    #     trajectory to the last angle found.

                for imn in range (im+1,parmList.inor+1):
                    cof[imn+1]=cof[im+1]
                break
                
            else:
                if im> imm:
                    imm=im  
    #     Add the contributions from the i-th trajectory to the collision 
    #     cross sections of all orders.
        for im in range(1,parmList.inor+1):
            crof[im-1]=crof[im-1]+pls*cof[im+1]*cof[im+1]
        if kp==1:
            crb=crb+pls
    # END of Monte Carlo
    u2=0.5*float(parmList.inum)
    #Normalize all collision cross-sections and the projection.
    crof=crof/u2
    crb=crb/float(parmList.inum)
    pacs=crb
    pamob=mconst/(np.sqrt(ion.temperature)*pacs)
    ehscs=crof[imm-1]
    ehsmob=mconst/(np.sqrt(ion.temperature)*ehscs)

    #Output the results of calculation.
    if switch.im4==0:
        print('\nmobility calculation by MOBIL4 (HS scattering)',file=fout)
        print('\nnumber of Monte Carlo trajectories =',parmList.inum,file=fout)
        print('\nmaximum number of reflections allowed =',parmList.inor,file=fout)
    print('\naverage PA cross section =',pamob,file=fout)
    print('\ninverse average PA mobility =',1.0/pamob,file=fout)
    print('\naverage PA cross section =',pacs*1.e20,file=fout)
    print('\nmaximum number of reflections encountered =',imm,file=fout)
    

    if switch.im4==0:
        print('\norder \t cross section',file=fout)
        for im in range(imm):
            print(im,crof[im],file=fout)

    print('\naverage EHS cross section =',ehsmob,file=fout)
    print('\ninverse average EHS mobility =',1.0/ehsmob,file=fout)
    print('average EHS cross section =',ehscs*1.e20,file=fout)
    

    return ehscs,ehsmob,pacs,pamob,imm
    

        
#############################    
##           MAIN          ##
#############################
#check if file is actually there
try:
    fin=open('mobcal.run',"r")
except IOError as e:
    print("Couldn't open file (%s)." % e)
    #return 0

tmp=fin.readline()
x=tmp.split()
inFileName=x[0]

tmp=fin.readline()
x=tmp.split()
outFileName=x[0]

tmp=fin.readline()
if tmp=='':
    Rseed=None
else:
    x=tmp.split()
    if x[0]=='':
        Rseed=None
    else:
        Rseed=x[0]

# inFileName="a10A1.mfj"
# outFileName="a10A1.out"
# Rseed=None

ConstList=constants()
printDebug=printSwitch()
#id is used to associated atom type to atom "mass" in input file. In mobcal atoms are identified by their mass and not symbol
AtomList=[
    atom('Hydrogen',1.008,0.65e-03*ConstList.xe,2.38e-10,2.2e-10,1),
    atom('Carbon',12.01,1.34e-3*ConstList.xe,3.043e-10,2.7e-10,12),
    atom('Nitrogen',14.01,1.34e-3*ConstList.xe,3.043e-10,2.7e-10,14),
    atom('Oxygen',16.00,1.34e-3*ConstList.xe,3.043e-10,2.7e-10,16),
    atom('Sodium',22.99,0.0278e-3*ConstList.xe,(3.97/1.12246)*1e-10,2.853e-10,23),
    atom('Silicon',28.09,1.35e-3*ConstList.xe,3.5e-10,2.95e-10,28),
    atom('Sulfur',32.06,1.35e-3*ConstList.xe,3.5e-10,3.5e-10,32),
    atom('Iron',55.85,1.35e-3*ConstList.xe,3.5e-10,3.5e-10,56)
]

# LJ scaling parameters
atoms,ion,InputData,tmp_asymp,lj,hs,fin,fout=fcoord(inFileName,outFileName,ConstList,AtomList,printDebug)

print('Lennard-Jones scaling parameters: eo= ',ConstList.eo/ConstList.xe,'ro= ',ConstList.ro*1e+10,file=fout)


    #  Constant for ion-induced dipole potential

    #  Polarizability of helium = 0.204956d-30 m3
    #  xeo is permitivity of vacuum, 8.854187817d-12 F.m-1

dipole=ConstList.xe**2*0.204956e-30/(8*ConstList.pi*ConstList.xeo)
ion.SetDipole(dipole)
print('dipole constant =',ion.dipole,file=fout)

    # Mass constants
ion.SetMassGas(4.0026,ConstList.xn*1000)

    # Mobility constant

mconst=np.sqrt(18*ConstList.pi)/16
mconst=mconst*np.sqrt(ConstList.xn*1000)*np.sqrt((1/ion.mGas)+(1/ion.TotalMass))
mconst=mconst*ConstList.xe/np.sqrt(ConstList.xk)
dens=ConstList.xn/ConstList.xmv
mconst=mconst/dens
ion.SetMobility(mconst)
print('mobility constant =',ion.mobility,file=fout)
print('temperature =',ion.temperature,file=fout)


#DIFF: In the original code follows an initialization of the random number generator.
# Here I am using numpy, so that section is unnecessary

np.random.RandomState(seed=Rseed)

tmm=np.zeros(InputData.icoord)
tmc=np.zeros(InputData.icoord)
asympp=np.zeros(InputData.icoord)
ehsc=np.zeros(InputData.icoord)
ehsm=np.zeros(InputData.icoord)
pac=np.zeros(InputData.icoord)
pam=np.zeros(InputData.icoord)
asympp[0]=tmp_asymp

####
imm=0
tmmob=0
immmin=0
immmax=0
###

InputData.icoord=2
for iic in range(InputData.icoord):
    print("",file=fout)
    if InputData.icoord > 0:
        print('\ncoordinate set =',iic+1,file=fout)
    print('\nstructural asymmetry parameter =',asympp[iic],file=fout)

    tmm[iic],tmc[iic],sdevpc=mobil2(atoms.Coords,atoms.Charges,fout,ConstList,printDebug,lj,ion.dipole,ion,mconst,iic)
    ehsc[iic],ehsm[iic],pac[iic],pam[iic],imm=mobil4(atoms.Coords,atoms.Charges,fout,ConstList,printDebug,hs,ion,mconst) 

    if(imm<immmin):
        immmin=imm
    if(imm>immmax):
        immmax=imm  
   
    if iic<InputData.icoord-1:
        Oxyz,asympp[iic+1],fin,fout=ncoord(fin,fout,ConstList,AtomList,printDebug,ion,InputData)
        atoms.UpdatePosition(Oxyz)
    printDebug.im2=1
    printDebug.im4=1

if not fin.closed:
    fin.close()
#     print out summary
print('\n\n\n SUMMARY\n',file=fout)
print('program version = 0.1',file=fout)
print('input file name =',inFileName,file=fout)
print('input file label = ',InputData.label,file=fout)
if tmmob!=0:
    if InputData.charges=='equal':
        fout.write("using a uniform charge distribution\n")
    elif InputData.charges=='calc':
        fout.write("using a calculated (non-uniform) charge distribution\n")
    elif InputData.charges=='none':
        fout.write("using no charge - only LJ interactions\n")
print('temperature = ',ion.temperature,file=fout)
print('using NUMPY random number generator with seed  =',Rseed,file=fout)

if InputData.icoord==1:
    print('structural asymmetry parameter =',asympp[0],file=fout)
    if ehsm[0]!=0:
        print('\nmobility calculation by MOBIL4 (HS scattering)',file=fout)
        print('\nnumber of Monte Carlo trajectories =',ConstList.inum,file=fout)
        print('\nmaximum number of reflections encountered =',imm,file=fout)
        print('\ninverse average PA mobility =',1.0/pam[0],file=fout)
        print('\naverage PA cross section =',pam[0]*1.e20,file=fout)
        print('\ninverse average EHS mobility =',1.0/ehsm[0],file=fout)
        print('\naverage EHS cross section =',ehsm[0]*1e20,file=fout)
    else:
        if tmm[0]==0:
            fout.close()
            sys.exit(0)
        print('\nmobility calculation by MOBIL2 (trajectory method)',file=fout)
        print('\ntrajectory parameters',file=fout)
        print('sw1 =',ConstList.sw1,'\tsw2 =',ConstList.sw2,file=fout)
        print('dtsf1 =',ConstList.dtsf1,'\tdtsf2 =',ConstList.dtsf2,file=fout)
        print('\n\nnumber of complete cycles (itn) =',ConstList.itn,file=fout)
        print('\nnumber of velocity points (inp) =',ConstList.inp,file=fout)
        print('\nnumber of random points (imp) =',ConstList.imp,file=fout)
        print('\ntotal number of points =',ConstList.itn*ConstList.inp*ConstList.imp,file=fout)
        print('\ninverse average (second order) TM mobility =',1.0/tmm[0],file=fout)
        print('\naverage TM cross section =',tmc[0]*1.e20,file=fout)
        print('\nstandard deviation (percent) =',sdevpc,file=fout)
        print('\nnumber of failed trajectories =',ConstList.ifailc,file=fout)

#     multiple coordinate sets
else:
    if ehsm[0]!=0:
        print('\nmobility calculation by MOBIL4 (HS scattering)',file=fout)
        print('\nnumber of Monte Carlo trajectories =',ConstList.inum,file=fout)
        print('\nminimum and maximum number of reflections =',immmin,immmax,file=fout)
        print('set \t PA CS \t\t PA MOB^-1 \t\t EHSS CS \t\t EHSS MOB^-1 \t ASYMP',file=fout)
        for i in range(InputData.icoord):
            print(i,pac[i]*1.0e20,1.0/pam[i],ehsc[i]*1.0e20,1.0/ehsm[i],asympp[i]/10.0,file=fout)
        pacs=np.sum(pac)/float(InputData.icoord)
        pamob=np.sum(pam)/float(InputData.icoord)
        ehscs=np.sum(ehsc)/float(InputData.icoord)
        ehsmob=np.sum(ehsm)/float(InputData.icoord)
        aasymp=np.sum(asympp)/float(InputData.icoord)
        print('\nAVGE',pacs*1.0e20,1.0/pamob,ehscs*1.0e20,1.0/ehsmob,aasymp/10.0,file=fout)
    else:
        if tmm[0]==0:
            fout.close()
            sys.exit(0)
        print('\nmobility calculation by MOBIL2 (trajectory method)',file=fout)
        print('\ntrajectory parameters',file=fout)
        print('sw1 =',ConstList.sw1,'\tsw2 =',ConstList.sw2,file=fout)
        print('dtsf1 =',ConstList.dtsf1,'\tdtsf2 =',ConstList.dtsf2,file=fout)
        print('\n\nnumber of complete cycles (itn) =',ConstList.itn,file=fout)
        print('\nnumber of velocity points (inp) =',ConstList.inp,file=fout)
        print('\nnumber of random points (imp) =',ConstList.imp,file=fout)
        print('\ntotal number of points =',ConstList.itn*ConstList.inp*ConstList.imp,file=fout)
        print('\nnumber of failed trajectories =',ConstList.ifailc,file=fout)
        print('\n\nset  \t\t TM CS \t\t\t TM MOB^-1')
        for i in range(InputData.icoord):
            print(i,tmc[i]*1.0e20,1.0/tmm[i])
        tmcs=np.sum(tmc)/float(InputData.icoord)
        tmmob=np.sum(tmm)/float(InputData.icoord)
        print('\nAVGE',tmcs*1.0e20,1.0/tmmob)

fout.close()
