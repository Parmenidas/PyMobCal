'''
Test of the function fcoord
Note that in this code we do *not* call the fortran function,
But instead we run both codes independently and compare outputs
'''

import math
import numpy as np
import scipy
from scipy.spatial.transform import Rotation as R
import sys

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
        self.inp=40   
        self.itn=10 
        self.imp=25
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
        self.v=2309.9              #
        self.b=1.5654e-10          #
        self.ntheta=33.3
        self.nphi=64.250
        self.ngamma=134.30
        self.ehsm=0
        self.tmm=0
        self.immmax=0
        self.immmin=self.inor
        self.traj_lost=30000        #DIFF: This was hard-coded into the function. I believe this is the right place
        self.vec_dim=100            #DIF: this is what is used within the FORTRAN function as dimension of arrays 


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
        self.ro2lj=np.square(rolj)
        self.ro6lj=np.power(self.ro2lj,3)
        self.ro12lj=np.square(self.ro6lj)
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
    for element in AtomList:
        if element.id==id:
            return element.Unpack()
    else:
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

def fcoord(fileIN,fileOUT,parmList,AtomList,switch):
    '''
    Read the input file. The format is the original mobcal format
    '''
    #TODO: fix error message/behaviour. Maybe raise exceptions?
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

    print(" input file name = {:30s}".format(fileIN),file=fout)
    tmp=fin.readline()
    inList["label"]=tmp.rstrip()
    print(" input file label = {:30s}".format(inList["label"]),file=fout)
    tmp=fin.readline()
    inList["icoord"]=int(tmp)
    print(" number of coordinate sets ={:5d}".format(inList["icoord"]),file=fout)
    tmp=fin.readline()
    inList["atoms"]=int(tmp)
    print(" number of atoms ={:4d}".format(inList["atoms"]),file=fout)
    
    tmp=fin.readline()
    x=tmp.split()
    inList["units"]=x[0]
    if x[0]=='au':
        print(" coordinates in atomic units",file=fout)
    elif x[0]=='ang':
        print(" coordinates in angstroms",file=fout)
    else:
        print(" coordinates selected not supported",file=fout)
        #return 0

    # Here we write the correction factor *before* charges.
    # Not sure why, but this is done in the original code
    # and I am retaining compatibility
    tmp=fin.readline()
    x=tmp.split()
    inList["charges"]=x[0]

    tmp=fin.readline()
    inList["cfact"]=float(tmp)
    print(" correction factor for coordinates ={: .4E}".format(inList["cfact"]),file=fout)

    if x[0]=='equal':
        print(" using a uniform charge distribution",file=fout)
    elif x[0]=='calc':
        print(" using a calculated (non-uniform) charge distribution",file=fout)
    elif x[0]=='none':
        print(" using no charge - only LJ interactions",file=fout)
    else:
        print("Charges model selected not supported",file=fout)
        fin.close()
        fout.close()
        sys.exit(0)
        #return 0

    #Setting up some arrays before reading coords. Try to be consistent with original code as well as optimize a bit
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

    #Reading atom coordinates starts here:
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
    print(" total charge ={: .4E}\n total absolute charge ={: .4E}".format(tcharge,acharge),file=fout)
    ion.SetCharges(tcharge,acharge)

    for i in range(inList["atoms"]):
        xmass[i],eolj[i],rolj[i],rhs[i]=GetAtom(AtomList,parmList,imass[i])
    tMass=np.sum(xmass)
    # In the next line the replacement of E with D is a hack for compatibility.
    # Briefly, fortran uses D instead of E to point out that you are using
    # the scientific notation to represent an integer number and not a float
    print(" mass of ion ={: .4E}".format(tMass).replace('E','D'),file=fout)
    ion.SetMass(tMass)
    
    lj=LennardJones(eolj,rolj)
    hs=HardSphere(rhs)

    if switch.iu1==1:
        x=' '
        print('\n'+9*x+'initial coordinates'+9*x+'mass'+3*x+'charge'+9*x+'LJ parameters\n',file=fout)
    
    #The center of mass is used as origin of the reference frame
    fxo=np.sum(np.multiply(fx,xmass))/tMass
    fyo=np.sum(np.multiply(fy,xmass))/tMass
    fzo=np.sum(np.multiply(fz,xmass))/tMass
    print(' center of mass coordinates = {: .4E},{: .4E},{: .4E}'.format(fxo,fyo,fzo),file=fout)

    Oxyz=np.zeros([inList["atoms"],3])   # original coordinates vector
    for i in range(inList["atoms"]):
        Oxyz[i][0]=(fx[i]-fxo)*1.e-10*inList["cfact"]
        Oxyz[i][1]=(fy[i]-fyo)*1.e-10*inList["cfact"]
        Oxyz[i][2]=(fz[i]-fzo)*1.e-10*inList["cfact"]

    if switch.iu1==1:
        for i in range(inList["atoms"]):
            print(' {: .4E} {: .4E} {: .4E} {:3d} {: .4E} {: .4E} {: .4E}'.format(Oxyz[i][0],Oxyz[i][1],Oxyz[i][2],imass[i],pcharges[i],eolj[i]/parmList.xe,rolj[i]*1e+10),file=fout)
        print("\n",file=fout)
    
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
#############################    
##           MAIN          ##
#############################
inFileName="a10A1.mfj"
outFileName="a10A1_py.out"
Rseed=43

ConstList=constants()
printDebug=printSwitch()
#id is used to associated atom type to atom "mass" in input file. In mobcal atoms are identified by their mass and not symbol
AtomList=(
    atom('Hydrogen',1.008,0.65e-03*ConstList.xe,2.38e-10,2.2e-10,1),
    atom('Carbon',12.01,1.34e-3*ConstList.xe,3.043e-10,2.7e-10,12),
    atom('Nitrogen',14.01,1.34e-3*ConstList.xe,3.043e-10,2.7e-10,14),
    atom('Oxygen',16.00,1.34e-3*ConstList.xe,3.043e-10,2.7e-10,16),
    atom('Sodium',22.99,0.0278e-3*ConstList.xe,(3.97/1.12246)*1e-10,2.853e-10,23),
    atom('Silicon',28.09,1.35e-3*ConstList.xe,3.5e-10,2.95e-10,28),
    atom('Sulfur',32.06,1.35e-3*ConstList.xe,3.5e-10,3.5e-10,32),
    atom('Iron',55.85,1.35e-3*ConstList.xe,3.5e-10,3.5e-10,56)
)

printDebug.iu1=1

atoms,ion,InputData,tmp_asymp,lj,hs,fin,fout=fcoord(inFileName,outFileName,ConstList,AtomList,printDebug)

if not fin.closed:
    fin.close()

fout.close()
