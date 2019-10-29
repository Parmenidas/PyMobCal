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


class AtomData:
    '''
    Contains information about each atom within the molecule
    * Attributes
    * ----------
    * coord : coordinates
    * charge : charges
    * mass : masses

    '''
    def __init__(self,mass=[],charge=[],coord=[]):
        self.Mass=mass
        self.Charges=charge
        self.Coords=coord
    
    def UpdatePosition(self,coord=[]):
        self.Coords=coord

    def UpdateCharges(self,charges=[]):
        self.Charges=charges

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

class infile:
    '''
    * This structure contains some general information obtained from the
    * input file. Ideally, no matter what the actually file format looks like,
    * the function reading it should populate this structure.
    * 
    * Attributes
    * ----------
    * label : a user-provided label for the simulation. Limit to MAX_LINE
    * n_atoms : number of atoms per molecules
    * n_coords : Sets of coordinates. Ideally, different structures of the same molecules. 
    *            *not* different orientation as the program takes care of that.
    * units : the units used for the coordinates of the atoms. AU=atomic units, ANG=angstroms
    * chargeUnits : How the charges are generated. NONE = no charges, EQUAL = all charges equal 1/n_Atoms, 
    *               CALC = read from file (an extra column must be present in the file)
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
    * Structure containing information about a specific atomtype
    * 
    * Attributes
    * ----------
    * name : name of the atom
    * mass : atomic mass
    * epsilonLj,sigmaLJ : Lennard-Jones parameters
    * rhs : parameters for hard sphere scattering
    * id : atomic number
    * 
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
        self.inum=250000
        self.inor=30
        
        # self.ipr=2
        # self.inp=5
        # self.itn=2 
        # self.imp=2
        # self.inum=10
        # self.inor=10

        self.cmin=0.0005
        self.sw1=0.00005
        self.sw2=0.005 
        self.dtsf1=0.5
        self.dtsf2=0.1
        self.inwr=1
        self.ifail=100
        self.ifailc=0              #This should probably be brought outside. It is changed by the code, so not a constant
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
        self.ibst_max=500           #DIFF: This was hard-coded into the function. I believe this is the right place

class printSwitch:
    '''    
    Some switches used to turn on/off some debugging features.
    Unless you are modifying the code, keep defaults

c     print switches ip=1  print scattering angles
c                    it=1  print trajectory
c                    iu1=1 print initial coordinates from FCOORD
c                    iu2=1 print angles and coordinates from ROTATE 
c                    iu3=1 print angles from ROTATE
c                    iv=1  print all potentials from POTENT
c                    im2=1 print brief information from MOBIL2
c                    im4=1 print brief information from MOBIL4
c                    igs=1 print out angles, v, and b from MOBIL4 into 
c                          a temporary file called hold
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
    /* 
    * Structure containing information for calculating Lennard Jones interactions
    * The code defines the potential as U = 4 * epsilon *[ (sigma/r)^12 - (sigma/r)^6 ]
    * 
    * Attributes
    * ----------
    * eolj : list of epsilon (one for each atom)
    * rolj : list of sigmas (one for each atom)
    * n : number of atoms
    * romax : largest value in rolj
    * 
    * The following attributes are derived from the previous.
    * As they are used over and over, we calculate them once and reuse.
    * 
    * eox4 : 4 * eolj (one for each atom)
    * ro2lj : square of rolj (one for each atom)
    * ro6lj : rolj^6 (one for each atom)
    * ro12lj : rolj^12 (one for each atom)
    * dro6 : 6 * ro6lj (one for each atom)
    * dro12 : 12 * ro12lj (one for each atom)    
    */    

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
    * Structure containing information for calculating hard sphere scattering
    * 
    * Attributes
    * ----------
    * rhs : list of radii  (one for each atom)
    * rhs2 : square of rhs. It is used over and over, so we compute and reuse to speed up calculations
    '''
    def __init__(self,rhs):
        self.rhs=rhs
        self.rhs2=np.multiply(rhs,rhs)

def GetAtom(AtomList,ConstList,id):
    '''
    * This function assign atom properties based on atom atomic number
    This function comapares the mass found in teh file with the database and set all
    the internal variable for that atom.
    If the atom is NOT in the database it is defaulted to a carbon atom and a warning is issued.
    Whether using a carbon atom as default is a good idea or not, I will let you decide

    * 
    * Arguments
    * ---------
    * id : atomic number
    * ConstList : list of important constants
    * AtomList : a list with the atoms recognized by mobcal and their proeprties
    * 
    * Returns
    * -------
    * mass : atom mass
    * epsilonLJ, sigmaLJ : Lennard-Jones parameters
    * rhs : hard shpere radius 

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
    * 
    * Arguments
    * ---------
    * vector : the vector(s) to  be rotated
    * theta,phi,gamma : the angles used to rotate vector
    * ConstList : List with the constant used in the simulations
    * switch : list with printing switches
    * fout : file for logging
    * 
    * Returns
    * -------
    * newVector : the vector(s) after the rotation

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

    Parameters
    ----------
    fileIN : name of the input file with the information to be read
    fileOUT : name of the file for logging
    parmList : list of important constants
    AtomList : a list with the atoms recognized by mobcal and their proeprties
    switch : switches for printing verbosity

    Returns
    -------
    atoms : structure with position, charges and masses for *each* atom
    ion : structure with information about the *whole* molecule
    InputData : structure with parameters set by the input file
    asymp : structural asymmetry parameter
    lj : structure with Lennard-Jones parameters
    hs : structure with hard-sphere parameters
    fin : handle for the input file 
    fout : handle for the output file

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

def ncoord(fin,fout,parmList,AtomList,switch,ion,inputList):
    '''
    Read coordinates after the first set (if any). The format is the original mobcal format

    Parameters
    ----------
    fin : input file with the information to be read
    fout : file for logging
    parmList : list of important constants
    AtomList : a list with the atoms recognized by mobcal and their proeprties
    switch : switches for printing verbosity
    ion : structure with information about the *whole* molecule
    InputData : structure with parameters set by the fcoord()

    Returns
    -------
    Oxyz,pcharges,asymp,fin,fout
    Oxyz : new position of the atoms
    pcharges : new charges of each atom
    asymp : structural asymmetry parameter
    fin : handle for the input file 
    fout : handle for the output file

    '''
    
    fx=np.zeros(inputList.atoms)
    fy=np.zeros(inputList.atoms)
    fz=np.zeros(inputList.atoms)
    imass=np.zeros(inputList.atoms)
    xmass=np.zeros(inputList.atoms)
    pcharges=np.zeros(inputList.atoms)

    tmp=fin.readline()   #unused line

    #Reading atoms coordinates starts here:
    for i in range(inputList.atoms):
        tmp=fin.readline()
        x=tmp.split()
        fx[i]=float(x[0])
        fy[i]=float(x[1])
        fz[i]=float(x[2])
        imass[i]=int(np.rint(float(x[3])))
        if inputList.charges=='calc':
            pcharges[i]=float(x[4])


    if inputList.units=='au':
        fx=fx*0.52917706
        fy=fy*0.52917706
        fz=fz*0.52917706

    for i in range(inputList.atoms):
        xmass[i],eolj,rolj,rhs=GetAtom(AtomList,parmList,imass[i])
    
    tMass=np.sum(xmass)
    if tMass!=ion.TotalMass:
        print(' masses do not add up',file=fout)
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

    return Oxyz,pcharges,asymp,fin,fout


def dljpot(position,vector,lj,charge,dipol):
    '''
    This function computes the force and energy between one particle
    of the gas and the molecule.

    In the TM method, one particle is launched against the molecule and its
    scattering is recorded. This function computes the force and energy at each
    time step and drives the motion of the particle.

    Parameters
    ----------
    position: position of the particle. it is a vector with three values for x,y,z
    vector: this vector contains the coordinates of the molecule
    lj: this is a structure that contains the Lennard-Jones parameters of each atom of the molecules
    charge: vector containing the charges of each atom in the molecule
    dipol: total dipole of the molecule

    Returns
    -------
    pot: Total potential energy
    dpotx,dpoty,dpotz: derivatives of the potential with respect x,y and z
    dmax: maximum distance between particle and molecule
    '''
    dmax=2.*lj.romax
    
    # Compute distance between molecule and gas atom
    # data is a vector containing the distance along *each* coordinate 
    # of *each* atom in the molecule with the gas molecules
    data=np.subtract(position,vector)
    data2=np.square(data)

    # r is the distace 
    r2=np.sum(data2,axis=1)
    r=np.sqrt(r2)

    # Define some quantities to be used later on for computing the actual
    # forces and potentials. We compute once and reuse
    r_min=np.amin(r)
    dmax=r_min if (r_min<dmax) else dmax
    
    r3=np.multiply(r2,r)
    r5=np.multiply(r3,r2)
    r6=np.multiply(r5,r)
    r8=np.multiply(r5,r3)
    r12=np.multiply(r6,r6)
    r14=np.multiply(r12,r2)

    # Lennard Jones potential
    e00=np.sum(lj.eox4*((lj.ro12lj/r12)-(lj.ro6lj/r6)))
    
    # Lennard Jones forces
    de00=lj.eox4*((lj.dro6/r8)-(lj.dro12/r14))
    de00x=np.sum(np.multiply(de00,data[:,0]))
    de00y=np.sum(np.multiply(de00,data[:,1]))
    de00z=np.sum(np.multiply(de00,data[:,2]))

    # Total potential energy from Lennard Jones
    pot=e00
    #  Derivatives of the potential with respect coords = - derivatives of the  momenta  with respect time
    dpotx=de00x
    dpoty=de00y
    dpotz=de00z

    # dipolar interactions (if any)
    if np.sum(charge): # TODO bring in total absolute charge. No reasons to recalculate it over and over
        qr3=charge/r3
        qr5=-3.*charge/r5
        rx=np.sum(np.multiply(qr3,data[:,0]))
        ry=np.sum(np.multiply(qr3,data[:,1]))
        rz=np.sum(np.multiply(qr3,data[:,2]))

        sum1=np.sum(qr3+(data2[:,0]*qr5))
        sum2=np.sum(data[:,0]*data[:,1]*qr5)
        sum3=np.sum(data[:,0]*data[:,2]*qr5)
        sum4=np.sum(qr3+(data2[:,1]*qr5))
        sum5=np.sum(data[:,1]*data[:,2]*qr5)
        sum6=np.sum(qr3+(data2[:,2]*qr5))

        # Total potential energy
        pot=e00-(dipol*((rx*rx)+(ry*ry)+(rz*rz)))
        #  Derivatives of the potential with respect coords = derivatives of the  momenta  with respect time
        dpotx=de00x-(dipol*((2.0*rx*sum1)+(2.*ry*sum2)+(2.*rz*sum3)))
        dpoty=de00y-(dipol*((2.0*rx*sum2)+(2.*ry*sum4)+(2.*rz*sum5)))
        dpotz=de00z-(dipol*((2.0*rx*sum3)+(2.*ry*sum5)+(2.*rz*sum6)))
    
    return pot,dpotx,dpoty,dpotz,dmax


def diffeq(tim,dt,position,momenta,parmList,ion,coord,lj,chargeList,dipole,Dposition,Dmomenta):
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
    pot,dpotx,dpoty,dpotz,dmax=dljpot(position,coord,lj,chargeList,dipole)

    Dmomenta=np.array([-dpotx,-dpoty,-dpotz])
 
    #finish step velocity
    momenta=momenta+0.5*dt*Dmomenta
    Dposition=momenta*ion.invMu

    tim=tim+dt
 
    return position,momenta,Dposition,Dmomenta,dt,tim,pot,dpotx,dpoty,dpotz,dmax

def deriv(position,momenta,ion,coord,lj,chargeList,dipole):
    '''
    /*
    * This function computes the derivative with respect time of position and momenta
    * 
    * Arguments
    * ---------
    * position : vector with the position of the gas particle. *must* be of dimesion = 3
    * momenta : vector with the linear  momentum (that is, velocity * mass) of the gas particle. *must* be of dimesion = 3
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
    pot,dpotx,dpoty,dpotz,dmax=dljpot(position,coord,lj,chargeList,dipole)

    # Update derivative of the momenta with respect to time
    Dmomenta=np.array([-dpotx,-dpoty,-dpotz])

    return Dposition,Dmomenta,pot,dpotx,dpoty,dpotz,dmax

def gsang(coord,theta,phi,gamma,chargeList,dipole,lj,v,b,switch,fout,ion,parmList):
    '''
    /*
    * This function computes the scattering angle using the trajectory method
    * 
    * Arguments
    * ---------
    * coord : coordinates of the atoms within the molecule
    * theta, phi, gamma : angle of the rotation
    * chargeList : charges of each atom in the molecules
    * dipole : molecular dipole
    * lj : lennard-jones parameters
    * v : velocity of a gas particle
    * b : impact parameter
    * switch : switches for controlling verbosity of logfile
    * fout : file for logging
    * ion : information about the molecule
    * parmList : list of important constants
    * 
    * Returns
    * -------
    * ang : scattering angle
    * erat : energy ration (kinetic/potential)
    * (ifailc : keeps track of how many attempted trajectories are rejected (this is kind of a global variable) TODO: fix)
    * d1,isteps : not used. TODO: remove
    * 
    */    

    '''
    velocity=np.array([0.,-v,0.])
    d1=0.0
    istep=0
    if switch.it:
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
    if switch.it:
        print(' time steps, dt1 ={0: .4E} dt2 ={1: .4E}'.format(dt1,dt2),file=fout)
    
#     determine trajectory start position
    e0=0.5*ion.mu*v*v
    ymax=np.amax(coord[:,1])*1e10
    ymin=np.amin(coord[:,1])*1e10
    iymin=np.trunc(ymin)-1
    iymax=np.trunc(ymax)+1
    iy=iymax
    position=np.array([b,np.multiply(iy,1.0e-10),0])
    #print(iymin,iymax,end=" ")
    
    pot,dpotx,dpoty,dpotz,dmax=dljpot(position,coord,lj,chargeList,dipole)
    #print(pot,dipole,lj.eolj[0],position[1]);

    if np.abs(pot/e0)<=parmList.sw1:
        #print(iy,"first loop")
        for iy in np.arange(iymax-1,iymin,-1):
            position=np.array([b,np.multiply(iy,1.0e-10),0])
            pot,dpotx,dpoty,dpotz,dmax=dljpot(position,coord,lj,chargeList,dipole)
            if np.abs(pot/e0)>=parmList.sw1:
                    break
        else:
            if switch.it:
                print(' trajectory not started - potential too small',file=fout)
            ang=0.0
            erat=1.0
            return ang,erat,d1,istep

    else:
        #print(iy,"second loop")
        while np.abs(pot/e0)>parmList.sw1:
            iy=iy+10
            position=np.array([b,np.multiply(iy,1.0e-10),0])
            pot,dpotx,dpoty,dpotz,dmax=dljpot(position,coord,lj,chargeList,dipole)
        
        while np.abs(pot/e0)<parmList.sw1:
            iy=iy-1
            position=np.array([b,np.multiply(iy,1.0e-10),0])
            pot,dpotx,dpoty,dpotz,dmax=dljpot(position,coord,lj,chargeList,dipole)

    etot=e0+pot
    if switch.it:
        print(' trajectory start position ={0: .4E}'.format(iy),file=fout)
    d1=iy*1.0e-10

#     initial coordinates and momenta
    momenta=np.multiply(velocity,ion.mu)
    tim=0.0
    if switch.it:
        x=' '
        print('\n\n trajectory ns, x,  y,  z,  kin e, dt,    tot e',file=fout)
        print(16*x+'vx, vy, vz, pot e, pot/e0\n',file=fout)    
    
    #     initialize the time derivatives of the coordinates and momenta
    Dposition,Dmomenta,pot,dpotx,dpoty,dpotz,dmax=deriv(position,momenta,ion,coord,lj,chargeList,dipole)

    ns=0
    ang=0.
    erat=0.0
    # if switch.it:
    #     e=np.dot(momenta,momenta)*(0.5*ion.invMu)
    #     x=' '
    #     print(" {:5d} {: .4E} {: .4E} {: .4E} {: .4E} {: .4E} {: .4E}".format(ns,position[0],position[1],position[2],e,dt,pot+e),file=fout)
    #     print(7*x+"{: .4E} {: .4E} {: .4E} {: .4E} {: .4E}".format(Dposition[0],Dposition[1],Dposition[2],pot,np.abs(pot/e0)),file=fout)

    while True:
        while True:
            while True:
                for nw in np.arange(1,parmList.inwr+1):
                    position,momenta,Dposition,Dmomenta,dt,tim,pot,dpotx,dpoty,dpotz,dmax=diffeq(tim,dt,position,momenta,parmList,ion,coord,lj,chargeList,dipole,Dposition,Dmomenta)
                ns=ns+parmList.inwr

            #     print out the trajectory coordinates and velocities
                if switch.it:
                    e=np.dot(momenta,momenta)*(0.5*ion.invMu)
                    x=' '
                    print(" {:5d} {: .4E} {: .4E} {: .4E} {: .4E} {: .4E} {: .4E}".format(ns,position[0],position[1],position[2],e,dt,pot+e),file=fout)
                    print(7*x+"{: .4E} {: .4E} {: .4E} {: .4E} {: .4E}".format(Dposition[0],Dposition[1],Dposition[2],pot,np.abs(pot/e0)),file=fout)

            #     check if trajectory has become "lost" (too many steps)
                if ns>parmList.traj_lost:
                    print(' trajectory lost: b ={0: .4E} v={1: .4E}'.format(b,v),file=fout)
                    ang=parmList.pi/2.0
                    e=0.5*ion.mu*np.dot(Dposition,Dposition)
                    erat=(e+pot)/etot
                    istep=ns
                    return ang,erat,d1,istep

            #     check if the trajectory is finished
                if dmax>=lj.romax:
                    break
            if np.abs(pot/e0)>parmList.sw2 and dt==dt1:
                dt=dt2
                
            if np.abs(pot/e0)<parmList.sw2 and dt==dt2:
                dt=dt1
                
            if np.abs(pot/e0)<=parmList.sw1:
                break
        if ns>=50:
            break
    istep=ns

    ang=np.sign(Dposition[0])*np.arccos(-Dposition[1]/np.linalg.norm(Dposition))
    
    #     check for energy conservation
    e=0.5*ion.mu*np.dot(Dposition,Dposition)
    erat=(e+pot)/etot
    if erat<1.01 and erat>0.99:
        return ang,erat,d1,istep  #ang
    print('\n energy not conserved: e ratio ={0: .4E} v ={1: .4E} b ={2: .4E}'.format(erat,v,b),file=fout)
    print(' gst2 ={0: .4E} theta ={1: .4E} phi ={2: .4E} gamma ={3: .4E}'.format(0.5*ion.mu*v*v/parmList.eo,theta*parmList.cang,phi*parmList.cang,gamma*parmList.cang),file=fout)
    if switch.ip==1:
        print('',file=fout)
    parmList.ifailc=parmList.ifailc+1   ## TODO: make ifailc external to list
    if parmList.ifailc==parmList.ifail:
        fout.close()
        print('Energy not conserved')
        sys.exit(0)
    return ang,erat,d1,istep

def rantate(vector,parmList,switch,fout,RngGen):
    '''
    Rotates the vector(s) provided by a random orientation
    
    Arguments
    ---------
    vector : the vector(s) to  be rotated
    parmList : List with the constant used in the simulations
    switch : list with printing switches
    fout : file for logging
    RngGen : random number generator
    
    Returns
    -------
    newVector : the vector(s) after the rotation
    thate,phi,gamma : the angle of the (random) rotation

    '''

    rnt=RngGen.rand()
    rnp=RngGen.rand()
    rng=RngGen.rand()

    theta=rnt*2.0*parmList.pi
    phi=np.arcsin((rnp*2.0)-1.0)+(parmList.pi/2.0)
    gamma=rng*2.0*parmList.pi
    newVector=rotate(vector,theta,phi,gamma,parmList,switch,fout)
    return newVector,theta,phi,gamma


def mobil2(coordList,chargeList,fout,parmList,switch,lj,dipole,ion,mconst,iic,RngGen):
    '''
    Computes the cross section using the TM algorithm
    
    Arguments
    ---------
    coordList : coordinates of the atoms
    chargeList : list of charge of *each* atom
    fout : file for logging
    parmList : List with the constant used in the simulations
    switch : list with printing switches
    lj : structure with the Lennard-Jones parameters
    dipole : molecular dipole
    ion : structure with information about the *whole* molecule
    mconst : mobility constant
    iic : set of data
    RngGen : random number generator
    
    Returns
    -------
    mob : mobility
    cs : cross-section
    sdevpc : standard deviation

    '''
    if switch.im2==0:
        print('\n mobility calculation by MOBIL2 (trajectory method)\n',file=fout)
        print(' global trajectory parameters\n\n sw1 ={: .4E}       sw2 ={: .4E}'.format(parmList.sw1,parmList.sw2),file=fout)
        print(' dtsf1 ={: .4E}     dtsf2 ={: .4E}'.format(parmList.dtsf1,parmList.dtsf2),file=fout)
        print(' inwr ={: 3d}              ifail ={: 5d}'.format(parmList.inwr,parmList.ifail),file=fout)

    # I can set those with the global variables
    switch.SetIt(0)
    switch.SetIu2(0)

    #     determine maximum extent and orient along x axis
    
    if switch.im2==0:
        print('\n maximum extent orientated along x axis',file=fout)

    rmax=0
    ihold=0
    for c,ivec in enumerate(coordList):
        tmp=np.linalg.norm(ivec)
        if tmp>rmax:
            rmax=tmp
            vector=ivec
            ihold=c
        
    rzy=np.sqrt(vector[2]**2+vector[1]**2)
    phi=np.arccos(vector[2]/rzy)
    phi=phi+(0.5*parmList.pi)  # I want to make the z component =0, thus shift to y-axis
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
        gamma=(2.*parmList.pi)-gamma
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
        print('\n Problem orientating along x axis\n',file=fout)
        for ivec in newVector:
            print(ivec[0],ivec[1],ivec[2],np.linalg.norm(ivec),file=fout)

    
    #     determine rmax, emax, and r00 along x, y, and z directions    
    if switch.ip==1:
        print("\n",file=fout)
    irn=1000
    ddd=(rmax+lj.romax)/float(irn)

    y=0.0
    z=0.0
    emaxx=0.0
    for ir in range(irn):
        x=rmax+lj.romax-(float(ir)*ddd)
        pot,dpotx,dpoty,dpotz,dmax=dljpot(np.array([x,y,z]),newVector,lj,chargeList,dipole)
        if pot>0:
            break
        r00x=x
        if  pot<emaxx:
            rmaxx=x
            emaxx=pot
    if switch.im2==0:
        print(' along x axis emax ={: .4E}eV rmax ={: .4E}A r00 ={: .4E}A'.format(emaxx/parmList.xe,rmaxx*1.0e10,r00x*1.0e10),file=fout)
    
    x=0.0
    z=0.0
    emaxy=0.0
    for ir in range(irn):
        y=rmax+lj.romax-(float(ir)*ddd)
        pot,dpotx,dpoty,dpotz,dmax=dljpot(np.array([x,y,z]),newVector,lj,chargeList,dipole)
        if pot>0:
            break
        r00y=y
        if  pot<emaxy:
            rmaxy=y
            emaxy=pot
    if switch.im2==0:
        print(' along y axis emax ={: .4E}eV rmax ={: .4E}A r00 ={: .4E}A'.format(emaxy/parmList.xe,rmaxy*1.0e10,r00y*1.0e10),file=fout)

    x=0.0
    y=0.0
    emaxz=0.0
    for ir in range(irn):
        z=rmax+lj.romax-(float(ir)*ddd)
        pot,dpotx,dpoty,dpotz,dmax=dljpot(np.array([x,y,z]),newVector,lj,chargeList,dipole)
        if pot>0:
            break
        r00z=z
        if  pot<emaxz:
            rmaxz=z
            emaxz=pot
    if switch.im2==0:
        print(' along z axis emax ={: .4E}eV rmax ={: .4E}A r00 ={: .4E}A\n'.format(emaxz/parmList.xe,rmaxz*1.0e10,r00z*1.0e10),file=fout)

#     set-up integration over gst
    tst=parmList.xk*ion.temperature/parmList.eo
    if switch.im2==0:
        print('\n t*={: .4E}'.format(tst),file=fout)
    tst3=tst*tst*tst
    dgst=5e-7*6*np.sqrt(tst)
    gst=dgst
    sum=0.0
    
    sequenceNumber=np.sqrt(np.arange(1,parmList.inp+1))
    sum1=np.sum(sequenceNumber)
    sum2=np.cumsum(np.sqrt(np.arange(0,parmList.inp)))
    if switch.im2==0:
        x=' '
        print('\n\n set-up gst integration - integration over velocity',file=fout)
        print('\n     pgst'+8*x+'wgst'+9*x+'v'+9*x+'ke/kt'+7*x+'gst^5*'+5*x+'frac of\n'+48*x+'exp(gst^2/tst)'+3*x+'sum\n',file=fout)
    
    wgst=sequenceNumber/sum1
    pgst=np.zeros(parmList.inp+1)
    gstt=tst3*(sum2+(0.5*sequenceNumber))/sum1

    for i in np.arange(1,parmList.inp+1):
        while True:
            sum=sum+(np.exp(-gst*gst/tst)*gst*gst*gst*gst*gst*dgst)
            gst=gst+dgst
            if sum>gstt[i-1]:
                pgst[i-1]=gst-(0.5*dgst)
                break
            if sum==gstt[i-1]:
                break

        if switch.im2==0:
            hold1=np.sqrt((pgst[i-1]*pgst[i-1]*parmList.eo)/(0.5*ion.mu))
            hold2=0.5*ion.mu*hold1*hold1/(parmList.xk*ion.temperature)
            hold3=np.exp(-pgst[i-1]*pgst[i-1]/tst)*pgst[i-1]**5
            print(' {: .4E} {: .4E} {: .4E} {: .4E} {: .4E} {: .4E}'.format(pgst[i-1],wgst[i-1],hold1,hold2,hold3,sum/tst3),file=fout)            

#     determine b2max
    b2max=np.zeros(parmList.inp)
    cosx=np.zeros(parmList.ibst_max)
    dbst2=1.0
    dbst22=0.1*dbst2
    cmin=0.0005
    if switch.im2==0:
        print('\n\n set up b2 integration - integration over impact parameter\n\n minimum value of (1-cosX) ={: .4E}\n'.format(cmin),file=fout)
    gst2=np.square(pgst)
    vList=np.sqrt((gst2*parmList.eo)/(0.5*ion.mu))
    for ig in range(parmList.inp-1,-1,-1):
        ibst=np.trunc(rmaxx/parmList.ro)-6
        if ig<parmList.inp-1:
            ibst=int(np.trunc(b2max[ig+1]/dbst2))-6
        if ibst<0:
            ibst=0
        if switch.ip==1:
            x=' '
            print('\n gst2 ={: .4E} v ={: .4E}'.format(gst2[ig],vList[ig])+'\n'+6*x+'b'+10*x+'bst2'+7*x+'X ang'+7*x+'cos(X)'+6*x+'e ratio',file=fout)
        while True:        
            bst2=dbst2*float(ibst)
            b=parmList.ro*np.sqrt(bst2)
            
            ang,erat,d1,istep=gsang(newVector,theta,phi,gamma,chargeList,dipole,lj,vList[ig],b,switch,fout,ion,parmList)
            #print (ibst,b,ang,np.cos(ang))
            cosx[ibst]=1.0-np.cos(ang)
            if(switch.ip==1):
                print(' {: .4E} {: .4E} {: .4E} {: .4E} {: .4E}'.format(b,bst2,ang,cosx[ibst],erat),file=fout)

            if ibst<5 or ( ibst>=5 and (not (cosx[ibst]<cmin and cosx[ibst-1]<cmin and cosx[ibst-2] < cmin and cosx[ibst-3]<cmin and cosx[ibst-4]<cmin))):
                ibst=ibst+1
                if ibst>=parmList.ibst_max:
                    print(' ibst greater than ',parmList.ibst_max)
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
        x=' '
        print('\n'+5*x+'gst'+11*x+'b2max/ro2'+9*x+'b/A\n',file=fout)
        for ig in range(parmList.inp):
            print(' {: .4E}     {: .4E}     {: .4E}'.format(pgst[ig],b2max[ig],parmList.ro*np.sqrt(b2max[ig])*1.0e10),file=fout)
            

#     Calculate Omega(1,1)*, Omega(1,2)*, Omega(1,3)*, and Omega(2,2)*
#     by integrating Q(1)* or Q(2)* over all orientations, and initial 
#     relative velocities.
    if(switch.im2==0):
        print('\n\n number of complete cycles (itn) ={:6d}'.format(parmList.itn),file=fout)
        print(' number of velocity points (inp) ={:6d}'.format(parmList.inp),file=fout)
        print(' number of random points (imp) ={:6d}'.format(parmList.imp),file=fout)
        print(' total number of points ={:7d}\n'.format(parmList.itn*parmList.inp*parmList.imp),file=fout)
    if switch.ip==1:
        print(' start mobility calculation',file=fout)
    
    q1st=np.zeros(parmList.inp)
    q2st=np.zeros(parmList.inp)
    om11st=np.zeros(parmList.itn)
    om12st=np.zeros(parmList.itn)
    om13st=np.zeros(parmList.itn)
    om22st=np.zeros(parmList.itn)

    for ic in range(parmList.itn):
        if switch.ip==1:
            print('\n cycle number, ic ={:3d}'.format(ic+1),file=fout)
        # DIFF : next 2 lines unecessary in original code? The are only useful if switch.ip=1. So I added those there
        #   gst2=pgst(ig)*pgst(ig)
        #   v=dsqrt((gst2*eo)/(0.5d0*mu))
        for ig in range(parmList.inp):
            if switch.ip==1:
                # gst2=pgst[ig]*pgst[ig]
                # v=np.sqrt((gst2*parmList.eo)/(0.5*ion.mu))
                print('\n ic ={:3d} ig ={:4d} gst2 ={: .4E} v ={: .4E}'.format(ic+1,ig+1,gst2[ig],vList[ig]),file=fout)
            temp1=0.0
            temp2=0.0
            for im in range(parmList.imp):
                if switch.ip==1 and im==0:
                    x=' '
                    print('\n'+5*x+'b/A'+8*x+'ang'+6*x+'(1-cosX)'+4*x+'e ratio'+4*x+'theta'+7*x+'phi'+7*x+'gamma',file=fout)
                rnb=RngGen.rand()
                newVector,theta,phi,gamma=rantate(coordList,parmList,switch,fout,RngGen)
                bst2=rnb*b2max[ig]
                b=parmList.ro*np.sqrt(bst2)
                if switch.igs==1:
                    ftemp=open('hold_py',"w")
                    print(iic,ic+1,ig+1,im+1,file=ftemp)
                    print(vList[ig],b,file=ftemp)
                    print(theta*parmList.cang,phi*parmList.cang,gamma*parmList.cang,file=ftemp)
                    ftemp.close()
                ang,erat,d1,istep=gsang(newVector,theta,phi,gamma,chargeList,dipole,lj,vList[ig],b,switch,fout,ion,parmList)
                #print("ang:",ang)
                hold1=1.0-np.cos(ang)
                hold2=np.sin(ang)*np.sin(ang)
                if switch.ip==1:
                    print(' {: .4E}{: .4E}{: .4E}{: .4E}{: .4E}{: .4E}{: .4E}'.format(b*1.e10,ang*parmList.cang,hold1,erat,theta*parmList.cang,phi*parmList.cang,gamma*parmList.cang),file=fout)
                temp1=temp1+(hold1*b2max[ig]/float(parmList.imp))
                temp2=temp2+(1.5*hold2*b2max[ig]/float(parmList.imp))
            om11st[ic]=om11st[ic]+(temp1*wgst[ig])
            om12st[ic]=om12st[ic]+(temp1*pgst[ig]*pgst[ig]*wgst[ig]*(1.0/(3.0*tst)))
            om13st[ic]=om13st[ic]+(temp1*(pgst[ig]**4)*wgst[ig]*(1.0/(12.0*tst*tst)))
            om22st[ic]=om22st[ic]+(temp2*pgst[ig]*pgst[ig]*wgst[ig]*(1.0/(3.0*tst)))
            q1st[ig]=q1st[ig]+temp1
            q2st[ig]=q2st[ig]+temp2
            if switch.ip==1:
                print('\n v ={: .4E}     q1st ={: .4E}\n'.format(vList[ig],q1st[ig]),file=fout)
        if switch.ip==1:
            print('\n OMEGA(1,1)*={: .4E}\n'.format(om11st[ic]),file=fout)

    #     calculate running averages
    temp=1.0/(mconst/(np.sqrt(ion.temperature)*om11st*parmList.pi*parmList.ro2))
    hold1=np.cumsum(om11st)
    hold2=np.cumsum(temp)

    if switch.im2==0:
        x=' '
        print('\n summary of mobility calculations\n\n cycle'+5*x+'cs/A^2'+6*x+'avge cs/A^2'+8*x+'Ko^-1'+7*x+'avge Ko^-1',file=fout)
        for icc in range(1,parmList.itn+1):
            print(' {: 3d}    {: .4E}    {: .4E}    {: .4E}    {: .4E}'.format(icc,om11st[icc-1]*parmList.pi*parmList.ro2*1.e20,hold1[icc-1]*parmList.pi*parmList.ro2*1e20/float(icc),temp[icc-1],hold2[icc-1]/float(icc)),file=fout)
        print('\n\n average values for q1st\n\n'+5*x+'gst2'+8*x+'wgst'+8*x+'q1st',file=fout)
        for ig in range(parmList.inp):
            print(' {: .4E} {: .4E} {: .4E}'.format(pgst[ig]*pgst[ig],wgst[ig],q1st[ig]/float(parmList.inp)),file=fout)
    
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
        print('\n\n mean OMEGA*(1,1) ={: .4E}\n standard deviation ={: .4E}\n standard error of mean ={: .4E}'.format(mom11st,sdom11st,sterr),file=fout)
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
        x=' '
        print('\n\n f value for second order correction={: .4E}'.format(f),file=fout)
        print(' (integrations for second order correction are not',file=fout)
        print(' accurate, check if correction becomes significant)',file=fout)
        print('\n omega*12 ={: .4E}  omega*13 ={: .4E}  omega*22 ={: .4E}'.format(mom12st,mom13st,mom22st),file=fout)
        print(7*x+'u2 ={: .4E}         w ={: .4E}     delta ={: .4E}'.format(u2,w,delta),file=fout)
    mob=(mconst*f)/(np.sqrt(ion.temperature)*cs)
    print('\n\n average (second order) TM mobility ={: .4E}'.format(mob),file=fout)
    print(' inverse average (second order) TM mobility ={: .4E}'.format(1.0/mob),file=fout)
    print(' average TM cross section ={: .4E}'.format(cs*1.e20),file=fout)

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
            item[1]=item[1]-yr
            item[2]=item[2]-zr
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
        # NOTE: 997 format(15x,1pe12.5,1x,e12.5,1x,e12.5) --Unnecessary?
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

def mobil4(coordList,chargeList,fout,parmList,switch,lj,ion,mconst,RngGen):
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
    imm=2
    #newVector=coordList

# Monte Carlo starts here

    for i in range(parmList.inum):
        newVector,theta,phi,gamma=rantate(coordList,parmList,switch,fout,RngGen)
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
        yr=ymin+ydi*RngGen.rand()
        zr=zmin+zdi*RngGen.rand()
        # Write the direction cosines of the initial incidence vector 
        # (collinear to the x axis).
        # DIFF: The original code added the versor as an extra coordinate.

        versor=np.array([1.,0.,0.])
        kp=0
        cof=np.zeros(parmList.inor+2)
        for im in range(2,parmList.inor+2):
            cof[im],newVector,versor,yr,zr,kp=che(im-1,newVector,lj,cof[im-1],versor,parmList,yr,zr,kp)
        # c     If the inclusion of the next order of collision had not changed
        # c     incidence/reflection angle (i.e., no collision occurred) stop 
        # c     following this trajectory.
            if cof[im]==cof[im-1]:
                if im-1 >imm:
                    imm=im-1
                    #     Set all further incidence/reflection angle cosines along this
                    #     trajectory to the last angle found.

                for imn in range (im+1,parmList.inor+2):
                    cof[imn]=cof[im]
                break
                
            else:
                if im> imm:
                    imm=im  
    #     Add the contributions from the i-th trajectory to the collision 
    #     cross sections of all orders.
        for im in range(2,parmList.inor+2):
            crof[im-2]=crof[im-2]+pls*cof[im]*cof[im]
        if kp==1:
            crb=crb+pls
    # END of Monte Carlo
    u2=0.5*float(parmList.inum)
    #Normalize all collision cross-sections and the projection.
    crof=crof/u2
    crb=crb/float(parmList.inum)
    pacs=crb
    pamob=mconst/(np.sqrt(ion.temperature)*pacs)
    ehscs=crof[imm-2]
    #print(imm,crof)
    ehsmob=mconst/(np.sqrt(ion.temperature)*ehscs)

    #Output the results of calculation.
    if switch.im4==0:
        print('\n\n mobility calculation by MOBIL4 (HS scattering)',file=fout)
        print('\n number of Monte Carlo trajectories ={:7d}'.format(parmList.inum),file=fout)
        print(' maximum number of reflections allowed ={:3d}'.format(parmList.inor),file=fout)
    print('\n average PA mobility ={: .4E}'.format(pamob),file=fout)
    print(' inverse average PA mobility ={: .4E}'.format(1.0/pamob),file=fout)
    print(' average PA cross section ={: .4E}'.format(pacs*1.e20),file=fout)
    print('\n maximum number of reflections encountered ={:3d}'.format(imm-1),file=fout)
    
    if switch.im4==0:
        print('\n     order   cross section',file=fout)
        for im in range(1,imm):
            print('     {:3d}     {: .5E}'.format(im,crof[im-1]),file=fout)

    print('\n average EHS mobility ={: .4E}'.format(ehsmob),file=fout)
    print(' inverse average EHS mobility ={: .4E}'.format(1.0/ehsmob),file=fout)
    print(' average EHS cross section ={: .4E}'.format(ehscs*1.e20),file=fout)
    
    return ehscs,ehsmob,pacs,pamob,imm-1

        
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
        Rseed=int(x[0])

# inFileName="a10A1.mfj"
# outFileName="a10A1.out"
# Rseed=None

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

# LJ scaling parameters
atoms,ion,InputData,tmp_asymp,lj,hs,fin,fout=fcoord(inFileName,outFileName,ConstList,AtomList,printDebug)

print(' Lennard-Jones scaling parameters: eo={: .4E} ro={: .4E}'.format(ConstList.eo/ConstList.xe,ConstList.ro*1e+10),file=fout)


    #  Constant for ion-induced dipole potential

    #  Polarizability of helium = 0.204956d-30 m3
    #  xeo is permitivity of vacuum, 8.854187817d-12 F.m-1

dipole=ConstList.xe**2*0.204956e-30/(8*ConstList.pi*ConstList.xeo)
ion.SetDipole(dipole)
print(' dipole constant ={: .4E}'.format(ion.dipole),file=fout)

    # Mass constants
ion.SetMassGas(4.0026,ConstList.xn*1000)

    # Mobility constant

mconst=np.sqrt(18*ConstList.pi)/16
mconst=mconst*np.sqrt(ConstList.xn*1000)*np.sqrt((1/ion.mGas)+(1/ion.TotalMass))
mconst=mconst*ConstList.xe/np.sqrt(ConstList.xk)
dens=ConstList.xn/ConstList.xmv
mconst=mconst/dens
ion.SetMobility(mconst)
print(' mobility constant ={: .4E}'.format(ion.mobility),file=fout)
print(' temperature ={: .4E}'.format(ion.temperature),file=fout)


#DIFF: In the original code follows an initialization of the random number generator.
# Here I am using numpy, so that section is unnecessary

RngGen=np.random.RandomState(seed=Rseed)

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
immmin=ConstList.inor
immmax=0
###

for iic in range(InputData.icoord):
    print('\n',file=fout)
    if InputData.icoord > 1:
        print('\n coordinate set ={:5d}'.format(iic+1),file=fout)
    print('\n structural asymmetry parameter ={: 8.4f}'.format(asympp[iic]),file=fout)

    ehsc[iic],ehsm[iic],pac[iic],pam[iic],imm=mobil4(atoms.Coords,atoms.Charges,fout,ConstList,printDebug,hs,ion,mconst,RngGen) 
    tmm[iic],tmc[iic],sdevpc=mobil2(atoms.Coords,atoms.Charges,fout,ConstList,printDebug,lj,ion.dipole,ion,mconst,iic,RngGen)

    if(imm<immmin):
        immmin=imm
    if(imm>immmax):
        immmax=imm  
   
    if iic<InputData.icoord-1:
        Oxyz,charges,asympp[iic+1],fin,fout=ncoord(fin,fout,ConstList,AtomList,printDebug,ion,InputData)
        atoms.UpdatePosition(Oxyz)
        atoms.UpdateCharges(charges)
    printDebug.im2=1
    printDebug.im4=1

if not fin.closed:
    fin.close()

#     print out summary
print('\n\n\n SUMMARY\n',file=fout)
print(' program version = 0.3',file=fout)
print(' input file name = {:30s}'.format(inFileName),file=fout)
print(' input file label = {:30s}'.format(InputData.label),file=fout)
if tmmob!=0:
    if InputData.charges=='equal':
        fout.write(" using a uniform charge distribution\n")
    elif InputData.charges=='calc':
        fout.write(" using a calculated (non-uniform) charge distribution\n")
    elif InputData.charges=='none':
        fout.write(" using no charge - only LJ interactions\n")
print(' temperature ={: .4E}'.format(ion.temperature),file=fout)
print(' using NUMPY random number generator with seed  =',Rseed,file=fout)

if InputData.icoord==1:
    print(' structural asymmetry parameter ={: .4E}'.format(asympp[0]),file=fout)
    if ehsm[0]!=0:
        print('\n mobility calculation by MOBIL4 (HS scattering)',file=fout)
        print('\n number of Monte Carlo trajectories ={:7d}'.format(ConstList.inum),file=fout)
        print(' maximum number of reflections encountered ={:3d}'.format(imm),file=fout)
        print('\n inverse average PA mobility ={: .4E}'.format(1.0/pam[0]),file=fout)
        print(' average PA cross section ={: .4E}'.format(pac[0]*1.e20),file=fout)
        print('\n inverse average EHS mobility ={: .4E}'.format(1.0/ehsm[0]),file=fout)
        print(' average EHS cross section ={: .4E}'.format(ehsc[0]*1e20),file=fout)
    
    if tmm[0]!=0:
        print('\n mobility calculation by MOBIL2 (trajectory method)',file=fout)
        print('\n trajectory parameters',file=fout)
        print(' sw1 ={: .4E}       sw2 ={: .4E}'.format(ConstList.sw1,ConstList.sw2),file=fout)
        print(' dtsf1 ={: .4E}     dtsf2 ={: .4E}'.format(ConstList.dtsf1,ConstList.dtsf2),file=fout)
        print('\n number of complete cycles (itn) ={:6d}'.format(ConstList.itn),file=fout)
        print(' number of velocity points (inp) ={:6d}'.format(ConstList.inp),file=fout)
        print(' number of random points (imp) ={:6d}'.format(ConstList.imp),file=fout)
        print(' total number of points ={:7d}'.format(ConstList.itn*ConstList.inp*ConstList.imp),file=fout)
        print('\n inverse average (second order) TM mobility ={: .4E}'.format(1.0/tmm[0]),file=fout)
        print(' average TM cross section ={: .4E}'.format(tmc[0]*1.e20),file=fout)
        print(' standard deviation (percent) ={: .4E}'.format(sdevpc),file=fout)
        print(' number of failed trajectories ={:4d}'.format(ConstList.ifailc),file=fout)

#     multiple coordinate sets
else:
    if ehsm[0]!=0:
        x=' '
        print('\n mobility calculation by MOBIL4 (HS scattering)',file=fout)
        print('\n number of Monte Carlo trajectories ={:7d}'.format(ConstList.inum),file=fout)
        print('\n minimum and maximum number of reflections ={:3d}  {:3d}'.format(immmin,immmax),file=fout)
        print('\n\n'+4*x+'set'+5*x+'PA CS'+7*x+'PA MOB^-1'+6*x+'EHSS CS'+6*x+'EHSS MOB^-1'+4*x+'ASYMP',file=fout)
        #DIFF: in the following printing the last element asympp is divided by 10 in the fortran code, but not here.
        # That is only a formatting hack, not needed here.
        for i in range(InputData.icoord):
            print(' {:5d}   {: .4E}   {: .4E}   {: .4E}   {: .4E}   {: 8.4f}'.format(i+1,pac[i]*1.0e20,1.0/pam[i],ehsc[i]*1.0e20,1.0/ehsm[i],asympp[i]),file=fout)
        pacs=np.sum(pac)/float(InputData.icoord)
        pamob=np.sum(pam)/float(InputData.icoord)
        ehscs=np.sum(ehsc)/float(InputData.icoord)
        ehsmob=np.sum(ehsm)/float(InputData.icoord)
        aasymp=np.sum(asympp)/float(InputData.icoord)
        #DIFF: in the following printing the last element aasymp is divided by 10 in the fortran code, but not here.
        # That is only a formatting hack, not needed here.

        print('\n   AVGE  {: .4E}   {: .4E}   {: .4E}   {: .4E}   {: 8.4f}'.format(pacs*1.0e20,1.0/pamob,ehscs*1.0e20,1.0/ehsmob,aasymp),file=fout)
    if tmm[0]!=0:
        x=' '
        print('\n mobility calculation by MOBIL2 (trajectory method)',file=fout)
        print('\n trajectory parameters',file=fout)
        print(' sw1 ={: .4E}       sw2 ={: .4E}'.format(ConstList.sw1,ConstList.sw2),file=fout)
        print(' dtsf1 ={: .4E}     dtsf2 ={: .4E}'.format(ConstList.dtsf1,ConstList.dtsf2),file=fout)
        print('\n number of complete cycles (itn) ={:6d}'.format(ConstList.itn),file=fout)
        print(' number of velocity points (inp) ={:6d}'.format(ConstList.inp),file=fout)
        print(' number of random points (imp) ={:6d}'.format(ConstList.imp),file=fout)
        print(' total number of points ={:7d}'.format(ConstList.itn*ConstList.inp*ConstList.imp),file=fout)
        print('\n total number of failed trajectories ={:4d}'.format(ConstList.ifailc),file=fout)
        print('\n\n'+4*x+'set'+5*x+'TM CS'+7*x+'TM MOB^-1',file=fout)
        for i in range(InputData.icoord):
            print(' {:5d}   {: .4E}   {: .4E}'.format(i+1,tmc[i]*1.0e20,1.0/tmm[i]),file=fout)
        tmcs=np.sum(tmc)/float(InputData.icoord)
        tmmob=np.sum(tmm)/float(InputData.icoord)
        print('\n   AVGE  {: .4E}   {: .4E}'.format(tmcs*1.0e20,1.0/tmmob),file=fout)

fout.close()
