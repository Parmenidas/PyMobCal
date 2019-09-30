'''
Test of mobil 4
'''

import math
import numpy as np
import scipy
from scipy.spatial.transform import Rotation as R
import sys
import mobil4_f
    
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


def rantate(vector,parmList,switch,fout,RngGen):
    '''
c     Rotates the cluster/molecule to a random orientation.  
    NOTE: here I have used rando_py generated by the script.
    '''

    rnt=RngGen.rand()
    rnp=RngGen.rand()
    rng=RngGen.rand()

    theta=rnt*2.0*parmList.pi
    phi=np.arcsin((rnp*2.0)-1.0)+(parmList.pi/2.0)
    gamma=rng*2.0*parmList.pi
    newVector=rotate(vector,theta,phi,gamma,parmList,switch,fout)
    return newVector,theta,phi,gamma

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
        # DIFF: The original code added this vector as an extra coordinate.

        versor=[1,0,0]
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

def MapCoords(f_pack,coords):
    for c,i in enumerate(coords):
        f_pack.coordinates.ox[c]=i[0]
        f_pack.coordinates.oy[c]=i[1]
        f_pack.coordinates.oz[c]=i[2]

def MapCoordsF(f_pack,coords):
    for c,i in enumerate(coords):
        f_pack.coordinates.fx[c]=i[0]
        f_pack.coordinates.fy[c]=i[1]
        f_pack.coordinates.fz[c]=i[2]

def MapCharges(f_pack,charges):
    for c,i in enumerate(charges):
        f_pack.charge.pcharge[c]=i

def MapHardSphere(f_pack,hs):
    for c,i in enumerate(hs.rhs):
        f_pack.hsparameters.rhs[c]=hs.rhs[c]
        f_pack.hsparameters.rhs2[c]=hs.rhs2[c]

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

#tmm=np.zeros(InputData.icoord)
#tmc=np.zeros(InputData.icoord)
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
###########################
## Set Python parameters ##
###########################
InputData.icoord=1
ConstList.ipr=2
ConstList.inp=2
ConstList.itn=2
ConstList.imp=2
fout_t_py=open("out.python.dat","w") # file with the result for python
ConstList.inum=10 
ConstList.inor=10
###########################
## Set Fortran parameters ##
###########################
rando_f=np.random.RandomState(seed=43)
mobil4_f.xrand=rando_f.rand
MapCommonsConstants(mobil4_f,printDebug,ConstList,ion.mu,ion.TotalMass,dipole,ion.mobility,InputData.cfact,InputData.atoms,InputData.icoord,1,lj)
MapLJConstants(mobil4_f,lj)
MapCoords(mobil4_f,atoms.Coords)
MapCharges(mobil4_f,atoms.Charges)
MapHardSphere(mobil4_f,hs)
#tmm_f=np.zeros(InputData.icoord)
#tmc_f=np.zeros(InputData.icoord)
ehsc_f=np.zeros(InputData.icoord)
ehsm_f=np.zeros(InputData.icoord)
pac_f=np.zeros(InputData.icoord)
pam_f=np.zeros(InputData.icoord)


fout_py=open("python_file.dat","w") # file with the result for python
fout_f=open("fortran_file.dat","w") # file with the result for fortran
fout_diff=open("diff_file.dat","w") # file with the difference of the two files above

for iic in range(InputData.icoord):
    print("",file=fout)
    if InputData.icoord > 0:
        print('\n coordinate set ={:5d}'.format(iic+1),file=fout)
    print('\n structural asymmetry parameter ={: .4E}'.format(asympp[iic]),file=fout)

    mobil4_f.constants.iic=iic
    ehsc[iic],ehsm[iic],pac[iic],pam[iic],imm=mobil4(atoms.Coords,atoms.Charges,fout_t_py,ConstList,printDebug,hs,ion,mconst,RngGen)
    #mobil4(atoms.Coords,atoms.Charges,fout_t_py,ConstList,printDebug,hs,ion,mconst,RngGen)
    ehsc_f[iic],ehsm_f[iic],pac_f[iic],pam_f[iic],imm_f=mobil4_f.mobil4(ion.temperature,ConstList.inum,ConstList.inor)
    #mobil4_f.mobil4(ion.temperature,ConstList.inum,ConstList.inor)

    print(ehsc[iic],ehsm[iic],pac[iic],pam[iic],imm,file=fout_py)
    print(ehsc_f[iic],ehsm_f[iic],pac_f[iic],pam_f[iic],imm_f,file=fout_f)
    print(ehsc_f[iic]-ehsc[iic],ehsm_f[iic]-ehsm[iic],pac_f[iic]-pac[iic],pam_f[iic]-pam[iic],imm_f-imm,file=fout_diff)

    if(imm<immmin):
        immmin=imm
    if(imm>immmax):
        immmax=imm  
   
    if iic<InputData.icoord-1:
        Oxyz,asympp[iic+1],fin,fout=ncoord(fin,fout,ConstList,AtomList,printDebug,ion,InputData)
        atoms.UpdatePosition(Oxyz)
        MapCoords(mobil4_f,atoms.Coords)
    printDebug.im2=1
    printDebug.im4=1

if not fin.closed:
    fin.close()

fout.close()
fout_t_py.close()
fout_py.close()
fout_f.close()
fout_diff.close()
