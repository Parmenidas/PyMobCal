'''
Test of function gsang. the new version of gsang use velocity
verlet algorithm for the integration of the equations of motions

note that the fortran version of gsang has been changed adding 
open/close the file and a check for ns>30000 has been replaced with 
a check for ns>1000 to speed up the test. The function diffeq has an
extra common /testing/ for local variables.

**IMPORTANT**: These function create their own output files that
should be checked for differences (out.fortrant.dat and out.pyhton.dat).
This is in *addition* to diff_file.dat

The last test appears to show that the velocity verlet better conserves the energy
of the system, while the original algorithm fails to conserve energy and stops 
the calculation.

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
ion.invMu=1./mu

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
