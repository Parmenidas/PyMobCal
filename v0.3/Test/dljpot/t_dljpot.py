'''
Here I test a new version of dljpot. 
This version is completely vectorized.
'''

import dljpot_f
import numpy as np

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
        self.vec_dim=100            #DIFF: this is what is used within the FORTRAN function as dimension of arrays 




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
pot,dpotx,dpoty,dpotz,dmax=0.,0.,0.,0.,10 # NOTE: they are output 
eolj=1.
eox4=4*eolj
rolj=1.
ro2lj=np.square(rolj)
ro6lj=np.power(ro2lj,3)
ro12lj=np.square(ro6lj)
dro6=6*ro6lj
dro12=12*ro12lj
romax=np.amax(rolj)

ConstList=constants()
#dipole=ConstList.xe**2*0.204956e-30/(8*ConstList.pi*ConstList.xeo)
dipole=0.

###############################################################
# Set Commons (FORTRAN). Not all variables in commons are actually used. 
# We only set the ones used.
###############################################################
dljpot_f.constants.inatom=1
dljpot_f.constants.romax=romax
dljpot_f.ljparameters.eox4[0]=eox4
dljpot_f.ljparameters.ro6lj[0]=ro6lj
dljpot_f.ljparameters.ro12lj[0]=ro12lj
dljpot_f.ljparameters.dro6[0]=dro6
dljpot_f.ljparameters.dro12[0]=dro12
#######################################
## Python code
#######################################
lj=LennardJones(np.array([eolj]),np.array([rolj]))

##########################
## Move Particle along x
##########################

fout_f=open("fortran_file.dat","w")
fout_py=open("python_file.dat","w")
fout_diff=open("diff_file.dat","w")
for dipole,tag_d in zip((0,1),("","d_")):
    for pcharge,tag_c in zip((0.,1.),("","q_")):
        fx,fy,fz=0,0,0
        tag=tag_d+tag_c
        header=tag+"x"
        print(header+"\nPos \t Pot \t dpotx \t dpoty \t dpotz \t dmax",file=fout_f)
        print(header+"\nPos \t Pot \t dpotx \t dpoty \t dpotz \t dmax",file=fout_py)
        print(header+"\nPos \t Pot \t dpotx \t dpoty \t dpotz \t dmax",file=fout_diff)

        for fx in np.arange(0.9,3,0.1):
            ###############################################################
            # Set Commons (FORTRAN). Not all variables in commons are actually used. 
            # We only set the ones used.
            ###############################################################
            dljpot_f.constants.dipol=dipole
            dljpot_f.coordinates.fx[0]=fx
            dljpot_f.coordinates.fy[0]=fy
            dljpot_f.coordinates.fz[0]=fz
            dljpot_f.charge.pcharge[0]=pcharge
            pot,dpotx,dpoty,dpotz,dmax=dljpot_f.dljpot(x,y,z)
            print (round(fx,1),"\t",pot,"\t",dpotx,"\t",dpoty,"\t",dpotz,"\t",dmax,file=fout_f)

            #######################################
            ## Python code
            #######################################
            vector=np.array([[fx,fy,fz]])
            charge=np.array([pcharge])
            dipol=dipole
            pot2,dpotx2,dpoty2,dpotz2,dmax2=dljpot(np.array([x,y,z]),vector,lj,charge,dipol)
            print (round(fx,1),"\t",pot2,"\t",dpotx2,"\t",dpoty2,"\t",dpotz2,"\t",dmax2,file=fout_py)
            # Difference   
            print (round(fx,1),"\t",pot2-pot,"\t",dpotx2-dpotx,"\t",dpoty2-dpoty,"\t",dpotz2-dpotz,"\t",dmax2-dmax,file=fout_diff)

        header=tag+"y"
        print("\n\n"+header+"\nPos \t Pot \t dpotx \t dpoty \t dpotz \t dmax",file=fout_f)
        print("\n\n"+header+"\nPos \t Pot \t dpotx \t dpoty \t dpotz \t dmax",file=fout_py)
        print("\n\n"+header+"\nPos \t Pot \t dpotx \t dpoty \t dpotz \t dmax",file=fout_diff)

        fx=0
        for fy in np.arange(0.9,3,0.1):
            ###############################################################
            # Set Commons (FORTRAN). Not all variables in commons are actually used. 
            # We only set the ones used.
            ###############################################################
            dljpot_f.constants.dipol=dipole
            dljpot_f.coordinates.fx[0]=fx
            dljpot_f.coordinates.fy[0]=fy
            dljpot_f.coordinates.fz[0]=fz
            dljpot_f.charge.pcharge[0]=pcharge
            pot,dpotx,dpoty,dpotz,dmax=dljpot_f.dljpot(x,y,z)
            print (round(fy,1),"\t",pot,"\t",dpotx,"\t",dpoty,"\t",dpotz,"\t",dmax,file=fout_f)

            #######################################
            ## Python code
            #######################################
            vector=np.array([[fx,fy,fz]])
            charge=np.array([pcharge])
            dipol=dipole
            pot2,dpotx2,dpoty2,dpotz2,dmax2=dljpot(np.array([x,y,z]),vector,lj,charge,dipol)
            print (round(fy,1),"\t",pot2,"\t",dpotx2,"\t",dpoty2,"\t",dpotz2,"\t",dmax2,file=fout_py)
            # Difference   
            print (round(fy,1),"\t",pot2-pot,"\t",dpotx2-dpotx,"\t",dpoty2-dpoty,"\t",dpotz2-dpotz,"\t",dmax2-dmax,file=fout_diff)

        header=tag+"z"
        print("\n\n"+header+"\nPos \t Pot \t dpotx \t dpoty \t dpotz \t dmax",file=fout_f)
        print("\n\n"+header+"\nPos \t Pot \t dpotx \t dpoty \t dpotz \t dmax",file=fout_py)
        print("\n\n"+header+"\nPos \t Pot \t dpotx \t dpoty \t dpotz \t dmax",file=fout_diff)
        fy=0
        for fz in np.arange(0.9,3,0.1):
            ###############################################################
            # Set Commons (FORTRAN). Not all variables in commons are actually used. 
            # We only set the ones used.
            ###############################################################

            dljpot_f.constants.dipol=dipole
            dljpot_f.coordinates.fx[0]=fx
            dljpot_f.coordinates.fy[0]=fy
            dljpot_f.coordinates.fz[0]=fz
            dljpot_f.charge.pcharge[0]=pcharge
            pot,dpotx,dpoty,dpotz,dmax=dljpot_f.dljpot(x,y,z)
            print (round(fz,1),"\t",pot,"\t",dpotx,"\t",dpoty,"\t",dpotz,"\t",dmax,file=fout_f)

            #######################################
            ## Python code
            #######################################
            vector=np.array([[fx,fy,fz]])
            charge=np.array([pcharge])
            dipol=dipole
            pot2,dpotx2,dpoty2,dpotz2,dmax2=dljpot(np.array([x,y,z]),vector,lj,charge,dipol)
            print (round(fz,1),"\t",pot2,"\t",dpotx2,"\t",dpoty2,"\t",dpotz2,"\t",dmax2,file=fout_py)
            # Difference   
            print (round(fz,1),"\t",pot2-pot,"\t",dpotx2-dpotx,"\t",dpoty2-dpoty,"\t",dpotz2-dpotz,"\t",dmax2-dmax,file=fout_diff)
        
        print("\n\n",file=fout_f)
        print("\n\n",file=fout_py)
        print("\n\n",file=fout_diff)

fout_f.close()
fout_py.close()
fout_diff.close()




