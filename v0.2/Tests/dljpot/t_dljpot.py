'''

 TEST of function dljpot  

We compute and plot the potential as the two particles moves away from each other
No integration of the equations of motion.

'''

import dljpot
import numpy as np

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


def dljpot_py(x,y,z,vector,lj,charge,dipol):
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
dljpot.constants.inatom=1
dljpot.constants.romax=romax
dljpot.ljparameters.eox4[0]=eox4
dljpot.ljparameters.ro6lj[0]=ro6lj
dljpot.ljparameters.ro12lj[0]=ro12lj
dljpot.ljparameters.dro6[0]=dro6
dljpot.ljparameters.dro12[0]=dro12
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
            dljpot.constants.dipol=dipole
            dljpot.coordinates.fx[0]=fx
            dljpot.coordinates.fy[0]=fy
            dljpot.coordinates.fz[0]=fz
            dljpot.charge.pcharge[0]=pcharge
            pot,dpotx,dpoty,dpotz,dmax=dljpot.dljpot(x,y,z)
            print (round(fx,1),"\t",pot,"\t",dpotx,"\t",dpoty,"\t",dpotz,"\t",dmax,file=fout_f)

            #######################################
            ## Python code
            #######################################
            vector=np.array([[fx,fy,fz]])
            charge=np.array([pcharge])
            dipol=dipole
            pot2,dpotx2,dpoty2,dpotz2,dmax2=dljpot_py(x,y,z,vector,lj,charge,dipol)
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
            dljpot.constants.dipol=dipole
            dljpot.coordinates.fx[0]=fx
            dljpot.coordinates.fy[0]=fy
            dljpot.coordinates.fz[0]=fz
            dljpot.charge.pcharge[0]=pcharge
            pot,dpotx,dpoty,dpotz,dmax=dljpot.dljpot(x,y,z)
            print (round(fy,1),"\t",pot,"\t",dpotx,"\t",dpoty,"\t",dpotz,"\t",dmax,file=fout_f)

            #######################################
            ## Python code
            #######################################
            vector=np.array([[fx,fy,fz]])
            charge=np.array([pcharge])
            dipol=dipole
            pot2,dpotx2,dpoty2,dpotz2,dmax2=dljpot_py(x,y,z,vector,lj,charge,dipol)
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

            dljpot.constants.dipol=dipole
            dljpot.coordinates.fx[0]=fx
            dljpot.coordinates.fy[0]=fy
            dljpot.coordinates.fz[0]=fz
            dljpot.charge.pcharge[0]=pcharge
            pot,dpotx,dpoty,dpotz,dmax=dljpot.dljpot(x,y,z)
            print (round(fz,1),"\t",pot,"\t",dpotx,"\t",dpoty,"\t",dpotz,"\t",dmax,file=fout_f)

            #######################################
            ## Python code
            #######################################
            vector=np.array([[fx,fy,fz]])
            charge=np.array([pcharge])
            dipol=dipole
            pot2,dpotx2,dpoty2,dpotz2,dmax2=dljpot_py(x,y,z,vector,lj,charge,dipol)
            print (round(fz,1),"\t",pot2,"\t",dpotx2,"\t",dpoty2,"\t",dpotz2,"\t",dmax2,file=fout_py)
            # Difference   
            print (round(fz,1),"\t",pot2-pot,"\t",dpotx2-dpotx,"\t",dpoty2-dpoty,"\t",dpotz2-dpotz,"\t",dmax2-dmax,file=fout_diff)
        
        print("\n\n",file=fout_f)
        print("\n\n",file=fout_py)
        print("\n\n",file=fout_diff)





fout_f.close()
fout_py.close()
fout_diff.close()




