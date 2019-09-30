c      Test for the full blown output to file
c      Here I only test for the first structure in the file.
c      icoord=1 is hardcoded in this test

c###################################################
c  MODIFIED
c  
c  * I made xrand() in rantate an external funcion
c    In this way I can use the python numdom number
c    generator
c  * We added a common /testing/ to diffeq
c  * mobcal is written to be a module for python. If 
c    we want to compare the two codes numerically
c    I must use the same rng
c###################################################


c
c     PROGRAM MOBCAL
c
c     Program to Calculate Mobilities
c
c     See: M. F. Mesleh, J. M. Hunter, A. A. Shvartsburg, G. C. Schatz,
c     and M. F. Jarrold, "Structural Information from Ion Mobility
c     Measurements: Effects of the Long Range Potential" J. Phys. Chem.
c     1996, 100, 16082-16086; A. A. Shvartsburg and M. F. Jarrold,
c     "An Exact Hard Spheres Scattering Model for the Mobilities of
c     Polyatomic Ions" Chem. Phys. Lett. 1996, 261, 86.
c
c     Version of 06/15/00 junko.f
c     Modified to include Na+ parameters and how masses are read
c     Version of 09/14/98
c     Corrected center of mass calculation in fcoord and ncoord and
c     changed how structural asymmetry parameter is calculated
c     Version of 08/15/98
c     Modified to calculate structural asymmetry parameter
c     Version of 02/20/98
c     Modified GSANG: removed bug that caused trajectory start position
c     to go into endless loop for long structures. Also fixed TESTXRAND
c     and summary print out.
c     Version of 01/29/98
c     Modified MOBIL2: removed a bug that messed up b2max calculation
c     for large molecules
c     Version of 12/24/97
c     Modified to permit averaging over multiple conformations 
c     Version of 12/08/97
c     Extensively modified to allow for a molecule consisting of
c     an unlimited variety of different atoms. Exact hard sphere
c     scattering incorporated. RANLUX random number generator 
c     incorporated. Several obsolete subroutines removed.
c     Version of 10/16/97
c     Modified to allow uniform and non-uniform charge distributions
c     Version of 08/23/97
c     Modified to allow calculations with a non-uniform charge distribution
c     Version of 04/17/96
c
c     MOBIL2 calculates the mobility using a Lennard-Jones plus
c     ion-induced dipole potential. MOBIL2 uses Monte Carlo integrations
c     over orientation and impact parameter and a numerical integration
c     over g* (the reduced velocity). MOBIL2 includes second order
c     corrections, though these don't appear to be important.
c
c     GSANG/DERIV/DIFFEQ are subroutines derived from code provided by
c     George Schatz. They use Runge-Kutta-Gill and Adams-Moulton predictor
c     -corrector integration methods. These subroutines are now set-up to
c     work in three dimensions.
c
c     DLJPOT calculates the potential and the derivatives of the potential.
c     The potential is given by a sum of 6-12 two body interactions and
c     ion-induced dipole interactions. The charge can be distributed over 
c     all atoms equally, a specific charge distribution can be employed, or
c     the charge can be set to zero (only the 6-12 part of the potential
c     is used). Each atom in the cluster or molecule can have different 
c     Lennard-Jones parameters.
c     
c     MOBIL4 calculates the exact hard sphere scattering mobility and
c     the projection approximation mobility. Adapted from code written
c     by Alexandre Shvartsburg.
c
c     XRAND/RANLUX are random number generators. RANLUX is the standard
c     ranlux number generator of M. Luscher (code from F. James). XRAND
c     is a function that calls either RANLUX or RAND (the standard F77
c     random number generator).
c
c     FCOORD reads in the coordinates, integer masses, and partial charges.  
c
c     RANTATE/ROTATE rotates the coordinates of the cluster.  
c
c     TRAJ calculates a series of trajectories for closer examination.
c
c     TRAJONE calculates one trajectory for a given velocity, impact
c     parameter, and angles (theta, phi, and gamma in degrees).
c
c     POTENT determines the average potential. (Only for near spherical
c     molecules or clusters).
c     
c
c     ***************************************************************
c
c
      subroutine Mobcal

      implicit double precision (a-h,m-z)
      dimension tmc(100),tmm(100),ehsc(100),ehsm(100),
     ?pac(100),pam(100),asympp(100)
      character*30 filen1,filen2,unit,dchar,xlabel
      parameter (len=1000)
      common/printswitch/ip,it,iu1,iu2,iu3,iv,im2,im4,igs
      common/constants/mu,ro,eo,pi,cang,ro2,dipol,emax,m1,m2,
     ?xe,xk,xn,mconst,correct,romax,inatom,icoord,iic
      common/charge/pcharge(len)
      common/coordinates/fx(len),fy(len),fz(len),
     ?ox(len),oy(len),oz(len)
      common/ljparameters/eolj(len),rolj(len),eox4(len),
     ?ro6lj(len),ro12lj(len),dro6(len),dro12(len)
      common/hsparameters/rhs(len),rhs2(len)
      common/trajectory/sw1,sw2,dtsf1,dtsf2,cmin,ifail,ifailc,inwr
      common/angles/theta,phi,gamma
      common/xrandom/i1,i2,i3,i4,i5,i6

c
c     ***************************************************************
c
c     To make it easier to run batches of programs the following
c     information is now entered from a file
c
c     1. file for input
c     2. file for output
c     3. random number seed
c
c     filen1='gly3h1'
c     filen2='temp.out'
c     i2=96
c
c     From a file
c
      filen1='a10A1.mfj'
      open (8,file='a10A1_f.out')
c      i2=43
c
c     ***************************************************************
c
c      open (8,file=filen2)
c
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
c
      ip=1
      it=1
      iu1=1
      iu2=0
      iu3=0
      iv=0
      im2=0
      im4=0
      igs=1
c
c     Constants from Handbook of Chemistry and Physics, 70th Edition
c     
      pi=3.14159265358979323846d0
      cang=180.d0/pi
      xe=1.60217733d-19
      xk=1.380658d-23
      xn=6.0221367d23
      xeo=8.854187817d-12
      xmv=0.02241410d0
c
c     Read in coordinates: The input file should have the following
c     structure:
c
c     line 1: label
c     line 2: number of coordinate sets  
c     line 3: total number of atoms
c     line 4: ang (to specify coordinates in angstroms) or 
c             au (for atomic units)
c     line 5: equal (to specify uniform charge distribution) or
c             none (to specify no charge) or
c             calc (to specify a non-uniform charge distribution)
c     line 6: correction factor for coordinates (usually 1.0000)
c     line 7: x, y, z, integer mass, and partial charge (if calc is set)
c     line 8: as above for all atoms. if multiple conformations are 
c             included the conformations are separated by a blank line
c
c
      call fcoord(filen1,unit,dchar,xlabel,asympp(1))
c
c     Lennard-Jones scaling parameters
c
      eo=1.34d-03*xe
      ro=3.043d0*1.0d-10
      ro2=ro*ro
      write(8,600) eo/xe,ro/1.0d-10
  600 format(1x,'Lennard-Jones scaling parameters: eo=',
     ?1pe11.4,1x,'ro=',e11.4)
c
c     Constant for ion-induced dipole potential
c
c     Polarizability of helium = 0.204956d-30 m3
c     xeo is permitivity of vacuum, 8.854187817d-12 F.m-1
c
      dipol=0.204956d-30/(2.d0*4.d0*pi*xeo)
      dipol=dipol*xe*xe 
      write(8,601) dipol
  601 format(1x,'dipole constant =',1pe11.4)
c
c     Mass constants
c
      m1=4.0026d0
      mu=((m1*m2)/(m1+m2))/(xn*1.0d3)
c
c     Mobility constant
c
      mconst=dsqrt(18.d0*pi)/16.d0
      mconst=mconst*dsqrt(xn*1.0d3)*dsqrt((1.d0/m1)+(1.d0/m2))
      mconst=mconst*xe/dsqrt(xk)
      dens=xn/xmv
      mconst=mconst/dens
      write(8,602) mconst
  602 format(1x,'mobility constant =',1pe11.4)
c
c     Temperature
c
      t=298.d0
      write(8,603) t
  603 format(1x,'temperature =',1pe11.4)
c
c     Define parameters for random number generator
c
c     If i5=1 RANLUX is used otherwise RAND is used. If RAND is used 
c     i2 is the seed integer. If RANLUX is used i2, i3, and i4 are seed
c     integers. i3 and i4 (which are used to start RANLUX in a particular
c     place) are usually set to zero. i1 contain the "luxury level" - how
c     good the random number generator is. Values of 0 to 4 can be used
c     and the default value in RANLUX is 3. At this level "any 
c     theoretically possible correlations have a very small chance of
c     being observed".
c
!       i1=3
!       i3=0
!       i4=0
!       i5=1
!       if(i5.ne.1) then 
!       write(8,604) i2
!   604 format(1x,'using RAND with seed integer =',i8)
!       else
!       write(8,619) 
!   619 format(/1x,'using RANLUX',/)
!       endif
! c     initialize the random number generators
!       call rluxgo(i1,i2,i3,i4)
!       call srand(i2)
!       hold=xrand()
c
c     Define parameters for POTENT
c
c     Number of rotations in average potential calculation. If ipr=0
c     an unrotated configuration is used. Otherwise ipr random rotations
c     are employed.
c
      ipr=2
c
c     Define parameters for MOBIL2
c
c     Number of complete cycles for average mobility calculation in 
c     MOBIL2. Default value is 10.
c
      itn=2
c
c     Number of points in velocity integration in MOBIL2. Default 
c     value is 20.
c
      inp=2    
c
c     Number of points in Monte Carlo integrations of impact parameter
c     and orientation in MOBIL2. Default value is 500.
c
      imp=2
c
c     Minimum value of (1-cosX). This quantity determines the maximum
c     impact parameter at each velocity. Default value is 0.0005.
c
      cmin=0.0005
c
c     Define some parameters for trajectories: sw1 defines the potential
c     energy where the trajectory starts and dtsf1 is related to the 
c     time step at the start of the trajectory. As the trajectory comes
c     close to a collision the time step is reduced. sw2 defines the 
c     potential energy where the time step is reduced and dtsf2 defines
c     the reduced time step. Default values are: 
c     sw1 = 0.00005   dtsf1=0.5
c     sw2 = 0.0025    dtsf2=0.05
c     inwr is the number of integration steps before the program tests
c     to see if the trajectory is done or lost. ifail is the number of
c     failed trajectories that are permitted (a failed trajectory is one
c     that does not conserve energy to within 1%. Default values are: 
c     inwr = 1        ifail = 100
c
      sw1=0.00005d0
      sw2=0.005d0 
      dtsf1=0.5d0
      dtsf2=0.1d0
      inwr=1
      ifail=100
      ifailc=0
c
c     Define parameters for MOBIL4
c
c     inum is the number of Monte Carlo trajectories employed and inor
c     is the maximum number of successive reflections followed. Default
c     values inum=600000 and inor=30
c
      inum=10
      inor=10
c
c     Define parameters for TRAJONE
c
      v=2309.9d0
      b=1.5654d-10
      ntheta=33.300d0
      nphi=64.250d0
      ngamma=134.30
c
c
c
      ehsm(1)=0.d0
      tmm(1)=0.d0
      immmax=0
      immmin=inor
c
c     ***************************************************************
c
c     call rantate
c     call potent(ipr) 
c     call traj(t)
c     call trajone(v,b,ntheta,nphi,ngamma)
c     call testxrand
c     goto 9999
c
      icoord=1
      do 1000 iic=1,icoord
      write(8,635)
  635 format(/)
      if(icoord.gt.1) write(8,621) iic
  621 format(/1x,'coordinate set =',i5)
      write(8,637) asympp(iic)
  637 format(/1x,'structural asymmetry parameter =',f8.4)
c

      call mobil4(t,inum,inor,ehsc(iic),ehsm(iic),pac(iic),pam(iic),imm)
      call mobil2(t,itn,inp,imp,tmm(iic),tmc(iic),sdevpc)


      if(imm.lt.immmin) immmin=imm
      if(imm.gt.immmax) immmax=imm  
      if(iic.ne.icoord) call ncoord(unit,dchar,asympp(iic+1))

      im2=1
      im4=1
 1000 continue
c
c     print out summary
c


      write(8,605) filen1,xlabel
  605 format(///1x,'SUMMARY',//1x,'program version = junkn.f',
     ?/1x,'input file name = ',a30,/1x,'input file label = ',a30)
      if(dchar.eq.'equal'.and.tmmob.ne.0.d0) write(8,606) 
  606 format(1x,'using a uniform charge distribution')
      if(dchar.eq.'calc'.and.tmmob.ne.0.d0) write(8,607)
  607 format(1x,'using a calculated (non-uniform) charge distribution')
      if(dchar.eq.'none'.and.tmmob.ne.0.d0) write(8,608)
  608 format(1x,'using no charge - only LJ interactions')
      write(8,616) t
  616 format(1x,'temperature =',1pe11.4)
      if(i5.ne.1) then 
      write(8,617) i2
  617 format(1x,'using RAND with seed integer =',i8)
      else
      write(8,620) i2 
  620 format(1x,'using RANLUX with seed integer =',i8)
      endif 
c

      if(icoord.gt.1) goto 9997
      write(8,636) asympp(1)
  636 format(1x,'structural asymmetry parameter =',1pe11.4) 
      if(ehsm(1).eq.0.d0) goto 9998
      write(8,609)
  609 format(/1x,'mobility calculation by MOBIL4 (HS scattering)')
      write(8,610) inum,imm
  610 format(/1x,'number of Monte Carlo trajectories =',i7,
     ?/1x,'maximum number of reflections encountered =',i3)
      write(8,611) 1.d0/pam(1),pac(1)*1.d20
  611 format(/1x,'inverse average PA mobility =',1pe11.4,
     ?/1x,'average PA cross section =',e11.4)
      write(8,612) 1.d0/ehsm(1),ehsc(1)*1.d20
  612 format(/1x,'inverse average EHS mobility =',1pe11.4,
     ?/1x,'average EHS cross section =',e11.4)
c

 9998 if(tmm(1).eq.0.d0) goto 9999
      write(8,613)
  613 format(/1x,'mobility calculation by MOBIL2 (trajectory method)')
      write(8,614) sw1,sw2,dtsf1,dtsf2
  614 format(/1x,'trajectory parameters',/1x,'sw1 =',1pe11.4,7x,
     ?'sw2 =',e11.4,/1x,'dtsf1 =',e11.4,5x,'dtsf2 =',e11.4)
      write(8,615) itn,inp,imp,itn*inp*imp
  615 format(/1x,'number of complete cycles (itn) =',i6,/1x,
     ?'number of velocity points (inp) =',i6,/1x,
     ?'number of random points (imp) =',i6,/1x,
     ?'total number of points =',i7)
      write(8,618) 1.d0/tmm(1),tmc(1)*1.d20,sdevpc,ifailc
  618 format(/1x,'inverse average (second order) TM mobility =',
     ?1pe11.4,/1x,'average TM cross section =',e11.4,/1x,
     ?'standard deviation (percent) =',e11.4,/1x,
     ?'number of failed trajectories =',i4)
      if(icoord.eq.1) goto 9999
c
c     multiple coordinate sets
c

 9997 if(ehsm(1).eq.0.d0) goto 9996
      write(8,625)
  625 format(/1x,'mobility calculation by MOBIL4 (HS scattering)')
      write(8,626) inum
  626 format(/1x,'number of Monte Carlo trajectories =',i7)
      write(8,634) immmin,immmax
  634 format(/1x,'minimum and maximum number of reflections =',
     ?i3,2x,i3)
      write(8,622)
  622 format(//4x,'set',5x,'PA CS',7x,'PA MOB^-1',6x,'EHSS CS',
     ?6x,'EHSS MOB^-1',4x,'ASYMP')  
      do 1010 i=1,icoord
 1010 write(8,623) i,pac(i)*1.0d20,1.d0/pam(i),
     ?ehsc(i)*1.0d20,1.d0/ehsm(i),asympp(i)/10.d0
  623 format(1x,i5,3x,1pe11.4,3x,e11.4,3x,e11.4,3x,e11.4,3x,f8.4) 
      pacs=0.d0
      pamob=0.d0
      ehscs=0.d0
      ehsmob=0.d0
      aasymp=0.d0
      do 1020 i=1,icoord
      pacs=pacs+pac(i)
      pamob=pamob+pam(i)
      ehscs=ehscs+ehsc(i)
      aasymp=aasymp+asympp(i)
 1020 ehsmob=ehsmob+ehsm(i)  
      pacs=pacs/dfloat(icoord)
      pamob=pamob/dfloat(icoord)
      ehscs=ehscs/dfloat(icoord)
      ehsmob=ehsmob/dfloat(icoord)
      aasymp=aasymp/dfloat(icoord)
      write(8,624) pacs*1.0d20,1.d0/pamob,ehscs*1.0d20,1.d0/ehsmob,
     ?aasymp/10.d0
  624 format(/3x,'AVGE',2x,1pe11.4,3x,e11.4,3x,e11.4,3x,e11.4,3x,f8.4)
c
 9996 if(tmm(1).eq.0.d0) goto 9999
      write(8,627)
  627 format(/1x,'mobility calculation by MOBIL2 (trajectory method)')
      write(8,628) sw1,sw2,dtsf1,dtsf2
  628 format(/1x,'trajectory parameters',/1x,'sw1 =',1pe11.4,7x,
     ?'sw2 =',e11.4,/1x,'dtsf1 =',e11.4,5x,'dtsf2 =',e11.4)
      write(8,629) itn,inp,imp,itn*inp*imp
  629 format(/1x,'number of complete cycles (itn) =',i6,/1x,
     ?'number of velocity points (inp) =',i6,/1x,
     ?'number of random points (imp) =',i6,/1x,
     ?'total number of points =',i7)
      write(8,633) ifailc
  633 format(/1x,'total number of failed trajectories =',i4)
      write(8,630)
  630 format(//4x,'set',5x,'TM CS',7x,'TM MOB^-1')
      do 1030 i=1,icoord
 1030 write(8,631) i,tmc(i)*1.0d20,1.d0/tmm(i)
  631 format(1x,i5,3x,1pe11.4,3x,e11.4)
      tmcs=0.d0
      tmmob=0.d0
      do 1040 i=1,icoord
      tmcs=tmcs+tmc(i)
 1040 tmmob=tmmob+tmm(i)  
      tmcs=tmcs/dfloat(icoord)
      tmmob=tmmob/dfloat(icoord)
      write(8,632) tmcs*1.0d20,1.d0/tmmob
  632 format(/3x,'AVGE',2x,1pe11.4,3x,e11.4)
c
c
c
 9999 continue
      close (8)
      stop
      end
c
c     ***************************************************************
c
      subroutine fcoord(filen1,unit,dchar,xlabel,asymp)
c
c     Reads in coordinates and other parameters.
c
      implicit double precision (a-h,m-z)
      character*30 filen1,unit,dchar,xlabel
      parameter (len=1000)
      dimension imass(len),xmass(len) 
      common/printswitch/ip,it,iu1,iu2,iu3,iv,im2,im4,igs
      common/constants/mu,ro,eo,pi,cang,ro2,dipol,emax,m1,m2,
     ?xe,xk,xn,mconst,correct,romax,inatom,icoord,iic
      common/charge/pcharge(len)
      common/coordinates/fx(len),fy(len),fz(len),
     ?ox(len),oy(len),oz(len)
      common/ljparameters/eolj(len),rolj(len),eox4(len),
     ?ro6lj(len),ro12lj(len),dro6(len),dro12(len)
      common/hsparameters/rhs(len),rhs2(len)
      common/trajectory/sw1,sw2,dtsf1,dtsf2,cmin,ifail,ifailc,inwr
      common/angles/theta,phi,gamma
      common/xrandom/i1,i2,i3,i4,i5,i6
c
      write(8,603) filen1
  603 format(1x,'input file name = ',a30)
      open (9,file=filen1)
      read(9,'(a30)') xlabel
      write(8,601) xlabel
  601 format(1x,'input file label = ',a30)
      read(9,*) icoord
      write(8,650) icoord
  650 format(1x,'number of coordinate sets =',i5)
c 
      read(9,*) inatom
      write(8,612) inatom
  612 format(1x,'number of atoms =',i4)
      read(9,'(a30)') unit
      read(9,'(a30)') dchar
c
      if(unit.eq.'au') write(8,611)
  611 format(1x,'coordinates in atomic units')
      if(unit.eq.'ang') write(8,614)
  614 format(1x,'coordinates in angstroms')
      if(unit.ne.'au'.and.unit.ne.'ang') then
      write(8,610)
  610 format(1x,'units not specified')
      close (8)
      stop
      endif
c
      read(9,*) correct
      write(8,613) correct
  613 format(1x,'correction factor for coordinates =',1pe11.4)
c
      if(dchar.eq.'equal') write(8,630) 
  630 format(1x,'using a uniform charge distribution')
      if(dchar.eq.'calc') write(8,631)
  631 format(1x,'using a calculated (non-uniform) charge distribution')
      if(dchar.eq.'none') write(8,633)
  633 format(1x,'using no charge - only LJ interactions')
      if(dchar.ne.'equal'.and.dchar.ne.'calc'.and.dchar
     ?.ne.'none') then
      write (8,632)
  632 format(1x,'charge distribution not specified')
      close (8)
      stop
      endif
c
      tcharge=0.d0
      acharge=0.d0
      do 2000 iatom=1,inatom
      if(dchar.eq.'calc') then
      read(9,*) fx(iatom),fy(iatom),fz(iatom),
     ?ximass,pcharge(iatom)
      tcharge=tcharge+pcharge(iatom)
      acharge=acharge+dabs(pcharge(iatom))
      else
      read(9,*) fx(iatom),fy(iatom),fz(iatom),ximass
      endif
      imass(iatom)=nint(ximass)
      if(unit.eq.'au') then
      fx(iatom)=fx(iatom)*0.52917706d0
      fy(iatom)=fy(iatom)*0.52917706d0
      fz(iatom)=fz(iatom)*0.52917706d0
      endif
 2000 continue
c
      if(dchar.eq.'equal') then
      do 2011 iatom=1,inatom
 2011 pcharge(iatom)=1.d0/dfloat(inatom)
      endif
c
      if(dchar.eq.'none') then
      do 2012 iatom=1,inatom
 2012 pcharge(iatom)=0.d0
      endif
c
      if(dchar.eq.'calc') write(8,615) tcharge,acharge
  615 format(1x,'total charge =',1pe11.4,/1x,
     ?'total absolute charge =',e11.4)
c
      do 2020 iatom=1,inatom
      itest=0
c             
c     hydrogen (average value of eo from ab initio calculations
c     and ro from fitting mobilities of C6H6 and others) 
c
      if(imass(iatom).eq.1) then
      itest=1
      xmass(iatom)=1.008d0
      eolj(iatom)=0.65d-03*xe
      rolj(iatom)=2.38d0*1.0d-10
      rhs(iatom)=2.2d0*1.0d-10
      endif
c
c     carbon (from fitting C60 mobility)
c
      if(imass(iatom).eq.12) then
      itest=1
      xmass(iatom)=12.01d0
      eolj(iatom)=1.34d-3*xe
      rolj(iatom)=3.043d0*1.0d-10
      rhs(iatom)=2.7d0*1.0d-10
      endif
c
c     nitrogen (same as carbon)
c
      if(imass(iatom).eq.14) then
      itest=1
      xmass(iatom)=14.01d0
      eolj(iatom)=1.34d-3*xe
      rolj(iatom)=3.043d0*1.0d-10
      rhs(iatom)=2.7d0*1.0d-10
      endif
c     
c     oxygen (same as carbon)
c
      if(imass(iatom).eq.16) then
      itest=1
      xmass(iatom)=16.00d0
      eolj(iatom)=1.34d-3*xe
      rolj(iatom)=3.043d0*1.0d-10
      rhs(iatom)=2.7d0*1.0d-10
      endif
c
c     sodium - Na+ from fitting potential derived by Viehland 
c     (Chem. Phys. 85 (1984) 291) for Na+ + He interactions from mobility
c     data. Note 12-6-4 doesn't provide a very good fit to Viehland's 
c     potential. Viehland's potential is flatter at small r. We didn't
c     include Viehland's first three points in the fit.
c     Hard sphere radius from fitting 300 K Na+ low field mobility in He.
c
      if(imass(iatom).eq.23) then
      itest=1
      xmass(iatom)=22.99d0
      eolj(iatom)=0.0278d-3*xe
      rolj(iatom)=(3.97d0/1.12246d0)*1.0d-10
      rhs(iatom)=2.853d0*1.0d-10
      endif
c
c     silicon (from fitting mobilities of small silicon clusters)
c
      if(imass(iatom).eq.28) then
      itest=1
      xmass(iatom)=28.09d0
      eolj(iatom)=1.35d-3*xe
      rolj(iatom)=3.5d0*1.0d-10
      rhs(iatom)=2.95d0*1.0d-10
      endif
c
c     sulfur (same as silicon)
c
      if(imass(iatom).eq.32) then
      itest=1
      xmass(iatom)=32.06d0
      eolj(iatom)=1.35d-3*xe
      rolj(iatom)=3.5d0*1.0d-10
      rhs(iatom)=3.5d0*1.0d-10
      endif
c
c     iron (same as silicon)
c
      if(imass(iatom).eq.56) then
      itest=1
      xmass(iatom)=55.85d0
      eolj(iatom)=1.35d-3*xe
      rolj(iatom)=3.5d0*1.0d-10
      rhs(iatom)=3.5d0*1.0d-10
      endif
c
      if(itest.eq.0) then
      write(8,602) iatom
  602 format(1x,'type not defined for atom number',i3)
      close (8)
      stop
      endif
c
 2020 continue  
      m2=0.d0
      do 2021 iatom=1,inatom
 2021 m2=m2+xmass(iatom)
      write(8,604) m2
  604 format(1x,'mass of ion =',1pd11.4)
      do 2030 iatom=1,inatom
      rhs2(iatom)=rhs(iatom)*rhs(iatom)
      eox4(iatom)=4.d0*eolj(iatom)
      ro2lj=rolj(iatom)*rolj(iatom)
      ro6lj(iatom)=ro2lj*ro2lj*ro2lj
      ro12lj(iatom)=ro6lj(iatom)*ro6lj(iatom)
      dro6(iatom)=6.d0*ro6lj(iatom)
 2030 dro12(iatom)=12.d0*ro12lj(iatom)
c
c
      if(iu1.eq.1) write(8,620)
  620 format(/9x,'initial coordinates',9x,'mass',3x,'charge',
     ?9x,'LJ parameters',/)
c
      fxo=0.d0
      fyo=0.d0
      fzo=0.d0
      do 2009 iatom=1,inatom
      fxo=fxo+(fx(iatom)*xmass(iatom))
      fyo=fyo+(fy(iatom)*xmass(iatom))
 2009 fzo=fzo+(fz(iatom)*xmass(iatom))
      fxo=fxo/m2
      fyo=fyo/m2
      fzo=fzo/m2
      write(8,623) fxo,fyo,fzo
  623 format(1x,'center of mass coordinates = ',1pe11.4,',',e11.4,
     ?',',e11.4)
      do 2010 iatom=1,inatom
      fx(iatom)=(fx(iatom)-fxo)*1.d-10*correct
      fy(iatom)=(fy(iatom)-fyo)*1.d-10*correct
      fz(iatom)=(fz(iatom)-fzo)*1.d-10*correct
      if(iu1.eq.1) write(8,600) fx(iatom),fy(iatom),fz(iatom),
     ?imass(iatom),pcharge(iatom),eolj(iatom)/xe,rolj(iatom)*1.0d10
  600 format(1x,1pe11.4,1x,e11.4,1x,e11.4,1x,i3,1x,
     ?e11.4,1x,e11.4,1x,e11.4)
 2010 continue
      if(iu1.eq.1) write(8,621)
  621 format(/)
c
      if(icoord.eq.1) close (9)
c
      do 3000 iatom=1,inatom
      ox(iatom)=fx(iatom)
      oy(iatom)=fy(iatom)
 3000 oz(iatom)=fz(iatom)
c
      romax=0.d0
      do 3001 iatom=1,inatom
 3001 if(rolj(iatom).gt.romax) romax=rolj(iatom)
c
c     determine structural asymmetry parameter
c
      theta=0.d0
      asymp=0.d0
      do 5000 igamma=0,360,2
      do 5000 iphi=0,180,2
      gamma=dfloat(igamma)/cang
      phi=dfloat(iphi)/cang
      call rotate
      xyzsum=0.d0
      yzsum=0.d0
      do 5005 iatom=1,inatom
      xyz=dsqrt(fx(iatom)**2+fy(iatom)**2+fz(iatom)**2)
      yz=dsqrt(fy(iatom)**2+fz(iatom)**2)
      xyzsum=xyzsum+xyz
      yzsum=yzsum+yz
 5005 continue
      hold=((pi/4.d0)*xyzsum)/yzsum
      if(hold.gt.asymp) asymp=hold
 5000 continue
c
      return
      end
c
c     ***************************************************************
c
      subroutine rantate
c
c     Rotates the cluster/molecule to a random orientation.  
c
      implicit double precision (a-h,m-z)
      parameter (len=1000)
      common/printswitch/ip,it,iu1,iu2,iu3,iv,im2,im4,igs
      common/constants/mu,ro,eo,pi,cang,ro2,dipol,emax,m1,m2,
     ?xe,xk,xn,mconst,correct,romax,inatom,icoord,iic
      common/charge/pcharge(len)
      common/coordinates/fx(len),fy(len),fz(len),
     ?ox(len),oy(len),oz(len)
      common/ljparameters/eolj(len),rolj(len),eox4(len),
     ?ro6lj(len),ro12lj(len),dro6(len),dro12(len)
      common/hsparameters/rhs(len),rhs2(len)
      common/trajectory/sw1,sw2,dtsf1,dtsf2,cmin,ifail,ifailc,inwr
      common/angles/theta,phi,gamma
      common/xrandom/i1,i2,i3,i4,i5,i6
c
c#######################################
cf2py intent(callback, hide) xrand
      external xrand
c#######################################

      rnt=xrand()
      rnp=xrand()
      rng=xrand()
      theta=rnt*2.d0*pi
      phi=dasin((rnp*2.d0)-1.d0)+(pi/2.d0)
      gamma=rng*2.d0*pi
      call rotate
c
      return
      end
c
c     ***************************************************************
c
      subroutine rotate
c
c     Rotates the cluster/molecule.  
c
      implicit double precision (a-h,m-z)
      parameter (len=1000)
      common/printswitch/ip,it,iu1,iu2,iu3,iv,im2,im4,igs
      common/constants/mu,ro,eo,pi,cang,ro2,dipol,emax,m1,m2,
     ?xe,xk,xn,mconst,correct,romax,inatom,icoord,iic
      common/charge/pcharge(len)
      common/coordinates/fx(len),fy(len),fz(len),
     ?ox(len),oy(len),oz(len)
      common/ljparameters/eolj(len),rolj(len),eox4(len),
     ?ro6lj(len),ro12lj(len),dro6(len),dro12(len)
      common/hsparameters/rhs(len),rhs2(len)
      common/trajectory/sw1,sw2,dtsf1,dtsf2,cmin,ifail,ifailc,inwr
      common/angles/theta,phi,gamma
      common/xrandom/i1,i2,i3,i4,i5,i6
c
      if(iu2.eq.1.or.iu3.eq.1) write(8,610) theta*cang,phi*cang,
     ?gamma*cang
  610 format(//1x,'coordinates rotated by ROTATE',//1x,
     ?'theta=',1pe11.4,1x,'phi=',e11.4,1x,'gamma=',1pe11.4,/)
c
      do 1000 iatom=1,inatom
      rxy=dsqrt(ox(iatom)*ox(iatom)+(oy(iatom)*oy(iatom)))
      if(rxy.eq.0.d0) goto 1010
      otheta=dacos(ox(iatom)/rxy)
      if(oy(iatom).lt.0.d0) otheta=(2.d0*pi)-otheta
      ntheta=otheta+theta
 1010 fx(iatom)=dcos(ntheta)*rxy
 1000 fy(iatom)=dsin(ntheta)*rxy
c
      do 2000 iatom=1,inatom
      rzy=dsqrt(oz(iatom)*oz(iatom)+(fy(iatom)*fy(iatom)))
      if(rzy.eq.0.d0) goto 2010
      ophi=dacos(oz(iatom)/rzy)
      if(fy(iatom).lt.0.d0) ophi=(2.d0*pi)-ophi
      nphi=ophi+phi
 2010 fz(iatom)=dcos(nphi)*rzy
 2000 fy(iatom)=dsin(nphi)*rzy
c
      do 3000 iatom=1,inatom
      rxy=dsqrt(fx(iatom)*fx(iatom)+(fy(iatom)*fy(iatom)))
      if(rxy.eq.0.d0) goto 3010
      ogamma=dacos(fx(iatom)/rxy)
      if(fy(iatom).lt.0.d0) ogamma=(2.d0*pi)-ogamma
      ngamma=ogamma+gamma
 3010 fx(iatom)=dcos(ngamma)*rxy
 3000 fy(iatom)=dsin(ngamma)*rxy
c
      if(iu2.eq.0) goto 4000
      write(8,620)
  620 format(9x,'initial coordinates',24x,'new coordinates',/)
      do 4020 iatom=1,inatom
 4020 write(8,600) ox(iatom),oy(iatom),oz(iatom),fx(iatom),
     ?fy(iatom),fz(iatom)
  600 format(1x,1pe11.4,2(1x,e11.4),5x,3(1x,e11.4))
      close (10)
c
 4000 return
      end
c
c     *****************************************************
c
      subroutine dljpot(x,y,z,pot,dpotx,dpoty,dpotz,dmax)
c
c     Subroutine to calculate L-J + ion-dipole potential.
c
      implicit double precision (a-h,m-z)
      parameter (len=1000)
      common/printswitch/ip,it,iu1,iu2,iu3,iv,im2,im4,igs
      common/constants/mu,ro,eo,pi,cang,ro2,dipol,emax,m1,m2,
     ?xe,xk,xn,mconst,correct,romax,inatom,icoord,iic
      common/charge/pcharge(len)
      common/coordinates/fx(len),fy(len),fz(len),
     ?ox(len),oy(len),oz(len)
      common/ljparameters/eolj(len),rolj(len),eox4(len),
     ?ro6lj(len),ro12lj(len),dro6(len),dro12(len)
      common/hsparameters/rhs(len),rhs2(len)
      common/trajectory/sw1,sw2,dtsf1,dtsf2,cmin,ifail,ifailc,inwr
      common/angles/theta,phi,gamma
      common/xrandom/i1,i2,i3,i4,i5,i6
c
      rx=0.d0
      ry=0.d0
      rz=0.d0
      e00=0.d0
      de00x=0.d0
      de00y=0.d0
      de00z=0.d0
      sum1=0.d0
      sum2=0.d0
      sum3=0.d0
      sum4=0.d0
      sum5=0.d0
      sum6=0.d0
      dmax=2.d0*romax
c
      do 1100 iatom=1,inatom
      xx=x-fx(iatom)
      xx2=xx*xx
      yy=y-fy(iatom)
      yy2=yy*yy
      zz=z-fz(iatom)
      zz2=zz*zz
      rxyz2=xx2+yy2+zz2
      rxyz=dsqrt(rxyz2)
      if(rxyz.lt.dmax) dmax=rxyz
      rxyz3=rxyz2*rxyz
      rxyz5=rxyz3*rxyz2
      rxyz6=rxyz5*rxyz
      rxyz8=rxyz5*rxyz3
      rxyz12=rxyz6*rxyz6
      rxyz14=rxyz12*rxyz2
c     LJ potential 
      e00=e00+(eox4(iatom)*((ro12lj(iatom)/rxyz12)-
     ?(ro6lj(iatom)/rxyz6)))
c     LJ derivative
      de00=eox4(iatom)*((dro6(iatom)/rxyz8)-
     ?(dro12(iatom)/rxyz14))
      de00x=de00x+(de00*xx)
      de00y=de00y+(de00*yy)
      de00z=de00z+(de00*zz)
c     ion-induced dipole potential
      if(pcharge(iatom).eq.0.d0) goto 1100
      rxyz3i=pcharge(iatom)/rxyz3
      rxyz5i=-3.d0*pcharge(iatom)/rxyz5
      rx=rx+(xx*rxyz3i)
      ry=ry+(yy*rxyz3i)
      rz=rz+(zz*rxyz3i)
c     ion-induced dipole derivative
      sum1=sum1+(rxyz3i+(xx2*rxyz5i))
      sum2=sum2+(xx*yy*rxyz5i)
      sum3=sum3+(xx*zz*rxyz5i)
      sum4=sum4+(rxyz3i+(yy2*rxyz5i))
      sum5=sum5+(yy*zz*rxyz5i)
      sum6=sum6+(rxyz3i+(zz2*rxyz5i))
c       
 1100 continue
c
      pot=e00-(dipol*((rx*rx)+(ry*ry)+(rz*rz)))
      dpotx=de00x-(dipol*((2.d0*rx*sum1)+(2.d0*ry*sum2)
     ?+(2.d0*rz*sum3)))
      dpoty=de00y-(dipol*((2.d0*rx*sum2)+(2.d0*ry*sum4)
     ?+(2.d0*rz*sum5)))
      dpotz=de00z-(dipol*((2.d0*rx*sum3)+(2.d0*ry*sum5)
     ?+(2.d0*rz*sum6)))
c
      return
      end
c
c     *******************************************************
c
      subroutine potent (ipr)
c
c     Evaluates average potential. Use of this subroutine is
c     only appropriate if all the surface atoms lie an equal
c     distance from the center of mass
c
      implicit double precision (a-h,m-z)
      dimension pot(3000),sdx(200),sdy(200),sdz(200)
      parameter (len=1000)
      common/printswitch/ip,it,iu1,iu2,iu3,iv,im2,im4,igs
      common/constants/mu,ro,eo,pi,cang,ro2,dipol,emax,m1,m2,
     ?xe,xk,xn,mconst,correct,romax,inatom,icoord,iic
      common/charge/pcharge(len)
      common/coordinates/fx(len),fy(len),fz(len),
     ?ox(len),oy(len),oz(len)
      common/ljparameters/eolj(len),rolj(len),eox4(len),
     ?ro6lj(len),ro12lj(len),dro6(len),dro12(len)
      common/hsparameters/rhs(len),rhs2(len)
      common/trajectory/sw1,sw2,dtsf1,dtsf2,cmin,ifail,ifailc,inwr
      common/angles/theta,phi,gamma
      common/xrandom/i1,i2,i3,i4,i5,i6
c
      write(8,680) ipr
  680 format(/1x,'evaluate average potential',//1x,
     ?'number of rotations =',i6)
c
c     evaluate average potential
c
      npp=200
      z=0.d0
      y=0.d0
      do 301 id=1,npp
      sdx(id)=0.d0
      sdy(id)=0.d0
      sdz(id)=0.d0
  301 pot(id)=0.d0
      iprr=1
      if(ipr.eq.0) then
      iprr=0
      ipr=1
      endif
      do 300 i=1,ipr
      if(iprr.ne.0) call rantate
      x=0.d0
      do 300 id=1,npp
      x=dfloat(id)*0.1d-10
      call dljpot(x,y,z,hold,dx,dy,dz,dmax)
      if(iv.eq.1) write(8,611) x,hold,dx,dy,dz
  611 format(1x,1pe11.4,5x,e11.4,5x,e11.4,1x,e11.4,1x,e11.4)
      sdx(id)=sdx(id)+dx
      sdy(id)=sdy(id)+dy
      sdz(id)=sdz(id)+dz
  300 pot(id)=pot(id)+hold
c
c     determine average distance to center of mass
c      
      hold=0.d0
      do 401 i=1,inatom
      temp=dsqrt((fx(i)*fx(i))+(fy(i)*fy(i))+(fz(i)*fz(i)))
  401 hold=hold+temp
      dist=hold/dfloat(inatom)
      write(8,612) dist
  612 format(1x,'average distance of atoms from COM =',1pe11.4)
c
c     write out average potential determine approximate
c     values for emax, rmax, and r00
c
      x=0.d0
      write(8,613)
  613 format(/1x,'Average Potential',/)
      ex=0.d0
      do 400 id=1,npp
      x=dfloat(id)*0.1d-10
      if(pot(id).lt.ex.and.x.gt.dist) then
      ex=pot(id)
      rx=x
      endif
      if(pot(id).gt.0.d0.and.x.gt.dist) r00=x
  400 write(8,611) x/1.0d-10,pot(id)/(xe*dfloat(ipr)),
     ?sdx(id)/dfloat(ipr),sdy(id)/dfloat(ipr),sdz(id)/dfloat(ipr) 
c
c     evaluate potential for narrow range around rmax
c
      do 402 i=1,3000
  402 pot(i)=0.d0
      do 403 i=1,ipr
      if(iprr.ne.0) call rantate
      x=rx-0.2d-10
      do 403 id=1,3000
      x=x+1.333333d-14
      call dljpot(x,y,z,hold,dx,dy,dz,dmax)
  403 pot(id)=pot(id)+hold
c
c     determine accurate values for emax and rmax
c
      x=rx-0.2d-10
      ex=0.d0
      do 404 i=1,3000
      x=x+1.333333d-14
      if(pot(i).lt.ex) then
      ex=pot(i)
      rx=x
      endif
  404 continue
      write(8,614) rx/1.0d-10,ex/(xe*dfloat(ipr))
  614 format(/1x,'rmax =',1pe11.4,3x,'emax =',e11.4)
c
c     evaluate potential for narrow range around r00
c
      do 405 i=1,3000
  405 pot(i)=0.d0
      do 406 i=1,ipr
      if(iprr.ne.0) call rantate
      x=r00-0.2d-10
      do 406 id=1,3000
      x=x+1.333333d-14
      call dljpot(x,y,z,hold,dx,dy,dz,dmax)
  406 pot(id)=pot(id)+hold      
c
c     determine accurate value for r00
c
      x=r00-0.2d-10
      do 407 i=1,3000
      x=x+1.333333d-14
  407 if(pot(i).gt.0.d0) r00=x
      r00=(r00+r00+1.333333d-14)/2.d0
      write(8,615) r00/1.0d-10
  615 format(/1x,'r00 =',1pe11.4)
c
c
      return
      end
c      
c     ***************************************************************
c
      subroutine traj (t)
c
c     Runs a small set of trajectories at a collision energy of kT.
c
      implicit double precision (a-h,m-z)
      parameter (len=1000)
      common/printswitch/ip,it,iu1,iu2,iu3,iv,im2,im4,igs
      common/constants/mu,ro,eo,pi,cang,ro2,dipol,emax,m1,m2,
     ?xe,xk,xn,mconst,correct,romax,inatom,icoord,iic
      common/charge/pcharge(len)
      common/coordinates/fx(len),fy(len),fz(len),
     ?ox(len),oy(len),oz(len)
      common/ljparameters/eolj(len),rolj(len),eox4(len),
     ?ro6lj(len),ro12lj(len),dro6(len),dro12(len) 
      common/hsparameters/rhs(len),rhs2(len)
      common/trajectory/sw1,sw2,dtsf1,dtsf2,cmin,ifail,ifailc,inwr
      common/angles/theta,phi,gamma
      common/xrandom/i1,i2,i3,i4,i5,i6
c
      write(8,603)
  603 format(/1x,'some sample trajectories run by TRAJ',/)
      write(8,604) sw1,sw2,dtsf1,dtsf2,inwr,ifail
  604 format(1x,'global trajectory parameters',//1x,'sw1 =',1pe11.4,7x,
     ?'sw2 =',e11.4,/1x,'dtsf1 =',e11.4,5x,'dtsf2 =',e11.4,/1x,
     ?'inwr =',i3,14x,'ifail =',i5,/)
c
      gst2=0.25d0*xk*t/eo
      v=dsqrt((gst2*eo)/(0.5d0*mu))
      if(it.ne.1) write(8,600) gst2,v
  600 format(/1x,'collision energy = 0.25kT',/1x,'gst2 =',1pe11.4,
     ?5x,'v =',e11.4)
      if(it.eq.1) write(8,630) gst2
  630 format(/1x,'collision energy = 0.25kT',/1x,'gst2 =',1pe11.4)
      if(it.ne.1) write(8,605)
  605 format(/6x,'b',13x,'ang',10x,'e ratio',10x,'d1',
     ?10x,'istep',/)
      bst1=0.5d0
      cosa=1.d0
 1000 b=ro*bst1
      call gsang(v,b,erat,ang,d1,istep)
      if(it.ne.1) write(8,610) b,ang*cang,erat,d1,istep
  610 format(1x,1pe11.4,4x,e11.4,4x,e11.4,4x,e11.4,4x,i5)
      if(it.eq.1) write(8,620) ang*cang,erat
  620 format(/1x,'ang =',1pe11.4,4x,'e ratio =',e11.4)
      cosao=cosa
      cosa=1.d0-dcos(ang)
      if(cosa.gt.1.0d-3.or.cosao.gt.1.0d-3) then
      bst1=bst1+0.5d0
      goto 1000
      endif
c
      gst2=xk*t/eo
      v=dsqrt((gst2*eo)/(0.5d0*mu))
      if(it.ne.1) write(8,601) gst2,v
  601 format(/1x,'collision energy = kT',/1x,'gst2 =',1pe11.4,
     ?5x,'v =',e11.4)
      if(it.eq.1) write(8,631) gst2
  631 format(/1x,'collision energy = kT',/1x,'gst2 =',1pe11.4)
      if(it.ne.1) write(8,606)
  606 format(/6x,'b',13x,'ang',10x,'e ratio',10x,'d1',
     ?10x,'istep',/)
      bst1=0.5d0
      cosa=1.d0
 1001 b=ro*bst1
      call gsang(v,b,erat,ang,d1,istep)
      if(it.ne.1) write(8,611) b,ang*cang,erat,d1,istep
  611 format(1x,1pe11.4,4x,e11.4,4x,e11.4,4x,e11.4,4x,i5)
      if(it.eq.1) write(8,621) ang*cang,erat
  621 format(/1x,'ang =',1pe11.4,4x,'e ratio =',e11.4)
      cosao=cosa
      cosa=1.d0-dcos(ang)
      if(cosa.gt.1.0d-3.or.cosao.gt.1.0d-3) then
      bst1=bst1+0.5d0
      goto 1001
      endif
c
      gst2=8.d0*xk*t/eo
      v=dsqrt((gst2*eo)/(0.5d0*mu))
      if(it.ne.1) write(8,602) gst2,v
  602 format(/1x,'collision energy = 8kT',/1x,'gst2 =',1pe11.4,
     ?5x,'v =',e11.4)
      if(it.eq.1) write(8,632) gst2
  632 format(/1x,'collision energy = 8kT',/1x,'gst2 =',1pe11.4)
      if(it.ne.1) write(8,607)
  607 format(/6x,'b',13x,'ang',10x,'e ratio',10x,'d1',
     ?10x,'istep',/)
      bst1=0.5d0
      cosa=1.d0
 1002 b=ro*bst1
      call gsang(v,b,erat,ang,d1,istep)
      if(it.ne.1) write(8,612) b,ang*cang,erat,d1,istep
  612 format(1x,1pe11.4,4x,e11.4,4x,e11.4,4x,e11.4,4x,i5)
      if(it.eq.1) write(8,622) ang*cang,erat
  622 format(/1x,'ang =',1pe11.4,4x,'e ratio =',e11.4)
      cosao=cosa
      cosa=1.d0-dcos(ang)
      if(cosa.gt.1.0d-3.or.cosao.gt.1.0d-3) then
      bst1=bst1+0.5d0
      goto 1002
      endif
c
      return
      end
c      
c     ***************************************************************
c
      subroutine trajone (v,b,ntheta,nphi,ngamma)
c
c     Runs one trajectory
c
      implicit double precision (a-h,m-z)
      parameter (len=1000)
      common/printswitch/ip,it,iu1,iu2,iu3,iv,im2,im4,igs
      common/constants/mu,ro,eo,pi,cang,ro2,dipol,emax,m1,m2,
     ?xe,xk,xn,mconst,correct,romax,inatom,icoord,iic
      common/charge/pcharge(len)
      common/coordinates/fx(len),fy(len),fz(len),
     ?ox(len),oy(len),oz(len)
      common/ljparameters/eolj(len),rolj(len),eox4(len),
     ?ro6lj(len),ro12lj(len),dro6(len),dro12(len) 
      common/hsparameters/rhs(len),rhs2(len)
      common/trajectory/sw1,sw2,dtsf1,dtsf2,cmin,ifail,ifailc,inwr
      common/angles/theta,phi,gamma
      common/xrandom/i1,i2,i3,i4,i5,i6
c
      write(8,600) v,b
  600 format(//1x,'run one trajectory by TRAJONE ',
     ?//1x,'v =',1pe11.4,4x,'b =',e11.4,/)

      write(8,604) sw1,sw2,dtsf1,dtsf2,inwr
  604 format(1x,'global trajectory parameters',//1x,'sw1 =',1pe11.4,7x,
     ?'sw2 =',e11.4,/1x,'dtsf1 =',e11.4,5x,'dtsf2 =',e11.4,/1x,
     ?'inwr =',i3)
c
      iu2=1
      theta=ntheta/cang
      phi=nphi/cang
      gamma=ngamma/cang
      call rotate
      iu2=0
      it=1
      call gsang(v,b,erat,ang,d1,istep)
      it=0
      write(8,610) ang*cang,erat
  610 format(/1x,'ang =',1pe11.4,4x,'e ratio =',e11.4)
c
      return
      end
c
c     ****************************************************************
c
      subroutine gsang(v,b,erat,ang,d1,istep)
c
c     Calculates trajectory. Adapted from code written by George Schatz
c     A fixed step integrator is used which combines a runge-kutta-gill
c     initiator with an adams-moulton predictor-corrector propagator.
c     
      implicit double precision (a-h,m-z)
      integer ns,nw,l
      dimension w(6),dw(6)
      parameter (len=1000)
      common/printswitch/ip,it,iu1,iu2,iu3,iv,im2,im4,igs
      common/constants/mu,ro,eo,pi,cang,ro2,dipol,emax,m1,m2,
     ?xe,xk,xn,mconst,correct,romax,inatom,icoord,iic
      common/charge/pcharge(len)
      common/coordinates/fx(len),fy(len),fz(len),
     ?ox(len),oy(len),oz(len)
      common/ljparameters/eolj(len),rolj(len),eox4(len),
     ?ro6lj(len),ro12lj(len),dro6(len),dro12(len)
      common/hsparameters/rhs(len),rhs2(len)
      common/trajectory/sw1,sw2,dtsf1,dtsf2,cmin,ifail,ifailc,inwr
      common/angles/theta,phi,gamma
      common/xrandom/i1,i2,i3,i4,i5,i6
c
      vy=-v
      vx=0.d0
      vz=0.d0
      vxyz=dabs(vy)
      if(it.eq.1) write(8,108) v,b
  108 format(/1x,'specific trajectory parameters',//1x,
     ?'v =',1pe11.4,4x,'b =',e11.4)
c
c     determine time step
c
      top=(v/95.2381d0)-0.5d0
      if(v.ge.1000.d0) top=10.d0
      if(v.ge.2000.d0) top=10.0d0-((v-2000.d0)*7.5d-3)
      if(v.ge.3000.d0) top=2.5d0
c
      dt1=top*dtsf1*1.0d-11/v
      dt2=dt1*dtsf2
      dt=dt1
      if(it.eq.1) write (8,107) dt1,dt2
  107 format(1x,'time steps, dt1 =',1pe11.4,' dt2 =',e11.4)
c
c     determine trajectory start position
c
      e0=0.5d0*mu*v*v
      x=b
      z=0.d0     
c
      ymin=0.d0
      ymax=0.d0
      do 200 i=1,inatom
      if(fy(i).gt.ymax) ymax=fy(i)
  200 if(fy(i).lt.ymin) ymin=fy(i)
      ymax=ymax/1.0d-10
      ymin=ymin/1.0d-10
      iymin=dint(ymin)-1
      iymax=dint(ymax)+1
c
      id2=iymax
      y=dfloat(id2)*1.0d-10
      call dljpot(x,y,z,pot,dpotx,dpoty,dpotz,dmax)
      if(dabs(pot/e0).gt.sw1) goto 302 
  300 id2=id2-1
      y=dfloat(id2)*1.0d-10
      call dljpot(x,y,z,pot,dpotx,dpoty,dpotz,dmax)
      if(id2.lt.iymin) then
      if(it.eq.1) write(8,621)
  621 format(1x,'trajectory not started - potential too small')
      ang=0.d0
      erat=1.000d0
      return
      endif
      y=dfloat(id2)*1.0d-10
      call dljpot(x,y,z,pot,dpotx,dpoty,dpotz,dmax)
      if(dabs(pot/e0).lt.sw1) goto 300
      goto 304
c
  302 id2=id2+10
      y=dfloat(id2)*1.0d-10
      call dljpot(x,y,z,pot,dpotx,dpoty,dpotz,dmax)
      if(dabs(pot/e0).gt.sw1) goto 302
  301 id2=id2-1
      y=dfloat(id2)*1.0d-10
      call dljpot(x,y,z,pot,dpotx,dpoty,dpotz,dmax)
      if(dabs(pot/e0).lt.sw1) goto 301
c
  304 y=dfloat(id2)*1.0d-10
      etot=e0+pot
      if(it.eq.1) write(8,622) y*1.0d10
  622 format(1x,'trajectory start position =',1pe11.4)
      d1=y
c
c     initial coordinates and momenta
c
      w(1)=x
      w(2)=vx*mu
      w(3)=y
      w(4)=vy*mu
      w(5)=z
      w(6)=vz*mu
      tim=0.d0
      if(it.eq.1) write(8,623)  
  623 format(//1x,'trajectory ns, x,  y,  z,  kin e, dt,    tot e',/1x,
     ?'               vx, vy, vz, pot e, pot/e0'/)    
c
c     initialize the time derivatives of the coordinates and momenta
c
      call deriv(w,dw,pot,dpotx,dpoty,dpotz,dmax)
      ns=0
      nw=0
      l=0
   15 call diffeq(l,tim,dt,w,dw,pot,dmax)
      nw=nw+1
      if (nw.ne.inwr) goto 15
      ns=ns+nw
      nw=0
c
c     print out the trajectory coordinates and velocities
c
      if(it.eq.1) then
      e=w(2)**2/(2.d0*mu)+w(4)**2/(2.d0*mu)+w(6)**2/(2.d0*mu)
      write(8,103) ns,w(1),w(3),w(5),e,dt,pot+e
  103 format(1x,i5,1x,1pe11.4,5(1x,e11.4))
      write(8,106) dw(1),dw(3),dw(5),pot,dabs(pot/e0)
  106 format(7x,1pe11.4,4(1x,e11.4))
      endif
c
c     check if trajectory has become "lost" (too many steps)
c
      if(ns.gt.30000) then
      write(8,105) b,v
  105 format(1x,'trajectory lost: b =',1pe11.4,' v =',e11.4)
      ang=pi/2.d0
      e=0.5d0*mu*(dw(1)**2+dw(3)**2+dw(5)**2)
      erat=(e+pot)/etot
      istep=ns
      return
      endif
c
c     check if the trajectory is finished
c
      if(dmax.lt.romax) goto 15
      if(dabs(pot/e0).gt.sw2.and.dt.eq.dt1) then
      dt=dt2
      l=0
      endif
      if(dabs(pot/e0).lt.sw2.and.dt.eq.dt2) then
      dt=dt1 
      l=0
      endif
      if(dabs(pot/e0).gt.sw1) goto 15
      if(ns.lt.50) goto 15
      istep=ns
c
c     determine scattering angle 
c
      if(dw(1).gt.0.d0) then
      num=dw(3)*(-v)
      den=v*dsqrt(dw(1)**2+dw(3)**2+dw(5)**2)
      ang=dacos(num/den)
      endif
      if(dw(1).lt.0.d0) then
      num=dw(3)*(-v)
      den=v*dsqrt(dw(1)**2+dw(3)**2+dw(5)**2)
      ang=(-dacos(num/den))
      endif
c
c     check for energy conservation
c
      e=0.5d0*mu*(dw(1)**2+dw(3)**2+dw(5)**2)
      erat=(e+pot)/etot
      if(erat.lt.1.01d0.and.erat.gt.0.99d0) return
      write(8,104) erat,v,b
  104 format(/1x,'energy not conserved: e ratio =',1pe11.4,
     ?' v =',e11.4,' b =',e11.4)
      write(8,109) 0.5d0*mu*v*v/eo,theta*cang,phi*cang,
     ?gamma*cang
  109 format(1x,'gst2 =',1pe11.4,' theta =',e11.4,' phi =',
     ?e11.4,' gamma =',e11.4)
      if(ip.eq.1) write(8,102)
  102 format()
      ifailc=ifailc+1
      if(ifailc.eq.ifail) stop
      return
c       
c
      end
c
c     *******************************************************
c
      subroutine diffeq(l,tim,dt,w,dw,pot,dmax)
c
c     Integration subroutine - uses 5th order runge-kutta-gill to 
c     initiate and 5th order adams-moulton predictor-corrector to 
c     propagate. Parameter l is initially set to zero and then 
c     incremented to tell the subroutine when to switch between 
c     integration methods. DIFFEQ calls subroutine DERIV to define 
c     the equations of motion to be integrated.
c
      implicit double precision (a-h,m-z)
      dimension w(6),dw(6),a(4),b(4),c(4),ampc(5),amcc(4),array(6,40),
     ?savw(40),savdw(40),q(40)
      parameter (len=1000) 
      common/printswitch/ip,it,iu1,iu2,iu3,iv,im2,im4,igs
      common/constants/mu,ro,eo,pi,cang,ro2,dipol,emax,m1,m2,
     ?xe,xk,xn,mconst,correct,romax,inatom,icoord,iic
      common/charge/pcharge(len)
      common/coordinates/fx(len),fy(len),fz(len),
     ?ox(len),oy(len),oz(len)
      common/ljparameters/eolj(len),rolj(len),eox4(len), 
     ?ro6lj(len),ro12lj(len),dro6(len),dro12(len)
      common/hsparameters/rhs(len),rhs2(len)
      common/trajectory/sw1,sw2,dtsf1,dtsf2,cmin,ifail,ifailc,inwr
      common/angles/theta,phi,gamma
      common/xrandom/i1,i2,i3,i4,i5,i6
      common/testing/hvar,hcvar,q,array
      data a/0.50d0,0.292893218814d0,1.70710678118d0,0.1666666666667d0/
      data b/2.0d0,1.0d0,1.0d0,2.0d0/
      data c/-0.5d0,-0.292893218814d0,-1.70710678118d0,-0.5d0/
      data ampc/-0.111059153612d0,0.672667757774d0,-1.70633621697d0,
     ?2.33387888707d0,-1.8524668225d0/
      data amcc/0.0189208128941d0,-0.121233356692d0,0.337771548703d0,
     ?-0.55921513665d0/
      data var,cvar,acst/2.97013888888d0,0.990972222222d0,
     ?0.332866152768d0/
c
      if (l) 4,1,3
    1 do 2 j=1,6
    2 q(j)=0.0
      hvar=dt*var
      hcvar=dt*cvar
      dt=0.5*dt
    3 l=l+1
c
c     This is the runge-kutta-gill part...the steps are broken up into
c     half steps to improve accuracy.
c
      k=0
   15 do 7 j=1,4
      if ((-1)**j.gt.0) tim=tim+0.5*dt
      call deriv(w,dw,pot,dpotx,dpoty,dpotz,dmax)
      do 7 i=1,6
      dw(i)=dt*dw(i)
      r=a(j)*(dw(i)-b(j)*q(i))
      w(i)=w(i)+r
    7 q(i)=q(i)+3.0*r+c(j)*dw(i)
      call deriv(w,dw,pot,dpotx,dpoty,dpotz,dmax)
      if (k) 13,13,14
   13 k=1
      goto 15
   14 if (l-6) 5,6,6
    6 l=-1
      dt=2.0*dt
      return
    5 do 8 j=1,6
    8 array(l,j)=dw(j)
      return
C
c     This is the adams-moulton predictor-corrector part.
c
    4 do 10 j=1,6
      savw(j)=w(j)
      savdw(j)=dw(j)
      array(6,j)=savdw(j)
      do 9 i=1,5
    9 array(6,j)=array(6,j)+ampc(i)*array(i,j)
   10 w(j)=array(6,j)*hvar+w(j)
      tim=tim+dt
      call deriv(w,dw,pot,dpotx,dpoty,dpotz,dmax)
      do 12 j=1,6
      array(6,j)=acst*dw(j)
      do 11 i=1,4
      array(i,j)=array(i+1,j)
   11 array(6,j)=array(i,j)*amcc(i)+array(6,j)
      array(5,j)=savdw(j)
   12 w(j)=savw(j)+hcvar*(array(5,j)+array(6,j))
      call deriv(w,dw,pot,dpotx,dpoty,dpotz,dmax)
      return
      end
c
c     ***************************************************************
c
      subroutine deriv(w,dw,pot,dpotx,dpoty,dpotz,dmax)
c
c     Defines Hamilton's equations of motion as the time derivatives 
c     of the coordinates and momenta.
c
      implicit double precision (a-h,m-z)
      dimension w(6),dw(6)
      parameter (len=1000)
      common/printswitch/ip,it,iu1,iu2,iu3,iv,im2,im4,igs
      common/constants/mu,ro,eo,pi,cang,ro2,dipol,emax,m1,m2,
     ?xe,xk,xn,mconst,correct,romax,inatom,icoord,iic
      common/charge/pcharge(len)
      common/coordinates/fx(len),fy(len),fz(len),
     ?ox(len),oy(len),oz(len)
      common/ljparameters/eolj(len),rolj(len),eox4(len),
     ?ro6lj(len),ro12lj(len),dro6(len),dro12(len)
      common/hsparameters/rhs(len),rhs2(len)
      common/trajectory/sw1,sw2,dtsf1,dtsf2,cmin,ifail,ifailc,inwr
      common/angles/theta,phi,gamma
      common/xrandom/i1,i2,i3,i4,i5,i6
c
c     From Hamilton's equations, the time derivatives of the coordinates
c     are the conjugates divided by the mass.
c
      dw(1)=w(2)/mu
      dw(3)=w(4)/mu
      dw(5)=w(6)/mu
c
c     Hamilton's equations for the time derivatives of the momenta
c     evaluated by using the coordinate derivatives together with the
c     chain rule.
c
      x=w(1)
      y=w(3)
      z=w(5)
c
c    These are analytical derivatives.
c
      call dljpot(x,y,z,pot,dpotx,dpoty,dpotz,dmax)
      dw(2)=-(dpotx)
      dw(4)=-(dpoty)
      dw(6)=-(dpotz)
c
c
      return
      end
c
c     ***************************************************************
c
      subroutine mobil2 (t,itn,inp,imp,mob,cs,sdevpc)
c
c     Subroutine to determine average mobility by trajectory method. 
c     All integrations Monte Carlo (except over velocity).
c
      implicit double precision (a-h,m-z)
      dimension pgst(100),wgst(100),b2max(100)
      dimension q1st(100),q2st(100),cosx(0:500)
      dimension om11st(100),om12st(100),om13st(100),om22st(100)
      parameter (len=1000)
      common/printswitch/ip,it,iu1,iu2,iu3,iv,im2,im4,igs
      common/constants/mu,ro,eo,pi,cang,ro2,dipol,emax,m1,m2,
     ?xe,xk,xn,mconst,correct,romax,inatom,icoord,iic
      common/charge/pcharge(len)
      common/coordinates/fx(len),fy(len),fz(len),
     ?ox(len),oy(len),oz(len)
      common/ljparameters/eolj(len),rolj(len),eox4(len),
     ?ro6lj(len),ro12lj(len),dro6(len),dro12(len)
      common/hsparameters/rhs(len),rhs2(len)
      common/trajectory/sw1,sw2,dtsf1,dtsf2,cmin,ifail,ifailc,inwr
      common/angles/theta,phi,gamma
      common/xrandom/i1,i2,i3,i4,i5,i6
c
c#######################################
cf2py intent(callback, hide) xrand
      external xrand
c#######################################

      if(im2.eq.0) write(8,631)
  631 format(/1x,'mobility calculation by MOBIL2 (trajectory method)',/)
      if(im2.eq.0) write(8,603) sw1,sw2,dtsf1,dtsf2,inwr,ifail
  603 format(1x,'global trajectory parameters',//1x,'sw1 =',1pe11.4,7x,
     ?'sw2 =',e11.4,/1x,'dtsf1 =',e11.4,5x,'dtsf2 =',e11.4,/1x,
     ?'inwr =',i3,14x,'ifail =',i5)
      it=0
      iu2=0
c
c     determine maximum extent and orientate along x axis
c
      if(im2.eq.0) write(8,632)
  632 format(/1x,'maximum extent orientated along x axis')
      rmax=0.d0
      do 1000 iatom=1,inatom
      r=dsqrt((ox(iatom)*ox(iatom))+(oy(iatom)*oy(iatom))+
     ?(oz(iatom)*oz(iatom)))
      if(r.gt.rmax) then
      rmax=r
      ihold=iatom
      endif
 1000 continue
c
      rzy=dsqrt((oz(ihold)*oz(ihold))+(oy(ihold)*oy(ihold)))
      phi=dacos(oz(ihold)/rzy)
      phi=phi+(pi/2.d0)
      if(oy(ihold).lt.0.d0) phi=(2.d0*pi)-phi
      phi=(2.d0*pi)-phi
      theta=0.d0 
      gamma=0.d0
      call rotate
      rxy=dsqrt((fx(ihold)*fx(ihold))+(fy(ihold)*fy(ihold)))
      gamma=dacos(fx(ihold)/rxy)
      if(fy(ihold).lt.0.d0) gamma=(2.d0*pi)-gamma
      gamma=(2.d0*pi)-gamma
      if(im2.eq.0) iu3=1
      if(ip.eq.1) iu2=1
      call rotate
      iu3=0
      iu2=0
      hold=fx(ihold)/rmax
      if(hold.lt.0.9999999999d0.or.hold.gt.1.0000000001d0.or.
     ?fy(ihold).gt.1.0d-20.or.fz(ihold).gt.1.0d-20.or.
     ?fy(ihold).lt.-1.0d-20.or.fz(ihold).lt.-1.0d-20) then 
      write(8,601)
  601 format(/1x,'Problem orientating along x axis',/)
      do 1001 iatom=1,inatom
      hold=dsqrt(fx(iatom)**2+fy(iatom)**2+fz(iatom)**2)
      write(8,602) iatom,fx(iatom),fy(iatom),fz(iatom),hold
  602 format(1x,i4,5x,1pe11.4,3(5x,e11.4))
 1001 continue     
      close (8)
      stop
      endif
c
c     determine rmax, emax, and r00 along x, y, and z directions
c
      if(ip.eq.1) write(8,689)
  689 format(/)
      irn=1000
      ddd=(rmax+romax)/dfloat(irn)
c
      y=0.d0
      z=0.d0
      emaxx=0.d0
      do 1101 ir=1,irn
      x=rmax+romax-(dfloat(ir)*ddd) 
      call dljpot(x,y,z,pot,dpotx,dpoty,dpotz,dmax)
      if(pot.gt.0.d0) goto 1111
      r00x=x
      if(pot.lt.emaxx) then
      rmaxx=x
      emaxx=pot
      endif
 1101 continue
 1111 continue
      if(im2.eq.0) write(8,614) emaxx/xe,rmaxx*1.0d10,r00x*1.0d10        
  614 format(1x,'along x axis emax =',1pe11.4,'eV rmax =',
     ?e11.4,'A r00 =',e11.4,'A')
c
      x=0.d0
      z=0.d0
      emaxy=0.d0
      do 1100 ir=1,irn
      y=rmax+romax-(dfloat(ir)*ddd)
      call dljpot(x,y,z,pot,dpotx,dpoty,dpotz,dmax)
      if(pot.gt.0.d0) goto 1110
      r00y=y
      if(pot.lt.emaxy) then
      rmaxy=y
      emaxy=pot
      endif
 1100 continue
 1110 continue
      if(im2.eq.0) write(8,613) emaxy/xe,rmaxy*1.0d10,r00y*1.0d10
  613 format(1x,'along y axis emax =',1pe11.4,'eV rmax =',
     ?e11.4,'A r00 =',e11.4,'A')
c
      x=0.d0
      y=0.d0
      emaxz=0.d0
      do 1102 ir=1,irn
      z=rmax+romax-(dfloat(ir)*ddd)
      call dljpot(x,y,z,pot,dpotx,dpoty,dpotz,dmax)
      if(pot.gt.0.d0) goto 1112
      r00z=z
      if(pot.lt.emaxz) then
      rmaxz=z
      emaxz=pot
      endif
 1102 continue
 1112 continue
      if(im2.eq.0) write(8,615) emaxz/xe,rmaxz*1.0d10,r00z*1.0d10
  615 format(1x,'along z axis emax =',1pe11.4,'eV rmax =',
     ?e11.4,'A r00 =',e11.4,'A',/)
c
c     set-up integration over gst
c     
      tst=xk*t/eo
      if(im2.eq.0) write(8,600) tst
  600 format(/1x,'t*=',1pe11.4)
      tst3=tst*tst*tst
c
      dgst=5.0d-7*6.d0*dsqrt(tst)
      gst=dgst
      sum=0.d0
      sum1=0.d0
      sum2=0.d0
c
      do 2020 i=1,inp
 2020 sum1=sum1+dsqrt(dfloat(i))
c
      if(im2.eq.0) write(8,611)
  611 format(//1x,'set-up gst integration - integration over velocity',
     ?//5x,'pgst',8x,'wgst',9x,'v',9x,'ke/kt',7x,'gst^5*',5x,'frac of',
     ?/48x,'exp(gst^2/tst)',3x,'sum',/)
      do 2000 i=1,inp
      hold1=dsqrt(dfloat(i))
      hold2=dsqrt(dfloat(i-1))
      sum2=sum2+hold2
      wgst(i)=hold1/sum1
      gstt=tst3*(sum2+(hold1/2.d0))/sum1
c
 2010 sum=sum+(dexp(-gst*gst/tst)*gst*gst*gst*gst*gst*dgst)
      gst=gst+dgst
      if(sum.gt.gstt) pgst(i)=gst-(dgst/2.d0)
      if(sum.lt.gstt) goto 2010
c
      hold1=dsqrt((pgst(i)*pgst(i)*eo)/(0.5d0*mu))
      hold2=0.5d0*mu*hold1*hold1/(xk*t)
      hold3=dexp(-pgst(i)*pgst(i)/tst)*pgst(i)**5.d0
      if(im2.eq.0) write(8,610) pgst(i),wgst(i),hold1,hold2,
     ?hold3,sum/tst3
  610 format(1x,1pe11.4,5(1x,e11.4))
 2000 continue
c
c     determine b2max
c
      dbst2=1.d0
      dbst22=dbst2/10.d0
      cmin=0.0005
      if(im2.eq.0) write(8,652) cmin
  652 format(//1x,'set up b2 integration - integration over',
     ?' impact parameter',//1x,
     ?'minimum value of (1-cosX) =',1pe11.4,/)
c
      do 3030 ig=inp,1,-1
      gst2=pgst(ig)*pgst(ig)
      v=dsqrt((gst2*eo)/(0.5d0*mu))
      ibst=dint(rmaxx/ro)-6
      if(ig.lt.inp)ibst=dint(b2max(ig+1)/dbst2)-6
      if(ibst.lt.0) ibst=0
      if(ip.eq.1) write(8,650) gst2,v
  650 format(/1x,'gst2 =',1pe11.4,1x,'v =',e11.4,/6x,'b',
     ?10x,'bst2',7x,'X ang',7x,'cos(X)',6x,'e ratio')
 3000 bst2=dbst2*dfloat(ibst)
      b=ro*dsqrt(bst2)
      call gsang(v,b,erat,ang,d1,istep)
      cosx(ibst)=1.d0-dcos(ang)
      if(ip.eq.1) write(8,651) b,bst2,ang,cosx(ibst),erat
  651 format(1x,1pe11.4,6(1x,e11.4))
      if(ibst.lt.5) goto 3010
      if(cosx(ibst).lt.cmin.and.cosx(ibst-1).lt.cmin.and.
     ?cosx(ibst-2).lt.cmin.and.cosx(ibst-3).lt.cmin.and.
     ?cosx(ibst-4).lt.cmin) goto 3020
 3010 ibst=ibst+1
      if(ibst.gt.500) then
      write(8,653)
  653 format(1x,'ibst greater than 500')
      close (8)
      stop     
      endif
      goto 3000
 3020 b2max(ig)=dfloat(ibst-5)*dbst2
 3040 b2max(ig)=b2max(ig)+dbst22
      b=ro*dsqrt(b2max(ig))
      call gsang(v,b,erat,ang,d1,istep)
      if(1.d0-dcos(ang).gt.cmin) goto 3040
 3030 continue
      if(im2.eq.0) then
      write(8,637) 
  637 format(/5x,'gst',11x,'b2max/ro2',9x,'b/A',/)
      do 3050 ig=1,inp
 3050 write(8,630) pgst(ig),b2max(ig),ro*dsqrt(b2max(ig))*1.0d10
  630 format(1x,1pe11.4,5x,e11.4,5x,e11.4)
      endif
c
c     Calculate Omega(1,1)*, Omega(1,2)*, Omega(1,3)*, and Omega(2,2)*
c     by integrating Q(1)* or Q(2)* over all orientations, and initial 
c     relative velocities.
c
      if(im2.eq.0) write(8,672) itn,inp,imp,itn*inp*imp
  672 format(//1x,'number of complete cycles (itn) =',i6,/1x,
     ?'number of velocity points (inp) =',i6,/1x,
     ?'number of random points (imp) =',i6,/1x,
     ?'total number of points =',i7,/)
c
      if(ip.eq.1) write(8,680)
  680 format(1x,'start mobility calculation')
      do 4011 ig=1,inp
      q1st(ig)=0.d0
      q2st(ig)=0.d0
 4011 continue
c
      do 4040 ic=1,itn
      if(ip.eq.1) write(8,681) ic
  681 format(/1x,'cycle number, ic =',i3)
      om11st(ic)=0.d0
      om12st(ic)=0.d0
      om13st(ic)=0.d0
      om22st(ic)=0.d0
      do 4010 ig=1,inp
      gst2=pgst(ig)*pgst(ig)
      v=dsqrt((gst2*eo)/(0.5d0*mu))
      if(ip.eq.1) write(8,682) ic,ig,gst2,v
  682 format(/1x,'ic =',i3,1x,'ig =',i4,1x,'gst2 =',1pe11.4,
     ?1x,'v =',e11.4)
      temp1=0.d0
      temp2=0.d0
c
      do 4000 im=1,imp
      if(ip.eq.1.and.im.eq.1) write(8,683)
  683 format(/5x,'b/A',8x,'ang',6x,'(1-cosX)',4x,'e ratio',4x,'theta',
     ?7x,'phi',7x,'gamma') 
      rnb=xrand()
      call rantate
      bst2=rnb*b2max(ig)
      b=ro*dsqrt(bst2)
      if(igs.eq.1) then
      open(15,file='hold')
      write(15,*) iic,ic,ig,im
      write(15,*) v,b
      write(15,*) theta*cang,phi*cang,gamma*cang
      close(15)
      endif
      call gsang(v,b,erat,ang,d1,istep)
      hold1=1.d0-dcos(ang)
      hold2=dsin(ang)*dsin(ang)
      if(ip.eq.1) write(8,684) b*1.d10,ang*cang,hold1,erat,
     ?theta*cang,phi*cang,gamma*cang
  684 format(1x,1pe11.4,7(e11.4))
      temp1=temp1+(hold1*b2max(ig)/dfloat(imp))
      temp2=temp2+(1.5d0*hold2*b2max(ig)/dfloat(imp))
 4000 continue
c
      om11st(ic)=om11st(ic)+(temp1*wgst(ig))
      om12st(ic)=om12st(ic)+(temp1*pgst(ig)*pgst(ig)*wgst(ig)*
     ?(1.d0/(3.d0*tst)))
      om13st(ic)=om13st(ic)+(temp1*(pgst(ig)**4)*wgst(ig)*
     ?(1.d0/(12.d0*tst*tst)))
      om22st(ic)=om22st(ic)+(temp2*pgst(ig)*pgst(ig)*wgst(ig)*
     ?(1.d0/(3.d0*tst)))
      q1st(ig)=q1st(ig)+temp1
      q2st(ig)=q2st(ig)+temp2
      if(ip.eq.1) write(8,670) v,q1st(ig) 
  670 format(/1x,'v =',1pe11.4,5x,'q1st =',e11.4,/)
 4010 continue
c
      if(ip.eq.1) write(8,620) om11st(ic)
  620 format(/1x,'OMEGA(1,1)*=',1pe11.4,/)
 4040 continue
c
c     calculate running averages
c
      hold1=0.d0
      hold2=0.d0
      if(im2.eq.0) write(8,685)
  685 format(/1x,'summary of mobility calculations',//1x,'cycle',
     ?5x,'cs/A^2',6x,'avge cs/A^2',8x,'Ko^-1',7x,'avge Ko^-1')
      do 4041 icc=1,itn
      temp=1.d0/(mconst/(dsqrt(t)*om11st(icc)*pi*ro*ro))
      hold1=hold1+om11st(icc)
      hold2=hold2+temp
 4041 if(im2.eq.0) write(8,622) icc,om11st(icc)*pi*ro*ro*1.d20,
     ?hold1*pi*ro*ro*1.d20/dfloat(icc),temp,hold2/dfloat(icc)
  622 format(1x,i3,4x,1pe11.4,4x,e11.4,4x,e11.4,4x,e11.4)
c
      if(im2.eq.0) then
      write(8,675)
  675 format(//1x,'average values for q1st',//5x,
     ?'gst2',8x,'wgst',8x,'q1st')
      do 4012 ig=1,inp
 4012 write(8,676) pgst(ig)*pgst(ig),wgst(ig),q1st(ig)/dfloat(inp)
  676 format(1x,1pe11.4,1x,e11.4,1x,e11.4)
      endif
c
      mom11st=0.d0
      mom12st=0.d0
      mom13st=0.d0
      mom22st=0.d0
      do 4050 ic=1,itn
      mom11st=mom11st+om11st(ic)
      mom12st=mom12st+om12st(ic)
      mom13st=mom13st+om13st(ic)
      mom22st=mom22st+om22st(ic)
 4050 continue
      mom11st=mom11st/dfloat(itn)
      mom12st=mom12st/dfloat(itn)
      mom13st=mom13st/dfloat(itn)
      mom22st=mom22st/dfloat(itn)
      sdom11st=0.d0
      do 4060 ic=1,itn
      hold=mom11st-om11st(ic)
 4060 sdom11st=sdom11st+(hold*hold)
      sdom11st=dsqrt(sdom11st/dfloat(itn))
      sterr=sdom11st/dsqrt(dfloat(itn))
      if(im2.eq.0) write(8,674) mom11st,sdom11st,sterr
  674 format(//1x,'mean OMEGA*(1,1) =',1pe11.4,/1x,
     ?'standard deviation =',
     ?e11.4,/1x,'standard error of mean =',e11.4)
      cs=mom11st*pi*ro*ro
      sdevpc=100.d0*sdom11st/mom11st
c
c     Use omegas to obtain higher order correction factor to mobility
c
      ayst=mom22st/mom11st
      best=((5.d0*mom12st)-(4.d0*mom13st))/mom11st
      cest=mom12st/mom11st
      term=((4.d0*ayst)/(15.d0))+(.5d0*((m2-m1)**2.d0)/(m1*m2))
      u2=term-(.08333d0*(2.4d0*best+1.d0)*(m1/m2))
      w=(m1/m2)
      delta=((((6.d0*cest)-5.d0)**2.d0)*w)/(60.d0*(1.d0+u2))
      f=1.d0/(1.d0-delta)
      if(im2.eq.0) write(8,673) f
  673 format(//1x,'f value for second order correction=',1pe11.4,/1x,
     ?'(integrations for second order correction are not',/1x,
     ?'accurate, check if correction becomes significant)')
      if(im2.eq.0) write(8,677) mom12st,mom13st,mom22st,u2,w,delta
  677 format(/1x,'omega*12 =',1pe11.4,2x,'omega*13 =',e11.4,2x,
     ?'omega*22 =',e11.4,/1x,'      u2 =',e11.4,2x,'       w =',
     ?e11.4,2x,'   delta =',e11.4)
      mob=(mconst*f)/(dsqrt(t)*cs)
      write(8,671) mob,1.d0/mob,cs*1.d20
  671 format(//1x,'average (second order) TM mobility =',1pe11.4,
     ?/1x,'inverse average (second order) TM mobility =',e11.4,
     ?/1x,'average TM cross section =',e11.4)
c
      return
      end
c
c     ***************************************************************
c
      subroutine mobil4 (t,inum,inor,ehscs,ehsmob,pacs,pamob,imm)
c
c     Calculates the collision cross-section using the projection 
c     approximation and the exact hard spheres scattering model. Adapted
c     from code written by Alexandre Shvartsburg.
c
      implicit double precision (a-h,m-z)
      dimension cof(-1:100),crof(100)
      parameter (len=1000)
      common/printswitch/ip,it,iu1,iu2,iu3,iv,im2,im4,igs
      common/constants/mu,ro,eo,pi,cang,ro2,dipol,emax,m1,m2,
     ?xe,xk,xn,mconst,correct,romax,inatom,icoord,iic
      common/charge/pcharge(len)
      common/coordinates/fx(len),fy(len),fz(len),
     ?ox(len),oy(len),oz(len)
      common/ljparameters/eolj(len),rolj(len),eox4(len),
     ?ro6lj(len),ro12lj(len),dro6(len),dro12(len)
      common/hsparameters/rhs(len),rhs2(len)
      common/trajectory/sw1,sw2,dtsf1,dtsf2,cmin,ifail,ifailc,inwr
      common/angles/theta,phi,gamma
      common/xrandom/i1,i2,i3,i4,i5,i6
c#######################################
cf2py intent(callback, hide) xrand
      external xrand
c#######################################

c
c     inum is the number of Monte Carlo trajectories employed.
c     inor is the maximum number of successive reflections followed 
c     for each trajectory. cof stores the cosines of halves of angles
c     between the incidence and reflected vectors (for successive 
c     collisions along a single trajectory). The corresponding collision
c     cross sections are accumulated in crof.
c
c     Initialize the cross-section arrays to zero.
      do 1018 im=1,inor
 1018 crof(im)=0.d0
c     crb is the projection (i.e. "hard-sphere cross-section").
      crb=0.d0
c     imm is the highest collision order encountered for any trajectory
      imm=1
c
c     start Monte Carlo integration.
c
      do 1100 i=1,inum
      call rantate
c     determine the y and z extremities to create a box
      ymin=0.d0
      ymax=0.d0
      zmin=0.d0
      zmax=0.d0
      do 1005 in=1,inatom
      y1=fy(in)+rhs(in)
      y2=fy(in)-rhs(in)
      z1=fz(in)+rhs(in)
      z2=fz(in)-rhs(in)
      if(y1.gt.ymax) ymax=y1
      if(y2.lt.ymin) ymin=y2
      if(z1.gt.zmax) zmax=z1
 1005 if(z2.lt.zmin) zmin=z2
c     ydi and zdi are the lengths of the box sides.
      ydi=ymax-ymin
      zdi=zmax-zmin
c     pls is the area of the box.
      pls=ydi*zdi
c     yr and zr are the random coordinates inside the box.
      yr=ymin+ydi*xrand()
      zr=zmin+zdi*xrand()
c     Write the direction cosines of the initial incidence vector 
c     (collinear to the x axis) into the inatom+1 coordinate slot.
      fx(inatom+1)=1.d0
      fy(inatom+1)=0.d0
      fz(inatom+1)=0.d0
      kp=0
      do 1110 im=1,inor  
      call che(im,cof(im),cof(im-1),yr,zr,kp)
c     If the inclusion of the next order of collision had not changed
c     incidence/reflection angle (i.e., no collision occurred) stop 
c     following this trajectory.
      if(cof(im).eq.cof(im-1)) then
c     If necessary, update imm.
      if(im-1.gt.imm) imm=im-1
c     Set all further incidence/reflection angle cosines along this
c     trajectory to the last angle found.
      do 1080 imn=im+1,inor
 1080 cof(imn)=cof(im)
      goto 1085
      else
      if(im.gt.imm) imm=Im
      endif
 1110 continue
c
c     Add the contributions from the i-th trajectory to the collision 
c     cross sections of all orders.
c
 1085 do 1095 im=1,inor
 1095 crof(im)=crof(im)+pls*cof(im)*cof(im)
      if(kp.eq.1) crb=crb+pls
 1100 continue
c
c     End Monte Carlo integration.
c
      u2=0.5d0*dfloat(inum)
c     Normalize all collision cross-sections and the projection.
      do 2110 im=1,imm
 2110 crof(im)=crof(im)/u2
      crb=crb/dfloat(inum)
      pacs=crb
      pamob=mconst/(dsqrt(t)*pacs)
      ehscs=crof(imm)
      ehsmob=mconst/(dsqrt(t)*ehscs)
c
c     Output the results of calculation.
c
      if(im4.eq.0) write(8,600)
  600 format(//1x,'mobility calculation by MOBIL4 (HS scattering)')
      if(im4.eq.0) write(8,601) inum,inor
  601 format(/1x,'number of Monte Carlo trajectories =',i7,
     ?/1x,'maximum number of reflections allowed =',i3)
      write(8,602) pamob,1.d0/pamob,pacs*1.d20
  602 format(/1x,'average PA mobility =',1pe11.4,
     ?/1x,'inverse average PA mobility =',e11.4,
     ?/1x,'average PA cross section =',e11.4)
      write(8,603) imm
  603 format(/1x,'maximum number of reflections encountered =',i3)
      if(im4.eq.0) then
      write(8,606)
  606 format(/5x,'order',3x,'cross section')
      do 2125 im=1,imm
 2125 write(8,604) im,crof(im)
  604 format(5x,i3,5x,1pe12.5)
      endif
      write(8,605) ehsmob,1.d0/ehsmob,ehscs*1.d20
  605 format(/1x,'average EHS mobility =',1pe11.4,
     ?/1x,'inverse average EHS mobility =',e11.4,
     ?/1x,'average EHS cross section =',e11.4)
c
      return
      end
c
c     ***************************************************************
c
      subroutine che (im,cof,cop,yr,zr,kp)
c
c     Guides hard sphere scattering trajectory. Adapted from code 
c     written by Alexandre Shvartsburg.
c 
      implicit double precision (a-h,m-z)
      parameter (len=1000)
      common/printswitch/ip,it,iu1,iu2,iu3,iv,im2,im4,igs
      common/constants/mu,ro,eo,pi,cang,ro2,dipol,emax,m1,m2,
     ?xe,xk,xn,mconst,correct,romax,inatom,icoord,iic
      common/charge/pcharge(len)
      common/coordinates/fx(len),fy(len),fz(len),
     ?ox(len),oy(len),oz(len)
      common/ljparameters/eolj(len),rolj(len),eox4(len),
     ?ro6lj(len),ro12lj(len),dro6(len),dro12(len)
      common/hsparameters/rhs(len),rhs2(len)
      common/trajectory/sw1,sw2,dtsf1,dtsf2,cmin,ifail,ifailc,inwr
      common/angles/theta,phi,gamma
      common/xrandom/i1,i2,i3,i4,i5,i6
c
c     If this is a secondary collision, the incident vector lies on the
c     x axis in transformed coordinates. 
c
      if(im.ne.1) then
      yr=0.d0
      zr=0.d0
      endif
c     xl is a large number greater than any cluster dimension.
      xl=1.0d6
      ki=0
      do 1000 in=1,inatom
      if(fx(in).gt.1.d-16.or.im.eq.1) then
c     yd and zd are the coordinates of the impact points for atom in
c     with respect to its own coordinates (if such a point exists).
c     dev is the impact parameter.   
      yd=yr-fy(in)
      zd=zr-fz(in)
      ras=(yd*yd)+(zd*zd)
      dev=dsqrt(ras)
c     If the collision with in-th atom occurs, then
      if(dev.le.rhs(in)) then
c     Find xc - the x coordinate of collision point with the in-th atom.
      xc=fx(in)-dsqrt(rhs2(in)-ras)
c     Find the smallest xc (earliest collision).
      if(xc.lt.xl) then
      xl=xc
      ki=in
      endif
      endif
      endif
 1000 continue
c     If a collisions took place, then xv, yv, and zv are the vectors
c     going from the center of ki-th atom to the collision point.
      if(ki.ne.0) then
      kp=1
      xv=xl-fx(ki)
      yv=yr-fy(ki)
      zv=zr-fz(ki)
c     Transform all coordinates by a parallel move such that the 
c     collision point is at (0,0,0) in new coordinates.
      do 1020 in=1,inatom
      fx(in)=fx(in)-xl
      fy(in)=fy(in)-yr
 1020 fz(in)=fz(in)-zr
c
c     Transform the coordinates of all atoms such that the direction
c     vector of the reflected particle is collinear to axis x. (Note
c     that the number of such transformations is infinite, thus the
c     direction of the y or z axis is arbitrary.) The direction cosines
c     of the incoming ray are also transformed. xve1, xve2, and xve3 
c     are the direction vectors of reflected ray in the coordinate 
c     system of incoming ray. 
c
c     Evaluate the transformation matrix elements.
      xxv=2.d0*xv*xv
      xyv=2.d0*xv*yv
      xzv=2.d0*xv*zv
      xyz=xyv*xzv
      rad2=rhs2(ki)-xxv
      ad1=(rad2*rad2)+(xyv*xyv)
      adr1=dsqrt(ad1)
      adr2=dsqrt(ad1*ad1+xyz*xyz+xzv*xzv*rad2*rad2)
      xve1=1.d0-2.d0*xv*xv/rhs2(ki)
      xve2=-xyv/rhs2(ki)
      xve3=-xzv/rhs2(ki)
      yve1=xyv/adr1
      yve2=rad2/adr1
      yve3=0.d0
      zve1=rad2*xzv/adr2
      zve2=-xyz/adr2
      zve3=ad1/adr2
  997 format(15x,1pe12.5,1x,e12.5,1x,e12.5)
c     Transform the coordinates and direction cosines of incoming ray.
      do 1010 in=1,inatom+1
      xne=xve1*fx(in)+xve2*fy(in)+xve3*fz(in)
      yne=yve1*fx(in)+yve2*fy(in)+yve3*fz(in)
      fz(in)=zve1*fx(in)+zve2*fy(in)+zve3*fz(in)
      fx(in)=xne
 1010 fy(in)=yne
c     Calculate cof - the cosine of the angle between the incident ray
c     and the normal to an imaginary plane, the reflection from which
c     would be equivalent to the actual multibody reflection. 
      cof=dcos(0.5d0*(pi-dacos(fx(inatom+1))))
      else
c     If no collision occurred, this cosine has not changed.
      cof=cop
      endif
      return
      end
c
c
c     ***************************************************************
c
      subroutine ncoord(unit,dchar,asymp)
c
c     Reads in a new set of coordinates
c
      implicit double precision (a-h,m-z)
      character*30 unit,dchar,dummy
      parameter (len=1000)
      dimension imass(len),xmass(len) 
      common/printswitch/ip,it,iu1,iu2,iu3,iv,im2,im4,igs
      common/constants/mu,ro,eo,pi,cang,ro2,dipol,emax,m1,m2,
     ?xe,xk,xn,mconst,correct,romax,inatom,icoord,iic
      common/charge/pcharge(len)
      common/coordinates/fx(len),fy(len),fz(len),
     ?ox(len),oy(len),oz(len)
      common/ljparameters/eolj(len),rolj(len),eox4(len),
     ?ro6lj(len),ro12lj(len),dro6(len),dro12(len)
      common/hsparameters/rhs(len),rhs2(len)
      common/trajectory/sw1,sw2,dtsf1,dtsf2,cmin,ifail,ifailc,inwr
      common/angles/theta,phi,gamma
      common/xrandom/i1,i2,i3,i4,i5,i6
c
      read(9,'(a30)',end=100) dummy
  100 continue
c
      do 2000 iatom=1,inatom
      if(dchar.eq.'calc') then
      read(9,*) fx(iatom),fy(iatom),fz(iatom),
     ?ximass,pcharge(iatom)
      else
      read(9,*) fx(iatom),fy(iatom),fz(iatom),ximass
      endif
      imass(iatom)=nint(ximass)
      if(unit.eq.'au') then
      fx(iatom)=fx(iatom)*0.52917706d0
      fy(iatom)=fy(iatom)*0.52917706d0
      fz(iatom)=fz(iatom)*0.52917706d0
      endif
 2000 continue
c
      do 2020 iatom=1,inatom
      itest=0
c
c     hydrogen (average value of eo from ab initio calculations
c     and ro from fitting mobilities of C6H6 and others) 
c
      if(imass(iatom).eq.1) then
      itest=1
      xmass(iatom)=1.008d0
      eolj(iatom)=0.65d-03*xe
      rolj(iatom)=2.38d0*1.0d-10
      rhs(iatom)=2.2d0*1.0d-10
      endif
c
c     carbon (from fitting C60 mobility)
c
      if(imass(iatom).eq.12) then
      itest=1
      xmass(iatom)=12.01d0
      eolj(iatom)=1.34d-3*xe
      rolj(iatom)=3.043d0*1.0d-10
      rhs(iatom)=2.7d0*1.0d-10
      endif
c
c     nitrogen (same as carbon)
c
      if(imass(iatom).eq.14) then
      itest=1
      xmass(iatom)=14.01d0
      eolj(iatom)=1.34d-3*xe
      rolj(iatom)=3.043d0*1.0d-10
      rhs(iatom)=2.7d0*1.0d-10
      endif
c     
c     oxygen (same as carbon)
c
      if(imass(iatom).eq.16) then
      itest=1
      xmass(iatom)=16.00d0
      eolj(iatom)=1.34d-3*xe
      rolj(iatom)=3.043d0*1.0d-10
      rhs(iatom)=2.7d0*1.0d-10
      endif
c
c     sodium - Na+ from fitting potential derived by Viehland 
c     (Chem. Phys. 85 (1984) 291) for Na+ + He interactions from mobility
c     data. Note 12-6-4 doesn't provide a very good fit to Viehland's 
c     potential. Viehland's potential is flatter at small r. We didn't
c     include Viehland's first three points in the fit.
c     Hard sphere radius from fitting 300 K Na+ low field mobility in He.
c
      if(imass(iatom).eq.23) then
      itest=1
      xmass(iatom)=22.99d0
      eolj(iatom)=0.0278d-3*xe
      rolj(iatom)=(3.97d0/1.12246d0)*1.0d-10
      rhs(iatom)=2.853d0*1.0d-10
      endif
c
c     silicon (from fitting mobilities of small silicon clusters)
c
      if(imass(iatom).eq.28) then
      itest=1
      xmass(iatom)=28.09d0
      eolj(iatom)=1.35d-3*xe
      rolj(iatom)=3.5d0*1.0d-10
      rhs(iatom)=2.95d0*1.0d-10
      endif
c
c     sulfur (same as silicon)
c
      if(imass(iatom).eq.32) then
      itest=1
      xmass(iatom)=32.06d0
      eolj(iatom)=1.35d-3*xe
      rolj(iatom)=3.5d0*1.0d-10
      rhs(iatom)=3.5d0*1.0d-10
      endif
c
c     iron (same as silicon)
c
      if(imass(iatom).eq.56) then
      itest=1
      xmass(iatom)=55.85d0
      eolj(iatom)=1.35d-3*xe
      rolj(iatom)=3.5d0*1.0d-10
      rhs(iatom)=3.5d0*1.0d-10
      endif
c
      if(itest.eq.0) then
      write(8,602) iatom
  602 format(1x,'type not defined for atom number',i3)
      close (8)
      stop
      endif
c
 2020 continue  
c
      mx=0.d0
      do 2021 iatom=1,inatom
 2021 mx=mx+xmass(iatom)
c
      if(mx.ne.m2) then
      write(8,624)
  624 format(1x,'masses do not add up')
      close (8)
      stop
      endif
c
      do 2030 iatom=1,inatom
      rhs2(iatom)=rhs(iatom)*rhs(iatom)
      eox4(iatom)=4.d0*eolj(iatom)
      ro2lj=rolj(iatom)*rolj(iatom)
      ro6lj(iatom)=ro2lj*ro2lj*ro2lj
      ro12lj(iatom)=ro6lj(iatom)*ro6lj(iatom)
      dro6(iatom)=6.d0*ro6lj(iatom)
 2030 dro12(iatom)=12.d0*ro12lj(iatom)
c
      fxo=0.d0
      fyo=0.d0
      fzo=0.d0
      do 2009 iatom=1,inatom
      fxo=fxo+(fx(iatom)*xmass(iatom))
      fyo=fyo+(fy(iatom)*xmass(iatom))
 2009 fzo=fzo+(fz(iatom)*xmass(iatom))
      fxo=fxo/m2
      fyo=fyo/m2
      fzo=fzo/m2
      do 2010 iatom=1,inatom
      fx(iatom)=(fx(iatom)-fxo)*1.d-10*correct
      fy(iatom)=(fy(iatom)-fyo)*1.d-10*correct
      fz(iatom)=(fz(iatom)-fzo)*1.d-10*correct
 2010 continue
c
      do 3000 iatom=1,inatom
      ox(iatom)=fx(iatom)
      oy(iatom)=fy(iatom)
 3000 oz(iatom)=fz(iatom)
c
c     determine structural asymmetry parameter
c
      theta=0.d0
      asymp=0.d0
      do 5000 igamma=0,360,2
      do 5000 iphi=0,180,2
      gamma=dfloat(igamma)/cang
      phi=dfloat(iphi)/cang
      call rotate
      xyzsum=0.d0
      yzsum=0.d0
      do 5005 iatom=1,inatom
      xyz=dsqrt(fx(iatom)**2+fy(iatom)**2+fz(iatom)**2)
      yz=dsqrt(fy(iatom)**2+fz(iatom)**2)
      xyzsum=xyzsum+xyz
      yzsum=yzsum+yz
 5005 continue
      hold=((pi/4.d0)*xyzsum)/yzsum
      if(hold.gt.asymp) asymp=hold
 5000 continue
c
c
      return
      end
c
c     ***************************************************************
c
