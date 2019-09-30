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

      filen1='a10A1.mfj'
      open (8,file='a10A1_f.out')
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
      ip=0
      it=0
      iu1=1
      iu2=0
      iu3=0
      iv=0
      im2=0
      im4=0
      igs=0
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

      close(8)
      stop
      end
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
