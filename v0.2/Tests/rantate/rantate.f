c     ****************************************************************
c     MODIFIED: Here, We open/close the file within rotate.
c     The rest is unchanged 
c     *****************************************************************

!       double precision function xrand()
!       implicit double precision (a-h,o-z)
!       dimension rvec(10)
!       common/xrandom/i1,i2,i3,i4,i5,i6
! c
! c     XRAND is a random number generator that uses RANLUX if i5=1
! c     otherwise it uses the standard RAND subroutine available in 
! c     FORTRAN 77. If RANLUX is used i1 contains the luxury level
! c     (1-4, 4 is the highest, 3 is default in RANLUX). i2, i3, and
! c     i4 are seed integers. i3 and i4 will normally be zero. If the
! c     standard RAND subroutine is to be employed, i2 contains the 
! c     seed integer. RANLUX was downloaded from http://kumo.swcp.
! c     com/fortran/random2.f90. 
! c 
!       i6=i6+1
! c
!       if (i5.eq.1) then
!       call ranlux(rvec,1)
!       xrand=rvec(1)
!       return
!       else
!       xrand=rand()
!       return
!       endif
! c
!       end
! c
c     ***************************************************************


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
      open (8,file='out.fortran.dat')
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
      close (8)
 4000 return
      end
