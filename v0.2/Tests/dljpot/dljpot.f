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
