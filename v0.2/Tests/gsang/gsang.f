c     ****************************************************************
c     MODIFIED: Here, We open/close the file within the function.
c     Also, we replaced a check for ns>30000 with ns>1000 to speed up
c     the test
c     We added a common /testing/ to diffeq
c     The rest is unchanged 
c     *****************************************************************

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
      open (8,file='out.fortran.out')
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
      if(ns.gt.1000) then
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
      if(ifailc.eq.ifail) then 
      close(8)
      stop
      endif
      return
c       
      close(8)
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
