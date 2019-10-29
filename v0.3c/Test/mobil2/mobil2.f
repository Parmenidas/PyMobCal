c###################################################
c  MODIFIED
c  
c  * I made xrand() in rantate an external funcion
c    In this way I can use the python numdom number
c    generator
c  * Mobil2 open/closes its own file
c  * We added a common /testing/ to diffeq
c###################################################

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
c#######################################
cf2py intent(callback, hide) xrand
      external xrand
c#######################################
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
c#####################################
      open (8,file='out.fortran.dat')
c#####################################
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
c#####################################
      close(8)
c#####################################

      return
      end
c
