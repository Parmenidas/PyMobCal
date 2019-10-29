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
c     "FAKE" subroutine that computes harmonic potential
c     instead of lennard-jones
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
      k=1
      dmax=0.d0

      pot=-0.5*k*(x**2+y**2+z**2)
      dpotx=k*x
      dpoty=k*y
      dpotz=k*z     
      return
      end
