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
