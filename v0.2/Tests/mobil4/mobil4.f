
c###################################################
c  MODIFIED
c  
c  * I made xrand() in rantate an external funcion
c    In this way I can use the python numdom number
c    generator
c  * Mobil4 open/closes its own file
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
c#####################################
      open (8,file='out.fortran.dat')
c#####################################

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
c#####################################
      close(8)
c#####################################

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
c     *******************************************************
