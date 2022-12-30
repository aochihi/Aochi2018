c program for strain-constrained stress model for strike-slip fault 
c (30/12/2022)
c reference: Aochi and Tsuda (GJI, 2022)
c See and read "layer5rev_distrib.f"
c
c This program provides stress field for striek-slip fault presented 
c in Figure 4b.
c
       IMPLICIT NONE
       INTEGER nn, nmax, izmax
       PARAMETER (nn=10, izmax=150) 
       INTEGER:: iz, i, nlayer
       REAL(8):: dep(nn), vs(nn), rho(nn), rid(nn)
       REAL(8):: depsmin(nn)
       REAL(8):: zdep, sz, deps, dz, rh2o
       REAL(8):: pi, g, fs, fd , tvalue, cohesive, dsg, phi, tb
       REAL(8):: s1, tau0, sn0, tp, tr, s3
	CHARACTER*60:: infile

c model parameters
       infile = "1d.txt"
       fs = 0.60d0
       fd = 0.45d0
       tvalue = 1.0d0
       cohesive = 5.0d0

       pi = 3.1415920d0 
       g = 9.81d0
       rh2o = 1000.d0
c dz in km
       dz = 0.1

       open(51, file=infile, status="old", err=99) 
       read(51, *) nmax
       if (nmax.gt.nn) then
 	 write(*,*) "too many layers. To increase nn", nmax, nn
	 stop
       endif
       write(*,*) "Top Depth (km): Vs (m/s): Rigidity (GPa)"
       do i=1, nmax
          read(51,*) dep(i), vs(i), rho(i)
          rid(i) = vs(i)*vs(i)*rho(i)
	  write(*, '(f12.2, f12.0,f12.1)') dep(i), vs(i), rid(i)/1.0e9
          depsmin(i) = 1.0d0
       enddo
       close(51)

c "before constrained" model
       open(53, file="ss_before_constrained.dat")
       do iz=0, izmax
         zdep = iz*dz
         i = nmax
         nlayer = -1
         do i = nmax, 1, -1
           if(nlayer.lt.0.and.zdep.gt.dep(i)) then
             nlayer = i
           endif
         enddo
         if (nlayer.lt.0) nlayer = 0

         sz = 0.0
         do i = nlayer, 1, -1
           if (i.eq.nlayer) then
             sz = sz + (zdep-dep(i))*rho(i)
           else
             sz = sz + (dep(i+1)-dep(i))*rho(i)
           endif
         enddo
         sz = (sz-zdep*rh2o)*g*0.001

         call general4ss(sz,fs,fd,cohesive,tvalue,phi,dsg,tb)

         s1 = sz + dsg
         s3 = sz - dsg
         if (s3.lt.0.0) then
            s3 = 0.0
            call general4r(s3, fs, fd, cohesive, tvalue, phi, dsg, tb)
            s1 = s3 + 2.*dsg
         endif

         deps = dsg/rid(nlayer)
         if (deps.lt.depsmin(nlayer)) depsmin(nlayer) = deps
         write(53, '(f8.2,4f10.3, 2e11.3)') zdep, s3, sz, s1, 
     &          dsg, rid(nlayer), deps
       enddo
       close(53)

       do i=nmax-1, 1, -1
         if(depsmin(i).gt.depsmin(i+1)) then
           depsmin(i) = depsmin(i+1)
         endif
       enddo
c 2nd turn
       open(391, file="ss_constrained.dat")
       write(391, '(a1, 4f12.4)') ">", fs, fd, tvalue, cohesive
       write(391, '(a1,i5,f12.4)') ">", izmax, dz
       do iz=0, izmax
         zdep = iz*dz
         i = nmax
         nlayer = -1
         do i = nmax, 1, -1
           if(nlayer.lt.0.and.zdep.gt.dep(i)) then
             nlayer = i
           endif
         enddo
         if (nlayer.lt.0) nlayer = 0

         sz = 0.0
         do i = nlayer, 1, -1
           if (i.eq.nlayer) then
             sz = sz + (zdep-dep(i))*rho(i)
           else
             sz = sz + (dep(i+1)-dep(i))*rho(i)
           endif
         enddo
         sz = (sz-zdep*rh2o)*g*0.001

         call general4ss(sz,fs,fd,cohesive,tvalue,phi,dsg,tb)
         s1 = sz + dsg
         s3 = sz - dsg
         if (s3.lt.0.0) then
            s3 = 0.0
            call general4r(s3, fs, fd, cohesive, tvalue, phi, dsg, tb)
            s1 = s3 + 2.*dsg
         endif

         deps = dsg/rid(nlayer)

         if(nlayer.ne.nmax.and.deps.gt.depsmin(nlayer+1)) then
            deps = depsmin(nlayer+1)
            dsg = rid(nlayer)*deps
         endif

c         s1 = sz + 2.*dsg
         s1 = sz + dsg
         s3 = sz - dsg
         if (s3.lt.0.0) then
c NEW DESCRIPTION
           s3 = 0.0d0
           s1 = s3 + 2.*dsg
         endif
         if (deps.lt.depsmin(nlayer)) depsmin(nlayer) = deps
         write(391, '(f8.2,4f10.3, 2e11.3)') zdep, s3, sz, s1, 
     &          dsg, rid(nlayer), deps
       enddo

 99    continue
       stop
       END
c program for reverse faulting 12/12/2017, created from general4.f
c given "p" is minimum principal stress "s_3"
c fs, fd, and tvalue fixed. 
c
	SUBROUTINE general4ss(s2, fs, fd, cohs, tval, phi, dsg, tb) 
	IMPLICIT NONE
	REAL(8) :: p, fs, tb, excess, phi, dsg, fd
	REAL(8) :: pi, phi2, s1, s2, s3, smean, theta, shdum, tndum, 
     &          sh0, tn0, tp, tr, cohs, tval, tdum, dsg0, dsg1, theta0
	
	pi = 3.1415920d0
	theta = atan(fs)
        theta0 = atan(fd)

c reliable direction
	phi = 0.25*pi - 0.5*theta
        phi2 = 2.*phi

c principal streses
c        s1 = p

c deviate stress : SS
        dsg0 = s2*sin(theta0)
        dsg1 = (cohs/tan(theta) + s2)*sin(theta)
	dsg = dsg0 + tval*(dsg1 - dsg0)


c Mohr's circle
        s1 = s2 + dsg
        s3 = s2 - dsg
        if ( s3.lt.0.0 ) then
          s3 = 0.0d0
          dsg = (s1 - s3)/2.
        endif
        smean = (s1 + s3)/2.
        
c initial stress field for optimal plane
	sh0 = dsg*sin(phi2)
	tn0 = smean - dsg*cos(phi2)

c peak strength
	tp = cohs + tn0*fs

c residual strength
	tr = tn0*fd
	tb = tp - tr

	tdum = (sh0-tr)/tb
c	write(*,'(8f8.2)') s1, s2, s3, tp, tr, sh0, tn0, tdum

	RETURN
	END
c
        SUBROUTINE general4r(s3, fs, fd, cohs, tval, phi, dsg, tb)
        IMPLICIT NONE
        REAL(8) :: p, fs, tb, excess, phi, dsg, fd
        REAL(8) :: pi, phi2, s1, s3, smean, theta, shdum, tndum,
     &          sh0, tn0, tp, tr, cohs, tval, tdum, dsg2

        pi = 3.1415920d0
        theta = atan(fs)

c reliable direction
        phi = 0.25*pi - 0.5*theta
        phi2 = 2.*phi

c principal streses
c        s1 = p

c deviate stress
        dsg = cohs*tval + (tval*(fs - fd) + fd)*s3
        dsg2 = sin(phi2) - (tval*(fs-fd) + fd)*(1.-cos(phi2))

        dsg = dsg/dsg2

c Mohr's circle
        s1 = s3 + 2.*dsg
        smean =  (s1 + s3)/2.

c initial stress field for optimal plane
        sh0 = dsg*sin(phi2)
        tn0 = smean - dsg*cos(phi2)

c peak strength
        tp = cohs + tn0*fs

c residual strength
        tr = tn0*fd
        tb = tp - tr

        tdum = (sh0-tr)/tb
c       write(*,'(7f8.2)') s1, s3, tp, tr, sh0, tn0, tdum

        RETURN
        END

