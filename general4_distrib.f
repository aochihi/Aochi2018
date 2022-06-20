c program for depth-variaion 30/05/2017 
c reference: Aochi (2018). 
c Dynamic asymmetry of normal and reverse faults 
c due to constrained depth-dependent stress accumulation, 
c Geophys. J. Int., 215, 2134-3243, doi:10.1093/ggy407.
c
c This is the program giving the intiial condition for a normal fault
c according to the formulation of Aochi (2018).
c
c given "p" is maximum principal stress "s_1"
c fs, fd, and tvalue fixed. 
c 
c NOTE:  THE PARAMETERS MODIFIABLES ARE GIVEN IN TEST PROGRAM. 
c Update: 2020/06/20 (H. Aochi)
c 
	SUBROUTINE general4n(s1, fs, fd, cohs, tval, phi, dsg, tb) 
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
	dsg = cohs*tval + tval*(fs - fd)*s1 + fd*s1
	dsg2 = (tval*(fs-fd)+fd)*(1.+cos(phi2)) + sin(phi2)

	dsg = dsg/dsg2

c Mohr's circle
        s3 = s1 - 2.*dsg
        smean =  (s1 + s3)/2.

        if (s3.lt.0.0) then
          s3 = 0.0d0
          dsg = (s1 - s3)/2.
          smean = (s1 + s3)/2.
        endif
        
c initial stress field for optimal plane
	sh0 = dsg*sin(phi2)
	tn0 = smean - dsg*cos(phi2)

c peak strength
	tp = cohs + tn0*fs

c residual strength
	tr = tn0*fd
	tb = tp - tr

	tdum = (sh0-tr)/tb
c	write(*,'(7f8.2)') s1, s3, tp, tr, sh0, tn0, tdum

	RETURN
	END
c
c test program
c
	IMPLICIT NONE
        REAL(8) :: p, fs, tb, excess, phi, dsg, fd, pi
        REAL(8) :: tvalue, cohesive, tp, tr
        REAL(8) :: dip, tau0, sn0, s1, s3, zdep, g, rho
        CHARACTER*60 :: outfile
	INTEGER :: i, iz
ccc
ccc PRAMETERS
ccc
	rho = 2700. - 1000.

        fs = 0.60d0
	fd = 0.48d0
	tvalue = 1.00d0
	cohesive = 5.0d0

        dip = 45.d0
        outfile = "test.out"
ccc        
ccc DO NOTE MODIFY BELOW
ccc
	pi = 3.1415920d0
	g = 9.81
        dip = pi*(90.0-dip)/180.
        open(42, file=outfile)

	do iz=0, 24
	  zdep = iz*0.5
	  s1 = rho*g*zdep*10.0**(-3)

	  call general4n(s1, fs, fd, cohesive, tvalue, phi, dsg, tb) 

	  s3 = s1 - 2.*dsg 
c added 2017/12/14 avoiding negative s3 -> in subroutine
          if (s3.lt.0.0d0) then
            s3 = 0.0d0
            dsg = (s1 - s3)/2.
          endif
          tau0 = dsg*sin(2.d0*dip)
          sn0  = s1 - dsg*(1.0 + cos(2.d0*dip))

	  tp = cohesive + fs*sn0
          tr = fd*sn0

	  write(42,'(8f10.3)') 
     &	    zdep, s1, s3, fd, tau0, tp, tr, sn0

	enddo
c repeating for figure
c        s1 = rho*g*5.0*10.0**(-3)
c	  call general4(s1, fs, fd, cohesive, tvalue, phi, dsg, tb) 
c	  s3 = s1 - 2.*dsg 
c        write(42, '(a5,2f10.3)') "p", (s1+s3)/2., 0.0
c        write(42, '(a5,2f10.3)') "s1", s1, 0.0
c        write(42, '(a5,2f10.3)') "s3", s3, 0.0
c        write(42, '(a5,2f10.3)') "fs", 0.0, cohesive+fs*0.0
c        write(42, '(a5,2f10.3)') "fs", 100.0, cohesive+fs*100.0
c        write(42, '(a5,2f10.3)') "fs", 200.0, cohesive+fs*200.0
c        write(42, '(a5,2f10.3)') "fd", 0.0, fd*0.0
c        write(42, '(a5,2f10.3)') "fd", 100.0, fd*100.0
c        write(42, '(a5,2f10.3)') "fd", 200.0, fd*200.0
c        do i=0,180
c          tau0 = dsg*sin(pi*real(i)/180.) 
c          sn0 = (s1+s3)/2. - dsg*cos(pi*real(i)/180.)
c          write(42, '(a5,2f10.3)') "mohr", sn0, tau0
c        enddo
      
c        tau0 = dsg*sin(2.d0*phi) 
c        sn0 = (s1+s3)/2. - dsg*cos(2.0*phi)
c        write(42, '(a5,2f10.3)') "phi", phi/pi*180.
c        write(42, '(a5,2f10.3)') "tpop", sn0, cohesive+sn0*fs
c        write(42, '(a5,2f10.3)') "t0op", sn0, tau0
c        write(42, '(a5,2f10.3)') "trop", sn0, sn0*fd
      
c        tau0 = dsg*sin(2.d0*dip) 
c        sn0 = (s1+s3)/2. - dsg*cos(2.0*dip)
c        write(42, '(a5,2f10.3)') "dip", (pi/2.-dip)/pi*180.
c        write(42, '(a5,2f10.3)') "tp", sn0, cohesive+sn0*fs
c        write(42, '(a5,2f10.3)') "t0", sn0, tau0
c        write(42, '(a5,2f10.3)') "tr", sn0, sn0*fd

        close(42)
        END

