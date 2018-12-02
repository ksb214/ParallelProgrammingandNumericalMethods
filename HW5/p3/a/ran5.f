      real*8 function ran5(idum,v,y)
	implicit none
      integer idum

      integer ia, im, iq, ir, ntab   
      real*8    am, atab
      parameter (ia=16807, im=2147483647, am=1.0/im)
      parameter (iq=127773, ir=2836, ntab=32, atab=ntab-1)
c-
      integer j, k
      real*8 v(ntab), y

c**********************************************************
c***  Contributed by Kursat Bekar, PSU, February 26, 2007
c**********************************************************
c
c  This double-precision fortran function generate a
c  sequence of pseudo-random numbers uniformly distributed
c  over the unit interval.
c
c  Compile with:  ifort -c ran5.f
c  then link the resulting ran5.o with your fortran/C/C++
c  compiled code in the standard way.
c
c  The program calling ran5 must declare as real*8:
c    ran5   the returned random number
c    v(32)  array used inside ran5 intialized as all 0.0
c    y      variable used inside ran5 initialized to 0.5
c
c  The argument idum is the integer seed, set only in first
c  call to ran5 and not changed in between calls.
c
c  An example of declarations in a fortran code:
c
c      integer seed
c      REAL*8  ran5,v(32),u
c      data v/32*0.0/, u/0.5/
c
c  An example of using ran5:
c      seed=1
c      do n=1,nmax
c         x=ran5(seed,v,u)
c         {code section that utilizes x, does not change
c          seed, v, or u}
c      enddo
c
c  For the purposes of achieving good parallel performance
c  you will need to declare all variables appearing in
c  ran5 private. This will ensure correct performance
c  since different threads will not overwrite intermediate
c  variabes in ran5, and will reduce contention due to
c  reading/writing to the same locations in the shared
c  memory. An easy way to do this is to declare a
c  Default(private) in the section where ran5 is called,
c  then explicitly declare as Shared those variables that
c  need be shared across threads. None of the variables
c  in ran5 need be Shared. Word of caution: By declaring
c  a private default, any variables you want the children
c  to know from the parent thread must be passed via a
c  Firstprivate clause.
c
c**********************************************************

      if (idum.le.0) then
          idum = max(-idum,1)
          do 12 j=ntab,1,-1
              k = idum/iq
              idum = ia*(idum-k*iq) - ir*k
              if (idum.lt.0) idum = idum+im
              v(j) = am*idum
   12     continue
          y = v(1)
      end if
    1 continue
      k = idum/iq
      idum = ia*(idum-k*iq) - ir*k
      if (idum.lt.0) idum = idum+im
      j = 1 + int(atab*y)
      y = v(j)
      ran5 = y
      v(j) = am*idum
      if (ran5.eq.0.0 .or. ran5.eq.1.0) goto 1

         return
         end
