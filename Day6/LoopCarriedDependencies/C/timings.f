      real function timings()
c     next three lines needed for hitachi
c     real*8 t
c     call xclock(t,5)
c     timings= t
c     next two lines are required for ncsa super cluster...
c     interface to real*4 function second [c]
c     end
c     real :: tmp(2)   ! needed for the SGI O200 definition.
c     integer :: tmp(4)   ! needed for the hitachi definition.
c------------------------------------------------------------------------------
C     IBM DEFINITION FOR TIMING
c------------------------------------------------------------------------------
c     timings = etime(tmp) ! SGI definition..
c     timings = times(tmp)  ! hitiachi definition..
      timings = second()  ! linux and CRAY definition..
c       timings = mclock()/100. ! IBM definition...
c     may be used someday to generate a system independant timing
c     routine..
c      call system_clock(count=clock1)
c      t = (clock1 - clock0)
c      atime = real( t / hz)
c      print *,atime,t,hz,clock1,clock0
      return
      end
