
      program run

      use const_def
      use ctrls_io 
      use mechanics_def
      use mechanics_lib

      implicit none

      integer(8) :: time0, time1, clock_rate
      real(dp) :: runtime

      integer :: id
      integer :: ierr
      type(orbit_info), pointer :: o

      logical, parameter :: dbg = .false.

      include 'include/formats.inc'

      ierr = 0

      call system_clock(time0, clock_rate)
      

      id = alloc_orbit(ierr)
      if (ierr /= 0) then
         write(*,'(a)') 'failed in alloc_binary'
         return
      end if

      call orbit_ptr(id, o, ierr)
      if (ierr /= 0) then
         write(*,'(a)') 'failed in orbit_ptr'
         return
      end if

      call do_one_setup(o, 'orbit_nml', ierr)
      if (ierr /= 0) then
         write(*,'(a)') 'failed in do_one_setup'
         stop
      end if
      
      write(*,'(a)')
      write(*,1) 'm1', o% m1 / Msun
      write(*,1) 'm2', o% m2 / Msun
      write(*,1) 'separation', o% separation / Rsun
      write(*,1) 'period', o% period / (24d0*60d0*60d0)
      write(*,1) 'eccentricity', o% eccentricity
      write(*,'(a)')

      call do_full_orbit(id, ierr)
      if (ierr /= 0) then
         write(*,'(a)') 'failed in do_full_orbit'
         stop
      end if

      call system_clock(time1, clock_rate)
      runtime = real(time1 - time0, dp) / clock_rate / 60
      if (dbg) write(*,1) 'runtime (minutes)', runtime
      
      end program run 
