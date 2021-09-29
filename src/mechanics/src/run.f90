
      program run

      use const_def
      use mechanics_def

      use mechanics_lib

      implicit none

      integer(8) :: time0, time1, clock_rate
      real(dp) :: runtime

      integer :: id
      integer :: ierr
      type(orbit_info), pointer :: o

      logical, parameter :: dbg = .true.

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

      ! now we can set values on orbit_info
      if (dbg) write(*,'(a)') 'setting init conditions'
      o% m1 = 1d0 * Msun
      o% m2 = 1d1 * Msun
      o% period = 2d0 * (24d0*60d0*60d0)
      o% separation = &
         ((standard_cgrav*(o% m1+o% m2))*(o% period/(2*pi))**2)**one_third
      o% eccentricity = 0d0

      o% initial_dt = 1d-3
      o% stop_after_n_period = 1d0

      o% implicit_scheme_tolerance = 1d-7

      ! initial position is the periastron
      if (dbg) write(*,'(a)') 'starting position at periastron'
      o% x0 = 0d0
      o% y0 = o% separation * (1 - o% eccentricity)
      o% vx0 = sqrt((standard_cgrav*o% m2/o% separation) * &
         (1 + o% eccentricity) / (1 - o% eccentricity))
      o% vy0 = 0d0

      ! begin loop at init conditions
      if (dbg) write(*,'(a)') 'begin loop'
      o% x = o% x0
      o% y = o% y0
      o% vx = o% vx0
      o% vy = o% vy0

      if (dbg) then
         write(*,'(a)')
         write(*,1) 'm1', o% m1 / Msun
         write(*,1) 'm2', o% m2 / Msun
         write(*,1) 'separation', o% separation / Rsun
         write(*,1) 'period', o% period / (24d0*60d0*60d0)
         write(*,1) 'eccentricity', o% eccentricity
         write(*,1) 'x0, y0', o% x0 / Rsun, o% y0 / Rsun
         write(*,1) 'xv0, vy0', o% vx0 / 1d5, o% vy0 / 1d5
         write(*,1) 'x, y', o% x / Rsun, o% y / Rsun
         write(*,1) 'xv, vy', o% vx / 1d5, o% vy / 1d5
         write(*,'(a)')
      end if


      call do_one_orbit(id, ierr)
      if (ierr /= 0) then
         write(*,'(a)') 'failed in do_one_orbit'
         return
      end if


      call system_clock(time1, clock_rate)
      runtime = real(time1 - time0, dp) / clock_rate / 60
      if (dbg) write(*,1) 'runtime (minutes)', runtime
      
      end program run 
