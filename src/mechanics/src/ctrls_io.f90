
      module ctrls_io

      use const_def
      use mechanics_def

      implicit none

      logical, parameter :: ctrls_io_dbg = .false.

      real(dp) :: &
         m1, m2, &
         separation, period, eccentricity, &
         initial_dt, &
         stop_after_this_many_periods, &
         x0, y0, &
         vx0, vy0, &
         x, y, &
         vx, vy, &
         implicit_scheme_tolerance
      
      logical :: verbose
      logical :: do_history

      character (len=strlen) :: history_fname

      namelist /orbit_controls/ &
         m1, &
         m2, &
         separation, &
         period, &
         eccentricity, &
         initial_dt, &
         stop_after_this_many_periods, &
         verbose, &
         do_history, &
         history_fname, &
         implicit_scheme_tolerance

      contains

      subroutine do_one_setup(o, filename, ierr)
         type(orbit_info), pointer :: o
         character (len=*), intent(in) :: filename
         integer, intent(out) :: ierr

         call set_defaults

         call read_controls(o, filename, ierr)
         if (ierr /= 0) then
            write(*,'(a)') 'failed in read_controls'
            return
         end if

      end subroutine do_one_setup


      subroutine read_controls(o, filename, ierr)
         type(orbit_info), pointer :: o
         character (len=*) :: filename
         integer, intent(out) :: ierr
         integer :: iounit

         if (len_trim(filename) > 0) then

            open(newunit=iounit, file=trim(filename), action='read', delim='quote', status='old', iostat=ierr)

            if (ierr /= 0) then
               write(*, *) 'failed to open controls namelist file ', trim(filename)
               return
            end if
            read(iounit, nml=orbit_controls, iostat=ierr)
            close(iounit)
            if (ierr /= 0) then
               write(*, *)
               write(*, *)
               write(*, *)
               write(*, *)
               write(*, '(a)') &
                  'failed while trying to read controls namelist file: ' // trim(filename)
               write(*, '(a)') &
                  'perhaps the following runtime error message will help you find the problem.'
               write(*, *)
               open(newunit=iounit, file=trim(filename), action='read', delim='quote', status='old', iostat=ierr)
               read(iounit, nml=orbit_controls)
               close(iounit)
               return
            end if
         end if
         
         call store_controls(o, ierr)
         if (ierr /= 0) then
            write(*,'(a)') 'failed in store_controls'
            return
         end if
         
      end subroutine read_controls 


      subroutine set_defaults
         m1 = 1d0
         m2 = 8d0
         period = 2d0
         separation = -1d0
         eccentricity = 0d0
         initial_dt = 1d-3
         stop_after_this_many_periods = 1.05
         verbose = .true.
         do_history = .true.
         history_fname = 'orbit_history.data'
         implicit_scheme_tolerance = 1d-6

      end subroutine set_defaults


      subroutine store_controls(o, ierr)
         type(orbit_info), pointer :: o
         integer, intent(out) :: ierr

         ierr = 0

         ! now we can set values on orbit_info
         if (ctrls_io_dbg) write(*,'(a)') 'setting init conditions'
         o% m1 = m1 * Msun
         o% m2 = m2 * Msun
         o% period = period * (24d0*60d0*60d0)
         o% separation = &
            ((standard_cgrav*(o% m1+o% m2))*(o% period/(2*pi))**2)**one_third
         o% eccentricity = eccentricity

         o% initial_dt = initial_dt
         o% stop_after_this_many_periods = stop_after_this_many_periods

         o% implicit_scheme_tolerance = implicit_scheme_tolerance

         o% verbose = verbose
         o% do_history = do_history
         o% history_fname = history_fname

         ! initial position is the periastron
         if (ctrls_io_dbg) write(*,'(a)') 'starting position at periastron'
         o% x0 = 0d0
         o% y0 = o% separation * (1 - o% eccentricity)
         o% vx0 = sqrt((standard_cgrav*o% m2/o% separation) * &
            (1 + o% eccentricity) / (1 - o% eccentricity))
         o% vy0 = 0d0

         ! begin loop at init conditions
         o% x = o% x0
         o% y = o% y0
         o% vx = o% vx0
         o% vy = o% vy0
         
      end subroutine store_controls

      end module ctrls_io 
