
      module mechanics_lib

      use const_def
      use mechanics_def

      implicit none

      logical, parameter :: mechanics_lib_dbg = .true.


      contains

      integer function rhs(x, y, vx, vy, ax, ay, m, ierr)
         real(dp), intent(inout) :: x, y, vx, vy
         real(dp) :: ax, ay
         real(dp), intent(in) :: m
         integer, intent(out) :: ierr

         real(dp) :: r

         ierr =  0

         rhs = keep_going

         if (mechanics_lib_dbg) write(*,'(a)') 'computing rhs'

         r = sqrt(x*x + y*y)

         ax = - standard_cgrav * m * x / (r*r*r)
         ay = - standard_cgrav * m * y / (r*r*r)

      end function rhs


      subroutine do_one_orbit(dt, ierr)

         
      end subroutine do_one_orbit 


      subroutine do_full_orbit(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type(orbit_info), pointer :: o
         real(dp) :: rel_error
         real(dp) :: step_dt, dt

         ierr = 0
         call orbit_ptr(id, o, ierr)
         if (ierr /= 0) then
            write(*,'(a)') 'failed in orbit_ptr'
            return
         end if

         ! begin with initial dt
         dt = o% initial_dt
         o% dt = dt

         evolve_loop: do while (o% t < o% period * o% stop_after_n_periods)

            if (mechanics_lib_dbg) write(*,'(a)') 'doing evolve loop'

            rel_error = 1d99

            step_loop: do while (rel_err > o% implicit_scheme_tolerance)

               step_dt = dt

               ! avoid going beyond max period
               if (o% t + step_dt > o% period * o% stop_after_n_periods) &
                  step_dt = o% period * o% stop_after_n_periods - o% t


            end do step_loop

            ! accept solution, update parameters
            o% dt = step_dt
            o% t = o% t + o% dt
            
         end do evolve_loop  
         
      end subroutine do_full_orbit 

      end module mechanics_lib
