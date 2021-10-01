
      module mechanics_def

      use const_def

      implicit none
      
      type orbit_info

         integer :: id

         logical :: in_use

         real(dp) :: &
            ! masses of binary system
            m1, m2, &
            ! separation, period & eccentricity
            separation, period, eccentricity, &
            ! initial timestep
            initial_dt, &
            ! fraction of periods to compute orbit
            stop_after_this_many_periods, &
            ! initial position of m2, assume m1 fixed in CM
            x0, y0, &
            ! initial velocity of m2
            vx0, vy0, &
            ! position of m2, assume m1 fixed in CM
            x, y, &
            ! velocity of m2
            vx, vy
         
         ! verbosity during evolution loop
         logical :: verbose

         ! flag to write output to file
         logical :: do_history

         ! fname of history file
         character (len=strlen) :: history_fname

         ! tolerance for Runge-Kutta solver
         real(dp) :: implicit_scheme_tolerance

      end type orbit_info

      logical :: have_initialized_orbit_handles = .false.
      integer, parameter :: max_orbit_handles = 10
      type(orbit_info), target, save :: orbit_handles(max_orbit_handles)

      logical, parameter :: mechanics_def_dbg = .false.

      contains
     

      integer function alloc_orbit(ierr)
         integer, intent(out) :: ierr
         integer :: i
         type (orbit_info), pointer :: o

         ierr = 0
         alloc_orbit = -1
!$omp critical (orbit_handle)
         if (.not. have_initialized_orbit_handles) then
            do i = 1, max_orbit_handles
               orbit_handles(i)% id = i
               orbit_handles(i)% in_use = .false.
            end do
            have_initialized_orbit_handles = .true.
         end if
         do i = 1, max_orbit_handles
            if (.not. orbit_handles(i)% in_use) then
               orbit_handles(i)% in_use = .true.
               alloc_orbit = i
               exit
            end if
         end do
!$omp end critical (orbit_handle)
         if (alloc_orbit == -1) then
            ierr = -1
            return
         end if
         if (orbit_handles(alloc_orbit)% id /= alloc_orbit) then
            ierr = -1
            return
         end if
         o => orbit_handles(alloc_orbit)

      end function alloc_orbit


      subroutine orbit_ptr(id, o, ierr)
         integer :: id
         type(orbit_info), pointer, intent(inout) :: o
         integer, intent(out) :: ierr

         if (mechanics_def_dbg) write(*,'(a)') 'setting orbit_ptr'

         if (id < 1 .or. id > max_orbit_handles) then
            ierr = -1
            return
         end if

         o => orbit_handles(id)
         ierr = 0

      end subroutine orbit_ptr


      subroutine free_orbit(o)
         type (orbit_info), pointer :: o

         orbit_handles(o% id)% in_use = .false.

      end subroutine free_orbit

      end module mechanics_def 
