
      module mechanics_lib

      use const_def
      use mechanics_def
      use num_def
      use num_lib

      implicit none

      logical, parameter :: mechanics_lib_dbg = .false.


      contains
      

      ! RHS of ODE system for the two body problem
      subroutine two_body_derivs(n, x, h, y, f, lrpar, rpar, lipar, ipar, ierr)
         integer, intent(in) :: n, lrpar, lipar
         real(dp), intent(in) :: x,h
         real(dp), intent(inout) :: y(:) ! (n)
         real(dp), intent(inout) :: f(:) ! (n)
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         integer, intent(out) :: ierr ! nonzero means retry with smaller timestep.
         real(dp) :: mass

         include 'include/formats.inc'

         mass = rpar(3)

         ierr = 0
         f(1) = y(3)
         f(2) = y(4)
         f(3) = - standard_cgrav * mass * y(1) / (y(1)**2 + y(2)**2)**(3d0/2d0)
         f(4) = - standard_cgrav * mass * y(2) / (y(1)**2 + y(2)**2)**(3d0/2d0)
         
      end subroutine two_body_derivs
      

      subroutine write_header_orbit_history(id, io, ierr)
         integer, intent(in) :: id
         integer, intent(in) :: io
         integer, intent(out) :: ierr

         type(orbit_info), pointer :: o
         character (len=strlen) :: fname

         ierr = 0
         call orbit_ptr(id, o, ierr)
         if (ierr /= 0) then
            write(*,'(a)') 'failed in orbit_ptr'
            return
         end if
         
         if (mechanics_lib_dbg) write(*,'(a)') 'writting header to file'

         fname = trim(o% history_fname)
         !inquire(file=fname, exist=history_file_exists)
         open(unit=io, file=trim(fname), action='write', iostat=ierr)
         if (ierr /= 0) then
            write(*,*) 'failed to open ' // trim(fname)
            return
         end if

         write(io, fmt='(a32,1x)', advance='no') trim('step_number')
         write(io, fmt='(a32,1x)', advance='no') trim('time')
         write(io, fmt='(a32,1x)', advance='no') trim('timestep')
         write(io, fmt='(a32,1x)', advance='no') trim('X')
         write(io, fmt='(a32,1x)', advance='no') trim('Y')
         write(io, fmt='(a32,1x)', advance='no') trim('vx')
         write(io, fmt='(a32,1x)', advance='yes') trim('vY')

      end subroutine write_header_orbit_history


      subroutine write_header_in_terminal(io)
         integer, intent(in) :: io

         write(io,'(a)') &
            '_______________________________________________________________________' // &
            '____________________'
         write(io,*)
         write(io,'(a)') &
            'step_number         time     timestep            X            Y        ' // &
            '   vX           vY'
         write(io,'(a)') &
            '_______________________________________________________________________' // &
            '____________________'
         write(io,*)

      end subroutine write_header_in_terminal


      ! print accepted solution during Runge-Kutta step loop
      subroutine solout(nr,xold,x,n,y,rwork,iwork,interp_y,lrpar,rpar,lipar,ipar,irtrn)
         ! nr is the step number.
         ! x is the current x value; xold is the previous x value.
         ! y is the current y value.
         ! irtrn negative means terminate integration.
         ! rwork and iwork hold info for
         integer, intent(in) :: nr, n, lrpar, lipar
         real(dp), intent(in) :: xold, x
         real(dp), intent(inout) :: y(:) ! (n)
         real(dp), intent(inout), target :: rwork(*)
         integer, intent(inout), target :: iwork(*)
         integer, intent(inout), pointer :: ipar(:) ! (lipar)
         real(dp), intent(inout), pointer :: rpar(:) ! (lrpar)
         interface
            ! this subroutine can be called from your solout routine.
            ! it computes interpolated values for y components during the just completed step.
            real(dp) function interp_y(i,s,rwork,iwork,ierr)
               use const_def, only: dp
               integer, intent(in) :: i ! result is interpolated approximation of y(i) at x=s.
               real(dp), intent(in) :: s ! interpolation x value (between xold and x).
               real(dp), intent(inout), target :: rwork(*)
               integer, intent(inout), target :: iwork(*)
               integer, intent(out) :: ierr
            end function interp_y
         end interface
         integer, intent(out) :: irtrn
         character (len=90) :: fmt, fmt1, fmt2, fmt3

         irtrn = 1
         
         ! dump to file
         if (ipar(2) >  0) then
            write(ipar(2), fmt='(i32,1x)', advance='no') nr  ! step number
            write(ipar(2), fmt='(1pes32.16e3, 1x)', advance='no') x / (24d0*60d0*60d0)  ! time
            write(ipar(2), fmt='(1pes32.16e3, 1x)', advance='no') (x-xold) / (24d0*60d0*60d0)  ! timestep
            write(ipar(2), fmt='(1pes32.16e3, 1x)', advance='no') y(1) / Rsun  ! X
            write(ipar(2), fmt='(1pes32.16e3, 1x)', advance='no') y(2) / Rsun  ! Y
            write(ipar(2), fmt='(1pes32.16e3, 1x)', advance='no') y(3) / 1d5   ! vX
            write(ipar(2), fmt='(1pes32.16e3, 1x)', advance='yes') y(4) / 1d5   ! vY
         end if

         ! no screen output
         if (ipar(1) /= 1) return

         if (mod(nr-1, 10) == 0) call write_header_in_terminal(6)


         ! format for nr and time
         if (x/(24d0*60d0*60d0) > 1d2) then
            fmt1 = '(3x,i8,2x,1pe11.3,0p,'
         else
            fmt1 = '(3x,i8,2x,f11.6,'
         end if
         ! format for timestep
         if ((x-xold)/(24d0*60d0*60d0) > 1d2) then
            fmt2 = '2x,1pe11.3,0p,'
         else
            fmt2 = '2x,f11.6,'
         end if
         ! format for X, Y, vX, vY
         fmt3 = '4(2X,1pe11.3,0p))'

         fmt = trim(fmt1) // trim(fmt2) // trim(fmt3)
         write(6,fmt=fmt) &
            nr-1, &
            x/(24d0*60d0*60d0), &
            (x-xold)/(24d0*60d0*60d0), &
            y(1)/Rsun, &
            y(2)/Rsun, &
            y(3)/Rsun, &
            y(4)/Rsun
         
      end subroutine solout


      subroutine solve_two_body_problem(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr

         integer, parameter :: nv = 4  ! the number of variables in the two body system of ODEs
         real(dp), parameter :: eps = 1d-3 ! stiffness parameter for van der Pol
         real(dp) :: rtol(1) ! relative error tolerance(s)
         real(dp) :: atol(1) ! absolute error tolerance(s)
         real(dp) :: x ! starting value for the interval of integration
         real(dp) :: xend ! ending value for the interval of integration
         integer, parameter :: lrpar = 3, lipar = 2
         real(dp) :: h, max_step_size
         integer :: lout, iout, idid, itol, j
         integer :: liwork, lwork, max_steps
         real(dp), pointer :: work(:)
         integer, pointer :: iwork(:)
         real(dp), target :: y_ary(nv)
         real(dp), pointer :: y(:)
         real(dp), target :: rpar_ary(lrpar)
         integer, target :: ipar_ary(lipar)
         integer, pointer :: ipar(:) ! (lipar)
         real(dp), pointer :: rpar(:) ! (lrpar)
         type(orbit_info), pointer :: o
         integer :: iounit

         include 'include/formats.inc'
         
         ierr = 0
         call orbit_ptr(id, o, ierr)
         if (ierr /= 0) then
            write(*,'(a)') 'failed in orbit_ptr'
            return
         end if

         ipar => ipar_ary
         rpar => rpar_ary

         if (mechanics_lib_dbg) write(*,'(a)') 'two body problem with the cash_karp method'

         y => y_ary

         ! time range of applicability is between [x, xend]
         x = 0
         xend = o% period * o% stop_after_this_many_periods


         ! initial conditions
         y(1) = o% x0
         y(2) = o% y0
         y(3) = o% vx0
         y(4) = o% vy0

         lout = 6
         max_steps = 10000
         max_step_size = o% period / 200
         itol = 0 ! scalar tolerances
         iout = 1 ! use 0 for no intermediate output
         rtol(1) = 1d-4
         atol(1) = 1d-4
         h = o% initial_dt
         rpar(1) = eps
         rpar(2) = 0
         rpar(3) = o% m2

         if (o% verbose) then
            ipar(1) = 1
         else
            ipar(1) = 0
         end if

         ! > 0 for output into a file with unit
         ipar(2) = 0
         if (o% do_history) then
            iounit = 22
            ipar(2) = iounit
            call write_header_orbit_history(id, iounit, ierr)
            if (ierr /= 0) return
         end if

         call cash_karp_work_sizes(nv,liwork,lwork)
         allocate(work(lwork), iwork(liwork))

         iwork = 0
         work = 0

         call cash_karp( &
               nv,two_body_derivs,x,y,xend, &
               h,max_step_size,max_steps, &
               rtol,atol,itol, &
               solout,iout,work,lwork,iwork,liwork, &
               lrpar,rpar,lipar,ipar,lout,idid)

         if (idid /= 1) then ! trouble
            write(*,*) 'idid', idid
            ierr = - 1
            return
         end if

         close(iounit)

         deallocate(work, iwork)
      
         write(*,'(a)')
         write(*,1) &
            'reached t limit: t_limit, t_final', &
            xend/(24d0*60d0*60d0), x/(24d0*60d0*60d0)
         write(*,'(a)')

      end subroutine solve_two_body_problem


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

         call solve_two_body_problem(id, ierr)
         if (ierr /= 0) then
            write(*,'(a)') 'failed in solve_two_body_problem'
            return
         end if
         
      end subroutine do_full_orbit 

      end module mechanics_lib
