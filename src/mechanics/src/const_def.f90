
      module const_def

      implicit none

      ! float precision
      integer, parameter :: dp = selected_real_kind(p=15)

      integer, parameter :: strlen = 256 ! for character (len=strlen)
      
      ! logging
      integer, parameter :: keep_going = 0
      integer, parameter :: terminate = 1

      ! math
      real(dp), parameter :: pi = 3.1415926535897932384626433832795028841971693993751d0
      real(dp), parameter :: one_third = 1d0 / 3d0

      real(dp), parameter :: secyer = 24*60*60*365.25d0 ! secs per year

      ! astro
      real(dp), parameter :: standard_cgrav = 6.67430d-8 ! gravitational constant (g^-1 cm^3 s^-2)
      real(dp), parameter :: mu_sun = 1.3271244d26
      real(dp), parameter :: Msun = mu_sun / standard_cgrav ! solar mass (g)
      real(dp), parameter :: Rsun = 6.957d10 ! solar radius (cm)

      ! debug
      logical, parameter :: const_def_dbg = .true.

      contains

      subroutine do_const_init(ierr)
         integer, intent(out) :: ierr

         ierr = keep_going

         if (const_def_dbg) write(*,'(a)') 'init const_def module'
         
      end subroutine do_const_init

      end module const_def
