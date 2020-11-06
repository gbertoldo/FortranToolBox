
program main

   ! Importing the module
   use mod_class_mcs

   implicit none

   ! Temp. variables
   integer i

   ! Creating an instance of the mcs class
   type(class_mcs) :: spline

   ! Defining a set of data
   integer, parameter :: N = 6
   real(8) :: xa = 0.0d0
   real(8) :: xb = 3.141592654d0
   real(8) :: x(N)
   real(8) :: y(N)


   do i = 1, N

      x(i) = (xb-xa) * dble(i-1)/dble(N-1) + xa

      y(i) = cos(x(i))

      write(11,*) x(i), y(i)

   end do

   ! Initializing the spline
   call spline%init(x, y)

   ! Print spline
   call print_spline(10, 10000)

contains

   !> \brief Prints the spline
   subroutine print_spline(unt, np)
      implicit none
      integer, intent(in) :: unt !< Unit to print
      integer, intent(in) :: np  !< Number of partitions

      integer :: i, n
      real(8) :: dx, xl

      n = size(x)

      dx = (x(n)-x(1))/dble(np)

      do i = 0, np

         xl = x(1) + dx * i

         write(unt,*) xl, spline%eval(xl)

      end do

   end subroutine

end program

