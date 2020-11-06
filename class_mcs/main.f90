
program main

   ! Importing the module
   use mod_class_mcs

   implicit none

   ! Creating an instance of the mcs class
   type(class_mcs) :: spline

   ! Defining a set of data
   integer              :: opt = 2 ! 1=without fictitious points, 2=with fictitious points
   integer              :: N   = 6
   real(8), allocatable :: x(:)
   real(8), allocatable :: y(:)

   ! Getting points for interpolation
   call get_dataset(opt, N, x, y)

   ! Initializing the spline
   call spline%init(x, y)

   ! Print spline
   call print_spline(10, 10000)

contains

   !> \brief Data set
   subroutine get_dataset(opt, N, x, y)
      implicit none
      integer,              intent(in)    :: opt
      integer,              intent(inout) :: N
      real(8), allocatable, intent(out)   :: x(:)
      real(8), allocatable, intent(out)   :: y(:)

      ! Temp. variables
      integer i

      select case (opt)
          case (1) ! Only boundary and internal points
            call get_dataset1(N, x, y)
          case (2) ! Includes fictitious points
            N = N + 2
            call get_dataset2(N, x, y)
          case default

      end select

      do i = 1, N

         write(11,*) x(i), y(i)

      end do

   end subroutine



   !> \brief Data set 1
   subroutine get_dataset1(N, x, y)
      implicit none
      integer,              intent(in)  :: N
      real(8), allocatable, intent(out) :: x(:)
      real(8), allocatable, intent(out) :: y(:)

      ! Temp. variables
      integer i
      real(8) :: xa = 0.d0
      real(8) :: xb = acos(-1.d0)

      ! Allocating memory
      allocate(x(N))
      allocate(y(N))

      do i = 1, N

         x(i) = (xb-xa) * dble(i-1)/dble(N-1) + xa

         y(i) = cos(x(i))

      end do

   end subroutine


   !> \brief Data set 2
   subroutine get_dataset2(N, x, y)
      implicit none
      integer,              intent(in)  :: N
      real(8), allocatable, intent(out) :: x(:)
      real(8), allocatable, intent(out) :: y(:)

      ! Temp. variables
      real(8) :: xb = acos(-1.d0)
      real(8), allocatable :: x1(:)
      real(8), allocatable :: y1(:)

      ! Allocating memory
      allocate(x(N))
      allocate(y(N))
      allocate(x1(N-2))
      allocate(y1(N-2))

      call get_dataset1(N-2,x1,y1)

      x(2:N-1) = x1
      y(2:N-1) = y1

      x(1) = -0.05 * xb
      x(N) =  1.05 * xb

      y(1) =  1.d0
      y(N) = -1.d0

      deallocate( x1, y1)

   end subroutine


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

