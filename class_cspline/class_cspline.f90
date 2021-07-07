!
!    Copyright (C) 2021 by Guilherme Bertoldo
!
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!    Contact:
!
!          Guilherme Bertoldo
!                 E-mail: glbertoldo@gmail.com

!    Institution
!          Federal University of Technology - Paraná - UTFPR
!          Linha Santa Bárbara, s/n, Francisco Beltrão, Paraná, Brazil
!          Zip Code 85601-970
!

!>
!! \brief Class class_cspline construct a natural cubic spline from a set of
!!        data (x,y). After constructed, y(x) may be calculated within [xa, xb],
!!        where xa and xb are the lower and upper limits of the data set used to
!!        construct the spline.
!!        For technical details, see:
!!          - Faires, J. D. & Burden, R. Thomson/Brooks/Cole (Ed.) Numerical methods.
!!
module mod_class_cspline
   use mod_tdma

   implicit none

   ! Makes everything private, except otherwise stated
   private

   !> \brief Class class_cspline for natural cubic spline construction
   !!
   !!        y(x) = yi + b (x-xi) + c (x-xi)^2 + d (x-xi)^3, for xi <= x <= xi+1
   !!
   type, public :: class_cspline
      integer              :: n     !< Number of points
      real(8), allocatable :: b(:)  !< Coefficients
      real(8), allocatable :: c(:)  !< Coefficients
      real(8), allocatable :: d(:)  !< Coefficients
      real(8), allocatable :: xv(:) !< Vector of discrete values of x
      real(8), allocatable :: yv(:) !< Vector of discrete values of y

   contains
      procedure, public, pass :: init !< Constructor
      procedure, public, pass :: eval !< Evaluates y(x)
   end type

contains

   !> \brief Initializes spline object
   subroutine init(this, x, y)
      implicit none
      class(class_cspline)              :: this !< A reference to this object
      real(8),               intent(in) :: x(:) !< Set of independent coordinates
      real(8),               intent(in) :: y(:) !< Set of dependent coordinates

      ! Seting the number of points
      this%n = size(x)

      if (allocated(this%xv)) then

         deallocate(this%b)
         deallocate(this%c)
         deallocate(this%d)
         deallocate(this%xv)
         deallocate(this%yv)

      end if

      ! Saving the points
      allocate(this%b(this%n) )
      allocate(this%c(this%n) )
      allocate(this%d(this%n) )
      allocate(this%xv(this%n))
      allocate(this%yv(this%n))

      this%xv = x
      this%yv = y

      ! Calculates the coefficients of the spline
      call get_csplines_coeff(this%n, x, y, this%b, this%c, this%d) ! Output: last three

   end subroutine


   !> \brief Calculates the natural cubic splines coefficients
   subroutine get_csplines_coeff(n, x, y, b, c, d) ! Output: last three
      implicit none
      integer, intent(in)  :: n !< Number of intervals
      real(8), intent(in)  :: x(0:n) !< x-coordinates
      real(8), intent(in)  :: y(0:n) !< y-coordinates
      real(8), intent(out) :: b(0:n) !< b coefficients
      real(8), intent(out) :: c(0:n) !< c coefficients
      real(8), intent(out) :: d(0:n) !< d coefficients

      ! Inner variables

      integer :: i
      real(8) :: h(0:n-1)
      real(8) :: bp(0:n)
      real(8) :: ap(0:n,3)


      ! Calculating h

      do i = 0, n-1

         h(i) = x(i+1)-x(i)

      end do


      ! Setting west boundary condition

      i = 0

      ap(i,1) = 0.d0
      ap(i,2) = 1.d0
      ap(i,3) = 0.d0
      bp(i) = 0.d0

      do i = 1, n-1

         ap(i,1) = h(i-1)
         ap(i,2) = 2.d0 * ( h(i) + h(i-1) )
         ap(i,3) = h(i)
         bp(i) = 3.d0 * ( (y(i+1)-y(i)) / h(i) - (y(i)-y(i-1)) / h(i-1) )

      end do

      ! Setting east boundary condition

      i = n

      ap(i,1) = 0.d0
      ap(i,2) = 1.d0
      ap(i,3) = 0.d0
      bp(i) = 0.d0

      ! Solving the linear system

      call tdma(n+1, ap, bp, c)

      ! Calculating the other coefficients

      do i = 0, n-1

         b(i) = ( y(i+1)-y(i) ) / h(i) - h(i) * ( c(i+1) + 2.d0 * c(i) ) / 3.d0

         d(i) = ( c(i+1)-c(i) ) / ( 3.d0 * h(i) )

      end do

      b(n) = b(n-1) + h(n-1) * ( c(n) + c(n-1) )

      d(n) = 0.d0

   end subroutine get_csplines_coeff


   !> \brief Evaluates y(x)
   real(8) function eval(this, x)
      implicit none
      class(class_cspline) :: this !< A reference to this object
      real(8),  intent(in) :: x    !< Independent variable

      ! Inner variables
      integer :: i

      associate (           &
            n  => this%n,   &
            xv => this%xv,  &
            yv => this%yv,  &
            b  => this%b,   &
            c  => this%c,   &
            d  => this%d    &
            )

         ! Searching for an interval that contains x
         do i = 2, n

            if ( x <= xv(i) ) exit

         end do

         i = i-1

         eval =  yv(i)                     &
            + b(i) * ( x - xv(i) )         &
            + c(i) * ( x - xv(i) ) ** 2.d0 &
            + d(i) * ( x - xv(i) ) ** 3.d0

      end associate

   end function

end module
