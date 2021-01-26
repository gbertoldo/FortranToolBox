!
!    Copyright (C) 2020 by Guilherme Bertoldo
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
!! \brief Class class_mcs construct a monotonic cubic spline from a set of
!!        data (x,y). After constructed, y(x) may be calculated within [xa, xb],
!!        where xa and xb are the lower and upper limits of the data set used to
!!        construct the spline.
!!        For technical details, see:
!!          - https://en.wikipedia.org/wiki/Monotone_cubic_interpolation
!!          - Fritsch, F. N.; Carlson, R. E. (1980). "Monotone Piecewise Cubic Interpolation".
!!            SIAM Journal on Numerical Analysis. SIAM. 17 (2): 238–246. doi:10.1137/0717021
!!
module mod_class_mcs
   implicit none

   ! Makes everything private, except otherwise stated
   private

   !> \brief Class class_mcs for monotonic cubic spline construction
   type, public :: class_mcs
      integer              :: n     !< Number of points
      real(8), allocatable :: m(:)  !< Derivatives dy/dx
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
      class(class_mcs)              :: this !< A reference to this object
      real(8),           intent(in) :: x(:) !< Set of independent coordinates
      real(8),           intent(in) :: y(:) !< Set of dependent coordinates

      ! Seting the number of points
      this%n = size(x)

      if (allocated(this%xv)) then

         deallocate(this%xv)
         deallocate(this%yv)
         deallocate(this%m)

      end if

      ! Saving the points
      allocate(this%xv(this%n))
      allocate(this%yv(this%n))
      allocate(this%m(this%n))

      this%xv = x
      this%yv = y

      ! Calculates the spline derivatives in order to ensure monotonicity
      call calc_spline_derivatives(x, y, this%m)

   end subroutine


   !> \brief Calculates the spline derivatives in order to ensure monotonicity
   subroutine calc_spline_derivatives(x, y, m)
      implicit none
      real(8),           intent(in)  :: x(:) !< Set of independent coordinates
      real(8),           intent(in)  :: y(:) !< Set of dependent coordinates
      real(8),           intent(out) :: m(:) !< Derivatives dy/dx

      ! Parameters
      real(8), parameter :: eps = 10.d0 * epsilon(1.d0)

      ! Inner variables
      integer :: i
      integer :: n
      real(8), allocatable :: alpha(:)
      real(8), allocatable :: beta (:)
      real(8), allocatable :: delta(:)

      ! Seting the number of points
      n = size(x)

      ! Allocating memory
      allocate(alpha(n-1))
      allocate(beta (n-1))
      allocate(delta(n-1))


      ! --------------------------------------------------
      ! STEP 1
      !
      ! Calculating the secants (delta)

      do i = 1, n-1

         delta(i)=(y(i+1)-y(i))/(x(i+1)-x(i))

      end do

      ! --------------------------------------------------
      ! STEP 2
      !
      ! Calculating the provisional tangents (m)

      m(1) = delta(1)

      do i = 2, n-1

         if (delta(i)*delta(i-1)<0.d0) then

            m(i) = 0.d0

         else

            m(i) = (delta(i)+delta(i-1))/2

         end if

      end do

      m(n) = delta(n-1)


      ! --------------------------------------------------
      ! STEP 3
      !
      ! If two successive y are equal (delta=0), set m=0 on
      ! the two boundaries of the partition
      do i = 1, n-1

         if ( abs(delta(i)) <= eps ) then

            m(i)   = 0.d0

            m(i+1) = 0.d0

         end if

      end do

      ! --------------------------------------------------
      ! STEP 4
      !
      ! Calculating alpha and beta (only for partitions where |delta| > 0)

      alpha = 0.d0
      beta  = 0.d0

      do i = 1, n-1

         if ( abs(delta(i)) > eps ) then

            alpha(i) = m(i) / delta(i)

            beta(i)  = m(i+1) / delta(i)

            if ( alpha(i) < 0.d0 ) then

               m(i) = 0.d0

               alpha(i) = 0.d0

            end if

            if ( beta(i) < 0.d0 ) then

               m(i+1) = 0.d0

               beta(i) = 0.d0

            end if

         end if

      end do


      ! --------------------------------------------------
      ! STEP 5
      !
      ! Adjusting monotonicity (only for partitions where |delta| > 0)

      do i = 1, n-1

         if ( abs(delta(i)) > eps ) then

            if ( alpha(i)**2+beta(i)**2 > 9.d0 ) then

               m(i)   = alpha(i) * delta(i) * 3.d0 / sqrt(alpha(i)**2+beta(i)**2)

               m(i+1) = beta(i)  * delta(i) * 3.d0 / sqrt(alpha(i)**2+beta(i)**2)

               !print*, "-->", atan(m(i))*180/3.14, atan(m(i+1))*180/3.14

            end if

         end if

         !print*, i, atan(m(i))*180/3.14, atan(m(i+1))*180/3.14, atan(delta(i+1))*180/3.14

      end do

      ! Cleaning memory
      deallocate(alpha)
      deallocate(beta)
      deallocate(delta)

   end subroutine


   !> \brief Hermite function H00
   real(8) function h00(t)
      implicit none
      real(8), intent(in) :: t

      h00 = (1.d0+2*t)*(1.d0-t)**2

   end function


   !> \brief Hermite function H01
   real(8) function h01(t)
      implicit none
      real(8), intent(in) :: t

      h01 = t**2*(3.d0-2*t)

   end function

   !> \brief Hermite function H10
   real(8) function h10(t)
      implicit none
      real(8), intent(in) :: t

      h10 = t*(1.d0-t)**2

   end function

   !> \brief Hermite function H11
   real(8) function h11(t)
      implicit none
      real(8), intent(in) :: t

      h11 = t**2*(t-1.d0)

   end function


   !> \brief Evaluates y(x)
   real(8) function eval(this, x)
      implicit none
      class(class_mcs)    :: this !< A reference to this object
      real(8), intent(in) :: x    !< Independent variable

      ! Inner variables
      integer :: i
      real(8) :: h
      real(8) :: t
      real(8) :: xlower
      real(8) :: ylower
      real(8) :: mlower
      real(8) :: xupper
      real(8) :: yupper
      real(8) :: mupper

      associate (          &
            n  => this%n,  &
            xv => this%xv, &
            yv => this%yv, &
            m  => this%m   &
            )

         ! Searching for an interval that contains x
         do i = 2, n

            if ( x <= xv(i) ) exit

         end do

         xlower = xv(i-1)
         ylower = yv(i-1)
         mlower = m(i-1)
         xupper = xv(i)
         yupper = yv(i)
         mupper = m(i)

         h = xupper - xlower

         t = (x-xlower) / h

         eval = ylower * h00(t) + h * mlower * h10(t) + yupper * h01(t) + h * mupper * h11(t)

      end associate

   end function

end module
