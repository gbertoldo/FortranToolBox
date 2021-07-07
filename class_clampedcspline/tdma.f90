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
!! \brief TriDiagonal Matrix Algorithm (TDMA) solves tridiagonal linear systems a x = b.
!!

module mod_tdma
  implicit none

contains

   subroutine tdma(n, a, b, x)
      implicit none
      integer, intent(in) :: n ! Number unknowns
      real(8), dimension(n,3), intent(in)  :: a ! Tri-diagonal matrix
      real(8), dimension(n),   intent(in)  :: b ! Source
      real(8), dimension(n),   intent(out) :: x ! Solution
      !
      ! a(i,1) = west coefficients
      ! a(i,2) = central coefficients
      ! a(i,3) = east coefficients
      !
      ! Auxiliary variables
      integer :: i
      real(8), dimension(n) :: P
      real(8), dimension(n) :: Q

      i = 1

      P(i) = - a(i,3) / a(i,2)

      Q(i) = b(i) / a(i,2)

      do i = 2, n

       P(i) = - a(i,3) / ( a(i,2) + a(i,1) * P(i-1) )

       Q(i) = ( b(i) - a(i,1) * Q(i-1) ) / ( a(i,2) + a(i,1) * P(i-1) )

      end do

      i = n
      x(i) = Q(i)

      do i = n-1, 1, -1
       x(i) = x(i+1) * P(i) + Q(i)
      end do

   end subroutine tdma

end module
