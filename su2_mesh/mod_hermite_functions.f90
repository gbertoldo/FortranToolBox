!> \brief Hermite functions
module mod_hermite_functions

contains

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

end module
