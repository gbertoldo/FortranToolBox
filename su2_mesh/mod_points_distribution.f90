!> \brief Provides procedures to distribute points along a line

module mod_points_distribution

    use mod_class_mcs

    implicit none

contains

    !> \brief Calculates a normalized uniform distribution t
    !! such that t(1)=0, t(N)=1.
    subroutine get_uniform_distribution(N, t) ! Output: last one
        implicit none
        integer, intent(in)  :: N          !< Number of points in the interval
        real(8), intent(out) :: t(1:N)     !< Distribution

        integer :: j ! Dummy index

        t(1) = 0.d0

        do j = 2, N-1

            t(j) = (dble(j-1)/dble(N-1))

        end do

        t(N) = 1.d0

    end subroutine


    !> \brief Calculates a normalized geometric progression distribution t
    !! t(1)=0, t(N)=1. Points are concentrated near t=0 by default, except
    !! if isReversed is true.
    subroutine get_gp_distribution(isReversed, N, r, t) ! Output: last one
        implicit none
        logical, intent(in)  :: isReversed !< Reverses the discretization distribution if true
        integer, intent(in)  :: N          !< Number of points in the interval
        real(8), intent(in)  :: r          !< Ratio of length to first partition width
        real(8), intent(out) :: t(1:N)     !< Distribution

        integer :: j ! Dummy index
        real(8) :: q ! GP ratio

        call get_GP_ratio(N-1,r,q) ! Last one is output

        t(1) = 0.d0

        do j = 2, N-1

            t(j) = t(j-1) + q**(j-2) / r

        end do

        t(N) = 1.d0

        ! Should be reversed?
        if (isReversed) call reverse(t)

    end subroutine


    !> \brief Calculates a normalized Monotonic Cubic Spline (MCS) distribution t
    !! t(1)=0, t(N)=1. Points are concentrated near t=0 and near t=1.
    !!   fa = fraction of the a uniform partition for left side
    !!   fb = fraction of the a uniform partition for left side
    !!   Example: consider a uniform partition of partition width 0.1.
    !!            If fa = 0.2, then the first partition on the left will
    !!            be fa * 0.1 = 0.02.
    subroutine get_mcs_distribution(N, fa, fb, t) ! Output: last one
        implicit none
        integer, intent(in)  :: N          !< Number of points in the entire interval
        real(8), intent(in)  :: fa         !< Fraction of the a uniform partition for left side
        real(8), intent(in)  :: fb         !< Fraction of the a uniform partition for left side
        real(8), intent(out) :: t(1:N)     !< Distribution

        integer :: i           ! Dummy index
        real(8) :: zv(4)       ! Values of the index
        real(8) :: tv(4)       ! Values of t
        type(class_mcs) :: mcs ! Creating a monotonic cubic spline interpolator

        ! Creating the points to be interpolated
        zv(1) = 0.d0
        zv(2) = (1.d0)/dble(N-1)
        zv(3) = dble(N-2)/dble(N-1)
        zv(4) = 1.d0

        tv(1) = 0.d0
        tv(2) = 1.d0 / dble(N) * fa
        tv(3) = 1.d0- 1.d0 / dble(N) * fb
        tv(4) = 1.d0

        ! Initializing the interpolator
        call mcs%init(zv, tv)

        ! Evaluating the spline
        do i = 1, N

            t(i) = mcs%eval(dble(i-1)/dble(N-1))

        end do

    end subroutine


    !> \brief Calculates a normalized geometric power law distribution t
    !! t(1)=0, t(N)=1. Points are concentrated near t=0 by default, except
    !! if isReversed is true.
    subroutine get_power_law_distribution(isReversed, N, r, t) ! Output: last one
        implicit none
        logical, intent(in)  :: isReversed !< Reverses the discretization distribution if true
        integer, intent(in)  :: N          !< Number of points in the interval
        real(8), intent(in)  :: r          !< Ratio of length to first partition width
        real(8), intent(out) :: t(1:N)     !< Distribution

        integer :: j ! Dummy index
        real(8) :: a ! Exponent of the power law

        if ( r > 1.d0 ) then

            a = log( r ) / log(dble(N-1))

        else

            a = 1.d0

        end if

        t(1) = 0.d0

        do j = 2, N-1

            t(j) = (dble(j-1)/dble(N-1))**a

        end do

        t(N) = 1.d0

        ! Should be reversed?
        if (isReversed) call reverse(t)

    end subroutine


    !> \brief Reverses a partition distribution
    subroutine reverse(t)
        implicit none
        real(8), intent(inout) :: t(:) !< Distribution of points

        integer :: i, N
        real(8) :: dt(size(t)-1)

        ! Number of partitions
        N = size(t)-1

        ! Copying the partition distribution in reversed order
        do i = 1, N
            dt(N-i+1) = t(i+1)-t(i)
        end do

        ! Reassigning t
        do i = 2, N
            t(i)=t(i-1)+dt(i-1)
        end do

    end subroutine


    !> \brief Calculates the geometric progression ratio
    subroutine get_gp_ratio(n, r, q)
        implicit none
        integer, intent(in)  ::   n !< number of partitions
        real(8), intent(in)  ::   r !< l/a1
        real(8), intent(out) ::   q !< q

        ! Parameters

        integer :: nit = 1000   ! Maximum number of iteractions
        real(8) :: tol = 1.d-15 ! Tolerance

        ! Inner variables

        integer ::   i ! Dummy index
        real(8) ::  qi ! inital value of q
        real(8) ::  qf ! final value of q
        real(8) ::  qm ! mean value of q

        if ( r < n ) then

            qi = 0.1d0

            qf = 1.d0 - 1.d-15

        else

            qi = 1.d0 + 1.d-15

            qf = 10.d0

        end if

        do i = 1, nit

            qm = 0.5d0 * qi + 0.5d0 * qf

            if ( 0.d0 < f(qi) * f(qm) ) then

                qi = qm

            else

                qf = qm

            end if

            if ( abs(qf-qi) < tol ) exit

        end do


        if ( i == nit ) then

            write(*,*) "get_gp_ratio: Maximum number of iteractions was exceeded."

            stop

        end if

        q = qm

    contains

        real(8) function f(q)
            implicit none
            real(8), intent(in) :: q

            f = q ** n + r * ( 1.d0 - q ) - 1.d0

        end function f

    end subroutine get_gp_ratio


end module
