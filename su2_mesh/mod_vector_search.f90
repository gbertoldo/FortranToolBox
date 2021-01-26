!> \brief Provides methods to search elements in vectors
module mod_vector_search
    implicit none

contains

    !>
    !! \brief Applies an algorithm similar to the binary search
    !! to find the interval within element e is in vector vec.
    !! If e is out of range, this function returns -1.
    !!
    !! Caution: vector vec must be sorted!
    !!
    integer function binary_search(e, vec)
        implicit none
        real(8),               intent(in) :: e
        real(8), dimension(:), intent(in) :: vec

        integer :: N
        integer :: kb, km, ke
        real(8) :: eb, em, ee

        N = size(vec)

        kb = 1
        ke = N
        eb = vec(kb)
        ee = vec(ke)

        ! Checking if the element is out of range
        if ( e < eb .or. e > ee ) then
            binary_search = -1
        end if

        do while ( ke-kb > 1 )

            km = (kb+ke)/2
            em = vec(km)

            if ( e <= em ) then
                ke = km
                ee = em
            else
                kb = km
                eb = em
            end if

        end do

        binary_search = kb

    end function


    !> \brief Provides a verification of the binary_search function
    subroutine binary_search_verification()
        implicit none
        integer, parameter :: N = 100000
        real(8) :: x(N)
        logical :: failure = .false.
        real(8) :: e
        integer :: k
        integer :: i

        do i = 1, N
            x(i) = i**0.5d0
        end do

        write(*,*) "Starting binary_search verification..."

        do i = 1, N

            e = x(i)

            k = binary_search(e, x)

            if ( .not. ( x(k) <= e .and. e <= x(k+1) ) ) then
                write(*,*) "Fail at index ", i
                failure = .true.
            end if

        end do

        do i = 1, N-1

            e = 0.5d0 * ( x(i) + x(i+1) )

            k = binary_search(e, x)

            if ( .not. ( x(k) <= e .and. e <= x(k+1) ) ) then
                write(*,*) "Fail at index ", i
                failure = .true.
            end if

        end do

        if ( .not. failure ) write(*,*) "Verification test concluded with success."

    end subroutine

end module
