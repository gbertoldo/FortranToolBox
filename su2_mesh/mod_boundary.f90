module mod_boundary
    implicit none
    
    type, public :: boundary
        integer :: sz
        real(8), allocatable, dimension(:) :: x
        real(8), allocatable, dimension(:) :: y       
    contains 
        procedure, public, pass :: init
        
    end type
    
contains

    subroutine init(this, sz)
        implicit none
        class(boundary) :: this
        integer, intent(in) :: sz
        
        this%sz = sz
        
        if ( allocated(this%x) ) deallocate(this%x)
        if ( allocated(this%y) ) deallocate(this%y)
        
        allocate( this%x(0:sz) )
        allocate( this%y(0:sz) )
        
        this%x = 0.d0
        this%y = 0.d0
        
    end subroutine        


    subroutine print_boundary(bound)
        implicit none
        type(boundary), intent(in) :: bound 
        
        integer :: i
        
        do i = 0, bound%sz
        
            write(*,"(I10,2X,ES23.16,2X,ES23.16)") i, bound%x(i), bound%y(i) 
        
        end do
        
    end subroutine
end module
