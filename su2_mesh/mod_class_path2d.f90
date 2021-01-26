!> \brief Defines a class for handling 2D paths

module mod_class_path2d

    use mod_vector_search

    implicit none

    ! Makes everything private, except otherwise stated
    private

    ! Public procedures
    public :: class_path2d_verification

    !> \brief Defines a 2D cartesian point
    type, public :: class_point2d
        real(8) :: x !< x coordinate
        real(8) :: y !< y coordinate
    end type

    !>
    !! \brief Defines a class to parameterize a path (x,y) as a function of t,
    !!        where t is dimensionless length of the path, e.g., t = integrate(ds, 0, smax)/smax.
    !!
    type, public :: class_path2d
        real(8), dimension(:), allocatable :: tv     !< Vector or t
        real(8), dimension(:), allocatable :: xv     !< Vector or x
        real(8), dimension(:), allocatable :: yv     !< Vector of y
        integer                            :: sz     !< Number of points of the discrete path
        real(8)                            :: length !< Returns the curve length

    contains

        procedure, public, pass :: init       !< Initializes the path
        procedure, public, pass :: point      !< Returns the point (x,y) for a given t
        procedure, public, pass :: x          !< Returns x for a given t
        procedure, public, pass :: y          !< Returns y for a given t
        final                   :: destructor !< Class destructor

    end type

contains

    !> \brief Class initializer
    subroutine init(this, x, y)
        implicit none
        class(class_path2d)               :: this !< A reference to this object
        real(8), dimension(:), intent(in) ::    x !< x coordinates of the discrete path
        real(8), dimension(:), intent(in) ::    y !< y coordinates of the discrete path

        ! Inner variables
        integer :: i

        ! x and y must have the same size
        if ( size(x) /= size(y) ) then
            write(*,*) "class_path2d%init: x and y must have the same size. Stopping..."
            stop
            return
        end if

        ! If memory was allocated before, deallocate ite
        if ( allocated(this%tv) ) deallocate(this%tv)
        if ( allocated(this%xv) ) deallocate(this%xv)
        if ( allocated(this%yv) ) deallocate(this%yv)

        ! Allocating memory
        allocate( this%tv(size(x)) )
        allocate( this%xv(size(x)) )
        allocate( this%yv(size(x)) )

        ! Initializing data
        this%tv = 0.d0
        this%xv = x
        this%yv = y

        ! Calculating the parameter t
        this%tv(1) = 0.d0

        do i = 2, size(x)

            this%tv(i) = this%tv(i-1) + sqrt( (x(i)-x(i-1))**2 + (y(i)-y(i-1))**2 )

        end do

        this%length = this%tv(size(x))

        ! Normalizing t
        this%tv = this%tv / this%length

        ! Defining the number of points of the discrete path
        this%sz = size(x)

    end subroutine



    !> \brief Calculates the point (x,y) for a given t
    type(class_point2d) function point(this, t)
        implicit none
        class(class_path2d) :: this !< A reference to this object
        real(8), intent(in) :: t    !< Parameter t [0,1]
        type(class_point2d) :: p    !< (x(t), y(t))

        integer :: k

        k = binary_search(t, this%tv)

        if ( k < 0 ) then
            write(*,*) "class_path2d%point: t must be in the interval [0,1]"
            stop
        end if

        p%x = this%xv(k) + (this%xv(k+1)-this%xv(k))/(this%tv(k+1)-this%tv(k))*(t-this%tv(k))
        p%y = this%yv(k) + (this%yv(k+1)-this%yv(k))/(this%tv(k+1)-this%tv(k))*(t-this%tv(k))

        point = p

     end function


    !> \brief Calculates the coordinate x for a given t
    real(8) function x(this, t)
        implicit none
        class(class_path2d) :: this !< A reference to this object
        real(8), intent(in) :: t    !< Parameter t [0,1]

        integer :: k

        k = binary_search(t, this%tv)

        if ( k < 0 ) then
            write(*,*) "class_path2d%x: t must be in the interval [0,1]"
            stop
        end if

        x = this%xv(k) + (this%xv(k+1)-this%xv(k))/(this%tv(k+1)-this%tv(k))*(t-this%tv(k))

     end function



    !> \brief Calculates the coordinate y for a given t
    real(8) function y(this, t)
        implicit none
        class(class_path2d) :: this !< A reference to this object
        real(8), intent(in) :: t    !< Parameter t [0,1]

        integer :: k

        k = binary_search(t, this%tv)

        if ( k < 0 ) then
            write(*,*) "class_path2d%y: t must be in the interval [0,1]"
            stop
        end if

        y = this%yv(k) + (this%yv(k+1)-this%yv(k))/(this%tv(k+1)-this%tv(k))*(t-this%tv(k))

    end function


    !> \brief Destructor
    subroutine destructor(this)
        implicit none
        type(class_path2d) :: this !< A reference to this object

        if ( allocated(this%tv) ) deallocate(this%tv)
        if ( allocated(this%xv) ) deallocate(this%xv)
        if ( allocated(this%yv) ) deallocate(this%yv)

    end subroutine


    !> \brief Verification subroutine
    !! TODO: Finish this procedure
    subroutine class_path2d_verification()
        implicit none

        integer, parameter :: N = 80
        real(8) :: x(N)
        real(8) :: y(N)
        type(class_path2d) :: path

        real(8) :: t
        integer :: i

        ! Creating the coordinates for the path = (x,y=x**2)
        do i = 1, N
            x(i) = dble(i-1)/dble(N-1)
            y(i) = x(i)**2
        end do

        ! Initializing the path
        call path%init(x, y)

        do i = 1, N

            t = dble(i-1)/dble(N-1)

            write(*,*) path%x(t), path%y(t)-path%x(t)**2

        end do

    end subroutine

end module
