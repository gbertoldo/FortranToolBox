!> \brief This program creates a Gmsh script that is used to generate a mesh for SU2 aerodynamics shape optimization.
program main
    use mod_boundary
    use mod_class_ifile
    use mod_class_path2d
    use mod_points_distribution
    implicit none

    character(len=100) :: version = "su2_mesh - 1.0.0. Last update 26/01/2021. Guilherme Bertoldo."

    ! Parameters to generate the mesh
    type params
        integer :: km
        integer :: kib
        integer :: nx
        integer :: ny
        real(8) :: lr
        real(8) :: rb
        real(8) :: n
        real(8) :: lb
        real(8) :: lf
        real(8) :: nib
        integer :: ko
        character(len=200) :: ogfile
        integer :: kpog
        real(8) :: faog
        real(8) :: fbog
        integer :: kpin
        real(8) :: fain
        real(8) :: fbin
    end type

    type(params) :: par

    write(*,*)
    write(*,*) version
    write(*,*)

    call load_params()

    call gmsh_script()

contains

    !> \brief Loads the parameters from the input file
    subroutine load_params()
        implicit none

        type(class_ifile) :: ifile

        call ifile%init("input.txt", "&")

        call ifile%load()

        call ifile%get_value(par%kib, "kib")
        call ifile%get_value( par%km,  "km")
        call ifile%get_value( par%nx,  "nx")
        call ifile%get_value( par%ny,  "ny")
        call ifile%get_value( par%lr,  "lr")
        call ifile%get_value( par%rb,  "rb")
        call ifile%get_value(  par%n,   "n")
        call ifile%get_value( par%lb,  "lb")
        call ifile%get_value( par%lf,  "lf")
        call ifile%get_value(par%nib, "nib")
        call ifile%get_value(par%ko, "ko")
        call ifile%get_value(par%ogfile, "ogfile")
        call ifile%get_value(par%kpog, "kpog")
        call ifile%get_value(par%faog, "faog")
        call ifile%get_value(par%fbog, "fbog")
        call ifile%get_value(par%kpin, "kpin")
        call ifile%get_value(par%fain, "fain")
        call ifile%get_value(par%fbin, "fbin")

    end subroutine


    !> \brief Generates the gmsh script
    subroutine gmsh_script()
        implicit none

        integer :: i, ib, ie
        integer :: Pcount = 0
        integer :: Lcount = 0

        type(boundary) :: inlet
        type(boundary) :: outlet
        type(boundary) :: symmetry
        type(boundary) :: ogive

        call inlet%init(   par%nx)
        call outlet%init(  par%ny)
        call symmetry%init(par%ny)
        call ogive%init(   par%nx)

        call generate_ogive(ogive)
        call generate_inlet(inlet)
        call generate_outlet(outlet)
        call generate_symmetry(symmetry)

        open(10, file="mesh.geo")

        ! Printing the points of the ogive
        write(10,*) "// Ogive points"
        do i = 0, ogive%sz
            Pcount = Pcount + 1
            write(10,"(A,I4,A,ES23.16,A,ES23.16,A)") "Point(", Pcount, ") = {", ogive%x(i), ", ", ogive%y(i), ", 0.0, 1.0};"
        end do

        ! Printing the points of the outlet
        write(10,*) "// Outlet points"
        do i = 1, outlet%sz-1
            Pcount = Pcount + 1
            write(10,"(A,I4,A,ES23.16,A,ES23.16,A)") "Point(", Pcount, ") = {", outlet%x(i), ", ", outlet%y(i), ", 0.0, 1.0};"
        end do

        ! Printing the points of the inlet
        write(10,*) "// Inlet points"
        do i = inlet%sz, 0, -1
            Pcount = Pcount + 1
            write(10,"(A,I4,A,ES23.16,A,ES23.16,A)") "Point(", Pcount, ") = {", inlet%x(i), ", ", inlet%y(i), ", 0.0, 1.0};"
        end do

        ! Printing the points of the symmetry
        write(10,*) "// Symmetry points"
        do i = 1, symmetry%sz-1
            Pcount = Pcount + 1
            write(10,"(A,I4,A,ES23.16,A,ES23.16,A)") "Point(", Pcount, ") = {", symmetry%x(i), ", ", symmetry%y(i), ", 0.0, 1.0};"
        end do

        ! Printing the lines of the ogive
        write(10,*) "// Ogive lines"
        do i = 0, ogive%sz-1
            Lcount = Lcount + 1
            write(10,"(A,I4,A,I4,A,I4,A)") "Line(", Lcount, ") = {",Lcount,", ", Lcount+1, "};"
        end do

        ! Printing the lines of the outlet
        write(10,*) "// Outlet lines"
        do i = 0, outlet%sz-1
            Lcount = Lcount + 1
            write(10,"(A,I4,A,I4,A,I4,A)") "Line(", Lcount, ") = {",Lcount,", ", Lcount+1, "};"
        end do

        ! Printing the lines of the inlet
        write(10,*) "// Inlet lines"
        do i = 0, inlet%sz-1
            Lcount = Lcount + 1
            write(10,"(A,I4,A,I4,A,I4,A)") "Line(", Lcount, ") = {",Lcount,", ", Lcount+1, "};"
        end do

        ! Printing the lines of the symmetry
        write(10,*) "// Symmetry lines"
        do i = 0, symmetry%sz-2
            Lcount = Lcount + 1
            write(10,"(A,I4,A,I4,A,I4,A)") "Line(", Lcount, ") = {",Lcount,", ", Lcount+1, "};"
        end do
        Lcount = Lcount + 1
        write(10,"(A,I4,A,I4,A,I4,A)") "Line(", Lcount, ") = {",Lcount,", ", 1, "};"

        ! Printing the curve loop
        write(10,"(A)",ADVANCE="NO") "Curve Loop(1) = {"
        do i = 1, Lcount-1
            write(10,"(I4,A)",ADVANCE="NO") i, ","
        end do
        write(10,"(I4,A)") Lcount, "};"

        ! Printing plane surface
        write(10,*) "Plane Surface(1) = {1};"

        ! Printing the physical curve "ogive1"
        write(10,"(A)",ADVANCE="NO") 'Physical Curve("ogive1") = {'
        ib = 1
        ie = ogive%sz-2
        do i = ib, ie
            write(10,"(I4,A)",ADVANCE="NO") i, ","
        end do
        write(10,"(I4,A)") ie+1,"};"

        ! Printing the physical curve "ogive2"
        write(10,"(A,I4,A)") 'Physical Curve("ogive2") = {', ie+2,"};"


        ! Printing the physical curve "outlet"
        write(10,"(A)",ADVANCE="NO") 'Physical Curve("outlet") = {'
        ib = ogive%sz+1
        ie = ogive%sz+outlet%sz-1
        do i = ib, ie
            write(10,"(I4,A)",ADVANCE="NO") i, ","
        end do
        write(10,"(I4,A)") ie+1,"};"


        ! Printing the physical curve "inlet"
        write(10,"(A)",ADVANCE="NO") 'Physical Curve("inlet") = {'
        ib = ogive%sz+outlet%sz+1
        ie = ib+inlet%sz-2
        do i = ib, ie
            write(10,"(I4,A)",ADVANCE="NO") i, ","
        end do
        write(10,"(I4,A)") ie+1,"};"

        ! Printing the physical curve "symmetry"
        write(10,"(A)",ADVANCE="NO") 'Physical Curve("symmetry") = {'
        ib = ie+2
        ie = ib+symmetry%sz-2
        do i = ib, ie
            write(10,"(I4,A)",ADVANCE="NO") i, ","
        end do
        write(10,"(I4,A)") ie+1,"};"

        ! Printing the domain
        write(10,*) 'Physical Surface("domain") = {1};'


        ! If km=0, the mesh is structured

        if ( par%km == 0 ) then

            write(10,"(A)",ADVANCE="NO") "Transfinite Surface {1} = {"
            write(10,"(I4,A)",ADVANCE="NO")        1, ","
            write(10,"(I4,A)",ADVANCE="NO") ogive%sz+1, ","
            write(10,"(I4,A)",ADVANCE="NO") ogive%sz+outlet%sz+1, ","
            write(10,"(I4,A)")              ogive%sz+outlet%sz+inlet%sz+1, "};"

            write(10,"(A)",ADVANCE="NO") "Transfinite Curve {"
            do i = 1, ogive%sz
                write(10,"(I4,A)",ADVANCE="NO") i, ","
            end do
            ib = ogive%sz+outlet%sz+1
            ie = ogive%sz+outlet%sz+inlet%sz-1
            do i = ib, ie
                write(10,"(I4,A)",ADVANCE="NO") i, ","
            end do
            write(10,"(I4,A)",ADVANCE="NO") ie+1, "} = "
            write(10,"(I4,A)") 1, " Using Progression 1;"


            write(10,"(A)",ADVANCE="NO") "Transfinite Curve {"
            ib = ogive%sz+1
            ie = ib+outlet%sz-1
            do i = ib, ie
                write(10,"(I4,A)",ADVANCE="NO") i, ","
            end do
            ib = ogive%sz+outlet%sz+inlet%sz+1
            ie = ib+symmetry%sz-2
            do i = ib, ie
                write(10,"(I4,A)",ADVANCE="NO") i, ","
            end do
            write(10,"(I4,A)",ADVANCE="NO") ie+1, "} = "
            write(10,"(I4,A)") 1, " Using Progression 1;"

            write(10,*) "Recombine Surface {1};"

        end if

    end subroutine


    !> \brief Generates the ogive boundary
    subroutine generate_ogive(bound)
        implicit none
        type(boundary), intent(inout) :: bound

        integer :: i
        integer :: iaux
        integer :: IO
        integer :: nlines
        integer :: np
        real(8) :: lr
        real(8) :: rb
        real(8) :: n
        real(8) :: aks
        real(8) :: raux
        real(8), allocatable :: xv(:)
        real(8), allocatable :: yv(:)

        type(class_path2d) :: path

        ! King of ogive (0=power-law, 1=load from file)
        if ( par%ko == 0 ) then

            lr = par%lr
            rb = par%rb
            n  = par%n

            aks = 2.d0 - 2.d0 * (n-0.5d0)

            ! Number of partitions
            np = 2000

            allocate( xv(0:np) )
            allocate( yv(0:np) )

            do i = 0, np

                xv(i) = (dble(i)/dble(np))**aks * lr
                yv(i) = ( xv(i) / lr ) ** n * rb

            end do

            call path%init(xv, yv)

            deallocate( xv )
            deallocate( yv )

        else ! Loads the ogive profile from file

            nlines = 0
            open(10, file=trim(par%ogfile))

            ! Skipping header
            read(10,*,IOStat=IO)
            do
                read(10,*,IOStat=IO) iaux, raux, raux
                if ( IO /= 0 ) exit
                nlines = nlines + 1
            end do

            if ( nlines < 3 ) then
                write(*,*) "generate_ogive: inapropriated ogive profile."
                stop
            end if

            allocate( xv(nlines) )
            allocate( yv(nlines) )

            rewind(10)

            ! Skipping header
            read(10,*,IOStat=IO)

            do i = 1, nlines

                read(10,*,IOStat=IO) iaux, xv(i), yv(i)

            end do

            call path%init(xv, yv)

            deallocate( xv )
            deallocate( yv )

        end if


        ! Distributing the points along the path to generate the discrete boundary
        call distribute_points(par%kpog, par%faog, par%fbog, path, bound)

    end subroutine


    !> \brief Generates the inlet boundary
    subroutine generate_inlet(bound)
        implicit none
        type(boundary), intent(inout) :: bound

        integer :: i
        real(8) :: la
        real(8) :: akn
        type(class_path2d) :: path

        la = par%lr + par%lf

        if ( par%kib == 0 ) then
            akn = 2.d0
        else
            akn = 2.d0 - 2.d0 * (par%nib-0.5d0)
        end if

        do i = 0, bound%sz

            bound%x(i) = (dble(i)/dble(bound%sz))**akn * la - par%lf

        end do

        if ( par%kib == 0 ) then

            do i = 0, bound%sz

                bound%y(i) = par%lb * sqrt(1.d0-((bound%x(i)-par%lr)/la)**2)

            end do
        else

            do i = 0, bound%sz

                bound%y(i) = par%lb * ((bound%x(i)+par%lf)/la) ** par%nib

            end do
        end if

        ! Initializing a parametric path
        call path%init(bound%x, bound%y)

        ! Distributing points along the discrete boundary
        call distribute_points(par%kpin, par%fain, par%fbin, path, bound)

    end subroutine


    !> \brief Generates the outlet boundary
    subroutine generate_outlet(bound)
        implicit none
        type(boundary), intent(inout) :: bound

        integer :: i

        do i = 0, bound%sz

            bound%x(i) = par%lr

            bound%y(i) = (dble(i)/dble(bound%sz))*(par%lb-par%rb)+par%rb

        end do

    end subroutine


    !> \brief Generates the symmetry boundary
    subroutine generate_symmetry(bound)
        implicit none
        type(boundary), intent(inout) :: bound

        integer :: i

        do i = 0, bound%sz

            bound%x(i) = (dble(i)/dble(bound%sz))*par%lf-par%lf

            bound%y(i) = 0.d0

        end do

    end subroutine


    !> \brief Distribute points according to a distribution
    subroutine distribute_points(kp, fa, fb, path, bound)
        implicit none
        integer,            intent(in)  :: kp ! Kind of partitioning (0=uniform, 1=monotonic cubic spline)
        real(8),            intent(in)  :: fa ! Factor for concentration of points near the left
        real(8),            intent(in)  :: fb ! Factor for concentration of points near the right
        type(class_path2d), intent(in)  :: path
        type(boundary),   intent(inout) :: bound

        integer :: i
        real(8), allocatable :: t(:)

        ! Distributing the points along the path to generate the discrete boundary

        allocate( t(0:bound%sz) )

        ! Uniform partioning
        if ( kp == 0 ) then
            call get_uniform_distribution(bound%sz+1, t)
        ! MCS partitioning
        else
            call get_mcs_distribution(bound%sz+1, fa, fb, t)
        end if

        do i = 0, bound%sz

            bound%x(i) = path%x(t(i))
            bound%y(i) = path%y(t(i))

        end do

        deallocate( t )

    end subroutine

end program
