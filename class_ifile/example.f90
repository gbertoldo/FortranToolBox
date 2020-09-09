!    Class ifile usage example.
!
!    Copyright (C) 2020 by Guilherme Bertoldo 
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
!          Guilherme Bertoldo (a)
!                 E-mail: glbertoldo@gmail.com
!
!    Institution
!          (a) Federal University of Technology - Paraná - UTFPR
!              Linha Santa Bárbara, s/n, Francisco Beltrão, Paraná, Brazil
!              Zip Code 85601-970
!

program main
    use mod_class_ifile

    implicit none
    
    ! Create an instance of the class_ifile
    type(class_ifile) :: ifile
    
    integer           :: ivar
    real(8)           :: rvar
    real(8)           :: rvvar(2)
    character(len=50) :: cvar
    character(len=50) :: caux
    
    ! Initializing the object with the name of the file, "input.txt", and the field delimiter caracter, &.
    call ifile%init("input.txt","&")
    
    
    ! Loading the content of the file
    call ifile%load()
    
    
    ! Now we are ready to get the values of interest
    call ifile%get_value(ivar, "DAY")    
    write(*,"(I23,A)") ivar, " : integer"
    
    call ifile%get_value(rvar, "PI")
    write(*,"(F23.3,A)") rvar, " : real(8)"

    call ifile%get_value(cvar, "LABEL")
    write(*,"(A23,A)") trim(cvar), " : string"

    
    ! This is an workaround to read vectors
    call ifile%get_value(caux, "VEC")
    read(caux,*) rvvar
    write(*,"(ES11.3,ES12.3,A)") rvvar, " : vector of real(8)"


    write(*,*)
    write(*,*) "You may print all data that you have used (choose the unit to print):"
    call ifile%print_pickedup_data(6)    


    write(*,*)
    write(*,*) "You may print all data that was loaded from the input file (choose the unit to print):"
    call ifile%print_read_data(6)    

end program
