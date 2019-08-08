program fabric_test

    use, intrinsic :: iso_c_binding
    implicit none

!    interface
!        subroutine fabric_init(argc,argv) bind(c, name="Fabric_Init")
!            import :: c_int, c_char
!            integer(kind=c_int) :: argc
!            character(kind=c_char), dimension(0:argc) :: argv
!        end subroutine
!    end interface

    interface
        subroutine fabric_init(in_name) bind(c, name="Fabric_Init")
            import :: c_char
            character(kind=c_char) :: in_name
        end subroutine
    end interface

!    integer(kind=c_int), parameter :: argc = 5
!    character(kind=c_char,len=20), dimension(0:argc) :: argv
!
!    argv = [character(kind=c_char,len=20) :: &
!        "-d","3","-i","in-drumC","-o","out-drumC"] 

    !print *, argv
    
    character(kind=c_char,len=20), parameter :: in_name = &
        "ftinit.txt"//c_null_char

    call fabric_init(in_name)
    !call fabric_init(argc,argv)


end program fabric_test
