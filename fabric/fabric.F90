program fabric_test

    use, intrinsic :: iso_c_binding
    implicit none

    interface
        subroutine fabric_init(argc,argv) bind(c, name="Fabric_Init")
            import :: c_int, c_char
            integer(kind=c_int) ::argc
            character(kind=c_char), dimension(0:6) :: argv
        end subroutine
    end interface


    integer(kind=c_int) :: argc = 6
    character(kind=c_char), dimension(0:6) :: argv

    argv = ["-d","3","-i","in-drumC","-o","out-drumC"] 

    print *, argv


end program fabric_test
