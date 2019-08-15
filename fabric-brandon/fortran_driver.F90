program fortran_driver

    use, intrinsic :: iso_c_binding
    implicit none

    interface
        subroutine fabric_test(in_name) bind(c, name="Fabric_test")
            import :: c_char
            character(kind=c_char) :: in_name
        end subroutine
    end interface

    
    character(kind=c_char,len=20), parameter :: in_name = &
        "ftinit.txt"//c_null_char

    call fabric_test(in_name)


end program fortran_driver
