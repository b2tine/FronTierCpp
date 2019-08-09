program fabric_test

    use, intrinsic :: iso_c_binding
    implicit none

    interface
        subroutine fabric_init(in_name) bind(c, name="Fabric_Init")
            import :: c_char
            character(kind=c_char) :: in_name
        end subroutine
    end interface

    
    character(kind=c_char,len=20), parameter :: in_name = &
        "ftinit.txt"//c_null_char

    call fabric_init(in_name)


end program fabric_test
