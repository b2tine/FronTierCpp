program fortran_driver

    use, intrinsic :: iso_c_binding
    implicit none

    interface
        subroutine smm_init(args) bind(c,name="SMM_Init")
            import :: c_char
            character(kind=c_char), dimension(*), intent(in) :: args
        end subroutine

        subroutine smm_init_modules() bind(c,name="SMM_InitModules")
        end subroutine

        subroutine smm_init_propagator() bind(c,name="SMM_InitPropagator")
        end subroutine

        subroutine smm_init_sm_params() bind(c,name="SMM_InitSpringMassParams")
        end subroutine

        subroutine smm_init_vel_func() bind(c,name="SMM_InitTestVelFunc")
        end subroutine

        subroutine smm_init_time_control() bind(c,name="SMM_InitTestTimeContrl")
        end subroutine

        subroutine smm_plot() bind(c, name="SMM_Plot")
        end subroutine

        subroutine smm_save() bind(c, name="SMM_Save")
        end subroutine

        subroutine smm_test_driver() bind(c, name="SMM_TestDriver")
        end subroutine

        subroutine smm_clean_up(exitcode) bind(c, name="SMM_CleanUp")
            import :: c_int
            integer(kind=c_int), optional, intent(in) :: exitcode
        end subroutine
    end interface


    integer :: i, j, k, l, m, n, argc, length, stat
    character(len=50) :: arg
    character(kind=c_char) :: d = c_null_char
    character(kind=c_char,len=50) :: ifile = c_null_char
    character(kind=c_char,len=50) :: ofile = c_null_char

    character(kind=c_char,len=256) :: args

    j = 0
    k = 0
    l = 0
    m = 0
    n = 0

    argc = command_argument_count()
    do i = 1, argc
        if (i == j .or. i == k .or. i == l) cycle
        call get_command_argument(i,arg,length=length,status=stat)

        select case(arg)

            case ("-d")
                j = i+1

            case ("-i")
                k = i+1

            case ("-o")
                l = i+1

            case ("-r")
                m = i+1

            case ("-t")
                n = i+1

        end select
    end do

    if (j /= 0) then
        call get_command_argument(j,arg,length=length,status=stat)
        read(unit=arg,fmt=*) d
        d = d//c_null_char
    end if

    if (k /= 0) then
        call get_command_argument(k,arg,length=length,status=stat)
        read(unit=arg,fmt=*) ifile
        ifile = ifile//c_null_char
    end if

    if (l /= 0) then
        call get_command_argument(l,arg,length=length,status=stat)
        read(unit=arg,fmt=*) ofile
        ofile = ofile//c_null_char
    end if


    !TODO: read in restart options -r and -t


    args = "-d "//d//" -i "//ifile//" -o "//ofile//c_null_char
    

    call smm_init(args)
    call smm_init_modules()
    call smm_init_sm_params()
    call smm_init_propagator()
    call smm_init_vel_func()
    call smm_init_time_control()

    call smm_plot()

    call smm_test_driver()
    call smm_clean_up()

end program fortran_driver
