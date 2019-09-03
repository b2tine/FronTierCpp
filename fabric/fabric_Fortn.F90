program fortran_driver

    use, intrinsic :: iso_c_binding
    implicit none

    interface
        subroutine smm_init(args) bind(c,name="SMM_Init")
            import :: c_char
            character(kind=c_char), dimension(*), intent(in) :: args
        end subroutine smm_init

        subroutine smm_init_modules() bind(c,name="SMM_InitModules")
        end subroutine smm_init_modules

        subroutine smm_init_propagator() bind(c,name="SMM_InitPropagator")
        end subroutine smm_init_propagator

        subroutine smm_init_sm_params() bind(c,name="SMM_InitSpringMassParams")
        end subroutine smm_init_sm_params

        subroutine smm_init_vel_func() bind(c,name="SMM_InitTestVelFunc")
        end subroutine smm_init_vel_func

        subroutine smm_init_time_control() bind(c,name="SMM_InitTestTimeContrl")
        end subroutine smm_init_time_control

        subroutine smm_plot() bind(c, name="SMM_Plot")
        end subroutine smm_plot

        subroutine smm_save() bind(c, name="SMM_Save")
        end subroutine smm_save

        subroutine smm_test_driver() bind(c, name="SMM_TestDriver")
        end subroutine smm_test_driver

        subroutine smm_clean_up(exitcode) bind(c, name="SMM_CleanUp")
            import :: c_int
            integer(kind=c_int), optional, intent(in) :: exitcode
        end subroutine smm_clean_up
    end interface


    character(len=50) :: arg
    integer :: i, k, argc, length, stat
    
    character(kind=c_char,len=10) :: d = c_null_char
    character(kind=c_char,len=50) :: ifile = c_null_char
    character(kind=c_char,len=50) :: ofile = c_null_char
    character(kind=c_char,len=50) :: rfile = c_null_char
    character(kind=c_char,len=10) :: rstep = c_null_char

    character(kind=c_char,len=:), allocatable :: args


    k = 0
    argc = command_argument_count()

    do i = 1, argc
        if (i == k) cycle
        call get_command_argument(i,arg,length=length,status=stat)

        select case(arg)

            case ("-d")
                k = i+1
                call get_command_argument(k,arg,length=length,status=stat)
                read(unit=arg,fmt="(a)") d
                d = "-d "//d//c_null_char

            case ("-i")
                k = i+1
                call get_command_argument(k,arg,length=length,status=stat)
                read(unit=arg,fmt="(a)") ifile
                ifile = "-i "//ifile//c_null_char

            case ("-o")
                k = i+1
                call get_command_argument(k,arg,length=length,status=stat)
                read(unit=arg,fmt="(a)") ofile
                ofile = "-o "//ofile//c_null_char

            case ("-r")
                k = i+1
                call get_command_argument(k,arg,length=length,status=stat)
                read(unit=arg,fmt="(a)") rfile
                rfile = "-r "//rfile//c_null_char

            case ("-t")
                k = i+1
                call get_command_argument(k,arg,length=length,status=stat)
                read(unit=arg,fmt="(a)") rstep
                rstep = "-t "//rstep//c_null_char

        end select
    end do


    args = d//" "//ifile//" "//ofile
    args = args//" "//rfile//" "//rstep
    

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
