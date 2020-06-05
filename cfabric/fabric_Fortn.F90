program fortran_driver

    use fortran_api
    use, intrinsic :: iso_c_binding
    implicit none

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

    call smm_test_driver_no_fluid()
    call smm_clean_up()

end program fortran_driver
