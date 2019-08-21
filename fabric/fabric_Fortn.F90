program fortran_driver

    use, intrinsic :: iso_c_binding
    implicit none

    interface
        subroutine smm_init(in_name) bind(c,name="SMM_Init")
            import :: c_char
            character(kind=c_char) :: in_name
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
    end interface

    
    character(kind=c_char,len=20), parameter :: in_name = &
        "ftinit.txt"//c_null_char

    call smm_init(in_name)
    call smm_init_modules()
    call smm_init_sm_params()
    call smm_init_propagator()
    call smm_init_vel_func()
    call smm_init_time_control()

    call smm_plot()

    call smm_test_driver()

end program fortran_driver
