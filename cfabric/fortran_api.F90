module fortran_api
    
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

        subroutine smm_init_time_control() bind(c,name="SMM_InitTestTimeControl")
        end subroutine smm_init_time_control

        subroutine smm_plot() bind(c, name="SMM_Plot")
        end subroutine smm_plot

        subroutine smm_save() bind(c, name="SMM_Save")
        end subroutine smm_save

        subroutine smm_driver() bind(c, name="SMM_Driver")
        end subroutine smm_driver

        subroutine smm_test_driver_no_fluid() bind(c, name="SMM_TestDriverNoFluid")
        end subroutine smm_test_driver_no_fluid

        subroutine smm_clean_up(exitcode) bind(c, name="SMM_CleanUp")
            import :: c_int
            integer(kind=c_int), optional, intent(in) :: exitcode
        end subroutine smm_clean_up
    end interface

end module fortran_api
