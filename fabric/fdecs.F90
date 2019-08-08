module fdecs

    use, intrinsic :: iso_c_binding
    implicit none
    private

    !boolean
    enum, bind(c)
        enumerator :: false = 0
        enumerator :: no = false
        enumerator :: function_failed = false
        enumerator :: true = 1
        enumerator :: yes = true
        enumerator :: function_succeeded = true
    end enum

    public false, no, function_failed, true, yes, function_succeeded

    !GEOMETRY_REMAP
    enum, bind(c)
        enumerator :: invalid_remap = 0
        enumerator :: identity_remap = 1
        enumerator :: cylindrical_remap = 2
        enumerator :: spherical_remap = 3
    end enum


    !F_BASIC_DATA
    type, public, bind(c) :: f_basic_data
        
        integer(kind=c_int) :: dim
        character(LEN=*,kind=c_char) :: in_name
        character(LEN=*,kind=c_char) :: out_name
        integer(kind=c_int), dimension(0:2) :: subdomains

        integer(kind(true)) :: ReadFromInput
        integer(kind(false)) :: RestartRun
        integer(kind(false)) :: ResetTime
        
        integer(kind=c_int) :: RestartStep
        character(LEN=*,kind=c_char) :: restart_name
        character(LEN=*,kind=c_char) :: restart_state_name

        real(kind=c_double), dimension(0:2,0:2) :: L, U
        integer(kind=c_int), dimension(0:2) :: gmax
        real(kind=c_double), dimension(0:2,0:1) :: boundary

        integer(kind=c_size_t) :: size_of_intfc_state
        integer(kind(identity_remap)) :: ReadFromInput
    
    end type f_basic_data

end module fdecs
