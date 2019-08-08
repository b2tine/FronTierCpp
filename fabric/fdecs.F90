module fdecs

    use, intrinsic :: iso_c_binding
    implicit none
    private

    type, public, bind(c) :: f_basic_data
        
        integer(kind=c_int) :: dim
        character(LEN=*,kind=c_char) :: in_name
        character(LEN=*,kind=c_char) :: out_name
        integer(kind=c_int), dimension(0:2) :: subdomains

        !TODO: need to deal with enum types

        !integer(kind=c_int) :: ReadFromInput !boolean
        !integer(kind=c_int) :: RestartRun    !boolean
        !integer(kind=c_int) :: ResetTime     !boolean
        
        integer(kind=c_int) :: RestartStep
        character(LEN=*,kind=c_char) :: restart_name
        character(LEN=*,kind=c_char) :: restart_state_name

        real(kind=c_double), dimension(0:2,0:2) :: L, U
        integer(kind=c_int), dimension(0:2) :: gmax
        real(kind=c_double), dimension(0:2,0:1) :: boundary

        integer(kind=c_size_t) :: size_of_intfc_state
        !GEOMETRY_REMAP enum
    
    end type f_basic_data

end module fdecs
