!>
!! @file m_compress.f90
!! @brief Contains the Fortran bindings for compression.
!! @author H. Le Berre


module m_compress

    use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_size_t, c_loc, c_int, c_double

    use nvtx

    implicit none

    type, bind(c) :: t_compress_state
        logical :: bActive = .false.

        type(c_ptr)    :: pDoubles = c_null_ptr
        integer(c_int) :: nDoubles = 0
        real(c_double) :: rate     = 0
        integer(c_int) :: from     = 0
        integer(c_int) :: to       = 0

        type(c_ptr)       :: pBytes = c_null_ptr
        integer(c_size_t) :: nBytes = 0

        type(c_ptr), private :: pInternal = c_null_ptr
    end type t_compress_state

    interface
        
        function c_compress_init(pBytes, nBytes, rate, from, to) result(state) bind(c, name='c_compress_init')
            import
            
            type(c_ptr),       value :: pBytes
            integer(c_size_t), value :: nBytes
            real(c_double),    value :: rate
            integer(c_int),    value :: from
            integer(c_int),    value :: to

            type(t_compress_state) :: state
        end function c_compress_init

        function c_compress(pState) result(offset) bind(c, name='c_compress')
            import

            type(c_ptr), value :: pState
            
            integer(c_size_t) :: offset
        end function c_compress

        function c_decompress(pState) result(offset) bind(c, name='c_decompress')
            import

            type(c_ptr), value :: pState
            
            integer(c_size_t) :: offset
        end function c_decompress

        subroutine c_compress_finalize(pState) bind(c, name='c_compress_finalize')
            import

            type(c_ptr), value :: pState
        end subroutine c_compress_finalize

        subroutine c_bench() bind(c, name="c_bench")
            import
        end subroutine c_bench

    end interface

    real :: s_time, e_time
    real :: compress_time, mpi_time, decompress_time
    integer :: nCalls_time = 0

    contains
        
        function f_compress_init(pDoubles, nDoubles, rate, from, to) result(state)
            !real*8, intent(in)  :: doubles(*)
            type(c_ptr), intent(in) :: pDoubles
            integer,     intent(in) :: nDoubles
            real,        intent(in) :: rate
            integer,     intent(in) :: from
            integer,     intent(in) :: to
            
            type(t_compress_state) :: state

            state = c_compress_init(pDoubles, int(nDoubles, c_size_t), real(rate, c_double), int(from, c_int), int(to, c_int))
        end function f_compress_init

        function f_compress(state) result(offset)
            type(t_compress_state), target, intent(in) :: state
            
            integer :: offset
            
            call nvtxStartRange("ZFP")

            call cpu_time(s_time)
            offset = c_compress(c_loc(state))
            call cpu_time(e_time)

            compress_time = compress_time + (e_time - s_time)
            nCalls_time = nCalls_time + 1

            print*, "f_compress"

            call nvtxEndRange()
        end function f_compress

        function f_decompress(state) result(offset)
            type(t_compress_state), target, intent(in) :: state
            
            integer :: offset

            call nvtxStartRange("ZFP")
            
            call cpu_time(s_time)
            offset = c_decompress(c_loc(state))
            call cpu_time(e_time)
            
            decompress_time = decompress_time + (e_time - s_time)
            
            call nvtxEndRange()
        end function f_decompress

        subroutine s_compress_finalize(state)
            type(t_compress_state), intent(inout) :: state

            call c_compress_finalize(c_loc(state))
        end subroutine s_compress_finalize

        subroutine s_compress_recap()
            print*, "m_compress timings:"
            print*, " - nCalls_time  ", nCalls_time
            print*, " - s_compress   ", (compress_time   / nCalls_time), "s"
            print*, " - mpi_sendrecv ", (mpi_time        / nCalls_time), "s"
            print*, " - s_decompress ", (decompress_time / nCalls_time), "s"
        end subroutine s_compress_recap

        subroutine s_bench()

            call c_bench()

        end subroutine s_bench

end module m_compress
