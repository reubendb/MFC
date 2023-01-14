!>
!! @file m_compress.f90
!! @brief Contains the Fortran bindings for compression.
!! @author Henry Le Berre


module m_compress

    use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_size_t, c_loc, c_int, c_double

    use m_nvtx

    implicit none

    type, bind(c) :: t_compress_config
        type(c_ptr)    :: pBuff = c_null_ptr
        integer(c_int) :: count = 0
        real(c_double) :: rate  = 0
        integer(c_int) :: from  = 0
        integer(c_int) :: to    = 0
    end type t_compress_config

    type, bind(c) :: t_compress_staging
        type(c_ptr)       :: pBuff = c_null_ptr
        integer(c_size_t) :: bytes = 0
    end type t_compress_staging

    type, bind(c) :: t_compress_state
        type(t_compress_config)  :: config
        type(c_ptr), private     :: pImpl   = c_null_ptr
        type(t_compress_staging) :: staging
    end type t_compress_state

    interface

        function c_compress_init(pConfig, pState) result(success) bind(c, name='c_compress_init')
            import

            type(c_ptr) :: pConfig
            type(c_ptr) :: pState

            logical(c_bool) :: success
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

        subroutine c_compress_bench() bind(c, name="c_compress_bench")
            import
        end subroutine c_compress_bench

    end interface

    real :: s_time, e_time
    real :: compress_time, mpi_time, decompress_time
    integer :: nCalls_time = 0

    contains

        function f_compress_init(config, state) result(success)
            type(t_compress_config), target :: config
            type(t_compress_state),  target :: state

            logical :: success

            success = c_compress_init(c_loc(config), c_loc(state))
        end function f_compress_init

        function f_compress(state) result(offset)
            type(t_compress_state), target, intent(in) :: state

            integer :: offset

            call nvtxStartRange("ZFP")

            call cpu_time(s_time)
            offset = c_compress(c_loc(state))
            call cpu_time(e_time)

            compress_time = compress_time + (e_time - s_time)
            nCalls_time   = nCalls_time + 1

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
            type(t_compress_state), target, intent(inout) :: state

            call c_compress_finalize(c_loc(state))
        end subroutine s_compress_finalize

        subroutine s_compress_recap()
            print*, "m_compress timings:"
            print*, " - nCalls_time  ", nCalls_time
            print*, " - s_compress   ", (compress_time   / nCalls_time), "s"
            print*, " - mpi_sendrecv ", (mpi_time        / nCalls_time), "s"
            print*, " - s_decompress ", (decompress_time / nCalls_time), "s"
        end subroutine s_compress_recap

        subroutine s_compress_bench()

            call c_compress_bench()

        end subroutine s_compress_bench

end module m_compress
