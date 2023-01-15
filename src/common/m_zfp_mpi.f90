module m_zfp_mpi

    use, intrinsic :: iso_c_binding, only: c_ptr, c_signed_char, c_null_ptr, &
                                           c_size_t, c_loc, c_int, c_double

    use m_nvtx

    implicit none

    type, bind(c) :: t_zfp_mpi_config
        integer(c_signed_char) :: bitfield
        type   (c_ptr)         :: pDoubles
        integer(c_size_t)      :: nDoubles
        real   (c_double)      :: rate
    end type t_zfp_mpi_config

    type, bind(c) :: t_zfp_mpi_staging
        type   (c_ptr)    :: pBytes
        integer(c_size_t) :: nBytes
    end type t_zfp_mpi_staging

    type, bind(c) :: t_zfp_mpi_worker
        type(c_ptr)          :: pConfig   = c_null_ptr
        type(c_ptr), private :: pInternal = c_null_ptr
    end type t_zfp_mpi_worker

    interface

        subroutine i_zfp_mpi_init(pConfig, pWorker) bind(c, name='zfp_mpi_init')
            import

            type(c_ptr), value :: pConfig
            type(c_ptr), value :: pWorker
        end subroutine i_zfp_mpi_init
    
    end interface

    contains

        subroutine s_zfp_mpi_init(config, worker)
            type(t_zfp_mpi_config), intent(in)    :: config
            type(t_zfp_mpi_worker), intent(inout) :: worker

            call i_zfp_mpi_init(c_loc(config), c_loc(worker))
        end subroutine s_zfp_mpi_init

end module m_zfp_mpi