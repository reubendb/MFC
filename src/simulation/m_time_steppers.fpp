!>
!! @file m_time_steppers.f90
!! @brief Contains module m_time_steppers

#:include 'macros.fpp'

!> @brief The following module features a variety of time-stepping schemes.
!!              Currently, it includes the following Runge-Kutta (RK) algorithms:
!!                   1) 1st Order TVD RK
!!                   2) 2nd Order TVD RK
!!                   3) 3rd Order TVD RK
!!              where TVD designates a total-variation-diminishing time-stepper.
module m_time_steppers

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_rhs                  !< Right-hand-side (RHS) evaluation procedures

    use m_data_output          !< Run-time info & solution data output procedures

    use m_bubbles              !< Bubble dynamics routines

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_fftw

    use m_nvtx
    ! ==========================================================================

    implicit none

#ifdef CRAY_ACC_WAR
    @:CRAY_DECLARE_GLOBAL(type(vector_field), dimension(:), q_cons_ts)
    !! Cell-average conservative variables at each time-stage (TS)

    @:CRAY_DECLARE_GLOBAL(type(scalar_field), dimension(:), q_prim_vf)
    !! Cell-average primitive variables at the current time-stage

    @:CRAY_DECLARE_GLOBAL(type(scalar_field), dimension(:), rhs_vf)
    !! Cell-average RHS variables at the current time-stage

    @:CRAY_DECLARE_GLOBAL(type(vector_field), dimension(:), q_prim_ts)
    !! Cell-average primitive variables at consecutive TIMESTEPS

    integer, private :: num_ts !<
    !! Number of time stages in the time-stepping scheme
#else
    type(vector_field), allocatable, dimension(:) :: q_cons_ts !<
    !! Cell-average conservative variables at each time-stage (TS)

    type(scalar_field), allocatable, dimension(:) :: q_prim_vf !<
    !! Cell-average primitive variables at the current time-stage

    type(scalar_field), allocatable, dimension(:) :: rhs_vf !<
    !! Cell-average RHS variables at the current time-stage

    type(vector_field), allocatable, dimension(:) :: q_prim_ts !<
    !! Cell-average primitive variables at consecutive TIMESTEPS

    integer, private :: num_ts !<
    !! Number of time stages in the time-stepping scheme
#endif

!$acc declare link(q_cons_ts,q_prim_vf,rhs_vf,q_prim_ts)

contains

    !> The computation of parameters, the allocation of memory,
        !!      the association of pointers and/or the execution of any
        !!      other procedures that are necessary to setup the module.
    subroutine s_initialize_time_steppers_module() ! -----------------------

        type(int_bounds_info) :: ix_t, iy_t, iz_t !<
            !! Indical bounds in the x-, y- and z-directions

        integer :: i, j !< Generic loop iterators

        ! Setting number of time-stages for selected time-stepping scheme
        if (time_stepper == 1) then
            num_ts = 1
        elseif (any(time_stepper == (/2, 3/))) then
            num_ts = 2
        end if

        ! Setting the indical bounds in the x-, y- and z-directions
        ix_t%beg = -buff_size; ix_t%end = m + buff_size

        if (n > 0) then
            iy_t%beg = -buff_size; iy_t%end = n + buff_size

            if (p > 0) then
                iz_t%beg = -buff_size; iz_t%end = p + buff_size
            else
                iz_t%beg = 0; iz_t%end = 0
            end if
        else
            iy_t%beg = 0; iy_t%end = 0
            iz_t%beg = 0; iz_t%end = 0
        end if

        ! Allocating the cell-average conservative variables
        @:ALLOCATE_GLOBAL(q_cons_ts(1:num_ts))

        do i = 1, num_ts
            allocate(q_cons_ts(i)%vf(1:sys_size))
        end do

        do i = 1, num_ts
            do j = 1, sys_size
                allocate(q_cons_ts(i)%vf(j)%sf(ix_t%beg:ix_t%end, &
                                                iy_t%beg:iy_t%end, &
                                                iz_t%beg:iz_t%end))
            end do
            
            @:ACC_SETUP_VFs(q_cons_ts(i))
        end do

        ! Allocating the cell-average primitive ts variables
        if (probe_wrt) then
            @:ALLOCATE_GLOBAL(q_prim_ts(0:3))

            do i = 0, 3
                allocate(q_prim_ts(i)%vf(1:sys_size))
            end do

            do i = 0, 3
                do j = 1, sys_size
                    allocate(q_prim_ts(i)%vf(j)%sf(ix_t%beg:ix_t%end, &
                                                   iy_t%beg:iy_t%end, &
                                                   iz_t%beg:iz_t%end))
                end do
            end do

            do i = 0, 3
                @:ACC_SETUP_VFs(q_prim_ts(i))
            end do
        end if

        ! Allocating the cell-average primitive variables
        @:ALLOCATE_GLOBAL(q_prim_vf(1:sys_size))
        
        do i = 1, adv_idx%end
            @:ALLOCATE(q_prim_vf(i)%sf(ix_t%beg:ix_t%end, &
                                      iy_t%beg:iy_t%end, &
                                      iz_t%beg:iz_t%end))
        end do

        if (bubbles) then
            do i = bub_idx%beg, bub_idx%end
                @:ALLOCATE(q_prim_vf(i)%sf(ix_t%beg:ix_t%end, &
                                          iy_t%beg:iy_t%end, &
                                          iz_t%beg:iz_t%end))
            end do
        end if

        if (hypoelasticity) then

            do i = stress_idx%beg, stress_idx%end
                @:ALLOCATE(q_prim_vf(i)%sf(ix_t%beg:ix_t%end, &
                                          iy_t%beg:iy_t%end, &
                                          iz_t%beg:iz_t%end))
            end do
        end if

        if (model_eqns == 3) then
            do i = internalEnergies_idx%beg, internalEnergies_idx%end
                @:ALLOCATE(q_prim_vf(i)%sf(ix_t%beg:ix_t%end, &
                                          iy_t%beg:iy_t%end, &
                                          iz_t%beg:iz_t%end))
            end do
        end if

        ! Allocating the cell-average RHS variables
        @:ALLOCATE_GLOBAL(rhs_vf(1:sys_size))

        do i = 1, sys_size
            @:ALLOCATE(rhs_vf(i)%sf(0:m, 0:n, 0:p))
        end do

        ! Opening and writing the header of the run-time information file
        if (proc_rank == 0 .and. run_time_info) then
            call s_open_run_time_information_file()
        end if

    end subroutine s_initialize_time_steppers_module ! ---------------------

    !> 1st order TVD RK time-stepping algorithm
        !! @param t_step Current time step
    subroutine s_1st_order_tvd_rk(t_step, time_avg) ! --------------------------------

        integer, intent(IN) :: t_step
        real(kind(0d0)), intent(INOUT) :: time_avg

        integer :: i, j, k, l !< Generic loop iterator
        real(kind(0d0)) :: start, finish

        ! Stage 1 of 1 =====================================================

        call cpu_time(start)

        call nvtxStartRange("Time_Step")

        call s_compute_rhs(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, t_step)

#ifdef DEBUG
        print *, 'got rhs'
#endif

        if (run_time_info) then
            call s_write_run_time_information(q_prim_vf, t_step)
        end if

#ifdef DEBUG
        print *, 'wrote runtime info'
#endif

        if (probe_wrt) then
            call s_time_step_cycling(t_step)
        end if

        if (t_step == t_step_stop) return

!$acc parallel loop collapse(4) gang vector default(present)
        do i = 1, sys_size
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_ts(1)%vf(i)%sf(j, k, l) = &
                            q_cons_ts(1)%vf(i)%sf(j, k, l) &
                            + dt*rhs_vf(i)%sf(j, k, l)
                    end do
                end do
            end do
        end do

        !print *, q_cons_ts(1)%vf(cont_idx%beg)%sf(102,0,0)
        !print *, q_cons_ts(1)%vf(E_idx)%sf(102,0,0)
        !print *, q_cons_ts(1)%vf(adv_idx%end)%sf(102,0,0)
        !print *, q_cons_ts(1)%vf(mom_idx%beg)%sf(102,0,0)

        if (grid_geometry == 3) call s_apply_fourier_filter(q_cons_ts(1)%vf)

        if (model_eqns == 3) call s_pressure_relaxation_procedure(q_cons_ts(1)%vf)

        call nvtxEndRange

        call cpu_time(finish)

        if (t_step >= 4) then
            time_avg = (abs(finish - start) + (t_step - 4)*time_avg)/(t_step - 3)
        else
            time_avg = 0d0
        end if

        ! ==================================================================

    end subroutine s_1st_order_tvd_rk ! ------------------------------------

    !> 2nd order TVD RK time-stepping algorithm
        !! @param t_step Current time-step
    subroutine s_2nd_order_tvd_rk(t_step, time_avg) ! --------------------------------

        integer, intent(IN) :: t_step
        real(kind(0d0)), intent(INOUT) :: time_avg

        integer :: i, j, k, l !< Generic loop iterator
        real(kind(0d0)) :: start, finish

        ! Stage 1 of 2 =====================================================

        call cpu_time(start)

        call nvtxStartRange("Time_Step")

        call s_compute_rhs(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, t_step)

        if (run_time_info) then
            call s_write_run_time_information(q_prim_vf, t_step)
        end if

        if (probe_wrt) then
            call s_time_step_cycling(t_step)
        end if

        if (t_step == t_step_stop) return

!$acc parallel loop collapse(4) gang vector default(present)
        do i = 1, sys_size
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_ts(2)%vf(i)%sf(j, k, l) = &
                            q_cons_ts(1)%vf(i)%sf(j, k, l) &
                            + dt*rhs_vf(i)%sf(j, k, l)
                    end do
                end do
            end do
        end do

        if (grid_geometry == 3) call s_apply_fourier_filter(q_cons_ts(2)%vf)

        if (model_eqns == 3) call s_pressure_relaxation_procedure(q_cons_ts(2)%vf)
        ! ==================================================================

        ! Stage 2 of 2 =====================================================

        call s_compute_rhs(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, t_step)

!$acc parallel loop collapse(4) gang vector default(present)
        do i = 1, sys_size
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_ts(1)%vf(i)%sf(j, k, l) = &
                            (q_cons_ts(1)%vf(i)%sf(j, k, l) &
                             + q_cons_ts(2)%vf(i)%sf(j, k, l) &
                             + dt*rhs_vf(i)%sf(j, k, l))/2d0
                    end do
                end do
            end do
        end do

        if (grid_geometry == 3) call s_apply_fourier_filter(q_cons_ts(1)%vf)

        if (model_eqns == 3) call s_pressure_relaxation_procedure(q_cons_ts(1)%vf)

        call nvtxEndRange

        call cpu_time(finish)

        if (t_step >= 4) then
            time_avg = (abs(finish - start) + (t_step - 4)*time_avg)/(t_step - 3)
        else
            time_avg = 0d0
        end if

        ! ==================================================================

    end subroutine s_2nd_order_tvd_rk ! ------------------------------------

    !> 3rd order TVD RK time-stepping algorithm
        !! @param t_step Current time-step
    subroutine s_3rd_order_tvd_rk(t_step, time_avg) ! --------------------------------

        integer, intent(IN) :: t_step
        real(kind(0d0)), intent(INOUT) :: time_avg

        integer :: i, j, k, l !< Generic loop iterator
        real(kind(0d0)) :: start, finish

        ! Stage 1 of 3 =====================================================

        call cpu_time(start)

        call nvtxStartRange("Time_Step")
#ifdef CRAY_PRINT_DEBUG
print*, "step -> rhs"
#endif
        call s_compute_rhs(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, t_step)
!        call s_compute_rhs_full(q_cons_ts(1)%vf, q_prim_vf, rhs_vf, t_step)

#ifdef CRAY_PRINT_DEBUG
print*, "step -< rti"
#endif
        if (run_time_info) then
            call s_write_run_time_information(q_prim_vf, t_step)
        end if
#ifdef CRAY_PRINT_DEBUG
print*, "step -> time_step_cycling"
#endif
        if (probe_wrt) then
            call s_time_step_cycling(t_step)
        end if

        if (t_step == t_step_stop) return

        !$acc parallel loop collapse(4) gang vector default(present)
        do i = 1, sys_size
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_ts(2)%vf(i)%sf(j, k, l) = &
                            q_cons_ts(1)%vf(i)%sf(j, k, l) &
                            + dt*rhs_vf(i)%sf(j, k, l)
                    end do
                end do
            end do
        end do
#if CRAY_PRINT_DEBUG
print*, "step -> s_apply_fourier_filter"
#endif
        if (grid_geometry == 3) call s_apply_fourier_filter(q_cons_ts(2)%vf)
#if CRAY_PRINT_DEBUG
print*, "step -> s_pressure_relaxation_procedure"
#endif
        if (model_eqns == 3) call s_pressure_relaxation_procedure(q_cons_ts(2)%vf)

        ! ==================================================================

        ! Stage 2 of 3 =====================================================
#ifdef CRAY_PRINT_DEBUG
print*, "step -> s_compute_rhs"
#endif
        call s_compute_rhs(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, t_step)
!        call s_compute_rhs_full(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, t_step)

        !$acc parallel loop collapse(4) gang vector default(present)
        do i = 1, sys_size
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_ts(2)%vf(i)%sf(j, k, l) = &
                            (3d0*q_cons_ts(1)%vf(i)%sf(j, k, l) &
                             + q_cons_ts(2)%vf(i)%sf(j, k, l) &
                             + dt*rhs_vf(i)%sf(j, k, l))/4d0
                    end do
                end do
            end do
        end do

        if (grid_geometry == 3) call s_apply_fourier_filter(q_cons_ts(2)%vf)

        if (model_eqns == 3) call s_pressure_relaxation_procedure(q_cons_ts(2)%vf)

        ! ==================================================================

        ! Stage 3 of 3 =====================================================
        call s_compute_rhs(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, t_step)
!        call s_compute_rhs_full(q_cons_ts(2)%vf, q_prim_vf, rhs_vf, t_step)

!$acc parallel loop collapse(4) gang vector default(present)
        do i = 1, sys_size
            do l = 0, p
                do k = 0, n
                    do j = 0, m
                        q_cons_ts(1)%vf(i)%sf(j, k, l) = &
                            (q_cons_ts(1)%vf(i)%sf(j, k, l) &
                             + 2d0*q_cons_ts(2)%vf(i)%sf(j, k, l) &
                             + 2d0*dt*rhs_vf(i)%sf(j, k, l))/3d0
                    end do
                end do
            end do
        end do

        if (grid_geometry == 3) call s_apply_fourier_filter(q_cons_ts(1)%vf)

        if (model_eqns == 3) call s_pressure_relaxation_procedure(q_cons_ts(1)%vf)

        call nvtxEndRange

        call cpu_time(finish)

        time = time + (finish - start)

        if (t_step >= 4) then
            time_avg = (abs(finish - start) + (t_step - 4)*time_avg)/(t_step - 3)
        else
            time_avg = 0d0
        end if

        ! ==================================================================

    end subroutine s_3rd_order_tvd_rk ! ------------------------------------

    !> This subroutine saves the temporary q_prim_vf vector
        !!      into the q_prim_ts vector that is then used in p_main
        !! @param t_step current time-step
    subroutine s_time_step_cycling(t_step) ! ----------------------------

        integer, intent(IN) :: t_step

        integer :: i !< Generic loop iterator

        do i = 1, sys_size
!$acc update host(q_prim_vf(i)%sf)
        end do

        if (t_step == t_step_start) then
            do i = 1, sys_size
                q_prim_ts(3)%vf(i)%sf(:, :, :) = q_prim_vf(i)%sf(:, :, :)
            end do
        elseif (t_step == t_step_start + 1) then
            do i = 1, sys_size
                q_prim_ts(2)%vf(i)%sf(:, :, :) = q_prim_vf(i)%sf(:, :, :)
            end do
        elseif (t_step == t_step_start + 2) then
            do i = 1, sys_size
                q_prim_ts(1)%vf(i)%sf(:, :, :) = q_prim_vf(i)%sf(:, :, :)
            end do
        elseif (t_step == t_step_start + 3) then
            do i = 1, sys_size
                q_prim_ts(0)%vf(i)%sf(:, :, :) = q_prim_vf(i)%sf(:, :, :)
            end do
        else ! All other timesteps
            do i = 1, sys_size
                q_prim_ts(3)%vf(i)%sf(:, :, :) = q_prim_ts(2)%vf(i)%sf(:, :, :)
                q_prim_ts(2)%vf(i)%sf(:, :, :) = q_prim_ts(1)%vf(i)%sf(:, :, :)
                q_prim_ts(1)%vf(i)%sf(:, :, :) = q_prim_ts(0)%vf(i)%sf(:, :, :)
                q_prim_ts(0)%vf(i)%sf(:, :, :) = q_prim_vf(i)%sf(:, :, :)
            end do
        end if

    end subroutine s_time_step_cycling ! -----------------------------------

    !> Module deallocation and/or disassociation procedures
    subroutine s_finalize_time_steppers_module() ! -------------------------

        integer :: i, j !< Generic loop iterators

        ! Deallocating the cell-average conservative variables
        do i = 1, num_ts

            do j = 1, sys_size
                @:DEALLOCATE(q_cons_ts(i)%vf(j)%sf)
            end do

            @:DEALLOCATE(q_cons_ts(i)%vf)

        end do

        @:DEALLOCATE_GLOBAL(q_cons_ts)

        ! Deallocating the cell-average primitive ts variables
        if (probe_wrt) then
            do i = 0, 3
                do j = 1, sys_size
                    @:DEALLOCATE(q_prim_ts(i)%vf(j)%sf)
                end do
                @:DEALLOCATE(q_prim_ts(i)%vf)
            end do
            @:DEALLOCATE_GLOBAL(q_prim_ts)
        end if

        ! Deallocating the cell-average primitive variables
        do i = 1, adv_idx%end
            @:DEALLOCATE(q_prim_vf(i)%sf)
        end do

        if (hypoelasticity) then
            do i = stress_idx%beg, stress_idx%end
                @:DEALLOCATE(q_prim_vf(i)%sf)
            end do
        end if

        if (bubbles) then
            do i = bub_idx%beg, bub_idx%end
                @:DEALLOCATE(q_prim_vf(i)%sf)
            end do
        end if

        if (model_eqns == 3) then
            do i = internalEnergies_idx%beg, internalEnergies_idx%end
                @:DEALLOCATE(q_prim_vf(i)%sf)
            end do
        end if

        @:DEALLOCATE_GLOBAL(q_prim_vf)

        ! Deallocating the cell-average RHS variables
        do i = 1, sys_size
            @:DEALLOCATE(rhs_vf(i)%sf)
        end do

        @:DEALLOCATE_GLOBAL(rhs_vf)

        ! Writing the footer of and closing the run-time information file
        if (proc_rank == 0 .and. run_time_info) then
            call s_close_run_time_information_file()
        end if

    end subroutine s_finalize_time_steppers_module ! -----------------------

end module m_time_steppers