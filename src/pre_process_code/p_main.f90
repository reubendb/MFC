!!       __  _______________
!!      /  |/  / ____/ ____/
!!     / /|_/ / /_  / /     
!!    / /  / / __/ / /___   
!!   /_/  /_/_/    \____/   
!!                       
!!  This file is part of MFC.
!!
!!  MFC is the legal property of its developers, whose names 
!!  are listed in the copyright file included with this source 
!!  distribution.
!!
!!  MFC is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published 
!!  by the Free Software Foundation, either version 3 of the license 
!!  or any later version.
!!
!!  MFC is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!!  GNU General Public License for more details.
!!  
!!  You should have received a copy of the GNU General Public License
!!  along with MFC (LICENSE).  
!!  If not, see <http://www.gnu.org/licenses/>.

!>
!! @file p_main.f90
!! @brief Contains program p_main
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

!> @brief This program takes care of setting up the initial condition and
!!              grid data for the multicomponent flow code.
PROGRAM p_main
    
    
    ! Dependencies =============================================================
    USE m_derived_types         !< Definitions of the derived types
    
    USE m_global_parameters     !< Global parameters for the code
    
    USE m_mpi_proxy             !< Message passing interface (MPI) module proxy
    
    USE m_variables_conversion  !< Subroutines to change the state variables from
                                !! one form to another
    
    USE m_start_up              !< Procedures to read in and check consistency of
                                !! the user provided inputs and data
    
    USE m_grid                  !< Procedures to generate (non-)uniform grids
    
    USE m_initial_condition     !< Procedures to generate initial condition
    
    USE m_data_output           !< Procedures to write the grid data and the
                                !! conservative variables to files
    ! ==========================================================================
    
    
    IMPLICIT NONE
    

    ! Initialization of the MPI environment
    CALL s_mpi_initialize()
    
    
    ! Rank 0 processor assigns default values to user inputs prior to reading
    ! those in from the input file. Next, the user inputs are read in and their
    ! consistency is checked. The detection of any inconsistencies automatically
    ! leads to the termination of the pre-process.
    IF (proc_rank == 0) THEN
        CALL s_assign_default_values_to_user_inputs()
        CALL s_read_input_file()
        CALL s_check_input_file()
    END IF
   

    ! Broadcasting the user inputs to all of the processors and performing the
    ! parallel computational domain decomposition. Neither procedure has to be
    ! carried out if pre-process is in fact not truly executed in parallel.
    CALL s_mpi_bcast_user_inputs()
    CALL s_initialize_parallel_io()
    CALL s_mpi_decompose_computational_domain()
    
    ! Computation of parameters, allocation procedures, and/or any other tasks
    ! needed to properly setup the modules
    CALL s_initialize_global_parameters_module()
    CALL s_initialize_data_output_module()
    CALL s_initialize_variables_conversion_module()
    CALL s_initialize_start_up_module()
    CALL s_initialize_grid_module()
    CALL s_initialize_initial_condition_module()

    ! Associate pointers for serial or parallel I/O
    IF (parallel_io .NEQV. .TRUE.) THEN
        s_generate_grid => s_generate_serial_grid
        s_read_grid_data_files => s_read_serial_grid_data_files
        s_read_ic_data_files => s_read_serial_ic_data_files
        s_write_data_files => s_write_serial_data_files
    ELSE
        s_generate_grid => s_generate_parallel_grid
        s_read_grid_data_files => s_read_parallel_grid_data_files
        s_read_ic_data_files => s_read_parallel_ic_data_files
        s_write_data_files => s_write_parallel_data_files
    END IF

    ! Setting up the grid and the initial condition. If the grid is read in from
    ! preexisting grid data files, it is checked for consistency. If the grid is
    ! not read in, it is generated from scratch according to the inputs provided
    ! by the user. The initial condition may also be read in. It in turn is not
    ! checked for consistency since it WILL further be edited by the pre-process
    ! and also because it may be incomplete at the time it is read in. Finally,
    ! when the grid and initial condition are completely setup, they are written
    ! to their respective data files.
    
    ! Setting up grid and initial condition
    
    IF(old_grid) THEN
        CALL s_read_grid_data_files(dflt_int)
        CALL s_check_grid_data_files()
    ELSE
        IF (parallel_io .NEQV. .TRUE.) THEN        
                CALL s_generate_grid(dflt_int)
        ELSE
            IF (proc_rank == 0) CALL s_generate_grid(dflt_int)
            CALL s_mpi_barrier()
            CALL s_read_grid_data_files(dflt_int)
            CALL s_check_grid_data_files()
        END IF
    END IF
        
    IF(old_ic) CALL s_read_ic_data_files(q_cons_vf)

    CALL s_generate_initial_condition()
    CALL s_write_data_files(q_cons_vf)

    ! Disassociate pointers for serial and parallel I/O
    s_generate_grid => NULL()
    s_read_grid_data_files => NULL()
    s_read_ic_data_files => NULL()
    s_write_data_files => NULL()
    
    ! Deallocation procedures for the modules
    CALL s_finalize_initial_condition_module()
    CALL s_finalize_grid_module()
    CALL s_finalize_start_up_module()
    CALL s_finalize_variables_conversion_module()
    CALL s_finalize_data_output_module()
    CALL s_finalize_global_parameters_module()    

    ! Finalization of the MPI environment
    CALL s_mpi_finalize()
    
END PROGRAM p_main
