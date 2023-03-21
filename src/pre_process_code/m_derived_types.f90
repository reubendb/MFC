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
!! @file m_derived_types.f90
!! @brief Contains module m_derived_types
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

!> @brief This file contains the definitions of all of the custom-defined
!!              types used in the pre-process code.
MODULE m_derived_types
    
    IMPLICIT NONE
    
    INTEGER, PARAMETER :: num_patches_max   = 10 !<
    !! Maximum number of patches allowed
    
    INTEGER, parameter :: num_fluids_max    = 10 !<
    !! Maximum number of fluids allowed
     
    !> Derived type annexing a scalar field (SF)
    TYPE scalar_field
        REAL(KIND(0d0)), POINTER, DIMENSION(:,:,:) :: sf => NULL()
    END TYPE scalar_field

    TYPE mpi_io_var
        INTEGER, ALLOCATABLE, DIMENSION(:) :: view
        TYPE(scalar_field), ALLOCATABLE, DIMENSION(:) :: var
    END TYPE  mpi_io_var

    !> Integer boounds for variables
    TYPE int_bounds_info
        INTEGER :: beg
        INTEGER :: end
    END TYPE int_bounds_info
    
    !> Derived type adding beginning (beg) and end bounds info as attributes
    TYPE bounds_info
        REAL(KIND(0d0)) :: beg
        REAL(KIND(0d0)) :: end
    END TYPE bounds_info
  
    !> bounds for the bubble dynamic variables
    TYPE bub_bounds_info
        INTEGER :: beg
        INTEGER :: end
        INTEGER, DIMENSION(:), ALLOCATABLE :: rs
        INTEGER, DIMENSION(:), ALLOCATABLE :: vs
        INTEGER, DIMENSION(:), ALLOCATABLE :: ps
        INTEGER, DIMENSION(:), ALLOCATABLE :: ms
        INTEGER, DIMENSION(:,:), ALLOCATABLE :: moms !< Moment indices for qbmm
        INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: fullmom !< Moment indices for qbmm
    END TYPE bub_bounds_info
    
    !> Derived type adding initial condition (ic) patch parameters as attributes
    !! NOTE: The requirements for the specification of the above parameters
    !! are strongly dependent on both the choice of the multicomponent flow
    !! model as well as the choice of the patch geometry.
    TYPE ic_patch_parameters
        
        INTEGER :: geometry !< Type of geometry for the patch
        
        REAL(KIND(0d0)) :: x_centroid, y_centroid, z_centroid !<
        !! Location of the geometric center, i.e. the centroid, of the patch. It
        !! is specified through its x-, y- and z-coordinates, respectively.
        
        REAL(KIND(0d0)) :: length_x, length_y, length_z !< Dimensions of the patch. x,y,z Lengths.
        REAL(KIND(0d0)) :: radius !< Dimensions of the patch. radius.

        REAL(KIND(0d0)), DIMENSION(3) :: radii !< 
        !! Vector indicating the various radii for the elliptical and ellipsoidal
        !! patch geometries. It is specified through its x-, y-, and z-components
        !! respectively.
 
        REAL(KIND(0d0)) :: epsilon, beta !< 
        !! The isentropic vortex parameters administrating, respectively, both
        !! the amplitude of the disturbance as well as its domain of influence.
        
        REAL(KIND(0d0)), DIMENSION(3) :: normal !<
        !! Normal vector indicating the orientation of the patch. It is specified
        !! through its x-, y- and z-components, respectively.
       
        LOGICAL, DIMENSION(0:num_patches_max-1) :: alter_patch !<
        !! List of permissions that indicate to the current patch which preceding
        !! patches it is allowed to overwrite when it is in process of being laid
        !! out in the domain
        
        LOGICAL :: smoothen !<
        !! Permission indicating to the current patch whether its boundaries will
        !! be smoothed out across a few cells or whether they are to remain sharp
        
        INTEGER :: smooth_patch_id !<
        !! Identity (id) of the patch with which current patch is to get smoothed
        
        REAL(KIND(0d0)) :: smooth_coeff !<
        !! Smoothing coefficient (coeff) adminstrating the size of the stencil of
        !! cells across which boundaries of the current patch will be smeared out
 
        REAL(KIND(0d0)), DIMENSION(num_fluids_max) :: alpha_rho
        REAL(KIND(0d0))                            :: rho
        REAL(KIND(0d0)), DIMENSION(3)              :: vel
        REAL(KIND(0d0))                            :: pres
        REAL(KIND(0d0)), DIMENSION(num_fluids_max) :: alpha
        REAL(KIND(0d0))                            :: gamma
        REAL(KIND(0d0))                            :: pi_inf !<
        !! Primitive variables associated with the patch. In order, these include
        !! the partial densities, density, velocity, pressure, volume fractions,
        !! specific heat ratio function and the liquid stiffness function.

        REAL(KIND(0d0)), DIMENSION(6)              :: tau_e
        !! Elastic stresses added to primitive variables if hypoelasticity = True
        
        REAL(KIND(0d0))    :: R0 !< Bubble size
        REAL(KIND(0d0))    :: V0 !< Bubble velocity

        REAL(KIND(0d0))    :: p0 !< Bubble size
        REAL(KIND(0d0))    :: m0 !< Bubble velocity
        
        
    END TYPE ic_patch_parameters
    
    !> Derived type annexing the physical parameters (PP) of the fluids. These
    !! include the specific heat ratio function and liquid stiffness function.
    TYPE physical_parameters
        REAL(KIND(0d0)) :: gamma
        REAL(KIND(0d0)) :: pi_inf
        REAL(KIND(0d0)) :: qv
        REAL(KIND(0d0)) :: mul0
        REAL(KIND(0d0)) :: ss
        REAL(KIND(0d0)) :: pv
        REAL(KIND(0d0)) :: gamma_v
        REAL(KIND(0d0)) :: M_v
        REAL(KIND(0d0)) :: mu_v
        REAL(KIND(0d0)) :: k_v
        REAL(KIND(0d0)) :: G
    END TYPE physical_parameters
    
    
END MODULE m_derived_types
