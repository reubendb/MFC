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
!! @file m_phasechange.f90
!! @brief Contains module m_phasechange
!! @author M. Rodriguez
!! @version 1.1
!! @date November 10, 2021

! TODO 
! 1. improve iteration method to include the Riemann method of Chiapolino
! 2. Add Pelanti's finite relaxation approach to the code
! 3. Add changes to the GPU branch using Spencer's suggestion of git diff of
! start of branch with latest commit

!> @brief This module is used to compute phase relaxation for pressure,
!         temperature and chemical interfacial relaxation
MODULE m_phasechange

    ! Dependencies =============================================================

    USE m_derived_types        !< Definitions of the derived types
    USE m_global_parameters    !< Definitions of the global parameters
    USE m_mpi_proxy            !< Message passing interface (MPI) module proxy
    USE m_variables_conversion !< State variables type conversion procedures

    ! ==========================================================================

    IMPLICIT NONE

    PRIVATE; PUBLIC :: s_initialize_phasechange_module, &
                       s_relaxation_solver,             & 
                       s_relaxation_finite_solver,      & 
                       s_finite_ptg_relaxation,         &
                       s_infinite_p_relaxation,         &
                       s_infinite_p_relaxation_k,       &
                       s_infinite_pt_relaxation,        &
                       s_infinite_pt_relaxation_k,      &
                       s_infinite_ptg_relaxation,       &
                       s_infinite_ptg_relaxation_k

    !> @name Abstract interface for creating function pointers
    !> @{
    ABSTRACT INTERFACE

        !> @name Abstract subroutine for the finite relaxation solver 
        !> @{
        SUBROUTINE s_abstract_relaxation_solver(q_cons_vf) ! -------
            IMPORT :: scalar_field, sys_size          
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(INOUT) :: q_cons_vf
        END SUBROUTINE

        !> @name Abstract subroutine for the finite relaxation solver 
        !> @{
        SUBROUTINE s_abstract_relaxation_finite_solver(q_cons_vf, rhs_vf) ! -------
            IMPORT :: scalar_field, sys_size
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN) :: q_cons_vf
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(INOUT) :: rhs_vf
        END SUBROUTINE

    END INTERFACE

    !> @name Parameters for the phase change part of the code
    !> @{
    INTEGER,         PARAMETER :: newton_iter       = 50        !< p_relaxk \alpha iter,                set to 25
    REAL(KIND(0d0)), PARAMETER :: pknewton_eps      = 1.d-15    !< p_relaxk \alpha threshold,           set to 1E-15
    REAL(KIND(0d0)), PARAMETER :: pTsatnewton_eps   = 1.d-10    !< Saturation temperature tol,          set to 1E-10
    REAL(KIND(0d0)), PARAMETER :: ptgnewton_eps     = 1.d-8     !< Saturation p-T-mu tolerance,         set to 1.d-10
    REAL(KIND(0d0)), PARAMETER :: pres_crit         = 22.06d6   !< Critical water pressure              set to 22.06d6
    REAL(KIND(0d0)), PARAMETER :: T_crit            = 648.d0    !< Critical water temperature           set to 648
    REAL(KIND(0d0)), PARAMETER :: TsatHv            = 1000.d0   !< Saturation temperature threshold,    set to 900
    REAL(KIND(0d0)), PARAMETER :: TsatLv            = 250.d0    !< Saturation temperature threshold,    set to 250
    !> @}

    !> @name Gibbs free energy phase change parameters
    !> @{
    REAL(KIND(0d0)) :: n1, n2, pinf1, pinf2
    REAL(KIND(0d0)) :: gibbsA, gibbsB, gibbsC, gibbsD
    !> @}

    PROCEDURE(s_abstract_relaxation_solver), & 
    POINTER :: s_relaxation_solver => NULL()

    PROCEDURE(s_abstract_relaxation_finite_solver), &
    POINTER :: s_relaxation_finite_solver => NULL()


    CONTAINS

        !>     The purpose of this subroutine is to determine the saturation
        !!         temperature by using a Newton-Raphson method from the provided
        !!         equilibrium pressure and EoS of the binary phase system.
        !!     @param q_cons_vf Cell-average conservative variables
        !!     @param p_star equilibrium pressure at the interface    
        SUBROUTINE s_initialize_phasechange_module()
            n1    = 1.d0/fluid_pp(1)%gamma + 1.d0
            pinf1 = fluid_pp(1)%pi_inf/(1.d0 + fluid_pp(1)%gamma)
            n2    = 1.d0/fluid_pp(2)%gamma + 1.d0
            pinf2 = fluid_pp(2)%pi_inf/(1.d0 + fluid_pp(2)%gamma)
            gibbsA = (n1*fluid_pp(1)%cv - n2*fluid_pp(2)%cv +  & 
                     fluid_pp(2)%qvp - fluid_pp(1)%qvp) / &
                     (n2*fluid_pp(2)%cv - fluid_pp(2)%cv)
            gibbsB = (fluid_pp(1)%qv - fluid_pp(2)%qv) / &
                     (n2*fluid_pp(2)%cv - fluid_pp(2)%cv)
            gibbsC = (n2*fluid_pp(2)%cv - n1*fluid_pp(1)%cv) / &
                     (n2*fluid_pp(2)%cv - fluid_pp(2)%cv)
            gibbsD = (n1*fluid_pp(1)%cv - fluid_pp(1)%cv) / & 
                     (n2*fluid_pp(2)%cv - fluid_pp(2)%cv)

            ! Associating procedural pointer to the subroutine that will be
            ! utilized to calculate the solution of a given Riemann problem
        
            ! Opening and writing the header of the run-time information file
            IF(relax_model == 0 .OR. relax_model == 1) THEN
                s_relaxation_solver => s_infinite_p_relaxation
            ELSEIF (relax_model == 2) THEN
                s_relaxation_solver => s_infinite_pt_relaxation
            ELSEIF (relax_model == 3) THEN
                s_relaxation_solver => s_infinite_ptg_relaxation
            ELSEIF (relax_model == 4) THEN
                s_relaxation_solver => s_infinite_p_relaxation_k
            ELSEIF (relax_model == 5) THEN
                s_relaxation_solver => s_infinite_pt_relaxation_k
            ELSEIF (relax_model == 6) THEN
                s_relaxation_solver => s_infinite_ptg_relaxation_k      
            ELSE
                PRINT '(A)', 'relaxation solver was not set!'
                CALL s_mpi_abort()
            END IF

            !IF (relax_model == 1) THEN
            !    s_relaxation_finite_solver => s_finite_ptg_relaxation
            !END IF      

        END SUBROUTINE s_initialize_phasechange_module !-------------------------------

        !> The purpose of this procedure is to employ the inputted
        !!      cell-average conservative variables in order to compute
        !!      the cell-average RHS variables of the semidiscrete form
        !!      of the governing equations by utilizing the appropriate
        !!      Riemann solver.        
        !!  @param q_cons_vf Cell-average conservative variables
        !!  @param q_prim_vf Cell-average primitive variables
        !!  @param rhs_vf Cell-average RHS variables
        SUBROUTINE s_finite_ptg_relaxation(q_cons_vf, rhs_vf) ! -------

            TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN) :: q_cons_vf
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(INOUT) :: rhs_vf

            REAL(KIND(0d0))                                   ::  sum_alpha, Tsat
            REAL(KIND(0d0)), DIMENSION(num_fluids)            ::  p_k, T_k, g_k, Z_k
            REAL(KIND(0d0))                                   ::  n_k, pinf_k
            REAL(KIND(0d0))                                   ::  rho_k, rhoeq_k, rhoe
            REAL(KIND(0d0))                                   ::  e_k, phi, psi
            REAL(KIND(0d0))                                   ::  f1, f2, f3, f4
            REAL(KIND(0d0))                                   ::  A1, B1, A2, B2
            REAL(KIND(0d0))                                   ::  rho_I, Q, kappa
            REAL(KIND(0d0))                                   ::  deltap, p_I, mu 
            REAL(KIND(0d0))                                   ::  e_I, mdot, nu, theta
            REAL(KIND(0d0))                                   ::  mdotalpha, mdotrhoe
            REAL(KIND(0d0)), DIMENSION(2)                     ::          Re
            REAL(KIND(0d0)), DIMENSION(num_fluids,num_fluids) ::          We

            INTEGER :: i,j,k,l,r !< Generic loop iterators

            DO j = 0, m
                DO k = 0, n
                    DO l = 0, p
                       ! Numerical correction of the volume fractions
                       IF (mpp_lim) THEN
                       END IF

                       IF (q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) .GT. palpha_eps .OR. &
                           q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) .LT. 1.d0-palpha_eps) THEN
       
                       ! computing n_k, pinf_k, p_k, T_k, and g_k for finite relaxation
                       phi = 0.d0; psi = 0.d0; f1 = 0.d0; f2 = 0.d0; f3 = 0.d0; f4 = 0.d0;
                       !rhoe = 0.d0;

                       DO i = 1, num_fluids
                          n_k    = 1.d0/fluid_pp(i)%gamma + 1.d0
                          pinf_k = fluid_pp(i)%pi_inf/(1.d0 + fluid_pp(i)%gamma)
                          rho_k  = q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l) &
                                  /q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l)
                          rhoeq_k = (q_cons_vf(i+internalEnergies_idx%beg-1)%sf(j,k,l) & 
                              -q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)*fluid_pp(i)%qv) &
                              /q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) 
                          !rhoe = rhoe + q_cons_vf(i+internalEnergies_idx%beg-1)%sf(j,k,l)
                          e_k = q_cons_vf(i+internalEnergies_idx%beg-1)%sf(j,k,l) &
                               /q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)

                          p_k(i) = (rhoeq_k-fluid_pp(i)%pi_inf)/fluid_pp(i)%gamma
                          T_k(i) = (rhoeq_k-fluid_pp(i)%pi_inf/(1.d0+fluid_pp(i)%gamma)) &
                              /(rho_k*fluid_pp(i)%cv)
                          g_k(i) = (n_k*fluid_pp(i)%cv-fluid_pp(i)%qvp)*T_k(i) &
                              -fluid_pp(i)%cv*T_k(i)*log(T_k(i)**(n_k) &
                              /((p_k(i)+pinf_k)**(n_k-1.d0)))+fluid_pp(i)%qv
                          Z_k(i) = n_k*(p_k(i)+pinf_k);

                          phi = phi + 1.d0/(q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)*fluid_pp(i)%cv)
                          psi = psi + 1.d0/(fluid_pp(i)%gamma*q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l))

                          f1 = f1 + (p_k(i)+n_k*pinf_k)/q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l)
                          f2 = f2 + pinf_k/(q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)*fluid_pp(i)%cv)
                          f3 = f3 + fluid_pp(i)%qv/(fluid_pp(i)%gamma & 
                                    *q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l))
                          f4 = f4 + (e_k-pinf_k/rho_k) &
                                    /(q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)*fluid_pp(i)%cv)
                       END DO
                       !TODO IMPROVE THIS APPROACH
                       Tsat = f_Tsat(p_k(1))
                       !PRINT *,'prelax =',p_k(1),'Trelax = ',T_k(1),', Tsat = ',Tsat

                       IF ( ISNAN(p_k(1)) .OR. p_k(1) < 0.d0 ) THEN 
                           PRINT *, 'code crashed' 
                           CALL s_mpi_abort()
                       END IF

                       kappa = f1/psi
                       rho_I = (phi*f1-psi*f2)/(psi*f4-phi*f3)
                       p_I = (Z_k(2)*p_k(1)+Z_k(1)*p_k(2))/(Z_k(1)+Z_k(2));
                       e_I = f4/phi + f2/(rho_I*phi)

                       !mu = 1.d8
                       mu = 0.d0;
                       theta = 1.d8 
                       nu = 0.d0
                       IF (T_k(1) .GT. Tsat) nu = 1d-3
  
                       deltap = mu*(p_k(1)-p_k(2))
                       Q = theta*(T_k(2)-T_k(1))
                       mdot = nu*(g_k(2)-g_k(1))
                       mdotalpha = mdot/rho_I
                       mdotrhoe = mdot*e_I

                       rhs_vf(1+adv_idx%beg-1)%sf(j,k,l) = & 
                              rhs_vf(1+adv_idx%beg-1)%sf(j,k,l) + deltap + Q/kappa + mdotalpha
                       rhs_vf(2+adv_idx%beg-1)%sf(j,k,l) = & 
                              rhs_vf(2+adv_idx%beg-1)%sf(j,k,l) - deltap - Q/kappa - mdotalpha
                       rhs_vf(1+cont_idx%beg-1)%sf(j,k,l) = & 
                              rhs_vf(1+cont_idx%beg-1)%sf(j,k,l) + mdot
                       rhs_vf(2+cont_idx%beg-1)%sf(j,k,l) = & 
                              rhs_vf(2+cont_idx%beg-1)%sf(j,k,l) - mdot
                       rhs_vf(1+internalEnergies_idx%beg-1)%sf(j,k,l) = &
                              rhs_vf(1+internalEnergies_idx%beg-1)%sf(j,k,l) - p_I*deltap + Q + mdotrhoe
                       rhs_vf(2+internalEnergies_idx%beg-1)%sf(j,k,l) = &
                              rhs_vf(2+internalEnergies_idx%beg-1)%sf(j,k,l) + p_I*deltap - Q - mdotrhoe
                       END IF
                    END DO
                END DO
            END DO
        END SUBROUTINE s_finite_ptg_relaxation ! --------------------------------------

        !> Description: The purpose of this procedure is to infinitely relax
        !!              the pressures from the internal-energy equations to a
        !!              unique pressure, from which the corresponding volume
        !!              fraction of each phase are recomputed. For conservation
        !!              purpose, this pressure is finally corrected using the
        !!              mixture-total-energy equation.
        SUBROUTINE s_infinite_p_relaxation_k(q_cons_vf) ! ----------------        

            ! Cell-average conservative variables
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(INOUT) :: q_cons_vf
         
            ! Relaxed pressure, initial partial pressures, function f(p) and its partial
            ! derivative df(p), isentropic partial density, sum of volume fractions,
            ! mixture density, dynamic pressure, surface energy, specific heat ratio
            ! function, liquid stiffness function (two variations of the last two
            ! ones), shear and volume Reynolds numbers and the Weber numbers
            REAL(KIND(0d0))                                   ::  pres_relax
            REAL(KIND(0d0)), DIMENSION(num_fluids)            :: pres_K_init
            REAL(KIND(0d0))                                   ::   numerator
            REAL(KIND(0d0))                                   :: denominator
            REAL(KIND(0d0))                                   ::      drhodp
            REAL(KIND(0d0))                                   ::      f_pres
            REAL(KIND(0d0))                                   ::     df_pres
            REAL(KIND(0d0)), DIMENSION(num_fluids)            ::     rho_K_s
            REAL(KIND(0d0))                                   ::         rho
            REAL(KIND(0d0))                                   ::    dyn_pres
            REAL(KIND(0d0))                                   ::        E_We
            REAL(KIND(0d0))                                   ::       gamma
            REAL(KIND(0d0))                                   ::      pi_inf
            REAL(KIND(0d0)), DIMENSION(num_fluids)            ::   gamma_min
            REAL(KIND(0d0)), DIMENSION(num_fluids)            ::    pres_inf
            REAL(KIND(0d0)), DIMENSION(2)                     ::          Re
            REAL(KIND(0d0)), DIMENSION(num_fluids,num_fluids) ::          We

            ! Generic loop iterators
            INTEGER :: i,j,k,l,iter

            ! Relaxation procedure determination variable
            LOGICAL :: relax
            
            DO i = 1, num_fluids
                gamma_min(i) = 1d0/fluid_pp(i)%gamma + 1d0
                pres_inf(i)  = fluid_pp(i)%pi_inf / (1d0+fluid_pp(i)%gamma)
            END DO

            DO j = 0, m
                DO k = 0, n
                    DO l = 0, p
                        ! Numerical correction of the volume fractions
                        IF (mpp_lim) THEN
                            CALL s_mixture_volume_fraction_correction(q_cons_vf, j, k, l )
                        END IF
                        ! Pressures relaxation procedure ===================================
                        ! Is the pressure relaxation procedure necessary?
                        relax = .TRUE.
                        DO i = 1, num_fluids
                            IF (q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) .GT. (1d0-sgm_eps)) relax = .FALSE.
                        END DO
                        IF (relax) THEN
                            ! Initial state
                            pres_relax = 0.d0
                            DO i = 1, num_fluids
                                IF (q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) .GT. sgm_eps) THEN
                                    pres_K_init(i) = &
                                            ((q_cons_vf(i+internalEnergies_idx%beg-1)%sf(j,k,l) & 
                                            -q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)*fluid_pp(i)%qv) &
                                            /q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) &
                                            - fluid_pp(i)%pi_inf) & 
                                            /fluid_pp(i)%gamma
                                    IF (pres_K_init(i) .LE. -(1d0 - 1d-8)*pres_inf(i) + 1d-8) & 
                                        pres_K_init(i) = -(1d0 - 1d-8)*pres_inf(i) + 1d-0
                                    !TODO BE VERY CAREFUL ABOUT THE LIMIT HERE as it may be 1d0 instead
                                ELSE
                                    pres_K_init(i) = 0d0
                                END IF
                                pres_relax = pres_relax + q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l)*pres_K_init(i)
                            END DO

                            ! Iterative process for relaxed pressure determination
                            ! Convergence?
                            iter    = 0; f_pres  = 1.d-9; df_pres = 1d9
                            DO i = 1, num_fluids
                                rho_K_s(i) = 0d0
                            END DO
                            DO WHILE (DABS(f_pres) .GT. pnewtonk_eps)
                                pres_relax = pres_relax - f_pres/df_pres
                                iter = iter + 1
                                IF ( iter == pnewtonk_iter ) THEN
                                    PRINT '(A)', 'Pressure relaxation procedure failed to converge to a solution. Exiting ...'
                                    CALL s_mpi_abort()
                                END IF
                                ! Physical pressure?
                                DO i = 1, num_fluids
                                    IF (pres_relax .LE. -(1d0 - 1d-8)*pres_inf(i) + 1d-8) THEN
                                        pres_relax = -(1d0 - 1d-8)*pres_inf(i) + 1d0
                                    END IF
                                END DO
                                ! Newton-Raphson method
                                f_pres  = -1d0
                                df_pres = 0d0
                                DO i = 1, num_fluids
                                    IF (q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) .GT. sgm_eps) THEN
                                        numerator   = gamma_min(i)*(pres_relax+pres_inf(i))
                                        denominator = numerator + pres_K_init(i)-pres_relax
                                        rho_K_s(i)  = q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)/&
                                            MAX(q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l),sgm_eps)*& 
                                                      numerator/denominator
                                        drhodp      = q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l) / & 
                                            MAX(q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l),sgm_eps) * & 
                                            gamma_min(i)*(pres_K_init(i)+pres_inf(i)) / (denominator*denominator)
                                        f_pres      = f_pres  + q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l) / rho_K_s(i)
                                        df_pres     = df_pres - q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l) * & 
                                                      drhodp / (rho_K_s(i)*rho_K_s(i))
                                        !rho_K_s(i) = q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l) / &
                                        !            MAX(q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l),sgm_eps) &
                                        !            * ((pres_relax+pres_inf(i)) / (pres_K_init(i) + &
                                        !            pres_inf(i)))**(1d0/gamma_min(i))
                                        !f_pres   = f_pres + q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)/rho_K_s(i)
                                        !df_pres  = df_pres - q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l) &
                                        !           /(gamma_min(i)*rho_K_s(i)*(pres_relax+pres_inf(i)))
                                    END IF
                                END DO
                            END DO
                            ! Cell update of the volume fraction
                            DO i = 1, num_fluids
                                IF (q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) .GT. sgm_eps) & 
                                    q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) = & 
                                    q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l) / rho_K_s(i)
                            END DO
                        END IF
                        CALL s_mixture_total_energy_correction(q_cons_vf, j, k, l )
                    END DO
                END DO
            END DO

        END SUBROUTINE s_infinite_p_relaxation_k ! -----------------------

        !>  The purpose of this procedure is to infinitely relax
        !!      the pressures from the internal-energy equations to a
        !!      unique pressure, from which the corresponding volume
        !!      fraction of each phase are recomputed. For conservation
        !!      purpose, this pressure is finally corrected using the
        !!      mixture-total-energy equation.
        !!  @param q_cons_vf Cell-average conservative variables
        SUBROUTINE s_infinite_p_relaxation(q_cons_vf) ! ----------------

            TYPE(scalar_field), DIMENSION(sys_size), INTENT(INOUT) :: q_cons_vf 
            !> @name Relaxed pressure, initial partial pressures, function f(p) and its partial
            !! derivative df(p), isentropic partial density, sum of volume fractions,
            !! mixture density, dynamic pressure, surface energy, specific heat ratio
            !! function, liquid stiffness function (two variations of the last two
            !! ones), shear and volume Reynolds numbers and the Weber numbers
            !> @{
            REAL(KIND(0d0))                                   ::  pres_relax
            REAL(KIND(0d0))                                   :: rho, rhoeq_k
            REAL(KIND(0d0))                                   ::          a1
            REAL(KIND(0d0))                                   ::    dyn_pres
            REAL(KIND(0d0))                                   ::E_We, p_infk
            REAL(KIND(0d0))                                   ::gamma, pi_inf
            REAL(KIND(0d0)), DIMENSION(num_fluids)            ::p_k, alpha_k
            REAL(KIND(0d0)), DIMENSION(2)                     ::          Re
            REAL(KIND(0d0)), DIMENSION(num_fluids,num_fluids) ::          We
            !> @}
            INTEGER :: i,j,k,l           !< Generic loop iterators
            LOGICAL :: relax             !< Relaxation procedure determination variable
            DO j = 0, m
                DO k = 0, n
                    DO l = 0, p
                        ! P RELAXATION ==================================
                        relax = .FALSE.
                        IF (mpp_lim) THEN
                            CALL s_mixture_volume_fraction_correction(q_cons_vf, j, k, l )
                        END IF
                        IF ( (q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) .GT. palpha_eps ) .AND. &
                              q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) .LT. 1.d0-palpha_eps ) relax = .TRUE.
                        IF (relax) THEN
                            DO i = 1, num_fluids
                                 alpha_k(i) = q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) 
                                 rhoeq_k = (q_cons_vf(i+internalEnergies_idx%beg-1)%sf(j,k,l) & 
                                          -q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)*fluid_pp(i)%qv) &
                                          /q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) 
                                 p_k(i) = (rhoeq_k-fluid_pp(i)%pi_inf)/fluid_pp(i)%gamma
                            END DO
                            a1 = f_alpha1_prelax(p_k,alpha_k)
                            ! Cell update of the volume fraction
                            q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) = a1
                            q_cons_vf(2+adv_idx%beg-1)%sf(j,k,l) = 1.d0 - a1
                            CALL s_mixture_total_energy_correction(q_cons_vf, j, k, l )
                        END IF
                    END DO
                END DO
            END DO

        END SUBROUTINE s_infinite_p_relaxation ! ----------------

        !> Description: The purpose of this procedure is to infinitely relax
        !!              the pressures from the internal-energy equations to a
        !!              unique pressure, from which the corresponding volume
        !!              fraction of each phase are recomputed. For conservation
        !!              purpose, this pressure is finally corrected using the
        !!              mixture-total-energy equation.
        SUBROUTINE s_infinite_p_relaxation_k(q_cons_vf) ! ----------------        
            ! Relaxed pressure, initial partial pressures, function f(p) and its partial
            ! derivative df(p), isentropic partial density, sum of volume fractions,
            ! mixture density, dynamic pressure, surface energy, specific heat ratio
            ! function, liquid stiffness function (two variations of the last two
            ! ones), shear and volume Reynolds numbers and the Weber numbers
            ! Cell-average conservative variables
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(INOUT) :: q_cons_vf         
            REAL(KIND(0d0)), DIMENSION(num_fluids) ::     rho_K_s, pres_K_init
            REAL(KIND(0d0)), DIMENSION(num_fluids) ::   gamma_min, pres_inf
            ! Generic loop iterators
            INTEGER :: i,j,k,l
            ! Relaxation procedure determination variable
            LOGICAL :: relax
            DO i = 1, num_fluids
                gamma_min(i) = 1d0/fluid_pp(i)%gamma + 1d0
                pres_inf(i)  = fluid_pp(i)%pi_inf / (1d0+fluid_pp(i)%gamma)
            END DO
            DO j = 0, m
                DO k = 0, n
                    DO l = 0, p
                        ! Numerical correction of the volume fractions
                        IF (mpp_lim) THEN
                            CALL s_mixture_volume_fraction_correction(q_cons_vf, j, k, l )
                        END IF
                        ! Pressures relaxation procedure ===================================
                        relax = .TRUE.
                        DO i = 1, num_fluids
                            IF (q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) .GT. (1.d0-palpha_eps)) relax = .FALSE.
                        END DO
                        IF (relax) THEN
                            ! Calculating the initial pressure
                            DO i = 1, num_fluids
                              pres_K_init(i) = ((q_cons_vf(i+internalEnergies_idx%beg-1)%sf(j,k,l) & 
                                 -q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)*fluid_pp(i)%qv) &
                                 /q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) &
                                 -fluid_pp(i)%pi_inf)/fluid_pp(i)%gamma
                              IF (pres_K_init(i) .LE. -(1d0 - 1d-8)*pres_inf(i) + 1d-8) & 
                                    pres_K_init(i) = -(1d0 - 1d-8)*pres_inf(i) + 1d-0
                            END DO
                            CALL s_compute_p_relax_k(rho_K_s,gamma_min,pres_inf,pres_K_init,q_cons_vf,j,k,l)
                            ! Cell update of the volume fraction
                            DO i = 1, num_fluids
                                IF (q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) .GT. palpha_eps) & 
                                    q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) = & 
                                    q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l) / rho_K_s(i)
                            END DO
                        END IF
                        CALL s_mixture_total_energy_correction(q_cons_vf, j, k, l )
                    END DO
                END DO
            END DO
        END SUBROUTINE s_infinite_p_relaxation_k ! -----------------------

        !>  The purpose of this procedure is to infinitely relax
        !!      the pressures from the internal-energy equations to a
        !!      unique pressure, from which the corresponding volume
        !!      fraction of each phase are recomputed. For conservation
        !!      purpose, this pressure is finally corrected using the
        !!      mixture-total-energy equation.
        !!  @param q_cons_vf Cell-average conservative variables
        SUBROUTINE s_infinite_pt_relaxation(q_cons_vf) ! ----------------
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(INOUT) :: q_cons_vf 
            !> @name Relaxed pressure, initial partial pressures, function f(p) and its partial
            !! derivative df(p), isentropic partial density, sum of volume fractions,
            !! mixture density, dynamic pressure, surface energy, specific heat ratio
            !! function, liquid stiffness function (two variations of the last two
            !! ones), shear and volume Reynolds numbers and the Weber numbers
            !> @{
            REAL(KIND(0d0))                                   :: pres_relax, Trelax
            REAL(KIND(0d0)), DIMENSION(num_fluids)            :: p_k, alpha_k, Tk
            REAL(KIND(0d0))                                   :: rhoalpha1, rhoalpha2
            REAL(KIND(0d0))                                   :: rho, rhoe, rhoeq_k
            REAL(KIND(0d0))                                   :: rho1, rho2
            REAL(KIND(0d0))                                   :: a1, a2
            REAL(KIND(0d0))                                   :: dyn_pres
            REAL(KIND(0d0))                                   :: E_We
            REAL(KIND(0d0))                                   :: gamma
            REAL(KIND(0d0))                                   :: pi_inf, p_infk
            REAL(KIND(0d0))                                   :: pres_sat, Tsat
            REAL(KIND(0d0))                                   ::  A, B, C, D
            REAL(KIND(0d0)), DIMENSION(2)                     ::          Re
            REAL(KIND(0d0)), DIMENSION(num_fluids,num_fluids) ::          We
            !> @}
            INTEGER :: i, j, k, l        !< Generic loop iterators
            LOGICAL :: relax             !< Relaxation procedure determination variable
            !< Computing the constant saturation properties 

            DO j = 0, m
                DO k = 0, n
                    DO l = 0, p
                        ! P RELAXATION ==================================
                        relax = .FALSE.
                        IF (mpp_lim) THEN
                            CALL s_mixture_volume_fraction_correction(q_cons_vf, j, k, l )
                        END IF
                        IF ( (q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) .GT. palpha_eps ) .AND. &
                              q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) .LT. 1.d0-palpha_eps ) relax = .TRUE.
                        IF (relax) THEN
                            DO i = 1, num_fluids
                                 alpha_k(i) = q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) 
                                 rhoeq_k = (q_cons_vf(i+internalEnergies_idx%beg-1)%sf(j,k,l) & 
                                          -q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)*fluid_pp(i)%qv) &
                                          /q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) 
                                 p_k(i) = (rhoeq_k-fluid_pp(i)%pi_inf)/fluid_pp(i)%gamma
                            END DO
                            a1 = f_alpha1_prelax(p_k,alpha_k)
                            ! Cell update of the volume fraction
                            q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) = a1
                            q_cons_vf(2+adv_idx%beg-1)%sf(j,k,l) = 1.d0 - a1
                        END IF
                        CALL s_mixture_total_energy_correction(q_cons_vf, j, k, l )

                        ! PT RELAXATION ==================================
                        rhoe = 0.d0
                        relax = .FALSE.
                        IF (mpp_lim) THEN
                            CALL s_mixture_volume_fraction_correction(q_cons_vf, j, k, l )
                        END IF
                        IF ( (q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) .GT. palpha_eps ) .AND. &
                              q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) .LT. 1.d0-palpha_eps ) relax = .TRUE.
                        IF (relax) THEN
                            rhoalpha1 = q_cons_vf(cont_idx%beg)%sf(j,k,l)
                            rhoalpha2 = q_cons_vf(1+cont_idx%beg)%sf(j,k,l)
                            DO i = 1, num_fluids
                                rhoe = rhoe + q_cons_vf(i+internalEnergies_idx%beg-1)%sf(j,k,l)
                            END DO
                            a1 = f_alpha1_ptrelax(rhoalpha1,rhoalpha2,rhoe)
                            ! Cell update of the volume fraction
                            q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l)  = a1
                            q_cons_vf(2+adv_idx%beg-1)%sf(j,k,l)  = 1.d0 - a1
                        END IF
                        CALL s_mixture_total_energy_correction(q_cons_vf, j, k, l )
                    END DO
                END DO
            END DO
        END SUBROUTINE s_infinite_pt_relaxation ! -----------------------

        ! modification from Jose Rodolfo Chreim, this subroutine is based on the work of Flatten
        ! Morin, and Munkejord: On Solutions to Equilibrium Problems for Systems of Stiffened Gases
        ! routine to calculate pressure and temperature equilibrium
        SUBROUTINE s_infinite_pt_relaxation_k(q_cons_vf) ! ----------------           
            ! initializing variables
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(INOUT)  ::  q_cons_vf         
            REAL(KIND(0.0D0)), DIMENSION(num_fluids)                ::  g_min, p_inf, pkHS
            REAL(KIND(0.0D0)), DIMENSION(num_fluids)                ::  ark, ekS, rhokS
            REAL(KIND(0.0D0)), DIMENSION(num_fluids)                ::  akT, arkT
            REAL(KIND(0.0D0))                                       ::  gp, gpp, hp, pO, pS, TS
            REAL(KIND(0.0D0))                                       ::  rhoe, Ewe, mCP, mQ, TvF
            REAL(KIND(0.0D0))                                       ::  rho
            REAL(KIND(0.0D0))                                       ::  TvFT
            REAL(KIND(0.0D0))                                       ::  dynP
            ! Generic loop iterators
            INTEGER :: i, j, k, l, ns

            ! gamma min for the SG EoS
            g_min( 1 : num_fluids ) = 1.0D0 / fluid_pp( 1 : num_fluids )%gamma + 1.0D0

            ! p_infty for the SG EoS
            p_inf( 1 : num_fluids )  = fluid_pp( 1 : num_fluids )%pi_inf &
            / ( 1.0D0 + fluid_pp( 1 : num_fluids )%gamma )

            DO j = 0, m
				DO k = 0, n
					DO l = 0, p
						
                        ! correcting negative volume fraction values
                        IF (mpp_lim) THEN
                            CALL s_mixture_volume_fraction_correction(q_cons_vf, j, k, l )
                        END IF

                        ! zeroing energy so this quantity can be updated at each (j,k,l) location
                        mQ = 0.0D0
                        mCP = 0.0D0
                        pO = 0.0D0
                        TvF = 0.0D0
                        pS = 0.0D0
                        ! Elastic energy. This is here just for future reference, in case the model needs to be extended to 
                        ! account for elasticity
                        EWe = 0.0D0
                        rho = 0.0D0

                        DO i = 1, num_fluids

                            ! check if all alpha(i), and alpha(i)*rho(i) are positive, as this one of the condition for the solver to
                            ! find a unique solution.

                            ! constraint on ith fluid volume fraction. Look for stronger condition, if possible
                            IF ( q_cons_vf( i + adv_idx%beg - 1 )%sf( j, k, l ) .LT. 0.0D0 ) THEN
                                
                                PRINT *, 'Constraint on individual volume fractions not satisfied &
                                i.e. alpha <= 0. m_phase_change, s_infinite_pt_relaxation_k. &
                                Aborting.'
                               
                                PRINT *, 'i, j, k, l', i, j, k, l
                                
                                PRINT *, 'alpha_i', q_cons_vf( i + adv_idx%beg - 1 )%sf( j, k, l )
                                
                                CALL s_mpi_abort()

                            END IF

                            ! constraint on individual masses
                            IF ( q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l ) .LT. 0.0D0 ) THEN
                                
                                PRINT *, 'Constraint on individual partial densities not satisfied &
                                i.e. alpha*rho <= 0. m_phase_change, s_infinite_pt_relaxation_k'
                                
                                PRINT *, 'i, j, k, l', i, j, k, l
                                
                                PRINT *, 'alpha_rho_i', q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l )

                                CALL s_mpi_abort()

                            END IF

                            ! pressure that comes from the homogeneous solver at that cell
                            pkHS( i ) = ( ( q_cons_vf( i + internalEnergies_idx%beg - 1 )%sf( j, k, l ) & 
							    	- q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l ) * fluid_pp( i )%qv ) &	
								    / q_cons_vf( i + adv_idx%beg - 1 )%sf( j, k, l ) &
								    - fluid_pp( i )%pi_inf ) / fluid_pp( i )%gamma
                
                            ! sum of the total alpha*rho*q of the system
                            mQ = mQ + q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l ) * fluid_pp( i )%qv
                
                            ! sum of the total alpha*rho*cp of the system
                            mCP = mCP + q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l ) &
                            * fluid_pp( i )%cv * g_min( i )

                            ! Checking value of alpha*rho
                            ark( i ) = q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l )

                            ! Mixture density
                            rho = rho + q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l )

                            ! checking volume fraction alpha
                            akT( i ) = q_cons_vf( i + adv_idx%beg - 1 )%sf( j, k, l )

                            ! Total Volume Fraction
                            TvF = TvF + akT( i )

    					END DO

                        ! dynamic pressure so as to calculate the total internal energy as the total energy minus the 
                        ! kinetic energy
                        dynP = 0.0D0                  
                        DO i = mom_idx%beg, mom_idx%end
                            dynP = dynP + 5.0D-1 * q_cons_vf( i )%sf( j, k, l ) * & 
                            q_cons_vf( i )%sf( j, k, l ) / rho
                        END DO

                        ! calculating the total energy that MUST be preserved throughout the pT-relaxation procedure at each of the cells
                        rhoe = q_cons_vf( E_idx )%sf( j, k, l ) - dynP

                        ! Checking the constraints that must be satisfied for the solver to work
                        ! constraint on energy
                        IF ( ( rhoe - mQ - MINVAL( p_inf ) ) .LT. 0.0D0 ) THEN
                            PRINT *, 'Constraint on total energy not satisfied &
                            i.e. rhoe - mQ - MINVAL( p_inf ) <= 0. m_phase_change &
                            s_infinite_pt_relaxation_k. Aborting'

                            PRINT *, 'j, k, l', j, k,l

                            PRINT *, 'rhoe - mQ - MINVAL( p_inf )', rhoe - mQ - MINVAL( p_inf )

                            CALL s_mpi_abort()

                        END IF

                        ! calculating initial estimate for pressure in the pT-relaxation procedure. I will also use this variable to 
                        ! iterate over the Newton's solver
                        pO = 0
                        
                        ns = 0

                        ! Maybe improve this condition afterwards. As long as the initial guess is in between -min(p_inf)
                        ! and infinity, a solution should be able to be found.
                        pS = MAX( MAXVAL( pkHS ), MINVAL( p_inf ) ) + 1.0D4

                        ! Maybe this test is not necessary once I defined pS as above
                        IF ( pS .LE. -1.0D0 * MINVAL( p_inf ) ) THEN
                            PRINT *, 'proposed IC is not a good one. Trying a new one &
                            (m_phase_change, s_infinite_pt_relaxation_k). Aborting'

                            CALL s_mpi_abort()

                        END IF

                        ! Newton Solver for the pT-equilibrium
                        DO WHILE ( DABS( pS - pO ) .GT. palpha_eps .AND. DABS( ( pS - pO ) / pO ) .GT. palpha_eps )
    
                            ! updating old pressure 
                            pO = pS

                            ! updating functions used in the Newton's solver
                            gpp = 0.0D0
                            gp = 0.0D0
                            hp = 0.0D0
                            DO i = 1, num_fluids
                                    gp =   gp + ( g_min( i ) - 1.0D0 ) & 
                                    * q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l ) * fluid_pp( i )%cv &
                                    * ( rhoe + pS - mQ ) / ( mCP * ( pS + p_inf( i ) ) )

                                    gpp = gpp + ( g_min( i ) - 1.0D0 ) & 
                                    * q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l ) * fluid_pp( i )%cv &
                                    * ( p_inf( i ) - rhoe + mQ ) / ( mCP * ( pS + p_inf( i ) ) ** 2 )
                            END DO

                            ! fix this line, I don't think it works
                            hp = 1.0D0 / ( rhoe + pS - mQ ) + 1.0D0 / ( pS + MINVAL( p_inf ) )

                            ! updating common pressure for the newton solver
                            pS = pO + ( ( 1.0D0 - gp ) / gpp ) / ( 1.0D0 - ( 1.0D0 - gp + DABS( 1.0D0 - gp ) ) &
                            / ( 2.0D0 * gpp ) * hp )

                            ns = ns + 1
                            IF ( ( ns .GT. 100 ) .OR. ( pS .LE. -1.0D0 * MINVAL( p_inf ) ) ) THEN
                                    PRINT *, 'Newton Solver for the pT-relaxation solver failed &
                                    (m_phase_change, s_infinite_pt_relaxation_k). Aborting'

                                    CALL s_mpi_abort()

                            END IF

                        END DO

                        ! UPDATING Variables after the newton solver. A priori, the following variables must be updated: alpha, alpha*rho, and alpha*rho*e
                        ! now, given alpha*rho remains constant in the process and we do not need to know rho separately, alpha*rho does not need to be updated
                        ! I do need, however, to update alpha, and alpha*rho*e. So, in the end I update alpha, e, and then alpha*rho*e
                        ! rhoeT = 0.0D0
                        TvFT = 0.0D0

                        ! GmT2 = SUM( q_cons_vf( adv_idx%beg : num_fluids + adv_idx%beg - 1 )%sf( j, k, l ) & 
                        ! * fluid_pp( 1 : num_fluids )%gamma )

                        ! common temperature
                        TS = ( rhoe + pS - mQ ) / mCP

                        ! updating variables
                        ! densities
                        rhokS( 1 : num_fluids ) = ( pS + p_inf( 1 : num_fluids ) ) &
                        / ( ( g_min( 1 : num_fluids ) - 1 ) * fluid_pp( 1 : num_fluids )%cv * TS )

                        ! internal energy
                        ekS( 1 : num_fluids ) = ( pS + g_min( 1 : num_fluids ) & 
                        * p_inf( 1 : num_fluids ) ) / ( pS + p_inf( 1 : num_fluids ) ) & 
                        * fluid_pp( 1 : num_fluids )%cv * TS + fluid_pp( 1 : num_fluids )%qv

                        DO i = 1, num_fluids
                                                       
                            ! volume fractions
                            q_cons_vf( i + adv_idx%beg - 1 )%sf( j, k, l ) = &
                            q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l ) / rhokS( i ) 

                            ! alpha*rho*e
                            q_cons_vf( i + internalEnergies_idx%beg - 1 )%sf( j, k, l ) = & 
                            q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l ) * ekS( i )
                            
                            ! Total volume Fraction Test
                            TvFT = TvFT + q_cons_vf( i + adv_idx%beg - 1 )%sf( j, k, l )

                        END DO

                        ! testing if hypotheses of the pT-relazation model are satisfied. Uncomment this part to 
                        ! test. Note that neither the "Total Thermodynamic Energy test" or the "indiviual masses test" 
                        ! is necessary as none of the variables used to calculate those quantities are updated.

                        ! sum of volume fractions test
                        IF( DABS( TvF - TvFT ) .GT. 1.0D-2 ) THEN

                            PRINT *, 'Total Volume Fraction Residual after pT-equilibrium', TvF - TvFT
                            
                            PRINT *, 'j, k, l ', j, k, l
                            
                        END IF

                    END DO
                END DO
            END DO

        END SUBROUTINE s_infinite_pt_relaxation_k ! -----------------------

        !>  The purpose of this procedure is to infinitely relax
        !!      the pressures from the internal-energy equations to a
        !!      unique pressure, from which the corresponding volume
        !!      fraction of each phase are recomputed. For conservation
        !!      purpose, this pressure is finally corrected using the
        !!      mixture-total-energy equation.
        !!  @param q_cons_vf Cell-average conservative variables
        SUBROUTINE s_infinite_ptg_relaxation(q_cons_vf) ! ----------------
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(INOUT) :: q_cons_vf 
            !> @name Relaxed pressure, initial partial pressures, function f(p) and its partial
            !! function, liquid stiffness function (two variations of the last two
            !! ones), shear and volume Reynolds numbers and the Weber numbers
            !> @{
            REAL(KIND(0d0))                                   :: pres_relax, Trelax
            REAL(KIND(0d0)), DIMENSION(num_fluids)            :: p_k, alpha_k, Tk
            REAL(KIND(0d0))                                   :: rhoalpha1, rhoalpha2
            REAL(KIND(0d0))                                   :: rho, rhoe, rhoeq_k
            REAL(KIND(0d0))                                   :: rho1, rho2
            REAL(KIND(0d0))                                   :: a1, a2
            REAL(KIND(0d0))                                   :: gamma, pi_inf, p_infk
            REAL(KIND(0d0))                                   :: Tsat
            REAL(KIND(0d0)), DIMENSION(2)                     :: Re
            REAL(KIND(0d0)), DIMENSION(num_fluids,num_fluids) :: We
            !> @}
            INTEGER :: i, j, k, l        !< Generic loop iterators
            LOGICAL :: relax             !< Relaxation procedure determination variable
            !< Computing the constant saturation properties 

            DO j = 0, m
                DO k = 0, n
                    DO l = 0, p
                        ! P RELAXATION ========================================
                        relax = .FALSE.
                        IF (mpp_lim) THEN
                            CALL s_mixture_volume_fraction_correction(q_cons_vf, j, k, l )
                        END IF
                        CALL s_convert_to_mixture_variables( q_cons_vf, rho, &
                                                             gamma, pi_inf,  &
                                                             Re, We, j, k, l )
                        IF ( (q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) .GT. palpha_eps ) .AND. &
                              q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) .LT. 1.d0-palpha_eps ) relax = .TRUE.
                        IF (relax) THEN
                            DO i = 1, num_fluids
                                 alpha_k(i) = q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) 
                                 rhoeq_k = (q_cons_vf(i+internalEnergies_idx%beg-1)%sf(j,k,l) & 
                                          -q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)*fluid_pp(i)%qv) &
                                          /q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) 
                                 p_k(i) = (rhoeq_k-fluid_pp(i)%pi_inf)/fluid_pp(i)%gamma
                            END DO
                            a1 = f_alpha1_prelax(p_k,alpha_k)
                            ! Cell update of the volume fraction
                            q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) = a1
                            q_cons_vf(2+adv_idx%beg-1)%sf(j,k,l) = 1.d0 - a1
                        END IF
                        CALL s_mixture_total_energy_correction(q_cons_vf, j, k, l )
                        ! PT RELAXATION ==========================================
                        rhoe = 0.d0
                        relax = .FALSE.
                        IF (mpp_lim) THEN
                            CALL s_mixture_volume_fraction_correction(q_cons_vf, j, k, l )
                        END IF
                        CALL s_convert_to_mixture_variables( q_cons_vf, rho, &
                                                             gamma, pi_inf,  &
                                                             Re, We, j, k, l )
                        IF ( (q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) .GT. palpha_eps ) .AND. &
                              q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) .LT. 1.d0-palpha_eps ) relax = .TRUE.
                        IF (relax) THEN
                            rhoalpha1 = q_cons_vf(cont_idx%beg)%sf(j,k,l)
                            rhoalpha2 = q_cons_vf(1+cont_idx%beg)%sf(j,k,l)
                            DO i = 1, num_fluids
                                rhoe = rhoe + q_cons_vf(i+internalEnergies_idx%beg-1)%sf(j,k,l)
                            END DO
                            a1 = f_alpha1_ptrelax(rhoalpha1,rhoalpha2,rhoe)
                            ! Cell update of the volume fraction
                            q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l)  = a1
                            q_cons_vf(2+adv_idx%beg-1)%sf(j,k,l)  = 1.d0 - a1
                        END IF
                        CALL s_mixture_total_energy_correction(q_cons_vf, j, k, l )
                        ! CHECKING IF PTG RELAXATION IS NEEDED  =====================
                        rhoe = 0.d0
                        relax = .FALSE.
                        IF (mpp_lim) THEN
                            CALL s_mixture_volume_fraction_correction(q_cons_vf, j, k, l )
                        END IF
                        CALL s_convert_to_mixture_variables( q_cons_vf, rho, &
                                                             gamma, pi_inf,  &
                                                             Re, We, j, k, l )
                        IF ((q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) .GT. ptgalpha_eps ) .AND. &
                             q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) .LT. 1.d0-ptgalpha_eps ) relax = .TRUE.
                        IF (relax) THEN
                           DO i = 1, num_fluids
                               rhoe = rhoe + q_cons_vf(i+internalEnergies_idx%beg-1)%sf(j,k,l) 
                           END DO                   
                           pres_relax = (rhoe - pi_inf)/gamma
                           Tsat = f_Tsat(pres_relax)
                           DO i = 1, num_fluids
                             Tk(i) = ((q_cons_vf(i+internalEnergies_idx%beg-1)%sf(j,k,l) & 
                                    -q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)*fluid_pp(i)%qv) &
                                    /q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) &
                                    -fluid_pp(i)%pi_inf & 
                                    /(1.d0+fluid_pp(i)%gamma)) &
                                    /(q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)*fluid_pp(i)%cv &
                                    /q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l)) 
                           END DO
                           IF (Tk(1) .LT. Tsat) relax = .FALSE.
                        END IF
                        ! PTG RELAXATION PROCEDURE ===========================
                        IF (relax) THEN
                            CALL s_compute_ptg_pTrelax(pres_relax,Trelax,rho,rhoe)
                            p_infk = fluid_pp(1)%pi_inf/(1.d0+fluid_pp(1)%gamma)
                            rho1 = (pres_relax + p_infk)*fluid_pp(1)%gamma /& 
                                   (fluid_pp(1)%cv*Trelax)
                            p_infk = fluid_pp(2)%pi_inf/(1.d0+fluid_pp(2)%gamma)
                            rho2 = (pres_relax + p_infk)*fluid_pp(2)%gamma /& 
                                   (fluid_pp(2)%cv*Trelax)
                            ! Calculate vapor and liquid volume fractions
                            a1 = (rho-rho2)/(rho1-rho2)
                            a2 = 1.d0 - a1
                            ! Cell update of the volume fraction
                            q_cons_vf(cont_idx%beg)%sf(j,k,l)   = rho1*a1
                            q_cons_vf(1+cont_idx%beg)%sf(j,k,l) = rho2*a2
                            q_cons_vf(adv_idx%beg)%sf(j,k,l)    = a1
                            q_cons_vf(1+adv_idx%beg)%sf(j,k,l)  = a2
                        END IF
                        CALL s_mixture_total_energy_correction(q_cons_vf, j, k, l )
                    END DO
                END DO
            END DO
        END SUBROUTINE s_infinite_ptg_relaxation ! -----------------------

        ! modification from Jose Rodolfo Chreim, this subroutine is a modification of the work of Flatten
        ! Morin, and Munkejord: On Solutions to Equilibrium Problems for Systems of Stiffened Gases
        ! routine to calculate pressure and temperature equilibrium. Here, I extend their approach to 
        ! consider pTg equilibrium allowing m1 and p to vary
        SUBROUTINE s_infinite_ptg_relaxation_k(q_cons_vf) ! ----------------

            TYPE(scalar_field), DIMENSION(sys_size), INTENT(INOUT)  :: q_cons_vf 
            REAL(KIND(0.0D0)), DIMENSION(num_fluids)                :: g_min, p_inf
            REAL(KIND(0.0D0)), DIMENSION(num_fluids)                :: pkHS, TkHS, sk, hk, gk
            REAL(KIND(0.0D0)), DIMENSION(num_fluids)                :: ekS, rhokS
            REAL(KIND(0.0D0)), DIMENSION(2)                         :: R2D, DeltamP
            REAL(KIND(0.0D0)), DIMENSION(2,2)                       :: Jac, InvJac, TJac
            REAL(KIND(0.0D0))                                       :: A, B, C, D, R, Om
            REAL(KIND(0.0D0))                                       :: gpp, hp, pO, pS, TS
            REAL(KIND(0.0D0))                                       :: rhoe, Ewe, mCP, mQ, TvF
            REAL(KIND(0.0D0))                                       :: rho, mCVGP, mCVGP2, mCPD, mQD
            REAL(KIND(0.0D0))                                       :: TvFT
            REAL(KIND(0.0D0))                                       :: dynP
            !< Generic loop iterators
            INTEGER :: i, j, k, l
            INTEGER :: lp, vp, ns

            ! indices for the reacting liquid and gas
            lp = 1 
            vp = 2

            ! gamma min for the SG EoS
            g_min( 1 : num_fluids ) = 1.0D0 / fluid_pp( 1 : num_fluids )%gamma + 1.0D0

            ! p_infty for the SG EoS
            p_inf( 1 : num_fluids )  = fluid_pp( 1 : num_fluids )%pi_inf &
            / ( 1.0D0 + fluid_pp( 1 : num_fluids )%gamma )
            
            ! variables used in the calculation of the saturation curves for fluids 1 and 2
            A = ( g_min( lp ) * fluid_pp( lp )%cv  - g_min( vp ) * fluid_pp( vp )%cv &
            + fluid_pp( vp )%qvp - fluid_pp( lp )%qvp ) &
            / ( ( g_min( vp ) - 1.0D0 ) * fluid_pp( vp )%cv )
            
            B = ( fluid_pp( lp )%qv - fluid_pp( vp )%qv ) &
            / ( ( g_min( vp ) - 1.0D0 ) * fluid_pp( vp )%cv )
            
            C = ( g_min( vp ) * fluid_pp( vp )%cv - g_min( lp ) * fluid_pp( lp )%cv ) & 
            / ( ( g_min( vp ) - 1.0D0 ) * fluid_pp( vp )%cv )
            
            D = ( ( g_min( lp ) - 1.0D0 ) * fluid_pp( lp )%cv ) &
            / ( ( g_min( vp ) - 1.0D0 ) * fluid_pp( vp )%cv )

            DO j = 0, m
                DO k = 0, n
                    DO l = 0, p
                        ! Numerical correction of the volume fractions
                        IF (mpp_lim) THEN
                            CALL s_mixture_volume_fraction_correction( q_cons_vf, j, k, l )
                        END IF

                        ! calculating fluids' thermodynamic quantities after the homogeneous solver
                        rho = 0.0D0
                        TvF = 0.0D0
                        mCVGP = 0.0D0
                        mCVGP2 = 0.0D0
                        mCP = 0.0D0
                        mCPD = 0.0D0
                        mQ = 0.0D0
                        mQD = 0.0D0
                        DO i = 1, num_fluids

                            ! constraint on ith fluid volume fraction. Look for stronger condition, if possible
                            IF ( ( q_cons_vf( i + adv_idx%beg - 1 )%sf( j, k, l ) .LE. 0.0D0 ) .OR. &
                                 ( q_cons_vf( i + adv_idx%beg - 1 )%sf( j, k, l ) .GT. 1.0D0 )     ) THEN
                                    
                                PRINT *, 'Constraint on individual volume fractions not satisfied &
                                i.e. alpha <= 0, or alpha > 1. m_phase_change, s_infinite_ptg_relaxation_k. &
                                Aborting.'
                            
                                PRINT *, 'i, j, k, l', i, j, k, l

                                PRINT *, 'alpha_i', q_cons_vf( i + adv_idx%beg - 1 )%sf( j, k, l )
                                
                                PRINT *, 'alpha_rho_i', q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l )

                                CALL s_mpi_abort()

                            END IF

                            ! Temperature
                            TkHS( i ) = ( ( q_cons_vf( i + internalEnergies_idx%beg -1 )%sf( j, k ,l ) & 
                                    - q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l ) * fluid_pp( i )%qv ) &
                                    / q_cons_vf( i + adv_idx%beg - 1)%sf( j, k, l ) - fluid_pp( i )%pi_inf & 
                                    / ( 1.0D0 + fluid_pp( i )%gamma ) ) &
                                    / ( q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l ) * fluid_pp( i )%cv &
                                    / q_cons_vf( i + adv_idx%beg - 1 )%sf( j, k, l ) ) 
                            ! pressure
                            pkHS( i ) = ( ( q_cons_vf( i + internalEnergies_idx%beg - 1 )%sf( j, k, l ) & 
                                -q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l ) * fluid_pp(i)%qv ) &
                                / q_cons_vf( i + adv_idx%beg - 1 )%sf( j, k, l ) - fluid_pp( i )%pi_inf ) &
                                / fluid_pp( i )%gamma

                            ! Mixture density
                            rho = rho + q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l )

                            ! checking volume fraction alpha
                            ! akT( i ) = q_cons_vf( i + adv_idx%beg - 1 )%sf( j, k, l )

                            ! Total Volume Fraction
                            TvF = TvF + q_cons_vf( i + adv_idx%beg - 1 )%sf( j, k, l )

                            ! sum of the total alpha*rho*q of the system
                            mQ = mQ + q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l ) * fluid_pp( i )%qv
            
                            ! sum of the total alpha*rho*cp of the system
                            mCP = mCP + q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l ) &
                            * fluid_pp( i )%cv * g_min( i )

                            IF ( ( i .NE. lp ) .AND. ( i .NE. vp ) ) THEN

                                mCVGP = mCVGP + q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l ) &
                                * fluid_pp( i )%cv * ( g_min( i ) - 1 ) / ( pS + p_inf( i ) )

                                mCVGP2 = mCVGP2 + q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l ) &
                                * fluid_pp( i )%cv * ( g_min( i ) - 1 ) / ( ( pS + p_inf( i ) ) ** 2 )

                                mQD = mQD + q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l ) * fluid_pp( i )%qv

                                ! sum of the total alpha*rho*cp of the system
                                mCPD = mCPD + q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l ) * fluid_pp( i )%cv &
                                * g_min( i )

                            END IF

                        END DO

                        ! Maybe improve this condition afterwards.
                        pS = MAX( MINVAL( pkHS ), MINVAL( p_inf ) ) + 0.0D4

                        ! dynamic pressure so as to calculate the total internal energy as the total energy minus the 
                        ! kinetic energy. Note that in this process, the total internal energy is not expected to change,
                        ! therefore, at least for now, nor the kinetic energy (even though the 'masses' change). This is 
                        ! because E_thermo = E_total - E_kin, in this module
                        dynP = 0.0D0                  
                        DO i = mom_idx%beg, mom_idx%end

                            dynP = dynP + 5.0D-1 * q_cons_vf( i )%sf( j, k, l ) * & 
                            q_cons_vf( i )%sf( j, k, l ) / rho

                        END DO

                        ! total thermodynamic energy
                        rhoe = q_cons_vf( E_idx )%sf( j, k, l ) - dynP

                        ! checking condition on the total thermodynamic energy such that the pT-equilibrium (attention: not pTg)
                        ! can be satisfied
                        IF ( ( rhoe - mQ - MINVAL( p_inf ) ) .LT. 0.0D0 ) THEN
                            PRINT *, 'Constraint on total energy not satisfied &
                            i.e. rhoe - mQ - MINVAL( p_inf ) <= 0. m_phase_change &
                            s_infinite_pt_relaxation_k'

                            PRINT *, 'j, k, l', j, k,l

                            PRINT *, 'rhoe - mQ - MINVAL( p_inf )', rhoe - mQ - MINVAL( p_inf )

                        END IF

                        ! calculating residual, i.e. the difference between the Gibbs Free energy of the gas and the liquid
                        ! this will be the condition to which pTg will be activated or not.
                        ! Total Internal Energy, and Gibbs-Free energy
                        CALL s_compute_pTg_residue(A, B, C, D, g_min, j, k, l, lp, mCPD, mCVGP, mQD &
                        , p_inf, q_cons_vf, pS, rhoe, R2D, vp)

                        ! Check whether I need to use both absolute and relative values
                        ! for the residual, and how to do it adequately.
                        ! calculating pT-equilibrium as the IC for the pTg-equilibrium. Note that this step might not be necessary
                        ! and I am not even sure if this would somehow complicate finding the solution of the pTg-equilibrium.
                        ! IF ( DSQRT( R2D( 1 ) * R2D( 1 ) + R2D( 2 ) * R2D( 2 ) ) .GT. ptgalpha_eps ) THEN
                        !     CALL s_infinite_pt_relaxation_k( q_cons_vf )
                        ! END IF

                        ! pTg-equilibrium solution procedure
                        ! Newton Sover parameters
                        ! counter
                        ns = 0

                        ! Relaxation factor
                        Om = 1.0D-14

                        ! PRINT *, q_cons_vf( lp + cont_idx%beg - 1 )%sf( j, k, l ) &
                        ! , q_cons_vf( vp + cont_idx%beg - 1 )%sf( j, k, l ) 
                        ! PRINT *, q_cons_vf( lp + cont_idx%beg - 1 )%sf( j, k, l ) &
                        ! + q_cons_vf( vp + cont_idx%beg - 1 )%sf( j, k, l ) 

                        ! checking if pTg relaxation is needed.
                        DO WHILE ( DSQRT( R2D( 1 ) * R2D( 1 ) + R2D( 2 ) * R2D( 2 ) ) .GT. ptgalpha_eps )

                            ! Updating counter for the iterative procedure
                            ns = ns + 1

                            DO i = 1, num_fluids

                                ! constraint on individual masses, note that this is needed as the masses are updated,
                                ! and they can get negative
                                IF ( q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l ) .LT. 0.0D0 ) THEN
                                    
                                    PRINT *, 'Constraint on individual partial densities not satisfied &
                                    i.e. alpha*rho <= 0. m_phase_change, s_infinite_pt_relaxation_k'
                                    
                                    PRINT *, 'i, j, k, l', i, j, k, l
                                                                       
                                    PRINT *, 'alpha_rho_i', q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l )
                                    
                                    PRINT *, 'pkHS', pkHS( i )

                                    PRINT *, 'pS', pS

                                    PRINT *, 'TS', ( rhoe + pS - mQ ) / mCP

                                    PRINT *, 'rhoe', rhoe

                                    PRINT *, 'sum(alpha)', TvF

                                    PRINT *, 'ns', ns

                                    PRINT *, DSQRT( R2D( 1 ) * R2D( 1 ) + R2D( 2 ) * R2D( 2 ) )

                                    PRINT *, 'R2D(1)', R2D( 1 )

                                    PRINT *, 'R2D(2)', R2D( 2 )
                                    
                                    PRINT *, 'J', Jac, 'J-1', InvJac

                                    CALL s_mpi_abort()

                                END IF

                            END DO

                            ! Checking criteria for the solver to find a solution
                            ! constrain on the individual masses
                            IF ( ( ns .GT. 1E5 ) .OR. ( pS .LE. -1.0D0 * MINVAL( p_inf ) ) ) THEN

                                PRINT *, 'Newton Solver for the pTg-relaxation solver failed &
                                (m_phase_change, s_infinite_ptg_relaxation_k). Aborting'

                                PRINT *, 'absolute value of residual is'
                                
                                PRINT *, DSQRT( R2D( 1 ) * R2D( 1 ) + R2D( 2 ) * R2D( 2 ) )

                                PRINT *, 'R2D(1)', R2D( 1 )

                                PRINT *, 'R2D(2)', R2D( 2 )
                                
                                PRINT *, 'J', Jac, 'J-1', InvJac

                                PRINT *, 'j, k, l', j, k, l

                                PRINT *, 'ns', ns
                                
                                PRINT *, 'pS', pS

                                PRINT *, 'Dm', DeltamP

                                CALL s_mpi_abort()

                            END IF    

                            ! checking condition on the total thermodynamic energy such that the pT-equilibrium (attention: not pTg)
                            ! can be satisfied. Note that because the masses are update, so mQ and, therefore, this constraint need to be tested
                            ! at every iteration
                            IF ( ( rhoe - mQ - MINVAL( p_inf ) ) .LT. 0.0D0 ) THEN

                                PRINT *, 'Constraint on total energy not satisfied &
                                i.e. rhoe - mQ - MINVAL( p_inf ) <= 0. m_phase_change &
                                s_infinite_pt_relaxation_k'

                                PRINT *, 'j, k, l', j, k, l

                                PRINT *, 'rhoe - mQ - MINVAL( p_inf )', rhoe - mQ - MINVAL( p_inf )

                            END IF

                            ! Calculating Jacobian matrix and its inverse, analitically
                            CALL s_compute_jacobian_matrix( B, C, D, rhoe, g_min, InvJac, j, Jac, k, l, &
                                                lp, mCPD, mCVGP, mCVGP2, mQD, p_inf, pS, q_cons_vf, TJac, vp )

                            ! calculating correction array for Netown's method
                            ! DeltamP = - 1.0D0 * MATMUL( InvJac, R2D )
                            DeltamP = - 1.0D0 * MATMUL( TJac, R2D )

                            ! updating two reacting 'masses'. Recall that the other 'masses' do not change during the phase change
                            ! liquid
                            q_cons_vf( lp + cont_idx%beg - 1 )%sf( j, k, l ) = &
                            q_cons_vf( lp + cont_idx%beg - 1 )%sf( j, k, l ) + Om * DeltamP( 1 )

                            ! gas
                            q_cons_vf( vp + cont_idx%beg - 1  )%sf( j, k, l ) = &
                            q_cons_vf( vp + cont_idx%beg - 1 )%sf( j, k, l ) - Om * DeltamP( 1 )

                            ! updating pressure
                            pS = pS + Om * DeltamP( 2 ) 

                            ! Auxiliary variable to help in the calculation of the residue
                            mCVGP = 0.0D0
                            mCVGP2 = 0.0D0
                            mCP = 0.0D0
                            mCPD = 0.0D0
                            mQ = 0.0D0
                            mQD = 0.0D0

                            DO i = 1, num_fluids

                                ! sum of the total alpha*rho*q of the system
                                mQ = mQ + q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l ) * fluid_pp( i )%qv
                
                                ! sum of the total alpha*rho*cp of the system
                                mCP = mCP + q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l ) &
                                * fluid_pp( i )%cv * g_min( i )

                                IF ( ( i .NE. lp ) .AND. ( i .NE. vp ) ) THEN

                                    mCVGP = mCVGP + q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l ) &
                                    * fluid_pp( i )%cv * ( g_min( i ) - 1 ) / ( pS + p_inf( i ) )

                                    mCVGP2 = mCVGP2 + q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l ) &
                                    * fluid_pp( i )%cv * ( g_min( i ) - 1 ) / ( ( pS + p_inf( i ) ) ** 2 )

                                    mQD = mQD + q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l ) * fluid_pp( i )%qv

                                    ! sum of the total alpha*rho*cp of the system
                                    mCPD = mCPD + q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l ) * fluid_pp( i )%cv &
                                    * g_min( i )

                                END IF

                            END DO
                        
                            ! calculating residualof the pTg-equilibrium procedure
                            ! Total Internal Energy, and Gibbs-Free energy
                            CALL s_compute_pTg_residue(A, B, C, D, g_min, j, k, l, lp, mCPD, mCVGP, mQD &
                                                            , p_inf, q_cons_vf, pS, rhoe, R2D, vp)

                            PRINT *, 'jacobianT', TJac

                            PRINT *, 'liquid', q_cons_vf( lp + cont_idx%beg - 1 )%sf( j, k, l )

                            PRINT *, 'vapor', q_cons_vf( lp + cont_idx%beg - 1 )%sf( j, k, l )

                            PRINT *, 'pS', pS

                            PRINT *, DSQRT( R2D( 1 ) * R2D( 1 ) + R2D( 2 ) * R2D( 2 ) )

                        END DO

                        ! updating variables as in the pT-equilibrium procedure, once the pTg-equilibrium has been achieved
                        ! common temperature
                        TS = ( rhoe + pS - mQ ) / mCP
                        
                        ! entropy
                        sk( 1 : num_fluids ) = &
                        fluid_pp( 1 : num_fluids )%cv * DLOG( ( TS ** g_min( 1 : num_fluids ) ) &
                        / ( ( pS + p_inf( 1 : num_fluids ) ) ** ( g_min( 1 : num_fluids ) - 1.0D0 ) ) ) & 
                        + fluid_pp( 1 : num_fluids )%qvp

                        ! enthalpy
                        hk( 1 : num_fluids ) = g_min( 1 : num_fluids ) * fluid_pp( 1 : num_fluids )%cv &
                        * TS + fluid_pp( 1 : num_fluids )%qv

                        ! Gibbs-free energy
                        gk( 1 : num_fluids ) = hk( 1 : num_fluids ) - TS * sk( 1 : num_fluids )

                        ! densities
                        rhokS( 1 : num_fluids ) = ( pS + p_inf( 1 : num_fluids ) ) &
                        / ( ( g_min( 1 : num_fluids ) - 1 ) * fluid_pp( 1 : num_fluids )%cv * TS )

                        ! internal energy
                        ekS( 1 : num_fluids ) = ( pS + g_min( 1 : num_fluids ) & 
                        * p_inf( 1 : num_fluids ) ) / ( pS + p_inf( 1 : num_fluids ) ) & 
                        * fluid_pp( 1 : num_fluids )%cv * TS + fluid_pp( 1 : num_fluids )%qv

                        TvFT = 0.0D0
                        DO i = 1, num_fluids
                                                       
                            ! volume fractions
                            q_cons_vf( i + adv_idx%beg - 1 )%sf( j, k, l ) = &
                            q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l ) / rhokS( i ) 

                            q_cons_vf( i + adv_idx%beg - 1 )%sf( j, k, l ) = ( g_min( i ) - 1.0D0 ) & 
                            * q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l ) & 
                            * fluid_pp( i )%cv * ( rhoe + pS - mQ ) / ( mCP * ( pS + p_inf( i ) ) ) 

                            ! alpha * rho test
                            ! arkT( i ) = rhokS( i ) * q_cons_vf( i + adv_idx%beg - 1 )%sf( j, k, l )

                            ! alpha*rho*e
                            q_cons_vf( i + internalEnergies_idx%beg - 1 )%sf( j, k, l ) = & 
                            q_cons_vf( i + cont_idx%beg - 1 )%sf( j, k, l ) * ekS( i )
                                                        
                            ! Total volume Fraction Test
                            TvFT = TvFT + q_cons_vf( i + adv_idx%beg - 1 )%sf( j, k, l )

                        END DO

                        ! sum of volume fractions test
                        IF( DABS( TvF - TvFT ) .GT. 1.0D-2 ) THEN

                            PRINT *, 'Total Volume Fraction (initial)', TvF
                            PRINT *, 'Total Volume Fraction (final)', TvFT
                            PRINT *, 'Total Volume Fraction Residual after pTg-equilibrium', TvF - TvFT
                            
                            PRINT *, 'j, k, l ', j, k, l

                            PRINT *, q_cons_vf( lp + cont_idx%beg - 1 )%sf( j, k, l ) &
                            , q_cons_vf( vp + cont_idx%beg - 1 )%sf( j, k, l ) 
                            PRINT *, q_cons_vf( lp + cont_idx%beg - 1 )%sf( j, k, l ) &
                            + q_cons_vf( vp + cont_idx%beg - 1 )%sf( j, k, l ) 

                            PRINT *, 'pS', pS, 'TS', TS

                            PRINT *, 'rhoe', rhoe

                            PRINT *, R2D
                
                            PRINT *, gK(vp) - gK(lp)

                            PAUSE
                        END IF

                    END DO
                END DO
            END DO
        END SUBROUTINE s_infinite_ptg_relaxation_k ! -----------------------

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!! SUBROUTINES SUBROUTINES SUBROUTINES !!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !> @name Relaxed pressure, initial partial pressures, function f(p) and its partial
        !! derivative df(p), isentropic partial density, sum of volume fractions,
        !! mixture density, dynamic pressure, surface energy, specific heat ratio
        !! function, liquid stiffness function (two variations of the last two
        !! ones), shear and volume Reynolds numbers and the Weber numbers
        !> @{
        SUBROUTINE s_mixture_volume_fraction_correction(q_cons_vf, j, k, l )

            TYPE(scalar_field), DIMENSION(sys_size), INTENT(INOUT) :: q_cons_vf 
            INTEGER, INTENT(IN)                                    :: j, k, l
            REAL(KIND(0d0))                                        :: sum_alpha
            !> @}
            INTEGER :: i           !< Generic loop iterators
            sum_alpha = 0d0
            DO i = 1, num_fluids
                ! I changed this line to .LE. so the phase cannot 'disappear'.
                ! Think of more appropriate conditions, later on
               IF ((q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l) .LT. sgm_eps) .OR. &
                   (q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) .LT. sgm_eps)) THEN
                    q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l) = sgm_eps
                    q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l)  = sgm_eps
                    q_cons_vf(i+internalEnergies_idx%beg-1)%sf(j,k,l) = 0d0
               END IF
               IF (q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) .GT. 1d0) & 
                   q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) = 1d0
               sum_alpha = sum_alpha + q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l)
             END DO
             DO i = 1, num_fluids
                   q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) = & 
                   q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) / sum_alpha
             END DO

        END SUBROUTINE s_mixture_volume_fraction_correction

        !> @name Relaxed pressure, initial partial pressures, function f(p) and its partial
        !! derivative df(p), isentropic partial density, sum of volume fractions,
        !! mixture density, dynamic pressure, surface energy, specific heat ratio
        !! function, liquid stiffness function (two variations of the last two
        !! ones), shear and volume Reynolds numbers and the Weber numbers
        !> @{
        SUBROUTINE s_mixture_total_energy_correction(q_cons_vf, j, k, l )

            TYPE(scalar_field), DIMENSION(sys_size), INTENT(INOUT) :: q_cons_vf 
            INTEGER, INTENT(IN)                                    :: j, k, l
            !> @name Relaxed pressure, initial partial pressures, function f(p) and its partial
            !! derivative df(p), isentropic partial density, sum of volume fractions,
            !! mixture density, dynamic pressure, surface energy, specific heat ratio
            !! function, liquid stiffness function (two variations of the last two
            !! ones), shear and volume Reynolds numbers and the Weber numbers
            !> @{
            REAL(KIND(0d0))                                   :: rho, dyn_pres, E_We
            REAL(KIND(0d0))                                   :: gamma, pi_inf, pres_relax
            REAL(KIND(0d0))                                   :: prT, rhoeT
            REAL(KIND(0d0)), DIMENSION(2)                     :: Re
            REAL(KIND(0d0)), DIMENSION(num_fluids,num_fluids) :: We
            !> @}
            INTEGER :: i     !< Generic loop iterators
            CALL s_convert_to_mixture_variables( q_cons_vf, rho, &
                                                 gamma, pi_inf,  &
                                                 Re, We, j, k, l )
            dyn_pres = 0.d0
            DO i = mom_idx%beg, mom_idx%end
                ! PRINT *, 'i from s_mixture_tec', i
                dyn_pres = dyn_pres + 5d-1*q_cons_vf(i)%sf(j,k,l) * & 
                q_cons_vf(i)%sf(j,k,l) / MAX(rho,sgm_eps)
            END DO
            E_We = 0.d0
            pres_relax = (q_cons_vf(E_idx)%sf(j,k,l) - dyn_pres - pi_inf - E_We)/gamma
            
            rhoeT = 0.0D0
            DO i = 1, num_fluids
                rhoeT = rhoeT + q_cons_vf(i+internalEnergies_idx%beg-1)%sf(j,k,l)
            END DO
            
            pRT = (rhoeT - pi_inf) / gamma

            ! PRINT *, 'P_RELAX', pres_relax
            ! PRINT *, 'pRT', pRT
            DO i = 1, num_fluids
                ! PRINT *, i
                ! PRINT *, q_cons_vf(i+internalEnergies_idx%beg-1)%sf(j,k,l)
                q_cons_vf(i+internalEnergies_idx%beg-1)%sf(j,k,l) = & 
                q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) * & 
                (fluid_pp(i)%gamma*pres_relax + fluid_pp(i)%pi_inf) +  & 
                q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)*fluid_pp(i)%qv
                ! PRINT *, q_cons_vf(i+internalEnergies_idx%beg-1)%sf(j,k,l)
            END DO
            ! ==================================================================
        END SUBROUTINE s_mixture_total_energy_correction

        !>     The purpose of this subroutine is to determine the saturation
        !!         temperature by using a Newton-Raphson method from the provided
        !!         equilibrium pressure and EoS of the binary phase system.
        !!     @param q_cons_vf Cell-average conservative variables
        !!     @param p_star equilibrium pressure at the interface    
        FUNCTION f_alpha1_prelax(p_k,alpha_k)
            REAL(KIND(0d0))                                       ::  pstar, f_alpha1_prelax
            REAL(KIND(0d0)), DIMENSION(num_fluids), INTENT(IN)    ::  p_k, alpha_k
            !> @name In-subroutine variables: vapor and liquid material properties n, p_infinity
            !!       heat capacities, cv, reference energy per unit mass, q, coefficients for the
            !!       iteration procedure, A-D, and iteration variables, f and df
            !> @{
            REAL(KIND(0d0))                                       ::  Z1, Z2, pI, C1, C2
            REAL(KIND(0d0))                                       ::  ap, bp, dp

            ! Calculating coefficients, Eq. C.6, Pelanti 2014
            Z1 = n1*(p_k(1)+pinf1)
            Z2 = n2*(p_k(2)+pinf2)
            pI = (Z2*p_k(1)+Z1*p_k(2))/(Z1+Z2)
            C1 = 2.d0*n1*pinf1+(n1-1.d0)*p_k(1)
            C2 = 2.d0*n2*pinf2+(n2-1.d0)*p_k(2)
            ap = 1.d0 + n2*alpha_k(1) + n1*alpha_k(2)
            bp = C1*alpha_k(2)+C2*alpha_k(1)-(n2+1.d0)*alpha_k(1)*p_k(1)-(n1+1.d0)*alpha_k(2)*p_k(2)
            dp = -(C2*alpha_k(1)*p_k(1) + C1*alpha_k(2)*p_k(2))
            ! Calculating the Tstar temperature, Eq. C.7, Pelanti 2014
            pstar = (-bp + DSQRT(bp*bp - 4.d0*ap*dp))/(2.d0*ap)
            f_alpha1_prelax = alpha_k(1)*((n1-1.d0)*pstar + 2.d0*p_k(1) + C1)/((n1+1.d0)*pstar+C1)
        END FUNCTION f_alpha1_prelax

        !>     The purpose of this subroutine is to determine the saturation
        !!         temperature by using a Newton-Raphson method from the provided
        !!         equilibrium pressure and EoS of the binary phase system.
        !!     @param q_cons_vf Cell-average conservative variables
        !!     @param p_star equilibrium pressure at the interface    
        FUNCTION f_alpha1_ptrelax(rhoalpha1,rhoalpha2,E0)

            REAL(KIND(0d0))                :: pstar, f_alpha1_ptrelax
            REAL(KIND(0d0)), INTENT(IN)    :: rhoalpha1, rhoalpha2, E0
            !> @name In-subroutine variables: vapor and liquid material properties n, p_infinity
            !!       heat capacities, cv, reference energy per unit mass, q, coefficients for the
            !!       iteration procedure, A-D, and iteration variables, f and df
            !> @{
            REAL(KIND(0d0))                ::  cv1, cv2, q1, q2
            REAL(KIND(0d0))                ::        ap, bp, dp

            cv1 = fluid_pp(1)%cv; q1 = fluid_pp(1)%qv;
            cv2 = fluid_pp(2)%cv; q2 = fluid_pp(2)%qv;
            ! Calculating coefficients, Eq. C.6, Pelanti 2014
            ap = rhoalpha1*cv1 + rhoalpha2*cv2
            bp = q1*cv1*(n1-1.d0)*rhoalpha1*rhoalpha1 + q2*cv2*(n2-1.d0)*rhoalpha2*rhoalpha2 + &
                 rhoalpha1*cv1*(n1*pinf1+pinf2) + rhoalpha2*cv2*(n2*pinf2+pinf1) + &
                 rhoalpha1*rhoalpha2*(q1*cv2*(n2-1.d0)+q2*cv1*(n1-1.d0)) - &
                 E0*(cv1*(n1-1.d0)*rhoalpha1 + cv2*(n2-1.d0)*rhoalpha2)
            dp = q1*cv1*(n1-1.d0)*pinf2*rhoalpha1*rhoalpha1 + q2*cv2*(n2-1.d0)*pinf1*rhoalpha2*rhoalpha2 + &
                 pinf1*pinf2*(rhoalpha1*cv1*n1 + rhoalpha2*cv2*n2) + & 
                 rhoalpha1*rhoalpha2*(q1*cv2*(n2-1.d0)*pinf1 + q2*cv1*(n1-1.d0)*pinf2) - &
                 E0*(cv1*(n1-1.d0)*pinf2*rhoalpha1 + cv2*(n2-1.d0)*pinf1*rhoalpha2)
            ! Calculating the Tstar temperature, Eq. C.7, Pelanti 2014
            pstar = (-bp + DSQRT(bp*bp - 4.d0*ap*dp))/(2.d0*ap)
            f_alpha1_ptrelax = (cv1*(n1-1.d0)*(pstar+pinf2)*rhoalpha1)/&
                     (cv1*(n1-1.d0)*(pstar+pinf2)*rhoalpha1 + &
                      cv2*(n2-1.d0)*(pstar+pinf1)*rhoalpha2)

        END FUNCTION f_alpha1_ptrelax

        !>     The purpose of this subroutine is to determine the saturation
        !!         temperature by using a Newton-Raphson method from the provided
        !!         equilibrium pressure and EoS of the binary phase system.
        !!     @param q_cons_vf Cell-average conservative variables
        !!     @param p_star equilibrium pressure at the interface    
        SUBROUTINE s_compute_fdfTsat(fp,dfdp,pstar,Tstar)

            REAL(KIND(0d0)), INTENT(OUT)        :: fp, dfdp
            REAL(KIND(0d0)), INTENT(IN)         :: pstar, Tstar
            fp = gibbsA + gibbsB/Tstar + gibbsC*DLOG(Tstar) - & 
                 DLOG((pstar+pinf2)/(pstar+pinf1)**gibbsD)
            dfdp = -gibbsB/(Tstar*Tstar) + gibbsC/Tstar

        END SUBROUTINE s_compute_fdfTsat !-------------------------------

        !>     The purpose of this subroutine is to determine the bracket of 
        !!         the pressure by finding the pressure at which b^2=4*a*c
        !!     @param p_star equilibrium pressure at the interface    
        SUBROUTINE s_compute_Tsat_bracket(TstarA,TstarB,pressure)

            REAL(KIND(0d0)), INTENT(OUT)   :: TstarA, TstarB
            REAL(KIND(0d0)), INTENT(IN)    :: pressure
            !> @name In-subroutine variables: vapor and liquid material properties
            !> @{
            REAL(KIND(0d0))                :: fA, fB, dfdp, factor

            ! Finding lower bound, getting the bracket 
            factor = 20.d0
            TstarA = TsatLv
            TstarB = TsatLv+factor
            CALL s_compute_fdfTsat(fA,dfdp,pressure,TstarA)
            CALL s_compute_fdfTsat(fB,dfdp,pressure,TstarB)
            !PRINT *, 'fA :: ',fA,', fB :: ',fB
            DO WHILE ( fA*fB .GT. 0.d0 )
                  IF (TstarA .GT. TsatHv) THEN
                         PRINT *, 'Tsat bracketing failed to find lower bound'
                         PRINT *, 'TstarA :: ',TstarA,', pressure :: ',pressure
                         PRINT *, 'fA :: ',fA,', fB :: ',fB
                         CALL s_mpi_abort()
                  END IF
                  fA = fB
                  TstarA = TstarB
                  TstarB = TstarA+factor
                  CALL s_compute_fdfTsat(fB,dfdp,pressure,TstarB)
                  !PRINT *, 'fB :: ',fB,', TstarB :: ',TstarB,', p :: ',pressure
                  !PRINT *,'fB :: ',fB,', TstarB :: ',TstarB
                  IF( ISNAN(fB) ) THEN
                        fB = fA
                        TstarB = TstarA
                        factor = factor-10.d0
                  ELSE 
                        factor = 20.d0
                  END IF
            END DO
        END SUBROUTINE s_compute_Tsat_bracket

        !>     The purpose of this subroutine is to determine the saturation
        !!         temperature by using a Newton-Raphson method from the provided
        !!         equilibrium pressure and EoS of the binary phase system.
        !!     @param p_star equilibrium pressure at the interface    
        FUNCTION f_Tsat(pressure)

            REAL(KIND(0d0)), INTENT(IN) :: pressure
            REAL(KIND(0d0))             :: f_Tsat, Tstar
            !> @name In-subroutine variables: vapor and liquid material properties n, p_infinity
            !!       heat capacities, cv, reference energy per unit mass, q, coefficients for the
            !!       iteration procedure, A-D, and iteration variables, f and df
            !> @{
            REAL(KIND(0d0)), INTENT(IN)    :: pressure
            REAL(KIND(0d0))                :: f_Tsat, Tstar
            REAL(KIND(0d0))                ::  delta, delta_old, fp, dfdp
            REAL(KIND(0d0))                ::  fL, fH, TstarL, TstarH, TsatA, TsatB
            INTEGER :: iter                !< Generic loop iterators
            CALL s_compute_Tsat_bracket(TsatA,TsatB,pressure)
            ! Computing f at lower and higher end of the bracket
            CALL s_compute_fdfTsat(fL,dfdp,pressure,TsatA)
            CALL s_compute_fdfTsat(fH,dfdp,pressure,TsatB)
            ! Establishing the direction of the descent to find zero
            IF(fL < 0.d0) THEN
                TstarL  = TsatA; TstarH  = TsatB;
            ELSE
                TstarL  = TsatA; TstarH  = TsatB;
            END IF
            Tstar = 0.5d0*(TstarL+TstarH)
            delta_old = DABS(TstarH-TstarL)
            delta = delta_old
            CALL s_compute_fdfTsat(fp,dfdp,pressure,Tstar)
            ! Combining bisection and newton-raphson methods
            DO iter = 0, newton_iter
                IF ((((Tstar-TstarH)*dfdp-fp)*((Tstar-TstarL)*dfdp-fp) > 0.d0) & ! Bisect if Newton out of range,
                        .OR. (DABS(2.0*fp) > DABS(delta_old*dfdp))) THEN         ! or not decreasing fast enough.
                   delta_old = delta
                   delta = 0.5d0*(TstarH-TstarL)
                   Tstar = TstarL + delta
                   IF (delta .EQ. 0.d0) EXIT
                ELSE                    ! Newton step acceptable, take it
                   delta_old = delta
                   delta = fp/dfdp
                   Tstar = Tstar - delta
                   IF (delta .EQ. 0.d0) EXIT
                END IF
                IF (DABS(delta/Tstar) < pTsatnewton_eps) EXIT
                CALL s_compute_fdfTsat(fp,dfdp,pressure,Tstar)           
                IF (fp < 0.d0) THEN     !Maintain the bracket on the root
                   TstarL = Tstar
                ELSE
                   TstarH = Tstar
                END IF
                IF (iter .EQ. newton_iter) THEN
                    PRINT *, 'Tsat : ',Tstar,', iter : ',iter
                    PRINT *, 'Tsat did not converge, stopping code'
                    CALL s_mpi_abort()
                END IF                  
            END DO
            f_Tsat = Tstar
        END FUNCTION f_Tsat !-------------------------------

        !>     The purpose of this subroutine is to determine the saturation
        !!         temperature by using a Newton-Raphson method from the provided
        !!         equilibrium pressure and EoS of the binary phase system.
        !!     @param p_star equilibrium pressure at the interface    
        SUBROUTINE s_compute_ptg_fdf(fp,dfdp,pstar,Tstar,rho0,E0)

            REAL(KIND(0d0)), INTENT(IN)    :: pstar, rho0, E0
            REAL(KIND(0d0)), INTENT(OUT)   :: fp, dfdp, Tstar
            !> @name In-subroutine variables: vapor and liquid material properties n, p_infinity
            !!       heat capacities, cv, reference energy per unit mass, q, coefficients for the
            !!       iteration procedure, A-D, and iteration variables, f and df
            !> @{
            REAL(KIND(0d0))                ::  cv1, cv2, q1, q2
            REAL(KIND(0d0))                ::        ap, bp, dp
            REAL(KIND(0d0))                ::  dadp, dbdp, dddp
            REAL(KIND(0d0))                ::              dTdp

            cv1 = fluid_pp(1)%cv; q1 = fluid_pp(1)%qv;
            cv2 = fluid_pp(2)%cv; q2 = fluid_pp(2)%qv;
            ! Calculating coefficients, Eq. C.6, Pelanti 2014
            ap = rho0*cv1*cv2*((n2-1.d0)*(pstar+n1*pinf1)-(n1-1.d0)*(pstar+n2*pinf2))
            bp = E0*((n1-1.d0)*cv1*(pstar+pinf2) - (n2-1.d0)*cv2*(pstar+pinf1)) + &
                 rho0*((n2-1.d0)*cv2*q1*(pstar+pinf1) - (n1-1.d0)*cv1*q2*(pstar+pinf2)) + &
                 cv2*(pstar+pinf1)*(pstar+n2*pinf2) - cv1*(pstar+pinf2)*(pstar+n1*pinf1)
            dp = (q2-q1)*(pstar+pinf1)*(pstar+pinf2)
            ! Calculating the Tstar temperature, Eq. C.7, Pelanti 2014
            Tstar = (-bp + DSQRT(bp*bp - 4.d0*ap*dp))/(2.d0*ap)
            ! PRINT *, 'ap: ', ap, 'bp: ', bp, 'dp:', dp
            ! PRINT *, 'sqrt(D)', DSQRT(bp*bp - 4.d0*ap*dp)
            ! PRINT *, 'Tstar: ', Tstar
            ! PRINT *, 'Tstar T : ', (-bp - DSQRT(bp*bp - 4.d0*ap*dp))/(2.d0*ap)
            ! Calculating the derivatives wrt pressure of the coefficients
            dadp = rho0*cv1*cv2*((n2-1.d0)-(n1-1.d0))
            dbdp = E0*((n1-1.d0)*cv1 - (n2-1.d0)*cv2) + &
                   rho0*((n2-1.d0)*cv2*q1 - (n1-1.d0)*cv1*q2) + &
                   cv2*((pstar+pinf1)+(pstar+n2*pinf2)) - & 
                   cv1*((pstar+pinf2)+(pstar+n1*pinf1))
            dddp = (q2-q1)*((pstar+pinf1)+(pstar+pinf2))
            ! Derivative of the temperature wrt to pressure, needed for dfdp
            dTdp = (-dbdp + (0.5d0/DSQRT(bp*bp-4.d0*ap*dp))*(2.d0*bp*dbdp-&
                   4.d0*(ap*dddp+dp*dadp)))/(2.d0*ap) - (dadp/ap)*Tstar
            fp   = gibbsA + gibbsB/Tstar + gibbsC*DLOG(Tstar) + & 
                   gibbsD*DLOG(pstar+pinf1) - DLOG(pstar+pinf2)
            dfdp = -gibbsB/(Tstar*Tstar)*dTdp + gibbsC/Tstar*dTdp + & 
                   gibbsD/(pstar+pinf1) - 1.d0/(pstar+pinf2)

            ! PRINT *, Tstar
            ! PRINT *, dadp, dbdp, dddp, dTdp, fp, dfdp

        END SUBROUTINE s_compute_ptg_fdf !------------------------

        !>     The purpose of this subroutine is to determine the bracket of 
        !!         the pressure by finding the pressure at which b^2=4*a*c
        !!     @param p_star equilibrium pressure at the interface    
        SUBROUTINE s_compute_ptg_bracket(pstarA,pstarB,pstar,rho0,E0)

            REAL(KIND(0d0)), INTENT(OUT)   :: pstarA, pstarB
            REAL(KIND(0d0)), INTENT(IN)    :: pstar, rho0, E0
            !> @name In-subroutine variables: vapor and liquid material properties n, p_infinity
            !!       heat capacities, cv, reference energy per unit mass, q, coefficients for the
            !!       iteration procedure, A-D, and iteration variables, f and df
            !> @{
            REAL(KIND(0d0))                :: fA, fB, dfdp, Tstar
            REAL(KIND(0d0))                :: pS, factor

            pstarA = 10000.d0
            pstarB = pstar
            ! PRINT *, 's_compute_ptg_fdf A'
            CALL s_compute_ptg_fdf(fA,dfdp,pstarA,Tstar,rho0,E0)
            ! PRINT *, 's_compute_ptg_fdf B'
            CALL s_compute_ptg_fdf(fB,dfdp,pstarB,Tstar,rho0,E0)
            factor = 1.05d0
            DO WHILE ( fA*fB .GT. 0.d0 )
                ! PRINT *, 'fA', fA
                ! PRINT *, 'fB', fb
                  IF (pstarA .GT. 1.d13) THEN
                         PRINT *, 'ptg bracketing failed to find lower bound'
                         PRINT *, 'pstarA :: ',pstarA
                         CALL s_mpi_abort()
                  END IF
                  fA = fB
                  pstarA = pstarB
                  pstarB = pstarA*factor
                  CALL s_compute_ptg_fdf(fB,dfdp,pstarB,Tstar,rho0,E0)
                  IF( ISNAN(fB) ) THEN
                        fB = fA
                        pstarB = pstarA
                        factor = factor*0.95d0
                  ELSE 
                        factor = 1.05d0
                  END IF
            END DO

        END SUBROUTINE s_compute_ptg_bracket

        !>     The purpose of this subroutine is to determine the saturation
        !!         temperature by using a Newton-Raphson method from the provided
        !!         equilibrium pressure and EoS of the binary phase system.
        !!     @param p_star equilibrium pressure at the interface    
        SUBROUTINE s_compute_ptg_pTrelax(pstar,Tstar,rho0,E0)

            REAL(KIND(0d0)), INTENT(INOUT) :: pstar
            REAL(KIND(0d0)), INTENT(OUT)   :: Tstar
            REAL(KIND(0d0)), INTENT(IN)    :: rho0, E0
            !> @name In-subroutine variables: vapor and liquid material properties n, p_infinity
            !!       heat capacities, cv, reference energy per unit mass, q, coefficients for the
            !!       iteration procedure, A-D, and iteration variables, f and df
            !> @{
            REAL(KIND(0d0))                ::  pstarA, pstarB, dTstar
            REAL(KIND(0d0))                ::  delta, delta_old, fp, dfdp
            REAL(KIND(0d0))                ::  fL, fH, pstarL, pstarH
            INTEGER :: iter                !< Generic loop iterators
            ! Computing the bracket of the root solution
            CALL s_compute_ptg_bracket(pstarA,pstarB,pstar,rho0,E0)
            ! PRINT *, 'pstarA', pstarA, 'pstarB', pstarB 
            ! Computing f at lower and higher end of the bracket
            CALL s_compute_ptg_fdf(fL,dfdp,pstarA,dTstar,rho0,E0)
            ! PRINT *, 'fL', fL, 'dfdpL', dfdp, 'TstarL', dTstar
            ! CALL s_compute_ptg_fdf(fH,dfdp,pstarB,dTstar,rho0,E0)
            ! PRINT *, 'fH', fH, 'dfdpH', dfdp, 'TstarH', dTstar
            ! Establishing the direction of the descent to find zero
            IF(fL < 0.d0) THEN
                pstarL  = pstarA; pstarH  = pstarB;
            ELSE
                pstarL  = pstarB; pstarH  = pstarA;
            END IF
            pstar = 0.5d0*(pstarA+pstarB)
            delta_old = DABS(pstarB-pstarA)
            delta = delta_old
            CALL s_compute_ptg_fdf(fp,dfdp,pstar,dTstar,rho0,E0)
            ! Combining bisection and newton-raphson methods
            DO iter = 0, newton_iter
                IF ((((pstar-pstarH)*dfdp-fp)*((pstar-pstarL)*dfdp-fp) > 0.d0) & ! Bisect if Newton out of range,
                        .OR. (DABS(2.0*fp) > DABS(delta_old*dfdp))) THEN         ! or not decreasing fast enough.
                   delta_old = delta
                   delta = 0.5d0*(pstarH-pstarL)
                   pstar = pstarL + delta
                   IF (delta .EQ. 0.d0) EXIT                    ! Change in root is negligible
                ELSE                                              ! Newton step acceptable, take it
                   delta_old = delta
                   delta = fp/dfdp
                   pstar = pstar - delta
                   IF (delta .EQ. 0.d0) EXIT
                END IF
                IF (DABS(delta/pstar) < ptgnewton_eps) EXIT           ! Stopping criteria
                ! Updating to next iteration
                CALL s_compute_ptg_fdf(fp,dfdp,pstar,Tstar,rho0,E0)
                IF (fp < 0.d0) THEN !Maintain the bracket on the root
                   pstarL = pstar
                ELSE
                   pstarH = pstar
                END IF
            END DO
        END SUBROUTINE s_compute_ptg_pTrelax !-------------------------------

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!  K-PHASE P RELAXATION FUNCTIONS !!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !>     The purpose of this subroutine is to determine the saturation
        !!         temperature by using a Newton-Raphson method from the provided
        !!         equilibrium pressure and EoS of the binary phase system.
        !!     @param p_star equilibrium pressure at the interface    
        SUBROUTINE s_compute_pk_fdf(fp,dfdp,pstar,rho_K_s,gamma_min,pres_inf,pres_K_init,q_cons_vf,j,k,l)
            !> @name In-subroutine variables: vapor and liquid material properties n, p_infinity
            !!       heat capacities, cv, reference energy per unit mass, q, coefficients for the
            !!       iteration procedure, A-D, and iteration variables, f and df
            !> @{
            ! Cell-average conservative variables
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN)  :: q_cons_vf
            REAL(KIND(0d0)), INTENT(OUT)                         :: fp, dfdp
            REAL(KIND(0d0)), INTENT(IN)                          :: pstar
            REAL(KIND(0d0)), DIMENSION(num_fluids), INTENT(OUT)  :: rho_K_s
            REAL(KIND(0d0)), DIMENSION(num_fluids), INTENT(IN)   :: gamma_min, pres_inf, pres_K_init
            REAL(KIND(0d0))                                      :: numerator, denominator
            REAL(KIND(0d0))                                      :: drhodp
            INTEGER, INTENT(IN)                                  :: j, k, l
            INTEGER                                              :: i
            fp = -1.d0; dfdp = 0.d0;
            DO i = 1, num_fluids
                  numerator   = gamma_min(i)*(pstar+pres_inf(i))
                  denominator = numerator + pres_K_init(i)-pstar
                  rho_K_s(i)  = q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)/&
                     MAX(q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l),sgm_eps)*numerator/denominator
                  drhodp      = q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l) / & 
                     MAX(q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l),sgm_eps) * & 
                     gamma_min(i)*(pres_K_init(i)+pres_inf(i)) / (denominator*denominator)
                  fp          = fp  + q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l) / rho_K_s(i)
                  dfdp        = dfdp - q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l) * & 
                     drhodp / (rho_K_s(i)*rho_K_s(i))
            END DO

        END SUBROUTINE s_compute_pk_fdf !------------------------

        !>     The purpose of this subroutine is to determine the bracket of 
        !!         the pressure by finding the pressure at which b^2=4*a*c
        !!     @param p_star equilibrium pressure at the interface    
        SUBROUTINE s_compute_pk_bracket(pstarA,pstarB,rho_K_s,gamma_min,pres_inf,pres_K_init,q_cons_vf,j,k,l)
            !> @name In-subroutine variables: vapor and liquid material properties n, p_infinity
            !!       heat capacities, cv, reference energy per unit mass, q, coefficients for the
            !!       iteration procedure, A-D, and iteration variables, f and df
            !> @{
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN)  :: q_cons_vf
            REAL(KIND(0d0)), INTENT(OUT)   :: pstarA, pstarB
            REAL(KIND(0d0))                :: fA, fB, dfdp
            REAL(KIND(0d0))                :: factor
            REAL(KIND(0d0)), DIMENSION(num_fluids), INTENT(IN)   :: gamma_min, pres_inf, pres_K_init
            REAL(KIND(0d0)), DIMENSION(num_fluids), INTENT(OUT)  :: rho_K_s
            INTEGER, INTENT(IN)            :: j,k,l
            !pstarA = 100.d0
            !pstarB = 5.d2
            pstarA = 1.d-15
            pstarB = 1.d1
            CALL s_compute_pk_fdf(fA,dfdp,pstarA,rho_K_s,gamma_min,pres_inf,pres_K_init,q_cons_vf,j,k,l)
            CALL s_compute_pk_fdf(fB,dfdp,pstarB,rho_K_s,gamma_min,pres_inf,pres_K_init,q_cons_vf,j,k,l)
            factor = 10.d0
            DO WHILE ( fA*fB .GT. 0.d0 )
                  IF (pstarA .GT. 1.d13) THEN
                         PRINT *, 'P-K bracketing failed to find lower bound'
                         PRINT *, 'pstarA :: ',pstarA,', pstarB :: ',pstarB
                         PRINT *, 'alpha1 :: ',q_cons_vf(adv_idx%beg)%sf(j,k,l)
                         PRINT *, 'alpha2 :: ',q_cons_vf(adv_idx%beg)%sf(j,k,l)
                         CALL s_mpi_abort()
                  END IF
                  fA = fB
                  pstarA = pstarB
                  pstarB = pstarA*factor
                  CALL s_compute_pk_fdf(fB,dfdp,pstarB,rho_K_s,gamma_min,pres_inf,pres_K_init,q_cons_vf,j,k,l)
                  IF( ISNAN(fB) ) THEN
                        fB = fA
                        pstarB = pstarA
                        factor = factor*0.5d0
                  ELSE 
                        factor = 10.d0
                  END IF
            END DO
        END SUBROUTINE s_compute_pk_bracket

        !>     The purpose of this subroutine is to determine the saturation
        !!         temperature by using a Newton-Raphson method from the provided
        !!         equilibrium pressure and EoS of the binary phase system.
        !!     @param p_star equilibrium pressure at the interface    
        SUBROUTINE s_compute_p_relax_k(rho_K_s,gamma_min,pres_inf,pres_K_init,q_cons_vf,j,k,l)
            !> @name In-subroutine variables: vapor and liquid material properties n, p_infinity
            !!       heat capacities, cv, reference energy per unit mass, q, coefficients for the
            !!       iteration procedure, A-D, and iteration variables, f and df
            !> @{
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN)  :: q_cons_vf
            REAL(KIND(0d0))                                      :: pstar, pstarA, pstarB
            REAL(KIND(0d0))                                      :: delta, delta_old, fp, dfdp
            REAL(KIND(0d0))                                      :: fL, fH, pstarL, pstarH
            REAL(KIND(0d0)), DIMENSION(num_fluids), INTENT(IN)   :: gamma_min, pres_inf, pres_K_init
            REAL(KIND(0d0)), DIMENSION(num_fluids), INTENT(OUT)  :: rho_K_s 
            INTEGER, INTENT(IN)                                  :: j,k,l
            INTEGER :: iter                !< Generic loop iterators
            ! Computing the bracket of the root solution
            CALL s_compute_pk_bracket(pstarA,pstarB,rho_K_s,gamma_min,pres_inf,pres_K_init,q_cons_vf,j,k,l)
            ! Computing f at lower and higher end of the bracket
            CALL s_compute_pk_fdf(fL,dfdp,pstarA,rho_K_s,gamma_min,pres_inf,pres_K_init,q_cons_vf,j,k,l)
            CALL s_compute_pk_fdf(fH,dfdp,pstarB,rho_K_s,gamma_min,pres_inf,pres_K_init,q_cons_vf,j,k,l)
            ! Establishing the direction of the descent to find zero
            IF(fL < 0.d0) THEN
                pstarL  = pstarA; pstarH  = pstarB;
            ELSE
                pstarL  = pstarB; pstarH  = pstarA;
            END IF
            pstar = 0.5d0*(pstarA+pstarB)
            delta_old = DABS(pstarB-pstarA)
            delta = delta_old
            CALL s_compute_pk_fdf(fp,dfdp,pstar,rho_K_s,gamma_min,pres_inf,pres_K_init,q_cons_vf,j,k,l)
            ! Combining bisection and newton-raphson methods
            DO iter = 0, newton_iter
                IF ((((pstar-pstarH)*dfdp-fp)*((pstar-pstarL)*dfdp-fp) > 0.d0) & ! Bisect if Newton out of range,
                        .OR. (DABS(2.0*fp) > DABS(delta_old*dfdp))) THEN         ! or not decreasing fast enough.
                   delta_old = delta
                   delta = 0.5d0*(pstarH-pstarL)
                   pstar = pstarL + delta
                   IF (delta .EQ. 0.d0) EXIT                    ! Change in root is negligible
                ELSE                                              ! Newton step acceptable, take it
                   delta_old = delta
                   delta = fp/dfdp
                   pstar = pstar - delta
                   IF (delta .EQ. 0.d0) EXIT
                END IF
                IF (DABS(delta/pstar) < pknewton_eps) EXIT           ! Stopping criteria
                ! Updating to next iteration
                CALL s_compute_pk_fdf(fp,dfdp,pstar,rho_K_s,gamma_min,pres_inf,pres_K_init,q_cons_vf,j,k,l)
                IF (fp < 0.d0) THEN !Maintain the bracket on the root
                   pstarL = pstar
                ELSE
                   pstarH = pstar
                END IF
            END DO
        END SUBROUTINE s_compute_p_relax_k !-------------------------------

        ! This SUBROUTINE IS USED TO CALCULATE THE 2X2 JACOBIAN AND ITS INVERSE FOR THE PROPOSED
        ! pTg-Equilibrium procedure
        SUBROUTINE s_compute_jacobian_matrix( B, C, D, rhoe, g_min, InvJac, j, Jac, k, l, lp, mCPD, &
                                                mCVGP, mCVGP2, mQD, p_inf, pS, q_cons_vf, TJac, vp )
        
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN)        :: q_cons_vf
            REAL( KIND( 0.0D0 ) ), INTENT(IN)                          :: B, C, D
            REAL( KIND( 0.0D0 ) ), INTENT(IN)                          :: pS, rhoe, mCPD, mCVGP, mCVGP2, mQD
            REAL( KIND( 0.0D0 ) ), DIMENSION(num_fluids), INTENT(IN)   :: g_min, p_inf
            REAL( KIND( 0.0D0 ) ), DIMENSION(2,2), INTENT(OUT)         :: Jac, InvJac, TJac
            INTEGER, INTENT(IN)                                        :: j, k, l, lp, vp
            REAL( KIND( 0.0D0 ) )                                      :: ml, mT
            INTEGER                                                    :: i
            
            ! mass of the reactant liquid
            ml = q_cons_vf( lp + cont_idx%beg - 1 )%sf( j, k, l )
            
            ! mass of the two participating fluids
            mT = q_cons_vf( lp + cont_idx%beg - 1 )%sf( j, k, l ) &
            + q_cons_vf( vp + cont_idx%beg - 1 )%sf( j, k, l )         

            ! F = (F1,F2) is the fuction whose roots we are looking for
            ! x = (m1, p) are the independent variables. m1 = mass of the first participant fluid, p = pressure
            ! F1 = 0 is the gibbs free energy quality
            ! F2 = 0 is the enforcement of the thermodynamic (total - kinectic) energy

            ! dF1dm
            Jac(1,1) = ( fluid_pp( lp )%cv * ( g_min( lp ) - 1 ) / ( pS + p_inf( lp ) )            &
            -            fluid_pp( vp )%cv * ( g_min( vp ) - 1 ) / ( pS + p_inf( vp ) ) )          &
            * ( B - C / ( ml * ( fluid_pp( lp )%cv * ( g_min( lp ) - 1 ) / ( pS + p_inf( lp ) )    & 
                               - fluid_pp( vp )%cv * ( g_min( vp ) - 1 ) / ( pS + p_inf( vp ) ) )  & 
            + mT * fluid_pp( vp )%cv * ( g_min( vp ) - 1 ) / ( pS + p_inf( vp ) ) + mCVGP ) )

            ! dF1dp
            Jac(1,2) = ( ml * ( fluid_pp( lp )%cv * ( g_min( lp ) - 1 ) / ( pS + p_inf( lp ) ) ** 2            &
                            - fluid_pp( vp )%cv * ( g_min( vp ) - 1 ) / ( pS + p_inf( vp ) )   ** 2 )          &
            +            mT * fluid_pp( vp )%cv * ( g_min( vp ) - 1 ) / ( pS + p_inf( vp ) )   ** 2 + mCVGP2 ) &
            * ( - B + C / ( ml * ( fluid_pp( lp )%cv * ( g_min( lp ) - 1 ) / ( pS + p_inf( lp ) )              &
                                 - fluid_pp( vp )%cv * ( g_min( vp ) - 1 ) / ( pS + p_inf( vp ) ) )            &
            + mT * fluid_pp( vp )%cv * ( g_min( vp ) - 1 ) / ( pS + p_inf( vp ) ) + mCVGP ) )                  &
            + D / ( pS + p_inf( lp ) ) - 1 / ( pS + p_inf( vp ) )

            ! dF2dm
            Jac(2,1) = fluid_pp( vp )%qv - fluid_pp( lp )%qv                                                &
            + ( fluid_pp( vp )%cv * g_min( vp ) - fluid_pp( lp )%cv * g_min( lp ) )                         &
            / ( ml * ( fluid_pp( lp )%cv * ( g_min( lp ) - 1 ) / ( pS + p_inf( lp ) )                       &
                     - fluid_pp( vp )%cv * ( g_min( vp ) - 1 ) / ( pS + p_inf( vp ) ) )                     & 
              + mT *   fluid_pp( vp )%cv * ( g_min( vp ) - 1 ) / ( pS + p_inf( vp ) ) + mCVGP )             &
            - ( ml * ( fluid_pp( vp )%cv * g_min( vp ) - fluid_pp( lp )%cv * g_min( lp ) )                  &
            -   mT * fluid_pp( vp )%cv * g_min( vp ) - mCPD )                                               &
            * ( fluid_pp( lp )%cv * ( g_min( lp ) - 1 ) / ( pS + p_inf( lp ) )                              &
            -   fluid_pp( vp )%cv * ( g_min( vp ) - 1 ) / ( pS + p_inf( vp ) ) )                            &
            / ( ( ml * ( fluid_pp( lp )%cv * ( g_min( lp ) - 1 ) / ( pS + p_inf( lp ) )                     &
            -            fluid_pp( vp )%cv * ( g_min( vp ) - 1 ) / ( pS + p_inf( vp ) ) )                   &
            +     mT *   fluid_pp( vp )%cv * ( g_min( vp ) - 1 ) / ( pS + p_inf( vp ) ) + mCVGP ) ** 2 )

            ! dF2dp
            ! Jac(2,2) = 4
            Jac(2,2) = 1 + ( ml * ( fluid_pp( vp )%cv * g_min( vp ) - fluid_pp( lp )%cv * g_min( lp ) )     &
            -                mT *   fluid_pp( vp )%cv * g_min( vp ) - mCPD )                                & 
            * ( ml * ( fluid_pp( lp )%cv * ( g_min( lp ) - 1 ) / ( pS + p_inf( lp ) )  ** 2                 & 
                     - fluid_pp( vp )%cv * ( g_min( vp ) - 1 ) / ( pS + p_inf( vp ) )  ** 2 )               &
            +   mT *   fluid_pp( vp )%cv * ( g_min( vp ) - 1 ) / ( pS + p_inf( vp ) )  ** 2 + mCVGP2 )      &
            / ( ml * ( fluid_pp( lp )%cv * ( g_min( lp ) - 1 ) / ( pS + p_inf( lp ) )                       &   
                     - fluid_pp( vp )%cv * ( g_min( vp ) - 1 ) / ( pS + p_inf( vp ) ) )                     &
            +   mT *   fluid_pp( vp )%cv * ( g_min( vp ) - 1 ) / ( pS + p_inf( vp ) ) + mCVGP ) ** 2

            ! elements of J^{-1}
            InvJac(1,1) =  Jac(2,2)
            InvJac(1,2) = - 1.0D0 * Jac(1,2)
            InvJac(2,1) = - 1.0D0 * Jac(2,1)
            InvJac(2,2) = Jac(1,1)

            ! elements of J^{T}
            TJac(1,1) = Jac(1,1)
            TJac(1,2) = Jac(2,1)
            TJac(2,1) = Jac(1,2)
            TJac(2,2) = Jac(2,2)

            ! dividing by det(J)
            InvJac = InvJac / ( Jac(1,1) * Jac(2,2) - Jac(1,2) * Jac(2,1) )             

        END SUBROUTINE s_compute_jacobian_matrix

        SUBROUTINE s_compute_pTg_residue(A, B, C, D, g_min, j, k, l, lp, mCPD, mCVGP, mQD &
                                            , p_inf, q_cons_vf, pS, rhoe, R2D, vp)
            
            TYPE(scalar_field), DIMENSION(sys_size), INTENT(IN)         :: q_cons_vf
            REAL( KIND( 0.0D0 ) ), INTENT(IN)                           :: A, B, C, D
            REAL( KIND( 0.0D0 ) ), INTENT(IN)                           :: pS, rhoe, mCPD, mCVGP, mQD
            REAL( KIND( 0.0D0 ) ), DIMENSION(num_fluids), INTENT(IN)    :: g_min, p_inf
            REAL( KIND( 0.0D0 ) ), DIMENSION(2), INTENT(OUT)            :: R2D
            INTEGER, INTENT(IN)                                         :: j, k, l, lp, vp
            REAL( KIND( 0.0D0 ) )                                       :: ml, mT
            INTEGER                                                     :: i
            
            ! mass of the reactant liquid
            ml = q_cons_vf( lp + cont_idx%beg - 1 )%sf( j, k, l )
            
            ! mass of the two participating fluids
            mT = q_cons_vf( lp + cont_idx%beg - 1 )%sf( j, k, l ) &
            + q_cons_vf( vp + cont_idx%beg - 1 )%sf( j, k, l )            

            ! Gibbs Free Energy Equality condition
            R2D( 1 ) = A &
                + B * ( ml * ( fluid_pp( lp )%cv * ( g_min( lp ) - 1 ) / ( pS + p_inf( lp ) ) &
                             - fluid_pp( vp )%cv * ( g_min( vp ) - 1 ) / ( pS + p_inf( vp ) ) ) &
                        + mT * fluid_pp( vp )%cv * ( g_min( vp ) - 1 ) / ( pS + p_inf( vp ) ) & 
                        + mCVGP ) &
            - C * DLOG( ml * ( fluid_pp( lp )%cv * ( g_min( lp ) - 1 ) / ( pS + p_inf( lp ) ) &
                             - fluid_pp( vp )%cv * ( g_min( vp ) - 1 ) / ( pS + p_inf( vp ) ) ) &
                        + mT * fluid_pp( vp )%cv * ( g_min( vp ) - 1 ) / ( pS + p_inf( vp ) ) &
                        + mCVGP ) &
            + D * DLOG( pS + p_inf( lp ) ) - DLOG( pS + p_inf( vp ) )
            
            ! Constant Energy Process condition
            R2D( 2 ) = rhoe + pS &
            + ml * ( fluid_pp( vp )%qv - fluid_pp( lp )%qv ) - mT * fluid_pp( vp )%qv - mQD &
            + ( ml * ( g_min( vp ) * fluid_pp( vp )%cv - g_min( lp ) * fluid_pp( lp )%cv ) &
            - mT * g_min( vp ) * fluid_pp( vp )%cv - mCPD ) &
            / ( ml * ( fluid_pp( lp )%cv * ( g_min( lp ) - 1 ) / ( pS + p_inf( lp ) ) &
                     - fluid_pp( vp )%cv * ( g_min( vp ) - 1 ) / ( pS + p_inf( vp ) ) ) &
            + mT * fluid_pp( vp )%cv * ( g_min( vp ) - 1 ) / ( pS + p_inf( vp ) ) + mCVGP )

            ! IF ( j .EQ. 1503 ) THEN
                
            !     PRINT *, A, B, C, D

            !     PRINT *, ml, mT

            !     PRINT *, pS

            !     PRINT *, rhoe

            !     PRINT *, R2D

            !     PAUSE

            ! END IF

        END SUBROUTINE s_compute_pTg_residue

        SUBROUTINE s_finalize_relaxation_solver_module()

           s_relaxation_solver => NULL()

        END SUBROUTINE

END MODULE m_phase_change