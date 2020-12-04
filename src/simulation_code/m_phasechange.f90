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
!! @brief Contains module m_qbmm
!! @author M. Rodriguez
!! @version 1.0
!! @date DEC 3, 2020

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
                       s_finite_ptg_relaxation, &
                       s_infinite_p_relaxation, &
                       s_infinite_p_relaxation_k, &
                       s_infinite_ptg_relaxation

    !> @name Parameters for the phase change part of the code
    !> @{
    REAL(KIND(0d0)), PARAMETER :: pnewtonk_eps      = 1.d-10    !< p_relaxk \alpha threshold,           set to 1E-15
    INTEGER,         PARAMETER :: pnewtonk_iter     = 50        !< p_relaxk \alpha iter,                set to 25
    REAL(KIND(0d0)), PARAMETER :: pTsatnewton_eps   = 1.d-10    !< Saturation temperature tol,          set to 1E-12
    INTEGER,         PARAMETER :: pTsatnewton_iter  = 50        !< Saturation temperature iteration,    set to 25
    REAL(KIND(0d0)), PARAMETER :: TsatHv            = 2000.d0    !< Saturation temperature threshold,    set to 900
    REAL(KIND(0d0)), PARAMETER :: TsatLv            = 100.d0    !< Saturation temperature threshold,    set to 250
    REAL(KIND(0d0)), PARAMETER :: palpha_epsH       = 1.d-6     !< p_relax high \alpha tolerance,       set to 1.d-6
    REAL(KIND(0d0)), PARAMETER :: palpha_epsL       = 1.d-6     !< p_relax low \alpha tolerance,        set to 1.d-6
    REAL(KIND(0d0)), PARAMETER :: ptgalpha_epsH     = 1.d-6     !< Saturation p-T-mu alpha tolerance,   set to 1.d-6
    REAL(KIND(0d0)), PARAMETER :: ptgalpha_epsL     = 1.d-6     !< Saturation p-T-mu alpha tolerance,   set to 1.d-6
    REAL(KIND(0d0)), PARAMETER :: ptgnewton_eps     = 1.d-10    !< Saturation p-T-mu tolerance,         set to 1.d-10
    INTEGER,         PARAMETER :: ptgnewton_iter    = 20        !< Saturation p-T-mu iteration,         set to 50
    !> @}

    !> @name Gibbs free energy phase change parameters
    !> @{
    REAL(KIND(0d0)) :: n1, n2, pinf1, pinf2
    REAL(KIND(0d0)) :: gibbsA, gibbsB, gibbsC, gibbsD
    !> @}

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

        END SUBROUTINE s_initialize_phasechange_module !-------------------------------

        ! ==================================================================                     
        ! Mixture-total-energy correction ==================================
        ! ==================================================================                     
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
            REAL(KIND(0d0)), DIMENSION(2)                     ::          Re
            REAL(KIND(0d0)), DIMENSION(num_fluids,num_fluids) ::          We
            !> @}
            INTEGER :: i           !< Generic loop iterators

            CALL s_convert_to_mixture_variables( q_cons_vf, rho, &
                                                 gamma, pi_inf,  &
                                                 Re, We, j, k, l )
            dyn_pres = 0.d0
            DO i = mom_idx%beg, mom_idx%end
                 dyn_pres = dyn_pres + 5d-1*q_cons_vf(i)%sf(j,k,l) * & 
                   q_cons_vf(i)%sf(j,k,l) / MAX(rho,sgm_eps)
            END DO
            E_We = 0.d0
            pres_relax = (q_cons_vf(E_idx)%sf(j,k,l) - dyn_pres - pi_inf - E_We)/gamma
            DO i = 1, num_fluids
                  q_cons_vf(i+internalEnergies_idx%beg-1)%sf(j,k,l) = & 
                  q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) * & 
                  (fluid_pp(i)%gamma*pres_relax + fluid_pp(i)%pi_inf) +  & 
                  q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)*fluid_pp(i)%qv
            END DO
            ! ==================================================================
        END SUBROUTINE s_mixture_total_energy_correction

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
               IF ((q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l) .LT. 0d0) .OR. &
                   (q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) .LT. 0d0)) THEN
                    q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l) = sgm_eps
                    q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l)  = sgm_eps
                    q_cons_vf(i+internalEnergies_idx%beg-1)%sf(j,k,l)  = 0d0
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
            REAL(KIND(0d0))                :: fA, fB, dfdpA, dfdpB, factor

            !> @name In-subroutine variables: vapor and liquid material properties n, p_infinity
            !!       heat capacities, cv, reference energy per unit mass, q, coefficients for the
            !!       iteration procedure, A-D, and iteration variables, f and df
            !> @{

            ! Finding lower bound, getting the bracket within one order of magnitude
            TstarA = 1.d0
            TstarB = TstarA
            CALL s_compute_fdfTsat(fA,dfdpA,pressure,TstarA)
            fB = fA
            factor = 10.d0
            DO WHILE ( fA*fB .GT. 0.d0 )
                  PRINT *, 'fA: ',fA,', fB: ',fB,', TstarA: ',TstarA,', TstarB :',TstarB
                  IF (TstarA .GT. 1.d13) THEN
                         PRINT *, 'Tsat bracketing failed to find lower bound'
                         PRINT *, 'TstarA :: ',TstarA
                         CALL s_mpi_abort()
                  END IF
                  fA = fB
                  TstarA = TstarB
                  TstarB = TstarA + factor
                  dfdpA = dfdpB
                  CALL s_compute_fdfTsat(fB,dfdpB,pressure,TstarB)
                  IF( ISNAN(fB) ) THEN
                        fB = fA
                        factor = factor*0.5d0
                  ELSE 
                        factor = 10.d0
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
            REAL(KIND(0d0))                ::  delta, delta_old, fp, dfdp
            REAL(KIND(0d0))                ::  fL, fH, TstarL, TstarH, TsatA, TsatB
            INTEGER :: iter      !< Generic loop iterators
            PRINT *, 'dfdpA :',dfdp,', pressure: ',pressure
            CALL s_compute_Tsat_bracket(TsatA,TsatB,pressure)
            ! Computing f at lower and higher end of the bracket
            CALL s_compute_fdfTsat(fL,dfdp,pressure,TsatA)
            PRINT *, 'dfdpA :',dfdp,', pressure: ',pressure
            CALL s_compute_fdfTsat(fH,dfdp,pressure,TsatB)
            PRINT *, 'dfdpB :',dfdp

            ! Establishing the direction of the descent to find zero
            IF(fL < 0.d0) THEN
                TstarL  = TsatA; TstarH  = TsatB;
            ELSE
                TstarL  = TsatA; TstarH  = TsatB;
            END IF
            PRINT *, 'fL : ',fL,', fH : ',fH,', TL : ',TstarL,', TH : ',TstarH
            Tstar = 0.5d0*(TstarL+TstarH)
            delta_old = DABS(TstarH-TstarL)
            delta = delta_old
            CALL s_compute_fdfTsat(fp,dfdp,pressure,Tstar)
            ! Combining bisection and newton-raphson methods
            DO iter = 0, pTsatnewton_iter
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
                IF (DABS(delta) < pTsatnewton_eps) EXIT
                CALL s_compute_fdfTsat(fp,dfdp,pressure,Tstar)           
                PRINT *,'Tstar : ',Tstar,', fp : ',fp,', dfdp : ',dfdp
                IF (fp < 0.d0) THEN     !Maintain the bracket on the root
                   TstarL = Tstar
                ELSE
                   TstarH = Tstar
                END IF
                IF (iter .EQ. pTsatnewton_iter-1) THEN
                    PRINT *, 'Tsat : ',Tstar,', iter : ',iter
                    PRINT *, 'Tsat did not converge, stopping code'
                    CALL s_mpi_abort()
                END IF                  
            END DO
            f_Tsat = Tstar
            PRINT *, 'Tsat : ',f_Tsat,', iter : ',iter,', delta : ',delta
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

            ! Material 1
            cv1 = fluid_pp(1)%cv; q1 = fluid_pp(1)%qv;
            ! Material 2
            cv2 = fluid_pp(2)%cv; q2 = fluid_pp(2)%qv;
            ! Calculating coefficients, Eq. C.6, Pelanti 2014
            ap = rho0*cv1*cv2*((n2-1.d0)*(pstar+n1*pinf1)-(n1-1.d0)*(pstar+n2*pinf2))
            bp = E0*((n1-1.d0)*cv1*(pstar+pinf2) - (n2-1.d0)*cv2*(pstar+pinf1)) + &
                 rho0*((n2-1.d0)*cv2*q1*(pstar+pinf1) - (n1-1.d0)*cv1*q2*(pstar+pinf2)) + &
                 cv2*(pstar+pinf1)*(pstar+n2*pinf2) - cv1*(pstar+pinf2)*(pstar+n1*pinf1)
            dp = (q2-q1)*(pstar+pinf1)*(pstar+pinf2)
            ! Calculating the Tstar temperature, Eq. C.7, Pelanti 2014
            Tstar = (-bp + DSQRT(bp*bp - 4.d0*ap*dp))/(2.d0*ap)
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

        END SUBROUTINE s_compute_ptg_fdf !------------------------

        !>     The purpose of this subroutine is to determine the bracket of 
        !!         the pressure by finding the pressure at which b^2=4*a*c
        !!     @param p_star equilibrium pressure at the interface    
        SUBROUTINE s_compute_ptg_bracket(pstarA,pstarB,pstar,rho0,E0)

            REAL(KIND(0d0)), INTENT(OUT)   :: pstarA, pstarB
            REAL(KIND(0d0)), INTENT(IN)    :: pstar, rho0, E0
            REAL(KIND(0d0))                :: fA, fB, dfdp, Tstar
            REAL(KIND(0d0))                :: factor

            pstarA = 1.d0
            pstarB = pstarA
            CALL s_compute_ptg_fdf(fA,dfdp,pstarA,Tstar,rho0,E0)
            fB = fA
            factor = 10.d0
            DO WHILE ( fA*fB .GT. 0.d0 )
                  IF (ISNAN(Tstar)) RETURN
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
                        factor = factor*0.5d0
                  ELSE 
                        factor = 10.d0
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
            REAL(KIND(0d0))                ::  pstarA, pstarB
            REAL(KIND(0d0))                ::  delta, delta_old, fp, dfdp
            REAL(KIND(0d0))                ::  fL, fH, pstarL, pstarH
            INTEGER :: iter      !< Generic loop iterators
            ! Computing the bracket of the root solution
            CALL s_compute_ptg_bracket(pstarA,pstarB,pstar,rho0,E0)
            ! Computing f at lower and higher end of the bracket
            CALL s_compute_ptg_fdf(fL,dfdp,pstarA,Tstar,rho0,E0)
            CALL s_compute_ptg_fdf(fH,dfdp,pstarB,Tstar,rho0,E0)
            ! Establishing the direction of the descent to find zero
            IF(fL < 0.d0) THEN
                pstarL  = pstarA; pstarH  = pstarB;
            ELSE
                pstarL  = pstarB; pstarH  = pstarA;
            END IF
            pstar = 0.5d0*(pstarA+pstarB)
            delta_old = DABS(pstarB-pstarA)
            delta = delta_old
            CALL s_compute_ptg_fdf(fp,dfdp,pstar,Tstar,rho0,E0)
            ! Combining bisection and newton-raphson methods
            DO iter = 0, ptgnewton_iter
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
                IF (DABS(delta) < ptgnewton_eps) EXIT           ! Stopping criteria
                ! Updating to next iteration
                CALL s_compute_ptg_fdf(fp,dfdp,pstar,Tstar,rho0,E0)
                IF (fp < 0.d0) THEN !Maintain the bracket on the root
                   pstarL = pstar
                ELSE
                   pstarH = pstar
                END IF
            END DO

        END SUBROUTINE s_compute_ptg_pTrelax !-------------------------------

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

            !! terms  for K div(u)
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

                       IF (q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) .GT. palpha_epsL .OR. &
                           q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l) .LT. 1.d0-palpha_epsH) THEN
       
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
                                        pres_K_init(i) = -(1d0 - 1d-8)*pres_inf(i) + 1d-8
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
                                pres_relax = pres_relax - f_pres / df_pres
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
                                        !numerator   = gamma_min(i)*(pres_relax+pres_inf(i))
                                        !denominator = numerator + pres_K_init(i)-pres_relax
                                        !rho_K_s(i)  = q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l)/&
                                        !    MAX(q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l),sgm_eps)*& 
                                        !              numerator/denominator
                                        !drhodp      = q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l) / & 
                                        !    MAX(q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l),sgm_eps) * & 
                                        !    gamma_min(i)*(pres_K_init(i)+pres_inf(i)) / (denominator*denominator)
                                        !f_pres      = f_pres  + q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l) / rho_K_s(i)
                                        !df_pres     = df_pres - q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l) * & 
                                        !              drhodp / (rho_K_s(i)*rho_K_s(i))
                                        rho_K_s(i) = q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l) / &
                                                    MAX(q_cons_vf(i+adv_idx%beg-1)%sf(j,k,l),sgm_eps) &
                                                    * ((pres_relax+pres_inf(i)) / (pres_K_init(i) + &
                                                    pres_inf(i)))**(1d0/gamma_min(i))
                                        f_pres      = f_pres  + q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l) &
                                            / rho_K_s(i)
                                        df_pres     = df_pres - q_cons_vf(i+cont_idx%beg-1)%sf(j,k,l) &
                                            / (gamma_min(i)*rho_K_s(i)*(pres_relax+pres_inf(i)))
                                    END IF
                                END DO
                                !pres_relax = pres_relax - f_pres / df_pres
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
            REAL(KIND(0d0))                                   ::Z1, Z2, pI, C1, C2
            REAL(KIND(0d0))                                   ::        ap, bp, dp
            INTEGER :: iter      !< Generic loop iterators

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
                        ! Numerical correction of the volume fractions
                        IF (mpp_lim) THEN
                            CALL s_mixture_volume_fraction_correction(q_cons_vf, j, k, l )
                        END IF
                        ! Thermodynamic equilibrium relaxation procedure ================================
                        relax = .FALSE.
                        IF ( (q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) .GT. palpha_epsL ) .AND. &
                              q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) .LT. 1.d0-palpha_epsH ) relax = .TRUE.
                        !> ==============================================================================
                        !! STARTING THE RELAXATION PROCEDURE ============================================
                        !< ==============================================================================
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
                    END DO
                END DO
            END DO
        END SUBROUTINE s_infinite_p_relaxation ! ----------------

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
            REAL(KIND(0d0))                                   ::  cv1, cv2, q1, q2
            REAL(KIND(0d0))                                   ::        ap, bp, dp
            INTEGER :: iter      !< Generic loop iterators

            ! Material 1
            cv1 = fluid_pp(1)%cv; q1 = fluid_pp(1)%qv;
            ! Material 2
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
            LOGICAL :: relax, failed     !< Relaxation procedure determination variable
            !< Computing the constant saturation properties 

            DO j = 0, m
                DO k = 0, n
                    DO l = 0, p
                        ! Resetting the internal energy value and the relax and failed flags
                        rhoe = 0.d0
                        relax = .FALSE.
                        ! Numerical correction of the volume fractions
                        IF (mpp_lim) THEN
                            CALL s_mixture_volume_fraction_correction(q_cons_vf, j, k, l )
                        END IF

                        CALL s_convert_to_mixture_variables( q_cons_vf, rho, &
                                                             gamma, pi_inf,  &
                                                             Re, We, j, k, l )

                        IF ( (q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) .GT. palpha_epsL ) .AND. &
                              q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) .LT. 1.d0-palpha_epsH ) relax = .TRUE.

                        !> ==============================================================================
                        !! P RELAXATION =================================================================
                        !< ==============================================================================
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
                        relax = .FALSE.

                        IF ( (q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) .GT. palpha_epsL ) .AND. &
                              q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) .LT. 1.d0-palpha_epsH ) relax = .TRUE.
                        !> ==============================================================================
                        !! PT RELAXATION ================================================================
                        !< ==============================================================================
                        IF (relax) THEN
                            rhoalpha1 = q_cons_vf(cont_idx%beg)%sf(j,k,l)
                            rhoalpha2 = q_cons_vf(1+cont_idx%beg)%sf(j,k,l)
                            rhoe = 0.d0
                            DO i = 1, num_fluids
                                rhoe = rhoe + q_cons_vf(i+internalEnergies_idx%beg-1)%sf(j,k,l)
                            END DO
                            a1 = f_alpha1_ptrelax(rhoalpha1,rhoalpha2,rhoe)
                            ! Cell update of the volume fraction
                            q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l)  = a1
                            q_cons_vf(2+adv_idx%beg-1)%sf(j,k,l)  = 1.d0 - a1
                        END IF
                        CALL s_mixture_total_energy_correction(q_cons_vf, j, k, l )
                        relax = .FALSE.

                        IF ((q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) .GT. ptgalpha_epsL ) .AND. &
                             q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) .LT. 1.d0-ptgalpha_epsH ) relax = .TRUE.
                        !> ==============================================================================
                        !! CHECKING IF PTG RELAXATION IS NEEDED  ========================================
                        !< ==============================================================================
                        IF (relax) THEN
                           DO i = 1, num_fluids
                               rhoe = rhoe + q_cons_vf(i+internalEnergies_idx%beg-1)%sf(j,k,l) 
                           END DO                   
                           pres_relax = (rhoe - pi_inf)/gamma
                           IF ( ISNAN(pres_relax) ) THEN
                               PRINT *, 'pressure is NaN, stopping code'
                               CALL s_mpi_abort()
                           END IF
                           PRINT *,'a1 : ',q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l),',a2 :',q_cons_vf(2+adv_idx%beg-1)%sf(j,k,l)  
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
                           IF ( (q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) .LT. ptgalpha_epsL ) .AND. &
                              (Tk(1) .GT. Tsat) ) relax = .FALSE.
                           IF ( (q_cons_vf(1+adv_idx%beg-1)%sf(j,k,l) .GT. 1.d0-ptgalpha_epsH ) .AND. &
                              (Tk(1) .LT. Tsat) ) relax = .FALSE.
                           IF (Tk(1) .LT. Tsat) relax = .FALSE.
                        END IF
                        !> ==============================================================================
                        !! PTG RELAXATION PROCEDURE =====================================================
                        !< ==============================================================================
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
        !CALL s_mpi_abort()
        END SUBROUTINE s_infinite_ptg_relaxation ! -----------------------


END MODULE m_phasechange
