#!/usr/bin/env python3

import json

print(json.dumps({
    # Computational Domain Parameters ==========================================
    'x_domain%beg'                 : -1.E+00,
    'y_domain%beg'                 : -1.E+00,
    'z_domain%beg'                 : -1.E+00,
    'x_domain%end'                 : +1.E+00,
    'y_domain%end'                 : +1.E+00,
    'z_domain%end'                 : +1.E+00,
    'm'                            : 200,
    'n'                            : 200,
    'p'                            : 200,
    'dt'                           : 4.E-08,
    't_step_start'                 : 0,
    't_step_stop'                  : 1,
    't_step_save'                  : 1,
    # ==========================================================================

    # Simulation Algorithm Parameters ==========================================
    'num_patches'                  : 2,
    'model_eqns'                   : 2,
    'num_fluids'                   : 2,
    'time_stepper'                 : 3,
    'weno_vars'                    : 2,
    'weno_order'                   : 5,
    'weno_eps'                     : 1.E-16,
    'mapped_weno'                  : 'T',
    'null_weights'                 : 'F',
    'mp_weno'                      : 'F',
    'riemann_solver'               : 2,
    'wave_speeds'                  : 1,
    'avg_state'                    : 2,
    'bc_x%beg'                     : -1,
    'bc_y%beg'                     : -1,
    'bc_z%beg'                     : -1,
    'bc_x%end'                     : -1,
    'bc_y%end'                     : -1,
    'bc_z%end'                     : -1,
    # ==========================================================================

    # Formatted Database Files Structure Parameters ============================
    'format'                       : 1,
    'precision'                    : 2,
    'prim_vars_wrt'                :'T',
    'parallel_io'                  :'F',
    # ==========================================================================

    # Patch 1: High pressured water ============================================
    'patch_icpp(1)%geometry'       : 9,
    'patch_icpp(1)%x_centroid'     : 0.0E+00,
    'patch_icpp(1)%y_centroid'     : 0.0E+00,
    'patch_icpp(1)%z_centroid'     : 0.0E+00,
    'patch_icpp(1)%length_x'       : 2.E+00,
    'patch_icpp(1)%length_y'       : 2.E+00,
    'patch_icpp(1)%length_z'       : 2.E+00,
    'patch_icpp(1)%vel(1)'         : 0.E+00,
    'patch_icpp(1)%vel(2)'         : 0.E+00,
    'patch_icpp(1)%vel(3)'         : 0.E+00,
    'patch_icpp(1)%pres'           : 1.E+09,
    'patch_icpp(1)%alpha_rho(1)'   : 1000.E+00,
    'patch_icpp(1)%alpha_rho(2)'   : 0.,
    'patch_icpp(1)%alpha(1)'       : 1.E+00,
    'patch_icpp(1)%alpha(2)'       : 0.E+00,
    # ==========================================================================

    # Patch 2: STL =============================================================
    'patch_icpp(2)%geometry'       : 20,
    'patch_icpp(2)%stl%file'       : 'Annulus.stl',
    'patch_icpp(2)%stl%offset(1)'  : -0.25,
    'patch_icpp(2)%stl%offset(2)'  : -0.25,
    'patch_icpp(2)%stl%offset(3)'  : -0.25,
    'patch_icpp(2)%stl%scale(1)'   : 1.0/20.E+00,
    'patch_icpp(2)%stl%scale(2)'   : 1.0/20.E+00,
    'patch_icpp(2)%stl%scale(3)'   : 1.0/20.E+00,
    'patch_icpp(2)%vel(1)'         : 0.E+00,
    'patch_icpp(2)%vel(2)'         : 0.E+00,
    'patch_icpp(2)%vel(3)'         : 0.E+00,
    'patch_icpp(2)%pres'           : 1.E+05,
    'patch_icpp(2)%alpha_rho(1)'   : 0.,
    'patch_icpp(2)%alpha_rho(2)'   : 50.E+0,
    'patch_icpp(2)%alpha(1)'       : 0,
    'patch_icpp(2)%alpha(2)'       : 1.,
    # ==========================================================================

    # Fluids Physical Parameters ===============================================
    'fluid_pp(1)%gamma'            : 1.E+00/(4.4E+00-1.E+00),
    'fluid_pp(1)%pi_inf'           : 4.4E+00*6.E+08/(4.4E+00-1.E+00),
    'fluid_pp(2)%gamma'            : 1.E+00/(1.4E+00-1.E+00),
    'fluid_pp(2)%pi_inf'           : 0.E+00,
    # ==========================================================================
}))

# ==============================================================================
