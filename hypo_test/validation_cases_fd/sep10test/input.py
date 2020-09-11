#!/usr/bin/python
import math

#Numerical setup
#Nx      = 499
Nx      = 201
dx      = 1./(1.*(Nx+1))

Tend    = 240E-06
#Tend = 64e-06
Nt      = 1000
mydt    = Tend/(1.*Nt)

# Command to navigate between directories
from os import chdir

# Command to acquire directory path
from os.path import dirname

# Command to acquire script name and module search path
from sys import argv, path

# Navigating to script directory
if len(dirname(argv[0])) != 0: chdir(dirname(argv[0]))

# Adding master_scripts directory to module search path
mfc_dir = '../../../src'; path[:0] = [mfc_dir + '/master_scripts']

# Command to execute the MFC components
from m_python_proxy import f_execute_mfc_component

# Selecting MFC component
comp_name = argv[1].strip()

# Serial or parallel computational engine

engine = 'serial'
if (comp_name == 'pre_process'): engine = 'serial'
# ==============================================================================

# Case Analysis Configuration ==================================================


# Configuring case dictionary
case_dict =                                                                     \
    {                                                                           \
                    # Logistics ================================================
                    'case_dir'                     : '\'.\'',                   \
                    'run_time_info'                : 'T',                       \
                    'nodes'                        : 1,                         \
                    'ppn'                          : 1,                         \
                    'queue'                        : 'normal',                  \
                    'walltime'                     : '24:00:00',                \
                    'mail_list'                    : '',                        \
                    # ==========================================================
                                                                                \
                    # Computational Domain Parameters ==========================
                    'x_domain%beg'                 : 0.E+00,                    \
                    'x_domain%end'                 : 1.E+00,                    \
                    'm'                            : Nx,                        \
                    'n'                            : 0,                         \
                    'p'                            : 0,                         \
                    'dt'                           : mydt,                      \
                    't_step_start'                 : 0,                         \
                    't_step_stop'                  : int(Nt),                        \
#                    't_step_save'                  : int(math.ceil(Nt/100.)),    \
                    't_step_save'                  : int(math.ceil(Nt/100.)),    \
		    # ==========================================================
                                                                                \
                    # Simulation Algorithm Parameters ==========================
                    'num_patches'                  : 2,                        \
                    'model_eqns'                   : 2,                        \
                    'alt_soundspeed'               : 'F',                      \
                    'num_fluids'                   : 2,                        \
		    'adv_alphan'                   : 'T',                      \
		    'mpp_lim'                      : 'F',                      \
		    'mixture_err'                  : 'F',                      \
		    'time_stepper'                 : 3,                        \
                    'weno_vars'                    : 2,                        \
                    'weno_order'                   : 5,                        \
                    'weno_eps'                     : 1.E-16,                   \
                    'char_decomp'                  : 'F',                      \
                    'mapped_weno'                  : 'F',                      \
                    'null_weights'                 : 'F',                      \
                    'mp_weno'                      : 'F',                      \
		    'riemann_solver'               : 1,                        \
                    'wave_speeds'                  : 1,                        \
                    'avg_state'                    : 2,                        \
                    'commute_err'                  : 'F',                      \
                    'split_err'                    : 'F',                      \
                    'bc_x%beg'                     : -3,                       \
                    'bc_x%end'                     : -3,                       \

                    'weno_Re_flux'                 : 'F',                      \
                    # ==========================================================

                    # Turning on Hypoelasticity ================================
#                    'tvd_rhs_flux'                 : ,                      \
                    'hypoelasticity'               : 'T',                      \
                    # ==========================================================
                                                                               \
                    # Formatted Database Files Structure Parameters ============
                    'format'                       : 1,                        \
                    'precision'                    : 2,                        \
                    'prim_vars_wrt'                :'T',                       \
		    'parallel_io'                  :'F',                       \
		    # ==========================================================
                                                                                
		    # Patch 1 L ================================================
                    'patch_icpp(1)%geometry'       : 1,                        \
                    'patch_icpp(1)%x_centroid'     : 0.5,                     \
                    'patch_icpp(1)%length_x'       : 1.,                      \
                    'patch_icpp(1)%vel(1)'         : 0.0,                      \
#                    'patch_icpp(1)%pres'           : 1E+09,                    \
                    'patch_icpp(1)%pres'           : 1E+09,                     \
                    'patch_icpp(1)%alpha_rho(1)'   : 1000,                     \
                    'patch_icpp(1)%alpha_rho(2)'   : 0.,                        \
                    'patch_icpp(1)%alpha(1)'       : 1.0,                       \
                    'patch_icpp(1)%alpha(2)'       : 0.,                         \
                    'patch_icpp(1)%tau_e(1)'       : 0.0,                      \
                    # ==========================================================

                    # Patch 2 R ================================================
                    'patch_icpp(2)%geometry'       : 1,                     \
                    'patch_icpp(2)%x_centroid'     : 0.85,                   \
                    'patch_icpp(2)%length_x'       : 0.3,                    \
                    'patch_icpp(2)%alter_patch(1)' : 'T',                       \
                    'patch_icpp(2)%vel(1)'         : 0,                    \
                    'patch_icpp(2)%pres'           : 1E+05,                    \
                    'patch_icpp(2)%alpha_rho(1)'   : 0,                    \
                    'patch_icpp(2)%alpha_rho(2)'   : 50,               \
#                    'patch_icpp(2)%alpha_rho(2)'   : 1000,               \
                    'patch_icpp(2)%alpha(1)'       : 0,                \
                    'patch_icpp(2)%alpha(2)'       : 1.0,
                   'patch_icpp(2)%tau_e(1)'       : 0.0,                \

                    # ==========================================================

                    # Patch 3 adding a smoothing region ========================
#                    'patch_icpp(3)%geometry'       : 1,                        \
#                    'patch_icpp(3)%x_centroid'     : 0.695,                       \
#                    'patch_icpp(3)%length_x'       : 0.01,                       \
#                    'patch_icpp(3)%alter_patch(1)' : 'T',                       \
#                    'patch_icpp(3)%vel(1)'         : 0,                         \
#                    'patch_icpp(3)%pres'           : 1E+05,                     \
#                    'patch_icpp(3)%alpha_rho(1)'   : 500,                    \
#                    'patch_icpp(3)%alpha_rho(2)'   : 25,               \
#                    'patch_icpp(3)%alpha(1)'       : 0.5,                \
#                    'patch_icpp(3)%alpha(2)'       : 0.5,
#                    'patch_icpp(3)%tau_e(1)'       : 0.0,                \

                    # Fluids Physical Parameters ===============================
#                    'fluid_pp(1)%gamma'            : 1.E+00/(1.4E+00-1.E+00),   \
#                    'fluid_pp(1)%pi_inf'           : 0.,                          \

                    'fluid_pp(1)%gamma'            : 1.E+00/(4.4E+00-1.E+00),   \
                    'fluid_pp(1)%pi_inf'           : 4.4E+00*6.E+08/(4.4E+00 - 1.E+00), \
#                    'fluid_pp(1)%G'                : 0.,                        \
                    'fluid_pp(1)%G'                : 1.E+09,                       \
#                    'fluid_pp(1)%Re(1)'            : 1.E+10,                   \
#                    'fluid_pp(1)%Re(2)'            : 50.E-03,                       \

                    'fluid_pp(2)%gamma'            : 1.E+00/(1.4E+00-1.E+00),   \
                    'fluid_pp(2)%pi_inf'           : 0.,                          \
                    'fluid_pp(2)%G'                : 0.,                        \

#                    'fluid_pp(2)%gamma'            : 1.E+00/(4.4E+00-1.E+00),   \
#                    'fluid_pp(2)%pi_inf'           : 4.4E+00*6.E+08/(4.4E+00 - 1.E+00), \
#                    'fluid_pp(2)%G'                : 1.E+09,                        \

#                    'fluid_pp(2)%gamma'            : 1.E+00/(4.4E+00-1.E+00),   \
#                    'fluid_pp(2)%pi_inf'           : 4.4E+00*6.E+08/(4.4E+00 - 1.E+00), \
#                    'fluid_pp(1)%G'                : 0.,                        \
#                    'fluid_pp(2)%G'                : 1.E+09,                       \
#                    'fluid_pp(2)%Re(1)'            : 50.E-03,                   \
#                    'fluid_pp(2)%Re(2)'            : 50.E-03,                   \
	            # ==========================================================
    }

# Executing MFC component
f_execute_mfc_component(comp_name, case_dict, mfc_dir, engine)

# ==============================================================================
