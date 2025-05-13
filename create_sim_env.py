from elastica._calculus import _isnan_check
from elastica import *
from elastica.rod.ribbon1D import Ribbon1D
from elastica.rod.linear_ribbon1D import LinearRibbon1D


from elastica._linalg import _batch_norm
from elastica._linalg import _batch_cross


import pickle
import numpy as np
import plotly.graph_objects as go
import pandas as pd
from collections import defaultdict
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from elastica.src_plotting_ribbon import *





def create_env(initial_position,d2_initial, position_sleeve, d2_sleeve, 
               target_force = 1e-1, k_p_input_force = 0.9, ramp_up_time = 0.03,
               dt = 2e-8, final_time = 0.1, 
               base_length = 50, thickness = 0.1, width = 5, 
               density = 1.017e-7, nu = 27e1, E = 2.77e3, poisson_ratio = 0.34, num_lagrange = 1e2
              ):

    class Ribbon(BaseSystemCollection, Contact, Forcing, Constraints, CallBacks, Damping):
        pass
    
    Sim_env = Ribbon()

    n_elem = initial_position.shape[1] - 1    
    position_diff = initial_position[..., 1:] - initial_position[..., :-1]
    rest_lengths = _batch_norm(position_diff)

    d3 = position_diff / rest_lengths
    d1 = _batch_cross(d2_initial,d3)
    d1 = d1 / _batch_norm(d1)

    
    directors = np.stack([d1,d2_initial,d3],axis=0)
    
    
    shearable_rod = Ribbon1D.straight_ribbon(
        n_elem,
        np.zeros((3,)),
        d3[:,0],
        d1[:,0],
        base_length,
        thickness,
        width,
        density,
        youngs_modulus=E,
        shear_modulus=E*num_lagrange, # influence the shearability of the ribbon
        poisson_ratio=0.34,
        position = initial_position,
        directors = directors,
    )
    
    
    Sim_env.append(shearable_rod)
    
    
    """ Add damping """
    Sim_env.dampen(shearable_rod).using(
        AnalyticalLinearDamper,
        damping_constant=nu,
        time_step=dt,
    )
    
    
    # Impact steady state !!! Add error when comparing with auto
    #Rod.dampen(shearable_rod).using(
    #         LaplaceDissipationFilter,
    #         filter_order=2,   # order of the filter (geater order means less damping)
    #    )
    

    """ Set up boundary conditions """
    Sim_env.constrain(shearable_rod).using(
        GeneralConstraint,
        constrained_position_idx=(0,),
        constrained_director_idx=(0,),
        translational_constraint_selector= np.abs(d3[:,0]) != np.max(np.abs(d3[:,0])),  # to block only the 2 direction that not the direction of d3 at the base
        rotational_constraint_selector=np.array([True, True, True])  # Fix all rotations
    )
    
    
    Sim_env.constrain(shearable_rod).using(
        GeneralConstraint,
        constrained_position_idx=(-1,),
        constrained_director_idx=(-1,),
        translational_constraint_selector=np.array([True, True, True]),  
        rotational_constraint_selector=np.array([False, False, False]) 
    )
    
    
    Sim_env.add_forcing_to(shearable_rod).using(
        ControledPushForce, target_force, d3[:,0], ramp_up_time, k_p = k_p_input_force
    )
    
    
    """ Set contact """

    position_diff = position_sleeve[..., 1:] - position_sleeve[..., :-1]
    rest_lengths = _batch_norm(position_diff)
    d3 = position_diff / rest_lengths
    d1 = _batch_cross(d2_sleeve,d3)
    d1 = d1 / _batch_norm(d1)


    sleeve = Sleeve(position_sleeve, d1)
    Sim_env.append(sleeve)


    Sim_env.detect_contact_between(shearable_rod, sleeve).using(
        RibbonSleeveContact,
        k=3*1.16e-3, #k = E_0 = 3*G_0 
        nu=0,
        alpha = -16.6,
        poisson_ratio=0.5,
        N_quad = 3,
    )
    
    
    class RibbonCallBack(CallBackBaseClass):
        """
        Call back function for Ribbon object
        """
    
        def __init__(self, step_skip: int, callback_params: dict):
            CallBackBaseClass.__init__(self)
            self.every = step_skip
            self.callback_params = callback_params
    
        def make_callback(self, system, time, current_step: int):
    
            if current_step % self.every == 0:
    
                self.callback_params["time"].append(time)
                self.callback_params["step"].append(current_step)
                self.callback_params["position"].append(system.position_collection.copy())
                self.callback_params["velocity"].append(system.velocity_collection.copy())
                self.callback_params["curvature"].append(system.kappa.copy())
                self.callback_params["sigma"].append(system.sigma.copy())
                self.callback_params["internal_stress"].append(system.internal_stress.copy())
                self.callback_params["internal_forces"].append(system.internal_forces.copy())
                self.callback_params["external_forces"].append(system.external_forces.copy())
                self.callback_params["internal_couple"].append(system.internal_couple.copy())
                self.callback_params["directors"].append(system.director_collection.copy())
                self.callback_params["dilatation"].append(system.dilatation.copy())
                self.callback_params["tangents"].append(system.tangents.copy())
                
            
                
                return
                
    class SleeveCallBack(CallBackBaseClass):
        """
        Call back function for Sleeve object
        """
    
        def __init__(self, step_skip: int, callback_params: dict):
            CallBackBaseClass.__init__(self)
            self.every = step_skip
            self.callback_params = callback_params
    
        def make_callback(self, system, time, current_step: int):
    
            if current_step % self.every == 0:
    
                self.callback_params["time"].append(time)
                self.callback_params["step"].append(current_step)
                self.callback_params["response_force_sleeve"].append(system.response_force_sleeve.copy())
                self.callback_params["displacement_sleeve"].append(system.displacement_sleeve.copy())
                self.callback_params["response_couple_sleeve"].append(system.response_couple_sleeve.copy())
                self.callback_params["rotation_sleeve"].append(system.rotation_sleeve.copy())
                
                return
                
    
    step_skip= int(final_time/dt/1000)
    
    ribbon_output_list = defaultdict(list)
    sleeve_output_list = defaultdict(list)
    
    Sim_env.collect_diagnostics(shearable_rod).using(
        RibbonCallBack, step_skip=step_skip, callback_params=ribbon_output_list
    )
    
    Sim_env.collect_diagnostics(sleeve).using(
        SleeveCallBack, step_skip=step_skip, callback_params=sleeve_output_list
    )

    Sim_env.finalize()

    return ribbon_output_list, sleeve_output_list, Sim_env