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





def create_env(initial_position,d1_initial, position_sleeve, d1_sleeve, 
               target_force = 1e-1, k_p_input_force = 0.9, ramp_up_time = 0.03,
               dt = 2e-8, final_time = 0.1, 
               base_length = 50, thickness = 0.1, width = 5, 
               density = 1.017e-7, nu = 27e1, E = 2.77e3, poisson_ratio = 0.34, num_lagrange = 1e2
              ):
    """
    Create a PyElastica simulation environment for ribbon-sleeve contact dynamics.
    
    This function sets up a complete simulation environment for studying the interaction
    between a flexible ribbon and a constraining sleeve. The ribbon is modeled as a 
    1D Cosserat rod with finite width and thickness, while the sleeve provides contact
    constraints along the ribbon's path.
    
    Parameters
    ----------
    initial_position : np.ndarray, shape (3, n_nodes)
        Initial 3D positions of ribbon nodes. Each column represents a node position
        in Cartesian coordinates [x, y, z].
        
    d1_initial : np.ndarray, shape (3, n_elem)
        Initial director vectors d1 for each ribbon element. These vectors define
        the ribbon's normal direction (perpendicular to the ribbon surface).
        
    position_sleeve : np.ndarray
        3D positions defining the sleeve geometry/constraints.
        
    d1_sleeve : np.ndarray
        Director vectors for the sleeve orientation.
        
    target_force : float, default=1e-1
        Target force magnitude [N] applied to drive the ribbon motion.
        
    k_p_input_force : float, default=0.9
        Proportional gain for the controlled force application. Controls how
        aggressively the force tracks the target value.
        
    ramp_up_time : float, default=0.03
        Time duration [s] over which forces are gradually increased from zero
        to target values to avoid numerical instabilities.
        
    dt : float, default=2e-8
        Time step size [s] for numerical integration. Small values ensure
        stability but increase computational cost.
        
    final_time : float, default=0.1
        Total simulation time [s].
        
    base_length : float, default=50
        Total length [mm] of the undeformed ribbon.
        
    thickness : float, default=0.1
        Ribbon thickness [mm] in the d2 direction (out-of-plane).
        
    width : float, default=5
        Ribbon width [mm] in the d1 direction (in-plane normal).
        
    density : float, default=1.017e-7
        Material density [kg/mm³]. Default value corresponds to typical
        polymer materials.
        
    nu : float, default=27e1
        Damping coefficient [Pa·s] for analytical linear damping. Controls
        energy dissipation and steady-state behavior.
        
    E : float, default=2.77e3
        Young's modulus [Pa] of the ribbon material. Defines bending stiffness.
        
    poisson_ratio : float, default=0.34
        Poisson's ratio [-] of the material. Affects cross-sectional deformation.
        
    num_lagrange : float, default=1e2
        Lagrange multiplier scaling factor that influences ribbon shearability.
        Higher values make the ribbon less shearable.
        
    Returns
    -------
    ribbon_output_list : defaultdict
        Dictionary containing time-series data for ribbon properties:
        - 'time': simulation time points
        - 'position': node positions over time
        - 'velocity': node velocities
        - 'curvature': element curvatures (kappa)
        - 'sigma': shear strains
        - 'internal_stress': internal stress distribution
        - 'internal_forces': internal force vectors
        - 'external_forces': applied external forces
        - 'internal_couple': internal moments
        - 'directors': director field evolution
        - 'dilatation': element stretch/compression
        - 'tangents': tangent vectors along the ribbon
        
    sleeve_output_list : defaultdict
        Dictionary containing sleeve interaction data:
        - 'response_force_sleeve': contact forces from sleeve
        - 'displacement_sleeve': sleeve displacement
        - 'response_couple_sleeve': contact moments from sleeve
        - 'rotation_sleeve': sleeve rotation
        
    Sim_env : Ribbon (BaseSystemCollection)
        Configured PyElastica simulation environment ready for time integration.
        Call `integrate()` method to run the simulation.
        
    Notes
    -----
    The function creates a ribbon with the following boundary conditions:
    - Base (index 0): Constrained in directions perpendicular to d3, all rotations fixed
    - Tip (index -1): All translations constrained, rotations free
    
    The ribbon-sleeve contact uses a Ogden contact model with:
    - Contact stiffness k = 3*1.16e-3 Pa (3*G_0)
    - Contact damping nu = 0
    - Contact parameter alpha = -16.6
    - Number of quadrature points N_quad = 3
    """

    class Ribbon(BaseSystemCollection, Contact, Forcing, Constraints, CallBacks, Damping):
        pass
    
    Sim_env = Ribbon()

    n_elem = initial_position.shape[1] - 1    
    position_diff = initial_position[..., 1:] - initial_position[..., :-1]
    rest_lengths = _batch_norm(position_diff)
    d3 = position_diff / rest_lengths
    
    d2 = _batch_cross(d3,d1_initial)
    d2 = d2 / _batch_norm(d2)

    d1 = _batch_cross(d2,d3)
    d1 = d1 / _batch_norm(d1)    
    
    directors = np.stack([d1,d2,d3],axis=0)
    
    
    ribbon = Ribbon1D.straight_ribbon(
        n_elem,
        np.zeros((3,)),
        d3[:,0],
        d1[:,0],
        base_length,
        thickness,
        width,
        density,
        youngs_modulus=E,
        shear_modulus=E*num_lagrange, # influence the shearability and strechability of the ribbon
        poisson_ratio=0.34,
        position = initial_position,
        directors = directors,
    )
    
    
    Sim_env.append(ribbon)
    
    
    """ Add damping """
    Sim_env.dampen(ribbon).using(
        AnalyticalLinearDamper,
        damping_constant=nu,
        time_step=dt,
    )
    
    
    # When tested, it impact steady state !!! Add error when comparing with auto
    #Rod.dampen(ribbon).using(
    #         LaplaceDissipationFilter,
    #         filter_order=2,   # order of the filter (geater order means less damping)
    #    )
    

    """ Set up boundary conditions """
    Sim_env.constrain(ribbon).using(
        GeneralConstraint,
        constrained_position_idx=(0,),
        constrained_director_idx=(0,),
        translational_constraint_selector= np.abs(d3[:,0]) != np.max(np.abs(d3[:,0])),  # to block only the 2 direction that are not the direction of d3 at the base
        rotational_constraint_selector=np.array([True, True, True])  # Fix all rotations
    )
    
    
    Sim_env.constrain(ribbon).using(
        GeneralConstraint,
        constrained_position_idx=(-1,),
        constrained_director_idx=(-1,),
        translational_constraint_selector=np.array([True, True, True]),  
        rotational_constraint_selector=np.array([False, False, False]) 
    )
    
    
    Sim_env.add_forcing_to(ribbon).using(
        ControledPushForce, target_force, d3[:,0], ramp_up_time, k_p = k_p_input_force
    )
    
    
    """ Set contact """


    sleeve = Sleeve(position_sleeve, d1_sleeve)
    Sim_env.append(sleeve)


    Sim_env.detect_contact_between(ribbon, sleeve).using(
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
    
    Sim_env.collect_diagnostics(ribbon).using(
        RibbonCallBack, step_skip=step_skip, callback_params=ribbon_output_list
    )
    
    Sim_env.collect_diagnostics(sleeve).using(
        SleeveCallBack, step_skip=step_skip, callback_params=sleeve_output_list
    )

    Sim_env.finalize()

    return ribbon_output_list, sleeve_output_list, Sim_env



######################################################################################################

def create_env_torsional_buckling(initial_position,d1_initial,
               tip_force = np.array([0.0, 2.0e-2, 2.0e-4]), ramp_up_time = 0.05,
               dt = 2e-8, final_time = 0.1, 
               base_length = 50, thickness = 0.1, width = 5, 
               density = 1.017e-7, nu = 27e1, E = 2.77e3, poisson_ratio = 0.34, num_lagrange = 1e2
              ):
    """
    Create a PyElastica simulation environment for studying torsional buckling of ribbons.
    
    This function sets up a ribbon simulation specifically designed to perform torsional buckling simulations. 
    The ribbon is fully constrained at the base and
    subjected to combined forces at the tip, which can induce buckling instabilities.
    
    Parameters
    ----------
    initial_position : np.ndarray, shape (3, n_nodes)
        Initial 3D positions of ribbon nodes. Each column represents a node position
        in Cartesian coordinates [x, y, z]. Typically starts as a straight configuration.
        
    d1_initial : np.ndarray, shape (3, n_elem)
        Initial director vectors d1 for each ribbon element. These define the ribbon's
        normal direction and initial orientation field.
        
    tip_force : np.ndarray, shape (3,), default=[0.0, 2.0e-2, 2.0e-4]
        Force vector [N] applied at the ribbon tip. Components represent:
        - [0]: Force in x-direction
        - [1]: Force in y-direction (dominant lateral force)
        - [2]: Force in z-direction (small axial component)
        
    ramp_up_time : float, default=0.05
        Time duration [s] for gradual force application. Longer ramp-up times
        help capture quasi-static behavior and avoid dynamic effects.
        
    dt : float, default=2e-8
        Time step size [s] for numerical integration. 
        
    final_time : float, default=0.1
        Total simulation time [s]. Should be long enough to reach steady-state.
        
    base_length : float, default=50
        Total undeformed length [mm] of the ribbon.
        
    thickness : float, default=0.1
        Ribbon thickness [mm].
        
    width : float, default=5
        Ribbon width [mm].
        
    density : float, default=1.017e-7
        Material density [kg/mm³].
        
    nu : float, default=27e1
        Damping coefficient [-]. Numerical damping that stabilize the simulation
        
    E : float, default=2.77e3
        Young's modulus [Pa]. Determines bending stiffness.
        
    poisson_ratio : float, default=0.34
        Poisson's ratio [-]. Affects the relationship between bending and
        torsional stiffnesses.
        
    num_lagrange : float, default=1e2
        Lagrange multiplier scaling factor that influences ribbon shearability.
        Higher values make the ribbon less shearable.
        
    Returns
    -------
    ribbon_output_list : defaultdict
        Dictionary containing comprehensive time-series data:
        - 'time': simulation time points
        - 'step': integration step numbers
        - 'position': node positions (tracks buckling shape evolution)
        - 'velocity': node velocities
        - 'curvature': element curvatures (kappa) - important for buckling analysis
        - 'sigma': shear strains
        - 'internal_stress': stress distribution along ribbon
        - 'internal_forces': internal force vectors
        - 'external_forces': applied external forces
        - 'internal_couple': internal moments (critical for torsional buckling)
        - 'directors': director field evolution (shows twist development)
        - 'dilatation': axial strain distribution
        - 'tangents': tangent vectors (track centerline orientation)
        
    Sim_env : Ribbon (BaseSystemCollection)
        Configured simulation environment ready for integration. The ribbon
        has clamped boundary conditions at the base.
        
    Notes
    -----
    Boundary Conditions:
    - Base (index 0): Fully clamped (all translations and rotations constrained)
    - Tip (index -1): Free to move and rotate
    
    The force is applied using EndpointForces with gradual ramp-up for numerical stability.
    """
    class Ribbon(BaseSystemCollection, Contact, Forcing, Constraints, CallBacks, Damping):
        pass
    
    Sim_env = Ribbon()

    n_elem = initial_position.shape[1] - 1    
    position_diff = initial_position[..., 1:] - initial_position[..., :-1]
    rest_lengths = _batch_norm(position_diff)
    d3 = position_diff / rest_lengths
    
    d2 = _batch_cross(d3,d1_initial)
    d2 = d2 / _batch_norm(d2)

    d1 = _batch_cross(d2,d3)
    d1 = d1 / _batch_norm(d1)  

    
    directors = np.stack([d1,d2,d3],axis=0)
    
    
    ribbon = Ribbon1D.straight_ribbon(
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
    
    
    Sim_env.append(ribbon)
    
    
    """ Add damping """
    Sim_env.dampen(ribbon).using(
        AnalyticalLinearDamper,
        damping_constant=nu,
        time_step=dt,
    )
    
    
    # When tested, it impact steady state !!! Add error when comparing with auto
    #Rod.dampen(ribbon).using(
    #         LaplaceDissipationFilter,
    #         filter_order=2,   # order of the filter (geater order means less damping)
    #    )
    

    """ Set up boundary conditions """
    Sim_env.constrain(ribbon).using(
        GeneralConstraint,
        constrained_position_idx=(0,),
        constrained_director_idx=(0,),
        translational_constraint_selector= np.array([True, True, True]) , # Fix all displacement
        rotational_constraint_selector=np.array([True, True, True])  # Fix all rotations
    )
    
    
    
    origin_force = np.array([0.0,0.0,0.0])
    Sim_env.add_forcing_to(ribbon).using(
        EndpointForces, origin_force, tip_force, ramp_up_time=ramp_up_time
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
                
    
    step_skip= int(final_time/dt/1000)
    
    ribbon_output_list = defaultdict(list)
    
    Sim_env.collect_diagnostics(ribbon).using(
        RibbonCallBack, step_skip=step_skip, callback_params=ribbon_output_list
    )

    Sim_env.finalize()

    return ribbon_output_list, Sim_env