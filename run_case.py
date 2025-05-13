from elastica._calculus import _isnan_check
from elastica.timestepper import extend_stepper_interface
from elastica import *
from elastica.rod.ribbon1D import Ribbon1D
from elastica.rod.linear_ribbon1D import LinearRibbon1D


from elastica._linalg import _batch_norm
from elastica._linalg import _batch_cross

import os
from elastica._rotations import _get_rotation_matrix

from itertools import groupby

import pickle
import numpy as np
import plotly.graph_objects as go
import pandas as pd
from collections import defaultdict
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from elastica.src_plotting_ribbon import *




class LinearRod(BaseSystemCollection, Contact, Forcing, Constraints, CallBacks, Damping):
    pass

Rod = LinearRod()

n_elem = 30


start = np.array([0.0, 0.0, 0.0])
direction = np.array([1.0, 0.0, 0.0])
normal = np.array([0.0, 0.0, 1.0])
base_length = 50
thickness = 0.1
width = 5

base_area = width*thickness
density = 1.017e-7
nu = 3e0
E = 2.77e3
num_lagrange = 1e2
poisson_ratio = 0.34

w0 = (np.pi/base_length)**2*np.sqrt(E*width*thickness**3/12/(density*base_area))
print("we have 2*w0 =", 2*w0)
print("we have P = ", 2*np.pi/w0)
dt = 1.0e-8

s = np.linspace(0,1,n_elem+1)
s_bc = np.concatenate((np.array([0]),s[:-1]))
w = 2*np.pi
A = 0.1

#position = base_length*np.array([s-A*np.sin(w*s_bc),np.zeros(n_elem+1),A*np.sin(w*s_bc)])

position = base_length*np.array([s,np.zeros(n_elem+1),A*np.sin(w*s+np.pi/2)-A])

d2 = np.array([np.zeros(n_elem),np.ones(n_elem),np.zeros(n_elem)])

position_diff = position[..., 1:] - position[..., :-1]
rest_lengths = _batch_norm(position_diff)
d3 = position_diff / rest_lengths

d1 = _batch_cross(d2,d3)

d1 = d1 / _batch_norm(d1)


directors = np.stack([d1,d2,d3],axis=0)


origin_force = np.array([1.0e-5, 0.0, 0.0])*0
end_force = np.array([0.0, 2.0e-2*0, 2.0e-4])*5*0 # (x,y,z)
ramp_up_time_force = 0.05

origin_torque = np.array([0, 0.0, 0.0])
end_torque = np.array([0.0, 5.0e-4, 0.0])*0
ramp_up_time_torque = 1.0

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
    position = position,
    directors = directors,
    #nu_for_torques=damp_coefficient*((radius_mean/radius_base)**4),
)



Rod.append(shearable_rod)


""" Add damping """
Rod.dampen(shearable_rod).using(
    AnalyticalLinearDamper,
    damping_constant=nu,
    time_step=dt,
)


# Impact steady state !!! Add error when comparing with auto
Rod.dampen(shearable_rod).using(
         LaplaceDissipationFilter,
         filter_order=2,   # order of the filter (geater order means less damping)
    )




""" Set up boundary conditions """
Rod.constrain(shearable_rod).using(
    GeneralConstraint,
    constrained_position_idx=(0,),
    constrained_director_idx=(0,),
    #translational_constraint_selector=np.array([False, True, True]),  # Allow X movement, fix Y & Z
    translational_constraint_selector=np.array([True, True, True]),  # Allow X movement, fix Y & Z
    rotational_constraint_selector=np.array([True, True, True])  # Fix all rotations
)


#Rod.constrain(shearable_rod).using(
#    GeneralConstraint,
#    constrained_position_idx=(-1,),
#    constrained_director_idx=(-1,),
#    rod_frame_bool = True,
    #translational_constraint_selector=np.array([False, False, True]),  
#    translational_constraint_selector=np.array([False, False, False]),  
#    rotational_constraint_selector=np.array([False, False, False]) 
#)


Rod.add_forcing_to(shearable_rod).using(
    EndpointForces, origin_force, end_force, ramp_up_time=ramp_up_time_force
)

#Rod.add_forcing_to(shearable_rod).using(
#    EndpointTorques, origin_torque, end_torque, ramp_up_time=ramp_up_time_torque
#)


""" Set contact """
s = np.linspace(0,1,n_elem)
w = np.pi/2
#d1 = np.array([np.zeros(n_elem),np.sin(w*s),-np.cos(w*s)])

s = np.linspace(0,1,n_elem+1)
s_bc = np.concatenate((np.array([0]),s[:-1]))
w = 2*np.pi
A = 0.1

position = base_length*np.array([s_bc,np.zeros(n_elem+1),A*np.sin(w*s_bc+np.pi/2)-A])

sleeve = Sleeve(position, d1)


Rod.append(sleeve)

Rod.detect_contact_between(shearable_rod, sleeve).using(
    RibbonSleeveContact,
    k=3*1.16e-3, #k = E_0 = 3*G_0 
    nu=0,
    alpha = -16.6,
    poisson_ratio=0.5,
    width_elem = width,
    length_elem = base_length/n_elem,
    N_quad = 3,
)



class RodCallBack(CallBackBaseClass):
    """
    Call back function
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
            self.callback_params["internal_couple"].append(system.internal_couple.copy())
            self.callback_params["directors"].append(system.director_collection.copy())
            self.callback_params["dilatation"].append(system.dilatation.copy())
            self.callback_params["tangents"].append(system.tangents.copy())
            
            
            return
            
#######################################
final_time = 5.0
#######################################

step_skip= int(final_time / dt/1000)

pp_list = defaultdict(list)
Rod.collect_diagnostics(shearable_rod).using(
    RodCallBack, step_skip=step_skip, callback_params=pp_list
)


Rod.finalize()
print("System finalized")

total_steps = int(final_time / dt)
print("Total steps to take", total_steps)

timestepper = PositionVerlet()
integrate(timestepper, Rod, final_time, total_steps, time_display_up = 60)

positions_over_time = np.array(pp_list["position"])
if (np.isnan(positions_over_time)==False).all()==False:
    print("Simulation diverge. Try lowering time step !")

with open("output/simulation_data.pickle", 'wb') as handle:
    pickle.dump(pp_list, handle, protocol=pickle.HIGHEST_PROTOCOL)

with open("output/simulation_data.pickle", 'rb') as handle:
    pp_list_read = pickle.load(handle)

solution_1 = process_solution_elastica(pp_list_read, step_skip=step_skip, base_length=base_length)


plot_multiple_solutions(
    solution_dfs=[solution_1],
    labels=["Linear Ribbon"],
    indices=np.linspace(0,solution_1.Index_solution.max(),100,dtype=np.int64),
    save_path="figure/run1.png" 
)

#
#solutions_auto = pd.read_csv('solution_auto_coupling.csv')

#temp1 = solutions_auto.Y
#temp2 = solutions_auto.Z 
#solutions_auto.Y = temp2
#solutions_auto.Z = temp1

#temp1 = solutions_auto.d1y
#temp2 = solutions_auto.d1z 
#solutions_auto.d1y = temp2
#solutions_auto.d1z = temp1

#temp1 = solutions_auto.d2y
#temp2 = solutions_auto.d2z 
#solutions_auto.d2y = temp2
#solutions_auto.d2z = temp1

#temp1 = solutions_auto.d3y
#temp2 = solutions_auto.d3z 
#solutions_auto.d3y = temp2
#solutions_auto.d3z = temp1

#solution_auto_unique = solutions_auto[solutions_auto.Index_solution == solutions_auto.Index_solution.max()]
#solution_elastica_unique = solution_1[solution_1.time==solution_1.time.max()]

#fig, ax = plt.subplots()

#solution_auto_unique.plot(x='s', y='X', ax=ax, label='Auto Solution', title = 'centerline X coordinate')
#solution_elastica_unique.plot(x='s', y='X', ax=ax, label='Elastica Solution')

#save_path="figure/X.png"
#fig.savefig(save_path, bbox_inches='tight', dpi=300)
#print(f"Plot saved to: {save_path}")
#plt.close(fig) 

#fig, ax = plt.subplots()

#solution_auto_unique.plot(x='s', y='Y', ax=ax, label='Auto Solution', title = 'centerline Y coordinate')
#solution_elastica_unique.plot(x='s', y='Y', ax=ax, label='Elastica Solution')

#save_path="figure/Y.png"
#fig.savefig(save_path, bbox_inches='tight', dpi=300)
#print(f"Plot saved to: {save_path}")
#plt.close(fig) 

#fig, ax = plt.subplots()

#solution_auto_unique.plot(x='s', y='Z', ax=ax, label='Auto Solution', title = 'centerline Z coordinate')
#solution_elastica_unique.plot(x='s', y='Z', ax=ax, label='Elastica Solution')

#save_path="figure/Z.png"
#fig.savefig(save_path, bbox_inches='tight', dpi=300)
#print(f"Plot saved to: {save_path}")
#plt.close(fig) 


#print("relative error on X at the tip:", 100*abs(solution_auto_unique.X.iloc[-1]-solution_elastica_unique.X.iloc[-1])/solution_auto_unique.X.iloc[-1],"%")
#print("relative error on Y at the tip:", 100*abs(solution_auto_unique.Y.iloc[-1]-solution_elastica_unique.Y.iloc[-1])/solution_auto_unique.Y.iloc[-1],"%")
#print("relative error on Z at the tip:", 100*abs(solution_auto_unique.Z.iloc[-1]-solution_elastica_unique.Z.iloc[-1])/solution_auto_unique.Z.iloc[-1],"%")