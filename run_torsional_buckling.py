import os
import pickle
import numpy as np
import pandas as pd

from collections import defaultdict
from elastica.src_plotting_ribbon import *
from create_sim_env import *



Y = 2.77e3
a = 5
t = 0.5
L = 50


force_adimentional = np.array([0.0, 0.15, 15])
force_newton = force_adimentional * (Y * a * t**3 / 12 / L**2)

print(force_newton)


direction = np.array([1,0,0])
normal = np.array([0,1,0])
final_time = 1.75
all_n_elem = np.array([5,10,50,100,150,200,300,500])
dt = 1e-6
base_length = L


for idx, n_elem in enumerate(all_n_elem):

    
    d1_initial = np.repeat(normal.reshape(3, 1), n_elem, axis=1)
    initial_position = np.array([np.linspace(0,1,n_elem+1)*base_length,np.zeros(n_elem+1),np.zeros(n_elem+1)])
    

    running = True
    while running:
        ribbon_output_list, Ribbon_sim = create_env_torsional_buckling(initial_position, d1_initial, dt = dt, 
                                                                       tip_force = force_newton, ramp_up_time = 0.1, nu = 3e1,
                                                                        base_length = base_length, thickness = t, width = a, 
                                                                        E = Y, poisson_ratio = 0.34, density = 1.017e-7, 
                                                                        final_time = final_time)
        timestepper = PositionVerlet()
        total_steps = int(final_time / dt)
        nan_detected, end_time = integrate(timestepper, Ribbon_sim, final_time, total_steps, time_display_up = 60)
        
        if nan_detected:
            dt /=1.2
            print("New dt:",dt)
        else:
            running = False
            
        if dt<=1e-9:
            running = False
            print("dt set below tolerance. Break simulation")

        with open(f"result/torsional_buckling_data_{idx}.pickle", 'wb') as handle:
            pickle.dump(ribbon_output_list, handle, protocol=pickle.HIGHEST_PROTOCOL)
    