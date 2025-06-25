import os
import pickle
import numpy as np
import pandas as pd

from collections import defaultdict
from elastica.src_plotting_ribbon import *
from create_sim_env import *


N_points = 50

direction = np.array([0,-1,0])
normal = np.array([0,0,1])
path_plan = pd.read_csv("path_0_0__15_-20__20_-40.csv")
#path_plan['twist (rad/mm)'] = [0, 4*np.pi/180, 0, 4*np.pi/180]
#path_plan = pd.read_csv("path_complex.csv")
#path_plan['twist (rad/mm)'] = [np.pi/20, 0, 0,0,0, -np.pi/20,0]

M = 50 # number of step in a path
min_length = 5
final_time = 0.1
kp=1.2


paths_points, base_length, normals = build_incremental_discretized_paths(path_plan, direction, normal, N_points, M, min_length)
idx = 0


for stat, normal in zip(paths_points, normals):
    
    d1_initial = normal
    initial_position = stat 
    sleeve_position = initial_position
    d1_sleeve = d1_initial

    dt = 1e-7
    
    running = True
    while running:
        ribbon_output_list, sleeve_output_list, Ribbon_sim = create_env(initial_position, d1_initial, sleeve_position, d1_sleeve, 
                                                                        base_length = base_length[idx], dt = dt, 
                                                                        final_time = final_time, k_p_input_force = kp)
        timestepper = PositionVerlet()
        total_steps = int(final_time / dt)
        nan_detected, end_time = integrate(timestepper, Ribbon_sim, final_time, total_steps, time_display_up = 60)
        
        if nan_detected:
            if end_time > 1e4*dt:
                kp /=1.5
                print("New kp", kp)
            dt /=1.2
            print("New dt:",dt,"\n")
        else:
            running = False
            
        if dt<=1e-9:
            running = False
            print("dt set below tolerance. Break simulation")

        with open(f"result/simulation_data_{idx}.pickle", 'wb') as handle:
            pickle.dump(ribbon_output_list, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
        with open(f"result/simulation_data_sleeve_{idx}.pickle", 'wb') as handle:
            pickle.dump(sleeve_output_list, handle, protocol=pickle.HIGHEST_PROTOCOL)

    print(f"Solution {idx} saved \n ------------ \n")
    idx+=1
