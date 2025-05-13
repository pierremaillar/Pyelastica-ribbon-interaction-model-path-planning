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
M = 10 # number of step in a path
min_length = 5

final_time = 0.1


paths_points, base_length = build_incremental_discretized_paths(path_plan, direction, normal, N_points, M, min_length)
plot_trajectory(paths_points[-1])



for idx, stat in enumerate(paths_points):


    d2_initial = np.array([np.ones(stat.shape[1]-1),np.zeros(stat.shape[1]-1),np.zeros(stat.shape[1]-1)])
    initial_position = stat 
    sleeve_position = initial_position
    d2_sleeve = d2_initial

    dt = 1e-7
    
    running = True
    while running:
        ribbon_output_list, sleeve_output_list, Ribbon_sim = create_env(initial_position, d2_initial, sleeve_position, d2_sleeve, base_length = base_length[idx], dt = dt, final_time = final_time)
        timestepper = PositionVerlet()
        total_steps = int(final_time / dt)
        nan_detected, _ = integrate(timestepper, Ribbon_sim, final_time, total_steps, time_display_up = 60)
        
        if nan_detected:
            dt /=1.2
            print("New dt:",dt)
        else:
            running = False
            
        if dt<=1e-9:
            running = False
            print("dt set below tolerance. Break simulation")

    with open(f"result/simulation_data_{idx}.pickle", 'wb') as handle:
        pickle.dump(ribbon_output_list, handle, protocol=pickle.HIGHEST_PROTOCOL)

    with open(f"result/simulation_data_sleeve_{idx}.pickle", 'wb') as handle:
        pickle.dump(sleeve_output_list, handle, protocol=pickle.HIGHEST_PROTOCOL)
