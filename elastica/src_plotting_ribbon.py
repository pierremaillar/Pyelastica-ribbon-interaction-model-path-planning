import numpy as np
import plotly.graph_objects as go
import pandas as pd
from collections import defaultdict
from matplotlib import pyplot as plt
import matplotlib.cm as cm
import math
import os

import matplotlib.animation as animation
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d import Axes3D





def plot_3D_ribbons_from_process_solution(solution_df, solution_indices=None, n_points=20, half_width=0.1, n_arrows = 10, save_path=None):
    """
    Plots 3D ribbon structures from a given partitioned solution file using Plotly.

    Parameters:
    -----------
    solution_df : pd.DataFrame
        A DataFrame containing the solution data with columns ['Index_solution', 'X', 'Y', 'Z', 
        'd3x', 'd3y', 'd3z', 'd2x', 'd2y', 'd2z', 's'].
    solution_indices : list of int, optional
        A list of solution indices to be visualized. If None, defaults to [1, 2].
    n_points : int, optional (default=20)
        Number of points across the ribbon width for surface plotting.
    n_arrows : int, optional (default=10)
        Number of arrows representing the cosserat fram allong the centerline to display
    half_width : float, optional (default=0.1)
        Half-width of the ribbon for visualization.
    save_path : str, optional (default="figure/3D_ribbons.png")
        File path to save the output image.

    Returns:
    --------
    fig : plotly.graph_objects.Figure
        A Plotly figure object displaying the 3D ribbon structures.

    Notes:
    ------
    - The function extracts centerline coordinates and constructs a ribbon surface using normal vectors.
    - Direction vectors (d1, d2, d3) are used to define the ribbon width and orientations.
    - Colored ribbons are plotted based on the 's' value in the DataFrame.
    - Arrows representing d1, d2, and d3 are plotted at intervals.
    - The figure is saved as an image and displayed.
    """
    if solution_indices is None:
        print("Don't forget to add the list of index you want to plot.")
        solution_indices = [1,2]

    fig = go.Figure()

    for index in solution_indices:
        one_solution = solution_df[solution_df['Index_solution'] == index]

        X_surf, Y_surf, Z_surf, colors = [], [], [], []

        min_index = min(solution_indices)
        max_index = max(solution_indices)
        if len(solution_indices) == 1:
            opacity = 1
        else: opacity = 0.2 + 0.8 * (index - min_index) / (max_index - min_index)

        for i in range(len(one_solution)):
            x, y, z = one_solution.iloc[i][['X', 'Y', 'Z']]
            
            d3 = one_solution.iloc[i][['d3x', 'd3y', 'd3z']].values
            d2 = one_solution.iloc[i][['d2x', 'd2y', 'd2z']].values
            d1 = np.cross(d3, d2)  # d1 is perpendicular to d2 and d3

            # Normalize vectors
            d3 /= np.linalg.norm(d3)
            d2 /= np.linalg.norm(d2)
            d1 /= np.linalg.norm(d1)

            # Create points along the width of the ribbon
            width_points = []
            for j in range(n_points):
                t = (j - n_points / 2) / (n_points / 2)
                width_point = np.array([x, y, z]) + t * half_width * d2
                width_points.append(width_point)

            width_points = np.array(width_points)
            X_surf.append(width_points[:, 0])
            Y_surf.append(width_points[:, 1])
            Z_surf.append(width_points[:, 2])

            color_value = one_solution.iloc[i]['s']
            colors.append([color_value] * n_points)

            if i % int(len(one_solution)/n_arrows) == 0:
                arrow_scale = 0.09
                fig.add_trace(go.Cone(x=[x], y=[y], z=[z], u=[-d1[0]], v=[-d1[1]], w=[-d1[2]], 
                                      colorscale='Blues', anchor='tail', sizemode='absolute',
                                      showscale=False, opacity=opacity, sizeref=arrow_scale, name='d1'))
                fig.add_trace(go.Cone(x=[x], y=[y], z=[z], u=[-d2[0]], v=[-d2[1]], w=[-d2[2]], 
                                      colorscale='Reds', anchor='tail', sizemode='absolute',
                                      showscale=False, opacity=opacity, sizeref=arrow_scale, name='d2'))
                fig.add_trace(go.Cone(x=[x], y=[y], z=[z], u=[d3[0]], v=[d3[1]], w=[d3[2]], 
                                      colorscale='Greens', anchor='tail', sizemode='absolute',
                                      showscale=False, opacity=opacity, sizeref=arrow_scale, name='d3'))

        X_surf, Y_surf, Z_surf, colors = map(np.array, (X_surf, Y_surf, Z_surf, colors))

        #print(f"Indice of the solution: {one_solution['Index_solution'].max()}")
        #print(f"X_max: {np.max(np.abs(one_solution.X)):.3f}")
        #print(f"Y_max: {np.max(np.abs(one_solution.Y)):.3f}")
        #print(f"Z_max: {np.max(np.abs(one_solution.Z)):.3f}\n")

        fig.add_trace(go.Surface(x=X_surf, y=Y_surf, z=Z_surf, surfacecolor=colors,
                                 colorscale='Plasma', showscale=False, opacity=opacity, name=f'Ribbon {index}'))

    fig.update_layout(
        scene=dict(
            xaxis_title='X',
            yaxis_title='Y',
            zaxis_title='Z',
            xaxis=dict(range=[-0.5, 1]),
            yaxis=dict(range=[-1, 1]),
            zaxis=dict(range=[-1, 1]),
            aspectmode='manual',
            aspectratio=dict(x=1, y=1, z=1)
        )
    )
    
    if save_path is not None:
        fig.savefig(save_path, bbox_inches='tight', dpi=300)
        print(f"Plot saved to: {save_path}")
        plt.close(fig) 
    else:
        return fig


def process_solution_file_auto(path_s):

    # Read the header to extract meta-information
    fl = pd.read_table(path_s, nrows=0, sep='\s+')
    lst = list(fl)
    con = [float(x) for x in lst]
    con = [int(x) for x in con]
    
    ntpl = con[6]  # Points in the time interval [0,1]
    nrowpr = con[8]  # Number of lines printed following the identifying line


    # Read the .s file, skipping bad lines
    sol_raw = pd.read_csv(path_s, sep='\s+', header=None, on_bad_lines='skip', skipinitialspace=True)

    # Search the index of the first line of solutions
    index_no_nan = sol_raw.dropna().index + 1  
    index_no_nan = index_no_nan[sol_raw.iloc[index_no_nan - 1, :].sum(axis=1) > 50]

    # Identify the rows with exactly 1 non-NaN value, likely separator rows
    index_one_non_nan = sol_raw[sol_raw.notna().sum(axis=1) <= 2].index

    # Find sequences of consecutive rows with single non-NaN values
    consecutive_indices = [index_one_non_nan[i] for i in range(len(index_one_non_nan) - 2)
                           if index_one_non_nan[i] == index_one_non_nan[i + 1] - 1
                           and index_one_non_nan[i] == index_one_non_nan[i + 2] - 2]

    consecutive_indices = np.array(consecutive_indices)
    index_no_nan = consecutive_indices - 4*ntpl + 1
    extracted_solutions = []

    # Extract solutions between corresponding pairs of index_no_nan and consecutive_indices
    for i, j in zip(index_no_nan, consecutive_indices):
        extracted_solutions.append(sol_raw.iloc[i:j + 1])

    # Concatenate extracted solutions
    extracted_solutions_df = pd.concat(extracted_solutions)

    # Clean and limit columns in the extracted solutions
    extracted_solutions_cleaned_df = extracted_solutions_df.iloc[:, :7]

    # Flatten and remove NaN values
    flattened = extracted_solutions_cleaned_df.values.flatten()
    non_nan_values = flattened[~np.isnan(flattened)]

    # Define the number of columns (22)
    num_columns = 22

    # Reshape into a DataFrame
    solution_df = pd.DataFrame(non_nan_values.reshape(-1, num_columns))

    # Set column names
    solution_df.columns = ['Cont_par', 'X', 'Y', 'Z', 'd3x', 'd3y', 'd3z', 'd2x', 'd2y', 'd2z', 'd1x', 'd1y', 'd1z',
                           'R1', 'R2', 'R3', 'm1', 'm2', 'm3', 'k2', 'k3', 's']

    # Calculate the number of points per solution
    n = (consecutive_indices[0] - index_no_nan[0] + 1) // 4

    # Determine number of rows
    num_rows = solution_df.shape[0]

    # Create an 'Index_solution' column to label each solution
    index_column = np.repeat(np.arange(1, (num_rows // n) + 2), n)[:num_rows]
    solution_df['Index_solution'] = index_column

    print(f"\nPartition found {solution_df['Index_solution'].max()} solutions with {n} points each")

    return solution_df.astype('float64')



def plot_multiple_solutions(
    solution_dfs, 
    labels, 
    indices, 
    start_color_idx=0, 
    true_solution=None, 
    print_legend=False,
    variables=[
        ('X', 'X Position'),
        ('Y', 'Y Position'),
        ('Z', 'Z Position'),
        ('R1', 'R1'),
        ('R2', 'R2'),
        ('R3', 'R3'),
        ('VX', 'VX'),
        ('VY', 'VY'),
        ('VZ', 'VZ'),
    ],
    save_path=None
):
    """
    Plot multiple ribbon simulation solutions as subplots.

    This function generates subplots for each specified variable (`X`, `Y`, `Z`, `R1`, etc.)
    for one or more solution datasets. Each variable is plotted as a function of arc-length `s`
    for a set of selected simulation frames (`Index_solution` values). The plots can also include
    a reference "true" solution for visual comparison.

    Parameters
    ----------
    solution_dfs : list of pandas.DataFrame
        A list of DataFrames, each representing a full simulation solution. Each DataFrame should
        contain columns such as: 'Index_solution', 's', 'X', 'Y', 'Z', 'R1', 'R2', 'R3', 'V1', 'V2', 'V3'.

    labels : list of str
        Labels for each solution in `solution_dfs`, used in the plot legend.

    indices : list or array-like of int
        Frame indices (`Index_solution`) to plot from each solution DataFrame.

    start_color_idx : int, optional (default=0)
        Index to shift the base colormap for the solutions. Useful for differentiating plot sets
        when calling this function multiple times.

    true_solution : pandas.DataFrame, optional
        A reference solution DataFrame with the same variables. If provided, it is plotted as a
        dashed red line for comparison.

    print_legend : bool, optional (default=False)
        Whether to include a legend in the plot. Only shown in the top-right subplot.

    variables : list of tuple, optional
        A list of (column_name, display_name) pairs representing the variables to plot and their
        corresponding y-axis labels.

    save_path : str or None, optional
        If provided, saves the resulting plot to the specified file path. The file format is
        inferred from the extension (e.g., .png, .pdf). If None, the plot is displayed instead.

    Returns
    -------
    None
        Displays the plot or saves it to a file.
    """

    num_solutions = len(solution_dfs)

    num_vars = len(variables)
    ncols = math.ceil(num_vars/3)
    nrows = 3
    fig_width = ncols * 4
    fig_height = nrows * 3

    fig, axs = plt.subplots(nrows, ncols, figsize=(fig_width, fig_height))

    colormaps = [cm.Blues, cm.Oranges, cm.Greens, cm.Purples, cm.Reds, cm.Greys]
    colormap = colormaps[start_color_idx % len(colormaps)]
    colors = [colormaps[i % len(colormaps)](np.linspace(0.3, 1, len(indices))) for i in range(num_solutions)]

    for i, (var, title) in enumerate(variables):
        row, col = divmod(i, 3)
        axs[row, col].grid(True)

        for j, (df, label) in enumerate(zip(solution_dfs, labels)):
            for k, index_solution in enumerate(indices):
                selected_df = df[df['Index_solution'] == index_solution].reset_index()
                color = colors[j][k]
                axs[row, col].plot(
                    selected_df['s'], selected_df[var],
                    color=color, linestyle='-', label=label if k == 0 else ""
                )

        if true_solution is not None:
            axs[row, col].plot(
                true_solution['s'], true_solution[var],
                color='red', linestyle='--', linewidth=2, label='True Solution'
            )

        axs[row, col].set_ylabel(title, fontsize=12)
        axs[row, col].set_xlabel('s', fontsize=12)
        axs[row, col].set_title(title, fontsize=14)
        axs[row, col].tick_params(axis='both', labelsize=10)

    if print_legend:
        axs[0, -1].legend(title='Legend', fontsize=10, loc='upper left')

    plt.tight_layout()
    plt.subplots_adjust(top=0.90)

    # Save or show
    if save_path is not None:
        fig.savefig(save_path, bbox_inches='tight', dpi=300)
        print(f"Plot saved to: {save_path}")
        plt.close(fig) 
    else:
        plt.show()



def plot_multiple_solutions_sleeve(
    solution_dfs, 
    labels, 
    indices, 
    start_color_idx=0, 
    print_legend=False,
    variables=[
        ('displacement_X', 'Displacement on X'),
        ('displacement_Y', 'Displacement on Y'),
        ('displacement_Z', 'Displacement on Z'),
        ('response_force_X', 'Force on X'),
        ('response_force_Y', 'Force on Y'),
        ('response_force_Z', 'Force on Z'),
        ('rotation_angle_d1', 'Rotation angle d1'),
        ('rotation_angle_d2', 'Rotation angle d2'),
        ('rotation_angle_d3', 'Rotation angle d3'),
        ('response_couple_X', 'Couple on X'),
        ('response_couple_Y', 'Couple on Y'),
        ('response_couple_Z', 'Couple on Z'),
    ],
    save_path=None
):


    num_solutions = len(solution_dfs)

    num_vars = len(variables)

    ncols = 3
    nrows = math.ceil(num_vars/3)
    fig_width = ncols * 4
    fig_height = nrows * 3

    fig, axs = plt.subplots(nrows, ncols, figsize=(fig_width, fig_height))

    colormaps = [cm.Blues, cm.Oranges, cm.Greens, cm.Purples, cm.Reds, cm.Greys]
    colormap = colormaps[start_color_idx % len(colormaps)]
    colors = [colormaps[i % len(colormaps)](np.linspace(0.3, 1, len(indices))) for i in range(num_solutions)]

    for i, (var, title) in enumerate(variables):
        row, col = divmod(i, 3)
        axs[row, col].grid(True)

        for j, (df, label) in enumerate(zip(solution_dfs, labels)):
            for k, index_solution in enumerate(indices):
                selected_df = df[df['Index_solution'] == index_solution].reset_index()
                color = colors[j][k]
                axs[row, col].plot(
                    selected_df['s'], selected_df[var],
                    color=color, linestyle='-', label=label if k == 0 else ""
                )

        axs[row, col].set_ylabel(title, fontsize=12)
        axs[row, col].set_xlabel('s', fontsize=12)
        axs[row, col].set_title(title, fontsize=14)
        axs[row, col].tick_params(axis='both', labelsize=10)

    if print_legend:
        axs[0, -1].legend(title='Legend', fontsize=10, loc='upper left')

    plt.tight_layout()
    plt.subplots_adjust(top=0.90)

    # Save or show
    if save_path is not None:
        fig.savefig(save_path, bbox_inches='tight', dpi=300)
        print(f"Plot saved to: {save_path}")
        plt.close(fig) 
    else:
        plt.show()


def plot_3D_ribbons_from_process_solutions(solution_df1, solution_indices1, solution_df2, solution_indices2,  color_df1 = 'viridis', color_df2 = 'Plasma',
                                           n_points=20, half_width=0.1, n_arrows=10, save_path="figure/3D_ribbons.png"):
    """
    Plots 3D ribbon structures from two different partitioned solution files using Plotly.

    Parameters:
    -----------
    solution_df1 : pd.DataFrame
        First DataFrame containing the solution data with columns ['Index_solution', 'X', 'Y', 'Z', 
        'd3x', 'd3y', 'd3z', 'd2x', 'd2y', 'd2z', 's'].
    
    solution_indices1 : list of int
        List of solution indices to be visualized from solution_df1.

    solution_df2 : pd.DataFrame
        Second DataFrame containing another set of solution data.

    solution_indices2 : list of int
        List of solution indices to be visualized from solution_df2.

    n_points : int, optional (default=20)
        Number of points across the ribbon width for surface plotting.

    n_arrows : int, optional (default=10)
        Number of arrows representing the Cosserat frame along the centerline to display.

    half_width : float, optional (default=0.1)
        Half-width of the ribbon for visualization.

    save_path : str, optional (default="figure/3D_ribbons.png")
        File path to save the output image.

    Returns:
    --------
    fig : plotly.graph_objects.Figure
        A Plotly figure object displaying the 3D ribbon structures.

    Notes:
    ------
    - Two datasets are plotted with separate colors.
    - Opacity varies within each group.
    - Cosserat frame arrows are included.
    """
    fig = go.Figure()

    # Define colormaps and opacity scaling
    datasets = [(solution_df1, solution_indices1, color_df1), (solution_df2, solution_indices2, color_df2)]
    
    for solution_df, solution_indices, colormap in datasets:
        min_index = min(solution_indices)
        max_index = max(solution_indices)

        for index in solution_indices:
            one_solution = solution_df[solution_df['Index_solution'] == index]

            X_surf, Y_surf, Z_surf, colors = [], [], [], []
            
            if len(solution_indices) == 1:
                opacity = 1
            else: opacity = 0.2 + 0.8 * (index - min_index) / (max_index - min_index)

            for i in range(len(one_solution)):
                x, y, z = one_solution.iloc[i][['X', 'Y', 'Z']]
                
                d3 = one_solution.iloc[i][['d3x', 'd3y', 'd3z']].values
                d2 = one_solution.iloc[i][['d2x', 'd2y', 'd2z']].values
                d1 = np.cross(d3, d2)  # d1 is perpendicular to d2 and d3

                # Normalize vectors
                d3 /= np.linalg.norm(d3)
                d2 /= np.linalg.norm(d2)
                d1 /= np.linalg.norm(d1)

                # Create points along the width of the ribbon
                width_points = []
                for j in range(n_points):
                    t = (j - n_points / 2) / (n_points / 2)
                    width_point = np.array([x, y, z]) + t * half_width * d2
                    width_points.append(width_point)

                width_points = np.array(width_points)
                X_surf.append(width_points[:, 0])
                Y_surf.append(width_points[:, 1])
                Z_surf.append(width_points[:, 2])

                color_value = one_solution.iloc[i]['s']
                colors.append([color_value] * n_points)

                if i % int(len(one_solution) / n_arrows) == 0:
                    arrow_scale = 0.09
                    fig.add_trace(go.Cone(x=[x], y=[y], z=[z], u=[-d1[0]], v=[-d1[1]], w=[-d1[2]], 
                                          colorscale='Blues', anchor='tail', sizemode='absolute',
                                          showscale=False, opacity=opacity, sizeref=arrow_scale, name=f'd1 ({index})'))
                    fig.add_trace(go.Cone(x=[x], y=[y], z=[z], u=[-d2[0]], v=[-d2[1]], w=[-d2[2]], 
                                          colorscale='Reds', anchor='tail', sizemode='absolute',
                                          showscale=False, opacity=opacity, sizeref=arrow_scale, name=f'd2 ({index})'))
                    fig.add_trace(go.Cone(x=[x], y=[y], z=[z], u=[d3[0]], v=[d3[1]], w=[d3[2]], 
                                          colorscale='Greens', anchor='tail', sizemode='absolute',
                                          showscale=False, opacity=opacity, sizeref=arrow_scale, name=f'd3 ({index})'))

            X_surf, Y_surf, Z_surf, colors = map(np.array, (X_surf, Y_surf, Z_surf, colors))

            #print(f"Solution {index} from dataset ({colormap}):")
            #print(f"  X_max: {np.max(np.abs(one_solution.X)):.3f}")
            #print(f"  Y_max: {np.max(np.abs(one_solution.Y)):.3f}")
            #print(f"  Z_max: {np.max(np.abs(one_solution.Z)):.3f}\n")

            fig.add_trace(go.Surface(
                x=X_surf, y=Y_surf, z=Z_surf, surfacecolor=colors,
                colorscale=colormap, showscale=False,
                opacity=opacity, name=f'Ribbon {index}'
            ))

    fig.update_layout(
        #title="3D Ribbon Comparison from Two Datasets",
        scene=dict(
            xaxis_title='X',
            yaxis_title='Y',
            zaxis_title='Z',
            xaxis=dict(range=[-0.5, 1]),
            yaxis=dict(range=[-1, 1]),
            zaxis=dict(range=[-1, 1]),
            aspectmode='manual',
            aspectratio=dict(x=1, y=1, z=1)
        )
    )
    
    #fig.write_image(save_path)
    fig.show()
    return fig

def process_solution_elastica(pp_list_read, base_length):
    """
    Processes solution data and converts it into a structured DataFrame.

    Parameters:
    -----------
    pp_list_read : dict
        DictionaR2 containing solution data with keys: 'time', 'step', 'position', 
        'directors', 'internal_stress', 'internal_couple', and 'curvature'.
        
    base_length : float
        Reference length to normalize positional values.

    Returns:
    --------
    pandas.DataFrame
        Processed DataFrame containing the solution data.
    """
    rows = []
    j = 0

    for t, step, pos, director, stress, couple, forces_in, forces_ex, curvature, strains, dilatation, tangents, velocities in zip(
        pp_list_read["time"], pp_list_read["step"],
        pp_list_read["position"], pp_list_read["directors"],
        pp_list_read["internal_stress"], pp_list_read["internal_couple"], pp_list_read["internal_forces"], pp_list_read["external_forces"],
        pp_list_read["curvature"], pp_list_read["sigma"], pp_list_read["dilatation"], pp_list_read["tangents"], pp_list_read["velocity"]
    ):
        num_elements = pos.shape[1]  # Number of points
        # Extend director matrix
        last_element = director[:, :, -1][:, :, np.newaxis] 
        director_extended = np.concatenate((director, last_element), axis=2)

        # Extend tangents matrix
        last_element = tangents[:, -1][:, np.newaxis] 
        tangents_extended = np.concatenate((tangents, last_element), axis=1)

        # Extend couple
        last_element = couple[:, -1][:, np.newaxis] 
        couple_extended = np.concatenate((couple, last_element), axis=1)
        couple_extended = np.concatenate((couple_extended, last_element), axis=1)

        # Extend stress
        last_element = stress[:, -1][:, np.newaxis] 
        stress_extended = np.concatenate((stress, last_element), axis=1)

        # Extend forces_int
        last_element = forces_in[:, -1][:, np.newaxis] 
        forces_in_extended = np.concatenate((forces_in, last_element), axis=1)

        # Extend forces_ex
        last_element = forces_ex[:, -1][:, np.newaxis] 
        forces_ex_extended = np.concatenate((forces_ex, last_element), axis=1)

        # Extend curvature
        last_element = curvature[:, -1][:, np.newaxis] 
        curvature_extended = np.concatenate((curvature, last_element), axis=1)
        curvature_extended = np.concatenate((curvature_extended, last_element), axis=1)

        # Extend strains
        last_element = strains[:, -1][:, np.newaxis] 
        strains_extended = np.concatenate((strains, last_element), axis=1)   

        # Extend dilatation
        last_element = dilatation[-1]
        dilatation_extended = np.append(dilatation,last_element)
       

        # Populate rows
        for i in range(num_elements):
            rows.append({
                "Index_solution": j,
                "time": t,
                "step": step,
                "s": i / (num_elements - 1),
                "X": pos[0, i] / base_length,
                "Y": pos[1, i] / base_length,
                "Z": pos[2, i] / base_length,
                'd3x': director_extended[2, 0, i], 
                'd3y': director_extended[2, 1, i], 
                'd3z': director_extended[2, 2, i],
                'd2x': director_extended[1, 0, i], 
                'd2y': director_extended[1, 1, i], 
                'd2z': director_extended[1, 2, i],
                'R1': stress_extended[0, i],
                'R2': stress_extended[1, i],
                'R3': stress_extended[2, i],
                'RX': forces_in_extended[0, i],
                'RY': forces_in_extended[1, i],
                'RZ': forces_in_extended[2, i],
                'ReX': forces_ex_extended[0, i],
                'ReY': forces_ex_extended[1, i],
                'ReZ': forces_ex_extended[2, i],
                'mX': couple_extended[2, i],
                'mY': couple_extended[2, i],
                'mZ': couple_extended[2, i],
                'k1': curvature_extended[2, i],
                'k2': curvature_extended[2, i],
                'k3': curvature_extended[2, i],
                'e1': strains_extended[0,i],
                'e2': strains_extended[1,i],
                'e3': strains_extended[2,i],
                'dilatation': dilatation_extended[i],
                'tx': tangents_extended[0, i], 
                'ty': tangents_extended[1, i], 
                'tz': tangents_extended[2, i],
                'VX' :velocities[0,i],
                'VY' :velocities[1,i],
                'VZ' :velocities[2,i]
            })
        j+=1

    return pd.DataFrame(rows)




def process_solution_elastica_sleeve(pp_list_read):
    """
    Processes solution data and converts it into a structured DataFrame.

    Parameters:
    -----------
    pp_list_read : dict
        Dictionary containing solution data with keys: 'response_force_sleeve', 'displacement_sleeve', 'response_couple_sleeve', 
        'rotation_sleeve'
    
    step_skip : int
        Step increment to normalize the solution index.

    Returns:
    --------
    pandas.DataFrame
        Processed DataFrame containing the solution data.
    """
    rows = []
    j = 0

    for t, step, force, displacement, couple, rotation in zip(
        pp_list_read["time"], pp_list_read["step"],
        pp_list_read["response_force_sleeve"], pp_list_read["displacement_sleeve"],
        pp_list_read["response_couple_sleeve"], pp_list_read["rotation_sleeve"]
    ):

        num_elements = force.shape[1]
        
        # Populate rows
        for i in range(num_elements):
            rows.append({
                "Index_solution": j,
                "time": t,
                "step": step,
                "s": i / (num_elements - 1),
                "response_force_X": force[0, i],
                "response_force_Y": force[1, i],
                "response_force_Z": force[2, i],
                "displacement_X": displacement[0, i],
                "displacement_Y": displacement[1, i],
                "displacement_Z": displacement[2, i],
                "response_couple_X": couple[0, i],
                "response_couple_Y": couple[1, i],
                "response_couple_Z": couple[2, i],
                "rotation_angle_d1": np.arcsin(rotation[0, i]),
                "rotation_angle_d2": np.arcsin(rotation[1, i]),
                "rotation_angle_d3": np.arcsin(rotation[2, i]),
            })
        j+=1

    return pd.DataFrame(rows)



def sanity_check_plot(solution):

    # Extra fields
    solution["dilatation_error"] = solution["dilatation"] - 1
    a_dot_b = solution['d3x'] * solution['tx'] + solution['d3y'] * solution['ty'] + solution['d3z'] * solution['tz']
    norm_a = np.sqrt(solution['d3x']**2 + solution['d3y']**2 + solution['d3z']**2)
    norm_b = np.sqrt(solution['tx']**2 + solution['ty']**2 + solution['tz']**2)
    solution["sin_theta_txd3"] = a_dot_b / (norm_a * norm_b) - 1

    solution["norm_V"] = np.sqrt((solution['VX'])**2 + solution['VY']**2 + solution['VZ']**2)

    solution_mean = solution.groupby("time").mean()
    solution_max = solution.abs().groupby("time").max()
    solution_l2 = solution.groupby("time").apply(lambda df: np.sqrt((df**2).sum()))

    # Additional processing for length calculation
    solution[['dX', 'dY', 'dZ']] = solution.groupby('time')[['X', 'Y', 'Z']].diff()
    solution['ds'] = np.sqrt(solution['dX']**2 + solution['dY']**2 + solution['dZ']**2)
    length_by_time = abs(solution.groupby('time')['ds'].sum() - 1).reset_index(name='curve_length_minus1')

    # Prepare selections
    one_solution_diff = solution[solution['s'] == 0.0].copy()
    one_solution_final = solution[solution['step'] == solution['step'].max()].copy()

    # Prepare figure with subplots
    fig, axes = plt.subplots(4, 2, figsize=(18, 22))
    fig.suptitle("Sanity Check Simulation", fontsize=22, fontweight='bold')

    fontdict = {'fontsize': 14}

    # Subplot 1: Steady state
    axes[0, 0].plot(solution_max.index, (solution_l2['norm_V']))
    axes[0, 0].set_title("Steady State", **fontdict)
    axes[0, 0].set_xlabel("Time", fontsize=12)
    axes[0, 0].set_ylabel("l2(|V|)", fontsize=12)
    axes[0, 0].set_yscale('log')
    axes[0, 0].grid(True)

    # Subplot 2: Time Step
    one_solution_diff.time.diff().reset_index(drop=True).plot(ax=axes[0, 1])
    axes[0, 1].set_title("Time Step", **fontdict)
    axes[0, 1].set_xlabel("Index", fontsize=12)
    axes[0, 1].set_ylabel("dt", fontsize=12)
    axes[0, 1].grid(True)

    # Subplot 3: Dilatation vs s
    axes[1, 0].plot(one_solution_final['s'], one_solution_final['dilatation_error'])
    axes[1, 0].set_title("Dilatation vs s", **fontdict)
    axes[1, 0].set_xlim(0, 1)
    axes[1, 0].set_xlabel("s", fontsize=12)
    axes[1, 0].set_ylabel("Dilatation Error", fontsize=12)
    axes[1, 0].grid(True)

    # Subplot 4: Shear Strain Error
    e1 = solution_l2.e1
    e2 = solution_l2.e2
    e3 = solution_l2.e3
    (np.sqrt(e1**2 + e2**2 + e3**2) / np.sqrt(solution_l2.X**2 + solution_l2.Y**2 + solution_l2.Z**2)).plot(ax=axes[1, 1], logy=True)
    axes[1, 1].set_title("|Shear Strain| over |r| (L2)", **fontdict)
    axes[1, 1].set_xlabel("Time", fontsize=12)
    axes[1, 1].set_ylabel("Error", fontsize=12)
    axes[1, 1].grid(True, which="both")

    # Subplot 5: Dilatation L2
    solution_l2.dilatation_error.iloc[2:].plot(ax=axes[2, 0], logy=True)
    axes[2, 0].set_title("Dilatation (L2)", **fontdict)
    axes[2, 0].set_xlabel("Time", fontsize=12)
    axes[2, 0].set_ylabel("L2 Norm", fontsize=12)
    axes[2, 0].grid(True, which="both")

    # Subplot 6: sin(theta) tx.d3
    solution_l2.sin_theta_txd3.plot(ax=axes[2, 1], logy=True)
    axes[2, 1].set_title("sin(theta) tx.d3", **fontdict)
    axes[2, 1].set_xlabel("Time", fontsize=12)
    axes[2, 1].set_ylabel("Value", fontsize=12)
    axes[2, 1].grid(True, which="both")

    # Subplot 7: Curve Length -1
    axes[3, 0].plot(length_by_time['time'], length_by_time['curve_length_minus1'])
    axes[3, 0].set_title("Centerline Length -1", **fontdict)
    axes[3, 0].set_xlabel("Time", fontsize=12)
    axes[3, 0].set_ylabel("|Length - 1|", fontsize=12)
    axes[3, 0].set_yscale("log")
    axes[3, 0].grid(True, which="both")

    # Leave subplot (3,1) empty
    axes[3, 1].axis('off')

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()

    print(f"Final stress stat at s = 0:\nR1: {one_solution_final.R1.iloc[0]:.4e}\nR2: {one_solution_final.R2.iloc[0]:.4e}\nR3: {one_solution_final.R3.iloc[0]:.4e}")

    print(f"Final stress stat at s = 1:\nR1: {one_solution_final.R1.iloc[-1]:.4e}\nR2: {one_solution_final.R2.iloc[-1]:.4e}\nR3: {one_solution_final.R3.iloc[-1]:.4e}")


def control_law_plot(solution):
    solution_at_base = solution[solution.s == 0]
    solution_at_tip = solution[solution.s == 1]

    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle("BC check and control law check", fontsize=22, fontweight='bold')

    fontdict = {'fontsize': 14}

    # Subplot 1: X

    axes[0, 0].plot(solution_at_tip.time, solution_at_tip.R3, label="R3 (internal, tip)")
    axes[0, 0].plot(solution_at_base.time, solution_at_base.ReX, label="Rex (external, base)")
    axes[0, 0].set_title("Forces in X-direction", **fontdict)
    axes[0, 0].set_xlabel("Time [s]", fontsize=12)
    axes[0, 0].set_ylabel("Force [N]", fontsize=12)
    axes[0, 0].grid(True)
    axes[0, 0].legend()

    # Subplot 2: Y
    
    axes[0, 1].plot(solution_at_tip.time, solution_at_tip.R3, label="R3 (internal, tip)")
    axes[0, 1].plot(solution_at_base.time, solution_at_base.ReY, label="ReY (external, base)")
    axes[0, 1].set_title("Forces in Y-direction", **fontdict)
    axes[0, 1].set_xlabel("Time [s]", fontsize=12)
    axes[0, 1].set_ylabel("Force [N]", fontsize=12)
    axes[0, 1].grid(True)
    axes[0, 1].legend()

    # Subplot 3: Z

    axes[0, 2].plot(solution_at_tip.time, solution_at_tip.R3, label="R3 (internal, tip)")
    axes[0, 2].plot(solution_at_base.time, solution_at_base.ReZ, label="ReZ (external, base)")
    axes[0, 2].set_title("Forces in Z-direction", **fontdict)
    axes[0, 2].set_xlabel("Time [s]", fontsize=12)
    axes[0, 2].set_ylabel("Force [N]", fontsize=12)
    axes[0, 2].grid(True)
    axes[0, 2].legend()

    # Subplot 4: Internal forces at base
    axes[1, 0].plot(solution_at_base.time, solution_at_base.RX, label="RX (internal, base)")
    axes[1, 0].plot(solution_at_base.time, solution_at_base.RY, label="RY (internal, base)")
    axes[1, 0].plot(solution_at_base.time, solution_at_base.RZ, label="RZ (internal, base)")
    axes[1, 0].set_title("Internal Forces at the Base", **fontdict)
    axes[1, 0].set_xlabel("Time [s]", fontsize=12)
    axes[1, 0].set_ylabel("Force [N]", fontsize=12)
    axes[1, 0].grid(True)
    axes[1, 0].legend()

    # Subplot 5: Internal forces at tip
    axes[1, 1].plot(solution_at_tip.time, solution_at_tip.RX, label="RX (internal, tip)")
    axes[1, 1].plot(solution_at_tip.time, solution_at_tip.RY, label="RY (internal, tip)")
    axes[1, 1].plot(solution_at_tip.time, solution_at_tip.RZ, label="RZ (internal, tip)")
    axes[1, 1].set_title("Internal Forces at the Tip", **fontdict)
    axes[1, 1].set_xlabel("Time [s]", fontsize=12)
    axes[1, 1].set_ylabel("Force [N]", fontsize=12)
    axes[1, 1].grid(True)
    axes[1, 1].legend()

    # Subplot 6: Leave empty
    axes[1, 2].axis('off')

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()
    

def build_incremental_discretized_paths(path_df, direction, normal, n_points, M, min_length):
    """
    Build discretized paths with twist support.
    
    Parameters:
    -----------
    path_df : DataFrame
        DataFrame with columns 'type', 'length (mm)', 'curvature (1/mm)', and optionally 'twist (rad/mm)'
    direction : ndarray
        Initial direction vector
    normal : ndarray
        Initial normal vector
    n_points : int
        Number of points in the final path
    M : int
        Number of paths to generate with different lengths
    min_length : float
        Minimum path length
        
    Returns:
    --------
    all_paths : list of arrays
        List of path arrays, each with shape (3, n_points)
    lengths : ndarray
        Array of path lengths
    all_normals : list of arrays
        List of normal vectors along each path segment, each with shape (3, n_points-1)
    """
    direction = direction / np.linalg.norm(direction)
    normal = normal - direction * np.dot(direction, normal)
    normal /= np.linalg.norm(normal)
    binormal = np.cross(direction, normal)
    
    total_length = path_df['length (mm)'].sum()
    max_length = total_length
    lengths = np.linspace(min_length, max_length, M)
    
    all_paths = []
    all_normals = []  # Store normals between consecutive points
    
    for target_length in lengths:
        path_points = [np.zeros(3)]
        segment_normals = []  # Will store normals for segments between points
        
        pos = np.zeros(3)
        T = direction.copy()
        N = normal.copy()
        B = binormal.copy()
        
        accumulated_length = 0.0
        
        for idx, row in path_df.iterrows():
            seg_type = row['type']
            L = row['length (mm)']
            k = row['curvature (1/mm)']
            # Get twist rate if it exists, otherwise default to 0
            twist_rate = row.get('twist (rad/mm)', 0.0)
            
            remaining_length = target_length - accumulated_length
            if remaining_length < 0:
                break
                
            # Clip segment if needed
            seg_length = min(L, remaining_length)
            if seg_length < 1e-6:
                continue
                
            # Allocate points based on segment proportion
            # Important: for large n_points, this ensures proper distribution
            seg_points = max(2, int(np.ceil((seg_length / target_length) * n_points)))
            ds = seg_length / (seg_points - 1)
            
            if seg_type == 'line' or np.isclose(k, 0.0):
                # Handle straight line segment
                for i in range(1, seg_points):
                    # Apply twist around tangent
                    if not np.isclose(twist_rate, 0.0):
                        twist_angle = twist_rate * ds
                        twist_matrix = skew_rotation(T, twist_angle)
                        N = twist_matrix @ N
                        B = twist_matrix @ B
                        
                        # Ensure orthonormality
                        B = np.cross(T, N)
                        B /= np.linalg.norm(B)
                        N = np.cross(B, T)
                        N /= np.linalg.norm(N)
                    
                    # Store the normal for this segment before moving
                    segment_normals.append(N.copy())
                    
                    # Move along tangent
                    pos = pos + T * ds
                    path_points.append(pos.copy())
            else:  # arc
                R = 1 / k
                theta_total = seg_length * k
                dtheta = theta_total / (seg_points - 1)
                center = pos + N * R
                
                for i in range(1, seg_points):
                    # Store the normal for this segment before moving
                    segment_normals.append(N.copy())
                    
                    # Apply curvature
                    theta = dtheta
                    rot_matrix = skew_rotation(B, theta)
                    
                    # Update position and frame for curvature
                    pos = center + rot_matrix @ (pos - center)
                    T = rot_matrix @ T
                    N = rot_matrix @ N
                    B = rot_matrix @ B
                    
                    # Apply twist if needed
                    if not np.isclose(twist_rate, 0.0):
                        twist_angle = twist_rate * ds
                        twist_matrix = skew_rotation(T, twist_angle)
                        N = twist_matrix @ N
                        B = twist_matrix @ B
                    
                    # Ensure orthonormality
                    T /= np.linalg.norm(T)
                    B = np.cross(T, N)
                    B /= np.linalg.norm(B)
                    N = np.cross(B, T)
                    N /= np.linalg.norm(N)
                    
                    # Store current position
                    path_points.append(pos.copy())
                    
                    # Move center for next iteration to maintain proper arc shape
                    center = pos + N * R
            
            accumulated_length += seg_length
            if accumulated_length >= target_length:
                break
        
        # Interpolate to ensure exactly `n_points` points
        path_points = np.array(path_points)
        segment_normals = np.array(segment_normals)
        
        if len(path_points) > 1:
            # Calculate distances between consecutive points
            distances = np.linalg.norm(np.diff(path_points, axis=0), axis=1)
            arc_lengths = np.insert(np.cumsum(distances), 0, 0.0)
            target_arc = np.linspace(0, arc_lengths[-1], n_points)
            
            # Interpolate path points
            interp_path = np.zeros((n_points, 3))
            for j in range(3):
                interp_path[:, j] = np.interp(target_arc, arc_lengths, path_points[:, j])
            
            # For the segment normals, we need n_points-1 normals
            # Calculate the centers of segments in the interpolated path
            segment_centers = (target_arc[:-1] + target_arc[1:]) / 2
            
            # Interpolate normals at segment centers
            interp_normals = np.zeros((n_points-1, 3))
            
            if len(segment_normals) > 0:
                # Calculate the centers of original segments for interpolation reference
                original_segment_centers = (arc_lengths[:-1] + arc_lengths[1:]) / 2
                
                for j in range(3):
                    interp_normals[:, j] = np.interp(segment_centers, original_segment_centers, segment_normals[:, j])
                
                # Normalize interpolated normals
                norms = np.linalg.norm(interp_normals, axis=1, keepdims=True)
                interp_normals = interp_normals / norms
            else:
                # Handle edge case with no segments
                interp_normals = np.tile(normal, (n_points-1, 1))
            
            all_paths.append(interp_path.T)  # shape (3, n_points)
            all_normals.append(interp_normals.T)  # shape (3, n_points-1)
        else:
            # Handle edge case with a single point (no segments)
            all_paths.append(np.zeros((3, n_points)))
            all_normals.append(np.tile(normal.reshape(3, 1), n_points-1))  # Empty normal array
    
    return all_paths, lengths, all_normals

def skew(v):
    """Skew-symmetric matrix for cross product"""
    return np.array([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0]
    ])

def skew_rotation(axis, angle):
    """
    Create rotation matrix from axis and angle using Rodrigues' formula
    """
    axis = axis / np.linalg.norm(axis)
    return (
        np.cos(angle) * np.eye(3) +
        np.sin(angle) * skew(axis) +
        (1 - np.cos(angle)) * np.outer(axis, axis)
    )
def plot_trajectory(path_points, show=True, color='b', label='Trajectory'):
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')
    
    path_points = np.array(path_points).T
    ax.scatter(path_points[:, 0], path_points[:, 1], path_points[:, 2], color=color, s = 1, label=label)

    # Mark start and end
    ax.scatter(*path_points[0], color='green', s=20, label='Start')
    ax.scatter(*path_points[-1], color='red', s=20, label='End')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    limit = np.max(abs(path_points))
    
    ax.set_xlim([-limit, limit])
    ax.set_ylim([-limit, limit])
    ax.set_zlim([-limit, limit])
    
    ax.set_title('Discretized 3D Path')
    ax.legend()
    ax.grid(True)
    
    if show:
        plt.tight_layout()
        plt.show()


def color_function(f1, m3, n_points_width, a):

    stress = np.linspace(a/2,-a/2,n_points_width)*12*m3/a**2+f1
    return stress


def plot_ribbon_with_views(X_surf, Y_surf, Z_surf, path, colors, stress_cmap, min_value, max_value, axes=None):
    def plot_single_view(ax, elev, azim, title, show_colorbar=False):
        ax.clear()  # Clear previous contents
        
        # Apply surface colors
        normalized_colors = Normalize(vmin=min_value, vmax=max_value)(colors)
        rgba_colors = stress_cmap(normalized_colors)

    
        line = ax.plot(path[2], 
                      path[0], 
                      path[1], 
                      linestyle='--', 
                      color='red',
                      linewidth=2,
                      zorder = 10, 
                      solid_capstyle='round')  
        
        surf = ax.plot_surface(Z_surf, X_surf, Y_surf, 
                              facecolors=rgba_colors, 
                              alpha=1, 
                              cstride=1, 
                              rstride=1, 
                              zorder = 1, 
                              shade=False)
        
        if show_colorbar:
            mappable = plt.cm.ScalarMappable(cmap=stress_cmap)
            mappable.set_array(np.linspace(min_value, max_value, 100))
            if not hasattr(ax, 'colorbar'):
                ax.colorbar = plt.colorbar(mappable, ax=ax, shrink=0.45, aspect=10)
                ax.colorbar.set_label("External Stress [Pa]")
        
        ax.view_init(elev=elev, azim=azim)
        ax.xaxis.set_pane_color('gray')
        ax.yaxis.set_pane_color('gray')
        ax.zaxis.set_pane_color('gray')
        ax.set_xlabel('Z [mm]')
        ax.set_ylabel('X [mm]')
        ax.set_zlabel('Y [mm]')
        ax.set_xlim([-25, 25])
        ax.set_ylim([-25, 25])
        ax.set_zlim([-50, 0])
        ax.set_title(title)
    
    # If no axes are passed, create new ones
    if axes is None:
        fig = plt.figure(figsize=(18, 6))
        axes = [
            fig.add_subplot(131, projection='3d'),
            fig.add_subplot(132, projection='3d'),
            fig.add_subplot(133, projection='3d'),
        ]
    
    # Plot into the given axes
    plot_single_view(axes[0], elev=22.5, azim=-45, title="Isometric View")
    plot_single_view(axes[1], elev=5, azim=-90, title="Front View")
    plot_single_view(axes[2], elev=5, azim=0, title="Side View", show_colorbar=True)
    
    return axes



def Plot_stress_field(solutions, solutions_sleeve, path, L, a, n_points=20, axes=None, min_value=None, max_value=None):


    # Define colors for the stress colormap
    stress_colors = [
        (0.0, "blue"),    # Low value
        (0.25, "cyan"),   # Mid-low
        (0.5, "lime"),   # Middle
        (0.75, "yellow"), # Mid-high
        (1.0, "red")      # High value
    ]
    
    stress_colors = [
        (0.0, "blue"),    # Low value
        (0.5, "white"),   # Middle
        (1.0, "red")      # High value
    ]
    
    # Create the colormap
    stress_cmap = LinearSegmentedColormap.from_list("StressColormap", stress_colors)
    
    # Get the solution with the maximum index
    last_solution = solutions[solutions.Index_solution == solutions.Index_solution.max()]
    last_solution_sleeve = solutions_sleeve[solutions_sleeve.Index_solution == solutions_sleeve.Index_solution.max()]

    last_solution_sleeve = pd.concat([last_solution_sleeve[last_solution_sleeve.s ==0],
                                     last_solution_sleeve,
                                     last_solution_sleeve[last_solution_sleeve.s ==1]])

    # Set parameters for the ribbon width and quiver arrows
    half_width = a / 2  # Ribbon width
    arrow_scale = 0.1  # Scale for quiver arrows
    
    # Prepare arrays for surface plotting
    X_surf, Y_surf, Z_surf, colors = [], [], [], []

    # Iterate over all points in the solution
    for i in range(len(last_solution)):
        # Centerline point
        x, y, z = last_solution.iloc[i][['X', 'Y', 'Z']] * L


        # Direction vectors (normalize d1 and d2)
        d3 = last_solution.iloc[i][['d3x', 'd3y', 'd3z']].values
        d2 = last_solution.iloc[i][['d2x', 'd2y', 'd2z']].values
        d3 /= np.linalg.norm(d3)
        d2 /= np.linalg.norm(d2)
        d1 = np.cross(d2,d3)

        # Create points for the ribbon width
        width_points = []
        for j in range(n_points):
            # Interpolate points along the width of the ribbon
            m = (j - n_points / 2) / (n_points / 2)
            width_point = np.array([x, y, z]) + m * half_width * d2
            width_points.append(width_point)

        # Append the X, Y, Z coordinates for the ribbon width
        width_points = np.array(width_points)
        X_surf.append(width_points[:, 0])
        Y_surf.append(width_points[:, 1])
        Z_surf.append(width_points[:, 2])

        rx = (last_solution_sleeve.iloc[i+1].response_force_X + last_solution_sleeve.iloc[i].response_force_X)/2
        ry = (last_solution_sleeve.iloc[i+1].response_force_Y + last_solution_sleeve.iloc[i].response_force_Y)/2
        rz = (last_solution_sleeve.iloc[i+1].response_force_Z + last_solution_sleeve.iloc[i].response_force_Z)/2
        

        mx = (last_solution_sleeve.iloc[i+1].response_couple_X + last_solution_sleeve.iloc[i].response_couple_X)/2
        my = (last_solution_sleeve.iloc[i+1].response_couple_Y + last_solution_sleeve.iloc[i].response_couple_Y)/2
        mz = (last_solution_sleeve.iloc[i+1].response_couple_Z + last_solution_sleeve.iloc[i].response_couple_Z)/2

        
        f_n = rx*d1[0] + ry*d1[1] + rz*d1[2]
        m3 = mx*d3[0] + my*d3[1] + mz*d3[2]
        
        colors.append([color_function(f_n*1e6, m3*1e6,n_points, a)])

    # Convert to numpy arrays for surface plotting
    X_surf = np.array(X_surf)
    Y_surf = np.array(Y_surf)
    Z_surf = np.array(Z_surf)
    colors = np.concatenate(colors, axis=0)

    if axes is not None:
        axes = plot_ribbon_with_views(X_surf, Y_surf, Z_surf, path, colors, stress_cmap, min_value, max_value, axes)



def compute_global_stress_range(all_solutions, all_sleeves, a, lengths, n_points):
    all_colors = []

    for solutions, sleeve, L in zip(all_solutions, all_sleeves, lengths):
        last_solution = solutions[solutions.Index_solution == solutions.Index_solution.max()]
        last_solution_sleeve = sleeve[sleeve.Index_solution == sleeve.Index_solution.max()]
        last_solution_sleeve = pd.concat([
            last_solution_sleeve[last_solution_sleeve.s == 0],
            last_solution_sleeve,
            last_solution_sleeve[last_solution_sleeve.s == 1]
        ])

        for i in range(len(last_solution)):
            d3 = last_solution.iloc[i][['d3x', 'd3y', 'd3z']].values
            d2 = last_solution.iloc[i][['d2x', 'd2y', 'd2z']].values
            d3 /= np.linalg.norm(d3)
            d2 /= np.linalg.norm(d2)
            d1 = np.cross(d2, d3)

            rx = (last_solution_sleeve.iloc[i+1].response_force_X + last_solution_sleeve.iloc[i].response_force_X)/2
            ry = (last_solution_sleeve.iloc[i+1].response_force_Y + last_solution_sleeve.iloc[i].response_force_Y)/2
            rz = (last_solution_sleeve.iloc[i+1].response_force_Z + last_solution_sleeve.iloc[i].response_force_Z)/2
            mx = (last_solution_sleeve.iloc[i+1].response_couple_X + last_solution_sleeve.iloc[i].response_couple_X)/2
            my = (last_solution_sleeve.iloc[i+1].response_couple_Y + last_solution_sleeve.iloc[i].response_couple_Y)/2
            mz = (last_solution_sleeve.iloc[i+1].response_couple_Z + last_solution_sleeve.iloc[i].response_couple_Z)/2

            f_n = rx*d1[0] + ry*d1[1] + rz*d1[2]
            m3 = mx*d3[0] + my*d3[1] + mz*d3[2]

            stress = color_function(f_n*1e6, m3*1e6, n_points, a)
            all_colors.append(stress)

    all_colors = np.concatenate(all_colors)
    min_val = -abs(all_colors).max()
    max_val = abs(all_colors).max()
    if max_val - min_val < 1e-3:
        max_val += 1e-3
        min_val -= 1e-3
    return min_val, max_val


def animate_stress_field_3views(all_solutions, all_sleeves, all_path, lengths, a, max_frame, n_points=20, interval=200):
    min_val, max_val = compute_global_stress_range(all_solutions, all_sleeves, a, lengths, n_points)

    fig = plt.figure(figsize=(18, 6))
    ax1 = fig.add_subplot(131, projection='3d')
    ax2 = fig.add_subplot(132, projection='3d')
    ax3 = fig.add_subplot(133, projection='3d')
    axes = [ax1, ax2, ax3]

    def update(frame):
        solutions = all_solutions[frame]
        solutions_sleeve = all_sleeves[frame]
        path = all_path[frame]
        L = lengths[frame]

        Plot_stress_field(solutions, solutions_sleeve, path, L, a, n_points=n_points,
                          axes=axes, min_value=min_val, max_value=max_val)

        return axes

    ani = animation.FuncAnimation(fig, update, frames=range(max_frame + 1), interval=interval, blit=False)
    ani.save("stress_field_animation.gif", fps=5)
    plt.tight_layout()
    plt.show()