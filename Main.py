# Import necessary modules
import numpy as np
import os
import glob
from multiprocessing import Pool

import pandas as pd
from structure_generation.structgen import generate_am_structure, generate_se_structure, v_frac_after_combining
from percolation_analysis.percolating_clusters import extract_se_percolating_clusters_y_range, find_percolating_se_particles, extract_am_percolating_clusters_y_range, find_percolating_am_particles
from percolation_analysis.output_parameters import calculate_utilization_level, calculate_utilization_level_se, specific_surface_area_with_overlaps, active_interface_area3
from excel_saving_handler.excel_handler import get_workbook_and_sheet, save_data_to_excel, copy_sheet_to_final_excel


# Define the main Excel file for final data collection and a directory for temporary outputs.
excel_filename = 'Simulation_R1.xlsx'
temp_output_dir = 'temp_excel_outputs_R1'

# Active Material (AM) volume fractions to be simulated.
#am_fractions = [0.25,0.27,0.29,0.31,0.33,0.35,0.37,0.39,0.41,0.43,0.45,0.47,0.49,0.51,0.53,0.55,0.57,0.59,0.61,0.63,0.65]
#am_fractions = [0.35,0.40,0.45,0.50]
am_fractions = [0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65]


# AM particle parameters.
am_diameter = 5
am_mean_radius = am_diameter / 2  # Mean AM particle radius (e.g., 5 μm diameter).
overlap_fraction_am = 0.8

# Solid Electrolyte (SE) particle parameters.
se_diameter = 3
se_mean_radius = se_diameter / 2  # Mean SE particle radius (e.g., 3 μm diameter).
overlap_fraction_se = 0.6  # This value scales the sum of radii for overlap determination.

# Microstructure porosity.
porosity = 0.2

# Microstructure dimensions (x, y, z).
size = (80, 140, 80)

# Other simulation parameters.
# Maximum search range for potential overlaps around each particle.
max_radius = max(am_mean_radius,se_mean_radius)
max_range = 2 * max_radius

num_points = 1000  # Number of points for surface area calculations.

voxel_resolution = 0.2  # Side length of each voxel in simulation units.

# Initialize lists for storing particle coordinates (though managed within functions).
am_coordinates_list = []
se_coordinates_list = []
results = []
SEED = 1    # Random seed for reproducibility.

# Additional parameters, currently set to zero for fixed radius.
deviation_fraction = 0
am_std_radius = 0
se_std_radius = 0

# Thresholds for considering AM and SE fractions.
am_fraction_threshold = [0.0001]
se_fraction_threshold = [0.0001]

# Define how often particle coordinates should be saved to temporary Excel files during generation.
PARTICLE_SAVE_INTERVAL = 1000

# Compile all initial simulation parameters into a dictionary for easy access and logging.
common_simulation_parameters = {
    'AM Mean Radius': am_mean_radius,
    'Overlap Fraction AM': overlap_fraction_am,
    'SE Mean Radius': se_mean_radius,
    'Overlap Fraction SE': overlap_fraction_se,
    'Deviation Fraction': deviation_fraction,
    'AM Standard Deviation Radius': am_std_radius,
    'SE Standard Deviation Radius': se_std_radius,
    'Max Range': max_range,
    'Number of Points': num_points,
    'Microstructure Size': size,
    'porosity': porosity,
    'Voxel Resolution': voxel_resolution,
    'SEED': SEED,
    'AM fraction threshold': am_fraction_threshold,
    'SE fraction threshold': se_fraction_threshold,
    'Particle Save Interval': PARTICLE_SAVE_INTERVAL,
    'Temp Output Dir': temp_output_dir
}

# Convert the initial parameters dictionary into a DataFrame for structured saving.
initial_params_df = pd.DataFrame([common_simulation_parameters])


def work2(am_fraction_tuple):
    """
    Executes a single simulation run for a given AM fraction, generating
    microstructure, performing percolation analysis, and calculating key metrics.

    Args:
        am_fraction_tuple (tuple): A tuple containing:
            - am_fraction (float): The active material volume fraction for this run.
            - common_params (dict): A dictionary of common simulation parameters.

    Returns:
        dict: A dictionary containing the simulation results for the given AM fraction,
              including volume fractions, utilization levels, surface areas,
              interface areas, and generation times.
    """
    # Unpack the active material fraction and common simulation parameters.
    am_fraction, common_params = am_fraction_tuple

    # Extract individual parameters from the common_params dictionary.
    am_mean_radius = common_params['AM Mean Radius']
    overlap_fraction_am = common_params['Overlap Fraction AM']
    se_mean_radius = common_params['SE Mean Radius']
    overlap_fraction_se = common_params['Overlap Fraction SE']
    max_range = common_params['Max Range']
    num_points = common_params['Number of Points']
    size = common_params['Microstructure Size']
    porosity = common_params['porosity']
    voxel_resolution = common_params['Voxel Resolution']
    SEED = common_params['SEED']
    am_fraction_threshold = common_params['AM fraction threshold']
    se_fraction_threshold = common_params['SE fraction threshold']
    particle_save_interval = common_params['Particle Save Interval']
    temp_output_dir = common_params['Temp Output Dir']

    # Initialize output metrics.
    ul, ul_se, sa, sa2, aia, aia2, aia3 = 0, 0, 0, 0, 0, 0, 0
    
    # Calculate the Solid Electrolyte (SE) fraction based on porosity and AM fraction.
    #se_fraction = 1 - porosity - am_fraction
    se_fraction = 0
    elapsed_time_am = 0
    elapsed_time_se = 0

    # Define the path for the temporary Excel file for this specific run.
    temp_excel_path = os.path.join(temp_output_dir, f'{size}_temp_data_AM_{am_fraction:.2f}.xlsx')

    # Voxelization parameters: calculate grid dimensions and create an empty voxel grid.
    voxel_dimensions = np.ceil(np.array(size) / voxel_resolution).astype(int)
    voxel_grid = np.zeros(voxel_dimensions, dtype=int)

    print(f'Starting simulation for microstructure size: {size}')


    # --- Microstructure Generation (most computationally intensive section) ---

    # Generate the Active Material (AM) structure and populate the voxel grid.
    voxel_grid, am_coordinates, elapsed_time_am = generate_am_structure(
        am_mean_radius, size, am_fraction, max_range, overlap_fraction_am, SEED,
        am_fraction_threshold, voxel_resolution, voxel_dimensions, temp_excel_path, particle_save_interval
    )
    # Generate the Solid Electrolyte (SE) structure, incorporating existing AM particles.
    voxel_grid, se_coordinates, elapsed_time_se = generate_se_structure(
        se_mean_radius, size, se_fraction, voxel_grid, max_range, overlap_fraction_se,
        am_coordinates, SEED, se_fraction_threshold, am_fraction, voxel_resolution,
        voxel_dimensions, temp_excel_path, particle_save_interval
    )

    # --- Percolation Analysis and Metric Calculation ---

    # Calculate actual volume fractions after combining AM and SE structures.
    vf_am_AC, vf_se_AC = v_frac_after_combining(voxel_grid, size, voxel_dimensions)

    # Extract percolating clusters for SE and AM phases along the y-range.
    percolating_am_clusters_yrange, num_am_percolating_clusters_yrange, am_cluster_sizes_yrange = extract_am_percolating_clusters_y_range(voxel_grid, size, voxel_resolution, am_fraction)
    percolating_se_clusters_yrange, num_se_percolating_clusters_yrange, se_cluster_sizes_yrange = extract_se_percolating_clusters_y_range(voxel_grid, size, voxel_resolution, am_fraction)

    # Identify coordinates of percolating AM and SE particles.
    percolating_am_coordinates, _ = find_percolating_am_particles(am_coordinates, percolating_am_clusters_yrange, am_fraction, size, voxel_resolution)
    percolating_se_coordinates, _ = find_percolating_se_particles(se_coordinates, percolating_se_clusters_yrange, am_fraction, size, voxel_resolution)

    # Calculate utilization levels and active interface area if AM clusters percolate.
    if num_am_percolating_clusters_yrange:
        ul = calculate_utilization_level(am_cluster_sizes_yrange, voxel_grid)
        ul_se = calculate_utilization_level_se(se_cluster_sizes_yrange, voxel_grid)
        sa = 1
        sa2 = specific_surface_area_with_overlaps(percolating_am_coordinates, num_points, max_range, am_fraction)
        aia = 1
        aia2 = 1
        aia3 = active_interface_area3(percolating_am_coordinates, percolating_se_coordinates, num_points, max_range, am_fraction)

    # Compile all results into a dictionary for the current simulation run.
    result_dict = {
        'Size': size,
        'AM Fraction': am_fraction,
        'vf_am_AC': vf_am_AC,
        'vf_se_AC': vf_se_AC,
        'UL': ul,
        'UL_se': ul_se,
        'SA': sa,
        'SA2': sa2,
        'AIA': aia,
        'AIA2': aia2,
        'AIA3': aia3,
        'AM Gen time(s)': elapsed_time_am,
        'SE Gen time(s)': elapsed_time_se,
        'Temp_Excel_Path': temp_excel_path  # Return the path to the temporary file for later consolidation.
    }

    # Save the run summary to a temporary Excel file for individual runs.
    try:
        temp_results_df = pd.DataFrame([result_dict])
        temp_results_df = temp_results_df.drop(columns=['Temp_Excel_Path'], errors='ignore')    # Exclude 'Temp_Excel_Path' column when saving the summary to the temporary file itself.
        workbook_summary, sheet_summary = get_workbook_and_sheet(temp_excel_path, 'Run_Summary', list(temp_results_df.columns)) # Use utility functions to get the workbook/sheet and save data.
        save_data_to_excel(workbook_summary, sheet_summary, temp_results_df.values.tolist(), temp_excel_path)
    except Exception as e:
        print(f"Error saving run summary to temporary file {temp_excel_path}: {e}")

    print(f'Completed simulation for microstructure size: {size}')
    print(f'DONE AM_fraction: {am_fraction:.2f}')
    print()

    return result_dict


if __name__ == "__main__":
    # Ensure the directory for temporary Excel files exists.
    os.makedirs(temp_output_dir, exist_ok=True)

    # Clean up any existing temporary files from previous runs to ensure a clean slate.
    for f in glob.glob(os.path.join(temp_output_dir, '*_temp_data_AM_*.xlsx')):
        os.remove(f)

    # Initialize the main Excel file by saving the common simulation parameters.
    workbook_initial, sheet_initial = get_workbook_and_sheet(excel_filename, 'Initial Parameters', list(initial_params_df.columns))
    save_data_to_excel(workbook_initial, sheet_initial, initial_params_df.values.tolist(), excel_filename)

    # Prepare arguments for multiprocessing: each AM fraction gets its own set of parameters.
    # This creates a list of tuples: (am_fraction, common_params_for_pool).
    pool_args = [(af, common_simulation_parameters) for af in am_fractions]

    # Use a multiprocessing Pool to run simulations in parallel for different AM fractions.
    with Pool() as pool:
        results = pool.map(work2, pool_args)

    # Consolidate all simulation results into the main Excel file.
    # Load the main workbook for consolidation; 'Results' sheet will be created/updated.
    final_workbook, _ = get_workbook_and_sheet(excel_filename, 'Results')

    # Create a DataFrame from the collected results, excluding the temporary file path.
    results_df_summary = pd.DataFrame([
        {k: v for k, v in res.items() if k not in ('Temp_Excel_Path')} for res in results
    ], columns=[
        'Size', 'AM Fraction', 'vf_am_AC', 'vf_se_AC', 'UL', 'UL_se', 'SA', 'SA2',
        'AIA', 'AIA2', 'AIA3', 'AM Gen time(s)', 'SE Gen time(s)'
    ])

    # Save the consolidated summary results to the 'Results' sheet in the main Excel file.
    workbook_results, sheet_results = get_workbook_and_sheet(excel_filename, 'Results', list(results_df_summary.columns))
    save_data_to_excel(workbook_results, sheet_results, results_df_summary.values.tolist(), excel_filename)

    # Iterate through the results and copy detailed particle coordinate data (AM and SE)
    # from temporary files to dedicated sheets in the final Excel file.
    with pd.ExcelWriter(excel_filename, mode='a', engine='openpyxl', if_sheet_exists='overlay') as writer:
        for res in results:
            am_fraction = res['AM Fraction']
            temp_excel_path = res['Temp_Excel_Path']

            # Copy AM and SE particle data sheets from temporary files.
            copy_sheet_to_final_excel(temp_excel_path, f'AM_{am_fraction:.2f}', writer)
            copy_sheet_to_final_excel(temp_excel_path, f'SE_{am_fraction:.2f}', writer)

    print("Simulation complete. All results consolidated into", excel_filename)