# structgen.py
import numpy as np
from tqdm import tqdm
import time
from scipy.spatial import cKDTree

from excel_saving_handler.excel_handler import get_workbook_and_sheet, save_data_to_excel


# Overlap, Generate particle, Voxelisation functions 
def check_overlap2(coordinates, center, radius, max_range, overlap_fraction):
    if not coordinates:
        return True
    
    # Convert center to a NumPy array once
    center_np = np.array(center)

    tree = cKDTree([c[0] for c in coordinates])
    neighbors_indices = tree.query_ball_point(center_np, max_range + radius)

    if not neighbors_indices:
        return True # No neighbors within the extended search range

    # Extract existing centers and radii for relevant neighbors
    existing_data = [coordinates[i] for i in neighbors_indices]
    existing_centers = np.array([data[0] for data in existing_data])
    existing_radii = np.array([data[1] for data in existing_data])

    distances_sq = np.sum((existing_centers - center_np)**2, axis=1)
    distances = np.sqrt(distances_sq) 

    critical_distances = (overlap_fraction) * (existing_radii + radius)

    # Check for any overlap using vectorized comparison
    # We compare distances directly with critical_distances
    if np.any(distances < critical_distances):
        return False
    
    return True

def generate_particle(coordinates, mean_radius, size, max_range, overlap_fraction):
    particle_coordinate=[]
    while True:
        radius = mean_radius
        center = np.random.rand(3) * np.array(size)

        if check_overlap2(coordinates, center, radius, max_range, overlap_fraction):
            particle_coordinate.append((center, radius))
            return particle_coordinate

def update_voxel_grid(particle_coordinates, sizes, voxel_grid, voxel_resolution, voxel_dimensions, particle_type):
    for center, radius in particle_coordinates:
        center = np.array(center)
        min_bounds = np.floor((center - radius) / voxel_resolution).astype(int)
        max_bounds = np.ceil((center + radius) / voxel_resolution).astype(int)
        min_bounds = np.clip(min_bounds, 0, voxel_dimensions - 1)
        max_bounds = np.clip(max_bounds, 0, voxel_dimensions - 1)
        voxel_positions = np.mgrid[min_bounds[0]:max_bounds[0] + 1,
                                   min_bounds[1]:max_bounds[1] + 1,
                                   min_bounds[2]:max_bounds[2] + 1].reshape(3, -1).T

        distances = np.linalg.norm(center - (voxel_positions + 0.5) * voxel_resolution, axis=1)
        inside_sphere = distances <= radius

        if particle_type == 'a':
            voxel_grid[tuple(voxel_positions[inside_sphere].T)] = 1
        elif particle_type == 's':
            voxel_positions_inside = voxel_positions[inside_sphere]
            not_already_occupied = voxel_grid[tuple(voxel_positions_inside.T)] == 0
            voxel_grid[tuple(voxel_positions_inside.T[:, not_already_occupied])] = -1
    return voxel_grid


# Removed excel_maker function as it's replaced by get_workbook_and_sheet in excel_handler.py

def generate_am_structure(am_mean_radius, size, am_fraction, max_range, overlap_fraction_am, SEED, am_fraction_threshold, voxel_resolution, voxel_dimensions, temp_excel_path, save_interval):
    am_coordinates = []
    vf_am = 0
    voxel_grid_am = np.zeros(voxel_dimensions, dtype=int)
    particle_buffer = []

    progress_bar = tqdm(total=am_fraction, desc=f'AM Volume Fraction {am_fraction:.2f}', unit='VF')
    np.random.seed(SEED)

    sheet_name = f'AM_{am_fraction:.2f}'
    headers = ['AM Fraction', 'x', 'y', 'z', 'radius'] # Define headers here
    workbook, sheet_main = get_workbook_and_sheet(temp_excel_path, sheet_name, headers) # Use new function

    for am_fraction_threshold_1 in am_fraction_threshold:
        start_time = time.time()
        elapsed_time = 0

        while vf_am < am_fraction:
            if am_fraction - vf_am < am_fraction_threshold_1:
                print("Threshold reached, AM fraction = ",vf_am)
                break

            am_coordinate = generate_particle(am_coordinates, am_mean_radius, size, max_range, overlap_fraction_am)
            am_coordinates.extend(am_coordinate)
            particle_buffer.append([am_fraction] + am_coordinate[0][0].tolist() + [am_coordinate[0][1]])

            if len(particle_buffer) % save_interval == 0:
                save_data_to_excel(workbook, sheet_main, particle_buffer, temp_excel_path) # Use new function
                particle_buffer = []
            
            particle_type='a'
            voxel_grid_am = update_voxel_grid(am_coordinate, size, voxel_grid_am, voxel_resolution, voxel_dimensions, particle_type)
            vf_am = np.count_nonzero(voxel_grid_am==1) / np.prod(voxel_dimensions)
            progress_bar.update(vf_am - progress_bar.n)
  
        
        if particle_buffer:
            save_data_to_excel(workbook, sheet_main, particle_buffer, temp_excel_path) # Use new function
            particle_buffer = []

        elapsed_time = time.time() - start_time

    progress_bar.close()
    return voxel_grid_am, am_coordinates, elapsed_time


def generate_se_structure(se_mean_radius, size, se_fraction, voxel_grid_am, max_range, overlap_fraction_se, am_coordinates, SEED, se_fraction_threshold, am_fraction, voxel_resolution, voxel_dimensions, temp_excel_path, save_interval):
    se_coordinates = []
    vf_se = 0
    particle_buffer = []

    progress_bar = tqdm(total=se_fraction, desc=f'SE Volume Fraction {se_fraction:.2f} for {am_fraction:.2f}', unit='VF')
    np.random.seed(SEED)

    sheet_name = f'SE_{am_fraction:.2f}'
    headers = ['AM Fraction', 'x', 'y', 'z', 'radius'] # Define headers here
    workbook, sheet_main = get_workbook_and_sheet(temp_excel_path, sheet_name, headers) # Use new function

    for se_fraction_threshold_1 in se_fraction_threshold:
        start_time = time.time()
        elapsed_time = 0

        while vf_se < se_fraction:
            if se_fraction - vf_se < se_fraction_threshold_1:
                print("Threshold reached, SE fraction = ",vf_se)
                break

            se_coordinate = generate_particle(am_coordinates + se_coordinates, se_mean_radius, size, max_range, overlap_fraction_se)
            se_coordinates.extend(se_coordinate)
            particle_buffer.append([am_fraction] + se_coordinate[0][0].tolist() + [se_coordinate[0][1]])

            if len(particle_buffer) % save_interval == 0:
                save_data_to_excel(workbook, sheet_main, particle_buffer, temp_excel_path) # Use new function
                particle_buffer = []

            particle_type = 's'
            voxel_grid_am = update_voxel_grid(se_coordinate, size, voxel_grid_am, voxel_resolution, voxel_dimensions, particle_type)
            vf_se = np.count_nonzero(voxel_grid_am==-1) / np.prod(voxel_dimensions)
            progress_bar.update(vf_se - progress_bar.n)

        if particle_buffer:
            save_data_to_excel(workbook, sheet_main, particle_buffer, temp_excel_path) # Use new function
            particle_buffer = []

        elapsed_time = time.time() - start_time

    progress_bar.close()
    return voxel_grid_am, se_coordinates, elapsed_time


def v_frac_after_combining(voxel_grid, size, voxel_dimensions):
    v_frac_am_AC = np.count_nonzero(voxel_grid == 1) / np.prod(voxel_dimensions)
    v_frac_se_AC = (np.count_nonzero(voxel_grid == -1)) / np.prod(voxel_dimensions)
    return v_frac_am_AC, v_frac_se_AC