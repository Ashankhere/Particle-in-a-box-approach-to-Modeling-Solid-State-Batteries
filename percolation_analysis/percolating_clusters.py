import numpy as np
from tqdm import tqdm
from scipy.ndimage import label

def _find_percolating_particles_core(coordinates, percolating_clusters_yrange, am_fraction, size, particle_type_name, threshold, voxel_resolution):
    """
    Core logic for finding percolating particles (AM or SE).
    Optimized for speed using pure NumPy vectorization.
    """
    voxel_dimensions = np.ceil(np.array(size) / voxel_resolution).astype(int)

    percolating_output_coordinates = []
    sus = []

    # Handle empty coordinates list gracefully
    if not coordinates:
        return [], []

    try:
        # Flatten the list of (center, radius) tuples into a 2D numpy array
        flattened_coords = [list(c) + [r] for c, r in coordinates]
        centers_radii = np.array(flattened_coords, dtype=np.float64)
    except Exception as e:
        print(f"Error: Could not process particle coordinates. Ensure 'center' is a 3-element iterable and 'radius' is a scalar. Error: {e}")
        return [], []

    # Use progress_bar directly in the loop as shown in percolating_clusters_11.py
    for i in tqdm(range(len(coordinates)), desc=f'Finding Percolating {particle_type_name} Particles {am_fraction:.2f}', unit='Particle'):
        center = centers_radii[i, :3]
        radius = centers_radii[i, 3]
        radius_sq = radius ** 2 # Pre-calculate squared radius for efficiency

        # Calculate the integer voxel bounds for the sphere
        min_bounds = np.floor((center - radius) / voxel_resolution).astype(int)
        max_bounds = np.ceil((center + radius) / voxel_resolution).astype(int)

        # Ensure bounds are within the voxel grid dimensions
        min_bounds = np.clip(min_bounds, 0, voxel_dimensions - 1)
        max_bounds = np.clip(max_bounds, 0, voxel_dimensions - 1)

        # Create ranges for voxel coordinates
        x_coords = np.arange(min_bounds[0], max_bounds[0] + 1)
        y_coords = np.arange(min_bounds[1], max_bounds[1] + 1)
        z_coords = np.arange(min_bounds[2], max_bounds[2] + 1)

        if not (x_coords.size > 0 and y_coords.size > 0 and z_coords.size > 0):
            sus.append(0.0)
            continue

        # Pure NumPy vectorized way to generate voxel_positions_bbox
        num_x = x_coords.size
        num_y = y_coords.size
        num_z = z_coords.size

        total_voxels_in_bbox = num_x * num_y * num_z
        voxel_positions_bbox = np.empty((total_voxels_in_bbox, 3), dtype=np.int32)

        voxel_positions_bbox[:, 0] = np.repeat(x_coords, num_y * num_z)
        voxel_positions_bbox[:, 1] = np.tile(np.repeat(y_coords, num_z), num_x)
        voxel_positions_bbox[:, 2] = np.tile(z_coords, num_x * num_y)

        # Calculate actual (float) coordinates for these voxels' centers
        voxel_real_coords = (voxel_positions_bbox + 0.1) * voxel_resolution

        # Calculate squared distances from the sphere's center to all these voxel real coordinates
        distances_sq = np.sum((voxel_real_coords - center)**2, axis=1)

        # Filter for voxels that are actually inside the sphere using squared distance
        valid_positions = voxel_positions_bbox[distances_sq <= radius_sq]

        if len(valid_positions) == 0:
            sus.append(0.0)
            continue

        # Get the boolean status of these valid voxels from the percolating cluster grid
        percolating_voxel_status = percolating_clusters_yrange[valid_positions[:, 0], valid_positions[:, 1], valid_positions[:, 2]]

        # Sum up the number of filled voxels
        filled_voxels = np.sum(percolating_voxel_status)
        total_voxels_inside = len(valid_positions)

        # Avoid division by zero
        ratio = filled_voxels / total_voxels_inside if total_voxels_inside > 0 else 0.0
        sus.append(ratio)

        if ratio > threshold:
            percolating_output_coordinates.append(coordinates[i])

    return percolating_output_coordinates, sus


def _extract_percolating_clusters_y_range_common(voxel_grid, voxel_material_value, start_y_coord, end_y_coord, material_name, am_fraction):
    """
    Optimized common logic for extracting percolating clusters along the y-axis.
    A cluster is considered percolating if it connects both start_y_coord and end_y_coord.
    """
    # 1. Label clusters for the specified material value
    clusters, num_clusters = label((voxel_grid == voxel_material_value))

    # Handle case with no clusters
    if num_clusters == 0:
        return np.zeros(voxel_grid.shape, dtype=bool), 0, []

    # 2. Get unique cluster labels at the start and end Y coordinates
    # We only care about the cluster labels, not their specific (x, z) positions
    
    # Extract slices for the start and end y-coordinates
    # Use .ravel() to flatten to 1D array of labels, then np.unique
    start_y_labels = np.unique(clusters[:, start_y_coord, :][clusters[:, start_y_coord, :] != 0])
    end_y_labels = np.unique(clusters[:, end_y_coord, :][clusters[:, end_y_coord, :] != 0])

    # 3. Find percolating cluster labels (those present in both start and end y slices)
    percolating_cluster_labels = np.intersect1d(start_y_labels, end_y_labels)

    # Handle case with no percolating clusters
    if len(percolating_cluster_labels) == 0:
        return np.zeros(voxel_grid.shape, dtype=bool), 0, []

    # 4. Create the boolean mask for percolating clusters
    # Initialize the mask as False
    percolating_clusters_yrange = np.zeros(voxel_grid.shape, dtype=bool)

    # Use isin for vectorized checking of all voxels against the percolating labels
    # This is much faster than looping through each cluster
    percolating_clusters_yrange = np.isin(clusters, percolating_cluster_labels)

    # 5. Calculate cluster sizes more efficiently
    cluster_sizes = []
    # If cluster sizes are indeed needed, we can iterate over the identified
    # percolating_cluster_labels and use np.count_nonzero
    # This loop is much shorter than the original if only a few clusters percolate
    for label_val in tqdm(percolating_cluster_labels, desc=f'Extracting {material_name} clusters for {am_fraction:.2f}', unit='cluster'):
        cluster_sizes.append(np.count_nonzero(clusters == label_val))

    num_percolating_clusters_yrange = len(percolating_cluster_labels)

    return percolating_clusters_yrange, num_percolating_clusters_yrange, cluster_sizes


def extract_se_percolating_clusters_y_range(voxel_grid, size, voxel_resolution, am_fraction):
    start_y_voxel = 0 # Starting boundary for y
    end_y_voxel = np.ceil(size[1] / voxel_resolution).astype(int) - 1 # Ending boundary for y
    return _extract_percolating_clusters_y_range_common(voxel_grid, -1, start_y_voxel, end_y_voxel, 'SE', am_fraction)


def extract_am_percolating_clusters_y_range(voxel_grid, size, voxel_resolution, am_fraction):
    start_y_voxel = 0 # Starting boundary for y
    end_y_voxel = np.ceil(size[1] / voxel_resolution).astype(int) - 1 # Ending boundary for y
    return _extract_percolating_clusters_y_range_common(voxel_grid, 1, start_y_voxel, end_y_voxel, 'AM', am_fraction)


def find_percolating_se_particles(se_coordinates, percolating_clusters_se_yrange, am_fraction, size, voxel_resolution):
    # Reuse the core logic
    return _find_percolating_particles_core(se_coordinates, percolating_clusters_se_yrange, am_fraction, size, 'SE', 0, voxel_resolution)


def find_percolating_am_particles(am_coordinates, percolating_clusters_yrange, am_fraction, size, voxel_resolution):
    # Reuse the core logic
    return _find_percolating_particles_core(am_coordinates, percolating_clusters_yrange, am_fraction, size, 'AM', 0, voxel_resolution)