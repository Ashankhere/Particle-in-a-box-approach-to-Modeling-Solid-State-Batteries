import numpy as np
from tqdm import tqdm
from scipy.spatial import cKDTree


# ## Utilization Level

def calculate_utilization_level(am_cluster_sizes_yrange, voxel_grid):
    """
    Calculates the utilization level for AM based on the sum of AM cluster sizes
    and the total number of non-zero voxels marked as 1.
    """
    UL_am = np.sum(am_cluster_sizes_yrange) / (np.count_nonzero(voxel_grid == 1))
    return UL_am


def calculate_utilization_level_se(se_cluster_sizes_yrange, voxel_grid):
    """
    Calculates the utilization level for SE based on the sum of SE cluster sizes
    and the total number of non-zero voxels marked as -1.
    """
    UL = np.sum(se_cluster_sizes_yrange) / (np.count_nonzero(voxel_grid == -1))
    return UL


def generate_uniform_points_on_sphere(center, radius, num_points):
    """
    Generates uniformly distributed points on the surface of a sphere
    using the Fibonacci sphere algorithm.

    Args:
        center (tuple or list): The (x, y, z) coordinates of the sphere's center.
        radius (float): The radius of the sphere.
        num_points (int): The number of points to generate.

    Returns:
        numpy.ndarray: An array of shape (num_points, 3) representing the
                       coordinates of the points on the sphere.
    """
    indices = np.arange(0, num_points, dtype=float) + 0.5
    phi = np.arccos(1 - 2 * indices / num_points)
    theta = np.pi * (1 + 5**0.5) * indices

    x = center[0] + radius * np.sin(phi) * np.cos(theta)
    y = center[1] + radius * np.sin(phi) * np.sin(theta)
    z = center[2] + radius * np.cos(phi)

    points_on_sphere = np.column_stack((x, y, z))
    return points_on_sphere

def build_kd_tree(coordinates_with_types):
    """
    Builds a cKDTree from a list of coordinates, where each coordinate
    entry can also contain additional data (e.g., particle type).

    Args:
        coordinates_with_types (list): A list of tuples, where each tuple
                                       starts with (center, radius, ...)
                                       and center is (x, y, z).

    Returns:
        scipy.spatial.cKDTree: The KD-tree built from the centers.
    """
    # Extract the centers of the coordinates and convert to NumPy array
    existing_centers = np.array([c[0] for c in coordinates_with_types])
    kdtree = cKDTree(existing_centers)
    return kdtree


def specific_surface_area_with_overlaps(percolating_am_coordinates, num_points, max_range, am_fraction):
    """
    Calculates the specific surface area of AM particles, considering
    self-overlaps (i.e., exposed surface only).
    """
    # Consolidate all coordinates with their type for a single KDTree
    all_am_coords = [(*coords, 'am') for coords in percolating_am_coordinates]
    kdtree_all = build_kd_tree(all_am_coords)

    percolating_surface_area = 0

    for current_index, (center, radius) in enumerate(tqdm(percolating_am_coordinates, desc=f"Processing Coordinates for SA2 {am_fraction:.2f}")):
        points_on_sphere = generate_uniform_points_on_sphere(center, radius, num_points)
        neighbor_indices = kdtree_all.query_ball_point(center, max_range)

        am_neighbors_data = []
        for idx in neighbor_indices:
            if idx == current_index: # Exclude the current particle itself
                continue
            am_neighbors_data.append((all_am_coords[idx][0], all_am_coords[idx][1])) # (center, radius)

        if not am_neighbors_data: # If no AM neighbors, this sphere is fully exposed
            percolating_surface_area += (4 * np.pi * radius**2)
            continue

        am_neighbor_centers = np.array([n[0] for n in am_neighbors_data])
        am_neighbor_radii = np.array([n[1] for n in am_neighbors_data])

        distances_matrix = np.linalg.norm(points_on_sphere[:, np.newaxis, :] - am_neighbor_centers[np.newaxis, :, :], axis=2)
        is_inside_any_am_neighbor = np.any(distances_matrix < am_neighbor_radii[np.newaxis, :], axis=1)

        exposed_points_count = np.sum(~is_inside_any_am_neighbor)
        sphere_area = 4 * np.pi * radius**2
        percolating_surface_area += (exposed_points_count / num_points) * sphere_area

    return percolating_surface_area


def active_interface_area3(percolating_am_coordinates, percolating_se_coordinates, num_points, max_range, am_fraction):
    """
    Calculates the active interface area between AM and SE particles,
    considering AM-AM self-overlap (i.e., only exposed AM surface that touches SE).
    """
    all_coords = [(*coords, 'se') for coords in percolating_se_coordinates] + \
                 [(*coords, 'am') for coords in percolating_am_coordinates]
    kdtree_all = build_kd_tree(all_coords)

    total_active_interface_area = 0

    for current_am_center, current_am_radius in tqdm(percolating_am_coordinates, desc=f"Processing Coordinates for AIA3 {am_fraction:.2f} (Points:{num_points})"):

        points_on_sphere = generate_uniform_points_on_sphere(current_am_center, current_am_radius, num_points)
        neighbor_indices = kdtree_all.query_ball_point(current_am_center, max_range)

        am_neighbors_data = []
        se_neighbors_data = []
        for idx in neighbor_indices:
            neighbor_center, neighbor_radius, neighbor_type = all_coords[idx]
            # Check if this neighbor is the current particle
            if np.array_equal(neighbor_center, current_am_center) and neighbor_radius == current_am_radius:
                continue

            if neighbor_type == 'am':
                am_neighbors_data.append((neighbor_center, neighbor_radius))
            elif neighbor_type == 'se':
                se_neighbors_data.append((neighbor_center, neighbor_radius))

        point_status = np.zeros(num_points, dtype=int) # 0: untouched, -1: inside AM neighbor, 1: inside SE neighbor

        # Pass 1: Check for overlaps with AM neighbors (if any point is inside an AM neighbor, it's marked -1)
        if am_neighbors_data:
            am_neighbor_centers = np.array([n[0] for n in am_neighbors_data])
            am_neighbor_radii = np.array([n[1] for n in am_neighbors_data])
            distances_to_am = np.linalg.norm(points_on_sphere[:, np.newaxis, :] - am_neighbor_centers[np.newaxis, :, :], axis=2)
            is_inside_am = np.any(distances_to_am < am_neighbor_radii[np.newaxis, :], axis=1)
            point_status[is_inside_am] = -1

        # Pass 2: Check for overlaps with SE neighbors (only for points not already marked -1)
        if se_neighbors_data:
            se_neighbor_centers = np.array([n[0] for n in se_neighbors_data])
            se_neighbor_radii = np.array([n[1] for n in se_neighbors_data])

            unmarked_points_indices = np.where(point_status == 0)[0]
            if len(unmarked_points_indices) > 0:
                distances_to_se = np.linalg.norm(points_on_sphere[unmarked_points_indices, np.newaxis, :] - se_neighbor_centers[np.newaxis, :, :], axis=2)
                is_inside_se_unmarked = np.any(distances_to_se < se_neighbor_radii[np.newaxis, :], axis=1)
                point_status[unmarked_points_indices[is_inside_se_unmarked]] = 1

        contact_points_count = np.sum(point_status == 1)
        sphere_area = 4 * np.pi * current_am_radius**2
        total_active_interface_area += (contact_points_count / num_points) * sphere_area

    return total_active_interface_area