# Particle-in-a-box-approach-to-Modeling-Solid-State-Batteries

## Project Overview

This repository hosts a Python-based simulation framework designed for the microstructural modeling of composite cathodes in all-solid-state batteries (ASSBs). The core objective is to understand the intricate relationship between the microstructure of battery components and their resulting electrochemical properties and overall performance. Leveraging percolation theory, the model simulates the distribution and connectivity of active material (AM) and solid electrolyte (SE) particles within a defined 3D volume (a "particle-in-a-box" approach) to predict key metrics such as utilization level, specific surface area, and active interface area.

This work is an attempt to develop an efficient and accurate open-source tool for exploring the vast parameter space of SSB microstructures, thereby accelerating research and development in this promising energy storage technology.

## Motivation

Solid-state batteries offer significant advantages over traditional liquid-electrolyte batteries, including enhanced safety, higher energy density, and improved cycle life. However, their performance is heavily influenced by the microstructure of the composite cathode, particularly the high resistance at the solid-solid interfaces between AM and SE. Effective cathode design requires a deep understanding of how factors like particle size, composition, porosity, and electrode thickness impact ionic and electronic conduction pathways. This simulation framework aims to provide this understanding, guiding the optimization of cathode microstructures for improved SSB performance.
## Key Features

* **3D Microstructure Generation:** Simulates the random placement of spherical Active Material (AM) and Solid Electrolyte (SE) particles within a cubic simulation box.
* **Percolation Analysis:** Identifies connected clusters of AM and SE particles that span the entire simulation domain, crucial for assessing effective ionic and electronic conductivity.
* **Property Calculation:** Computes critical metrics for battery performance assessment:
    * **Utilization Level (UL):** Quantifies the efficiency of material usage within percolating clusters for both AM & SE.
    * **Specific Surface Area (SSA):** Calculates the exposed surface area of percolating AM particles, accounting for self-overlaps using an improved "Golden Spiral Method" for point distribution. 
    * **Active Interface Area (AIA):** Determines the contact area between AM and SE particles, vital for electrochemical reactions, with multiple calculation methods to account for complex overlaps.
* **Parallel Processing:** Utilizes Python's `multiprocessing` library to distribute simulation runs for different AM fractions across multiple CPU cores, significantly reducing overall runtime.
* **Excel Output:** Organizes and exports all simulation parameters, particle coordinates, and calculated results to structured Excel files for easy analysis and visualization.

## Code Structure

The project is modularized into several Python files, each responsible for a specific aspect of the simulation:

* `Main.py`: The primary script orchestrating the simulation workflow, including parameter setup, parallel execution, and results consolidation.
* `structgen.py`: Contains functions for generating the 3D microstructure, including particle placement and voxelization.
* `output_parameters.py`: Implements the algorithms for calculating utilization level, specific surface area, and active interface area.
* `percolating_clusters.py`: Focuses on identifying and extracting percolating clusters within the generated voxel grid.
* `excel_handler.py`: Manages all Excel-related operations, such as creating workbooks, sheets, writing data, and applying styles.

## How It Works (Simplified Workflow)

1.  **Initialization:** Define global simulation parameters like microstructure dimensions, particle sizes, porosity, and desired AM volume fractions.
2.  **Parallel Simulation Runs:** For each specified AM volume fraction, a separate simulation run is initiated in parallel.
3.  **Structure Generation:** Inside each run, AM and SE particles are randomly placed within the 3D box, ensuring adherence to defined overlap constraints. The space is then voxelized.
4.  **Percolation Analysis:** The voxelized structure is analyzed to identify connected pathways (percolating clusters) of both AM (for electronic conduction) and SE (for ionic conduction).
5.  **Property Calculation:** Based on the identified percolating clusters, metrics like utilization level, specific surface area, and active interface area are computed.
6.  **Data Export:** All raw particle data and calculated results are systematically saved to temporary Excel files.
7.  **Consolidation:** After all parallel runs are complete, a main script consolidates the data from all temporary files into a single, comprehensive Excel workbook for final analysis.

## Installation and Usage

### Prerequisites

* Python 3.x
* `numpy`
* `pandas`
* `openpyxl`
* `tqdm`
* `scipy`

You can install the required packages using pip:

```bash
pip install numpy pandas openpyxl tqdm scipy
