# Ribbon Interaction Simulation

A modified PyElastica framework for simulating elastic ribbon mechanics and interaction with hyperelastic media using 1D ribbon models.

## Overview

This project implements a numerical simulation framework for analyzing stresses and deformations in elastic ribbons navigating through hyperelastic media. It combines the 1D ribbon model by Audoly & Neukirch with PyElastica's time-integration capabilities to study ribbon-medium interactions.

## Key Features

- **1D Ribbon Model**: Implementation of the Audoly-Neukirch ribbon model for efficient simulation
- **Hyperelastic Contact**: Ogden-type contact forces for ribbon-medium interaction
- **Quasi-static Analysis**: Optimized for obtaining static equilibrium configurations
- **Cluster Computing**: Designed for high-performance computing environments

## Installation

**Required packages:**
- numpy 2.2.6
- scipy 1.15.3
- matplotlib 3.10.3
- numba 0.61.2
- plotly 6.1.1
- pandas 2.2.3
- tqdm 4.67.1
- kaleido 0.2.1

### Setup

1. Clone the repository:
```bash
git clone https://github.com/pierremaillar/Ribbon_interaction.git
cd Ribbon_interaction
```

2. Install dependencies:
```bash
pip install -U -r requirements.txt
```

3. Run simulations directly from source - no additional installation required.

## Usage

### Guides

Follow the Jupyter notebooks guide:

- **`Noto_guide_Auto.ipynb`**: How to generate solution using Auto-07P. These solutions are used to validate the model implemented in PyElastica.
- **`Noto_elastica_ribbon_model.ipynb`**: How to setup and use the Ribbon model for simple load cases.
- **`Noto_elastica_sleeve_interaction.ipynb`**: How to implement the media interaction using the `Sleeve` object. 
- **`Noto_elastica_full_trajectory.ipynb`**: How to perform a full trajectory analysis.
- 
### Basic Simulation

The framework provides several script types:

- **`create_sim_env.py`**: Initialize simulation environment
- **`run_*.py`**: Execute simulations with specific parameters
- **`sim_*.run`**: Sbatch scripts for cluster submission

### Example Workflow

1. Configure simulation parameters in a `run_*.py` script
2. Submit to cluster: `sbatch sim_trajectory.run`
3. Monitor progress: `squeue` and `sjob {job_id}`
4. Analyze results using provided visualization tools

## Key Components

### Ribbon1D Class
- Implementation of 1D ribbon mechanics
- Constitutive laws based on plate theory
- Compatible with PyElastica's time-stepping framework

### Sleeve Contact Model
- Hyperelastic contact forces using Ogden model
- Torque computation for bending and twisting resistance
- Quasi-static equilibrium solver

### Memory Block Architecture
- Optimized contiguous memory allocation
- Numba JIT compilation for performance
- Vectorized operations across multiple ribbons

## Performance Notes

- **Recommended**: 16 CPU cores per simulation
- **Runtime**: ~6 hours for complete trajectory analysis
- **Stability**: Requires small time steps due to explicit integration
- **Convergence**: Manual tuning of end time and damping required

## Cluster Usage

Compatible with SCITAS clusters (JED, Helvetios). See Section 7 of the user guide for detailed setup instructions.

### Quick Start on Cluster
```bash
# Connect to cluster
ssh username@jed.hpc.epfl.ch

# Create environment
python -m venv --system-site-packages venvs/ribbon-sim
source venvs/ribbon-sim/bin/activate

# Clone and setup
git clone https://github.com/pierremaillar/Ribbon_interaction.git
cd Ribbon_interaction
pip install -U -r requirements.txt

# Submit job
sbatch sim_trajectory.run
```

## Limitations

- **Contact Model**: Simplified Ogden-type indentation assumptions
- **Time Integration**: Explicit scheme requires very small time steps
- **Convergence**: Manual monitoring required for quasi-static equilibrium
- **Angles**: Large deformation accuracy limited for twist/bend > 45°

## Future Improvements

- FEM-based contact model for improved accuracy
- Adaptive time stepping with automatic convergence criteria
- GPU parallelization for enhanced performance
- Better load balancing for cluster computing

## Documentation

For detailed technical documentation, see the complete user guide (`User_guide___elastica___auto.pdf`).

## Citation

Based on the 1D ribbon model from:
> Audoly, B., & Neukirch, S. (2021). A one-dimensional model for elastic ribbons: A little stretching makes a big difference. *Journal of the Mechanics and Physics of Solids*, 153, 104457.

## Contact

**Author**: Pierre Maillard  
**Institution**: Master student in Mechanical Engineering  
**Supervisor**: Lorenzo Noseda  
**Date**: June 26, 2025

---

⚠️ **Note**: This is research-grade software. Validate results carefully and consider the documented limitations for your specific use case.
