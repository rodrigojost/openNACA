# OpenNACA — OpenFOAM NACA Airfoil Aerodynamic Parameter Study

This project implements a fully automated CFD pipeline using **OpenFOAM 13** and **Python** to study the aerodynamic performance of 4-digit NACA airfoils. A structured C-grid mesh is generated automatically from the airfoil geometry, incompressible turbulent flow is solved using the Spalart–Allmaras model, and lift and drag coefficients are extracted and plotted. We conduct a **grid convergence study** and a **camber parameter study** (NACA 0412 → 8412) at zero angle of attack.

---

## Installation & Environment Setup

OpenFOAM 13 is natively supported on Linux and macOS. On Windows, WSL2 is required.

### 1. Install OpenFOAM 13

**Linux / WSL2 (Ubuntu):**
```bash
wget -qO - https://dl.openfoam.org/gpg.key | sudo tee /etc/apt/trusted.gpg.d/openfoam.asc > /dev/null
sudo add-apt-repository http://dl.openfoam.org/ubuntu
sudo apt-get update && sudo apt-get install -y openfoam13
echo "source /opt/openfoam13/etc/bashrc" >> ~/.bashrc
source ~/.bashrc
```

**macOS:** Follow the [official OpenFOAM macOS instructions](https://openfoam.org/download/).

**Docker (cross-platform alternative):**
```bash
docker run -it --rm -v $(pwd):/work -w /work microfluidica/openfoam:13
```

### 2. Set Up the Python Environment

Create and activate a virtual environment and install the required dependencies:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install numpy matplotlib pandas
```

---

## Running the Studies

Ensure your virtual environment is activated (`source .venv/bin/activate`) and that OpenFOAM commands are on your PATH before running.

### 1. Run the Grid Convergence Study

Runs NACA 2412 at four mesh resolutions (20, 50, 100, 200 surface points), extracts steady-state $C_l$ and $C_d$ for each, and plots convergence:

```bash
python3 convergence_study.py
```

Outputs:
- `convergence_plot.png` — $C_l$ and $C_d$ vs. number of discretization points (dual Y axes)
- `force_coefficients_time.png` — force coefficients over simulation time (finest grid)
- `residuals_time.png` — solver residuals over time (log scale)
- `convergence_2412_20pts/`, `convergence_2412_50pts/`, ... — full OpenFOAM case directories

### 2. Run the Camber Parameter Study

Uses the grid resolution and simulation time determined from the convergence study. Sweeps camber M from 0 to 8 (NACA 0412 → 8412) at fixed P=4, TT=12:

```bash
python3 camber_study.py
```

Outputs:
- `camber_study_plot.png` — $C_l$ and $C_d$ vs. maximum camber M (dual Y axes)
- `camber_0412/`, `camber_1412/`, ..., `camber_8412/` — full OpenFOAM case directories

---

## Project Structure

```
OpenNACA/
├── README.md                          # Setup and run instructions (this file)
├── .gitignore                         # Ignores OpenFOAM case dirs, raw data, and plots
├── simulation_utils.py                # Core library: NACA geometry, OF dictionaries,
│                                      #   command runner, results parser
├── curiosityFluidsAirfoilMesher.py    # Structured C-grid blockMeshDict generator
│                                      #   (GPLv3, adapted from curiosityFluids)
├── convergence_study.py               # Task 1 & 2: grid + time convergence study
└── camber_study.py                    # Task 3: camber parameter sweep M=0..8
```

---

## Notes on the Symmetric Airfoil (NACA 0412)

The `curiosityFluidsAirfoilMesher` builds boundary-layer blocks by extruding surface normals outward. For a perfectly symmetric profile (M=0), the top and bottom trailing-edge points share Y=0, which causes the two BL blocks to collapse into each other, producing negative-volume cells. A microscopic trailing-edge offset (`±0.05% chord`) is applied automatically in `simulation_utils.py` to break this degeneracy without affecting the aerodynamic result. The solver is also configured with `nNonOrthogonalCorrectors 2` for additional robustness.

---

## Authors & Course Context

Developed as part of **Project 02 — OpenFOAM NACA Airfoil Simulations** for the course  
**11.00153 Modern Simulation Software Development** (Dr. Lambert Theisen & Dr. Georgii Oblapenko), Summer Semester 2026.
