# Project 02 - NACA Airfoil simulations with OpenFOAM

This repository contains the code to automatically generate meshes, run OpenFOAM simulations, and post-process results for NACA 4-digit airfoils.

## Scripts

1. `simulation_utils.py`: Contains core functions for generating the airfoil, generating OpenFOAM dictionaries, and running Docker OpenFOAM commands.
2. `convergence_study.py`: Script to perform a grid convergence study based on the number of discretization points.
3. `camber_study.py`: Script to perform a parameter study of the impact of camber on the lift and drag coefficients.

## Requirements
- Python 3
- numpy, matplotlib, pandas
- Docker (with the `microfluidica/openfoam:13` image)

## Usage

To run the convergence study:
```bash
python convergence_study.py
```

To run the camber study:
```bash
python camber_study.py
```
