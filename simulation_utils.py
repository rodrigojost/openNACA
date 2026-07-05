import os
import numpy as np
import math
from pathlib import Path
import subprocess
import pandas as pd
from curiosityFluidsAirfoilMesher import curiosityFluidsAirfoilMesher

def generate_naca4(digits, chord=1.0, points=100):
    """Generates 2D coordinates for a 4-digit NACA airfoil."""
    m = int(digits[0]) / 100.0  # Max camber
    p = int(digits[1]) / 10.0   # Position of max camber
    t = int(digits[2:]) / 100.0 # Max thickness

    # Use cosine spacing to cluster more points at the leading edge for a smoother curve
    beta = np.linspace(0, np.pi, points)
    x = (1 - np.cos(beta)) / 2

    # Thickness equation
    # (Note: -0.1015 leaves a tiny trailing edge gap, standard for NACA.)
    yt = 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x**2 + 0.2843 * x**3 - 0.1015 * x**4)

    # Camber line calculations
    yc = np.zeros_like(x)
    dyc_dx = np.zeros_like(x)
    
    if p > 0:
        # Ahead of max camber
        front = x <= p
        yc[front] = (m / p**2) * (2 * p * x[front] - x[front]**2)
        dyc_dx[front] = (2 * m / p**2) * (p - x[front])
        
        # Behind max camber
        back = x > p
        yc[back] = (m / (1 - p)**2) * ((1 - 2 * p) + 2 * p * x[back] - x[back]**2)
        dyc_dx[back] = (2 * m / (1 - p)**2) * (p - x[back])
        
    theta = np.arctan(dyc_dx)
    
    # Upper surface coordinates
    xu = x - yt * np.sin(theta)
    yu = yc + yt * np.cos(theta)
    
    # Lower surface coordinates
    xl = x + yt * np.sin(theta)
    yl = yc - yt * np.cos(theta)
    
    # Combine into one continuous loop (Trailing Edge -> Leading Edge -> Trailing Edge)
    x_coords = np.concatenate([xu[::-1], xl[1:]])
    y_coords = np.concatenate([yu[::-1], yl[1:]])

    # FIX for symmetric (M=0) airfoils:
    # The curiosityFluidsAirfoilMesher splits the profile into top/bottom by checking
    # whether consecutive X values increase or decrease. For a perfectly symmetric
    # airfoil the trailing-edge points (top and bottom) share Y=0, which makes the
    # boundary-layer blocks collapse and creates negative-volume cells.
    # It's applied a microscopic vertical offset (+/- 0.5% chord * thickness) to the
    # last point of each surface so the mesher sees a non-zero trailing-edge gap.
    # This is physically invisible but prevents the mesh singularity.
    if m == 0:
        te_offset = 0.0005 * chord * t
        # Upper TE (first point in concatenated array)
        y_coords[0] += te_offset
        # Lower TE (last point in concatenated array)
        y_coords[-1] -= te_offset
    
    # Scale by chord length
    return np.column_stack((x_coords, y_coords)) * chord

def setup_case_directories(case_dir):
    """Sets up the directory structure for an OpenFOAM case."""
    os.makedirs(os.path.join(case_dir, "constant"), exist_ok=True)
    os.makedirs(os.path.join(case_dir, "system"), exist_ok=True)
    # write empty file for use with paraview
    with open(os.path.join(case_dir, "case.foam"), "w") as f:
        f.write("")

def create_control_dict(case_dir, end_time=0.6, u_inf=30, rho_inf=1.225, z_length=1.0, write_interval=0.1):
    system_dir = os.path.join(case_dir, "system")
    os.makedirs(system_dir, exist_ok=True)
    area = 1.0 * z_length
    
    content = f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Website:  https://openfoam.org                  |
|   \\\\  /    A nd           | Version:  13                                    |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    format      ascii;
    class       dictionary;
    location    "system";
    object      controlDict;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
application     foamRun;
solver          incompressibleFluid;

startFrom       startTime;
startTime       0;

stopAt          endTime;
endTime         {end_time};

deltaT          0.0005; // Set conservative for initialization stability

writeControl    runTime;
writeInterval   {write_interval};

purgeWrite      0;
writeFormat     ascii;
writePrecision  6;
writeCompression off;
timeFormat      general;
timePrecision   6;
runTimeModifiable true;

functions
{{
    computeLiftDrag
    {{
        type            forceCoeffs;      // Must be an unquoted word
        libs            ("libforces.so"); // Shared libraries CAN be a string list
        writeControl    timeStep;
        writeInterval   1;

        patches         (airfoil);        // Name of your boundary patch in patches/polyMesh
        p               p;                // Pressure field name
        U               U;                // Velocity field name

        // For incompressible flows (using foamRun with incompressibleFluid module)
        rho             rhoInf;
        rhoInf          {rho_inf}; // reference density

        CofR            (0 0 0);          // Center of Rotation (for moment coefficients)

        // Direction vectors (Adjust these based on your coordinate system orientation!)
        liftDir         (0 1 0);          // If lift is pointing up the Y-axis
        dragDir         (1 0 0);          // If drag is pointing along the X-axis
        pitchAxis       (0 0 1);          // Pitch axis (Z-axis)

        // Reference values used to normalize forces into Cd and Cl
        magUInf         {u_inf};          // Freestream velocity magnitude
        lRef            1.0;              // Reference length (e.g., chord length of airfoil) 
        Aref            {area};           // Planform Area (Chord 1.0m * Span 0.1m)
    }}

    residuals
    {{
        type            residuals;
        libs            ("libutilityFunctionObjects.so");
        writeControl    timeStep;
        writeInterval   1;
        fields          (U p);
    }}
}}
// ************************************************************************* //
"""
    dict_path = os.path.join(system_dir, "controlDict")
    with open(dict_path, "w") as f:
        f.write(content)
    return dict_path

def create_momentum_transport_sa(case_dir):
    constant_dir = os.path.join(case_dir, "constant")
    os.makedirs(constant_dir, exist_ok=True)
    
    content = f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Website:  https://openfoam.org                  |
|   \\\\  /    A nd           | Version:  13                                    |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    format      ascii;
    class       dictionary;
    location    "constant";
    object      momentumTransport;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

simulationType RAS;

RAS
{{
    model           SpalartAllmaras;
    turbulence      on;
    printCoeffs     on;
}}
// ************************************************************************* //
"""
    dict_path = os.path.join(constant_dir, "momentumTransport")
    with open(dict_path, "w") as f:
        f.write(content)
    return dict_path

def generate_openfoam_fvschemes_fvsolution(case_dir, steady_state=True):
    system_dir = os.path.join(case_dir, "system")
    
    ddt_scheme = "steadyState" if steady_state else "Euler"
    
    fv_schemes = f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Website:  https://openfoam.org                  |
|   \\\\  /    A nd           | Version:  13                                    |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSchemes;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

ddtSchemes
{{
    default         {ddt_scheme};
}}

gradSchemes
{{
    default         Gauss linear;
    grad(U)         Gauss linear;
}}

divSchemes
{{
    default         none;
    div(phi,U)      Gauss linearUpwind grad(U);
    div(phi,nuTilda) Gauss upwind;
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}}

laplacianSchemes
{{
    default         Gauss linear corrected;
}}

interpolationSchemes
{{
    default         linear;
}}

snGradSchemes
{{
    default         corrected;
}}

wallDist
{{
    method          meshWave;
}}
// ************************************************************************* //
"""
    n_outer_correctors = 1 if steady_state else 2
    relaxation_factors = """relaxationFactors
{
    fields
    {
        p               0.3;
    }
    equations
    {
        U               0.7;
        nuTilda         0.7;
    }
}""" if steady_state else ""

    fv_solution = f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Website:  https://openfoam.org                  |
|   \\\\  /    A nd           | Version:  13                                    |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      fvSolution;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

solvers
{{
    p
    {{
        solver          GAMG;
        tolerance       1e-06;
        relTol          0.01;
        smoother        GaussSeidel;
    }}
    pFinal
    {{
        $p;
        tolerance       1e-06;
        relTol          0;
    }}
    "(U|nuTilda)"
    {{
        solver          PBiCGStab;
        preconditioner  DILU;
        tolerance       1e-05;
        relTol          0.1;
    }}
    "(U|nuTilda)Final"
    {{
        solver          PBiCGStab;
        preconditioner  DILU;
        tolerance       1e-05;
        relTol          0;
    }}
}}

PIMPLE
{{
    nOuterCorrectors {n_outer_correctors};
    nCorrectors      2;
    nNonOrthogonalCorrectors 2;
}}

{relaxation_factors}
// ************************************************************************* //
"""
    with open(os.path.join(system_dir, "fvSchemes"), "w") as f:
        f.write(fv_schemes)
    with open(os.path.join(system_dir, "fvSolution"), "w") as f:
        f.write(fv_solution)

def generate_physical_properties(case_dir, nu_value=1.5e-5):
    constant_dir = os.path.join(case_dir, "constant")
    
    file_content = f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Website:  https://openfoam.org                  |
|   \\\\  /    A nd           | Version:  13                                    |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    format      ascii;
    class       dictionary;
    location    "constant";
    object      physicalProperties;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

viscosityModel  constant;

nu              [0 2 -1 0 0 0 0] {nu_value};
// ************************************************************************* //
"""
    file_path = os.path.join(constant_dir, "physicalProperties")
    with open(file_path, "w") as f:
        f.write(file_content)

def generate_0_directory(case_dir, u_freestream, rho_freestream, nu=1.5e-5):
    zero_dir = os.path.join(case_dir, "0")
    os.makedirs(zero_dir, exist_ok=True)
    
    nu_tilda_freestream = 3.0 * nu
    cv1 = 7.1
    chi = nu_tilda_freestream / nu
    fv1 = (chi**3) / (chi**3 + cv1**3)
    nut_freestream = nu_tilda_freestream * fv1
    
    u_content = f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Website:  https://openfoam.org                  |
|   \\\\  /    A nd           | Version:  13                                    |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       volVectorField;
    location    "0";
    object      U;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 1 -1 0 0 0 0];

internalField   uniform ({u_freestream} 0 0);

boundaryField
{{
    farfield
    {{
        type            freestream;
        freestreamValue uniform ({u_freestream} 0 0);
    }}
    airfoil
    {{
        type            noSlip;
    }}
    frontAndBack
    {{
        type            empty;
    }}
}}
// ************************************************************************* //
"""

    p_content = f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Website:  https://openfoam.org                  |
|   \\\\  /    A nd           | Version:  13                                    |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      p;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 2 -2 0 0 0 0];

internalField   uniform 0;

boundaryField
{{
    farfield
    {{
        type            freestreamPressure;
        freestreamValue uniform 0;
    }}
    airfoil
    {{
        type            zeroGradient;
    }}
    frontAndBack
    {{
        type            empty;
    }}
}}
// ************************************************************************* //
"""

    nutilda_content = f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Website:  https://openfoam.org                  |
|   \\\\  /    A nd           | Version:  13                                    |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      nuTilda;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 2 -1 0 0 0 0];

internalField   uniform {nu_tilda_freestream};

boundaryField
{{
    farfield
    {{
        type            freestream;
        freestreamValue uniform {nu_tilda_freestream};
    }}
    airfoil
    {{
        type            fixedValue;
        value           uniform 0;
    }}
    frontAndBack
    {{
        type            empty;
    }}
}}
// ************************************************************************* //
"""

    nut_content = f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Website:  https://openfoam.org                  |
|   \\\\  /    A nd           | Version:  13                                    |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      nut;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 2 -1 0 0 0 0];

internalField   uniform {nut_freestream};

boundaryField
{{
    farfield
    {{
        type            calculated;
        value           uniform {nut_freestream};
    }}
    airfoil
    {{
        type            nutUSpaldingWallFunction;
        value           uniform 0;
    }}
    frontAndBack
    {{
        type            empty;
    }}
}}
// ************************************************************************* //
"""
    fields = {"U": u_content, "p": p_content, "nuTilda": nutilda_content, "nut": nut_content}
    for field_name, script_content in fields.items():
        file_path = os.path.join(zero_dir, field_name)
        with open(file_path, "w") as f:
            f.write(script_content)

OPENFOAM_IMAGE = "microfluidica/openfoam:13"
PROJECT_DIR = Path.cwd().resolve()

def foam_cmd(*args, case_dir="openfoam_naca", workdir="/work"):
    # Runs the OpenFOAM command directly instead of using Docker.
    # This requires running the script in an environment where OpenFOAM is installed (like WSL).
    return [
        *args,
        "-case",
        case_dir
    ]

def parse_force_coeffs(case_dir="openfoam_naca"):
    forcecoeffs_file = os.path.join(case_dir, "postProcessing", "computeLiftDrag", "0", "forceCoeffs.dat")
    if not os.path.exists(forcecoeffs_file):
        return None
    # Read the space-delimited file, ignoring the OpenFOAM header lines starting with '#'
    # Column 0 is Time, Column 1 is Cd, Column 3 is Cl
    df = pd.read_csv(forcecoeffs_file, sep=r'\s+', comment='#', header=None)
    time = df[0].values
    cd = df[2].values
    cl = df[3].values
    return time, cd, cl

def parse_residuals(case_dir="openfoam_naca"):
    residuals_file = os.path.join(case_dir, "postProcessing", "residuals", "0", "residuals.dat")
    if not os.path.exists(residuals_file):
        return None
    df = pd.read_csv(residuals_file, sep=r'\s+', comment='#', header=None)
    time = df[0].values
    r_ux = df[1].values
    r_uy = df[2].values
    r_p = df[3].values
    return time, r_ux, r_uy, r_p

def setup_and_run_simulation(points, airfoil_code="2412", end_time=0.6, case_dir="openfoam_naca"):
    """Complete workflow for a single simulation run."""
    chord_length = 1.0
    rho_freestream = 1.225
    u_freestream = 30
    
    # 1. Generate Airfoil
    discretized_surface_points = generate_naca4(airfoil_code, chord=chord_length, points=points)
    
    # 2. Setup OpenFOAM Case Directory
    setup_case_directories(case_dir)
    
    # 3. Meshing
    curiosityFluidsAirfoilMesher(discretized_surface_points, case_dir=case_dir)
    
    # 4. Generate OpenFOAM Dictionaries
    create_control_dict(case_dir, end_time=end_time, u_inf=u_freestream, rho_inf=rho_freestream)
    create_momentum_transport_sa(case_dir)
    generate_openfoam_fvschemes_fvsolution(case_dir)
    generate_physical_properties(case_dir)
    generate_0_directory(case_dir, u_freestream, rho_freestream)
    
    # 5. Run OpenFOAM Commands
    subprocess.run(foam_cmd("blockMesh", case_dir=case_dir), check=True)
    subprocess.run(foam_cmd("checkMesh", case_dir=case_dir), check=False) # eventual mesh warnings are printed but foamRun is still attempted.
    subprocess.run(foam_cmd("foamRun", case_dir=case_dir), check=True)
