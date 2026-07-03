import os
import matplotlib.pyplot as plt
from simulation_utils import setup_and_run_simulation, parse_force_coeffs, parse_residuals

def main():
    points_list = [20, 50, 100, 200]
    airfoil_code = "2412"
    end_time = 0.6
    
    steady_cd_list = []
    steady_cl_list = []
    
    print("Starting Grid Convergence Study...")
    
    for pts in points_list:
        case_dir = f"openfoam_naca_{pts}"
        print(f"\n--- Running simulation with {pts} points ---")
        
        setup_and_run_simulation(points=pts, airfoil_code=airfoil_code, end_time=end_time, case_dir=case_dir)
        
        time, cd, cl = parse_force_coeffs(case_dir=case_dir)
        
        # Get the final steady state values
        final_cd = cd[-1]
        final_cl = cl[-1]
        
        steady_cd_list.append(final_cd)
        steady_cl_list.append(final_cl)
        
        print(f"Results for {pts} points -> Cd: {final_cd:.4f}, Cl: {final_cl:.4f}")
        
    # Plotting Grid Convergence
    plt.figure(figsize=(10, 5))
    plt.plot(points_list, steady_cl_list, 'o-', color='blue', label='Lift Coefficient (Cl)')
    plt.title(f"Grid Convergence Study (NACA {airfoil_code})")
    plt.xlabel("Number of Discretization Points")
    plt.ylabel("Steady-State Lift Coefficient")
    plt.grid(True, linestyle="--", alpha=0.7)
    plt.legend()
    plt.tight_layout()
    plt.savefig("convergence_plot_cl.png")
    
    plt.figure(figsize=(10, 5))
    plt.plot(points_list, steady_cd_list, 'o-', color='red', label='Drag Coefficient (Cd)')
    plt.title(f"Grid Convergence Study (NACA {airfoil_code})")
    plt.xlabel("Number of Discretization Points")
    plt.ylabel("Steady-State Drag Coefficient")
    plt.grid(True, linestyle="--", alpha=0.7)
    plt.legend()
    plt.tight_layout()
    plt.savefig("convergence_plot_cd.png")
    
    print("\nConvergence study completed. Plots saved as 'convergence_plot_cl.png' and 'convergence_plot_cd.png'.")
    
    # Task 2: Plotting residuals and force coeffs over time for the finest grid (to verify simulation time)
    best_pts = points_list[-1]
    best_case_dir = f"openfoam_naca_{best_pts}"
    time, cd, cl = parse_force_coeffs(case_dir=best_case_dir)
    time_r, r_ux, r_uy, r_p = parse_residuals(case_dir=best_case_dir)
    
    # Force Coefficients Over Time
    plt.figure(figsize=(10, 5))
    plt.plot(time, cd, label="Drag Coefficient (Cd)", color="red", linewidth=2)
    plt.plot(time, cl, label="Lift Coefficient (Cl)", color="blue", linewidth=2)
    plt.title("Aerodynamic Force Coefficients Over Time (Finest Grid)")
    plt.xlabel("Time (s)")
    plt.ylabel("Coefficient Value")
    plt.grid(True, linestyle="--", alpha=0.7)
    plt.legend()
    plt.tight_layout()
    plt.savefig("force_coefficients_time.png")
    
    # Residuals Over Time
    plt.figure(figsize=(10, 6))
    plt.plot(time_r, r_p, label='Pressure', color='black', alpha=0.8)
    plt.plot(time_r, r_ux, label='Velocity Ux', color='blue', alpha=0.8)
    plt.plot(time_r, r_uy, label='Velocity Uy', color='green', alpha=0.8)
    plt.yscale('log')
    plt.title("Solver Residuals Over Time (Finest Grid)")
    plt.xlabel("Time (s)")
    plt.ylabel("Residual (Log Scale)")
    plt.grid(True, which="both", linestyle="--", alpha=0.5)
    plt.legend()
    plt.tight_layout()
    plt.savefig("residuals_time.png")
    
    print("Time convergence plots saved as 'force_coefficients_time.png' and 'residuals_time.png'.")

if __name__ == "__main__":
    main()
