import os
import matplotlib.pyplot as plt
from simulation_utils import setup_and_run_simulation, parse_force_coeffs

def main():
    # Assuming from the convergence study that 100 points and 0.6s is sufficient
    # (Users can adjust these values based on the convergence study results)
    best_points = 100
    end_time = 0.6
    
    # Camber study from M=0 to M=8, P=4, TT=12 -> 0412 to 8412
    camber_values = list(range(0, 9))
    airfoil_codes = [f"{m}412" for m in camber_values]
    
    steady_cd_list = []
    steady_cl_list = []
    
    print("Starting Camber Parameter Study...")
    
    for code, m in zip(airfoil_codes, camber_values):
        case_dir = f"openfoam_naca_{code}"
        print(f"\n--- Running simulation for NACA {code} ---")
        
        setup_and_run_simulation(points=best_points, airfoil_code=code, end_time=end_time, case_dir=case_dir)
        
        time, cd, cl = parse_force_coeffs(case_dir=case_dir)
        
        final_cd = cd[-1]
        final_cl = cl[-1]
        
        steady_cd_list.append(final_cd)
        steady_cl_list.append(final_cl)
        
        print(f"Results for NACA {code} -> Cd: {final_cd:.4f}, Cl: {final_cl:.4f}")
        
    # Plotting Camber Study
    plt.figure(figsize=(10, 5))
    plt.plot(camber_values, steady_cl_list, 'o-', color='blue', label='Lift Coefficient (Cl)')
    plt.plot(camber_values, steady_cd_list, 's-', color='red', label='Drag Coefficient (Cd)')
    plt.title("Impact of Camber on Lift and Drag Coefficients")
    plt.xlabel("Maximum Camber (%)")
    plt.ylabel("Steady-State Coefficient Value")
    plt.grid(True, linestyle="--", alpha=0.7)
    plt.legend()
    plt.tight_layout()
    plt.savefig("camber_study_plot.png")
    
    print("\nCamber parameter study completed. Plot saved as 'camber_study_plot.png'.")

if __name__ == "__main__":
    main()
