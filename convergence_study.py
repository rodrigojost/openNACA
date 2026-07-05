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
        # The mesh and results are preserved for further inspection.
        case_dir = f"convergence_{airfoil_code}_{pts}pts"
        print(f"\n--- Running simulation with {pts} points -> {case_dir} ---")

        setup_and_run_simulation(points=pts, airfoil_code=airfoil_code,
                                 end_time=end_time, case_dir=case_dir)

        time, cd, cl = parse_force_coeffs(case_dir=case_dir)

        final_cd = cd[-1]
        final_cl = cl[-1]

        steady_cd_list.append(final_cd)
        steady_cl_list.append(final_cl)

        print(f"Results for {pts} pts -> Cd: {final_cd:.4f}, Cl: {final_cl:.4f}")

    # ── Plot 1: Cl and Cd vs. number of points (dual Y axes) ──────────────────
    fig, ax1 = plt.subplots(figsize=(10, 5))
    color_cl = "steelblue"
    color_cd = "crimson"

    ax1.set_xlabel("Number of Discretization Points")
    ax1.set_ylabel("Lift Coefficient  $C_l$", color=color_cl)
    ax1.plot(points_list, steady_cl_list, 'o-', color=color_cl, linewidth=2,
             label="$C_l$")
    ax1.tick_params(axis='y', labelcolor=color_cl)

    ax2 = ax1.twinx()
    ax2.set_ylabel("Drag Coefficient  $C_d$", color=color_cd)
    ax2.plot(points_list, steady_cd_list, 's--', color=color_cd, linewidth=2,
             label="$C_d$")
    ax2.tick_params(axis='y', labelcolor=color_cd)

    # Combined legend
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="center right")

    fig.suptitle(f"Grid Convergence Study  (NACA {airfoil_code})")
    ax1.grid(True, linestyle="--", alpha=0.5)
    fig.tight_layout()
    fig.savefig("convergence_plot.png", dpi=150)
    print("\nConvergence plot saved as 'convergence_plot.png'.")

    # ── Plot 2 & 3: Force coefficients and residuals over time (finest grid) ──
    best_pts     = points_list[-1]
    best_case    = f"convergence_{airfoil_code}_{best_pts}pts"
    time, cd, cl = parse_force_coeffs(case_dir=best_case)
    time_r, r_ux, r_uy, r_p = parse_residuals(case_dir=best_case)

    # Force coefficients over time – dual Y axes
    fig, ax1 = plt.subplots(figsize=(10, 5))
    ax1.set_xlabel("Time (s)")
    ax1.set_ylabel("Lift Coefficient  $C_l$", color=color_cl)
    ax1.plot(time, cl, color=color_cl, linewidth=2, label="$C_l$")
    ax1.tick_params(axis='y', labelcolor=color_cl)

    ax2 = ax1.twinx()
    ax2.set_ylabel("Drag Coefficient  $C_d$", color=color_cd)
    ax2.plot(time, cd, color=color_cd, linewidth=2, linestyle="--", label="$C_d$")
    ax2.tick_params(axis='y', labelcolor=color_cd)

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper right")

    fig.suptitle(f"Force Coefficients over Time  (NACA {airfoil_code}, {best_pts} pts)")
    ax1.grid(True, linestyle="--", alpha=0.5)
    fig.tight_layout()
    fig.savefig("force_coefficients_time.png", dpi=150)

    # Residuals over time (log scale, single axis is fine)
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.semilogy(time_r, r_p,  color='black', alpha=0.85, label='Pressure $p$')
    ax.semilogy(time_r, r_ux, color='steelblue', alpha=0.85, label='$U_x$')
    ax.semilogy(time_r, r_uy, color='seagreen',  alpha=0.85, label='$U_y$')
    ax.set_title(f"Solver Residuals  (NACA {airfoil_code}, {best_pts} pts)")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Initial Residual  (log scale)")
    ax.grid(True, which="both", linestyle="--", alpha=0.5)
    ax.legend()
    fig.tight_layout()
    fig.savefig("residuals_time.png", dpi=150)

    print("Time convergence plots saved as 'force_coefficients_time.png'"
          " and 'residuals_time.png'.")

if __name__ == "__main__":
    main()
