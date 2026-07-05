import os
import matplotlib.pyplot as plt
from simulation_utils import setup_and_run_simulation, parse_force_coeffs

def main():
    # Values chosen from the convergence study results
    best_points = 100
    end_time = 0.6

    # M = 0..8, NACA 0412 … 8412
    camber_values = list(range(0, 9))
    airfoil_codes = [f"{m}412" for m in camber_values]

    valid_cambers = []
    steady_cd_list = []
    steady_cl_list = []

    print("Starting Camber Parameter Study...")

    for code, m in zip(airfoil_codes, camber_values):
        # Case directory is named after the airfoil code so the mesh is kept
        case_dir = f"camber_{code}"
        print(f"\n--- Running simulation for NACA {code} ---")

        try:
            setup_and_run_simulation(points=best_points, airfoil_code=code,
                                     end_time=end_time, case_dir=case_dir)

            result = parse_force_coeffs(case_dir=case_dir)
            if result is None:
                raise RuntimeError("No forceCoeffs.dat found – foamRun may have failed silently.")

            time, cd, cl = result
            final_cd = cd[-1]
            final_cl = cl[-1]

            valid_cambers.append(m)
            steady_cd_list.append(final_cd)
            steady_cl_list.append(final_cl)

            print(f"Results for NACA {code} -> Cd: {final_cd:.4f}, Cl: {final_cl:.4f}")

        except Exception as e:
            print(f"[WARNING] NACA {code} failed: {e}. Skipping.")

    if not valid_cambers:
        print("No successful simulations – nothing to plot.")
        return

    # ── Camber study plot: dual Y axes ────────────────────────────────────────
    color_cl = "steelblue"
    color_cd = "crimson"

    fig, ax1 = plt.subplots(figsize=(10, 5))
    ax1.set_xlabel("Maximum Camber  M (%)")
    ax1.set_ylabel("Lift Coefficient  $C_l$", color=color_cl)
    ax1.plot(valid_cambers, steady_cl_list, 'o-', color=color_cl, linewidth=2,
             label="$C_l$")
    ax1.tick_params(axis='y', labelcolor=color_cl)

    ax2 = ax1.twinx()
    ax2.set_ylabel("Drag Coefficient  $C_d$", color=color_cd)
    ax2.plot(valid_cambers, steady_cd_list, 's--', color=color_cd, linewidth=2,
             label="$C_d$")
    ax2.tick_params(axis='y', labelcolor=color_cd)

    # Combined legend
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper left")

    ax1.set_xticks(valid_cambers)
    fig.suptitle("Impact of Camber on Lift and Drag Coefficients  (NACA M412, AoA=0°)")
    ax1.grid(True, linestyle="--", alpha=0.5)
    fig.tight_layout()
    fig.savefig("camber_study_plot.png", dpi=150)

    print(f"\nCamber study completed ({len(valid_cambers)}/{len(camber_values)} cases).")
    print("Plot saved as 'camber_study_plot.png'.")

if __name__ == "__main__":
    main()
