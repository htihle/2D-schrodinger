#include <iostream>
#include <armadillo>
#include <fstream>
#include <cmath>

using namespace arma;
using namespace std;

// Function to normalize the wavefunction
void normalize(cx_mat& psi, double dx) {
    double norm = sqrt(accu(abs(psi) % abs(psi)) * dx * dx);
    if (norm > 1e-10) {
        psi = psi / norm;
    }
}

// Function to calculate kinetic energy expectation value
double kinetic_energy(const cx_mat& psi, double hbar, double m, double dx) {
    int n = psi.n_rows;
    double KE = 0.0;

    for (int i = 1; i < n-1; i++) {
        for (int j = 1; j < n-1; j++) {
            complex<double> laplacian = (psi(i+1,j) + psi(i-1,j) + psi(i,j+1) + psi(i,j-1) - 4.0*psi(i,j)) / (dx*dx);
            KE += real(conj(psi(i,j)) * laplacian);
        }
    }

    return -0.5 * hbar * hbar / m * KE * dx * dx;
}

// Function to calculate potential energy expectation value
double potential_energy(const cx_mat& psi, const mat& V_original, double dx, double m, double hbar) {
    int n = psi.n_rows;
    double PE = 0.0;

    for (int i = 1; i < n-1; i++) {
        for (int j = 1; j < n-1; j++) {
            // V_original is already the physical potential (no conversion needed)
            double V_phys = V_original(i,j);
            PE += V_phys * real(conj(psi(i,j)) * psi(i,j));
        }
    }

    return PE * dx * dx;
}

int main()
{
    // Physical parameters
    double m = 1.0;
    double hbar = 1.0;
    double omega = 30.0;
    double omega2 = omega * omega;
    double x0 = 0.5;  // Center of harmonic oscillator

    // Numerical parameters
    int n = 100;      // Grid points (increased from 30)
    int nt = 2000;    // Time steps
    int max_iter = 1000;  // Max iterations for convergence
    double tol = 1e-8;    // Convergence tolerance

    double dx = 1.0 / (n - 1);
    double dt = 0.00005;  // Smaller time step for better accuracy

    // Complex unit
    complex<double> i_unit(0, 1.0);
    complex<double> alpha = i_unit * hbar * dt / (2.0 * m * dx * dx);

    // Initialize grids
    cx_mat psi = cx_mat(n, n, fill::zeros);
    cx_mat psi_new = cx_mat(n, n, fill::zeros);
    cx_mat psi_old = cx_mat(n, n, fill::zeros);
    mat V(n, n, fill::zeros);
    mat V_original(n, n, fill::zeros);  // Store original V for energy calculation

    vec xx = linspace(0, 1, n);
    vec yy = linspace(0, 1, n);

    // Set up potential (2D harmonic oscillator)
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double V_phys = 0.5 * m * omega2 * ((xx(i) - x0) * (xx(i) - x0) +
                                                 (yy(j) - x0) * (yy(j) - x0));
            V_original(i, j) = V_phys;
            // Dimensionless potential for numerical scheme
            V(i, j) = 2 * m * dx * dx / (hbar * hbar) * V_phys;
        }
    }

    // Initial condition: Ground state of 2D harmonic oscillator
    // ψ₀(x,y) = (mω/πℏ)^(1/2) * exp(-mω((x-x₀)² + (y-y₀)²)/(2ℏ))
    double sigma = sqrt(hbar / (m * omega));
    double norm_factor = 1.0 / (sigma * sqrt(M_PI));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double r2 = (xx(i) - x0) * (xx(i) - x0) + (yy(j) - x0) * (yy(j) - x0);
            psi(i, j) = norm_factor * exp(-r2 / (2.0 * sigma * sigma));
        }
    }

    // Normalize initial state
    normalize(psi, dx);

    // Set up output files
    ofstream outfile("wavefunction.dat");
    ofstream diagfile("diagnostics.dat");

    // Calculate and output initial energy
    double E_kin = kinetic_energy(psi, hbar, m, dx);
    double E_pot = potential_energy(psi, V_original, dx, m, hbar);
    double E_total = E_kin + E_pot;
    double E_initial = E_total;
    double E_theory = hbar * omega;  // Ground state: E = ℏω(nx+1/2) + ℏω(ny+1/2) = ℏω for nx=ny=0

    // Print setup info
    printf("═══════════════════════════════════════════\n");
    printf("  2D Schrödinger Equation Solver\n");
    printf("═══════════════════════════════════════════\n");
    printf("Grid parameters:\n");
    printf("  n   = %d × %d\n", n, n);
    printf("  dx  = %.6f\n", dx);
    printf("  dt  = %.6f\n", dt);
    printf("───────────────────────────────────────────\n");
    printf("Physical parameters:\n");
    printf("  ℏ   = %.2f\n", hbar);
    printf("  m   = %.2f\n", m);
    printf("  ω   = %.2f\n", omega);
    printf("───────────────────────────────────────────\n");
    printf("Initial state: Ground state\n");
    printf("  E_kinetic    = %.6f\n", E_kin);
    printf("  E_potential  = %.6f\n", E_pot);
    printf("  E_total      = %.6f\n", E_total);
    printf("  E_theory     = %.6f (ground state)\n", E_theory);
    printf("═══════════════════════════════════════════\n\n");

    // Write header for diagnostics
    diagfile << "# time  probability  energy  E_kinetic  E_potential  rel_error\n";

    psi_old = psi;

    // Time evolution loop
    int output_interval = 10;  // Output every N steps

    for (int k = 0; k < nt; k++) {
        double t = k * dt;

        // Implicit solver with Gauss-Seidel iteration
        psi_new = psi;  // Initial guess
        int iter;

        for (iter = 0; iter < max_iter; iter++) {
            double max_diff = 0.0;

            // Gauss-Seidel: update in-place, using latest values
            for (int i = 1; i < n-1; i++) {
                for (int j = 1; j < n-1; j++) {
                    complex<double> old_val = psi_new(i, j);

                    // Use updated values for i-1, j-1 (Gauss-Seidel)
                    complex<double> neighbor_sum = psi_new(i+1, j) + psi_new(i-1, j) +
                                                    psi_new(i, j+1) + psi_new(i, j-1);

                    psi_new(i, j) = (psi_old(i, j) + alpha * neighbor_sum) /
                                     (1.0 + (4.0 + V(i, j)) * alpha);

                    double diff = abs(psi_new(i, j) - old_val);
                    if (diff > max_diff) max_diff = diff;
                }
            }

            // Check convergence
            if (max_diff < tol) {
                break;
            }
        }

        psi = psi_new;
        psi_old = psi;

        // Normalize
        normalize(psi, dx);

        // Calculate diagnostics
        if (k % output_interval == 0) {
            double prob = accu(abs(psi) % abs(psi)) * dx * dx;
            E_kin = kinetic_energy(psi, hbar, m, dx);
            E_pot = potential_energy(psi, V_original, dx, m, hbar);
            E_total = E_kin + E_pot;
            double rel_error = abs(E_total - E_initial) / abs(E_initial);

            // Write diagnostics
            diagfile << t << "  " << prob << "  " << E_total << "  "
                     << E_kin << "  " << E_pot << "  " << rel_error << "\n";

            // Write wavefunction
            outfile << "# Time = " << t << "\n";
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    outfile << xx(i) << " " << yy(j) << " " << abs(psi(i,j)) * abs(psi(i,j)) << "\n";
                }
                outfile << "\n";
            }
            outfile << "\n";

            // Print progress
            printf("Progress: %6.2f%% | t = %.4f | Iter: %3d | E = %.6f | ΔE/E = %.2e\r",
                   100.0 * k / nt, t, iter, E_total, rel_error);
            fflush(stdout);
        }
    }

    printf("\n\n");

    // Final diagnostics
    double prob_final = accu(abs(psi) % abs(psi)) * dx * dx;
    E_kin = kinetic_energy(psi, hbar, m, dx);
    E_pot = potential_energy(psi, V_original, dx, m, hbar);
    E_total = E_kin + E_pot;

    printf("═══════════════════════════════════════════\n");
    printf("  Final Results\n");
    printf("═══════════════════════════════════════════\n");
    printf("  Final time       = %.4f\n", nt * dt);
    printf("  Total probability = %.8f (should be 1.0)\n", prob_final);
    printf("  Final energy     = %.6f\n", E_total);
    printf("  Energy drift     = %.6f\n", E_total - E_initial);
    printf("  Relative error   = %.2e\n", abs(E_total - E_initial) / abs(E_initial));
    printf("═══════════════════════════════════════════\n");

    outfile.close();
    diagfile.close();

    cout << "\nOutput files created:\n";
    cout << "  - wavefunction.dat (probability density over time)\n";
    cout << "  - diagnostics.dat (energy and probability conservation)\n";
    cout << "\nRun 'python3 visualize.py' to see the animation!\n";

    return 0;
}
