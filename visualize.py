#!/usr/bin/env python3
"""
Visualization script for 2D Schrödinger equation simulation
Creates animated plots of the wavefunction probability density
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from mpl_toolkits.mplot3d import Axes3D
import sys

def read_wavefunction_data(filename):
    """Read wavefunction data from file"""
    frames = []
    times = []
    current_frame = []
    current_time = 0

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('# Time'):
                # Save previous frame if it exists
                if current_frame:
                    frames.append(np.array(current_frame))
                current_frame = []
                # Extract time
                current_time = float(line.split('=')[1])
                times.append(current_time)
            elif line and not line.startswith('#'):
                parts = line.split()
                if len(parts) == 3:
                    x, y, prob = map(float, parts)
                    current_frame.append([x, y, prob])

    # Add last frame
    if current_frame:
        frames.append(np.array(current_frame))

    return frames, times

def read_diagnostics(filename):
    """Read energy and probability diagnostics"""
    data = np.loadtxt(filename)
    return {
        'time': data[:, 0],
        'probability': data[:, 1],
        'energy': data[:, 2],
        'E_kinetic': data[:, 3],
        'E_potential': data[:, 4],
        'rel_error': data[:, 5]
    }

def create_visualization():
    """Create comprehensive visualization of the simulation"""

    print("Loading data...")
    frames, times = read_wavefunction_data('wavefunction.dat')
    diag = read_diagnostics('diagnostics.dat')

    n = int(np.sqrt(len(frames[0])))
    print(f"Grid size: {n}×{n}")
    print(f"Number of frames: {len(frames)}")

    # Create figure with subplots
    fig = plt.figure(figsize=(16, 10))

    # 1. 3D surface plot of probability density
    ax1 = fig.add_subplot(2, 3, 1, projection='3d')

    # 2. 2D heatmap
    ax2 = fig.add_subplot(2, 3, 2)

    # 3. Cross-section through center
    ax3 = fig.add_subplot(2, 3, 3)

    # 4. Energy conservation
    ax4 = fig.add_subplot(2, 3, 4)
    ax4.plot(diag['time'], diag['energy'], 'b-', linewidth=2, label='Total Energy')
    ax4.plot(diag['time'], diag['E_kinetic'], 'r--', linewidth=1.5, label='Kinetic')
    ax4.plot(diag['time'], diag['E_potential'], 'g--', linewidth=1.5, label='Potential')
    ax4.set_xlabel('Time', fontsize=12)
    ax4.set_ylabel('Energy', fontsize=12)
    ax4.set_title('Energy Conservation', fontsize=14, fontweight='bold')
    ax4.legend()
    ax4.grid(True, alpha=0.3)

    # 5. Probability conservation
    ax5 = fig.add_subplot(2, 3, 5)
    ax5.plot(diag['time'], diag['probability'], 'b-', linewidth=2)
    ax5.axhline(y=1.0, color='r', linestyle='--', linewidth=1, label='Expected (1.0)')
    ax5.set_xlabel('Time', fontsize=12)
    ax5.set_ylabel('Total Probability', fontsize=12)
    ax5.set_title('Probability Conservation', fontsize=14, fontweight='bold')
    ax5.legend()
    ax5.grid(True, alpha=0.3)
    ax5.set_ylim([0.999, 1.001])

    # 6. Relative energy error
    ax6 = fig.add_subplot(2, 3, 6)
    ax6.semilogy(diag['time'], diag['rel_error'], 'b-', linewidth=2)
    ax6.set_xlabel('Time', fontsize=12)
    ax6.set_ylabel('Relative Energy Error', fontsize=12)
    ax6.set_title('Energy Error (log scale)', fontsize=14, fontweight='bold')
    ax6.grid(True, alpha=0.3, which='both')

    # Prepare first frame
    frame_data = frames[0]
    X = frame_data[:, 0].reshape(n, n)
    Y = frame_data[:, 1].reshape(n, n)
    Z = frame_data[:, 2].reshape(n, n)

    # 3D surface
    surf = ax1.plot_surface(X, Y, Z, cmap='viridis', edgecolor='none', alpha=0.9)
    ax1.set_xlabel('x', fontsize=10)
    ax1.set_ylabel('y', fontsize=10)
    ax1.set_zlabel('|ψ|²', fontsize=10)
    ax1.set_title('Probability Density (3D)', fontsize=14, fontweight='bold')
    ax1.view_init(elev=30, azim=45)

    # Set fixed z-axis limits
    max_prob = max([frame[:, 2].max() for frame in frames])
    ax1.set_zlim([0, max_prob * 1.1])

    # 2D heatmap
    im = ax2.imshow(Z, extent=[0, 1, 0, 1], origin='lower', cmap='hot',
                    vmin=0, vmax=max_prob * 1.1, aspect='auto')
    ax2.set_xlabel('x', fontsize=12)
    ax2.set_ylabel('y', fontsize=12)
    ax2.set_title('Probability Density (2D)', fontsize=14, fontweight='bold')
    plt.colorbar(im, ax=ax2, label='|ψ|²')

    # Cross-section
    center_idx = n // 2
    line_x, = ax3.plot(X[center_idx, :], Z[center_idx, :], 'b-', linewidth=2, label='y=0.5')
    line_y, = ax3.plot(Y[:, center_idx], Z[:, center_idx], 'r--', linewidth=2, label='x=0.5')
    ax3.set_xlabel('Position', fontsize=12)
    ax3.set_ylabel('|ψ|²', fontsize=12)
    ax3.set_title('Cross-sections', fontsize=14, fontweight='bold')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    ax3.set_ylim([0, max_prob * 1.1])

    # Time indicator on energy plot
    time_line = ax4.axvline(x=times[0], color='k', linestyle=':', linewidth=1.5, alpha=0.5)

    # Add overall title
    title = fig.suptitle(f'2D Schrödinger Equation: Ground State Evolution (t = {times[0]:.4f})',
                         fontsize=16, fontweight='bold')

    plt.tight_layout(rect=[0, 0, 1, 0.96])

    def update(frame_idx):
        """Update function for animation"""
        frame_data = frames[frame_idx]
        t = times[frame_idx]

        # Reshape data
        Z = frame_data[:, 2].reshape(n, n)

        # Update 3D surface
        ax1.clear()
        ax1.plot_surface(X, Y, Z, cmap='viridis', edgecolor='none', alpha=0.9)
        ax1.set_xlabel('x', fontsize=10)
        ax1.set_ylabel('y', fontsize=10)
        ax1.set_zlabel('|ψ|²', fontsize=10)
        ax1.set_title('Probability Density (3D)', fontsize=14, fontweight='bold')
        ax1.set_zlim([0, max_prob * 1.1])
        ax1.view_init(elev=30, azim=45)

        # Update heatmap
        im.set_data(Z)

        # Update cross-sections
        line_x.set_ydata(Z[center_idx, :])
        line_y.set_ydata(Z[:, center_idx])

        # Update time indicator
        time_line.set_xdata([t, t])

        # Update title
        title.set_text(f'2D Schrödinger Equation: Ground State Evolution (t = {t:.4f})')

        return surf, im, line_x, line_y, time_line, title

    # Create animation
    print("Creating animation...")
    anim = FuncAnimation(fig, update, frames=len(frames), interval=100, blit=False)

    # Save animation
    print("Saving animation as 'schrodinger_2d.gif'...")
    writer = PillowWriter(fps=10)
    anim.save('schrodinger_2d.gif', writer=writer, dpi=80)
    print("Animation saved!")

    # Save final frame as static image
    print("Saving final frame as 'schrodinger_2d_final.png'...")
    update(len(frames) - 1)
    plt.savefig('schrodinger_2d_final.png', dpi=150, bbox_inches='tight')
    print("Final frame saved!")

    # Show interactive plot
    print("\nShowing interactive plot (close window to exit)...")
    plt.show()

def create_summary_plots():
    """Create summary diagnostic plots"""
    print("Creating summary diagnostic plots...")

    diag = read_diagnostics('diagnostics.dat')

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Energy components over time
    ax = axes[0, 0]
    ax.plot(diag['time'], diag['energy'], 'b-', linewidth=2, label='Total')
    ax.plot(diag['time'], diag['E_kinetic'], 'r--', linewidth=1.5, label='Kinetic')
    ax.plot(diag['time'], diag['E_potential'], 'g--', linewidth=1.5, label='Potential')
    ax.set_xlabel('Time', fontsize=12)
    ax.set_ylabel('Energy', fontsize=12)
    ax.set_title('Energy Components', fontsize=14, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

    # Energy conservation (zoomed)
    ax = axes[0, 1]
    E_mean = np.mean(diag['energy'])
    ax.plot(diag['time'], diag['energy'] - E_mean, 'b-', linewidth=2)
    ax.set_xlabel('Time', fontsize=12)
    ax.set_ylabel('Energy - Mean(Energy)', fontsize=12)
    ax.set_title('Energy Drift from Mean', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.ticklabel_format(axis='y', style='scientific', scilimits=(-3, 3))

    # Probability conservation
    ax = axes[1, 0]
    ax.plot(diag['time'], diag['probability'], 'b-', linewidth=2)
    ax.axhline(y=1.0, color='r', linestyle='--', linewidth=1, label='Expected')
    ax.set_xlabel('Time', fontsize=12)
    ax.set_ylabel('Total Probability ∫|ψ|²dxdy', fontsize=12)
    ax.set_title('Probability Conservation', fontsize=14, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

    # Relative energy error (log scale)
    ax = axes[1, 1]
    ax.semilogy(diag['time'], diag['rel_error'], 'b-', linewidth=2)
    ax.set_xlabel('Time', fontsize=12)
    ax.set_ylabel('|ΔE| / |E₀|', fontsize=12)
    ax.set_title('Relative Energy Error', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3, which='both')

    plt.suptitle('2D Schrödinger Equation: Diagnostics', fontsize=16, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    plt.savefig('diagnostics.png', dpi=150, bbox_inches='tight')
    print("Diagnostic plots saved as 'diagnostics.png'")
    plt.show()

if __name__ == '__main__':
    try:
        print("=" * 60)
        print("  2D Schrödinger Equation Visualization")
        print("=" * 60)

        # First create diagnostic plots
        create_summary_plots()

        # Then create main animation
        create_visualization()

        print("\n" + "=" * 60)
        print("Visualization complete!")
        print("Generated files:")
        print("  - schrodinger_2d.gif (animation)")
        print("  - schrodinger_2d_final.png (final state)")
        print("  - diagnostics.png (energy and probability plots)")
        print("=" * 60)

    except FileNotFoundError as e:
        print(f"Error: Could not find data file. Make sure to run the simulation first!")
        print(f"Details: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error during visualization: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
