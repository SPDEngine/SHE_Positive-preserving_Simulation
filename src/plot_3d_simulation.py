"""plot_3d_simulation.py

Visualise the msSupDxDt results in three dimensions.

This script loads ``msSupDxDt.npz`` containing ``dt_num``, ``h_num`` and
``maxLT`` arrays and produces a 3-D curve where the axes correspond to the
time-step sizes ``dt_num``, spatial mesh widths ``h_num`` and the square root
of the maximum mean-square error ``maxLT``.

Usage
-----
python plot_3d_simulation.py
"""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  # for 3-D plotting
from tqdm import tqdm


def main() -> None:
    """Load experiment data and display a 3-D plot."""
    data_path = Path("msSupDxDt.npz")
    if not data_path.exists():
        raise FileNotFoundError(f"{data_path} not found")

    data = np.load(data_path)
    print(data.files)
    for key in data.files:
        arr = data[key]
        print(f"{key}: shape={arr.shape}, dtype={arr.dtype}")
        print(arr if arr.size < 10 else arr[:10])  # Preview first 10 elements if large
        print("-" * 40)
    dt_num = data["dt_num"]
    h_num = data["h_num"]
    # maxLT = data["maxLT"]
    u_exLT = data["u_exLT"]

    # Create time and space grids based on u_exLT shape (space, time)
    num_spatial_points, num_time_steps = u_exLT.shape
    print(f"Number of time steps: {num_time_steps}, Computed steps: {0.5/dt_num}")
    print(f"Number of spatial points: {num_spatial_points}, Computed points: {1/h_num}")

    time_grid = np.linspace(0, 0.5, num_time_steps)
    print(f"Time grid: {time_grid}")
    space_grid = np.linspace(0, 1, num_spatial_points)
    print(f"Space grid: {space_grid}")

    # Downsample grids for faster 3D plotting and animation
    max_pts = 100
    s_space = max(1, num_spatial_points // max_pts)
    s_time = max(1, num_time_steps // max_pts)
    space_sub = space_grid[::s_space]
    time_sub = time_grid[::s_time]
    u_sub = u_exLT[::s_space, ::s_time]

    # 3D Plot
    fig_3d = plt.figure()
    ax_3d = fig_3d.add_subplot(111, projection='3d')
    X_3d, T_3d = np.meshgrid(space_sub, time_sub, indexing='ij')
    surf = ax_3d.plot_surface(X_3d, T_3d, u_sub, cmap=cm.viridis, linewidth=0, antialiased=False)
    fig_3d.colorbar(surf, shrink=0.5, aspect=5)

    ax_3d.set_xlabel('Space')
    ax_3d.set_ylabel('Time')
    ax_3d.set_zlabel('u_exLT')
    ax_3d.set_title('3D Mesh Plot of u_exLT vs. Time and Space')
    plt.tight_layout()
    fig_3d.savefig('3d_simulation.png')
    print("3D plot saved as 3d_simulation.png")

    # 2D Animation
    fig_2d, ax_2d = plt.subplots()
    line, = ax_2d.plot(space_grid, u_exLT[:, 0], color='blue')

    ax_2d.set_xlabel('Space')
    ax_2d.set_ylabel('u_exLT')
    ax_2d.set_title(f'2D Plot of u_exLT at Time = {time_grid[0]:.3f}')
    ax_2d.set_ylim(u_exLT.min(), u_exLT.max()) # Set fixed y-limits
    plt.tight_layout()
    plt.show()

    # Create the animation
    def animate(i):
        line.set_ydata(u_exLT[:, i])
        ax_2d.set_title(f'2D Plot of u_exLT at Time = {time_grid[i]:.3f}')
        return line,

    # Subsample frames for faster animation
    frame_idx = np.arange(0, num_time_steps, s_time)
    ani = animation.FuncAnimation(fig_2d, animate, frames=frame_idx, interval=50, blit=True)

    # Save the animation
    print("Saving animation... This may take a while.")
    try:
        with tqdm(total=num_time_steps, desc="Saving animation") as pbar:
            def update_progress(i, n):
                pbar.update(1)

            ani.save('2d_simulation.mp4', writer='ffmpeg', fps=10, progress_callback=update_progress)

        print()
        print("Animation saved as 2d_simulation.mp4")
    except ValueError as e:
        print(f"Error saving animation: {e}")
        print("Please ensure ffmpeg is installed and accessible in your PATH.")


if __name__ == "__main__":
    main()
