from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from matplotlib.collections import LineCollection
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  # for 3-D plotting
from tqdm import tqdm

# Set up Matplotlib for LaTeX rendering
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Computer Modern Roman']


def main(fps: int = 10, output_prefix: str = "msSupDxDt", M: int = 4) -> None:
    """Load experiment data and display a 3-D plot."""
    data_path = Path(f"{output_prefix}.npz")
    if not data_path.exists():
        raise FileNotFoundError(f"{data_path} not found")

    data = np.load(data_path)
    all_u_exLT = data["u_exLT"]
    dt_num = data["dt_num"]
    h_num = data["h_num"]

    for m_idx in range(M):
        u_exLT = all_u_exLT[m_idx, :, :]

        # Create time and space grids based on u_exLT shape (space, time)
        num_spatial_points, num_time_steps = u_exLT.shape
        print(f"Number of time steps: {num_time_steps}, Computed steps: {0.5/dt_num}")
        print(f"Number of spatial points: {num_spatial_points}, Computed points: {1/h_num}")

        time_grid = np.linspace(0, 0.5, num_time_steps)
        # print(f"Time grid: {time_grid}")
        space_grid = np.linspace(0, 1, num_spatial_points)
        # print(f"Space grid: {space_grid}")

        # Use full resolution for plotting
        s_space = 1
        s_time = 1
        space_sub = space_grid
        time_sub = time_grid
        u_sub = u_exLT

        # 3D Plot
        fig_3d = plt.figure()
        ax_3d = fig_3d.add_subplot(111, projection='3d')
        X_3d, T_3d = np.meshgrid(space_sub, time_sub, indexing='ij')
        surf = ax_3d.plot_surface(X_3d, T_3d, u_sub, cmap=cm.viridis, linewidth=0, antialiased=False)
        # fig_3d.colorbar(surf, shrink=0.5, aspect=5)

        ax_3d.set_xlabel('Space')
        ax_3d.set_ylabel('Time')
        ax_3d.set_zlabel(r'$u(t,x)$')
        ax_3d.set_title(r'3D Mesh Plot of $u(t,x)$ vs. Time and Space')
        plt.tight_layout()
        fig_3d.savefig(f'{output_prefix}_{M}-{m_idx+1}_3d.png')
        print(f"3D plot saved as {output_prefix}_{M}-{m_idx+1}_3d.png")

        # 2D Animation
        fig_2d, ax_2d = plt.subplots()
        # Set up colormap and normalization for the 2D plot
        norm = colors.Normalize(vmin=u_exLT.min(), vmax=u_exLT.max())
        cmap = cm.viridis
        # Initial plot with LineCollection for gradient line
        points = np.array([space_grid, u_exLT[:, 0]]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        line_collection = LineCollection(segments, cmap=cmap, norm=norm)
        line_collection.set_array(u_exLT[:, 0])
        line_collection.set_linewidth(2) # Adjust linewidth as desired
        ax_2d.add_collection(line_collection)
        # Add a thin black line on top for better contour visibility
        line_contour, = ax_2d.plot(space_grid, u_exLT[:, 0], color='black', linewidth=0.5)

        ax_2d.set_xlabel('Space')
        ax_2d.set_ylabel(r'$u(t,x)$')
        ax_2d.set_title(f'$u(t,x)$ at $t$= {time_grid[0]:.3f}')
        ax_2d.set_ylim(u_exLT.min(), u_exLT.max()) # Set fixed y-limits
        plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax_2d)
        plt.tight_layout()

        # Create the animation
        def animate(i):
            # Update LineCollection with new data and colors
            points = np.array([space_grid, u_exLT[:, i]]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            line_collection.set_segments(segments)
            line_collection.set_array(u_exLT[:, i])
            line_contour.set_ydata(u_exLT[:, i])
            ax_2d.set_title(f'$u(t,x)$ at $t$ = {time_grid[i]:.3f}')
            return line_collection, line_contour,

        # Subsample frames for faster animation
        frame_idx = np.arange(0, num_time_steps, s_time)
        ani = animation.FuncAnimation(fig_2d, animate, frames=frame_idx, interval=50, blit=False)

        # Save the animation
        print("Saving animation... This may take a while.")
        try:
            with tqdm(total=len(frame_idx), desc="Saving animation") as pbar:
                def update_progress(i, n):
                    pbar.update(1)

                ani.save(f'{output_prefix}_{M}-{m_idx+1}_2d.mp4', writer='ffmpeg', fps=fps, progress_callback=update_progress)

            print()
            print(f"Animation saved as {output_prefix}_{M}-{m_idx+1}_2d.mp4")
        except ValueError as e:
            print(f"Error saving animation: {e}")
            print("Please ensure ffmpeg is installed and accessible in your PATH.")
        plt.close(fig_3d)
        plt.close(fig_2d)
