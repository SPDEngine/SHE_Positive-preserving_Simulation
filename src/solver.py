"""Module for running the 'msSupDxDt' stochastic heat equation experiments in Python.

This module provides:
- run_sup_dxdt: Function to execute the splitting scheme experiment, compute
  mean-square errors, and save results to a .npz file and a PNG plot.

The simulation is performed on the spatial domain (0,1) with homogeneous Dirichlet boundary conditions.
"""
import argparse
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from .laplacian import make_laplacian
from .noise import generate_noise
from tqdm.auto import trange


def run_sup_dxdt(seed: int, out_name: str, M: int = 4, alphaN: float = 1.0, T: float = 0.5, discretization_level: int = 10) -> None:
    """Run the splitting scheme experiment for the stochastic heat equation.

    Parameters
    ----------
    seed : int
        Random seed for reproducibility.
    out_name : str
        Prefix for output files ('.npz' and '.png' will be appended).
    M : int, optional
        Number of Monte Carlo samples (default is 4).
    alphaN : float, optional
        Noise strength parameter (default is 1.0).
    T : float, optional
        End time of the simulation (default is 0.5).
    discretization_level : int, optional
        The level of spatial discretization (default is 10). The spatial step size `h` is calculated as `2**(-discretization_level)`.

    Returns
    -------
    None

    Side Effects
    ------------
    Saves a '.npz' file containing 'dt_num', 'h_num', and 'maxLT', and a
    log-log error plot as a PNG image.

    The simulation is performed on the spatial domain (0,1) with homogeneous Dirichlet boundary conditions.
    """

    print(f"run_sup_dxdt: seed={seed}, output={out_name}, M={M}")
    rng = np.random.default_rng(seed)

    
    x_0, x_end = 0.0, 1.0
    t_0 = 0.5 * 0

    u_0 = lambda x: np.cos(np.pi * (x - 0.5))

    h = 2.0 ** -discretization_level
    dt = h ** 2

    Nt = int((T - t_0) / dt)
    x = np.arange(x_0 + h, x_end, h)
    Nx = len(x) + 1

    dt_num = np.array([dt])
    h_num = np.array([h])
    Nt_num = np.array([Nt])
    Nx_num = np.array([Nx])

    # For compatibility with existing code that expects _ex variables
    h_ex = h
    dt_ex = dt
    Nt_ex = Nt
    x_ex = x
    Nx_ex = Nx

    # Precompute fine-level operator and noise factors
    D_ex = make_laplacian(Nx_ex)
    expD_ex = expm(dt_ex * D_ex)
    f_factor = alphaN
    sqrt_factor = np.sqrt(Nx_ex - 1) * f_factor
    const_factor = -0.5 * (f_factor ** 2) * (Nx_ex - 1) * dt_ex

    # storage for all samples of u_exLT
    all_u_exLT = np.empty((M, Nx_ex - 1, Nt_ex + 1))

    for m in trange(M, desc="samples"):
        # generate fine noise increments
        W = rng.standard_normal((Nt_ex, Nx_ex - 1)) * np.sqrt(dt_ex)
        R = np.exp(const_factor + sqrt_factor * W)

        # Reference fine solution
        u_exLT = np.empty((Nx_ex - 1, Nt_ex + 1))
        u_exLT[:, 0] = u_0(x_ex)
        # Reference fine solution time-stepping with progress bar
        for t in trange(Nt_ex, desc="reference", leave=False):
            u_exLT[:, t + 1] = expD_ex @ (u_exLT[:, t] * R[t, :])
        all_u_exLT[m, :, :] = u_exLT

    out = Path(f"{out_name}.npz")
    np.savez(out, dt_num=dt_num, h_num=h_num, u_exLT=all_u_exLT)
    print(f"Results stored in {out}")

