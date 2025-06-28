"""Module for running the 'msSupDxDt' stochastic heat equation experiments in Python.

This module provides:
- run_sup_dxdt: Function to execute the splitting scheme experiment, compute
  mean-square errors, and save results to a .npz file and a PNG plot.
- main: Command-line interface for running experiments.

Usage
-----
python solver.py --seed 93 --output msSupDxDt --samples 18
"""
import argparse
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from .laplacian import make_laplacian
from .noise import generate_noise
from tqdm.auto import trange


def run_sup_dxdt(seed: int, out_name: str, M: int = 18) -> None:
    """Run the msSupDxDt experiment replicating the MATLAB implementation.

    Parameters
    ----------
    seed : int
        Random seed for reproducibility.
    out_name : str
        Prefix for output files ('.npz' and '.png' will be appended).
    M : int, optional
        Number of Monte Carlo samples (default is 18).

    Returns
    -------
    None

    Side Effects
    ------------
    Saves a '.npz' file containing 'dt_num', 'h_num', and 'maxLT', and a
    log-log error plot as a PNG image.
    """

    print(f"run_sup_dxdt: seed={seed}, output={out_name}, M={M}")
    rng = np.random.default_rng(seed)

    alphaN = 1.0
    x_0, x_end = 0.0, 1.0
    t_0, t_end = 0.5 * 0, 0.5  # t_0=0; t_end=.5 in the scripts

    u_0 = lambda x: np.cos(np.pi * (x - 0.5))

    g = lambda u: alphaN * u
    f = lambda u: alphaN

    hh = 2.0 ** np.arange(-10, -10, -1)
    ddt = hh ** 2
    h_ex = 2.0 ** -10
    dt_ex = h_ex ** 2

    NNt = ((t_end - t_0) / ddt).astype(int)
    Nt_ex = int((t_end - t_0) / dt_ex)

    Nx = []
    for h in hh:
        xx = np.arange(x_0 + h, x_end, h)
        Nx.append(len(xx) + 1)
    x_ex = np.arange(x_0 + h_ex, x_end, h_ex)
    Nx_ex = len(x_ex) + 1

    dt_num = np.concatenate([ddt, [dt_ex]])
    h_num = np.concatenate([hh, [h_ex]])
    Nt_num = np.concatenate([NNt, [Nt_ex]])
    Nx_num = np.concatenate([Nx, [Nx_ex]])

    # storage for mean-square errors
    errorLT = [[None] * M for _ in range(len(dt_num))]

    # Precompute fine-level operator and noise factors
    D_ex = make_laplacian(Nx_ex)
    expD_ex = expm(dt_ex * D_ex)
    f_factor = alphaN
    sqrt_factor = np.sqrt(Nx_ex - 1) * f_factor
    const_factor = -0.5 * (f_factor ** 2) * (Nx_ex - 1) * dt_ex

    # Precompute coarse-level operators and grouping factors
    expD_list = []
    k_list = []
    kx_list = []
    levels = len(dt_num)
    for l in range(levels):
        dt = dt_num[l]
        h = h_num[l]
        k = int(round(dt / dt_ex))
        kx = int(round(h / h_ex))
        k_list.append(k)
        kx_list.append(kx)
        Nx_l = int(Nx_num[l])
        D_l = make_laplacian(Nx_l)
        expD_list.append(expm(dt * D_l))

    for m in trange(M, desc="samples"):
        # generate fine noise increments and precompute multiplicative factors
        W = rng.standard_normal((Nt_ex, Nx_ex - 1)) * np.sqrt(dt_ex)
        R = np.exp(const_factor + sqrt_factor * W)

        # Reference fine solution
        u_exLT = np.empty((Nx_ex - 1, Nt_ex + 1))
        u_exLT[:, 0] = u_0(x_ex)
        # Reference fine solution time-stepping with progress bar
        for t in trange(Nt_ex, desc="reference", leave=False):
            u_exLT[:, t + 1] = expD_ex @ (u_exLT[:, t] * R[t, :])

        # for l in reversed(range(len(dt_num))):
        #     Nt = int(Nt_num[l])
        #     Nx_l = int(Nx_num[l])
        #     dt = dt_num[l]
        #     h = h_num[l]
        #     k = k_list[l]
        #     kx = kx_list[l]
        #     expD = expD_list[l]

        #     xx = np.arange(x_0 + h, x_end, h)


        #     err = np.zeros((Nx_l - 1, Nt + 1))

        #     # Coarse noise coupling: vectorized block sums
        #     if k == 1 and kx == 1:
        #         V = W
        #     else:
        #         # reshape W of shape (Nt_ex, Nx_ex-1) to (Nt, k, Nx_l-1, kx)
        #         V = W.reshape(Nt, k, Nx_l - 1, kx).sum(axis=(1, 3)) / np.sqrt(kx)

        #     uLT = u_0(xx)
        #     for t in trange(Nt, desc="scheme", leave=False):
        #         tempo = uLT * np.exp(
        #             (-f(uLT) ** 2 / 2 * (Nx_l - 1)) * dt
        #             + np.sqrt(Nx_l - 1) * f(uLT) * V[t, :]
        #         )
        #         uLT = expD @ tempo
        #         for j in range(Nx_l - 1):
        #             err[j, t + 1] = (
        #                 u_exLT[kx * j, k * t + 1] - uLT[j]
        #             ) ** 2

        #     errorLT[l][m] = err

    # maxLT = np.zeros(len(dt_num))
    # for l in range(len(dt_num)):
        # tmp = sum(errorLT[l]) / M
        # maxLT[l] = tmp.max()

    out = Path(f"{out_name}.npz")
    # np.savez(out, dt_num=dt_num, h_num=h_num, maxLT=maxLT, u_exLT=u_exLT)
    np.savez(out, dt_num=dt_num, h_num=h_num, u_exLT=u_exLT)
    print(f"Results stored in {out}")

    # # simple plot
    # import matplotlib.pyplot as plt

    # plt.loglog(dt_num, np.sqrt(maxLT), "ms-", label="Error LT")
    # plt.loglog(dt_num, dt_num ** 0.25, "r--", label="Slope 1/4")
    # plt.xlabel("$\\Delta t$")
    # plt.ylabel("Error")
    # plt.legend(loc="lower right")
    # plt.title(f"$M_s={M}$")

    # fig_path = Path(f"{out_name}.png")
    # plt.savefig(fig_path, dpi=150)
    # plt.close()
    # print(f"Results stored in {out} and {fig_path}")



