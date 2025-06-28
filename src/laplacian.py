"""Module for constructing finite-difference Laplacian matrices with homogeneous
Dirichlet boundary conditions on a unit interval.

This module assembles a second-order central difference Laplacian matrix,
scaled appropriately for a unit domain discretization.

Functions
---------
make_laplacian(Nx)
    Assemble a scaled (Nx-1)x(Nx-1) finite-difference Laplacian matrix.
"""
import numpy as np


def make_laplacian(Nx: int) -> np.ndarray:
    """Assemble finite difference Laplacian matrix.

    Parameters
    ----------
    Nx : int
        Number of grid points including the boundaries.

    Returns
    -------
    ndarray
        The ``(Nx-1) x (Nx-1)`` Laplacian matrix for Dirichlet conditions
        scaled by ``Nx**2``.
    """
    D = np.zeros((Nx - 1, Nx - 1))
    for j in range(Nx - 2):
        D[j, j] = -2.0
        D[j, j + 1] = 1.0
        D[j + 1, j] = 1.0
    D[Nx - 2, Nx - 2] = -2.0
    D *= Nx ** 2
    return D
