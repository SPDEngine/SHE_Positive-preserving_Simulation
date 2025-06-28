"""Module for generating space-time white noise increments.

This module provides functions to generate realizations of space-time
white noise increments for stochastic heat equation simulations.

Functions
---------
generate_noise(Nt, Nx, dt)
    Generate independent Gaussian increments scaled by sqrt(dt).
"""
import numpy as np


def generate_noise(Nt: int, Nx: int, dt: float) -> np.ndarray:
    """Generate space-time white noise increments.

    Parameters
    ----------
    Nt : int
        Number of time steps.
    Nx : int
        Number of spatial grid points including boundaries.
    dt : float
        Time step size.

    Returns
    -------
    ndarray
        Array of shape ``(Nt, Nx-1)`` with independent normal
        variables scaled by ``sqrt(dt)``.
    """
    return np.sqrt(dt) * np.random.randn(Nt, Nx - 1)
