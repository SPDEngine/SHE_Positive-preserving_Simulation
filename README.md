# Simulation of SHE with positivity-preserving splitting scheme

The goal of this repository is to simulate the Stochastic Heat Equation (SHE)
using a positivity-preserving splitting scheme, with outputs presented as 3D
plots and MP4 animations.

## Reference

These MATLAB scripts support the results published in

Bréhier, Charles-Edouard; Cohen, David; Ulander, Johan.
"Analysis of a positivity-preserving splitting scheme for some semilinear
stochastic heat equations." *ESAIM: Mathematical Modelling and Numerical
Analysis* 58(4):1317–1346, 2024. DOI:
<https://doi.org/10.1051/m2an/2024032>.

Some functions were translated from the original MATLAB code available at [Zenodo](https://zenodo.org/records/10300733) to Python.

Further background and a brief podcast discussion can be found in the
official reference entry at
<https://spdes-bib.readthedocs.io/en/latest/bib_entries/brehier.cohen.ea:24:analysis.html>.

## Usage

To install and run this project, follow these steps:

1.  **Clone the repository:**

    ```bash
    git clone https://github.com/lzc0090/SHE_Positive-preserving_Simulation.git
    cd SHE_Positive-preserving_Simulation
    ```

2.  **Install the package:**

    ```bash
    pip install .
    ```

3.  **Run the simulation:**

    To run the simulation with default parameters and generate plots, use the following command:

    ```bash
    she-sim
    ```

    You can also customize the simulation with the following options:

    ```bash
    she-sim --seed <seed> --samples <samples> --output <output_prefix> --alphaN <alphaN> --no-plot
    ```

    -   `--seed`: Random seed for the simulation (default: 93).
    -   `--samples`: Number of Monte Carlo samples (default: 18).
    -   `--output`: Prefix for the output file names (default: `msSupDxDt`).
    -   `--alphaN`: Noise strength parameter (default: 1.0).
    -   `--no-plot`: Suppress the generation of plots.

## License

[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)

The code is translated from Matlab to Python by Le Chen. Le Chen would like to thank David Cohen for his help in understanding the original MATLAB code.

