
import argparse
import sys
from .solver import run_sup_dxdt
from .plot_3d_simulation import main as plot_main

def main():
    parser = argparse.ArgumentParser(description='Run SHE simulation and plot results.')
    parser.add_argument('--seed',    type=int,            default=93,          help='Random seed for the simulation.')
    parser.add_argument('--output',  type=str,            default='msSupDxDt', help='Output file name prefix.')
    parser.add_argument('--samples', type=int,            default=18,          help='Number of Monte Carlo samples.')
    parser.add_argument('--alphaN',  type=float,          default=1.0,         help='Noise strength parameter.')
    parser.add_argument('--no-plot', action='store_true',                      help='Suppress plotting the results.')

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    run_sup_dxdt(args.seed, args.output, args.samples, args.alphaN)

    if not args.no_plot:
        plot_main()

if __name__ == '__main__':
    main()
