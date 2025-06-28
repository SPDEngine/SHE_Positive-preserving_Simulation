
import argparse
import sys
import yaml
from .solver import run_sup_dxdt
from .plot_3d_simulation import main as plot_main

def main():
    parser = argparse.ArgumentParser(description='Run SHE simulation and plot results.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--config', type=str, help='Path to a YAML configuration file.')
    parser.add_argument('--seed',                 type=int,            default=93,          help='Random seed for the simulation.')
    parser.add_argument('--output',               type=str,            default='msSupDxDt', help='Output file name prefix.')
    parser.add_argument('--samples',              type=int,            default=4,           help='Number of Monte Carlo samples.')
    parser.add_argument('--alphaN',               type=float,          default=1.0,         help='Noise strength parameter.')
    parser.add_argument('--T',                    type=float,          default=0.5,         help='End time of the simulation.')
    parser.add_argument('--discretization_level', type=int,            default=10,          help='Level of spatial discretization.')
    parser.add_argument('--fps',                  type=int,            default=10,          help='Frames per second for the animation.')
    parser.add_argument('--max_animation_frames', type=int,            default=120,         help='Maximum number of frames for the animation.')
    parser.add_argument('--no-plot',              action='store_true', help='Suppress plotting the results.')

    # Parse arguments from the command line first
    # This allows --config to be processed before other arguments are finalized
    args, remaining_argv = parser.parse_known_args()

    # If a config file is specified, load it and update defaults
    if args.config:
        with open(args.config, 'r') as f:
            config = yaml.safe_load(f)
        # Update parser defaults with values from the config file
        parser.set_defaults(**config)

    # Re-parse all arguments, now with updated defaults from config file
    # Command-line arguments will override config file values
    args = parser.parse_args(remaining_argv)

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    output_file_path = Path(f"{args.output}.npz")

    if output_file_path.exists():
        print(f"Skipping simulation: {output_file_path} already exists.")
    else:
        run_sup_dxdt(args.seed, args.output, args.samples, args.alphaN, args.T, args.discretization_level)

    if not args.no_plot:
        plot_main(args.fps, args.output, args.samples, args.max_animation_frames)

if __name__ == '__main__':
    main()
