from features import *
import argparse

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("--input_dir", default="GPSExample", type=str)
	parser.add_argument("--output_file", default="output.tsv", type=str)	
	parser.add_argument("--wtype", default='GLR', type=str)
	parser.add_argument("--center_rad", default=200, type=int)
	parser.add_argument("--interval", default=10, type=int)
	parser.add_argument("--accuracy", default=51.0, type=float)
	parser.add_argument("--n_reps", default=1, type=int)
	parser.add_argument("--min_pause_dur", default=300, type=int)
	parser.add_argument("--min_pause_dist", default=60, type=int)
	parser.add_argument("--r", default=None, type=float)
	parser.add_argument("--w", default=None, type=float)

	args = parser.parse_args()
	input_dir = args.input_dir
	output_file = args.output_file
	wtype = args.wtype
	center_rad = args.center_rad
	interval = args.interval
	acc_threshold = args.accuracy
	n_reps = args.n_reps
	min_pause_dur = args.min_pause_dur
	min_pause_dist = args.min_pause_dist

	print("Parameters: ")
	print("Input directory: ", input_dir)
	print("Output file: ", output_file)
	print("W type: ", wtype)
	print("Center radius: ", center_rad)
	print("Interval: ", interval)
	print("Accuracy threshold: ", acc_threshold)
	print("N reps: ", n_reps)
	print("Minimum pause duration: ", min_pause_dur)
	print("Minimum pause distance: ", min_pause_dist)

	run_barnett_features(input_dir, output_file, wtype="GLR", spread_pars=[10,1], tz="", center_rad=center_rad, interval=interval, acc_threshold=acc_threshold, n_reps=n_reps, min_pause_dur=min_pause_dur, min_pause_dist=min_pause_dist, r=None, w=None, tint_m=None, tint_k=None)

main()
