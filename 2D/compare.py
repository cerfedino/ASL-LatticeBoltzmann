import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib   
import sys
import argparse


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'



def valid_directory(path):
    if os.path.isdir(path):
        return path
    else:
        raise argparse.ArgumentTypeError(f"The path '{path}' does not exist or is not a directory.")

parser = argparse.ArgumentParser(description='Compare two result sets for the 3D Lattice Boltzmann equation.')
parser.add_argument('--baseline', type=valid_directory, required=True, help='The folder containing the reference results')
parser.add_argument('--output', type=valid_directory, required=True, help='The folder containing the output results')

ARGS = parser.parse_args()

reference_folder = ARGS.baseline
output_folder = ARGS.output

# List files ending with .csv in ref folder
reference_output_files = [f for f in os.listdir(reference_folder) if f.endswith(".npy")]
output_files = [f for f in os.listdir(output_folder) if f.endswith(".npy")]


if len(reference_output_files) != len(output_files):
    print(f"{bcolors.FAIL}[X] The number of output files does not match the number of reference files{bcolors.ENDC}")
    exit(1)


TOLERANCE = 0.01


def is_equal(name):
    output = ""
    equal = True
    
    reference_file = os.path.join(f"{reference_folder}/{name}")
    output_file = os.path.join(f"{output_folder}/{name}")

    output+=f"\nComparing {output_file} with reference"
    ref_arr = np.load(reference_file)
    arr = np.load(output_file)

    if ref_arr.shape == arr.shape:
        output+=f"\n\t[+] Matrix shape matches {ref_arr.shape}"
    else:
        output+=f"\n\t{bcolors.FAIL}[ERROR] Matrix shape does not match! {arr.shape} should be {ref_arr.shape}{bcolors.ENDC}"
        equal = False

    if np.all(ref_arr == arr):
        output+=f"\n\t[+] Matrices are equal"
    else:
        diff = np.std(ref_arr - arr)
        if diff < TOLERANCE:
            output+=f"\n\t[+] Matrices are equal within tolerance {TOLERANCE}, diff: {diff}"
        else:
            output+=f"\n\t{bcolors.FAIL}[ERROR] Matrices are not equal within tolerance {TOLERANCE}{bcolors.ENDC}, diff: {diff}"
            equal = False

    #if not equal:
    print(output)
    return equal


error = False
for i, file in enumerate(reference_output_files):
    print(f"[-] Checking file {i+1}/{len(reference_output_files)}\r", end="")
    error = not is_equal(os.path.basename(file)) or error


print(f"[-] Checked the first {len(reference_output_files)} timesteps")    
if error:
    print(f"{bcolors.FAIL}[X] The implementation does not follow the baseline{bcolors.ENDC}")
    exit(1)
else:
    print(f"{bcolors.OKGREEN}[+] The implementation follows the baseline{bcolors.ENDC}")
    exit(0)