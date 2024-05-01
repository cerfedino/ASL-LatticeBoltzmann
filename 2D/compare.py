import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib   


TOLERANCE = 0.01

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


def is_equal(folder, name):
    output = ""
    error = False


    output+=f"\nComparing {folder}/{name} with reference"
    ref_arr = np.load("output/reference/" + name + ".npy")
    arr = np.load(f"{folder}/{name}.npy")

    if ref_arr.shape == arr.shape:
        output+=f"\n\t[+] Matrix shape matches {ref_arr.shape}"
    else:
        output+=f"\n\t{bcolors.FAIL}[ERROR] Matrix shape does not match! {arr.shape} should be {ref_arr.shape}{bcolors.ENDC}"
        error = True

    if np.all(ref_arr == arr):
        output+=f"\n\t[+] Matrices are equal"
    else:
        diff = np.std(ref_arr - arr)
        if diff < TOLERANCE:
            output+=f"\n\t[+] Matrices are equal within tolerance {TOLERANCE}"
        else:
            output+=f"\n\t{bcolors.FAIL}[ERROR] Matrices are not equal within tolerance {TOLERANCE}{bcolors.ENDC}"
            error = True

        ##print(f"Reference: {ref_arr[:5]}")
        ##print(f"Output: {arr[:5]}")
    if error:
        print(output)
    return error


# reference files in output/reference
output_refernce_folders = [x[0] for x in os.walk("output/reference")][1:]

# get folders in output that start with 2024
# output_folders = [x[0] for x in os.walk("output/")][1:]

# output_folder = output_folders[0]
output_folder = "output/00_latest"

print(f"Selected folder: {output_folder}")


# compare npy file F_start.npy
# is_equal(output_folder, "F_start")

# # compare npy file X.npy
# is_equal(output_folder, "X_start")

# # compare npy file Y.npy
# is_equal(output_folder, "X_start")

# # compare npy file F_3_setup.npy
# is_equal(output_folder, "F_3_setup")

# # compare npy file rho_start.npy
# is_equal(output_folder, "rho_start")

# # compare npy file F_normalized_start.npy
# is_equal(output_folder, "F_normalized_start")

# # compare npy file cylinder_start.npy
# is_equal(output_folder, "cylinder_start")

N = 30

error = False
for i in range(N):
    # is_equal(output_folder, f"F_before_drift_{i}")
    # is_equal(output_folder, f"F_after_drift_{i}")

    # is_equal(output_folder, f"bndryF_start_{i}")
    # is_equal(output_folder, f"bndryF_reordered_{i}")

    # is_equal(output_folder, f"rho_loop_{i}")
    # is_equal(output_folder, f"ux_loop_{i}")
    # is_equal(output_folder, f"uy_loop_{i}")

    # is_equal(output_folder, f"Feq_{i}")

    error = is_equal(output_folder, f"vorticity_{i:05}") or error 

print(f"[-] Checked the first {N} timesteps")    
if error:
    print(f"{bcolors.FAIL}[X] The implementation does not follow the baseline{bcolors.ENDC}")
    exit(1)
else:
    print(f"{bcolors.OKGREEN}[+] The implementation follows the baseline{bcolors.ENDC}")
    exit(0)