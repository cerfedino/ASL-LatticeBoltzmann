import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib   

#matplotlib.use("Qt5agg")


def is_equal(folder, name):
    print("--------------------------------------------")
    print(f"Comparing {name} in {folder} with reference")
    ref_arr = np.load("output/reference/" + name + ".npy")
    arr = np.load(f"{folder}/{name}.npy")

    print(f"    Shape of reference: {ref_arr.shape}")
    print(f"    Shape of output: {arr.shape}")

    print(f"  Is equal: {np.all(ref_arr == arr)}")

    # if not equal, print the difference
    if not np.all(ref_arr == arr):
        print(f"Difference: {np.sum(ref_arr - arr)}")

        # if difference is more than 1000 absolute
#        if np.sum(np.abs(ref_arr - arr)) > 1000:
#            print(f"Reference: {ref_arr[:5]}")
#            print(f"Output: {arr[:5]}")

        ##print(f"Reference: {ref_arr[:5]}")
        ##print(f"Output: {arr[:5]}")


# reference files in output/reference
output_refernce_folders = [x[0] for x in os.walk("output/reference")][1:]

# get folders in output that start with 2024
output_folders = [x[0] for x in os.walk("output/")][1:]

output_folder = output_folders[0]

print(f"Selected folder: {output_folder}")


# compare npy file F_start.npy
is_equal(output_folder, "F_start")

# compare npy file X.npy
is_equal(output_folder, "X_start")

# compare npy file Y.npy
is_equal(output_folder, "X_start")

# compare npy file F_3_setup.npy
is_equal(output_folder, "F_3_setup")

# compare npy file rho_start.npy
is_equal(output_folder, "rho_start")

# compare npy file F_normalized_start.npy
is_equal(output_folder, "F_normalized_start")

# compare npy file cylinder_start.npy
is_equal(output_folder, "cylinder_start")

N = 5

for i in range(N):
    print("############################################")
    print("Iteration: ", i)
    print("############################################")

    is_equal(output_folder, f"F_before_drift_{i}")
    is_equal(output_folder, f"F_after_drift_{i}")

    is_equal(output_folder, f"bndryF_start_{i}")
    is_equal(output_folder, f"bndryF_reordered_{i}")

    is_equal(output_folder, f"rho_loop_{i}")
    is_equal(output_folder, f"ux_loop_{i}")
    is_equal(output_folder, f"uy_loop_{i}")

    is_equal(output_folder, f"Feq_{i}")


    # with 5 leading zeros
    is_equal(output_folder, f"vorticity_{i:05}")

