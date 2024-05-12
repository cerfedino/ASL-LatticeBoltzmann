import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import matplotlib   
import json
import subprocess


TOLERANCE = 0.01
step_size = 5

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
    #gotta fix it up from here on
    reference_location = "build_original/output/" + name +".csv"
    comp_location = folder +"/"+ name +".csv"

    ref_data = pd.read_csv(reference_location)
    comp_data = pd.read_csv(comp_location)

    #print(ref_data.shape)
    if ref_data.shape == comp_data.shape:
        output += f"\n\t[+] Matrix shape matches {ref_data.shape}"
    else:
        output+=f"\n\t{bcolors.FAIL}[ERROR] Matrix shape does not match! {comp_data.shape} should be {ref_data.shape}{bcolors.ENDC}"
        error = True
        
   
    if((ref_data.values == comp_data.values).all()):
         output+=f"\n\t[+] Matrices are equal"

    else:
        diff = ref_data.values - comp_data.values
        if(abs(diff) <= TOLERANCE).all():
            output+=f"\n\t[+] Matrices are equal within tolerance {TOLERANCE}"
        else:
            output+=f"\n\t{bcolors.FAIL}[ERROR] Matrices are not equal within tolerance {TOLERANCE}{bcolors.ENDC}"
            error = True

    if(error):
        print(output)

    return error

# Function to update a value in the options.json for both build and build_original
def update_json_value(key, new_value):
    build_json = 'build/options.json'
    # Read the JSON file and parse its contents
    with open(build_json, 'r') as file:
        data = json.load(file)

    # Update the value corresponding to the given key
    data[key] = new_value
    
    # Write the updated data back to the JSON file
    with open(build_json, 'w') as file:
        json.dump(data, file, indent=4)

    build_original_json = 'build_original/options.json'
    with open(build_original_json, 'r') as file:
        data = json.load(file)

    # Update the value corresponding to the given key
    data[key] = new_value
    
    # Write the updated data back to the JSON file
    with open(build_original_json, 'w') as file:
        json.dump(data, file, indent=4)
    
    
    

def test():
    # reference files in output/reference
    output_refernce_folders = [x[0] for x in os.walk("output/reference")][1:]

    # output_folder = output_folders[0]
    output_folder = "build/output"

    print(f"Selected folder: {output_folder}")


    N = 30*step_size

    error = False
    for i in range(0,N,5):

        error = is_equal(output_folder, f"{i}") or error 

    print(f"[-] Checked the first {N} timesteps with step size of {step_size}")    
    if error:
        print(f"{bcolors.FAIL}[X] The implementation does not follow the baseline{bcolors.ENDC}")
        exit(1)
    else:
        print(f"{bcolors.OKGREEN}[+] The implementation follows the baseline{bcolors.ENDC}")
    

grid_sizes = [
    ('grid_size', '[1,2,3]'),
    ('grid_size', '[3,10,300]'),
    ('grid_size', '[100,200,300]')
]

# Loop through the list of values and update each value in the JSON file, then compare the values
for grid, grid_size in grid_sizes:
    print("attempt with "+ grid_size)
    update_json_value(grid, grid_size)
    subprocess.check_output('build/main_cmdline')
    process = subprocess.Popen('build_original/main_cmdline', stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate(input=b'1\n')
    test()