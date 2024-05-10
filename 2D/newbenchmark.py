import os
import subprocess
import re
import numpy as np
from dotenv import load_dotenv
import matplotlib.pyplot as plt


PROJECT_ROOT = "../"
PREVIOUS_VERSIONS_PATH = f"./previous_versions/"
OUTPUT_FOLDER = f"./artifacts"



# Read MEM_BW and PEAK_FLOPS from environment variables, otherwise from .env file, otherwise display error message
try:
    print("Reading environment variables")
    MEM_BW = os.environ["MEM_BW"]
    PEAK_FLOPS = os.environ["PEAK_FLOPS"]
    SIMD_LEN_BITS = os.environ["SIMD_LEN_BITS"]
except KeyError:
    print("Environment variables not found, reading from .env file")
    # make sure keys are there AND they hgave avalue set
    try:
        load_dotenv()
        MEM_BW = int(os.getenv("MEM_BW"))
        PEAK_FLOPS = float(os.getenv("PEAK_FLOPS"))
        SIMD_LEN_BITS = int(os.getenv("SIMD_LEN_BITS"))
    except FileNotFoundError:
        print("Error: Environment variables not found and .env file not found")
        exit(1)
    except KeyError:
        print("Error: Environment variables not found and .env file does not contain the required keys")
        exit(1)

print("MEM_BW: ", MEM_BW)
print("PEAK_FLOPS: ", PEAK_FLOPS)
print("SIMD_LEN_BITS: ", SIMD_LEN_BITS)
    


def make_roofline_plot(PEAK_SCALAR, MEM_BW, SIMD_LEN):
    XMIN = 0.03125
    XMAX = 2**4
    
    PLT_FACECOLOR="#E4E4E4"
    PLT_BOUND_COLOR="#55575C"
    
    DOUBLE = 64

    PEAK_simd = PEAK_SCALAR * (SIMD_LEN/DOUBLE)
    
    fig = plt.figure(figsize=(20, 12), facecolor="white")
    plt.gca().set_facecolor(PLT_FACECOLOR)
    
    
    #  Memory bound
    x = np.linspace(XMIN, XMAX, 100)
    y_mem_scalar = [(MEM_BW * i) for i in x ]
    plt.plot(x, y_mem_scalar, color=PLT_BOUND_COLOR, linestyle='--', linewidth=1.4)
    
    ### No SIMD bounds
    # Compute bound
    y_comp_scalar = [PEAK_SCALAR for i in x]
    plt.plot(x, y_comp_scalar, color=PLT_BOUND_COLOR, linestyle='--', linewidth=1.4)
    plt.text(0.04,y_comp_scalar[0], f"Peak $\pi$ scalar ({y_comp_scalar[0]} iops/cycle)", fontsize=12, va='bottom', ha='left', color=PLT_BOUND_COLOR) 

    plt.text(0.035, 0.6, f"Read $\\beta$ ({MEM_BW} Bytes/Cycle)", fontsize=12, va='bottom', ha='left', color=PLT_BOUND_COLOR, rotation=33)

    # plt.plot(x, [min(y_mem_scalar[i], y_comp_scalar[i]) for i in range(len(x))], color='black', linestyle='-', linewidth=1.4)

    ### SIMD bounds
    # Compute bound
    y_comp_simd = [PEAK_simd for i in x]
    plt.plot(x, y_comp_simd, color=PLT_BOUND_COLOR, linestyle='--', linewidth=1.4)
    plt.text(0.04,y_comp_simd[0], f"Peak $\pi$ SIMD ({y_comp_simd[0]} iops/cycle)", fontsize=12, va='bottom', ha='left', color=PLT_BOUND_COLOR)

    plt.gca().yaxis.grid(True, which='major', color='w', linestyle='-', linewidth=1.5)
    plt.xlabel("Operational Intensity I(n) [iops/byte]", fontsize=14, rotation=0, labelpad=5, ha='center')
    plt.ylabel("Performance P(n) [iops/cycle]", fontsize=14, rotation=0, ha='left', labelpad=-15, position=(1,1))

    # Plot bounds
    plt.ylim(np.min(y_mem_scalar), 2**6)
    plt.xlim(np.min(x), XMAX)
    plt.xscale("log", base=2)
    plt.yscale("log", base=2)
    
    
    return plt


def run_executable_and_get_output(executable_path, *args):
    print(f"Running {executable_path} {' '.join(map(str, args))}")
    result = subprocess.run([executable_path] + list(args), stdout=subprocess.PIPE, text=True)
    return result.stdout



def main():
    previous_versions = []
    # Scan for all previous version that we plan to benchmark
    for folder in os.listdir(PREVIOUS_VERSIONS_PATH):
        # Check if the subfolder is a directory
        if not os.path.isdir(PREVIOUS_VERSIONS_PATH + folder):
            continue
        previous_versions.append(PREVIOUS_VERSIONS_PATH + folder)
    
    
    
    # Init empty roofline plot
    plt = make_roofline_plot(PEAK_FLOPS, MEM_BW, SIMD_LEN_BITS)
    plt.savefig(f"{OUTPUT_FOLDER}/roofline_plot.out.pdf")
    plt.show()
    
    for path in previous_versions:
        print("=========")
        env = {}
        with open(path + "/.env") as file:
            for line in file:
                key, value = line.strip().split("=")
                env[key] = value
                
        print("\nVERSION: ", env["VERSION"])
        print("Path: ", path)
        print("TITLE: ", env["TITLE"])
        print("DESCRIPTION: ", env["DESC"])
        print("COMMIT: ", env["COMMIT"])
        
        print("")
        # Execute executable file in subprocess and read stdout
        output = run_executable_and_get_output(f"{path}/main_optimized.old.o", '400', '200', '100')
        flops_per_cycle = np.sum([float(f) for f in re.findall(r"(\d+\.\d+) Flops/Cycle", output)])
        print("Recorded total flops/cycle : ", flops_per_cycle)
        
main()