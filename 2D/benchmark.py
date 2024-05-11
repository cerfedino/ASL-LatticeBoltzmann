import os
import subprocess
import re
import numpy as np
from dotenv import load_dotenv
import matplotlib.pyplot as plt
from sympy import symbols, evalf, sympify


PREVIOUS_VERSIONS_PATH = f"{os.path.dirname(__file__)}/previous_versions/"
OUTPUT_FOLDER = f"{os.path.dirname(__file__)}/artifacts/"

Nx, Ny, Nt, NL = symbols('Nx Ny Nt NL')

PARAMS= {Nx: 400, Ny: 200, NL: 9, Nt:100}
# Read MEM_BW and PEAK_SCALAR from environment variables, otherwise from .env file, otherwise display error message
try:
    print("Reading environment variables")
    MEM_BW = int(os.environ["MEM_BW"])
    PEAK_SCALAR = float(os.environ["PEAK_SCALAR"])
    SIMD_LEN_BITS = int(os.environ["SIMD_LEN_BITS"])
except KeyError:
    print("Environment variables not found, reading from .env file")
    # make sure keys are there AND they hgave avalue set
    try:
        load_dotenv()
        MEM_BW = int(os.getenv("MEM_BW"))
        PEAK_SCALAR = float(os.getenv("PEAK_SCALAR"))
        SIMD_LEN_BITS = int(os.getenv("SIMD_LEN_BITS"))
    except FileNotFoundError:
        print("Error: Environment variables not found and .env file not found")
        exit(1)
    except KeyError:
        print("Error: Environment variables not found and .env file does not contain the required keys")
        exit(1)

print("MEM_BW: ", MEM_BW)
print("PEAK_SCALAR: ", PEAK_SCALAR)
print("SIMD_LEN_BITS: ", SIMD_LEN_BITS)
PEAK_simd = PEAK_SCALAR * (SIMD_LEN_BITS/64)


def make_roofline_plot(PEAK_SCALAR, MEM_BW, SIMD_LEN_BITS):
    XMIN = 0.03125
    XMAX = 2**4
    
    PLT_FACECOLOR="#E4E4E4"
    PLT_BOUND_COLOR="#55575C"
    
    
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

def plot_roofline(plt, label, perf, intensity, simd=False):
    plt.plot(intensity, perf, 'ro', label=label)
    plt.text(intensity, perf, f"{label}\n ({perf:.2f} iops/cycle)", fontsize=12, va='top', ha='center', color='red')
    return plt



def run_executable_and_get_output(executable_path, args):
    print(f"Running {executable_path} {' '.join(args)}")
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
    plt = make_roofline_plot(PEAK_SCALAR, MEM_BW, SIMD_LEN_BITS)
    
    for path in previous_versions:
        print("=========")
        
        load_dotenv(dotenv_path=path+"/.env")
        VERSION = os.environ["VERSION"]
        PATH = path
        TITLE = os.environ["TITLE"]
        DESC = os.environ["DESC"]
        COMMIT = os.environ["COMMIT"]
        WORK = os.environ["WORK"]
        DATA_MOVEMENT = os.environ["DATA_MOV"]

        print("\nVERSION: ", VERSION)
        print("Path: ", PATH)
        print("TITLE: ", TITLE)
        print("DESCRIPTION: ", DESC)
        print("COMMIT: ", COMMIT)
        print("WORK: ", WORK)
        print("DATA_MOV: ", DATA_MOVEMENT)
        
        # Execute executable file in subprocess and read stdout
        output = run_executable_and_get_output(f"{path}/main_optimized.o", [f"{PARAMS[Nx]}", f"{PARAMS[Ny]}", f"{PARAMS[Nt]}"])
        
        print(output)
        runs = int(re.search(r"Run (\d+)/\d+ done\nProfiling results:", output).group(1))
        cycles = int(re.search(r"Cycles taken: (\d+)", output).group(1))
        print(f"Total Number of Cycles: {int(cycles/runs)}")
        
        WORK_eval = int(sympify(WORK).subs(PARAMS).evalf())
        DATA_eval = int(sympify(DATA_MOVEMENT).subs(PARAMS).evalf())
        INTENSITY = WORK_eval / DATA_eval
        PERFORMANCE = WORK_eval / int(cycles/runs)
        
        print("")
        print("Work: ", WORK, "= " + str(WORK_eval))
        print("Data Movement: ", DATA_MOVEMENT, "= " + str(DATA_eval))
        print("Operational Intensity: ", INTENSITY)
        print("Performance: ", PERFORMANCE)
        
        
        plt = plot_roofline(plt, TITLE, PERFORMANCE, INTENSITY, simd=True)
        plt.savefig(f"{OUTPUT_FOLDER}/roofline_plot.out.pdf")
        plt.show()

        
main()