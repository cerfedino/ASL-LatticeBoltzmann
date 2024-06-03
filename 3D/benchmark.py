import os
import subprocess
import re
import math
import numpy as np
import pandas as pd
from dotenv import load_dotenv, dotenv_values
import matplotlib.pyplot as plt
from sympy import symbols, evalf, sympify
import sys


PREVIOUS_VERSIONS_PATH = os.path.join(os.path.dirname(__file__), 'previous_versions/')
OUTPUT_FOLDER = os.path.join(os.path.dirname(__file__),'artifacts/')


if len(sys.argv) != 5:
    print("Usage: python3 benchmark.py <Nx> <Ny> <Nz> <Nt>")
    exit(1)

MNx = int(sys.argv[1])
MNy = int(sys.argv[2])
MNz = int(sys.argv[3])
MNt = int(sys.argv[4])

Nx, Ny, Nz, Nt, NL = symbols('Nx Ny Nz Nt NL')
NX, NY, NZ, NT, NL = symbols('NX NY NZ NT NL')

PARAMS= {Nx: MNx, Ny: MNy, Nz: MNz, NL: 15, Nt: MNt}
PARAMS2= {NX: MNx, NY: MNy, NZ: MNz, NL: 15, NT: MNt}

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

fig_index = 0
def make_roofline_plot(PEAK_SCALAR, MEM_BW, SIMD_LEN_BITS):
    global fig_index
    XMIN = 1E-5
    XMAX = 2**4
    
    PLT_FACECOLOR="#E4E4E4"
    PLT_BOUND_COLOR="#55575C"
    
    fig = plt.figure(fig_index, figsize=(20, 12), facecolor="white")
    fig_index+=1
    plt.gca().set_facecolor(PLT_FACECOLOR)
    
    plt.xscale("log", base=2)
    plt.yscale("log", base=2)
    
    #  Memory bound
    x = np.linspace(XMIN, XMAX, 120)
    memory_line, = plt.plot((XMIN, XMAX), (MEM_BW*XMIN, MEM_BW*XMAX), color=PLT_BOUND_COLOR, linestyle='--', linewidth=1.4)
    
    ### No SIMD bounds
    # Compute bound
    y_comp_scalar = [PEAK_SCALAR for i in x]
    plt.plot(x, y_comp_scalar, color=PLT_BOUND_COLOR, linestyle='--', linewidth=1.4)
    plt.text(0.04,y_comp_scalar[0], f"Peak $\\pi$ scalar ({y_comp_scalar[0]} iops/cycle)", fontweight="bold", fontsize=12, va='bottom', ha='left', color=PLT_BOUND_COLOR) 

    ax = plt.gca()
    memory_coords = memory_line.get_data()
    p1 = ax.transData.transform_point(memory_coords[0])
    p2 = ax.transData.transform_point(memory_coords[-1])
    text_rot = 90 - np.degrees(np.arctan2(p2[1] - p1[1], p2[0] - p1[0]))
    text = plt.text(np.min(x) * 1.5, MEM_BW * np.min(x) * 1.5 * 1.5, f"Read $\\beta$ ({MEM_BW} Bytes/Cycle)", fontweight="bold", fontsize=12, va='bottom', ha='left', color=PLT_BOUND_COLOR, rotation=text_rot)

    ### SIMD bounds
    # Compute bound
    y_comp_simd = [PEAK_SCALAR * SIMD_LEN_BITS / 64 for i in x]
    plt.plot(x, y_comp_simd, color=PLT_BOUND_COLOR, linestyle='--', linewidth=1.4)
    plt.text(0.04,y_comp_simd[0], f"Peak $\\pi$ SIMD ({y_comp_simd[0]} iops/cycle)", fontweight="bold", fontsize=12, va='bottom', ha='left', color=PLT_BOUND_COLOR)

    plt.gca().yaxis.grid(True, which='major', color='w', linestyle='-', linewidth=1.5)
    plt.xlabel("Operational Intensity I(n) [iops/byte]", fontsize=14, rotation=0, labelpad=5, ha='center')
    plt.ylabel("Performance P(n) [iops/cycle]", fontsize=14, rotation=0, ha='left', labelpad=-15, position=(1,1))

    # Plot bounds
    plt.xlim(np.min(x), XMAX)
    
    
    return fig

def plot_roofline(plt, label, perf, intensity, simd=False):
    plt.gca().plot(intensity, perf, 'ro', label=label)
    plt.gca().text(intensity * 1.15, perf, f"{label} ({perf:.2f} iops/cycle)", fontsize=12, ha='left', color='red')
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

    if len(previous_versions) == 0:
        print("No previous versions found, stopping...")
        exit(0)
    
    previous_versions.sort()

    # Init empty roofline plot
    plt_all = make_roofline_plot(PEAK_SCALAR, MEM_BW, SIMD_LEN_BITS)
  
    '''
    
    We have 14*4 additions, 1 division and 3 multiplications
    Total flops - 14*4+1+3
    Ports 2 and 3 have to process 14*4 additions, which would take 14*2 cycles
    Thus we would have (14*4+1+3)/28 = 2.14 flops/cycle
    '''
    plt_density_momentum  = make_roofline_plot(2.14, MEM_BW, SIMD_LEN_BITS); plt_density_momentum.suptitle("Compute density momentum moment")
 
    '''
    We got only 2 FMAs... Its kinda sad
    '''
    plt_stream  = make_roofline_plot(4, MEM_BW, SIMD_LEN_BITS); plt_stream.suptitle("Stream")
    
    '''
    We can do 2 FMAs and 2 Adds, and those are the only ops we do in the whole for-loops basically, so in the best case we would start one at every cycle at every port
    '''
    plt_collision   = make_roofline_plot(6, MEM_BW, SIMD_LEN_BITS); plt_collision.suptitle("Collision")

    min_performance = math.inf

    min_performance_density_momentum  = math.inf
    min_performance_collision  = math.inf
    min_performance_stream    = math.inf

    version_labels = []
    loop_cycles = {
        "density_momentum": [],
        "stream": [],
        "collision": [],
    }
  
    for path in previous_versions:
        print("=========")
        
        conf = dotenv_values(path+"/.env")
        VERSION = conf["VERSION"]
        PATH = path
        TITLE = conf["TITLE"]
        DESC = conf["DESC"]
        COMMIT = conf["COMMIT"]
        WORK = conf["WORK"]
        DATA_MOVEMENT = conf["DATA_MOV"]

        version_labels.append(f"{TITLE}\n{DESC}")

        print("\nVERSION: ", VERSION)
        print("Path: ", PATH)
        print("TITLE: ", TITLE)
        print("DESCRIPTION: ", DESC)
        print("COMMIT: ", COMMIT)
        print("WORK: ", WORK)
        print("DATA_MOV: ", DATA_MOVEMENT)
        
        # Execute executable file in subprocess and read stdout
        # If there is a specific executable for the current configuration, use it instead
        if os.path.exists(f"{path}/src/cmdline/main_{MNx}_{MNy}_{MNz}_{MNt}_profiler.o"):
            executable = f"{path}/src/cmdline/main_{MNx}_{MNy}_{MNz}_{MNt}_profiler.o"
        else:
            print(f"Didn't find profiler for this path: {PATH}")
            exit(1)
            executable = f"{path}/src/main_optimized_profile.o"
            
        output = run_executable_and_get_output(executable, [f"{PARAMS[Nx]}", f"{PARAMS[Ny]}", f"{PARAMS[Nz]}", f"{PARAMS[Nt]}"])
        print("done")
        try:
            runs = int(re.search(r"Run (\d+)/\d+ done\nProfiling results:", output).group(1))

            stats_density_momentum  = re.search(r".Compute density\s*Calculation: ([\d|\\.]*) Flops/Cycle.*Arithmetic intensity: ([\d|\\.]*)", output).groups()
            performance_density_momentum, intensity_density_momentum = float(stats_density_momentum[0]), float(stats_density_momentum[1])
            stats_collision  = re.search(r".*Compute collision\s*Calculation: ([\d|\\.]*) Flops/Cycle.*Arithmetic intensity: ([\d|\\.]*)", output).groups()
            performance_collision, intensity_collision = float(stats_collision[0]), float(stats_collision[1])
            stats_stream    = re.search(r".*Compute stream\s*Calculation: ([\d|\\.]*) Flops/Cycle.*Arithmetic intensity: ([\d|\\.]*)", output).groups()
            performance_stream, intensity_stream = float(stats_stream[0]), float(stats_stream[1])

            cycles_density_momentum  = int(re.search(r".*Compute density\s*Calculation: .* (\d+) cycles.*", output).group(1))
            cycles_collision  = int(re.search(r".*Compute collision\s*Calculation: .* (\d+) cycles.*", output).group(1))
            cycles_stream    = int(re.search(r".*Compute stream\s*Calculation: .* (\d+) cycles.*", output).group(1))
            
            cycles = cycles_density_momentum + cycles_collision + cycles_stream

            loop_cycles["density_momentum"].append(cycles_density_momentum); loop_cycles["collision"].append(cycles_collision); loop_cycles["stream"].append(cycles_stream); 
        except AttributeError as e:
            print("Couldn't match regex to output:\n")
            print(output)
            print("\n\nException: ", e)
        print()
        print(f"Compute density cycles: {cycles_density_momentum }")
        print(f"Compute collision cycles: {cycles_collision }")
        print(f"Compute stream cycles: {cycles_stream }")
        print(f"Total Number of Cycles: {cycles}")
        
        work_string = sympify(WORK).subs(PARAMS).subs(PARAMS2).evalf()
        print(work_string)
        WORK_eval = int(work_string)

        DATA_eval = int(sympify(DATA_MOVEMENT).subs(PARAMS).subs(PARAMS2).evalf())
        INTENSITY = WORK_eval / DATA_eval
        PERFORMANCE = WORK_eval / cycles
        
        min_performance = min(min_performance, PERFORMANCE)

        min_performance_density_momentum  = min(min_performance_density_momentum , performance_density_momentum )
        min_performance_collision  = min(min_performance_collision , performance_collision )
        min_performance_stream    = min(min_performance_stream   , performance_stream   )
        
        print()
        print("Work: ", WORK, "= " + str(WORK_eval))
        print("Data Movement: ", DATA_MOVEMENT, "= " + str(DATA_eval))
        print("Operational Intensity: ", INTENSITY)
        print("Performance: ", PERFORMANCE)
        
        plt_all = plot_roofline(plt_all, TITLE, PERFORMANCE, INTENSITY, simd=True)

        plt_density_momentum  = plot_roofline(plt_density_momentum , TITLE, performance_density_momentum , intensity_density_momentum , simd=True)
        plt_collision  = plot_roofline(plt_collision , TITLE, performance_collision , intensity_collision , simd=True)
        plt_stream    = plot_roofline(plt_stream   , TITLE, performance_stream   , intensity_stream   , simd=True)

    plt_all.gca().set_ylim(min(min_performance, PEAK_SCALAR)/2.5, 2.5 * PEAK_simd)

    plt_density_momentum.gca().set_ylim(min(min_performance_density_momentum, PEAK_SCALAR)/2.5, 2.5 * PEAK_simd)
    plt_collision.gca().set_ylim(min(min_performance_collision, PEAK_SCALAR)/2.5, 2.5 * PEAK_simd)
    plt_stream.gca().set_ylim(min(min_performance_stream, PEAK_SCALAR)/2.5, 2.5 * PEAK_simd)
    print(version_labels)

    bar_fig = plt.figure("bar", figsize=(20, 12), facecolor="white")
    df = pd.DataFrame.from_dict(loop_cycles, columns=version_labels, orient='index')
    normalized_df = df/df.sum()
    normalized_df = normalized_df.transpose()
    bottom = np.zeros(len(version_labels))
    for col in normalized_df.columns:
        p = bar_fig.gca().bar(version_labels, normalized_df[col], 0.5, label=col, bottom=bottom)
        bar_fig.gca().bar_label(p, fmt=col.upper(), label_type='center', size='x-large')
        bottom += normalized_df[col]
    bar_fig.gca().legend(loc="upper right", fontsize='large')
    bar_fig.gca().tick_params(labelsize='x-large')

    plt_all.savefig(f"{OUTPUT_FOLDER}/roofline_plot.out.pdf")
    plt_density_momentum.savefig(f"{OUTPUT_FOLDER}/roofline_plot_density_momentum.out.pdf")
    plt_collision.savefig(f"{OUTPUT_FOLDER}/roofline_plot_collision.out.pdf")
    #plt_stream.savefig(f"{OUTPUT_FOLDER}/roofline_plot_stream.out.pdf")
    bar_fig.savefig(f"{OUTPUT_FOLDER}/bar_plot.out.pdf")

    plt.show()

        
main()