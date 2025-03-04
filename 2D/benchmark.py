import os
import subprocess
import re
import math
import numpy as np
import pandas as pd
from dotenv import load_dotenv, dotenv_values
import matplotlib.pyplot as plt
from matplotlib.colors import hsv_to_rgb
from sympy import symbols, evalf, sympify
import sys


PREVIOUS_VERSIONS_PATH = os.path.join(os.path.dirname(__file__), 'previous_versions/')
OUTPUT_FOLDER = os.path.join(os.path.dirname(__file__),'artifacts/')

if len(sys.argv) != 4:
    print("Usage: python3 benchmark.py <Nx> <Ny> <Nt>")
    exit(1)

MNx = int(sys.argv[1])
MNy = int(sys.argv[2])
MNt = int(sys.argv[3])


Nx, Ny, Nt, NL = symbols('Nx Ny Nt NL')

PARAMS= {Nx: MNx, Ny: MNy, NL: 9, Nt: MNt}
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
    XMAX = 2**6
    
    PLT_FACECOLOR="#E4E4E4"
    PLT_BOUND_COLOR="#55575C"
    
    fig = plt.figure(fig_index, figsize=(20, 12), facecolor="white")
    fig_index+=1
    plt.gca().set_facecolor(PLT_FACECOLOR)
    
    plt.xscale("log", base=2)
    plt.yscale("log", base=2)
    
    #  Memory bound
    x = np.linspace(XMIN, XMAX, 100)
    memory_line, = plt.plot((XMIN, XMAX), (MEM_BW*XMIN, MEM_BW*XMAX), color=PLT_BOUND_COLOR, linestyle='--', linewidth=1.4)
    
    ### No SIMD bounds
    # Compute bound
    y_comp_scalar = [PEAK_SCALAR for i in x]
    plt.plot(x, y_comp_scalar, color=PLT_BOUND_COLOR, linestyle='--', linewidth=1.4)
    plt.text(0.04,y_comp_scalar[0], f"Peak $\\pi$ scalar ({y_comp_scalar[0]} flops/cycle)", fontweight="bold", fontsize=12, va='bottom', ha='left', color=PLT_BOUND_COLOR) 

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
    plt.text(0.04,y_comp_simd[0], f"Peak $\\pi$ SIMD ({y_comp_simd[0]} flops/cycle)", fontweight="bold", fontsize=12, va='bottom', ha='left', color=PLT_BOUND_COLOR)

    plt.gca().yaxis.grid(True, which='major', color='w', linestyle='-', linewidth=1.5)
    plt.xlabel("Operational Intensity I(n) [flops/byte]", fontsize=14, rotation=0, labelpad=5, ha='center')
    plt.ylabel("Performance P(n) [flops/cycle]", fontsize=14, rotation=0, ha='left', labelpad=-15, position=(1,1))

    # Plot bounds
    plt.xlim(np.min(x), XMAX)
    
    
    return fig

def plot_roofline(plt, label, perf, intensity, color, name=None):
    print(f"Plotting {name} with performance {perf} and intensity {intensity}")
    plt.gca().plot(intensity, perf, 'ro', label=f"{label} ({perf:.2f} flops/cycle)", color=color)
    texts = plt.gca().text(intensity * 1.05, perf, label, fontsize=12, va='center', ha='left', color=color)
    #adjust_text(texts)
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
        
    # Get sorted indexes of previous_versions
    previous_versions.sort(key=lambda x: int(os.path.basename(x).lstrip("_").lstrip("v")))
        
    # Init empty roofline plot
    plt_all = make_roofline_plot(PEAK_SCALAR, MEM_BW, SIMD_LEN_BITS)

    '''
    Port 1 is the bottleneck. We have N divisions and 20N mults. N divisions are scheduled on Port 1, and the remaining 19N multiplications on both ports 0 and 1.
    Thus, port 1 has to process N + 19N/2 = 10.5N ops. With a throughput of 2 ops/cycle, this results in a run time of 5.25N cycles.
    Total flops are 48N flops.
    Max performance is thus 48N/5.25N =~ 9.15 flops/cycle
    '''
    plt_rho  = make_roofline_plot(9.15, MEM_BW, SIMD_LEN_BITS); plt_rho.suptitle("Rho")
    '''
    Ports 0/1 are the bottlenecks for this calculation
    These perform 4 ops per cycle together.
    In total, they have to process 10*N ops, which thus take 10*N/4 cycles = 2.5*N cycles.
    The calculation has a total of 15*N ops, meaning we have a max performance of 15/2.5 = 6 flops/cycle
    '''
    plt_feq  = make_roofline_plot(6, MEM_BW, SIMD_LEN_BITS); plt_feq.suptitle("FEQ")
    # One FMA per cycle, thus we can only utilize the FMA ports, of which we have 2, which perform two ops per cycle
    plt_f    = make_roofline_plot(8, MEM_BW, SIMD_LEN_BITS); plt_f.suptitle("F")
    # Only has subtractions, so we can make use of the two ADD/SUB ports, which perform two ops per cycle
    plt_vort = make_roofline_plot(4, MEM_BW, SIMD_LEN_BITS); plt_vort.suptitle("Vort")

    min_performance = math.inf

    min_performance_rho  = math.inf
    min_performance_feq  = math.inf
    min_performance_f    = math.inf
    min_performance_vort = math.inf

  
    colorrange = (0, 90)
    value_factor = 0.7  # Adjust this factor to make yellows darker

    version_labels = []
    loop_cycles = {
        "rho": [],
        "feq": [],
        "f": [],
        "vort": [],
    }
    
    for path in previous_versions:
        print("=========")
        v_idx = previous_versions.index(path)
        
        if not os.path.exists(path+"/.env"):
            print(f"The .env file for {path} does not exist, skipping")
            exit(1)
        
        conf = dotenv_values(path+"/.env")
        VERSION = conf["VERSION"]
        PATH = path
        TITLE = conf["TITLE"]
        DESC = conf["DESC"]
        COMMIT = conf["COMMIT"]
        WORK = conf["WORK"]
        DATA_MOVEMENT = conf["DATA_MOV"]

        version_labels.append(VERSION)

        print("\nVERSION: ", VERSION)
        print("Path: ", PATH)
        print("TITLE: ", TITLE)
        print("DESCRIPTION: ", DESC)
        print("COMMIT: ", COMMIT)
        print("WORK: ", WORK)
        print("DATA_MOV: ", DATA_MOVEMENT)
        
        # Run the profiling executable in order to read cycles measurements
        executable = f"{path}/src/main_optimized_profile_{MNx}_{MNy}_{MNt}.o"
            
        output = run_executable_and_get_output(executable, [f"{PARAMS[Nx]}", f"{PARAMS[Ny]}", f"{PARAMS[Nt]}"])
        
        try:
            cycles_drift = re.search(r".*Drift\s*Calculation: .* (\d+) cycles.*", output)
            cycles_drift = int(cycles_drift.group(1)) if cycles_drift is not None and int(cycles_drift.group(1)) > 1000 else 1
            
            cycles_rho  = re.search(r".*Rho\s*Calculation: .* (\d+) cycles.*", output)
            cycles_rho  = int(cycles_rho.group(1)) if cycles_rho is not None and int(cycles_rho.group(1)) > 1000 else 1
            
            cycles_feq  = re.search(r".*FEQ\s*Calculation: .* (\d+) cycles.*", output)
            cycles_feq  = int(cycles_feq.group(1)) if cycles_feq is not None and int(cycles_feq.group(1)) > 1000 else 1
            
            cycles_f    = re.search(r".*F\s*Calculation: .* (\d+) cycles.*", output)
            cycles_f    = int(cycles_f.group(1)) if cycles_f is not None and int(cycles_f.group(1)) > 1000 else 1
            
            cycles_vort = re.search(r".*Vort\s*Calculation: .* (\d+) cycles.*", output)
            cycles_vort = int(cycles_vort.group(1)) if cycles_vort is not None and int(cycles_vort.group(1)) > 1000 else 1
            
            cycles = cycles_drift + cycles_rho + cycles_feq + cycles_f + cycles_vort

            loop_cycles["rho"].append(cycles_rho); loop_cycles["feq"].append(cycles_feq); loop_cycles["f"].append(cycles_f); loop_cycles["vort"].append(cycles_vort)
        except AttributeError as e:
            print("Couldn't match regex to output:\n")
            print(output)
            print("\n\nException: ", e)
            exit(1)

        print()
        if cycles_drift > 1:
          print(f"Drift cycles: {cycles_drift}")
        if cycles_rho > 1:
          print(f"Rho  cycles: {cycles_rho }")
        if cycles_feq > 1:
          print(f"FEQ  cycles: {cycles_feq }")
        if cycles_f > 1:
          print(f"F    cycles: {cycles_f   }")
        if cycles_vort > 1:
          print(f"Vort cycles: {cycles_vort}")
          print(f"Total Number of Cycles: {cycles}")
        
        # Run the benchmark executable in order to read PAPI measurements
        executable = f"{path}/src/main_optimized_benchmark_{MNx}_{MNy}_{MNt}.o"
            
        output = run_executable_and_get_output(executable, [f"{PARAMS[Nx]}", f"{PARAMS[Ny]}", f"{PARAMS[Nt]}"])
        
        try:
            stats_drift = re.search(r".*Drift Mem Transfer: ([\d|\.]*).*Drift Floating point operations: ([\d|\.]*).*Drift Arithmetic Intensity: ([\d|\.]*).*", output, re.DOTALL)
            if stats_drift is not None and stats_drift.groups()[2] != "":
              stats_drift = stats_drift.groups()
              mem_drift, flops_drift, intensity_drift = float(stats_drift[0]), float(stats_drift[1]), float(stats_drift[2])
            else:
              mem_drift, flops_drift, intensity_drift = 1, 1, 1
            
            stats_rho = re.search(r".*Rho Mem Transfer: ([\d|\.]*).*Rho Floating point operations: ([\d|\.]*).*Rho Arithmetic Intensity: ([\d|\.]*).*", output, re.DOTALL)
            if stats_rho is not None and stats_rho.groups()[2] != "":
              stats_rho = stats_rho.groups()
              mem_rho, flops_rho, intensity_rho = float(stats_rho[0]), float(stats_rho[1]), float(stats_rho[2])
            else:
              mem_rho, flops_rho, intensity_rho = 1, 1, 1
            
            stats_feq = re.search(r".*FEQ Mem Transfer: ([\d|\.]*).*FEQ Floating point operations: ([\d|\.]*).*FEQ Arithmetic Intensity: ([\d|\.]*).*", output, re.DOTALL)
            if stats_feq is not None and stats_feq.groups()[2] != "":
              stats_feq = stats_feq.groups()
              mem_feq, flops_feq, intensity_feq = float(stats_feq[0]), float(stats_feq[1]), float(stats_feq[2])
            else:
              mem_feq, flops_feq, intensity_feq = 1, 1, 1
            stats_f    = re.search(r".*F Mem Transfer: ([\d|\.]*).*F Floating point operations: ([\d|\.]*).*F Arithmetic Intensity: ([\d|\.]*).*", output, re.DOTALL)
            if stats_f is not None and stats_f.groups()[2] != "":
              stats_f = stats_f.groups()
              mem_f, flops_f, intensity_f = float(stats_f[0]), float(stats_f[1]), float(stats_f[2])
            else:
              mem_f, flops_f, intensity_f = 1, 1, 1
            
            stats_vort = re.search(r".*Vort Mem Transfer: ([\d|\.]*).*Vort Floating point operations: ([\d|\.]*).*Vort Arithmetic Intensity: ([\d|\.]*).*", output, re.DOTALL)
            if stats_vort is not None and stats_vort.groups()[2] != "":
              stats_vort = stats_vort.groups()
              mem_vort, flops_vort, intensity_vort = float(stats_vort[0]), float(stats_vort[1]), float(stats_vort[2])
            else:
              mem_vort, flops_vort, intensity_vort = 1, 1, 1

        except AttributeError as e:
            print("Couldn't match regex to output:\n")
            print(output)
            print("\n\nException: ", e)
            exit(1)
        
        WORK_eval = flops_rho + flops_feq + flops_f + flops_vort
        DATA_eval = mem_rho + mem_feq + mem_f + mem_vort
        INTENSITY = WORK_eval / DATA_eval
        PERFORMANCE = WORK_eval / cycles

        min_performance = min(min_performance, PERFORMANCE)
        
        performance_rho = flops_rho / cycles_rho
        performance_feq = flops_feq / cycles_feq
        performance_f   = flops_f   / cycles_f
        performance_vort= flops_vort/ cycles_vort
        

        min_performance_rho  = min(min_performance_rho , performance_rho )
        min_performance_feq  = min(min_performance_feq , performance_feq )
        min_performance_f    = min(min_performance_f   , performance_f   )
        min_performance_vort = min(min_performance_vort, performance_vort)
        
        print()
        print("Work: ", WORK, "= " + str(WORK_eval))
        print("Data Movement: ", DATA_MOVEMENT, "= " + str(DATA_eval))
        print("Operational Intensity: ", INTENSITY)
        print("Performance: ", PERFORMANCE)
        
        curcolor = hsv_to_rgb([[(colorrange[0]+((colorrange[1]-colorrange[0])*((v_idx+1)/len(previous_versions))))/360,  1, value_factor]])[0]

        if not os.path.basename(path).startswith("_"):
          plt_all = plot_roofline(plt_all, TITLE, PERFORMANCE, INTENSITY, curcolor)

        plt_rho  = plot_roofline(plt_rho , TITLE, performance_rho , intensity_rho, curcolor)
        plt_feq  = plot_roofline(plt_feq , TITLE, performance_feq , intensity_feq, curcolor)
        plt_f    = plot_roofline(plt_f   , TITLE, performance_f   , intensity_f, curcolor)
        plt_vort = plot_roofline(plt_vort, TITLE, performance_vort, intensity_vort, curcolor)


    plt_all.gca().set_ylim(min(min_performance, PEAK_SCALAR)/2.5, 2.5 * PEAK_simd)

    plt_rho.gca().set_ylim(min(min_performance_rho, PEAK_SCALAR)/2.5, 2.5 * PEAK_simd)
    plt_feq.gca().set_ylim(min(min_performance_feq, PEAK_SCALAR)/2.5, 2.5 * PEAK_simd)
    plt_f.gca().set_ylim(min(min_performance_f, PEAK_SCALAR)/2.5, 2.5 * PEAK_simd)
    plt_vort.gca().set_ylim(min(min_performance_vort, PEAK_SCALAR)/2.5, 2.5 * PEAK_simd)

    bar_fig = plt.figure("bar", figsize=(40, 12), facecolor="white")
    df = pd.DataFrame.from_dict(loop_cycles, columns=version_labels, orient='index')
    normalized_df = df/df.sum()
    normalized_df = normalized_df.transpose()
    bottom = np.zeros(len(version_labels))
    for col in normalized_df.columns:
        p = bar_fig.gca().bar(['_'+i for i in version_labels], normalized_df[col], 0.5, label=col, bottom=bottom)
        # bar_fig.gca().bar_label(p, fmt=col.upper(), label_type='center', size='x-large')
        bottom += normalized_df[col]
    bar_fig.gca().legend(loc="upper right", fontsize='large') 
    bar_fig.gca().tick_params(labelsize='large')
    
    x_ticks = range(len(version_labels))
    x_labels = ['v'+i for i in version_labels]
    plt.xticks(x_ticks, x_labels, rotation=0, horizontalalignment='center')

    plt_all.legend(loc="upper right", fontsize='large')
    plt_all.savefig(f"{OUTPUT_FOLDER}/roofline_plot.out.pdf", bbox_inches='tight', pad_inches=0)
    plt_rho.savefig(f"{OUTPUT_FOLDER}/roofline_plot_rho.out.pdf", bbox_inches='tight', pad_inches=0)
    plt_feq.savefig(f"{OUTPUT_FOLDER}/roofline_plot_feq.out.pdf", bbox_inches='tight', pad_inches=0)
    plt_f.savefig(f"{OUTPUT_FOLDER}/roofline_plot_f.out.pdf", bbox_inches='tight', pad_inches=0)
    plt_vort.savefig(f"{OUTPUT_FOLDER}/roofline_plot_vort.out.pdf", bbox_inches='tight', pad_inches=0)
    bar_fig.savefig(f"{OUTPUT_FOLDER}/bar_plot.out.pdf", bbox_inches='tight', pad_inches=0)

    plt.show()

        
main()
