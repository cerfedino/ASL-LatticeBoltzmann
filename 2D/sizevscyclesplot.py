import os
import subprocess
import re
import numpy as np
from dotenv import load_dotenv, dotenv_values
import matplotlib.pyplot as plt
from matplotlib.colors import hsv_to_rgb
import sys
import argparse


PREVIOUS_VERSIONS_PATH = f"{os.path.dirname(__file__)}/previous_versions/"
OUTPUT_FOLDER = f"{os.path.dirname(__file__)}/artifacts/"


if len(sys.argv) <= 5:
    print("Usage: python3 sizevscyclesplot.py <Nx> <Ny> <Nt> < i0, i1 .. values")
    exit(1)

Nx = int(sys.argv[1])
Ny = int(sys.argv[2])
Nt = int(sys.argv[3])
i = [int(x) for x in sys.argv[4:]]


print(f"Nx: {Nx}, Ny: {Ny}, Nt: {Nt}, i: {i}")

# Read from environment variables, otherwise from .env file, otherwise display error message
try:
    print("Reading environment variables")
    L1_SIZE_BYTES = int(os.environ["L1_SIZE_BITS"])/8
    L2_SIZE_BYTES = int(os.environ["L2_SIZE_BITS"])/8
    L3_SIZE_BYTES = int(os.environ["L3_SIZE_BITS"])/8
except KeyError:
    print("Environment variables not found, reading from .env file")
    # make sure keys are there AND they hgave avalue set
    try:
        load_dotenv()
        L1_SIZE_BYTES = int(os.getenv("L1_SIZE_BITS"))/8
        L2_SIZE_BYTES = int(os.getenv("L2_SIZE_BITS"))/8
        L3_SIZE_BYTES = int(os.getenv("L3_SIZE_BITS"))/8
    except FileNotFoundError:
        print("Error: Environment variables not found and .env file not found")
        exit(1)
    except KeyError:
        print("Error: Environment variables not found and .env file does not contain the required keys")
        exit(1)

print("L1_SIZE_BYTES: ", L1_SIZE_BYTES)
print("L2_SIZE_BYTES: ", L2_SIZE_BYTES)
print("L3_SIZE_BYTES: ", L3_SIZE_BYTES)



def init_plot():

    PLT_FACECOLOR="#E4E4E4"
    PLT_BOUND_COLOR="#55575C"

    fig = plt.figure(figsize=(20, 12), facecolor="white")
    plt.gca().set_facecolor(PLT_FACECOLOR)

    plt.gca().yaxis.grid(True, which='major', color='w', linestyle='-', linewidth=1.5)
    plt.xlabel("Input size [bytes]", fontsize=14, rotation=0, labelpad=5, ha='center')
    plt.ylabel("Cycles [log scale]", fontsize=14, rotation=0, ha='left', labelpad=-15, position=(1,1))

    return fig


def run_executable_and_get_output(executable_path, args):
    try:
        print(f"Running {executable_path} {' '.join(args)}")
        result = subprocess.run([executable_path] + list(args), stdout=subprocess.PIPE, text=True)
        result.check_returncode() # Ensure we catch non-zero exit statuses
        return result.stdout
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Error running {executable_path}: {e}") from e

def main():
  XMAX=L1_SIZE_BYTES+L2_SIZE_BYTES+L3_SIZE_BYTES*1.2


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
  
  fig = init_plot()
  


  YMAX = 0
  for v_idx,path in enumerate(previous_versions):
    print("=========")
    
    conf = dotenv_values(path+"/.env")
    VERSION = conf["VERSION"]
    PATH = path
    TITLE = conf["TITLE"]
    DESC = conf["DESC"]
    COMMIT = conf["COMMIT"]
    WORK = conf["WORK"]
    DATA_MOVEMENT = conf["DATA_MOV"]

    print("\nVERSION: ", VERSION)
    print("Path: ", PATH)
    print("TITLE: ", TITLE)
    print("DESCRIPTION: ", DESC)
    print("COMMIT: ", COMMIT)
    print("WORK: ", WORK)
    print("DATA_MOV: ", DATA_MOVEMENT)
    
    print("")
    input_sizes  = []
    total_cycles = []
    for idx,multiplier in enumerate(i):
      if len(input_sizes) > 0 and input_sizes[-1] > XMAX:
        continue
      
      MNx = Nx * multiplier
      MNy = Ny * multiplier
      MNt = Nt
      
      print(f"[{idx+1}/{len(i)}] Running {VERSION} with {MNx}x{MNy} ...")
      
      memsize_bytes = MNy * MNx * 9 * 8 +  MNy * MNx * 9 * 8 + MNy * MNx * 8 + MNy * MNx * 8 + MNy * MNx * 8 + MNy * MNx * 8 + MNy * MNx * 8; + MNy * MNx * 8 + MNy * MNx * 4 + MNy * MNx * 4;
    
      # Execute executable file in subprocess and read stdout
      # If there is a specific executable for the current configuration, use it instead
      if os.path.exists(f"{path}/src/main_optimized_profile_{MNx}_{MNy}_{MNt}.o"):
        executable = f"{path}/src/main_optimized_profile_{MNx}_{MNy}_{MNt}.o"
      else:
        executable = f"{path}/src/main_optimized_profile.o"
          
      output = run_executable_and_get_output(executable, [str(MNx), str(MNy), str(MNt)])
      # print(output)
      
      try:
        runs = int(re.search(r"Run (\d+)/\d+ done\nProfiling results:", output).group(1))

        cycles_rho  = int(re.search(r".*Rho\s*Calculation: .* (\d+) cycles.*", output).group(1))
        cycles_feq  = int(re.search(r".*FEQ\s*Calculation: .* (\d+) cycles.*", output).group(1))
        cycles_f    = int(re.search(r".*F\s*Calculation: .* (\d+) cycles.*", output).group(1))
        cycles_vort = int(re.search(r".*Vort\s*Calculation: .* (\d+) cycles.*", output).group(1))
        cycles = cycles_rho + cycles_feq + cycles_f + cycles_vort
      except AttributeError as e:
        print("Couldn't match regex to output:\n")
        print(output)
        print("\n\nException: ", e)

      print(f"Rho  cycles: {cycles_rho }")
      print(f"FEQ  cycles: {cycles_feq }")
      print(f"F    cycles: {cycles_f   }")
      print(f"Vort cycles: {cycles_vort}")
      print(f"Total Number of Cycles: {cycles}")
      
      input_sizes.append(memsize_bytes)
      total_cycles.append(cycles)
    
    colorrange = (0, 90)
    curcolor = hsv_to_rgb([[(colorrange[0]+((colorrange[1]-colorrange[0])*((v_idx+1)/len(previous_versions))))/360,  1, 1]])[0]
    plt.gca().plot(input_sizes, total_cycles, '-ro', label=TITLE, color=curcolor)
    
    YMAX = max(YMAX, max(total_cycles))
    
  YMAX*=1.1

  plt.axvline(x=(L1_SIZE_BYTES), color='gray', linestyle='--', label=f"L1 Cache Size: {int(L1_SIZE_BYTES)} bytes")
  plt.text((L1_SIZE_BYTES), YMAX*0.99, f"L1", va='top', ha='right')

  plt.axvline(x=(L1_SIZE_BYTES+L2_SIZE_BYTES), color='gray', linestyle='--', label=f"L2 Cache Size: {int(L2_SIZE_BYTES)} bytes")
  plt.text((L1_SIZE_BYTES+L2_SIZE_BYTES), YMAX*0.99, f"L2", va='top', ha='right')

  plt.axvline(x=(L1_SIZE_BYTES+L2_SIZE_BYTES+L3_SIZE_BYTES), color='gray', linestyle='--', label=f"L3 Cache Size: {int(L3_SIZE_BYTES)} bytes")
  plt.text((L1_SIZE_BYTES+L2_SIZE_BYTES+L3_SIZE_BYTES), YMAX*0.99, f"L3", va='top', ha='right')

  plt.gca().get_yaxis().get_major_formatter().set_scientific(False)  # Disable scientific notation
  plt.yscale('log')
  plt.xlim(0, XMAX)
  plt.ylim(0, YMAX)
  fig.legend()
  fig.savefig(f"{OUTPUT_FOLDER}/sizevscycles.out.pdf")

  fig.show()

        
main()