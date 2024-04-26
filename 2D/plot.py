import matplotlib.pyplot as plt
import numpy as np
import os

output_folders = [x[0] for x in os.walk("output/")][1:]

# if the output folder is empty, remove it
for folder in output_folders:
    if len(os.listdir(folder)) == 0:
        output_folders.remove(folder)
        os.rmdir(folder)
    
output_folders.sort(reverse=True)
    
# TODO show may last x ones
print("Select a folder to plot:")
for i, folder in enumerate(output_folders):
    print(f"[{i}] {folder}")
    
# select a folder
folder = output_folders[int(input("Enter the number of the folder to plot: "))]
print(f"Selected folder: {folder}")

files = os.listdir(folder)
npy_files = [file for file in files if file.endswith(".npy")]

# take first file for now
file = npy_files[0]
print(f"Selected file: {file}")

data = np.load(f"{folder}/{file}")

print(f"Shape of data: {data.shape}")
print(f"Data: {data[:5]}")



# get all .npy that start with step_
npy_files = [file for file in files if file.startswith("step_")]

# sort the files
npy_files.sort()


cylinder = np.load("cylinder.npy") # todo store cylinder numpy array in the output folder


# plot the files
count = 0

for file in npy_files:
    
    print(f"Plotting file: {file}")
    
    # load the file
    vorticity = np.load(f"{folder}/{file}")
    
    plt.cla()

    vorticity = np.ma.array(vorticity, mask=cylinder)
    plt.imshow(vorticity, cmap='bwr')
    plt.imshow(~cylinder, cmap='gray', alpha=0.3)
    plt.clim(-.1, .1)
    ax = plt.gca()
    ax.invert_yaxis()
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)    
    ax.set_aspect('equal')    
    plt.pause(0.001)