import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib   
matplotlib.use("Qt5agg")

output_folders = [x[0] for x in os.walk("output/")][1:]

# if the output folder is empty, remove it
for folder in output_folders:
    if len(os.listdir(folder)) == 0:
        output_folders.remove(folder)
        os.rmdir(folder)
    
# remove folders that end in images or video
output_folders = [folder for folder in output_folders if not folder.endswith("images") and not folder.endswith("video")]

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
npy_files = [file for file in files if file.startswith("vorticity_")]

# sort the files
npy_files.sort()


cylinder = np.load(folder + "/cylinder.npy") # todo store cylinder numpy array in the output folder
print(f"Shape of cylinder: {cylinder.shape}")


#convert cylinder to int values
cylinder = cylinder.astype(bool)



Nx = cylinder.shape[1]
Ny = cylinder.shape[0]

X, Y = np.meshgrid(range(Nx), range(Ny))
cylinder_reference = (X - Nx/4)**2 + (Y - Ny/2)**2 < (Ny/4)**2


print(f"cylinder == cylinder_reference: {np.all(cylinder == cylinder_reference)}")



# write cyliner_reference to txt file
np.savetxt(f"{folder}/00_cylinder_reference.txt", cylinder_reference, fmt="%d")

# write cyliner to txt file
np.savetxt(f"{folder}/00_cylinder.txt", cylinder, fmt="%d")

# create images folder
if not os.path.exists(f"{folder}/images"):
    os.mkdir(f"{folder}/images")

# plot the files
count = 0

fig = plt.figure(figsize=(8,4), dpi=80)



for file in npy_files:
    print(f"Plotting file: {file}")
    # save the 50th file to txt
    if count == 50:
        np.savetxt(f"{folder}/50_vorticity.txt", np.load(f"{folder}/{file}"), fmt="%.6f")

    #print(f"Plotting file: {file}")
    
    # load the file
    vorticity = np.load(f"{folder}/{file}")
    #print(f"Shape of vorticity: {vorticity.shape}")

    
    plt.cla()
    #vorticity[cylinder] = np.nan
    vorticity = np.ma.array(vorticity, mask=cylinder)
    plt.imshow(vorticity, cmap='bwr')
    plt.imshow(~cylinder, cmap='gray', alpha=0.3)

    # put filename as title
    plt.title(file)

    plt.clim(-.1, .1)
    ax = plt.gca()
    ax.invert_yaxis()
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)    
    ax.set_aspect('equal')    
    #plt.pause(0.02)
    # save image in the same folder under /images
    plt.savefig(f"{folder}/images/{count}.png")
    count += 1


if not os.path.exists(f"{folder}/video"):
    os.mkdir(f"{folder}/video")

# combine the images into a video
# use codec H264
os.system(f"ffmpeg -r 30 -i {folder}/images/%01d.png -vcodec libx264 -y {folder}/video/output.mp4")